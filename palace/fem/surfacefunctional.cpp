// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfacefunctional.hpp"

#include <cstdlib>
#include <map>
#include <memory>
#include <mutex>
#include "fem/coefficient.hpp"
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "fem/libceed/basis.hpp"
#include "fem/libceed/coefficient.hpp"
#include "fem/libceed/functional.hpp"
#include "fem/libceed/integrator.hpp"
#include "fem/libceed/restriction.hpp"
#include "fem/mesh.hpp"
#include "linalg/vector.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/diagnostic.hpp"

PalacePragmaDiagnosticPush
PalacePragmaDiagnosticDisableUnused

#include "fem/qfunctions/32/surf_32_qf.h"
#include "fem/qfunctions/33/eval_33_qf.h"

PalacePragmaDiagnosticPop

namespace palace
{

namespace
{

// Key uniquely identifying the map(s) from boundary element reference coordinates to
// attached volume element reference coordinates (and thus the tabulated bases): the
// geometries involved, the orientation flag, and the quantized mapped quadrature point
// coordinates themselves. Keying directly on the mapped points is robust to the
// underlying mfem face-orientation encodings.
using FaceConfigKey = std::vector<long long>;

constexpr double QUANTIZE_SCALE = 1.0e10;

// Registry of mapped face integration rules. The IntegrationRule objects must have
// application lifetime: mfem::FiniteElement::GetDofToQuad caches tabulations keyed by
// the IntegrationRule pointer inside the (global, shared) FiniteElement objects, so
// destroying an IntegrationRule which was used for tabulation would leave a dangling
// cache entry.
const mfem::IntegrationRule *
GetRegisteredMappedIr(const FaceConfigKey &key,
                      const std::vector<mfem::IntegrationPoint> &pts)
{
  static std::map<FaceConfigKey, std::unique_ptr<mfem::IntegrationRule>> registry;
  static std::mutex registry_mutex;
  std::lock_guard<std::mutex> lock(registry_mutex);
  auto it = registry.find(key);
  if (it == registry.end())
  {
    auto ir = std::make_unique<mfem::IntegrationRule>(static_cast<int>(pts.size()));
    for (std::size_t q = 0; q < pts.size(); q++)
    {
      ir->IntPoint(static_cast<int>(q)) = pts[q];
    }
    it = registry.emplace(key, std::move(ir)).first;
  }
  return it->second.get();
}

// Evaluation plan for a single marked boundary element: which volume element(s) the
// field is evaluated from, and the face quadrature point positions mapped into the
// volume element reference space(s). elem_b >= 0 indicates two-sided evaluation.
struct ElemPlan
{
  int bdr_elem;
  int elem_a = -1, elem_b = -1;
  std::vector<mfem::IntegrationPoint> pts_a, pts_b;
  bool flip = false;  // Boundary element and face have the same orientation (o % 2 == 0),
                      // legacy coefficients invert the boundary element normal
};

// Data for one group of boundary elements sharing the same face configuration.
struct FaceGroup
{
  mfem::Geometry::Type bdr_geom;
  mfem::Geometry::Type vol_geom_a = mfem::Geometry::INVALID;
  mfem::Geometry::Type vol_geom_b = mfem::Geometry::INVALID;
  bool flip_normal = false;
  const mfem::IntegrationRule *mapped_ir_a = nullptr;  // Registry, application lifetime
  const mfem::IntegrationRule *mapped_ir_b = nullptr;
  std::vector<int> bdr_indices, vol_indices_a, vol_indices_b;
  std::vector<int> out_slots;  // Output vector slot for each boundary element
};

void AppendPoints(FaceConfigKey &key, const std::vector<mfem::IntegrationPoint> &pts)
{
  for (const auto &ip : pts)
  {
    key.push_back(std::llround(ip.x * QUANTIZE_SCALE));
    key.push_back(std::llround(ip.y * QUANTIZE_SCALE));
    key.push_back(std::llround(ip.z * QUANTIZE_SCALE));
  }
}

// Re-point the passive field inputs of each group operator at the given source vectors
// and accumulate into the output vector with CeedOperatorApplyAdd.
void ApplyAddGroups(const std::vector<fem::CeedGroupOperator> &groups,
                    const std::array<const Vector *, 4> &srcs, const Vector &out)
{
  for (const auto &[ceed, op, field_sources] : groups)
  {
    for (const auto &[name, source] : field_sources)
    {
      MFEM_ASSERT(srcs[source], "Missing source vector for libCEED field input!");
      CeedOperatorField field;
      CeedVector field_vec;
      PalaceCeedCall(ceed, CeedOperatorGetFieldByName(op, name.c_str(), &field));
      PalaceCeedCall(ceed, CeedOperatorFieldGetVector(field, &field_vec));
      ceed::InitCeedVector(*srcs[source], ceed, &field_vec, false);
    }
    CeedVector out_vec;
    ceed::InitCeedVector(out, ceed, &out_vec);
    PalaceCeedCall(
        ceed, CeedOperatorApplyAdd(op, CEED_VECTOR_NONE, out_vec, CEED_REQUEST_IMMEDIATE));
    PalaceCeedCall(ceed, CeedVectorDestroy(&out_vec));
  }
}

// Holds libCEED object references created during operator assembly for destruction once
// the assembled operator owns them.
struct CeedAssemblyScratch
{
  Ceed ceed;
  std::vector<CeedVector> vecs;
  std::vector<CeedElemRestriction> restrs;
  std::vector<CeedBasis> bases;

  CeedAssemblyScratch(Ceed ceed) : ceed(ceed) {}
  ~CeedAssemblyScratch()
  {
    for (auto &v : vecs)
    {
      PalaceCeedCall(ceed, CeedVectorDestroy(&v));
    }
    for (auto &r : restrs)
    {
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&r));
    }
    for (auto &b : bases)
    {
      PalaceCeedCall(ceed, CeedBasisDestroy(&b));
    }
  }
};

}  // namespace

SurfaceFunctional::SurfaceFunctional(Kind kind, const Mesh &mesh,
                                     const mfem::Array<int> &bdr_attr_marker,
                                     const mfem::ParFiniteElementSpace *fespace)
  : kind(kind), fespace_e(fespace), fespace_b(nullptr), mat_op(nullptr),
    comm(mesh.GetComm())
{
  MFEM_VERIFY(kind == Kind::AREA || kind == Kind::HCURL_NORM2,
              "Invalid SurfaceFunctional constructor for the requested functional kind!");
  MFEM_VERIFY(kind == Kind::AREA || fespace,
              "SurfaceFunctional requires a field finite element space for functionals "
              "with field inputs!");
  Assemble(mesh, bdr_attr_marker);
}

SurfaceFunctional::SurfaceFunctional(Kind kind, const Mesh &mesh,
                                     const mfem::Array<int> &bdr_attr_marker,
                                     const mfem::ParFiniteElementSpace &fespace, int lod)
  : kind(kind), fespace_e(kind == Kind::BDR_FIELD_E ? &fespace : nullptr),
    fespace_b(kind == Kind::BDR_FIELD_B ? &fespace : nullptr), mat_op(nullptr),
    comm(mesh.GetComm()), viz_lod(lod)
{
  MFEM_VERIFY(kind == Kind::BDR_FIELD_E || kind == Kind::BDR_FIELD_B,
              "Invalid SurfaceFunctional constructor for the requested functional kind!");
  Assemble(mesh, bdr_attr_marker);
}

SurfaceFunctional::SurfaceFunctional(Kind kind, const Mesh &mesh,
                                     const mfem::Array<int> &bdr_attr_marker,
                                     const mfem::ParFiniteElementSpace &fespace,
                                     const MaterialOperator &mat_op, int lod,
                                     double scaling)
  : kind(kind),
    fespace_e((kind == Kind::BDR_FLUX_Q || kind == Kind::BDR_ENERGY_E) ? &fespace
                                                                       : nullptr),
    fespace_b((kind == Kind::BDR_CURRENT_J || kind == Kind::BDR_ENERGY_M) ? &fespace
                                                                          : nullptr),
    mat_op(&mat_op), comm(mesh.GetComm()), viz_lod(lod), viz_scaling(scaling)
{
  MFEM_VERIFY(kind == Kind::BDR_FLUX_Q || kind == Kind::BDR_CURRENT_J ||
                  kind == Kind::BDR_ENERGY_E || kind == Kind::BDR_ENERGY_M,
              "Invalid SurfaceFunctional constructor for the requested functional kind!");
  Assemble(mesh, bdr_attr_marker);
}

SurfaceFunctional::SurfaceFunctional(const Mesh &mesh,
                                     const mfem::Array<int> &bdr_attr_marker,
                                     const mfem::ParFiniteElementSpace &nd_fespace,
                                     const MaterialOperator &mat_op,
                                     InterfaceDielectric type, double t_i, double epsilon_i)
  : kind(Kind::INTERFACE_EPR), epr_type(type), epr_t(t_i), epr_epsilon(epsilon_i),
    fespace_e(&nd_fespace), fespace_b(nullptr), mat_op(&mat_op), comm(mesh.GetComm())
{
  Assemble(mesh, bdr_attr_marker);
}

SurfaceFunctional::SurfaceFunctional(const Mesh &mesh,
                                     const mfem::Array<int> &bdr_attr_marker,
                                     const mfem::ParFiniteElementSpace *nd_fespace,
                                     const mfem::ParFiniteElementSpace *rt_fespace,
                                     const MaterialOperator &mat_op, SurfaceFlux type,
                                     bool two_sided, const mfem::Vector &x0)
  : kind(Kind::SURFACE_FLUX), flux_type(type), flux_two_sided(two_sided), flux_x0(x0),
    fespace_e(nd_fespace), fespace_b(rt_fespace), mat_op(&mat_op), comm(mesh.GetComm())
{
  MFEM_VERIFY(
      (nd_fespace || (type != SurfaceFlux::ELECTRIC && type != SurfaceFlux::POWER)) &&
          (rt_fespace || (type != SurfaceFlux::MAGNETIC && type != SurfaceFlux::POWER)),
      "Missing finite element space for surface flux functional!");
  Assemble(mesh, bdr_attr_marker);
}

SurfaceFunctional::SurfaceFunctional(const Mesh &mesh,
                                     const mfem::Array<int> &bdr_attr_marker,
                                     const mfem::ParFiniteElementSpace &nd_fespace,
                                     const mfem::ParFiniteElementSpace &rt_fespace,
                                     const MaterialOperator &mat_op,
                                     const std::vector<std::array<double, 3>> &r_naughts)
  : kind(Kind::FARFIELD), farfield_dirs(r_naughts), farfield_mesh(&mesh),
    fespace_e(&nd_fespace), fespace_b(&rt_fespace), mat_op(&mat_op), comm(mesh.GetComm())
{
  farfield_marker = bdr_attr_marker;
  Assemble(mesh, bdr_attr_marker);
}

SurfaceFunctional::~SurfaceFunctional()
{
  for (auto &group : groups)
  {
    PalaceCeedCall(group.ceed, CeedOperatorDestroy(&group.op));
  }
}

bool SurfaceFunctional::Enabled()
{
  static const bool enabled = !std::getenv("PALACE_LEGACY_SURFACE_POSTPRO");
  return enabled;
}

void SurfaceFunctional::Assemble(const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker)
{
  AssembleLocal(mesh, bdr_attr_marker);

  // The validity decision must be globally consistent: a rank may locally require an
  // unsupported configuration (e.g. two-sided evaluation across a process boundary)
  // while others do not. All ranks must agree on whether the libCEED path or the
  // legacy fallback is used so that the collective evaluation calls match.
  bool global_valid = valid;
  Mpi::GlobalAnd(1, &global_valid, comm);
  if (!global_valid && valid)
  {
    // Discard the locally assembled operators; the legacy path will be used.
    for (auto &group : groups)
    {
      PalaceCeedCall(group.ceed, CeedOperatorDestroy(&group.op));
    }
    groups.clear();
    valid = false;
  }
}

void SurfaceFunctional::AssembleLocal(const Mesh &mesh,
                                      const mfem::Array<int> &bdr_attr_marker)
{
  if (mesh.Dimension() != 3 || mesh.SpaceDimension() != 3)
  {
    // Not yet supported (2D solver meshes): callers fall back to the legacy paths.
    valid = false;
    return;
  }
  const mfem::ParMesh &pmesh = mesh.Get();
  MFEM_VERIFY(pmesh.GetNodes(), "The mesh has no nodal FE space!");
  const mfem::FiniteElementSpace &mesh_fespace = *pmesh.GetNodes()->FESpace();
  const bool need_field = (kind != Kind::AREA);

  // Plan the evaluation for each marked boundary element and group the elements by
  // their face configuration. All elements in a group share the tabulated bases and are
  // processed by a single libCEED operator.
  std::map<FaceConfigKey, FaceGroup> face_groups;
  int num_marked = 0;
  {
    constexpr double threshold = 1.0 - 1.0e-6;
    mfem::FaceElementTransformations FET;
    mfem::IsoparametricTransformation T1, T2;
    for (int i = 0; i < pmesh.GetNBE(); i++)
    {
      const int attr = pmesh.GetBdrAttribute(i);
      if (!bdr_attr_marker[attr - 1])
      {
        continue;
      }

      // Get the face and attached element transformations following the same
      // conventions as the legacy BdrGridFunctionCoefficient evaluation path
      // (FET.Elem1 always exists and is local, FET.Elem2 exists for interior
      // boundaries and may correspond to a face neighbor element on another process
      // for boundaries on parallel interfaces).
      const bool flip = BdrGridFunctionCoefficient::GetBdrElementNeighborTransformations(
          i, pmesh, FET, T1, T2);
      const bool has_elem2 = (FET.Elem2 != nullptr);
      const bool elem2_local = has_elem2 && (FET.Elem2No < pmesh.GetNE());

      // Decide which side(s) the field is evaluated from following the conventions of
      // the legacy coefficients (see InterfaceDielectricCoefficient and
      // BdrSurfaceFluxCoefficient).
      ElemPlan plan;
      plan.bdr_elem = i;
      plan.flip = flip;
      if (!need_field)
      {
        // No field inputs (side selection does not apply).
      }
      else if (kind == Kind::HCURL_NORM2)
      {
        plan.elem_a = FET.Elem1No;
      }
      else if (kind == Kind::FARFIELD)
      {
        if (has_elem2)
        {
          // Far-field computations are only supported on external boundaries (the
          // legacy path errors in this case as well).
          valid = false;
          return;
        }
        plan.elem_a = FET.Elem1No;
      }
      else if (kind == Kind::SURFACE_FLUX || IsBufferKind(kind) ||
               (kind == Kind::INTERFACE_EPR && epr_type == InterfaceDielectric::DEFAULT))
      {
        plan.elem_a = FET.Elem1No;
        if (has_elem2)
        {
          plan.elem_b = FET.Elem2No;
        }
      }
      else  // INTERFACE_EPR with MA, MS, or SA
      {
        // Single-sided evaluation on the vacuum (MA, SA) or substrate (MS) side, with
        // averaging if both sides qualify, skipping the element if neither does.
        const bool vacuum_side = (epr_type != InterfaceDielectric::MS);
        auto OnSide = [&](int elem_attr)
        {
          const double ls = mat_op->GetLightSpeedMax(elem_attr);
          return vacuum_side ? (ls >= threshold) : (ls < threshold);
        };
        const bool use_elem1 = OnSide(FET.Elem1->Attribute);
        const bool use_elem2 = has_elem2 && OnSide(FET.Elem2->Attribute);
        if (use_elem1)
        {
          plan.elem_a = FET.Elem1No;
          if (use_elem2)
          {
            plan.elem_b = FET.Elem2No;
          }
        }
        else if (use_elem2)
        {
          plan.elem_a = FET.Elem2No;
        }
        else
        {
          // Neither side qualifies: the element contributes zero, skip it.
          continue;
        }
      }

      // TODO: Support face neighbor (other process) elements for two-sided interior
      // boundaries on parallel interfaces, as handled by the legacy coefficients via
      // ParMesh::ExchangeFaceNbrData. Requires element restrictions indexing into the
      // face neighbor data vector. Until then, fall back to the legacy paths.
      if ((plan.elem_a >= 0 && plan.elem_a >= pmesh.GetNE()) ||
          (plan.elem_b >= 0 && !elem2_local))
      {
        valid = false;
        return;
      }

      // Map the face quadrature points to the volume element reference space(s):
      // boundary element reference coordinates -> face reference coordinates
      // (TransformBdrElementToFace with the boundary element to face orientation) ->
      // element 1/2 reference coordinates (FET.Loc1/Loc2).
      const auto bdr_geom = pmesh.GetBdrElementGeometry(i);
      const bool buffer_kind = IsBufferKind(kind);
      const mfem::IntegrationRule &face_ir =
          buffer_kind ? mfem::GlobGeometryRefiner.Refine(bdr_geom, viz_lod, 1)->RefPts
                      : mfem::IntRules.Get(
                            bdr_geom, fem::DefaultIntegrationOrder::Get(pmesh, bdr_geom));
      const int nq = face_ir.GetNPoints();
      int f, o;
      pmesh.GetBdrElementFace(i, &f, &o);
      auto MapPoints = [&](mfem::IntegrationPointTransformation &loc,
                           std::vector<mfem::IntegrationPoint> &pts)
      {
        pts.resize(nq);
        for (int q = 0; q < nq; q++)
        {
          mfem::IntegrationPoint fip = mfem::Mesh::TransformBdrElementToFace(
              FET.GetGeometryType(), o, face_ir.IntPoint(q));
          loc.Transform(fip, pts[q]);
          pts[q].weight = face_ir.IntPoint(q).weight;
        }
      };
      if (plan.elem_a >= 0)
      {
        MapPoints(plan.elem_a == FET.Elem1No ? FET.Loc1 : FET.Loc2, plan.pts_a);
      }
      if (plan.elem_b >= 0)
      {
        MapPoints(FET.Loc2, plan.pts_b);
      }

      // Build the group key from the geometries, orientation flag, and quantized
      // mapped quadrature point coordinates.
      const auto vol_geom_a = (plan.elem_a >= 0) ? pmesh.GetElementGeometry(plan.elem_a)
                                                 : mfem::Geometry::INVALID;
      const auto vol_geom_b = (plan.elem_b >= 0) ? pmesh.GetElementGeometry(plan.elem_b)
                                                 : mfem::Geometry::INVALID;
      FaceConfigKey key;
      key.reserve(5 + 3 * (plan.pts_a.size() + plan.pts_b.size()));
      key.push_back(static_cast<long long>(bdr_geom));
      key.push_back(static_cast<long long>(vol_geom_a));
      key.push_back(static_cast<long long>(vol_geom_b));
      key.push_back(static_cast<long long>(plan.flip));
      key.push_back(static_cast<long long>(nq));
      AppendPoints(key, plan.pts_a);
      AppendPoints(key, plan.pts_b);

      auto it = face_groups.find(key);
      if (it == face_groups.end())
      {
        FaceGroup group;
        group.bdr_geom = bdr_geom;
        group.vol_geom_a = vol_geom_a;
        group.vol_geom_b = vol_geom_b;
        group.flip_normal = plan.flip;
        if (plan.elem_a >= 0)
        {
          FaceConfigKey key_a = key;
          key_a.push_back(0);  // Side tag
          group.mapped_ir_a = GetRegisteredMappedIr(key_a, plan.pts_a);
        }
        if (plan.elem_b >= 0)
        {
          FaceConfigKey key_b = key;
          key_b.push_back(1);  // Side tag
          group.mapped_ir_b = GetRegisteredMappedIr(key_b, plan.pts_b);
        }
        it = face_groups.emplace(key, std::move(group)).first;
      }
      it->second.bdr_indices.push_back(i);
      if (plan.elem_a >= 0)
      {
        it->second.vol_indices_a.push_back(plan.elem_a);
      }
      if (plan.elem_b >= 0)
      {
        it->second.vol_indices_b.push_back(plan.elem_b);
      }
      if (IsBufferKind(kind))
      {
        if (buffer_bases.empty())
        {
          buffer_bases.resize(pmesh.GetNBE(), -1);
        }
        buffer_bases[i] = buffer_size;
        it->second.out_slots.push_back(buffer_size);
        buffer_size += nq * BufferNumComp(kind);
        num_marked++;
      }
      else
      {
        it->second.out_slots.push_back(num_marked++);
      }
    }
  }

  // Initialize the local output vector and field staging vector. Far-field operators
  // produce 6 values (Re/Im of a 3-vector) per direction per element.
  const int num_out =
      (kind == Kind::FARFIELD) ? 6 * static_cast<int>(farfield_dirs.size()) : 1;
  local_out.SetSize(
      (kind == Kind::BDR_FIELD_E || kind == Kind::BDR_FIELD_B) ? 0 : num_marked * num_out);
  local_out.UseDevice(true);
  if (need_field)
  {
    const int max_vsize = std::max(fespace_e ? fespace_e->GetVSize() : 0,
                                   fespace_b ? fespace_b->GetVSize() : 0);
    field_staging.SetSize(max_vsize);
    field_staging.UseDevice(true);
    field_staging = 0.0;
  }

  // Build the (group independent part of the) QFunction context for the integrand.
  std::vector<CeedIntScalar> base_ctx;
  if (kind == Kind::INTERFACE_EPR)
  {
    base_ctx.resize(2);
    base_ctx[0].second = 0.0;
    base_ctx[1].second = 0.0;
    switch (epr_type)
    {
      case InterfaceDielectric::DEFAULT:
        base_ctx[0].second = 0.5 * epr_t * epr_epsilon;
        break;
      case InterfaceDielectric::MA:
        base_ctx[0].second = 0.5 * epr_t / epr_epsilon;
        break;
      case InterfaceDielectric::MS:
        {
          base_ctx[0].second = 0.5 * epr_t / epr_epsilon;
          MaterialPropertyCoefficient epsilon_func(mat_op->GetAttributeToMaterial(),
                                                   mat_op->GetPermittivityReal());
          auto mat_ctx = ceed::PopulateCoefficientContext(3, &epsilon_func);
          base_ctx.insert(base_ctx.end(), mat_ctx.begin(), mat_ctx.end());
        }
        break;
      case InterfaceDielectric::SA:
        base_ctx[0].second = 0.5 * epr_t * epr_epsilon;
        base_ctx[1].second = 0.5 * epr_t / epr_epsilon;
        break;
    }
  }
  else if (kind == Kind::SURFACE_FLUX)
  {
    base_ctx.resize(5);
    base_ctx[0].second = 1.0;  // Normal sign, set per group
    base_ctx[1].first = flux_two_sided;
    for (int d = 0; d < 3; d++)
    {
      base_ctx[2 + d].second = (flux_x0.Size() > d) ? flux_x0(d) : 0.0;
    }
    if (flux_type == SurfaceFlux::ELECTRIC)
    {
      MaterialPropertyCoefficient epsilon_func(mat_op->GetAttributeToMaterial(),
                                               mat_op->GetPermittivityReal());
      auto mat_ctx = ceed::PopulateCoefficientContext(3, &epsilon_func);
      base_ctx.insert(base_ctx.end(), mat_ctx.begin(), mat_ctx.end());
    }
    else if (flux_type == SurfaceFlux::POWER)
    {
      MaterialPropertyCoefficient invmu_func(mat_op->GetAttributeToMaterial(),
                                             mat_op->GetInvPermeability());
      auto mat_ctx = ceed::PopulateCoefficientContext(3, &invmu_func);
      base_ctx.insert(base_ctx.end(), mat_ctx.begin(), mat_ctx.end());
    }
  }
  else if (kind == Kind::BDR_FLUX_Q || kind == Kind::BDR_CURRENT_J ||
           kind == Kind::BDR_ENERGY_E || kind == Kind::BDR_ENERGY_M)
  {
    base_ctx.resize(2);
    base_ctx[0].second = 1.0;  // Normal sign, set per group
    base_ctx[1].second = viz_scaling;
    MaterialPropertyCoefficient coeff_func(
        mat_op->GetAttributeToMaterial(),
        (kind == Kind::BDR_FLUX_Q || kind == Kind::BDR_ENERGY_E)
            ? mat_op->GetPermittivityReal()
            : mat_op->GetInvPermeability());
    auto mat_ctx = ceed::PopulateCoefficientContext(3, &coeff_func);
    base_ctx.insert(base_ctx.end(), mat_ctx.begin(), mat_ctx.end());
  }
  else if (kind == Kind::FARFIELD)
  {
    const int N = static_cast<int>(farfield_dirs.size());
    base_ctx.resize(4 + 3 * N);
    base_ctx[0].second = 1.0;  // Normal sign, set per group
    base_ctx[1].second = farfield_omega_re;
    base_ctx[2].second = farfield_omega_im;
    base_ctx[3].first = N;
    for (int d = 0; d < N; d++)
    {
      for (int c = 0; c < 3; c++)
      {
        base_ctx[4 + 3 * d + c].second = farfield_dirs[d][c];
      }
    }
    MaterialPropertyCoefficient c0_func(mat_op->GetAttributeToMaterial(),
                                        mat_op->GetLightSpeed());
    auto mat_ctx = ceed::PopulateCoefficientContext(3, &c0_func);
    base_ctx.insert(base_ctx.end(), mat_ctx.begin(), mat_ctx.end());
  }
  else
  {
    base_ctx.resize(2);
    base_ctx[0].second = 0.0;
    base_ctx[1].second = 0.0;
  }

  // Assemble a libCEED operator for each group. For now, all operators are constructed
  // on a single Ceed context (no OpenMP parallel assembly or application; correctness
  // first, this can be extended with the thread partitioning of fem/mesh.cpp later).
  Ceed ceed = ceed::internal::GetCeedObjects()[0];
  for (auto &[key, group] : face_groups)
  {
    const std::size_t num_elem = group.bdr_indices.size();
    const bool has_b = !group.vol_indices_b.empty();
    const bool buffer_kind = IsBufferKind(kind);
    const mfem::IntegrationRule &face_ir =
        buffer_kind
            ? mfem::GlobGeometryRefiner.Refine(group.bdr_geom, viz_lod, 1)->RefPts
            : mfem::IntRules.Get(group.bdr_geom,
                                 fem::DefaultIntegrationOrder::Get(pmesh, group.bdr_geom));

    // Objects are owned by the assembled operator; scratch destroys our references.
    CeedAssemblyScratch scratch(ceed);

    // Assemble the inputs in the order expected by the QFunctions: quadrature weights
    // and boundary element mesh node gradients (surface measure and normal), per-side
    // volume element attributes (material lookup) and mesh node gradients at the mapped
    // points (Piola transformations), coordinates (SURFACE_FLUX only), then the field
    // inputs (side 1 fields, then side 2 fields). The geometry is computed on the fly
    // in the QFunctions (no stored geometry factor data).
    std::vector<ceed::CeedFunctionalFieldInput> inputs;
    std::vector<std::pair<std::string, int>> field_sources;
    const mfem::FiniteElement *face_mesh_fe =
        mesh_fespace.FEColl()->FiniteElementForGeometry(group.bdr_geom);
    if (!face_mesh_fe)
    {
      face_mesh_fe = mesh_fespace.FEColl()->TraceFiniteElementForGeometry(group.bdr_geom);
    }
    const bool field_kind = (kind == Kind::BDR_FIELD_E || kind == Kind::BDR_FIELD_B);
    if (!field_kind)
    {
      CeedBasis face_mesh_basis;
      ceed::InitBasisAtPoints(*face_mesh_fe, face_ir, mesh_fespace.GetVDim(), ceed,
                              &face_mesh_basis);
      CeedElemRestriction face_mesh_restr = FiniteElementSpace::BuildCeedElemRestriction(
          mesh_fespace, ceed, group.bdr_geom, group.bdr_indices, /*is_interp*/ true);
      CeedVector mesh_nodes_vec;
      ceed::InitCeedVector(*mesh_fespace.GetMesh()->GetNodes(), ceed, &mesh_nodes_vec);
      if (!buffer_kind)
      {
        inputs.push_back({"qw", nullptr, nullptr, face_mesh_basis, ceed::EvalMode::Weight});
      }
      inputs.push_back(
          {"x_f", mesh_nodes_vec, face_mesh_restr, face_mesh_basis, ceed::EvalMode::Grad});
      scratch.vecs.push_back(mesh_nodes_vec);
      scratch.restrs.push_back(face_mesh_restr);
      scratch.bases.push_back(face_mesh_basis);
    }
    auto AddVolGeomInputs = [&](const std::string &suffix, const std::vector<int> &indices,
                                mfem::Geometry::Type geom, const mfem::IntegrationRule &ir)
    {
      const int num_pts = ir.GetNPoints();
      auto &elem_attr = elem_attrs.emplace_back(indices.size());
      const auto &loc_attr = mesh.GetCeedAttributes();
      for (std::size_t k = 0; k < indices.size(); k++)
      {
        elem_attr[k] = loc_attr.at(pmesh.GetAttribute(indices[k]));
      }
      CeedElemRestriction attr_restr;
      CeedBasis attr_basis;
      CeedVector attr_vec;
      PalaceCeedCall(
          ceed, CeedElemRestrictionCreateStrided(ceed, indices.size(), 1, 1, indices.size(),
                                                 CEED_STRIDES_BACKEND, &attr_restr));
      {
        // Note: ceed::GetCeedTopology(CEED_TOPOLOGY_LINE) == 1.
        mfem::Vector Bt(num_pts), Gt(num_pts), qX(num_pts), qW(num_pts);
        Bt = 1.0;
        Gt = 0.0;
        qX = 0.0;
        qW = 0.0;
        PalaceCeedCall(ceed, CeedBasisCreateH1(ceed, CEED_TOPOLOGY_LINE, 1, 1, num_pts,
                                               Bt.GetData(), Gt.GetData(), qX.GetData(),
                                               qW.GetData(), &attr_basis));
      }
      ceed::InitCeedVector(elem_attr, ceed, &attr_vec);
      inputs.push_back(
          {"attr_" + suffix, attr_vec, attr_restr, attr_basis, ceed::EvalMode::Interp});
      CeedElemRestriction mesh_restr = FiniteElementSpace::BuildCeedElemRestriction(
          mesh_fespace, ceed, geom, indices, /*is_interp*/ true);
      const mfem::FiniteElement *mesh_fe =
          mesh_fespace.FEColl()->FiniteElementForGeometry(geom);
      CeedBasis mesh_basis;
      ceed::InitBasisAtPoints(*mesh_fe, ir, mesh_fespace.GetVDim(), ceed, &mesh_basis);
      CeedVector mesh_nodes_vec;
      ceed::InitCeedVector(*mesh_fespace.GetMesh()->GetNodes(), ceed, &mesh_nodes_vec);
      inputs.push_back(
          {"x_" + suffix, mesh_nodes_vec, mesh_restr, mesh_basis, ceed::EvalMode::Grad});
      scratch.vecs.push_back(attr_vec);
      scratch.vecs.push_back(mesh_nodes_vec);
      scratch.restrs.push_back(attr_restr);
      scratch.restrs.push_back(mesh_restr);
      scratch.bases.push_back(attr_basis);
      scratch.bases.push_back(mesh_basis);
    };
    if (!group.vol_indices_a.empty())
    {
      AddVolGeomInputs("1", group.vol_indices_a, group.vol_geom_a, *group.mapped_ir_a);
    }
    if (has_b)
    {
      AddVolGeomInputs("2", group.vol_indices_b, group.vol_geom_b, *group.mapped_ir_b);
    }
    if (kind == Kind::SURFACE_FLUX || kind == Kind::FARFIELD)
    {
      // Coordinates at the face quadrature points, interpolated from the mesh nodes
      // with the boundary element basis (constant input, never re-pointed).
      CeedElemRestriction x_restr = FiniteElementSpace::BuildCeedElemRestriction(
          mesh_fespace, ceed, group.bdr_geom, group.bdr_indices, /*is_interp*/ true);
      CeedBasis x_basis;
      ceed::InitBasisAtPoints(*face_mesh_fe, face_ir, mesh_fespace.GetVDim(), ceed,
                              &x_basis);
      CeedVector x_vec;
      ceed::InitCeedVector(*mesh_fespace.GetMesh()->GetNodes(), ceed, &x_vec);
      inputs.push_back({"x", x_vec, x_restr, x_basis, ceed::EvalMode::Interp});
      scratch.vecs.push_back(x_vec);
      scratch.restrs.push_back(x_restr);
      scratch.bases.push_back(x_basis);
    }
    auto AddFieldInput = [&](const std::string &name, int source,
                             const mfem::ParFiniteElementSpace &fespace,
                             const std::vector<int> &indices, mfem::Geometry::Type geom,
                             const mfem::IntegrationRule &ir)
    {
      CeedElemRestriction restr;
      CeedBasis basis;
      CeedVector vec;
      ceed::InitRestriction(fespace, indices, false, /*is_interp*/ true, false, ceed,
                            &restr);
      const mfem::FiniteElement *fe = fespace.FEColl()->FiniteElementForGeometry(geom);
      MFEM_VERIFY(fe, "Unable to get field finite element for surface functional!");
      ceed::InitBasisAtPoints(*fe, ir, fespace.GetVDim(), ceed, &basis);
      ceed::InitCeedVector(field_staging, ceed, &vec);
      inputs.push_back({name, vec, restr, basis, ceed::EvalMode::Interp});
      field_sources.emplace_back(name, source);
      scratch.vecs.push_back(vec);
      scratch.restrs.push_back(restr);
      scratch.bases.push_back(basis);
    };
    if (kind == Kind::HCURL_NORM2 || kind == Kind::INTERFACE_EPR || buffer_kind)
    {
      const auto &field_fespace = fespace_b ? *fespace_b : *fespace_e;
      AddFieldInput("u_1", 0, field_fespace, group.vol_indices_a, group.vol_geom_a,
                    *group.mapped_ir_a);
      if (has_b)
      {
        AddFieldInput("u_2", 0, field_fespace, group.vol_indices_b, group.vol_geom_b,
                      *group.mapped_ir_b);
      }
    }
    else if (kind == Kind::FARFIELD)
    {
      AddFieldInput("u_1", 0, *fespace_e, group.vol_indices_a, group.vol_geom_a,
                    *group.mapped_ir_a);
      AddFieldInput("u_2", 1, *fespace_e, group.vol_indices_a, group.vol_geom_a,
                    *group.mapped_ir_a);
      AddFieldInput("u_3", 2, *fespace_b, group.vol_indices_a, group.vol_geom_a,
                    *group.mapped_ir_a);
      AddFieldInput("u_4", 3, *fespace_b, group.vol_indices_a, group.vol_geom_a,
                    *group.mapped_ir_a);
    }
    else if (kind == Kind::SURFACE_FLUX)
    {
      int count = 0;
      auto AddSide = [&](const std::vector<int> &indices, mfem::Geometry::Type geom,
                         const mfem::IntegrationRule &ir)
      {
        if (flux_type == SurfaceFlux::ELECTRIC || flux_type == SurfaceFlux::POWER)
        {
          AddFieldInput("u_" + std::to_string(++count), 0, *fespace_e, indices, geom, ir);
        }
        if (flux_type == SurfaceFlux::MAGNETIC || flux_type == SurfaceFlux::POWER)
        {
          AddFieldInput("u_" + std::to_string(++count), 1, *fespace_b, indices, geom, ir);
        }
      };
      AddSide(group.vol_indices_a, group.vol_geom_a, *group.mapped_ir_a);
      if (has_b)
      {
        AddSide(group.vol_indices_b, group.vol_geom_b, *group.mapped_ir_b);
      }
    }

    // Output restriction: for integral kinds, num_out slots per boundary element in
    // the local output vector (component stride num_marked); for the boundary
    // visualization field kinds, 3 components per lattice point scattering into the
    // output buffer at the per-element base offsets.
    CeedElemRestriction out_restr;
    if (buffer_kind)
    {
      const int nq = face_ir.GetNPoints();
      const int nc = BufferNumComp(kind);
      std::vector<CeedInt> offsets(num_elem * nq);
      for (std::size_t e = 0; e < num_elem; e++)
      {
        for (int j = 0; j < nq; j++)
        {
          offsets[e * nq + j] = group.out_slots[e] + nc * j;
        }
      }
      PalaceCeedCall(ceed, CeedElemRestrictionCreate(ceed, static_cast<CeedInt>(num_elem),
                                                     nq, nc, 1, (CeedSize)buffer_size,
                                                     CEED_MEM_HOST, CEED_COPY_VALUES,
                                                     offsets.data(), &out_restr));
    }
    else
    {
      PalaceCeedCall(ceed, CeedElemRestrictionCreate(
                               ceed, static_cast<CeedInt>(num_elem), 1, num_out, num_marked,
                               (CeedSize)num_marked * num_out, CEED_MEM_HOST,
                               CEED_COPY_VALUES, group.out_slots.data(), &out_restr));
    }
    scratch.restrs.push_back(out_restr);

    // Select the QFunction and finalize the (group dependent) context.
    std::vector<CeedIntScalar> ctx = base_ctx;
    ceed::CeedQFunctionInfo info;
    switch (kind)
    {
      case Kind::AREA:
        info.apply_qf = f_integ_surf_area_32;
        info.apply_qf_path = PalaceQFunctionRelativePath(f_integ_surf_area_32_loc);
        break;
      case Kind::HCURL_NORM2:
        info.apply_qf = f_integ_surf_hcurl_norm2_32;
        info.apply_qf_path = PalaceQFunctionRelativePath(f_integ_surf_hcurl_norm2_32_loc);
        break;
      case Kind::INTERFACE_EPR:
        switch (epr_type)
        {
          case InterfaceDielectric::DEFAULT:
            info.apply_qf = has_b ? f_integ_surf_epr_def_2_32 : f_integ_surf_epr_def_1_32;
            info.apply_qf_path = PalaceQFunctionRelativePath(
                has_b ? f_integ_surf_epr_def_2_32_loc : f_integ_surf_epr_def_1_32_loc);
            break;
          case InterfaceDielectric::MA:
            info.apply_qf = has_b ? f_integ_surf_epr_ma_2_32 : f_integ_surf_epr_ma_1_32;
            info.apply_qf_path = PalaceQFunctionRelativePath(
                has_b ? f_integ_surf_epr_ma_2_32_loc : f_integ_surf_epr_ma_1_32_loc);
            break;
          case InterfaceDielectric::MS:
            info.apply_qf = has_b ? f_integ_surf_epr_ms_2_32 : f_integ_surf_epr_ms_1_32;
            info.apply_qf_path = PalaceQFunctionRelativePath(
                has_b ? f_integ_surf_epr_ms_2_32_loc : f_integ_surf_epr_ms_1_32_loc);
            break;
          case InterfaceDielectric::SA:
            info.apply_qf = has_b ? f_integ_surf_epr_sa_2_32 : f_integ_surf_epr_sa_1_32;
            info.apply_qf_path = PalaceQFunctionRelativePath(
                has_b ? f_integ_surf_epr_sa_2_32_loc : f_integ_surf_epr_sa_1_32_loc);
            break;
        }
        break;
      case Kind::BDR_FIELD_E:
        info.apply_qf = has_b ? f_eval_bdr_hcurl_2_32 : f_eval_bdr_hcurl_1_32;
        info.apply_qf_path = PalaceQFunctionRelativePath(has_b ? f_eval_bdr_hcurl_2_32_loc
                                                               : f_eval_bdr_hcurl_1_32_loc);
        break;
      case Kind::BDR_FIELD_B:
        info.apply_qf = has_b ? f_eval_bdr_hdiv_2_32 : f_eval_bdr_hdiv_1_32;
        info.apply_qf_path = PalaceQFunctionRelativePath(has_b ? f_eval_bdr_hdiv_2_32_loc
                                                               : f_eval_bdr_hdiv_1_32_loc);
        break;
      case Kind::BDR_FLUX_Q:
        ctx[0].second = group.flip_normal ? -1.0 : 1.0;
        info.apply_qf = has_b ? f_eval_bdr_flux_q_2_32 : f_eval_bdr_flux_q_1_32;
        info.apply_qf_path = PalaceQFunctionRelativePath(
            has_b ? f_eval_bdr_flux_q_2_32_loc : f_eval_bdr_flux_q_1_32_loc);
        break;
      case Kind::BDR_CURRENT_J:
        ctx[0].second = group.flip_normal ? -1.0 : 1.0;
        info.apply_qf = has_b ? f_eval_bdr_current_j_2_32 : f_eval_bdr_current_j_1_32;
        info.apply_qf_path = PalaceQFunctionRelativePath(
            has_b ? f_eval_bdr_current_j_2_32_loc : f_eval_bdr_current_j_1_32_loc);
        break;
      case Kind::BDR_ENERGY_E:
        info.apply_qf = has_b ? f_eval_bdr_energy_e_2_32 : f_eval_bdr_energy_e_1_32;
        info.apply_qf_path = PalaceQFunctionRelativePath(
            has_b ? f_eval_bdr_energy_e_2_32_loc : f_eval_bdr_energy_e_1_32_loc);
        break;
      case Kind::BDR_ENERGY_M:
        info.apply_qf = has_b ? f_eval_bdr_energy_m_2_32 : f_eval_bdr_energy_m_1_32;
        info.apply_qf_path = PalaceQFunctionRelativePath(
            has_b ? f_eval_bdr_energy_m_2_32_loc : f_eval_bdr_energy_m_1_32_loc);
        break;
      case Kind::FARFIELD:
        ctx[0].second = group.flip_normal ? -1.0 : 1.0;
        info.apply_qf = f_integ_surf_farfield_32;
        info.apply_qf_path = PalaceQFunctionRelativePath(f_integ_surf_farfield_32_loc);
        break;
      case Kind::SURFACE_FLUX:
        ctx[0].second = group.flip_normal ? -1.0 : 1.0;
        switch (flux_type)
        {
          case SurfaceFlux::ELECTRIC:
            info.apply_qf = has_b ? f_integ_surf_flux_e_2_32 : f_integ_surf_flux_e_1_32;
            info.apply_qf_path = PalaceQFunctionRelativePath(
                has_b ? f_integ_surf_flux_e_2_32_loc : f_integ_surf_flux_e_1_32_loc);
            break;
          case SurfaceFlux::MAGNETIC:
            info.apply_qf = has_b ? f_integ_surf_flux_m_2_32 : f_integ_surf_flux_m_1_32;
            info.apply_qf_path = PalaceQFunctionRelativePath(
                has_b ? f_integ_surf_flux_m_2_32_loc : f_integ_surf_flux_m_1_32_loc);
            break;
          case SurfaceFlux::POWER:
            info.apply_qf = has_b ? f_integ_surf_flux_p_2_32 : f_integ_surf_flux_p_1_32;
            info.apply_qf_path = PalaceQFunctionRelativePath(
                has_b ? f_integ_surf_flux_p_2_32_loc : f_integ_surf_flux_p_1_32_loc);
            break;
        }
        break;
    }

    // Assemble the operator.
    CeedOperator op;
    if (buffer_kind)
    {
      ceed::AssembleCeedPointEvaluator(info, ctx.data(), ctx.size() * sizeof(CeedIntScalar),
                                       ceed, inputs, BufferNumComp(kind), out_restr, &op);
    }
    else
    {
      ceed::AssembleCeedSurfaceFunctional(info, ctx.data(),
                                          ctx.size() * sizeof(CeedIntScalar), ceed, inputs,
                                          num_out, out_restr, &op);
    }
    groups.push_back({ceed, op, std::move(field_sources)});
  }
}

void SurfaceFunctional::ApplyAdd(const std::array<const Vector *, 4> &srcs) const
{
  ApplyAddGroups(groups, srcs, local_out);
}

double SurfaceFunctional::EvalLocal(const std::array<const Vector *, 4> &srcs) const
{
  MFEM_VERIFY(valid, "Eval called on an invalid (unassembled) SurfaceFunctional!");
  if (local_out.Size() == 0)
  {
    return 0.0;
  }
  local_out = 0.0;
  ApplyAdd(srcs);
  return linalg::LocalSum(local_out);
}

double SurfaceFunctional::Eval(const Vector *u) const
{
  MFEM_VERIFY(kind == Kind::AREA || u,
              "SurfaceFunctional::Eval requires a field vector for functionals with "
              "field inputs!");
  double dot = EvalLocal({u, nullptr});
  Mpi::GlobalSum(1, &dot, comm);
  return dot;
}

double SurfaceFunctional::Eval(const GridFunction &u) const
{
  MFEM_VERIFY(valid, "Eval called on an invalid (unassembled) SurfaceFunctional!");
  MFEM_VERIFY(kind == Kind::HCURL_NORM2 || kind == Kind::INTERFACE_EPR,
              "SurfaceFunctional::Eval with a grid function is only valid for functionals "
              "quadratic in a single field!");
  double dot = 0.0;
  if (local_out.Size() > 0)
  {
    local_out = 0.0;
    ApplyAdd({&u.Real(), nullptr});
    if (u.HasImag())
    {
      ApplyAdd({&u.Imag(), nullptr});
    }
    dot = linalg::LocalSum(local_out);
  }
  Mpi::GlobalSum(1, &dot, comm);
  return dot;
}

std::complex<double> SurfaceFunctional::EvalFlux(const GridFunction *E,
                                                 const GridFunction *B) const
{
  MFEM_VERIFY(kind == Kind::SURFACE_FLUX,
              "SurfaceFunctional::EvalFlux is only valid for surface flux functionals!");
  MFEM_VERIFY(
      (E || (flux_type != SurfaceFlux::ELECTRIC && flux_type != SurfaceFlux::POWER)) &&
          (B || (flux_type != SurfaceFlux::MAGNETIC && flux_type != SurfaceFlux::POWER)),
      "Missing E or B field grid function for surface flux evaluation!");

  // For complex-valued fields, output the separate real and imaginary parts for the
  // time-harmonic quantity. For power flux (Poynting vector), output only the
  // stationary real part and not the part which has double the frequency (the real and
  // imaginary part contributions add).
  const bool has_imag = E ? E->HasImag() : B->HasImag();
  std::complex<double> dot(EvalLocal({E ? &E->Real() : nullptr, B ? &B->Real() : nullptr}),
                           0.0);
  if (has_imag)
  {
    const double doti = EvalLocal({E ? &E->Imag() : nullptr, B ? &B->Imag() : nullptr});
    if (flux_type == SurfaceFlux::POWER)
    {
      dot += doti;
    }
    else
    {
      dot.imag(doti);
    }
  }
  Mpi::GlobalSum(1, &dot, comm);
  return dot;
}

void SurfaceFunctional::EvalBuffer(const Vector &u, Vector &buffer) const
{
  MFEM_VERIFY(valid && IsBufferKind(kind),
              "EvalBuffer requires a valid boundary visualization field functional!");
  MFEM_ASSERT(buffer.Size() == buffer_size, "Invalid buffer size for EvalBuffer!");
  buffer = 0.0;
  ApplyAddGroups(groups, {&u}, buffer);
}

void SurfaceFunctional::EvalBuffer(const GridFunction &u, Vector &buffer) const
{
  MFEM_VERIFY(valid && IsBufferKind(kind),
              "EvalBuffer requires a valid boundary visualization field functional!");
  MFEM_ASSERT(buffer.Size() == buffer_size, "Invalid buffer size for EvalBuffer!");
  buffer = 0.0;
  ApplyAddGroups(groups, {&u.Real()}, buffer);
  if (u.HasImag())
  {
    ApplyAddGroups(groups, {&u.Imag()}, buffer);
  }
}

std::vector<std::array<std::complex<double>, 3>>
SurfaceFunctional::EvalFarField(const GridFunction &E, const GridFunction &B,
                                double omega_re, double omega_im)
{
  MFEM_VERIFY(kind == Kind::FARFIELD && E.HasImag() && B.HasImag(),
              "SurfaceFunctional::EvalFarField requires a far-field functional and "
              "complex-valued fields!");
  MFEM_VERIFY(valid, "EvalFarField called on an invalid (unassembled) SurfaceFunctional!");

  // The frequency enters the QFunction context: reassemble when it changes.
  if (omega_re != farfield_omega_re || omega_im != farfield_omega_im)
  {
    for (auto &group : groups)
    {
      PalaceCeedCall(group.ceed, CeedOperatorDestroy(&group.op));
    }
    groups.clear();
    elem_attrs.clear();
    farfield_omega_re = omega_re;
    farfield_omega_im = omega_im;
    Assemble(*farfield_mesh, farfield_marker);
  }

  // Integrate, reduce each component over the local elements and all processes, and
  // apply the final cross products (following GetFarFieldrE).
  const int N = static_cast<int>(farfield_dirs.size());
  const int num_marked = local_out.Size() / std::max(6 * N, 1);
  std::vector<double> integrals(6 * N, 0.0);
  if (local_out.Size() > 0)
  {
    local_out = 0.0;
    ApplyAddGroups(groups, {&E.Real(), &E.Imag(), &B.Real(), &B.Imag()}, local_out);
    Vector slice;
    slice.UseDevice(true);
    for (int c = 0; c < 6 * N; c++)
    {
      slice.MakeRef(local_out, c * num_marked, num_marked);
      integrals[c] = linalg::LocalSum(slice);
    }
  }
  Mpi::GlobalSum(6 * N, integrals.data(), comm);

  std::vector<std::array<std::complex<double>, 3>> result(N);
  for (int d = 0; d < N; d++)
  {
    const auto &r = farfield_dirs[d];
    const double *Ir = integrals.data() + 6 * d, *Ii = integrals.data() + 6 * d + 3;
    const double cr[3] = {r[1] * Ir[2] - r[2] * Ir[1], r[2] * Ir[0] - r[0] * Ir[2],
                          r[0] * Ir[1] - r[1] * Ir[0]};
    const double ci[3] = {r[1] * Ii[2] - r[2] * Ii[1], r[2] * Ii[0] - r[0] * Ii[2],
                          r[0] * Ii[1] - r[1] * Ii[0]};
    for (int c = 0; c < 3; c++)
    {
      result[d][c] = {cr[c], ci[c]};
    }
  }
  return result;
}

std::complex<double> SurfaceFunctional::EvalComplexPower(const GridFunction &E,
                                                         const GridFunction &B) const
{
  MFEM_VERIFY(kind == Kind::SURFACE_FLUX && flux_type == SurfaceFlux::POWER &&
                  flux_two_sided,
              "SurfaceFunctional::EvalComplexPower is only valid for two-sided POWER "
              "flux functionals!");
  MFEM_VERIFY(E.HasImag() == B.HasImag(),
              "Mismatch between real- and complex-valued E and B fields in port power "
              "calculation!");

  // Following LumpedPortData::GetPower: P = ∫ E ⋅ (n x H) dS with H = μ⁻¹ B and n the
  // normal oriented into element 1 (contributions from both sides of an interior
  // boundary add). With S(e, b) = ∫ (e x μ⁻¹ b) ⋅ n dS (the two-sided POWER flux
  // functional), E ⋅ (n x H) = -(E x H) ⋅ n gives
  //   Re{P} = S(E_re, B_re) + S(E_im, B_im)
  //   Im{P} = S(E_im, B_re) - S(E_re, B_im) .
  const bool has_imag = E.HasImag();
  std::complex<double> dot(EvalLocal({&E.Real(), &B.Real()}), 0.0);
  if (has_imag)
  {
    dot += EvalLocal({&E.Imag(), &B.Imag()});
    dot.imag(EvalLocal({&E.Imag(), &B.Real()}) - EvalLocal({&E.Real(), &B.Imag()}));
  }
  Mpi::GlobalSum(1, &dot, comm);
  return dot;
}

DomainFieldEvaluator::DomainFieldEvaluator(
    Kind kind, const Mesh &mesh, const MaterialOperator &mat_op,
    const mfem::ParFiniteElementSpace *nd_fespace,
    const mfem::ParFiniteElementSpace *rt_fespace,
    const mfem::ParFiniteElementSpace &target_fespace, double scaling)
  : kind(kind), fespace_e(nd_fespace), fespace_b(rt_fespace)
{
  MFEM_VERIFY((nd_fespace || kind == Kind::ENERGY_M) &&
                  (rt_fespace || kind == Kind::ENERGY_E),
              "Missing finite element space for domain field evaluator!");
  Assemble(mesh, mat_op, target_fespace, scaling);
}

DomainFieldEvaluator::~DomainFieldEvaluator()
{
  for (auto &group : groups)
  {
    PalaceCeedCall(group.ceed, CeedOperatorDestroy(&group.op));
  }
}

void DomainFieldEvaluator::Assemble(const Mesh &mesh, const MaterialOperator &mat_op,
                                    const mfem::ParFiniteElementSpace &target_fespace,
                                    double scaling)
{
  if (mesh.Dimension() != 3 || mesh.SpaceDimension() != 3)
  {
    valid = false;
    return;
  }
  const mfem::ParMesh &pmesh = mesh.Get();
  const mfem::FiniteElementSpace &mesh_fespace = *pmesh.GetNodes()->FESpace();

  // Group the elements by geometry type.
  std::map<mfem::Geometry::Type, std::vector<int>> geom_elems;
  for (int e = 0; e < pmesh.GetNE(); e++)
  {
    geom_elems[pmesh.GetElementGeometry(e)].push_back(e);
  }

  // QFunction context: scaling factor followed by the material property table.
  std::vector<CeedIntScalar> ctx(1);
  ctx[0].second = scaling;
  {
    MaterialPropertyCoefficient coeff_func(mat_op.GetAttributeToMaterial(),
                                           kind == Kind::ENERGY_E
                                               ? mat_op.GetPermittivityReal()
                                               : mat_op.GetInvPermeability());
    auto mat_ctx = ceed::PopulateCoefficientContext(3, &coeff_func);
    ctx.insert(ctx.end(), mat_ctx.begin(), mat_ctx.end());
  }

  field_staging.SetSize(std::max(fespace_e ? fespace_e->GetVSize() : 0,
                                 fespace_b ? fespace_b->GetVSize() : 0));
  field_staging.UseDevice(true);
  field_staging = 0.0;

  Ceed ceed = ceed::internal::GetCeedObjects()[0];
  for (const auto &[geom, indices] : geom_elems)
  {
    CeedAssemblyScratch scratch(ceed);

    // Evaluation points are the nodal points of the (interpolatory) target space.
    const mfem::FiniteElement *target_fe =
        target_fespace.FEColl()->FiniteElementForGeometry(geom);
    MFEM_VERIFY(target_fe, "Unable to get target finite element for field evaluator!");
    const mfem::IntegrationRule &nodes_ir = target_fe->GetNodes();
    const int num_pts = nodes_ir.GetNPoints();

    // Element attributes (libCEED local, for material lookup) and mesh node gradients
    // (for on the fly geometry evaluation at the points).
    std::vector<ceed::CeedFunctionalFieldInput> inputs;
    std::vector<std::pair<std::string, int>> field_sources;
    {
      auto &elem_attr = elem_attrs.emplace_back(indices.size());
      const auto &loc_attr = mesh.GetCeedAttributes();
      for (std::size_t k = 0; k < indices.size(); k++)
      {
        elem_attr[k] = loc_attr.at(pmesh.GetAttribute(indices[k]));
      }
      CeedElemRestriction attr_restr;
      CeedBasis attr_basis;
      CeedVector attr_vec;
      PalaceCeedCall(
          ceed, CeedElemRestrictionCreateStrided(ceed, indices.size(), 1, 1, indices.size(),
                                                 CEED_STRIDES_BACKEND, &attr_restr));
      {
        // Note: ceed::GetCeedTopology(CEED_TOPOLOGY_LINE) == 1.
        mfem::Vector Bt(num_pts), Gt(num_pts), qX(num_pts), qW(num_pts);
        Bt = 1.0;
        Gt = 0.0;
        qX = 0.0;
        qW = 0.0;
        PalaceCeedCall(ceed, CeedBasisCreateH1(ceed, CEED_TOPOLOGY_LINE, 1, 1, num_pts,
                                               Bt.GetData(), Gt.GetData(), qX.GetData(),
                                               qW.GetData(), &attr_basis));
      }
      ceed::InitCeedVector(elem_attr, ceed, &attr_vec);
      inputs.push_back({"attr", attr_vec, attr_restr, attr_basis, ceed::EvalMode::Interp});
      scratch.vecs.push_back(attr_vec);
      scratch.restrs.push_back(attr_restr);
      scratch.bases.push_back(attr_basis);
    }
    {
      CeedElemRestriction mesh_restr = FiniteElementSpace::BuildCeedElemRestriction(
          mesh_fespace, ceed, geom, indices, /*is_interp*/ true);
      const mfem::FiniteElement *mesh_fe =
          mesh_fespace.FEColl()->FiniteElementForGeometry(geom);
      CeedBasis mesh_basis;
      ceed::InitBasisAtPoints(*mesh_fe, nodes_ir, mesh_fespace.GetVDim(), ceed,
                              &mesh_basis);
      CeedVector mesh_nodes_vec;
      ceed::InitCeedVector(*mesh_fespace.GetMesh()->GetNodes(), ceed, &mesh_nodes_vec);
      inputs.push_back({"x", mesh_nodes_vec, mesh_restr, mesh_basis, ceed::EvalMode::Grad});
      scratch.vecs.push_back(mesh_nodes_vec);
      scratch.restrs.push_back(mesh_restr);
      scratch.bases.push_back(mesh_basis);
    }

    // Field inputs evaluated at the nodal points.
    auto AddFieldInput =
        [&](const std::string &name, int source, const mfem::ParFiniteElementSpace &fespace)
    {
      CeedElemRestriction restr;
      CeedBasis basis;
      CeedVector vec;
      ceed::InitRestriction(fespace, indices, false, /*is_interp*/ true, false, ceed,
                            &restr);
      const mfem::FiniteElement *fe = fespace.FEColl()->FiniteElementForGeometry(geom);
      ceed::InitBasisAtPoints(*fe, nodes_ir, fespace.GetVDim(), ceed, &basis);
      ceed::InitCeedVector(field_staging, ceed, &vec);
      inputs.push_back({name, vec, restr, basis, ceed::EvalMode::Interp});
      field_sources.emplace_back(name, source);
      scratch.vecs.push_back(vec);
      scratch.restrs.push_back(restr);
      scratch.bases.push_back(basis);
    };
    if (kind == Kind::ENERGY_E || kind == Kind::POYNTING)
    {
      AddFieldInput("u_1", 0, *fespace_e);
    }
    if (kind == Kind::ENERGY_M || kind == Kind::POYNTING)
    {
      AddFieldInput(kind == Kind::POYNTING ? "u_2" : "u_1", 1, *fespace_b);
    }

    // Output restriction scattering nodal values into the target grid function.
    CeedElemRestriction out_restr;
    ceed::InitRestriction(target_fespace, indices, false, /*is_interp*/ true, false, ceed,
                          &out_restr);
    scratch.restrs.push_back(out_restr);

    ceed::CeedQFunctionInfo info;
    switch (kind)
    {
      case Kind::ENERGY_E:
        info.apply_qf = f_eval_energy_e_33;
        info.apply_qf_path = PalaceQFunctionRelativePath(f_eval_energy_e_33_loc);
        break;
      case Kind::ENERGY_M:
        info.apply_qf = f_eval_energy_m_33;
        info.apply_qf_path = PalaceQFunctionRelativePath(f_eval_energy_m_33_loc);
        break;
      case Kind::POYNTING:
        info.apply_qf = f_eval_poynting_33;
        info.apply_qf_path = PalaceQFunctionRelativePath(f_eval_poynting_33_loc);
        break;
    }

    CeedOperator op;
    ceed::AssembleCeedPointEvaluator(info, ctx.data(), ctx.size() * sizeof(CeedIntScalar),
                                     ceed, inputs, target_fespace.GetVDim(), out_restr,
                                     &op);
    groups.push_back({ceed, op, std::move(field_sources)});
  }
}

void DomainFieldEvaluator::Eval(const GridFunction *E, const GridFunction *B,
                                Vector &out) const
{
  MFEM_VERIFY(valid, "Eval called on an invalid (unassembled) DomainFieldEvaluator!");
  MFEM_VERIFY((E || kind == Kind::ENERGY_M) && (B || kind == Kind::ENERGY_E),
              "Missing field grid function for domain field evaluator!");
  out = 0.0;
  ApplyAddGroups(groups, {E ? &E->Real() : nullptr, B ? &B->Real() : nullptr}, out);
  if (E ? E->HasImag() : B->HasImag())
  {
    ApplyAddGroups(groups, {E ? &E->Imag() : nullptr, B ? &B->Imag() : nullptr}, out);
  }
}

}  // namespace palace
