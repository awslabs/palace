// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "output_functionals.hpp"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <map>
#include <memory>
#include <mutex>
#include "fem/coefficient.hpp"
#include "fem/facenbrexchange.hpp"
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
  // A side whose volume element is a face neighbor (ghost) on another process: its field
  // is not in the local vector, so it is pulled via FaceNbrFieldExchange and fed to the
  // operator at the mapped points (only Elem2 can be a ghost, so at most one side).
  bool ghost_a = false, ghost_b = false;
  int face_nbr = -1;    // Face neighbor element index (Elem2No - ParMesh::GetNE())
  int ghost_attr = 0;   // Ghost element mesh attribute (material lookup on the requester)
  mfem::Geometry::Type ghost_geom = mfem::Geometry::INVALID;  // Ghost element geometry
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
  // Ghost (face neighbor) side, if any: at most one of ghost_a/ghost_b is set (only
  // Elem2 can be a ghost). For the ghost side, vol_indices are empty; face_nbr and
  // ghost_attr hold the per-element face neighbor index and mesh attribute, and req_idx
  // the FaceNbrFieldExchange request index (filled when the exchange is built).
  bool ghost_a = false, ghost_b = false;
  bool at_points = false;
  // Split two-sided boundary-visualization groups into one-sided AtPoints operators.
  // ApplyAdd accumulation combines the side contributions in the shared output buffer;
  // side_scale handles averages, normal_scale handles signed jump quantities.
  double side_scale = 1.0;
  double normal_scale = 1.0;
  std::vector<int> face_nbr, ghost_attr, req_idx;
  // For AtPoints groups, the mapped volume reference coordinates vary by boundary
  // element and are stored in boundary-element order, nq entries per element. Local
  // split SURFACE_FLUX groups may also carry a per-entry normal scale folded into the
  // precomputed face Jacobian so opposite sides can share one AtPoints operator.
  std::vector<mfem::IntegrationPoint> mapped_pts_a, mapped_pts_b;
  std::vector<double> normal_scales;
};

void FoldNormalScaleIntoFaceJacobian(const mfem::DenseMatrix &J, double normal_scale,
                                     double *Jf)
{
  // Surface flux kernels recover the unnormalized physical normal from the columns of
  // the face Jacobian. Multiplying one column by s multiplies both the surface measure
  // and the oriented normal by s, including flipping the normal direction for s < 0.
  // Split local two-sided SURFACE_FLUX AtPoints groups use this to carry the per-entry
  // signed jump/average convention in geometry while sharing one operator across sides.
  for (int d = 0; d < 2; d++)
  {
    for (int c = 0; c < 3; c++)
    {
      Jf[c + 3 * d] = (d == 0 ? normal_scale : 1.0) * J(c, d);
    }
  }
}

void AppendPoints(FaceConfigKey &key, const std::vector<mfem::IntegrationPoint> &pts)
{
  for (const auto &ip : pts)
  {
    key.push_back(std::llround(ip.x * QUANTIZE_SCALE));
    key.push_back(std::llround(ip.y * QUANTIZE_SCALE));
    key.push_back(std::llround(ip.z * QUANTIZE_SCALE));
  }
}

const char *KindName(SurfaceFunctional::Kind kind)
{
  switch (kind)
  {
    case SurfaceFunctional::Kind::AREA:
      return "AREA";
    case SurfaceFunctional::Kind::HCURL_NORM2:
      return "HCURL_NORM2";
    case SurfaceFunctional::Kind::INTERFACE_EPR:
      return "INTERFACE_EPR";
    case SurfaceFunctional::Kind::SURFACE_FLUX:
      return "SURFACE_FLUX";
    case SurfaceFunctional::Kind::FARFIELD:
      return "FARFIELD";
    case SurfaceFunctional::Kind::BDR_FIELD_E:
      return "BDR_FIELD_E";
    case SurfaceFunctional::Kind::BDR_FIELD_B:
      return "BDR_FIELD_B";
    case SurfaceFunctional::Kind::BDR_FLUX_Q:
      return "BDR_FLUX_Q";
    case SurfaceFunctional::Kind::BDR_CURRENT_J:
      return "BDR_CURRENT_J";
    case SurfaceFunctional::Kind::BDR_ENERGY_E:
      return "BDR_ENERGY_E";
    case SurfaceFunctional::Kind::BDR_ENERGY_M:
      return "BDR_ENERGY_M";
    case SurfaceFunctional::Kind::BDR_POYNTING:
      return "BDR_POYNTING";
  }
  return "UNKNOWN";
}

bool CeedSupportsNonTensorAtPoints(Ceed ceed)
{
  // This proof path uses the non-tensor simplex AtPoints kernels added for CUDA ref
  // and MAGMA. Other backends keep the existing mapped-point tabulation path.
  const char *resource;
  PalaceCeedCall(ceed, CeedGetResource(ceed, &resource));
  return std::strstr(resource, "/gpu/cuda/ref") ||
         std::strstr(resource, "/gpu/cuda/magma");
}

int TetNumModes(int degree)
{
  return (degree + 1) * (degree + 2) * (degree + 3) / 6;
}

mfem::IntegrationRule MakeTetLatticeRule(int degree)
{
  mfem::IntegrationRule ir(TetNumModes(degree));
  int q = 0;
  if (degree == 0)
  {
    ir.IntPoint(q).Set3(0.25, 0.25, 0.25);
    ir.IntPoint(q++).weight = 1.0;
  }
  else
  {
    for (int total = 0; total <= degree; total++)
    {
      for (int i = 0; i <= total; i++)
      {
        for (int j = 0; j <= total - i; j++)
        {
          const int k = total - i - j;
          ir.IntPoint(q).Set3(static_cast<double>(i) / degree,
                              static_cast<double>(j) / degree,
                              static_cast<double>(k) / degree);
          ir.IntPoint(q++).weight = 1.0;
        }
      }
    }
  }
  return ir;
}

void InitTetBasisForAtPoints(const mfem::FiniteElement &fe, bool grad_only,
                             CeedInt num_comp, Ceed ceed, CeedBasis *basis)
{
  MFEM_VERIFY(fe.GetGeomType() == mfem::Geometry::TETRAHEDRON,
              "AtPoints surface output proof currently supports tetrahedral volume "
              "elements only!");
  const int degree = std::max(0, fe.GetOrder() - (grad_only ? 1 : 0));
  const mfem::IntegrationRule ir = MakeTetLatticeRule(degree);
  ceed::InitBasisAtPoints(fe, ir, num_comp, ceed, basis);
}

void CreateSequentialPointRestriction(Ceed ceed, std::size_t num_elem, int nq, int num_comp,
                                      CeedSize l_size, CeedElemRestriction *restr)
{
  const CeedInt total_pts = static_cast<CeedInt>(num_elem * nq);
  std::vector<CeedInt> offsets(num_elem + 1 + total_pts);
  for (std::size_t e = 0; e <= num_elem; e++)
  {
    offsets[e] = static_cast<CeedInt>(num_elem + 1 + e * nq);
  }
  for (CeedInt i = 0; i < total_pts; i++)
  {
    offsets[num_elem + 1 + i] = i;
  }
  PalaceCeedCall(ceed, CeedElemRestrictionCreateAtPoints(
                           ceed, static_cast<CeedInt>(num_elem), total_pts, num_comp,
                           l_size, CEED_MEM_HOST, CEED_COPY_VALUES, offsets.data(), restr));
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

void fem::ApplyAddGroupOperators(const std::vector<fem::CeedGroupOperator> &groups,
                                 const std::array<const Vector *, 4> &srcs,
                                 const Vector &out, const Vector *imported)
{
  for (const auto &group : groups)
  {
    for (const auto &[name, source] : group.field_sources)
    {
      // Source index 4 selects the imported face neighbor field values (see
      // SurfaceFunctional::face_nbr_exchange); the operator's restriction slices and
      // transposes the shared vector to the per-element layout.
      const Vector *sv = (source < 4) ? srcs[source] : imported;
      MFEM_ASSERT(sv, "Missing source vector for libCEED field input!");
      CeedOperatorField field;
      CeedVector field_vec;
      PalaceCeedCall(group.ceed,
                     CeedOperatorGetFieldByName(group.op, name.c_str(), &field));
      PalaceCeedCall(group.ceed, CeedOperatorFieldGetVector(field, &field_vec));
      ceed::InitCeedVector(*sv, group.ceed, &field_vec, false);
    }
    CeedMemType out_mem;
    PalaceCeedCall(group.ceed, CeedGetPreferredMemType(group.ceed, &out_mem));
    if (!mfem::Device::Allows(mfem::Backend::DEVICE_MASK) && out_mem == CEED_MEM_DEVICE)
    {
      out_mem = CEED_MEM_HOST;
    }
    auto *out_data = const_cast<Vector &>(out).ReadWrite(out_mem == CEED_MEM_DEVICE);
    const CeedSize out_size = out.Size();
    if (!group.out_vec || group.out_size != out_size)
    {
      if (group.out_vec)
      {
        PalaceCeedCall(group.ceed, CeedVectorDestroy(&group.out_vec));
      }
      PalaceCeedCall(group.ceed, CeedVectorCreate(group.ceed, out_size, &group.out_vec));
      group.out_size = out_size;
    }
    PalaceCeedCall(group.ceed,
                   CeedVectorSetArray(group.out_vec, out_mem, CEED_USE_POINTER, out_data));
    PalaceCeedCall(group.ceed, CeedOperatorApplyAdd(group.op, CEED_VECTOR_NONE,
                                                    group.out_vec,
                                                    CEED_REQUEST_IMMEDIATE));
    PalaceCeedCall(group.ceed, CeedVectorTakeArray(group.out_vec, out_mem, nullptr));
  }
}

SurfaceFunctional::SurfaceFunctional(Kind kind, const Mesh &mesh,
                                     const mfem::Array<int> &bdr_attr_marker,
                                     const mfem::ParFiniteElementSpace *fespace)
  : kind(kind), nd_fespace(fespace), rt_fespace(nullptr), mat_op(nullptr),
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
  : kind(kind), nd_fespace(kind == Kind::BDR_FIELD_E ? &fespace : nullptr),
    rt_fespace(kind == Kind::BDR_FIELD_B ? &fespace : nullptr), mat_op(nullptr),
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
    nd_fespace((kind == Kind::BDR_FLUX_Q || kind == Kind::BDR_ENERGY_E) ? &fespace
                                                                       : nullptr),
    rt_fespace((kind == Kind::BDR_CURRENT_J || kind == Kind::BDR_ENERGY_M) ? &fespace
                                                                          : nullptr),
    mat_op(&mat_op), comm(mesh.GetComm()), viz_lod(lod), viz_scaling(scaling)
{
  MFEM_VERIFY(kind == Kind::BDR_FLUX_Q || kind == Kind::BDR_CURRENT_J ||
                  kind == Kind::BDR_ENERGY_E || kind == Kind::BDR_ENERGY_M,
              "Invalid SurfaceFunctional constructor for the requested functional kind!");
  Assemble(mesh, bdr_attr_marker);
}

SurfaceFunctional::SurfaceFunctional(Kind kind, const Mesh &mesh,
                                     const mfem::Array<int> &bdr_attr_marker,
                                     const mfem::ParFiniteElementSpace &nd_fespace,
                                     const mfem::ParFiniteElementSpace &rt_fespace,
                                     const MaterialOperator &mat_op, int lod,
                                     double scaling)
  : kind(kind), nd_fespace(&nd_fespace), rt_fespace(&rt_fespace), mat_op(&mat_op),
    comm(mesh.GetComm()), viz_lod(lod), viz_scaling(scaling)
{
  MFEM_VERIFY(kind == Kind::BDR_POYNTING,
              "Invalid SurfaceFunctional constructor for the requested functional kind!");
  Assemble(mesh, bdr_attr_marker);
}

SurfaceFunctional::SurfaceFunctional(const Mesh &mesh,
                                     const mfem::Array<int> &bdr_attr_marker,
                                     const mfem::ParFiniteElementSpace &nd_fespace,
                                     const MaterialOperator &mat_op,
                                     InterfaceDielectric type, double t_i, double epsilon_i)
  : kind(Kind::INTERFACE_EPR), epr_type(type), epr_t(t_i), epr_epsilon(epsilon_i),
    nd_fespace(&nd_fespace), rt_fespace(nullptr), mat_op(&mat_op), comm(mesh.GetComm())
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
    nd_fespace(nd_fespace), rt_fespace(rt_fespace), mat_op(&mat_op), comm(mesh.GetComm())
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
  : kind(Kind::FARFIELD), farfield_dirs(r_naughts), nd_fespace(&nd_fespace),
    rt_fespace(&rt_fespace), mat_op(&mat_op), comm(mesh.GetComm())
{
  Assemble(mesh, bdr_attr_marker);
}

SurfaceFunctional::~SurfaceFunctional()
{
  for (auto &group : groups)
  {
    PalaceCeedCall(group.ceed, CeedOperatorDestroy(&group.op));
    if (group.ctx)
    {
      PalaceCeedCall(group.ceed, CeedQFunctionContextDestroy(&group.ctx));
    }
    if (group.out_vec)
    {
      PalaceCeedCall(group.ceed, CeedVectorDestroy(&group.out_vec));
    }
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
      if (group.ctx)
      {
        PalaceCeedCall(group.ceed, CeedQFunctionContextDestroy(&group.ctx));
      }
      if (group.out_vec)
      {
        PalaceCeedCall(group.ceed, CeedVectorDestroy(&group.out_vec));
      }
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
  Ceed ceed = ceed::internal::GetCeedObjects()[0];
  const bool use_at_points = CeedSupportsNonTensorAtPoints(ceed);

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

      // A two-sided interior boundary on a parallel interface: the Elem2 side is a face
      // neighbor (ghost) element on another process. Its field is not in the local
      // vector; it is pulled via FaceNbrFieldExchange (built below from the mapped
      // points) and fed to the ghost side of the operator, so the libCEED path handles
      // it without falling back to the legacy coefficients. Only Elem2 can be a ghost.
      if (plan.elem_a >= 0 && plan.elem_a >= pmesh.GetNE())
      {
        plan.ghost_a = true;
      }
      if (plan.elem_b >= 0 && !elem2_local)
      {
        plan.ghost_b = true;
      }
      if (plan.ghost_a || plan.ghost_b)
      {
        plan.face_nbr = FET.Elem2No - pmesh.GetNE();
        plan.ghost_attr = FET.Elem2->Attribute;
        // The ghost element index is out of range for pmesh.GetElementGeometry; capture
        // its geometry here while the face neighbor transformation is in scope.
        plan.ghost_geom = FET.Elem2->GetGeometryType();
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

      // Build the group key. For AtPoints-capable boundary visualization groups, the
      // mapped volume reference coordinates are runtime data rather than part of the
      // basis/JIT key; older paths still key on the mapped coordinates themselves.
      const auto vol_geom_a = plan.ghost_a ? plan.ghost_geom
                              : (plan.elem_a >= 0) ? pmesh.GetElementGeometry(plan.elem_a)
                                                   : mfem::Geometry::INVALID;
      const auto vol_geom_b = plan.ghost_b ? plan.ghost_geom
                              : (plan.elem_b >= 0) ? pmesh.GetElementGeometry(plan.elem_b)
                                                   : mfem::Geometry::INVALID;
      const bool can_surface_flux_at_points_a =
          use_at_points && kind == Kind::SURFACE_FLUX && plan.elem_a >= 0 &&
          !plan.ghost_a && vol_geom_a == mfem::Geometry::TETRAHEDRON;
      const bool can_surface_flux_at_points_b =
          use_at_points && kind == Kind::SURFACE_FLUX && plan.elem_b >= 0 &&
          !plan.ghost_b && vol_geom_b == mfem::Geometry::TETRAHEDRON;
      const bool can_epr_at_points =
          use_at_points && kind == Kind::INTERFACE_EPR && plan.elem_a >= 0 &&
          plan.elem_b < 0 && !plan.ghost_a &&
          vol_geom_a == mfem::Geometry::TETRAHEDRON;
      const bool can_farfield_at_points =
          use_at_points && kind == Kind::FARFIELD && plan.elem_a >= 0 &&
          plan.elem_b < 0 && !plan.ghost_a &&
          vol_geom_a == mfem::Geometry::TETRAHEDRON;
      const bool can_at_points_a =
          (use_at_points && buffer_kind && plan.elem_a >= 0 && !plan.ghost_a &&
           vol_geom_a == mfem::Geometry::TETRAHEDRON) ||
          can_surface_flux_at_points_a || can_epr_at_points || can_farfield_at_points;
      const bool can_at_points_b =
          (use_at_points && buffer_kind && plan.elem_b >= 0 && !plan.ghost_b &&
           vol_geom_b == mfem::Geometry::TETRAHEDRON) ||
          can_surface_flux_at_points_b;

      int out_slot;
      if (IsBufferKind(kind))
      {
        if (buffer_bases.empty())
        {
          buffer_bases.resize(pmesh.GetNBE(), -1);
        }
        out_slot = buffer_size;
        buffer_bases[i] = out_slot;
        buffer_size += nq * BufferNumComp(kind);
        num_marked++;
      }
      else
      {
        out_slot = num_marked++;
        local_out_attrs.push_back(attr);
      }

      auto AddGroup = [&](int elem_a, bool ghost_a, mfem::Geometry::Type geom_a,
                          const std::vector<mfem::IntegrationPoint> &pts_a, int elem_b,
                          bool ghost_b, mfem::Geometry::Type geom_b,
                          const std::vector<mfem::IntegrationPoint> &pts_b,
                          bool at_points_group, double side_scale, double normal_scale)
      {
        FaceConfigKey key;
        key.reserve(10 + (at_points_group ? 0 : 3 * (pts_a.size() + pts_b.size())));
        key.push_back(static_cast<long long>(bdr_geom));
        key.push_back(static_cast<long long>(geom_a));
        key.push_back(static_cast<long long>(geom_b));
        key.push_back(static_cast<long long>(plan.flip));
        key.push_back(static_cast<long long>(ghost_a));
        key.push_back(static_cast<long long>(ghost_b));
        key.push_back(static_cast<long long>(nq));
        key.push_back(static_cast<long long>(at_points_group));
        key.push_back(static_cast<long long>(std::llround(side_scale * QUANTIZE_SCALE)));
        key.push_back(static_cast<long long>(
            (at_points_group && kind == Kind::SURFACE_FLUX)
                ? 0
                : std::llround(normal_scale * QUANTIZE_SCALE)));
        if (!at_points_group)
        {
          AppendPoints(key, pts_a);
          AppendPoints(key, pts_b);
        }

        auto it = face_groups.find(key);
        if (it == face_groups.end())
        {
          FaceGroup group;
          group.bdr_geom = bdr_geom;
          group.vol_geom_a = geom_a;
          group.vol_geom_b = geom_b;
          group.flip_normal = plan.flip;
          group.ghost_a = ghost_a;
          group.ghost_b = ghost_b;
          group.at_points = at_points_group;
          group.side_scale = side_scale;
          group.normal_scale = normal_scale;
          if (!at_points_group && elem_a >= 0)
          {
            FaceConfigKey key_a = key;
            key_a.push_back(0);  // Side tag
            group.mapped_ir_a = GetRegisteredMappedIr(key_a, pts_a);
          }
          if (!at_points_group && elem_b >= 0)
          {
            FaceConfigKey key_b = key;
            key_b.push_back(1);  // Side tag
            group.mapped_ir_b = GetRegisteredMappedIr(key_b, pts_b);
          }
          it = face_groups.emplace(key, std::move(group)).first;
        }
        it->second.bdr_indices.push_back(i);
        if (at_points_group)
        {
          it->second.mapped_pts_a.insert(it->second.mapped_pts_a.end(), pts_a.begin(),
                                         pts_a.end());
          if (kind == Kind::SURFACE_FLUX)
          {
            it->second.normal_scales.push_back(normal_scale);
          }
        }
        if (elem_a >= 0)
        {
          if (ghost_a)
          {
            it->second.face_nbr.push_back(plan.face_nbr);
            it->second.ghost_attr.push_back(plan.ghost_attr);
          }
          else
          {
            it->second.vol_indices_a.push_back(elem_a);
          }
        }
        if (elem_b >= 0)
        {
          if (ghost_b)
          {
            it->second.face_nbr.push_back(plan.face_nbr);
            it->second.ghost_attr.push_back(plan.ghost_attr);
          }
          else
          {
            it->second.vol_indices_b.push_back(elem_b);
          }
        }
        it->second.out_slots.push_back(out_slot);
      };

      if (kind == Kind::SURFACE_FLUX && can_surface_flux_at_points_a &&
          can_surface_flux_at_points_b)
      {
        // Split local two-sided flux into one-sided AtPoints operators so each side's
        // mapped volume reference coordinates are runtime data. The shared output slot
        // preserves the original two-sided difference or non-two-sided average.
        const double scale_a = flux_two_sided ? 1.0 : 0.5;
        const double scale_b = flux_two_sided ? -1.0 : 0.5;
        AddGroup(plan.elem_a, false, vol_geom_a, plan.pts_a, -1, false,
                 mfem::Geometry::INVALID, {}, true, 1.0, scale_a);
        AddGroup(plan.elem_b, false, vol_geom_b, plan.pts_b, -1, false,
                 mfem::Geometry::INVALID, {}, true, 1.0, scale_b);
      }
      else if (buffer_kind && can_at_points_a && can_at_points_b)
      {
        // Split two-sided boundary visualization into two one-sided AtPoints operators
        // that accumulate into the same output slot. One AtPoints operator can only own
        // one point-coordinate set; the two sides generally use different mapped
        // reference coordinates, so preserving the two-sided kernels would reintroduce
        // mapped-point specialization.
        const bool average_quantity =
            kind == Kind::BDR_FIELD_E || kind == Kind::BDR_FIELD_B ||
            kind == Kind::BDR_ENERGY_E || kind == Kind::BDR_ENERGY_M ||
            kind == Kind::BDR_POYNTING;
        const double side_weight = average_quantity ? 0.5 : 1.0;
        AddGroup(plan.elem_a, false, vol_geom_a, plan.pts_a, -1, false,
                 mfem::Geometry::INVALID, {}, true, side_weight, 1.0);
        AddGroup(plan.elem_b, false, vol_geom_b, plan.pts_b, -1, false,
                 mfem::Geometry::INVALID, {}, true, side_weight, -1.0);
      }
      else
      {
        const bool at_points_group = can_at_points_a && plan.elem_b < 0;
        AddGroup(plan.elem_a, plan.ghost_a, vol_geom_a, plan.pts_a, plan.elem_b,
                 plan.ghost_b, vol_geom_b, plan.pts_b, at_points_group, 1.0, 1.0);
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
    const int max_vsize = std::max(nd_fespace ? nd_fespace->GetVSize() : 0,
                                   rt_fespace ? rt_fespace->GetVSize() : 0);
    field_staging.SetSize(max_vsize);
    field_staging.UseDevice(true);
    field_staging = 0.0;
  }

  // Build the (group independent part of the) QFunction context for the integrand.
  std::vector<CeedIntScalar> base_ctx;
  if (kind == Kind::INTERFACE_EPR)
  {
    // CeedIntScalar is a union, so the runtime integrand selector (epr_type) needs its
    // own slot: [0].first = epr_type (0 = DEFAULT, 1 = MA, 2 = MS, 3 = SA), then
    // [1].second = scale0, [2].second = scale1, then (MS only) the material context. The
    // shared kernel passes ctx + 1 to the per-type helpers so their relative layout
    // (scale0, scale1, material) is unchanged.
    base_ctx.resize(3);
    base_ctx[1].second = 0.0;
    base_ctx[2].second = 0.0;
    switch (epr_type)
    {
      case InterfaceDielectric::DEFAULT:
        base_ctx[0].first = 0;
        base_ctx[1].second = 0.5 * epr_t * epr_epsilon;
        break;
      case InterfaceDielectric::MA:
        base_ctx[0].first = 1;
        base_ctx[1].second = 0.5 * epr_t / epr_epsilon;
        break;
      case InterfaceDielectric::MS:
        {
          base_ctx[0].first = 2;
          base_ctx[1].second = 0.5 * epr_t / epr_epsilon;
          MaterialPropertyCoefficient epsilon_func(mat_op->GetAttributeToMaterial(),
                                                   mat_op->GetPermittivityReal());
          auto mat_ctx = ceed::PopulateCoefficientContext(3, &epsilon_func);
          base_ctx.insert(base_ctx.end(), mat_ctx.begin(), mat_ctx.end());
        }
        break;
      case InterfaceDielectric::SA:
        base_ctx[0].first = 3;
        base_ctx[1].second = 0.5 * epr_t * epr_epsilon;
        base_ctx[2].second = 0.5 * epr_t / epr_epsilon;
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
  else if (kind == Kind::BDR_FLUX_Q || kind == Kind::BDR_CURRENT_J)
  {
    base_ctx.resize(2);
    base_ctx[0].second = 1.0;  // Normal sign, set per group
    base_ctx[1].second = viz_scaling;
    MaterialPropertyCoefficient coeff_func(mat_op->GetAttributeToMaterial(),
                                           (kind == Kind::BDR_FLUX_Q)
                                               ? mat_op->GetPermittivityReal()
                                               : mat_op->GetInvPermeability());
    auto mat_ctx = ceed::PopulateCoefficientContext(3, &coeff_func);
    base_ctx.insert(base_ctx.end(), mat_ctx.begin(), mat_ctx.end());
  }
  else if (kind == Kind::BDR_POYNTING)
  {
    // Pointwise S = scale * E x (mu^-1 B). For split two-sided AtPoints groups,
    // group.side_scale folds in the legacy 1/2 side average.
    base_ctx.resize(1);
    base_ctx[0].second = viz_scaling;
    MaterialPropertyCoefficient invmu_func(mat_op->GetAttributeToMaterial(),
                                           mat_op->GetInvPermeability());
    auto mat_ctx = ceed::PopulateCoefficientContext(3, &invmu_func);
    base_ctx.insert(base_ctx.end(), mat_ctx.begin(), mat_ctx.end());
  }
  else if (kind == Kind::BDR_ENERGY_E || kind == Kind::BDR_ENERGY_M)
  {
    // Shared energy kernel: [0].first = piola (0 = ND/E, 1 = RT/B), [1].second = scaling,
    // material table at +2.
    base_ctx.resize(2);
    base_ctx[0].first = (kind == Kind::BDR_ENERGY_M) ? 1 : 0;
    base_ctx[1].second = viz_scaling;
    MaterialPropertyCoefficient coeff_func(mat_op->GetAttributeToMaterial(),
                                           (kind == Kind::BDR_ENERGY_E)
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
    base_ctx[1].second = farfield_omega.real();
    base_ctx[2].second = farfield_omega.imag();
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
    base_ctx[1].second = 1.0;
    if (kind == Kind::BDR_FIELD_E || kind == Kind::BDR_FIELD_B)
    {
      // Shared field kernel: [0].first = piola (0 = ND/E, 1 = RT/B),
      // [1].second = output scale (used when split AtPoints sides accumulate).
      base_ctx[0].first = (kind == Kind::BDR_FIELD_B) ? 1 : 0;
    }
  }

  // Build the face neighbor field exchange for any ghost (face neighbor) sides of
  // two-sided interior boundaries on parallel interfaces. Collective: all processes
  // participate (those without ghost faces pose no requests), so the decision is reduced
  // globally. Each ghost face requests the neighbor's volume field at the mapped face
  // quadrature points; the returned physical-space values feed the ghost side below.
  {
    int any_ghost = 0;
    for (const auto &[k, g] : face_groups)
    {
      if (g.ghost_a || g.ghost_b)
      {
        any_ghost = 1;
        break;
      }
    }
    Mpi::GlobalSum(1, &any_ghost, comm);
    if (any_ghost > 0)
    {
      // Source slots match the field inputs: SURFACE_FLUX uses slot 0 = E (nd_fespace)
      // and slot 1 = B (rt_fespace); all other kinds carry a single field at slot 0.
      std::array<const mfem::ParFiniteElementSpace *, FaceNbrFieldExchange::MaxSources>
          ex_fes = {nullptr, nullptr, nullptr, nullptr};
      unsigned int source_mask;
      if (kind == Kind::SURFACE_FLUX)
      {
        ex_fes[0] = nd_fespace;
        ex_fes[1] = rt_fespace;
        source_mask = (flux_type == SurfaceFlux::ELECTRIC)   ? 0b01u
                      : (flux_type == SurfaceFlux::MAGNETIC) ? 0b10u
                                                             : 0b11u;
      }
      else if (kind == Kind::BDR_POYNTING)
      {
        ex_fes[0] = nd_fespace;
        ex_fes[1] = rt_fespace;
        source_mask = 0b11u;
      }
      else
      {
        ex_fes[0] = rt_fespace ? rt_fespace : nd_fespace;
        source_mask = 0b01u;
      }
      std::vector<FaceNbrFieldExchange::Request> requests;
      for (auto &[k, g] : face_groups)
      {
        if (!(g.ghost_a || g.ghost_b))
        {
          continue;
        }
        const mfem::IntegrationRule &ir = g.ghost_a ? *g.mapped_ir_a : *g.mapped_ir_b;
        std::vector<mfem::IntegrationPoint> pts(ir.GetNPoints());
        for (int q = 0; q < ir.GetNPoints(); q++)
        {
          pts[q] = ir.IntPoint(q);
        }
        g.req_idx.resize(g.face_nbr.size());
        for (std::size_t e = 0; e < g.face_nbr.size(); e++)
        {
          g.req_idx[e] = static_cast<int>(requests.size());
          auto &req = requests.emplace_back();
          req.face_nbr_elem = g.face_nbr[e];
          req.source_mask = source_mask;
          req.pts = pts;
        }
      }
      face_nbr_exchange = std::make_unique<FaceNbrFieldExchange>(mesh, ex_fes, requests);
    }
  }

  // Assemble a libCEED operator for each group. For now, all operators are constructed
  // on a single Ceed context (no OpenMP parallel assembly or application; correctness
  // first, this can be extended with the thread partitioning of fem/mesh.cpp later).
  const bool profile = (std::getenv("PALACE_SURFACE_PROFILE") != nullptr);
  long long profile_counts[4] = {0, 0, 0, 0};  // groups, elems, AtPoints groups, AtPoints elems
  for (auto &face_group : face_groups)
  {
    auto &group = face_group.second;
    const std::size_t num_elem = group.bdr_indices.size();
    if (profile)
    {
      profile_counts[0]++;
      profile_counts[1] += static_cast<long long>(num_elem);
      if (group.at_points)
      {
        profile_counts[2]++;
        profile_counts[3] += static_cast<long long>(num_elem);
      }
    }
    const bool has_b = !group.vol_indices_b.empty() || group.ghost_b;
    const bool buffer_kind = IsBufferKind(kind);
    const mfem::IntegrationRule &face_ir =
        buffer_kind
            ? mfem::GlobGeometryRefiner.Refine(group.bdr_geom, viz_lod, 1)->RefPts
            : mfem::IntRules.Get(group.bdr_geom,
                                 fem::DefaultIntegrationOrder::Get(pmesh, group.bdr_geom));

    // Objects are owned by the assembled operator; scratch destroys our references.
    CeedAssemblyScratch scratch(ceed);
    CeedElemRestriction points_restr = nullptr;
    CeedVector points_vec = nullptr;
    if (group.at_points)
    {
      const int nq = face_ir.GetNPoints();
      MFEM_VERIFY(group.mapped_pts_a.size() == num_elem * static_cast<std::size_t>(nq),
                  "Invalid AtPoints coordinate data for surface output group!");
      auto &points = elem_attrs.emplace_back(3 * num_elem * nq);
      for (std::size_t e = 0; e < num_elem; e++)
      {
        for (int q = 0; q < nq; q++)
        {
          const auto &ip = group.mapped_pts_a[e * nq + q];
          const std::size_t off = 3 * (e * nq + q);
          points[off + 0] = ip.x;
          points[off + 1] = ip.y;
          points[off + 2] = ip.z;
        }
      }
      CreateSequentialPointRestriction(ceed, num_elem, nq, 3, points.Size(), &points_restr);
      ceed::InitCeedVector(points, ceed, &points_vec);
      scratch.restrs.push_back(points_restr);
      scratch.vecs.push_back(points_vec);
    }

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
    auto AddSequentialPointInput = [&](const std::string &name, Vector &data,
                                       int num_comp)
    {
      const int nq = face_ir.GetNPoints();
      CeedElemRestriction restr;
      CeedVector vec;
      CreateSequentialPointRestriction(ceed, num_elem, nq, num_comp, data.Size(), &restr);
      ceed::InitCeedVector(data, ceed, &vec);
      inputs.push_back({name, vec, restr, nullptr, ceed::EvalMode::None});
      scratch.vecs.push_back(vec);
      scratch.restrs.push_back(restr);
    };
    auto AddFaceGeomInputAtPoints = [&]()
    {
      const int nq = face_ir.GetNPoints();
      if (!buffer_kind)
      {
        auto &qw = elem_attrs.emplace_back(num_elem * nq);
        for (std::size_t e = 0; e < num_elem; e++)
        {
          for (int q = 0; q < nq; q++)
          {
            qw[e * nq + q] = face_ir.IntPoint(q).weight;
          }
        }
        AddSequentialPointInput("qw", qw, 1);
      }
      if (kind == Kind::SURFACE_FLUX && !group.normal_scales.empty())
      {
        MFEM_VERIFY(group.normal_scales.size() == num_elem,
                    "SURFACE_FLUX AtPoints normal scales must be entry-indexed!");
      }
      auto &face_geom = elem_attrs.emplace_back(num_elem * nq * 6);
      face_geom = 0.0;
      for (std::size_t e = 0; e < num_elem; e++)
      {
        mfem::IsoparametricTransformation T;
        pmesh.GetBdrElementTransformation(group.bdr_indices[e], &T);
        for (int q = 0; q < nq; q++)
        {
          const mfem::IntegrationPoint &ip = face_ir.IntPoint(q);
          T.SetIntPoint(&ip);
          const mfem::DenseMatrix &J = T.Jacobian();
          const double normal_scale =
              (kind == Kind::SURFACE_FLUX && !group.normal_scales.empty())
                  ? group.normal_scales[e]
                  : 1.0;
          const std::size_t off = 6 * (e * nq + q);
          FoldNormalScaleIntoFaceJacobian(J, normal_scale,
                                          face_geom.HostReadWrite() + off);
        }
      }
      AddSequentialPointInput("grad_x_f", face_geom, 6);
    };
    const bool field_kind = (kind == Kind::BDR_FIELD_E || kind == Kind::BDR_FIELD_B);
    if (!field_kind)
    {
      if (group.at_points)
      {
        AddFaceGeomInputAtPoints();
      }
      else
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
          inputs.push_back(
              {"qw", nullptr, nullptr, face_mesh_basis, ceed::EvalMode::Weight});
        }
        inputs.push_back({"x_f", mesh_nodes_vec, face_mesh_restr, face_mesh_basis,
                          ceed::EvalMode::Grad});
        scratch.vecs.push_back(mesh_nodes_vec);
        scratch.restrs.push_back(face_mesh_restr);
        scratch.bases.push_back(face_mesh_basis);
      }
    }
    auto AddVolGeomInputs = [&](const std::string &suffix, const std::vector<int> &indices,
                                mfem::Geometry::Type geom, const mfem::IntegrationRule &ir)
    {
      const int num_pts = group.at_points ? face_ir.GetNPoints() : ir.GetNPoints();
      const auto &loc_attr = mesh.GetCeedAttributes();
      if (group.at_points)
      {
        auto &elem_attr_pts = elem_attrs.emplace_back(indices.size() * num_pts);
        for (std::size_t k = 0; k < indices.size(); k++)
        {
          const double attr = loc_attr.at(pmesh.GetAttribute(indices[k]));
          for (int q = 0; q < num_pts; q++)
          {
            elem_attr_pts[k * num_pts + q] = attr;
          }
        }
        AddSequentialPointInput("attr_" + suffix, elem_attr_pts, 1);

        CeedElemRestriction mesh_restr = FiniteElementSpace::BuildCeedElemRestriction(
            mesh_fespace, ceed, geom, indices, /*is_interp*/ true);
        const mfem::FiniteElement *mesh_fe =
            mesh_fespace.FEColl()->FiniteElementForGeometry(geom);
        CeedBasis mesh_basis;
        InitTetBasisForAtPoints(*mesh_fe, /*grad_only*/ true, mesh_fespace.GetVDim(), ceed,
                                &mesh_basis);
        CeedVector mesh_nodes_vec;
        ceed::InitCeedVector(*mesh_fespace.GetMesh()->GetNodes(), ceed, &mesh_nodes_vec);
        inputs.push_back(
            {"x_" + suffix, mesh_nodes_vec, mesh_restr, mesh_basis, ceed::EvalMode::Grad});
        scratch.vecs.push_back(mesh_nodes_vec);
        scratch.restrs.push_back(mesh_restr);
        scratch.bases.push_back(mesh_basis);
        return;
      }

      auto &elem_attr = elem_attrs.emplace_back(indices.size());
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
    // Ghost (face neighbor) side variant: the volume element lives on another process,
    // so its attributes come from group.ghost_attr and the Piola Jacobian is replaced by
    // a constant identity. The exchanged field values are already physical (the owning
    // process applied the Piola map), so Piola(I) = I in the shared kernels reproduces
    // them. Same input names/positions as AddVolGeomInputs ("grad_x_<suffix>" via
    // EVAL_NONE matches the EVAL_GRAD naming), so the kernels are unchanged.
    auto AddVolGeomInputsGhost = [&](const std::string &suffix,
                                     const mfem::IntegrationRule &ir)
    {
      const int num_pts = ir.GetNPoints();
      auto &elem_attr = elem_attrs.emplace_back(group.ghost_attr.size());
      const auto &loc_attr = mesh.GetCeedAttributes();
      for (std::size_t k = 0; k < group.ghost_attr.size(); k++)
      {
        elem_attr[k] = loc_attr.at(group.ghost_attr[k]);
      }
      CeedElemRestriction attr_restr;
      CeedBasis attr_basis;
      CeedVector attr_vec;
      PalaceCeedCall(ceed, CeedElemRestrictionCreateStrided(
                               ceed, group.ghost_attr.size(), 1, 1,
                               group.ghost_attr.size(), CEED_STRIDES_BACKEND, &attr_restr));
      {
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
      // Constant 3x3 identity Jacobian (9 components, component-major [elem][comp][pt]),
      // passed directly to the kernel's grad_x_<suffix> input (EVAL_NONE).
      auto &ident = elem_attrs.emplace_back(num_elem * 9 * num_pts);
      ident = 0.0;
      for (std::size_t e = 0; e < num_elem; e++)
      {
        for (int c : {0, 4, 8})
        {
          for (int i = 0; i < num_pts; i++)
          {
            ident[e * 9 * num_pts + c * num_pts + i] = 1.0;
          }
        }
      }
      CeedElemRestriction ident_restr;
      const CeedInt strides[3] = {1, num_pts, 9 * num_pts};
      PalaceCeedCall(ceed, CeedElemRestrictionCreateStrided(
                               ceed, static_cast<CeedInt>(num_elem), num_pts, 9,
                               (CeedSize)num_elem * 9 * num_pts, strides, &ident_restr));
      CeedVector ident_vec;
      ceed::InitCeedVector(ident, ceed, &ident_vec);
      inputs.push_back(
          {"grad_x_" + suffix, ident_vec, ident_restr, nullptr, ceed::EvalMode::None});
      scratch.vecs.push_back(attr_vec);
      scratch.vecs.push_back(ident_vec);
      scratch.restrs.push_back(attr_restr);
      scratch.restrs.push_back(ident_restr);
      scratch.bases.push_back(attr_basis);
    };
    if (group.ghost_a)
    {
      AddVolGeomInputsGhost("1", *group.mapped_ir_a);
    }
    else if (!group.vol_indices_a.empty())
    {
      AddVolGeomInputs("1", group.vol_indices_a, group.vol_geom_a,
                       group.at_points ? face_ir : *group.mapped_ir_a);
    }
    if (group.ghost_b)
    {
      AddVolGeomInputsGhost("2", *group.mapped_ir_b);
    }
    else if (!group.vol_indices_b.empty())
    {
      AddVolGeomInputs("2", group.vol_indices_b, group.vol_geom_b,
                       group.at_points ? face_ir : *group.mapped_ir_b);
    }
    if (kind == Kind::SURFACE_FLUX || kind == Kind::FARFIELD)
    {
      if (group.at_points)
      {
        // Physical coordinates are surface-only data. Keep them on the boundary
        // quadrature rule rather than evaluating the boundary basis at the volume
        // reference points used by the AtPoints trace operator. Two-sided surface flux
        // kernels never use x (their normal orientation is fixed by side convention), so
        // avoid the host boundary-coordinate transform for those large trace surfaces.
        const int nq = face_ir.GetNPoints();
        auto &x_pts = elem_attrs.emplace_back(num_elem * nq * 3);
        x_pts = 0.0;
        if (!(kind == Kind::SURFACE_FLUX && flux_two_sided))
        {
          for (std::size_t e = 0; e < num_elem; e++)
          {
            mfem::IsoparametricTransformation T;
            pmesh.GetBdrElementTransformation(group.bdr_indices[e], &T);
            for (int q = 0; q < nq; q++)
            {
              const mfem::IntegrationPoint &ip = face_ir.IntPoint(q);
              mfem::Vector x_loc(3);
              T.Transform(ip, x_loc);
              const std::size_t off = 3 * (e * nq + q);
              x_pts[off + 0] = x_loc(0);
              x_pts[off + 1] = x_loc(1);
              x_pts[off + 2] = x_loc(2);
            }
          }
        }
        AddSequentialPointInput("x", x_pts, 3);
      }
      else
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
      if (group.at_points)
      {
        InitTetBasisForAtPoints(*fe, /*grad_only*/ false, fespace.GetVDim(), ceed, &basis);
      }
      else
      {
        ceed::InitBasisAtPoints(*fe, ir, fespace.GetVDim(), ceed, &basis);
      }
      ceed::InitCeedVector(field_staging, ceed, &vec);
      inputs.push_back({name, vec, restr, basis, ceed::EvalMode::Interp});
      field_sources.emplace_back(name, source);
      scratch.vecs.push_back(vec);
      scratch.restrs.push_back(restr);
      scratch.bases.push_back(basis);
    };
    // Ghost (face neighbor) field input: the neighbor's physical field values at the
    // mapped points, imported via FaceNbrFieldExchange. The indexed restriction slices
    // the shared imported vector per element (base = ImportOffset) and transposes its
    // point-major layout ([pt][comp]) to the component-major layout the kernel reads
    // (offset base + 3*pt, comp_stride 1). Re-pointed at the imported vector on each
    // apply via source index 4 (see fem::ApplyAddGroupOperators).
    auto AddFieldInputGhost = [&](const std::string &name, int slot,
                                  const mfem::IntegrationRule &ir)
    {
      const int nq = ir.GetNPoints();
      std::vector<CeedInt> offsets(num_elem * nq);
      for (std::size_t e = 0; e < num_elem; e++)
      {
        const int base = face_nbr_exchange->ImportOffset(group.req_idx[e], slot);
        for (int i = 0; i < nq; i++)
        {
          offsets[e * nq + i] = base + 3 * i;
        }
      }
      CeedElemRestriction restr;
      PalaceCeedCall(ceed, CeedElemRestrictionCreate(
                               ceed, static_cast<CeedInt>(num_elem), nq, 3, 1,
                               (CeedSize)face_nbr_exchange->ImportSize(), CEED_MEM_HOST,
                               CEED_COPY_VALUES, offsets.data(), &restr));
      CeedVector vec;
      ceed::InitCeedVector(face_nbr_exchange->Imported(), ceed, &vec);
      inputs.push_back({name, vec, restr, nullptr, ceed::EvalMode::None});
      field_sources.emplace_back(name, 4);  // 4 -> imported (see ApplyAddGroupOperators)
      scratch.vecs.push_back(vec);
      scratch.restrs.push_back(restr);
    };
    if (kind == Kind::BDR_POYNTING)
    {
      auto AddPoyntingSide = [&](const std::string &suffix,
                                 const std::vector<int> &indices,
                                 mfem::Geometry::Type geom,
                                 const mfem::IntegrationRule &ir, bool ghost)
      {
        const std::string e_name = "u_e_" + suffix;
        const std::string b_name = "u_b_" + suffix;
        if (ghost)
        {
          AddFieldInputGhost(e_name, 0, ir);
          AddFieldInputGhost(b_name, 1, ir);
        }
        else
        {
          AddFieldInput(e_name, 0, *nd_fespace, indices, geom, ir);
          AddFieldInput(b_name, 1, *rt_fespace, indices, geom, ir);
        }
      };
      AddPoyntingSide("1", group.vol_indices_a, group.vol_geom_a,
                      group.at_points ? face_ir : *group.mapped_ir_a, group.ghost_a);
      if (has_b)
      {
        AddPoyntingSide("2", group.vol_indices_b, group.vol_geom_b,
                        group.at_points ? face_ir : *group.mapped_ir_b, group.ghost_b);
      }
    }
    else if (kind == Kind::HCURL_NORM2 || kind == Kind::INTERFACE_EPR || buffer_kind)
    {
      const auto &field_fespace = rt_fespace ? *rt_fespace : *nd_fespace;
      if (group.ghost_a)
      {
        AddFieldInputGhost("u_1", 0, *group.mapped_ir_a);
      }
      else
      {
        AddFieldInput("u_1", 0, field_fespace, group.vol_indices_a, group.vol_geom_a,
                      group.at_points ? face_ir : *group.mapped_ir_a);
      }
      if (has_b)
      {
        if (group.ghost_b)
        {
          AddFieldInputGhost("u_2", 0, *group.mapped_ir_b);
        }
        else
        {
          AddFieldInput("u_2", 0, field_fespace, group.vol_indices_b, group.vol_geom_b,
                        group.at_points ? face_ir : *group.mapped_ir_b);
        }
      }
    }
    else if (kind == Kind::FARFIELD)
    {
      const mfem::IntegrationRule &ir = group.at_points ? face_ir : *group.mapped_ir_a;
      AddFieldInput("u_1", 0, *nd_fespace, group.vol_indices_a, group.vol_geom_a, ir);
      AddFieldInput("u_2", 1, *nd_fespace, group.vol_indices_a, group.vol_geom_a, ir);
      AddFieldInput("u_3", 2, *rt_fespace, group.vol_indices_a, group.vol_geom_a, ir);
      AddFieldInput("u_4", 3, *rt_fespace, group.vol_indices_a, group.vol_geom_a, ir);
    }
    else if (kind == Kind::SURFACE_FLUX)
    {
      int count = 0;
      auto AddSide = [&](const std::vector<int> &indices, mfem::Geometry::Type geom,
                         const mfem::IntegrationRule &ir, bool ghost)
      {
        if (flux_type == SurfaceFlux::ELECTRIC || flux_type == SurfaceFlux::POWER)
        {
          const std::string nm = "u_" + std::to_string(++count);
          if (ghost)
          {
            AddFieldInputGhost(nm, 0, ir);
          }
          else
          {
            AddFieldInput(nm, 0, *nd_fespace, indices, geom, ir);
          }
        }
        if (flux_type == SurfaceFlux::MAGNETIC || flux_type == SurfaceFlux::POWER)
        {
          const std::string nm = "u_" + std::to_string(++count);
          if (ghost)
          {
            AddFieldInputGhost(nm, 1, ir);
          }
          else
          {
            AddFieldInput(nm, 1, *rt_fespace, indices, geom, ir);
          }
        }
      };
      AddSide(group.vol_indices_a, group.vol_geom_a,
              group.at_points ? face_ir : *group.mapped_ir_a, group.ghost_a);
      if (has_b)
      {
        AddSide(group.vol_indices_b, group.vol_geom_b,
                group.at_points ? face_ir : *group.mapped_ir_b, group.ghost_b);
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
      // Even for AtPoints operators, keep the output as an ordinary EVAL_NONE
      // restriction. libCEED requires all AtPoints restrictions on the same operator to
      // use identical point-offset layouts; the output buffer is intentionally scattered
      // by boundary-element slot and should not constrain the point-coordinate layout.
      PalaceCeedCall(ceed, CeedElemRestrictionCreate(ceed, static_cast<CeedInt>(num_elem),
                                                     nq, nc, 1, (CeedSize)buffer_size,
                                                     CEED_MEM_HOST, CEED_COPY_VALUES,
                                                     offsets.data(), &out_restr));
    }
    else if (group.at_points)
    {
      const int nq = face_ir.GetNPoints();
      std::vector<CeedInt> offsets(num_elem * nq);
      for (std::size_t e = 0; e < num_elem; e++)
      {
        for (int j = 0; j < nq; j++)
        {
          offsets[e * nq + j] = group.out_slots[e];
        }
      }
      // Each quadrature-point contribution scatters to the element's output slot;
      // ApplyAdd performs the reduction over the duplicate offsets.
      PalaceCeedCall(ceed, CeedElemRestrictionCreate(
                               ceed, static_cast<CeedInt>(num_elem), nq, num_out,
                               num_marked, (CeedSize)num_marked * num_out,
                               CEED_MEM_HOST, CEED_COPY_VALUES, offsets.data(),
                               &out_restr));
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
        // All four interface types share one kernel; epr_type is set in base_ctx[0].first.
        info.apply_qf = has_b ? f_integ_surf_epr_2_32 : f_integ_surf_epr_1_32;
        info.apply_qf_path = PalaceQFunctionRelativePath(has_b ? f_integ_surf_epr_2_32_loc
                                                               : f_integ_surf_epr_1_32_loc);
        break;
        break;
      case Kind::BDR_FIELD_E:
      case Kind::BDR_FIELD_B:
        // Shared kernel; the ND/RT Piola is set in base_ctx[0].first.
        ctx[1].second *= group.side_scale;
        info.apply_qf = has_b ? f_eval_bdr_field_2_32 : f_eval_bdr_field_1_32;
        info.apply_qf_path = PalaceQFunctionRelativePath(has_b ? f_eval_bdr_field_2_32_loc
                                                               : f_eval_bdr_field_1_32_loc);
        break;
      case Kind::BDR_FLUX_Q:
        ctx[0].second = (group.flip_normal ? -1.0 : 1.0) * group.normal_scale;
        info.apply_qf = has_b ? f_eval_bdr_flux_q_2_32 : f_eval_bdr_flux_q_1_32;
        info.apply_qf_path = PalaceQFunctionRelativePath(
            has_b ? f_eval_bdr_flux_q_2_32_loc : f_eval_bdr_flux_q_1_32_loc);
        break;
      case Kind::BDR_CURRENT_J:
        ctx[0].second = (group.flip_normal ? -1.0 : 1.0) * group.normal_scale;
        info.apply_qf = has_b ? f_eval_bdr_current_j_2_32 : f_eval_bdr_current_j_1_32;
        info.apply_qf_path = PalaceQFunctionRelativePath(
            has_b ? f_eval_bdr_current_j_2_32_loc : f_eval_bdr_current_j_1_32_loc);
        break;
      case Kind::BDR_ENERGY_E:
      case Kind::BDR_ENERGY_M:
        // Shared kernel; the ND/RT Piola is set in base_ctx[0].first.
        ctx[1].second *= group.side_scale;
        info.apply_qf = has_b ? f_eval_bdr_energy_2_32 : f_eval_bdr_energy_1_32;
        info.apply_qf_path = PalaceQFunctionRelativePath(
            has_b ? f_eval_bdr_energy_2_32_loc : f_eval_bdr_energy_1_32_loc);
        break;
      case Kind::BDR_POYNTING:
        ctx[0].second *= group.side_scale;
        info.apply_qf = has_b ? f_eval_bdr_poynting_2_32 : f_eval_bdr_poynting_1_32;
        info.apply_qf_path = PalaceQFunctionRelativePath(
            has_b ? f_eval_bdr_poynting_2_32_loc : f_eval_bdr_poynting_1_32_loc);
        break;
      case Kind::FARFIELD:
        ctx[0].second = group.flip_normal ? -1.0 : 1.0;
        info.apply_qf = f_integ_surf_farfield_32;
        info.apply_qf_path = PalaceQFunctionRelativePath(f_integ_surf_farfield_32_loc);
        break;
      case Kind::SURFACE_FLUX:
        ctx[0].second =
            (group.flip_normal ? -1.0 : 1.0) *
            ((group.at_points && !group.normal_scales.empty())
                 ? 1.0  // Per-entry normal_scale is folded into face geometry.
                 : group.normal_scale);
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
    CeedQFunctionContext op_ctx = nullptr;
    if (buffer_kind)
    {
      if (group.at_points)
      {
        ceed::AssembleCeedPointEvaluatorAtPoints(
            info, ctx.data(), ctx.size() * sizeof(CeedIntScalar), ceed, inputs,
            points_restr, points_vec, BufferNumComp(kind), out_restr, &op);
      }
      else
      {
        ceed::AssembleCeedPointEvaluator(info, ctx.data(),
                                         ctx.size() * sizeof(CeedIntScalar), ceed, inputs,
                                         BufferNumComp(kind), out_restr, &op);
      }
    }
    else if (group.at_points)
    {
      ceed::AssembleCeedPointEvaluatorAtPoints(
          info, ctx.data(), ctx.size() * sizeof(CeedIntScalar), ceed, inputs,
          points_restr, points_vec, num_out, out_restr, &op,
          (kind == Kind::FARFIELD) ? &op_ctx : nullptr);
    }
    else
    {
      ceed::AssembleCeedSurfaceFunctional(info, ctx.data(),
                                          ctx.size() * sizeof(CeedIntScalar), ceed, inputs,
                                          num_out, out_restr, &op, &op_ctx);
    }
    groups.push_back({ceed, op, std::move(field_sources), op_ctx});
  }
  if (profile)
  {
    Mpi::GlobalSum(4, profile_counts, comm);
    Mpi::Print(comm,
               "SurfaceFunctional profile kind={} groups={} elems={} at_points_groups={} "
               "at_points_elems={}\n",
               KindName(kind), profile_counts[0], profile_counts[1], profile_counts[2],
               profile_counts[3]);
  }
}

void SurfaceFunctional::ApplyAdd(const std::array<const Vector *, 4> &srcs) const
{
  if (face_nbr_exchange)
  {
    face_nbr_exchange->Exchange(srcs);
    fem::ApplyAddGroupOperators(groups, srcs, local_out, &face_nbr_exchange->Imported());
  }
  else
  {
    fem::ApplyAddGroupOperators(groups, srcs, local_out);
  }
}

double SurfaceFunctional::EvalLocal(const std::array<const Vector *, 4> &srcs) const
{
  MFEM_VERIFY(valid, "Eval called on an invalid (unassembled) SurfaceFunctional!");
  // A process holding the face neighbor exchange must still apply (the exchange is
  // collective: it may need to export its local field for a neighbor's ghost request
  // even with no local marked elements). ApplyAdd is a no-op on the empty group set.
  if (local_out.Size() == 0 && !face_nbr_exchange)
  {
    return 0.0;
  }
  if (local_out.Size() > 0)
  {
    local_out = 0.0;
  }
  ApplyAdd(srcs);
  return (local_out.Size() > 0) ? linalg::LocalSum(local_out) : 0.0;
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
  if (local_out.Size() > 0 || face_nbr_exchange)
  {
    if (local_out.Size() > 0)
    {
      local_out = 0.0;
    }
    ApplyAdd({&u.Real(), nullptr});
    if (u.HasImag())
    {
      ApplyAdd({&u.Imag(), nullptr});
    }
    if (local_out.Size() > 0)
    {
      dot = linalg::LocalSum(local_out);
    }
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
  MFEM_VERIFY(valid && IsBufferKind(kind) && kind != Kind::BDR_POYNTING,
              "EvalBuffer requires a valid single-field boundary visualization functional!");
  MFEM_ASSERT(buffer.Size() == buffer_size, "Invalid buffer size for EvalBuffer!");
  buffer = 0.0;
  if (face_nbr_exchange)
  {
    face_nbr_exchange->Exchange({&u});
    fem::ApplyAddGroupOperators(groups, {&u}, buffer, &face_nbr_exchange->Imported());
  }
  else
  {
    fem::ApplyAddGroupOperators(groups, {&u}, buffer);
  }
}

void SurfaceFunctional::EvalBuffer(const GridFunction &u, Vector &buffer) const
{
  MFEM_VERIFY(valid && IsBufferKind(kind) && kind != Kind::BDR_POYNTING,
              "EvalBuffer requires a valid single-field boundary visualization functional!");
  MFEM_ASSERT(buffer.Size() == buffer_size, "Invalid buffer size for EvalBuffer!");
  buffer = 0.0;
  auto Apply = [&](const Vector &v)
  {
    if (face_nbr_exchange)
    {
      face_nbr_exchange->Exchange({&v});
      fem::ApplyAddGroupOperators(groups, {&v}, buffer, &face_nbr_exchange->Imported());
    }
    else
    {
      fem::ApplyAddGroupOperators(groups, {&v}, buffer);
    }
  };
  Apply(u.Real());
  if (u.HasImag())
  {
    Apply(u.Imag());
  }
}

void SurfaceFunctional::EvalBuffer(const GridFunction &E, const GridFunction &B,
                                   Vector &buffer) const
{
  MFEM_VERIFY(valid && kind == Kind::BDR_POYNTING,
              "EvalBuffer requires a valid boundary Poynting visualization functional!");
  MFEM_VERIFY(E.HasImag() == B.HasImag(),
              "Mismatch between real- and complex-valued E and B fields in boundary "
              "Poynting visualization!");
  MFEM_ASSERT(buffer.Size() == buffer_size, "Invalid buffer size for EvalBuffer!");
  buffer = 0.0;
  auto Apply = [&](const Vector &e, const Vector &b)
  {
    if (face_nbr_exchange)
    {
      face_nbr_exchange->Exchange({&e, &b});
      fem::ApplyAddGroupOperators(groups, {&e, &b}, buffer,
                                  &face_nbr_exchange->Imported());
    }
    else
    {
      fem::ApplyAddGroupOperators(groups, {&e, &b}, buffer);
    }
  };
  Apply(E.Real(), B.Real());
  if (E.HasImag())
  {
    Apply(E.Imag(), B.Imag());
  }
}

std::vector<std::array<std::complex<double>, 3>>
SurfaceFunctional::EvalFarField(const GridFunction &E, const GridFunction &B,
                                std::complex<double> omega)
{
  MFEM_VERIFY(kind == Kind::FARFIELD && E.HasImag() && B.HasImag(),
              "SurfaceFunctional::EvalFarField requires a far-field functional and "
              "complex-valued fields!");
  MFEM_VERIFY(valid, "EvalFarField called on an invalid (unassembled) SurfaceFunctional!");

  // The frequency enters only the QFunction context (the omega slots); update it in
  // place rather than reassembling the operators. Reassembly would rebuild the bases,
  // restrictions, and on-the-fly geometry inputs and re-JIT the (expensive) far-field
  // kernel on every frequency -- none of which depend on omega. FARFIELD context layout
  // (see Assemble): [0] normal sign, [1] omega_re, [2] omega_im, [3] N, [4..] directions,
  // then the material context.
  if (omega != farfield_omega)
  {
    farfield_omega = omega;
    for (auto &group : groups)
    {
      if (!group.ctx)
      {
        continue;
      }
      CeedIntScalar *data;
      PalaceCeedCall(group.ceed,
                     CeedQFunctionContextGetData(group.ctx, CEED_MEM_HOST, &data));
      data[1].second = omega.real();
      data[2].second = omega.imag();
      PalaceCeedCall(group.ceed, CeedQFunctionContextRestoreData(group.ctx, &data));
    }
  }

  // Integrate, reduce each component over the local elements and all processes, and
  // apply the final cross products (following GetFarFieldrE).
  const int N = static_cast<int>(farfield_dirs.size());
  const int num_marked = local_out.Size() / std::max(6 * N, 1);
  std::vector<double> integrals(6 * N, 0.0);
  if (local_out.Size() > 0)
  {
    local_out = 0.0;
    fem::ApplyAddGroupOperators(groups, {&E.Real(), &E.Imag(), &B.Real(), &B.Imag()},
                                local_out);
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

std::vector<std::complex<double>> SurfaceFunctional::EvalComplexPowerByAttribute(
    const GridFunction &E, const GridFunction &B, const mfem::Array<int> &attr_to_bin,
    int num_bins) const
{
  MFEM_VERIFY(kind == Kind::SURFACE_FLUX && flux_type == SurfaceFlux::POWER &&
                  flux_two_sided,
              "SurfaceFunctional::EvalComplexPowerByAttribute is only valid for "
              "two-sided POWER flux functionals!");
  MFEM_VERIFY(E.HasImag() == B.HasImag(),
              "Mismatch between real- and complex-valued E and B fields in batched "
              "port power calculation!");
  MFEM_VERIFY(num_bins >= 0, "Invalid number of output bins!");
  MFEM_VERIFY(local_out_attrs.size() == static_cast<std::size_t>(local_out.Size()),
              "SurfaceFunctional attribute bins require one output slot per element!");

  auto AccumulateBins = [&](const std::array<const Vector *, 4> &srcs,
                            std::vector<double> &bins, double scale)
  {
    if (local_out.Size() > 0)
    {
      local_out = 0.0;
    }
    // Keep the apply collective even on ranks with no local marked elements, since a
    // face-neighbor exchange may need this rank to export field data for another rank's
    // processor-boundary side.
    ApplyAdd(srcs);
    if (local_out.Size() == 0)
    {
      return;
    }
    const double *vals = local_out.HostRead();
    for (int i = 0; i < local_out.Size(); i++)
    {
      const int attr = local_out_attrs[i];
      const int bin = (attr > 0 && attr <= attr_to_bin.Size()) ? attr_to_bin[attr - 1] : -1;
      if (bin >= 0)
      {
        MFEM_VERIFY(bin < num_bins, "SurfaceFunctional attribute bin out of range!");
        bins[bin] += scale * vals[i];
      }
    }
  };

  std::vector<double> real(num_bins, 0.0), imag(num_bins, 0.0);
  AccumulateBins({&E.Real(), &B.Real()}, real, 1.0);
  if (E.HasImag())
  {
    AccumulateBins({&E.Imag(), &B.Imag()}, real, 1.0);
    AccumulateBins({&E.Imag(), &B.Real()}, imag, 1.0);
    AccumulateBins({&E.Real(), &B.Imag()}, imag, -1.0);
  }

  std::vector<double> packed(2 * num_bins, 0.0);
  for (int i = 0; i < num_bins; i++)
  {
    packed[2 * i + 0] = real[i];
    packed[2 * i + 1] = imag[i];
  }
  Mpi::GlobalSum(static_cast<int>(packed.size()), packed.data(), comm);

  std::vector<std::complex<double>> result(num_bins);
  for (int i = 0; i < num_bins; i++)
  {
    result[i] = {packed[2 * i + 0], packed[2 * i + 1]};
  }
  return result;
}

}  // namespace palace
