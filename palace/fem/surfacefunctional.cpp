// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfacefunctional.hpp"

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

// Build the libCEED geometry factor quadrature data for the given elements with the
// given integration rule (mirrors the geometry data assembly in fem/mesh.cpp, but with
// caller-provided integration rules so that face and volume data are point-consistent).
// The attribute channel is filled with the provided per-element attribute values.
// Returns the geometry data vector and restriction (owned by the caller).
void BuildGeometryData(Ceed ceed, const mfem::FiniteElementSpace &mesh_fespace,
                       mfem::Geometry::Type geom, const std::vector<int> &indices,
                       const mfem::IntegrationRule &ir, const Vector &elem_attr,
                       CeedVector *geom_data, CeedElemRestriction *geom_data_restr)
{
  const std::size_t num_elem = indices.size();
  const int dim = mfem::Geometry::Dimension[geom];
  const int space_dim = mesh_fespace.GetMesh()->SpaceDimension();

  // Construct mesh node element restriction and basis (basis tabulated at the provided
  // integration rule points, which may lie on a face of the reference element for
  // volume geometry data evaluated at face quadrature points). Native dof ordering is
  // required to pair with the full (non-tensor) basis tabulation.
  CeedElemRestriction mesh_restr = FiniteElementSpace::BuildCeedElemRestriction(
      mesh_fespace, ceed, geom, indices, /*is_interp*/ true);
  CeedBasis mesh_basis;
  {
    const mfem::FiniteElement *fe = mesh_fespace.FEColl()->FiniteElementForGeometry(geom);
    if (!fe)
    {
      fe = mesh_fespace.FEColl()->TraceFiniteElementForGeometry(geom);
    }
    MFEM_VERIFY(fe, "Unable to get mesh nodal finite element for geometry data!");
    ceed::InitBasisAtPoints(*fe, ir, mesh_fespace.GetVDim(), ceed, &mesh_basis);
  }
  CeedVector mesh_nodes_vec;
  ceed::InitCeedVector(*mesh_fespace.GetMesh()->GetNodes(), ceed, &mesh_nodes_vec);
  CeedInt num_qpts;
  PalaceCeedCall(ceed, CeedBasisGetNumQuadraturePoints(mesh_basis, &num_qpts));

  // Construct element attribute element restriction and basis (all-ones basis
  // broadcasting the per-element attribute to quadrature points).
  CeedElemRestriction attr_restr;
  CeedBasis attr_basis;
  PalaceCeedCall(ceed, CeedElemRestrictionCreateStrided(ceed, num_elem, 1, 1, num_elem,
                                                        CEED_STRIDES_BACKEND, &attr_restr));
  {
    // Note: ceed::GetCeedTopology(CEED_TOPOLOGY_LINE) == 1.
    mfem::Vector Bt(num_qpts), Gt(num_qpts), qX(num_qpts), qW(num_qpts);
    Bt = 1.0;
    Gt = 0.0;
    qX = 0.0;
    qW = 0.0;
    PalaceCeedCall(ceed, CeedBasisCreateH1(ceed, CEED_TOPOLOGY_LINE, 1, 1, num_qpts,
                                           Bt.GetData(), Gt.GetData(), qX.GetData(),
                                           qW.GetData(), &attr_basis));
  }
  CeedVector elem_attr_vec;
  ceed::InitCeedVector(elem_attr, ceed, &elem_attr_vec);

  // Allocate storage for geometry factor data.
  CeedInt geom_data_size = 2 + space_dim * dim;
  PalaceCeedCall(
      ceed,
      CeedVectorCreate(ceed, (CeedSize)num_elem * num_qpts * geom_data_size, geom_data));
  PalaceCeedCall(
      ceed, CeedElemRestrictionCreateStrided(ceed, num_elem, num_qpts, geom_data_size,
                                             (CeedSize)num_elem * num_qpts * geom_data_size,
                                             CEED_STRIDES_BACKEND, geom_data_restr));

  // Compute the required geometry factors at quadrature points.
  ceed::AssembleCeedGeometryData(ceed, mesh_restr, mesh_basis, mesh_nodes_vec, attr_restr,
                                 attr_basis, elem_attr_vec, *geom_data, *geom_data_restr);
  PalaceCeedCall(ceed, CeedVectorDestroy(&mesh_nodes_vec));
  PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&mesh_restr));
  PalaceCeedCall(ceed, CeedBasisDestroy(&mesh_basis));
  PalaceCeedCall(ceed, CeedVectorDestroy(&elem_attr_vec));
  PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&attr_restr));
  PalaceCeedCall(ceed, CeedBasisDestroy(&attr_basis));
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

}  // namespace

SurfaceFunctional::SurfaceFunctional(Kind kind, const Mesh &mesh,
                                     const mfem::Array<int> &bdr_attr_marker,
                                     const FiniteElementSpace *fespace)
  : kind(kind), fespace_e(fespace), fespace_b(nullptr), mat_op(nullptr),
    comm(mesh.GetComm())
{
  MFEM_VERIFY(kind == Kind::AREA || kind == Kind::HCURL_NORM2,
              "Invalid SurfaceFunctional constructor for the requested functional kind!");
  MFEM_VERIFY(kind == Kind::AREA || fespace,
              "SurfaceFunctional requires a field finite element space for functionals "
              "with field inputs!");
  MFEM_VERIFY(mesh.Dimension() == 3 && mesh.SpaceDimension() == 3,
              "SurfaceFunctional is only implemented for 3D meshes!");
  Assemble(mesh, bdr_attr_marker);
}

SurfaceFunctional::SurfaceFunctional(const Mesh &mesh,
                                     const mfem::Array<int> &bdr_attr_marker,
                                     const FiniteElementSpace &nd_fespace,
                                     const MaterialOperator &mat_op,
                                     InterfaceDielectric type, double t_i, double epsilon_i)
  : kind(Kind::INTERFACE_EPR), epr_type(type), epr_t(t_i), epr_epsilon(epsilon_i),
    fespace_e(&nd_fespace), fespace_b(nullptr), mat_op(&mat_op), comm(mesh.GetComm())
{
  MFEM_VERIFY(mesh.Dimension() == 3 && mesh.SpaceDimension() == 3,
              "SurfaceFunctional is only implemented for 3D meshes!");
  Assemble(mesh, bdr_attr_marker);
}

SurfaceFunctional::SurfaceFunctional(const Mesh &mesh,
                                     const mfem::Array<int> &bdr_attr_marker,
                                     const FiniteElementSpace *nd_fespace,
                                     const FiniteElementSpace *rt_fespace,
                                     const MaterialOperator &mat_op, SurfaceFlux type,
                                     bool two_sided, const mfem::Vector &x0)
  : kind(Kind::SURFACE_FLUX), flux_type(type), flux_two_sided(two_sided), flux_x0(x0),
    fespace_e(nd_fespace), fespace_b(rt_fespace), mat_op(&mat_op), comm(mesh.GetComm())
{
  MFEM_VERIFY(
      (nd_fespace || (type != SurfaceFlux::ELECTRIC && type != SurfaceFlux::POWER)) &&
          (rt_fespace || (type != SurfaceFlux::MAGNETIC && type != SurfaceFlux::POWER)),
      "Missing finite element space for surface flux functional!");
  MFEM_VERIFY(mesh.Dimension() == 3 && mesh.SpaceDimension() == 3,
              "SurfaceFunctional is only implemented for 3D meshes!");
  Assemble(mesh, bdr_attr_marker);
}

SurfaceFunctional::~SurfaceFunctional()
{
  for (auto &group : groups)
  {
    PalaceCeedCall(group.ceed, CeedOperatorDestroy(&group.op));
  }
}

void SurfaceFunctional::Assemble(const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker)
{
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
      else if (kind == Kind::SURFACE_FLUX ||
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
      // face neighbor data vector.
      MFEM_VERIFY((plan.elem_a < 0 || plan.elem_a < pmesh.GetNE()) &&
                      (plan.elem_b < 0 || elem2_local),
                  "SurfaceFunctional does not yet support two-sided evaluation on "
                  "process-boundary interior surfaces (boundary element "
                      << i << ")!");

      // Map the face quadrature points to the volume element reference space(s):
      // boundary element reference coordinates -> face reference coordinates
      // (TransformBdrElementToFace with the boundary element to face orientation) ->
      // element 1/2 reference coordinates (FET.Loc1/Loc2).
      const auto bdr_geom = pmesh.GetBdrElementGeometry(i);
      const int q_order = fem::DefaultIntegrationOrder::Get(pmesh, bdr_geom);
      const mfem::IntegrationRule &face_ir = mfem::IntRules.Get(bdr_geom, q_order);
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
      it->second.out_slots.push_back(num_marked++);
    }
  }

  // Initialize the local output vector and field staging vector.
  local_out.SetSize(num_marked);
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
    const mfem::IntegrationRule &face_ir = mfem::IntRules.Get(
        group.bdr_geom, fem::DefaultIntegrationOrder::Get(pmesh, group.bdr_geom));

    // Handles to destroy after operator assembly (the operator holds references).
    std::vector<CeedVector> tmp_vecs;
    std::vector<CeedElemRestriction> tmp_restrs;
    std::vector<CeedBasis> tmp_bases;

    // Face geometry data (boundary element Jacobians at the face quadrature points).
    CeedVector face_geom_data;
    CeedElemRestriction face_geom_data_restr;
    {
      Vector elem_attr(num_elem);
      for (std::size_t k = 0; k < num_elem; k++)
      {
        elem_attr[k] = pmesh.GetBdrAttribute(group.bdr_indices[k]);
      }
      BuildGeometryData(ceed, mesh_fespace, group.bdr_geom, group.bdr_indices, face_ir,
                        elem_attr, &face_geom_data, &face_geom_data_restr);
      tmp_vecs.push_back(face_geom_data);
      tmp_restrs.push_back(face_geom_data_restr);
    }

    // Volume geometry data evaluated at the mapped face quadrature points (for the
    // Piola transformations of the field inputs and the material property lookups with
    // the local libCEED attribute).
    auto GetCeedElemAttr = [&](const std::vector<int> &indices)
    {
      Vector elem_attr(indices.size());
      const auto &loc_attr = mesh.GetCeedAttributes();
      for (std::size_t k = 0; k < indices.size(); k++)
      {
        elem_attr[k] = loc_attr.at(pmesh.GetAttribute(indices[k]));
      }
      return elem_attr;
    };
    CeedVector vol_geom_data = nullptr, vol_geom_data_b = nullptr;
    CeedElemRestriction vol_geom_data_restr = nullptr, vol_geom_data_b_restr = nullptr;
    if (!group.vol_indices_a.empty())
    {
      Vector elem_attr = GetCeedElemAttr(group.vol_indices_a);
      BuildGeometryData(ceed, mesh_fespace, group.vol_geom_a, group.vol_indices_a,
                        *group.mapped_ir_a, elem_attr, &vol_geom_data,
                        &vol_geom_data_restr);
      tmp_vecs.push_back(vol_geom_data);
      tmp_restrs.push_back(vol_geom_data_restr);
    }
    if (has_b)
    {
      Vector elem_attr = GetCeedElemAttr(group.vol_indices_b);
      BuildGeometryData(ceed, mesh_fespace, group.vol_geom_b, group.vol_indices_b,
                        *group.mapped_ir_b, elem_attr, &vol_geom_data_b,
                        &vol_geom_data_b_restr);
      tmp_vecs.push_back(vol_geom_data_b);
      tmp_restrs.push_back(vol_geom_data_b_restr);
    }

    // Assemble the inputs in the order expected by the QFunctions: side b volume
    // geometry data (EvalMode::None), coordinates (SURFACE_FLUX only), then the field
    // inputs (side a fields, then side b fields).
    std::vector<ceed::CeedFunctionalFieldInput> inputs;
    std::vector<std::pair<std::string, int>> field_sources;
    if (has_b)
    {
      inputs.push_back({"vol_geom_data_b", vol_geom_data_b, vol_geom_data_b_restr, nullptr,
                        ceed::EvalMode::None});
    }
    if (kind == Kind::SURFACE_FLUX)
    {
      // Coordinates at the face quadrature points, interpolated from the mesh nodes
      // with the boundary element basis (constant input, never re-pointed).
      CeedElemRestriction x_restr = FiniteElementSpace::BuildCeedElemRestriction(
          mesh_fespace, ceed, group.bdr_geom, group.bdr_indices, /*is_interp*/ true);
      const mfem::FiniteElement *fe =
          mesh_fespace.FEColl()->FiniteElementForGeometry(group.bdr_geom);
      if (!fe)
      {
        fe = mesh_fespace.FEColl()->TraceFiniteElementForGeometry(group.bdr_geom);
      }
      CeedBasis x_basis;
      ceed::InitBasisAtPoints(*fe, face_ir, mesh_fespace.GetVDim(), ceed, &x_basis);
      CeedVector x_vec;
      ceed::InitCeedVector(*mesh_fespace.GetMesh()->GetNodes(), ceed, &x_vec);
      inputs.push_back({"x", x_vec, x_restr, x_basis, ceed::EvalMode::Interp});
      tmp_vecs.push_back(x_vec);
      tmp_restrs.push_back(x_restr);
      tmp_bases.push_back(x_basis);
    }
    auto AddFieldInput = [&](const std::string &name, int source,
                             const FiniteElementSpace &fespace,
                             const std::vector<int> &indices, mfem::Geometry::Type geom,
                             const mfem::IntegrationRule &ir)
    {
      CeedElemRestriction restr;
      CeedBasis basis;
      CeedVector vec;
      ceed::InitRestriction(fespace.Get(), indices, false, /*is_interp*/ true, false, ceed,
                            &restr);
      const mfem::FiniteElement *fe = fespace.GetFEColl().FiniteElementForGeometry(geom);
      MFEM_VERIFY(fe, "Unable to get field finite element for surface functional!");
      ceed::InitBasisAtPoints(*fe, ir, fespace.GetVDim(), ceed, &basis);
      ceed::InitCeedVector(field_staging, ceed, &vec);
      inputs.push_back({name, vec, restr, basis, ceed::EvalMode::Interp});
      field_sources.emplace_back(name, source);
      tmp_vecs.push_back(vec);
      tmp_restrs.push_back(restr);
      tmp_bases.push_back(basis);
    };
    if (kind == Kind::HCURL_NORM2 || kind == Kind::INTERFACE_EPR)
    {
      AddFieldInput("u_1", 0, *fespace_e, group.vol_indices_a, group.vol_geom_a,
                    *group.mapped_ir_a);
      if (has_b)
      {
        AddFieldInput("u_2", 0, *fespace_e, group.vol_indices_b, group.vol_geom_b,
                      *group.mapped_ir_b);
      }
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

    // Output restriction: one slot per boundary element in the local output vector.
    CeedElemRestriction out_restr;
    PalaceCeedCall(ceed, CeedElemRestrictionCreate(ceed, static_cast<CeedInt>(num_elem), 1,
                                                   1, num_marked, num_marked, CEED_MEM_HOST,
                                                   CEED_COPY_VALUES, group.out_slots.data(),
                                                   &out_restr));
    tmp_restrs.push_back(out_restr);

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
    ceed::AssembleCeedSurfaceFunctional(
        info, ctx.data(), ctx.size() * sizeof(CeedIntScalar), ceed, inputs, face_geom_data,
        face_geom_data_restr, vol_geom_data, vol_geom_data_restr, 1, out_restr, &op);
    groups.push_back({ceed, op, std::move(field_sources)});

    // Cleanup (objects are now owned through the operator's references).
    for (auto &v : tmp_vecs)
    {
      PalaceCeedCall(ceed, CeedVectorDestroy(&v));
    }
    for (auto &r : tmp_restrs)
    {
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&r));
    }
    for (auto &b : tmp_bases)
    {
      PalaceCeedCall(ceed, CeedBasisDestroy(&b));
    }
  }
}

void SurfaceFunctional::ApplyAdd(const std::array<const Vector *, 2> &srcs) const
{
  for (const auto &group : groups)
  {
    Ceed ceed = group.ceed;

    // Re-point the passive field inputs at the caller's data (the data or its location
    // may have changed since the last call).
    for (const auto &[name, source] : group.field_sources)
    {
      MFEM_ASSERT(srcs[source], "Missing source vector for SurfaceFunctional field input!");
      CeedOperatorField field;
      CeedVector field_vec;
      PalaceCeedCall(ceed, CeedOperatorGetFieldByName(group.op, name.c_str(), &field));
      PalaceCeedCall(ceed, CeedOperatorFieldGetVector(field, &field_vec));
      ceed::InitCeedVector(*srcs[source], ceed, &field_vec, false);
    }

    // Accumulate the per-element integrals into the local output vector.
    CeedVector out_vec;
    ceed::InitCeedVector(local_out, ceed, &out_vec);
    PalaceCeedCall(ceed, CeedOperatorApplyAdd(group.op, CEED_VECTOR_NONE, out_vec,
                                              CEED_REQUEST_IMMEDIATE));
    PalaceCeedCall(ceed, CeedVectorDestroy(&out_vec));
  }
}

double SurfaceFunctional::EvalLocal(const std::array<const Vector *, 2> &srcs) const
{
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

}  // namespace palace
