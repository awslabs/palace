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
// volume element reference space(s). elem_b >= 0 indicates two-sided evaluation with
// averaging of the fields from both sides.
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
  : kind(kind), fespace(fespace), mat_op(nullptr), comm(mesh.GetComm())
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
    fespace(&nd_fespace), mat_op(&mat_op), comm(mesh.GetComm())
{
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
      // the legacy coefficients (see InterfaceDielectricCoefficient).
      ElemPlan plan;
      plan.bdr_elem = i;
      plan.flip = flip;
      if (!need_field)
      {
        // No field inputs (side selection does not apply).
      }
      else if (kind == Kind::HCURL_NORM2 ||
               (kind == Kind::INTERFACE_EPR && epr_type == InterfaceDielectric::DEFAULT &&
                !has_elem2))
      {
        plan.elem_a = FET.Elem1No;
      }
      else if (kind == Kind::INTERFACE_EPR && epr_type == InterfaceDielectric::DEFAULT)
      {
        plan.elem_a = FET.Elem1No;
        plan.elem_b = FET.Elem2No;
      }
      else if (kind == Kind::INTERFACE_EPR)
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
  if (need_field && fespace)
  {
    field_staging.SetSize(fespace->GetVSize());
    field_staging.UseDevice(true);
    field_staging = 0.0;
  }

  // Build the QFunction context for the integrand.
  std::vector<CeedIntScalar> ctx(2);
  ctx[0].second = 0.0;
  ctx[1].second = 0.0;
  if (kind == Kind::INTERFACE_EPR)
  {
    switch (epr_type)
    {
      case InterfaceDielectric::DEFAULT:
        ctx[0].second = 0.5 * epr_t * epr_epsilon;
        break;
      case InterfaceDielectric::MA:
        ctx[0].second = 0.5 * epr_t / epr_epsilon;
        break;
      case InterfaceDielectric::MS:
        {
          ctx[0].second = 0.5 * epr_t / epr_epsilon;
          MaterialPropertyCoefficient epsilon_func(mat_op->GetAttributeToMaterial(),
                                                   mat_op->GetPermittivityReal());
          auto mat_ctx = ceed::PopulateCoefficientContext(3, &epsilon_func);
          ctx.insert(ctx.end(), mat_ctx.begin(), mat_ctx.end());
        }
        break;
      case InterfaceDielectric::SA:
        ctx[0].second = 0.5 * epr_t * epr_epsilon;
        ctx[1].second = 0.5 * epr_t / epr_epsilon;
        break;
    }
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
    }

    // Volume geometry data evaluated at the mapped face quadrature points (for the
    // Piola transformations of the field inputs and, for INTERFACE_EPR MS, the material
    // property lookup with the local libCEED attribute).
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
    }
    if (has_b)
    {
      Vector elem_attr = GetCeedElemAttr(group.vol_indices_b);
      BuildGeometryData(ceed, mesh_fespace, group.vol_geom_b, group.vol_indices_b,
                        *group.mapped_ir_b, elem_attr, &vol_geom_data_b,
                        &vol_geom_data_b_restr);
    }

    // Field input restrictions (volume element dofs) and bases (volume element basis
    // tabulated at the mapped face quadrature points). The side b volume geometry data
    // is passed as an additional EvalMode::None field input ahead of the fields.
    std::vector<ceed::CeedFunctionalFieldInput> inputs;
    int num_fields = 0;
    CeedElemRestriction field_restr = nullptr, field_b_restr = nullptr;
    CeedBasis field_basis = nullptr, field_b_basis = nullptr;
    CeedVector field_vec = nullptr, field_b_vec = nullptr;
    if (has_b)
    {
      inputs.push_back({"vol_geom_data_b", vol_geom_data_b, vol_geom_data_b_restr, nullptr,
                        ceed::EvalMode::None});
    }
    auto AddFieldInput = [&](const std::string &name, const std::vector<int> &indices,
                             mfem::Geometry::Type geom, const mfem::IntegrationRule &ir,
                             CeedElemRestriction *restr, CeedBasis *basis, CeedVector *vec)
    {
      ceed::InitRestriction(fespace->Get(), indices, false, /*is_interp*/ true, false, ceed,
                            restr);
      const mfem::FiniteElement *fe = fespace->GetFEColl().FiniteElementForGeometry(geom);
      MFEM_VERIFY(fe, "Unable to get field finite element for surface functional!");
      ceed::InitBasisAtPoints(*fe, ir, fespace->GetVDim(), ceed, basis);
      ceed::InitCeedVector(field_staging, ceed, vec);
      inputs.push_back({name, *vec, *restr, *basis, ceed::EvalMode::Interp});
      num_fields++;
    };
    if (need_field)
    {
      AddFieldInput("u_1", group.vol_indices_a, group.vol_geom_a, *group.mapped_ir_a,
                    &field_restr, &field_basis, &field_vec);
      if (has_b)
      {
        AddFieldInput("u_2", group.vol_indices_b, group.vol_geom_b, *group.mapped_ir_b,
                      &field_b_restr, &field_b_basis, &field_b_vec);
      }
    }

    // Output restriction: one slot per boundary element in the local output vector.
    CeedElemRestriction out_restr;
    PalaceCeedCall(ceed, CeedElemRestrictionCreate(ceed, static_cast<CeedInt>(num_elem), 1,
                                                   1, num_marked, num_marked, CEED_MEM_HOST,
                                                   CEED_COPY_VALUES, group.out_slots.data(),
                                                   &out_restr));

    // Select the QFunction.
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
    }

    // Assemble the operator.
    CeedOperator op;
    ceed::AssembleCeedSurfaceFunctional(
        info, ctx.data(), ctx.size() * sizeof(CeedIntScalar), ceed, inputs, face_geom_data,
        face_geom_data_restr, vol_geom_data, vol_geom_data_restr, 1, out_restr, &op);
    groups.push_back({ceed, op, num_fields});

    // Cleanup (objects are now owned through the operator's references).
    PalaceCeedCall(ceed, CeedVectorDestroy(&face_geom_data));
    PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&face_geom_data_restr));
    if (vol_geom_data)
    {
      PalaceCeedCall(ceed, CeedVectorDestroy(&vol_geom_data));
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&vol_geom_data_restr));
    }
    if (vol_geom_data_b)
    {
      PalaceCeedCall(ceed, CeedVectorDestroy(&vol_geom_data_b));
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&vol_geom_data_b_restr));
    }
    if (field_restr)
    {
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&field_restr));
      PalaceCeedCall(ceed, CeedBasisDestroy(&field_basis));
      PalaceCeedCall(ceed, CeedVectorDestroy(&field_vec));
    }
    if (field_b_restr)
    {
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&field_b_restr));
      PalaceCeedCall(ceed, CeedBasisDestroy(&field_b_basis));
      PalaceCeedCall(ceed, CeedVectorDestroy(&field_b_vec));
    }
    PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&out_restr));
  }
}

void SurfaceFunctional::ApplyAdd(const Vector *u) const
{
  for (const auto &group : groups)
  {
    Ceed ceed = group.ceed;

    // Re-point the passive field inputs at the caller's data (the data or its location
    // may have changed since the last call).
    if (u)
    {
      static const char *field_names[] = {"u_1", "u_2"};
      for (int k = 0; k < group.num_fields; k++)
      {
        CeedOperatorField field;
        CeedVector field_vec;
        PalaceCeedCall(ceed, CeedOperatorGetFieldByName(group.op, field_names[k], &field));
        PalaceCeedCall(ceed, CeedOperatorFieldGetVector(field, &field_vec));
        ceed::InitCeedVector(*u, ceed, &field_vec, false);
      }
    }

    // Accumulate the per-element integrals into the local output vector.
    CeedVector out_vec;
    ceed::InitCeedVector(local_out, ceed, &out_vec);
    PalaceCeedCall(ceed, CeedOperatorApplyAdd(group.op, CEED_VECTOR_NONE, out_vec,
                                              CEED_REQUEST_IMMEDIATE));
    PalaceCeedCall(ceed, CeedVectorDestroy(&out_vec));
  }
}

double SurfaceFunctional::Eval(const Vector *u) const
{
  MFEM_VERIFY(kind == Kind::AREA || u,
              "SurfaceFunctional::Eval requires a field vector for functionals with "
              "field inputs!");
  MFEM_VERIFY(!u || !fespace || u->Size() == fespace->GetVSize(),
              "Invalid field vector size for SurfaceFunctional::Eval ("
                  << (u ? u->Size() : 0) << " vs. " << fespace->GetVSize() << ")!");

  double dot = 0.0;
  if (local_out.Size() > 0)
  {
    local_out = 0.0;
    ApplyAdd(u);
    dot = linalg::LocalSum(local_out);
  }
  Mpi::GlobalSum(1, &dot, comm);
  return dot;
}

double SurfaceFunctional::Eval(const GridFunction &u) const
{
  MFEM_VERIFY(kind != Kind::AREA, "SurfaceFunctional::Eval with a grid function is only "
                                  "valid for functionals with field inputs!");
  double dot = 0.0;
  if (local_out.Size() > 0)
  {
    local_out = 0.0;
    ApplyAdd(&u.Real());
    if (u.HasImag())
    {
      ApplyAdd(&u.Imag());
    }
    dot = linalg::LocalSum(local_out);
  }
  Mpi::GlobalSum(1, &dot, comm);
  return dot;
}

}  // namespace palace
