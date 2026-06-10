// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfacefunctional.hpp"

#include <map>
#include <memory>
#include <mutex>
#include "fem/coefficient.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "fem/libceed/basis.hpp"
#include "fem/libceed/functional.hpp"
#include "fem/libceed/integrator.hpp"
#include "fem/libceed/restriction.hpp"
#include "fem/mesh.hpp"
#include "linalg/vector.hpp"
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

// Key uniquely identifying the map from boundary element reference coordinates to
// attached volume element reference coordinates (and thus the tabulated basis): the
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
std::map<FaceConfigKey, std::unique_ptr<mfem::IntegrationRule>> &MappedIrRegistry()
{
  static std::map<FaceConfigKey, std::unique_ptr<mfem::IntegrationRule>> registry;
  return registry;
}

std::mutex &MappedIrRegistryMutex()
{
  static std::mutex m;
  return m;
}

// Data for one group of boundary elements sharing the same face configuration (mapped
// quadrature point positions in the volume element reference space).
struct FaceGroup
{
  mfem::Geometry::Type bdr_geom, vol_geom;
  bool flip_normal;  // Boundary element and face have the same orientation (o % 2 == 0),
                     // legacy coefficients invert the boundary element normal
  const mfem::IntegrationRule *mapped_ir;  // From the registry, application lifetime
  std::vector<int> bdr_indices, vol_indices;
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

}  // namespace

SurfaceFunctional::SurfaceFunctional(Kind kind, const Mesh &mesh,
                                     const mfem::Array<int> &bdr_attr_marker,
                                     const FiniteElementSpace *fespace)
  : kind(kind), fespace(fespace), comm(mesh.GetComm())
{
  MFEM_VERIFY(kind == Kind::AREA || fespace,
              "SurfaceFunctional requires a field finite element space for functionals "
              "with field inputs!");
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

  // Group the marked boundary elements by their face configuration: the mapped
  // positions of the face quadrature points in the attached volume element's reference
  // space (and the geometries and normal orientation flag). All elements in a group
  // share the tabulated basis and are processed by a single libCEED operator.
  std::map<FaceConfigKey, FaceGroup> face_groups;
  int num_marked = 0;
  {
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
      // (FET.Elem1 always exists and is local).
      const bool flip = BdrGridFunctionCoefficient::GetBdrElementNeighborTransformations(
          i, pmesh, FET, T1, T2);
      const auto bdr_geom = pmesh.GetBdrElementGeometry(i);
      const auto vol_geom = pmesh.GetElementGeometry(FET.Elem1No);
      int f, o;
      pmesh.GetBdrElementFace(i, &f, &o);

      // Map the face quadrature points to the volume element reference space:
      // boundary element reference coordinates -> face reference coordinates
      // (TransformBdrElementToFace with the boundary element to face orientation) ->
      // element 1 reference coordinates (FET.Loc1).
      const int q_order = fem::DefaultIntegrationOrder::Get(pmesh, bdr_geom);
      const mfem::IntegrationRule &face_ir = mfem::IntRules.Get(bdr_geom, q_order);
      const int nq = face_ir.GetNPoints();
      std::vector<mfem::IntegrationPoint> mapped_ips(nq);
      for (int q = 0; q < nq; q++)
      {
        mfem::IntegrationPoint fip = mfem::Mesh::TransformBdrElementToFace(
            FET.GetGeometryType(), o, face_ir.IntPoint(q));
        FET.Loc1.Transform(fip, mapped_ips[q]);
        mapped_ips[q].weight = face_ir.IntPoint(q).weight;
      }

      // Build the group key from the geometries, orientation flag, and quantized
      // mapped quadrature point coordinates.
      FaceConfigKey key;
      key.reserve(4 + 3 * nq);
      key.push_back(static_cast<long long>(vol_geom));
      key.push_back(static_cast<long long>(bdr_geom));
      key.push_back(static_cast<long long>(flip));
      key.push_back(static_cast<long long>(nq));
      for (int q = 0; q < nq; q++)
      {
        const auto &ip = mapped_ips[q];
        key.push_back(std::llround(ip.x * QUANTIZE_SCALE));
        key.push_back(std::llround(ip.y * QUANTIZE_SCALE));
        key.push_back(std::llround(ip.z * QUANTIZE_SCALE));
      }

      auto it = face_groups.find(key);
      if (it == face_groups.end())
      {
        FaceGroup group;
        group.bdr_geom = bdr_geom;
        group.vol_geom = vol_geom;
        group.flip_normal = flip;
        {
          // Insert the mapped integration rule into the (application lifetime)
          // registry, or reuse an existing entry.
          std::lock_guard<std::mutex> lock(MappedIrRegistryMutex());
          auto &registry = MappedIrRegistry();
          auto ir_it = registry.find(key);
          if (ir_it == registry.end())
          {
            auto ir = std::make_unique<mfem::IntegrationRule>(nq);
            for (int q = 0; q < nq; q++)
            {
              ir->IntPoint(q) = mapped_ips[q];
            }
            ir_it = registry.emplace(key, std::move(ir)).first;
          }
          group.mapped_ir = ir_it->second.get();
        }
        it = face_groups.emplace(key, std::move(group)).first;
      }
      it->second.bdr_indices.push_back(i);
      it->second.vol_indices.push_back(FET.Elem1No);
      it->second.out_slots.push_back(num_marked++);
    }
  }

  // Initialize the local output vector and field staging vector.
  local_out.SetSize(num_marked);
  local_out.UseDevice(true);
  if (fespace)
  {
    field_staging.SetSize(fespace->GetVSize());
    field_staging.UseDevice(true);
    field_staging = 0.0;
  }

  // Assemble a libCEED operator for each group. For now, all operators are constructed
  // on a single Ceed context (no OpenMP parallel assembly or application; correctness
  // first, this can be extended with the thread partitioning of fem/mesh.cpp later).
  Ceed ceed = ceed::internal::GetCeedObjects()[0];
  for (auto &[key, group] : face_groups)
  {
    const std::size_t num_elem = group.bdr_indices.size();
    const mfem::IntegrationRule &face_ir = mfem::IntRules.Get(
        group.bdr_geom, fem::DefaultIntegrationOrder::Get(pmesh, group.bdr_geom));

    // Face geometry data (boundary element Jacobians at the face quadrature points).
    // The attribute channel currently stores the raw mfem boundary attribute (material
    // lookups via libCEED attribute indirection are wired up with the specific
    // integrand kinds that require them).
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
    // Piola transformations of the field inputs). Not needed for field-less
    // functionals.
    CeedVector vol_geom_data = nullptr;
    CeedElemRestriction vol_geom_data_restr = nullptr;
    const bool need_vol_geom = (kind != Kind::AREA);
    if (need_vol_geom)
    {
      Vector elem_attr(num_elem);
      for (std::size_t k = 0; k < num_elem; k++)
      {
        elem_attr[k] = pmesh.GetAttribute(group.vol_indices[k]);
      }
      BuildGeometryData(ceed, mesh_fespace, group.vol_geom, group.vol_indices,
                        *group.mapped_ir, elem_attr, &vol_geom_data, &vol_geom_data_restr);
    }

    // Field input restriction (volume element dofs) and basis (volume element basis
    // tabulated at the mapped face quadrature points).
    std::vector<ceed::CeedFunctionalFieldInput> inputs;
    CeedElemRestriction field_restr = nullptr;
    CeedBasis field_basis = nullptr;
    CeedVector field_vec = nullptr;
    if (kind == Kind::HCURL_NORM2)
    {
      ceed::InitRestriction(fespace->Get(), group.vol_indices, false, /*is_interp*/ true,
                            false, ceed, &field_restr);
      const mfem::FiniteElement *fe =
          fespace->GetFEColl().FiniteElementForGeometry(group.vol_geom);
      MFEM_VERIFY(fe, "Unable to get field finite element for surface functional!");
      ceed::InitBasisAtPoints(*fe, *group.mapped_ir, fespace->GetVDim(), ceed,
                              &field_basis);
      ceed::InitCeedVector(field_staging, ceed, &field_vec);
      inputs.push_back(
          {"u_1", field_vec, field_restr, field_basis, ceed::EvalMode::Interp});
    }

    // Output restriction: one slot per boundary element in the local output vector.
    CeedElemRestriction out_restr;
    PalaceCeedCall(ceed, CeedElemRestrictionCreate(ceed, static_cast<CeedInt>(num_elem), 1,
                                                   1, num_marked, num_marked, CEED_MEM_HOST,
                                                   CEED_COPY_VALUES, group.out_slots.data(),
                                                   &out_restr));

    // Assemble the operator.
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
    }

    CeedOperator op;
    ceed::AssembleCeedSurfaceFunctional(info, nullptr, 0, ceed, inputs, face_geom_data,
                                        face_geom_data_restr, vol_geom_data,
                                        vol_geom_data_restr, 1, out_restr, &op);
    groups.push_back({ceed, op});

    // Cleanup (objects are now owned through the operator's references).
    PalaceCeedCall(ceed, CeedVectorDestroy(&face_geom_data));
    PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&face_geom_data_restr));
    if (vol_geom_data)
    {
      PalaceCeedCall(ceed, CeedVectorDestroy(&vol_geom_data));
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&vol_geom_data_restr));
    }
    if (field_restr)
    {
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&field_restr));
      PalaceCeedCall(ceed, CeedBasisDestroy(&field_basis));
      PalaceCeedCall(ceed, CeedVectorDestroy(&field_vec));
    }
    PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&out_restr));
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
    for (const auto &group : groups)
    {
      Ceed ceed = group.ceed;

      // Re-point the passive field input at the caller's data (the data or its
      // location may have changed since the last call).
      if (u)
      {
        CeedOperatorField field;
        CeedVector field_vec;
        PalaceCeedCall(ceed, CeedOperatorGetFieldByName(group.op, "u_1", &field));
        PalaceCeedCall(ceed, CeedOperatorFieldGetVector(field, &field_vec));
        ceed::InitCeedVector(*u, ceed, &field_vec, false);
      }

      // Accumulate the per-element integrals into the local output vector.
      CeedVector out_vec;
      ceed::InitCeedVector(local_out, ceed, &out_vec);
      PalaceCeedCall(ceed, CeedOperatorApplyAdd(group.op, CEED_VECTOR_NONE, out_vec,
                                                CEED_REQUEST_IMMEDIATE));
      PalaceCeedCall(ceed, CeedVectorDestroy(&out_vec));
    }
    dot = linalg::LocalSum(local_out);
  }
  Mpi::GlobalSum(1, &dot, comm);
  return dot;
}

}  // namespace palace
