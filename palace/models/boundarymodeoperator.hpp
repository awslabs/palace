// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_BOUNDARY_MODE_OPERATOR_HPP
#define PALACE_MODELS_BOUNDARY_MODE_OPERATOR_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
#include "models/materialoperator.hpp"
#include "utils/iodata.hpp"

namespace palace
{

class FarfieldBoundaryOperator;
class SurfaceConductivityOperator;
class SurfaceImpedanceOperator;

//
// Top-level operator for 2D boundary mode analysis, analogous to SpaceOperator for 3D
// driven/eigenmode problems. Owns the mesh, FE spaces, material operator, and boundary
// operators. Constructed from IoData and a mesh (2D directly or 3D with boundary
// attributes for submesh extraction). Does not own the eigenvalue solver or the error
// estimator — the driver constructs those on top of this operator's FE context.
//
class BoundaryModeOperator
{
public:
  BoundaryModeOperator(const IoData &iodata,
                       const std::vector<std::unique_ptr<Mesh>> &mesh);

  // Access FE spaces.
  FiniteElementSpace &GetNDSpace() { return nd_fespaces.GetFinestFESpace(); }
  FiniteElementSpace &GetH1Space() { return h1_fespaces.GetFinestFESpace(); }
  FiniteElementSpace &GetRTSpace() { return rt_fespaces.GetFinestFESpace(); }
  FiniteElementSpace &GetCurlSpace() { return *l2_curl_fespace; }
  const FiniteElementSpace &GetNDSpace() const { return nd_fespaces.GetFinestFESpace(); }
  const FiniteElementSpace &GetH1Space() const { return h1_fespaces.GetFinestFESpace(); }
  const FiniteElementSpace &GetRTSpace() const { return rt_fespaces.GetFinestFESpace(); }
  const FiniteElementSpace &GetCurlSpace() const { return *l2_curl_fespace; }

  // Access space hierarchies (for error estimation and multigrid).
  FiniteElementSpaceHierarchy &GetNDSpaceHierarchy() { return nd_fespaces; }
  FiniteElementSpaceHierarchy &GetH1SpaceHierarchy() { return h1_fespaces; }
  FiniteElementSpaceHierarchy &GetH1AuxSpaceHierarchy() { return h1_aux_fespaces; }
  FiniteElementSpaceHierarchy &GetRTSpaceHierarchy() { return rt_fespaces; }
  const FiniteElementSpaceHierarchy &GetNDSpaceHierarchy() const { return nd_fespaces; }
  const FiniteElementSpaceHierarchy &GetH1SpaceHierarchy() const { return h1_fespaces; }
  const FiniteElementSpaceHierarchy &GetH1AuxSpaceHierarchy() const
  {
    return h1_aux_fespaces;
  }
  const FiniteElementSpaceHierarchy &GetRTSpaceHierarchy() const { return rt_fespaces; }

  // Access per-level essential BC true DOF lists for the block system.
  std::vector<mfem::Array<int>> &GetNDDbcTDofLists() { return nd_dbc_tdof_lists; }
  std::vector<mfem::Array<int>> &GetH1DbcTDofLists() { return h1_dbc_tdof_lists; }
  std::vector<mfem::Array<int>> &GetH1AuxDbcTDofLists() { return h1_aux_dbc_tdof_lists; }

  // Access material and boundary operators.
  const MaterialOperator &GetMaterialOp() const { return *mat_op; }
  SurfaceImpedanceOperator *GetSurfZOp() { return surf_z_op.get(); }
  FarfieldBoundaryOperator *GetFarfieldOp() { return farfield_op.get(); }
  SurfaceConductivityOperator *GetSurfSigmaOp() { return surf_sigma_op.get(); }

  // Access mesh.
  Mesh &GetMesh() { return *solve_mesh; }
  const Mesh &GetMesh() const { return *solve_mesh; }
  MPI_Comm GetComm() const { return solve_mesh->GetComm(); }

  // Access solver order.
  int GetSolverOrder() const { return solver_order; }

  // Submesh projection data (for coordinate transforms and material tensor rotation).
  // For direct 2D, normal is nullptr.
  bool IsFromSubmesh() const { return use_submesh; }
  const mfem::Vector *GetSubmeshNormal() const
  {
    return use_submesh ? &submesh_normal : nullptr;
  }
  const mfem::Vector &GetSubmeshCentroid() const { return submesh_centroid; }
  const mfem::Vector &GetSubmeshE1() const { return submesh_e1; }
  const mfem::Vector &GetSubmeshE2() const { return submesh_e2; }

  // True vector sizes.
  int GetNDTrueVSize() const { return nd_fespaces.GetFinestFESpace().GetTrueVSize(); }
  int GetH1TrueVSize() const { return h1_fespaces.GetFinestFESpace().GetTrueVSize(); }

private:
  const IoData &iodata;
  int solver_order;
  bool use_submesh;

  // Mesh (owned for submesh case, non-owning pointer for direct 2D).
  std::unique_ptr<Mesh> owned_mesh;
  Mesh *solve_mesh;

  // Submesh projection geometry.
  mfem::Vector submesh_centroid, submesh_e1, submesh_e2, submesh_normal;

  // FE collections and space hierarchies.
  std::vector<std::unique_ptr<mfem::ND_FECollection>> nd_fecs;
  std::vector<std::unique_ptr<mfem::H1_FECollection>> h1_fecs;
  std::vector<std::unique_ptr<mfem::H1_FECollection>> h1_aux_fecs;
  std::vector<std::unique_ptr<mfem::RT_FECollection>> rt_fecs;
  FiniteElementSpaceHierarchy nd_fespaces, h1_fespaces, h1_aux_fespaces, rt_fespaces;
  std::vector<mfem::Array<int>> nd_dbc_tdof_lists, h1_dbc_tdof_lists, h1_aux_dbc_tdof_lists;

  // L2 curl space for B-field in 2D.
  std::unique_ptr<mfem::L2_FECollection> l2_curl_fec;
  std::unique_ptr<FiniteElementSpace> l2_curl_fespace;

  // Material and boundary operators.
  std::unique_ptr<MaterialOperator> mat_op;
  std::unique_ptr<SurfaceImpedanceOperator> surf_z_op;
  std::unique_ptr<FarfieldBoundaryOperator> farfield_op;
  std::unique_ptr<SurfaceConductivityOperator> surf_sigma_op;

  // DBC attributes.
  mfem::Array<int> dbc_bcs;

  // Setup helpers.
  void SetUpMesh(const std::vector<std::unique_ptr<Mesh>> &mesh);
  void SetUpFESpaces(const std::vector<std::unique_ptr<Mesh>> &mesh);
};

}  // namespace palace

#endif  // PALACE_MODELS_BOUNDARY_MODE_OPERATOR_HPP
