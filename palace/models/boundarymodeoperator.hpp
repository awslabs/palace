// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_BOUNDARY_MODE_OPERATOR_HPP
#define PALACE_MODELS_BOUNDARY_MODE_OPERATOR_HPP

#include <complex>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/vector.hpp"
#include "models/materialoperator.hpp"
#include "models/waveportoperator.hpp"
#include "utils/iodata.hpp"

namespace palace
{

class FarfieldBoundaryOperator;
class SurfaceConductivityOperator;
class SurfaceImpedanceOperator;

//
// Top-level operator for 2D boundary mode analysis, analogous to SpaceOperator for 3D
// driven/eigenmode problems. Owns the mesh, FE spaces, material operator, boundary
// operators, and the eigenvalue solver. Constructed from IoData and a mesh (2D directly
// or 3D with boundary attributes for submesh extraction).
//
class BoundaryModeOperator
{
public:
  // Result of an eigenvalue solve.
  struct SolveResult
  {
    int num_converged;
    double sigma;
  };

  BoundaryModeOperator(const IoData &iodata,
                       const std::vector<std::unique_ptr<Mesh>> &mesh);

  // Solve the eigenvalue problem at the given frequency.
  SolveResult Solve(double omega, double kn_target);

  // Access converged eigenvalues and eigenvectors.
  std::complex<double> GetEigenvalue(int i) const;
  void GetEigenvector(int i, ComplexVector &x) const;
  double GetError(int i, EigenvalueSolver::ErrorType type) const;

  // Access FE spaces.
  FiniteElementSpace &GetNDSpace() { return nd_fespaces.GetFinestFESpace(); }
  FiniteElementSpace &GetH1Space() { return h1_fespaces.GetFinestFESpace(); }
  FiniteElementSpace &GetCurlSpace() { return *l2_curl_fespace; }
  const FiniteElementSpace &GetNDSpace() const { return nd_fespaces.GetFinestFESpace(); }
  const FiniteElementSpace &GetH1Space() const { return h1_fespaces.GetFinestFESpace(); }
  const FiniteElementSpace &GetCurlSpace() const { return *l2_curl_fespace; }

  // Access space hierarchies (for error estimation).
  const FiniteElementSpaceHierarchy &GetNDSpaceHierarchy() const { return nd_fespaces; }
  const FiniteElementSpaceHierarchy &GetH1SpaceHierarchy() const { return h1_fespaces; }

  // Access FE collections.
  const mfem::FiniteElementCollection *GetNDFEColl() const { return nd_fecs.back().get(); }
  const mfem::FiniteElementCollection *GetH1FEColl() const { return h1_fecs.back().get(); }

  // Error estimation: add error indicator for a converged mode.
  void AddErrorIndicator(const ComplexVector &et, const ComplexVector &bz,
                         double total_domain_energy, ErrorIndicator &indicator);

  // Access material and mesh.
  const MaterialOperator &GetMaterialOp() const { return *mat_op; }
  Mesh &GetMesh() { return *solve_mesh; }
  const Mesh &GetMesh() const { return *solve_mesh; }
  MPI_Comm GetComm() const { return solve_mesh->GetComm(); }

  // Access the assembled Btt and Atn matrices (for power normalization).
  const mfem::HypreParMatrix *GetBtt() const;
  const mfem::HypreParMatrix *GetAtnr() const;
  const mfem::HypreParMatrix *GetAtni() const;

  // Access solver order.
  int GetSolverOrder() const { return solver_order; }

  // Get submesh projection data (for coordinate transforms).
  bool IsFromSubmesh() const { return use_submesh; }
  const mfem::Vector &GetSubmeshCentroid() const { return submesh_centroid; }
  const mfem::Vector &GetSubmeshE1() const { return submesh_e1; }
  const mfem::Vector &GetSubmeshE2() const { return submesh_e2; }

  // Access the linear solver (for metadata reporting).
  const ComplexKspSolver *GetLinearSolver() const;

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
  FiniteElementSpaceHierarchy nd_fespaces, h1_fespaces, h1_aux_fespaces;
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

  // Multigrid configuration for the eigenvalue solver (must outlive port_data).
  std::unique_ptr<ModeEigenSolverMultigridConfig> mg_config;

  // Wave port data (owns the eigenvalue solver internally).
  std::unique_ptr<WavePortData> port_data;

  // Error estimator (owns its FE spaces for flux recovery).
  std::unique_ptr<mfem::RT_FECollection> rt_fec_est;
  std::unique_ptr<FiniteElementSpaceHierarchy> nd_fespaces_est, rt_fespaces_est,
      h1_fespaces_est;
  std::unique_ptr<BoundaryModeFluxErrorEstimator<ComplexVector>> estimator;

  // Setup helpers.
  void SetUpMesh(const std::vector<std::unique_ptr<Mesh>> &mesh);
  void SetUpFESpaces(const std::vector<std::unique_ptr<Mesh>> &mesh);
  void SetUpEigenSolver();
};

}  // namespace palace

#endif  // PALACE_MODELS_BOUNDARY_MODE_OPERATOR_HPP
