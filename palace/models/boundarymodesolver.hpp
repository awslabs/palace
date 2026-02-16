// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_BOUNDARY_MODE_SOLVER_HPP
#define PALACE_MODELS_BOUNDARY_MODE_SOLVER_HPP

#include <complex>
#include <memory>
#include <mfem.hpp>
#include "linalg/eps.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "utils/configfile.hpp"
#include "utils/labels.hpp"

namespace palace
{

class FiniteElementSpace;

//
// Configuration for the BoundaryModeSolver, parameterizing the differences between the
// 2D mode analysis (ModeAnalysisSolver) and 3D wave port (WavePortOperator) eigenvalue
// problems.
//
struct BoundaryModeSolverConfig
{
  // Material property mappings. For mode analysis, these come from
  // GetAttributeToMaterial() (volume attrs). For wave ports, from
  // GetBdrAttributeToMaterial() (boundary attrs).
  const mfem::Array<int> *attr_to_material;

  // Inverse permeability tensor. Used for Atn, Ant, and Btt bilinear forms.
  const mfem::DenseTensor *inv_permeability;

  // Curl-curl inverse permeability. In 2D mode analysis this is the scalar mu_inv
  // (GetCurlCurlInvPermeability); for 3D wave ports it is the full tensor
  // (GetInvPermeability), with subsequent normal projection.
  const mfem::DenseTensor *curlcurl_inv_permeability;

  // Real part of permittivity tensor. Used for Att, Ant bilinear forms.
  const mfem::DenseTensor *permittivity_real;

  // Scalar permittivity for the Ann mass matrix. In 2D mode analysis this is
  // GetPermittivityScalar(); for 3D wave ports it is GetPermittivityReal() with
  // subsequent normal projection.
  const mfem::DenseTensor *permittivity_scalar;

  // Imaginary part of permittivity (loss tangent contribution). Set to nullptr when
  // there is no loss tangent.
  const mfem::DenseTensor *permittivity_imag;

  // Surface normal vector for 3D wave port boundaries. Set to nullptr for 2D domain
  // meshes (mode analysis), which do not need normal projection.
  const mfem::Vector *normal;

  // Whether the material operator reports a nonzero loss tangent.
  bool has_loss_tangent;

  // Eigenvalue solver configuration.
  int num_modes;
  int num_vec;
  double eig_tol;
  EigenvalueSolver::WhichType which_eig;

  // Linear solver configuration.
  const config::LinearSolverData *linear;
  EigenSolverBackend eigen_backend;

  // Verbosity level for solvers.
  int verbose;
};

//
// Shared eigenvalue problem solver for 2D boundary mode computation. Assembles and solves
// the block eigenvalue problem:
//
//   [Att  Atn] [et]          [Btt  0] [et]
//   [Ant  Ann] [en] = lambda  [0    0] [en]
//
// with shift-and-invert spectral transformation. Used by both ModeAnalysisSolver (2D
// waveguide cross-section) and WavePortOperator (3D boundary submesh).
//
class BoundaryModeSolver
{
public:
  struct SolveResult
  {
    int num_converged;
    double sigma;
  };

  // The constructor assembles frequency-independent matrices (Atn, Ant, Ann, Btt, block B)
  // and configures linear and eigenvalue solvers. Matrix assembly uses the FE space
  // communicator (all processes). If solver_comm != MPI_COMM_NULL, the linear and
  // eigenvalue solvers are configured on that communicator (subset of processes, e.g. for
  // wave port boundaries). If solver_comm == MPI_COMM_NULL, the FE space communicator is
  // used for solvers as well.
  BoundaryModeSolver(const BoundaryModeSolverConfig &config,
                     const FiniteElementSpace &nd_fespace,
                     const FiniteElementSpace &h1_fespace,
                     const mfem::Array<int> &dbc_tdof_list,
                     MPI_Comm solver_comm = MPI_COMM_NULL);

  ~BoundaryModeSolver();

  // Assemble frequency-dependent Att matrix and solve the eigenvalue problem. The shift
  // sigma = -kn_target^2 is applied. An optional initial space vector can be provided
  // for eigenvalue solver warm-starting. Use Solve() when all processes have solvers
  // (ModeAnalysis). Use SolveSplit() when only a subset of processes have solvers
  // (WavePort): assembly runs on all processes, solve only on those with has_solver=true.
  SolveResult Solve(double omega, double sigma,
                    const ComplexVector *initial_space = nullptr);
  SolveResult SolveSplit(double omega, double sigma, bool has_solver,
                         const ComplexVector *initial_space = nullptr);

  // Access converged eigenvalues and eigenvectors.
  std::complex<double> GetEigenvalue(int i) const;
  void GetEigenvector(int i, ComplexVector &x) const;

  // Access the assembled Btt matrix (needed for impedance postprocessing).
  const mfem::HypreParMatrix *GetBtt() const { return Bttr.get(); }

  // Get the true vector sizes for the ND and H1 FE spaces.
  int GetNDTrueVSize() const { return nd_size; }
  int GetH1TrueVSize() const { return h1_size; }

  // Access the linear solver (for metadata reporting). Returns nullptr if this process
  // does not have a solver configured (non-port process in wave port mode).
  const ComplexKspSolver *GetLinearSolver() const { return ksp.get(); }

private:
  // Configuration (stored for Solve-time assembly).
  BoundaryModeSolverConfig config;

  // References to FE spaces (not owned).
  const FiniteElementSpace &nd_fespace;
  const FiniteElementSpace &h1_fespace;

  // Essential boundary condition true DOF list for the combined block system.
  mfem::Array<int> dbc_tdof_list;

  // Cached FE space sizes.
  int nd_size, h1_size;

  // Frequency-independent assembled matrices.
  std::unique_ptr<mfem::HypreParMatrix> Atnr, Atni, Antr, Anti, Annr, Anni;
  std::unique_ptr<mfem::HypreParMatrix> Bttr;
  std::unique_ptr<ComplexOperator> opB;

  // Frequency-dependent block A operator (rebuilt each Solve).
  std::unique_ptr<ComplexOperator> opA;

  // Eigenvalue and linear solvers (null on processes without solver_comm).
  std::unique_ptr<EigenvalueSolver> eigen;
  std::unique_ptr<ComplexKspSolver> ksp;

  // Assemble frequency-dependent Att and build block A (MPI collective on FE space comm).
  void AssembleFrequencyDependent(double omega, double sigma);

  // Private helper methods for bilinear form assembly.
  using ComplexHypreParMatrix = std::tuple<std::unique_ptr<mfem::HypreParMatrix>,
                                           std::unique_ptr<mfem::HypreParMatrix>>;

  ComplexHypreParMatrix AssembleAtt(double omega, double sigma) const;
  ComplexHypreParMatrix AssembleAtn() const;
  ComplexHypreParMatrix AssembleAnt() const;
  ComplexHypreParMatrix AssembleAnn() const;
  ComplexHypreParMatrix AssembleBtt() const;

  ComplexHypreParMatrix
  BuildSystemMatrixA(const mfem::HypreParMatrix *Attr, const mfem::HypreParMatrix *Atti,
                     const mfem::HypreParMatrix *Atnr, const mfem::HypreParMatrix *Atni,
                     const mfem::HypreParMatrix *Antr, const mfem::HypreParMatrix *Anti,
                     const mfem::HypreParMatrix *Annr,
                     const mfem::HypreParMatrix *Anni) const;

  ComplexHypreParMatrix BuildSystemMatrixB(const mfem::HypreParMatrix *Bttr,
                                           const mfem::HypreParMatrix *Btti,
                                           const mfem::HypreParMatrix *Dnn) const;

  // Set up the linear solver (GMRES + sparse direct preconditioner).
  void SetUpLinearSolver(MPI_Comm comm);

  // Set up the eigenvalue solver (SLEPc or ARPACK).
  void SetUpEigenSolver(MPI_Comm comm);
};

}  // namespace palace

#endif  // PALACE_MODELS_BOUNDARY_MODE_SOLVER_HPP
