// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_NEW_BOUNDARY_MODE_SOLVER_HPP
#define PALACE_MODELS_NEW_BOUNDARY_MODE_SOLVER_HPP

#include <complex>
#include <memory>
#include <mfem.hpp>
#include "linalg/eps.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/boundarymodesolver.hpp"

namespace palace
{

class FiniteElementSpace;

//
// Linear eigenvalue solver for 2D boundary mode computation using a direct linearization
// of the transverse curl-curl (Equation 1) and normal curl-curl (Equation 2) equations
// with the Vardapetyan-Demkowicz substitution e_n_tilde = i*kn*E_n, e_t = E_t.
//
// Assembles and solves the generalized eigenvalue problem:
//
//   [Att  Atn] [et]          [Btt  Btn] [et]
//   [ 0   Ann] [en] = lambda  [ 0    0 ] [en]
//
// where:
//   Att = mu_cc^{-1}(curl_t u, curl_t v) - omega^2(eps u, v) + BC-t  [same as BMS]
//   Atn = -(mu^{-1} grad_t u, v)                                      [same as BMS]
//   Ann = -(mu^{-1} grad u, grad v) + omega^2(eps u, v) + BC-n        [NEW: H1 stiffness
//          + mass + impedance boundary]
//   Btt = -(mu^{-1} u, v)                                              [same as BMS]
//   Btn = -(mu^{-1} u, grad v) = Atn^T                                 [NEW: transpose]
//
// The A matrix is upper block-triangular (Ant = 0). The B matrix has Btn coupling from
// Equation 2 and Bnn = 0. This is a standard GEP with shift-and-invert.
//
// The eigenvector contains [e_t; e_n_tilde] where e_n_tilde = i*kn*E_n. Recovery of
// E_n = e_n_tilde / (i*kn) is performed in the driver, not here.
//
class NewBoundaryModeSolver
{
public:
  struct SolveResult
  {
    int num_converged;
    double sigma;
  };

  // The constructor assembles frequency-independent matrices (Atn, Btn = Atn^T, Btt)
  // and configures linear and eigenvalue solvers. Frequency-dependent matrices (Att, Ann)
  // are assembled at solve time via AssembleFrequencyDependent(). Matrix assembly uses
  // the FE space communicator (all processes). If solver_comm != MPI_COMM_NULL, the
  // linear and eigenvalue solvers are configured on that communicator.
  NewBoundaryModeSolver(const BoundaryModeSolverConfig &config,
                        const FiniteElementSpace &nd_fespace,
                        const FiniteElementSpace &h1_fespace,
                        const mfem::Array<int> &dbc_tdof_list,
                        MPI_Comm solver_comm = MPI_COMM_NULL);

  ~NewBoundaryModeSolver();

  // Assemble frequency-dependent matrices and solve the eigenvalue problem. The shift
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
  // Atn: gradient coupling -(mu^{-1} grad_t u, v).
  // Btn: transpose of Atn, coupling from Equation 2: -(mu^{-1} u, grad_t v).
  // Bttr: positive copy of Btt for external use (postprocessing).
  std::unique_ptr<mfem::HypreParMatrix> Atnr, Atni;
  std::unique_ptr<mfem::HypreParMatrix> Btnr, Btni;
  std::unique_ptr<mfem::HypreParMatrix> Bttr;
  std::unique_ptr<ComplexOperator> opB;

  // Frequency-dependent block A operator (rebuilt each Solve).
  std::unique_ptr<ComplexOperator> opA;

  // Eigenvalue and linear solvers (null on processes without solver_comm).
  std::unique_ptr<EigenvalueSolver> eigen;
  std::unique_ptr<ComplexKspSolver> ksp;

  // Assemble frequency-dependent Att and Ann, then build block A (MPI collective on FE
  // space comm).
  void AssembleFrequencyDependent(double omega, double sigma);

  // Private helper methods for bilinear form assembly.
  using ComplexHypreParMatrix = std::tuple<std::unique_ptr<mfem::HypreParMatrix>,
                                           std::unique_ptr<mfem::HypreParMatrix>>;

  // Att: same as BoundaryModeSolver (curl-curl + mass + BC-t).
  ComplexHypreParMatrix AssembleAtt(double omega, double sigma) const;

  // Atn: gradient coupling -(mu^{-1} grad_t u, v), same as BoundaryModeSolver.
  ComplexHypreParMatrix AssembleAtn() const;

  // Ann: H1 stiffness + mass + BC-n (NEW formulation).
  // Ann = -(mu^{-1} grad u, grad v) + omega^2(eps u, v) + BC-n impedance.
  ComplexHypreParMatrix AssembleAnn(double omega) const;

  // Btt: ND mass -(mu^{-1} u, v), same as BoundaryModeSolver.
  ComplexHypreParMatrix AssembleBtt() const;

  // Build the 2x2 block A matrix. The (1,0) block is -sigma * Btn from the
  // shift-and-invert transformation (nullptr when sigma = 0).
  ComplexHypreParMatrix
  BuildSystemMatrixA(const mfem::HypreParMatrix *Attr, const mfem::HypreParMatrix *Atti,
                     const mfem::HypreParMatrix *Atnr, const mfem::HypreParMatrix *Atni,
                     const mfem::HypreParMatrix *Annr, const mfem::HypreParMatrix *Anni,
                     const mfem::HypreParMatrix *shifted_Btnr = nullptr) const;

  // Build the 2x2 block B matrix: [Btt, 0; Btn, 0].
  ComplexHypreParMatrix BuildSystemMatrixB(const mfem::HypreParMatrix *Bttr,
                                           const mfem::HypreParMatrix *Btti,
                                           const mfem::HypreParMatrix *Btnr,
                                           const mfem::HypreParMatrix *Btni,
                                           const mfem::HypreParMatrix *Dnn) const;

  // Set up the linear solver (GMRES + sparse direct preconditioner).
  void SetUpLinearSolver(MPI_Comm comm);

  // Set up the eigenvalue solver (SLEPc or ARPACK).
  void SetUpEigenSolver(MPI_Comm comm);
};

}  // namespace palace

#endif  // PALACE_MODELS_NEW_BOUNDARY_MODE_SOLVER_HPP
