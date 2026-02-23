// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_QUAD_BOUNDARY_MODE_SOLVER_HPP
#define PALACE_MODELS_QUAD_BOUNDARY_MODE_SOLVER_HPP

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
// Quadratic eigenvalue problem (QEP) solver for 2D boundary mode computation with
// impedance boundary conditions. Derives the eigenvalue problem from the full 3D
// curl-curl weak form (not the divergence constraint), giving:
//
//   (M0 + lambda * M1 + lambda^2 * M2) [et; en] = 0,  lambda = ikn
//
// where impedance BCs enter M0 as boundary mass terms on BOTH et (ND tangential) and
// en (H1 normal), with no kn-dependence in any coefficient matrix. This avoids the
// limitation of the Vardapetyan-Demkowicz linearization (en_code = ikn * en_phys) which
// cannot capture the en component of the impedance BC.
//
// Uses SLEPc's PEP (Polynomial Eigenvalue Problem) solver.
//
class QuadBoundaryModeSolver
{
public:
  struct SolveResult
  {
    int num_converged;
    std::complex<double> lambda_target;
  };

  // Assemble frequency-dependent M0, frequency-independent M1 and M2, and configure the
  // PEP solver. Boundary conditions: PEC attributes get Dirichlet on both ND and H1.
  // Impedance attributes get Robin BCs (boundary mass in M0).
  QuadBoundaryModeSolver(const BoundaryModeSolverConfig &config,
                         const FiniteElementSpace &nd_fespace,
                         const FiniteElementSpace &h1_fespace,
                         const mfem::Array<int> &dbc_tdof_list, double omega,
                         MPI_Comm solver_comm = MPI_COMM_NULL);

  ~QuadBoundaryModeSolver();

  // Solve the QEP for eigenvalues lambda = ikn near the target kn_target.
  SolveResult Solve(double kn_target);

  // Get the i-th eigenvalue lambda = ikn.
  std::complex<double> GetEigenvalue(int i) const;

  // Get the i-th eigenvector [et; en] (size nd_size + h1_size).
  void GetEigenvector(int i, ComplexVector &x) const;

  // Convenience: get kn = -i * lambda from the i-th eigenvalue.
  std::complex<double> GetKn(int i) const;

  // Get the eigenpair error for the i-th converged eigenvalue.
  double GetError(int i, EigenvalueSolver::ErrorType type) const;

  // Get the true vector sizes for the ND and H1 FE spaces.
  int GetNDTrueVSize() const { return nd_size; }
  int GetH1TrueVSize() const { return h1_size; }

  // Access the Btt matrix (M2(1,1), negated) for impedance postprocessing.
  const mfem::HypreParMatrix *GetBtt() const { return Bttr.get(); }

  // Access the linear solver for metadata reporting.
  const ComplexKspSolver *GetLinearSolver() const { return ksp.get(); }

private:
  // Configuration.
  BoundaryModeSolverConfig config;

  // FE space references.
  const FiniteElementSpace &nd_fespace;
  const FiniteElementSpace &h1_fespace;

  // Dirichlet DOF list for the block system.
  mfem::Array<int> dbc_tdof_list;

  // FE space sizes.
  int nd_size, h1_size;

  // Assembled block operators: M0 + lambda * M1 + lambda^2 * M2 = 0.
  std::unique_ptr<ComplexOperator> opM0, opM1, opM2;

  // Btt (for postprocessing).
  std::unique_ptr<mfem::HypreParMatrix> Bttr;

  // Eigenvalue and linear solvers.
  std::unique_ptr<EigenvalueSolver> eigen;
  std::unique_ptr<ComplexKspSolver> ksp;

  // Shifted polynomial operator P(sigma) for the KSP (kept alive for the solve duration).
  std::unique_ptr<ComplexOperator> opP_shifted;

  // Assembly helpers.
  using HyprePtr = std::unique_ptr<mfem::HypreParMatrix>;
  using ComplexHypreParMatrix = std::tuple<HyprePtr, HyprePtr>;

  // Assemble M0: curl-curl + stiffness + mass + impedance boundary (frequency-dependent).
  ComplexHypreParMatrix AssembleM0(double omega) const;

  // Assemble M1: gradient coupling (frequency-independent).
  ComplexHypreParMatrix AssembleM1() const;

  // Assemble M2: ND mass (frequency-independent). Also saves Bttr for postprocessing.
  ComplexHypreParMatrix AssembleM2();

  // Build a 2x2 block HypreParMatrix from individual blocks and apply Dirichlet BCs.
  ComplexHypreParMatrix BuildBlockMatrix(
      const mfem::HypreParMatrix *block00r, const mfem::HypreParMatrix *block01r,
      const mfem::HypreParMatrix *block10r, const mfem::HypreParMatrix *block11r,
      const mfem::HypreParMatrix *block00i, const mfem::HypreParMatrix *block01i,
      const mfem::HypreParMatrix *block10i, const mfem::HypreParMatrix *block11i,
      Operator::DiagonalPolicy diag_policy) const;

  // Set up solvers.
  void SetUpLinearSolver(MPI_Comm comm);
  void SetUpEigenSolver(MPI_Comm comm);
};

}  // namespace palace

#endif  // PALACE_MODELS_QUAD_BOUNDARY_MODE_SOLVER_HPP
