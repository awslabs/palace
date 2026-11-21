// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_KSP_HPP
#define PALACE_LINALG_KSP_HPP

#include <memory>
#include <type_traits>
#include <vector>
#include "linalg/iterative.hpp"
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "utils/labels.hpp"

namespace palace
{

class FiniteElementSpaceHierarchy;
class IoData;

namespace config
{

struct LinearSolverData;

}  // namespace config

// Derive the preconditioner matrix symmetry (used by sparse direct solvers) from the
// problem type, shifted-preconditioner flag, and periodicity. This is a derived runtime
// property, not a configuration field, so it is computed here rather than stored on
// LinearSolverData.
MatrixSymmetry GetPreconditionerMatrixSymmetry(const IoData &iodata);

// Whether the coarse-level solver of a complex-valued multigrid preconditioner is applied
// to the real-valued Br + Bi of its complex-valued operator (see MfemWrapperSolver), rather
// than to the complex-valued operator: AMS, BoomerAMG, and the sparse direct solvers
// without ComplexCoarseSolve.
bool UseRealCoarseApproximation(const config::LinearSolverData &linear);

//
// Linear solver class composing an iterative solver and preconditioner object.
//
template <OperatorType OperType>
class BaseKspSolver
{
  using VecType =
      std::conditional_t<std::is_same_v<OperType, ComplexOperator>, ComplexVector, Vector>;

protected:
  // The actual solver and preconditioner objects.
  std::unique_ptr<IterativeSolver<OperType>> ksp;
  std::unique_ptr<Solver<OperType>> pc;

  // Counters for number of calls to Mult method for linear solves, and cumulative number
  // of iterations.
  mutable int ksp_mult, ksp_mult_it;

  // Enable timer contribution for Timer::KSP_PRECONDITIONER.
  bool use_timer;

public:
  BaseKspSolver(const config::LinearSolverData &linear, MatrixSymmetry pc_mat_sym,
                int verbose, FiniteElementSpaceHierarchy &fespaces,
                FiniteElementSpaceHierarchy *aux_fespaces = nullptr);
  BaseKspSolver(const IoData &iodata, FiniteElementSpaceHierarchy &fespaces,
                FiniteElementSpaceHierarchy *aux_fespaces = nullptr);
  BaseKspSolver(std::unique_ptr<IterativeSolver<OperType>> &&ksp,
                std::unique_ptr<Solver<OperType>> &&pc);

  int NumTotalMult() const { return ksp_mult; }
  int NumTotalMultIterations() const { return ksp_mult_it; }

  // Whether the most recent Mult converged to the requested tolerance.
  bool GetConverged() const { return ksp->GetConverged(); }

  // Forward tolerance access to the underlying iterative solver.
  double GetRelTol() const { return ksp->GetRelTol(); }
  double GetAbsTol() const { return ksp->GetAbsTol(); }
  void SetRelTol(double tol) { ksp->SetRelTol(tol); }
  void SetAbsTol(double tol) { ksp->SetAbsTol(tol); }

  // The least-squares convergence test of the CG solver: its tolerance and ‖d‖² (see
  // CgSolver::SetLeastSquaresTol), and the decrease of the objective over the iterations of
  // the most recent Mult (see CgSolver::GetLeastSquaresDecrease). The iterative solver must
  // be a CgSolver.
  void SetLeastSquaresTol(double tol);
  void SetLeastSquaresNorm(double norm2);
  std::vector<double> GetLeastSquaresDecrease() const;

  // Number of iterations of the most recent Mult.
  int GetNumIterations() const { return ksp->GetNumIterations(); }

  void SetOperators(const OperType &op, const OperType &pc_op);

  // Whether the iterative solver uses the given solution vector as initial guess.
  void SetInitialGuess(bool guess) { ksp->SetInitialGuess(guess); }

  // Stop the solve as soon as the correction to a nonzero initial guess exceeds tol times
  // the solution norm, when the Krylov solver can monitor it (GMRES with left or no
  // preconditioning; see GmresSolver::SetMaxRelativeCorrection). No limit if tol = 0.
  void SetMaxRelativeCorrection(double tol);

  // Whether the most recent Mult was stopped by the correction limit.
  bool GetCorrectionExceeded() const;

  void Mult(const VecType &x, VecType &y) const;
};

using KspSolver = BaseKspSolver<Operator>;
using ComplexKspSolver = BaseKspSolver<ComplexOperator>;

}  // namespace palace

#endif  // PALACE_LINALG_KSP_HPP
