// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_ITERATIVE_HPP
#define PALACE_LINALG_ITERATIVE_HPP

#include <type_traits>
#include <vector>
#include <mfem.hpp>
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "linalg/vector.hpp"
#include "utils/labels.hpp"

namespace palace
{

//
// Iterative solvers based on Krylov subspace methods with optional preconditioning, for
// real- or complex-valued systems.
//

// Base class for iterative solvers based on Krylov subspace methods with optional
// preconditioning.
template <OperatorType OperType>
class IterativeSolver : public Solver<OperType>
{
protected:
  using RealType = double;
  using ScalarType = std::conditional_t<std::is_same_v<OperType, ComplexOperator>,
                                        std::complex<RealType>, RealType>;

  // MPI communicator associated with the solver.
  MPI_Comm comm;

  // Control level of printing during solves.
  mfem::IterativeSolver::PrintLevel print_opts;
  int int_width, tab_width;

  // Relative and absolute tolerances.
  double rel_tol, abs_tol;

  // Limit for the number of solver iterations.
  int max_it;

  // Operator and (optional) preconditioner associated with the iterative solver (not
  // owned).
  const OperType *A;
  const Solver<OperType> *B;

  // Variables set during solve to capture solve statistics.
  mutable bool converged;
  mutable double initial_res, final_res;
  mutable int final_it;

  // Enable timer contribution for Timer::PRECONDITIONER.
  bool use_timer;

public:
  IterativeSolver(MPI_Comm comm, int print);

  // Set an indentation for all log printing.
  void SetTabWidth(int width) { tab_width = width; }

  // Set the relative convergence tolerance.
  void SetTol(double tol) { SetRelTol(tol); }
  void SetRelTol(double tol) { rel_tol = tol; }

  // Get the relative convergence tolerance.
  double GetRelTol() const { return rel_tol; }

  // Get the absolute convergence tolerance.
  double GetAbsTol() const { return abs_tol; }

  // Set the absolute convergence tolerance.
  void SetAbsTol(double tol) { abs_tol = tol; }

  // Set the maximum number of iterations.
  void SetMaxIter(int its)
  {
    max_it = its;
    int_width = 1 + static_cast<int>(std::log10(its));
  }

  // Set the operator for the solver.
  void SetOperator(const OperType &op) override
  {
    A = &op;
    this->height = op.Height();
    this->width = op.Width();
  }

  // Set the preconditioner for the solver.
  void SetPreconditioner(const Solver<OperType> &pc) { B = &pc; }

  // Returns if the previous solve converged or not.
  bool GetConverged() const { return converged && (rel_tol > 0.0 || abs_tol > 0.0); }

  // Returns the initial (absolute) residual for the previous solve.
  double GetInitialRes() const { return initial_res; }

  // Returns the final (absolute) residual for the previous solve, which may be an estimate
  // to the true residual.
  double GetFinalRes() const { return final_res; }

  // Returns the number of iterations for the previous solve.
  int GetNumIterations() const { return final_it; }

  // Get the associated MPI communicator.
  MPI_Comm GetComm() const { return comm; }

  // Activate preconditioner timing during solves.
  void EnableTimer() { use_timer = true; }
};

// Preconditioned Conjugate Gradient (CG) method for SPD linear systems.
template <OperatorType OperType>
class CgSolver : public IterativeSolver<OperType>
{
protected:
  using VecType = typename Solver<OperType>::VecType;
  using RealType = typename IterativeSolver<OperType>::RealType;
  using ScalarType = typename IterativeSolver<OperType>::ScalarType;

  using IterativeSolver<OperType>::comm;
  using IterativeSolver<OperType>::print_opts;
  using IterativeSolver<OperType>::int_width;
  using IterativeSolver<OperType>::tab_width;

  using IterativeSolver<OperType>::rel_tol;
  using IterativeSolver<OperType>::abs_tol;
  using IterativeSolver<OperType>::max_it;

  using IterativeSolver<OperType>::A;
  using IterativeSolver<OperType>::B;

  using IterativeSolver<OperType>::converged;
  using IterativeSolver<OperType>::initial_res;
  using IterativeSolver<OperType>::final_res;
  using IterativeSolver<OperType>::final_it;

  // Temporary workspace for solve.
  mutable VecType r, z, p;

public:
  CgSolver(MPI_Comm comm, int print) : IterativeSolver<OperType>(comm, print) {}

  void Mult(const VecType &b, VecType &x) const override;
};

// Preconditioned Generalized Minimum Residual Method (GMRES) for general nonsymmetric
// linear systems.
template <OperatorType OperType>
class GmresSolver : public IterativeSolver<OperType>
{
protected:
  using VecType = typename Solver<OperType>::VecType;
  using RealType = typename IterativeSolver<OperType>::RealType;
  using ScalarType = typename IterativeSolver<OperType>::ScalarType;

  using IterativeSolver<OperType>::comm;
  using IterativeSolver<OperType>::print_opts;
  using IterativeSolver<OperType>::int_width;
  using IterativeSolver<OperType>::tab_width;

  using IterativeSolver<OperType>::rel_tol;
  using IterativeSolver<OperType>::abs_tol;
  using IterativeSolver<OperType>::max_it;

  using IterativeSolver<OperType>::A;
  using IterativeSolver<OperType>::B;

  using IterativeSolver<OperType>::converged;
  using IterativeSolver<OperType>::initial_res;
  using IterativeSolver<OperType>::final_res;
  using IterativeSolver<OperType>::final_it;

  // Maximum subspace dimension for restarted GMRES.
  mutable int max_dim;

  // Orthogonalization method for orthonormalizing a newly computed vector against a basis
  // at each iteration.
  Orthogonalization gs_orthog;

  // Use left or right preconditioning.
  PreconditionerSide pc_side;

  // Temporary workspace for solve.
  mutable std::vector<VecType> V;
  mutable VecType r;
  mutable std::vector<ScalarType> H;
  mutable std::vector<ScalarType> s, sn;
  mutable std::vector<RealType> cs;

  // Stop as soon as the correction to a nonzero initial guess exceeds this fraction of the
  // solution norm (no limit if zero), and whether the previous solve stopped for it.
  RealType max_rel_correction;
  mutable bool correction_exceeded;

  // With CGS2 (see Arnoldi), the Hessenberg matrix before the Givens rotations, the upper
  // triangular U such that each V[j] before its reorthogonalization and normalization is
  // Σₗ U[l, j] V[l] in the final basis, s[j] before the latest rotation, and the power of
  // two by which V[j + 1] was divided.
  mutable std::vector<ScalarType> Hraw, U;
  mutable ScalarType s_prev = 0.0;
  mutable RealType tau = 1.0;

  // Allocate storage for solve.
  virtual void Initialize() const;
  virtual void Update(int j) const;

  // Norm of the correction x_j - x₀ = V_j y_j in the current Krylov cycle, with y_j solving
  // the rotated least-squares system of iteration j. V_j is orthonormal (with CGS2, V[j] is
  // completed by Arnoldi(j)), so this is ‖y_j‖.
  RealType CorrectionNorm(int j) const;

  // Orthogonalize V[j + 1] = M V[j] against V[0], ..., V[j] and normalize it, with the
  // coefficients in column j of the Hessenberg matrix, then apply the Givens rotations to
  // it and to s. With CGS2, classical Gram-Schmidt with delayed reorthogonalization (DCGS2,
  // Świrydowicz et al., Numer. Linear Algebra Appl. 28 (2021), e2343) takes one global
  // reduction per iteration instead of two: V[j] and column j - 1 are completed at
  // iteration j, and V[j + 1] and column j after one pass give the residual norm.
  void Arnoldi(int j) const;

public:
  GmresSolver(MPI_Comm comm, int print)
    : IterativeSolver<OperType>(comm, print), max_dim(-1),
      gs_orthog(Orthogonalization::MGS), pc_side(PreconditionerSide::LEFT),
      max_rel_correction(0.0), correction_exceeded(false)
  {
  }

  // Set the dimension for restart.
  void SetRestartDim(int dim) { max_dim = dim; }

  // Set the orthogonalization method.
  void SetOrthogonalization(Orthogonalization orthog) { gs_orthog = orthog; }

  // Set the side for preconditioning.
  virtual void SetPreconditionerSide(PreconditionerSide side) { pc_side = side; }

  // Stop the solve, unconverged, as soon as the norm of the correction to a nonzero initial
  // guess x₀ exceeds tol times the norm of the iterate: ‖x_j - x₀‖ > tol (‖x₀‖ + ‖x_j -
  // x₀‖) (no limit if tol = 0). The correction norm is exact and free for left (or no)
  // preconditioning, where it is ‖y_j‖; it is not monitored with right preconditioning, by
  // FGMRES, or after a restart.
  void SetMaxRelativeCorrection(double tol) { max_rel_correction = tol; }

  // Whether the previous solve was stopped by the correction limit.
  bool GetCorrectionExceeded() const { return correction_exceeded; }

  void Mult(const VecType &b, VecType &x) const override;
};

// Preconditioned Flexible Generalized Minimum Residual Method (FGMRES) for general
// nonsymmetric linear systems with a non-constant preconditioner.
template <OperatorType OperType>
class FgmresSolver : public GmresSolver<OperType>
{
protected:
  using VecType = typename GmresSolver<OperType>::VecType;
  using RealType = typename GmresSolver<OperType>::RealType;
  using ScalarType = typename GmresSolver<OperType>::ScalarType;

  using GmresSolver<OperType>::comm;
  using GmresSolver<OperType>::print_opts;
  using GmresSolver<OperType>::int_width;
  using GmresSolver<OperType>::tab_width;

  using GmresSolver<OperType>::rel_tol;
  using GmresSolver<OperType>::abs_tol;
  using GmresSolver<OperType>::max_it;

  using GmresSolver<OperType>::A;
  using GmresSolver<OperType>::B;

  using GmresSolver<OperType>::converged;
  using GmresSolver<OperType>::initial_res;
  using GmresSolver<OperType>::final_res;
  using GmresSolver<OperType>::final_it;

  using GmresSolver<OperType>::max_dim;
  using GmresSolver<OperType>::gs_orthog;
  using GmresSolver<OperType>::pc_side;
  using GmresSolver<OperType>::V;
  using GmresSolver<OperType>::H;
  using GmresSolver<OperType>::s;
  using GmresSolver<OperType>::sn;
  using GmresSolver<OperType>::cs;
  using GmresSolver<OperType>::U;

  // Temporary workspace for solve.
  mutable std::vector<VecType> Z;

  // Allocate storage for solve.
  void Initialize() const override;
  void Update(int j) const override;

public:
  FgmresSolver(MPI_Comm comm, int print) : GmresSolver<OperType>(comm, print)
  {
    pc_side = PreconditionerSide::RIGHT;
  }

  void SetPreconditionerSide(const PreconditionerSide side) override
  {
    MFEM_VERIFY(side == PreconditionerSide::RIGHT,
                "FGMRES solver only supports right preconditioning!");
  }

  void Mult(const VecType &b, VecType &x) const override;
};

}  // namespace palace

#endif  // PALACE_LINALG_ITERATIVE_HPP
