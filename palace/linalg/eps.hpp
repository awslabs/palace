// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_EPS_HPP
#define PALACE_LINALG_EPS_HPP

#include <complex>
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

template <typename VecType>
class DivFreeSolver;

//
// Pure abstract base class for solving generalized linear eigenvalue problems problems or
// quadratic polynomial eigenvalue problems.
//
class EigenvalueSolver
{
public:
  enum class ScaleType
  {
    NONE,
    NORM_2
  };

  enum class WhichType
  {
    LARGEST_MAGNITUDE,
    SMALLEST_MAGNITUDE,
    LARGEST_REAL,
    SMALLEST_REAL,
    LARGEST_IMAGINARY,
    SMALLEST_IMAGINARY,
    TARGET_MAGNITUDE,
    TARGET_REAL,
    TARGET_IMAGINARY
  };

  enum class ErrorType
  {
    ABSOLUTE,
    RELATIVE,
    BACKWARD
  };

public:
  EigenvalueSolver() = default;
  virtual ~EigenvalueSolver() = default;

  // Set operators for the generalized eigenvalue problem, quadratic polynomial
  // eigenvalue problem, or nonlinear eigenvalue problem.
  virtual void SetOperators(const ComplexOperator &K, const ComplexOperator &M,
                            ScaleType type)
  {
    MFEM_ABORT("SetOperators not defined!");
  }

  virtual void SetOperators(const ComplexOperator &K, const ComplexOperator &C,
                            const ComplexOperator &M, ScaleType type)
  {
    MFEM_ABORT("SetOperators not defined!");
  }

  virtual void SetOperators(const ComplexOperator &K, const ComplexOperator &M,
                            std::function<const ComplexOperator &(std::complex<double>)> A2,
                            ScaleType type)
  {
    MFEM_ABORT("SetOperators not defined!");
  }

  // Set the frequency-dependent A2 matrix function A2(λ) for the nonlinear eigensolvers.
  virtual void SetExtraSystemMatrix(
      std::function<std::unique_ptr<ComplexOperator>(std::complex<double>)>)
  {
    MFEM_ABORT("SetExtraSystemMatrix not defined!");
  }

  virtual void SetPreconditionerUpdate(std::function<std::unique_ptr<ComplexOperator>(
                                           std::complex<double>, std::complex<double>,
                                           std::complex<double>, std::complex<double>)>)
  {
    MFEM_ABORT("SetPreconditionerUpdate not defined!");
  }

  // For the linear generalized case, the linear solver should be configured to compute the
  // action of M⁻¹ (with no spectral transformation) or (K - σ M)⁻¹. For the quadratic
  // case, the linear solver should be configured to compute the action of M⁻¹ (with no
  // spectral transformation) or P(σ)⁻¹.
  virtual void SetLinearSolver(ComplexKspSolver &ksp) = 0;

  // Set the projection operator for enforcing the divergence-free constraint.
  virtual void SetDivFreeProjector(const DivFreeSolver<ComplexVector> &divfree) = 0;

  // Set optional B matrix used for weighted inner products. This must be set explicitly
  // even for generalized problems, otherwise the identity will be used.
  virtual void SetBMat(const Operator &B) = 0;

  // Get scaling factors used by the solver.
  virtual double GetScalingGamma() const = 0;
  virtual double GetScalingDelta() const = 0;

  // Set the number of required eigenmodes.
  virtual void SetNumModes(int num_eig, int num_vec = 0) = 0;

  // Set solver tolerance.
  virtual void SetTol(double tol) = 0;

  // Set maximum number of Arnoldi update iterations.
  virtual void SetMaxIter(int max_it) = 0;

  // Set target spectrum for the eigensolver. When a spectral transformation is used, this
  // applies to the spectrum of the shifted operator.
  virtual void SetWhichEigenpairs(WhichType type) = 0;

  // Set shift-and-invert spectral transformation.
  virtual void SetShiftInvert(std::complex<double> s, bool precond = false) = 0;

  // Set an initial vector for the solution subspace.
  virtual void SetInitialSpace(const ComplexVector &v) = 0;

  // Set an initial vector for a quadratic eigenvalue problem solved by linearization, with
  // eigenvectors [x; λ x]: v combines eigenvector approximations xᵢ, and v_lambda the same
  // combination of λᵢ xᵢ. Other solvers use v only.
  virtual void SetInitialSpaceBlocks(const ComplexVector &v, const ComplexVector &v_lambda)
  {
    SetInitialSpace(v);
  }

  // Relative tolerance of the linear solves of the first Krylov cycle, used if positive and
  // smaller than that of the linear solver. Not all solvers support it.
  virtual void SetFirstCycleLinearSolverTol(double tol) {}

  // Whether the initial vector of the next solve approximates the wanted eigenvectors (warm
  // start). Not all solvers use it.
  virtual void SetWarmStart(bool warm) {}

  // Number of iterations (restarts) of the last solve, or -1 if not available.
  virtual int GetNumIterations() const { return -1; }

  // Restrict the search to modes with quality factor |Q| = |ω| / (2 |Im{ω}|) >= q_min
  // (eigenvalues λ = iω of a quadratic problem), so that strongly damped or growing Ritz
  // values cannot stall the convergence of the others. An eigenpair then converges only
  // once its Ritz value was also present at the previous restart, except in the first
  // Krylov cycle of a solve which is not warm-started. A q_min <= 0.5 disables the
  // restriction. Not all solvers support it.
  virtual void SetMinQualityFactor(double q_min) {}

  // Solve the eigenvalue problem. Returns the number of converged eigenvalues.
  virtual int Solve() = 0;

  // Get the corresponding eigenvalue.
  virtual std::complex<double> GetEigenvalue(int i) const = 0;

  // Get the corresponding eigenvector. Eigenvectors are normalized such that ||x||₂ = 1,
  // unless the B-matrix is set for weighted inner products.
  virtual void GetEigenvector(int i, ComplexVector &x) const = 0;

  // Get the corresponding eigenpair error.
  virtual double GetError(int i, ErrorType type) const = 0;

  // Re-normalize the given number of eigenvectors, for example if the matrix B for weighted
  // inner products has changed. This does not perform re-orthogonalization with respect to
  // the new matrix, only normalization.
  virtual void RescaleEigenvectors(int num_eig) = 0;
};

}  // namespace palace

#endif  // PALACE_LINALG_EPS_HPP
