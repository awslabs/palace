// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_EPS_HPP
#define PALACE_LINALG_EPS_HPP

#include <algorithm>
#include <complex>
#include <vector>
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

  // Set the target of the eigenvalue search (sorting criterion and search region) apart
  // from the shift of the shift-and-invert spectral transformation (SetShiftInvert sets
  // both). Only solvers for which SupportsTarget() is true support a target different from
  // the shift.
  virtual bool SupportsTarget() const { return false; }
  virtual void SetTarget(std::complex<double> t) {}

  // Set the initial space of the last call to SetInitialSpace or SetInitialSpaceBlocks
  // again, to repeat a solve. Not all solvers support it.
  virtual void ResetInitialSpace() {}

  // Set an initial vector for the solution subspace.
  virtual void SetInitialSpace(const ComplexVector &v) = 0;

  // Set an initial vector for the solution subspace of a quadratic eigenvalue problem
  // solved by linearization, whose eigenvectors are [x; λ x]: v is a combination of
  // eigenvector approximations xᵢ and v_lambda the same combination of λᵢ xᵢ. Solvers
  // which do not work on a linearization use v only.
  virtual void SetInitialSpaceBlocks(const ComplexVector &v, const ComplexVector &v_lambda)
  {
    SetInitialSpace(v);
  }

  // Relax the relative tolerance of the linear solves of the spectral transformation as
  // the wanted eigenpairs converge (inexact shift-and-invert), with the constant c of
  // GetRelaxedLinearSolverTol. A c <= 0 disables the relaxation. Not all solvers support
  // it.
  virtual void SetLinearSolverRelaxation(double c) {}

  // Relative tolerance of the linear solves of the first Krylov cycle of the next solve (a
  // tol <= 0 keeps the linear solver tolerance). With a start vector which already
  // approximates the wanted eigenvectors closely (warm start), the eigenpairs which
  // converge in the first cycle have Ritz values set by the first linear solves, whose
  // error is then of the order of their tolerance. Not all solvers support it.
  virtual void SetFirstCycleLinearSolverTol(double tol) {}

  // Whether the start vector of the next solve already approximates the wanted
  // eigenvectors (warm start): the eigenpairs which converge in its first Krylov cycle are
  // then accepted only after one more restart (see the stable Ritz value requirement of the
  // convergence test). Not all solvers support it.
  virtual void SetWarmStart(bool warm) {}

  // Test convergence every k steps of the Krylov factorization of each cycle, and end the
  // cycle once the requested eigenpairs have converged and their Ritz values have settled
  // (k = 0: at the end of each cycle only). This saves operator applications when they are
  // expensive, as the linear solves of shift-and-invert with iterative solvers. Not all
  // solvers support it.
  virtual void SetConvergenceCheckInterval(int k) {}

  // Estimates of the wanted eigenvalues of the next solve (as returned by GetEigenvalue),
  // for example those of the previous adaptation iteration. With shift-and-invert, they
  // set the starting tolerance of the divergence-free projection of the linear solve
  // results. Not all solvers use them.
  virtual void SetEigenvalueEstimates(const std::vector<std::complex<double>> &lambda) {}

  // Number of iterations (restarts) of the last solve, or -1 if not available.
  virtual int GetNumIterations() const { return -1; }

  // Restrict the eigenvalues searched for to quality factors Q = |ω| / (2 Im{ω}) >= q_min
  // for damped modes, and |Q| >= 10 for growing ones (eigenvalue λ = iω of a quadratic or
  // nonlinear problem with damping), so that strongly damped Ritz values cannot stall the
  // convergence of the wanted eigenpairs. A q_min <= 0.5 disables the restriction. Not all
  // solvers support it.
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

// Relative tolerance for the linear solves of an inexact shift-and-invert eigenvalue
// solver, for a wanted unconverged eigenpair with relative error estimate r (r <= 0 when
// none is left), eigenvalue solver tolerance eig_tol, and linear solver tolerance lin_tol:
// a new basis vector enters the wanted Ritz vector with a weight of the order of its
// residual divided by the separation of its Ritz value from the others, so the error of
// the linear solves can grow as the residual decreases (Simoncini, SINUM 2005; Freitag and
// Spence, SIMAX 2009). The error of a linear solve enters the residual of an eigenpair
// with transformed eigenvalue θ amplified by about θ_max / |θ| (theta_ratio, for the
// largest transformed eigenvalue θ_max) and divided by the relative separation gap of θ
// from the other Ritz values, so the tolerance is c eig_tol / r with c reduced to
// 0.01 gap / theta_ratio, bounded by lin_tol from below and by max(lin_tol, 1e-4) from
// above. Taking the minimum over the wanted unconverged eigenpairs keeps their residuals
// within a small fraction of eig_tol of their error estimates over the restarts.
inline double GetRelaxedLinearSolverTol(double c, double eig_tol, double lin_tol, double r,
                                        double theta_ratio = 1.0, double gap = 1.0)
{
  constexpr double max_tol = 1.0e-4, max_amplified_error = 1.0e-2;
  if (c <= 0.0 || r <= 0.0)
  {
    return lin_tol;
  }
  c = std::min(c, max_amplified_error * std::min(gap, 1.0) / std::max(theta_ratio, 1.0));
  return std::max(lin_tol, std::min(max_tol, c * eig_tol / r));
}

// Shift s of the shift-and-invert spectral transformation for the search of the n
// eigenvalues λ = iω above the target τ > 0 (Im{λ} >= τ, nearest to τ first), from
// estimates λᵢ of them (the n wanted ones first, further ones may follow), for example the
// converged eigenvalues of the previous adaptation iteration. With ωₙ the largest Im{λᵢ}
// of the wanted estimates and s <= (τ + ωₙ) / 2, every mode with τ <= ω <= ωₙ is within
// ωₙ - s of s, so the wanted modes are the dominant ones of the transformed operator as
// with s = τ, while the modes below the target are farther from s than all of them. The
// shift is s = τ + f (ωₙ - τ), f <= 1/2, kept clear of the estimates of modes with Q >= 10
// (and of τ): within a quarter of the interval between consecutive ones from one of them,
// it moves down to the middle of that interval, or of the one below, and it moves further
// down to the first interval below at least half as wide as the mean spacing (ωₙ - τ) / n
// of the wanted modes (clusters of close modes make the shifted linear systems harder).
// Returns τ without estimates of every wanted mode, or without such an interval.
inline double GetEigenvalueShift(double target,
                                 const std::vector<std::complex<double>> &lambda, int n,
                                 double f = 0.5)
{
  if (target <= 0.0 || n <= 0 || static_cast<int>(lambda.size()) < n)
  {
    return target;
  }
  double w_n = target;
  for (int i = 0; i < n; i++)
  {
    w_n = std::max(w_n, lambda[i].imag());
  }
  std::vector<double> w = {target};
  for (const auto &l : lambda)
  {
    // Q = |λ| / (-2 Re{λ}) for λ = iω.
    if (l.imag() > target && std::abs(l) >= -20.0 * l.real())
    {
      w.push_back(l.imag());
    }
  }
  std::sort(w.begin(), w.end());
  w.erase(std::unique(w.begin(), w.end(),
                      [](double x, double y) { return y - x <= 1.0e-9 * std::abs(y); }),
          w.end());  // Degenerate modes
  double s = target + std::clamp(f, 0.0, 0.5) * (w_n - target);
  const auto it = std::upper_bound(w.begin(), w.end(), s);  // First estimate above s
  if (it == w.end() || it == w.begin())
  {
    return target;
  }
  // Interval [w[j], w[j + 1]] of the shift.
  auto j = static_cast<std::size_t>(it - w.begin()) - 1;
  if (w[j + 1] - s < 0.25 * (w[j + 1] - w[j]))
  {
    s = 0.5 * (w[j] + w[j + 1]);
  }
  else if (s - w[j] < 0.25 * (w[j + 1] - w[j]))
  {
    if (j == 0)
    {
      return target;
    }
    j--;
    s = 0.5 * (w[j] + w[j + 1]);
  }
  const double min_width = 0.5 * (w_n - target) / n;
  while (w[j + 1] - w[j] < min_width)
  {
    if (j == 0)
    {
      return target;
    }
    j--;
    s = 0.5 * (w[j] + w[j + 1]);
  }
  return s;
}

// Check of an eigenvalue solve with the shift s apart from the target τ (see
// GetEigenvalueShift), from its eigenvalues λᵢ (the n requested ones first): at least n
// converged, and s at most in the middle of τ and the farthest of the n requested ones, so
// that every mode nearer to the target was at least as dominant in the transformed operator
// as the accepted ones, as with the shift at the target.
inline bool CheckEigenvalueShift(double target, double s,
                                 const std::vector<std::complex<double>> &lambda, int n)
{
  if (static_cast<int>(lambda.size()) < n || n <= 0)
  {
    return false;
  }
  double w_n = target;
  for (int i = 0; i < n; i++)
  {
    w_n = std::max(w_n, lambda[i].imag());
  }
  return s <= 0.5 * (target + w_n);
}

}  // namespace palace

#endif  // PALACE_LINALG_EPS_HPP
