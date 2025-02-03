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

  // Set operators for the generalized eigenvalue problem or for the quadratic polynomial
  // eigenvalue problem.
  virtual void SetOperators(const ComplexOperator &K, const ComplexOperator &M,
                            ScaleType type) = 0;
  virtual void SetOperators(const ComplexOperator &K, const ComplexOperator &C,
                            const ComplexOperator &M, ScaleType type) = 0;
  // For the linear generalized case, the linear solver should be configured to compute the
  // action of M⁻¹ (with no spectral transformation) or (K - σ M)⁻¹. For the quadratic
  // case, the linear solver should be configured to compute the action of M⁻¹ (with no
  // spectral transformation) or P(σ)⁻¹.
  virtual void SetLinearSolver(const ComplexKspSolver &ksp) = 0;

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
