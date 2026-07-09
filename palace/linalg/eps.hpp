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

namespace config
{

struct EigenSolverData;

}  // namespace config

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

// Construct and configure an eigenvalue solver from configuration: backend selection
// (with compile-time availability guards), problem-structure dispatch (linear EPS,
// quadratic PEP or linearized PEP, or NEP for the SLP nonlinear path), problem type,
// Gram-Schmidt orthogonalization variant, number of modes, tolerance, and maximum
// iterations. Does NOT call SetOperators or SetLinearSolver — the caller is responsible
// for configuring the pencil operators and the spectral-transformation solver, since
// these legitimately differ between consumers (the eigenmode solver driver may extend
// K/C/M with interpolated frequency-nonlinear A2 terms, while the driven-solver PROM
// enrichment uses its already-assembled bare K/C/M).
std::unique_ptr<EigenvalueSolver>
BuildEigenvalueSolver(MPI_Comm comm, const config::EigenSolverData &eigenmode,
                      Orthogonalization gs_orthog, int verbose, bool quadratic,
                      bool nonlinear_slp = false);

// Configure the shift-and-invert spectral transformation about the target frequency σ,
// with the eigenvalue conventions of the two pencils: λ = iω for the quadratic (and SLP
// nonlinear) problem, so the shift is iσ; μ = ω² for the linear problem, so the shift is
// σ². Also selects the backend-appropriate WhichType for eigenvalues near the target.
void SetEigenSolverShiftInvert(EigenvalueSolver &eigen, EigenSolverBackend backend,
                               double target, bool quadratic, bool nonlinear_slp = false);

// Convert a converged eigenvalue to complex angular frequency ω for either pencil
// convention (λ = iω quadratic, μ = ω² linear).
std::complex<double> EigenvalueToOmega(std::complex<double> lambda, bool quadratic);

}  // namespace palace

#endif  // PALACE_LINALG_EPS_HPP
