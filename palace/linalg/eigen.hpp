// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_EIGEN_HPP
#define PALACE_LINALG_EIGEN_HPP

namespace palace
{

class DivFreeSolver;
class KspSolver;

namespace petsc
{

class PetscParMatrix;
class PetscParVector;

}  // namespace petsc

//
// Pure abstract base class for solving generalized linear eigenvalue problems problems or
// quadratic polynomial eigenvalue problems.
//
class EigenSolverBase
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

  EigenSolverBase() = default;
  virtual ~EigenSolverBase() = default;

  // Set operators for the generalized eigenvalue problem or for the quadratic polynomial
  // eigenvalue problem.
  virtual void SetOperators(const petsc::PetscParMatrix &K, const petsc::PetscParMatrix &M,
                            ScaleType type) = 0;
  virtual void SetOperators(const petsc::PetscParMatrix &K, const petsc::PetscParMatrix &C,
                            const petsc::PetscParMatrix &M, ScaleType type) = 0;

  // For the linear generalized case, the linear solver should be configured to compute the
  // action of M⁻¹ (with no spectral transformation) or (K - σ M)⁻¹. For the quadratic
  // case, the linear solver should be configured to compute the action of M⁻¹ (with no
  // spectral transformation) or P(σ)⁻¹.
  virtual void SetLinearSolver(const KspSolver &ksp) = 0;

  // Set the projection operator for the divergence-free constraint.
  virtual void SetProjector(const DivFreeSolver &divfree) = 0;

  // Get scaling factors used by the solver.
  virtual double GetScalingGamma() const = 0;
  virtual double GetScalingDelta() const = 0;

  // Set the number of required eigenmodes.
  virtual void SetNumModes(int numeig, int numvec = 0) = 0;

  // Set solver tolerance.
  virtual void SetTol(double tol) = 0;

  // Set maximum number of Arnoldi update iterations.
  virtual void SetMaxIter(int maxits) = 0;

  // Set target spectrum for the eigensolver. When a spectral transformation is used, this
  // applies to the spectrum of the shifted operator.
  virtual void SetWhichEigenpairs(WhichType type) = 0;

  // Set shift-and-invert spectral transformation.
  virtual void SetShiftInvert(double tr, double ti, bool precond = false) = 0;

  // Set optional B matrix used for weighted inner products. This must be set explicitly
  // even for generalized problems, otherwise the identity will be used.
  virtual void SetBMat(const petsc::PetscParMatrix &B) = 0;

  // Set an initial vector for the solution subspace.
  virtual void SetInitialSpace(const petsc::PetscParVector &v) = 0;

  // Solve the eigenvalue problem. Returns the number of converged eigenvalues.
  virtual int Solve() = 0;

  // Get the corresponding eigenvalue.
  virtual void GetEigenvalue(int i, double &eigr, double &eigi) const = 0;

  // Get the corresponding eigenvector.
  virtual void GetEigenvector(int i, petsc::PetscParVector &x) const = 0;

  // Get the corresponding eigenpair error.
  virtual void GetError(int i, ErrorType type, double &err) const = 0;
};

}  // namespace palace

#endif  // PALACE_LINALG_EIGEN_HPP
