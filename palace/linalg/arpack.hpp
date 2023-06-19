// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_ARPACK_HPP
#define PALACE_LINALG_ARPACK_HPP

#if defined(PALACE_WITH_ARPACK)

#include "linalg/petsc.hpp"

#if !defined(PETSC_USE_COMPLEX)
#error "ARPACK interface requires PETSc built with complex scalars!"
#endif

#include <parpack.hpp>
#include "linalg/eigen.hpp"

namespace palace
{

class DivFreeSolver;
class KspSolver;

namespace arpack
{

//
// A wrapper for the ARPACK/PARPACK library for generalized linear eigenvalue problems or
// quadratic polynomial eigenvalue problems. Shift-and-invert spectral transformations are
// used to compute interior eigenvalues. Currently only implemented for complex scalar
// interface.
//
class ArpackEigenSolver : public EigenSolverBase
{
protected:
  // Control print level for debugging.
  int print;

  // Status variable for ARPACK.
  int info;

  // Number eigenvalues to be computed, and dimension.
  int nev, ncv;

  // Relative eigenvalue error convergence tolerance for the solver.
  double rtol;

  // Maximum number of Arnoldi update iterations.
  int max_it;

  // Specifies which part of the spectrum to search for.
  ::arpack::which which_option;

  // Variables for scaling, from Higham et al., IJNME 2008.
  double gamma, delta;

  // Parameters defining the spectral transformation.
  PetscScalar sigma;
  bool sinvert;

  // Storage for computed eigenvalues.
  PetscScalar *eig;
  int *perm;

  // Storage for Arnoldi basis vectors.
  petsc::PetscDenseMatrix *V;

  // Storage for computed residual norms.
  mutable double *res;

  // On input used to define optional initial guess, on output stores final residual
  // vector.
  mutable petsc::PetscParVector *r;

  // Reference to linear solver used for operator action for M⁻¹ (with no spectral
  // transformation) or (K - σ M)⁻¹ (generalized EVP with shift-and- invert) or P(σ)⁻¹
  // (polynomial with shift-and-invert) (not owned).
  const KspSolver *opInv;

  // Reference to solver for projecting an intermediate vector onto a divergence-free space
  // (not owned).
  const DivFreeSolver *opProj;

  // Reference to matrix used for weighted inner products (not owned). May be nullptr, in
  // which case identity is used.
  const petsc::PetscParMatrix *opB;

  // Perform the ARPACK RCI loop.
  int SolveInternal(petsc::PetscParVector &r_, petsc::PetscDenseMatrix &V_,
                    PetscScalar *eig_, int *perm_);

  // Helper routines for parameter checking.
  void CheckParameters() const;
  void CheckInfoAUPD(int info) const;
  void CheckInfoEUPD(int info) const;

  // Helper routines for ARPACK RCI.
  virtual void ApplyOp(const petsc::PetscParVector &x, petsc::PetscParVector &y) const = 0;
  virtual void ApplyOpB(const petsc::PetscParVector &x, petsc::PetscParVector &y) const = 0;

  // Helper routine for computing the eigenpair residual.
  virtual void GetResidual(PetscScalar l, const petsc::PetscParVector &x,
                           petsc::PetscParVector &r) const = 0;

  // Helper routine for computing the backward error.
  virtual double GetBackwardScaling(PetscScalar l) const = 0;

  // Return problem type name.
  virtual const char *GetName() const = 0;

public:
  ArpackEigenSolver(int print_lvl);
  ~ArpackEigenSolver() override;

  // Set operators for the generalized eigenvalue problem or for the quadratic polynomial
  // eigenvalue problem.
  void SetOperators(const petsc::PetscParMatrix &K, const petsc::PetscParMatrix &M,
                    ScaleType type) override;
  void SetOperators(const petsc::PetscParMatrix &K, const petsc::PetscParMatrix &C,
                    const petsc::PetscParMatrix &M, ScaleType type) override;

  // For the linear generalized case, the linear solver should be configured to compute the
  // action of M⁻¹ (with no spectral transformation) or (K - σ M)⁻¹. For the quadratic
  // case, the linear solver should be configured to compute the action of M⁻¹ (with no
  // spectral transformation) or P(σ)⁻¹.
  void SetLinearSolver(const KspSolver &ksp) override;

  // Set the projection operator for the divergence-free constraint.
  void SetProjector(const DivFreeSolver &divfree) override;

  // Set optional B matrix used for weighted inner products. This must be set explicitly
  // even for generalized problems, otherwise the identity will be used.
  void SetBMat(const petsc::PetscParMatrix &B) override;

  // Get scaling factors used by the solver.
  double GetScalingGamma() const override { return gamma; }
  double GetScalingDelta() const override { return delta; }

  // Set the number of required eigenmodes.
  void SetNumModes(int numeig, int numvec = 0) override;

  // Set solver tolerance.
  void SetTol(double tol) override;

  // Set maximum number of Arnoldi update iterations.
  void SetMaxIter(int maxits) override;

  // Set target spectrum for the eigensolver. When a spectral transformation is used, this
  // applies to the spectrum of the shifted operator.
  void SetWhichEigenpairs(WhichType type) override;

  // Set shift-and-invert spectral transformation.
  void SetShiftInvert(double tr, double ti, bool precond = false) override;

  // Set an initial vector for the solution subspace.
  void SetInitialSpace(const petsc::PetscParVector &v) override;

  // Solve the eigenvalue problem. Returns the number of converged eigenvalues.
  int Solve() override = 0;

  // Get the corresponding eigenvalue.
  void GetEigenvalue(int i, double &eigr, double &eigi) const override;

  // Get the corresponding eigenvector.
  void GetEigenvector(int i, petsc::PetscParVector &x) const override;

  // Get the corresponding eigenpair error.
  void GetError(int i, ErrorType type, double &err) const override;
};

// Generalized eigenvalue problem solver: K x = λ M x .
class ArpackEPSSolver : public ArpackEigenSolver
{
private:
  // References to matrices defining the generalized eigenvalue problem (not owned).
  const petsc::PetscParMatrix *opK, *opM;

  // Operator norms for scaling.
  mutable double normK, normM;

  // Workspace vector for operator applications.
  mutable petsc::PetscParVector *z;

protected:
  // Helper routines for ARPACK RCI interface.
  void ApplyOp(const petsc::PetscParVector &x, petsc::PetscParVector &y) const override;
  void ApplyOpB(const petsc::PetscParVector &x, petsc::PetscParVector &y) const override;

  // Helper routine for computing the eigenpair residual: r = (K - λ M) x .
  void GetResidual(PetscScalar l, const petsc::PetscParVector &x,
                   petsc::PetscParVector &r) const override;

  // Helper routine for computing the backward error.
  double GetBackwardScaling(PetscScalar l) const override;

  // Return problem type name.
  const char *GetName() const override { return "EPS"; }

public:
  ArpackEPSSolver(int print_lvl);
  ~ArpackEPSSolver() override;

  // Set operators for the generalized eigenvalue problem.
  void SetOperators(const petsc::PetscParMatrix &K, const petsc::PetscParMatrix &M,
                    ScaleType type) override;

  // Solve the eigenvalue problem. Returns the number of converged eigenvalues.
  int Solve() override;
};

// Quadratic eigenvalue problem solver: P(λ) x = (K + λ C + λ² M) x = 0 .
class ArpackPEPSolver : public ArpackEigenSolver
{
private:
  // References to matrices defining the quadratic polynomial eigenvalue problem
  // (not owned).
  const petsc::PetscParMatrix *opK, *opC, *opM;

  // Operator norms for scaling.
  mutable double normK, normC, normM;

  // Workspace vectors for operator applications.
  mutable petsc::PetscParVector *x1, *x2, *y1, *y2, *z;

  // Do eigenvector extraction from the linearized problem to the actual eigenvectors.
  void ExtractEigenvector(PetscScalar l, petsc::PetscParVector &y,
                          petsc::PetscParVector &x);

  // Helper methods for splitting a block vector from the linearized problem into its into
  // two parts.
  PetscScalar *GetBlocks(petsc::PetscParVector &v, petsc::PetscParVector &v1,
                         petsc::PetscParVector &v2) const;
  void RestoreBlocks(PetscScalar *pv, petsc::PetscParVector &v, petsc::PetscParVector &v1,
                     petsc::PetscParVector &v2) const;

protected:
  // Helper routines for ARPACK RCI interface.
  void ApplyOp(const petsc::PetscParVector &x, petsc::PetscParVector &y) const override;
  void ApplyOpB(const petsc::PetscParVector &x, petsc::PetscParVector &y) const override;

  // Helper routine for computing the eigenpair residual: r = P(λ) x .
  void GetResidual(PetscScalar l, const petsc::PetscParVector &x,
                   petsc::PetscParVector &r) const override;

  // Helper routine for computing the backward error.
  double GetBackwardScaling(PetscScalar l) const override;

  // Return problem type name.
  const char *GetName() const override { return "PEP"; }

public:
  ArpackPEPSolver(int print_lvl);
  ~ArpackPEPSolver() override;

  // Set operators for the quadratic polynomial eigenvalue problem.
  void SetOperators(const petsc::PetscParMatrix &K, const petsc::PetscParMatrix &C,
                    const petsc::PetscParMatrix &M, ScaleType type) override;

  // Solve the eigenvalue problem. Returns the number of converged eigenvalues.
  int Solve() override;
};

}  // namespace arpack

}  // namespace palace

#endif

#endif  // PALACE_LINALG_ARPACK_HPP
