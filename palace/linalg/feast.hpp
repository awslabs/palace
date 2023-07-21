// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_FEAST_HPP
#define PALACE_LINALG_FEAST_HPP

#if defined(PALACE_WITH_SLEPC)

#include "linalg/petsc.hpp"

#if !defined(PETSC_USE_COMPLEX)
#error "FEAST eigenvalue solver requires PETSc built with complex scalars!"
#endif

#include <mfem.hpp>
#include "linalg/eigen.hpp"

// Forward declarations of SLEPc objects.
typedef struct _p_RG *RG;

namespace palace
{

class DivFreeSolver;
class IoData;
class KspSolver;
class SpaceOperator;

namespace feast
{

namespace internal
{

class FeastLinearSolver;

}  // namespace internal

//
// A native implementation of the FEAST eigenvalue solver, with residual- inverse iteration
// for linear and quadratic eigenvalue problems with complex- symmetric matrices. Depends
// on SLEPc for some functionality like quadrature and solving projected the eigenvalue
// problem.
//
class FeastEigenSolver : public EigenSolverBase
{
protected:
  // Control print level for debugging.
  int print;

  // Status variable.
  int info;

  // Number eigenvalues to be computed. Also the subspace and projected system dimensions.
  PetscInt nev, m0, mQ;

  // Number of moments to consider for subspace construction.
  PetscInt N_moments;

  // Relative eigenvalue error convergence tolerance for the solver.
  PetscReal rtol;

  // Maximum number of FEAST iterations.
  PetscInt max_it;

  // Variables for scaling, from Higham et al., IJNME 2008.
  PetscReal gamma, delta;

  // Parameters defining the integration contour.
  PetscScalar bl, tr;
  bool real_threshold, imag_threshold;

  // Storage for computed eigenvalues.
  PetscScalar *eig;
  PetscInt *perm;

  // Storage for computed eigenvectors.
  petsc::PetscDenseMatrix *X;

  // Storage for computed residual norms.
  mutable PetscReal *res;

  // Workspace vector for initial space and residual calculations.
  mutable petsc::PetscParVector *r0;

  // Solvers for Linear systems for the different quadrature points along the contour.
  std::vector<internal::FeastLinearSolver> opInv;

  // Reference to solver for projecting an intermediate vector onto a divergence-free space
  // (not owned).
  const DivFreeSolver *opProj;

  // Reference to matrix used for weighted inner products (not owned). May be nullptr, in
  // which case identity is used.
  const petsc::PetscParMatrix *opB;

  // Perform the FEAST solve.
  int SolveInternal(RG rg);

  // Helper routine for parameter checking.
  void CheckParameters();

  // Helper routine to construct the integration contour.
  RG ConfigureRG(PetscScalar *&z, PetscScalar *&w);

  // Helper routine for sorting eigenvalues of the projected problem.
  PetscInt *SortEigenvalues(const PetscScalar *eig_, PetscInt m) const;

  // Helper routine for computing the Qᴴ A Q matrix product.
  void BVMatProjectInternal(const petsc::PetscDenseMatrix &Q,
                            const petsc::PetscParMatrix &A, petsc::PetscDenseMatrix &Ar,
                            petsc::PetscDenseMatrix &R, PetscReal scale) const;

  // Helper routine for solving the projected eigenvalue problem.
  virtual void SolveProjectedProblem(const petsc::PetscDenseMatrix &Q_,
                                     petsc::PetscDenseMatrix &R_,
                                     petsc::PetscDenseMatrix &X_, PetscScalar *eig_) = 0;

  // Helper routine for computing the eigenpair residual.
  virtual void GetResidual(PetscScalar eig_, const petsc::PetscParVector &x_,
                           petsc::PetscParVector &r_) const = 0;

  // Helper routine for computing the backward error.
  virtual PetscReal GetBackwardScaling(PetscScalar eig_) const = 0;

  // Return problem type name.
  virtual const char *GetName() const = 0;

public:
  FeastEigenSolver(MPI_Comm comm, const IoData &iodata, SpaceOperator &spaceop, int np,
                   int print_lvl);
  ~FeastEigenSolver() override;

  // Set operators for the generalized eigenvalue problem or for the quadratic polynomial
  // eigenvalue problem.
  void SetOperators(const petsc::PetscParMatrix &K, const petsc::PetscParMatrix &M,
                    ScaleType type) override;
  void SetOperators(const petsc::PetscParMatrix &K, const petsc::PetscParMatrix &C,
                    const petsc::PetscParMatrix &M, ScaleType type) override;

  // Set the projection operator for the divergence-free constraint.
  void SetProjector(const DivFreeSolver &divfree) override;

  // Get scaling factors used by the solver.
  double GetScalingGamma() const override { return (double)gamma; }
  double GetScalingDelta() const override { return (double)delta; }

  // Set the number of required eigenmodes.
  void SetNumModes(int numeig, int numvec = 0) override;

  // Set solver tolerance.
  void SetTol(double tol) override;

  // Set maximum number of FEAST iterations.
  void SetMaxIter(int maxits) override;

  // Set up region for contour integration. Region is defined by the bottom- left and
  // top-right points in the complex plane.
  void SetContour(double blr, double bli, double trr, double tri,
                  bool filter_small_real = false, bool filter_small_imag = false);

  // Set optional B matrix used for weighted inner products. This must be set explicitly
  // even for generalized problems, otherwise the identity will be used.
  void SetBMat(const petsc::PetscParMatrix &B) override;

  // Set an initial vector for the solution subspace.
  void SetInitialSpace(const petsc::PetscParVector &v) override;

  // Solve the eigenvalue problem. Returns the number of converged eigenvalues.
  int Solve() override = 0;

  // Return number of linear solves and linear solver iterations performed during the FEAST
  // solve.
  int GetTotalKspMult() const;
  int GetTotalKspIter() const;

  // Get the corresponding eigenvalue.
  void GetEigenvalue(int i, double &eigr, double &eigi) const override;

  // Get the corresponding eigenvector.
  void GetEigenvector(int i, petsc::PetscParVector &v) const override;

  // Get the corresponding eigenpair error.
  void GetError(int i, ErrorType type, double &err) const override;

  // Methods not relevant to the FEAST eigenvalue solver.
  void SetLinearSolver(const KspSolver &ksp) override
  {
    MFEM_ABORT("SetLinearSolver not defined for FeastEigenSolver!");
  }
  void SetWhichEigenpairs(WhichType type) override
  {
    MFEM_ABORT("SetWhichEigenpairs not defined for FeastEigenSolver!");
  }
  void SetShiftInvert(double tr_in, double ti_in, bool precond = false) override
  {
    MFEM_ABORT("SetShiftInvert not defined for FeastEigenSolver!");
  }
};

// Generalized eigenvalue problem solver: K x = λ M x .
class FeastEPSSolver : public FeastEigenSolver
{
private:
  // References to matrices defining the generalized eigenvalue problem (not owned).
  const petsc::PetscParMatrix *opK, *opM;

  // Operator norms for scaling.
  mutable PetscReal normK, normM;

  // Sequential workspace matrices for projected problem.
  petsc::PetscDenseMatrix *AQ, *BQ, *XQ, *XQ0;

protected:
  // Helper routine for solving the projected eigenvalue problem.
  void SolveProjectedProblem(const petsc::PetscDenseMatrix &Q_, petsc::PetscDenseMatrix &R_,
                             petsc::PetscDenseMatrix &X_, PetscScalar *eig_) override;

  // Helper routine for computing the eigenpair residuals: R = K X - M X Λ .
  void GetResidual(PetscScalar eig_, const petsc::PetscParVector &x_,
                   petsc::PetscParVector &r_) const override;

  // Helper routine for computing the backward error.
  PetscReal GetBackwardScaling(PetscScalar eig_) const override;

  // Return problem type name.
  const char *GetName() const override { return "EPS"; }

public:
  FeastEPSSolver(MPI_Comm comm, const IoData &iodata, SpaceOperator &spaceop, int np,
                 int print_lvl);

  // Set operators for the generalized eigenvalue problem.
  void SetOperators(const petsc::PetscParMatrix &K, const petsc::PetscParMatrix &M,
                    ScaleType type) override;

  // Solve the eigenvalue problem. Returns the number of converged eigenvalues.
  int Solve() override;
};

// Quadratic eigenvalue problem solver: P(λ) x = (K + λ C + λ² M) x = 0 .
class FeastPEPSolver : public FeastEigenSolver
{
private:
  // References to matrices defining the quadratic eigenvalue problem (not owned).
  const petsc::PetscParMatrix *opK, *opC, *opM;

  // Operator norms for scaling.
  mutable PetscReal normK, normC, normM;

  // Sequential workspace matrices for projected problem.
  petsc::PetscDenseMatrix *AQ, *BQ, *AQ0, *XQ, *XQ0;

protected:
  // Helper routine for solving the projected eigenvalue problem.
  void SolveProjectedProblem(const petsc::PetscDenseMatrix &Q_, petsc::PetscDenseMatrix &R_,
                             petsc::PetscDenseMatrix &X_, PetscScalar *eig_) override;

  // Helper routine for computing the eigenpair residuals: R = P(Λ, X) .
  void GetResidual(PetscScalar eig_, const petsc::PetscParVector &x_,
                   petsc::PetscParVector &r_) const override;

  // Helper routine for computing the backward error.
  PetscReal GetBackwardScaling(PetscScalar eig_) const override;

  // Return problem type name.
  const char *GetName() const override { return "PEP"; }

public:
  FeastPEPSolver(MPI_Comm comm, const IoData &iodata, SpaceOperator &spaceop, int np,
                 int print_lvl);

  // Set operators for the quadratic polynomial eigenvalue problem.
  void SetOperators(const petsc::PetscParMatrix &K, const petsc::PetscParMatrix &C,
                    const petsc::PetscParMatrix &M, ScaleType type) override;

  // Solve the eigenvalue problem. Returns the number of converged eigenvalues.
  int Solve() override;
};

}  // namespace feast

}  // namespace palace

#endif

#endif  // PALACE_LINALG_FEAST_HPP
