// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_SLEPC_HPP
#define PALACE_LINALG_SLEPC_HPP

#if defined(PALACE_WITH_SLEPC)

#include "linalg/petsc.hpp"

#if !defined(PETSC_USE_COMPLEX)
#error "SLEPc interface requires PETSc compiled with complex scalars!"
#endif

#include <complex>
#include <memory>
#include <string>
#include <mpi.h>
#include "linalg/eps.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

// Forward declarations of SLEPc objects.
typedef struct _p_EPS *EPS;
typedef struct _p_PEP *PEP;
typedef struct _p_BV *BV;
typedef struct _p_ST *ST;
typedef struct _p_RG *RG;

namespace palace
{

namespace slepc
{

// Wrappers for SlepcInitialize/SlepcInitializeNoArguments/SlepcFinalize.
void Initialize(int &argc, char **&argv, const char rc_file[], const char help[]);
void Initialize();
void Finalize();

// Compute and return the maximum singular value of the given operator, σₙ² = λₙ(Aᴴ A) .
PetscReal GetMaxSingularValue(MPI_Comm comm, const ComplexOperator &A, bool herm = false,
                              PetscReal tol = PETSC_DEFAULT,
                              PetscInt max_it = PETSC_DEFAULT);

//
// A wrapper for the SLEPc library for generalized linear eigenvalue problems or quadratic
// polynomial eigenvalue problems. Shift-and-invert spectral transformations can be used to
// compute interior eigenvalues.
//
class SlepcEigenvalueSolver : public EigenvalueSolver
{
public:
  enum class ProblemType
  {
    HERMITIAN,
    NON_HERMITIAN,
    GEN_HERMITIAN,
    GEN_INDEFINITE,
    GEN_NON_HERMITIAN,
    HYPERBOLIC,
    GYROSCOPIC
  };

  enum class Type
  {
    KRYLOVSCHUR,
    POWER,
    SUBSPACE,
    TOAR,
    STOAR,
    QARNOLDI,
    JD
  };

  // Workspace vector for operator applications.
  mutable ComplexVector x1, y1;

protected:
  // Control print level for debugging.
  int print;

  // Variables for scaling, from Higham et al., IJNME 2008.
  PetscReal gamma, delta;

  // Parameters defining the spectral transformation.
  PetscScalar sigma;
  bool sinvert, region;

  // Storage for computed residual norms and eigenvector normalizations.
  std::unique_ptr<PetscReal[]> res, xscale;

  // Reference to linear solver used for operator action for M⁻¹ (with no spectral
  // transformation) or (K - σ M)⁻¹ (generalized EVP with shift-and- invert) or P(σ)⁻¹
  // (polynomial with shift-and-invert) (not owned).
  const ComplexKspSolver *opInv;

  // Reference to solver for projecting an intermediate vector onto a divergence-free space
  // (not owned).
  const DivFreeSolver<ComplexVector> *opProj;

  // Reference to matrix used for weighted inner products (not owned). May be nullptr, in
  // which case identity is used.
  const Operator *opB;

  // Workspace objects for eigenvalue calculations.
  Mat B0;
  Vec v0;

  // Boolean to handle SetFromOptions calls.
  mutable bool cl_custom;

  // Customize object with command line options set.
  virtual void Customize();

  // Helper routine for computing the eigenvector normalization.
  PetscReal GetEigenvectorNorm(const ComplexVector &x, ComplexVector &Bx) const;

  // Helper routine for computing the eigenpair residual.
  virtual PetscReal GetResidualNorm(PetscScalar l, const ComplexVector &x,
                                    ComplexVector &r) const = 0;

  // Helper routine for computing the backward error.
  virtual PetscReal GetBackwardScaling(PetscScalar l) const = 0;

public:
  SlepcEigenvalueSolver(int print);
  ~SlepcEigenvalueSolver() override;

  // Set operators for the generalized eigenvalue problem or for the quadratic polynomial
  // eigenvalue problem.
  void SetOperators(const ComplexOperator &K, const ComplexOperator &M,
                    ScaleType type) override;
  void SetOperators(const ComplexOperator &K, const ComplexOperator &C,
                    const ComplexOperator &M, ScaleType type) override;
  // For the linear generalized case, the linear solver should be configured to compute the
  // action of M⁻¹ (with no spectral transformation) or (K - σ M)⁻¹. For the quadratic
  // case, the linear solver should be configured to compute the action of M⁻¹ (with no
  // spectral transformation) or P(σ)⁻¹.
  void SetLinearSolver(const ComplexKspSolver &ksp) override;

  // Set the projection operator for enforcing the divergence-free constraint.
  void SetDivFreeProjector(const DivFreeSolver<ComplexVector> &divfree) override;

  // Set optional B matrix used for weighted inner products. This must be set explicitly
  // even for generalized problems, otherwise the identity will be used.
  void SetBMat(const Operator &B) override;

  // Get scaling factors used by the solver.
  PetscReal GetScalingGamma() const override { return gamma; }
  PetscReal GetScalingDelta() const override { return delta; }

  // Set shift-and-invert spectral transformation.
  void SetShiftInvert(std::complex<double> s, bool precond = false) override;

  // Set problem type.
  virtual void SetProblemType(ProblemType type) = 0;

  // Set eigenvalue solver.
  virtual void SetType(Type type) = 0;

  // Configure the basis vectors object associated with the eigenvalue solver.
  void SetOrthogonalization(bool mgs, bool cgs2);

  // Get the corresponding eigenpair error.
  PetscReal GetError(int i, ErrorType type) const override;

  // Re-normalize the given number of eigenvectors, for example if the matrix B for weighted
  // inner products has changed. This does not perform re-orthogonalization with respect to
  // the new matrix, only normalization.
  void RescaleEigenvectors(int num_eig) override;

  // Get the basis vectors object.
  virtual BV GetBV() const = 0;

  // Get the spectral transformation object.
  virtual ST GetST() const = 0;

  // Get the filtering region object.
  virtual RG GetRG() const = 0;

  // Get the associated MPI communicator.
  virtual MPI_Comm GetComm() const = 0;

  // Conversion function to PetscObject.
  virtual operator PetscObject() const = 0;
};

// Base class for SLEPc's EPS problem type.
class SlepcEPSSolverBase : public SlepcEigenvalueSolver
{
protected:
  // SLEPc eigensolver object. Polynomial problems are handled using linearization.
  EPS eps;

  // Shell matrices for the generalized eigenvalue problem.
  Mat A0, A1;

  void Customize() override;

public:
  // Calls SLEPc's EPSCreate. Expects SLEPc to be initialized/finalized externally.
  SlepcEPSSolverBase(MPI_Comm comm, int print, const std::string &prefix = std::string());

  // Call's SLEPc's EPSDestroy.
  ~SlepcEPSSolverBase() override;

  // Conversion function to SLEPc's EPS type.
  operator EPS() const { return eps; }

  void SetNumModes(int num_eig, int num_vec = 0) override;

  void SetTol(PetscReal tol) override;

  void SetMaxIter(int max_it) override;

  void SetWhichEigenpairs(WhichType type) override;

  void SetProblemType(ProblemType type) override;

  void SetType(Type type) override;

  void SetInitialSpace(const ComplexVector &v) override;

  int Solve() override;

  std::complex<double> GetEigenvalue(int i) const override;

  void GetEigenvector(int i, ComplexVector &x) const override;

  BV GetBV() const override;

  ST GetST() const override;

  RG GetRG() const override;

  MPI_Comm GetComm() const override
  {
    return eps ? PetscObjectComm(reinterpret_cast<PetscObject>(eps)) : MPI_COMM_NULL;
  }

  operator PetscObject() const override { return reinterpret_cast<PetscObject>(eps); };
};

// Generalized eigenvalue problem solver: K x = λ M x .
class SlepcEPSSolver : public SlepcEPSSolverBase
{
public:
  using SlepcEigenvalueSolver::delta;
  using SlepcEigenvalueSolver::gamma;
  using SlepcEigenvalueSolver::opB;
  using SlepcEigenvalueSolver::opInv;
  using SlepcEigenvalueSolver::opProj;
  using SlepcEigenvalueSolver::sigma;
  using SlepcEigenvalueSolver::sinvert;

  // References to matrices defining the generalized eigenvalue problem (not owned).
  const ComplexOperator *opK, *opM;

private:
  // Operator norms for scaling.
  mutable PetscReal normK, normM;

protected:
  PetscReal GetResidualNorm(PetscScalar l, const ComplexVector &x,
                            ComplexVector &r) const override;

  PetscReal GetBackwardScaling(PetscScalar l) const override;

public:
  SlepcEPSSolver(MPI_Comm comm, int print, const std::string &prefix = std::string());

  using SlepcEigenvalueSolver::SetOperators;
  void SetOperators(const ComplexOperator &K, const ComplexOperator &M,
                    ScaleType type) override;

  void SetBMat(const Operator &B) override;
};

// Quadratic eigenvalue problem solver: P(λ) x = (K + λ C + λ² M) x = 0 , solved via
// linearization: L₀ y = λ L₁ y .
class SlepcPEPLinearSolver : public SlepcEPSSolverBase
{
public:
  using SlepcEigenvalueSolver::delta;
  using SlepcEigenvalueSolver::gamma;
  using SlepcEigenvalueSolver::opB;
  using SlepcEigenvalueSolver::opInv;
  using SlepcEigenvalueSolver::opProj;
  using SlepcEigenvalueSolver::sigma;
  using SlepcEigenvalueSolver::sinvert;

  // References to matrices defining the quadratic polynomial eigenvalue problem
  // (not owned).
  const ComplexOperator *opK, *opC, *opM;

  // Workspace vectors for operator applications.
  mutable ComplexVector x2, y2;

private:
  // Operator norms for scaling.
  mutable PetscReal normK, normC, normM;

protected:
  PetscReal GetResidualNorm(PetscScalar l, const ComplexVector &x,
                            ComplexVector &r) const override;

  PetscReal GetBackwardScaling(PetscScalar l) const override;

public:
  SlepcPEPLinearSolver(MPI_Comm comm, int print, const std::string &prefix = std::string());

  using SlepcEigenvalueSolver::SetOperators;
  void SetOperators(const ComplexOperator &K, const ComplexOperator &C,
                    const ComplexOperator &M, ScaleType type) override;

  void SetBMat(const Operator &B) override;

  void SetInitialSpace(const ComplexVector &v) override;

  void GetEigenvector(int i, ComplexVector &x) const override;
};

// Base class for SLEPc's PEP problem type.
class SlepcPEPSolverBase : public SlepcEigenvalueSolver
{
protected:
  // SLEPc eigensolver object.
  PEP pep;

  // Shell matrices for the quadratic polynomial eigenvalue problem
  Mat A0, A1, A2;

  void Customize() override;

public:
  // Calls SLEPc's PEPCreate. Expects SLEPc to be initialized/finalized externally.
  SlepcPEPSolverBase(MPI_Comm comm, int print, const std::string &prefix = std::string());

  // Call's SLEPc's PEPDestroy.
  ~SlepcPEPSolverBase() override;

  // Conversion function to SLEPc's PEP type.
  operator PEP() const { return pep; }

  void SetNumModes(int num_eig, int num_vec = 0) override;

  void SetTol(PetscReal tol) override;

  void SetMaxIter(int max_it) override;

  void SetWhichEigenpairs(WhichType type) override;

  void SetProblemType(ProblemType type) override;

  void SetType(Type type) override;

  void SetInitialSpace(const ComplexVector &v) override;

  int Solve() override;

  std::complex<double> GetEigenvalue(int i) const override;

  void GetEigenvector(int i, ComplexVector &x) const override;

  BV GetBV() const override;

  ST GetST() const override;

  RG GetRG() const override;

  MPI_Comm GetComm() const override
  {
    return pep ? PetscObjectComm(reinterpret_cast<PetscObject>(pep)) : MPI_COMM_NULL;
  }

  operator PetscObject() const override { return reinterpret_cast<PetscObject>(pep); };
};

// Quadratic eigenvalue problem solver: P(λ) x = (K + λ C + λ² M) x = 0 .
class SlepcPEPSolver : public SlepcPEPSolverBase
{
public:
  using SlepcEigenvalueSolver::delta;
  using SlepcEigenvalueSolver::gamma;
  using SlepcEigenvalueSolver::opB;
  using SlepcEigenvalueSolver::opInv;
  using SlepcEigenvalueSolver::opProj;
  using SlepcEigenvalueSolver::sigma;
  using SlepcEigenvalueSolver::sinvert;

  // References to matrices defining the quadratic polynomial eigenvalue problem
  // (not owned).
  const ComplexOperator *opK, *opC, *opM;

private:
  // Operator norms for scaling.
  mutable PetscReal normK, normC, normM;

protected:
  PetscReal GetResidualNorm(PetscScalar l, const ComplexVector &x,
                            ComplexVector &r) const override;

  PetscReal GetBackwardScaling(PetscScalar l) const override;

public:
  SlepcPEPSolver(MPI_Comm comm, int print, const std::string &prefix = std::string());

  using SlepcEigenvalueSolver::SetOperators;
  void SetOperators(const ComplexOperator &K, const ComplexOperator &C,
                    const ComplexOperator &M, ScaleType type) override;
  void SetBMat(const Operator &B) override;
};

}  // namespace slepc

}  // namespace palace

#endif

#endif  // PALACE_LINALG_SLEPC_HPP
