// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_SLEPC_HPP
#define PALACE_SLEPC_HPP

#if defined(PALACE_WITH_SLEPC)

#include "linalg/petsc.hpp"

#if !defined(PETSC_USE_COMPLEX)
#error "SLEPc interface requires PETSc built with complex scalars!"
#endif

#include <string>
#include "linalg/eigen.hpp"

// Forward declarations of SLEPc objects.
typedef struct _p_EPS *EPS;
typedef struct _p_PEP *PEP;
typedef struct _p_BV *BV;
typedef struct _p_ST *ST;
typedef struct _p_RG *RG;

namespace palace
{

class DivFreeSolver;
class KspSolver;

namespace petsc
{

class PetscParMatrix;
class PetscParVector;

}  // namespace petsc

namespace slepc
{

// Wrappers for SlepcInitialize/SlepcInitializeNoArguments/SlepcFinalize.
void Initialize(int &argc, char **&argv, const char rc_file[], const char help[]);
void Initialize();
void Finalize();

// Compute and return the maximum singular value of the given operator, σₙ² = λₙ(Aᴴ A) .
PetscReal GetMaxSingularValue(const petsc::PetscParMatrix &A, PetscReal tol = PETSC_DEFAULT,
                              PetscInt maxits = PETSC_DEFAULT);

//
// A wrapper for the SLEPc library for generalized linear eigenvalue problems or quadratic
// polynomial eigenvalue problems. Shift-and-invert spectral transformations can be used to
// compute interior eigenvalues.
//
class SlepcEigenSolver : public EigenSolverBase
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
    ARPACK,
    POWER,
    SUBSPACE,
    TOAR,
    STOAR,
    QARNOLDI,
    JD
  };

protected:
  // Boolean to handle SetFromOptions calls.
  mutable bool clcustom;

  // Control print level for debugging.
  int print;

  // Variables for scaling, from Higham et al., IJNME 2008.
  PetscReal gamma, delta;

  // Parameters defining the spectral transformation.
  PetscScalar sigma;
  bool sinvert, region;

  // Storage for computed residual norms.
  mutable PetscReal *res;

  // Workspace vectors for initial space and residual calculations.
  mutable petsc::PetscParVector *v0, *r0;

  // Reference to linear solver used for operator action for M⁻¹ (with no spectral
  // transformation) or (K - σ M)⁻¹ (generalized EVP with shift-and- invert) or P(σ)⁻¹
  // (polynomial with shift-and-invert) (not owned).
  const KspSolver *opInv;

  // Reference to solver for projecting an intermediate vector onto a divergence-free space
  // (not owned).
  const DivFreeSolver *opProj;

  // Customize object with command line options set.
  virtual void Customize();

  // Configure KSP object associated with the spectral transformation.
  void SetPCShell(void *ctx, PetscErrorCode (*__pc_apply)(PC, Vec, Vec));

  // Specify rectangular region of the complex plane, bounded by[rminr, rmaxr] x
  // [rmini, rmaxi] in which to constrain eigenvalue search.
  void SetRegion(PetscReal rminr, PetscReal rmaxr, PetscReal rmini, PetscReal rmaxi,
                 bool complement = false);

  // Perform the back-transformation from the spectrally transformed eigenvalue back to the
  // original problem.
  void GetBackTransform(PetscScalar eig, PetscReal &eigr, PetscReal &eigi) const;

  // Helper routine for computing the eigenpair residual.
  virtual void GetResidual(PetscScalar eig, const petsc::PetscParVector &v,
                           petsc::PetscParVector &r) const = 0;

  // Helper routine for computing the backward error.
  virtual PetscReal GetBackwardScaling(PetscScalar eig) const = 0;

public:
  SlepcEigenSolver(int print_lvl);
  ~SlepcEigenSolver() override;

  SlepcEigenSolver(const SlepcEigenSolver &) = delete;
  SlepcEigenSolver(SlepcEigenSolver &&) = delete;
  SlepcEigenSolver &operator=(const SlepcEigenSolver &) = delete;
  SlepcEigenSolver &operator=(SlepcEigenSolver &&) = delete;

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

  // Set the projection operator or operators for the divergence-free constraint.
  void SetProjector(const DivFreeSolver &divfree) override;

  // Set optional B matrix used for weighted inner products. This must be set explicitly
  // even for generalized problems, otherwise the identity will be used.
  void SetBMat(const petsc::PetscParMatrix &B) override;

  // Get spectral transformation target used by the solver.
  bool IsShiftInvert() const { return sinvert; }
  PetscScalar GetTarget() const { return sigma; }

  // Get scaling factors used by the solver.
  double GetScalingGamma() const override { return (double)gamma; }
  double GetScalingDelta() const override { return (double)delta; }

  // Set shift-and-invert spectral transformation.
  void SetShiftInvert(double tr, double ti, bool precond = false) override;

  // Configure the basis vectors object associated with the eigenvalue solver.
  void SetOrthogonalization(bool mgs, bool cgs2);

  // Set the number of required eigenmodes.
  void SetNumModes(int numeig, int numvec = 0) override = 0;

  // Set solver tolerance.
  void SetTol(double tol) override = 0;

  // Set maximum number of iterations.
  void SetMaxIter(int maxits) override = 0;

  // Set target spectrum for the eigensolver. When a spectral transformation is used, this
  // applies to the spectrum of the shifted operator.
  void SetWhichEigenpairs(WhichType type) override = 0;

  // Set problem type.
  virtual void SetProblemType(ProblemType type) = 0;

  // Set eigenvalue solver.
  virtual void SetType(Type type) = 0;

  // Set an initial vector for the solution subspace.
  void SetInitialSpace(const petsc::PetscParVector &v) override = 0;

  // Solve the eigenvalue problem. Returns the number of converged eigenvalues.
  int Solve() override = 0;

  // Get the corresponding eigenvalue.
  void GetEigenvalue(int i, double &eigr, double &eigi) const override = 0;

  // Get the corresponding eigenvector.
  void GetEigenvector(int i, petsc::PetscParVector &v) const override = 0;

  // Get the corresponding eigenpair error.
  void GetError(int i, ErrorType type, double &err) const override;

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

  // Access solver object for callback functions.
  const KspSolver *GetKspSolver() const { return opInv; }

  // Access solver object for callback functions.
  const DivFreeSolver *GetDivFreeSolver() const { return opProj; }
};

// Base class for SLEPc's EPS problem type.
class SlepcEPSSolverBase : public SlepcEigenSolver
{
protected:
  // SLEPc eigensolver object. Polynomial problems are handled using linearization.
  EPS eps;

  // Shell matrices for the generalized eigenvalue problem.
  petsc::PetscParMatrix *A, *B;

  // Customize object with command line options set.
  void Customize() override;

public:
  // Calls SLEPc's EPSCreate. Expects SLEPc to be initialized/finalized externally.
  SlepcEPSSolverBase(MPI_Comm comm, int print_lvl,
                     const std::string &prefix = std::string());

  // Call's SLEPc's EPSDestroy.
  ~SlepcEPSSolverBase() override;
  SlepcEPSSolverBase(const SlepcEPSSolverBase &) = delete;
  SlepcEPSSolverBase(SlepcEPSSolverBase &&) = delete;
  SlepcEPSSolverBase &operator=(const SlepcEPSSolverBase &) = delete;
  SlepcEPSSolverBase &operator=(SlepcEPSSolverBase &&) = delete;

  // Set the number of required eigenmodes.
  void SetNumModes(int numeig, int numvec = 0) override;

  // Set solver tolerance.
  void SetTol(double tol) override;

  // Set maximum number of iterations.
  void SetMaxIter(int maxits) override;

  // Set target spectrum for the eigensolver. When a spectral transformation is used, this
  // applies to the spectrum of the shifted operator.
  void SetWhichEigenpairs(WhichType type) override;

  // Set problem type.
  void SetProblemType(ProblemType type) override;

  // Set eigenvalue solver.
  void SetType(Type type) override;

  // Set an initial vector for the solution subspace.
  void SetInitialSpace(const petsc::PetscParVector &v) override;

  // Solve the eigenvalue problem. Returns the number of converged eigenvalues.
  int Solve() override;

  // Get the corresponding eigenvalue.
  void GetEigenvalue(int i, double &eigr, double &eigi) const override;

  // Get the corresponding eigenvector.
  void GetEigenvector(int i, petsc::PetscParVector &v) const override;

  // Get the basis vectors object.
  BV GetBV() const override;

  // Get the spectral transformation object.
  ST GetST() const override;

  // Get the filtering region object.
  RG GetRG() const override;

  // Get the associated MPI communicator.
  MPI_Comm GetComm() const override;

  // Conversion function to SLEPc's EPS type.
  operator EPS() const { return eps; }

  // Conversion function to PetscObject.
  operator PetscObject() const override { return reinterpret_cast<PetscObject>(eps); };
};

// Generalized eigenvalue problem solver: K x = λ M x .
class SlepcEPSSolver : public SlepcEPSSolverBase
{
private:
  // References to matrices defining the generalized eigenvalue problem (not owned).
  const petsc::PetscParMatrix *opK, *opM;

  // Operator norms for scaling.
  mutable PetscReal normK, normM;

protected:
  // Helper routine for computing the eigenpair residual.
  void GetResidual(PetscScalar eig, const petsc::PetscParVector &v,
                   petsc::PetscParVector &r) const override;

  // Helper routine for computing the backward error.
  PetscReal GetBackwardScaling(PetscScalar eig) const override;

public:
  SlepcEPSSolver(MPI_Comm comm, int print_lvl, const std::string &prefix = std::string());

  // Set operators for the generalized eigenvalue problem.
  void SetOperators(const petsc::PetscParMatrix &K, const petsc::PetscParMatrix &M,
                    ScaleType type) override;

  // Access methods for operator application.
  const petsc::PetscParMatrix *GetOpK() { return opK; }
  const petsc::PetscParMatrix *GetOpM() { return opM; }
};

// Quadratic eigenvalue problem solver: P(λ) x = (K + λ C + λ² M) x = 0 , solved via
// linearization: L₀ y = λ L₁ y .
class SlepcPEPLinearSolver : public SlepcEPSSolverBase
{
private:
  // References to matrices defining the quadratic polynomial eigenvalue problem
  // (not owned).
  const petsc::PetscParMatrix *opK, *opC, *opM;

  // Operator norms for scaling.
  mutable PetscReal normK, normC, normM;

  // Shell matrix used for weighted inner products. May be nullptr, in which case identity
  // is used. Also a reference to the original passed in matrix.
  petsc::PetscParMatrix *B0;
  const petsc::PetscParMatrix *opB;

  // Workspace vectors for operator applications.
  mutable petsc::PetscParVector *x1, *x2, *y1, *y2, *z;

protected:
  // Helper routine for computing the eigenpair residual.
  void GetResidual(PetscScalar eig, const petsc::PetscParVector &v,
                   petsc::PetscParVector &r) const override;

  // Helper routine for computing the backward error.
  PetscReal GetBackwardScaling(PetscScalar eig) const override;

public:
  SlepcPEPLinearSolver(MPI_Comm comm, int print_lvl,
                       const std::string &prefix = std::string());
  ~SlepcPEPLinearSolver() override;
  SlepcPEPLinearSolver(const SlepcPEPLinearSolver &) = delete;
  SlepcPEPLinearSolver(SlepcPEPLinearSolver &&) = delete;
  SlepcPEPLinearSolver &operator=(const SlepcPEPLinearSolver &) = delete;
  SlepcPEPLinearSolver &operator=(SlepcPEPLinearSolver &&) = delete;

  // Set operators for the quadratic polynomial eigenvalue problem.
  void SetOperators(const petsc::PetscParMatrix &K, const petsc::PetscParMatrix &C,
                    const petsc::PetscParMatrix &M, ScaleType type) override;

  // Configure the basis vectors object associated with the eigenvalue solver.
  void SetBMat(const petsc::PetscParMatrix &B) override;

  // Set an initial vector for the solution subspace.
  void SetInitialSpace(const petsc::PetscParVector &v) override;

  // Get the corresponding eigenvector.
  void GetEigenvector(int i, petsc::PetscParVector &v) const override;

  // Helper methods for splitting a block vector from the linearized problem into its into
  // two parts.
  PetscScalar *GetBlocks(petsc::PetscParVector &v, petsc::PetscParVector &v1,
                         petsc::PetscParVector &v2) const;
  const PetscScalar *GetBlocksRead(const petsc::PetscParVector &v,
                                   petsc::PetscParVector &v1,
                                   petsc::PetscParVector &v2) const;
  void RestoreBlocks(PetscScalar *pv, petsc::PetscParVector &v, petsc::PetscParVector &v1,
                     petsc::PetscParVector &v2) const;
  void RestoreBlocksRead(const PetscScalar *pv, const petsc::PetscParVector &v,
                         petsc::PetscParVector &v1, petsc::PetscParVector &v2) const;

  // Access methods for operator application.
  const petsc::PetscParMatrix *GetOpK() { return opK; }
  const petsc::PetscParMatrix *GetOpC() { return opC; }
  const petsc::PetscParMatrix *GetOpM() { return opM; }
  const petsc::PetscParMatrix *GetOpB() { return opB; }
  petsc::PetscParVector *GetX1() { return x1; }
  petsc::PetscParVector *GetX2() { return x2; }
  petsc::PetscParVector *GetY1() { return y1; }
  petsc::PetscParVector *GetY2() { return y2; }
};

// Base class for SLEPc's PEP problem type.
class SlepcPEPSolverBase : public SlepcEigenSolver
{
protected:
  // SLEPc eigensolver object.
  PEP pep;

  // Shell matrices for the quadratic polynomial eigenvalue problem
  petsc::PetscParMatrix *A0, *A1, *A2;

  // Customize object with command line options set.
  void Customize() override;

public:
  // Calls SLEPc's PEPCreate. Expects SLEPc to be initialized/finalized externally.
  SlepcPEPSolverBase(MPI_Comm comm, int print_lvl,
                     const std::string &prefix = std::string());

  // Call's SLEPc's PEPDestroy.
  ~SlepcPEPSolverBase() override;
  SlepcPEPSolverBase(const SlepcPEPSolverBase &) = delete;
  SlepcPEPSolverBase(SlepcPEPSolverBase &&) = delete;
  SlepcPEPSolverBase &operator=(const SlepcPEPSolverBase &) = delete;
  SlepcPEPSolverBase &operator=(SlepcPEPSolverBase &&) = delete;

  // Set the number of required eigenmodes.
  void SetNumModes(int numeig, int numvec = 0) override;

  // Set solver tolerance.
  void SetTol(double tol) override;

  // Set maximum number of iterations.
  void SetMaxIter(int maxits) override;

  // Set target spectrum for the eigensolver. When a spectral transformation is used, this
  // applies to the spectrum of the shifted operator.
  void SetWhichEigenpairs(WhichType type) override;

  // Set problem type.
  void SetProblemType(ProblemType type) override;

  // Set eigenvalue solver.
  void SetType(Type type) override;

  // Set an initial vector for the solution subspace.
  void SetInitialSpace(const petsc::PetscParVector &v) override;

  // Solve the eigenvalue problem. Returns the number of converged eigenvalues.
  int Solve() override;

  // Get the corresponding eigenvalue.
  void GetEigenvalue(int i, double &eigr, double &eigi) const override;

  // Get the corresponding eigenvector.
  void GetEigenvector(int i, petsc::PetscParVector &v) const override;

  // Get the basis vectors object.
  BV GetBV() const override;

  // Get the spectral transformation object.
  ST GetST() const override;

  // Get the filtering region object.
  RG GetRG() const override;

  // Get the associated MPI communicator.
  MPI_Comm GetComm() const override;

  // Conversion function to SLEPc's PEP type.
  operator PEP() const { return pep; }

  // Conversion function to PetscObject.
  operator PetscObject() const override { return reinterpret_cast<PetscObject>(pep); };
};

// Quadratic eigenvalue problem solver: P(λ) x = (K + λ C + λ² M) x = 0 .
class SlepcPEPSolver : public SlepcPEPSolverBase
{
private:
  // References to matrices defining the quadratic polynomial eigenvalue problem
  // (not owned).
  const petsc::PetscParMatrix *opK, *opC, *opM;

  // Operator norms for scaling.
  mutable PetscReal normK, normC, normM;

protected:
  // Helper routine for computing the eigenpair residual.
  void GetResidual(PetscScalar eig, const petsc::PetscParVector &v,
                   petsc::PetscParVector &r) const override;

  // Helper routine for computing the backward error.
  PetscReal GetBackwardScaling(PetscScalar eig) const override;

public:
  SlepcPEPSolver(MPI_Comm comm, int print_lvl, const std::string &prefix = std::string());

  // Set operators for the quadratic polynomial eigenvalue problem.
  void SetOperators(const petsc::PetscParMatrix &K, const petsc::PetscParMatrix &C,
                    const petsc::PetscParMatrix &M, ScaleType type) override;

  // Access methods for operator application.
  const petsc::PetscParMatrix *GetOpK() { return opK; }
  const petsc::PetscParMatrix *GetOpC() { return opC; }
  const petsc::PetscParMatrix *GetOpM() { return opM; }
};

}  // namespace slepc

}  // namespace palace

#endif

#endif  // PALACE_SLEPC_HPP
