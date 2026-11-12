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
#include <optional>
#include <string>
#include <vector>
#include <mpi.h>
#include "linalg/eps.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

// Forward declarations of SLEPc objects.
typedef struct _p_EPS *EPS;
typedef struct _p_PEP *PEP;
typedef struct _p_NEP *NEP;
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
    GYROSCOPIC,
    GENERAL,
  };

  enum class Type
  {
    KRYLOVSCHUR,
    POWER,
    SUBSPACE,
    TOAR,
    STOAR,
    QARNOLDI,
    JD,
    SLP,
    NLEIGS
  };

  // Workspace vector for operator applications.
  mutable ComplexVector x1, y1;

  // References to matrices defining the (possibly quadratic) eigenvalue problem
  // (not owned).
  const ComplexOperator *opK, *opC, *opM;

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
  ComplexKspSolver *opInv;

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

  // Set operators for the generalized eigenvalue problem, the quadratic polynomial
  // eigenvalue problem, or the nonlinear eigenvalue problem.
  void SetOperators(const ComplexOperator &K, const ComplexOperator &M,
                    ScaleType type) override;
  void SetOperators(const ComplexOperator &K, const ComplexOperator &C,
                    const ComplexOperator &M, ScaleType type) override;

  //  For the linear generalized case, the linear solver should be configured to compute the
  //  action of M⁻¹ (with no spectral transformation) or (K - σ M)⁻¹. For the quadratic
  //  case, the linear solver should be configured to compute the action of M⁻¹ (with no
  //  spectral transformation) or P(σ)⁻¹.
  void SetLinearSolver(ComplexKspSolver &ksp) override;

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
// linearization: L₀ y = λ L₁ y . With rational terms (AddPoleTerms), the problem is
//   T(λ) x = (P(λ) + Σₜ (Σⱼ cₜⱼ λ² / (λ - pₜⱼ)) Bₜ) x = 0 ,
// and the linearization is still exact: each pole adds the auxiliary block
// wₜⱼ = λ² x / (λ - pₜⱼ), which satisfies pₜⱼ wₜⱼ = λ (wₜⱼ - λ x) (Su and Bai, Solving
// rational eigenvalue problems via linearization, SIAM J. Matrix Anal. Appl. 2011).
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

  // Workspace vectors for operator applications.
  mutable ComplexVector x2, y2;

  // A rational term (Σⱼ cⱼ λ² / (λ - pⱼ)) B of the eigenvalue problem (B is not owned).
  // The auxiliary unknown of each pole exists only on the true dofs where B has a nonzero
  // row (mask is 1 there, 0 elsewhere). A full-size auxiliary block would give the
  // linearization the eigenvalue pⱼ once per true dof outside the support of B, with a
  // zero physical part: these are not eigenvalues of T.
  struct PoleTerm
  {
    const ComplexOperator *B;
    std::vector<std::complex<double>> poles, coeffs;
    Vector mask;
  };
  std::vector<PoleTerm> pole_terms;

  // Workspace vectors for the auxiliary blocks, in the order of pole_terms and their poles.
  mutable std::vector<ComplexVector> xa, ya;
  mutable ComplexVector za;

  // y += Σₜ Bₜ (Σⱼ sₜⱼ aₜⱼ) for the auxiliary blocks a: sₜⱼ = cₜⱼ for the operator L₀, or
  // sₜⱼ = cₜⱼ σ γ / (σ - pₜⱼ) for the right-hand side of the shift-and-invert solve.
  void AddPoleTermsMult(const std::vector<ComplexVector> &a, bool shifted,
                        ComplexVector &y) const;

  // Restrict an auxiliary block of the given rational term to the support of its B.
  void MaskPoleBlock(const PoleTerm &term, ComplexVector &a) const;

private:
  // Operator norms for scaling.
  mutable PetscReal normK, normC, normM;

  void SetOperatorsImpl(const ComplexOperator &K, const ComplexOperator *C,
                        const ComplexOperator &M, ScaleType type);

protected:
  PetscReal GetResidualNorm(PetscScalar l, const ComplexVector &x,
                            ComplexVector &r) const override;

  PetscReal GetBackwardScaling(PetscScalar l) const override;

public:
  SlepcPEPLinearSolver(MPI_Comm comm, int print, const std::string &prefix = std::string());

  // Without C, the problem is still solved through the linearization (for rational terms).
  void SetOperators(const ComplexOperator &K, const ComplexOperator &M,
                    ScaleType type) override;
  void SetOperators(const ComplexOperator &K, const ComplexOperator &C,
                    const ComplexOperator &M, ScaleType type) override;

  // Add the rational term (Σⱼ cⱼ λ² / (λ - pⱼ)) B to the problem, with distinct nonzero
  // poles pⱼ. Must be called on every process before SetOperators. The linear solver
  // (SetLinearSolver) must then solve with T(σ), which includes every rational term at the
  // shift σ, and σ must not be a pole.
  void AddPoleTerms(const ComplexOperator &B,
                    const std::vector<std::complex<double>> &poles,
                    const std::vector<std::complex<double>> &coeffs);

  void SetBMat(const Operator &B) override;

  void SetInitialSpace(const ComplexVector &v) override;

  // A second start for Solve, which then solves the same problem again from it. A Krylov
  // solve from one start finds only one direction of an exactly degenerate eigenspace (in
  // exact arithmetic), so one of a degenerate pair can be missing. With a second start the
  // eigenpairs are those converged in both solves, the first solve's first, each in SLEPc's
  // order: a pool of initial guesses for a nonlinear refinement, which may hold the same
  // eigenpair twice, not a set of distinct eigenpairs. It must be set on every process
  // (with the local part of the start, which may be empty) before SetOperators; the solver
  // then takes one SetOperators and one Solve.
  void SetSecondInitialSpace(const ComplexVector &v);

  int Solve() override;

  std::complex<double> GetEigenvalue(int i) const override;

  void GetEigenvector(int i, ComplexVector &x) const override;

private:
  ComplexVector v_second;
  bool second_start = false, solved = false;

  // Eigenpairs of both solves (the vectors before the normalization by xscale), used by
  // the getters once both solves are done.
  std::vector<std::complex<double>> pool_eig;
  std::vector<ComplexVector> pool_vec;
  bool pooled = false;

  // The physical block of the i-th eigenvector of the last SLEPc solve, unnormalized.
  void GetSolveEigenvector(int i, ComplexVector &x) const;
};

// Base class for SLEPc's PEP problem type.
class SlepcPEPSolverBase : public SlepcEigenvalueSolver
{
protected:
  // SLEPc eigensolver object.
  PEP pep;

  // Shell matrices for the quadratic polynomial eigenvalue problem.
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

// Base class for SLEPc's NEP problem type
class SlepcNEPSolverBase : public SlepcEigenvalueSolver
{
protected:
  // SLEPc eigensolver object.
  NEP nep;

  // Shell matrices for the nonlinear eigenvalue problem.
  Mat A, J;

  // Order of sorted eigenvalues.
  std::unique_ptr<int[]> perm;

  void Customize() override;

public:
  // Calls SLEPc's NEPCreate. Expects SLEPc to be initialized/finalized externally.
  SlepcNEPSolverBase(MPI_Comm comm, int print, const std::string &prefix = std::string());

  // Call's SLEPc's NEPDestroy.
  ~SlepcNEPSolverBase() override;

  // Conversion function to SLEPc's PEP type.
  operator NEP() const { return nep; }

  void SetNumModes(int num_eig, int num_vec = 0) override;

  void SetTol(PetscReal tol) override;

  void SetMaxIter(int max_it) override;

  void SetWhichEigenpairs(WhichType type) override;

  void SetShiftInvert(std::complex<double> s, bool precond = false) override;

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
    return nep ? PetscObjectComm(reinterpret_cast<PetscObject>(nep)) : MPI_COMM_NULL;
  }

  operator PetscObject() const override { return reinterpret_cast<PetscObject>(nep); };
};

// Nonlinear eigenvalue problem solver: T(λ) x = (K + λ C + λ² M + A2(λ)) x = 0.
class SlepcNEPSolver : public SlepcNEPSolverBase
{
public:
  using SlepcEigenvalueSolver::delta;
  using SlepcEigenvalueSolver::gamma;
  using SlepcEigenvalueSolver::opB;
  using SlepcEigenvalueSolver::opInv;
  using SlepcEigenvalueSolver::opProj;
  using SlepcEigenvalueSolver::sigma;
  using SlepcEigenvalueSolver::sinvert;

  // Operators for the nonlinear eigenvalue problem.
  // Function and Jacobian shell matrices must retain independent A2 owners: SLEPc may
  // apply either shell after invoking the other callback. The top-level opA/opJ sums hold
  // non-owning references to these subordinate operators.
  std::unique_ptr<ComplexOperator> opA2, opA;                   // Function T(λ)
  std::unique_ptr<ComplexOperator> opA2_jac, opA2p, opAJ, opJ;  // Jacobian T′(λ)
  std::unique_ptr<ComplexOperator> opA2_pc, opA_pc, opP_pc;     // Preconditioner

  // Known poles of A2(λ): the Jacobian shell applies T′ + L(λ) T, the derivative of the
  // pole-cleared d T (see nleps::PoleTerms). opA_jac is T at the Jacobian's λ.
  nleps::PoleTerms poles;
  std::unique_ptr<ComplexOperator> opA_jac;
  std::complex<double> logd_J = 0.0;

  // Function to compute the A2 operator.
  std::optional<std::function<std::unique_ptr<ComplexOperator>(std::complex<double>)>>
      funcA2;

  // Function to compute the preconditioner matrix.
  std::optional<std::function<std::unique_ptr<ComplexOperator>(
      std::complex<double>, std::complex<double>, std::complex<double>,
      std::complex<double>)>>
      funcP;

  // Eigenvalue estimate at current iteration.
  PetscScalar lambda;

  // Boolean flag to identify new λ estimate requiring a preconditioner update.
  bool new_lambda = true;

  // Boolean flag to avoid modifying an unused preconditioner.
  bool first_pc = true;

private:
  // Operator norms for scaling.
  mutable PetscReal normK, normC, normM;

  // Shared by the SetOperators overloads, with C = nullptr when there is no damping.
  void SetOperatorsImpl(const ComplexOperator &K, const ComplexOperator *C,
                        const ComplexOperator &M, ScaleType type);

protected:
  PetscReal GetResidualNorm(PetscScalar l, const ComplexVector &x,
                            ComplexVector &r) const override;

  PetscReal GetBackwardScaling(PetscScalar l) const override;

public:
  SlepcNEPSolver(MPI_Comm comm, int print, const std::string &prefix = std::string());

  using SlepcEigenvalueSolver::SetOperators;
  void SetOperators(const ComplexOperator &K, const ComplexOperator &M,
                    ScaleType type) override;
  void SetOperators(const ComplexOperator &K, const ComplexOperator &C,
                    const ComplexOperator &M, ScaleType type) override;
  void SetBMat(const Operator &B) override;

  // Set the frequency-dependent A2 matrix function.
  void SetExtraSystemMatrix(
      std::function<std::unique_ptr<ComplexOperator>(std::complex<double>)>) override;

  // Set the preconditioner update function.
  void SetPreconditionerUpdate(std::function<std::unique_ptr<ComplexOperator>(
                                   std::complex<double>, std::complex<double>,
                                   std::complex<double>, std::complex<double>)>) override;

  // Set the known poles of A2(λ).
  void SetExtraSystemPoles(const std::vector<std::complex<double>> &p) override
  {
    poles.poles = p;
  }
};

}  // namespace slepc

}  // namespace palace

#endif

#endif  // PALACE_LINALG_SLEPC_HPP
