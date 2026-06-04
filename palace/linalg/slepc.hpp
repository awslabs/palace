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
#include <slepcfn.h>
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
  std::unique_ptr<ComplexOperator> opA2, opA2p, opJ, opA, opAJ, opA2_pc, opA_pc, opP_pc;

  // Function to compute the A2 operator at real ω = |Im λ|.
  std::optional<std::function<std::unique_ptr<ComplexOperator>(double)>> funcA2;

  // Optional A2(λ) evaluated at complex λ (analytic continuation). When set, the NEP
  // residual and Jacobian use this directly so the operator becomes holomorphic
  // (required by NLEIGS, beneficial for SLP). When unset, BuildA2(λ) falls back to
  // funcA2(|Im λ|) — the legacy real-ω stamping.
  std::optional<std::function<std::unique_ptr<ComplexOperator>(std::complex<double>)>>
      funcA2Complex;

  // Toggle between A2(λ) (when true) and A2(|Im λ|) (when false) for the production
  // residual / Jacobian path. Default false.
  bool use_complex_a2 = false;

  // Helper: prefer funcA2Complex(λ) when set AND use_complex_a2 is true, else fall
  // back to funcA2(|Im λ|).
  std::unique_ptr<ComplexOperator> BuildA2(std::complex<double> lambda) const;

  // BC frequency ω passed to funcP for the eigenvalue estimate λ, matching the A2
  // baked into the SLP system matrix: complex ω = -i·λ when the complex-λ A2 path is
  // active (so the preconditioner's wave-port / ABC / surf-σ terms match), else the
  // real frequency |Im λ|. SLP is a Newton iteration whose λ leaves the imaginary axis,
  // so this matters exactly as in the HYBRID path.
  std::complex<double> PreconditionerBCFreq(std::complex<double> lambda) const
  {
    return use_complex_a2 ? lambda / std::complex<double>(0.0, 1.0)
                          : std::complex<double>(std::abs(lambda.imag()), 0.0);
  }

  // Function to compute the preconditioner matrix. 4th arg (a3) is the complex BC
  // frequency ω = -i·λ; SLP is a Newton iteration so λ leaves the imaginary axis, and
  // passing the full complex ω makes the preconditioner's wave-port / ABC / surf-σ
  // terms match the exact complex system matrix.
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

  // Set the frequency-dependent A2 matrix function (real-ω stamping).
  void SetExtraSystemMatrix(
      std::function<std::unique_ptr<ComplexOperator>(double)>) override;

  // Set the holomorphic A2(λ) builder. Always installable; whether the solver uses
  // it (vs A2(|Im λ|)) is governed by SetUseComplexA2.
  void SetExtraSystemMatrixComplex(
      std::function<std::unique_ptr<ComplexOperator>(std::complex<double>)>) override;

  // Toggle the production path between A2(λ) and A2(|Im λ|).
  void SetUseComplexA2(bool b) override { use_complex_a2 = b; }

  // Set the preconditioner update function.
  void SetPreconditionerUpdate(std::function<std::unique_ptr<ComplexOperator>(
                                   std::complex<double>, std::complex<double>,
                                   std::complex<double>, std::complex<double>)>) override;
};

// SLEPc NEP NLEIGS solver in *split-form* mode. The nonlinear operator is exposed as
//   T(λ) = Σ_j  A_j · f_j(λ)
// with each A_j a constant PETSc Mat and each f_j a scalar SLEPc FN object. NLEIGS
// builds a rational interpolant of T(λ) on a target region Σ ⊂ ℂ via Leja–Bagby
// nodes; the rational degree adapts to user-supplied tolerance.
//
// The split-form path is fundamentally different from the SLP / matrix-callback path
// (which goes through NEPSetFunction with a shell A) — see comments below for the
// reason a separate class exists.
class SlepcNEPNLEIGSSolver : public SlepcNEPSolverBase
{
public:
  using SlepcEigenvalueSolver::delta;
  using SlepcEigenvalueSolver::gamma;
  using SlepcEigenvalueSolver::opB;
  using SlepcEigenvalueSolver::opInv;
  using SlepcEigenvalueSolver::opProj;
  using SlepcEigenvalueSolver::sigma;
  using SlepcEigenvalueSolver::sinvert;

  SlepcNEPNLEIGSSolver(MPI_Comm comm, int print, const std::string &prefix = std::string());
  ~SlepcNEPNLEIGSSolver() override;

  using SlepcEigenvalueSolver::SetOperators;
  void SetOperators(const ComplexOperator &K, const ComplexOperator &C,
                    const ComplexOperator &M, ScaleType type) override;
  void SetBMat(const Operator &B) override;

  // Append a (constant matrix, scalar FN) pair to the split-form operator. Must be
  // called *after* SetOperators (which populates K, C, M) and *before* Solve. The
  // ComplexOperator pointed to by `op` is NOT owned — the caller must keep it alive
  // for the lifetime of the solver. The FN passed in is moved into the solver and
  // destroyed by it.
  void AddSplitFormTerm(const ComplexOperator *op, FN fn);

  // Tuning knobs. SetFullBasis: NEPNLEIGSSetFullBasis (default false, i.e., SLEPc
  // default; the matrix-callback path was sensitive to TRUE in PR #467, but the
  // split-form path uses a different basis-conditioning code path so the same
  // sensitivity is not expected). SetInterpolation: tol/deg for the rational
  // interpolant. SetSingularities: the singularity set Ξ (poles of the rational
  // approximation). SetRegion: target region Σ as a flat 4-tuple [Re_min, Re_max,
  // Im_min, Im_max] (RGINTERVAL); empty → auto-derive from sigma at solve time.
  void SetNLEIGSFullBasis(bool full_basis) { full_basis_ = full_basis; }
  void SetNLEIGSInterpolation(double tol, int deg)
  {
    interp_tol_ = tol;
    interp_deg_ = deg;
  }
  void SetNLEIGSSingularities(const std::vector<std::complex<double>> &xi)
  {
    singularities_ = xi;
  }
  void SetNLEIGSRegion(const std::vector<double> &region) { region_ = region; }
  // Communicate the eigensolver's TargetUpper so the auto-derived RG bracket
  // matches [target, target_upper] on the imaginary axis. If unset (≤ 0), defaults
  // to 3·target at solve time.
  void SetTargetUpper(double target_upper) { target_upper_ = target_upper; }

  // The preconditioner builder funcP(a0, a1, a2, ω) returns Palace's preconditioner
  // for the operator a0·K + a1·C + a2·M + i·ω·(extra). NLEIGS calls KSPSolve at
  // multiple rational interpolation shifts σ_i; for each shift we need a PC that
  // approximates (T(σ_i))⁻¹. This setter installs the same builder as the SLP path
  // uses; the PCSHELL constructed inside Customize() invokes it lazily per shift.
  void SetPreconditionerBuilder(
      std::function<std::unique_ptr<ComplexOperator>(std::complex<double>,
                                                     std::complex<double>,
                                                     std::complex<double>,
                                                     std::complex<double>)>
          P)
  {
    funcP_ = std::move(P);
  }

  // SLP-only setters that are not meaningful for NLEIGS — provided so generic call
  // sites can configure both solvers uniformly. They are no-ops here.
  void
  SetExtraSystemMatrix(std::function<std::unique_ptr<ComplexOperator>(double)>) override
  {
  }
  void SetExtraSystemMatrixComplex(
      std::function<std::unique_ptr<ComplexOperator>(std::complex<double>)>) override
  {
  }
  void SetExtraSystemMatrixDerivative(
      std::function<std::unique_ptr<ComplexOperator>(double)>) override
  {
  }
  void SetExtraSystemMatrixDerivativeComplex(
      std::function<std::unique_ptr<ComplexOperator>(std::complex<double>)>) override
  {
  }
  void SetUseComplexA2(bool) override {}
  void SetPreconditionerUpdate(std::function<std::unique_ptr<ComplexOperator>(
                                   std::complex<double>, std::complex<double>,
                                   std::complex<double>, std::complex<double>)>) override
  {
  }

  // Override Solve to skip the base's `A && J && opInv` check (NLEIGS has no shell
  // A/J — split-form Mats are stored in term_mats_ instead). Also override
  // SetInitialSpace for the same reason.
  int Solve() override;
  void SetInitialSpace(const ComplexVector &v) override;

protected:
  PetscReal GetResidualNorm(PetscScalar l, const ComplexVector &x,
                            ComplexVector &r) const override;
  PetscReal GetBackwardScaling(PetscScalar l) const override;

public:
  // Per-term shell-matrix context. Must be public so the static shell-MULT callback
  // can access fields without a friend declaration.
  struct TermCtx
  {
    const ComplexOperator *op;
    mutable ComplexVector x, y;
  };

private:
  // Apply NLEIGS-specific tuning (FullBasis, Interpolation, Singularities, RG) and
  // call NEPSetSplitOperator. Invoked from Solve() after all terms are registered.
  void Customize();

  // Heap-allocated so MatShellGetContext() returns a stable pointer for NEP's
  // internal lifetime. TermCtx is declared public above so the static shell-MULT
  // callback can access it.
  std::vector<std::unique_ptr<TermCtx>> term_ctxs_;
  std::vector<Mat> term_mats_;
  std::vector<FN> term_fns_;

  // Tuning state.
  bool full_basis_ = false;
  double interp_tol_ = 1.0e-7;
  int interp_deg_ = 80;
  double target_upper_ = -1.0;
  std::vector<std::complex<double>> singularities_;
  std::vector<double> region_;

  // Preconditioner builder (installed by the eigensolver driver) — same signature as
  // SLP's funcP. Used by the per-shift PCSHELL that wraps each NLEIGS inner KSP.
  std::optional<std::function<std::unique_ptr<ComplexOperator>(
      std::complex<double>, std::complex<double>, std::complex<double>,
      std::complex<double>)>>
      funcP_;

public:
  // Per-inner-KSP PCSHELL context. NLEIGS holds nshifts inner KSPs; we attach one
  // KSPCtx to each via PCShellSetContext. The PCSHELL apply lazily builds T(σ_i)
  // and a preconditioner for it on first call, then forwards to the solver-level
  // opInv (Palace's ComplexKspSolver) reconfigured for that shift.
  struct KSPCtx
  {
    SlepcNEPNLEIGSSolver *solver;
    std::complex<double> sigma_shift;
    std::unique_ptr<ComplexOperator> opT;
    std::unique_ptr<ComplexOperator> opP;
    bool built = false;
  };

  // Helpers used by the file-scope PCSHELL apply.
  bool HasFuncP() const { return funcP_.has_value(); }
  std::unique_ptr<ComplexOperator> BuildTAtShift(std::complex<double> sigma) const;
  std::unique_ptr<ComplexOperator> BuildPAtShift(std::complex<double> sigma) const;
  ComplexKspSolver *GetOpInv() const { return opInv; }
  const DivFreeSolver<ComplexVector> *GetOpProj() const { return opProj; }

private:
  std::vector<std::unique_ptr<KSPCtx>> ksp_ctxs_;
  // Latched problem size from SetOperators, needed to construct shell matrices for
  // additional split-form terms.
  PetscInt problem_size_ = 0;

  // Operator norms for residual / backward-error scaling.
  mutable PetscReal normK = 0.0, normC = 0.0, normM = 0.0;
  // Per-A2-term spectral norms (lazy; -1 means not yet computed). Aligned with
  // term_ctxs_; slots 0..2 unused (those are K, C, M tracked in normK/C/M).
  mutable std::vector<PetscReal> term_norms_;
  const ComplexOperator *opK_ref = nullptr;
  const ComplexOperator *opC_ref = nullptr;
  const ComplexOperator *opM_ref = nullptr;
};

}  // namespace slepc

}  // namespace palace

#endif

#endif  // PALACE_LINALG_SLEPC_HPP
