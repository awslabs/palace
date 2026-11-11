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
#include <utility>
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

// Configure the region of the eigenvalues λ = iω (scaled) of a quadratic eigenvalue problem
// searched for above the shift-and-invert target λ = i s, s > 0: Im{λ} >= s and quality
// factors Q = |ω| / (2 Im{ω}) >= q_min > 0.5 for damped modes (Re{λ} < 0) and |Q| >= 10 for
// growing ones (Re{λ} > 0).
void ConfigureQualityFactorRegion(RG rg, PetscReal s, double q_min);

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

  // Parameters defining the spectral transformation: the shift σ, and the target τ of the
  // eigenvalue search (sorting criterion and search region), which is σ unless set apart.
  PetscScalar sigma, tau;
  bool sinvert, region, has_tau, st_shift_set;

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

  // Minimum quality factor of the eigenvalues searched for (see SetMinQualityFactor),
  // whether the search region for it is configured, and, with printing, the Ritz values of
  // the last iteration: their number, those outside of the region, and those outside with
  // an error estimate below the tolerance (see CountRegionExcluded).
  double q_min;
  bool q_region, q_region_monitor;
  int q_region_count[3];

  // Print the Ritz values counted outside of the quality-factor region.
  void PrintRegionExcluded() const;

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

  // Set the target of the eigenvalue search apart from the shift (after SetShiftInvert).
  void SetTarget(std::complex<double> t) override;

  // Factor converting the relative error estimate of a Ritz value θ of the transformed
  // problem, a relative error of λ - σ = 1 / θ, into one for the eigenvalue λ, with a
  // target apart from the shift: max(1, |λ - σ| / |λ|). An accepted eigenvalue then has a
  // relative error estimate below the tolerance, as with σ = τ (|λ - τ| < |λ| for all λ in
  // the search region).
  PetscReal GetTargetErrorScale(PetscScalar theta) const;

  // Set problem type.
  virtual void SetProblemType(ProblemType type) = 0;

  // Set eigenvalue solver.
  virtual void SetType(Type type) = 0;

  // Configure the basis vectors object associated with the eigenvalue solver.
  void SetOrthogonalization(bool mgs, bool cgs2);

  void SetMinQualityFactor(double q) override { q_min = q; }

  // Count the nest current Ritz values of the transformed problem outside of the search
  // region, and those among them with an error estimate below tol (called by the
  // eigensolver monitor at every iteration).
  void CountRegionExcluded(PetscReal tol, const PetscScalar *eigr, const PetscScalar *eigi,
                           const PetscReal *errest, int nest);

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

  // Constant for the relaxation of the linear solver tolerance (see
  // SetLinearSolverRelaxation), the linear solver tolerance outside of the relaxation, the
  // largest transformed eigenvalue estimate of the solve, and whether the monitor updating
  // the tolerance is installed.
  double relax_c, relax_lin_tol, relax_theta_max;
  bool relax_monitor;

  // Largest linear solver tolerance of the last solve, whether the relaxation is suspended
  // (strict tolerance, for a repeated solve), and a copy of the initial space vector.
  double relax_max_used = 0.0;
  bool relax_suspend = false;

  // Linear solver tolerance of the first Krylov cycle (see SetFirstCycleLinearSolverTol).
  double first_cycle_tol = 0.0;

  // For the convergence test (see CheckRitzValue): the Ritz values of the transformed
  // problem at the previous restart, at the last completed convergence check and at the
  // current one (a Krylov-Schur cycle can test convergence several times, see
  // EndRitzValueCheck), each with the number of consecutive restarts at which it was kept
  // unconverged, the Ritz values accepted although they did not settle, the iteration of
  // the current restart, whether an unconverged Ritz value was already tested at the
  // current check, and the numbers of convergences delayed and of eigenpairs accepted with
  // a changing Ritz value in the last solve.
  std::vector<std::pair<PetscScalar, int>> ritz_prev, ritz_last, ritz_curr;
  std::vector<PetscScalar> ritz_unsettled;
  PetscInt ritz_its = -1;
  bool ritz_blocked = false;
  int ritz_num_delayed = 0, ritz_num_unsettled = 0;

  // Whether the start vector of the next solve is a warm start (see SetWarmStart).
  bool warm_start = false;

  // Krylov-Schur steps between convergence checks (see SetConvergenceCheckInterval).
  int check_interval = 0;
  Vec v_init = nullptr;
  bool has_init = false;

  // Relative tolerance of the divergence-free projection of the shift-and-invert solve
  // results (see CheckDivFreeTol), the ratio |λ - σ| / |σ| of the wanted Ritz values which
  // stops the current solve to repeat it with a smaller tolerance (0 while it continues),
  // the largest such ratio of the current solve and that of its last restart, the number
  // of solves repeated for it, whether the current solve is the last one allowed (it is
  // not stopped), whether the monitor checking it is installed, and the eigenvalue
  // estimates for the next solve (see SetEigenvalueEstimates).
  double divfree_tol = 1.0e-3, divfree_stop = 0.0, divfree_ratio = 0.0,
         divfree_ratio_last = 0.0;
  int divfree_repeats = 0;
  bool divfree_final = false, divfree_monitor = false;
  std::vector<std::complex<double>> divfree_estimates;

  // Keep a copy of the initial space vector, for a repeated solve.
  void SaveInitialSpace(Vec v);

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

  bool SupportsTarget() const override { return true; }

  void SetProblemType(ProblemType type) override;

  void SetType(Type type) override;

  void SetInitialSpace(const ComplexVector &v) override;

  void ResetInitialSpace() override;

  void SetLinearSolverRelaxation(double c) override { relax_c = c; }

  void SetFirstCycleLinearSolverTol(double tol) override { first_cycle_tol = tol; }

  void SetWarmStart(bool warm) override { warm_start = warm; }

  void SetConvergenceCheckInterval(int k) override { check_interval = k; }

  void SetEigenvalueEstimates(const std::vector<std::complex<double>> &lambda) override
  {
    divfree_estimates = lambda;
  }

  // Update the error estimate of the Ritz value theta of the transformed problem at the
  // current restart (called by the convergence test, after SLEPc's relative estimate).
  void CheckRitzValue(PetscScalar theta, PetscReal *errest);

  // End a convergence check (called by the stopping test after the convergence tests of a
  // check): the Ritz values of the last check of a Krylov-Schur cycle are the ones present
  // at its restart.
  void EndRitzValueCheck();

  // Print the number of convergences delayed for Ritz values not settled in the last solve.
  void PrintRitzValueChecks() const;

  int GetNumIterations() const override;

  // Update the linear solver tolerance from the nest current Ritz values of the transformed
  // problem and their error estimates, of which the first nconv have converged and the
  // first nev are wanted (called by the eigensolver monitor after each restart).
  void UpdateLinearSolverTol(int nconv, int nev, const PetscScalar *eigr,
                             const PetscScalar *eigi, const PetscReal *errest, int nest);

  // Check the tolerance of the divergence-free projection of the solve results against the
  // nest current Ritz values of the transformed problem, of which the first nev are wanted
  // (called by the eigensolver monitor after each restart).
  void CheckDivFreeTol(int nev, const PetscScalar *eigr, const PetscScalar *eigi, int nest);

  // Whether the solve stops at the end of the current Krylov cycle, to be repeated with a
  // smaller tolerance of the divergence-free projection (see CheckDivFreeTol).
  bool StopForDivFreeTol() const { return divfree_stop > 0.0; }

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
  using SlepcEPSSolverBase::divfree_stop;
  using SlepcEPSSolverBase::divfree_tol;

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
  using SlepcEPSSolverBase::divfree_stop;
  using SlepcEPSSolverBase::divfree_tol;

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

  void SetInitialSpaceBlocks(const ComplexVector &v,
                             const ComplexVector &v_lambda) override;

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
  // Function and Jacobian shell matrices must retain independent A2 owners: SLEPc may
  // apply either shell after invoking the other callback. The top-level opA/opJ sums hold
  // non-owning references to these subordinate operators.
  std::unique_ptr<ComplexOperator> opA2, opA;                   // Function T(λ)
  std::unique_ptr<ComplexOperator> opA2_jac, opA2p, opAJ, opJ;  // Jacobian T′(λ)
  std::unique_ptr<ComplexOperator> opA2_pc, opA_pc, opP_pc;     // Preconditioner

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
};

}  // namespace slepc

}  // namespace palace

#endif

#endif  // PALACE_LINALG_SLEPC_HPP
