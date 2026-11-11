// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_KSP_HPP
#define PALACE_LINALG_KSP_HPP

#include <memory>
#include <string>
#include <string_view>
#include <type_traits>
#include "linalg/iterative.hpp"
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "linalg/solverstats.hpp"
#include "utils/labels.hpp"

namespace palace
{

class DeflatedPreconditioner;
class DeflationLevels;
class FiniteElementSpaceHierarchy;
class IoData;

namespace config
{

struct LinearSolverData;

}  // namespace config

// Derive the preconditioner matrix symmetry (used by sparse direct solvers) from the
// problem type, shifted-preconditioner flag, and periodicity. This is a derived runtime
// property, not a configuration field, so it is computed here rather than stored on
// LinearSolverData.
MatrixSymmetry GetPreconditionerMatrixSymmetry(const IoData &iodata);

//
// Linear solver class composing an iterative solver and preconditioner object.
//
template <OperatorType OperType>
class BaseKspSolver
{
  using VecType =
      std::conditional_t<std::is_same_v<OperType, ComplexOperator>, ComplexVector, Vector>;

protected:
  // The actual solver and preconditioner objects.
  std::unique_ptr<IterativeSolver<OperType>> ksp;
  std::unique_ptr<Solver<OperType>> pc;

  // Counters for number of calls to Mult method for linear solves, and cumulative number
  // of iterations.
  mutable int ksp_mult, ksp_mult_it;

  // Enable timer contribution for Timer::KSP_PRECONDITIONER.
  bool use_timer;

  // Solver role for the performance report, and its statistics.
  std::string role;
  int perf_solve_tally, perf_setup_tally, perf_basis;
  perf::KspStats perf_stats;

  // Preconditioner wrapper recording an inclusive tally for nested solvers.
  std::unique_ptr<Solver<OperType>> pc_tally;

  // Relative tolerance of a repeated solve after a false convergence (see
  // SetFallbackRelTol), 0 for none.
  double fallback_rel_tol = 0.0;

  // Whether a false convergence is returned to the caller instead of being an error (see
  // SetAllowFalseConvergence), and whether the last solve converged falsely (the reference
  // solves of the automatic deflation do not change it).
  bool allow_false_convergence = false;
  mutable bool false_converged = false;

  // Deflation of the solves after the first one with the same operator (see SetDeflation):
  // requested number of vectors (< 0 for the automatic choice), deflated preconditioner and
  // the nested sets of vectors of the automatic choice, operator, whether the next solve
  // (without deflation) provides the basis, whether the deflation stopped, whether the
  // model predicts a large gain, and whether the next solve confirms a verification, the
  // deflated solves since the last check and the count of the next verification, the
  // relative cost per iteration of the deflation per vector and per application, and the
  // iterations of the last deflated solves.
  int defl_size = 0;
  std::unique_ptr<DeflatedPreconditioner> defl;
  std::unique_ptr<DeflationLevels> defl_levels;
  const OperType *defl_op = nullptr;
  mutable bool defl_build = false, defl_off = false, defl_confident = false,
               defl_confirm = false;
  mutable int defl_solves = 0, defl_next_check = 0;
  mutable double defl_cost_vec = 0.0, defl_cost_red = 0.0;
  mutable std::vector<int> defl_its;

  // The preconditioner of the iterative solver without deflation.
  const Solver<OperType> &GetBasePreconditioner() const;

  // Analyze the Arnoldi relation of the solve just completed without deflation (t_it
  // seconds per iteration) and build the deflation basis, choosing its size for the
  // automatic mode.
  void BuildDeflation(double t_it) const;

  // Relative cost per iteration of the solves deflated with k vectors.
  double DeflationCost(int k) const;

  // Automatic deflation: from the solve of b just completed with its iterations (without
  // deflation for a verification, with deflation for its confirmation and with the set
  // just activated for a comparison), decide whether the current set of vectors pays,
  // against no deflation or against the previous set.
  enum class DeflationCheckKind
  {
    VERIFICATION,
    CONFIRMATION,
    COMPARISON
  };
  void CheckDeflation(const VecType &b, int its, DeflationCheckKind kind) const;

  // Iterations of the solve of b with the first k deflation vectors (k = 0: without
  // deflation), stopped after limit iterations (limit + 1 if not converged by then).
  int ReferenceSolve(const VecType &b, int k, int limit) const;

public:
  BaseKspSolver(const config::LinearSolverData &linear, MatrixSymmetry pc_mat_sym,
                int verbose, FiniteElementSpaceHierarchy &fespaces,
                FiniteElementSpaceHierarchy *aux_fespaces = nullptr);
  BaseKspSolver(const IoData &iodata, FiniteElementSpaceHierarchy &fespaces,
                FiniteElementSpaceHierarchy *aux_fespaces = nullptr);
  // An empty role means the role active at construction (perf::CurrentRole()).
  BaseKspSolver(std::unique_ptr<IterativeSolver<OperType>> &&ksp,
                std::unique_ptr<Solver<OperType>> &&pc, std::string_view role = {});
  ~BaseKspSolver();

  int NumTotalMult() const { return ksp_mult; }
  int NumTotalMultIterations() const { return ksp_mult_it; }

  const std::string &GetRole() const { return role; }

  // The preconditioner, without deflation (nullptr for none).
  const Solver<OperType> *GetPreconditioner() const { return pc.get(); }
  Solver<OperType> *GetPreconditioner() { return pc.get(); }

  // Whether the most recent Mult converged to the requested tolerance.
  bool GetConverged() const { return ksp->GetConverged(); }

  // Forward tolerance access to the underlying iterative solver.
  double GetRelTol() const { return ksp->GetRelTol(); }
  double GetAbsTol() const { return ksp->GetAbsTol(); }
  void SetRelTol(double tol) { ksp->SetRelTol(tol); }
  void SetAbsTol(double tol) { ksp->SetAbsTol(tol); }

  // A solve at a relative tolerance above tol which converges falsely (true residual not
  // reduced, see IterativeSolver::SetTrueResidualCheck) is repeated with the relative
  // tolerance tol before the error, for relaxed solves whose true residual check is loose
  // for nearly singular systems. A tol <= 0 disables the repetition.
  void SetFallbackRelTol(double tol) { fallback_rel_tol = tol; }

  // Whether the most recent Mult converged falsely (true residual not reduced, see
  // IterativeSolver::SetTrueResidualCheck). A false convergence is an error unless allowed,
  // in which case the caller must not use the solution.
  bool GetFalseConverged() const { return false_converged; }
  void SetAllowFalseConvergence(bool allow) { allow_false_convergence = allow; }

  int GetMaxIter() const { return ksp->GetMaxIter(); }
  void SetMaxIter(int its) { ksp->SetMaxIter(its); }
  void SetInitialGuess(bool guess) { ksp->SetInitialGuess(guess); }

  // Start a solve from zero when its initial guess has a true residual no smaller than tol
  // times that of zero (GMRES, see GmresSolver::SetInitialGuessCheck).
  void SetInitialGuessCheck(double tol);

  // With reuse_coarse, a multigrid preconditioner keeps the coarse-level solver
  // (factorization) of the previous operators.
  void SetOperators(const OperType &op, const OperType &pc_op, bool reuse_coarse = false);

  // Whether the preconditioner is geometric multigrid, whose coarse-level solver
  // SetOperators can keep.
  bool HasMultigridPreconditioner() const;

  // Deflate the solves after the first one with the same operators with the k harmonic Ritz
  // vectors of the preconditioned operator with the eigenvalues closest to zero, from the
  // Krylov subspace of the first solve (complex GMRES with left preconditioning only, see
  // DeflatedPreconditioner): they use the A-DEF1 deflated preconditioner. k = 0 disables
  // the deflation; k < 0 chooses the vectors from the harmonic Ritz values: only those of
  // magnitude below 0.25 at first (at most 64, possibly none), if GMRES for the projected
  // operator predicts a gain of at least 5% at the measured cost of the deflation per
  // iteration; a solve without deflation after 8 and then every 32 deflated solves verifies
  // the gain (the deflation stops, or goes back to fewer vectors, if it finds none), and
  // while the deflation pays the vectors of the values below 0.5 and then 1, from the same
  // first solve, are tried when the model predicts at least 10% fewer iterations with them
  // (and kept if they save more than 3% on the same right-hand side). The basis is
  // computed again after each SetOperators.
  void SetDeflation(int k);

  // Number of vectors of the current deflation basis.
  int GetDeflationSize() const;

  // Solve, warning if not converged when warn_unconverged.
  void Mult(const VecType &x, VecType &y, bool warn_unconverged = true) const;
};

using KspSolver = BaseKspSolver<Operator>;
using ComplexKspSolver = BaseKspSolver<ComplexOperator>;

// Reuse of the coarse-level factorization of a multigrid preconditioner across the
// frequencies of a sweep, as an approximate coarse solve for nearby frequencies (see
// config["Solver"]["Driven"]["CoarseReuseTol"]).
class CoarseReuse
{
private:
  // Relative frequency distance for reuse (0 disables it), and allowed relative growth of
  // the number of iterations with a reused factorization.
  double tol, growth;

  // Frequency of the current coarse factorization, and the largest number of iterations of
  // the solves with a new factorization.
  double omega_ref = 0.0;
  int its_ref = 0;
  bool valid = false;

public:
  CoarseReuse(double tol, double growth = 0.5) : tol(tol), growth(growth) {}

  // Set the operators at frequency omega and solve. The coarse factorization of an earlier
  // frequency within the tolerance is kept, and the solve is then limited to (1 + growth)
  // times the largest number of iterations of the solves with a new factorization: if it
  // does not converge, the coarse level is set up at omega and the solve continues from the
  // current iterate.
  void Solve(ComplexKspSolver &ksp, const ComplexOperator &A, const ComplexOperator &P,
             double omega, const ComplexVector &b, ComplexVector &x, bool initial_guess);
};

}  // namespace palace

#endif  // PALACE_LINALG_KSP_HPP
