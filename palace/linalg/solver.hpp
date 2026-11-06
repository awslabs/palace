// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_SOLVER_HPP
#define PALACE_LINALG_SOLVER_HPP

#include <string>
#include <type_traits>
#include <mfem.hpp>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "utils/timer.hpp"

namespace palace
{

//
// The base Solver<OperType> class is a templated version of mfem::Solver for operation with
// real- or complex-valued operators.
//

// Abstract base class for real-valued or complex-valued solvers.
template <OperatorType OperType>
class Solver : public OperType
{
protected:
  using VecType =
      std::conditional_t<std::is_same_v<OperType, ComplexOperator>, ComplexVector, Vector>;

  // Whether or not to use the second argument of Mult() as an initial guess.
  bool initial_guess;

public:
  Solver(bool initial_guess = false) : OperType(), initial_guess(initial_guess) {}
  virtual ~Solver() = default;

  // Configure whether or not to use an initial guess when applying the solver.
  virtual void SetInitialGuess(bool guess) { initial_guess = guess; }

  // Set the operator associated with the solver, or update it if called repeatedly.
  virtual void SetOperator(const OperType &op) = 0;

  // Register inclusive performance tallies named under prefix, measured within the given
  // exclusive timer row. No-op by default.
  virtual void SetPerfTallies(const std::string &prefix, Timer::Index parent) {}

  // Apply the solver for the transpose problem.
  void MultTranspose(const VecType &x, VecType &y) const override
  {
    MFEM_ABORT("MultTranspose() is not implemented for base class Solver<OperType>!");
  }

  // Apply the solver with a preallocated temporary storage vector.
  virtual void Mult2(const VecType &x, VecType &y, VecType &r) const
  {
    MFEM_ABORT("Mult2() with temporary storage vector is not implemented for base class "
               "Solver<OperType>!");
  }

  // Apply the solver for the transpose problem with a preallocated temporary storage
  //  vector.
  virtual void MultTranspose2(const VecType &x, VecType &y, VecType &r) const
  {
    MFEM_ABORT("MultTranspose2() with temporary storage vector is not implemented for base "
               "class Solver<OperType>!");
  }
};

// This solver wraps a real-valued mfem::Solver for application to complex-valued problems
// as a preconditioner inside of a Solver<OperType> or for assembling the matrix-free
// preconditioner operator as an mfem::HypreParMatrix.
template <OperatorType OperType>
class MfemWrapperSolver : public Solver<OperType>
{
  using VecType = typename Solver<OperType>::VecType;

private:
  // The actual mfem::Solver.
  std::unique_ptr<mfem::Solver> pc;

  // System matrix A in parallel assembled form.
  std::unique_ptr<mfem::HypreParMatrix> A;

  // Whether or not to save the parallel assembled matrix after calling
  // mfem::Solver::SetOperator (some solvers copy their input).
  bool save_assembled;

  // Whether to use the exact complex-valued system matrix or the real-valued
  // approximation A = Ar + Ai.
  bool complex_matrix = true;

  // Whether to drop small entries (< ε) in the sparse system matrix.
  bool drop_small_entries = true;

  // Whether to reuse the column reordering of previous factorizations.
  bool reorder_reuse = true;

  // Number of small entries dropped by the most recent DropSmallEntries() call.
  int num_dropped_entries = 0;

  // Drop small entries.
  void DropSmallEntries();

  // Performance reporting: the first application after SetOperator (where lazy solvers
  // factorize or set up) is timed as Timer::KSP_COARSE_SETUP, later ones as apply_timer.
  // Solvers which factorize in SetOperator (MUMPS) time that as Timer::KSP_COARSE_SETUP.
  bool time_setup = false, eager_setup = false;
  Timer::Index apply_timer = Timer::NUM_TIMINGS;
  mutable bool setup_pending = false;
  int num_setups = 0;
  int perf_setups = -1, perf_refact, perf_applies, perf_rows, perf_rows_rank;

  // Call the wrapped solver's SetOperator.
  void SetInnerOperator(const mfem::HypreParMatrix &op);

public:
  MfemWrapperSolver(std::unique_ptr<mfem::Solver> &&pc, bool save_assembled = true,
                    bool complex_matrix = true, bool drop_small_entries = true,
                    bool reorder_reuse = true)
    : Solver<OperType>(pc->iterative_mode), pc(std::move(pc)),
      save_assembled(save_assembled), complex_matrix(complex_matrix),
      drop_small_entries(drop_small_entries), reorder_reuse(reorder_reuse)
  {
    InitEagerSetup();
  }

  // Enable timing (and coarse solver counters, keyed by the current solver role), see
  // above. Applications after the first are timed as apply_timer (NUM_TIMINGS: not timed).
  void EnableTimer(bool setup, Timer::Index apply = Timer::NUM_TIMINGS);
  void InitEagerSetup();

  // Access the underlying solver.
  const mfem::Solver &GetSolver() { return *pc; }

  // Configure whether or not to save the assembled operator.
  void SetSaveAssembled(bool save) { save_assembled = save; }

  // Configure whether or not to drop small entries in the system matrix.
  void SetDropSmallEntries(bool drop) { drop_small_entries = drop; }

  // Configure whether or not to use the exact complex-valued system matrix
  // or the real-valued approximation.
  void SetComplexMatrix(bool use_complex) { complex_matrix = use_complex; }

  void SetInitialGuess(bool guess) override
  {
    Solver<OperType>::SetInitialGuess(guess);
    pc->iterative_mode = guess;
  }

  void SetOperator(const OperType &op) override;

  void Mult(const VecType &x, VecType &y) const override;
};

}  // namespace palace

#endif  // PALACE_LINALG_SOLVER_HPP
