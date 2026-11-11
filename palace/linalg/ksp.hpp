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

  // Verification of the error propagator of the preconditioner (main solver only): the
  // operator to check after the next solve with at least one iteration, for the first
  // setup which is solved (once per solver object, which is created anew in each adaptive
  // mesh refinement iteration).
  mutable bool check_pc = true;
  mutable const OperType *check_op = nullptr;

public:
  BaseKspSolver(const config::LinearSolverData &linear, MatrixSymmetry pc_mat_sym,
                int verbose, FiniteElementSpaceHierarchy &fespaces,
                FiniteElementSpaceHierarchy *aux_fespaces = nullptr);
  BaseKspSolver(const IoData &iodata, FiniteElementSpaceHierarchy &fespaces,
                FiniteElementSpaceHierarchy *aux_fespaces = nullptr);
  // The constructors above configure the solver of the role active at construction
  // (perf::CurrentRole(), "Main" by default). Other solvers name their role, under which
  // the performance report keeps their statistics.
  BaseKspSolver(std::unique_ptr<IterativeSolver<OperType>> &&ksp,
                std::unique_ptr<Solver<OperType>> &&pc, std::string_view role);

  int NumTotalMult() const { return ksp_mult; }
  int NumTotalMultIterations() const { return ksp_mult_it; }

  const std::string &GetRole() const { return role; }

  // Whether the most recent Mult converged to the requested tolerance.
  bool GetConverged() const { return ksp->GetConverged(); }

  // Forward tolerance access to the underlying iterative solver.
  double GetRelTol() const { return ksp->GetRelTol(); }
  double GetAbsTol() const { return ksp->GetAbsTol(); }
  void SetRelTol(double tol) { ksp->SetRelTol(tol); }
  void SetAbsTol(double tol) { ksp->SetAbsTol(tol); }

  void SetOperators(const OperType &op, const OperType &pc_op);

  void Mult(const VecType &x, VecType &y) const;
};

using KspSolver = BaseKspSolver<Operator>;
using ComplexKspSolver = BaseKspSolver<ComplexOperator>;

}  // namespace palace

#endif  // PALACE_LINALG_KSP_HPP
