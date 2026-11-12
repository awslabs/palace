// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_SOLVER_STATS_HPP
#define PALACE_LINALG_SOLVER_STATS_HPP

#include <array>
#include <string>
#include <string_view>
#include "utils/timer.hpp"

namespace palace::perf
{

// Linear solve statistics for one solver role, written under LinearSolver/<role> (see
// perf::LinearSolverPath). All values are uniform over the processes which run the role's
// solves (all of them, except for solvers on a subcommunicator, such as the mode solver of
// each wave port, which has a role of its own): they are reduced as such.
class KspStats
{
public:
  // Iteration histogram buckets: 0-0, 1-1, 2-3, 4-7, ..., 512-1023, 1024+.
  static constexpr int NUM_BUCKETS = 12;
  static int Bucket(int its);
  static std::string BucketName(int b);

  explicit KspStats(std::string_view role);

  // Per solve: iterations, convergence, restarts, and the initial and final residual norms
  // (the reduction factor of a solve with n > 0 iterations is (final / initial)^(1 / n),
  // the geometric mean of the residual reduction per iteration).
  void Record(int its, bool converged, int restarts, double initial_res,
              double final_res) const;

private:
  int solves, its, max_its, not_converged, not_converged_its, zero_its, restarts,
      mean_reduction, min_reduction, max_reduction, max_rel_res;
  std::array<int, NUM_BUCKETS> histogram;
};

// Eigenvalue solve statistics under EigenSolver (main role) or EigenSolver/<role>: outer
// iterations (restarts for Krylov-Schur and ARPACK), converged and requested eigenpairs.
void RecordEigenSolve(std::string_view role, int outer_its, int num_conv, int num_req);

// Inclusive tally for a smoother's largest eigenvalue estimates, keyed by the current role:
// register it at the construction of the smoother, whose estimates may run outside the
// role's scope.
int LambdaMaxTally();

}  // namespace palace::perf

#endif  // PALACE_LINALG_SOLVER_STATS_HPP
