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
// perf::LinearSolverPath). All values are rank-uniform.
class KspStats
{
public:
  // Iteration histogram buckets: 0, 1, 2-3, 4-7, ..., 512-1023, 1024+.
  static constexpr int NUM_BUCKETS = 12;
  static int Bucket(int its);
  static std::string BucketName(int b);

  explicit KspStats(std::string_view role);

  void Record(int its, bool converged, int restarts, double initial_res,
              double final_res) const;

private:
  int solves, its, max_its, not_converged, not_converged_its, zero_its, restarts,
      mean_reduction, max_rel_res;
  std::array<int, NUM_BUCKETS> histogram;
};

// Eigenvalue solve statistics under EigenSolver (main role) or EigenSolver/<role>: outer
// iterations (restarts for Krylov-Schur and ARPACK), converged and requested eigenpairs.
void RecordEigenSolve(std::string_view role, int outer_its, int num_conv, int num_req);

// Inclusive tally for a smoother's largest eigenvalue estimate, keyed by the current role.
int LambdaMaxTally();

// Exclusive timer row a nested solver role runs within (parent for its tallies).
Timer::Index RoleParent(std::string_view role);

}  // namespace palace::perf

#endif  // PALACE_LINALG_SOLVER_STATS_HPP
