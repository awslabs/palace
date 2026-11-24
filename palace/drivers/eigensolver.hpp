// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_EIGEN_SOLVER_HPP
#define PALACE_DRIVERS_EIGEN_SOLVER_HPP

#include <algorithm>
#include <memory>
#include <vector>
#include "drivers/basesolver.hpp"

namespace palace
{

class ErrorIndicator;
class Mesh;
class VectorFieldTransfer;

//
// Decides whether the eigenvalue solve of the next adaptation iteration starts from the
// eigenvectors of the current one (warm start). The more accurate first Krylov cycle of a
// warm start can cost more than it saves on a solve of two cycles, which it shortens by one
// cycle at most: the warm start follows the cold solves of three cycles or more. With the
// quality factor region, which shortens the cold solves that stall behind strongly damped
// Ritz values to two cycles, it also follows those of two cycles, as long as the
// warm-started solves after them converge in a single cycle.
//
class WarmStartGate
{
private:
  // Number of Krylov cycles of the last solve without a warm start (-1 before the first),
  // and whether a warm start still follows the cold solves of two cycles.
  int last_cold_its = -1;
  bool after_two_cycles = true;

public:
  // Whether the next solve is warm-started, with or without the quality factor region.
  bool Next(bool q_region) const
  {
    return last_cold_its >= 3 || (q_region && last_cold_its == 2 && after_two_cycles);
  }

  // Records a solve: whether it was warm-started, and its number of Krylov cycles.
  void Update(bool warm, int its)
  {
    if (!warm)
    {
      last_cold_its = its;
    }
    else if (last_cold_its < 3 && its > 1)
    {
      // The warm start after a cold solve of two cycles saved no cycle.
      after_two_cycles = false;
    }
  }

  int GetLastColdIterations() const { return last_cold_its; }
};

//
// Decides the Krylov subspace dimension of the eigenvalue solves of an adaptive run. On a
// dense cluster of modes, the solves with the default dimension, max(2 N, N + 15), can
// stall for tens of Krylov cycles, and the quality factor region changes for how long,
// while a few more basis vectors resolve the cluster in a third of the cycles. The solves
// which converge in a few cycles pay for more vectors, since every cycle does more linear
// solves: with the region, the solves after one which took more than 10 cycles use 6 more
// vectors than the default, unless another dimension is set.
//
class SubspaceGrowth
{
private:
  bool grow = false;

public:
  static constexpr int max_cycles = 10, extra_vectors = 6;

  // The subspace dimension of the next solve for num_modes modes, from the configured one
  // (the configuration resolves an unset dimension to the default), at most the size of
  // the problem.
  int GetMaxSize(int num_modes, int max_size, long long int size) const
  {
    const int n = std::max(2 * num_modes, num_modes + 15);
    if (!grow || max_size != n)
    {
      return max_size;
    }
    return static_cast<int>(std::min(static_cast<long long int>(n + extra_vectors), size));
  }

  // Records a solve: whether it used the quality factor region, and its number of Krylov
  // cycles.
  void Update(bool q_region, int its) { grow = grow || (q_region && its > max_cycles); }
};

//
// Driver class for eigenmode simulations.
//
class EigenSolver : public BaseSolver
{
private:
  // Combination of the eigenvectors of the previous adaptation iteration (and of λ x for
  // quadratic problems), kept through the refinement and rebalancing of the finest mesh to
  // start the next eigenvalue solve, and the number of eigenvectors in it.
  mutable std::unique_ptr<VectorFieldTransfer> prev_eigenvectors;
  mutable int prev_num_eigenvectors = 0;

  // Warm starts and subspace dimensions of the next solves, and number of solves so far
  // (the initial one, then one per adaptation iteration).
  mutable WarmStartGate warm_start_gate;
  mutable SubspaceGrowth subspace_growth;
  mutable int num_solves = 0;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

  void UpdateAfterMeshChange(Mesh &mesh) const override;

public:
  EigenSolver(const IoData &iodata, bool root, int size = 0, int num_thread = 0,
              const char *git_tag = nullptr);
  ~EigenSolver() override;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_EIGEN_SOLVER_HPP
