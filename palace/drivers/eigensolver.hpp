// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_EIGEN_SOLVER_HPP
#define PALACE_DRIVERS_EIGEN_SOLVER_HPP

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

  // Warm starts of the next solves, and number of solves so far (the initial one, then one
  // per adaptation iteration).
  mutable WarmStartGate warm_start_gate;
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
