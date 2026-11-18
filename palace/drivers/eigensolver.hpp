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

  // Number of Krylov cycles of the last solve without a warm start, and number of solves so
  // far (the initial one, then one per adaptation iteration).
  mutable int last_cold_its = -1;
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
