// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_EIGEN_SOLVER_HPP
#define PALACE_DRIVERS_EIGEN_SOLVER_HPP

#include <complex>
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
  // Eigenpairs of the previous adaptation iteration, kept through the refinement and
  // rebalancing of the finest mesh for the starting vector of the next eigenvalue solve,
  // and the number of iterations (restarts) of that solve.
  mutable std::unique_ptr<VectorFieldTransfer> prev_eigenvectors;
  mutable std::vector<std::complex<double>> prev_eigenvalues;
  mutable int prev_num_its = -1;

  // Number of restarts of the last solve which was not warm-started.
  mutable int cold_num_its = -1;

  // Number of solves so far: the initial one, then one per adaptation iteration.
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
