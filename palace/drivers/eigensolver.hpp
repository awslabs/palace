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

struct CoarseSolveChoice;
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
  // rebalancing of the finest mesh for the starting vector of the next eigenvalue solve.
  mutable std::unique_ptr<VectorFieldTransfer> prev_eigenvectors;
  mutable std::vector<std::complex<double>> prev_eigenvalues;

  // Number of eigenvalue solver iterations (restarts) of the previous adaptation iteration.
  mutable int prev_num_its = -1;

  // Converged eigenvalues λ = iω of the previous adaptation iteration (the requested ones
  // first), for the shift of the spectral transformation of the next one, the number of
  // iterations (restarts) of its solve, and whether its shift was apart from the target.
  mutable std::vector<std::complex<double>> prev_shift_eigenvalues;
  mutable int prev_shift_its = -1;
  mutable bool prev_shifted = false;

  // Whether the last solve used the tolerances of intermediate adaptation iterations.
  mutable bool intermediate_solve = false;

  // Choice of the coarse solve of the multigrid preconditioner made for an earlier shift
  // or adaptation iteration, reused by the next ones (see GeometricMultigridSolver).
  mutable std::shared_ptr<CoarseSolveChoice> coarse_choice;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

  void UpdateAfterMeshChange(Mesh &mesh) const override;

  bool SolvedWithIntermediateTolerances() const override { return intermediate_solve; }

public:
  EigenSolver(const IoData &iodata, bool root, int size = 0, int num_thread = 0,
              const char *git_tag = nullptr);
  ~EigenSolver() override;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_EIGEN_SOLVER_HPP
