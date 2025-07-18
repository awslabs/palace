// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_DRIVEN_SOLVER_HPP
#define PALACE_DRIVERS_DRIVEN_SOLVER_HPP

#include <memory>
#include <vector>
#include "drivers/basesolver.hpp"
#include "utils/configfile.hpp"

namespace palace
{

class ErrorIndicator;
class Mesh;
template <ProblemType>
class PostOperator;
class SpaceOperator;

// Mini helper class that stored indexing information for printing and restart.
struct DrivenSolverIndexing
{
  std::size_t nr_total_samples;

  // Restart (1-Based indexing from config file)
  std::size_t restart = 1;

  // Offsets from restart (0-Based)
  std::size_t excitation_n0 = 0;
  std::size_t omega_n0 = 0;

  DrivenSolverIndexing(std::size_t nr_port_excitations, std::size_t nr_freq_samples,
                       std::size_t restart_)
    : nr_total_samples(nr_port_excitations * nr_freq_samples), restart(restart_),
      excitation_n0(std::size_t(restart - 1) / nr_freq_samples),
      omega_n0(std::size_t(restart - 1) % nr_freq_samples)
  {
  }
};

//
// Driver class for driven terminal simulations.
//
class DrivenSolver : public BaseSolver
{
private:
  ErrorIndicator SweepUniform(SpaceOperator &space_op) const;

  ErrorIndicator SweepAdaptive(SpaceOperator &space_op) const;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

public:
  using BaseSolver::BaseSolver;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_DRIVEN_SOLVER_HPP
