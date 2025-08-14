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
