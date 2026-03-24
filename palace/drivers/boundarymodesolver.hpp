// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_BOUNDARY_MODE_SOLVER_HPP
#define PALACE_DRIVERS_BOUNDARY_MODE_SOLVER_HPP

#include "drivers/basesolver.hpp"

namespace palace
{

//
// Driver class for 2D waveguide mode analysis. Constructs a BoundaryModeOperator (which
// owns the mesh, FE spaces, material operator, and eigenvalue solver), then orchestrates
// the solve and postprocessing.
//
class BoundaryModeSolver : public BaseSolver
{
public:
  using BaseSolver::BaseSolver;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_BOUNDARY_MODE_SOLVER_HPP
