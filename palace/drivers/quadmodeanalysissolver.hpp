// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_QUAD_MODE_ANALYSIS_SOLVER_HPP
#define PALACE_DRIVERS_QUAD_MODE_ANALYSIS_SOLVER_HPP

#include "drivers/basesolver.hpp"

namespace palace
{

//
// Driver class for 2D waveguide mode analysis using the quadratic eigenvalue problem (QEP)
// formulation. Solves for propagation modes with full impedance BC support on both et and
// en components. Uses the same JSON config as ModeAnalysisSolver (Problem.Type =
// "ModeAnalysis") but internally constructs a QuadBoundaryModeSolver instead.
//
class QuadModeAnalysisSolver : public BaseSolver
{
public:
  using BaseSolver::BaseSolver;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_QUAD_MODE_ANALYSIS_SOLVER_HPP
