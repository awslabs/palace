// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_MODE_ANALYSIS_SOLVER_HPP
#define PALACE_DRIVERS_MODE_ANALYSIS_SOLVER_HPP

#include "drivers/basesolver.hpp"

namespace palace
{

//
// Driver class for 2D waveguide mode analysis using a linear eigenvalue formulation
// (Eq 1 + Eq 2 with VD substitution). This formulation supports full impedance BC
// handling (both BC-t and BC-n) while maintaining a standard generalized eigenvalue
// problem in kn^2 (no quadratic linearization).
//
class ModeAnalysisSolver : public BaseSolver
{
public:
  using BaseSolver::BaseSolver;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_MODE_ANALYSIS_SOLVER_HPP
