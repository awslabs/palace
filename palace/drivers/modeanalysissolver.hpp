// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_MODE_ANALYSIS_SOLVER_HPP
#define PALACE_DRIVERS_MODE_ANALYSIS_SOLVER_HPP

#include "drivers/basesolver.hpp"

namespace palace
{

//
// Driver class for 2D waveguide mode analysis. Solves for propagation modes of a
// waveguide cross-section at fixed frequency, computing propagation constant (kn),
// effective index (n_eff), characteristic impedance (Z0), and per-unit-length
// inductance and capacitance.
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
