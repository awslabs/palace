// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_MAGNETOSTATIC_SOLVER_HPP
#define PALACE_DRIVERS_MAGNETOSTATIC_SOLVER_HPP

#include <memory>
#include <vector>
#include "drivers/basesolver.hpp"
#include "linalg/vector.hpp"
#include "utils/configfile.hpp"

namespace palace
{

class ErrorIndicator;
class Mesh;
template <ProblemType>
class PostOperator;
class SurfaceCurrentOperator;

//
// Driver class for magnetostatic simulations.
//
class MagnetostaticSolver : public BaseSolver
{
private:
  void PostprocessTerminals(PostOperator<ProblemType::MAGNETOSTATIC> &post_op,
                            const SurfaceCurrentOperator &surf_j_op,
                            const std::vector<Vector> &A,
                            const std::vector<double> &I_inc) const;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

public:
  using BaseSolver::BaseSolver;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_MAGNETOSTATIC_SOLVER_HPP
