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
class SurfacePostGeometry;
class SurfaceResponseGeometry;

//
// Driver class for eigenmode simulations.
//
class EigenSolver : public BaseSolver
{
private:
  mutable std::shared_ptr<const SurfacePostGeometry> surface_post_geometry;
  mutable std::shared_ptr<const SurfaceResponseGeometry> response_geometry;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

public:
  using BaseSolver::BaseSolver;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_EIGEN_SOLVER_HPP
