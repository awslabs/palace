// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_FLUX_ERROR_ESTIMATOR_HPP
#define PALACE_DRIVERS_FLUX_ERROR_ESTIMATOR_HPP

#include <memory>
#include <type_traits>
#include "linalg/errorestimator.hpp"
#include "models/spaceoperator.hpp"
#include "utils/iodata.hpp"

namespace palace
{

// Flux error estimator for the problem dimension: boundary mode estimator in 2D, time
// dependent (3D) estimator otherwise.
template <typename VecType>
class FluxErrorEstimator
{
  static_assert(std::is_same_v<VecType, Vector> || std::is_same_v<VecType, ComplexVector>,
                "FluxErrorEstimator can only be defined for VecType = Vector or "
                "ComplexVector!");

  std::unique_ptr<TimeDependentFluxErrorEstimator<VecType>> estimator_3d;
  std::unique_ptr<BoundaryModeFluxErrorEstimator<VecType>> estimator_2d;

public:
  FluxErrorEstimator(SpaceOperator &space_op, const IoData &iodata)
  {
    const bool is_2d = (space_op.GetNDSpace().Dimension() < 3);
    if (is_2d)
    {
      estimator_2d = std::make_unique<BoundaryModeFluxErrorEstimator<VecType>>(
          space_op.GetMaterialOp(), space_op.GetNDSpaces(), space_op.GetRTSpaces(),
          space_op.GetCurlSpace(), space_op.GetH1Spaces(),
          iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
          iodata.solver.linear.estimator_mg);
    }
    else
    {
      estimator_3d = std::make_unique<TimeDependentFluxErrorEstimator<VecType>>(
          space_op.GetMaterialOp(), space_op.GetNDSpaces(), space_op.GetRTSpaces(),
          iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
          iodata.solver.linear.estimator_mg);
    }
  }

  void AddErrorIndicator(const VecType &E, const VecType &B, double Et,
                         ErrorIndicator &indicator) const
  {
    if (estimator_2d)
    {
      estimator_2d->AddErrorIndicator(E, B, Et, indicator);
    }
    else
    {
      estimator_3d->AddErrorIndicator(E, B, Et, indicator);
    }
  }
};

}  // namespace palace

#endif  // PALACE_DRIVERS_FLUX_ERROR_ESTIMATOR_HPP
