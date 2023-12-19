// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_TRANSIENT_SOLVER_HPP
#define PALACE_DRIVERS_TRANSIENT_SOLVER_HPP

#include <functional>
#include <memory>
#include <vector>
#include "drivers/basesolver.hpp"

namespace palace
{

class ErrorIndicator;
class IoData;
class LumpedPortOperator;
class Mesh;
class PostOperator;
class SurfaceCurrentOperator;
class Timer;

//
// Driver class for time-dependent driven terminal simulations.
//
class TransientSolver : public BaseSolver
{
private:
  std::function<double(double)> GetTimeExcitation(bool dot) const;

  int GetNumSteps(double start, double end, double delta) const;

  void Postprocess(const PostOperator &postop, const LumpedPortOperator &lumped_port_op,
                   const SurfaceCurrentOperator &surf_j_op, int step, double t,
                   double J_coef, double E_elec, double E_mag, bool full,
                   const ErrorIndicator *indicator) const;

  void PostprocessCurrents(const PostOperator &postop,
                           const SurfaceCurrentOperator &surf_j_op, int step, double t,
                           double J_coef) const;

  void PostprocessPorts(const PostOperator &postop,
                        const LumpedPortOperator &lumped_port_op, int step, double t,
                        double J_coef) const;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

public:
  using BaseSolver::BaseSolver;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_TRANSIENT_SOLVER_HPP
