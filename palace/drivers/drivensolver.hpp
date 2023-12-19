// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_DRIVEN_SOLVER_HPP
#define PALACE_DRIVERS_DRIVEN_SOLVER_HPP

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
class SpaceOperator;
class SurfaceCurrentOperator;
class Timer;
class WavePortOperator;

//
// Driver class for driven terminal simulations.
//
class DrivenSolver : public BaseSolver
{
private:
  int GetNumSteps(double start, double end, double delta) const;

  ErrorIndicator SweepUniform(SpaceOperator &spaceop, PostOperator &postop, int nstep,
                              int step0, double omega0, double delta_omega) const;

  ErrorIndicator SweepAdaptive(SpaceOperator &spaceop, PostOperator &postop, int nstep,
                               int step0, double omega0, double delta_omega) const;

  void Postprocess(const PostOperator &postop, const LumpedPortOperator &lumped_port_op,
                   const WavePortOperator &wave_port_op,
                   const SurfaceCurrentOperator &surf_j_op, int step, double omega,
                   double E_elec, double E_mag, bool full,
                   const ErrorIndicator *indicator) const;

  void PostprocessCurrents(const PostOperator &postop,
                           const SurfaceCurrentOperator &surf_j_op, int step,
                           double omega) const;

  void PostprocessPorts(const PostOperator &postop,
                        const LumpedPortOperator &lumped_port_op, int step,
                        double omega) const;

  void PostprocessSParameters(const PostOperator &postop,
                              const LumpedPortOperator &lumped_port_op,
                              const WavePortOperator &wave_port_op, int step,
                              double omega) const;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

public:
  using BaseSolver::BaseSolver;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_DRIVEN_SOLVER_HPP
