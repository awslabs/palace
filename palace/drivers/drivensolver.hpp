// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_DRIVEN_SOLVER_HPP
#define PALACE_DRIVERS_DRIVEN_SOLVER_HPP

#include <memory>
#include <vector>
#include "drivers/basesolver.hpp"

namespace mfem
{

class ParMesh;

}  // namespace mfem

namespace palace
{

class CurlFluxErrorEstimator;
class ErrorIndicator;
class IoData;
class LumpedPortOperator;
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

  ErrorIndicator SweepUniform(SpaceOperator &spaceop, PostOperator &postop,
                              CurlFluxErrorEstimator &estimator, int nstep, int step0,
                              double omega0, double delta_omega) const;
  ErrorIndicator SweepAdaptive(SpaceOperator &spaceop, PostOperator &postop,
                               CurlFluxErrorEstimator &estimator, int nstep, int step0,
                               double omega0, double delta_omega) const;

  void Postprocess(const PostOperator &postop, const LumpedPortOperator &lumped_port_op,
                   const WavePortOperator &wave_port_op,
                   const SurfaceCurrentOperator &surf_j_op, int step, double omega,
                   double E_elec, double E_mag, bool full) const;

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

public:
  using BaseSolver::BaseSolver;

  ErrorIndicator Solve(const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh) const final;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_DRIVEN_SOLVER_HPP
