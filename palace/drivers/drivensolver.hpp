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

  void SweepUniform(SpaceOperator &spaceop, PostOperator &postop, int nstep, int step0,
                    double omega0, double delta_omega, Timer &timer) const;
  void SweepAdaptive(SpaceOperator &spaceop, PostOperator &postop, int nstep, int step0,
                     double omega0, double delta_omega, Timer &timer) const;

  void Postprocess(const PostOperator &postop, const LumpedPortOperator &lumped_port_op,
                   const WavePortOperator &wave_port_op,
                   const SurfaceCurrentOperator &surf_j_op, int step, double omega,
                   double E_elec, double E_mag, bool full, Timer &timer) const;

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
  DrivenSolver(const IoData &iodata, bool root, int size = 0, int num_thread = 0,
               const char *git_tag = nullptr)
    : BaseSolver(iodata, root, size, num_thread, git_tag)
  {
  }

  void Solve(std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
             Timer &timer) const override;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_DRIVEN_SOLVER_HPP
