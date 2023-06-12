// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_TRANSIENT_SOLVER_HPP
#define PALACE_DRIVERS_TRANSIENT_SOLVER_HPP

#include <functional>
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
                   double J_coef, double E_elec, double E_mag, bool ful,
                   Timer &timer) const;

  void PostprocessCurrents(const PostOperator &postop,
                           const SurfaceCurrentOperator &surf_j_op, int step, double t,
                           double J_coef) const;
  void PostprocessPorts(const PostOperator &postop,
                        const LumpedPortOperator &lumped_port_op, int step, double t,
                        double J_coef) const;

public:
  TransientSolver(const IoData &iodata, bool root, int size = 0, int num_thread = 0,
                  const char *git_tag = nullptr)
    : BaseSolver(iodata, root, size, num_thread, git_tag)
  {
  }

  void Solve(std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
             Timer &timer) const override;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_TRANSIENT_SOLVER_HPP
