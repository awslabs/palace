// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_EIGEN_SOLVER_HPP
#define PALACE_DRIVERS_EIGEN_SOLVER_HPP

#include <complex>
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
class Timer;

//
// Driver class for eigenmode simulations.
//
class EigenSolver : public BaseSolver
{
private:
  void Postprocess(const PostOperator &postop, const LumpedPortOperator &lumped_port_op,
                   int i, std::complex<double> omega, double error1, double error2,
                   int num_conv, Timer &timer) const;

  void PostprocessEigen(int i, std::complex<double> omega, double error1, double error2,
                        int num_conv) const;
  void PostprocessEPR(const PostOperator &postop, const LumpedPortOperator &lumped_port_op,
                      int i, std::complex<double> omega, double Em) const;

public:
  EigenSolver(const IoData &iodata, bool root, int size = 0, int num_thread = 0,
              const char *git_tag = nullptr)
    : BaseSolver(iodata, root, size, num_thread, git_tag)
  {
  }

  void Solve(std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
             Timer &timer) const override;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_EIGEN_SOLVER_HPP
