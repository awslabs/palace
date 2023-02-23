// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_EIGEN_SOLVER_HPP
#define PALACE_EIGEN_SOLVER_HPP

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
  void Postprocess(const std::string &post_dir, const PostOperator &postop,
                   const LumpedPortOperator &lumped_port_op, int i,
                   std::complex<double> omega, double error1, double error2, int num_conv,
                   Timer &timer) const;

  void PostprocessEigen(const std::string &post_dir, int i, std::complex<double> omega,
                        double error1, double error2, int num_conv) const;
  void PostprocessEPR(const std::string &post_dir, const PostOperator &postop,
                      const LumpedPortOperator &lumped_port_op, int i,
                      std::complex<double> omega, double Em) const;

public:
  using BaseSolver::BaseSolver;

  BaseSolver::ErrorIndicators Solve(std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
                                    Timer &timer, int iter) const override;
};

}  // namespace palace

#endif  // PALACE_EIGEN_SOLVER_HPP
