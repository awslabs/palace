// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_EIGEN_SOLVER_HPP
#define PALACE_DRIVERS_EIGEN_SOLVER_HPP

#include <complex>
#include <memory>
#include <vector>
#include "drivers/basesolver.hpp"

namespace palace
{

class ErrorIndicator;
class LumpedPortOperator;
class Mesh;
class PostOperator;

//
// Driver class for eigenmode simulations.
//
class EigenSolver : public BaseSolver
{
private:
  void Postprocess(const PostOperator &post_op, const LumpedPortOperator &lumped_port_op,
                   int i, std::complex<double> omega, double error_bkwd, double error_abs,
                   int num_conv, double E_elec, double E_mag,
                   const ErrorIndicator *indicator) const;

  void PostprocessEigen(int i, std::complex<double> omega, double error_bkwd,
                        double error_abs, int num_conv) const;

  void PostprocessPorts(const PostOperator &post_op,
                        const LumpedPortOperator &lumped_port_op, int i) const;

  void PostprocessEPR(const PostOperator &post_op, const LumpedPortOperator &lumped_port_op,
                      int i, std::complex<double> omega, double E_m) const;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

public:
  using BaseSolver::BaseSolver;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_EIGEN_SOLVER_HPP
