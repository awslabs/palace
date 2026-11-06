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

class ComplexOperator;
class ComplexVector;
template <typename VecType>
class DivFreeSolver;
class EigenvalueSolver;
class ErrorIndicator;
class Mesh;
class SpaceOperator;

//
// Driver class for eigenmode simulations.
//
class EigenSolver : public BaseSolver
{
private:
  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

  // Compute and write the derivatives of the frequencies and quality factors of the
  // configured objectives, for the modes of angular frequencies omega.
  void PostprocessSensitivities(SpaceOperator &space_op, const EigenvalueSolver &eigen,
                                const ComplexOperator *C, const ComplexOperator &M,
                                bool has_A2, const DivFreeSolver<ComplexVector> *divfree,
                                const std::vector<std::complex<double>> &omega) const;

public:
  using BaseSolver::BaseSolver;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_EIGEN_SOLVER_HPP
