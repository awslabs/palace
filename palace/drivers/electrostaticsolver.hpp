// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_ELECTROSTATIC_SOLVER_HPP
#define PALACE_DRIVERS_ELECTROSTATIC_SOLVER_HPP

#include <map>
#include <memory>
#include <vector>
#include "drivers/basesolver.hpp"
#include "linalg/vector.hpp"

namespace mfem
{

template <typename T>
class Array;
class DenseMatrix;

}  // namespace mfem

namespace palace
{

class ErrorIndicator;
class IoData;
class LaplaceOperator;
class Mesh;
class PostOperator;
class Timer;

//
// Driver class for electrostatic simulations.
//
class ElectrostaticSolver : public BaseSolver
{
private:
  ErrorIndicator Postprocess(LaplaceOperator &laplaceop, PostOperator &postop,
                             const std::vector<Vector> &V) const;

  void PostprocessTerminals(const std::map<int, mfem::Array<int>> &terminal_sources,
                            const mfem::DenseMatrix &C, const mfem::DenseMatrix &Cinv,
                            const mfem::DenseMatrix &Cm) const;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

public:
  using BaseSolver::BaseSolver;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_ELECTROSTATIC_SOLVER_HPP
