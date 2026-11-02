// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_ELECTROSTATIC_SOLVER_HPP
#define PALACE_DRIVERS_ELECTROSTATIC_SOLVER_HPP

#include <map>
#include <memory>
#include <vector>
#include "drivers/basesolver.hpp"
#include "linalg/vector.hpp"
#include "utils/configfile.hpp"

namespace mfem
{

template <typename T>
class Array;
class DenseMatrix;

}  // namespace mfem

namespace palace
{

class ErrorIndicator;
class LaplaceOperator;
class Mesh;
template <ProblemType>
class PostOperator;

//
// Driver class for electrostatic simulations.
//
class ElectrostaticSolver : public BaseSolver
{
private:
  // Compute and write the capacitance matrix, returned in nondimensional units.
  mfem::DenseMatrix
  PostprocessTerminals(PostOperator<ProblemType::ELECTROSTATIC> &post_op,
                       const std::map<int, mfem::Array<int>> &terminal_sources,
                       const std::vector<Vector> &V) const;

  // Compute and write the derivatives of the capacitances of the configured objectives.
  void PostprocessSensitivities(const LaplaceOperator &laplace_op,
                                const std::map<int, mfem::Array<int>> &terminal_sources,
                                const std::vector<Vector> &V,
                                const mfem::DenseMatrix &C) const;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

public:
  using BaseSolver::BaseSolver;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_ELECTROSTATIC_SOLVER_HPP
