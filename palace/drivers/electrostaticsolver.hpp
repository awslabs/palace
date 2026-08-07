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
protected:
  void PostprocessTerminals(PostOperator<ProblemType::ELECTROSTATIC> &post_op,
                            const std::map<int, mfem::Array<int>> &terminal_sources,
                            const std::vector<Vector> &V) const;

  // Inner solve: fills V with the per-terminal potentials for a prebuilt operator, so
  // tests/MMS can inspect the solved field. No postprocessing or error estimation.
  void Solve(std::vector<Vector> &V, LaplaceOperator &laplace_op) const;

  // Dispatched entry point: builds the operator, runs the inner solve, then postprocesses.
  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

public:
  using BaseSolver::BaseSolver;
};

// Test-only subclass: re-exposes the protected Solve overloads as public for unit tests.
class ExposedElectrostaticSolver : public ElectrostaticSolver
{
public:
  using ElectrostaticSolver::ElectrostaticSolver;
  using ElectrostaticSolver::Solve;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_ELECTROSTATIC_SOLVER_HPP
