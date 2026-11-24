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
private:
  void PostprocessTerminals(PostOperator<ProblemType::ELECTROSTATIC> &post_op,
                            const std::map<int, mfem::Array<int>> &terminal_sources,
                            const std::vector<Vector> &V) const;

  // The equilibrated error estimate of the terminal solves (see
  // EquilibratedFluxErrorEstimator), its guaranteed bounds of the relative errors and its
  // guaranteed brackets of the entries of the capacitance and mutual capacitance matrices
  // (log, and terminal-C-lower.csv, terminal-C-upper.csv, terminal-Cm-lower.csv and
  // terminal-Cm-upper.csv).
  void EstimateEquilibrated(LaplaceOperator &laplace_op, const std::vector<Vector> &V,
                            const std::vector<Vector> &G,
                            const std::vector<double> &energies,
                            ErrorIndicator &indicator) const;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

public:
  using BaseSolver::BaseSolver;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_ELECTROSTATIC_SOLVER_HPP
