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

  // The Maxwell, inverse and mutual capacitance matrices from C (of the given terminals),
  // written to terminal-C.csv, terminal-Cinv.csv and terminal-Cm.csv with the terminal
  // excitations.
  void PostprocessCapacitance(const std::vector<int> &terminal_idx,
                              const mfem::DenseMatrix &C) const;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

public:
  using BaseSolver::BaseSolver;

  // An online substructured run takes the reference length of its saved model, so that a
  // mesh of the region alone is nondimensionalized as the model was.
  void Preprocess(IoData &iodata, std::unique_ptr<mfem::Mesh> &smesh,
                  MPI_Comm comm) const override;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_ELECTROSTATIC_SOLVER_HPP
