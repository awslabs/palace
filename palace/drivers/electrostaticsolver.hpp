// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_ELECTROSTATIC_SOLVER_HPP
#define PALACE_ELECTROSTATIC_SOLVER_HPP

#include <map>
#include <memory>
#include <vector>
#include "drivers/basesolver.hpp"

namespace mfem
{

template <typename T>
class Array;
class DenseMatrix;
class ParMesh;
class Vector;

}  // namespace mfem

namespace palace
{

class IoData;
class LaplaceOperator;
class PostOperator;
class Timer;

//
// Driver class for electrostatic simulations.
//
class ElectrostaticSolver : public BaseSolver
{
private:
  void Postprocess(LaplaceOperator &laplaceop, PostOperator &postop,
                   const std::vector<mfem::Vector> &V, Timer &timer) const;

  void PostprocessTerminals(const std::map<int, mfem::Array<int>> &terminal_sources,
                            const mfem::DenseMatrix &C, const mfem::DenseMatrix &Cinv,
                            const mfem::DenseMatrix &Cm) const;

public:
  using BaseSolver::BaseSolver;

  BaseSolver::SolveOutput Solve(std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
                                Timer &timer) const override;
};

}  // namespace palace

#endif  // PALACE_ELECTROSTATIC_SOLVER_HPP
