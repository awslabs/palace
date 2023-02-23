// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MAGNETOSTATIC_SOLVER_HPP
#define PALACE_MAGNETOSTATIC_SOLVER_HPP

#include <memory>
#include <vector>
#include "drivers/basesolver.hpp"

namespace mfem
{

class DenseMatrix;
class ParMesh;
class Vector;

}  // namespace mfem

namespace palace
{

class CurlCurlOperator;
class IoData;
class PostOperator;
class SurfaceCurrentOperator;
class Timer;

//
// Driver class for magnetostatic simulations.
//
class MagnetostaticSolver : public BaseSolver
{
private:
  void Postprocess(const std::string &post_dir, CurlCurlOperator &curlcurlop,
                   PostOperator &postop, const std::vector<mfem::Vector> &A,
                   Timer &timer) const;

  void PostprocessTerminals(const std::string &post_dir,
                            const SurfaceCurrentOperator &surf_j_op,
                            const mfem::DenseMatrix &M, const mfem::DenseMatrix &Minv,
                            const mfem::DenseMatrix &Mm) const;

public:
  using BaseSolver::BaseSolver;

  BaseSolver::SolveOutput Solve(std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
                                Timer &timer, int iter) const override;
};

}  // namespace palace

#endif  // PALACE_MAGNETOSTATIC_SOLVER_HPP
