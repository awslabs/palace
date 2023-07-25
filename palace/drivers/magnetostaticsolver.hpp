// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_MAGNETOSTATIC_SOLVER_HPP
#define PALACE_DRIVERS_MAGNETOSTATIC_SOLVER_HPP

#include <memory>
#include <vector>
#include "drivers/basesolver.hpp"
#include "linalg/vector.hpp"

namespace mfem
{

class DenseMatrix;
class ParMesh;

}  // namespace mfem

namespace palace
{

class CurlCurlOperator;
struct ErrorIndicators;
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
  void Postprocess(CurlCurlOperator &curlcurlop, PostOperator &postop,
                   const std::vector<Vector> &A, Timer &timer) const;

  void PostprocessTerminals(const SurfaceCurrentOperator &surf_j_op,
                            const mfem::DenseMatrix &M, const mfem::DenseMatrix &Minv,
                            const mfem::DenseMatrix &Mm) const;

public:
  using BaseSolver::BaseSolver;

  ErrorIndicators Solve(const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
                        Timer &timer) const final;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_MAGNETOSTATIC_SOLVER_HPP
