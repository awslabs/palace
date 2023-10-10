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
class ErrorIndicator;
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
  ErrorIndicator Postprocess(CurlCurlOperator &curlcurlop, PostOperator &postop,
                             const std::vector<Vector> &A) const;

  void PostprocessTerminals(const SurfaceCurrentOperator &surf_j_op,
                            const mfem::DenseMatrix &M, const mfem::DenseMatrix &Minv,
                            const mfem::DenseMatrix &Mm) const;

public:
  using BaseSolver::BaseSolver;

  ErrorIndicator
  Solve(const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh) const override;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_MAGNETOSTATIC_SOLVER_HPP
