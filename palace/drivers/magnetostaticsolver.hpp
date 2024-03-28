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

}  // namespace mfem

namespace palace
{

class CurlCurlOperator;
class ErrorIndicator;
class IoData;
class Mesh;
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
                   const std::vector<Vector> &A, const ErrorIndicator &indicator) const;

  void PostprocessTerminals(const SurfaceCurrentOperator &surf_j_op,
                            const mfem::DenseMatrix &M, const mfem::DenseMatrix &Minv,
                            const mfem::DenseMatrix &Mm) const;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

public:
  using BaseSolver::BaseSolver;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_MAGNETOSTATIC_SOLVER_HPP
