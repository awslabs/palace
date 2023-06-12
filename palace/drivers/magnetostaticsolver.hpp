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
  MagnetostaticSolver(const IoData &iodata, bool root, int size = 0, int num_thread = 0,
                      const char *git_tag = nullptr)
    : BaseSolver(iodata, root, size, num_thread, git_tag)
  {
  }

  void Solve(std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
             Timer &timer) const override;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_MAGNETOSTATIC_SOLVER_HPP
