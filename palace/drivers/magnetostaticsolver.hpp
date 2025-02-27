// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_MAGNETOSTATIC_SOLVER_HPP
#define PALACE_DRIVERS_MAGNETOSTATIC_SOLVER_HPP

#include <memory>
#include <vector>
#include "drivers/basesolver.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class ErrorIndicator;
class Mesh;
class PostOperator;
class SurfaceCurrentOperator;

//
// Driver class for magnetostatic simulations.
//
class MagnetostaticSolver : public BaseSolver
{
private:
  struct PostprocessPrintResults
  {
    bool write_paraview_fields = false;
    int n_post = 0;

    DomainsPostPrinter domains;
    SurfacesPostPrinter surfaces;
    ProbePostPrinter probes;

    ErrorIndicatorPostPrinter error_indicator;

    PostprocessPrintResults(bool is_mpi_root, const std::string &post_dir,
                            const PostOperator &post_op, int n_post_);
    void PostprocessStep(const IoData &iodata, const PostOperator &post_op, int step,
                         int idx, double E_mag);
    void PostprocessFinal(const PostOperator &post_op, const ErrorIndicator &indicator);
  };

  void PostprocessTerminals(PostOperator &post_op, const SurfaceCurrentOperator &surf_j_op,
                            const std::vector<Vector> &A,
                            const std::vector<double> &I_inc) const;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

public:
  using BaseSolver::BaseSolver;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_MAGNETOSTATIC_SOLVER_HPP
