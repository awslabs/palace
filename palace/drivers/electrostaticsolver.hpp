// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_ELECTROSTATIC_SOLVER_HPP
#define PALACE_DRIVERS_ELECTROSTATIC_SOLVER_HPP

#include <map>
#include <memory>
#include <vector>
#include "drivers/basesolver.hpp"
#include "linalg/vector.hpp"

namespace mfem
{

template <typename T>
class Array;

}  // namespace mfem

namespace palace
{

class ErrorIndicator;
class Mesh;
class PostOperator;

//
// Driver class for electrostatic simulations.
//
class ElectrostaticSolver : public BaseSolver
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

    PostprocessPrintResults(bool is_mpi_root, const fs::path &post_dir,
                            const PostOperator &post_op, int n_post_);
    void PostprocessStep(const IoData &iodata, const PostOperator &post_op, int step,
                         int idx);
    void PostprocessFinal(const PostOperator &post_op, const ErrorIndicator &indicator);
  };

  void PostprocessTerminals(PostOperator &post_op,
                            const std::map<int, mfem::Array<int>> &terminal_sources,
                            const std::vector<Vector> &V) const;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

public:
  using BaseSolver::BaseSolver;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_ELECTROSTATIC_SOLVER_HPP
