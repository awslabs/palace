// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_EIGEN_SOLVER_HPP
#define PALACE_DRIVERS_EIGEN_SOLVER_HPP

#include <complex>
#include <memory>
#include <vector>
#include "drivers/basesolver.hpp"

namespace palace
{

class ErrorIndicator;
class LumpedPortOperator;
class Mesh;
class PostOperator;
class SpaceOperator;

//
// Driver class for eigenmode simulations.
//
class EigenSolver : public BaseSolver
{
private:
  struct EigenPostPrinter
  {
    bool root_ = false;
    TableWithCSVFile eig;

    // Print data to stdout with custom table formatting
    void PrintStdoutHeader();
    void PrintStdoutRow(size_t j);

  public:
    int stdout_int_print_width = 0;

    EigenPostPrinter(const fs::path &post_dir, const SpaceOperator &space_op, int n_post);
    void AddMeasurement(int eigen_print_idx, const PostOperator &post_op, double error_bkwd,
                        double error_abs, const IoData &iodata);
  };

  class PortsPostPrinter
  {
    TableWithCSVFile port_V;
    TableWithCSVFile port_I;

  public:
    PortsPostPrinter(const fs::path &post_dir, const SpaceOperator &space_op,
                     int n_expected_rows);
    void AddMeasurement(int eigen_print_idx, const PostOperator &post_op,
                        const LumpedPortOperator &lumped_port_op, const IoData &iodata);
  };

  class EPRPostPrinter
  {
    TableWithCSVFile port_EPR;
    TableWithCSVFile port_Q;

    std::vector<int> ports_with_L;
    std::vector<int> ports_with_R;
    void AddMeasurementEPR(double eigen_print_idx, const PostOperator &post_op,
                           const LumpedPortOperator &lumped_port_op, const IoData &iodata);
    void AddMeasurementQ(double eigen_print_idx, const PostOperator &post_op,
                         const LumpedPortOperator &lumped_port_op, const IoData &iodata);

  public:
    EPRPostPrinter(const fs::path &post_dir, const SpaceOperator &space_op,
                   int n_expected_rows);

    void AddMeasurement(double eigen_print_idx, const PostOperator &post_op,
                        const LumpedPortOperator &lumped_port_op, const IoData &iodata);
  };

  struct PostprocessPrintResults
  {
    bool write_paraview_fields = false;
    int n_post = 0;

    DomainsPostPrinter domains;
    SurfacesPostPrinter surfaces;
    ProbePostPrinter probes;
    EigenPostPrinter eigen;
    EPRPostPrinter epr;

    ErrorIndicatorPostPrinter error_indicator;

    PostprocessPrintResults(const fs::path &post_dir, const PostOperator &post_op,
                            const SpaceOperator &space_op, int n_post_);
    void PostprocessStep(const IoData &iodata, const PostOperator &post_op,
                         const SpaceOperator &space_op, int step, double error_abs,
                         double error_bkward);
    void PostprocessFinal(const PostOperator &post_op, const ErrorIndicator &indicator);
  };

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

public:
  using BaseSolver::BaseSolver;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_EIGEN_SOLVER_HPP
