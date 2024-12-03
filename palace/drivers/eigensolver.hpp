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
    bool do_measurement_ = false;
    TableWithCSVFile eig = {};

    // Print data to stdout with custom table formatting
    void PrintStdoutHeader();
    void PrintStdoutRow(size_t j);

  public:
    int stdout_int_print_width = 0;

    EigenPostPrinter() = default;
    EigenPostPrinter(bool do_measurement, bool root, const fs::path &post_dir, int n_post);
    void AddMeasurement(int eigen_print_idx, std::complex<double> omega, double error_bkwd,
                        double error_abs, const IoData &iodata);
  };

  class PortsPostPrinter
  {
    bool root_ = false;
    bool do_measurement_ = false;
    TableWithCSVFile port_V = {};
    TableWithCSVFile port_I = {};

  public:
    PortsPostPrinter() = default;
    PortsPostPrinter(bool do_measurement, bool root, const fs::path &post_dir,
                     const LumpedPortOperator &lumped_port_op, int n_expected_rows);
    void AddMeasurement(int eigen_print_idx, const PostOperator &post_op,
                        const LumpedPortOperator &lumped_port_op, const IoData &iodata);
  };

  // Common domain postprocessing for all simulation types.
  class EPRPostPrinter
  {
    bool root_ = false;
    bool do_measurement_EPR_ = false;
    bool do_measurement_Q_ = false;
    TableWithCSVFile port_EPR = {};
    TableWithCSVFile port_Q = {};

    std::vector<int> ports_with_L;
    std::vector<int> ports_with_R;

  public:
    EPRPostPrinter() = default;
    EPRPostPrinter(bool do_measurement, bool root, const fs::path &post_dir,
                   const LumpedPortOperator &lumped_port_op, int n_expected_rows);

    void AddMeasurementEPR(double eigen_print_idx, const PostOperator &post_op,
                           const LumpedPortOperator &lumped_port_op, double E_m,
                           const IoData &iodata);
    void AddMeasurementQ(double eigen_print_idx, const PostOperator &post_op,
                         const LumpedPortOperator &lumped_port_op,
                         std::complex<double> omega, double E_m, const IoData &iodata);

    void AddMeasurement(double eigen_print_idx, const PostOperator &post_op,
                        const LumpedPortOperator &lumped_port_op,
                        std::complex<double> omega, double E_m, const IoData &iodata);
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

    PostprocessPrintResults(bool is_mpi_root, const fs::path &post_dir,
                            const PostOperator &post_op, const SpaceOperator &space_op,
                            int n_post_);
    void PostprocessStep(const IoData &iodata, const PostOperator &post_op,
                         const SpaceOperator &space_op, int step,
                         std::complex<double> omega, double E_elec, double E_mag,
                         double error_abs, double error_bkward);
    void PostprocessFinal(const PostOperator &post_op, const ErrorIndicator &indicator);
  };

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

public:
  using BaseSolver::BaseSolver;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_EIGEN_SOLVER_HPP
