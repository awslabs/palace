// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_DRIVEN_SOLVER_HPP
#define PALACE_DRIVERS_DRIVEN_SOLVER_HPP

#include <memory>
#include <vector>
#include "drivers/basesolver.hpp"
#include "utils/tablecsv.hpp"

namespace palace
{

class ErrorIndicator;
class LumpedPortOperator;
class Mesh;
class PostOperator;
class SpaceOperator;
class SurfaceCurrentOperator;
class WavePortOperator;

//
// Driver class for driven terminal simulations.
//
class DrivenSolver : public BaseSolver
{
private:
  int GetNumSteps(double start, double end, double delta) const;

  // Printers for storing and printing postprocessing mesurements

  class CurrentsPostPrinter
  {
    TableWithCSVFile surface_I;
  public:
    CurrentsPostPrinter() = default;
    CurrentsPostPrinter(const fs::path &post_dir, const SpaceOperator &space_op, int n_expected_rows);
    void AddMeasurement(double freq, const SurfaceCurrentOperator &surf_j_op, const IoData &iodata);
  };

  class PortsPostPrinter
  {
    TableWithCSVFile port_V;
    TableWithCSVFile port_I;
  public:
    PortsPostPrinter() = default;
    PortsPostPrinter(const fs::path &post_dir, const SpaceOperator &space_op, int n_expected_rows);
    void AddMeasurement(double freq, const PostOperator &post_op,
                        const LumpedPortOperator &lumped_port_op, const IoData &iodata);
  };

  class SParametersPostPrinter
  {
    // Postprocess S-parameters. This computes a column of the S matrix corresponding to the
    // excited port index specified in the configuration file, storing |S_ij| and arg
    // (S_ij) in dB and degrees, respectively. S-parameter output is only available for a
    // single lumped or wave port excitation.

    bool root_ = false;
    bool do_measurement_ = false;
    TableWithCSVFile port_S;

    // Currently can't mix lumped and surface ports for s-matrix
    bool src_lumped_port = true;
    int source_idx = -1;

  public:
    SParametersPostPrinter() = default;
    SParametersPostPrinter(const fs::path &post_dir, const SpaceOperator &space_op, int n_expected_rows);
    void AddMeasurement(double freq, const PostOperator &post_op,
                        const LumpedPortOperator &lumped_port_op,
                        const WavePortOperator &wave_port_op, const IoData &iodata);
  };

  struct PostprocessPrintResults
  {
    bool write_paraview_fields = false;
    int delta_post = 0;

    DomainsPostPrinter domains;
    SurfacesPostPrinter surfaces;
    CurrentsPostPrinter currents;
    ProbePostPrinter probes;
    PortsPostPrinter ports;
    SParametersPostPrinter s_parameters;

    ErrorIndicatorPostPrinter error_indicator;

    PostprocessPrintResults(bool is_mpi_root, const fs::path &post_dir,
                            const PostOperator &post_op, const SpaceOperator &space_op,
                            int n_expected_rows, int delta_post);
    void PostprocessStep(const IoData &iodata, const PostOperator &post_op,
                         const SpaceOperator &space_op, int step);
    void PostprocessFinal(const PostOperator &post_op, const ErrorIndicator &indicator);
  };

  ErrorIndicator SweepUniform(SpaceOperator &space_op, PostOperator &post_op,
                              PostprocessPrintResults &post_results, int n_step, int step0,
                              double omega0, double delta_omega) const;

  ErrorIndicator SweepAdaptive(SpaceOperator &space_op, PostOperator &post_op,
                               PostprocessPrintResults &post_results, int n_step, int step0,
                               double omega0, double delta_omega) const;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

public:
  using BaseSolver::BaseSolver;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_DRIVEN_SOLVER_HPP
