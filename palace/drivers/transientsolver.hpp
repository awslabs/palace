// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_TRANSIENT_SOLVER_HPP
#define PALACE_DRIVERS_TRANSIENT_SOLVER_HPP

#include <functional>
#include <memory>
#include <vector>
#include "drivers/basesolver.hpp"

namespace palace
{

class ErrorIndicator;
class LumpedPortOperator;
class Mesh;
class PostOperator;
class SurfaceCurrentOperator;
class SpaceOperator;

//
// Driver class for time-dependent driven terminal simulations.
//
class TransientSolver : public BaseSolver
{
private:
  std::function<double(double)> GetTimeExcitation(bool dot) const;

  int GetNumSteps(double start, double end, double delta) const;

  class CurrentsPostPrinter
  {
    bool root_ = false;
    bool do_measurement_ = false;
    TableWithCSVFile surface_I = {};

  public:
    CurrentsPostPrinter() = default;
    CurrentsPostPrinter(bool do_measurement, bool root, const fs::path &post_dir,
                        const SurfaceCurrentOperator &surf_j_op, int n_expected_rows);
    void AddMeasurement(double t, double J_coef, const SurfaceCurrentOperator &surf_j_op,
                        const IoData &iodata);
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
    void AddMeasurement(double t, double J_coef, const PostOperator &post_op,
                        const LumpedPortOperator &lumped_port_op, const IoData &iodata);
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

    ErrorIndicatorPostPrinter error_indicator;

    PostprocessPrintResults(bool is_mpi_root, const fs::path &post_dir,
                            const PostOperator &post_op, const SpaceOperator &space_op,
                            int n_expected_rows, int delta_post);
    void PostprocessStep(const IoData &iodata, const PostOperator &post_op,
                         const SpaceOperator &space_op, int step, double t, double J_coef);
    void PostprocessFinal(const PostOperator &post_op, const ErrorIndicator &indicator);
  };

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

public:
  using BaseSolver::BaseSolver;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_TRANSIENT_SOLVER_HPP
