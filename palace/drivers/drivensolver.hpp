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
class PortExcitationHelper;

class PostprocessCurrentsHelper
{
  bool root_ = false;
  bool do_measurement_ = false;
  TableWithCSVFile surface_I = {};

public:
  PostprocessCurrentsHelper() = default;
  PostprocessCurrentsHelper(bool do_measurement, bool root, const std::string &post_dir,
                            const SurfaceCurrentOperator &surf_j_op, int n_expected_rows);
  void AddMeasurement(double omega, const SurfaceCurrentOperator &surf_j_op,
                      const IoData &iodata);
};

class PostprocessPortsHelper
{
  bool root_ = false;
  bool do_measurement_ = false;
  TableWithCSVFile port_V = {};
  TableWithCSVFile port_I = {};

public:
  PostprocessPortsHelper() = default;
  PostprocessPortsHelper(bool do_measurement, bool root, const std::string &post_dir,
                         const LumpedPortOperator &lumped_port_op, int n_expected_rows);
  void AddMeasurement(double omega, const PostOperator &post_op,
                      const LumpedPortOperator &lumped_port_op, const IoData &iodata);
};

class PostprocessSParametersHelper
{
  // Postprocess S-parameters. This computes a column of the S matrix corresponding to the
  // excited port index specified in the configuration file, storing |S_ij| and arg
  // (S_ij) in dB and degrees, respectively. S-parameter output is only available for a
  // single lumped or wave port excitation.

  bool root_ = false;
  bool do_measurement_ = false;
  TableWithCSVFile port_S = {};

public:
  PostprocessSParametersHelper() = default;
  PostprocessSParametersHelper(bool do_measurement, bool root, const std::string &post_dir,
                               const PortExcitationHelper &port_excitations,
                               const LumpedPortOperator &lumped_port_op,
                               const WavePortOperator &wave_port_op, int n_expected_rows);
  void AddMeasurement(double omega, int excitation_idx, const PostOperator &post_op,
                      const PortExcitationHelper &port_excitations,
                      const LumpedPortOperator &lumped_port_op,
                      const WavePortOperator &wave_port_op, const IoData &iodata);
};

//
// Driver class for driven terminal simulations.
//
class DrivenSolver : public BaseSolver
{
private:
  int GetNumSteps(double start, double end, double delta) const;

  struct PostprocessResults
  {
    PostprocessCurrentsHelper currents;
    PostprocessPortsHelper ports;
    PostprocessSParametersHelper s_parameters;

    PostprocessResults(bool is_mpi_root, const std::string &post_dir,
                       const SpaceOperator &space_op,
                       const PortExcitationHelper &port_excitations, int n_expected_rows);
    void PostprocessStep(const IoData &iodata, const PostOperator &post_op,
                         const SpaceOperator &space_op,
                         const PortExcitationHelper &port_excitations, int step,
                         double omega, int excitation_idx, double E_elec, double E_mag,
                         const ErrorIndicator *indicator);
  };

  ErrorIndicator SweepUniform(SpaceOperator &space_op, PostOperator &post_op,
                              PostprocessResults &post_results, int n_step, int step0,
                              double omega0, double delta_omega) const;

  ErrorIndicator SweepAdaptive(SpaceOperator &space_op, PostOperator &post_op,
                               PostprocessResults &post_results, int n_step, int step0,
                               double omega0, double delta_omega) const;

  void Postprocess(const PostOperator &post_op, const LumpedPortOperator &lumped_port_op,
                   const WavePortOperator &wave_port_op,
                   const SurfaceCurrentOperator &surf_j_op, int step, double omega,
                   double E_elec, double E_mag, const ErrorIndicator *indicator) const;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

public:
  using BaseSolver::BaseSolver;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_DRIVEN_SOLVER_HPP
