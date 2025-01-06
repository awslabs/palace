// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_BASE_SOLVER_HPP
#define PALACE_DRIVERS_BASE_SOLVER_HPP

#include <memory>
#include <string>
#include <vector>
#include <fmt/os.h>
#include "fem/errorindicator.hpp"
#include "utils/filesystem.hpp"
#include "utils/tablecsv.hpp"

namespace palace
{

class DomainPostOperator;
class FiniteElementSpaceHierarchy;
class IoData;
class Mesh;
class PostOperator;
class Timer;

//
// Base driver class for all simulation types.
//
class BaseSolver
{
protected:
  // Reference to configuration file data (not owned).
  const IoData &iodata;

  // Parameters for writing postprocessing outputs.
  fs::path post_dir;

  // Common domain postprocessing for all simulation types.
  class DomainsPostPrinter
  {
    TableWithCSVFile domain_E;
  public:
    DomainsPostPrinter() = default;
    DomainsPostPrinter(const fs::path &post_dir, const PostOperator &post_op, const std::string &idx_col_name, int n_expected_rows);
    void AddMeasurement(double idx_value_dimensionful, const PostOperator &post_op, const IoData &iodata);
  };

  // Common surface postprocessing for all simulation types.
  class SurfacesPostPrinter
  {
    TableWithCSVFile surface_F;
    TableWithCSVFile surface_Q;
    void AddMeasurementFlux(double idx_value_dimensionful, const PostOperator &post_op,
                            const IoData &iodata);
    void AddMeasurementEps(double idx_value_dimensionful, const PostOperator &post_op,
                           const IoData &iodata);
  public:
    SurfacesPostPrinter() = default;
    SurfacesPostPrinter(const fs::path &post_dir, const PostOperator &post_op, const std::string &idx_col_name, int n_expected_rows);
    void AddMeasurement(double idx_value_dimensionful, const PostOperator &post_op, const IoData &iodata);
  };

  // Common probe postprocessing for all simulation types.
  class ProbePostPrinter
  {
    TableWithCSVFile probe_E;
    TableWithCSVFile probe_B;

    void AddMeasurementE(double idx_value_dimensionful, const PostOperator &post_op,
                         const IoData &iodata);
    void AddMeasurementB(double idx_value_dimensionful, const PostOperator &post_op,
                         const IoData &iodata);
  public:
    ProbePostPrinter() = default;
    ProbePostPrinter(const fs::path &post_dir, const PostOperator &post_op, const std::string &idx_col_name, int n_expected_rows);

    void AddMeasurement(double idx_value_dimensionful, const PostOperator &post_op, const IoData &iodata);
  };

  // Common error indicator postprocessing for all simulation types. //
  // This is trivial since data is only added at the end of the solve, rather after each
  // step (time / frequency / eigenvector).
  class ErrorIndicatorPostPrinter
  {
    bool root_ = false;
    bool do_measurement_ = false;
    TableWithCSVFile error_indicator;

  public:
    ErrorIndicatorPostPrinter() = default;
    ErrorIndicatorPostPrinter(const fs::path &post_dir);

    void PrintIndicatorStatistics(const PostOperator &post_op, const ErrorIndicator::SummaryStatistics &indicator_stats);
  };

  // Performs a solve using the mesh sequence, then reports error indicators and the number
  // of global true dofs.
  virtual std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const = 0;

public:
  BaseSolver(const IoData &iodata, bool root, int size = 0, int num_thread = 0,
             const char *git_tag = nullptr);
  virtual ~BaseSolver() = default;

  // Performs adaptive mesh refinement using the solve-estimate-mark-refine paradigm.
  // Dispatches to the Solve method for the driver specific calculations.
  void SolveEstimateMarkRefine(std::vector<std::unique_ptr<Mesh>> &mesh) const;

  // These methods write different simulation metadata to a JSON file in post_dir.
  void SaveMetadata(const FiniteElementSpaceHierarchy &fespaces) const;
  template <typename SolverType>
  void SaveMetadata(const SolverType &ksp) const;
  void SaveMetadata(const Timer &timer) const;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_BASE_SOLVER_HPP
