// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_PERFREPORT_HPP
#define PALACE_UTILS_PERFREPORT_HPP

#include <string_view>
#include <vector>
#include <mpi.h>
#include <nlohmann/json_fwd.hpp>
#include "utils/perfmetrics.hpp"

namespace palace::perf
{

inline constexpr int METADATA_VERSION = 2;

// Reduced values of one period (this iteration, or the run so far).
struct PeriodResults
{
  std::vector<Stats> time, count, mem;  // Per timer row (mem: not TOTAL)
  std::vector<Stats> node_mem;          // Per timer row, statistics over nodes
  Stats peak, baseline, unattributed;
  Stats node_peak, node_baseline, node_unattributed;
  std::vector<Stats> metric;  // Per GlobalMetrics() entry
};

struct ReportResults
{
  int iteration = 1;
  PeriodResults iter, cum;

  // Statistics of a metric, or nullptr if it was set nowhere in the period.
  const Stats *Metric(std::string_view path, bool cumulative) const;
};

// Reduce the report of the current iteration. Collective on comm (Mpi::World()).
const ReportResults &Report(MPI_Comm comm);

// Write the last report into palace.json (global root only), replacing what the previous
// report wrote and keeping the keys of other writers.
void WriteJson(nlohmann::json &meta);

// Print the log tables of the last report.
void Print(MPI_Comm comm);

}  // namespace palace::perf

#endif  // PALACE_UTILS_PERFREPORT_HPP
