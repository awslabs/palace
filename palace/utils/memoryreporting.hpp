// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_MEMORY_REPORTING_HPP
#define PALACE_UTILS_MEMORY_REPORTING_HPP

#include <string>
#include <mfem.hpp>

namespace palace::memory_reporting
{

// Per-process memory usage in bytes.
long GetCurrentMemory();
long GetPeakMemory();

// Aggregated memory statistics, in bytes.
struct MemoryStats
{
  std::string label;
  long min;
  long max;
  long sum;
  double avg;
};

// Per-rank memory statistics (aggregated across all MPI ranks).
MemoryStats GetCurrentMemoryStats(MPI_Comm comm);
MemoryStats GetPeakMemoryStats(MPI_Comm comm);

// Per-node memory statistics (processes sharing physical memory), in bytes.
// Uses MPI_Comm_split_type to identify shared memory groups, sums memory per node,
// and reports statistics across nodes. Assumes shared memory is negligible.
MemoryStats GetCurrentNodeMemoryStats(MPI_Comm comm);
MemoryStats GetPeakNodeMemoryStats(MPI_Comm comm);

// Print memory usage summary.
void PrintMemoryUsage(MPI_Comm comm, const MemoryStats &stats);

}  // namespace palace::memory_reporting

#endif  // PALACE_UTILS_MEMORY_REPORTING_HPP
