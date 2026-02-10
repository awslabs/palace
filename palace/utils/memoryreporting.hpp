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

// Aggregated memory statistics across all MPI processes, in bytes.
struct MemoryStats
{
  std::string label;
  long min;
  long max;
  long sum;
  double avg;
};
MemoryStats GetCurrentMemoryStats(MPI_Comm comm);
MemoryStats GetPeakMemoryStats(MPI_Comm comm);

// Print memory usage summary.
void PrintMemoryUsage(MPI_Comm comm, const MemoryStats &stats);

}  // namespace palace::memory_reporting

#endif  // PALACE_UTILS_MEMORY_REPORTING_HPP
