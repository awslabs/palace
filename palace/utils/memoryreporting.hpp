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

// Return the free memory of the heap to the operating system, where supported (glibc
// malloc_trim), after large data structures have been released: freed memory scattered
// among small allocations otherwise stays resident.
void TrimHeap();

// Reset the peak RSS to the current RSS, where supported (Linux /proc/self/clear_refs).
// Returns whether the reset succeeded.
bool ResetPeakMemory();

// Memory available to the process on its node, in bytes: MemAvailable of /proc/meminfo,
// or less if a memory control group (version 1 or 2) of the process or of one of its
// ancestors has less headroom (its limit minus its usage without the inactive file
// cache). Returns -1 if unknown: off Linux, or if a file it needs cannot be read.
long GetAvailableMemory();

// Make GetAvailableMemory return bytes instead (for tests; negative restores it).
void SetAvailableMemoryForTesting(long bytes);

// Shared-memory (node) communicator and node-leader communicator (MPI_COMM_NULL on
// non-leaders) for comm, created once and cached.
struct NodeComms
{
  MPI_Comm node = MPI_COMM_NULL, leaders = MPI_COMM_NULL;
  int node_rank = 0;
};
const NodeComms &GetNodeComms(MPI_Comm comm);

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

// Compute per-node statistics for an arbitrary per-rank value (bytes). Sums across
// processes sharing a node and reports min/max/avg across nodes.
MemoryStats ComputeNodeMemoryStats(std::string label, long local_value, MPI_Comm comm);

// Format a byte count as a human-readable string (e.g., "1.5G", "228.0M").
std::string FormatBytes(double bytes);

// Print memory usage summary.
void PrintMemoryUsage(MPI_Comm comm, const MemoryStats &stats);

}  // namespace palace::memory_reporting

#endif  // PALACE_UTILS_MEMORY_REPORTING_HPP
