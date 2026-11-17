// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_PERFMETRICS_HPP
#define PALACE_UTILS_PERFMETRICS_HPP

#include <vector>
#include <mpi.h>

namespace palace::perf
{

// Statistics of a value over the processes which contribute to it.
struct Stats
{
  double min = 0.0, max = 0.0, sum = 0.0;
  int max_rank = -1;  // Lowest rank holding the maximum
  int ranks = 0;      // Number of contributing processes (0: set nowhere)
  double Avg() const { return (ranks > 0) ? sum / ranks : 0.0; }
};

// Reduce each value over the processes of comm; value i contributes on a process if
// set[i]. Collective: three MPI_Allreduce.
std::vector<Stats> ReduceStats(const std::vector<double> &value,
                               const std::vector<char> &set, MPI_Comm comm);

}  // namespace palace::perf

#endif  // PALACE_UTILS_PERFMETRICS_HPP
