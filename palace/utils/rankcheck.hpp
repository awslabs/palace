// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_RANKCHECK_HPP
#define PALACE_UTILS_RANKCHECK_HPP

// CI diagnostic only: report when a result that every rank computes redundantly differs
// across the ranks of a communicator. All ranks must call it at the same point.

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <vector>
#include <mpi.h>

namespace palace::rankcheck
{

inline void Check(MPI_Comm comm, const char *what, const double *data, int n)
{
  if (comm == MPI_COMM_NULL)
  {
    return;
  }
  int size = 1, rank = 0;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);
  if (size == 1)
  {
    return;
  }
  int nminmax[2] = {-n, n};
  MPI_Allreduce(MPI_IN_PLACE, nminmax, 2, MPI_INT, MPI_MAX, comm);
  if (-nminmax[0] != nminmax[1])
  {
    if (rank == 0)
    {
      std::fprintf(stderr, "RANKDIFF %s: sizes differ (min %d, max %d)\n", what,
                   -nminmax[0], nminmax[1]);
    }
    return;
  }
  std::vector<double> root(data, data + n);
  MPI_Bcast(root.data(), n, MPI_DOUBLE, 0, comm);
  double stats[2] = {0.0, 0.0};  // Number of differing entries, max relative difference
  for (int i = 0; i < n; i++)
  {
    if (root[i] != data[i])
    {
      stats[0] += 1.0;
      stats[1] = std::max(stats[1], std::abs(root[i] - data[i]) /
                                        std::max(std::abs(root[i]), 1.0e-300));
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, stats, 2, MPI_DOUBLE, MPI_MAX, comm);
  if (rank == 0)
  {
    std::fprintf(stderr, "RANKCHECK %s n=%d %s ndiff=%.0f maxrel=%.3e\n", what, n,
                 stats[0] > 0.0 ? "DIFF" : "same", stats[0], stats[1]);
  }
}

inline void Check(MPI_Comm comm, const char *what, const std::complex<double> *data, int n)
{
  Check(comm, what, reinterpret_cast<const double *>(data), 2 * n);
}

inline void Check(MPI_Comm comm, const char *what, const std::vector<int> &data)
{
  std::vector<double> d(data.begin(), data.end());
  Check(comm, what, d.data(), static_cast<int>(d.size()));
}

}  // namespace palace::rankcheck

#endif  // PALACE_UTILS_RANKCHECK_HPP
