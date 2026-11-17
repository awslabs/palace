// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "perfmetrics.hpp"

#include <limits>
#include <mfem.hpp>
#include "utils/communication.hpp"

namespace palace::perf
{

std::vector<Stats> ReduceStats(const std::vector<double> &value,
                               const std::vector<char> &set, MPI_Comm comm)
{
  MFEM_VERIFY(value.size() == set.size(), "ReduceStats: size mismatch!");
  const int n = static_cast<int>(value.size());
  const int rank = Mpi::Rank(comm);
  constexpr double inf = std::numeric_limits<double>::infinity();
  struct DoubleInt
  {
    double v;
    int r;
  };
  std::vector<double> mn(n), sm(2 * n);
  std::vector<DoubleInt> mx(n);
  for (int i = 0; i < n; i++)
  {
    mn[i] = set[i] ? value[i] : inf;
    mx[i] = {set[i] ? value[i] : -inf, rank};
    sm[i] = set[i] ? value[i] : 0.0;
    sm[n + i] = set[i] ? 1.0 : 0.0;
  }
  // Raw MPI (not Mpi::Global*), so that future tallies of Palace collectives do not
  // measure the report itself.
  MPI_Allreduce(MPI_IN_PLACE, mn.data(), n, MPI_DOUBLE, MPI_MIN, comm);
  MPI_Allreduce(MPI_IN_PLACE, mx.data(), n, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
  MPI_Allreduce(MPI_IN_PLACE, sm.data(), 2 * n, MPI_DOUBLE, MPI_SUM, comm);
  std::vector<Stats> out(n);
  for (int i = 0; i < n; i++)
  {
    const int ranks = static_cast<int>(sm[n + i]);
    if (ranks > 0)
    {
      out[i] = {mn[i], mx[i].v, sm[i], mx[i].r, ranks};
    }
  }
  return out;
}

}  // namespace palace::perf
