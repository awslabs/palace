// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "schursolver.hpp"

#include <algorithm>
#include <complex>
#include <utility>
#include "utils/communication.hpp"

namespace palace
{

namespace
{

HYPRE_BigInt FirstRow(MPI_Comm comm, int n_loc)
{
  HYPRE_BigInt n = n_loc, first = 0;
  MPI_Exscan(&n, &first, 1, HYPRE_MPI_BIG_INT, MPI_SUM, comm);
  return Mpi::Root(comm) ? 0 : first;
}

}  // namespace

SchurRowGather::SchurRowGather(MPI_Comm comm, HYPRE_BigInt row0, int n_loc,
                               const std::vector<HYPRE_BigInt> &schur_vars)
  : comm(comm), n_schur(static_cast<int>(schur_vars.size()))
{
  for (int p = 0; p < n_schur; p++)
  {
    if (schur_vars[p] >= row0 && schur_vars[p] < row0 + n_loc)
    {
      row.push_back(static_cast<int>(schur_vars[p] - row0));
      pos.push_back(p);
    }
  }
  const int m = static_cast<int>(row.size()), nranks = Mpi::Size(comm);
  const bool root = Mpi::Root(comm);
  pos_cnt.assign(root ? nranks : 0, 0);
  pos_dsp.assign(root ? nranks : 0, 0);
  MPI_Gather(&m, 1, MPI_INT, pos_cnt.data(), 1, MPI_INT, 0, comm);
  for (int r = 1; r < static_cast<int>(pos_cnt.size()); r++)
  {
    pos_dsp[r] = pos_dsp[r - 1] + pos_cnt[r - 1];
  }
  all_pos.resize(root ? n_schur : 0);
  MPI_Gatherv(pos.data(), m, MPI_INT, all_pos.data(), pos_cnt.data(), pos_dsp.data(),
              MPI_INT, 0, comm);
}

void SchurRowGather::Gather(const std::vector<double> &loc, int nb, double *g) const
{
  const bool root = Mpi::Root(comm);
  const int m = static_cast<int>(row.size());
  std::vector<int> cnt(pos_cnt.size()), dsp(pos_dsp.size());
  for (std::size_t q = 0; q < cnt.size(); q++)
  {
    cnt[q] = nb * pos_cnt[q];
    dsp[q] = nb * pos_dsp[q];
  }
  std::vector<double> all(root ? static_cast<std::size_t>(n_schur) * nb : 0);
  MPI_Gatherv(loc.data(), nb * m, MPI_DOUBLE, all.data(), cnt.data(), dsp.data(),
              MPI_DOUBLE, 0, comm);
  for (std::size_t q = 0; q < cnt.size(); q++)
  {
    for (int k = 0; k < nb; k++)
    {
      for (int r = 0; r < pos_cnt[q]; r++)
      {
        g[all_pos[pos_dsp[q] + r] + static_cast<std::size_t>(n_schur) * k] =
            all[dsp[q] + static_cast<std::size_t>(k) * pos_cnt[q] + r];
      }
    }
  }
}

SchurRowGather::SchurRowGather(MPI_Comm comm, int n_loc,
                               const std::vector<HYPRE_BigInt> &schur_vars)
  : SchurRowGather(comm, FirstRow(comm, n_loc), n_loc, schur_vars)
{
}

template <typename T>
std::vector<T> SchurRowGather::GatherSchur(const std::vector<T> &rows) const
{
  const int m = static_cast<int>(row.size());
  std::vector<T> full;
  if (n_schur == 0)
  {
    return full;
  }
  const RowTypeT<T> type(n_schur);
  if (!Mpi::Root(comm))
  {
    MPI_Send(rows.data(), m, type.type, 0, 0, comm);
    return full;
  }
  full.resize(static_cast<std::size_t>(n_schur) * n_schur);
  std::vector<T> buf;
  for (std::size_t q = 0; q < pos_cnt.size(); q++)
  {
    const T *x = rows.data();
    if (q > 0)
    {
      buf.resize(static_cast<std::size_t>(pos_cnt[q]) * n_schur);
      MPI_Recv(buf.data(), pos_cnt[q], type.type, static_cast<int>(q), 0, comm,
               MPI_STATUS_IGNORE);
      x = buf.data();
    }
    for (int r = 0; r < pos_cnt[q]; r++)
    {
      const std::size_t p = all_pos[pos_dsp[q] + r];
      for (int c = 0; c < n_schur; c++)
      {
        full[p + static_cast<std::size_t>(n_schur) * c] =
            x[static_cast<std::size_t>(r) * n_schur + c];
      }
    }
  }
  return full;
}

template <typename T>
std::vector<T> SchurRowGather::ScatterSchur(const std::vector<T> &full) const
{
  const int m = static_cast<int>(row.size());
  std::vector<T> rows(static_cast<std::size_t>(m) * n_schur);
  if (n_schur == 0)
  {
    return rows;
  }
  const RowTypeT<T> type(n_schur);
  if (!Mpi::Root(comm))
  {
    MPI_Recv(rows.data(), m, type.type, 0, 0, comm, MPI_STATUS_IGNORE);
    return rows;
  }
  std::vector<T> buf;
  for (std::size_t q = 0; q < pos_cnt.size(); q++)
  {
    if (q > 0)
    {
      buf.resize(static_cast<std::size_t>(pos_cnt[q]) * n_schur);
    }
    T *x = (q > 0) ? buf.data() : rows.data();
    for (int r = 0; r < pos_cnt[q]; r++)
    {
      const std::size_t p = all_pos[pos_dsp[q] + r];
      for (int c = 0; c < n_schur; c++)
      {
        x[static_cast<std::size_t>(r) * n_schur + c] =
            full[p + static_cast<std::size_t>(n_schur) * c];
      }
    }
    if (q > 0)
    {
      MPI_Send(x, pos_cnt[q], type.type, static_cast<int>(q), 0, comm);
    }
  }
  return rows;
}

template std::vector<double> SchurRowGather::GatherSchur(const std::vector<double> &) const;
template std::vector<std::complex<double>>
SchurRowGather::GatherSchur(const std::vector<std::complex<double>> &) const;
template std::vector<double>
SchurRowGather::ScatterSchur(const std::vector<double> &) const;
template std::vector<std::complex<double>>
SchurRowGather::ScatterSchur(const std::vector<std::complex<double>> &) const;

}  // namespace palace
