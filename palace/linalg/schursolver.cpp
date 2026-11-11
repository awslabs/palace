// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "schursolver.hpp"

#include <algorithm>
#include <utility>
#include "utils/communication.hpp"

namespace palace
{

BackSolveSchurSolver::BackSolveSchurSolver(const mfem::HypreParMatrix &A,
                                           const std::vector<HYPRE_BigInt> &schur_vars,
                                           InternalSolver solve, int block)
  : A(A), solve(std::move(solve)), n_schur(static_cast<int>(schur_vars.size())),
    block(block)
{
  MPI_Comm comm = A.GetComm();
  const HYPRE_BigInt row0 = A.GetRowStarts()[0];
  for (int p = 0; p < n_schur; p++)
  {
    if (schur_vars[p] >= row0 && schur_vars[p] < row0 + A.Height())
    {
      schur_row.push_back(static_cast<int>(schur_vars[p] - row0));
      schur_pos.push_back(p);
    }
  }
  const int m = static_cast<int>(schur_row.size()), nranks = Mpi::Size(comm);
  const bool root = Mpi::Root(comm);
  pos_cnt.assign(root ? nranks : 0, 0);
  pos_dsp.assign(root ? nranks : 0, 0);
  MPI_Gather(&m, 1, MPI_INT, pos_cnt.data(), 1, MPI_INT, 0, comm);
  for (int r = 1; r < static_cast<int>(pos_cnt.size()); r++)
  {
    pos_dsp[r] = pos_dsp[r - 1] + pos_cnt[r - 1];
  }
  all_pos.resize(root ? n_schur : 0);
  MPI_Gatherv(schur_pos.data(), m, MPI_INT, all_pos.data(), pos_cnt.data(), pos_dsp.data(),
              MPI_INT, 0, comm);
}

std::vector<double> BackSolveSchurSolver::ReleaseSchur()
{
  MPI_Comm comm = A.GetComm();
  const bool root = Mpi::Root(comm);
  const int n_loc = A.Height(), m = static_cast<int>(schur_row.size());
  // The local row of each Schur variable on this rank, or -1.
  std::vector<int> local(n_schur, -1);
  for (int r = 0; r < m; r++)
  {
    local[schur_pos[r]] = schur_row[r];
  }
  std::vector<double> S(root ? static_cast<std::size_t>(n_schur) * n_schur : 0);
  std::vector<mfem::Vector> t(std::min(block, n_schur), mfem::Vector(n_loc));
  mfem::Vector z(n_loc);
  // The entries of the local Schur rows of the columns of a batch, gathered on rank 0.
  std::vector<double> loc(static_cast<std::size_t>(m) * block),
      all(root ? static_cast<std::size_t>(n_schur) * block : 0);
  std::vector<int> cnt(pos_cnt.size()), dsp(pos_dsp.size());
  for (int c0 = 0; c0 < n_schur; c0 += block)
  {
    const int nb = std::min(block, n_schur - c0);
    std::vector<const mfem::Vector *> X(nb);
    std::vector<mfem::Vector *> Y(nb);
    for (int k = 0; k < nb; k++)
    {
      // Column c of A: A_IS e_c on the internal variables, A_SS e_c on the Schur ones.
      z = 0.0;
      if (local[c0 + k] >= 0)
      {
        z.HostReadWrite()[local[c0 + k]] = 1.0;
      }
      A.Mult(z, t[k]);
      const double *tk = t[k].HostRead();
      for (int r = 0; r < m; r++)
      {
        loc[static_cast<std::size_t>(k) * m + r] = tk[schur_row[r]];
      }
      X[k] = &t[k];
      Y[k] = &t[k];
    }
    solve(X, Y);
    for (int k = 0; k < nb; k++)
    {
      A.Mult(t[k], z);  // A_SI A_II^-1 A_IS e_c on the Schur variables
      const double *zz = z.HostRead();
      for (int r = 0; r < m; r++)
      {
        loc[static_cast<std::size_t>(k) * m + r] -= zz[schur_row[r]];
      }
    }
    for (std::size_t q = 0; q < cnt.size(); q++)
    {
      cnt[q] = nb * pos_cnt[q];
      dsp[q] = nb * pos_dsp[q];
    }
    MPI_Gatherv(loc.data(), nb * m, MPI_DOUBLE, all.data(), cnt.data(), dsp.data(),
                MPI_DOUBLE, 0, comm);
    for (std::size_t q = 0; q < cnt.size(); q++)
    {
      for (int k = 0; k < nb; k++)
      {
        for (int r = 0; r < pos_cnt[q]; r++)
        {
          S[all_pos[pos_dsp[q] + r] + static_cast<std::size_t>(n_schur) * (c0 + k)] =
              all[dsp[q] + static_cast<std::size_t>(k) * pos_cnt[q] + r];
        }
      }
    }
  }
  return S;
}

void BackSolveSchurSolver::SolveInternal(const std::vector<const mfem::Vector *> &X,
                                         const std::vector<mfem::Vector *> &Y)
{
  solve(X, Y);
}

}  // namespace palace
