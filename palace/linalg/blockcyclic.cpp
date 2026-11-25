// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "blockcyclic.hpp"

#include <cmath>
#include <complex>
#include <type_traits>
#include <mfem.hpp>
#include "utils/communication.hpp"

extern "C"
{
  void dpotrf_(const char *, const int *, double *, const int *, int *);
  void dpotrs_(const char *, const int *, const int *, const double *, const int *,
               double *, const int *, int *);
  void zgetrf_(const int *, const int *, std::complex<double> *, const int *, int *, int *);
  void zgetrs_(const char *, const int *, const int *, const std::complex<double> *,
               const int *, const int *, std::complex<double> *, const int *, int *);
#if defined(PALACE_WITH_SCALAPACK)
  int Csys2blacs_handle(MPI_Comm);
  void Cfree_blacs_system_handle(int);
  void Cblacs_gridinit(int *, const char *, int, int);
  void Cblacs_gridexit(int);
  void descinit_(int *, const int *, const int *, const int *, const int *, const int *,
                 const int *, const int *, const int *, int *);
  void pdpotrf_(const char *, const int *, double *, const int *, const int *, const int *,
                int *);
  void pdpotrs_(const char *, const int *, const int *, const double *, const int *,
                const int *, const int *, double *, const int *, const int *, const int *,
                int *);
  void pzgetrf_(const int *, const int *, std::complex<double> *, const int *, const int *,
                const int *, int *, int *);
  void pzgetrs_(const char *, const int *, const int *, const std::complex<double> *,
                const int *, const int *, const int *, const int *, std::complex<double> *,
                const int *, const int *, const int *, int *);
#endif
}

namespace palace
{

namespace
{

// The failure of a LAPACK or ScaLAPACK call on any process: an illegal argument (info < 0)
// before a singular, or not positive definite, matrix (info > 0). Collective.
int Failure(int info, MPI_Comm comm)
{
  int codes[2] = {-std::min(info, 0), std::max(info, 0)};
  Mpi::GlobalMax(2, codes, comm);
  return codes[0] ? -codes[0] : codes[1];
}

}  // namespace

int BlockCyclicGrid::LocalCount(int n, int b, int p, int np)
{
  const int blocks = n / b, extra = blocks % np;
  return (blocks / np) * b + ((p < extra) ? b : (p == extra) ? n % b : 0);
}

BlockCyclicGrid::BlockCyclicGrid(MPI_Comm comm, int n, int nproc, int mb)
  : comm(comm), n(n), mb(mb)
{
  const int nranks = Mpi::Size(comm), rank = Mpi::Rank(comm);
#if defined(PALACE_WITH_SCALAPACK)
  if (nproc <= 0)
  {
    const double side = static_cast<double>(n) / kMinLocal;
    nproc = static_cast<int>(std::min<double>(nranks, std::max(1.0, side * side)));
  }
  nproc = std::clamp(nproc, 1, nranks);
#else
  nproc = 1;
#endif
  nprow = static_cast<int>(std::sqrt(static_cast<double>(nproc)));
  npcol = nproc / nprow;
  if (rank < nprow * npcol)
  {
    myrow = rank / npcol;
    mycol = rank % npcol;
    mloc = LocalCount(n, mb, myrow, nprow);
    nloc = LocalCount(n, mb, mycol, npcol);
  }
}

template <typename T>
std::vector<T> BlockCyclicGrid::ScatterFromRoot(const std::vector<T> &full,
                                                bool row_major) const
{
  // Each local part as columns of its leading dimension, so that the counts are columns.
  const auto type = mpi::DataType<T>();
  auto send = [&](const std::vector<T> &part, int lld, int ncol, int to)
  {
    MPI_Datatype column;
    MPI_Type_contiguous(lld, type, &column);
    MPI_Type_commit(&column);
    MPI_Send(part.data(), std::max(1, ncol), column, to, 0, comm);
    MPI_Type_free(&column);
  };
  std::vector<T> local(LocalSize(), T(0.0));
  if (Mpi::Root(comm))
  {
    MFEM_VERIFY(full.size() == static_cast<std::size_t>(n) * n,
                "Wrong size of a matrix to distribute!");
    std::vector<T> block;
    for (int q = nprow * npcol - 1; q >= 0; q--)
    {
      const int pr = q / npcol, pc = q % npcol, mq = LocalCount(n, mb, pr, nprow),
                nq = LocalCount(n, mb, pc, npcol), lld = std::max(1, mq);
      auto &part = (q == 0) ? local : block;
      part.assign(static_cast<std::size_t>(lld) * std::max(1, nq), T(0.0));
      for (int jl = 0; jl < nq; jl++)
      {
        const auto j = static_cast<std::size_t>(GlobalIndex(jl, mb, pc, npcol));
        for (int il = 0; il < mq; il++)
        {
          const auto i = static_cast<std::size_t>(GlobalIndex(il, mb, pr, nprow));
          part[il + static_cast<std::size_t>(lld) * jl] =
              row_major ? full[i * n + j] : full[i + j * n];
        }
      }
      if (q > 0)
      {
        send(block, lld, nq, q);
      }
    }
  }
  else if (InGrid())
  {
    MPI_Datatype column;
    MPI_Type_contiguous(LeadingDim(), type, &column);
    MPI_Type_commit(&column);
    MPI_Recv(local.data(), std::max(1, nloc), column, 0, 0, comm, MPI_STATUS_IGNORE);
    MPI_Type_free(&column);
  }
  return local;
}

void BlockCyclicGrid::AddFromRows(const std::vector<double> &rows, int row0, int nrows,
                                  std::vector<double> &local, int chunk) const
{
  // Each process sends the entries of its rows to the processes of the grid which hold
  // them, in chunks of columns. The sender and the receiver of the entries of a pair of
  // processes traverse them by rows, then columns, so that only their values are sent.
  const int nranks = Mpi::Size(comm);
  std::vector<int> off(nranks), cnt(nranks);
  MPI_Allgather(&row0, 1, MPI_INT, off.data(), 1, MPI_INT, comm);
  MPI_Allgather(&nrows, 1, MPI_INT, cnt.data(), 1, MPI_INT, comm);
  // The local rows of the grid in the rows of each process, [l0, l1).
  std::vector<int> l0(nranks, 0), l1(nranks, 0);
  for (int s = 0, l = 0; s < nranks && InGrid(); s++)
  {
    while (l < mloc && GlobalRow(l) < off[s])
    {
      l++;
    }
    l0[s] = l;
    while (l < mloc && GlobalRow(l) < off[s] + cnt[s])
    {
      l++;
    }
    l1[s] = l;
  }
  const int most =
      std::max({1, *std::ranges::max_element(cnt), LocalCount(n, mb, 0, nprow)});
  const int cols = std::max(1, chunk / most);
  const std::size_t lld = LeadingDim();
  std::vector<int> scnt(nranks), sdsp(nranks), rcnt(nranks), rdsp(nranks), pos;
  std::vector<double> send, recv;
  for (int c0 = 0; c0 < n; c0 += cols)
  {
    const int c1 = std::min(n, c0 + cols);
    std::ranges::fill(scnt, 0);
    for (int r = 0; r < nrows; r++)
    {
      for (int c = c0; c < c1; c++)
      {
        scnt[Owner(row0 + r, c)]++;
      }
    }
    // The local columns of the grid in [c0, c1), [m0, m1).
    int m0 = 0;
    while (InGrid() && m0 < nloc && GlobalCol(m0) < c0)
    {
      m0++;
    }
    int m1 = m0;
    while (InGrid() && m1 < nloc && GlobalCol(m1) < c1)
    {
      m1++;
    }
    for (int s = 0; s < nranks; s++)
    {
      rcnt[s] = (l1[s] - l0[s]) * (m1 - m0);
      sdsp[s] = (s == 0) ? 0 : sdsp[s - 1] + scnt[s - 1];
      rdsp[s] = (s == 0) ? 0 : rdsp[s - 1] + rcnt[s - 1];
    }
    send.resize(static_cast<std::size_t>(sdsp.back()) + scnt.back());
    recv.resize(static_cast<std::size_t>(rdsp.back()) + rcnt.back());
    pos = sdsp;
    for (int r = 0; r < nrows; r++)
    {
      const double *row = &rows[static_cast<std::size_t>(r) * n];
      for (int c = c0; c < c1; c++)
      {
        send[pos[Owner(row0 + r, c)]++] = row[c];
      }
    }
    MPI_Alltoallv(send.data(), scnt.data(), sdsp.data(), MPI_DOUBLE, recv.data(),
                  rcnt.data(), rdsp.data(), MPI_DOUBLE, comm);
    for (int s = 0; s < nranks; s++)
    {
      const double *v = recv.data() + rdsp[s];
      for (int l = l0[s]; l < l1[s]; l++)
      {
        for (int m = m0; m < m1; m++)
        {
          local[l + lld * m] += *v++;
        }
      }
    }
  }
}

template <typename T>
BlockCyclicFactorization<T>::BlockCyclicFactorization(const BlockCyclicGrid &grid)
  : grid(grid)
{
#if defined(PALACE_WITH_SCALAPACK)
  system = Csys2blacs_handle(grid.GetComm());
  ctxt = system;
  Cblacs_gridinit(&ctxt, "Row", grid.NumProcessRows(), grid.NumProcessCols());
#endif
}

template <typename T>
BlockCyclicFactorization<T>::~BlockCyclicFactorization()
{
#if defined(PALACE_WITH_SCALAPACK)
  if (grid.InGrid())
  {
    Cblacs_gridexit(ctxt);
  }
  Cfree_blacs_system_handle(system);
#endif
}

template <typename T>
void BlockCyclicFactorization<T>::Factor(std::vector<T> &A)
{
  constexpr bool real = std::is_same_v<T, double>;
  const int n = grid.Size();
  MFEM_VERIFY(A.size() == grid.LocalSize(), "Wrong size of a block cyclic matrix!");
  factors = &A;
  int info = 0;
  if (grid.InGrid() && n > 0)
  {
#if defined(PALACE_WITH_SCALAPACK)
    const int mb = grid.BlockSize(), lld = grid.LeadingDim(), izero = 0, ione = 1;
    descinit_(desc.data(), &n, &n, &mb, &mb, &izero, &izero, &ctxt, &lld, &info);
    if (info == 0)
    {
      if constexpr (real)
      {
        pdpotrf_("L", &n, A.data(), &ione, &ione, desc.data(), &info);
      }
      else
      {
        ipiv.assign(grid.LocalRows() + mb, 0);
        pzgetrf_(&n, &n, A.data(), &ione, &ione, desc.data(), ipiv.data(), &info);
      }
    }
#else
    if constexpr (real)
    {
      dpotrf_("L", &n, A.data(), &n, &info);
    }
    else
    {
      ipiv.assign(n, 0);
      zgetrf_(&n, &n, A.data(), &n, ipiv.data(), &info);
    }
#endif
  }
  info = Failure(info, grid.GetComm());
  MFEM_VERIFY(info == 0, "Dense " << (real ? "Cholesky" : "LU")
                                  << " factorization failed (info = " << info << ")!");
}

template <typename T>
void BlockCyclicFactorization<T>::Solve(std::vector<T> &b, int nb) const
{
  constexpr bool real = std::is_same_v<T, double>;
  const int n = grid.Size();
  if (n == 0 || nb == 0)
  {
    return;
  }
  MPI_Comm comm = grid.GetComm();
  const bool root = Mpi::Root(comm);
#if defined(PALACE_WITH_SCALAPACK)
  // The right-hand sides go to the first column of the grid, which holds the rows of B
  // (one block column), and the solutions back to rank 0.
  const auto type = mpi::DataType<T>();
  const std::size_t size = static_cast<std::size_t>(n) * nb;
  std::vector<T> all(size);
  if (root)
  {
    std::copy(b.begin(), b.begin() + static_cast<std::ptrdiff_t>(size), all.begin());
  }
  MPI_Bcast(all.data(), static_cast<int>(size), type, 0, comm);
  const bool first_col = grid.InGrid() && grid.MyCol() == 0;
  const int mloc = grid.LocalRows(), lld = grid.LeadingDim();
  std::vector<T> B(first_col ? static_cast<std::size_t>(lld) * nb : 1);
  if (first_col)
  {
    for (int k = 0; k < nb; k++)
    {
      for (int il = 0; il < mloc; il++)
      {
        B[il + static_cast<std::size_t>(lld) * k] =
            all[grid.GlobalRow(il) + static_cast<std::size_t>(n) * k];
      }
    }
  }
  int info = 0;
  if (grid.InGrid())
  {
    const int mb = grid.BlockSize(), izero = 0, ione = 1;
    std::array<int, 9> descb{};
    descinit_(descb.data(), &n, &nb, &mb, &nb, &izero, &izero, &ctxt, &lld, &info);
    if (info == 0)
    {
      if constexpr (real)
      {
        pdpotrs_("L", &n, &nb, factors->data(), &ione, &ione, desc.data(), B.data(), &ione,
                 &ione, descb.data(), &info);
      }
      else
      {
        pzgetrs_("N", &n, &nb, factors->data(), &ione, &ione, desc.data(), ipiv.data(),
                 B.data(), &ione, &ione, descb.data(), &info);
      }
    }
  }
  info = Failure(info, comm);
  MFEM_VERIFY(info == 0, "Dense solve failed (info = " << info << ")!");
  std::fill(all.begin(), all.end(), T(0.0));
  if (first_col)
  {
    for (int k = 0; k < nb; k++)
    {
      for (int il = 0; il < mloc; il++)
      {
        all[grid.GlobalRow(il) + static_cast<std::size_t>(n) * k] =
            B[il + static_cast<std::size_t>(lld) * k];
      }
    }
  }
  MPI_Reduce(root ? MPI_IN_PLACE : all.data(), all.data(), static_cast<int>(size), type,
             MPI_SUM, 0, comm);
  if (root)
  {
    std::copy(all.begin(), all.end(), b.begin());
  }
#else
  if (root)
  {
    int info = 0;
    if constexpr (real)
    {
      dpotrs_("L", &n, &nb, factors->data(), &n, b.data(), &n, &info);
    }
    else
    {
      zgetrs_("N", &n, &nb, factors->data(), &n, ipiv.data(), b.data(), &n, &info);
    }
    MFEM_VERIFY(info == 0, "Dense solve failed (info = " << info << ")!");
  }
#endif
}

template std::vector<double> BlockCyclicGrid::ScatterFromRoot(const std::vector<double> &,
                                                              bool) const;
template std::vector<std::complex<double>>
BlockCyclicGrid::ScatterFromRoot(const std::vector<std::complex<double>> &, bool) const;
template class BlockCyclicFactorization<double>;
template class BlockCyclicFactorization<std::complex<double>>;

}  // namespace palace
