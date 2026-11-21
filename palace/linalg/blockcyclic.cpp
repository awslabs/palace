// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "blockcyclic.hpp"

#include <cmath>
#include <mfem.hpp>
#include "utils/communication.hpp"

extern "C"
{
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
  void pzgetrf_(const int *, const int *, std::complex<double> *, const int *, const int *,
                const int *, int *, int *);
  void pzgetrs_(const char *, const int *, const int *, const std::complex<double> *,
                const int *, const int *, const int *, const int *, std::complex<double> *,
                const int *, const int *, const int *, int *);
#endif
}

namespace palace
{

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

std::vector<std::complex<double>>
BlockCyclicGrid::ScatterFromRoot(const std::vector<std::complex<double>> &full,
                                 bool row_major) const
{
  // Each local part as columns of its leading dimension, so that the counts are columns.
  const auto type = mpi::DataType<std::complex<double>>();
  auto send = [&](const std::vector<std::complex<double>> &part, int lld, int ncol, int to)
  {
    MPI_Datatype column;
    MPI_Type_contiguous(lld, type, &column);
    MPI_Type_commit(&column);
    MPI_Send(part.data(), std::max(1, ncol), column, to, 0, comm);
    MPI_Type_free(&column);
  };
  std::vector<std::complex<double>> local(LocalSize(), 0.0);
  if (Mpi::Root(comm))
  {
    MFEM_VERIFY(full.size() == static_cast<std::size_t>(n) * n,
                "Wrong size of a matrix to distribute!");
    std::vector<std::complex<double>> block;
    for (int q = nprow * npcol - 1; q >= 0; q--)
    {
      const int pr = q / npcol, pc = q % npcol, mq = LocalCount(n, mb, pr, nprow),
                nq = LocalCount(n, mb, pc, npcol), lld = std::max(1, mq);
      auto &part = (q == 0) ? local : block;
      part.assign(static_cast<std::size_t>(lld) * std::max(1, nq), 0.0);
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

BlockCyclicLU::BlockCyclicLU(const BlockCyclicGrid &grid) : grid(grid)
{
#if defined(PALACE_WITH_SCALAPACK)
  system = Csys2blacs_handle(grid.GetComm());
  ctxt = system;
  Cblacs_gridinit(&ctxt, "Row", grid.NumProcessRows(), grid.NumProcessCols());
#endif
}

BlockCyclicLU::~BlockCyclicLU()
{
#if defined(PALACE_WITH_SCALAPACK)
  if (grid.InGrid())
  {
    Cblacs_gridexit(ctxt);
  }
  Cfree_blacs_system_handle(system);
#endif
}

void BlockCyclicLU::Factor(std::vector<std::complex<double>> &A)
{
  const int n = grid.Size();
  MFEM_VERIFY(A.size() == grid.LocalSize(), "Wrong size of a block cyclic matrix!");
  factors = &A;
  int info = 0;
  if (grid.InGrid() && n > 0)
  {
#if defined(PALACE_WITH_SCALAPACK)
    const int mb = grid.BlockSize(), lld = grid.LeadingDim(), izero = 0, ione = 1;
    descinit_(desc.data(), &n, &n, &mb, &mb, &izero, &izero, &ctxt, &lld, &info);
    ipiv.assign(grid.LocalRows() + mb, 0);
    if (info == 0)
    {
      pzgetrf_(&n, &n, A.data(), &ione, &ione, desc.data(), ipiv.data(), &info);
    }
#else
    ipiv.assign(n, 0);
    zgetrf_(&n, &n, A.data(), &n, ipiv.data(), &info);
#endif
  }
  Mpi::GlobalMax(1, &info, grid.GetComm());
  MFEM_VERIFY(info == 0, "Dense LU factorization failed (info = " << info << ")!");
}

void BlockCyclicLU::Solve(std::vector<std::complex<double>> &b, int nb) const
{
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
  const auto type = mpi::DataType<std::complex<double>>();
  const std::size_t size = static_cast<std::size_t>(n) * nb;
  std::vector<std::complex<double>> all(size);
  if (root)
  {
    std::copy(b.begin(), b.begin() + static_cast<std::ptrdiff_t>(size), all.begin());
  }
  MPI_Bcast(all.data(), static_cast<int>(size), type, 0, comm);
  const bool first_col = grid.InGrid() && grid.MyCol() == 0;
  const int mloc = grid.LocalRows(), lld = grid.LeadingDim();
  std::vector<std::complex<double>> B(first_col ? static_cast<std::size_t>(lld) * nb : 1);
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
      pzgetrs_("N", &n, &nb, factors->data(), &ione, &ione, desc.data(), ipiv.data(),
               B.data(), &ione, &ione, descb.data(), &info);
    }
  }
  Mpi::GlobalMax(1, &info, comm);
  MFEM_VERIFY(info == 0, "Dense LU solve failed (info = " << info << ")!");
  std::fill(all.begin(), all.end(), 0.0);
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
    zgetrs_("N", &n, &nb, factors->data(), &n, ipiv.data(), b.data(), &n, &info);
    MFEM_VERIFY(info == 0, "Dense LU solve failed (info = " << info << ")!");
  }
#endif
}

}  // namespace palace
