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
  void dpotrf_(const char *, const int *, double *, const int *, int *);
  void dsyevr_(const char *, const char *, const char *, const int *, double *, const int *,
               const double *, const double *, const int *, const int *, const double *,
               int *, double *, double *, const int *, int *, double *, const int *, int *,
               const int *, int *);
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
  void pdpotrf_(const char *, const int *, double *, const int *, const int *, const int *,
                int *);
  void pdsyevx_(const char *, const char *, const char *, const int *, double *,
                const int *, const int *, const int *, const double *, const double *,
                const int *, const int *, const double *, int *, int *, double *,
                const double *, double *, const int *, const int *, const int *, double *,
                const int *, int *, const int *, int *, int *, double *, int *);
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

namespace
{

#if defined(PALACE_WITH_SCALAPACK)
// The BLACS context of a grid, while in scope, and the ScaLAPACK descriptor of an n x n
// matrix on it. Collective.
class BlacsGrid
{
public:
  explicit BlacsGrid(const BlockCyclicGrid &grid) : grid(grid)
  {
    system = Csys2blacs_handle(grid.GetComm());
    ctxt = system;
    Cblacs_gridinit(&ctxt, "Row", grid.NumProcessRows(), grid.NumProcessCols());
  }
  ~BlacsGrid()
  {
    if (grid.InGrid())
    {
      Cblacs_gridexit(ctxt);
    }
    Cfree_blacs_system_handle(system);
  }
  BlacsGrid(const BlacsGrid &) = delete;
  BlacsGrid &operator=(const BlacsGrid &) = delete;

  std::array<int, 9> Descriptor(int &info) const
  {
    const int n = grid.Size(), mb = grid.BlockSize(), lld = grid.LeadingDim(), izero = 0;
    std::array<int, 9> desc{};
    descinit_(desc.data(), &n, &n, &mb, &mb, &izero, &izero, &ctxt, &lld, &info);
    return desc;
  }

private:
  const BlockCyclicGrid &grid;
  int system = -1, ctxt = -1;
};
#endif

}  // namespace

bool PositiveDefinite(const BlockCyclicGrid &grid, std::vector<double> &A)
{
  const int n = grid.Size();
  MFEM_VERIFY(A.size() == grid.LocalSize(), "Wrong size of a block cyclic matrix!");
  int info = 0;
#if defined(PALACE_WITH_SCALAPACK)
  const BlacsGrid blacs(grid);
  if (grid.InGrid() && n > 0)
  {
    const auto desc = blacs.Descriptor(info);
    const int ione = 1;
    if (info == 0)
    {
      pdpotrf_("L", &n, A.data(), &ione, &ione, desc.data(), &info);
    }
  }
#else
  if (grid.InGrid() && n > 0)
  {
    dpotrf_("L", &n, A.data(), &n, &info);
  }
#endif
  // info > 0: a leading minor is not positive definite.
  int failed = (info < 0);
  Mpi::GlobalMax(1, &failed, grid.GetComm());
  MFEM_VERIFY(!failed, "Dense Cholesky factorization failed (info = " << info << ")!");
  Mpi::GlobalMax(1, &info, grid.GetComm());
  return info == 0;
}

double LeastEigenvalue(const BlockCyclicGrid &grid, std::vector<double> &A)
{
  const int n = grid.Size();
  MFEM_VERIFY(A.size() == grid.LocalSize(), "Wrong size of a block cyclic matrix!");
  const int one = 1;
  const double zero = 0.0;
  int found = 1, info = 0, lwork = -1, liwork = -1;
  double z = 0.0;
  // W needs n entries even for one eigenvalue (the bisection stores the whole bracket).
  std::vector<double> w(std::max(n, 1), 0.0), work(3);
  std::vector<int> iwork(1);
#if defined(PALACE_WITH_SCALAPACK)
  const BlacsGrid blacs(grid);
  if (grid.InGrid() && n > 0)
  {
    // The least eigenvalue only, without eigenvectors (Z is not referenced). The workspace
    // query also writes WORK(1:3), the tolerance and bounds that rank 0 broadcasts.
    const auto desc = blacs.Descriptor(info);
    const int nprocs = grid.NumProcessRows() * grid.NumProcessCols();
    const double orfac = 1.0e-3;
    int nz = 0;
    std::vector<int> ifail(n), iclustr(2 * nprocs);
    std::vector<double> gap(nprocs);
    for (int pass = 0; pass < 2 && info == 0; pass++)
    {
      pdsyevx_("N", "I", "L", &n, A.data(), &one, &one, desc.data(), &zero, &zero, &one,
               &one, &zero, &found, &nz, w.data(), &orfac, &z, &one, &one, desc.data(),
               work.data(), &lwork, iwork.data(), &liwork, ifail.data(), iclustr.data(),
               gap.data(), &info);
      if (pass == 0)
      {
        lwork = static_cast<int>(work[0]);
        liwork = iwork[0];
        work.resize(lwork);
        iwork.resize(liwork);
      }
    }
  }
#else
  if (grid.InGrid() && n > 0)
  {
    std::vector<int> isuppz(2);
    for (int pass = 0; pass < 2 && info == 0; pass++)
    {
      dsyevr_("N", "I", "L", &n, A.data(), &n, &zero, &zero, &one, &one, &zero, &found,
              w.data(), &z, &one, isuppz.data(), work.data(), &lwork, iwork.data(), &liwork,
              &info);
      if (pass == 0)
      {
        lwork = static_cast<int>(work[0]);
        liwork = iwork[0];
        work.resize(lwork);
        iwork.resize(liwork);
      }
    }
  }
#endif
  int failed = (info != 0 || (n > 0 && found != 1));
  Mpi::GlobalMax(1, &failed, grid.GetComm());
  MFEM_VERIFY(!failed,
              "Dense symmetric eigenvalue computation failed (info = " << info << ")!");
  // The grid includes rank 0, which holds the eigenvalue.
  Mpi::Broadcast(1, w.data(), 0, grid.GetComm());
  return (n > 0) ? w[0] : 0.0;
}

}  // namespace palace
