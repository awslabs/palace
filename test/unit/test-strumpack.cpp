// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "linalg/strumpack.hpp"
#include "utils/communication.hpp"

#if defined(MFEM_USE_STRUMPACK)

namespace palace
{

namespace
{

// Number of stored entries (i, j) of a square hypre_ParCSRMatrix whose transpose (j, i) is
// not stored.
HYPRE_BigInt CountEntriesWithoutTranspose(const mfem::HypreParMatrix &A)
{
  MPI_Comm comm = A.GetComm();
  const int P = Mpi::Size(comm);
  auto *parcsr = (hypre_ParCSRMatrix *)const_cast<mfem::HypreParMatrix &>(A);
  hypre_CSRMatrix *csr = hypre_MergeDiagAndOffd(parcsr);
  hypre_CSRMatrixMigrate(csr, HYPRE_MEMORY_HOST);
  const HYPRE_Int n = hypre_CSRMatrixNumRows(csr);
  const HYPRE_Int *I = hypre_CSRMatrixI(csr);
  const HYPRE_BigInt *J = hypre_CSRMatrixBigJ(csr);
  const HYPRE_BigInt first = hypre_ParCSRMatrixFirstRowIndex(parcsr);
  std::vector<HYPRE_BigInt> starts(P + 1);
  MPI_Allgather(&first, 1, HYPRE_MPI_BIG_INT, starts.data(), 1, HYPRE_MPI_BIG_INT, comm);
  starts[P] = hypre_ParCSRMatrixGlobalNumRows(parcsr);

  // Send every (j, i) to the owner of row j, which checks whether it is stored.
  std::vector<std::vector<HYPRE_BigInt>> send(P);
  for (HYPRE_Int r = 0; r < n; r++)
  {
    for (HYPRE_Int k = I[r]; k < I[r + 1]; k++)
    {
      const int q = static_cast<int>(std::upper_bound(starts.begin(), starts.end(), J[k]) -
                                     starts.begin() - 1);
      send[q].push_back(J[k]);
      send[q].push_back(first + r);
    }
  }
  std::vector<int> scnt(P), rcnt(P), sdsp(P, 0), rdsp(P, 0);
  for (int q = 0; q < P; q++)
  {
    scnt[q] = static_cast<int>(send[q].size());
  }
  MPI_Alltoall(scnt.data(), 1, MPI_INT, rcnt.data(), 1, MPI_INT, comm);
  std::vector<HYPRE_BigInt> sbuf, rbuf;
  for (int q = 0; q < P; q++)
  {
    sdsp[q] = static_cast<int>(sbuf.size());
    sbuf.insert(sbuf.end(), send[q].begin(), send[q].end());
    rdsp[q] = (q > 0) ? rdsp[q - 1] + rcnt[q - 1] : 0;
  }
  rbuf.resize(rdsp[P - 1] + rcnt[P - 1]);
  MPI_Alltoallv(sbuf.data(), scnt.data(), sdsp.data(), HYPRE_MPI_BIG_INT, rbuf.data(),
                rcnt.data(), rdsp.data(), HYPRE_MPI_BIG_INT, comm);
  HYPRE_BigInt missing = 0;
  for (std::size_t t = 0; t < rbuf.size(); t += 2)
  {
    const HYPRE_Int r = static_cast<HYPRE_Int>(rbuf[t] - first);
    missing += (std::find(J + I[r], J + I[r + 1], rbuf[t + 1]) == J + I[r + 1]);
  }
  hypre_CSRMatrixDestroy(csr);
  Mpi::GlobalSum(1, &missing, comm);
  return missing;
}

}  // namespace

// STRUMPACK builds its elimination tree from the sparsity pattern. The parallel RAP can
// store explicit zeros without their transposes (for example hypre's host matrix product
// inserts zero "diagonal" entries into the off-diagonal block of a process whose number of
// rows equals its number of off-diagonal columns), so the pattern passed to STRUMPACK is
// not guaranteed to be symmetric. Declaring such a pattern symmetric gave an invalid
// elimination tree, and factorization or solve failures.
TEST_CASE("StrumpackSolver unsymmetric sparsity pattern", "[strumpack][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const int P = Mpi::Size(comm), rank = Mpi::Rank(comm);

  // Five-point Laplacian on an m x m grid, rows distributed in contiguous blocks, plus
  // explicit zeros which couple distant unknowns in one direction only.
  const int m = 24;
  const HYPRE_BigInt N = m * m;
  std::vector<HYPRE_BigInt> starts(P + 1);
  for (int q = 0; q <= P; q++)
  {
    starts[q] = (N * q) / P;
  }
  const HYPRE_BigInt first = starts[rank], last = starts[rank + 1];
  const int n_loc = static_cast<int>(last - first);
  std::vector<int> I(n_loc + 1, 0);
  std::vector<HYPRE_BigInt> J;
  std::vector<double> data;
  for (HYPRE_BigInt i = first; i < last; i++)
  {
    const HYPRE_BigInt ix = i % m, iy = i / m;
    std::vector<std::pair<HYPRE_BigInt, double>> row;
    row.emplace_back(i, 4.0);
    if (ix > 0)
    {
      row.emplace_back(i - 1, -1.0);
    }
    if (ix < m - 1)
    {
      row.emplace_back(i + 1, -1.0);
    }
    if (iy > 0)
    {
      row.emplace_back(i - m, -1.0);
    }
    if (iy < m - 1)
    {
      row.emplace_back(i + m, -1.0);
    }
    // One-sided explicit zeros: (i, N - 1 - i) for some i in the first half of the grid,
    // and (i, i - N / 2) for some i in the second half. Neither has its transpose.
    if (i < N / 2 && i % 5 == 0)
    {
      row.emplace_back(N - 1 - i, 0.0);
    }
    if (i >= N / 2 && i % 7 == 3 && (i - N / 2) % 5 != 0)
    {
      row.emplace_back(i - N / 2, 0.0);
    }
    std::sort(row.begin(), row.end());
    for (const auto &[j, v] : row)
    {
      J.push_back(j);
      data.push_back(v);
    }
    I[i - first + 1] = static_cast<int>(J.size());
  }
  HYPRE_BigInt row_starts[2] = {first, last};
  mfem::HypreParMatrix A(comm, n_loc, N, N, I.data(), J.data(), data.data(), row_starts,
                         row_starts);
  REQUIRE(CountEntriesWithoutTranspose(A) > 0);

  // Manufactured solution.
  mfem::HypreParVector x_ex(A), b(A), x(A);
  for (int i = 0; i < n_loc; i++)
  {
    const double t = static_cast<double>(first + i) / static_cast<double>(N);
    x_ex(i) = std::sin(7.0 * t) + t;
  }
  A.Mult(x_ex, b);

  StrumpackSolver strumpack(comm, SymbolicFactorization::DEFAULT, SparseCompression::NONE,
                            0.0, 0, 0, false, true, 0);
  strumpack.SetOperator(A);
  x = 0.0;
  strumpack.Mult(b, x);

  // Solve a second time with the same factorization, and after refactorizing the same
  // matrix.
  mfem::HypreParVector x2(A);
  x2 = 0.0;
  strumpack.Mult(b, x2);
  strumpack.SetOperator(A);
  mfem::HypreParVector x3(A);
  x3 = 0.0;
  strumpack.Mult(b, x3);

  double err[3] = {0.0, 0.0, 0.0}, nrm = 0.0;
  for (int i = 0; i < n_loc; i++)
  {
    err[0] = std::max(err[0], std::abs(x(i) - x_ex(i)));
    err[1] = std::max(err[1], std::abs(x2(i) - x_ex(i)));
    err[2] = std::max(err[2], std::abs(x3(i) - x_ex(i)));
    nrm = std::max(nrm, std::abs(x_ex(i)));
  }
  Mpi::GlobalMax(3, err, comm);
  Mpi::GlobalMax(1, &nrm, comm);
  CHECK(err[0] <= 1.0e-12 * nrm);
  CHECK(err[1] <= 1.0e-12 * nrm);
  CHECK(err[2] <= 1.0e-12 * nrm);
}

// A small matrix on many processes, as the coarse matrix of a multigrid hierarchy built by
// uniform refinement of a small mesh: the separator tree of the column ordering has fewer
// leaves than processes, so only some processes get a subgraph and a local tree. The two
// assignments must agree for every ordering (the sequential METIS ordering of a matrix with
// fewer than P / 2 leaves failed in the symbolic factorization). Exercised on 4 or more
// processes.
TEST_CASE("StrumpackSolver small matrices on many processes",
          "[strumpack][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const int P = Mpi::Size(comm), rank = Mpi::Rank(comm);
  const auto reorder =
      GENERATE(SymbolicFactorization::METIS, SymbolicFactorization::AMD,
               SymbolicFactorization::PTSCOTCH, SymbolicFactorization::DEFAULT);
  for (int m = 2; m <= 20; m++)
  {
    // Five-point Laplacian on an m x m grid, at least one row on every process.
    const HYPRE_BigInt N = m * m;
    if (N < P)
    {
      continue;
    }
    INFO("Grid " << m << " x " << m << ", ordering " << static_cast<int>(reorder));
    const HYPRE_BigInt first = (N * rank) / P, last = (N * (rank + 1)) / P;
    const int n_loc = static_cast<int>(last - first);
    std::vector<int> I(n_loc + 1, 0);
    std::vector<HYPRE_BigInt> J;
    std::vector<double> data;
    for (HYPRE_BigInt i = first; i < last; i++)
    {
      const HYPRE_BigInt ix = i % m, iy = i / m;
      const HYPRE_BigInt cols[5] = {i - m, i - 1, i, i + 1, i + m};
      const bool stored[5] = {iy > 0, ix > 0, true, ix < m - 1, iy < m - 1};
      for (int k = 0; k < 5; k++)
      {
        if (stored[k])
        {
          J.push_back(cols[k]);
          data.push_back((k == 2) ? 4.0 : -1.0);
        }
      }
      I[i - first + 1] = static_cast<int>(J.size());
    }
    HYPRE_BigInt row_starts[2] = {first, last};
    mfem::HypreParMatrix A(comm, n_loc, N, N, I.data(), J.data(), data.data(), row_starts,
                           row_starts);

    mfem::HypreParVector x_ex(A), b(A), x(A);
    for (int i = 0; i < n_loc; i++)
    {
      const double t = static_cast<double>(first + i) / static_cast<double>(N);
      x_ex(i) = std::cos(5.0 * t) + t;
    }
    A.Mult(x_ex, b);
    StrumpackSolver strumpack(comm, reorder, SparseCompression::NONE, 0.0, 0, 0, false,
                              true, 0);
    strumpack.SetOperator(A);
    x = 0.0;
    strumpack.Mult(b, x);

    double err = 0.0, nrm = 0.0;
    for (int i = 0; i < n_loc; i++)
    {
      err = std::max(err, std::abs(x(i) - x_ex(i)));
      nrm = std::max(nrm, std::abs(x_ex(i)));
    }
    Mpi::GlobalMax(1, &err, comm);
    Mpi::GlobalMax(1, &nrm, comm);
    CHECK(err <= 1.0e-12 * nrm);
  }
}

}  // namespace palace

#endif
