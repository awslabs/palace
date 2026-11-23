// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "linalg/strumpack.hpp"
#include "test-helpers.hpp"
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

  // The solver stores the permutation once per node, when STRUMPACK is built with that
  // option. The tests run on one node, where the permutation is not sent to other nodes:
  // with STRUMPACK_SHARED_PERM_GROUP=1, every process acts as a node, so that on several
  // processes the root broadcasts the permutation to the other node leaders as on several
  // nodes. STRUMPACK reads the variable on every process at each reordering, so it is set
  // before the solver is constructed and restored after the solver is destroyed.
  const bool process_as_node = GENERATE(false, true);
  CAPTURE(process_as_node);
  ScopedEnvironmentVariable group("STRUMPACK_SHARED_PERM_GROUP",
                                  process_as_node ? "1" : nullptr);

  // With and without the selective inversion of the distributed fronts, which has no
  // effect on one process or when STRUMPACK is built without that option.
  const bool selective_inversion = GENERATE(false, true);
  CAPTURE(selective_inversion);
  StrumpackSolver strumpack(comm, SymbolicFactorization::DEFAULT, SparseCompression::NONE,
                            0.0, 0, 0, false, selective_inversion, 0);
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

// STRUMPACK 8.1.0 did not apply the row interchanges of a diagonal tile of a distributed
// BLR front to the tiles of L on its left, while the solve applies all the interchanges to
// the right-hand side first. With BLR compression on two or more processes, the
// factorization of an indefinite matrix was then inaccurate whatever the compression
// tolerance (fixed by extern/patch/STRUMPACK/strumpack_blr_pivoting.diff). On one process
// no front is distributed and this checks the BLR factorization of the sequential code,
// which did not have the defect.
TEST_CASE("StrumpackSolver BLR compression of an indefinite matrix",
          "[strumpack][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const int P = Mpi::Size(comm), rank = Mpi::Rank(comm);

  // K - 104 M, where K and M are the stiffness and mass matrices of the Laplacian on the
  // unit cube with natural boundary conditions, for n^3 hexahedral elements of order 2
  // (tensor products of the 1D matrices), with the rows distributed in contiguous blocks.
  // Its largest separator, a plane of (2 n + 1)^2 = 625 unknowns, is above the 512 unknowns
  // from which STRUMPACK compresses a front, and is the only compressed front. The matrix
  // is indefinite, so that the factorization interchanges rows. The shift is between two
  // eigenvalues of the discrete problem, 98.7 and 108.6 (10 pi^2 and 11 pi^2 for the
  // Laplacian), with 35 eigenvalues below it: the factorization of the compressed front
  // then does 16 to 29 row interchanges on 1 to 8 processes. A small shift is a weaker
  // test: with 20, just above the eigenvalue 2 pi^2, it does 1 to 3 on most process counts.
  const int n = 12, m = 2 * n + 1;
  const double h = 1.0 / n, shift = 104.0;
  const double K_e[3][3] = {{7.0, -8.0, 1.0}, {-8.0, 16.0, -8.0}, {1.0, -8.0, 7.0}};
  const double M_e[3][3] = {{4.0, 2.0, -1.0}, {2.0, 16.0, 2.0}, {-1.0, 2.0, 4.0}};
  std::vector<double> K1(m * m, 0.0), M1(m * m, 0.0);
  for (int e = 0; e < n; e++)
  {
    for (int a = 0; a < 3; a++)
    {
      for (int b = 0; b < 3; b++)
      {
        K1[(2 * e + a) * m + 2 * e + b] += K_e[a][b] / (3.0 * h);
        M1[(2 * e + a) * m + 2 * e + b] += M_e[a][b] * h / 30.0;
      }
    }
  }
  const HYPRE_BigInt N = m * m * m;
  const HYPRE_BigInt first = (N * rank) / P, last = (N * (rank + 1)) / P;
  const int n_loc = static_cast<int>(last - first);
  std::vector<int> I(n_loc + 1, 0);
  std::vector<HYPRE_BigInt> J;
  std::vector<double> data;
  for (HYPRE_BigInt g = first; g < last; g++)
  {
    const int i = static_cast<int>(g % m), j = static_cast<int>((g / m) % m),
              k = static_cast<int>(g / (m * m));
    for (int k2 = std::max(0, k - 2); k2 <= std::min(m - 1, k + 2); k2++)
    {
      for (int j2 = std::max(0, j - 2); j2 <= std::min(m - 1, j + 2); j2++)
      {
        for (int i2 = std::max(0, i - 2); i2 <= std::min(m - 1, i + 2); i2++)
        {
          // An entry of the 1D mass matrix is nonzero exactly when its two nodes belong to
          // one element.
          const double Mx = M1[i * m + i2], My = M1[j * m + j2], Mz = M1[k * m + k2];
          if (Mx == 0.0 || My == 0.0 || Mz == 0.0)
          {
            continue;
          }
          const double Kx = K1[i * m + i2], Ky = K1[j * m + j2], Kz = K1[k * m + k2];
          J.push_back(i2 + m * (j2 + m * k2));
          data.push_back(Kx * My * Mz + Mx * Ky * Mz + Mx * My * Kz - shift * Mx * My * Mz);
        }
      }
    }
    I[g - first + 1] = static_cast<int>(J.size());
  }
  HYPRE_BigInt row_starts[2] = {first, last};
  mfem::HypreParMatrix A(comm, n_loc, N, N, I.data(), J.data(), data.data(), row_starts,
                         row_starts);

  mfem::HypreParVector b(A), x(A), r(A);
  for (int i = 0; i < n_loc; i++)
  {
    b(i) = std::sin(static_cast<double>((first + i) * (first + i)));
  }

  // The residual of one solve is 1e-10 to 3e-10 on 1 to 8 processes at this compression
  // tolerance. With the defect it was 0.23 on 2 processes, 0.44 on 4 and 0.46 on 8, and at
  // least 0.045 on 2 to 8 processes.
  StrumpackSolver strumpack(comm, SymbolicFactorization::DEFAULT, SparseCompression::BLR,
                            1.0e-9, 0, 0, false, false, 0);
  strumpack.SetOperator(A);
  x = 0.0;
  strumpack.Mult(b, x);
  A.Mult(x, r);
  r -= b;
  const double res =
      std::sqrt(mfem::InnerProduct(comm, r, r) / mfem::InnerProduct(comm, b, b));
  CHECK(res <= 1.0e-6);
}

}  // namespace palace

#endif
