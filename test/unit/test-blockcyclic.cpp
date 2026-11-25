// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <complex>
#include <type_traits>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include "linalg/blockcyclic.hpp"
#include "utils/communication.hpp"

namespace palace
{

namespace
{

using cplx = std::complex<double>;

// A pseudo-random value in [-1, 1] x [-1, 1] i from (i, j).
cplx Value(int i, int j)
{
  const double a = std::sin(12.9898 * i + 78.233 * j + 1.0) * 43758.5453,
               b = std::sin(39.346 * i + 11.135 * j + 2.0) * 23421.631;
  return {2.0 * (a - std::floor(a)) - 1.0, 2.0 * (b - std::floor(b)) - 1.0};
}

template <typename T>
T Entry(cplx v)
{
  if constexpr (std::is_same_v<T, double>)
  {
    return v.real();
  }
  else
  {
    return v;
  }
}

}  // namespace

TEST_CASE("Block cyclic grids hold each entry once", "[blockcyclic][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const int nranks = Mpi::Size(comm), rank = Mpi::Rank(comm);
  const int n = GENERATE(1, 15, 16, 17, 100);
  const int nproc = GENERATE_COPY(Catch::Generators::range(1, nranks + 1));
  const bool row_major = GENERATE(false, true);
  BlockCyclicGrid grid(comm, n, nproc, 16);
  CHECK(grid.NumProcessRows() * grid.NumProcessCols() <= nproc);

  // The local parts of the matrix with entries (i, j), from rank 0.
  std::vector<cplx> full;
  if (Mpi::Root(comm))
  {
    full.resize(static_cast<std::size_t>(n) * n);
    for (int j = 0; j < n; j++)
    {
      for (int i = 0; i < n; i++)
      {
        full[row_major ? i * n + j : i + j * n] = {double(i), double(j)};
      }
    }
  }
  const auto local = grid.ScatterFromRoot(full, row_major);
  REQUIRE(local.size() == grid.LocalSize());
  // The entries held, those misplaced, and those which Owner and LocalIndex do not find.
  long long count[3] = {0, 0, 0};
  for (int jl = 0; jl < grid.LocalCols(); jl++)
  {
    for (int il = 0; il < grid.LocalRows(); il++)
    {
      const std::size_t q = il + static_cast<std::size_t>(grid.LeadingDim()) * jl;
      const int g = grid.GlobalRow(il), h = grid.GlobalCol(jl);
      count[0]++;
      count[1] += (local[q] != cplx(g, h));
      count[2] += (grid.Owner(g, h) != rank || grid.LocalIndex(g, h) != q);
    }
  }
  Mpi::GlobalSum(3, count, comm);
  CHECK(count[0] == static_cast<long long>(n) * n);
  CHECK(count[1] == 0);
  CHECK(count[2] == 0);
}

TEST_CASE("Block cyclic grids add a matrix given by rows",
          "[blockcyclic][Serial][Parallel]")
{
  // The rows of each process (the last holds none), added to the local parts of grids of
  // one process to all, in chunks of a few entries or of all.
  MPI_Comm comm = Mpi::World();
  const int nranks = Mpi::Size(comm), rank = Mpi::Rank(comm);
  const int n = GENERATE(1, 17, 100);
  const int nproc = GENERATE_COPY(Catch::Generators::range(1, nranks + 1));
  const int chunk = GENERATE(7, 1 << 24);
  BlockCyclicGrid grid(comm, n, nproc, 4);
  std::vector<int> start(nranks + 1, 0);
  for (int r = 0; r < nranks; r++)
  {
    const int weight = (r + 1 < nranks || nranks == 1) ? r + 1 : 0;
    start[r + 1] = start[r] + weight;
  }
  for (int r = 0; r <= nranks; r++)
  {
    start[r] = static_cast<int>(static_cast<long long>(start[r]) * n / start[nranks]);
  }
  const int row0 = start[rank], nrows = start[rank + 1] - start[rank];
  auto f = [](int i, int j) { return 1000.0 * i + j; };
  std::vector<double> rows(static_cast<std::size_t>(nrows) * n);
  for (int r = 0; r < nrows; r++)
  {
    for (int j = 0; j < n; j++)
    {
      rows[static_cast<std::size_t>(r) * n + j] = f(row0 + r, j);
    }
  }
  std::vector<double> local(grid.LocalSize(), 1.0);
  grid.AddFromRows(rows, row0, nrows, local, chunk);
  long long wrong = 0;
  for (int jl = 0; jl < grid.LocalCols(); jl++)
  {
    for (int il = 0; il < grid.LocalRows(); il++)
    {
      wrong += (local[il + static_cast<std::size_t>(grid.LeadingDim()) * jl] !=
                1.0 + f(grid.GlobalRow(il), grid.GlobalCol(jl)));
    }
  }
  Mpi::GlobalSum(1, &wrong, comm);
  CHECK(wrong == 0);
}

TEMPLATE_TEST_CASE("Block cyclic factorizations solve dense systems",
                   "[blockcyclic][Serial][Parallel]", double, std::complex<double>)
{
  // Cholesky's factorization of a real symmetric positive definite matrix, and LU's of a
  // complex one: diagonally dominant, or a cyclic shift with a small diagonal (which needs
  // pivoting).
  using T = TestType;
  constexpr bool real = std::is_same_v<T, double>;
  MPI_Comm comm = Mpi::World();
  const int nranks = Mpi::Size(comm);
  const int n = GENERATE(1, 40, 101);
  const int nproc = GENERATE_COPY(Catch::Generators::range(1, nranks + 1));
  bool pivoting = false;
  if constexpr (!real)
  {
    pivoting = GENERATE(false, true);
  }
  BlockCyclicGrid grid(comm, n, nproc, 16);
  const bool root = Mpi::Root(comm);
  std::vector<T> A;
  if (root)
  {
    A.resize(static_cast<std::size_t>(n) * n);
    for (int j = 0; j < n; j++)
    {
      for (int i = 0; i < n; i++)
      {
        T v = Entry<T>(0.05 * (Value(i, j) + Value(j, i)));
        if (pivoting)
        {
          v += (i == (j + 1) % n) ? 1.0 : 0.0;
          v += (i == j) ? 1.0e-3 : 0.0;
        }
        else
        {
          v += (i == j) ? double(n) : 0.0;
        }
        A[i + static_cast<std::size_t>(n) * j] = v;
      }
    }
  }
  auto local = grid.ScatterFromRoot(A);
  BlockCyclicFactorization<T> factors(grid);
  factors.Factor(local);
  constexpr int nb = 2;
  std::vector<T> b, x;
  if (root)
  {
    b.resize(static_cast<std::size_t>(n) * nb);
    for (std::size_t q = 0; q < b.size(); q++)
    {
      b[q] = Entry<T>(Value(static_cast<int>(q), -1));
    }
    x = b;
  }
  factors.Solve(x, nb);
  double res[2] = {0.0, 0.0};
  if (root)
  {
    for (int k = 0; k < nb; k++)
    {
      for (int i = 0; i < n; i++)
      {
        T r = b[i + static_cast<std::size_t>(n) * k];
        for (int j = 0; j < n; j++)
        {
          r -= A[i + static_cast<std::size_t>(n) * j] *
               x[j + static_cast<std::size_t>(n) * k];
        }
        res[0] += std::norm(r);
        res[1] += std::norm(b[i + static_cast<std::size_t>(n) * k]);
      }
    }
  }
  Mpi::Broadcast(2, res, 0, comm);
  CHECK(std::sqrt(res[0] / res[1]) < 1.0e-12);
}

}  // namespace palace
