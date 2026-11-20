// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <complex>
#include <vector>
#include <mfem.hpp>
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

}  // namespace

TEST_CASE("Block cyclic grids hold each entry once", "[blockcyclic][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const int nranks = Mpi::Size(comm);
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
  long long count[2] = {0, 0};  // the entries held, and those misplaced
  for (int jl = 0; jl < grid.LocalCols(); jl++)
  {
    for (int il = 0; il < grid.LocalRows(); il++)
    {
      count[0]++;
      const cplx v = local[il + static_cast<std::size_t>(grid.LeadingDim()) * jl];
      count[1] += (v != cplx(grid.GlobalRow(il), grid.GlobalCol(jl)));
    }
  }
  Mpi::GlobalSum(2, count, comm);
  CHECK(count[0] == static_cast<long long>(n) * n);
  CHECK(count[1] == 0);
}

TEST_CASE("Block cyclic LU solves dense systems", "[blockcyclic][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const int nranks = Mpi::Size(comm);
  const int n = GENERATE(1, 40, 101);
  const int nproc = GENERATE_COPY(Catch::Generators::range(1, nranks + 1));
  // Diagonally dominant, or a cyclic shift with a small diagonal (which needs pivoting).
  const bool pivoting = GENERATE(false, true);
  BlockCyclicGrid grid(comm, n, nproc, 16);
  const bool root = Mpi::Root(comm);
  std::vector<cplx> A;
  if (root)
  {
    A.resize(static_cast<std::size_t>(n) * n);
    for (int j = 0; j < n; j++)
    {
      for (int i = 0; i < n; i++)
      {
        cplx v = 0.1 * Value(i, j);
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
  BlockCyclicLU lu(grid);
  lu.Factor(local);
  constexpr int nb = 2;
  std::vector<cplx> b, x;
  if (root)
  {
    b.resize(static_cast<std::size_t>(n) * nb);
    for (std::size_t q = 0; q < b.size(); q++)
    {
      b[q] = Value(static_cast<int>(q), -1);
    }
    x = b;
  }
  lu.Solve(x, nb);
  double res[2] = {0.0, 0.0};
  if (root)
  {
    for (int k = 0; k < nb; k++)
    {
      for (int i = 0; i < n; i++)
      {
        cplx r = b[i + static_cast<std::size_t>(n) * k];
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
