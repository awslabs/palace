// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "linalg/superlu.hpp"
#include "utils/communication.hpp"
#include "utils/labels.hpp"

#if defined(MFEM_USE_SUPERLU)

using namespace palace;

namespace
{

// A 54 x 54 coarse multigrid matrix of an eigenmode solve, as captured at the input of
// SuperLU (test/data/superlu): 966 stored entries, 808 of them exact zeros. This process
// takes its rows [first, first + num), in the stored order.
struct CoarseMatrix
{
  std::vector<int> I;
  std::vector<HYPRE_BigInt> J;
  std::vector<double> a;
  std::unique_ptr<mfem::HypreParMatrix> A;

  CoarseMatrix(const std::string &name, HYPRE_BigInt first, HYPRE_BigInt num)
  {
    std::ifstream is(std::string(PALACE_TEST_DATA_DIR) + "/superlu/" + name);
    REQUIRE(is.good());
    I.assign(num + 1, 0);
    for (std::string line; std::getline(is, line);)
    {
      if (line.empty() || line[0] == '#')
      {
        continue;
      }
      HYPRE_BigInt i, j;
      double v;
      std::istringstream(line) >> i >> j >> v;
      if (i >= first && i < first + num)  // rows are stored in increasing order
      {
        I[i - first + 1]++;
        J.push_back(j);
        a.push_back(v);
      }
    }
    for (HYPRE_BigInt i = 0; i < num; i++)
    {
      I[i + 1] += I[i];
    }
    HYPRE_BigInt rows[2] = {first, first + num};
    A = std::make_unique<mfem::HypreParMatrix>(Mpi::World(), static_cast<int>(num), 54, 54,
                                               I.data(), J.data(), a.data(), rows, rows);
    REQUIRE(A->NNZ() == 966);  // the stored zeros are kept (808 of the 966 entries)
  }
};

// The rows [first, first + num) of this process: all on the first process (the others
// have none), or split between the processes.
std::pair<HYPRE_BigInt, HYPRE_BigInt> Rows(bool split)
{
  const HYPRE_BigInt np = Mpi::Size(Mpi::World()), rank = Mpi::Rank(Mpi::World());
  if (!split)
  {
    return {rank == 0 ? 0 : 54, rank == 0 ? 54 : 0};
  }
  return {rank * (54 / np) + std::min(rank, 54 % np), 54 / np + (rank < 54 % np)};
}

// ||A x_k - b_k|| / ||b_k|| for the solutions x_k of the solver for the right-hand sides
// b_k,i = 1 + sin((1 + k)(1 + i)), k < num_rhs (one Mult, or one ArrayMult), with the
// original matrix (its stored zeros do not change its action).
std::vector<double> SolveResiduals(const SuperLUSolver &solver,
                                   const mfem::HypreParMatrix &A, HYPRE_BigInt first,
                                   int num_rhs)
{
  std::vector<mfem::Vector> b(num_rhs, mfem::Vector(A.Height())), x(b);
  mfem::Array<const mfem::Vector *> B(num_rhs);
  mfem::Array<mfem::Vector *> X(num_rhs);
  for (int k = 0; k < num_rhs; k++)
  {
    for (int i = 0; i < A.Height(); i++)
    {
      b[k](i) = 1.0 + std::sin((1.0 + k) * (1.0 + static_cast<double>(first + i)));
    }
    B[k] = &b[k];
    X[k] = &x[k];
  }
  if (num_rhs == 1)
  {
    solver.Mult(b[0], x[0]);
  }
  else
  {
    solver.ArrayMult(B, X);
  }
  std::vector<double> res;
  for (int k = 0; k < num_rhs; k++)
  {
    mfem::Vector r(A.Height());
    A.Mult(x[k], r);
    r -= b[k];
    double norms[2] = {r * r, b[k] * b[k]};
    Mpi::GlobalSum(2, norms, Mpi::World());
    res.push_back(std::sqrt(norms[0] / norms[1]));
  }
  return res;
}

double SolveResidual(const SuperLUSolver &solver, const mfem::HypreParMatrix &A,
                     HYPRE_BigInt first)
{
  return SolveResiduals(solver, A, first, 1)[0];
}

std::unique_ptr<SuperLUSolver> MakeSolver()
{
  return std::make_unique<SuperLUSolver>(Mpi::World(), SymbolicFactorization::DEFAULT,
                                         false, true, 0);
}

}  // namespace

// The stored zeros of the coarse matrix made the factorization inaccurate (relative
// residual 0.55 for a matrix with condition number 33), and reusing the row permutation of
// the first matrix for the second gave an exactly zero pivot, U(54, 54), although the
// second has condition number 92.
TEST_CASE("SuperLU factors coarse eigenmode matrices with stored zeros",
          "[superlu][Serial][Parallel]")
{
  const bool split = GENERATE(false, true);
  CAPTURE(split);
  const auto [first, num] = Rows(split);
  const CoarseMatrix seed("coarse-seed.txt", first, num),
      target("coarse-target.txt", first, num);

  auto solver = MakeSolver();
  solver->SetOperator(*seed.A);
  CHECK(SolveResidual(*solver, *seed.A, first) < 1.0e-12);
  solver->SetOperator(*target.A);
  CHECK(SolveResidual(*solver, *target.A, first) < 1.0e-12);

  // A second operator before any solve replaces the first.
  auto other = MakeSolver();
  other->SetOperator(*seed.A);
  other->SetOperator(*target.A);
  CHECK(SolveResidual(*other, *target.A, first) < 1.0e-12);

  // A solver given an operator but destroyed without a solve.
  MakeSolver()->SetOperator(*seed.A);

  // Another distribution of the rows is an error, on all processes at once.
  if (Mpi::Size(Mpi::World()) > 1)
  {
    const auto [other_first, other_num] = Rows(!split);
    const CoarseMatrix moved("coarse-target.txt", other_first, other_num);
    CHECK_THROWS(solver->SetOperator(*moved.A));
  }
}

// The number of right-hand sides is fixed by the first solve: the native solver checks it
// only after a new operator has released the previous factors, so it is checked first, and
// a rejected solve leaves the operator to the next one.
TEST_CASE("SuperLU solves with a fixed number of right-hand sides",
          "[superlu][Serial][Parallel]")
{
  const bool split = GENERATE(false, true);
  CAPTURE(split);
  const auto [first, num] = Rows(split);
  const CoarseMatrix seed("coarse-seed.txt", first, num),
      target("coarse-target.txt", first, num);

  // Two right-hand sides for every operator.
  auto solver = MakeSolver();
  solver->SetOperator(*seed.A);
  for (const double res : SolveResiduals(*solver, *seed.A, first, 2))
  {
    CHECK(res < 1.0e-12);
  }
  solver->SetOperator(*target.A);
  for (const double res : SolveResiduals(*solver, *target.A, first, 2))
  {
    CHECK(res < 1.0e-12);
  }

  // No right-hand sides is an error before the operator is used (the solver is then
  // destroyed without a solve).
  {
    auto empty = MakeSolver();
    empty->SetOperator(*seed.A);
    mfem::Array<const mfem::Vector *> B;
    mfem::Array<mfem::Vector *> X;
    CHECK_THROWS(empty->ArrayMult(B, X));
    CHECK_THROWS(empty->ArrayMultTranspose(B, X));
  }

  // Another number of right-hand sides is an error, and the operator is still solved next.
  auto single = MakeSolver();
  single->SetOperator(*seed.A);
  CHECK(SolveResidual(*single, *seed.A, first) < 1.0e-12);
  single->SetOperator(*target.A);
  CHECK_THROWS(SolveResiduals(*single, *target.A, first, 2));
  CHECK(SolveResidual(*single, *target.A, first) < 1.0e-12);
}

#endif
