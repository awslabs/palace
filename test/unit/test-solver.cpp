// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

namespace palace
{

namespace
{

// Stands for a parallel direct solver, and counts how a process applies it.
class CountingSolver : public mfem::Solver
{
public:
  mutable std::array<int, 2> calls = {0, 0};  // Mult, ArrayMult

  void SetOperator(const mfem::Operator &op) override
  {
    height = op.Height();
    width = op.Width();
  }

  void Mult(const mfem::Vector &x, mfem::Vector &y) const override
  {
    calls[0]++;
    y = x;
  }

  void ArrayMult(const mfem::Array<const mfem::Vector *> &X,
                 mfem::Array<mfem::Vector *> &Y) const override
  {
    calls[1]++;
    for (int i = 0; i < X.Size(); i++)
    {
      *Y[i] = *X[i];
    }
  }
};

// The identity matrix of size n with all of its rows on the root process.
std::unique_ptr<mfem::HypreParMatrix> RootIdentity(MPI_Comm comm, HYPRE_BigInt n)
{
  const int n_loc = Mpi::Root(comm) ? static_cast<int>(n) : 0;
  std::vector<int> I(n_loc + 1);
  std::vector<HYPRE_BigInt> J(n_loc);
  std::vector<double> data(n_loc, 1.0);
  for (int i = 0; i < n_loc; i++)
  {
    I[i] = i;
    J[i] = i;
  }
  I[n_loc] = n_loc;
  HYPRE_BigInt row_starts[2] = {Mpi::Root(comm) ? 0 : n, n};
  return std::make_unique<mfem::HypreParMatrix>(comm, n_loc, n, n, I.data(), J.data(),
                                                data.data(), row_starts, row_starts);
}

}  // namespace

// The processes of a parallel solver must all make the same calls to it. On the coarsest
// multigrid level of a small mesh on many processes, some processes own no unknowns, and
// their local sizes do not tell the real-valued block form of a complex system, of twice
// the size of the vectors, from a system of the size of the vectors.
TEST_CASE("MfemWrapperSolver applies a complex system alike on all processes",
          "[solver][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const HYPRE_BigInt N = 4;
  const int n_loc = Mpi::Root(comm) ? static_cast<int>(N) : 0;
  const auto Ar = RootIdentity(comm, N), Ai = RootIdentity(comm, N),
             A_block = RootIdentity(comm, 2 * N);

  // The calls of this process to the solver for one application to a vector.
  auto Apply = [n_loc](const ComplexOperator &A, bool complex_matrix, bool caller = false)
  {
    auto counting = std::make_unique<CountingSolver>();
    const auto *counts = counting.get();
    MfemWrapperSolver<ComplexOperator> solver(std::move(counting), true, complex_matrix);
    solver.SetCallerBlockForm(caller);
    solver.SetOperator(A);
    ComplexVector x(n_loc), y(n_loc);
    x = 1.0;
    y = 0.0;
    solver.Mult(x, y);
    return counts->calls;
  };
  auto OnAllProcesses = [comm](std::array<int, 2> calls, std::array<int, 2> expected)
  {
    int same = (calls == expected);
    Mpi::GlobalMin(1, &same, comm);
    return same == 1;
  };

  // The block form assembled by the solver, and the real-valued approximation.
  const ComplexWrapperOperator A(Ar.get(), Ai.get());
  CHECK(OnAllProcesses(Apply(A, true), {1, 0}));
  CHECK(OnAllProcesses(Apply(A, false), {0, 1}));

  // A block form assembled by the caller, as the real part of the operator.
  const ComplexWrapperOperator A_caller(A_block.get(), nullptr);
  CHECK(OnAllProcesses(Apply(A_caller, false, true), {1, 0}));
}

}  // namespace palace
