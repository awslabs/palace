// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <memory>
#include <utility>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "linalg/mumps.hpp"
#include "linalg/mumpscomplex.hpp"
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "linalg/strumpack.hpp"
#include "linalg/strumpackcomplex.hpp"
#include "linalg/superlu.hpp"
#include "linalg/superlucomplex.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

#if defined(MFEM_USE_STRUMPACK) || defined(MFEM_USE_SUPERLU) || defined(MFEM_USE_MUMPS)

namespace palace
{

namespace
{

// Complex-valued indefinite matrix K - k² (1 - i τ(x)) M (+ i C with a convection matrix C
// for an unsymmetric matrix) of an H1 space on a hexahedral mesh of the unit cube, with a
// loss τ(x) in part of the domain, as a ComplexOperator with HypreParMatrix real and
// imaginary parts.
std::unique_ptr<ComplexOperator> AssembleTestMatrix(mfem::ParFiniteElementSpace &fespace,
                                                    double imag_scale = 1.0,
                                                    bool unsymmetric = false)
{
  constexpr double k2 = 20.0;
  mfem::ConstantCoefficient one(1.0), mk2(-k2);
  mfem::FunctionCoefficient tau([](const mfem::Vector &x)
                                { return (x(2) > 0.5) ? 5.0 * (x(2) - 0.5) : 0.0; });
  mfem::ProductCoefficient k2tau(k2 * imag_scale, tau);
  mfem::Vector b(3);
  b(0) = 1.0;
  b(1) = 2.0;
  b(2) = 3.0;
  mfem::VectorConstantCoefficient beta(b);
  mfem::ParBilinearForm ar(&fespace), ai(&fespace);
  ar.AddDomainIntegrator(new mfem::DiffusionIntegrator(one));
  ar.AddDomainIntegrator(new mfem::MassIntegrator(mk2));
  ai.AddDomainIntegrator(new mfem::MassIntegrator(k2tau));
  if (unsymmetric)
  {
    ai.AddDomainIntegrator(new mfem::ConvectionIntegrator(beta));
  }
  ar.Assemble();
  ar.Finalize();
  ai.Assemble();
  ai.Finalize();
  return std::make_unique<ComplexWrapperOperator>(
      std::unique_ptr<Operator>(ar.ParallelAssemble()),
      std::unique_ptr<Operator>(ai.ParallelAssemble()));
}

// Diagonal matrix (2 + i) M_L with the lumped mass matrix M_L, which has a different
// sparsity pattern than the matrices of AssembleTestMatrix.
std::unique_ptr<ComplexOperator>
AssembleDiagonalTestMatrix(mfem::ParFiniteElementSpace &fespace)
{
  mfem::ConstantCoefficient two(2.0), one(1.0);
  mfem::ParBilinearForm ar(&fespace), ai(&fespace);
  ar.AddDomainIntegrator(new mfem::LumpedIntegrator(new mfem::MassIntegrator(two)));
  ai.AddDomainIntegrator(new mfem::LumpedIntegrator(new mfem::MassIntegrator(one)));
  ar.Assemble();
  ar.Finalize();
  ai.Assemble();
  ai.Finalize();
  return std::make_unique<ComplexWrapperOperator>(
      std::unique_ptr<Operator>(ar.ParallelAssemble()),
      std::unique_ptr<Operator>(ai.ParallelAssemble()));
}

// The solution of the complex arithmetic solver matches the solution of the solver of the
// equivalent real-valued system, also for repeated factorizations of matrices with the same
// sparsity pattern and with a different sparsity pattern.
void TestComplexDirectSolver(Solver<ComplexOperator> &solver, Solver<ComplexOperator> &ref,
                             double tol, bool unsymmetric = false)
{
  MPI_Comm comm = Mpi::World();
  mfem::Mesh serial_mesh =
      mfem::Mesh::MakeCartesian3D(6, 6, 6, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0);
  mfem::ParMesh mesh(comm, serial_mesh);
  mfem::H1_FECollection fec(2, 3);
  mfem::ParFiniteElementSpace fespace(&mesh, &fec);
  auto A = AssembleTestMatrix(fespace, 1.0, unsymmetric);

  ComplexVector b(A->Height()), x(A->Height()), x_ref(A->Height()), r(A->Height());
  b.Real().Randomize(1);
  b.Imag().Randomize(2);
  const double norm_b = linalg::Norml2(comm, b);
  auto Residual = [&](const ComplexOperator &M)
  {
    M.Mult(x, r);
    r -= b;
    return linalg::Norml2(comm, r) / norm_b;
  };

  solver.SetOperator(*A);
  solver.Mult(b, x);
  ref.SetOperator(*A);
  ref.Mult(b, x_ref);
  CHECK(Residual(*A) < tol);
  x -= x_ref;
  CHECK(linalg::Norml2(comm, x) / linalg::Norml2(comm, x_ref) < 10.0 * tol);

  // Same sparsity pattern, different values (reuse of the reordering).
  auto A2 = AssembleTestMatrix(fespace, 2.0, unsymmetric);
  solver.SetOperator(*A2);
  solver.Mult(b, x);
  CHECK(Residual(*A2) < tol);

  // Different sparsity pattern, and back.
  auto D = AssembleDiagonalTestMatrix(fespace);
  solver.SetOperator(*D);
  solver.Mult(b, x);
  CHECK(Residual(*D) < tol);
  solver.SetOperator(*A);
  solver.Mult(b, x);
  CHECK(Residual(*A) < tol);
}

}  // namespace

#if defined(MFEM_USE_STRUMPACK)
TEST_CASE("STRUMPACK complex arithmetic solver", "[strumpack][Serial][Parallel]")
{
  // Complex symmetric and unsymmetric matrices, with and without the reuse of the
  // reordering.
  MPI_Comm comm = Mpi::World();
  const bool unsymmetric = GENERATE(false, true);
  const bool reorder_reuse = GENERATE(true, false);
  StrumpackComplexSolver solver(comm, SymbolicFactorization::DEFAULT, reorder_reuse, 0);
  MfemWrapperSolver<ComplexOperator> ref(
      std::make_unique<StrumpackSolver>(comm, SymbolicFactorization::DEFAULT,
                                        SparseCompression::NONE, 0.0, 0, 0, true, 0),
      true, true, false, true);
  TestComplexDirectSolver(solver, ref, 1.0e-11, unsymmetric);
}

// As for StrumpackSolver (see test-strumpack.cpp), the parallel matrix products of hypre
// can store explicit zeros without their transposes, so the sparsity pattern given to
// STRUMPACK is not declared symmetric: a five-point Laplacian with a loss, plus explicit
// zeros which couple distant unknowns in one direction only.
TEST_CASE("STRUMPACK complex arithmetic solver unsymmetric sparsity pattern",
          "[strumpack][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const int P = Mpi::Size(comm), rank = Mpi::Rank(comm);
  const int m = 24;
  const HYPRE_BigInt N = m * m;
  const HYPRE_BigInt first = (N * rank) / P, last = (N * (rank + 1)) / P;
  const int n_loc = static_cast<int>(last - first);
  std::vector<int> Ir(n_loc + 1, 0), Ii(n_loc + 1, 0);
  std::vector<HYPRE_BigInt> Jr, Ji;
  std::vector<double> Vr, Vi;
  for (HYPRE_BigInt i = first; i < last; i++)
  {
    const HYPRE_BigInt ix = i % m, iy = i / m;
    std::vector<std::pair<HYPRE_BigInt, double>> row = {{i, 4.0}};
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
      Jr.push_back(j);
      Vr.push_back(v);
    }
    Ir[i - first + 1] = static_cast<int>(Jr.size());
    Ji.push_back(i);
    Vi.push_back(0.5);
    Ii[i - first + 1] = static_cast<int>(Ji.size());
  }
  HYPRE_BigInt row_starts[2] = {first, last};
  auto *Ar = new mfem::HypreParMatrix(comm, n_loc, N, N, Ir.data(), Jr.data(), Vr.data(),
                                      row_starts, row_starts);
  auto *Ai = new mfem::HypreParMatrix(comm, n_loc, N, N, Ii.data(), Ji.data(), Vi.data(),
                                      row_starts, row_starts);
  ComplexWrapperOperator A{std::unique_ptr<Operator>(Ar), std::unique_ptr<Operator>(Ai)};

  ComplexVector x_ex(n_loc), b(n_loc), x(n_loc), r(n_loc);
  for (int i = 0; i < n_loc; i++)
  {
    const double t = static_cast<double>(first + i) / static_cast<double>(N);
    x_ex.Real()[i] = std::sin(7.0 * t) + t;
    x_ex.Imag()[i] = std::cos(3.0 * t);
  }
  A.Mult(x_ex, b);
  StrumpackComplexSolver solver(comm, SymbolicFactorization::DEFAULT, true, 0);
  for (int k = 0; k < 2; k++)
  {
    // Factorization, and refactorization reusing the reordering.
    solver.SetOperator(A);
    x = 0.0;
    solver.Mult(b, x);
    r = x;
    r -= x_ex;
    CHECK(linalg::Norml2(comm, r) / linalg::Norml2(comm, x_ex) < 1.0e-12);
  }
}
#endif

#if defined(MFEM_USE_SUPERLU)
TEST_CASE("SuperLU complex arithmetic solver", "[superlu][Serial][Parallel]")
{
  // Complex symmetric and unsymmetric matrices, with and without the reuse of the row and
  // column permutations.
  MPI_Comm comm = Mpi::World();
  const bool unsymmetric = GENERATE(false, true);
  const bool reorder_reuse = GENERATE(true, false);
  SuperLUComplexSolver solver(comm, SymbolicFactorization::DEFAULT, reorder_reuse, 0);
  MfemWrapperSolver<ComplexOperator> ref(
      std::make_unique<SuperLUSolver>(comm, SymbolicFactorization::DEFAULT, false, true, 0),
      true, true, false, true);
  TestComplexDirectSolver(solver, ref, 1.0e-10, unsymmetric);
}
#endif

#if defined(MFEM_USE_MUMPS)
TEST_CASE("MUMPS complex arithmetic solver", "[mumps][Serial][Parallel]")
{
  // Complex symmetric and unsymmetric matrices.
  MPI_Comm comm = Mpi::World();
  const auto sym = GENERATE(MatrixSymmetry::SYMMETRIC, MatrixSymmetry::UNSYMMETRIC);
  MumpsComplexSolver solver(comm, sym, SymbolicFactorization::DEFAULT, 0.0, true, 0);
  MfemWrapperSolver<ComplexOperator> ref(
      std::make_unique<MumpsSolver>(comm, sym, SymbolicFactorization::DEFAULT, 0.0, true,
                                    0),
      true, true, false, true);
  TestComplexDirectSolver(solver, ref, 1.0e-10, sym == MatrixSymmetry::UNSYMMETRIC);
}
#endif

}  // namespace palace

#endif
