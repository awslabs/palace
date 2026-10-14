// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <memory>
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
  // Double precision, and single precision with iterative refinement.
  MPI_Comm comm = Mpi::World();
  const bool single = GENERATE(false, true);
  StrumpackComplexSolver solver(comm, SymbolicFactorization::DEFAULT, single,
                                single ? 2 : 0, 0);
  MfemWrapperSolver<ComplexOperator> ref(
      std::make_unique<StrumpackSolver>(comm, SymbolicFactorization::DEFAULT,
                                        SparseCompression::NONE, 0.0, 0, 0, true, 0),
      true, true, false, true);
  TestComplexDirectSolver(solver, ref, single ? 1.0e-8 : 1.0e-11);
}
#endif

#if defined(MFEM_USE_SUPERLU)
TEST_CASE("SuperLU complex arithmetic solver", "[superlu][Serial][Parallel]")
{
  // With and without the accuracy check of the factorizations.
  MPI_Comm comm = Mpi::World();
  const bool robust = GENERATE(false, true);
  SuperLUComplexSolver solver(comm, SymbolicFactorization::DEFAULT, true, 0, robust);
  MfemWrapperSolver<ComplexOperator> ref(
      std::make_unique<SuperLUSolver>(comm, SymbolicFactorization::DEFAULT, false, true, 0),
      true, true, false, true);
  TestComplexDirectSolver(solver, ref, 1.0e-10);
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
