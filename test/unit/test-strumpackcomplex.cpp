// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <memory>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "linalg/strumpack.hpp"
#include "linalg/strumpackcomplex.hpp"
#include "linalg/superlu.hpp"
#include "linalg/superlucomplex.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

#if defined(MFEM_USE_STRUMPACK) || defined(MFEM_USE_SUPERLU)

namespace palace
{

namespace
{

// Complex symmetric indefinite matrix K - k² (1 - i τ(x)) M of an H1 space on a hexahedral
// mesh of the unit cube, with a loss τ(x) in part of the domain, as a ComplexOperator with
// HypreParMatrix real and imaginary parts.
std::unique_ptr<ComplexOperator> AssembleTestMatrix(mfem::ParFiniteElementSpace &fespace)
{
  constexpr double k2 = 20.0;
  mfem::ConstantCoefficient one(1.0), mk2(-k2);
  mfem::FunctionCoefficient tau([](const mfem::Vector &x)
                                { return (x(2) > 0.5) ? 5.0 * (x(2) - 0.5) : 0.0; });
  mfem::ProductCoefficient k2tau(k2, tau);
  mfem::ParBilinearForm ar(&fespace), ai(&fespace);
  ar.AddDomainIntegrator(new mfem::DiffusionIntegrator(one));
  ar.AddDomainIntegrator(new mfem::MassIntegrator(mk2));
  ai.AddDomainIntegrator(new mfem::MassIntegrator(k2tau));
  ar.Assemble();
  ar.Finalize();
  ai.Assemble();
  ai.Finalize();
  return std::make_unique<ComplexWrapperOperator>(
      std::unique_ptr<Operator>(ar.ParallelAssemble()),
      std::unique_ptr<Operator>(ai.ParallelAssemble()));
}

}  // namespace

#if defined(MFEM_USE_STRUMPACK)
TEST_CASE("STRUMPACK complex arithmetic solver", "[strumpack][Serial][Parallel]")
{
  // The solution of the complex arithmetic solver matches the solution of the solver of the
  // equivalent real-valued system, in double precision and in single precision with
  // iterative refinement.
  MPI_Comm comm = Mpi::World();
  const bool single = GENERATE(false, true);
  mfem::Mesh serial_mesh =
      mfem::Mesh::MakeCartesian3D(6, 6, 6, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0);
  mfem::ParMesh mesh(comm, serial_mesh);
  mfem::H1_FECollection fec(2, 3);
  mfem::ParFiniteElementSpace fespace(&mesh, &fec);
  auto A = AssembleTestMatrix(fespace);

  ComplexVector b(A->Height()), x(A->Height()), x_ref(A->Height()), r(A->Height());
  b.Real().Randomize(1);
  b.Imag().Randomize(2);

  StrumpackComplexSolver solver(comm, SymbolicFactorization::DEFAULT, single,
                                single ? 2 : 0, 0);
  solver.SetOperator(*A);
  solver.Mult(b, x);

  MfemWrapperSolver<ComplexOperator> ref(
      std::make_unique<StrumpackSolver>(comm, SymbolicFactorization::DEFAULT,
                                        SparseCompression::NONE, 0.0, 0, 0, true, 0),
      true, true, false, true);
  ref.SetOperator(*A);
  ref.Mult(b, x_ref);

  // Relative residuals and difference.
  A->Mult(x, r);
  r -= b;
  const double norm_b = linalg::Norml2(comm, b);
  CHECK(linalg::Norml2(comm, r) / norm_b < (single ? 1.0e-8 : 1.0e-11));
  x -= x_ref;
  CHECK(linalg::Norml2(comm, x) / linalg::Norml2(comm, x_ref) <
        (single ? 1.0e-7 : 1.0e-10));

  // Factorization of a matrix with the same sparsity pattern and different values (reuse of
  // the reordering).
  auto A2 = AssembleTestMatrix(fespace);
  *const_cast<mfem::HypreParMatrix *>(
      dynamic_cast<const mfem::HypreParMatrix *>(A2->Imag())) *= 2.0;
  solver.SetOperator(*A2);
  solver.Mult(b, x);
  A2->Mult(x, r);
  r -= b;
  CHECK(linalg::Norml2(comm, r) / norm_b < (single ? 1.0e-8 : 1.0e-11));
}
#endif

#if defined(MFEM_USE_SUPERLU)
TEST_CASE("SuperLU complex arithmetic solver", "[superlu][Serial][Parallel]")
{
  // The solution of the complex arithmetic solver matches the solution of the solver of the
  // equivalent real-valued system, also for repeated factorizations with the same sparsity
  // pattern, with and without the accuracy check of the factorizations.
  MPI_Comm comm = Mpi::World();
  const bool robust = GENERATE(false, true);
  mfem::Mesh serial_mesh =
      mfem::Mesh::MakeCartesian3D(6, 6, 6, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0);
  mfem::ParMesh mesh(comm, serial_mesh);
  mfem::H1_FECollection fec(2, 3);
  mfem::ParFiniteElementSpace fespace(&mesh, &fec);
  auto A = AssembleTestMatrix(fespace);

  ComplexVector b(A->Height()), x(A->Height()), x_ref(A->Height()), r(A->Height());
  b.Real().Randomize(1);
  b.Imag().Randomize(2);

  SuperLUComplexSolver solver(comm, SymbolicFactorization::DEFAULT, true, 0, robust);
  solver.SetOperator(*A);
  solver.Mult(b, x);

  MfemWrapperSolver<ComplexOperator> ref(
      std::make_unique<SuperLUSolver>(comm, SymbolicFactorization::DEFAULT, false, true, 0),
      true, true, false, true);
  ref.SetOperator(*A);
  ref.Mult(b, x_ref);

  A->Mult(x, r);
  r -= b;
  const double norm_b = linalg::Norml2(comm, b);
  CHECK(linalg::Norml2(comm, r) / norm_b < 1.0e-10);
  x -= x_ref;
  CHECK(linalg::Norml2(comm, x) / linalg::Norml2(comm, x_ref) < 1.0e-9);

  auto A2 = AssembleTestMatrix(fespace);
  *const_cast<mfem::HypreParMatrix *>(
      dynamic_cast<const mfem::HypreParMatrix *>(A2->Imag())) *= 2.0;
  solver.SetOperator(*A2);
  solver.Mult(b, x);
  A2->Mult(x, r);
  r -= b;
  CHECK(linalg::Norml2(comm, r) / norm_b < 1.0e-10);
}
#endif

}  // namespace palace

#endif
