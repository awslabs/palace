// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <memory>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include "linalg/ras.hpp"
#include "utils/communication.hpp"

namespace palace
{

TEST_CASE("RAS preconditions a distributed diffusion problem", "[Serial][Parallel]")
{
  // Assemble a small SPD operator with interface coupling when run on multiple ranks. The
  // serial registration also exercises RasSolver's no-interface block-ILU fallback.
  mfem::Mesh serial_mesh = mfem::Mesh::MakeCartesian1D(16);
  mfem::ParMesh mesh(Mpi::World(), serial_mesh);
  mfem::H1_FECollection fec(1, mesh.Dimension());
  mfem::ParFiniteElementSpace fespace(&mesh, &fec);

  mfem::ParBilinearForm form(&fespace);
  form.AddDomainIntegrator(new mfem::DiffusionIntegrator);
  form.AddDomainIntegrator(new mfem::MassIntegrator);
  form.Assemble();
  form.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> A(form.ParallelAssemble());

  RasSolver ras(0);
  ras.SetOperator(*A);

  mfem::Vector exact(fespace.GetTrueVSize()), rhs(fespace.GetTrueVSize());
  mfem::Vector solution(fespace.GetTrueVSize()), error(fespace.GetTrueVSize());
  exact = 1.0;
  solution = 0.0;
  A->Mult(exact, rhs);

  mfem::GMRESSolver solver(Mpi::World());
  solver.SetOperator(*A);
  solver.SetPreconditioner(ras);
  solver.SetRelTol(1.0e-12);
  solver.SetAbsTol(0.0);
  solver.SetMaxIter(100);
  solver.SetKDim(30);
  solver.SetPrintLevel(0);
  solver.Mult(rhs, solution);

  error = solution;
  error -= exact;
  const double local_error = error * error;
  double global_error = 0.0;
  MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_SUM, Mpi::World());
  CHECK(solver.GetConverged());
  CHECK(std::sqrt(global_error) < 1.0e-9);
}

}  // namespace palace
