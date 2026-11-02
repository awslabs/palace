// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <memory>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include "linalg/amg.hpp"
#include "utils/communication.hpp"

namespace palace
{

namespace
{

// Max-norm of x - y relative to the max-norm of y, over all processes.
double RelMaxDiff(MPI_Comm comm, const mfem::Vector &x, const mfem::Vector &y)
{
  double d[2] = {0.0, 0.0};
  for (int i = 0; i < x.Size(); i++)
  {
    d[0] = std::max(d[0], std::abs(x[i] - y[i]));
    d[1] = std::max(d[1], std::abs(y[i]));
  }
  Mpi::GlobalMax(2, d, comm);
  return d[1] > 0.0 ? d[0] / d[1] : d[0];
}

}  // namespace

TEST_CASE("BoomerAMG zero initial guess", "[amg][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  auto smesh = mfem::Mesh::MakeCartesian3D(10, 8, 6, mfem::Element::TETRAHEDRON);
  mfem::ParMesh pmesh(comm, smesh);
  mfem::H1_FECollection fec(1, pmesh.Dimension());
  mfem::ParFiniteElementSpace fespace(&pmesh, &fec);
  mfem::Array<int> dbc_marker(pmesh.bdr_attributes.Max()), dbc_tdof_list;
  dbc_marker = 0;
  dbc_marker[0] = 1;
  fespace.GetEssentialTrueDofs(dbc_marker, dbc_tdof_list);

  mfem::ConstantCoefficient one(1.0);
  mfem::ParBilinearForm a(&fespace);
  a.AddDomainIntegrator(new mfem::DiffusionIntegrator(one));
  a.Assemble();
  a.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> A(a.ParallelAssemble());
  std::unique_ptr<mfem::HypreParMatrix> Ae(A->EliminateRowsCols(dbc_tdof_list));

  BoomerAmgSolver amg(1, 1, true, 0);
  amg.SetOperator(*A);
  const int n = A->Height();
  mfem::Vector b0(n), x0(n), z0(n);
  b0.Randomize(1 + Mpi::Rank(comm));

  // Zero initial guess marked as zero for the first relaxation.
  amg.Mult(b0, x0);

  // Same result when the zero initial guess is given explicitly (not marked).
  amg.iterative_mode = true;
  z0 = 0.0;
  amg.Mult(b0, z0);
  amg.iterative_mode = false;
  CHECK(RelMaxDiff(comm, z0, x0) == 0.0);
}

}  // namespace palace
