// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <memory>
#include <tuple>
#include <vector>
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

TEST_CASE("BoomerAMG zero initial guess and multiple right-hand sides",
          "[amg][Serial][Parallel]")
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

  const int n = A->Height();
  mfem::Vector b0(n), b1(n), x0(n), x1(n), y0(n), y1(n), z0(n);
  b0.Randomize(1 + Mpi::Rank(comm));
  b1.Randomize(101 + Mpi::Rank(comm));

  // One or several V-cycles and sweeps (the zero flag must not survive the first
  // relaxation), with zero, one, or two levels of aggressive coarsening.
  for (auto [cycle_it, smooth_it, agg_levels] :
       {std::tuple{1, 1, 1}, std::tuple{2, 2, 1}, std::tuple{3, 1, 0}, std::tuple{1, 1, 2}})
  {
    BoomerAmgSolver amg(cycle_it, smooth_it, agg_levels, 0);
    amg.SetOperator(*A);
    for (int rep = 0; rep < 2; rep++)
    {
      // Zero initial guess marked as zero for the first relaxation.
      x0 = 1.0;
      amg.Mult(b0, x0);

      // Same result when the zero initial guess is given explicitly (not marked).
      amg.iterative_mode = true;
      z0 = 0.0;
      amg.Mult(b0, z0);
      amg.iterative_mode = false;
      CHECK(RelMaxDiff(comm, z0, x0) == 0.0);
    }
  }
}

TEST_CASE("BoomerAMG multiple right-hand sides", "[amg][Serial][Parallel]")
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

  const int n = A->Height(), nv_max = 3;
  std::vector<mfem::Vector> b(nv_max, mfem::Vector(n)), x(nv_max, mfem::Vector(n)),
      y(nv_max, mfem::Vector(n));
  for (int k = 0; k < nv_max; k++)
  {
    b[k].Randomize(1 + 100 * k + Mpi::Rank(comm));
  }
  mfem::Vector z(n);

  for (auto [cycle_it, smooth_it, agg_levels] :
       {std::tuple{1, 1, 1}, std::tuple{2, 2, 1}, std::tuple{3, 1, 0}, std::tuple{1, 1, 2}})
  {
    BoomerAmgSolver amg(cycle_it, smooth_it, agg_levels, 0);
    amg.SetOperator(*A);
    for (int k = 0; k < nv_max; k++)
    {
      amg.Mult(b[k], x[k]);
    }

    // Several right-hand sides in one V-cycle, alternating with single applications (the
    // communication pattern switches between one and several components) and with a
    // different number of components (the multicomponent vectors are reallocated).
    for (int nv : {2, 3, 2, 1})
    {
      mfem::Array<const mfem::Vector *> B(nv);
      mfem::Array<mfem::Vector *> Y(nv);
      for (int k = 0; k < nv; k++)
      {
        B[k] = &b[k];
        Y[k] = &y[k];
        y[k] = 1.0;
      }
      amg.ArrayMult(B, Y);
      for (int k = 0; k < nv; k++)
      {
        CHECK(RelMaxDiff(comm, y[k], x[k]) <= 1.0e-12);
      }
      amg.Mult(b[nv - 1], z);
      CHECK(RelMaxDiff(comm, z, x[nv - 1]) == 0.0);
    }
  }
}

// The hybrid Gauss-Seidel relaxation skips the exchange of a vector flagged as zero. AMS
// relaxes its iterate on the Nedelec space with several sweeps of hybrid Gauss-Seidel
// directly (hypre_ParCSRRelax, not through hypre_BoomerAMGRelax), so the flag must be
// cleared once the iterate is updated, or later sweeps use zero off-processor values.
TEST_CASE("AMS relaxation of a zero-flagged initial guess", "[amg][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  auto smesh = mfem::Mesh::MakeCartesian3D(6, 5, 4, mfem::Element::TETRAHEDRON);
  mfem::ParMesh pmesh(comm, smesh);
  mfem::ND_FECollection fec(1, pmesh.Dimension());
  mfem::ParFiniteElementSpace fespace(&pmesh, &fec);

  mfem::ConstantCoefficient one(1.0);
  mfem::ParBilinearForm a(&fespace);
  a.AddDomainIntegrator(new mfem::CurlCurlIntegrator(one));
  a.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(one));
  a.Assemble();
  a.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> A(a.ParallelAssemble());

  mfem::HypreAMS ams(*A, &fespace);
  ams.SetPrintLevel(0);
  // l1-scaled hybrid symmetric Gauss-Seidel with several sweeps.
  HYPRE_AMSSetSmoothingOptions(ams, 2, 3, 1.0, 1.0);
  ams.iterative_mode = true;

  mfem::HypreParVector b(*A), x(*A), z(*A);
  b.Randomize(1 + Mpi::Rank(comm));

  // Explicit zero initial guess, not flagged.
  z = 0.0;
  hypre_ParVectorAllZeros((hypre_ParVector *)z) = 0;
  ams.Mult(b, z);

  // Same initial guess flagged as zero.
  x.HypreWrite();
  hypre_ParVectorSetZeros(x);
  ams.Mult(b, x);
  CHECK(RelMaxDiff(comm, x, z) == 0.0);
}

}  // namespace palace
