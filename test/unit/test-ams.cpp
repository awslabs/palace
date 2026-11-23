// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <tuple>
#include <utility>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "linalg/ams.hpp"
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

#if defined(HYPRE_AMS_AMG_REDIST)
// Check which processes the hypre patch for the coarse-level agglomeration gathers the
// scalar Pi-space hierarchies of AMS onto. The additive Pi-space corrections are solved
// concurrently when their hierarchies are gathered from the finest level on: every group
// of processes then has a different owner for each hierarchy, so that no process works on
// two of them, and the group sizes differ by at most one. A hierarchy which is not solved
// concurrently with the others is only gathered from the finest level on if that level is
// small.
void CheckAgglomeration(MPI_Comm comm, const HypreAmsSolver &ams, int target)
{
  const auto *data =
      reinterpret_cast<const hypre_AMSData *>(static_cast<HYPRE_Solver>(ams));
  int num_pi = 0, num_agglomerated = 0, num_concurrent = 0, num_planned = 0, owned = 0;
  for (auto [B, A] :
       {std::pair{data->B_Pix, data->A_Pix}, std::pair{data->B_Piy, data->A_Piy},
        std::pair{data->B_Piz, data->A_Piz}})
  {
    if (!A)
    {
      continue;
    }
    num_pi++;
    HYPRE_Int level, nconc, group_size;
    hypre_BoomerAMGRedistGetInfo(B, &level, &nconc, &group_size);
    num_agglomerated += (level >= 0);
    if (nconc == 1)
    {
      CHECK((level != 0 || hypre_ParCSRMatrixGlobalNumRows(A) < 2 * target));
      continue;
    }
    num_concurrent++;
    num_planned = nconc;
    owned += (group_size > 0);
    int group_min = (group_size > 0) ? group_size : std::numeric_limits<int>::max();
    int group_max = group_size;
    Mpi::GlobalMin(1, &group_min, comm);
    Mpi::GlobalMax(1, &group_max, comm);
    CHECK(level == 0);
    CHECK(group_min >= nconc);
    CHECK(group_max <= group_min + 1);
  }
  Mpi::GlobalMax(1, &owned, comm);
  CHECK(owned <= 1);
  CHECK(num_concurrent == num_planned);

  // hypre leaves a hierarchy distributed, without a message, if its cycle depends on the
  // data distribution (a hybrid Gauss-Seidel relaxation, for example). The hierarchies of
  // these tests are all agglomerated on two or more processes.
  if (Mpi::Size(comm) > 1)
  {
    CHECK(num_agglomerated == num_pi);
  }

  // With at least two processes for each hierarchy, the corrections are always concurrent
  // (if hypre can agglomerate all hierarchies: not the ones with a coarsest level of
  // several times the target, when their coarsening stalls).
  if (Mpi::Size(comm) >= 2 * num_pi && num_agglomerated == num_pi)
  {
    CHECK(num_concurrent == num_pi);
  }
}
#endif

// Apply HypreAmsSolver with the coarse levels of its auxiliary space AMG hierarchies
// agglomerated (config["Solver"]["Linear"]["AMSAgglomerationTarget"] > 0), and compare with
// the hierarchies kept distributed (0): the hierarchies and the cycles are the same, so the
// results agree up to the summation order.
void TestAgglomeration(mfem::Mesh &smesh)
{
  MPI_Comm comm = Mpi::World();
  Mesh mesh(comm, smesh);
  constexpr int order = 1;
  // The discrete gradient for AMS is assembled with libCEED, which needs the default
  // quadrature order (set by IoData in a simulation).
  fem::DefaultIntegrationOrder::p_trial = order;
  mfem::ND_FECollection nd_fec(order, mesh.Dimension());
  mfem::H1_FECollection h1_fec(order, mesh.Dimension());
  FiniteElementSpace nd_fespace(mesh, &nd_fec), h1_fespace(mesh, &h1_fec);
  mfem::Array<int> dbc_marker(mesh.Get().bdr_attributes.Max()), dbc_tdof_list;
  dbc_marker = 1;
  nd_fespace.Get().GetEssentialTrueDofs(dbc_marker, dbc_tdof_list);

  // With and without the mass term (singular for magnetostatics, without the G-space
  // hierarchy), and with one or several smoothing sweeps.
  for (const auto &params :
       {std::pair{true, 1}, std::pair{false, 1}, std::pair{true, 3}, std::pair{false, 3}})
  {
    // Not a structured binding: Clang cannot capture one in a lambda with OpenMP.
    const bool singular = params.first;
    const int smooth_it = params.second;
    mfem::ConstantCoefficient one(1.0);
    mfem::ParBilinearForm a(&nd_fespace.Get());
    a.AddDomainIntegrator(new mfem::CurlCurlIntegrator(one));
    if (!singular)
    {
      a.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(one));
    }
    a.Assemble();
    a.Finalize();
    std::unique_ptr<mfem::HypreParMatrix> A(a.ParallelAssemble());
    std::unique_ptr<mfem::HypreParMatrix> Ae(A->EliminateRowsCols(dbc_tdof_list));

    const int n = A->Height();
    mfem::Vector b(n), x_ref(n), x(n);
    b.Randomize(1 + Mpi::Rank(comm));
    auto Solver = [&](int target)
    {
      auto ams = std::make_unique<HypreAmsSolver>(
          nd_fespace, h1_fespace, /*cycle_it*/ 2, smooth_it, /*vector_interp*/ false,
          singular, /*agg_coarsen*/ singular, target, /*print*/ 0);
      ams->SetOperator(*A);
      return ams;
    };

    // Reference: hierarchies kept distributed during the solve.
    Solver(0)->Mult(b, x_ref);

    // The default target gathers these small hierarchies entirely, each onto one process,
    // and the small one gathers them onto several processes first.
    for (int target : {2048, 64})
    {
      auto ams = Solver(target);
      for (int rep = 0; rep < 2; rep++)
      {
        x = 1.0;
        ams->Mult(b, x);
        CAPTURE(singular, smooth_it, target, rep);
        CHECK(RelMaxDiff(comm, x, x_ref) <= 1.0e-10);
      }
#if defined(HYPRE_AMS_AMG_REDIST)
      CheckAgglomeration(comm, *ams, target);
#endif
    }
  }
}

}  // namespace

TEST_CASE("AMS with agglomerated auxiliary space coarse levels", "[ams][Serial][Parallel]")
{
  SECTION("3D")
  {
    auto smesh = mfem::Mesh::MakeCartesian3D(12, 10, 8, mfem::Element::TETRAHEDRON);
    TestAgglomeration(smesh);
  }
  SECTION("2D")
  {
    auto smesh = mfem::Mesh::MakeCartesian2D(40, 32, mfem::Element::TRIANGLE);
    TestAgglomeration(smesh);
  }
}

// The Pi-space hierarchies are only gathered onto groups of processes with a process for
// each hierarchy from 6 processes on (3D), which is more than the processes of the
// [Parallel] tests: run at the regression rank count too.
TEST_CASE("AMS with agglomerated auxiliary space coarse levels on more processes",
          "[ams][Parallel][Regression]")
{
  SECTION("3D")
  {
    auto smesh = mfem::Mesh::MakeCartesian3D(12, 10, 8, mfem::Element::TETRAHEDRON);
    TestAgglomeration(smesh);
  }
  SECTION("2D")
  {
    auto smesh = mfem::Mesh::MakeCartesian2D(40, 32, mfem::Element::TRIANGLE);
    TestAgglomeration(smesh);
  }
}

// A singular operator (curl-curl without a mass term or essential DOFs) on two separate
// parts, a box and a small column, like the inside of a via enclosed by the Dirichlet
// boundary of a broken space. The Galerkin product of the constant of the column can be an
// exactly zero row on a coarse level of the Pi-space hierarchies, which hypre's l1-Jacobi
// relaxation for all points divided by, 0 / 0, before
// https://github.com/hypre-space/hypre/pull/1676. Whether the row is exactly zero depends
// on rounding, so this only checks that the result is finite; "BoomerAMG l1-Jacobi with
// zero rows" in test-amg.cpp has zero rows on every platform.
TEST_CASE("AMS for a singular operator with a separate part", "[ams][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  auto box = mfem::Mesh::MakeCartesian3D(6, 6, 6, mfem::Element::TETRAHEDRON);
  auto column =
      mfem::Mesh::MakeCartesian3D(3, 3, 6, mfem::Element::TETRAHEDRON, 0.1, 0.1, 0.5);
  column.Transform(
      [](const mfem::Vector &x, mfem::Vector &y)
      {
        y = x;
        y(0) += 2.0;
      });
  mfem::Mesh *pieces[2] = {&box, &column};
  mfem::Mesh smesh(pieces, 2);
  Mesh mesh(comm, smesh);
  constexpr int order = 1;
  fem::DefaultIntegrationOrder::p_trial = order;
  mfem::ND_FECollection nd_fec(order, mesh.Dimension());
  mfem::H1_FECollection h1_fec(order, mesh.Dimension());
  FiniteElementSpace nd_fespace(mesh, &nd_fec), h1_fespace(mesh, &h1_fec);
  mfem::ConstantCoefficient one(1.0);
  mfem::ParBilinearForm a(&nd_fespace.Get());
  a.AddDomainIntegrator(new mfem::CurlCurlIntegrator(one));
  a.Assemble();
  a.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> A(a.ParallelAssemble());

  // A consistent right-hand side.
  const int n = A->Height();
  mfem::Vector u(n), b(n), x(n);
  u.Randomize(1 + Mpi::Rank(comm));
  A->Mult(u, b);

  // With the coarse levels of the hierarchies kept distributed, and agglomerated (with the
  // hypre patch for it in extern/patch/hypre; otherwise the same).
  for (int target : {0, 64, 2048})
  {
    HypreAmsSolver ams(nd_fespace, h1_fespace, /*cycle_it*/ 1, /*smooth_it*/ 1,
                       /*vector_interp*/ false, /*singular_op*/ true, /*agg_coarsen*/ true,
                       target, /*print*/ 0);
    ams.SetOperator(*A);
    x = 0.0;
    ams.Mult(b, x);
    int nonfinite = 0;
    for (int i = 0; i < n; i++)
    {
      nonfinite += !std::isfinite(x(i));
    }
    Mpi::GlobalSum(1, &nonfinite, comm);
    CAPTURE(target);
    CHECK(nonfinite == 0);
  }
}

}  // namespace palace
