// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <utility>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
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

// hypre AMS with the options of palace::HypreAmsSolver (cycle 14, l1-Jacobi BoomerAMG with
// Gaussian elimination on the coarsest level).
struct AmsSolver
{
  HYPRE_Solver ams = nullptr;
  AmsSolver(const mfem::HypreParMatrix &A, const mfem::HypreParMatrix &G,
            mfem::HypreParVector &x, mfem::HypreParVector &y, mfem::HypreParVector &z,
            bool singular, int maxit, int sweeps)
  {
    HYPRE_AMSCreate(&ams);
    HYPRE_AMSSetDimension(ams, 3);
    HYPRE_AMSSetCycleType(ams, 14);
    HYPRE_AMSSetPrintLevel(ams, 0);
    HYPRE_AMSSetMaxIter(ams, maxit);
    // A fixed number of iterations, as for a preconditioner: no convergence error is
    // flagged in hypre's global error state (which later solver setups would report).
    HYPRE_AMSSetTol(ams, 0.0);
    if (singular)
    {
      HYPRE_AMSSetBetaPoissonMatrix(ams, nullptr);
    }
    HYPRE_AMSSetSmoothingOptions(ams, 2, sweeps, 1.0, 1.0);
    HYPRE_AMSSetAlphaAMGOptions(ams, 10, singular ? 1 : 0, 18, 0.5, 6, 4);
    HYPRE_AMSSetBetaAMGOptions(ams, 10, singular ? 1 : 0, 18, 0.5, 6, 4);
    HYPRE_AMSSetAlphaAMGCoarseRelaxType(ams, 9);
    HYPRE_AMSSetBetaAMGCoarseRelaxType(ams, 9);
    HYPRE_AMSSetDiscreteGradient(ams, G);
    HYPRE_AMSSetCoordinateVectors(ams, x, y, z);
    mfem::HypreParVector b(A), u(A);
    HYPRE_AMSSetup(ams, A, b, u);
  }
  ~AmsSolver() { HYPRE_AMSDestroy(ams); }
  void Mult(const mfem::HypreParMatrix &A, mfem::HypreParVector &b, mfem::HypreParVector &u,
            bool flag_zero)
  {
    if (flag_zero)
    {
      hypre_ParVectorSetZeros(u);
    }
    else
    {
      u = 0.0;
    }
    HYPRE_AMSSolve(ams, A, b, u);
  }
};

}  // namespace

TEST_CASE("AMS with agglomerated auxiliary space coarse levels", "[ams][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  auto smesh = mfem::Mesh::MakeCartesian3D(12, 10, 8, mfem::Element::TETRAHEDRON);
  mfem::ParMesh pmesh(comm, smesh);
  mfem::ND_FECollection nd_fec(1, 3);
  mfem::H1_FECollection h1_fec(1, 3);
  mfem::ParFiniteElementSpace nd(&pmesh, &nd_fec), h1(&pmesh, &h1_fec);
  mfem::Array<int> dbc_marker(pmesh.bdr_attributes.Max()), dbc_tdof_list;
  dbc_marker = 1;  // PEC enclosure: well-posed also without the mass term
  nd.GetEssentialTrueDofs(dbc_marker, dbc_tdof_list);

  mfem::ParDiscreteLinearOperator grad(&h1, &nd);
  grad.AddDomainInterpolator(new mfem::GradientInterpolator);
  grad.Assemble();
  grad.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> G(grad.ParallelAssemble());

  mfem::ParGridFunction xc(&h1), yc(&h1), zc(&h1);
  mfem::FunctionCoefficient fx([](const mfem::Vector &p) { return p(0); }),
      fy([](const mfem::Vector &p) { return p(1); }),
      fz([](const mfem::Vector &p) { return p(2); });
  xc.ProjectCoefficient(fx);
  yc.ProjectCoefficient(fy);
  zc.ProjectCoefficient(fz);
  std::unique_ptr<mfem::HypreParVector> X(xc.ParallelProject()), Y(yc.ParallelProject()),
      Z(zc.ParallelProject());

  // One and several Nedelec smoothing sweeps per AMS iteration (more than one sweep checks
  // that the zero flag of the initial guess is cleared after the first sweep).
  for (auto [singular, sweeps] :
       {std::pair{true, 1}, std::pair{false, 1}, std::pair{true, 3}, std::pair{false, 3}})
  {
    mfem::ConstantCoefficient one(1.0);
    mfem::ParBilinearForm a(&nd);
    a.AddDomainIntegrator(new mfem::CurlCurlIntegrator(one));
    if (!singular)
    {
      a.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(one));
    }
    a.Assemble();
    a.Finalize();
    std::unique_ptr<mfem::HypreParMatrix> A(a.ParallelAssemble());
    std::unique_ptr<mfem::HypreParMatrix> Ae(A->EliminateRowsCols(dbc_tdof_list));

    mfem::HypreParVector b(*A), u_ref(*A), u(*A), u0(*A);
    b.Randomize(1 + Mpi::Rank(comm));

    // Reference: hierarchies kept distributed during the solve.
    setenv("HYPRE_AMG_REDIST", "0", 1);
    {
      AmsSolver ams(*A, *G, *X, *Y, *Z, singular, 2, sweeps);
      ams.Mult(*A, b, u_ref, false);
      // A zero initial guess flagged to hypre gives the same result.
      ams.Mult(*A, b, u0, true);
      CHECK(RelMaxDiff(comm, u0, u_ref) == 0.0);
    }

    // Coarse levels agglomerated (and concurrent Pi-space corrections), several targets.
    for (const char *target : {"2048", "16", "4"})
    {
      setenv("HYPRE_AMG_REDIST", target, 1);
      AmsSolver ams(*A, *G, *X, *Y, *Z, singular, 2, sweeps);
      for (int it = 0; it < 2; it++)
      {
        ams.Mult(*A, b, u, it == 0);
        CHECK(RelMaxDiff(comm, u, u_ref) <= 1.0e-10);
      }
    }
    unsetenv("HYPRE_AMG_REDIST");
  }
}

}  // namespace palace
