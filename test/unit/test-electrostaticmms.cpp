// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Method of Manufactured Solutions (MMS) verification for the electrostatic operator.
//
// We manufacture a smooth potential on the unit cube [0,1]³ that vanishes on the boundary:
//     V_mms(x,y,z) = sin(πx) sin(πy) sin(πz)
// Since ∇²(sin πx sin πy sin πz) = -3π² V_mms, the volumetric charge that makes V_mms an
// exact solution of  -∇·(ε ∇V) = ρ  (constant ε) is
//     ρ_mms = -ε ∇²V_mms = 3π²ε V_mms.
// We solve -∇·(ε ∇V) = ρ_mms with homogeneous Dirichlet BCs (all six faces grounded,
// matching V_mms = 0 on ∂Ω) and check that the discrete solution converges to V_mms in the
// L2 norm.

#include <cmath>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "drivers/electrostaticsolver.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/laplaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

using namespace palace;
using namespace Catch::Matchers;

namespace
{

// Constant relative permittivity on the domain. A single scalar ε is what lets the source
// be ρ_mms = -ε ∇²V_mms (ε pulls out of the divergence); a spatially-varying or tensor ε
// would add a ∇ε·∇V term.
constexpr double kEpsilonR = 1.0;

// Manufactured solution V_mms = sin(πx) sin(πy) sin(πz). Vanishes on all faces of [0,1]³
// (so it matches a grounded boundary) and is a Laplacian eigenfunction (∇²V_mms = -3π²
// V_mms), which makes the source a trivial multiple of the solution.
double VmmsFunction(const mfem::Vector &x)
{
  return std::sin(M_PI * x[0]) * std::sin(M_PI * x[1]) * std::sin(M_PI * x[2]);
}

// Manufactured source ρ_mms = -ε ∇²V_mms = 3π²ε · V_mms (the 3 is the dimension count).
double RhommsFunction(const mfem::Vector &x)
{
  return 3.0 * M_PI * M_PI * kEpsilonR * VmmsFunction(x);
}

// Solve the manufactured electrostatic problem on an N×N×N unit cube and return the L2
// error ‖V_h - V_mms‖ against the manufactured solution.
double SolveMmsL2Error(int n, int order)
{
  MPI_Comm comm = Mpi::World();

  // Build the IoData in-memory (no mesh file): constant permittivity on the single domain
  // attribute, all six cube faces grounded (Dirichlet).
  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.model.Lc = 1.0;
  iodata.problem.type = ProblemType::ELECTROSTATIC;
  auto &material = iodata.domains.materials.emplace_back();
  material.attributes = {1};
  material.epsilon_r.s = {kEpsilonR, kEpsilonR, kEpsilonR};
  iodata.boundaries.pec.attributes = {1, 2, 3, 4, 5, 6};
  iodata.solver.order = order;
  // Drive the iterative-solver (algebraic) error below the discretization error we measure.
  // The default 1e-6 leaves an error floor that the finest/highest-order discretization
  // error falls below, which would collapse the measured convergence rate.
  iodata.solver.linear.tol = 1.0e-9;

  // Structured unit-cube mesh; MakeCartesian3D assigns domain attribute 1 and face
  // attributes 1..6, so h = 1/n and refining is just changing n.
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(n, n, n, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0));
  iodata.NondimensionalizeInputs(serial_mesh);
  iodata.CheckConfiguration();  // also initializes DefaultIntegrationOrder used by
                                // integrators
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(std::move(par_mesh)));

  LaplaceOperator laplace_op(iodata, mesh);
  auto &h1_fespace = laplace_op.GetH1Space();

  // Set the manufactured source ρ_mms on the operator, then run the production solve path
  // via the exposed inner Solve: the operator assembles ∫ρφ into the RHS and applies
  // homogeneous Dirichlet BCs internally. root=false so BaseSolver's ctor writes no
  // metadata/output files.
  mfem::FunctionCoefficient rho_coef(RhommsFunction);
  laplace_op.SetRhsSource(rho_coef);
  ExposedElectrostaticSolver solver(iodata, /*root=*/false);
  std::vector<Vector> V;
  solver.Solve(V, laplace_op);

  // Recover the grid function and compute the L2 error against V_mms. Integrate two orders
  // above the FE order so quadrature error doesn't contaminate the measured discretization
  // error.
  mfem::ParGridFunction V_gf(&h1_fespace.Get());
  V_gf.SetFromTrueDofs(V[0]);

  const int order_quad = order + 2;
  const mfem::IntegrationRule *irs[mfem::Geometry::NumGeom];
  for (int i = 0; i < mfem::Geometry::NumGeom; i++)
  {
    irs[i] = &mfem::IntRules.Get(i, order_quad);
  }
  mfem::FunctionCoefficient v_exact(VmmsFunction);
  return V_gf.ComputeL2Error(v_exact, irs);
}

// Observed convergence rate between two (mesh size h, L2 error) points: the slope of the
// error in log-log space, log(e0/e1) / log(h0/h1).
double ConvergenceRate(double h0, double e0, double h1, double e1)
{
  return std::log(e0 / e1) / std::log(h0 / h1);
}

// Convergence sweep configuration. Extend these lists to add polynomial orders or
// refinement levels; the test logic below adapts to any number of entries.
const std::vector<int> kOrders = {1, 2, 3};
const std::vector<int> kResolutions = {8, 16,
                                       32};  // N per side; h = 1/N (must be ascending)

// Tolerance on the measured convergence rate vs. the theoretical p+1.
constexpr double kRateTol = 0.3;

}  // namespace

// Single-resolution sanity check: the manufactured problem solves and the L2 error is
// small.
TEST_CASE("Electrostatic MMS solves to a small L2 error",
          "[electrostaticmms][Serial][Parallel]")
{
  const double err = SolveMmsL2Error(/*n=*/16, /*order=*/2);
  CHECK(err < 1.0e-3);
}

// Convergence-rate verification: for each polynomial order p, the L2 error must decrease
// monotonically under mesh refinement and every consecutive-level rate must match the
// theoretical p+1 (to within kRateTol). Checking every consecutive pair confirms the rate
// is stable across refinement, not just at the finest level. Order is parametrized so a
// failure reports which order regressed.
TEST_CASE("Electrostatic MMS converges at the optimal rate",
          "[electrostaticmms][Serial][Parallel][Long]")
{
  const int order = GENERATE(from_range(kOrders));
  CAPTURE(order);

  // Solve on each refinement level and record (h, error).
  std::vector<double> h(kResolutions.size()), err(kResolutions.size());
  for (std::size_t i = 0; i < kResolutions.size(); i++)
  {
    h[i] = 1.0 / kResolutions[i];
    err[i] = SolveMmsL2Error(kResolutions[i], order);
  }

  // Every consecutive pair: error decreases and the observed rate matches p+1.
  const double expected_rate = order + 1.0;
  for (std::size_t i = 0; i + 1 < kResolutions.size(); i++)
  {
    CAPTURE(kResolutions[i], kResolutions[i + 1], err[i], err[i + 1]);
    CHECK(err[i + 1] < err[i]);
    const double rate = ConvergenceRate(h[i], err[i], h[i + 1], err[i + 1]);
    CAPTURE(rate, expected_rate);
    CHECK_THAT(rate, WithinAbs(expected_rate, kRateTol));
  }
}
