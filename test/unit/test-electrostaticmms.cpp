// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Method of Manufactured Solutions (MMS) verification for the electrostatic operator.
//
// We manufacture a potential V_mms on the unit cube [0,1]³, substitute it into
//     -∇·(ε ∇V) = ρ   (constant ε)
// to get the source ρ_mms = -ε ∇²V_mms that makes V_mms exact, solve with V = V_mms on ∂Ω,
// and compare the discrete solution to V_mms in the L2 norm. Several manufactured solutions
// (defined below) exercise homogeneous and non-homogeneous Dirichlet BCs, checked both for
// optimal convergence rate under mesh refinement and for exactness when V_mms lies in the
// FE space.

#include <cmath>
#include <memory>
#include <string>
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

// sin(πx)sin(πy)sin(πz): vanishes on ∂Ω (homogeneous Dirichlet). Laplacian eigenfunction,
// so ρ_mms = -ε ∇²V_mms = 3π²ε V_mms.
double VmmsSin(const mfem::Vector &x)
{
  return std::sin(M_PI * x[0]) * std::sin(M_PI * x[1]) * std::sin(M_PI * x[2]);
}
double RhoSin(const mfem::Vector &x)
{
  return 3.0 * M_PI * M_PI * kEpsilonR * VmmsSin(x);
}

// cos(πx)cos(πy)cos(πz): nonzero on ∂Ω (non-homogeneous Dirichlet). Also an eigenfunction.
double VmmsCos(const mfem::Vector &x)
{
  return std::cos(M_PI * x[0]) * std::cos(M_PI * x[1]) * std::cos(M_PI * x[2]);
}
double RhoCos(const mfem::Vector &x)
{
  return 3.0 * M_PI * M_PI * kEpsilonR * VmmsCos(x);
}

// x² + y² + z²: a degree-2 polynomial, nonzero on ∂Ω. ∇²V_mms = 6, so ρ_mms = -6ε
// (constant).
double VmmsPoly(const mfem::Vector &x)
{
  return x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
}
double RhoPoly(const mfem::Vector &)
{
  return -6.0 * kEpsilonR;
}

// A manufactured case: exact solution V_mms, its source ρ_mms = -ε ∇²V_mms, and whether
// V_mms is nonzero on the boundary (i.e. needs a non-homogeneous Dirichlet lift).
struct MmsCase
{
  double (*v_mms)(const mfem::Vector &);
  double (*rho_mms)(const mfem::Vector &);
  bool nonzero_boundary;
};

const MmsCase kHomogeneous{VmmsSin, RhoSin, false};
const MmsCase kNonHomogeneous{VmmsCos, RhoCos, true};
const MmsCase kPolynomial{VmmsPoly, RhoPoly, true};

// Solve the manufactured electrostatic problem for the given case on an N×N×N unit cube and
// return the L2 error ‖V_h - V_mms‖ against the manufactured solution.
double SolveMmsL2Error(const MmsCase &mms, int n, int order, double linear_tol = 1.0e-9)
{
  MPI_Comm comm = Mpi::World();

  // Build the IoData in-memory (no mesh file): constant permittivity on the single domain
  // attribute, all six cube faces marked Dirichlet.
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
  // error falls below, which would collapse the measured convergence rate. The
  // polynomial-exactness test needs an even tighter tol, since its discretization error is
  // zero and the algebraic error is then all that remains.
  iodata.solver.linear.tol = linear_tol;

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

  // Set the manufactured source ρ_mms and, for the non-homogeneous case, the prescribed
  // boundary values V_mms. Then run the production solve path via the exposed inner Solve
  // (root=false so BaseSolver's ctor writes no metadata/output files).
  mfem::FunctionCoefficient v_exact(mms.v_mms);
  mfem::FunctionCoefficient rho_coef(mms.rho_mms);
  laplace_op.SetRhsSource(rho_coef);
  if (mms.nonzero_boundary)
  {
    laplace_op.SetDbcCoefficient(v_exact);
  }
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

// Single-resolution sanity check: both the homogeneous and non-homogeneous manufactured
// problems solve and reach a small L2 error.
TEST_CASE("Electrostatic MMS solves to a small L2 error",
          "[electrostaticmms][Serial][Parallel]")
{
  CHECK(SolveMmsL2Error(kHomogeneous, /*n=*/16, /*order=*/2) < 1.0e-3);
  CHECK(SolveMmsL2Error(kNonHomogeneous, /*n=*/16, /*order=*/2) < 1.0e-3);
}

// A polynomial of degree ≤ p lies exactly in the order-p H1 space, so the FEM has zero
// discretization error and should reproduce V_mms to ~machine precision on any mesh — a
// sharper check than the convergence rate (a wrong sign/constant gives an O(1) error, not a
// bad slope). With discretization error zero, the iterative-solver error is all that
// remains, so use a tight linear tolerance to keep it below the threshold.
TEST_CASE("Electrostatic MMS is exact for a polynomial in the FE space",
          "[electrostaticmms][Serial][Parallel]")
{
  CHECK(SolveMmsL2Error(kPolynomial, /*n=*/4, /*order=*/2, /*linear_tol=*/1.0e-13) <
        1.0e-10);
}

// Convergence-rate verification: for each polynomial order p, the L2 error must decrease
// monotonically under mesh refinement and every consecutive-level rate must match the
// theoretical p+1 (to within kRateTol). Checking every consecutive pair confirms the rate
// is stable across refinement, not just at the finest level. Order is parametrized so a
// failure reports which order regressed.
TEST_CASE("Electrostatic MMS converges at the optimal rate",
          "[electrostaticmms][Serial][Parallel][Long]")
{
  // Run the sweep for both the homogeneous and non-homogeneous manufactured solutions, at
  // each polynomial order. GENERATE takes the Cartesian product, so each (case, order) is a
  // separate reported section — a failure identifies exactly which combination regressed.
  const auto [name, mms] = GENERATE(table<std::string, MmsCase>(
      {{"homogeneous", kHomogeneous}, {"non-homogeneous", kNonHomogeneous}}));
  const int order = GENERATE(from_range(kOrders));
  CAPTURE(name, order);

  // Solve on each refinement level and record (h, error).
  std::vector<double> h(kResolutions.size()), err(kResolutions.size());
  for (std::size_t i = 0; i < kResolutions.size(); i++)
  {
    h[i] = 1.0 / kResolutions[i];
    err[i] = SolveMmsL2Error(mms, kResolutions[i], order);
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
