// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Method of Manufactured Solutions (MMS) verification for the electrostatic operator.
//
// We manufacture potentials on affine unit-cube meshes and a curved cylindrical mesh,
// substitute them into
//     -∇·(ε ∇V) = ρ   (constant diagonal tensor ε)
// to get the source ρ_mms that makes V_mms exact, solve with V = V_mms on ∂Ω, and compare
// both the discrete potential and recovered electric field E = -∇V to their manufactured
// counterparts in the L2 norm. The smooth convergence cases use an isotropic material
// baseline, while polynomial exactness and Neumann cases exercise anisotropic permittivity
// as a Palace-specific extension.

#include <array>
#include <cmath>
#include <fstream>
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

// The smooth cases use the isotropic material baseline assumed by the paper. Focused
// Palace-specific cases retain distinct principal values to detect omitted, scalarized, or
// permuted material components.
constexpr std::array<double, 3> kIsotropicEpsilonR = {1.0, 1.0, 1.0};
constexpr std::array<double, 3> kAnisotropicEpsilonR = {2.0, 3.0, 5.0};

// sin(πx)sin(2πy)sin(πz): vanishes on ∂Ω (homogeneous Dirichlet). For ε = I,
// ρ_mms = 6π²V_mms.
double VmmsSin(const mfem::Vector &x)
{
  return std::sin(M_PI * x[0]) * std::sin(2.0 * M_PI * x[1]) * std::sin(M_PI * x[2]);
}
void EmmsSin(const mfem::Vector &x, mfem::Vector &E)
{
  E[0] =
      -M_PI * std::cos(M_PI * x[0]) * std::sin(2.0 * M_PI * x[1]) * std::sin(M_PI * x[2]);
  E[1] = -2.0 * M_PI * std::sin(M_PI * x[0]) * std::cos(2.0 * M_PI * x[1]) *
         std::sin(M_PI * x[2]);
  E[2] =
      -M_PI * std::sin(M_PI * x[0]) * std::sin(2.0 * M_PI * x[1]) * std::cos(M_PI * x[2]);
}
double RhoSin(const mfem::Vector &x)
{
  return M_PI * M_PI *
         (kIsotropicEpsilonR[0] + 4.0 * kIsotropicEpsilonR[1] + kIsotropicEpsilonR[2]) *
         VmmsSin(x);
}

// cos(πx)cos(2πy)cos(πz): nonzero on ∂Ω (non-homogeneous Dirichlet), with the
// same isotropic source factor as VmmsSin.
double VmmsCos(const mfem::Vector &x)
{
  return std::cos(M_PI * x[0]) * std::cos(2.0 * M_PI * x[1]) * std::cos(M_PI * x[2]);
}
void EmmsCos(const mfem::Vector &x, mfem::Vector &E)
{
  E[0] = M_PI * std::sin(M_PI * x[0]) * std::cos(2.0 * M_PI * x[1]) * std::cos(M_PI * x[2]);
  E[1] = 2.0 * M_PI * std::cos(M_PI * x[0]) * std::sin(2.0 * M_PI * x[1]) *
         std::cos(M_PI * x[2]);
  E[2] = M_PI * std::cos(M_PI * x[0]) * std::cos(2.0 * M_PI * x[1]) * std::sin(M_PI * x[2]);
}
double RhoCos(const mfem::Vector &x)
{
  return M_PI * M_PI *
         (kIsotropicEpsilonR[0] + 4.0 * kIsotropicEpsilonR[1] + kIsotropicEpsilonR[2]) *
         VmmsCos(x);
}

// x² + 2y² + 3z²: a degree-2 polynomial, nonzero on ∂Ω. Thus
// ρ_mms = -2(ε_x + 2ε_y + 3ε_z), and ∇V_mms = (2x, 4y, 6z).
double VmmsPoly(const mfem::Vector &x)
{
  return x[0] * x[0] + 2.0 * x[1] * x[1] + 3.0 * x[2] * x[2];
}
void EmmsPoly(const mfem::Vector &x, mfem::Vector &E)
{
  E[0] = -2.0 * x[0];
  E[1] = -4.0 * x[1];
  E[2] = -6.0 * x[2];
}
double RhoPoly(const mfem::Vector &)
{
  return -2.0 * (kAnisotropicEpsilonR[0] + 2.0 * kAnisotropicEpsilonR[1] +
                 3.0 * kAnisotropicEpsilonR[2]);
}
// Isotropic counterpart used on the curved cylinder to isolate the isoparametric mapping.
double RhoPolyIsotropic(const mfem::Vector &)
{
  return -2.0 * (kIsotropicEpsilonR[0] + 2.0 * kIsotropicEpsilonR[1] +
                 3.0 * kIsotropicEpsilonR[2]);
}
// Neumann flux g = n̂·ε∇V for the polynomial on the x = 1 face (outward normal +x̂):
// g = 2ε_x x, which equals 2ε_x there.
double NeumannPolyX1(const mfem::Vector &x)
{
  return 2.0 * kAnisotropicEpsilonR[0] * x[0];
}

// A manufactured case: exact solution V_mms, its source ρ_mms = -∇·(ε∇V_mms), constant
// diagonal relative permittivity, and whether V_mms needs a non-homogeneous Dirichlet lift.
struct MmsCase
{
  double (*v_mms)(const mfem::Vector &);
  void (*e_mms)(const mfem::Vector &, mfem::Vector &);
  double (*rho_mms)(const mfem::Vector &);
  std::array<double, 3> epsilon_r;
  bool nonzero_boundary;
};

const MmsCase kHomogeneous{VmmsSin, EmmsSin, RhoSin, kIsotropicEpsilonR, false};
const MmsCase kNonHomogeneous{VmmsCos, EmmsCos, RhoCos, kIsotropicEpsilonR, true};
const MmsCase kCurvedCylinder{VmmsPoly, EmmsPoly, RhoPolyIsotropic, kIsotropicEpsilonR,
                              true};
const MmsCase kPolynomial{VmmsPoly, EmmsPoly, RhoPoly, kAnisotropicEpsilonR, true};

// Solve the manufactured electrostatic problem for the given case on an N×N×N unit cube and
// return the L2 errors in the potential and electric field against the manufactured
// solution.
struct MmsError
{
  double potential;
  double electric;
};

MmsError ComputeMmsErrors(LaplaceOperator &laplace_op, const Vector &V,
                          mfem::Coefficient &v_exact, mfem::VectorCoefficient &e_exact)
{
  mfem::ParGridFunction V_gf(&laplace_op.GetH1Space().Get());
  V_gf.SetFromTrueDofs(V);

  Vector E(laplace_op.GetGradMatrix().Height());
  E = 0.0;
  laplace_op.GetGradMatrix().AddMult(V, E, -1.0);
  mfem::ParGridFunction E_gf(&laplace_op.GetNDSpace().Get());
  E_gf.SetFromTrueDofs(E);

  return {V_gf.ComputeL2Error(v_exact), E_gf.ComputeL2Error(e_exact)};
}

MmsError SolveMmsError(const MmsCase &mms, std::unique_ptr<mfem::Mesh> serial_mesh,
                       int order, double linear_tol = 1.0e-9)
{
  MPI_Comm comm = Mpi::World();

  // Build the IoData in-memory with Dirichlet conditions on every mesh boundary.
  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.model.Lc = 1.0;
  iodata.problem.type = ProblemType::ELECTROSTATIC;
  auto &material = iodata.domains.materials.emplace_back();
  material.attributes = {1};
  material.epsilon_r.s = mms.epsilon_r;
  for (int i = 0; i < serial_mesh->bdr_attributes.Size(); i++)
  {
    iodata.boundaries.pec.attributes.push_back(serial_mesh->bdr_attributes[i]);
  }
  iodata.solver.order = order;
  // Drive the iterative-solver (algebraic) error below the discretization error we measure.
  // The default 1e-6 leaves an error floor that the finest/highest-order discretization
  // error falls below, which would collapse the measured convergence rate. The
  // polynomial-exactness test needs an even tighter tol, since its discretization error is
  // zero and the algebraic error is then all that remains.
  iodata.solver.linear.tol = linear_tol;

  iodata.NondimensionalizeInputs(serial_mesh);
  iodata.CheckConfiguration();  // also initializes DefaultIntegrationOrder used by
                                // integrators
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(std::move(par_mesh)));

  LaplaceOperator laplace_op(iodata, mesh);

  // Set the manufactured source ρ_mms and, for the non-homogeneous case, the prescribed
  // boundary values V_mms. Then run the production solve path via the exposed inner Solve
  // (root=false so BaseSolver's ctor writes no metadata/output files).
  mfem::FunctionCoefficient v_exact(mms.v_mms);
  mfem::VectorFunctionCoefficient e_exact(/*vdim=*/3, mms.e_mms);
  mfem::FunctionCoefficient rho_coef(mms.rho_mms);
  laplace_op.SetRhsSource(rho_coef);
  if (mms.nonzero_boundary)
  {
    laplace_op.SetDbcCoefficient(v_exact);
  }
  ExposedElectrostaticSolver solver(iodata, /*root=*/false);
  std::vector<Vector> V;
  solver.Solve(V, laplace_op);

  // Compute errors in V and in the production electric-field recovery E = -∇V. MFEM selects
  // an error-integration rule based on each finite-element space's polynomial order.
  return ComputeMmsErrors(laplace_op, V[0], v_exact, e_exact);
}

MmsError
SolveCartesianMmsError(const MmsCase &mms, int n, int order, double linear_tol = 1.0e-9,
                       mfem::Element::Type element_type = mfem::Element::HEXAHEDRON)
{
  // MakeCartesian3D assigns domain attribute 1 and face attributes 1..6. The default
  // hexahedral convergence meshes have h = 1/n.
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(n, n, n, element_type, 1.0, 1.0, 1.0));
  return SolveMmsError(mms, std::move(serial_mesh), order, linear_tol);
}

std::unique_ptr<mfem::Mesh> LoadCurvedCylinderMesh(int refinement)
{
  const std::string mesh_path =
      std::string(PALACE_TEST_DATA_DIR) + "/mesh/cylinder-tet-p2.msh";
  std::ifstream mesh_file(mesh_path);
  REQUIRE(mesh_file.good());
  constexpr bool generate_edges = false, refine = true, fix_orientation = true;
  auto serial_mesh =
      std::make_unique<mfem::Mesh>(mesh_file, generate_edges, refine, fix_orientation);
  REQUIRE(serial_mesh->GetNodes());
  REQUIRE(serial_mesh->GetNodes()->FESpace()->GetMaxElementOrder() == 2);
  for (int l = 0; l < refinement; l++)
  {
    serial_mesh->UniformRefinement();
  }
  return serial_mesh;
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
const std::vector<int> kResolutions = {4, 8,
                                       16};  // N per side; h = 1/N (must be ascending)

// Tolerance on the measured convergence rates vs. the theoretical p+1 for V and p for E.
constexpr double kRateTol = 0.3;
// The coarsest curved cylinder mesh is intentionally nonuniform, so its asymptotic rate has
// slightly more variation than the structured Cartesian meshes.
constexpr double kCurvedRateTol = 0.35;

}  // namespace

// Single-resolution sanity check: both the homogeneous and non-homogeneous manufactured
// problems solve and reach small potential and electric-field L2 errors.
TEST_CASE("Electrostatic MMS solves to small potential and field errors",
          "[electrostaticmms][Serial][Parallel]")
{
  const auto homogeneous = SolveCartesianMmsError(kHomogeneous, /*n=*/16, /*order=*/2);
  CHECK(homogeneous.potential < 1.0e-3);
  CHECK(homogeneous.electric < 2.0e-2);

  const auto nonhomogeneous =
      SolveCartesianMmsError(kNonHomogeneous, /*n=*/16, /*order=*/2);
  CHECK(nonhomogeneous.potential < 1.0e-3);
  CHECK(nonhomogeneous.electric < 2.0e-2);
}

// On these affine Cartesian meshes, a polynomial of degree ≤ p lies exactly in the order-p
// H1 space, so the FEM has zero discretization error and should reproduce V_mms to ~machine
// precision. The straight-tetrahedron case mirrors the paper's polynomial exactness check.
// With discretization error zero, the iterative-solver error is all that remains, so use a
// tight linear tolerance to keep it below the threshold.
TEST_CASE("Electrostatic MMS is exact for a polynomial in the FE space",
          "[electrostaticmms][Serial][Parallel]")
{
  const auto [element_name, element_type] =
      GENERATE(table<std::string, mfem::Element::Type>(
          {{"hexahedron", mfem::Element::HEXAHEDRON},
           {"tetrahedron", mfem::Element::TETRAHEDRON}}));
  CAPTURE(element_name);

  const auto error = SolveCartesianMmsError(kPolynomial, /*n=*/4, /*order=*/2,
                                            /*linear_tol=*/1.0e-13, element_type);
  CHECK(error.potential < 1.0e-10);
  CHECK(error.electric < 1.0e-10);
}

// Mixed Dirichlet/Neumann: leave the x = 1 face (attribute 3) out of the grounded set so it
// becomes a natural boundary carrying the prescribed flux g = n̂·ε∇V, and impose V_mms on
// the other five faces. Verifies the Neumann boundary term ∮ g φ dS. The polynomial V_mms
// is FE-exact at order 2, so a correct flux assembly gives ~machine-precision error; a
// missing or wrong Neumann term leaves an O(1) error.
TEST_CASE("Electrostatic MMS handles a Neumann boundary",
          "[electrostaticmms][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.model.Lc = 1.0;
  iodata.problem.type = ProblemType::ELECTROSTATIC;
  auto &material = iodata.domains.materials.emplace_back();
  material.attributes = {1};
  material.epsilon_r.s = kAnisotropicEpsilonR;
  iodata.boundaries.pec.attributes = {1, 2, 4, 5,
                                      6};  // all faces except x = 1 (attribute 3)
  iodata.solver.order = 2;
  iodata.solver.linear.tol = 1.0e-13;

  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0));
  iodata.NondimensionalizeInputs(serial_mesh);
  iodata.CheckConfiguration();
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(std::move(par_mesh)));

  LaplaceOperator laplace_op(iodata, mesh);
  mfem::FunctionCoefficient v_exact(VmmsPoly);
  mfem::VectorFunctionCoefficient e_exact(/*vdim=*/3, EmmsPoly);
  mfem::FunctionCoefficient rho_coef(RhoPoly);
  mfem::FunctionCoefficient neumann_coef(NeumannPolyX1);
  laplace_op.SetRhsSource(rho_coef);
  laplace_op.SetDbcCoefficient(v_exact);
  laplace_op.SetNeumannCoefficient(neumann_coef, {3});

  ExposedElectrostaticSolver solver(iodata, /*root=*/false);
  std::vector<Vector> V;
  solver.Solve(V, laplace_op);

  const auto error = ComputeMmsErrors(laplace_op, V[0], v_exact, e_exact);
  CHECK(error.potential < 1.0e-10);
  CHECK(error.electric < 1.0e-10);
}

// Convergence-rate verification: for each polynomial order p, the potential and electric
// field L2 errors must decrease monotonically under mesh refinement, with consecutive-level
// rates matching the theoretical p+1 for V and p for E (to within kRateTol). Checking every
// consecutive pair confirms the rates are stable across refinement, not just at the finest
// level. Order is parametrized so a failure reports which order regressed.
TEST_CASE("Electrostatic MMS converges at the optimal rate",
          "[electrostaticmms][Serial][Parallel]")
{
  // Run the sweep for both the homogeneous and non-homogeneous manufactured solutions, at
  // each polynomial order. GENERATE takes the Cartesian product, so each (case, order) is a
  // separate reported section — a failure identifies exactly which combination regressed.
  const auto [name, mms] = GENERATE(table<std::string, MmsCase>(
      {{"homogeneous", kHomogeneous}, {"non-homogeneous", kNonHomogeneous}}));
  const int order = GENERATE(from_range(kOrders));
  CAPTURE(name, order);

  // Solve on each refinement level and record (h, error).
  std::vector<double> h(kResolutions.size());
  std::vector<MmsError> error(kResolutions.size());
  for (std::size_t i = 0; i < kResolutions.size(); i++)
  {
    h[i] = 1.0 / kResolutions[i];
    error[i] = SolveCartesianMmsError(mms, kResolutions[i], order);
  }

  // Every consecutive pair: both errors decrease and the observed rates match theory.
  for (std::size_t i = 0; i + 1 < kResolutions.size(); i++)
  {
    CAPTURE(kResolutions[i], kResolutions[i + 1], error[i].potential,
            error[i + 1].potential, error[i].electric, error[i + 1].electric);

    CHECK(error[i + 1].potential < error[i].potential);
    const double potential_rate =
        ConvergenceRate(h[i], error[i].potential, h[i + 1], error[i + 1].potential);
    CAPTURE(potential_rate);
    CHECK_THAT(potential_rate, WithinAbs(order + 1.0, kRateTol));

    CHECK(error[i + 1].electric < error[i].electric);
    const double electric_rate =
        ConvergenceRate(h[i], error[i].electric, h[i + 1], error[i + 1].electric);
    CAPTURE(electric_rate);
    CHECK_THAT(electric_rate, WithinAbs(static_cast<double>(order), kRateTol));
  }
}

// Reuse the polynomial that is represented exactly on affine order-2 tetrahedra. The
// example cylinder instead uses quadratic isoparametric tetrahedra, so the physical-space
// basis functions are not polynomial and the error is nonzero. Uniformly refining the
// quadratic map should recover the optimal rates, as in the paper's curved-cylinder test.
// This exercises curved coordinate transformations; it does not regenerate CAD-projected
// meshes at each level.
TEST_CASE("Electrostatic MMS converges on curved tetrahedra",
          "[electrostaticmms][Serial][Parallel]")
{
  constexpr int order = 2;
  const std::vector<int> refinement = {0, 1, 2};
  std::vector<MmsError> error(refinement.size());
  for (std::size_t i = 0; i < refinement.size(); i++)
  {
    error[i] = SolveMmsError(kCurvedCylinder, LoadCurvedCylinderMesh(refinement[i]), order);
  }

  for (std::size_t i = 0; i + 1 < refinement.size(); i++)
  {
    CAPTURE(refinement[i], refinement[i + 1], error[i].potential, error[i + 1].potential,
            error[i].electric, error[i + 1].electric);

    CHECK(error[i + 1].potential < error[i].potential);
    const double potential_rate =
        ConvergenceRate(1.0, error[i].potential, 0.5, error[i + 1].potential);
    CAPTURE(potential_rate);
    CHECK_THAT(potential_rate, WithinAbs(order + 1.0, kCurvedRateTol));

    CHECK(error[i + 1].electric < error[i].electric);
    const double electric_rate =
        ConvergenceRate(1.0, error[i].electric, 0.5, error[i + 1].electric);
    CAPTURE(electric_rate);
    CHECK_THAT(electric_rate, WithinAbs(static_cast<double>(order), kCurvedRateTol));
  }
}
