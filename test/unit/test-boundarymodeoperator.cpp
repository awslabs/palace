// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <functional>
#include <limits>
#include <numbers>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "drivers/boundarymodesolver.hpp"

#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
#include "fixtures.hpp"
#include "models/boundarymodeoperator.hpp"
#include "models/farfieldboundaryoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/modeoperatorassembly.hpp"
#include "models/surfaceconductivityoperator.hpp"
#include "models/surfaceimpedanceoperator.hpp"
#include "models/surfacerationalimpedanceoperator.hpp"
#include "models/waveportoperator.hpp"
#include "utils/communication.hpp"
#include "utils/constants.hpp"
#include "utils/filesystem.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

namespace palace
{
using namespace Catch::Matchers;

namespace
{

// Solve for modes of a 2D rectangular waveguide cross-section using
// ModeEigenSolver. Returns eigenvalues as complex kn values.
// MakeCartesian2D boundary attributes: bottom=1, right=2, top=3, left=4.
struct ModeResult
{
  std::vector<std::complex<double>> kn;
  std::vector<std::complex<double>> final_kn;
  std::vector<double> reduced_backward_errors;
  std::vector<double> final_reduced_backward_errors;
  int num_converged;
  ModeEigenSolver::ReducedModelStats reduced_stats;
  ModeEigenSolver::ReducedModelStats first_target_stats;
  std::size_t reduced_basis_size = 0;
  std::size_t first_target_basis_size = 0;
  double reduced_tol = 0.0;
  int complex_exact_converged = -1;
  // Shift-and-invert target of the in-band solve that produced kn.
  double kn_target = 0.0;
};

ModeResult SolveRectangularModes(double width, double height, double freq_ghz,
                                 double epsilon_r, int order, int num_modes,
                                 const std::function<void(IoData &)> &configure_bcs,
                                 bool exercise_reduced_model = false,
                                 int reduced_evaluations = 1,
                                 std::size_t reduced_training_capacity = 16,
                                 std::vector<double> reduced_training_factors = {0.9, 1.1},
                                 double adaptive_tol = 1.0e-3, double eig_tol = 1.0e-8)
{
  MPI_Comm comm = Mpi::World();
  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.model.Lc = 1.0;

  auto &material = iodata.domains.materials.emplace_back();
  material.attributes = {1};
  material.epsilon_r.s = {epsilon_r, epsilon_r, epsilon_r};

  // Default: PEC on all boundaries.
  iodata.boundaries.pec.attributes = {1, 2, 3, 4};

  // Let the caller configure specific BCs.
  configure_bcs(iodata);

  iodata.solver.order = order;
  iodata.solver.boundary_mode.freq = freq_ghz;
  iodata.solver.boundary_mode.n = num_modes;
  iodata.solver.boundary_mode.tol = eig_tol;
  iodata.solver.linear.tol = eig_tol;
  iodata.solver.linear.max_it = 200;

  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(10, 5, mfem::Element::TRIANGLE, false, width, height));
  iodata.NondimensionalizeInputs(serial_mesh);
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  iodata.CheckConfiguration();
  Mesh palace_mesh(std::move(par_mesh));

  auto nd_fec = std::make_unique<mfem::ND_FECollection>(order, palace_mesh.Dimension());
  auto h1_fec = std::make_unique<mfem::H1_FECollection>(order, palace_mesh.Dimension());
  FiniteElementSpace nd_fespace(palace_mesh, nd_fec.get());
  FiniteElementSpace h1_fespace(palace_mesh, h1_fec.get());
  MaterialOperator mat_op(iodata, palace_mesh);

  SurfaceImpedanceOperator surf_z_op(iodata, mat_op, palace_mesh.Get());
  FarfieldBoundaryOperator farfield_op(iodata, mat_op, palace_mesh.Get());
  SurfaceConductivityOperator surf_sigma_op(iodata, mat_op, palace_mesh.Get());
  SurfaceRationalImpedanceOperator surf_rz_op(iodata, mat_op, palace_mesh.Get());

  mfem::Array<int> nd_dbc_tdof_list, h1_dbc_tdof_list;
  {
    const auto &pmesh = palace_mesh.Get();
    int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
    auto dbc_marker = mesh::AttrToMarker(bdr_attr_max, iodata.boundaries.pec.attributes);
    nd_fespace.Get().GetEssentialTrueDofs(dbc_marker, nd_dbc_tdof_list);
    h1_fespace.Get().GetEssentialTrueDofs(dbc_marker, h1_dbc_tdof_list);
  }

  int nd_size = nd_fespace.GetTrueVSize();
  mfem::Array<int> dbc_tdof_list;
  dbc_tdof_list.Append(nd_dbc_tdof_list);
  for (int i = 0; i < h1_dbc_tdof_list.Size(); i++)
  {
    dbc_tdof_list.Append(nd_size + h1_dbc_tdof_list[i]);
  }

  double omega = 2.0 * std::numbers::pi *
                 iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(freq_ghz);

  // ModeEigenSolver requires a positive Krylov subspace size (num_vec). Mirror the
  // formula used by IoData::CheckConfiguration for eigenmode.max_size.
  const int num_vec = std::max(2 * num_modes, num_modes + 15);
  ModeEigenSolver mode_solver(mat_op, nullptr, surf_z_op, farfield_op, surf_sigma_op,
                              surf_rz_op, nd_fespace, h1_fespace, dbc_tdof_list, num_modes,
                              num_vec, eig_tol, EigenvalueSolver::WhichType::LARGEST_REAL,
                              iodata.solver.linear, iodata.solver.boundary_mode.type, 0,
                              nd_fespace.GetComm());

  auto target_at = [&](std::complex<double> w)
  { return w.real() * std::sqrt(1.1 * mat_op.GetMaxMuEpsilon()); };
  auto solve_at = [&](std::complex<double> w)
  {
    const double kn_target = target_at(w);
    const double sigma = -kn_target * kn_target;
    return std::make_pair(mode_solver.Solve(w, sigma), sigma);
  };

  if (exercise_reduced_model)
  {
    mode_solver.ConfigureReducedModelTraining(reduced_training_capacity);
    for (double factor : reduced_training_factors)
    {
      solve_at(factor * omega);
    }
    mode_solver.EnableReducedModel(adaptive_tol);
  }

  auto result = solve_at(omega).first;
  ModeResult out;
  out.num_converged = result.num_converged;
  out.kn_target = target_at(omega);
  for (int i = 0; i < result.num_converged; i++)
  {
    // Capture the first in-band evaluation, which is the reduced result under test.
    out.kn.push_back(mode_solver.GetPropagationConstant(i));
    if (exercise_reduced_model)
    {
      out.reduced_backward_errors.push_back(
          mode_solver.GetError(i, EigenvalueSolver::ErrorType::BACKWARD));
    }
  }

  if (exercise_reduced_model)
  {
    out.first_target_stats = mode_solver.GetReducedModelStats();
    out.first_target_basis_size = mode_solver.GetReducedBasisSize();
    for (int i = 1; i < reduced_evaluations; i++)
    {
      solve_at(omega);
    }
    for (int i = 0; i < result.num_converged; i++)
    {
      out.final_kn.push_back(mode_solver.GetPropagationConstant(i));
      out.final_reduced_backward_errors.push_back(
          mode_solver.GetError(i, EigenvalueSolver::ErrorType::BACKWARD));
    }
    // Complex-frequency queries must bypass the real-axis reduced model. Issue one complex
    // query solely to verify stats after retaining the reduced real-frequency result above.
    out.complex_exact_converged =
        solve_at(std::complex<double>(omega, 1.0e-3 * omega)).first.num_converged;
  }
  out.reduced_stats = mode_solver.GetReducedModelStats();
  out.reduced_basis_size = mode_solver.GetReducedBasisSize();
  out.reduced_tol = mode_solver.GetReducedTolerance();
  return out;
}

// Exact mode index of the TM0 mode of a parallel plate guide with a PEC plate at y = 0, a
// vacuum gap of height h, and a London slab of thickness d and penetration depth lambda,
// whose back face is PEC or free (PMC), at free-space wavenumber k0 (all lengths in um).
// With p^2 = k0^2 (n^2 - 1) and q^2 = 1/lambda^2 + p^2, H_x = cosh(p y) in the gap and cosh
// or sinh of q (h + d - y) in the slab, and continuity of H_x and of (1/eps) dH_x/dy at y =
// h, where eps = 1 - 1/(k0 lambda)^2 is the effective permittivity of the slab, give p
// tanh(p h) = -(q/eps) coth(q d) (free back face) or -(q/eps) tanh(q d) (PEC back face). In
// the quasi-static limit, n^2 = 1 + (lambda/h) coth(d/lambda) or tanh(d/lambda) (Swihart).
double ExactLondonSlabIndex(double h, double d, double lambda, double k0, bool pec_back)
{
  const double eps = 1.0 - 1.0 / (k0 * k0 * lambda * lambda);
  auto F = [&](double n)
  {
    const double p = k0 * std::sqrt(n * n - 1.0);
    const double q = std::sqrt(1.0 / (lambda * lambda) + p * p);
    const double g = pec_back ? std::tanh(q * d) : 1.0 / std::tanh(q * d);
    return p * std::tanh(p * h) + q * g / eps;
  };
  double lo = 1.0 + 1.0e-12, hi = 10.0;
  for (int it = 0; it < 200; it++)
  {
    const double mid = 0.5 * (lo + hi);
    ((F(mid) < 0.0) ? lo : hi) = mid;
  }
  return 0.5 * (lo + hi);
}

// Propagation constant normalized by the free-space wavenumber (n_eff) of a mode of a 2D
// cross-section of width w, in um: a PEC plate at y = 0 (attribute 1), a vacuum gap of
// height h, and a London superconductor slab of thickness d with penetration depth lambda
// above it, whose back face y = h + d (attribute 2) is PEC or left free (PMC). The sides x
// = 0, w (attribute 3) are PEC or free.
double SolveLondonSlabMode(double w, double h, double d, double lambda, bool pec_back,
                           bool pec_sides, double freq_ghz, double n_target, int order,
                           int ny_gap = 8)
{
  MPI_Comm comm = Mpi::World();
  Units units(1.0e-6, 1.0e-6);
  IoData iodata(units);
  iodata.model.Lc = 1.0;
  auto &vacuum = iodata.domains.materials.emplace_back();
  vacuum.attributes = {1};
  auto &london = iodata.domains.materials.emplace_back();
  london.attributes = {2};
  london.lambda_L = lambda;
  iodata.boundaries.pec.attributes = {1};
  if (pec_back)
  {
    iodata.boundaries.pec.attributes.push_back(2);
  }
  if (pec_sides)
  {
    iodata.boundaries.pec.attributes.push_back(3);
  }
  iodata.solver.order = order;
  iodata.solver.boundary_mode.freq = freq_ghz;
  iodata.solver.boundary_mode.n = 1;
  iodata.solver.boundary_mode.tol = 1.0e-12;
  iodata.solver.linear.tol = 1.0e-12;
  iodata.solver.linear.max_it = 200;

  // Structured quadrilaterals: 2 across the width, ny_gap across the gap, and enough across
  // the slab to resolve the penetration depth.
  constexpr int nx = 2;
  const int ny_slab = std::max(8, static_cast<int>(std::ceil(8.0 * d / lambda)));
  std::vector<double> y;
  for (int j = 0; j <= ny_gap; j++)
  {
    y.push_back(h * j / ny_gap);
  }
  for (int j = 1; j <= ny_slab; j++)
  {
    y.push_back(h + d * j / ny_slab);
  }
  const int ny = static_cast<int>(y.size()) - 1;
  auto serial_mesh =
      std::make_unique<mfem::Mesh>(2, (nx + 1) * (ny + 1), nx * ny, 2 * (nx + ny));
  for (int j = 0; j <= ny; j++)
  {
    for (int i = 0; i <= nx; i++)
    {
      serial_mesh->AddVertex(w * i / nx, y[j]);
    }
  }
  auto v = [&](int i, int j) { return i + (nx + 1) * j; };
  for (int j = 0; j < ny; j++)
  {
    for (int i = 0; i < nx; i++)
    {
      serial_mesh->AddQuad(v(i, j), v(i + 1, j), v(i + 1, j + 1), v(i, j + 1),
                           (j < ny_gap) ? 1 : 2);
    }
  }
  for (int i = 0; i < nx; i++)
  {
    serial_mesh->AddBdrSegment(v(i, 0), v(i + 1, 0), 1);
    serial_mesh->AddBdrSegment(v(i + 1, ny), v(i, ny), 2);
  }
  for (int j = 0; j < ny; j++)
  {
    serial_mesh->AddBdrSegment(v(0, j + 1), v(0, j), 3);
    serial_mesh->AddBdrSegment(v(nx, j), v(nx, j + 1), 3);
  }
  serial_mesh->FinalizeQuadMesh(1, 1, true);

  iodata.NondimensionalizeInputs(serial_mesh);
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  iodata.CheckConfiguration();
  Mesh palace_mesh(std::move(par_mesh));

  auto nd_fec = std::make_unique<mfem::ND_FECollection>(order, palace_mesh.Dimension());
  auto h1_fec = std::make_unique<mfem::H1_FECollection>(order, palace_mesh.Dimension());
  FiniteElementSpace nd_fespace(palace_mesh, nd_fec.get());
  FiniteElementSpace h1_fespace(palace_mesh, h1_fec.get());
  MaterialOperator mat_op(iodata, palace_mesh);
  SurfaceImpedanceOperator surf_z_op(iodata, mat_op, palace_mesh.Get());
  FarfieldBoundaryOperator farfield_op(iodata, mat_op, palace_mesh.Get());
  SurfaceConductivityOperator surf_sigma_op(iodata, mat_op, palace_mesh.Get());
  SurfaceRationalImpedanceOperator surf_rz_op(iodata, mat_op, palace_mesh.Get());

  mfem::Array<int> nd_dbc_tdof_list, h1_dbc_tdof_list;
  {
    const auto &pmesh = palace_mesh.Get();
    int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
    auto dbc_marker = mesh::AttrToMarker(bdr_attr_max, iodata.boundaries.pec.attributes);
    nd_fespace.Get().GetEssentialTrueDofs(dbc_marker, nd_dbc_tdof_list);
    h1_fespace.Get().GetEssentialTrueDofs(dbc_marker, h1_dbc_tdof_list);
  }
  const int nd_size = nd_fespace.GetTrueVSize();
  mfem::Array<int> dbc_tdof_list;
  dbc_tdof_list.Append(nd_dbc_tdof_list);
  for (int i = 0; i < h1_dbc_tdof_list.Size(); i++)
  {
    dbc_tdof_list.Append(nd_size + h1_dbc_tdof_list[i]);
  }

  const double omega =
      2.0 * std::numbers::pi *
      iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(freq_ghz);
  // Shift-and-invert about the target, as in BoundaryModeSolver with a target n_eff.
  const int num_vec = 16;
  ModeEigenSolver mode_solver(mat_op, nullptr, surf_z_op, farfield_op, surf_sigma_op,
                              surf_rz_op, nd_fespace, h1_fespace, dbc_tdof_list, 1, num_vec,
                              1.0e-12, EigenvalueSolver::WhichType::LARGEST_MAGNITUDE,
                              iodata.solver.linear, iodata.solver.boundary_mode.type, 0,
                              nd_fespace.GetComm());
  const double kn_target = omega * n_target;
  auto result = mode_solver.Solve(omega, -kn_target * kn_target);
  REQUIRE(result.num_converged >= 1);
  const auto kn = mode_solver.GetPropagationConstant(0);
  CHECK(std::abs(kn.imag()) < 1.0e-6 * std::abs(kn.real()));
  return kn.real() / omega;
}

}  // namespace

TEST_CASE("ModeEigenSolver PEC", "[boundarymodeoperator][Serial]")
{
  // Rectangular waveguide: 1000×500 μm (L0=1e-6), ε=4, f=500 GHz.
  // Analytical kn for TE10 mode:
  //   kc = π / (a * L0) = π / 1e-3 ≈ 3141.6 1/m
  //   ω = 2π * 500e9 ≈ 3.1416e12 rad/s
  //   kn = sqrt(ω²ε/c² - kc²) = sqrt(4*(π*1e12/c)² - (π/1e-3)²)
  //   In nondimensional units (Lc = a*L0 = 1e-3 m):
  //     kn_nd = kn * Lc
  auto result = SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 3, [](IoData &) {});

  REQUIRE(result.num_converged >= 1);

  double kn_real = result.kn[0].real();
  double kn_imag = result.kn[0].imag();
  CAPTURE(kn_real, kn_imag, result.num_converged);

  // First mode should be propagating (real kn, negligible imaginary part).
  CHECK(kn_real > 0.0);
  CHECK(std::abs(kn_imag) < 1.0e-6 * std::abs(kn_real));

  // Analytical kn for TE10 mode of rectangular waveguide with PEC walls:
  //   a = 1000 μm = 1e-3 m, ε_r = 4, f = 500 GHz
  //   kc = π / a = π / 1e-3 m
  //   kn = sqrt(ω²ε_r/c² - kc²)
  //   kn ≈ 20708 1/m → nondimensional (×Lc where Lc = 1e-6 m) ≈ 0.02071
  // Allow 5% tolerance for the coarse 10×5 mesh at order 2.
  CHECK_THAT(kn_real, WithinRel(0.02071, 0.05));
}

TEST_CASE("ModeEigenSolver guarded reduced real-frequency solve",
          "[boundarymodeoperator][Serial][Parallel]")
{
  constexpr int num_modes = 3;
  auto exact =
      SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, num_modes, [](IoData &) {});
  REQUIRE(exact.num_converged >= num_modes);
  CHECK(exact.reduced_basis_size == 0);
  CHECK(exact.reduced_stats.reduced_solves == 0);
  CHECK(exact.reduced_stats.exact_solves == 1);

  auto reduced = SolveRectangularModes(
      1000.0, 500.0, 500.0, 4.0, 2, num_modes, [](IoData &) {}, true, 21);
  REQUIRE(reduced.num_converged >= num_modes);
  REQUIRE(reduced.reduced_basis_size >= num_modes);
  CHECK(reduced.reduced_stats.reduced_solves == 21);
  CHECK(reduced.reduced_stats.exact_solves == 2);
  CHECK(reduced.reduced_stats.worst_residual <= reduced.reduced_tol);
  REQUIRE(reduced.reduced_backward_errors.size() == num_modes);
  for (double error : reduced.reduced_backward_errors)
  {
    CHECK(error <= reduced.reduced_tol);
  }
  CHECK(reduced.complex_exact_converged >= num_modes);
  CHECK(reduced.reduced_stats.offline_basis_rank >= num_modes);
  CHECK(reduced.reduced_stats.online_basis_cap >=
        reduced.reduced_stats.offline_basis_rank + 4 * num_modes);
  for (int i = 0; i < num_modes; i++)
  {
    CHECK_THAT(reduced.kn[i].real(), WithinRel(exact.kn[i].real(), 1.0e-6));
    CHECK_THAT(reduced.kn[i].imag(), WithinAbs(exact.kn[i].imag(), 1.0e-8));
  }
}

TEST_CASE("ModeEigenSolver reduced basis capacity lifecycle",
          "[boundarymodeoperator][Serial]")
{
  auto result =
      SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 1, [](IoData &) {}, true, 1, 1);
  REQUIRE(result.num_converged >= 1);
  CHECK(result.reduced_stats.offline_basis_rank == 1);
  CHECK(result.reduced_stats.online_basis_cap == 5);
  CHECK(result.reduced_stats.reduced_solves == 1);
}

TEST_CASE("ModeEigenSolver Impedance shifts kn", "[boundarymodeoperator][Serial]")
{
  auto pec_result = SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 3, [](IoData &) {});

  // Use a large enough inductance so the impedance shift is well above numerical
  // noise across different BLAS/LAPACK implementations.
  auto imp_result = SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 3,
                                          [](IoData &iodata)
                                          {
                                            iodata.boundaries.pec.attributes = {1, 3, 4};
                                            auto &imp =
                                                iodata.boundaries.impedance.emplace_back();
                                            imp.attributes = {2};
                                            imp.Ls = 1.0e-8;
                                          });

  REQUIRE(pec_result.num_converged >= 1);
  REQUIRE(imp_result.num_converged >= 1);

  CHECK(imp_result.kn[0].real() > pec_result.kn[0].real());
}

TEST_CASE("ModeEigenSolver rational impedance affine component",
          "[boundarymodeoperator][Serial]")
{
  auto configure_rational = [](IoData &iodata)
  {
    iodata.boundaries.pec.attributes = {1, 3, 4};
    auto &rz = iodata.boundaries.rational_impedance.emplace_back();
    rz.attributes = {2};
    // Series RL impedance Z(s) = R + sL gives the genuinely rational Robin coefficient
    // g(s) = s/(R+sL), rather than merely duplicating the linear resistive component.
    rz.num = {1.0e-12, 50.0};
    rz.den = {1.0};
  };
  auto exact = SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 1, configure_rational);
  auto reduced =
      SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 1, configure_rational, true);

  REQUIRE(exact.num_converged >= 1);
  REQUIRE(reduced.num_converged >= 1);
  CHECK(reduced.reduced_stats.reduced_solves == 1);
  CHECK(reduced.reduced_stats.worst_residual <= reduced.reduced_tol);
  CHECK_THAT(reduced.kn[0].real(), WithinRel(exact.kn[0].real(), 1.0e-6));
  CHECK_THAT(reduced.kn[0].imag(), WithinAbs(exact.kn[0].imag(), 1.0e-8));
}

TEST_CASE("ModeEigenSolver reduced rejection fallback enrichment",
          "[boundarymodeoperator][Serial][Parallel]")
{
  auto configure_rational = [](IoData &iodata)
  {
    iodata.boundaries.pec.attributes = {1, 3};
    auto &rz = iodata.boundaries.rational_impedance.emplace_back();
    rz.attributes = {2, 4};
    rz.num = {1.0e-12, 1.0};
    rz.den = {1.0};
  };
  auto result = SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 1, configure_rational,
                                      true, 2, 16, {0.1}, 1.0e-4, 1.0e-10);

  REQUIRE(result.num_converged >= 1);
  CHECK(result.first_target_stats.offline_basis_rank == 1);
  CHECK(result.first_target_stats.exact_solves == 2);
  CHECK(result.first_target_stats.fallbacks == 1);
  CHECK(result.first_target_stats.reduced_solves == 0);
  CHECK(result.first_target_basis_size > result.first_target_stats.offline_basis_rank);

  CHECK(result.reduced_stats.exact_solves == 2);
  CHECK(result.reduced_stats.fallbacks == 1);
  CHECK(result.reduced_stats.reduced_solves == 1);
  CHECK(result.reduced_basis_size == result.first_target_basis_size);
  CHECK(result.reduced_stats.worst_residual <= result.reduced_tol);
  REQUIRE(result.final_reduced_backward_errors.size() == 1);
  CHECK(result.final_reduced_backward_errors[0] <= result.reduced_tol);
  REQUIRE(result.final_kn.size() == 1);
  CHECK_THAT(result.final_kn[0].real(), WithinRel(result.kn[0].real(), 1.0e-6));
  CHECK_THAT(result.final_kn[0].imag(), WithinAbs(result.kn[0].imag(), 1.0e-8));
}

TEST_CASE("ModeEigenSolver Conductivity adds loss", "[boundarymodeoperator][Serial]")
{
  auto pec_result = SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 3, [](IoData &) {});

  auto configure_conductivity = [](IoData &iodata)
  {
    iodata.boundaries.pec.attributes = {1, 3, 4};
    auto &cond = iodata.boundaries.conductivity.emplace_back();
    cond.attributes = {2};
    cond.sigma = 5.0e7;
    cond.h = 0.001;
  };
  auto cond_result =
      SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 3, configure_conductivity);
  auto cond_reduced =
      SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 1, configure_conductivity, true);

  REQUIRE(pec_result.num_converged >= 1);
  REQUIRE(cond_result.num_converged >= 1);
  REQUIRE(cond_reduced.num_converged >= 1);

  CHECK(std::abs(cond_result.kn[0].imag()) > std::abs(pec_result.kn[0].imag()));
  CHECK(cond_reduced.reduced_stats.reduced_solves == 1);
  CHECK(cond_reduced.reduced_stats.worst_residual <= cond_reduced.reduced_tol);
  CHECK_THAT(cond_reduced.kn[0].real(), WithinRel(cond_result.kn[0].real(), 1.0e-6));
  CHECK_THAT(cond_reduced.kn[0].imag(), WithinAbs(cond_result.kn[0].imag(), 1.0e-8));
}

TEST_CASE("ModeEigenSolver London slab", "[boundarymodeoperator][Serial][Parallel]")
{
  // London superconductor slabs of thickness d, compared to exact solutions.
  constexpr double lambda = 0.1;
  SECTION("Out-of-plane current")
  {
    // TM0 mode of a parallel plate guide with gap h, with the slab carrying the current
    // along the propagation direction, compared to the exact dispersion relation. The
    // frequency is high enough for the propagation constant to be well conditioned.
    constexpr double h = 0.5, freq_ghz = 500.0;
    const double k0 =
        2.0 * std::numbers::pi * freq_ghz * 1.0e9 / electromagnetics::c0_ * 1.0e-6;
    for (double d : {0.01, 0.1, 1.0})
    {
      for (bool pec_back : {false, true})
      {
        const double n_exact = ExactLondonSlabIndex(h, d, lambda, k0, pec_back);
        const double n_eff =
            SolveLondonSlabMode(1.0, h, d, lambda, pec_back, false, freq_ghz, n_exact, 3);
        CAPTURE(d, pec_back, n_exact, n_eff);
        CHECK_THAT(n_eff, WithinRel(n_exact, 1.0e-8));
      }
    }
  }
  SECTION("In-plane current")
  {
    // TE1 mode between the PEC plate and the slab (PEC sides), with E_x = sin(k_y y) in the
    // gap and the current along the slab, perpendicular to the propagation direction. In
    // the slab, q^2 = 1/lambda^2 - k_y^2 and k_y cot(k_y h) = -q tanh(q d) (free back
    // face).
    constexpr double h = 500.0, d = 0.5, freq_ghz = 500.0;
    const double k0 =
        2.0 * std::numbers::pi * freq_ghz * 1.0e9 / electromagnetics::c0_ * 1.0e-6;
    auto F = [&](double ky)
    {
      const double q = std::sqrt(1.0 / (lambda * lambda) - ky * ky);
      return ky / std::tan(ky * h) + q * std::tanh(q * d);
    };
    double lo = 0.9 * std::numbers::pi / h, hi = (1.0 - 1.0e-12) * std::numbers::pi / h;
    for (int it = 0; it < 200; it++)
    {
      const double mid = 0.5 * (lo + hi);
      ((F(lo) * F(mid) <= 0.0) ? hi : lo) = mid;
    }
    const double ky = 0.5 * (lo + hi);
    const double n_exact = std::sqrt(1.0 - (ky / k0) * (ky / k0));
    const double n_pec =
        std::sqrt(1.0 - (std::numbers::pi / h / k0) * (std::numbers::pi / h / k0));
    const double n_eff =
        SolveLondonSlabMode(50.0, h, d, lambda, false, true, freq_ghz, n_exact, 3, 32);
    CAPTURE(n_exact, n_pec, n_eff);
    CHECK_THAT(n_eff - n_pec, WithinRel(n_exact - n_pec, 1.0e-5));
  }
}

TEST_CASE("BoundaryModeOperator interior PEC on a nonconforming mesh",
          "[boundarymodeoperator][Serial][Parallel]")
{
  // A PEC line through the interior of the cross-section at x = 4, and an interior boundary
  // without boundary condition at y = 2 for x > 4, with the elements of the quadrant x > 4,
  // y > 2 refined. The edges of both lines on the unrefined side are master edges without
  // boundary elements: those of the PEC line must be essential, and those of the other line
  // not.
  const auto type = GENERATE(mfem::Element::TRIANGLE, mfem::Element::QUADRILATERAL);
  const int order = GENERATE(1, 2);
  MPI_Comm comm = Mpi::World();
  constexpr double width = 8.0, height = 4.0, x_pec = 4.0, y_other = 2.0, eps = 1.0e-9;
  constexpr int pec_attr = 5, other_attr = 6;

  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(8, 4, type, false, width, height));
  {
    auto &smesh = *serial_mesh;
    mfem::Array<int> v;
    for (int f = 0; f < smesh.GetNumFaces(); f++)
    {
      int e1, e2;
      smesh.GetFaceElements(f, &e1, &e2);
      if (e1 < 0 || e2 < 0)
      {
        continue;
      }
      smesh.GetFaceVertices(f, v);
      const double *a = smesh.GetVertex(v[0]), *b = smesh.GetVertex(v[1]);
      if (std::abs(a[0] - x_pec) < eps && std::abs(b[0] - x_pec) < eps)
      {
        smesh.AddBdrSegment(v[0], v[1], pec_attr);
      }
      else if (std::abs(a[1] - y_other) < eps && std::abs(b[1] - y_other) < eps &&
               std::min(a[0], b[0]) > x_pec - eps)
      {
        smesh.AddBdrSegment(v[0], v[1], other_attr);
      }
    }
    smesh.FinalizeTopology();
    smesh.Finalize(true, true);
    smesh.EnsureNCMesh(true);
    mfem::Array<int> marked;
    mfem::Vector c;
    for (int e = 0; e < smesh.GetNE(); e++)
    {
      smesh.GetElementCenter(e, c);
      if (c(0) > x_pec && c(1) > y_other)
      {
        marked.Append(e);
      }
    }
    smesh.GeneralRefinement(marked, 1, 0);
  }

  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.model.Lc = 1.0;
  auto &material = iodata.domains.materials.emplace_back();
  material.attributes = {1};
  material.epsilon_r.s = {1.0, 1.0, 1.0};
  iodata.boundaries.pec.attributes = {pec_attr};
  iodata.solver.order = order;
  iodata.solver.boundary_mode.freq = 1.0;
  iodata.solver.boundary_mode.n = 1;
  iodata.NondimensionalizeInputs(serial_mesh);
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  iodata.CheckConfiguration();
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(std::move(par_mesh)));
  MaterialOperator mat_op(iodata, *mesh.back());
  BoundaryModeOperator mode_op(iodata, mesh, mat_op);

  auto CheckDbcTDofs = [&](const FiniteElementSpace &fespace,
                           const mfem::Array<int> &dbc_tdof_list, bool edge_dofs)
  {
    const auto &fes = fespace.Get();
    const auto &pmesh = *fes.GetParMesh();
    std::vector<bool> is_dbc(fes.GetTrueVSize(), false);
    for (auto t : dbc_tdof_list)
    {
      is_dbc[t] = true;
    }
    // Numbers of checked true DOFs on the PEC line, on the other line, and on the master
    // edges of the PEC line (which have DOFs if the space has edge DOFs).
    int counts[3] = {0, 0, 0};
    auto Check = [&](const mfem::Array<int> &dofs, bool pec, bool master)
    {
      for (auto d : dofs)
      {
        const int t = fes.GetLocalTDofNumber((d >= 0) ? d : -1 - d);
        if (t < 0)
        {
          continue;
        }
        CHECK(is_dbc[t] == pec);
        counts[pec ? 0 : 1]++;
        counts[2] += (pec && master);
      }
    };
    auto OnPEC = [&](const double *x) { return std::abs(x[0] - x_pec) < eps; };
    auto OnOther = [&](const double *x)
    { return std::abs(x[1] - y_other) < eps && x[0] > x_pec + eps && x[0] < width - eps; };
    const auto &edge_list = pmesh.ncmesh->GetEdgeList();
    mfem::Array<int> v, dofs;
    for (int e = 0; e < pmesh.GetNEdges(); e++)
    {
      pmesh.GetEdgeVertices(e, v);
      const double *a = pmesh.GetVertex(v[0]), *b = pmesh.GetVertex(v[1]);
      const bool pec = OnPEC(a) && OnPEC(b);
      const bool other = std::abs(a[1] - y_other) < eps && std::abs(b[1] - y_other) < eps &&
                         std::min(a[0], b[0]) > x_pec - eps;
      if (pec || other)
      {
        fes.GetEdgeInteriorDofs(e, dofs);
        Check(dofs, pec,
              edge_list.GetMeshIdAndType(e).type ==
                  mfem::NCMesh::NCList::MeshIdType::MASTER);
      }
    }
    for (int i = 0; i < pmesh.GetNV(); i++)
    {
      const double *x = pmesh.GetVertex(i);
      if (OnPEC(x) || OnOther(x))
      {
        fes.GetVertexDofs(i, dofs);
        Check(dofs, OnPEC(x), false);
      }
    }
    Mpi::GlobalSum(3, counts, comm);
    CHECK(counts[0] > 0);
    CHECK(counts[1] > 0);
    CHECK((counts[2] > 0) == edge_dofs);
  };
  CheckDbcTDofs(mode_op.GetNDSpace(), mode_op.GetNDDbcTDofLists().back(), true);
  CheckDbcTDofs(mode_op.GetH1Space(), mode_op.GetH1DbcTDofLists().back(), order > 1);
}

// BoundaryMode cross-section extracted from a nonconforming 3D mesh: the wave port of
// cpw_wave_2dmode on the uncracked mesh, with the elements next to the port refined on one
// side of the traces (as an adapted mesh can be). The 2D mesh must keep its boundary
// elements through partitioning and uniform refinement, and the DoFs on the PEC edges
// (master edges included) must be essential, and no others.
TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "BoundaryMode cross-section of a nonconforming 3D mesh",
                 "[boundarymodeoperator][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const auto cpw_dir = fs::path(PALACE_TEST_DATA_DIR) / "regression" / "input" / "cpw";
  auto buffer = PreprocessFile((cpw_dir / "cpw_wave_2dmode.json").string().c_str());
  auto setup = nlohmann::json::parse(buffer);
  setup["Model"]["Mesh"] = (cpw_dir / setup["Model"]["Mesh"].get<std::string>()).string();
  setup["Model"]["CrackInternalBoundaryElements"] = false;
  setup["Model"]["Refinement"] = {{"UniformLevels", 1}};
  setup["Problem"]["Output"] = temp_dir.string();
  IoData iodata(setup, /*print=*/false);
  const std::vector<int> pec = iodata.boundaries.pec.attributes;

  auto smesh = mesh::Load(iodata, comm);
  if (smesh)
  {
    // Refine the elements next to the cross-section on one side of the traces, the
    // interior PEC boundary elements.
    double z_trace = mfem::infinity();
    for (int be = 0; be < smesh->GetNBE(); be++)
    {
      int e1, e2;
      smesh->GetFaceElements(smesh->GetBdrElementFaceIndex(be), &e1, &e2);
      if (e2 >= 0 &&
          std::find(pec.begin(), pec.end(), smesh->GetBdrAttribute(be)) != pec.end())
      {
        mfem::Array<int> v;
        smesh->GetBdrElementVertices(be, v);
        z_trace = std::min(z_trace, smesh->GetVertex(v[0])[2]);
      }
    }
    REQUIRE(z_trace < mfem::infinity());
    const auto &attrs = iodata.solver.boundary_mode.attributes;
    smesh->EnsureNCMesh(true);
    mfem::Array<int> marked;
    mfem::Vector c;
    for (int be = 0; be < smesh->GetNBE(); be++)
    {
      if (std::find(attrs.begin(), attrs.end(), smesh->GetBdrAttribute(be)) == attrs.end())
      {
        continue;
      }
      int e1, e2;
      smesh->GetFaceElements(smesh->GetBdrElementFaceIndex(be), &e1, &e2);
      smesh->GetElementCenter(e1, c);
      if (c(2) > z_trace)
      {
        marked.Append(e1);
      }
    }
    marked.Sort();
    marked.Unique();
    REQUIRE(marked.Size() > 0);
    smesh->GeneralRefinement(marked, 1, 0);
  }
  BoundaryModeSolver solver(iodata, Mpi::Root(comm), Mpi::Size(comm), 1, nullptr);
  solver.Preprocess(iodata, smesh, comm);
  std::vector<std::unique_ptr<mfem::ParMesh>> mfem_mesh;
  mfem_mesh.push_back(mesh::Partition(iodata, std::move(smesh), comm));

  // Segments of the PEC boundary elements of the 2D mesh before its uniform refinement,
  // gathered from all ranks (the other wave ports are relabelled as PEC), and the numbers
  // of PEC boundary elements and of those in the interior.
  auto CountPEC = [&](const mfem::ParMesh &pmesh, std::vector<double> *segs)
  {
    std::array<int, 2> n = {0, 0};
    for (int be = 0; be < pmesh.GetNBE(); be++)
    {
      if (std::find(pec.begin(), pec.end(), pmesh.GetBdrAttribute(be)) == pec.end())
      {
        continue;
      }
      n[0]++;
      n[1] += pmesh.FaceIsInterior(pmesh.GetBdrElementFaceIndex(be));
      mfem::Array<int> v;
      pmesh.GetBdrElementVertices(be, v);
      for (int i : v)
      {
        double x[2];
        pmesh.GetNode(i, x);
        if (segs)
        {
          segs->insert(segs->end(), x, x + 2);
        }
      }
    }
    Mpi::GlobalSum(2, n.data(), comm);
    return n;
  };
  std::vector<double> segs;
  const auto ne = mfem_mesh.back()->GetGlobalNE();
  const auto nbe_pec = CountPEC(*mfem_mesh.back(), &segs);
  CHECK(nbe_pec[0] > 0);
  CHECK(nbe_pec[1] > 0);
  mesh::RefineMesh(iodata, mfem_mesh);
  std::vector<std::unique_ptr<Mesh>> mesh;
  for (auto &m : mfem_mesh)
  {
    mesh.push_back(std::make_unique<Mesh>(std::move(m)));
  }
  const auto &pmesh = mesh.back()->Get();
  REQUIRE(pmesh.Nonconforming());
  REQUIRE(pmesh.GetGlobalNE() == 4 * ne);
  CHECK(CountPEC(pmesh, nullptr) == std::array<int, 2>{2 * nbe_pec[0], 2 * nbe_pec[1]});
  MaterialOperator mat_op(iodata, *mesh.back());
  BoundaryModeOperator mode_op(iodata, mesh, mat_op);

  int nloc = static_cast<int>(segs.size()), nproc = Mpi::Size(comm);
  std::vector<int> sizes(nproc), displs(nproc, 0);
  MPI_Allgather(&nloc, 1, MPI_INT, sizes.data(), 1, MPI_INT, comm);
  for (int r = 1; r < nproc; r++)
  {
    displs[r] = displs[r - 1] + sizes[r - 1];
  }
  std::vector<double> all_segs(displs.back() + sizes.back());
  MPI_Allgatherv(segs.data(), nloc, MPI_DOUBLE, all_segs.data(), sizes.data(),
                 displs.data(), MPI_DOUBLE, comm);
  auto OnPEC = [&](const double *x)
  {
    for (std::size_t k = 0; k < all_segs.size(); k += 4)
    {
      const double *a = &all_segs[k], *b = &all_segs[k + 2];
      const double ab2 = (b[0] - a[0]) * (b[0] - a[0]) + (b[1] - a[1]) * (b[1] - a[1]);
      const double ax_ab = (x[0] - a[0]) * (b[0] - a[0]) + (x[1] - a[1]) * (b[1] - a[1]);
      const double ax2 = (x[0] - a[0]) * (x[0] - a[0]) + (x[1] - a[1]) * (x[1] - a[1]);
      const double t = ax_ab / ab2;
      if (t > -1.0e-9 && t < 1.0 + 1.0e-9 && ax2 - t * t * ab2 < 1.0e-12 * ab2)
      {
        return true;
      }
    }
    return false;
  };

  // Numbers of checked true DoFs on the PEC edges and off them, and on master edges.
  int counts[3] = {0, 0, 0};
  auto CheckDbcTDofs =
      [&](const FiniteElementSpace &fespace, const mfem::Array<int> &dbc_tdof_list)
  {
    const auto &fes = fespace.Get();
    std::vector<bool> is_dbc(fes.GetTrueVSize(), false);
    for (auto t : dbc_tdof_list)
    {
      is_dbc[t] = true;
    }
    auto Check = [&](const mfem::Array<int> &dofs, bool on_pec, bool master)
    {
      for (auto d : dofs)
      {
        const int t = fes.GetLocalTDofNumber((d >= 0) ? d : -1 - d);
        if (t < 0)
        {
          continue;
        }
        CHECK(is_dbc[t] == on_pec);
        counts[on_pec ? 0 : 1]++;
        counts[2] += (on_pec && master);
      }
    };
    const auto &edge_list = pmesh.ncmesh->GetEdgeList();
    mfem::Array<int> v, dofs;
    for (int e = 0; e < pmesh.GetNEdges(); e++)
    {
      double x0[2], x1[2];
      pmesh.GetEdgeVertices(e, v);
      pmesh.GetNode(v[0], x0);
      pmesh.GetNode(v[1], x1);
      const double mid[2] = {0.5 * (x0[0] + x1[0]), 0.5 * (x0[1] + x1[1])};
      fes.GetEdgeInteriorDofs(e, dofs);
      Check(dofs, OnPEC(mid),
            edge_list.GetMeshIdAndType(e).type == mfem::NCMesh::NCList::MeshIdType::MASTER);
    }
    for (int i = 0; i < pmesh.GetNV(); i++)
    {
      double x[2];
      pmesh.GetNode(i, x);
      fes.GetVertexDofs(i, dofs);
      Check(dofs, OnPEC(x), false);
    }
  };
  CheckDbcTDofs(mode_op.GetNDSpace(), mode_op.GetNDDbcTDofLists().back());
  CheckDbcTDofs(mode_op.GetH1Space(), mode_op.GetH1DbcTDofLists().back());
  Mpi::GlobalSum(3, counts, comm);
  CHECK(counts[0] > 0);
  CHECK(counts[1] > 0);
  CHECK(counts[2] > 0);
}

TEST_CASE("Mode target distance", "[boundarymodeoperator][Serial]")
{
  using mode_assembly::TargetDistance;

  // Reported wave-port selections of a strongly evanescent mode of a lossy cross-section
  // whose real part lies closer to the target than that of the guided mode, so ranking by
  // |Re{kn} - kn_target| picked it: a CPW with absorbing and conductivity boundaries at
  // 1 GHz (awslabs/palace#920, in 1/m) and examples/cpw/cpw_wave_uniform.json at 0.1 GHz
  // after a sweep to 50 GHz (awslabs/palace#996, nondimensional).
  struct ReportedSelection
  {
    int issue;
    double kn_target;
    std::complex<double> kn_guided, kn_evanescent;
  };
  const ReportedSelection reported[] = {
      {920, 60.2, {38.07, -4.445}, {51.03, 5677.0}},
      {996, 2.981703e-2, {1.988221e-2, -5.136255e-7}, {3.341085e-2, 1.574046e1}}};
  for (const auto &r : reported)
  {
    CAPTURE(r.issue);
    REQUIRE(std::abs(r.kn_evanescent.real() - r.kn_target) <
            std::abs(r.kn_guided.real() - r.kn_target));
    CHECK(TargetDistance(r.kn_guided, r.kn_target) <
          TargetDistance(r.kn_evanescent, r.kn_target));
  }
  const double kn_target = reported[1].kn_target;
  const std::complex<double> kn_cpw = reported[1].kn_guided;
  const std::complex<double> kn_evanescent = reported[1].kn_evanescent;

  // Lossless propagating modes below the target rank by descending kn, ahead of all
  // lossless evanescent modes, which rank by ascending attenuation.
  const double kt = 1.0;
  CHECK(TargetDistance(0.9 * kt, kt) < TargetDistance(0.5 * kt, kt));
  CHECK(TargetDistance(0.5 * kt, kt) < TargetDistance(1.0e-3 * kt, kt));
  CHECK(TargetDistance(1.0e-3 * kt, kt) <
        TargetDistance(std::complex<double>(0.0, 1.0e-3 * kt), kt));
  CHECK(TargetDistance(std::complex<double>(0.0, 1.0e-3 * kt), kt) <
        TargetDistance(std::complex<double>(0.0, 10.0 * kt), kt));

  // Invariant under conjugation: a near-cutoff mode may land on either side of the branch
  // cut of the principal square root.
  for (const auto kn : {kn_cpw, kn_evanescent, std::complex<double>(0.3, -0.7)})
  {
    CHECK(TargetDistance(std::conj(kn), kn_target) == TargetDistance(kn, kn_target));
  }

  // Non-finite values rank last, after any finite candidate.
  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  CHECK(TargetDistance(std::complex<double>(nan, 0.0), kt) == inf);
  CHECK(TargetDistance(std::complex<double>(0.0, nan), kt) == inf);
  CHECK(TargetDistance(std::complex<double>(nan, inf), kt) == inf);
  CHECK(TargetDistance(std::complex<double>(inf, 0.0), kt) == inf);
  CHECK(TargetDistance(kn_evanescent, kt) < inf);
}

TEST_CASE("ModeEigenSolver ranks evanescent modes below cutoff by attenuation",
          "[boundarymodeoperator][Serial][Parallel]")
{
  // Below the TE10 cutoff (75 GHz for this 1000 x 500 um guide with eps = 4) every mode of
  // a lossless guide is evanescent: the real parts of the computed kn are round-off and
  // carry no ranking information, and ranking on them returned an arbitrary evanescent mode
  // as mode 1 (e.g. |S11| > 0 dB for a lossless shorted guide driven below cutoff). Both
  // the exact and the reduced solve must rank the modes least-attenuated first, so mode 1
  // is the continuation of TE10 through cutoff.
  constexpr int num_modes = 3;
  constexpr double freq_ghz = 60.0, epsilon_r = 4.0;
  auto exact = SolveRectangularModes(1000.0, 500.0, freq_ghz, epsilon_r, 2, num_modes,
                                     [](IoData &) {});
  REQUIRE(exact.num_converged >= num_modes);
  for (int i = 0; i < exact.num_converged; i++)
  {
    CAPTURE(i, exact.kn[i]);
    REQUIRE(std::abs(exact.kn[i].real()) <= 1.0e-6 * std::abs(exact.kn[i].imag()));
    if (i > 0)
    {
      CHECK(std::abs(exact.kn[i - 1].imag()) <=
            (1.0 + 1.0e-6) * std::abs(exact.kn[i].imag()));
    }
  }

  // TE10 attenuation constant sqrt((pi/a)^2 - eps k0^2), nondimensionalized by 1 um.
  const double kc = std::numbers::pi / 1000.0;
  const double k = std::sqrt(epsilon_r) * 2.0 * std::numbers::pi * freq_ghz * 1.0e9 /
                   electromagnetics::c0_ * 1.0e-6;
  CHECK_THAT(std::abs(exact.kn[0].imag()), WithinRel(std::sqrt(kc * kc - k * k), 1.0e-3));

  // Train the reduced model below cutoff as well (0.9 and 1.1 times the frequency). Compare
  // kn^2, which does not depend on the sign of i|kn| returned for an evanescent mode.
  auto reduced = SolveRectangularModes(
      1000.0, 500.0, freq_ghz, epsilon_r, 2, num_modes, [](IoData &) {}, true);
  REQUIRE(reduced.reduced_stats.reduced_solves == 1);
  REQUIRE(reduced.kn.size() == num_modes);
  for (int i = 0; i < num_modes; i++)
  {
    CAPTURE(i, reduced.kn[i], exact.kn[i]);
    CHECK(std::abs(reduced.kn[i] * reduced.kn[i] - exact.kn[i] * exact.kn[i]) <=
          1.0e-6 * std::norm(exact.kn[i]));
  }
}

TEST_CASE("ModeEigenSolver ranks lossy modes by complex target distance",
          "[boundarymodeoperator][Serial][Parallel]")
{
  // A resistive top wall gives the modes attenuation comparable to the spacing of their
  // phase constants, so ranking by |Re{kn} - kn_target| and by |kn - kn_target| disagree.
  // Exact and reduced solves must both return the modes in ascending complex distance.
  auto configure_resistive_wall = [](IoData &iodata)
  {
    iodata.boundaries.pec.attributes = {1, 2, 4};
    auto &imp = iodata.boundaries.impedance.emplace_back();
    imp.attributes = {3};
    imp.Rs = 377.0;
  };
  constexpr int num_modes = 3;
  auto exact = SolveRectangularModes(1000.0, 500.0, 200.0, 4.0, 2, num_modes,
                                     configure_resistive_wall);
  REQUIRE(exact.num_converged >= num_modes);
  const double kn_target = exact.kn_target;
  for (int i = 1; i < exact.num_converged; i++)
  {
    CAPTURE(i, exact.kn[i - 1], exact.kn[i]);
    CHECK(mode_assembly::TargetDistance(exact.kn[i - 1], kn_target) <=
          mode_assembly::TargetDistance(exact.kn[i], kn_target));
  }

  // Guard the premise: some returned mode has a real part closer to the target than the
  // first-ranked mode, so a real-part ranking would select differently.
  bool real_part_ranking_differs = false;
  for (int i = 1; i < exact.num_converged; i++)
  {
    real_part_ranking_differs =
        real_part_ranking_differs ||
        std::abs(exact.kn[i].real() - kn_target) < std::abs(exact.kn[0].real() - kn_target);
  }
  REQUIRE(real_part_ranking_differs);

  // The reduced wave-port solve, trained with all ranked modes around the target, must
  // select the same modes in the same order.
  auto reduced = SolveRectangularModes(1000.0, 500.0, 200.0, 4.0, 2, num_modes,
                                       configure_resistive_wall, true, 1, 16,
                                       {0.95, 0.98, 1.02, 1.05}, 1.0e-2);
  REQUIRE(reduced.reduced_stats.reduced_solves == 1);
  REQUIRE(reduced.kn.size() == num_modes);
  for (int i = 0; i < num_modes; i++)
  {
    CAPTURE(i, reduced.kn[i], exact.kn[i]);
    CHECK(std::abs(reduced.kn[i] - exact.kn[i]) <= 1.0e-5 * std::abs(exact.kn[i]));
  }
}

}  // namespace palace
