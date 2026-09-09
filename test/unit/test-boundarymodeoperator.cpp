// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <complex>
#include <functional>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
#include "models/farfieldboundaryoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/surfaceconductivityoperator.hpp"
#include "models/surfaceimpedanceoperator.hpp"
#include "models/surfacerationalimpedanceoperator.hpp"
#include "models/waveportoperator.hpp"
#include "utils/communication.hpp"
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

  double omega =
      2.0 * M_PI * iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(freq_ghz);

  // ModeEigenSolver requires a positive Krylov subspace size (num_vec). Mirror the
  // formula used by IoData::CheckConfiguration for eigenmode.max_size.
  const int num_vec = std::max(2 * num_modes, num_modes + 15);
  ModeEigenSolver mode_solver(mat_op, nullptr, surf_z_op, farfield_op, surf_sigma_op,
                              surf_rz_op, nd_fespace, h1_fespace, dbc_tdof_list, num_modes,
                              num_vec, eig_tol, EigenvalueSolver::WhichType::LARGEST_REAL,
                              iodata.solver.linear, iodata.solver.boundary_mode.type, 0,
                              nd_fespace.GetComm());

  auto solve_at = [&](std::complex<double> w)
  {
    const double kn_target = w.real() * std::sqrt(1.1 * mat_op.GetMaxMuEpsilon());
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

}  // namespace palace
