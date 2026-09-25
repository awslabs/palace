// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
#include "linalg/vector.hpp"
#include "models/boundarymodeoperator.hpp"
#include "models/farfieldboundaryoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/modeeigensolver.hpp"
#include "models/modeoperatorassembly.hpp"
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

// Solve for the fundamental mode of the rectangular waveguide cross-section (at 100 GHz,
// where it is the only propagating mode) through the 2D BoundaryModeOperator path used by
// the BoundaryMode driver, which uses p-multigrid preconditioning when mg_max_levels > 1.
std::vector<std::complex<double>>
SolveRectangularModesMultigrid(int mg_max_levels,
                               const std::function<void(IoData &)> &configure_bcs)
{
  constexpr double width = 1000.0, height = 500.0, freq_ghz = 100.0, eig_tol = 1.0e-8;
  constexpr int order = 2, num_modes = 1;
  MPI_Comm comm = Mpi::World();
  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.problem.type = ProblemType::BOUNDARYMODE;
  iodata.model.Lc = 1.0;

  auto &material = iodata.domains.materials.emplace_back();
  material.attributes = {1};
  material.epsilon_r.s = {4.0, 4.0, 4.0};
  iodata.boundaries.pec.attributes = {1, 2, 3, 4};
  configure_bcs(iodata);

  iodata.solver.order = order;
  iodata.solver.boundary_mode.freq = freq_ghz;
  iodata.solver.boundary_mode.n = num_modes;
  iodata.solver.boundary_mode.tol = eig_tol;
  iodata.solver.linear.tol = 1.0e-10;
  iodata.solver.linear.max_it = 200;
  iodata.solver.linear.mg_max_levels = mg_max_levels;

  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(10, 5, mfem::Element::TRIANGLE, false, width, height));
  iodata.NondimensionalizeInputs(serial_mesh);
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(
      std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(comm, *serial_mesh)));
  iodata.CheckConfiguration();

  MaterialOperator mat_op(iodata, *mesh.back());
  BoundaryModeOperator mode_op(iodata, mesh, mat_op);
  REQUIRE(mode_op.GetNDSpaceHierarchy().GetNumLevels() ==
          static_cast<std::size_t>(mg_max_levels));

  const int nd_size = mode_op.GetNDTrueVSize();
  mfem::Array<int> dbc_tdof_list;
  dbc_tdof_list.Append(mode_op.GetNDDbcTDofLists().back());
  for (auto tdof : mode_op.GetH1DbcTDofLists().back())
  {
    dbc_tdof_list.Append(nd_size + tdof);
  }

  const int num_vec = std::max(2 * num_modes, num_modes + 15);
  ModeEigenSolver mode_solver(mode_op, dbc_tdof_list, num_modes, num_vec, eig_tol,
                              EigenvalueSolver::WhichType::LARGEST_REAL,
                              iodata.solver.linear, iodata.solver.boundary_mode.type, 0);
  const double omega =
      2.0 * M_PI * iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(freq_ghz);
  const double kn_target = omega * std::sqrt(1.1 * mat_op.GetMaxMuEpsilon());
  auto result = mode_solver.Solve(omega, -kn_target * kn_target);

  std::vector<std::complex<double>> kn;
  for (int i = 0; i < result.num_converged; i++)
  {
    kn.push_back(mode_solver.GetPropagationConstant(i));
  }
  return kn;
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

TEST_CASE("ModeEigenSolver p-multigrid preconditioning",
          "[boundarymodeoperator][Serial][Parallel]")
{
  // The H1 block of the multigrid preconditioner is negative definite (the diffusion term
  // keeps its sign from the integration by parts), and essential DOFs are eliminated with a
  // unit diagonal, so its Chebyshev smoothers see a diagonal of mixed sign (with essential
  // boundaries) or a negative diagonal (without). The multigrid preconditioned solve must
  // reproduce the mode of the sparse direct one.
  auto check = [](const std::function<void(IoData &)> &configure_bcs)
  {
    const auto direct = SolveRectangularModesMultigrid(1, configure_bcs);
    const auto multigrid = SolveRectangularModesMultigrid(2, configure_bcs);
    REQUIRE(direct.size() >= 1);
    REQUIRE(multigrid.size() >= 1);
    CAPTURE(direct[0], multigrid[0]);
    CHECK_THAT(multigrid[0].real(), WithinRel(direct[0].real(), 1.0e-6));
    CHECK_THAT(multigrid[0].imag(),
               WithinAbs(direct[0].imag(), 1.0e-6 * std::abs(direct[0])));
  };

  SECTION("PEC walls")
  {
    check([](IoData &) {});
  }

  SECTION("Impedance walls")
  {
    check(
        [](IoData &iodata)
        {
          iodata.boundaries.pec.attributes.clear();
          auto &imp = iodata.boundaries.impedance.emplace_back();
          imp.attributes = {1, 2, 3, 4};
          imp.Ls = 1.0e-8;
        });
  }
}

TEST_CASE("ModeOperatorModel farfield damping uses the neighboring material",
          "[boundarymodeoperator][Serial][Parallel]")
{
  // The first-order absorbing boundary condition adds iω/Z₀ times the boundary mass of the
  // in-plane (ND, tangential trace) and out-of-plane (H1, with the negative sign of the
  // H1 block) field components, where 1/Z₀ = √(ε/μ) is taken from the domain material
  // adjacent to each boundary element. The absorbing side x = W spans two materials, with
  // ε = 4 (1/Z₀ = 2) along y < H/4 and vacuum (1/Z₀ = 1) above, so a field tangential to
  // it with unit magnitude gives ±(2 H/4 + 3 H/4).
  constexpr double W = 2.0, H = 1.0;
  MPI_Comm comm = Mpi::World();
  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.problem.type = ProblemType::BOUNDARYMODE;
  iodata.model.Lc = 1.0;
  {
    auto &substrate = iodata.domains.materials.emplace_back();
    substrate.attributes = {1};
    substrate.epsilon_r.s = {4.0, 4.0, 4.0};
    auto &vacuum = iodata.domains.materials.emplace_back();
    vacuum.attributes = {2};
  }
  iodata.boundaries.farfield.attributes = {2};  // MakeCartesian2D: x = W
  iodata.solver.order = 2;
  iodata.solver.boundary_mode.freq = 1.0;
  iodata.solver.boundary_mode.n = 1;

  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(8, 4, mfem::Element::TRIANGLE, false, W, H));
  for (int i = 0; i < serial_mesh->GetNE(); i++)
  {
    mfem::Vector center;
    serial_mesh->GetElementCenter(i, center);
    serial_mesh->SetAttribute(i, (center(1) < 0.25 * H) ? 1 : 2);
  }
  serial_mesh->SetAttributes();
  iodata.NondimensionalizeInputs(serial_mesh);
  Mesh palace_mesh(std::make_unique<mfem::ParMesh>(comm, *serial_mesh));
  iodata.CheckConfiguration();

  mfem::ND_FECollection nd_fec(iodata.solver.order, palace_mesh.Dimension());
  mfem::H1_FECollection h1_fec(iodata.solver.order, palace_mesh.Dimension());
  FiniteElementSpace nd_fespace(palace_mesh, &nd_fec);
  FiniteElementSpace h1_fespace(palace_mesh, &h1_fec);
  MaterialOperator mat_op(iodata, palace_mesh);
  SurfaceImpedanceOperator surf_z_op(iodata, mat_op, palace_mesh.Get());
  FarfieldBoundaryOperator farfield_op(iodata, mat_op, palace_mesh.Get());
  SurfaceConductivityOperator surf_sigma_op(iodata, mat_op, palace_mesh.Get());
  SurfaceRationalImpedanceOperator surf_rz_op(iodata, mat_op, palace_mesh.Get());

  auto [Atnr, Atni] = mode_assembly::AssembleAtn(nd_fespace, h1_fespace, mat_op);
  std::unique_ptr<mfem::HypreParMatrix> Btnr(Atnr->Transpose());
  *Btnr *= -1.0;
  auto [Bttr, Btti] = mode_assembly::AssembleBtt(nd_fespace, mat_op);
  mfem::Array<int> dbc_tdof_list;  // No essential boundaries
  mode_assembly::ModeOperatorModel model(nd_fespace, h1_fespace, mat_op, nullptr, surf_z_op,
                                         farfield_op, surf_sigma_op, surf_rz_op, *Bttr,
                                         Atnr.get(), Atni.get(), Btnr.get(), dbc_tdof_list);
  const auto &components = model.GetComponents();
  auto omega_component =
      std::find_if(components.begin(), components.end(), [](const auto &component)
                   { return component.type == mode_assembly::CoefficientType::OMEGA; });
  REQUIRE(omega_component != components.end());

  // Apply the component to [Eₜ; Eₙ] and return the pairing with the input. The damping
  // only has an imaginary part (the frequency scalar iω is applied separately).
  const int nd_size = nd_fespace.GetTrueVSize(), h1_size = h1_fespace.GetTrueVSize();
  auto Pairing = [&](const Vector &et, double en)
  {
    ComplexVector x(nd_size + h1_size), y(nd_size + h1_size);
    x = 0.0;
    x.Real().SetVector(et, 0);
    for (int i = nd_size; i < nd_size + h1_size; i++)
    {
      x.Real()[i] = en;
    }
    omega_component->op->Mult(x, y);
    CHECK(linalg::Norml2(comm, y.Real()) == 0.0);
    return linalg::Dot(comm, x.Real(), y.Imag());
  };
  const double expected = 2.0 * 0.25 * H + 1.0 * 0.75 * H;

  // In-plane field ŷ, tangential to the absorbing boundary.
  Vector et(nd_size);
  {
    mfem::ParGridFunction E(&nd_fespace.Get());
    mfem::Vector yhat(2);
    yhat(0) = 0.0;
    yhat(1) = 1.0;
    mfem::VectorConstantCoefficient coeff(yhat);
    E.ProjectCoefficient(coeff);
    E.GetTrueDofs(et);
  }
  CHECK_THAT(Pairing(et, 0.0), WithinRel(expected, 1.0e-12));

  // Unit out-of-plane field.
  Vector zero(nd_size);
  zero = 0.0;
  CHECK_THAT(Pairing(zero, 1.0), WithinRel(-expected, 1.0e-12));
}

}  // namespace palace
