// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
#include "models/materialoperator.hpp"
#include "models/waveportoperator.hpp"
#include "utils/communication.hpp"
#include "utils/constants.hpp"
#include "utils/iodata.hpp"
#include "utils/labels.hpp"
#include "utils/units.hpp"

namespace palace
{
using namespace Catch::Matchers;

// Expose private WavePortOperator::Initialize so tests can drive a single port without
// needing the full simulation harness (and without const_cast around GetPort()).
struct ExposedWavePortOperator : public WavePortOperator
{
  using WavePortOperator::Initialize;
  using WavePortOperator::WavePortOperator;
};

// Validate WavePortData::GetCharacteristicImpedance against the analytical TE10
// rectangular-waveguide power-voltage impedance.
//
//   Z_PV(TE10) = 2 · (b/a) · Z_TE,    Z_TE = η₀ / √(1 − (fc/f)²),    fc = c / (2a)
//
// This is the standard convention for power-voltage impedance with the voltage
// path taken across the narrow dimension at the broadside center (where E is
// maximum). The factor of 2 over (b/a)·Z_TE comes from the integration of
// sin²(πy/a) over the broadside, contributing a/2 to the power integral.
//
// Reference: Pozar, "Microwave Engineering" (4th ed.), Sec. 3.3 / Table 3.1
// (equivalent voltage / current / impedance for waveguide modes).
TEST_CASE("WavePort TE10 Z_PV", "[waveportimpedance][Serial]")
{
  MPI_Comm comm = Mpi::World();

  // WR-90 rectangular waveguide, broadside a, narrow side b, length L. Vacuum-
  // filled. Frequency 10 GHz, well above the TE10 cutoff fc = c/(2a) ≈ 6.56 GHz.
  // Match the SolveRectangularModes pattern: mesh coords in meters, Units(1,1)
  // so the nondimensionalization is identity.
  const double a_m = 22.86e-3;  // broadside (along y)
  const double b_m = 10.16e-3;  // narrow side (along z)
  const double L_m = 10.0e-3;   // waveguide propagation length (along x)
  const double f_GHz = 10.0;

  // Build a 3D rectangular box. MakeCartesian3D face-attribute convention:
  //   1 = bottom (z=0), 2 = front (y=0), 3 = right (x=L),
  //   4 = back (y=ly),  5 = left (x=0),  6 = top (z=lz).
  // We place the wave port on face 5 (x=0); PEC on the other five faces.
  // Use a moderately fine mesh so the line integral converges within ~1%.
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(8, 8, 4, mfem::Element::TETRAHEDRON, L_m, a_m, b_m));

  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.model.L0 = 1.0;  // Mesh coordinates are in raw meters.
  iodata.model.Lc = 1.0;

  // Vacuum domain.
  auto &material = iodata.domains.materials.emplace_back();
  material.attributes = {1};
  material.epsilon_r.s = {1.0, 1.0, 1.0};
  material.mu_r.s = {1.0, 1.0, 1.0};

  // PEC walls (everything except the port face 5).
  iodata.boundaries.pec.attributes = {1, 2, 3, 4, 6};

  // Wave port at face 5 (x=0). Voltage path across the narrow dimension at the
  // broadside center: from (0, a/2, 0) to (0, a/2, b). At this y-location the
  // TE10 E-field amplitude E_y = E_0·sin(π/2) = E_0 (its peak), so the line
  // integral V = E_0·b directly. Coordinates in meters (matching mesh units).
  auto &wave = iodata.boundaries.waveport.try_emplace(1).first->second;
  wave.attributes = {5};
  wave.mode_idx = 1;
  wave.excitation = 0;
  wave.active = true;
  wave.voltage_path = {{0.0, a_m / 2.0, 0.0}, {0.0, a_m / 2.0, b_m}};
  wave.n_samples = 200;
  wave.eig_tol = 1.0e-8;
  wave.ksp_tol = 1.0e-8;
  wave.ksp_max_its = 100;
  // Mirror the formula used by IoData::CheckConfiguration for eigenmode.max_size.
  wave.max_size = std::max(2 * wave.mode_idx, wave.mode_idx + 15);

  iodata.solver.order = 2;
  iodata.solver.linear.tol = 1.0e-8;
  iodata.solver.linear.max_it = 200;
  iodata.problem.type = ProblemType::DRIVEN;

  // Nondimensionalize and build mesh / FE spaces.
  iodata.NondimensionalizeInputs(serial_mesh);
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  iodata.CheckConfiguration();
  Mesh palace_mesh(std::move(par_mesh));

  auto nd_fec =
      std::make_unique<mfem::ND_FECollection>(iodata.solver.order, palace_mesh.Dimension());
  auto h1_fec =
      std::make_unique<mfem::H1_FECollection>(iodata.solver.order, palace_mesh.Dimension());
  FiniteElementSpace nd_fespace_palace(palace_mesh, nd_fec.get());
  FiniteElementSpace h1_fespace_palace(palace_mesh, h1_fec.get());
  MaterialOperator mat_op(iodata, palace_mesh);

  // Construct the WavePortOperator (which builds the WavePortData internally).
  ExposedWavePortOperator wave_port_op(iodata, mat_op, nd_fespace_palace.Get(),
                                       h1_fespace_palace.Get());
  wave_port_op.SetSuppressOutput(true);

  // Initialize the modes at the test frequency. Internal units: ω·tc.
  const double omega_nondim =
      2.0 * M_PI * iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(f_GHz);
  wave_port_op.Initialize(omega_nondim);

  std::complex<double> Z_PV_nondim = wave_port_op.GetPort(1).GetCharacteristicImpedance();
  // Dimensionalize: Z_PV in physical Ω = Z_PV_nondim · Z₀ (Z₀ = 376.73 Ω is the
  // internal-to-physical impedance scale factor).
  const double Z_PV_ohm = Z_PV_nondim.real() * electromagnetics::Z0_;

  // Analytical TE10 Z_PV:
  const double c0 = electromagnetics::c0_;
  const double f_Hz = f_GHz * 1.0e9;
  const double fc_Hz = c0 / (2.0 * a_m);
  const double Z_TE = electromagnetics::Z0_ / std::sqrt(1.0 - std::pow(fc_Hz / f_Hz, 2));
  const double Z_PV_analytical = 2.0 * (b_m / a_m) * Z_TE;

  CAPTURE(Z_PV_ohm, Z_PV_analytical);

  CHECK_THAT(Z_PV_ohm, WithinRel(Z_PV_analytical, 1.0e-3));
}

// Validate WavePortData::GetModePolaritySign returns the expected sign for the TE10
// mode in a rectangular waveguide. We use distinct PEC attributes for the y=0 and
// y=a broadside walls. The convention is (high, low) = (signal terminal, ground
// terminal) — E should point from high to low. TE10 has E = E_0·sin(πy/a)·ŷ which
// is +y-aligned, so naming y=0 as "signal" (high) and y=a as "ground" (low) gives a
// (high → low) direction of +y, matching the field, sign = +1. Swapping the two
// attribute roles flips the expected sign to -1.
TEST_CASE("WavePort TE10 mode polarity sign", "[waveportimpedance][Serial]")
{
  MPI_Comm comm = Mpi::World();

  const double a_m = 22.86e-3, b_m = 10.16e-3, L_m = 10.0e-3;
  const double f_GHz = 10.0;

  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(8, 8, 4, mfem::Element::TETRAHEDRON, L_m, a_m, b_m));

  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.model.L0 = 1.0;
  iodata.model.Lc = 1.0;

  auto &material = iodata.domains.materials.emplace_back();
  material.attributes = {1};
  material.epsilon_r.s = {1.0, 1.0, 1.0};
  material.mu_r.s = {1.0, 1.0, 1.0};

  // PEC walls (everything except port at face 5). Same as the Z_PV test; we don't
  // need separate physical attributes here because GetModePolaritySign matches
  // against the parent boundary attribute, which is already distinct for each face
  // of MakeCartesian3D.
  iodata.boundaries.pec.attributes = {1, 2, 3, 4, 6};

  auto &wave = iodata.boundaries.waveport.try_emplace(1).first->second;
  wave.attributes = {5};
  wave.mode_idx = 1;
  wave.excitation = 0;
  wave.active = true;
  wave.eig_tol = 1.0e-8;
  wave.ksp_tol = 1.0e-8;
  wave.ksp_max_its = 100;
  // Mirror the formula used by IoData::CheckConfiguration for eigenmode.max_size.
  wave.max_size = std::max(2 * wave.mode_idx, wave.mode_idx + 15);

  iodata.solver.order = 2;
  iodata.solver.linear.tol = 1.0e-8;
  iodata.solver.linear.max_it = 200;
  iodata.problem.type = ProblemType::DRIVEN;

  iodata.NondimensionalizeInputs(serial_mesh);
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  iodata.CheckConfiguration();
  Mesh palace_mesh(std::move(par_mesh));

  auto nd_fec =
      std::make_unique<mfem::ND_FECollection>(iodata.solver.order, palace_mesh.Dimension());
  auto h1_fec =
      std::make_unique<mfem::H1_FECollection>(iodata.solver.order, palace_mesh.Dimension());
  FiniteElementSpace nd_fespace_palace(palace_mesh, nd_fec.get());
  FiniteElementSpace h1_fespace_palace(palace_mesh, h1_fec.get());
  MaterialOperator mat_op(iodata, palace_mesh);

  ExposedWavePortOperator wave_port_op(iodata, mat_op, nd_fespace_palace.Get(),
                                       h1_fespace_palace.Get());
  wave_port_op.SetSuppressOutput(true);

  const double omega_nondim =
      2.0 * M_PI * iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(f_GHz);
  wave_port_op.Initialize(omega_nondim);

  // MakeCartesian3D face attributes: 2 = y=0 (broadside), 4 = y=a (broadside). TE10
  // has E = E_0·sin(πy/a)·ŷ pointing in +y. With (high, low) = (y=0, y=a) the
  // expected (high → low) direction is +y, matching the field → sign = +1. Swapping
  // gives (high → low) = -y, opposite the field → sign = -1.
  const auto &port = wave_port_op.GetPort(1);
  CHECK(port.GetModePolaritySign(/*high_attr=*/2, /*low_attr=*/4) == +1);
  CHECK(port.GetModePolaritySign(/*high_attr=*/4, /*low_attr=*/2) == -1);
}

TEST_CASE("WavePort anisotropic mode ordering", "[waveportimpedance][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();

  // WR-28 guide propagating along z. The dominant TE10 mode has E along y and H along x.
  const double a_m = 7.112e-3;
  const double b_m = 3.556e-3;
  const double L_m = 8.0e-3;
  const double f_GHz = 20.0;

  auto solve_kn =
      [&](const std::array<double, 3> &mu_r, const std::array<double, 3> &epsilon_r)
  {
    auto serial_mesh = std::make_unique<mfem::Mesh>(
        mfem::Mesh::MakeCartesian3D(8, 4, 2, mfem::Element::TETRAHEDRON, a_m, b_m, L_m));

    Units units(1.0, 1.0);
    IoData iodata(units);
    iodata.model.L0 = 1.0;
    iodata.model.Lc = 1.0;

    auto &material = iodata.domains.materials.emplace_back();
    material.attributes = {1};
    material.mu_r.s = mu_r;
    material.epsilon_r.s = epsilon_r;

    // MakeCartesian3D face 1 is z=0; use it as the port and make all other faces PEC.
    iodata.boundaries.pec.attributes = {2, 3, 4, 5, 6};
    auto &wave = iodata.boundaries.waveport.try_emplace(1).first->second;
    wave.attributes = {1};
    wave.mode_idx = 1;
    wave.active = true;
    wave.eig_tol = 1.0e-9;
    wave.ksp_tol = 1.0e-9;
    wave.ksp_max_its = 200;
    wave.max_size = std::max(2 * wave.mode_idx, wave.mode_idx + 15);

    iodata.solver.order = 1;
    iodata.solver.linear.tol = 1.0e-9;
    iodata.solver.linear.max_it = 200;
    iodata.problem.type = ProblemType::DRIVEN;

    iodata.NondimensionalizeInputs(serial_mesh);
    auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
    iodata.CheckConfiguration();
    Mesh palace_mesh(std::move(par_mesh));

    auto nd_fec = std::make_unique<mfem::ND_FECollection>(iodata.solver.order,
                                                          palace_mesh.Dimension());
    auto h1_fec = std::make_unique<mfem::H1_FECollection>(iodata.solver.order,
                                                          palace_mesh.Dimension());
    FiniteElementSpace nd_fespace(palace_mesh, nd_fec.get());
    FiniteElementSpace h1_fespace(palace_mesh, h1_fec.get());
    MaterialOperator mat_op(iodata, palace_mesh);

    ExposedWavePortOperator wave_port_op(iodata, mat_op, nd_fespace.Get(),
                                         h1_fespace.Get());
    wave_port_op.SetSuppressOutput(true);
    const double omega =
        2.0 * M_PI * iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(f_GHz);
    wave_port_op.Initialize(omega);

    const double kn_scale =
        1.0 / iodata.units.Dimensionalize<Units::ValueType::LENGTH>(1.0);
    return wave_port_op.GetPort(1).kn0 * kn_scale;
  };

  const double k0 = 2.0 * M_PI * f_GHz * 1.0e9 / electromagnetics::c0_;
  const double kc = M_PI / a_m;
  auto te10_kn =
      [&](const std::array<double, 3> &mu_r, const std::array<double, 3> &epsilon_r)
  { return std::sqrt(mu_r[0] * epsilon_r[1] * k0 * k0 - (mu_r[0] / mu_r[2]) * kc * kc); };

  const std::array<double, 3> mu_r = {1.0, 1.0, 1.0};
  const std::array<double, 3> epsilon_r = {2.0, 4.0, 3.0};
  const auto kn = solve_kn(mu_r, epsilon_r);
  const double kn_ref = te10_kn(mu_r, epsilon_r);
  CAPTURE(kn, kn_ref);
  CHECK_THAT(kn.real(), WithinRel(kn_ref, 0.02));
  CHECK_THAT(kn.imag(), WithinAbs(0.0, 0.01 * kn_ref));
}

// Validate the COMPLEX-frequency wave-port cross-section solve
// (WavePortData::SolveKnComplex) against the closed-form TE10 propagation constant,
// evaluated at a genuinely complex frequency ω. For a lossless homogeneous guide k_n(ω) =
// √(εμ ω² − k_c²), k_c = π/a, is an EXACT analytic function of ω. The whole
// complex-wave-port path (used by the eigenmode nonlinear solve at ω = -i·λ) rests on
// SolveKnComplex reproducing this analytic continuation off the real axis — i.e. k_n(ω) for
// complex ω must equal √(εμ ω² − k_c²) to discretization error, NOT the real-ω value k_n(Re
// ω). This pins that down directly, independent of any cavity / eigenmode / Q
// considerations.
TEST_CASE("WavePortData TE10 at complex ω", "[waveportimpedance][Serial]")
{
  MPI_Comm comm = Mpi::World();

  // Same WR-90 vacuum guide as the Z_PV test. TE10 cutoff fc = c/(2a) ≈ 6.56 GHz.
  const double a_m = 22.86e-3;  // broadside (along y) — sets k_c = π/a
  const double b_m = 10.16e-3;  // narrow side (along z)
  const double L_m = 10.0e-3;   // propagation length (along x)

  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(8, 8, 4, mfem::Element::TETRAHEDRON, L_m, a_m, b_m));

  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.model.L0 = 1.0;  // Mesh coordinates are in raw meters → k_n_nondim is in rad/m.
  iodata.model.Lc = 1.0;

  auto &material = iodata.domains.materials.emplace_back();
  material.attributes = {1};
  material.epsilon_r.s = {1.0, 1.0, 1.0};
  material.mu_r.s = {1.0, 1.0, 1.0};

  iodata.boundaries.pec.attributes = {1, 2, 3, 4, 6};

  auto &wave = iodata.boundaries.waveport.try_emplace(1).first->second;
  wave.attributes = {5};
  wave.mode_idx = 1;
  wave.excitation = 0;
  wave.active = true;
  wave.n_samples = 200;
  wave.eig_tol = 1.0e-10;
  wave.ksp_tol = 1.0e-10;
  wave.ksp_max_its = 200;
  // Mirror the formula used by IoData::CheckConfiguration for eigenmode.max_size.
  wave.max_size = std::max(2 * wave.mode_idx, wave.mode_idx + 15);

  iodata.solver.order = 2;
  iodata.solver.linear.tol = 1.0e-10;
  iodata.solver.linear.max_it = 200;
  iodata.problem.type = ProblemType::DRIVEN;

  iodata.NondimensionalizeInputs(serial_mesh);
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  iodata.CheckConfiguration();
  Mesh palace_mesh(std::move(par_mesh));

  auto nd_fec =
      std::make_unique<mfem::ND_FECollection>(iodata.solver.order, palace_mesh.Dimension());
  auto h1_fec =
      std::make_unique<mfem::H1_FECollection>(iodata.solver.order, palace_mesh.Dimension());
  FiniteElementSpace nd_fespace_palace(palace_mesh, nd_fec.get());
  FiniteElementSpace h1_fespace_palace(palace_mesh, h1_fec.get());
  MaterialOperator mat_op(iodata, palace_mesh);

  WavePortOperator wave_port_op(iodata, mat_op, nd_fespace_palace.Get(),
                                h1_fespace_palace.Get());
  wave_port_op.SetSuppressOutput(true);
  auto &port = const_cast<WavePortData &>(wave_port_op.GetPort(1));

  // Closed-form TE10 propagation constant k_n(ω) = √(εμ ω² − k_c²). In this Units(1,1)
  // / Lc = 1 m setup the internal angular frequency is ω·tc with tc = Lc/c0, and k_n comes
  // back in rad/m. Build the analytic reference directly in physical (rad/m) units and
  // compare against the dimensionalized SolveKnComplex result k_n_phys = k_n_nondim ·
  // (1/Lc).
  const double c0 = electromagnetics::c0_;
  const double kc = M_PI / a_m;  // TE10 transverse cutoff wavenumber [rad/m]
  auto kn_closed_form = [&](std::complex<double> omega_rad_s) -> std::complex<double>
  {
    // k_n = √((ω/c)² − k_c²); principal branch gives Re(k_n) ≥ 0 (forward sheet).
    const std::complex<double> k0 = omega_rad_s / c0;
    return std::sqrt(k0 * k0 - kc * kc);
  };

  // Helper: nondimensional ω for a physical frequency f [GHz].
  auto omega_nd = [&](double f_GHz)
  {
    return 2.0 * M_PI * iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(f_GHz);
  };
  // k_n scale factor (nondim → rad/m): kc_len = 1/Lc(meters).
  const double kn_scale = 1.0 / iodata.units.Dimensionalize<Units::ValueType::LENGTH>(1.0);

  SECTION("real ω above cutoff")
  {
    const double f_GHz = 10.0;
    std::complex<double> kn_nd =
        port.SolveKnComplex(std::complex<double>(omega_nd(f_GHz), 0.0));
    const std::complex<double> kn_phys = kn_nd * kn_scale;
    const std::complex<double> kn_ref =
        kn_closed_form(2.0 * M_PI * f_GHz * 1.0e9);  // ≈ 158.24 rad/m, ~0 imag
    CAPTURE(kn_phys, kn_ref);
    CHECK_THAT(kn_phys.real(), WithinRel(kn_ref.real(), 1.0e-4));
    CHECK_THAT(kn_phys.imag(), WithinAbs(0.0, 1.0e-4 * kn_ref.real()));
  }

  SECTION("complex ω")
  {
    // Probe a genuinely complex frequency: f = 10 GHz with a 5% imaginary part (a Q≈10
    // quasinormal-mode-like point, ω = ω_r(1 + i/20)). The exact solve MUST track the
    // analytic continuation here — this is precisely where a real-ω evaluation (k_n at
    // Re ω) would differ by O(Im ω / Re ω) ~ 1/(2Q).
    const double f_r_GHz = 10.0;
    const std::complex<double> scale_c(1.0, 0.05);
    std::complex<double> kn_nd = port.SolveKnComplex(omega_nd(f_r_GHz) * scale_c);
    const std::complex<double> kn_phys = kn_nd * kn_scale;
    const std::complex<double> omega_rad_s = 2.0 * M_PI * f_r_GHz * 1.0e9 * scale_c;
    const std::complex<double> kn_ref = kn_closed_form(omega_rad_s);
    CAPTURE(kn_phys, kn_ref);
    // Both real and imaginary parts must match the closed-form analytic continuation.
    CHECK_THAT(kn_phys.real(), WithinRel(kn_ref.real(), 1.0e-4));
    CHECK_THAT(kn_phys.imag(), WithinRel(kn_ref.imag(), 1.0e-4));
    // And the imaginary part must be genuinely nonzero (this is the whole point — the
    // complex-ω solve sees the off-axis dispersion that real-ω evaluation cannot).
    CHECK(std::abs(kn_ref.imag()) > 1.0);
  }

  SECTION("complex ω near cutoff")
  {
    // f_r just above cutoff (≈6.56 GHz) with loss: this is the branch-sensitive regime.
    // The principal-branch √ must give Re(k_n) ≥ 0 (forward / decaying sheet), matching the
    // closed form — guards against the wrong-Riemann-sheet hazard.
    const double f_r_GHz = 7.0;
    const std::complex<double> scale_c(1.0, 0.08);
    std::complex<double> kn_nd = port.SolveKnComplex(omega_nd(f_r_GHz) * scale_c);
    const std::complex<double> kn_phys = kn_nd * kn_scale;
    const std::complex<double> kn_ref =
        kn_closed_form(2.0 * M_PI * f_r_GHz * 1.0e9 * scale_c);
    CAPTURE(kn_phys, kn_ref);
    CHECK(kn_phys.real() >= 0.0);  // forward / decaying sheet
    CHECK_THAT(kn_phys.real(), WithinRel(kn_ref.real(), 1.0e-4));
    CHECK_THAT(kn_phys.imag(), WithinRel(kn_ref.imag(), 1.0e-4));
  }
}

// Validate the wave-port cross-section solve for a guide with a lossy (finite domain
// conductivity σ) fill against the closed-form dispersion with complex permittivity,
//
//   k_n(ω) = √(ω²μ(ε − iσ/ω) − k_c²) = √(ω²με − iωμσ − k_c²),
//
// which holds for BOTH the TE and TM modes of a homogeneous rectangular guide. The
// conduction current enters the 2D mode solver's two blocks with opposite signs (+iωσ in
// the transverse Att block, −iωσ in the normal Ann block, following the ∓ω²ε_c sign
// convention). A TE mode has no normal field component (e_n = 0), so it exercises only
// Att; a TM mode has e_n ≠ 0 and requires the Ann σ terms to be assembled consistently
// with Att — with σ missing from Ann (the historical gap), the degenerate TE11/TM11 pair
// splits into one correct root and one root missing its conduction loss, while TE10 stays
// correct. Checking BOTH members of the pair against the same closed form therefore pins
// the Ann σ assembly and the Att/Ann consistency regardless of how the degenerate pair is
// ordered. The complex-ω section additionally covers the Im(ω)·σ cross terms on the
// opposite slots (Re(±iωσ) = ∓wi·σ).
TEST_CASE("WavePortData lossy fill at real and complex ω", "[waveportimpedance][Serial]")
{
  MPI_Comm comm = Mpi::World();

  // WR-90 guide, uniformly filled with a weakly conducting dielectric. At 20 GHz,
  // σ/(ωε₀) ≈ 0.09, so the conduction term is a ~9% effect — far above discretization
  // error but small enough that the mode structure is unchanged. Cutoffs: TE10 6.56,
  // TE20 13.11, TE01 14.75, TE11/TM11 16.15 GHz (degenerate pair, modes 4 and 5).
  const double a_m = 22.86e-3;
  const double b_m = 10.16e-3;
  const double L_m = 10.0e-3;
  const double sigma_Sm = 0.1;
  const double f_GHz = 20.0;

  // The TM11 transverse wavenumber k_c = π√(1/a² + 1/b²) ≈ 338 rad/m is 2.5x the TE10
  // one; refine the port face beyond the TE10 tests' 8x8x4 for comparable accuracy.
  auto build_kn = [&](int mode_idx, std::complex<double> omega_scale)
  {
    auto serial_mesh = std::make_unique<mfem::Mesh>(
        mfem::Mesh::MakeCartesian3D(4, 14, 8, mfem::Element::TETRAHEDRON, L_m, a_m, b_m));

    Units units(1.0, 1.0);
    IoData iodata(units);
    iodata.model.L0 = 1.0;
    iodata.model.Lc = 1.0;

    auto &material = iodata.domains.materials.emplace_back();
    material.attributes = {1};
    material.epsilon_r.s = {1.0, 1.0, 1.0};
    material.mu_r.s = {1.0, 1.0, 1.0};
    material.sigma.s = {sigma_Sm, sigma_Sm, sigma_Sm};

    iodata.boundaries.pec.attributes = {1, 2, 3, 4, 6};

    auto &wave = iodata.boundaries.waveport.try_emplace(1).first->second;
    wave.attributes = {5};
    wave.mode_idx = mode_idx;
    wave.excitation = 0;
    wave.active = true;
    wave.eig_tol = 1.0e-10;
    wave.ksp_tol = 1.0e-10;
    wave.ksp_max_its = 200;
    // Mirror the formula used by IoData::CheckConfiguration for eigenmode.max_size.
    wave.max_size = std::max(2 * wave.mode_idx, wave.mode_idx + 15);

    iodata.solver.order = 2;
    iodata.solver.linear.tol = 1.0e-10;
    iodata.solver.linear.max_it = 200;
    iodata.problem.type = ProblemType::DRIVEN;

    iodata.NondimensionalizeInputs(serial_mesh);
    auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
    iodata.CheckConfiguration();
    Mesh palace_mesh(std::move(par_mesh));

    auto nd_fec = std::make_unique<mfem::ND_FECollection>(iodata.solver.order,
                                                          palace_mesh.Dimension());
    auto h1_fec = std::make_unique<mfem::H1_FECollection>(iodata.solver.order,
                                                          palace_mesh.Dimension());
    FiniteElementSpace nd_fespace_palace(palace_mesh, nd_fec.get());
    FiniteElementSpace h1_fespace_palace(palace_mesh, h1_fec.get());
    MaterialOperator mat_op(iodata, palace_mesh);

    WavePortOperator wave_port_op(iodata, mat_op, nd_fespace_palace.Get(),
                                  h1_fespace_palace.Get());
    wave_port_op.SetSuppressOutput(true);
    auto &port = const_cast<WavePortData &>(wave_port_op.GetPort(1));

    const double omega_nd =
        2.0 * M_PI * iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(f_GHz);
    const std::complex<double> kn_nd = port.SolveKnComplex(omega_nd * omega_scale);
    const double kn_scale =
        1.0 / iodata.units.Dimensionalize<Units::ValueType::LENGTH>(1.0);
    return kn_nd * kn_scale;  // physical rad/m
  };

  const double c0 = electromagnetics::c0_;
  const double mu0 = 4.0e-7 * M_PI;
  auto kn_closed_form = [&](double kc, std::complex<double> omega_rad_s)
  {
    const std::complex<double> k0 = omega_rad_s / c0;
    return std::sqrt(
        k0 * k0 - std::complex<double>(0.0, 1.0) * omega_rad_s * mu0 * sigma_Sm - kc * kc);
  };
  const double kc_te10 = M_PI / a_m;
  const double kc_tm11 = M_PI * std::sqrt(1.0 / (a_m * a_m) + 1.0 / (b_m * b_m));
  const std::complex<double> omega_r(2.0 * M_PI * f_GHz * 1.0e9, 0.0);

  SECTION("TE10 (transverse only: validates Att +iωσ)")
  {
    const std::complex<double> kn = build_kn(1, {1.0, 0.0});
    const std::complex<double> kn_ref = kn_closed_form(kc_te10, omega_r);
    CAPTURE(kn, kn_ref);
    CHECK_THAT(kn.real(), WithinRel(kn_ref.real(), 1.0e-3));
    CHECK_THAT(kn.imag(), WithinRel(kn_ref.imag(), 1.0e-3));
    // The conduction loss must actually be there (Im(k_n) ≈ -19.9 rad/m at 20 GHz).
    CHECK(std::abs(kn.imag()) > 10.0);
  }

  SECTION("TE11/TM11 pair (normal component: validates Ann -iωσ and Att/Ann consistency)")
  {
    // Both members of the degenerate pair must satisfy the same closed form
    // (kn ≈ 249.4 - 31.7i rad/m); a σ term missing from Ann breaks exactly one of them.
    const std::complex<double> kn_ref = kn_closed_form(kc_tm11, omega_r);
    for (int mode_idx : {4, 5})
    {
      const std::complex<double> kn = build_kn(mode_idx, {1.0, 0.0});
      CAPTURE(mode_idx, kn, kn_ref);
      CHECK_THAT(kn.real(), WithinRel(kn_ref.real(), 5.0e-3));
      CHECK_THAT(kn.imag(), WithinRel(kn_ref.imag(), 5.0e-3));
    }
  }

  SECTION("TE11/TM11 pair at complex ω (validates the Im(ω)·σ cross terms)")
  {
    // ω = ω_r(1 + i/20): the σ cross terms Re(±iωσ) = ∓wi·σ land on the real slots of
    // Att/Ann and only contribute off the real-ω axis.
    const std::complex<double> scale_c(1.0, 0.05);
    const std::complex<double> kn_ref = kn_closed_form(kc_tm11, omega_r * scale_c);
    for (int mode_idx : {4, 5})
    {
      const std::complex<double> kn = build_kn(mode_idx, scale_c);
      CAPTURE(mode_idx, kn, kn_ref);
      CHECK_THAT(kn.real(), WithinRel(kn_ref.real(), 5.0e-3));
      CHECK_THAT(kn.imag(), WithinRel(kn_ref.imag(), 5.0e-3));
    }
  }
}

}  // namespace palace
