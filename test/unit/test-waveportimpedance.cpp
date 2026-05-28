// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

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
TEST_CASE("WavePortData TE10 Z_PV matches analytical formula",
          "[waveportimpedance][Serial]")
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
  WavePortOperator wave_port_op(iodata, mat_op, nd_fespace_palace.Get(),
                                h1_fespace_palace.Get());
  wave_port_op.SetSuppressOutput(true);

  // Initialize the modes at the test frequency. Internal units: ω·tc. Initialize
  // is per-port; const_cast around the GetPort() const is safe here because we
  // know the wave_port_op is non-const in this scope.
  const double omega_nondim =
      2.0 * M_PI * iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(f_GHz);
  auto &port = const_cast<WavePortData &>(wave_port_op.GetPort(1));
  port.Initialize(omega_nondim);

  std::complex<double> Z_PV_nondim = port.GetCharacteristicImpedance();
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

}  // namespace palace
