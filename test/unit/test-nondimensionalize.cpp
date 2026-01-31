// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/iodata.hpp"

namespace palace
{
using namespace Catch;

// Test that verifies NondimensionalizeInputs behavior before refactoring.
// After refactoring to free functions, this test ensures no behavior changed.
TEST_CASE("Nondimensionalize via IoData", "[nondimensionalize][Serial]")
{
  // Use non-trivial units: L0 = 1e-3 (mm), Lc = 2.0 (in mesh units, so 2mm)
  // This gives Lc_m = 2e-3 m
  constexpr double L0 = 1e-3;
  constexpr double Lc = 2.0;
  Units units(L0, Lc * L0);

  IoData iodata(units);
  iodata.model.L0 = L0;
  iodata.model.Lc = Lc;

  // Fill in values that will be nondimensionalized

  // Material: conductivity and London penetration depth
  auto &mat = iodata.domains.materials.emplace_back();
  mat.attributes = {1};
  mat.sigma = config::SymmetricMatrixData<3>(1e6);  // 1e6 S/m
  mat.lambda_L = 0.1;                               // 0.1 mm (in mesh units)

  // Probe location
  iodata.domains.postpro.probe[1].center = {1.0, 2.0, 3.0};  // in mesh units

  // Current dipole location
  iodata.domains.current_dipole[1].center = {4.0, 5.0, 6.0};

  // Conductivity boundary
  auto &cond = iodata.boundaries.conductivity.emplace_back();
  cond.attributes = {1};
  cond.sigma = 5e7;  // S/m
  cond.h = 0.01;     // thickness in mesh units

  // Impedance boundary
  auto &imp = iodata.boundaries.impedance.emplace_back();
  imp.attributes = {2};
  imp.Rs = 50.0;   // Ohms
  imp.Ls = 1e-9;   // H
  imp.Cs = 1e-12;  // F

  // Lumped port
  iodata.boundaries.lumpedport[1].R = 100.0;
  iodata.boundaries.lumpedport[1].L = 2e-9;
  iodata.boundaries.lumpedport[1].C = 2e-12;
  iodata.boundaries.lumpedport[1].Rs = 25.0;
  iodata.boundaries.lumpedport[1].Ls = 0.5e-9;
  iodata.boundaries.lumpedport[1].Cs = 0.5e-12;

  // Periodic boundary wave vector
  iodata.boundaries.periodic.wave_vector = {0.1, 0.2, 0.3};

  // Wave port offset
  iodata.boundaries.waveport[1].d_offset = 0.5;  // mesh units

  // Surface flux center
  iodata.boundaries.postpro.flux[1].center = {7.0, 8.0, 9.0};

  // Dielectric interface thickness
  iodata.boundaries.postpro.dielectric[1].t = 0.001;  // mesh units

  // Eigenmode frequencies
  iodata.solver.eigenmode.target = 5.0;         // GHz
  iodata.solver.eigenmode.target_upper = 10.0;  // GHz

  // Driven frequencies
  iodata.solver.driven.sample_f = {1.0, 2.0, 3.0};  // GHz

  // Transient parameters
  iodata.solver.transient.pulse_f = 2.5;    // GHz
  iodata.solver.transient.pulse_tau = 1.0;  // ns
  iodata.solver.transient.max_t = 10.0;     // ns
  iodata.solver.transient.delta_t = 0.1;    // ns

  // Create a simple mesh to pass to NondimensionalizeInputs
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::TETRAHEDRON, Lc, Lc, Lc));
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);

  // Call the original nondimensionalization
  iodata.NondimensionalizeInputs(*par_mesh);

  // Now verify all the nondimensionalized values
  const double mesh_scale = units.GetMeshLengthRelativeScale();
  const double sigma_scale = units.GetScaleFactor<Units::ValueType::CONDUCTIVITY>();
  const double Z_scale = units.GetScaleFactor<Units::ValueType::IMPEDANCE>();
  const double L_scale = units.GetScaleFactor<Units::ValueType::INDUCTANCE>();
  const double C_scale = units.GetScaleFactor<Units::ValueType::CAPACITANCE>();
  const double t_scale = units.GetScaleFactor<Units::ValueType::TIME>();
  const double f_scale = units.GetScaleFactor<Units::ValueType::FREQUENCY>();

  // Material
  CHECK(mat.sigma.s[0] == Approx(1e6 / sigma_scale));
  CHECK(mat.lambda_L == Approx(0.1 / mesh_scale));

  // Probe location (divided by mesh_scale)
  CHECK(iodata.domains.postpro.probe[1].center[0] == Approx(1.0 / mesh_scale));
  CHECK(iodata.domains.postpro.probe[1].center[1] == Approx(2.0 / mesh_scale));
  CHECK(iodata.domains.postpro.probe[1].center[2] == Approx(3.0 / mesh_scale));

  // Current dipole location
  CHECK(iodata.domains.current_dipole[1].center[0] == Approx(4.0 / mesh_scale));
  CHECK(iodata.domains.current_dipole[1].center[1] == Approx(5.0 / mesh_scale));
  CHECK(iodata.domains.current_dipole[1].center[2] == Approx(6.0 / mesh_scale));

  // Conductivity boundary
  CHECK(cond.sigma == Approx(5e7 / sigma_scale));
  CHECK(cond.h == Approx(0.01 / mesh_scale));

  // Impedance boundary
  CHECK(imp.Rs == Approx(50.0 / Z_scale));
  CHECK(imp.Ls == Approx(1e-9 / L_scale));
  CHECK(imp.Cs == Approx(1e-12 / C_scale));

  // Lumped port
  CHECK(iodata.boundaries.lumpedport[1].R == Approx(100.0 / Z_scale));
  CHECK(iodata.boundaries.lumpedport[1].L == Approx(2e-9 / L_scale));
  CHECK(iodata.boundaries.lumpedport[1].C == Approx(2e-12 / C_scale));
  CHECK(iodata.boundaries.lumpedport[1].Rs == Approx(25.0 / Z_scale));
  CHECK(iodata.boundaries.lumpedport[1].Ls == Approx(0.5e-9 / L_scale));
  CHECK(iodata.boundaries.lumpedport[1].Cs == Approx(0.5e-12 / C_scale));

  // Periodic wave vector (multiplied by mesh_scale)
  CHECK(iodata.boundaries.periodic.wave_vector[0] == Approx(0.1 * mesh_scale));
  CHECK(iodata.boundaries.periodic.wave_vector[1] == Approx(0.2 * mesh_scale));
  CHECK(iodata.boundaries.periodic.wave_vector[2] == Approx(0.3 * mesh_scale));

  // Wave port offset
  CHECK(iodata.boundaries.waveport[1].d_offset == Approx(0.5 / mesh_scale));

  // Surface flux center
  CHECK(iodata.boundaries.postpro.flux[1].center[0] == Approx(7.0 / mesh_scale));
  CHECK(iodata.boundaries.postpro.flux[1].center[1] == Approx(8.0 / mesh_scale));
  CHECK(iodata.boundaries.postpro.flux[1].center[2] == Approx(9.0 / mesh_scale));

  // Dielectric interface thickness
  CHECK(iodata.boundaries.postpro.dielectric[1].t == Approx(0.001 / mesh_scale));

  // Eigenmode frequencies (converted to angular frequency)
  CHECK(iodata.solver.eigenmode.target ==
        Approx(2 * M_PI * units.Nondimensionalize<Units::ValueType::FREQUENCY>(5.0)));
  CHECK(iodata.solver.eigenmode.target_upper ==
        Approx(2 * M_PI * units.Nondimensionalize<Units::ValueType::FREQUENCY>(10.0)));

  // Driven frequencies
  CHECK(iodata.solver.driven.sample_f[0] ==
        Approx(2 * M_PI * units.Nondimensionalize<Units::ValueType::FREQUENCY>(1.0)));
  CHECK(iodata.solver.driven.sample_f[1] ==
        Approx(2 * M_PI * units.Nondimensionalize<Units::ValueType::FREQUENCY>(2.0)));
  CHECK(iodata.solver.driven.sample_f[2] ==
        Approx(2 * M_PI * units.Nondimensionalize<Units::ValueType::FREQUENCY>(3.0)));

  // Transient parameters
  CHECK(iodata.solver.transient.pulse_f ==
        Approx(2 * M_PI * units.Nondimensionalize<Units::ValueType::FREQUENCY>(2.5)));
  CHECK(iodata.solver.transient.pulse_tau ==
        Approx(units.Nondimensionalize<Units::ValueType::TIME>(1.0)));
  CHECK(iodata.solver.transient.max_t ==
        Approx(units.Nondimensionalize<Units::ValueType::TIME>(10.0)));
  CHECK(iodata.solver.transient.delta_t ==
        Approx(units.Nondimensionalize<Units::ValueType::TIME>(0.1)));
}

// Test that free functions produce identical results to IoData::NondimensionalizeInputs
TEST_CASE("Nondimensionalize free functions match IoData", "[nondimensionalize][Serial]")
{
  constexpr double L0 = 1e-3;
  constexpr double Lc = 2.0;
  Units units(L0, Lc * L0);

  // Test MaterialData
  {
    config::MaterialData data;
    data.sigma = config::SymmetricMatrixData<3>(1e6);
    data.lambda_L = 0.1;

    config::Nondimensionalize(units, data);

    CHECK(data.sigma.s[0] ==
          Approx(1e6 / units.GetScaleFactor<Units::ValueType::CONDUCTIVITY>()));
    CHECK(data.lambda_L == Approx(0.1 / units.GetMeshLengthRelativeScale()));
  }

  // Test ConductivityData
  {
    config::ConductivityData data;
    data.sigma = 5e7;
    data.h = 0.01;

    config::Nondimensionalize(units, data);

    CHECK(data.sigma ==
          Approx(5e7 / units.GetScaleFactor<Units::ValueType::CONDUCTIVITY>()));
    CHECK(data.h == Approx(0.01 / units.GetMeshLengthRelativeScale()));
  }

  // Test ImpedanceData
  {
    config::ImpedanceData data;
    data.Rs = 50.0;
    data.Ls = 1e-9;
    data.Cs = 1e-12;

    config::Nondimensionalize(units, data);

    CHECK(data.Rs == Approx(50.0 / units.GetScaleFactor<Units::ValueType::IMPEDANCE>()));
    CHECK(data.Ls == Approx(1e-9 / units.GetScaleFactor<Units::ValueType::INDUCTANCE>()));
    CHECK(data.Cs == Approx(1e-12 / units.GetScaleFactor<Units::ValueType::CAPACITANCE>()));
  }

  // Test LumpedPortData
  {
    config::LumpedPortData data;
    data.R = 100.0;
    data.L = 2e-9;
    data.C = 2e-12;
    data.Rs = 25.0;
    data.Ls = 0.5e-9;
    data.Cs = 0.5e-12;

    config::Nondimensionalize(units, data);

    CHECK(data.R == Approx(100.0 / units.GetScaleFactor<Units::ValueType::IMPEDANCE>()));
    CHECK(data.L == Approx(2e-9 / units.GetScaleFactor<Units::ValueType::INDUCTANCE>()));
    CHECK(data.C == Approx(2e-12 / units.GetScaleFactor<Units::ValueType::CAPACITANCE>()));
    CHECK(data.Rs == Approx(25.0 / units.GetScaleFactor<Units::ValueType::IMPEDANCE>()));
    CHECK(data.Ls == Approx(0.5e-9 / units.GetScaleFactor<Units::ValueType::INDUCTANCE>()));
    CHECK(data.Cs ==
          Approx(0.5e-12 / units.GetScaleFactor<Units::ValueType::CAPACITANCE>()));
  }

  // Test PeriodicBoundaryData
  {
    config::PeriodicBoundaryData data;
    data.wave_vector = {0.1, 0.2, 0.3};

    config::Nondimensionalize(units, data);

    CHECK(data.wave_vector[0] == Approx(0.1 * units.GetMeshLengthRelativeScale()));
    CHECK(data.wave_vector[1] == Approx(0.2 * units.GetMeshLengthRelativeScale()));
    CHECK(data.wave_vector[2] == Approx(0.3 * units.GetMeshLengthRelativeScale()));
  }

  // Test EigenSolverData
  {
    config::EigenSolverData data;
    data.target = 5.0;
    data.target_upper = 10.0;

    config::Nondimensionalize(units, data);

    CHECK(data.target ==
          Approx(2 * M_PI * units.Nondimensionalize<Units::ValueType::FREQUENCY>(5.0)));
    CHECK(data.target_upper ==
          Approx(2 * M_PI * units.Nondimensionalize<Units::ValueType::FREQUENCY>(10.0)));
  }

  // Test TransientSolverData
  {
    config::TransientSolverData data;
    data.pulse_f = 2.5;
    data.pulse_tau = 1.0;
    data.max_t = 10.0;
    data.delta_t = 0.1;

    config::Nondimensionalize(units, data);

    CHECK(data.pulse_f ==
          Approx(2 * M_PI * units.Nondimensionalize<Units::ValueType::FREQUENCY>(2.5)));
    CHECK(data.pulse_tau == Approx(units.Nondimensionalize<Units::ValueType::TIME>(1.0)));
    CHECK(data.max_t == Approx(units.Nondimensionalize<Units::ValueType::TIME>(10.0)));
    CHECK(data.delta_t == Approx(units.Nondimensionalize<Units::ValueType::TIME>(0.1)));
  }

  // Test RefinementData
  {
    config::RefinementData data;
    config::BoxRefinementData box;
    box.bbmin = {1.0, 2.0, 3.0};
    box.bbmax = {4.0, 5.0, 6.0};
    box.ref_levels = 1;
    data.GetBoxes().push_back(box);

    config::SphereRefinementData sphere;
    sphere.center = {7.0, 8.0, 9.0};
    sphere.r = 0.5;
    sphere.ref_levels = 2;
    data.GetSpheres().push_back(sphere);

    config::Nondimensionalize(units, data);

    const double Lc0 = units.GetMeshLengthRelativeScale();
    CHECK(data.GetBoxes()[0].bbmin[0] == Approx(1.0 / Lc0));
    CHECK(data.GetBoxes()[0].bbmin[1] == Approx(2.0 / Lc0));
    CHECK(data.GetBoxes()[0].bbmin[2] == Approx(3.0 / Lc0));
    CHECK(data.GetBoxes()[0].bbmax[0] == Approx(4.0 / Lc0));
    CHECK(data.GetBoxes()[0].bbmax[1] == Approx(5.0 / Lc0));
    CHECK(data.GetBoxes()[0].bbmax[2] == Approx(6.0 / Lc0));
    CHECK(data.GetSpheres()[0].center[0] == Approx(7.0 / Lc0));
    CHECK(data.GetSpheres()[0].center[1] == Approx(8.0 / Lc0));
    CHECK(data.GetSpheres()[0].center[2] == Approx(9.0 / Lc0));
    CHECK(data.GetSpheres()[0].r == Approx(0.5 / Lc0));
  }

  // Test ProbeData
  {
    config::ProbeData data;
    data.center = {1.0, 2.0, 3.0};

    config::Nondimensionalize(units, data);

    const double Lc0 = units.GetMeshLengthRelativeScale();
    CHECK(data.center[0] == Approx(1.0 / Lc0));
    CHECK(data.center[1] == Approx(2.0 / Lc0));
    CHECK(data.center[2] == Approx(3.0 / Lc0));
  }

  // Test CurrentDipoleData
  {
    config::CurrentDipoleData data;
    data.center = {4.0, 5.0, 6.0};

    config::Nondimensionalize(units, data);

    const double Lc0 = units.GetMeshLengthRelativeScale();
    CHECK(data.center[0] == Approx(4.0 / Lc0));
    CHECK(data.center[1] == Approx(5.0 / Lc0));
    CHECK(data.center[2] == Approx(6.0 / Lc0));
  }

  // Test WavePortData
  {
    config::WavePortData data;
    data.d_offset = 0.5;

    config::Nondimensionalize(units, data);

    CHECK(data.d_offset == Approx(0.5 / units.GetMeshLengthRelativeScale()));
  }

  // Test SurfaceFluxData
  {
    config::SurfaceFluxData data;
    data.center = {7.0, 8.0, 9.0};

    config::Nondimensionalize(units, data);

    const double Lc0 = units.GetMeshLengthRelativeScale();
    CHECK(data.center[0] == Approx(7.0 / Lc0));
    CHECK(data.center[1] == Approx(8.0 / Lc0));
    CHECK(data.center[2] == Approx(9.0 / Lc0));
  }

  // Test InterfaceDielectricData
  {
    config::InterfaceDielectricData data;
    data.t = 0.001;

    config::Nondimensionalize(units, data);

    CHECK(data.t == Approx(0.001 / units.GetMeshLengthRelativeScale()));
  }

  // Test DrivenSolverData
  {
    config::DrivenSolverData data;
    data.sample_f = {1.0, 2.0, 3.0};

    config::Nondimensionalize(units, data);

    CHECK(data.sample_f[0] ==
          Approx(2 * M_PI * units.Nondimensionalize<Units::ValueType::FREQUENCY>(1.0)));
    CHECK(data.sample_f[1] ==
          Approx(2 * M_PI * units.Nondimensionalize<Units::ValueType::FREQUENCY>(2.0)));
    CHECK(data.sample_f[2] ==
          Approx(2 * M_PI * units.Nondimensionalize<Units::ValueType::FREQUENCY>(3.0)));
  }
}

}  // namespace palace
