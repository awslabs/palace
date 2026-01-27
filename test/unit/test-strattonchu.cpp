// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Test the near-to-far field transformations implemented in
// SurfacePostOperator.
//
// This test verifies the far-field computation by setting up a known dipole
// field and comparing the computed far-field against analytical expectations.
//
// More specifically, we set up a time-harmonic Hertzian dipole aligned on the z
// axis, with fields:
//
// E_r = (p₀/(4πε₀r³)) * 2cos(θ) * (1 - ikr) * e^(ikr)
// E_θ = (p₀/(4πε₀r³)) * sin(θ) * (1 - ikr - k²r²) * e^(ikr)
// E_φ = 0
//
// H_r = 0
// H_θ = 0
// H_φ = iω/(4πε₀r²) * sin(θ) * (1 - ikr) * e^(ikr)
//
// where k = 2π/λ = ω/c and p₀ is the dipole moment.

#include <complex>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <mfem/fem/coefficient.hpp>
#include <mfem/linalg/vector.hpp>
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "models/materialoperator.hpp"
#include "models/surfacepostoperator.hpp"
#include "utils/communication.hpp"
#include "utils/constants.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

namespace palace
{
using namespace Catch;
using namespace electromagnetics;
using namespace Catch::Matchers;

namespace
{

// Compute the Cartesian components electric field of a time-harmonic dipole
// aligned on the z axis. The returned field is non-dimensionalized according to
// the provided units.
//
// We want the non-dimensional E field to mock a Palace-computed field.
std::array<std::complex<double>, 3> ComputeDipoleENonDim(const mfem::Vector &x_nondim,
                                                         const Units &units, double p0,
                                                         double freq_Hz)
{
  double r_nondim = x_nondim.Norml2();
  if (r_nondim < 1e-12)
    return {0.0, 0.0, 0.0};

  // Convert to spherical coordinates.
  double theta = std::acos(x_nondim(2) / r_nondim);
  double phi = std::atan2(x_nondim(1), x_nondim(0));
  double r = units.Dimensionalize<Units::ValueType::LENGTH>(r_nondim);

  double k = 2.0 * M_PI * freq_Hz / c0_;  // Wave number k = ω/c = 2πf/c.
  double kr = k * r;

  double factor = p0 / (4.0 * M_PI * epsilon0_ * r * r * r);
  std::complex<double> jkr(0, kr);
  std::complex<double> exp_jkr = std::exp(jkr);

  std::complex<double> Er = factor * 2.0 * std::cos(theta) * (1.0 - jkr) * exp_jkr;
  std::complex<double> Etheta = factor * std::sin(theta) * (1.0 - jkr - kr * kr) * exp_jkr;

  Er = units.Nondimensionalize<Units::ValueType::FIELD_E>(Er);
  Etheta = units.Nondimensionalize<Units::ValueType::FIELD_E>(Etheta);

  return {{Er * std::sin(theta) * std::cos(phi) + Etheta * std::cos(theta) * std::cos(phi),
           Er * std::sin(theta) * std::sin(phi) + Etheta * std::cos(theta) * std::sin(phi),
           Er * std::cos(theta) - Etheta * std::sin(theta)}};
}

std::array<std::complex<double>, 3> ComputeDipoleBNonDim(const mfem::Vector &x_nondim,
                                                         const Units &units, double p0,
                                                         double freq_Hz)
{
  double r_nondim = x_nondim.Norml2();
  if (r_nondim < 1e-12)
    return {0.0, 0.0, 0.0};

  // Convert to spherical coordinates.
  double theta = std::acos(x_nondim(2) / r_nondim);
  double phi = std::atan2(x_nondim(1), x_nondim(0));

  double r = units.Dimensionalize<Units::ValueType::LENGTH>(r_nondim);

  double omega_rad_per_sec = 2 * M_PI * freq_Hz;
  double k = omega_rad_per_sec / c0_;
  double kr = k * r;

  // The magnetic field has only φ component for z-directed dipole.
  std::complex<double> factor(0, omega_rad_per_sec * mu0_ * p0 / (4.0 * M_PI * r * r));
  std::complex<double> jkr(0, kr);
  std::complex<double> exp_jkr = std::exp(jkr);
  std::complex<double> Bphi = factor * std::sin(theta) * (1.0 - jkr) * exp_jkr;

  Bphi = units.Nondimensionalize<Units::ValueType::FIELD_B>(Bphi);

  return {{-Bphi * std::sin(phi), Bphi * std::cos(phi), 0.0}};
}

// Take the analytic limit for kr >> 1 of the equations for E above and return
// rE / exp(ijk). Note, this is dimensional.
std::array<std::complex<double>, 3> ComputeAnalyticalFarFieldrE(double theta, double phi,
                                                                double p0, double freq_Hz)
{
  double k = 2.0 * M_PI * freq_Hz / c0_;  // Wave number k = ω/c = 2πf/c
  double factor = k * k * p0 / (4.0 * M_PI * epsilon0_);

  // The only component that survives is the transverse (as expected for a wave).
  std::complex<double> E_theta(factor * std::sin(theta), 0);

  // Return in Cartesian coordinates.
  return {{E_theta * std::cos(theta) * std::cos(phi),
           E_theta * std::cos(theta) * std::sin(phi), -E_theta * std::sin(theta)}};
}

// Generate uniformly distributed test points on unit sphere.
std::vector<std::pair<double, double>> GenerateSphericalTestPoints()
{
  // Test constants.
  constexpr int num_theta = 32;
  constexpr int num_phi = 32;

  std::vector<std::pair<double, double>> thetaphis;
  thetaphis.reserve(num_theta * num_phi);

  for (int i = 0; i < num_theta; ++i)
  {
    double theta = acos(1.0 - 2.0 * i / (num_theta - 1.0));
    for (int j = 0; j < num_phi; ++j)
    {
      double phi = 2.0 * M_PI * j / num_phi;
      thetaphis.emplace_back(theta, phi);
    }
  }
  return thetaphis;
}

// Compare the implementation in SurfacePostOperator with the analytic
// expectation in ComputeAnalyticalFarFieldrE. Note, the agreement has to be up
// to a phase, so we compare the magnitudes along each direction.
void runFarFieldTest(double freq_Hz, std::unique_ptr<mfem::Mesh> serial_mesh,
                     const std::vector<int> &attributes)
{
  // Test constants.
  constexpr double atol = 1e-4;
  constexpr double rtol = 5e-6;

  constexpr double p0 = 1e-9;  // Dipole moment [C⋅m]

  Units units(0.496, 1.453);  // Pick some arbitrary non-trivial units for testing

  IoData iodata = IoData(units);
  // We need to have at least one material. By default it's vacuum.
  iodata.domains.materials.emplace_back().attributes = {1};
  // We also need to have a non-empty farfield and a compatible problem type, or
  // the constructor for SurfacepostOperator will skip everything.
  iodata.boundaries.postpro.farfield.attributes = attributes;
  iodata.boundaries.postpro.farfield.thetaphis.emplace_back();
  iodata.problem.type = ProblemType::DRIVEN;

  auto comm = Mpi::World();

  // Read parallel mesh.
  const int dim = serial_mesh->Dimension();
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  iodata.NondimensionalizeInputs(*par_mesh);
  Mesh palace_mesh(std::move(par_mesh));

  // Set up complex fields on Nédélec element space.
  mfem::ND_FECollection nd_fec(3, dim);  // Third order elements
  FiniteElementSpace nd_fespace(palace_mesh, &nd_fec);

  // We need lambdas to set up the fields because of the mfem ProjectCoefficient
  // interface.
  auto E_real = [=](const mfem::Vector &x, mfem::Vector &E)
  {
    auto E_complex = ComputeDipoleENonDim(x, units, p0, freq_Hz);
    E(0) = std::real(E_complex[0]);
    E(1) = std::real(E_complex[1]);
    E(2) = std::real(E_complex[2]);
  };

  auto E_imag = [=](const mfem::Vector &x, mfem::Vector &E)
  {
    auto E_complex = ComputeDipoleENonDim(x, units, p0, freq_Hz);
    E(0) = std::imag(E_complex[0]);
    E(1) = std::imag(E_complex[1]);
    E(2) = std::imag(E_complex[2]);
  };

  auto B_real = [=](const mfem::Vector &x, mfem::Vector &B)
  {
    auto B_complex = ComputeDipoleBNonDim(x, units, p0, freq_Hz);
    B(0) = std::real(B_complex[0]);
    B(1) = std::real(B_complex[1]);
    B(2) = std::real(B_complex[2]);
  };

  auto B_imag = [=](const mfem::Vector &x, mfem::Vector &B)
  {
    auto B_complex = ComputeDipoleBNonDim(x, units, p0, freq_Hz);
    B(0) = std::imag(B_complex[0]);
    B(1) = std::imag(B_complex[1]);
    B(2) = std::imag(B_complex[2]);
  };

  bool complex = true;
  GridFunction E_field(nd_fespace, complex);
  mfem::VectorFunctionCoefficient Ereal(dim, E_real);
  mfem::VectorFunctionCoefficient Eimag(dim, E_imag);
  E_field.Real().ProjectCoefficient(Ereal);
  E_field.Imag().ProjectCoefficient(Eimag);

  GridFunction B_field(nd_fespace, complex);
  mfem::VectorFunctionCoefficient Breal(dim, B_real);
  mfem::VectorFunctionCoefficient Bimag(dim, B_imag);
  B_field.Real().ProjectCoefficient(Breal);
  B_field.Imag().ProjectCoefficient(Bimag);

  // Setup Palace operators for far-field computation.
  MaterialOperator mat_op(iodata, palace_mesh);
  SurfacePostOperator surf_post_op(iodata, mat_op, nd_fespace, nd_fespace);

  auto thetaphis = GenerateSphericalTestPoints();
  double omega_rad_per_time =
      2 * M_PI * units.Nondimensionalize<Units::ValueType::FREQUENCY>(freq_Hz / 1e9);
  constexpr double omega_im = 0.0;
  auto rE_computed =
      surf_post_op.GetFarFieldrE(thetaphis, E_field, B_field, omega_rad_per_time, omega_im);

  // Verify computed far-field against analytical solution
  for (size_t i = 0; i < thetaphis.size(); i++)
  {
    const auto &E_phys = units.Dimensionalize<Units::ValueType::VOLTAGE>(rE_computed[i]);
    const auto &[theta, phi] = thetaphis[i];
    auto rE_far = ComputeAnalyticalFarFieldrE(theta, phi, p0, freq_Hz);

    for (size_t j = 0; j < dim; j++)
    {
      // The agreement has to be up to a phase, so we compare the absolute values.
      CHECK_THAT(std::abs(E_phys[j]),
                 WithinRel(std::abs(rE_far[j]), rtol) || WithinAbs(0.0, atol));
    }
  }
}

TEST_CASE("Dipole field implementation", "[strattonchu][Serial]")
{
  // Test constants.
  constexpr double atol = 1e-4;
  constexpr double rtol = 1e-7;

  // Check that the ComputeAnalyticalFarFieldrE is the limit ComputeDipoleENonDim
  // by evaluating some points with large radii. The agreement has to be up to a
  // phase, so we check the absolute value.
  double p0 = 1e-9;
  double freq_Hz = 50e6;

  Units units(1.0, 1.0);
  double k = 2.0 * M_PI * freq_Hz / c0_;

  // Test at kr = 1000 (far-field regime) for various angles.
  double r = 1000.0 / k;
  std::vector<std::pair<double, double>> test_angles = {{M_PI / 6, 0.0},
                                                        {M_PI / 3, M_PI / 4},
                                                        {M_PI / 2, M_PI / 2},
                                                        {2 * M_PI / 3, 3 * M_PI / 4}};

  for (const auto &[theta, phi] : test_angles)
  {
    mfem::Vector x(3);
    x(0) = r * std::sin(theta) * std::cos(phi);
    x(1) = r * std::sin(theta) * std::sin(phi);
    x(2) = r * std::cos(theta);

    // Compare near-field and far-field with physical units.
    auto E_near = units.Dimensionalize<Units::ValueType::FIELD_E>(
        ComputeDipoleENonDim(x, units, p0, freq_Hz));
    auto rE_far = ComputeAnalyticalFarFieldrE(theta, phi, p0, freq_Hz);

    for (int i = 0; i < 3; i++)
    {
      std::complex<double> rE_near = r * E_near[i];
      CHECK_THAT(std::abs(rE_far[i]),
                 WithinRel(std::abs(rE_near), atol) || WithinAbs(0.0, rtol));
    }
  }
}

TEST_CASE("PostOperator", "[strattonchu][Serial]")
{
  // This test checks the Stratton-Chu code using a non-trivial mesh:
  // 1. The outer boundary has multiple attributes.
  // 2. The mesh is offset with respect to the source.
  // 2. The outer boundary is not a sphere.

  double freq_Hz = GENERATE(35e6, 50e6);
  std::vector<int> attributes = {1, 2, 3, 4, 5, 6};

  // Make mesh for a cube [0, 1] x [0, 1] x [0, 1].
  int resolution = 20;
  std::unique_ptr<mfem::Mesh> serial_mesh =
      std::make_unique<mfem::Mesh>(mfem::Mesh::MakeCartesian3D(
          resolution, resolution, resolution, mfem::Element::TETRAHEDRON));

  // Offset the cube a little bit.
  serial_mesh->Transform(
      [](const mfem::Vector &x, mfem::Vector &p)
      {
        p = x;
        p(0) -= 0.25;
        p(1) -= 0.25;
        p(2) -= 0.25;
      });

  runFarFieldTest(freq_Hz, std::move(serial_mesh), attributes);
}

TEST_CASE("FarField constructor fails with anisotropic materials", "[strattonchu][Serial]")
{
  Units units(0.496, 1.453);
  IoData iodata = IoData(units);

  auto &material = iodata.domains.materials.emplace_back();
  material.attributes = {1};
  material.mu_r.s[0] = 2;  // Make it anisotropic.

  iodata.boundaries.postpro.farfield.attributes = {1};
  iodata.boundaries.postpro.farfield.thetaphis.emplace_back();
  iodata.problem.type = ProblemType::DRIVEN;

  MPI_Comm comm = Mpi::World();
  std::unique_ptr<mfem::Mesh> serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::TETRAHEDRON));
  const int dim = serial_mesh->Dimension();
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);

  iodata.NondimensionalizeInputs(*par_mesh);
  Mesh palace_mesh(std::move(par_mesh));

  mfem::ND_FECollection nd_fec(1, dim);
  FiniteElementSpace nd_fespace(palace_mesh, &nd_fec);
  MaterialOperator mat_op(iodata, palace_mesh);

  CHECK_THROWS(SurfacePostOperator(iodata, mat_op, nd_fespace, nd_fespace));
}

}  // namespace
}  // namespace palace
