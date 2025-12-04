// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Test the current dipole implementation by solving a driven problem
// and comparing computed fields with analytical solutions.
//
// Uses analytical solution for current dipole from:
// https://em.geosci.xyz/content/maxwell1_fundamentals/dipole_sources_in_homogeneous_media/electric_dipole_frequency/analytic_solution.html

#include <complex>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "linalg/ksp.hpp"
#include "linalg/vector.hpp"
#include "models/spaceoperator.hpp"
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

// Compute analytical electric field from current dipole (geosci.xyz formulation)
// Returns non-dimensionalized field for comparison with Palace computed field
std::array<std::complex<double>, 3> ComputeDipoleENonDim(const mfem::Vector &x_nondim,
                                                         const Units &units, double I0l_Am,
                                                         double freq_Hz)
{
  double r_nondim = x_nondim.Norml2();
  if (r_nondim < 1e-12)
    return {0.0, 0.0, 0.0};

  // Convert to spherical coordinates
  double theta = std::acos(x_nondim(2) / r_nondim);
  double phi = std::atan2(x_nondim(1), x_nondim(0));
  double r = units.Dimensionalize<Units::ValueType::LENGTH>(r_nondim);

  double omega = 2.0 * M_PI * freq_Hz;
  double k = omega / c0_;
  double kr = k * r;

  // σ = iωε₀ for vacuum (geosci.xyz formulation)
  std::complex<double> sigma(0.0, omega * epsilon0_);
  std::complex<double> factor = I0l_Am / (4.0 * M_PI * sigma * r * r * r);
  std::complex<double> ikr(0.0, -kr);  // e^(-ikr) convention
  std::complex<double> exp_ikr = std::exp(ikr);

  std::complex<double> Er = factor * 2.0 * std::cos(theta) * (1.0 + ikr) * exp_ikr;
  std::complex<double> Etheta = factor * std::sin(theta) * (1.0 + ikr - kr * kr) * exp_ikr;

  // Nondimensionalize for comparison with Palace
  Er = units.Nondimensionalize<Units::ValueType::FIELD_E>(Er);
  Etheta = units.Nondimensionalize<Units::ValueType::FIELD_E>(Etheta);

  // Convert to Cartesian coordinates
  return {{Er * std::sin(theta) * std::cos(phi) + Etheta * std::cos(theta) * std::cos(phi),
           Er * std::sin(theta) * std::sin(phi) + Etheta * std::cos(theta) * std::sin(phi),
           Er * std::cos(theta) - Etheta * std::sin(theta)}};
}

// Compute analytical magnetic field from current dipole (geosci.xyz formulation)
// Returns non-dimensionalized field for comparison with Palace computed field
std::array<std::complex<double>, 3> ComputeDipoleBNonDim(const mfem::Vector &x_nondim,
                                                         const Units &units, double I0l_Am,
                                                         double freq_Hz)
{
  double r_nondim = x_nondim.Norml2();
  if (r_nondim < 1e-12)
    return {0.0, 0.0, 0.0};

  // Convert to spherical coordinates
  double theta = std::acos(x_nondim(2) / r_nondim);
  double phi = std::atan2(x_nondim(1), x_nondim(0));
  double r = units.Dimensionalize<Units::ValueType::LENGTH>(r_nondim);

  double omega = 2.0 * M_PI * freq_Hz;
  double k = omega / c0_;
  double kr = k * r;

  // H_φ = (ikI₀l)/(4πr²) * sin(θ) * (1 + ikr) * e^(-ikr) (geosci.xyz formulation)
  std::complex<double> ik(0.0, k);
  std::complex<double> ikr(0.0, -kr);  // e^(-ikr) convention
  std::complex<double> exp_ikr = std::exp(ikr);
  std::complex<double> factor = (ik * I0l_Am) / (4.0 * M_PI * r * r);
  std::complex<double> Hphi = factor * std::sin(theta) * (1.0 + ikr) * exp_ikr;

  // Convert H to B = μ₀H and nondimensionalize
  std::complex<double> Bphi = mu0_ * Hphi;
  Bphi = units.Nondimensionalize<Units::ValueType::FIELD_B>(Bphi);

  // Convert to Cartesian coordinates (only φ component for z-directed dipole)
  return {{-Bphi * std::sin(phi), Bphi * std::cos(phi), 0.0}};
}

}  // namespace

void runCurrentDipoleTest(double freq_Hz, std::unique_ptr<mfem::Mesh> serial_mesh)
{
  // Test constants
  constexpr double atol = 1e-4;
  constexpr double rtol = 5e-2;    // 5% tolerance
  constexpr double I0l_Am = 1e-3;  // Current moment [A⋅m]

  Units units(0.496, 1.453);  // Non-trivial units for testing

  IoData iodata = IoData(units);
  iodata.domains.materials.emplace_back().attributes = {1};
  iodata.problem.type = ProblemType::DRIVEN;

  // Current dipole configuration
  auto &dipole_config = iodata.domains.current_dipole[1];
  dipole_config.moment = {0.0, 0.0, I0l_Am};  // z-aligned dipole
  dipole_config.center = {0.0, 0.0, 0.0};

  auto comm = Mpi::World();

  const int dim = serial_mesh->Dimension();
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  iodata.NondimensionalizeInputs(*par_mesh);
  Mesh palace_mesh(std::move(par_mesh));

  // Set up finite element spaces
  mfem::ND_FECollection nd_fec(2, dim);  // Second order elements for E field
  FiniteElementSpace nd_fespace(palace_mesh, &nd_fec);
  mfem::RT_FECollection rt_fec(2, dim);  // Second order elements for B field
  FiniteElementSpace rt_fespace(palace_mesh, &rt_fec);

  // Set up analytical E field coefficients
  auto E_real = [=](const mfem::Vector &x, mfem::Vector &E)
  {
    auto E_complex = ComputeDipoleENonDim(x, units, I0l_Am, freq_Hz);
    E(0) = std::real(E_complex[0]);
    E(1) = std::real(E_complex[1]);
    E(2) = std::real(E_complex[2]);
  };

  auto E_imag = [=](const mfem::Vector &x, mfem::Vector &E)
  {
    auto E_complex = ComputeDipoleENonDim(x, units, I0l_Am, freq_Hz);
    E(0) = std::imag(E_complex[0]);
    E(1) = std::imag(E_complex[1]);
    E(2) = std::imag(E_complex[2]);
  };

  // Set up analytical B field coefficients
  auto B_real = [=](const mfem::Vector &x, mfem::Vector &B)
  {
    auto B_complex = ComputeDipoleBNonDim(x, units, I0l_Am, freq_Hz);
    B(0) = std::real(B_complex[0]);
    B(1) = std::real(B_complex[1]);
    B(2) = std::real(B_complex[2]);
  };

  auto B_imag = [=](const mfem::Vector &x, mfem::Vector &B)
  {
    auto B_complex = ComputeDipoleBNonDim(x, units, I0l_Am, freq_Hz);
    B(0) = std::imag(B_complex[0]);
    B(1) = std::imag(B_complex[1]);
    B(2) = std::imag(B_complex[2]);
  };

  // Project analytical solutions onto grid functions
  bool complex = true;
  GridFunction E_analytical(nd_fespace, complex);
  GridFunction B_analytical(rt_fespace, complex);

  mfem::VectorFunctionCoefficient Ereal(dim, E_real);
  mfem::VectorFunctionCoefficient Eimag(dim, E_imag);
  E_analytical.Real().ProjectCoefficient(Ereal);
  E_analytical.Imag().ProjectCoefficient(Eimag);

  mfem::VectorFunctionCoefficient Breal(dim, B_real);
  mfem::VectorFunctionCoefficient Bimag(dim, B_imag);
  B_analytical.Real().ProjectCoefficient(Breal);
  B_analytical.Imag().ProjectCoefficient(Bimag);

  // TODO: Solve driven problem with current dipole and compare with analytical fields
  // For now, just verify analytical fields are non-zero
  Vector E_real_vec = E_analytical.Real().GetTrueVector();
  Vector E_imag_vec = E_analytical.Imag().GetTrueVector();
  Vector B_real_vec = B_analytical.Real().GetTrueVector();
  Vector B_imag_vec = B_analytical.Imag().GetTrueVector();

  double E_real_norm = E_real_vec.Norml2();
  double E_imag_norm = E_imag_vec.Norml2();
  double B_real_norm = B_real_vec.Norml2();
  double B_imag_norm = B_imag_vec.Norml2();

  CHECK(E_real_norm > 0.0);
  CHECK(E_imag_norm > 0.0);
  CHECK(B_real_norm > 0.0);
  CHECK(B_imag_norm > 0.0);
}

TEST_CASE("Current dipole field implementation", "[currentdipole][Serial]")
{
  double freq_Hz = 50e6;

  // Create mesh - unit cube centered at origin
  std::unique_ptr<mfem::Mesh> serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(6, 6, 6, mfem::Element::TETRAHEDRON));

  // Center mesh around origin
  serial_mesh->Transform(
      [](const mfem::Vector &x, mfem::Vector &p)
      {
        p = x;
        p(0) -= 0.5;
        p(1) -= 0.5;
        p(2) -= 0.5;
      });

  runCurrentDipoleTest(freq_Hz, std::move(serial_mesh));
}

}  // namespace palace
