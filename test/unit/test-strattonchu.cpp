// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Test the near-to-far field transformations implemented in
// SurfacePostOperator.
//
// This test validates the far-field computation by setting up a known dipole
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
#include "drivers/drivensolver.hpp"
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/mesh.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/materialoperator.hpp"
#include "models/spaceoperator.hpp"
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
using namespace std::complex_literals;

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
  std::complex<double> exp_jkr = std::exp(-jkr);

  std::complex<double> Er = factor * 2.0 * std::cos(theta) * (1.0 + jkr) * exp_jkr;
  std::complex<double> Etheta = factor * std::sin(theta) * (1.0 + jkr - kr * kr) * exp_jkr;

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

  // Validate computed far-field against analytical solution
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

// Compare the implementation in CurrentDipoleOperator projected to the far-field, using
// SurfacePostOperator, with the analytic expectation in ComputeAnalyticalFarFieldrE.
void runCurrentDipoleTest(double freq_Hz, std::unique_ptr<mfem::Mesh> serial_mesh,
                          const std::vector<int> &attributes)
{
  constexpr double atol = 1e-4;
  constexpr double rtol = 0.05;  // 5% relative error for large domain with absorbing BC
  constexpr double Ids = 1.;


  Units units(1.0, 1.0);  // Pick some arbitrary non-trivial units for testing
  IoData iodata{units};
  iodata.domains.materials.emplace_back().attributes = {1};
  auto &dipole_config = iodata.domains.current_dipole[1];
  double omega = 2.0 * M_PI * units.Nondimensionalize<Units::ValueType::FREQUENCY>(freq_Hz / 1e9);
  double p0 = Ids / (2.0 * M_PI * freq_Hz); // Ids / ω

  dipole_config.moment = {0.0, 0.0, Ids};
  dipole_config.center = {0.0, 0.0, 0.0};
  iodata.boundaries.postpro.farfield.attributes = attributes;
  iodata.boundaries.postpro.farfield.thetaphis.emplace_back();
  iodata.problem.type = ProblemType::DRIVEN;
  iodata.CheckConfiguration();

  auto comm = Mpi::World();

  // Read parallel mesh.
  const int dim = serial_mesh->Dimension();
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  iodata.NondimensionalizeInputs(*par_mesh);
  Mesh palace_mesh(std::move(par_mesh));

  std::vector<std::unique_ptr<Mesh>> mesh_vec;
  mesh_vec.push_back(std::make_unique<Mesh>(palace_mesh));
  SpaceOperator space_op(iodata, mesh_vec);

  auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
  auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  const auto &Curl = space_op.GetCurlMatrix();
  ComplexKspSolver ksp(iodata, space_op.GetNDSpaces(), &space_op.GetH1Spaces());
  ComplexVector RHS(Curl.Width()), E(Curl.Width()), B(Curl.Height());
  E = 0.0;
  B = 0.0;

  auto A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO);
  auto A = space_op.GetSystemMatrix(1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i, K.get(),
                                    C.get(), M.get(), A2.get());
  auto P = space_op.GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, 1i * omega,
                                                             -omega * omega + 0.0i, omega);
  ksp.SetOperators(*A, *P);
  space_op.GetExcitationVector(1, omega, RHS);
  RHS *= (-1i * omega); // RHS = -iω RHS
  ksp.Mult(RHS, E);

  // Compute B = -1/(iω) ∇ x E on the true dofs.
  Curl.Mult(E.Real(), B.Real());
  Curl.Mult(E.Imag(), B.Imag());
  B *= -1.0 / (1i * omega);  // TODO: add error computation for B as well.

  // Recover complete GridFunction from constrained solver solution for error computation.
  GridFunction E_field(space_op.GetNDSpace(), true);
  E_field.Real().SetFromTrueDofs(E.Real());
  E_field.Imag().SetFromTrueDofs(E.Imag());

  // Create analytical solution coefficients.
  auto E_exact_real = [=](const mfem::Vector &x, mfem::Vector &E_out)
  {
    auto E_complex = ComputeDipoleENonDim(x, units, p0, freq_Hz);
    E_out(0) = std::real(E_complex[0]);
    E_out(1) = std::real(E_complex[1]);
    E_out(2) = std::real(E_complex[2]);
  };
  auto E_exact_imag = [=](const mfem::Vector &x, mfem::Vector &E_out)
  {
    auto E_complex = ComputeDipoleENonDim(x, units, p0, freq_Hz);
    E_out(0) = std::imag(E_complex[0]);
    E_out(1) = std::imag(E_complex[1]);
    E_out(2) = std::imag(E_complex[2]);
  };
  mfem::VectorFunctionCoefficient E_exact_real_coef(dim, E_exact_real);
  mfem::VectorFunctionCoefficient E_exact_imag_coef(dim, E_exact_imag);

  // Higher-order integration for accuracy
  int order = 3;
  int order_quad = std::max(2, 2 * order + 1);
  const mfem::IntegrationRule *irs[mfem::Geometry::NumGeom];
  for (int i = 0; i < mfem::Geometry::NumGeom; ++i)
  {
    irs[i] = &(mfem::IntRules.Get(i, order_quad));
  }

  // Compute L2 error
  double error_real = E_field.Real().ComputeL2Error(E_exact_real_coef, irs);
  double error_imag = E_field.Imag().ComputeL2Error(E_exact_imag_coef, irs);
  double total_error = std::sqrt(error_real * error_real + error_imag * error_imag);

  // Compute exact solution norm using zero field trick
  // When field=0, ComputeL2Error(exact_coef) returns ||0 - exact_coef||_L2 =
  // ||exact_coef||_L2
  GridFunction zero_field_for_norm_computation(space_op.GetNDSpace(), true);
  zero_field_for_norm_computation = 0.0;  // Initialize to zero
  double exact_norm_real =
      zero_field_for_norm_computation.Real().ComputeL2Error(E_exact_real_coef, irs);
  double exact_norm_imag =
      zero_field_for_norm_computation.Imag().ComputeL2Error(E_exact_imag_coef, irs);
  double total_exact_norm =
      std::sqrt(exact_norm_real * exact_norm_real + exact_norm_imag * exact_norm_imag);

  // Compute relative error
  double relative_error = total_error / total_exact_norm;

  // ------------------------Debugging--------------------------------
  // Debug analytical field values at key locations to understand the field pattern
  std::cout << "\n=== Analytical Field Debug at Key Points ===";
  std::vector<std::pair<mfem::Vector, std::string>> test_points = {
      {mfem::Vector({0.0, 0.0, 0.0}), "Origin"},
      {mfem::Vector({0.0, 0.0, 0.02}), "+0.02 origin on z-axis"},
      {mfem::Vector({0.0, 0.0, -0.1}), "-0.1m from origin on z-axis"},
      {mfem::Vector({-0.495, 0.0, 0.0}), "-0.495 field on x-axis"},
      {mfem::Vector({0.0, 0.495, 0.0}), "+0.495 on y-axis"},
      {mfem::Vector({0.0, 0.0, 0.495}), "+0.495m from origin on z-axis"},
      {mfem::Vector({0.5, 0.0, 0.0}), "-0.5m from origin on x-axis"},

  };

  for (const auto &[pt, description] : test_points)
  {
    double r_nondim = pt.Norml2();
    double r_dim = units.Dimensionalize<Units::ValueType::LENGTH>(r_nondim);

    // Evaluate nondimensional analytical solution
    auto E_analytical = ComputeDipoleENonDim(pt, units, p0, freq_Hz);

    double analytical_magnitude =
        std::sqrt(std::real(E_analytical[0] * std::conj(E_analytical[0])) +
                  std::real(E_analytical[1] * std::conj(E_analytical[1])) +
                  std::real(E_analytical[2] * std::conj(E_analytical[2])));

    // Try to evaluate Palace FEM solution at this point
    mfem::Vector E_palace_real(dim), E_palace_imag(dim);
    E_palace_real = 0.0;
    E_palace_imag = 0.0;
    bool palace_eval_success = false;

    try
    {
      // Find which element contains this point
      auto &mesh = space_op.GetNDSpace().GetParMesh();
      mfem::Array<int> elem_ids;
      mfem::Array<mfem::IntegrationPoint> ips;

      // Convert Vector to DenseMatrix for FindPoints
      mfem::DenseMatrix point_mat(dim, 1);
      for (int i = 0; i < dim; i++)
      {
        point_mat(i, 0) = pt(i);
      }

      mesh.FindPoints(point_mat, elem_ids, ips);

      if (elem_ids.Size() > 0 && elem_ids[0] >= 0)
      {
        // Evaluate the gridfunction at this point
        E_field.Real().GetVectorValue(elem_ids[0], ips[0], E_palace_real);
        E_field.Imag().GetVectorValue(elem_ids[0], ips[0], E_palace_imag);
        palace_eval_success = true;
      }
    }
    catch (...)
    {
      // Evaluation failed - just continue with zero values
    }

    double palace_magnitude = std::sqrt(
        E_palace_real(0) * E_palace_real(0) + E_palace_real(1) * E_palace_real(1) +
        E_palace_real(2) * E_palace_real(2) + E_palace_imag(0) * E_palace_imag(0) +
        E_palace_imag(1) * E_palace_imag(1) + E_palace_imag(2) * E_palace_imag(2));

    std::cout << "\n"
              << description << " r=" << r_dim << "m:" << " r_nondim=" << r_nondim << "m:";
    std::cout << "\n  PALACE FEM SOLUTION:";
    if (palace_eval_success)
    {
      std::cout << "\n    Palace |E|: " << palace_magnitude;
      std::cout << "\n    Palace Ex: "
                << std::complex<double>(E_palace_real(0), E_palace_imag(0));
      std::cout << "\n    Palace Ey: "
                << std::complex<double>(E_palace_real(1), E_palace_imag(1));
      std::cout << "\n    Palace Ez: "
                << std::complex<double>(E_palace_real(2), E_palace_imag(2));
    }
    else
    {
      std::cout << "\n    Palace evaluation failed at this point";
    }
    std::cout << "\n  ANALYTICAL SOLUTION:";
    std::cout << "\n    Analytical |E|: " << analytical_magnitude;
    std::cout << "\n    Analytical Ex: " << E_analytical[0];
    std::cout << "\n    Analytical Ey: " << E_analytical[1];
    std::cout << "\n    Analytical Ez: " << E_analytical[2];
    if (palace_eval_success && analytical_magnitude > 1e-12)
    {
      std::cout << "\n  COMPARISON:";
      std::cout << "\n    |E| ratio (Palace/Analytical): "
                << palace_magnitude / analytical_magnitude;
    }
    std::cout << "\n  Field decay: |E| ~ 1/r^"
              << (r_dim > 0.5 ? std::log(analytical_magnitude) / std::log(1.0 / r_dim)
                              : 0.0);
  }
  std::cout << "\n============================================\n";

  std::cout << "\n\n--- Absolute L2 Errors ---";
  std::cout << "\nReal part error: || E_h_Re - E_Re ||_L2 = " << error_real;
  std::cout << "\nImag part error: || E_h_Im - E_Im ||_L2 = " << error_imag;
  std::cout << "\nTotal error:     || E_h - E ||_L2 = " << total_error;

  std::cout << "\n\n--- Exact Solution Norms ---";
  std::cout << "\nReal part norm: || E_Re ||_L2 = " << exact_norm_real;
  std::cout << "\nImag part norm: || E_Im ||_L2 = " << exact_norm_imag;
  std::cout << "\nTotal norm:     || E ||_L2 = " << total_exact_norm;

  std::cout << "\n\n--- Relative L2 Errors ---";
  std::cout << "\nReal part: || E_h_Re - E_Re || / ||E_Re|| = "
            << error_real / exact_norm_real;
  std::cout << "\nImag part: || E_h_Im - E_Im || / ||E_Im|| = "
            << error_imag / exact_norm_imag;
  std::cout << "\nTotal:     || E_h - E || / ||E|| = " << relative_error;
  // -----------------------------------------------------------------------

  CHECK_THAT(relative_error, WithinAbs(0.0, rtol));
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

TEST_CASE("Electric Current Dipole implementation", "[currentdipole][strattonchu][Serial]")
{
  double freq_Hz = 35e6;
  std::vector<int> attributes = {1, 2, 3, 4, 5, 6};

  // Make mesh for a cube [0, 1] x [0, 1] x [0, 1]
  int resolution = 20;
  std::unique_ptr<mfem::Mesh> serial_mesh =
      std::make_unique<mfem::Mesh>(mfem::Mesh::MakeCartesian3D(
          resolution, resolution, resolution, mfem::Element::TETRAHEDRON));

  // Offset the cube to center the origin
  serial_mesh->Transform(
      [](const mfem::Vector &x, mfem::Vector &p)
      {
        p = x;
        p(0) -= 0.5;  // Transform [0,1] -> [-0.5,0.5]
        p(1) -= 0.5;
        p(2) -= 0.5;
      });

  runCurrentDipoleTest(freq_Hz, std::move(serial_mesh), attributes);
}

}  // namespace
}  // namespace palace
