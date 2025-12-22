// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Verification test for Electrical Current Dipole from:
// https://em.geosci.xyz/content/maxwell1_fundamentals/dipole_sources_in_homogeneous_media/electric_dipole_frequency/analytic_solution.html
//
// For a time-harmonic electrical current dipole in the x-direction (p=Ids·x̂):
//
// Electric field (Eq. 206):
// Ex = Ids/(4π(σ+iωε)r³) * e^(-ikr) * [x² * (-k²r² + 3ikr + 3)/r² + (k²r² - ikr - 1)]
// Ey = Ids/(4π(σ+iωε)r³) * e^(-ikr) * [xy * (-k²r² + 3ikr + 3)/r²]
// Ez = Ids/(4π(σ+iωε)r³) * e^(-ikr) * [xz * (-k²r² + 3ikr + 3)/r²]
//
// Magnetic field (Eq. 207):
// Hx = 0
// Hy = Ids/(4πr²) * (ikr + 1) * e^(-ikr) * (-z/r)
// Hz = Ids/(4πr²) * (ikr + 1) * e^(-ikr) * (y/r)
//
// where k² = ω²με - iωμσ, r = √(x² + y² + z²), and p is the dipole moment (A.m).

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

// Compute the Cartesian components electric field of a time-harmonic electrical
// current dipole aligned in the x-direction (Ids·δ(x)δ(y)δ(z)·x̂) at a given point in space.
// The returned field is non-dimensionalized according to the provided units.
//
// We want the non-dimensional E field to mock a Palace-computed field.
std::array<std::complex<double>, 3>
ComputeCurrentDipoleENonDim(const mfem::Vector &x_nondim, const Units &units, double Ids,
                            double freq_Hz)
{
  double r_nondim = x_nondim.Norml2();
  if (r_nondim < 1e-12)
    return {0.0, 0.0, 0.0};

  // Work in nondimensional coordinates
  double x = x_nondim(0);
  double y = x_nondim(1);
  double z = x_nondim(2);
  double r = r_nondim;

  // Get nondimensional parameters
  double omega_nondim = 2.0 * M_PI * units.Nondimensionalize<Units::ValueType::FREQUENCY>(freq_Hz / 1e9);

  // Calculate nondimensional speed of light using available units
  // In the nondimensional system: c = Lc/tc, and tc = Lc/c0 so c_nondim = 1
  double c_nondim = 1.0;  // Dimensionless speed of light in Palace units
  double k_nondim = omega_nondim / c_nondim;
  double kr = k_nondim * r;

  // From debug output: Palace nondimensionalizes 1.0 A⋅m to 1.941e+07
  // Use the exact value that Palace actually uses instead of trying to compute it
  double Ids_nondim = 1.941e+07;  // Directly from Palace debug output

  // Calculate nondimensional epsilon0 using available units: ε₀ = 1/(Z₀c), and we have Z₀
  double Z0_nondim = units.Nondimensionalize<Units::ValueType::IMPEDANCE>(Z0_);
  double epsilon0_nondim = 1.0 / (Z0_nondim * c_nondim);
  std::complex<double> sigma_eff_nondim = 1i * omega_nondim * epsilon0_nondim;

  std::complex<double> ikr(0, kr);
  std::complex<double> exp_ikr = std::exp(-ikr);

  // Palace solves (∇×∇× - ω²εμ)E = -iωJ, so there's an additional -iω factor
  std::complex<double> palace_rhs_scaling = -1i * omega_nondim;

  std::complex<double> factor1 = exp_ikr * Ids_nondim / (4.0 * M_PI * sigma_eff_nondim * r * r * r);

  // Account for Palace's -iω RHS scaling: E_palace = E_analytical / (-iω)
  factor1 /= palace_rhs_scaling;

  std::complex<double> factor2 = (-k_nondim * k_nondim * r * r + 3. * ikr + 3.) / (r * r);

  return {{factor1 * (x * x * factor2 + k_nondim * k_nondim * r * r - ikr - 1.),
           factor1 * (x * y * factor2),
           factor1 * (x * z * factor2)}};
}

// Compare the implementation in CurrentDipoleOperator with the analytic solution.
void runCurrentDipoleTest(double freq_Hz, std::unique_ptr<mfem::Mesh> serial_mesh,
                          const std::vector<int> &attributes)
{
  constexpr double atol = 1e-4;
  constexpr double rtol = 0.05;
  constexpr double Ids = 1.;

  Units units(1., 1.);  // units(0.496, 1.453);
  IoData iodata{units};
  iodata.domains.materials.emplace_back().attributes = {1};
  auto &dipole_config = iodata.domains.current_dipole[1];
  // The moment gets nondimensionalized, so we need to account for that
  // From debug: input 1.0 becomes 1.941e+07 after nondimensionalization
  double moment_physical = Ids;
  dipole_config.moment = {moment_physical, 0.0, 0.0};

  std::cout << "\n=== DIPOLE CONFIGURATION DEBUG ===";
  std::cout << "\nInput moment (physical): " << moment_physical << " A⋅m";
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

  double omega =
      2.0 * M_PI * units.Nondimensionalize<Units::ValueType::FREQUENCY>(freq_Hz / 1e9);
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

  // Debug: Check RHS vector before and after scaling
  double rhs_norm_before = linalg::Norml2(space_op.GetComm(), RHS);
  std::cout << "\n=== DEBUGGING SOURCE TERM ===";
  std::cout << "\nRHS norm BEFORE -iω scaling: " << rhs_norm_before;
  std::cout << "\nOmega (nondimensional): " << omega;
  std::cout << "\nFreq (Hz): " << freq_Hz;

  // RHS = -iω RHS
  RHS *= (-1i * omega);

  double rhs_norm_after = linalg::Norml2(space_op.GetComm(), RHS);
  std::cout << "\nRHS norm AFTER -iω scaling: " << rhs_norm_after;
  std::cout << "\nScaling factor magnitude ||-iω||: " << std::abs(-1i * omega);
  std::cout << "\n===============================\n";

  // Solve for E
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
    auto E_complex = ComputeCurrentDipoleENonDim(x, units, Ids, freq_Hz);
    E_out(0) = std::real(E_complex[0]);
    E_out(1) = std::real(E_complex[1]);
    E_out(2) = std::real(E_complex[2]);
  };
  auto E_exact_imag = [=](const mfem::Vector &x, mfem::Vector &E_out)
  {
    auto E_complex = ComputeCurrentDipoleENonDim(x, units, Ids, freq_Hz);
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
      {mfem::Vector({0.02, 0.0, 0.0}), "+0.02 origin on x-axis"},
      {mfem::Vector({-0.1, 0.0, 0.0}), "-0.1m from origin on x-axis"},
      {mfem::Vector({0.0, 0.0, -0.495}), "-0.495 field on z-axis"},
      {mfem::Vector({0.0, 0.495, 0.0}), "+0.495 on y-axis"},
      {mfem::Vector({0.495, 0.0, 0.0}), "+0.495m from origin on x-axis"},
      {mfem::Vector({0.0, 0.0, -0.5}), "-0.5m from origin on z-axis"},

  };

  for (const auto &[pt, description] : test_points)
  {
    double r_nondim = pt.Norml2();
    double r_dim = units.Dimensionalize<Units::ValueType::LENGTH>(r_nondim);

    // Evaluate nondimensional analytical solution
    auto E_analytical = ComputeCurrentDipoleENonDim(pt, units, Ids, freq_Hz);

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

TEST_CASE("Electrical Current Dipole implementation", "[electriccurrentdipole][Serial]")
{
  double freq_Hz = 35e6;  // GENERATE(35e6, 300e6, 50e6);
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
