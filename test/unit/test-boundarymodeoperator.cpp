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
#include "models/boundarymodeoperator.hpp"
#include "models/farfieldboundaryoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/surfaceconductivityoperator.hpp"
#include "models/surfaceimpedanceoperator.hpp"
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
// BoundaryModeOperator. Returns eigenvalues as complex kn values.
// MakeCartesian2D boundary attributes: bottom=1, right=2, top=3, left=4.
struct ModeResult
{
  std::vector<std::complex<double>> kn;
  int num_converged;
};

ModeResult SolveRectangularModes(double width, double height, double freq_ghz,
                                 double epsilon_r, int order, int num_modes,
                                 const std::function<void(IoData &)> &configure_bcs)
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
  iodata.solver.mode_analysis.freq = freq_ghz;
  iodata.solver.mode_analysis.n = num_modes;
  iodata.solver.mode_analysis.tol = 1.0e-8;
  iodata.solver.linear.tol = 1.0e-8;
  iodata.solver.linear.max_it = 200;

  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(10, 5, mfem::Element::TRIANGLE, false, width, height));
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  iodata.NondimensionalizeInputs(*par_mesh);
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
  double c_min = mat_op.GetLightSpeedMax().Min();
  Mpi::GlobalMin(1, &c_min, nd_fespace.GetComm());
  double kn_target = omega / c_min * std::sqrt(1.1);

  BoundaryModeOperatorConfig config;
  config.attr_to_material = &mat_op.GetAttributeToMaterial();
  config.inv_permeability = &mat_op.GetInvPermeability();
  config.curlcurl_inv_permeability = &mat_op.GetCurlCurlInvPermeability();
  config.permittivity_real = &mat_op.GetPermittivityReal();
  config.permittivity_scalar = &mat_op.GetPermittivityScalar();
  config.normal = nullptr;
  config.permittivity_imag =
      mat_op.HasLossTangent() ? &mat_op.GetPermittivityImag() : nullptr;
  config.permittivity_imag_scalar =
      mat_op.HasLossTangent() ? &mat_op.GetPermittivityImagScalar() : nullptr;
  config.has_loss_tangent = mat_op.HasLossTangent();
  config.conductivity = mat_op.HasConductivity() ? &mat_op.GetConductivity() : nullptr;
  config.has_conductivity = mat_op.HasConductivity();
  config.inv_london_depth = mat_op.HasLondonDepth() ? &mat_op.GetInvLondonDepth() : nullptr;
  config.has_london_depth = mat_op.HasLondonDepth();
  config.mat_op = &mat_op;
  config.surf_z_op = &surf_z_op;
  config.farfield_op = &farfield_op;
  config.surf_sigma_op = &surf_sigma_op;
  config.num_modes = num_modes;
  config.num_vec = -1;
  config.eig_tol = 1.0e-8;
  config.which_eig = EigenvalueSolver::WhichType::LARGEST_REAL;
  config.linear = &iodata.solver.linear;
  config.eigen_backend = iodata.solver.mode_analysis.type;
  config.verbose = 0;

  BoundaryModeOperator mode_solver(config, nd_fespace, h1_fespace, dbc_tdof_list,
                                   nd_fespace.GetComm());

  double sigma = -kn_target * kn_target;
  auto result = mode_solver.Solve(omega, sigma);

  ModeResult out;
  out.num_converged = result.num_converged;
  for (int i = 0; i < result.num_converged; i++)
  {
    auto lambda = mode_solver.GetEigenvalue(i);
    out.kn.push_back(std::sqrt(-sigma - 1.0 / lambda));
  }
  return out;
}

}  // namespace

TEST_CASE("BoundaryModeOperator PEC", "[boundarymodeoperator][Serial]")
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
  // Allow 2% tolerance for the coarse 10×5 mesh at order 2.
  CHECK_THAT(kn_real, WithinRel(0.02071, 0.02));
}

TEST_CASE("BoundaryModeOperator Impedance shifts kn", "[boundarymodeoperator][Serial]")
{
  auto pec_result = SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 3, [](IoData &) {});

  auto imp_result = SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 3,
                                          [](IoData &iodata)
                                          {
                                            iodata.boundaries.pec.attributes = {1, 3, 4};
                                            auto &imp =
                                                iodata.boundaries.impedance.emplace_back();
                                            imp.attributes = {2};
                                            imp.Ls = 1.0e-10;
                                          });

  REQUIRE(pec_result.num_converged >= 1);
  REQUIRE(imp_result.num_converged >= 1);

  CHECK(imp_result.kn[0].real() > pec_result.kn[0].real());
}

TEST_CASE("BoundaryModeOperator Conductivity adds loss", "[boundarymodeoperator][Serial]")
{
  auto pec_result = SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 3, [](IoData &) {});

  auto cond_result =
      SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 3,
                            [](IoData &iodata)
                            {
                              iodata.boundaries.pec.attributes = {1, 3, 4};
                              auto &cond = iodata.boundaries.conductivity.emplace_back();
                              cond.attributes = {2};
                              cond.sigma = 5.0e7;
                              cond.h = 0.001;
                            });

  REQUIRE(pec_result.num_converged >= 1);
  REQUIRE(cond_result.num_converged >= 1);

  CHECK(std::abs(cond_result.kn[0].imag()) > std::abs(pec_result.kn[0].imag()));
}

}  // namespace palace
