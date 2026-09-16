// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <complex>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "models/domainpostoperator.hpp"
#include "models/lumpedportoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/spaceoperator.hpp"
#include "models/surfacepostoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/constants.hpp"
#include "utils/filesystem.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

// Physical normalization of quadratic (power, energy) outputs for the peak-phasor time-
// harmonic convention. The expected values below are written directly from the physical
// definitions and do not use any Palace helper:
//
//   - For a complex peak phasor F, the time average of the physical product over one
//     period is <Re{F e^{iωt}} Re{G e^{iωt}}> = 1/2 Re{F G⋆}, so the time-averaged energy
//     densities and Poynting vector are U_e = 1/4 ε |E|², U_m = 1/4 μ |H|², and
//     S = 1/2 Re{E x H⋆}.
//   - For real-valued instantaneous fields, U_e = 1/2 ε E², U_m = 1/2 μ H², S = E x H.
//   - An excited lumped port carries unit time-averaged incident power in the frequency
//     domain: for R = 50 Ω, V_inc = 10 V and I_inc = 0.2 A (peak), and the incident field
//     itself projects to S = 1 and carries 1 W through the port.
//
// The reference fields are E = (1 + i, 0, 0) V/m and H = (0, 2 + 3i, 0) A/m in vacuum,
// so that |E|² = 2, |H|² = 13, and E x H⋆ = (1 + i)(2 - 3i) ẑ = (5 - i) ẑ: the
// time-averaged power density is (2.5 - 0.5i) W/m² along +z, U_e = ε₀/2 J/m³, and U_m = 13
// μ₀/4 J/m³.

namespace palace
{
using namespace Catch::Matchers;
using json = nlohmann::json;
using VT = Units::ValueType;

namespace
{

struct UniformFields
{
  // Cartesian components of the (possibly complex) uniform E [V/m] and H [A/m] fields.
  std::complex<double> Ex, Hy;
  bool complex;
};

// Project uniform E (ND) and B = μ₀ H (RT) fields, given in SI units, onto the grid
// functions in nondimensional units.
void ProjectUniformFields(const Units &units, const UniformFields &f, GridFunction &E,
                          GridFunction &B)
{
  auto Project = [&](mfem::ParGridFunction &gf, double vx, double vy, double vz)
  {
    mfem::Vector v(3);
    v(0) = vx;
    v(1) = vy;
    v(2) = vz;
    mfem::VectorConstantCoefficient coeff(v);
    gf.ProjectCoefficient(coeff);
  };
  const double mu0 = electromagnetics::mu0_;
  Project(E.Real(), units.Nondimensionalize<VT::FIELD_E>(f.Ex.real()), 0.0, 0.0);
  Project(B.Real(), 0.0, units.Nondimensionalize<VT::FIELD_B>(mu0 * f.Hy.real()), 0.0);
  if (f.complex)
  {
    Project(E.Imag(), units.Nondimensionalize<VT::FIELD_E>(f.Ex.imag()), 0.0, 0.0);
    Project(B.Imag(), 0.0, units.Nondimensionalize<VT::FIELD_B>(mu0 * f.Hy.imag()), 0.0);
  }
}

// Time-averaged (complex phasor) or instantaneous (real) physical products, written from
// the definitions.
double ExpectedElectricEnergyDensity(const UniformFields &f)
{
  return (f.complex ? 0.25 : 0.5) * electromagnetics::epsilon0_ * std::norm(f.Ex);
}

double ExpectedMagneticEnergyDensity(const UniformFields &f)
{
  return (f.complex ? 0.25 : 0.5) * electromagnetics::mu0_ * std::norm(f.Hy);
}

std::complex<double> ExpectedPoyntingZ(const UniformFields &f)
{
  // E x H⋆ = Ex Hy⋆ ẑ.
  return (f.complex ? 0.5 : 1.0) * f.Ex * std::conj(f.Hy);
}

}  // namespace

TEST_CASE("Time-averaged domain energies, surface power flux, and interface energy",
          "[timeaverage][Serial][Parallel]")
{
  constexpr int order = 1;
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  // Cases: complex phasor with nonzero real and imaginary parts, a global phase rotation of
  // the same phasor (all outputs invariant), a purely real complex phasor (still time-
  // averaged), and the real instantaneous control.
  const std::complex<double> rot = std::polar(1.0, 0.7);
  const std::vector<UniformFields> cases = {
      {std::complex<double>(1.0, 1.0), std::complex<double>(2.0, 3.0), true},
      {rot * std::complex<double>(1.0, 1.0), rot * std::complex<double>(2.0, 3.0), true},
      {std::complex<double>(1.0, 0.0), std::complex<double>(2.0, 0.0), true},
      {std::complex<double>(1.0, 0.0), std::complex<double>(2.0, 0.0), false}};
  const auto &fields = cases[GENERATE(0, 1, 2, 3)];
  CAPTURE(fields.Ex, fields.Hy, fields.complex);

  Units units(0.496, 1.453);
  MPI_Comm comm = Mpi::World();

  // Vacuum cube [0, sx] x [0, sy] x [0, sz], in nondimensional mesh units.
  constexpr double sx = 1.1, sy = 2.5, sz = 3.8;
  constexpr int resolution = 4;
  auto serial_mesh = std::make_unique<mfem::Mesh>(mfem::Mesh::MakeCartesian3D(
      resolution, resolution, resolution, mfem::Element::HEXAHEDRON, sx, sy, sz));
  const int dim = serial_mesh->Dimension();
  Mesh mesh(std::make_unique<mfem::ParMesh>(comm, *serial_mesh));
  serial_mesh.reset();

  config::MaterialData material;
  material.attributes = {1};
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({material}, periodic, ProblemType::DRIVEN, mesh);

  mfem::H1_FECollection h1_fec(order, dim);
  mfem::ND_FECollection nd_fec(order, dim);
  mfem::RT_FECollection rt_fec(order - 1, dim);
  FiniteElementSpace h1_fespace(mesh, &h1_fec), nd_fespace(mesh, &nd_fec),
      rt_fespace(mesh, &rt_fec);

  GridFunction E(nd_fespace, fields.complex), B(rt_fespace, fields.complex);
  ProjectUniformFields(units, fields, E, B);

  const double sx_m = units.Dimensionalize<VT::LENGTH>(sx);
  const double sy_m = units.Dimensionalize<VT::LENGTH>(sy);
  const double sz_m = units.Dimensionalize<VT::LENGTH>(sz);
  const double volume_m3 = sx_m * sy_m * sz_m;
  const double top_area_m2 = sx_m * sy_m;

  SECTION("Domain energies")
  {
    config::DomainPostData postpro;
    config::DomainEnergyData energy;
    energy.attributes = {1};
    postpro.energy.emplace(1, energy);
    DomainPostOperator dom_post_op(postpro, mat_op, nd_fespace, rt_fespace);

    const double U_e =
        units.Dimensionalize<VT::ENERGY>(dom_post_op.GetElectricFieldEnergy(E));
    const double U_m =
        units.Dimensionalize<VT::ENERGY>(dom_post_op.GetMagneticFieldEnergy(B));
    CHECK_THAT(U_e, WithinRel(ExpectedElectricEnergyDensity(fields) * volume_m3, 1.0e-10));
    CHECK_THAT(U_m, WithinRel(ExpectedMagneticEnergyDensity(fields) * volume_m3, 1.0e-10));

    // Subdomain reducers follow the same convention.
    const double U_e_1 =
        units.Dimensionalize<VT::ENERGY>(dom_post_op.GetDomainElectricFieldEnergy(1, E));
    const double U_m_1 =
        units.Dimensionalize<VT::ENERGY>(dom_post_op.GetDomainMagneticFieldEnergy(1, B));
    CHECK_THAT(U_e_1, WithinRel(U_e, 1.0e-12));
    CHECK_THAT(U_m_1, WithinRel(U_m, 1.0e-12));
  }

  SECTION("Surface power flux and interface energy")
  {
    // Power flux through the top face z = sz (boundary attribute 6 of MakeCartesian3D),
    // with the normal oriented outward from the cube center, so that the +z power density
    // gives a positive flux. The interface dielectric layer (thickness t, relative
    // permittivity ε_r, both nondimensional here) on the same face has energy
    // 1/2 t ε_r |E|² A (real) or 1/4 t ε_r |E|² A (complex peak phasor), in the internal
    // units where ε₀ = 1.
    constexpr int top_attr = 6;
    constexpr double t = 0.01, eps_r = 3.0;
    config::BoundaryPostData postpro;
    config::SurfaceFluxData flux;
    flux.type = SurfaceFlux::POWER;
    flux.two_sided = false;
    flux.center = {0.5 * sx, 0.5 * sy, 0.5 * sz};
    flux.no_center = false;
    flux.attributes = {top_attr};
    postpro.flux.emplace(1, flux);
    config::InterfaceDielectricData dielectric;
    dielectric.type = InterfaceDielectric::DEFAULT;
    dielectric.t = t;
    dielectric.epsilon_r = eps_r;
    dielectric.tandelta = 0.0;
    dielectric.attributes = {top_attr};
    postpro.dielectric.emplace(1, dielectric);
    SurfacePostOperator surf_post_op(postpro, ProblemType::DRIVEN, mat_op, h1_fespace,
                                     nd_fespace);

    const auto flux_nondim = surf_post_op.GetSurfaceFlux(1, &E, &B);
    const double P = units.Dimensionalize<VT::POWER>(flux_nondim.real());
    CHECK_THAT(P, WithinRel(ExpectedPoyntingZ(fields).real() * top_area_m2, 1.0e-10));
    // Only the stationary real part of the power flux is reported.
    CHECK_THAT(flux_nondim.imag(), WithinAbs(0.0, 1.0e-14));

    const double E2_nondim = std::norm(std::complex<double>(
        units.Nondimensionalize<VT::FIELD_E>(fields.Ex.real()),
        fields.complex ? units.Nondimensionalize<VT::FIELD_E>(fields.Ex.imag()) : 0.0));
    const double expected_interface =
        (fields.complex ? 0.25 : 0.5) * t * eps_r * E2_nondim * sx * sy;
    CHECK_THAT(surf_post_op.GetInterfaceElectricFieldEnergy(1, E),
               WithinRel(expected_interface, 1.0e-10));
  }
}

TEST_CASE("Lumped port unit time-averaged incident power",
          "[timeaverage][Serial][Parallel]")
{
  // A 50 Ω rectangular lumped port on the square (0,0,0)-(1,1,0) of the (3,2,1) box,
  // direction +Y, with the domain at z > 0 (see
  // LumpedPort_BasicTests_1ElementPort_Cube321). Build the incident wave of the port
  // explicitly, E = V_inc / l ŷ and H = -E / R_s x̂ (propagating into the domain, +z), from
  // the excitation amplitudes reported by the port, and check with the production port
  // functionals that the incident wave carries 1 W of time-averaged power, projects to S =
  // 1, and reproduces V_inc.
  MPI_Comm comm = Mpi::World();
  auto solver_order = GENERATE(1, 2);
  // Purely real phasor, and the same wave with a global phase.
  const std::complex<double> phase =
      GENERATE(std::complex<double>(1.0, 0.0), std::polar(1.0, -1.1));
  CAPTURE(solver_order, phase);

  const double L0 = 1.0e-6;
  const double Lc = 7.0;
  IoData iodata{Units(L0, Lc)};
  iodata.model.mesh =
      fs::path(PALACE_TEST_DATA_DIR) / "lumpedport_mesh/cube_mesh_3_2_1_hex.msh";
  iodata.model.L0 = L0;
  iodata.model.Lc = Lc;
  iodata.model.crack_bdr_elements = false;
  iodata.solver.order = solver_order;
  json domains_json = {
      {"Materials",
       json::array({json::object({{"Attributes", json::array({1, 2, 3, 4, 5, 6})},
                                  {"Permeability", 1.0},
                                  {"Permittivity", 1.0},
                                  {"LossTan", 0.0}})})}};
  iodata.domains = config::DomainData(domains_json);
  json boundaries_json = {
      {"LumpedPort", json::array({json::object({{"Index", 1},
                                                {"R", 50.0},
                                                {"Excitation", 1},
                                                {"Attributes", json::array({1})},
                                                {"Direction", "+Y"}})})}};
  iodata.boundaries = config::BoundaryData(boundaries_json);
  iodata.CheckConfiguration();

  std::vector<std::unique_ptr<Mesh>> mesh;
  {
    std::vector<std::unique_ptr<mfem::ParMesh>> mfem_mesh;
    auto smesh = mesh::Load(iodata, comm);
    iodata.NondimensionalizeInputs(smesh);
    mfem_mesh.push_back(mesh::Partition(iodata, std::move(smesh), comm));
    mesh::RefineMesh(iodata, mfem_mesh);
    for (auto &m : mfem_mesh)
    {
      mesh.push_back(std::make_unique<Mesh>(std::move(m)));
    }
  }
  SpaceOperator space_op(iodata, mesh);
  const auto &port = space_op.GetLumpedPortOp().GetPort(1);
  REQUIRE(port.HasExcitation());
  REQUIRE(port.elems.size() == 1);
  const auto &elem = *port.elems.front();

  // Excitation amplitudes: 10 V and 0.2 A peak for the 1 W time-averaged wave.
  const double V_inc = port.GetExcitationVoltage(true);
  const double I_inc = port.GetExcitationCurrent(true);
  CHECK_THAT(iodata.units.Dimensionalize<VT::VOLTAGE>(V_inc), WithinRel(10.0, 1.0e-12));
  CHECK_THAT(iodata.units.Dimensionalize<VT::CURRENT>(I_inc), WithinRel(0.2, 1.0e-12));
  CHECK_THAT(iodata.units.Dimensionalize<VT::POWER>(0.5 * V_inc * I_inc),
             WithinRel(1.0, 1.0e-12));

  // Incident wave fields from the port geometry: E_inc = V_inc / l along the port direction
  // and H_inc = E_inc / R_s with R_s = R w / l the sheet resistance, oriented so that
  // E x H points into the domain (+z).
  const double E_inc = V_inc / elem.GetGeometryLength();
  const double Rs = port.GetExcitationRefResistance() * port.GetToSquare(elem);
  const double H_inc = E_inc / Rs;
  GridFunction E(space_op.GetNDSpace(), true), B(space_op.GetRTSpace(), true);
  auto Project = [&](mfem::ParGridFunction &gf, double vx, double vy, double vz)
  {
    mfem::Vector v(3);
    v(0) = vx;
    v(1) = vy;
    v(2) = vz;
    mfem::VectorConstantCoefficient coeff(v);
    gf.ProjectCoefficient(coeff);
  };
  const std::complex<double> Ey = phase * E_inc, Hx = -phase * H_inc;
  Project(E.Real(), 0.0, Ey.real(), 0.0);
  Project(E.Imag(), 0.0, Ey.imag(), 0.0);
  Project(B.Real(), Hx.real(), 0.0, 0.0);  // μ_r = 1 in nondimensional units.
  Project(B.Imag(), Hx.imag(), 0.0, 0.0);

  // Port voltage and S-parameter of the incident wave itself.
  const auto V = port.GetVoltage(E);
  CHECK_THAT(V.real(), WithinRel((phase * V_inc).real(), 1.0e-10));
  CHECK_THAT(V.imag(), WithinAbs((phase * V_inc).imag(), 1.0e-10 * V_inc));
  const auto S = port.GetSParameter(E);
  CHECK_THAT(S.real(), WithinRel(phase.real(), 1.0e-10));
  CHECK_THAT(S.imag(), WithinAbs(phase.imag(), 1.0e-10));

  // Time-averaged power through the port, through both the per-port and batched paths,
  // and the independent legacy coefficient path.
  const auto P = port.GetPower(E, B);
  CHECK_THAT(iodata.units.Dimensionalize<VT::POWER>(P.real()), WithinRel(1.0, 1.0e-10));
  CHECK_THAT(P.imag(), WithinAbs(0.0, 1.0e-10));
  const auto P_legacy = port.GetPowerLegacy(E, B);
  CHECK_THAT(P_legacy.real(), WithinRel(P.real(), 1.0e-10));
  CHECK_THAT(P_legacy.imag(), WithinAbs(0.0, 1.0e-10));
  const auto powers = space_op.GetLumpedPortOp().GetPowers(E, B);
  REQUIRE(powers.count(1) == 1);
  CHECK_THAT(powers.at(1).real(), WithinRel(P.real(), 1.0e-10));
  CHECK_THAT(powers.at(1).imag(), WithinAbs(0.0, 1.0e-10));

  // The assembled frequency domain source term is the dual of 2 H_inc for this incident
  // wave: its action on the incident E field is 2 ∫ E_inc ⋅ H_inc dS = 2 V_inc I_inc = 4 W,
  // twice the full power overlap of the unit time-averaged power mode.
  ComplexVector RHS1;
  REQUIRE(space_op.GetExcitationVector1(1, RHS1));
  Vector e_inc(space_op.GetNDSpace().GetTrueVSize());
  E.Real().ParallelProject(e_inc);
  double rhs_dot_e = linalg::Dot(comm, RHS1.Real(), e_inc);
  CHECK_THAT(iodata.units.Dimensionalize<VT::POWER>(rhs_dot_e),
             WithinRel(4.0 * phase.real(), 1.0e-10));
}

}  // namespace palace
