// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <cmath>
#include <complex>
#include <memory>
#include <numbers>
#include <tuple>
#include <vector>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "fem/bilinearform.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/complexconductivitymodel.hpp"
#include "models/dispersiveconductivityoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/constants.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

namespace palace
{
using namespace Catch;
using namespace Catch::Matchers;
using namespace std::complex_literals;

namespace
{

// Complete elliptic integrals of the first and second kind K(k), E(k) of modulus k, by the
// arithmetic-geometric mean (libc++ lacks the C++17 special math functions).
std::pair<double, double> EllipticKE(double k)
{
  double a = 1.0, b = std::sqrt(1.0 - k * k), c = k;
  double sum = 0.5 * c * c, pow2 = 0.5;
  for (int n = 0; n < 50 && std::abs(c) > 1.0e-17; n++)
  {
    const double an = 0.5 * (a + b);
    c = 0.5 * (a - b);
    b = std::sqrt(a * b);
    a = an;
    pow2 *= 2.0;
    sum += pow2 * c * c;
  }
  const double K = std::numbers::pi / (2.0 * a);
  return {K, K * (1.0 - sum)};
}

// Nondimensional Mattis-Bardeen model data, with gap and temperature as angular
// frequencies.
config::ComplexConductivityData MattisBardeenData(double sigma_n, double omega_gap,
                                                  double omega_T, double omega_ref = 0.0)
{
  config::ComplexConductivityData data;
  data.type = ComplexConductivityType::MATTIS_BARDEEN;
  data.sigma_n = sigma_n;
  data.gap = omega_gap;
  data.temperature = omega_T;
  data.reference_freq = omega_ref;
  return data;
}

// Sets (and restores on destruction) the global quadrature settings which
// IoData::CheckConfiguration normally initializes.
struct IntegrationSettingsGuard
{
  int pa_order_threshold = BilinearForm::pa_order_threshold;
  int p_trial = fem::DefaultIntegrationOrder::p_trial;
  bool q_order_jac = fem::DefaultIntegrationOrder::q_order_jac;
  int q_order_extra_pk = fem::DefaultIntegrationOrder::q_order_extra_pk;
  int q_order_extra_qk = fem::DefaultIntegrationOrder::q_order_extra_qk;

  IntegrationSettingsGuard(const config::SolverData &solver)
  {
    BilinearForm::pa_order_threshold = solver.pa_order_threshold;
    fem::DefaultIntegrationOrder::p_trial = solver.order;
    fem::DefaultIntegrationOrder::q_order_jac = solver.q_order_jac;
    fem::DefaultIntegrationOrder::q_order_extra_pk = solver.q_order_extra;
    fem::DefaultIntegrationOrder::q_order_extra_qk = solver.q_order_extra;
  }
  ~IntegrationSettingsGuard()
  {
    BilinearForm::pa_order_threshold = pa_order_threshold;
    fem::DefaultIntegrationOrder::p_trial = p_trial;
    fem::DefaultIntegrationOrder::q_order_jac = q_order_jac;
    fem::DefaultIntegrationOrder::q_order_extra_pk = q_order_extra_pk;
    fem::DefaultIntegrationOrder::q_order_extra_qk = q_order_extra_qk;
  }
};

// Two-cell mesh with domain attribute 1 for x < 1/2 and 2 otherwise.
std::unique_ptr<mfem::Mesh> TwoAttributeMesh(MPI_Comm comm)
{
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(2, 1, 1, mfem::Element::TETRAHEDRON));
  while (serial_mesh->GetNE() < 4 * Mpi::Size(comm))
  {
    serial_mesh->UniformRefinement();
  }
  for (int i = 0; i < serial_mesh->GetNE(); i++)
  {
    mfem::Vector center;
    serial_mesh->GetElementCenter(i, center);
    serial_mesh->SetAttribute(i, center(0) < 0.5 ? 1 : 2);
  }
  serial_mesh->SetAttributes();
  return serial_mesh;
}

}  // namespace

TEST_CASE("Mattis-Bardeen conductivity", "[complexconductivity][Serial]")
{
  SECTION("Zero temperature closed forms")
  {
    // Mattis and Bardeen (1958), Eqs. (3.11)-(3.13), with x = ħω/Δ.
    for (double x : {0.2, 1.0, 1.8, 2.4, 4.0, 10.0})
    {
      const auto sigma = superconductor::MattisBardeenRatio(x, 0.0);
      const double k = std::abs(2.0 - x) / (2.0 + x), kp = std::sqrt(1.0 - k * k);
      const double sigma_1 = (x > 2.0) ? (1.0 + 2.0 / x) * EllipticKE(k).second -
                                             (4.0 / x) * EllipticKE(k).first
                                       : 0.0;
      const auto [Kp, Ep] = EllipticKE(kp);
      const double sigma_2 = 0.5 * (2.0 / x + 1.0) * Ep + 0.5 * (2.0 / x - 1.0) * Kp;
      CHECK_THAT(sigma.real(), WithinAbs(sigma_1, 1.0e-11));
      CHECK_THAT(-sigma.imag(), WithinRel(sigma_2, 1.0e-11));
    }
  }

  SECTION("Finite temperature reference values")
  {
    // Independent evaluation (SciPy adaptive quadrature, cross-checked against mpmath) at
    // kB T/Δ = 1/4.
    const std::vector<std::array<double, 3>> ref = {
        {0.01, 0.5926330884266798, 303.07535362217567},
        {0.3, 0.13492767097301378, 10.191942805497058},
        {1.8, 0.019043255413684488, 1.2974026759733266},
        {2.5, 0.3054135914933136, 0.5584558160090177},
        {6.0, 0.8340239838641357, 0.0885321361837766}};
    for (const auto &[x, sigma_1, sigma_2] : ref)
    {
      const auto sigma = superconductor::MattisBardeenRatio(x, 0.25);
      CHECK_THAT(sigma.real(), WithinRel(sigma_1, 1.0e-9));
      CHECK_THAT(-sigma.imag(), WithinRel(sigma_2, 1.0e-9));
    }
  }

  SECTION("London and normal-state limits")
  {
    // ħω σ2/σn → π Δ tanh(Δ/2kB T) as ω → 0, and σ → σn for ħω >> Δ.
    const double theta = 0.3;
    const double x = 1.0e-5;
    const auto sigma_low = superconductor::MattisBardeenRatio(x, theta);
    CHECK_THAT(-x * sigma_low.imag(),
               WithinRel(std::numbers::pi * std::tanh(0.5 / theta), 1.0e-5));
    const auto sigma_high = superconductor::MattisBardeenRatio(2.0e3, theta);
    CHECK_THAT(sigma_high.real(), WithinAbs(1.0, 1.0e-4));
    CHECK(-sigma_high.imag() > 0.0);
    CHECK(-sigma_high.imag() < 1.0e-3);
  }
}

TEST_CASE("BCS gap", "[complexconductivity][Serial]")
{
  CHECK(superconductor::BCSGapRatio(0.0) == 1.0);
  CHECK(superconductor::BCSGapRatio(1.0) == 0.0);
  CHECK(superconductor::BCSGapRatio(1.5) == 0.0);
  // Reference roots of the weak-coupling gap equation (independent SciPy evaluation).
  CHECK_THAT(superconductor::BCSGapRatio(0.3), WithinRel(0.9971237624901917, 1.0e-10));
  CHECK_THAT(superconductor::BCSGapRatio(0.5), WithinRel(0.9568847158748907, 1.0e-10));
  CHECK_THAT(superconductor::BCSGapRatio(0.9), WithinRel(0.5263419345511171, 1.0e-10));
  CHECK_THAT(superconductor::BCSGapRatio(0.99), WithinRel(0.17295385407970149, 1.0e-9));
  // Ginzburg-Landau limit Δ(T)/Δ(0) ≈ 1.74 √(1 - t) near Tc.
  CHECK_THAT(superconductor::BCSGapRatio(0.9999), WithinRel(1.74 * 0.01, 1.0e-2));
}

TEST_CASE("ComplexConductivityModel", "[complexconductivity][Serial]")
{
  SECTION("Mattis-Bardeen with a given gap")
  {
    const ComplexConductivityModel model(MattisBardeenData(5.0, 2.0, 0.5));
    for (double omega : {0.1, 1.0, 3.0, 6.0})
    {
      const auto ref = 5.0 * superconductor::MattisBardeenRatio(omega / 2.0, 0.25);
      CHECK(model(omega) == ref);
    }
  }

  SECTION("Mattis-Bardeen with the BCS gap")
  {
    auto data = MattisBardeenData(1.0, 0.0, 0.5);
    data.critical_temperature = 1.0;
    const ComplexConductivityModel model(data);
    const double omega_gap =
        superconductor::BCS_GAP_RATIO * superconductor::BCSGapRatio(0.5);
    CHECK_THAT(model.GetGapFrequency(), WithinRel(omega_gap, 1.0e-14));
    const auto ref = superconductor::MattisBardeenRatio(0.7 / omega_gap, 0.5 / omega_gap);
    CHECK_THAT(model(0.7).real(), WithinRel(ref.real(), 1.0e-12));
    CHECK_THAT(model(0.7).imag(), WithinRel(ref.imag(), 1.0e-12));
  }

  SECTION("Table interpolation and extrapolation")
  {
    // Two-fluid data: σ1 constant and ωσ2 constant are reproduced exactly between the
    // nodes, and held constant outside the table.
    config::ComplexConductivityData data;
    data.type = ComplexConductivityType::TABLE;
    data.table_f = {1.0, 2.0, 5.0, 10.0};
    data.table_sigma1 = {3.0, 3.0, 3.0, 3.0};
    data.table_sigma2 = {4.0, 2.0, 0.8, 0.4};
    const ComplexConductivityModel model(data);
    for (double omega : {1.0, 1.3, 2.0, 3.7, 9.9, 10.0})
    {
      CHECK_THAT(model(omega).real(), WithinRel(3.0, 1.0e-14));
      CHECK_THAT(-omega * model(omega).imag(), WithinRel(4.0, 1.0e-14));
    }
    CHECK_THAT(-0.5 * model(0.5).imag(), WithinRel(4.0, 1.0e-14));
    CHECK_THAT(-20.0 * model(20.0).imag(), WithinRel(4.0, 1.0e-14));

    // Nodes are interpolated, and the interpolant is monotone between monotone data.
    data.table_sigma1 = {1.0, 2.0, 2.5, 6.0};
    const ComplexConductivityModel model2(data);
    for (std::size_t i = 0; i < data.table_f.size(); i++)
    {
      CHECK_THAT(model2(data.table_f[i]).real(), WithinRel(data.table_sigma1[i], 1.0e-14));
    }
    double prev = model2(1.0).real();
    for (double omega = 1.05; omega <= 10.0; omega += 0.05)
    {
      const double s = model2(omega).real();
      CHECK(s >= prev - 1.0e-14);
      prev = s;
    }
  }
}

TEST_CASE("ComplexConductivity configuration", "[complexconductivity][Serial]")
{
  using json = nlohmann::json;

  SECTION("Mattis-Bardeen")
  {
    const auto j = json::parse(R"({"Attributes": [1], "ComplexConductivity":
        {"Model": "MattisBardeen", "NormalConductivity": 1.6e7, "Temperature": 4.2,
         "CriticalTemperature": 9.0, "ReferenceFrequency": 100.0}})");
    const config::MaterialData data(j);
    REQUIRE(data.complex_sigma);
    const auto &cs = *data.complex_sigma;
    CHECK(cs.IsDispersive());
    CHECK(cs.type == ComplexConductivityType::MATTIS_BARDEEN);
    CHECK(cs.sigma_n == 1.6e7);
    CHECK(cs.temperature == 4.2);
    CHECK(cs.critical_temperature == 9.0);
    CHECK(cs.reference_freq == 100.0);
    CHECK(data.sigma.s[0] == 0.0);
    CHECK(data.lambda_L == 0.0);
  }

  SECTION("Two-fluid models resolve into a conductivity and London depth")
  {
    const double sigma_n = 2.0e7, T = 4.5, Tc = 9.0, lambda_0 = 0.08, t = T / Tc;
    auto j = json::parse(R"({"Attributes": [1], "ComplexConductivity":
        {"Model": "TwoFluid", "TwoFluidModel": "Classical"}})");
    auto &cs = j["ComplexConductivity"];
    cs["NormalConductivity"] = sigma_n;
    cs["Temperature"] = T;
    cs["CriticalTemperature"] = Tc;
    cs["LondonDepth0"] = lambda_0;
    {
      const config::MaterialData data(j);
      REQUIRE(data.complex_sigma);
      CHECK_FALSE(data.complex_sigma->IsDispersive());
      CHECK_THAT(data.sigma.s[0], WithinRel(sigma_n * std::pow(t, 4), 1.0e-14));
      CHECK_THAT(data.sigma.s[2], WithinRel(sigma_n * std::pow(t, 4), 1.0e-14));
      CHECK_THAT(data.lambda_L,
                 WithinRel(lambda_0 / std::sqrt(1.0 - std::pow(t, 4)), 1.0e-14));
    }
    cs["TwoFluidModel"] = "Enhanced";
    {
      const config::MaterialData data(j);
      CHECK_THAT(data.sigma.s[0], WithinRel(sigma_n * std::sqrt(t), 1.0e-14));
      CHECK_THAT(data.lambda_L,
                 WithinRel(lambda_0 / std::sqrt(1.0 - std::pow(t, 1.5)), 1.0e-14));
    }
  }

  SECTION("Invalid input")
  {
    // Exclusive with the frequency-independent conductivity and London depth.
    CHECK_THROWS_WITH(config::MaterialData(json::parse(R"({"Attributes": [1],
        "Conductivity": 1.0, "ComplexConductivity": {"Model": "MattisBardeen",
        "NormalConductivity": 1.0, "Gap": 1.0e-3}})")),
                      Matchers::ContainsSubstring("cannot be combined"));
    // Mattis-Bardeen needs a gap or a critical temperature above the temperature.
    CHECK_THROWS(config::MaterialData(json::parse(R"({"Attributes": [1],
        "ComplexConductivity": {"Model": "MattisBardeen", "NormalConductivity": 1.0}})")));
    CHECK_THROWS(config::MaterialData(json::parse(R"({"Attributes": [1],
        "ComplexConductivity": {"Model": "MattisBardeen", "NormalConductivity": 1.0,
        "Temperature": 10.0, "CriticalTemperature": 9.0}})")));
    // Tables need strictly increasing frequencies.
    CHECK_THROWS(config::MaterialData(json::parse(R"({"Attributes": [1],
        "ComplexConductivity": {"Model": "Table", "Frequency": [2.0, 1.0],
        "Sigma1": [1.0, 1.0], "Sigma2": [1.0, 1.0]}})")));
  }
}

TEST_CASE("ComplexConductivity reference frequency", "[complexconductivity][Serial]")
{
  // IoData::CheckConfiguration resolves the reference frequency [GHz] from the problem
  // type, and NondimensionalizeInputs converts the frequencies and energies.
  constexpr double L0 = 1.0e-6, Lc = 100.0;
  const Units units(L0, Lc * L0);
  auto setup = [&](ProblemType type)
  {
    auto iodata = std::make_unique<IoData>(units);
    iodata->model.L0 = L0;
    iodata->model.Lc = Lc;
    iodata->problem.type = type;
    auto &mat = iodata->domains.materials.emplace_back();
    mat.attributes = {1};
    auto &cs = mat.complex_sigma.emplace();
    cs.sigma_n = 1.0e7;
    cs.temperature = 4.2;
    cs.gap = 1.4e-3;
    return iodata;
  };
  auto omega_nd = [&units](double f_GHz)
  {
    return 2 * std::numbers::pi *
           units.Nondimensionalize<Units::ValueType::FREQUENCY>(f_GHz);
  };
  std::unique_ptr<mfem::Mesh> no_mesh;

  SECTION("Driven: geometric mean of the sweep band, or the given value")
  {
    auto iodata = setup(ProblemType::DRIVEN);
    iodata->solver.driven.sample_f = {10.0, 250.0, 1000.0};
    iodata->CheckConfiguration();
    const auto &cs = *iodata->domains.materials[0].complex_sigma;
    CHECK_THAT(cs.reference_freq, WithinRel(100.0, 1.0e-14));
    iodata->NondimensionalizeInputs(no_mesh);
    CHECK_THAT(cs.reference_freq, WithinRel(omega_nd(100.0), 1.0e-14));
    CHECK_THAT(cs.gap, WithinRel(omega_nd(1.0e-9 * 1.4e-3 * electromagnetics::e_ /
                                          electromagnetics::h_),
                                 1.0e-14));
    CHECK_THAT(cs.temperature, WithinRel(omega_nd(1.0e-9 * 4.2 * electromagnetics::kB_ /
                                                  electromagnetics::h_),
                                         1.0e-14));
    CHECK_THAT(
        cs.sigma_n,
        WithinRel(1.0e7 / units.GetScaleFactor<Units::ValueType::CONDUCTIVITY>(), 1.0e-14));

    auto iodata2 = setup(ProblemType::DRIVEN);
    iodata2->solver.driven.sample_f = {10.0, 1000.0};
    iodata2->domains.materials[0].complex_sigma->reference_freq = 300.0;
    iodata2->CheckConfiguration();
    CHECK(iodata2->domains.materials[0].complex_sigma->reference_freq == 300.0);
  }

  SECTION("Boundary mode: the mode frequency, or the given value")
  {
    auto iodata = setup(ProblemType::BOUNDARYMODE);
    iodata->solver.boundary_mode.freq = 42.0;
    iodata->CheckConfiguration();
    CHECK(iodata->domains.materials[0].complex_sigma->reference_freq == 42.0);

    auto iodata2 = setup(ProblemType::BOUNDARYMODE);
    iodata2->solver.boundary_mode.freq = 42.0;
    iodata2->domains.materials[0].complex_sigma->reference_freq = 300.0;
    iodata2->CheckConfiguration();
    CHECK(iodata2->domains.materials[0].complex_sigma->reference_freq == 300.0);
  }

  SECTION("Unsupported problem types and table ranges")
  {
    for (auto type : {ProblemType::EIGENMODE, ProblemType::TRANSIENT})
    {
      CHECK_THROWS_WITH(setup(type)->CheckConfiguration(),
                        Matchers::ContainsSubstring("only supported for driven"));
    }
    auto synthesis = setup(ProblemType::DRIVEN);
    synthesis->solver.driven.sample_f = {10.0, 200.0};
    synthesis->solver.driven.adaptive_tol = 1.0e-6;
    synthesis->solver.driven.adaptive_circuit_synthesis = true;
    synthesis->boundaries.lumpedport[1];
    CHECK_THROWS_WITH(synthesis->CheckConfiguration(),
                      Matchers::ContainsSubstring("dispersive \"ComplexConductivity\""));
    auto table = [&](ProblemType type)
    {
      auto iodata = setup(type);
      auto &cs = *iodata->domains.materials[0].complex_sigma;
      cs.type = ComplexConductivityType::TABLE;
      cs.table_f = {1.0, 100.0};
      cs.table_sigma1 = {1.0, 1.0};
      cs.table_sigma2 = {1.0, 1.0};
      iodata->solver.driven.sample_f = {10.0, 200.0};
      iodata->solver.boundary_mode.freq = 200.0;
      return iodata;
    };
    for (auto type : {ProblemType::DRIVEN, ProblemType::BOUNDARYMODE})
    {
      CHECK_THROWS_WITH(
          table(type)->CheckConfiguration(),
          Matchers::ContainsSubstring("within the \"ComplexConductivity\" table"));
    }
  }
}

TEST_CASE("Dispersive conductivity anchor and residual",
          "[complexconductivity][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  auto serial_mesh = TwoAttributeMesh(comm);
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  Mesh palace_mesh(std::move(par_mesh));

  const double omega_ref = 0.7;
  config::MaterialData sc, dielectric;
  sc.attributes = {1};
  sc.complex_sigma = MattisBardeenData(5.0, 2.0, 0.4, omega_ref);
  dielectric.attributes = {2};
  dielectric.epsilon_r = config::SymmetricMatrixData<3>(4.0);
  const std::vector<config::MaterialData> materials = {sc, dielectric};
  config::PeriodicBoundaryData periodic;
  const ComplexConductivityModel model(*sc.complex_sigma);
  const auto sigma_ref = model(omega_ref);

  MaterialOperator mat_op(materials, periodic, ProblemType::DRIVEN, palace_mesh);
  CHECK(mat_op.HasConductivity());
  CHECK(mat_op.HasLondonDepth());
  if (palace_mesh.GetCeedAttributes().contains(1))
  {
    // Anchor: σ1(ω_ref) on the damping and ω_ref σ2(ω_ref) on the stiffness operator.
    const auto sigma = mat_op.GetConductivity(1);
    const auto london = mat_op.GetInvLondonDepth(1);
    for (int d = 0; d < 3; d++)
    {
      CHECK(sigma(d, d) == sigma_ref.real());
      CHECK(london(d, d) == -omega_ref * sigma_ref.imag());
    }
    CHECK(sigma(0, 1) == 0.0);
  }
  if (palace_mesh.GetCeedAttributes().contains(2))
  {
    CHECK(mat_op.GetConductivity(2).MaxMaxNorm() == 0.0);
    CHECK(mat_op.GetInvLondonDepth(2).MaxMaxNorm() == 0.0);
  }

  const Units units(1.0, 1.0);
  DispersiveConductivityOperator disp_op(materials, ProblemType::DRIVEN, units, mat_op);
  REQUIRE(disp_op.Size() == 1);
  // The residual vanishes exactly at the reference frequency.
  CHECK(disp_op.EvaluateScalar(0, omega_ref) == 0.0);
  // Elsewhere it is iω (σ(ω) - σ_a(ω)), σ_a(ω) = σ1(ω_ref) - i ω_ref σ2(ω_ref)/ω.
  for (double omega : {0.2, 1.5, 5.0})
  {
    const std::complex<double> sigma_a(sigma_ref.real(),
                                       omega_ref * sigma_ref.imag() / omega);
    const auto ref = 1i * omega * (model(omega) - sigma_a);
    const auto s = disp_op.EvaluateScalar(0, omega);
    CHECK_THAT(s.real(), WithinAbs(ref.real(), 1.0e-12 * std::abs(ref)));
    CHECK_THAT(s.imag(), WithinAbs(ref.imag(), 1.0e-12 * std::abs(ref)));
  }

  // Coefficients are added on the dispersive domain attributes only.
  const double omega = 1.5;
  MaterialPropertyCoefficient fr(mat_op.MaxCeedAttribute()), fi(mat_op.MaxCeedAttribute());
  disp_op.AddExtraSystemCoefficients(omega, fr, fi);
  const auto s = disp_op.EvaluateScalar(0, omega);
  for (auto [attr, value_r, value_i] :
       {std::tuple{1, s.real(), s.imag()}, std::tuple{2, 0.0, 0.0}})
  {
    const auto ceed_attrs = palace_mesh.GetCeedAttributes(attr);
    if (ceed_attrs.Size() == 0)
    {
      continue;
    }
    const int k_r = fr.GetAttributeToMaterial()[ceed_attrs[0] - 1];
    const int k_i = fi.GetAttributeToMaterial()[ceed_attrs[0] - 1];
    if (attr == 1)
    {
      REQUIRE(k_r >= 0);
      REQUIRE(k_i >= 0);
      CHECK(fr.GetMaterialProperties()(0, 0, k_r) == value_r);
      CHECK(fi.GetMaterialProperties()(0, 0, k_i) == value_i);
    }
    else
    {
      CHECK(k_r < 0);
      CHECK(k_i < 0);
    }
  }
}

TEST_CASE("SpaceOperator dispersive conductivity",
          "[complexconductivity][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  auto serial_mesh = TwoAttributeMesh(comm);
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(comm, *serial_mesh));

  const double omega_ref = 0.7;
  config::MaterialData sc, dielectric;
  sc.attributes = {1};
  sc.complex_sigma = MattisBardeenData(5.0, 2.0, 0.4, omega_ref);
  dielectric.attributes = {2};
  config::SolverData solver;
  solver.order = 2;
  IntegrationSettingsGuard settings_guard(solver);
  config::DomainData domains;
  domains.attributes = {1, 2};
  domains.materials = {sc, dielectric};
  config::BoundaryData boundaries;
  const Units units(1.0, 1.0);
  SpaceOperator space_op(solver, domains, boundaries, ProblemType::DRIVEN, units, mesh);
  REQUIRE(space_op.GetDispersiveConductivityOp().Size() == 1);

  ComplexVector x(space_op.GetNDSpace().GetTrueVSize()), y(x.Size()), z(x.Size());
  x.UseDevice(true);
  y.UseDevice(true);
  z.UseDevice(true);
  linalg::SetRandom(comm, x);

  // The residual vanishes at the reference frequency, but A2 is still returned so that
  // callers detect the frequency dependence.
  auto A2_ref =
      space_op.GetExtraSystemMatrix<ComplexOperator>(omega_ref, Operator::DIAG_ZERO);
  REQUIRE(A2_ref);
  A2_ref->Mult(x, y);
  CHECK(linalg::Norml2(comm, y) == 0.0);

  // Away from the reference frequency, A2(ω) = δy(ω) M_g with M_g the unit domain mass,
  // returned on the imaginary slot.
  auto M_g =
      space_op.GetDispersiveConductivityMassMatrix<ComplexOperator>(0, Operator::DIAG_ZERO);
  REQUIRE(M_g);
  M_g->Mult(x, z);
  const double norm_Mx = linalg::Norml2(comm, z);
  REQUIRE(norm_Mx > 0.0);
  for (double omega : {0.3, 1.5})
  {
    const auto s = space_op.GetDispersiveConductivityOp().EvaluateScalar(0, omega);
    auto A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO);
    REQUIRE(A2);
    A2->Mult(x, y);
    M_g->Mult(x, z);
    z *= s / 1i;
    y -= z;
    CHECK(linalg::Norml2(comm, y) < 1.0e-12 * std::abs(s) * norm_Mx);
  }

  // Preconditioner assembly with the residual term.
  CHECK(space_op.GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, 1.5i, -2.25 + 0.0i,
                                                          1.5) != nullptr);
}

}  // namespace palace
