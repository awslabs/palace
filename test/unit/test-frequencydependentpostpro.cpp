// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Consumers of the complete permittivity ε(iω) and the stored-energy permittivity W(ω) of
// frequency-dependent materials: material evaluation, postprocessing, and the flux error
// estimators. Expected values come from the model formulas, independent of the
// implementation (internal units, ε₀ = 1).

#include <array>
#include <cmath>
#include <complex>
#include <fstream>
#include <memory>
#include <numbers>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/coefficient.hpp"
#include "fem/errorindicator.hpp"
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "fem/output_functionals.hpp"
#include "fem/point_field_evaluator.hpp"
#include "fixtures.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/farfieldboundaryoperator.hpp"
#include "models/floquetportoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/romoperator.hpp"
#include "models/spaceoperator.hpp"
#include "models/surfacepostoperator.hpp"
#include "test-helpers.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/constants.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

namespace palace
{

using namespace Catch::Matchers;
using namespace std::complex_literals;

namespace
{

constexpr double eps_inf = 2.0;

// Model terms with their closed-form susceptibility χ(iω) and stored-energy
// susceptibility w(ω) (permittivitymodel.hpp: Re χ for a relaxation, plus the kinetic
// energy of an oscillator).
struct Model
{
  const char *name;
  std::vector<config::PermittivityTermData> terms;
  std::complex<double> (*chi)(double);
  double (*w)(double);
};

// Drude ωp = 3, γ = 0.5: χ = ωp² / (s (s + γ)), w = ωp² / (ω² + γ²). Its zero pole is
// folded into the conductivity of the solver, and Re ε < 0 below ω ≈ ωp / √ε∞.
std::complex<double> DrudeChi(double w)
{
  const std::complex<double> s(0.0, w);
  return 9.0 / (s * (s + 0.5));
}
double DrudeW(double w)
{
  return 9.0 / (w * w + 0.25);
}

// Lorentz Δε = 1.5, ω₀ = 1, γ = 0.2: χ = Δε ω₀² / (s² + γ s + ω₀²),
// w = Δε ω₀² (ω² + ω₀²) / ((ω₀² − ω²)² + γ² ω²).
std::complex<double> LorentzChi(double w)
{
  const std::complex<double> s(0.0, w);
  return 1.5 / (s * s + 0.2 * s + 1.0);
}
double LorentzW(double w)
{
  return 1.5 * (w * w + 1.0) / ((1.0 - w * w) * (1.0 - w * w) + 0.04 * w * w);
}

// Debye Δε = 3, τ = 0.7 plus Djordjevic-Sarkar (0.4, 0.1, 50): relaxations, w = Re χ.
std::complex<double> DebyeDSChi(double w)
{
  const std::complex<double> s(0.0, w);
  return 3.0 / (1.0 + 0.7 * s) + 0.4 * std::log((s + 50.0) / (s + 0.1));
}
double DebyeDSW(double w)
{
  return DebyeDSChi(w).real();
}

// A pure zero pole r / s (r = 2), a free-carrier conductivity of the model: no nonlinear
// term, no stored energy.
std::complex<double> ZeroPoleChi(double w)
{
  return 2.0 / std::complex<double>(0.0, w);
}
double ZeroPoleW(double)
{
  return 0.0;
}

std::vector<Model> Models()
{
  return {
      {"Drude", {config::DrudeData{3.0, 0.5}}, DrudeChi, DrudeW},
      {"Lorentz", {config::LorentzData{1.5, 1.0, 0.2}}, LorentzChi, LorentzW},
      {"DebyeDS",
       {config::DebyeData{3.0, 0.7}, config::DjordjevicSarkarData{0.4, 0.1, 50.0}},
       DebyeDSChi,
       DebyeDSW},
      {"ZeroPole", {config::PermittivityPoleData{0.0, 2.0}}, ZeroPoleChi, ZeroPoleW},
  };
}

config::MaterialData MakeMaterial(const std::vector<int> &attrs,
                                  const std::vector<config::PermittivityTermData> &terms)
{
  config::MaterialData material;
  material.attributes = attrs;
  material.epsilon_r.s = {eps_inf, eps_inf, eps_inf};
  material.permittivity_terms = terms;
  return material;
}

Mesh MakeBoxMesh(int dim, int n, mfem::Element::Type type, double sx = 1.0, double sy = 1.3,
                 double sz = 0.8)
{
  auto serial = std::make_unique<mfem::Mesh>(
      (dim == 3) ? mfem::Mesh::MakeCartesian3D(n, n, n, type, sx, sy, sz)
                 : mfem::Mesh::MakeCartesian2D(n, n, type, true, sx, sy));
  // Attribute 1 on the first element only: on several processes, its material has no
  // elements on most of them.
  for (int i = 0; i < serial->GetNE(); i++)
  {
    serial->SetAttribute(i, (i == 0) ? 1 : 2);
  }
  serial->SetAttributes();
  return Mesh(std::make_unique<mfem::ParMesh>(Mpi::World(), *serial));
}

}  // namespace

TEST_CASE("Complete frequency-dependent permittivity and stored energy",
          "[frequencydependent][Serial][Parallel]")
{
  auto mesh = MakeBoxMesh(3, 2, mfem::Element::HEXAHEDRON);
  config::PeriodicBoundaryData periodic;
  for (const auto &model : Models())
  {
    CAPTURE(model.name);
    config::MaterialData lossy;
    lossy.attributes = {2};
    lossy.epsilon_r.s = {3.0, 3.0, 3.0};
    lossy.tandelta.s = {0.1, 0.1, 0.1};
    MaterialOperator mat_op({MakeMaterial({1}, model.terms), lossy}, periodic,
                            ProblemType::DRIVEN, mesh);
    REQUIRE(mat_op.HasFrequencyDependentPermittivityModel(0));
    CHECK_FALSE(mat_op.HasFrequencyDependentPermittivityModel(1));
    // The pure zero pole has no nonlinear term, but its model is still configured.
    CHECK(mat_op.HasFrequencyDependentPermittivitySupport(0) ==
          (std::string(model.name) != "ZeroPole"));
    for (double omega : {0.3, 1.1, 4.0})
    {
      CAPTURE(omega);
      const std::complex<double> eps = eps_inf + model.chi(omega);
      const auto value = mat_op.EvaluateFrequencyDependentPermittivity(0, {0.0, omega});
      CHECK_THAT(value.real(), WithinRel(eps.real(), 1.0e-12) || WithinAbs(0.0, 1.0e-12));
      CHECK_THAT(value.imag(), WithinRel(eps.imag(), 1.0e-12) || WithinAbs(0.0, 1.0e-12));
      const double W = eps_inf + model.w(omega);
      CHECK_THAT(mat_op.GetFrequencyDependentStoredEnergyPermittivity(0, omega),
                 WithinRel(W, 1.0e-12));

      // Tables in the layout of the solver getters, which keep ε∞.
      const auto eps_r = mat_op.EvaluatePermittivityReal(omega);
      const auto eps_i = mat_op.EvaluatePermittivityImag(omega);
      const auto eps_w = mat_op.EvaluateStoredEnergyPermittivity(omega);
      if (mesh.GetCeedAttributes(1).Size() > 0)
      {
        const int k = mat_op.GetLocalMaterial(1);  // Global attribute, present here
        for (int d = 0; d < 3; d++)
        {
          CHECK_THAT(eps_r(d, d, k), WithinRel(eps.real(), 1.0e-12));
          CHECK_THAT(eps_i(d, d, k), WithinRel(eps.imag(), 1.0e-12));
          CHECK_THAT(eps_w(d, d, k), WithinRel(W, 1.0e-12));
          CHECK(mat_op.GetPermittivityReal()(d, d, k) == eps_inf);
        }
        CHECK(eps_r(0, 1, k) == 0.0);
      }
      if (mesh.GetCeedAttributes(2).Size() > 0)
      {
        const int k = mat_op.GetLocalMaterial(2);
        CHECK(eps_r(0, 0, k) == 3.0);
        CHECK(eps_i(0, 0, k) == Catch::Approx(-0.3));
        CHECK(eps_w(0, 0, k) == 3.0);
      }
    }
  }

  // A singular evaluation is an error, not a silent infinity.
  MaterialOperator lossless({MakeMaterial({1, 2}, {config::LorentzData{1.5, 1.0, 0.0}})},
                            periodic, ProblemType::DRIVEN, mesh);
  CHECK_THROWS_WITH(lossless.EvaluateFrequencyDependentPermittivity(0, {0.0, 1.0}),
                    ContainsSubstring("singular"));
}

TEST_CASE("Frequency-dependent material identity and exterior wave speed",
          "[frequencydependent][Serial][Parallel]")
{
  auto mesh = MakeBoxMesh(3, 2, mfem::Element::HEXAHEDRON);
  config::PeriodicBoundaryData periodic;
  // A frequency-dependent material with ε∞ = 1 is not vacuum for interface
  // postprocessing, at any frequency.
  auto unit = MakeMaterial({1}, {config::DrudeData{3.0, 0.5}});
  unit.epsilon_r.s = {1.0, 1.0, 1.0};
  config::MaterialData vacuum;
  vacuum.attributes = {2};
  MaterialOperator mat_op({unit, vacuum}, periodic, ProblemType::DRIVEN, mesh);
  for (int attr : {1, 2})
  {
    if (mesh.GetCeedAttributes(attr).Size() > 0)
    {
      CHECK(mat_op.IsVacuum(attr) == (attr == 2));
      CHECK(mat_op.GetLightSpeedMax(attr) == Catch::Approx(1.0));
    }
  }

  // The wave speed of a lossless exterior uses Re ε(iω), and requires it positive: on every
  // process, also where the material has no elements.
  const double omega = 4.0;
  const double eps = (1.0 + DrudeChi(omega)).real();
  REQUIRE(eps > 0.0);
  const auto c0 = mat_op.EvaluateLightSpeed(omega, {0});
  const auto c0_inf = mat_op.EvaluateLightSpeed(omega, {});
  if (mesh.GetCeedAttributes(1).Size() > 0)
  {
    const int k = mat_op.GetLocalMaterial(1);  // Global attribute, present here
    CHECK_THAT(c0(0, 0, k), WithinRel(1.0 / std::sqrt(eps), 1.0e-12));
    CHECK(c0_inf(0, 0, k) == 1.0);
  }
  CHECK_THROWS_WITH(mat_op.EvaluateLightSpeed(1.0, {0}), ContainsSubstring("Re ε(iω) > 0"));
}

TEST_CASE("Absorbing boundaries next to Drude and zero-pole materials",
          "[frequencydependent][Serial][Parallel]")
{
  // The absorbing boundary scale r(ω) = √(ε_w(iω)/ε∞) of the neighboring material uses its
  // wave permittivity ε_w = ε + σ/s: the complete permittivity, including the zero pole
  // folded into the conductivity (whose omission dropped the loss and reversed the sign of
  // the boundary loss), plus the user conductivity σ, once (not the conductivity matrix,
  // which also holds the folded zero pole). A configured model with only a zero pole has no
  // nonlinear term but still a frequency-dependent boundary. At complex ω, ε_w is continued
  // analytically. The displacement and stored energy keep ε and W, without σ.
  auto mesh = MakeBoxMesh(3, 2, mfem::Element::HEXAHEDRON);
  config::PeriodicBoundaryData periodic;
  for (const auto &model : {Models()[0], Models()[3]})
  {
    for (double sigma : {0.0, 1.0})
    {
      CAPTURE(model.name, sigma);
      auto material = MakeMaterial({1, 2}, model.terms);
      material.sigma.s = {sigma, sigma, sigma};
      MaterialOperator mat_op({material}, periodic, ProblemType::DRIVEN, mesh);
      config::FarfieldBoundaryData farfield;
      farfield.attributes = {1, 2, 3, 4, 5, 6};
      FarfieldBoundaryOperator farfield_op(farfield, ProblemType::DRIVEN, mat_op,
                                           mesh.Get());
      REQUIRE(farfield_op.NumDispersiveMaterials() == 1);
      for (std::complex<double> omega : {std::complex<double>(4.0, 0.0), {2.0, 0.3}})
      {
        CAPTURE(omega);
        const std::complex<double> s = 1i * omega;
        std::complex<double> chi = 0.0;
        if (std::string(model.name) == "Drude")
        {
          chi = 9.0 / (s * (s + 0.5));
        }
        else
        {
          chi = 2.0 / s;
        }
        const auto r = std::sqrt((eps_inf + chi + sigma / s) / eps_inf);
        const auto value = farfield_op.GetDispersiveScale(0, omega);
        CHECK_THAT(value.real(), WithinRel(r.real(), 1.0e-12));
        CHECK_THAT(value.imag(), WithinRel(r.imag(), 1.0e-12));
        const auto eps = mat_op.EvaluateFrequencyDependentPermittivity(0, s);
        CHECK(std::abs(eps - (eps_inf + chi)) < 1.0e-12 * std::abs(eps));
        // Passive at real ω: the boundary dissipates, Re{iω r} ≤ 0 is the loss sign of the
        // first-order term i ω r (here r has Im r < 0).
        if (omega.imag() == 0.0)
        {
          CHECK(value.imag() < 0.0);
        }
      }
    }
  }

  // Debye ε∞ = 2, Δε = 2, τ = 1 with σ = 1 at ω = 1: the displacement and stored energy use
  // ε = 3 − i and W = 3, the wave ε_w = 3 − 2i.
  auto debye = MakeMaterial({1, 2}, {config::DebyeData{2.0, 1.0}});
  debye.sigma.s = {1.0, 1.0, 1.0};
  MaterialOperator debye_op({debye}, periodic, ProblemType::DRIVEN, mesh);
  const std::complex<double> s = 1i;
  CHECK(std::abs(debye_op.EvaluateFrequencyDependentPermittivity(0, s) -
                 std::complex<double>(3.0, -1.0)) < 1.0e-14);
  CHECK(std::abs(debye_op.EvaluateFrequencyDependentWavePermittivity(0, s) -
                 std::complex<double>(3.0, -2.0)) < 1.0e-14);
  CHECK_THAT(debye_op.GetFrequencyDependentStoredEnergyPermittivity(0, 1.0),
             WithinRel(3.0, 1.0e-14));

  // An anisotropic user conductivity has no scalar wave impedance: the boundaries which
  // need it are an error, not a silent approximation; other consumers are unaffected.
  auto anisotropic = debye;
  anisotropic.sigma.s = {1.0, 2.0, 1.0};
  MaterialOperator anisotropic_op({anisotropic}, periodic, ProblemType::DRIVEN, mesh);
  CHECK(std::abs(anisotropic_op.EvaluateFrequencyDependentPermittivity(0, s) -
                 std::complex<double>(3.0, -1.0)) < 1.0e-14);
  CHECK_THROWS_WITH(anisotropic_op.EvaluateFrequencyDependentWavePermittivity(0, s),
                    ContainsSubstring("scalar Conductivity"));
}

TEST_CASE("Floquet propagation constants in lossy and gain media", "[floquetport][Serial]")
{
  using F = FloquetPortData;
  // Lossless propagating and evanescent orders keep their established roots.
  CHECK(F::OutgoingPropagationConstant({4.0, 0.0}) == std::complex<double>(2.0, 0.0));
  const auto evanescent = F::OutgoingPropagationConstant({-4.0, 0.0});
  CHECK(evanescent.real() == 0.0);
  CHECK(evanescent.imag() == -2.0);
  // A passive lossy medium decays in the outgoing direction.
  const auto lossy = F::OutgoingPropagationConstant({4.0, -1.0});
  CHECK(lossy.real() > 0.0);
  CHECK(lossy.imag() < 0.0);
  CHECK(std::abs(lossy * lossy - std::complex<double>(4.0, -1.0)) < 1.0e-14);
  // A gain medium keeps the outgoing direction (Re γ > 0) and grows.
  const auto gain = F::OutgoingPropagationConstant({4.0, 1.0});
  CHECK(gain.real() > 0.0);
  CHECK(gain.imag() > 0.0);
  CHECK(std::abs(gain * gain - std::complex<double>(4.0, 1.0)) < 1.0e-14);
  // DtN eigenvalues: TE γ, TM ω²με/γ, and their real parts are the power per amplitude.
  const double omega = 2.0;
  const std::complex<double> mu_eps(1.5, -0.2);
  CHECK(F::DtNEigenvalue(true, lossy, omega, mu_eps) == lossy);
  const auto tm = F::DtNEigenvalue(false, lossy, omega, mu_eps);
  CHECK(std::abs(tm - omega * omega * mu_eps / lossy) < 1.0e-14);
}

namespace
{

// Access to the factored-term self-check of the PROM online assembly.
class FloquetRomOperator : public RomOperator
{
public:
  using RomOperator::RomOperator;
  bool SelfChecked() const { return other_A2_self_checked; }
  bool FactoredOk() const { return other_A2_factored_ok; }
};

}  // namespace

TEST_CASE_METHOD(test::SharedTempDir, "Floquet ports in a lossy frequency-dependent medium",
                 "[floquetport][frequencydependent][Serial][Parallel]")
{
  // A plane wave at normal incidence through a short slab between two Floquet ports at
  // z = 0 (excited) and z = L in a homogeneous medium, periodic in x and y. The ports use
  // the complete ε(iω) of the medium in γ₀ = ω √(μ ε) (the outgoing root), in their Robin
  // term, and in the excitation, so the driven port is matched, S11 = 0, and the
  // transmission is S21 = e^{−iγ₀ L}. This holds for a lossless constant medium
  // (|S21| = 1), a lossy Drude medium with Re ε > 0, and a Drude medium below its plasma
  // frequency with Re ε < 0 and Re γ₀² < 0, whose (0,0) mode still carries power,
  // Re γ₀ > 0. With a user conductivity σ the wave uses ε_w = ε + σ/(iω) (the volume
  // includes σ, so the ports must too, once: Debye with ε = 3 − i and ε_w = 3 − 2i, and
  // Drude below its plasma frequency, whose folded zero pole must not be counted again).
  // The PROM, with the solution in its basis, reproduces it through its factored Floquet
  // Robin term.
  using json = nlohmann::json;
  enum class Medium
  {
    LOSSLESS,
    LOSSY,
    NEGATIVE,
    DEBYE_SIGMA,
    NEGATIVE_SIGMA
  };
  const auto medium = GENERATE(Medium::LOSSLESS, Medium::LOSSY, Medium::NEGATIVE,
                               Medium::DEBYE_SIGMA, Medium::NEGATIVE_SIGMA);
  CAPTURE(static_cast<int>(medium));
  MPI_Comm comm = Mpi::World();
  IoData iodata{Units(1.0e-3, 1.0e-3)};
  iodata.model.L0 = 1.0e-3;
  iodata.model.Lc = 1.0;
  iodata.problem.type = ProblemType::DRIVEN;
  iodata.solver.order = 3;
  iodata.solver.driven.sample_f = {100.0};  // GHz
  // Drude ε = 2 + (fp/f)² / (i (i + fγ/f)) at f = 100 GHz: fp = 100 GHz, fγ = 50 GHz gives
  // 1.2 − 0.4i, and fp = 300 GHz gives −5.2 − 3.6i.
  // Debye Δε = 2 with ωτ = 1 at 100 GHz gives ε = 3 − i. The user conductivity
  // σ = ω ε₀ adds −i to the wave permittivity.
  json permittivity = 2.0;
  std::complex<double> eps = 2.0;
  const bool with_sigma =
      (medium == Medium::DEBYE_SIGMA || medium == Medium::NEGATIVE_SIGMA);
  if (medium == Medium::DEBYE_SIGMA)
  {
    permittivity = {
        {"HighFrequency", 2.0},
        {"Terms",
         json::array({{{"Type", "Debye"},
                       {"DeltaPermittivity", 2.0},
                       {"RelaxationTime", 1.0e9 / (2.0 * std::numbers::pi * 1.0e11)}}})}};
    eps = {3.0, -1.0};
  }
  else if (medium != Medium::LOSSLESS)
  {
    const double fp = (medium == Medium::LOSSY) ? 100.0 : 300.0;
    permittivity = {{"HighFrequency", 2.0},
                    {"Terms", json::array({{{"Type", "Drude"},
                                            {"PlasmaFrequency", fp},
                                            {"CollisionFrequency", 50.0}}})}};
    eps = (medium == Medium::LOSSY) ? std::complex<double>(1.2, -0.4)
                                    : std::complex<double>(-5.2, -3.6);
  }
  json material = {{"Attributes", json::array({1})}, {"Permittivity", permittivity}};
  if (with_sigma)
  {
    material["Conductivity"] =
        2.0 * std::numbers::pi * 1.0e11 * electromagnetics::epsilon0_;
    eps -= 1i;
  }
  iodata.domains = config::DomainData(json{{"Materials", json::array({material})}});
  // MakeCartesian3D boundary attributes: 1 (z = 0), 6 (z = L), 5 and 3 (x faces), 2 and
  // 4 (y faces).
  iodata.boundaries = config::BoundaryData(json{
      {"Periodic",
       {{"BoundaryPairs", json::array({{{"DonorAttributes", json::array({5})},
                                        {"ReceiverAttributes", json::array({3})}},
                                       {{"DonorAttributes", json::array({2})},
                                        {"ReceiverAttributes", json::array({4})}}})}}},
      {"FloquetPort",
       json::array({{{"Index", 1}, {"Attributes", json::array({1})}, {"Excitation", 1}},
                    {{"Index", 2}, {"Attributes", json::array({6})}}})}});
  iodata.CheckConfiguration();
  constexpr double a = 0.5, L = 0.25;
  // The periodic mesh is made by the mesh loading path from the boundary pairs.
  iodata.model.mesh = (temp_dir / "floquet_box.mesh").string();
  if (Mpi::Root(comm))
  {
    auto box = mfem::Mesh::MakeCartesian3D(3, 3, 6, mfem::Element::HEXAHEDRON, a, a, L);
    std::ofstream ofs(iodata.model.mesh);
    box.Print(ofs);
  }
  Mpi::Barrier(comm);
  auto smesh = mesh::Load(iodata, comm);
  iodata.NondimensionalizeInputs(smesh);
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(mesh::Partition(iodata, std::move(smesh), comm)));
  SpaceOperator space_op(iodata, mesh);

  const double omega = iodata.solver.driven.sample_f.front();
  const std::complex<double> gamma = omega * std::sqrt(eps);
  REQUIRE(gamma.real() > 0.0);
  if (medium == Medium::NEGATIVE || medium == Medium::NEGATIVE_SIGMA)
  {
    REQUIRE((gamma * gamma).real() < 0.0);
  }
  auto &floquet_op = space_op.GetFloquetPortOp();
  floquet_op.Initialize(omega);
  const auto gamma0 = floquet_op.GetPort(1).GetGamma0();
  CHECK_THAT(gamma0.real(), WithinRel(gamma.real(), 1.0e-10));
  CHECK_THAT(gamma0.imag(), WithinAbs(gamma.imag(), 1.0e-10));

  const std::complex<double> S21_expected = std::exp(-1i * gamma * L);
  auto CheckS = [&](const ComplexVector &E, const char *path)
  {
    INFO(path);
    GridFunction E_gf(space_op.GetNDSpace(), true);
    E_gf.Real().SetFromTrueDofs(E.Real());
    E_gf.Imag().SetFromTrueDofs(E.Imag());
    const auto key = std::make_tuple(0, 0, true);
    const auto S1 = floquet_op.GetPort(1).GetAllSParameters(E_gf, true);
    const auto S2 = floquet_op.GetPort(2).GetAllSParameters(E_gf);
    REQUIRE(S1.contains(key));
    REQUIRE(S2.contains(key));
    const auto S11 = S1.at(key), S21 = S2.at(key);
    CAPTURE(S11, S21, S21_expected);
    CHECK(std::abs(S11) < 1.0e-3);
    CHECK(std::abs(S21 - S21_expected) < 1.0e-3);
  };

  // HDM: (K + iωC − ω²M + A2(ω)) E = RHS.
  ComplexVector E;
  {
    auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
    auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto A2 = space_op.GetExtraSystemOperator(omega, Operator::DIAG_ZERO);
    auto A = space_op.GetSystemMatrix(1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i,
                                      K.get(), C.get(), M.get(), A2.get());
    auto P = space_op.GetPreconditionerMatrix<ComplexOperator>(
        1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i, omega);
    ComplexKspSolver ksp(iodata, space_op.GetNDSpaces(), &space_op.GetH1Spaces());
    ksp.SetOperators(*A, *P);
    ComplexVector RHS;
    RHS.UseDevice(true);
    E.UseDevice(true);
    REQUIRE(space_op.GetExcitationVector(1, omega, RHS));
    E.SetSize(RHS.Size());
    E = 0.0;
    ksp.Mult(RHS, E);
  }
  CheckS(E, "HDM");
  if (medium != Medium::LOSSLESS)
  {
    CHECK(std::abs(S21_expected) < 0.95);
  }

  // PROM with the HDM solution in its basis.
  FloquetRomOperator prom_op(iodata, space_op, 4);
  prom_op.SetExcitationIndex(1);
  ComplexVector E_hdm(E.Size()), E_rom(E.Size());
  E_hdm = 0.0;
  E_rom = 0.0;
  E_hdm.UseDevice(true);
  E_rom.UseDevice(true);
  prom_op.SolveHDM(1, omega, E_hdm);
  prom_op.UpdatePROM(E_hdm, "sample");
  prom_op.SolvePROM(1, omega, E_rom);
  CHECK(prom_op.SelfChecked());
  CHECK(prom_op.FactoredOk());
  CheckS(E_rom, "PROM");
}

TEST_CASE("Electric surface flux of frequency-dependent materials",
          "[frequencydependent][Serial][Parallel][GPU]")
{
  // A uniform complex field E = E_n n through the face of the box: Φ = ε(iω) E_n A with
  // the complete permittivity, for the 3D libCEED and the 2D legacy paths, refreshed when
  // the frequency changes and back.
  const int dim = GENERATE(3, 2);
  CAPTURE(dim);
  test::IntegrationSettingsGuard guard;
  constexpr int order = 1;
  fem::DefaultIntegrationOrder::p_trial = order;
  auto mesh = MakeBoxMesh(
      dim, 3, (dim == 3) ? mfem::Element::TETRAHEDRON : mfem::Element::TRIANGLE);
  config::PeriodicBoundaryData periodic;
  const auto model = Models()[0];  // Drude
  MaterialOperator mat_op({MakeMaterial({1}, model.terms), MakeMaterial({2}, model.terms)},
                          periodic, ProblemType::DRIVEN, mesh);
  mfem::H1_FECollection h1_fec(order, dim);
  mfem::ND_FECollection nd_fec(order, dim);
  FiniteElementSpace h1_fespace(mesh, &h1_fec), nd_fespace(mesh, &nd_fec);
  GridFunction E(nd_fespace, true);
  const std::complex<double> En(1.0, 0.5);
  {
    mfem::Vector vr(dim), vi(dim);
    vr = 0.0;
    vi = 0.0;
    vr(dim - 1) = En.real();
    vi(dim - 1) = En.imag();
    mfem::VectorConstantCoefficient cr(vr), ci(vi);
    E.Real().ProjectCoefficient(cr);
    E.Imag().ProjectCoefficient(ci);
  }
  // Top face: z = 0.8 (attribute 6) in 3D, y = 1.3 (attribute 3) in 2D.
  config::BoundaryPostData postpro;
  config::SurfaceFluxData flux;
  flux.type = SurfaceFlux::ELECTRIC;
  flux.center = {0.5, 0.65, 0.4};
  flux.no_center = false;
  flux.attributes = {(dim == 3) ? 6 : 3};
  postpro.flux.emplace(1, flux);
  SurfacePostOperator surf_post_op(postpro, ProblemType::DRIVEN, mat_op, h1_fespace,
                                   nd_fespace);
  const double area = (dim == 3) ? 1.0 * 1.3 : 1.0;
  CHECK_THROWS(surf_post_op.GetSurfaceFlux(1, &E, nullptr));  // Requires a frequency
  for (double omega : {4.0, 1.1, 4.0})
  {
    CAPTURE(omega);
    const std::complex<double> expected = (eps_inf + model.chi(omega)) * En * area;
    const auto value = surf_post_op.GetSurfaceFlux(1, &E, nullptr, omega);
    CHECK_THAT(value.real(), WithinRel(expected.real(), 1.0e-10));
    CHECK_THAT(value.imag(), WithinRel(expected.imag(), 1.0e-10));
  }
}

TEST_CASE("Metal-substrate interface of a frequency-dependent substrate",
          "[frequencydependent][Serial][Parallel][GPU]")
{
  // The MS participation of a uniform normal field E = E_n n on the substrate side of the
  // face is 1/2 t / ε_MS |D_n|² A with |D_n| = |ε(iω)| |E_n|: the complete displacement,
  // not Re ε² of a lossy model. The substrate has ε∞ = 1, so it is identified as substrate
  // by its configured model, not by its light speed. The other side is vacuum.
  const int dim = GENERATE(3, 2);
  CAPTURE(dim);
  test::IntegrationSettingsGuard guard;
  constexpr int order = 1;
  fem::DefaultIntegrationOrder::p_trial = order;
  auto serial = std::make_unique<mfem::Mesh>(
      (dim == 3) ? mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::HEXAHEDRON)
                 : mfem::Mesh::MakeCartesian2D(2, 2, mfem::Element::QUADRILATERAL));
  // Substrate below the interface (last coordinate < 1/2, attribute 1), vacuum above, and
  // the interface boundary elements with a new attribute.
  for (int i = 0; i < serial->GetNE(); i++)
  {
    mfem::Vector c;
    serial->GetElementCenter(i, c);
    serial->SetAttribute(i, (c(dim - 1) < 0.5) ? 1 : 2);
  }
  const int interface_attr = serial->bdr_attributes.Max() + 1;
  for (int f = 0; f < serial->GetNumFaces(); f++)
  {
    int e1, e2;
    serial->GetFaceElements(f, &e1, &e2);
    if (e1 >= 0 && e2 >= 0 && serial->GetAttribute(e1) != serial->GetAttribute(e2))
    {
      auto *be = serial->GetFace(f)->Duplicate(serial.get());
      be->SetAttribute(interface_attr);
      serial->AddBdrElement(be);
    }
  }
  serial->FinalizeTopology();
  serial->Finalize();
  serial->SetAttributes();
  Mesh mesh(std::make_unique<mfem::ParMesh>(Mpi::World(), *serial));
  config::PeriodicBoundaryData periodic;
  auto substrate = MakeMaterial({1}, Models()[0].terms);
  substrate.epsilon_r.s = {1.0, 1.0, 1.0};
  config::MaterialData vacuum;
  vacuum.attributes = {2};
  MaterialOperator mat_op({substrate, vacuum}, periodic, ProblemType::DRIVEN, mesh);
  mfem::H1_FECollection h1_fec(order, dim);
  mfem::ND_FECollection nd_fec(order, dim);
  FiniteElementSpace h1_fespace(mesh, &h1_fec), nd_fespace(mesh, &nd_fec);
  GridFunction E(nd_fespace, true);
  const std::complex<double> En(0.8, -0.3);
  {
    mfem::Vector vr(dim), vi(dim);
    vr = 0.0;
    vi = 0.0;
    vr(dim - 1) = En.real();
    vi(dim - 1) = En.imag();
    mfem::VectorConstantCoefficient cr(vr), ci(vi);
    E.Real().ProjectCoefficient(cr);
    E.Imag().ProjectCoefficient(ci);
  }
  config::BoundaryPostData postpro;
  config::InterfaceDielectricData data;
  data.type = InterfaceDielectric::MS;
  data.t = 0.01;
  data.epsilon_r = 5.0;
  data.tandelta = 1.0e-3;
  data.attributes = {interface_attr};
  postpro.dielectric.emplace(1, data);
  SurfacePostOperator surf_post_op(postpro, ProblemType::DRIVEN, mat_op, h1_fespace,
                                   nd_fespace);
  for (double omega : {4.0, 1.1, 4.0})
  {
    CAPTURE(omega);
    const double eps_abs = std::abs(1.0 + DrudeChi(omega));
    const double expected =
        0.5 * data.t / data.epsilon_r * eps_abs * eps_abs * std::norm(En);
    CHECK_THAT(surf_post_op.GetInterfaceElectricFieldEnergy(1, E, omega),
               WithinRel(expected, 1.0e-10));
  }
}

TEST_CASE("Electric energy density of frequency-dependent materials",
          "[frequencydependent][Serial][Parallel][GPU]")
{
  // U_e = 1/2 W(ω) |E|² for a uniform complex field, from the libCEED domain evaluator and
  // the legacy coefficient, refreshed when the frequency changes and back.
  test::IntegrationSettingsGuard guard;
  constexpr int order = 1;
  fem::DefaultIntegrationOrder::p_trial = order;
  auto mesh = MakeBoxMesh(3, 2, mfem::Element::HEXAHEDRON);
  config::PeriodicBoundaryData periodic;
  const auto model = Models()[1];  // Lorentz
  MaterialOperator mat_op({MakeMaterial({1}, model.terms), MakeMaterial({2}, model.terms)},
                          periodic, ProblemType::DRIVEN, mesh);
  mfem::ND_FECollection nd_fec(order, 3);
  FiniteElementSpace nd_fespace(mesh, &nd_fec);
  GridFunction E(nd_fespace, true);
  {
    mfem::Vector vr({0.3, -0.2, 1.0}), vi({0.1, 0.4, -0.5});
    mfem::VectorConstantCoefficient cr(vr), ci(vi);
    E.Real().ProjectCoefficient(cr);
    E.Imag().ProjectCoefficient(ci);
  }
  const double E2 = 0.09 + 0.04 + 1.0 + 0.01 + 0.16 + 0.25;
  mfem::L2_FECollection viz_fec(order, 3);
  mfem::ParFiniteElementSpace viz(&mesh.Get(), &viz_fec);
  PointFieldEvaluator eval(PointFieldEvaluator::Kind::ENERGY_E, mesh, mat_op,
                           E.ParFESpace(), nullptr, viz, 1.0);
  REQUIRE(eval.IsValid());
  EnergyDensityCoefficient<EnergyDensityType::ELECTRIC> legacy(E, mat_op);
  mfem::ParGridFunction val(&viz), ref(&viz);
  CHECK_THROWS(eval.Eval(&E, nullptr, val));  // Requires a frequency
  // Below and above the resonance, where Re ε < 0 but W > 0.
  for (double omega : {0.5, 1.3, 0.5})
  {
    CAPTURE(omega);
    const double expected = 0.5 * (eps_inf + LorentzW(omega)) * E2;
    eval.SetFrequency(omega);
    eval.Eval(&E, nullptr, val);
    legacy.SetStoredEnergyPermittivity(mat_op.EvaluateStoredEnergyPermittivity(omega));
    ref.ProjectCoefficient(legacy);
    for (int i = 0; i < val.Size(); i++)
    {
      CHECK_THAT(val(i), WithinRel(expected, 1.0e-10));
      CHECK_THAT(ref(i), WithinRel(expected, 1.0e-10));
    }
  }
}

namespace
{

// The gradient flux estimator in a homogeneous frequency-dependent medium (two materials
// with the same model, one on a single element) against the same estimator in vacuum:
// the projection commutes with the constant ε(iω), so D_h − D = ε (E − P E) and
// η² = |ε|² / W η²_vacuum on every element, with the cross-coupled complex ε (not Re ε²)
// and the positive stored-energy weight W.
void CheckHomogeneousEstimator(int dim, mfem::Element::Type type)
{
  test::IntegrationSettingsGuard guard;
  constexpr int order = 2;
  fem::DefaultIntegrationOrder::p_trial = order;
  auto mesh = MakeBoxMesh(dim, 2, type);
  config::PeriodicBoundaryData periodic;
  mfem::ND_FECollection nd_fec(order, dim);
  mfem::RT_FECollection rt_fec(order - 1, dim);
  FiniteElementSpaceHierarchy rt_fespaces(
      std::make_unique<FiniteElementSpace>(mesh, &rt_fec));
  FiniteElementSpace nd_fespace(mesh, &nd_fec);
  constexpr double tol = 1.0e-14;
  constexpr int max_it = 1000;

  ComplexVector E(nd_fespace.GetTrueVSize());
  E.UseDevice(true);
  linalg::SetRandom(Mpi::World(), E);

  config::MaterialData vacuum;
  vacuum.attributes = {1, 2};
  MaterialOperator vacuum_op({vacuum}, periodic, ProblemType::DRIVEN, mesh);
  GradFluxErrorEstimator<ComplexVector> vacuum_est(vacuum_op, nd_fespace, rt_fespaces, tol,
                                                   max_it, 0, false);
  ErrorIndicator ref;
  vacuum_est.AddErrorIndicator(E, 1.0, ref);

  auto Check = [&](const ErrorIndicator &ind, double ratio)
  {
    const auto *h = ind.Local().HostRead();
    const auto *h_ref = ref.Local().HostRead();
    REQUIRE(ind.Local().Size() == ref.Local().Size());
    for (int i = 0; i < ind.Local().Size(); i++)
    {
      CHECK_THAT(h[i], WithinRel(std::sqrt(ratio) * h_ref[i], 1.0e-8));
    }
  };
  for (const auto &model : Models())
  {
    CAPTURE(model.name);
    MaterialOperator mat_op(
        {MakeMaterial({1}, model.terms), MakeMaterial({2}, model.terms)}, periodic,
        ProblemType::DRIVEN, mesh);
    GradFluxErrorEstimator<ComplexVector> estimator(mat_op, nd_fespace, rt_fespaces, tol,
                                                    max_it, 0, false);
    // The frequency is required, and each change refreshes the retained contexts and
    // projection weights.
    {
      ErrorIndicator ind;
      CHECK_THROWS(estimator.AddErrorIndicator(E, 1.0, ind));
    }
    for (double omega : {0.6, 2.5, 0.6})
    {
      CAPTURE(omega);
      const double eps_abs2 = std::norm(eps_inf + model.chi(omega));
      const double W = eps_inf + model.w(omega);
      ErrorIndicator ind;
      estimator.AddErrorIndicator(E, 1.0, omega, ind);
      Check(ind, eps_abs2 / W);
    }
  }

  // A general PoleResidue model with W(ω) ≤ 0 keeps a positive indicator with ε∞ in place
  // of W (and warns). A relaxation with a negative residue has W = ε∞ + a b / (ω² + a²).
  MaterialOperator general_op(
      {MakeMaterial({1, 2}, {config::PermittivityPoleData{-1.0, -10.0}})}, periodic,
      ProblemType::DRIVEN, mesh);
  GradFluxErrorEstimator<ComplexVector> general(general_op, nd_fespace, rt_fespaces, tol,
                                                max_it, 0, false);
  const double omega = 0.5;
  REQUIRE(general_op.GetFrequencyDependentStoredEnergyPermittivity(0, omega) < 0.0);
  ErrorIndicator ind;
  general.AddErrorIndicator(E, 1.0, omega, ind);
  const std::complex<double> eps =
      eps_inf - 10.0 / (std::complex<double>(0.0, omega) + 1.0);
  Check(ind, std::norm(eps) / eps_inf);
}

}  // namespace

TEST_CASE("Gradient flux estimator of frequency-dependent materials",
          "[frequencydependent][errorestimator][Serial][Parallel][GPU]")
{
  SECTION("3D")
  {
    CheckHomogeneousEstimator(3, mfem::Element::TETRAHEDRON);
  }
  SECTION("2D")
  {
    CheckHomogeneousEstimator(2, mfem::Element::QUADRILATERAL);
  }
}

TEST_CASE("Flux error estimator with a nonpositive energy normalization",
          "[frequencydependent][errorestimator][Serial][Parallel][GPU]")
{
  // A general model can make the total stored energy nonpositive. The indicator then uses
  // the energy norm of the recovered fluxes, 1/2 Σ (D*ᴴ W⁻¹ D* + H*ᴴ μ H*), not 1, NaN, or
  // a rejection: it equals the indicator normalized by that energy passed explicitly.
  test::IntegrationSettingsGuard guard;
  constexpr int order = 1;
  fem::DefaultIntegrationOrder::p_trial = order;
  auto mesh = MakeBoxMesh(3, 2, mfem::Element::TETRAHEDRON);
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op(
      {MakeMaterial({1}, Models()[0].terms), MakeMaterial({2}, Models()[2].terms)},
      periodic, ProblemType::DRIVEN, mesh);
  mfem::ND_FECollection nd_fec(order, 3);
  mfem::RT_FECollection rt_fec(order - 1, 3);
  FiniteElementSpaceHierarchy nd_fespaces(
      std::make_unique<FiniteElementSpace>(mesh, &nd_fec));
  FiniteElementSpaceHierarchy rt_fespaces(
      std::make_unique<FiniteElementSpace>(mesh, &rt_fec));
  TimeDependentFluxErrorEstimator<ComplexVector> estimator(mat_op, nd_fespaces, rt_fespaces,
                                                           1.0e-14, 1000, 0, false);
  ComplexVector E(nd_fespaces.GetFinestFESpace().GetTrueVSize()),
      B(rt_fespaces.GetFinestFESpace().GetTrueVSize());
  E.UseDevice(true);
  B.UseDevice(true);
  linalg::SetRandom(Mpi::World(), E);
  linalg::SetRandom(Mpi::World(), B);
  const double omega = 1.7;
  ErrorIndicator fallback, unit;
  estimator.AddErrorIndicator(E, B, -1.0, omega, fallback);
  estimator.AddErrorIndicator(E, B, 0.5, omega, unit);  // Scale 1: η² unnormalized
  const auto *h = fallback.Local().HostRead();
  const auto *h_unit = unit.Local().HostRead();
  // η_fallback² = η² · 0.5 / E_r, so Σ η_fallback² / Σ η² = 0.5 / E_r is constant, and
  // positive and finite.
  double ratio = -1.0;
  for (int i = 0; i < fallback.Local().Size(); i++)
  {
    REQUIRE(std::isfinite(h[i]));
    REQUIRE(h[i] > 0.0);
    const double r = (h[i] * h[i]) / (h_unit[i] * h_unit[i]);
    if (ratio < 0.0)
    {
      ratio = r;
    }
    CHECK_THAT(r, WithinRel(ratio, 1.0e-10));
  }
}

namespace
{

// The 2D magnetic recovery test setup: the curl space
// with the INTEGRAL map (B = B̂ / |J|) and its H1 recovery space.
struct Curl2D
{
  Mesh mesh;
  mfem::L2_FECollection l2_fec;
  mfem::H1_FECollection h1_fec;
  FiniteElementSpace l2_fespace;
  FiniteElementSpaceHierarchy h1_fespaces;
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op;

  static Mesh Make(mfem::Element::Type type, bool distort)
  {
    auto serial = std::make_unique<mfem::Mesh>(
        mfem::Mesh::MakeCartesian2D(3, 3, type, false, 1.0, 0.7));
    // Unequal element sizes, and for quadrilaterals non-parallelograms (straight sides).
    serial->Transform(
        [distort](const mfem::Vector &x, mfem::Vector &y)
        {
          y = x;
          y(0) = x(0) * x(0) * 0.6 + 0.4 * x(0);
          if (distort)
          {
            y(1) += 0.12 * std::sin(std::numbers::pi * x(0)) *
                    std::sin(std::numbers::pi * x(1) / 0.7);
          }
        });
    for (int i = 0; i < serial->GetNE(); i++)
    {
      serial->SetAttribute(i, 1);
    }
    serial->SetAttributes();
    return Mesh(std::make_unique<mfem::ParMesh>(Mpi::World(), *serial));
  }

  Curl2D(mfem::Element::Type type, bool distort, int order)
    : mesh(Make(type, distort)),
      l2_fec(order - 1, 2, mfem::BasisType::GaussLegendre, mfem::FiniteElement::INTEGRAL),
      h1_fec(order, 2), l2_fespace(mesh, &l2_fec),
      h1_fespaces(std::make_unique<FiniteElementSpace>(mesh, &h1_fec)),
      mat_op(
          [&]()
          {
            config::MaterialData m;
            m.attributes = {1};
            m.mu_r.s = {2.0, 2.0, 2.0};
            return std::vector<config::MaterialData>{m};
          }(),
          periodic, ProblemType::EIGENMODE, mesh)
  {
  }
};

}  // namespace

TEST_CASE("2D magnetic flux recovery uses the physical B of the INTEGRAL map",
          "[errorestimator][Serial][Parallel][GPU]")
{
  test::IntegrationSettingsGuard guard;
  constexpr int order = 2;
  fem::DefaultIntegrationOrder::p_trial = order;
  constexpr double mu = 2.0;

  SECTION("Constant B has zero indicator on unequal and distorted cells")
  {
    const auto type = GENERATE(mfem::Element::TRIANGLE, mfem::Element::QUADRILATERAL);
    const bool distort = GENERATE(false, true);
    CAPTURE(type, distort);
    Curl2D c(type, distort, order);
    GridFunction B(c.l2_fespace);
    mfem::ConstantCoefficient one(1.0);
    B.Real().ProjectCoefficient(one);
    CurlFluxErrorEstimator<Vector> estimator(c.mat_op, c.l2_fespace, c.h1_fespaces, 1.0e-14,
                                             1000, 0, false);
    ErrorIndicator ind;
    estimator.AddErrorIndicator(B.Real(), 0.5, ind);  // Scale 1
    const auto *h = ind.Local().HostRead();
    for (int i = 0; i < ind.Local().Size(); i++)
    {
      CHECK(h[i] < 1.0e-10);
    }
  }

  SECTION("Nonzero indicator matches physical quadrature")
  {
    // Discontinuous B on unequal triangles; H* is the H1 projection of μ⁻¹ B computed with
    // MFEM physical shape functions and η²_K = ∫_K |μ^(-1/2) B − μ^(1/2) H*|².
    Curl2D c(mfem::Element::TRIANGLE, false, order);
    auto &pmesh = c.mesh.Get();
    GridFunction B(c.l2_fespace);
    mfem::FunctionCoefficient f([](const mfem::Vector &x)
                                { return 1.0 + x(0) * x(1) + std::sin(3.0 * x(0)); });
    B.Real().ProjectCoefficient(f);
    auto &h1 = c.h1_fespaces.GetFinestFESpace();
    mfem::ParBilinearForm m(&h1.Get());
    m.AddDomainIntegrator(new mfem::MassIntegrator);
    m.Assemble();
    m.Finalize();
    std::unique_ptr<mfem::HypreParMatrix> M(m.ParallelAssemble());
    mfem::GridFunctionCoefficient B_coeff(&B.Real());
    mfem::ProductCoefficient rhs_coeff(1.0 / mu, B_coeff);
    mfem::ParLinearForm b(&h1.Get());
    const auto &lf_ir = mfem::IntRules.Get(mfem::Geometry::TRIANGLE, 10);
    b.AddDomainIntegrator(new mfem::DomainLFIntegrator(rhs_coeff, &lf_ir));
    b.Assemble();
    std::unique_ptr<mfem::HypreParVector> rhs(b.ParallelAssemble());
    mfem::HypreParVector Ht(*rhs);
    Ht = 0.0;
    mfem::CGSolver cg(pmesh.GetComm());
    cg.SetRelTol(1.0e-15);
    cg.SetAbsTol(0.0);
    cg.SetMaxIter(2000);
    cg.SetOperator(*M);
    cg.Mult(*rhs, Ht);
    mfem::ParGridFunction H(&h1.Get());
    H.SetFromTrueDofs(Ht);

    CurlFluxErrorEstimator<Vector> estimator(c.mat_op, c.l2_fespace, c.h1_fespaces, 1.0e-15,
                                             2000, 0, false);
    ErrorIndicator ind;
    estimator.AddErrorIndicator(B.Real(), 0.5, ind);  // Scale 1: η_K
    const auto *h = ind.Local().HostRead();
    double max_ref = 0.0;
    for (int e = 0; e < pmesh.GetNE(); e++)
    {
      auto *T = pmesh.GetElementTransformation(e);
      const auto &ir = mfem::IntRules.Get(pmesh.GetElementGeometry(e), 12);
      double eta2 = 0.0;
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
        const auto &ip = ir.IntPoint(q);
        T->SetIntPoint(&ip);
        const double d =
            B.Real().GetValue(*T, ip) / std::sqrt(mu) - std::sqrt(mu) * H.GetValue(*T, ip);
        eta2 += ip.weight * T->Weight() * d * d;
      }
      max_ref = std::max(max_ref, std::sqrt(eta2));
      CHECK_THAT(h[e], WithinRel(std::sqrt(eta2), 1.0e-8));
    }
    Mpi::GlobalMax(1, &max_ref, pmesh.GetComm());
    CHECK(max_ref > 1.0e-4);  // Not vacuous
  }
}

}  // namespace palace
