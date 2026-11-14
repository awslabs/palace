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
#include <limits>
#include <memory>
#include <numbers>
#include <optional>
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
#include "models/domainpostoperator.hpp"
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
#include "utils/permittivitymodel.hpp"
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
    // Typed terms, and the RAW zero pole (a current, which stores no energy alone).
    CHECK(mat_op.HasStoredElectricEnergy());
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

TEST_CASE("Stored energy of a frequency-dependent material requires a finite frequency",
          "[frequencydependent][Serial][Parallel]")
{
  // At ω = +∞ the Debye, Drude, and zero-pole energies are finite (zero), so the frequency
  // itself is checked, not only the result: the same error for every model, and for a
  // Debye alone.
  auto mesh = MakeBoxMesh(3, 2, mfem::Element::HEXAHEDRON);
  config::PeriodicBoundaryData periodic;
  const double inf = std::numeric_limits<double>::infinity();
  std::vector<std::vector<config::PermittivityTermData>> models = {
      {config::DebyeData{3.0, 0.7}}};
  for (const auto &model : Models())
  {
    models.push_back(model.terms);
  }
  for (std::size_t k = 0; k < models.size(); k++)
  {
    CAPTURE(k);
    MaterialOperator mat_op({MakeMaterial({1, 2}, models[k])}, periodic,
                            ProblemType::DRIVEN, mesh);
    REQUIRE(mat_op.HasStoredElectricEnergy(0));
    CHECK_THROWS_WITH(mat_op.GetFrequencyDependentStoredEnergySusceptibility(0, inf),
                      ContainsSubstring("finite positive frequency"));
    CHECK_THROWS_WITH(mat_op.GetFrequencyDependentStoredEnergyPermittivity(0, inf),
                      ContainsSubstring("finite positive frequency"));
  }
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
  // Lossless propagating and evanescent orders keep their established roots, exactly: a
  // real part of round-off size on an evanescent order would count it as propagating.
  CHECK(F::OutgoingPropagationConstant({4.0, 0.0}) == std::complex<double>(2.0, 0.0));
  CHECK(F::OutgoingPropagationConstant({0.0, 0.0}) == std::complex<double>(0.0, 0.0));
  for (const double zero : {0.0, -0.0})
  {
    const auto evanescent = F::OutgoingPropagationConstant({-4.0, zero});
    CHECK(evanescent.real() == 0.0);
    CHECK(evanescent.imag() == -2.0);
    // No power in either polarization of a lossless evanescent order.
    CHECK(F::DtNEigenvalue(true, evanescent, 2.0, 1.5).real() == 0.0);
    CHECK(F::DtNEigenvalue(false, evanescent, 2.0, 1.5).real() == 0.0);
    const auto far = F::OutgoingPropagationConstant({-3.0e7, zero});
    CHECK(far.real() == 0.0);
    CHECK(far.imag() == -std::sqrt(3.0e7));
  }
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

TEST_CASE_METHOD(test::SharedTempDir,
                 "Floquet port orders capped by the largest element stretch",
                 "[floquetport][Serial][Parallel]")
{
  // The diffraction orders of a Floquet port are capped at the Nyquist limit of the port
  // elements, |b| h_max < p π, with h_max the largest singular value of the 3 x 2 boundary
  // element Jacobians. The port faces of a periodic box of 3 x 4 elements over a x a are
  // affine rectangles with sides a / 3 and a / 4, so h_max = a / 3 and |b| = 2π / a: for
  // p = 3, the cap is ⌊1.5 · 3⌋ = 4 below MaxOrder = 5. The shorter side a / 4 gives no
  // cap, the diagonal (Frobenius norm) gives 3. In a lossy medium every order carries
  // power, so the port reports S-parameters for all (2 M + 1)² orders in both
  // polarizations.
  using json = nlohmann::json;
  MPI_Comm comm = Mpi::World();
  IoData iodata{Units(1.0e-3, 1.0e-3)};
  iodata.model.L0 = 1.0e-3;
  iodata.model.Lc = 1.0;
  iodata.problem.type = ProblemType::DRIVEN;
  constexpr int p = 3, max_order = 5, nx = 3, ny = 4;
  iodata.solver.order = p;
  iodata.solver.driven.sample_f = {100.0};  // GHz
  json permittivity = {{"HighFrequency", 2.0},
                       {"Terms", json::array({{{"Type", "Drude"},
                                               {"PlasmaFrequency", 100.0},
                                               {"CollisionFrequency", 50.0}}})}};
  iodata.domains = config::DomainData(json{
      {"Materials",
       json::array({{{"Attributes", json::array({1})}, {"Permittivity", permittivity}}})}});
  iodata.boundaries = config::BoundaryData(json{
      {"Periodic",
       {{"BoundaryPairs", json::array({{{"DonorAttributes", json::array({5})},
                                        {"ReceiverAttributes", json::array({3})}},
                                       {{"DonorAttributes", json::array({2})},
                                        {"ReceiverAttributes", json::array({4})}}})}}},
      {"FloquetPort",
       json::array(
           {{{"Index", 1},
             {"Attributes", json::array({1})},
             {"Excitation", 1},
             {"MaxOrder", max_order}},
            {{"Index", 2}, {"Attributes", json::array({6})}, {"MaxOrder", max_order}}})}});
  iodata.CheckConfiguration();
  constexpr double a = 0.5, L = 0.25;
  iodata.model.mesh = (temp_dir / "floquet_cap_box.mesh").string();
  if (Mpi::Root(comm))
  {
    auto box = mfem::Mesh::MakeCartesian3D(nx, ny, 2, mfem::Element::HEXAHEDRON, a, a, L);
    std::ofstream ofs(iodata.model.mesh);
    box.Print(ofs);
  }
  Mpi::Barrier(comm);
  auto smesh = mesh::Load(iodata, comm);
  iodata.NondimensionalizeInputs(smesh);
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(mesh::Partition(iodata, std::move(smesh), comm)));
  SpaceOperator space_op(iodata, mesh);

  auto cap = [&](double h)
  {
    const double b = 2.0 * std::numbers::pi / a;
    return std::min(max_order,
                    static_cast<int>(std::floor(p * std::numbers::pi / (b * h))));
  };
  const int M = cap(a / nx);
  REQUIRE(M == 4);
  REQUIRE(cap(a / ny) == max_order);
  REQUIRE(cap(std::hypot(a / nx, a / ny)) == 3);

  auto &floquet_op = space_op.GetFloquetPortOp();
  floquet_op.Initialize(iodata.solver.driven.sample_f.front());
  GridFunction E(space_op.GetNDSpace(), true);
  E.Real() = 0.0;
  E.Imag() = 0.0;
  for (int idx : {1, 2})
  {
    const auto S = floquet_op.GetPort(idx).GetAllSParameters(E);
    CAPTURE(idx, S.size());
    CHECK(S.size() == static_cast<std::size_t>(2 * (2 * M + 1) * (2 * M + 1)));
    CHECK(S.contains(std::make_tuple(M, -M, false)));
    CHECK_FALSE(S.contains(std::make_tuple(M + 1, 0, true)));
  }
}

// Whether the frequency-domain postprocessing of a SurfacePostOperator can be called with
// and without the frequency of the fields.
template <typename Op>
concept CallableWithFrequency = requires(const Op &op, const GridFunction *E) {
  op.GetSurfaceFlux(1, E, nullptr, 1.0);
} && requires(const Op &op, const GridFunction &E) {
  op.GetInterfaceElectricFieldEnergy(1, E, 1.0);
};
template <typename Op>
concept CallableWithoutFrequency = requires(const Op &op, const GridFunction *E) {
  op.GetSurfaceFlux(1, E, nullptr);
} || requires(const Op &op, const GridFunction &E) {
  op.GetInterfaceElectricFieldEnergy(1, E);
};

// Whether the gradient flux error indicator of fields of type V can be computed without
// their frequency.
template <typename Est, typename V = ComplexVector>
concept IndicatorWithoutFrequency =
    requires(const Est &est, const V &E, ErrorIndicator &ind) {
      est.AddErrorIndicator(E, 1.0, ind);
    };

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
  // The frequency of the fields is a required argument, and positive for a
  // frequency-dependent material.
  static_assert(CallableWithFrequency<SurfacePostOperator> &&
                !CallableWithoutFrequency<SurfacePostOperator>);
  CHECK_THROWS(surf_post_op.GetSurfaceFlux(1, &E, nullptr, 0.0));
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

    // The MFEM coefficient of the legacy path, with the same substrate table (and the
    // neighbor field values of interfaces on process boundaries).
    E.Real().ExchangeFaceNbrData();
    E.Imag().ExchangeFaceNbrData();
    const auto epsilon_s = MetalSubstratePermittivity(mat_op, omega);
    InterfaceDielectricCoefficient<InterfaceDielectric::MS> coeff(
        E, mat_op, data.t, data.epsilon_r, nullptr, &epsilon_s);
    const auto &pmesh = mesh.Get();
    auto marker = mesh::AttrToMarker(pmesh.bdr_attributes.Max(), data.attributes);
    mfem::LinearForm lf(&h1_fespace.Get());
    lf.AddBoundaryIntegrator(new BoundaryLFIntegrator(coeff), marker);
    lf.UseFastAssembly(false);
    lf.UseDevice(false);
    lf.Assemble();
    double legacy = linalg::LocalSum(lf);
    Mpi::GlobalSum(1, &legacy, Mpi::World());
    CHECK_THAT(legacy, WithinRel(expected, 1.0e-10));
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
    const double *v = val.HostRead();  // Eval writes on the device
    for (int i = 0; i < val.Size(); i++)
    {
      CHECK_THAT(v[i], WithinRel(expected, 1.0e-10));
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
  // Complex fields require their frequency. The constant reference does not depend on it,
  // so it is evaluated at one fixture frequency and used at every frequency below.
  constexpr double omega_ref = 1.0;
  ErrorIndicator ref;
  vacuum_est.AddErrorIndicator(E, 1.0, omega_ref, ref);

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
    // The frequency is required (at compile time), and each change refreshes the retained
    // contexts and projection weights.
    static_assert(!IndicatorWithoutFrequency<GradFluxErrorEstimator<ComplexVector>> &&
                  IndicatorWithoutFrequency<GradFluxErrorEstimator<Vector>, Vector>);
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

  // A material without a stored electric energy, a RAW −0.5/(s + 1) compensated by a
  // typed Debye 1/(1 + s) (PoleResidue terms specify none), weights its indicator
  // by ε∞, a numerical weight, and not by a termwise W (here Re ε). Next to a constant
  // lossy material with the same complex ε(iω), the medium is homogeneous in ε, and every
  // element has η² = c |ε|² / weight η²_vacuum with one c > 0: the normalization is the
  // numerical recovered-flux norm, as no total energy is passed.
  {
    const double omega = 0.5;
    const std::complex<double> s(0.0, omega);
    const std::complex<double> eps = eps_inf + 1.0 / (1.0 + s) - 0.5 / (s + 1.0);
    config::MaterialData lossy;
    lossy.attributes = {1};
    lossy.epsilon_r.s = {eps.real(), eps.real(), eps.real()};
    const double tandelta = -eps.imag() / eps.real();
    lossy.tandelta.s = {tandelta, tandelta, tandelta};
    MaterialOperator compensated_op(
        {lossy, MakeMaterial({2}, {config::DebyeData{1.0, 1.0},
                                   config::PermittivityPoleData{-1.0, -0.5}})},
        periodic, ProblemType::DRIVEN, mesh);
    REQUIRE_FALSE(compensated_op.HasStoredElectricEnergy(1));
    REQUIRE_FALSE(compensated_op.HasStoredElectricEnergy());
    GradFluxErrorEstimator<ComplexVector> compensated(compensated_op, nd_fespace,
                                                      rt_fespaces, tol, max_it, 0, false);
    ErrorIndicator ind;
    compensated.AddErrorIndicator(E, std::nullopt, omega, ind);
    // A total energy is not accepted without one.
    {
      ErrorIndicator unused;
      CHECK_THROWS(compensated.AddErrorIndicator(E, 1.0, omega, unused));
    }
    const auto *h = ind.Local().HostRead();
    const auto *h_ref = ref.Local().HostRead();
    REQUIRE(ind.Local().Size() == mesh.GetNE());
    double c_min = mfem::infinity(), c_max = 0.0;
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      REQUIRE(std::isfinite(h[e]));
      REQUIRE(h[e] > 0.0);
      const double weight = (mesh.Get().GetAttribute(e) == 1) ? eps.real() : eps_inf;
      const double c = (h[e] * h[e]) / (h_ref[e] * h_ref[e]) * weight / std::norm(eps);
      c_min = std::min(c_min, c);
      c_max = std::max(c_max, c);
    }
    Mpi::GlobalMin(1, &c_min, Mpi::World());
    Mpi::GlobalMax(1, &c_max, Mpi::World());
    CHECK_THAT(c_max, WithinRel(c_min, 1.0e-8));
    // An available material does not accept a missing total energy.
    MaterialOperator lorentz_op({MakeMaterial({1, 2}, Models()[1].terms)}, periodic,
                                ProblemType::DRIVEN, mesh);
    GradFluxErrorEstimator<ComplexVector> lorentz_est(lorentz_op, nd_fespace, rt_fespaces,
                                                      tol, max_it, 0, false);
    ErrorIndicator unused;
    CHECK_THROWS(lorentz_est.AddErrorIndicator(E, std::nullopt, omega, unused));
  }

  // A frequency-dependent material (attribute 2) next to a constant lossy one (attribute 1)
  // with the same complex ε at ω (Re ε, tan δ = −Im ε / Re ε): the medium is homogeneous
  // in ε only if the constant material also projects its loss, and each element is
  // weighted by its own W, the stored energy of the model or Re ε of the constant material.
  const auto lorentz = Models()[1];
  for (double w : {0.6, 2.5})
  {
    CAPTURE(w);
    const std::complex<double> eps_w = eps_inf + lorentz.chi(w);
    REQUIRE(eps_w.real() > 0.0);
    config::MaterialData lossy;
    lossy.attributes = {1};
    lossy.epsilon_r.s = {eps_w.real(), eps_w.real(), eps_w.real()};
    const double tandelta = -eps_w.imag() / eps_w.real();
    lossy.tandelta.s = {tandelta, tandelta, tandelta};
    MaterialOperator mixed_op({lossy, MakeMaterial({2}, lorentz.terms)}, periodic,
                              ProblemType::DRIVEN, mesh);
    GradFluxErrorEstimator<ComplexVector> mixed(mixed_op, nd_fespace, rt_fespaces, tol,
                                                max_it, 0, false);
    ErrorIndicator mixed_ind;
    mixed.AddErrorIndicator(E, 1.0, w, mixed_ind);
    const double W = eps_inf + lorentz.w(w);
    const auto *h = mixed_ind.Local().HostRead();
    const auto *h_ref = ref.Local().HostRead();
    REQUIRE(mixed_ind.Local().Size() == mesh.GetNE());
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      const double weight = (mesh.Get().GetAttribute(e) == 1) ? eps_w.real() : W;
      CHECK_THAT(h[e], WithinRel(std::sqrt(std::norm(eps_w) / weight) * h_ref[e], 1.0e-8));
    }
  }
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

TEST_CASE("Gradient flux estimator of a layered medium with continuous normal D",
          "[frequencydependent][errorestimator][Serial][Parallel][GPU]")
{
  // Two layers: a Drude medium (Re ε < 0 at ω) below a constant lossy dielectric, with the
  // piecewise constant field E = (c, 0, 1/ε) (2D: (c, 1/ε)), whose tangential part is
  // continuous. Its complete displacement D = ε E has a continuous normal component and is
  // in the recovery space, so every indicator vanishes: the complex constitutive flux with
  // the loss of both materials and its real-imaginary coupling. The same field with
  // E_n = 1 on both sides (a jump in D_n) has nonzero indicators.
  const int dim = GENERATE(3, 2);
  CAPTURE(dim);
  test::IntegrationSettingsGuard guard;
  constexpr int order = 1;
  fem::DefaultIntegrationOrder::p_trial = order;
  auto serial = std::make_unique<mfem::Mesh>(
      (dim == 3) ? mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::HEXAHEDRON)
                 : mfem::Mesh::MakeCartesian2D(2, 2, mfem::Element::QUADRILATERAL, true));
  for (int i = 0; i < serial->GetNE(); i++)
  {
    mfem::Vector c;
    serial->GetElementCenter(i, c);
    serial->SetAttribute(i, (c(dim - 1) < 0.5) ? 1 : 2);
  }
  serial->SetAttributes();
  Mesh mesh(std::make_unique<mfem::ParMesh>(Mpi::World(), *serial));
  const double omega = 1.0;
  const std::complex<double> eps1 = eps_inf + DrudeChi(omega), eps2(3.0, -0.6);
  REQUIRE(eps1.real() < 0.0);
  config::MaterialData lossy;
  lossy.attributes = {2};
  lossy.epsilon_r.s = {3.0, 3.0, 3.0};
  lossy.tandelta.s = {0.2, 0.2, 0.2};
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({MakeMaterial({1}, Models()[0].terms), lossy}, periodic,
                          ProblemType::DRIVEN, mesh);
  mfem::ND_FECollection nd_fec(order, dim);
  mfem::RT_FECollection rt_fec(order - 1, dim);
  FiniteElementSpace nd_fespace(mesh, &nd_fec);
  FiniteElementSpaceHierarchy rt_fespaces(
      std::make_unique<FiniteElementSpace>(mesh, &rt_fec));
  GradFluxErrorEstimator<ComplexVector> estimator(mat_op, nd_fespace, rt_fespaces, 1.0e-14,
                                                  1000, 0, false);

  auto Estimate = [&](bool continuous)
  {
    // Real and imaginary parts of E on each layer.
    std::array<mfem::Vector, 4> values;
    for (int k = 0; k < 2; k++)
    {
      const std::complex<double> En = continuous ? 1.0 / ((k == 0) ? eps1 : eps2) : 1.0;
      values[2 * k].SetSize(dim);
      values[2 * k + 1].SetSize(dim);
      values[2 * k] = 0.0;
      values[2 * k + 1] = 0.0;
      values[2 * k](0) = 0.3;
      values[2 * k](dim - 1) = En.real();
      values[2 * k + 1](dim - 1) = En.imag();
    }
    mfem::VectorConstantCoefficient r1(values[0]), i1(values[1]), r2(values[2]),
        i2(values[3]);
    mfem::Array<int> attrs({1, 2});
    mfem::Array<mfem::VectorCoefficient *> fr({&r1, &r2}), fi({&i1, &i2});
    mfem::PWVectorCoefficient cr(dim, attrs, fr), ci(dim, attrs, fi);
    GridFunction E_gf(nd_fespace, true);
    E_gf.Real().ProjectCoefficient(cr);
    E_gf.Imag().ProjectCoefficient(ci);
    ComplexVector E(nd_fespace.GetTrueVSize());
    E_gf.Real().GetTrueDofs(E.Real());
    E_gf.Imag().GetTrueDofs(E.Imag());
    ErrorIndicator ind;
    estimator.AddErrorIndicator(E, 1.0, omega, ind);
    double max = ind.Local().Size() ? ind.Local().Normlinf() : 0.0;
    Mpi::GlobalMax(1, &max, Mpi::World());
    return max;
  };
  const double jump = Estimate(false);
  REQUIRE(jump > 0.1);
  CHECK(Estimate(true) < 1.0e-8 * jump);
}

TEST_CASE("Flux error estimator without a total energy normalization",
          "[frequencydependent][errorestimator][Serial][Parallel][GPU]")
{
  // Without a stored electric energy (a RAW relaxation, compensated by a typed Debye
  // term) there is no total energy. The indicator then uses the energy norm of
  // the recovered fluxes, 1/2 Σ (D*ᴴ W⁻¹ D* + H*ᴴ μ H*), a numerical norm, not 1, NaN, or a
  // rejection: positive, finite, and invariant to the scale of the fields, unlike a fixed
  // normalization.
  test::IntegrationSettingsGuard guard;
  constexpr int order = 1;
  fem::DefaultIntegrationOrder::p_trial = order;
  auto mesh = MakeBoxMesh(3, 2, mfem::Element::TETRAHEDRON);
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({MakeMaterial({1}, {config::DebyeData{1.0, 1.0},
                                              config::PermittivityPoleData{-1.0, -0.5}}),
                           MakeMaterial({2}, Models()[2].terms)},
                          periodic, ProblemType::DRIVEN, mesh);
  REQUIRE_FALSE(mat_op.HasStoredElectricEnergy());
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
  ErrorIndicator fallback, scaled;
  estimator.AddErrorIndicator(E, B, std::nullopt, omega, fallback);
  {
    ComplexVector E2(E), B2(B);
    E2 *= 2.0;
    B2 *= 2.0;
    estimator.AddErrorIndicator(E2, B2, std::nullopt, omega, scaled);
  }
  const auto *h = fallback.Local().HostRead();
  const auto *h_scaled = scaled.Local().HostRead();
  for (int i = 0; i < fallback.Local().Size(); i++)
  {
    REQUIRE(std::isfinite(h[i]));
    REQUIRE(h[i] > 0.0);
    CHECK_THAT(h_scaled[i], WithinRel(h[i], 1.0e-8));
  }
  {
    ErrorIndicator unused;
    CHECK_THROWS(estimator.AddErrorIndicator(E, B, -1.0, omega, unused));
  }

  // Without a frequency-dependent material, the normalization is that of the constant
  // estimators also for a nonpositive energy (a zero field): scale 1, with no fallback to
  // the recovered fluxes.
  MaterialOperator constant_op({MakeMaterial({1}, {}), MakeMaterial({2}, {})}, periodic,
                               ProblemType::DRIVEN, mesh);
  TimeDependentFluxErrorEstimator<ComplexVector> constant_estimator(
      constant_op, nd_fespaces, rt_fespaces, 1.0e-14, 1000, 0, false);
  // The established scale for Et ≤ 0 is 1, that of Et = 0.5.
  ErrorIndicator established;
  constant_estimator.AddErrorIndicator(E, B, 0.5, omega, established);
  for (double Et : {0.0, -1.0})
  {
    CAPTURE(Et);
    ErrorIndicator with_omega;
    constant_estimator.AddErrorIndicator(E, B, Et, omega, with_omega);
    const auto *a = with_omega.Local().HostRead();
    const auto *b = established.Local().HostRead();
    for (int i = 0; i < with_omega.Local().Size(); i++)
    {
      // Two evaluations of the same indicators with the same scale, equal up to roundoff:
      // device reductions (CUDA/MAGMA) are not bitwise repeatable, and differ by 1-2 ULPs.
      CAPTURE(i, a[i], b[i]);
      CHECK_THAT(a[i], Catch::Matchers::WithinULP(b[i], 4));
    }
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

  static Mesh Make(mfem::Element::Type type, bool distort, bool curved)
  {
    auto serial = std::make_unique<mfem::Mesh>(
        mfem::Mesh::MakeCartesian2D(3, 3, type, false, 1.0, 0.7));
    if (curved)
    {
      // Quadratic geometry: the transform below then bends the element edges.
      serial->SetCurvature(2);
    }
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

  Curl2D(mfem::Element::Type type, bool distort, int order, bool curved = false)
    : mesh(Make(type, distort, curved)),
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

  SECTION("Constant B has zero indicator on curved cells")
  {
    // On quadratic triangles |J| is quadratic, so B̂ = |J| B for a constant B is exact in
    // the second-order curl space, and so is its recovery, independent of quadrature.
    Curl2D c(mfem::Element::TRIANGLE, true, 3, true);
    REQUIRE(c.mesh.Get().GetNodes());
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

namespace
{

// The terms as RAW PoleResidue terms with the same susceptibility.
std::vector<config::PermittivityTermData>
RawTerms(const std::vector<config::PermittivityTermData> &terms)
{
  std::vector<config::PermittivityTermData> raw;
  for (const auto &t : permittivity::PoleTerms(terms))
  {
    raw.push_back(t);
  }
  return raw;
}

// The energies of a field over the entire domain and the postprocessing domains, for
// comparison between materials.
struct EnergyConsumers
{
  Mesh &mesh;
  mfem::ND_FECollection nd_fec;
  mfem::RT_FECollection rt_fec;
  FiniteElementSpace nd_fespace, rt_fespace;
  config::DomainPostData postpro;
  GridFunction E;
  ComplexVector e;

  explicit EnergyConsumers(Mesh &mesh)
    : mesh(mesh), nd_fec(1, 3), rt_fec(0, 3), nd_fespace(mesh, &nd_fec),
      rt_fespace(mesh, &rt_fec), E(nd_fespace, true), e(nd_fespace.GetTrueVSize())
  {
    postpro.energy[1].attributes = {2};
    postpro.energy[2].attributes = {1, 2};
    e.UseDevice(true);
    linalg::SetRandom(Mpi::World(), e);
    E.Real().SetFromTrueDofs(e.Real());
    E.Imag().SetFromTrueDofs(e.Imag());
  }
};

}  // namespace

TEST_CASE("Single sections as typed or RAW terms", "[frequencydependent][Serial][Parallel]")
{
  // A single Debye, Drude, or under- or overdamped Lorentz section has the same response
  // as typed terms and as RAW PoleResidue terms, also with the RAW terms reversed and
  // split, and so the same field energy at Re ε in the material and the domains. Only the
  // typed terms specify a stored energy: the RAW terms specify a susceptibility, whose
  // stored energy is not evaluated (an error, not a termwise value). A Drude and an
  // overdamped Lorentz are RAW fits with a negative residue.
  test::IntegrationSettingsGuard guard;
  fem::DefaultIntegrationOrder::p_trial = 1;
  auto mesh = MakeBoxMesh(3, 2, mfem::Element::HEXAHEDRON);
  config::PeriodicBoundaryData periodic;
  EnergyConsumers c(mesh);
  using ElectricEnergy = DomainPostOperator::ElectricEnergy;
  const std::vector<config::PermittivityTermData> sections = {
      config::DebyeData{3.0, 0.7}, config::DrudeData{3.0, 0.5},
      config::LorentzData{1.5, 1.0, 0.2}, config::LorentzData{1.5, 1.0, 5.0}};
  for (std::size_t k = 0; k < sections.size(); k++)
  {
    CAPTURE(k);
    const std::vector<config::PermittivityTermData> typed = {sections[k]};
    const auto raw = RawTerms(typed);
    std::vector<config::PermittivityTermData> split;
    for (auto it = raw.rbegin(); it != raw.rend(); ++it)
    {
      auto half = std::get<config::PermittivityPoleData>(*it);
      half.residue *= 0.5;
      split.push_back(half);
      split.push_back(half);
    }
    std::vector<std::unique_ptr<MaterialOperator>> ops;
    std::vector<std::unique_ptr<DomainPostOperator>> posts;
    const std::array<const std::vector<config::PermittivityTermData> *, 3> models = {
        &typed, &raw, &split};
    for (const auto *terms : models)
    {
      ops.push_back(std::make_unique<MaterialOperator>(
          std::vector<config::MaterialData>{MakeMaterial({1}, *terms),
                                            MakeMaterial({2}, *terms)},
          periodic, ProblemType::DRIVEN, mesh));
      const bool is_typed = (terms == &typed);
      REQUIRE(ops.back()->HasStoredElectricEnergy() == is_typed);
      posts.push_back(std::make_unique<DomainPostOperator>(c.postpro, *ops.back(),
                                                           c.nd_fespace, c.rt_fespace));
      REQUIRE(posts.back()->HasStoredElectricFieldEnergy() == is_typed);
    }
    CHECK_THAT(ops[1]->GetStoredElectricEnergyReason(0),
               ContainsSubstring("no stored-energy realization"));
    for (double omega : {0.3, 1.1, 4.0})
    {
      CAPTURE(omega);
      const auto eps = ops[0]->EvaluateFrequencyDependentPermittivity(0, {0.0, omega});
      CHECK(ops[0]->GetFrequencyDependentStoredEnergyPermittivity(0, omega) > eps_inf);
      for (std::size_t i = 1; i < ops.size(); i++)
      {
        CAPTURE(i);
        const auto eps_i = ops[i]->EvaluateFrequencyDependentPermittivity(0, {0.0, omega});
        CHECK(std::abs(eps_i - eps) <= 1.0e-12 * std::abs(eps));
        CHECK_THROWS(ops[i]->GetFrequencyDependentStoredEnergyPermittivity(0, omega));
        CHECK_THROWS(posts[i]->GetElectricFieldEnergy(c.E, omega, ElectricEnergy::STORED));
        const auto type = ElectricEnergy::MATERIAL;
        CHECK_THAT(posts[i]->GetElectricFieldEnergy(c.E, omega, type),
                   WithinRel(posts[0]->GetElectricFieldEnergy(c.E, omega, type), 1.0e-10));
        for (int idx : {1, 2})
        {
          CHECK_THAT(
              posts[i]->GetDomainElectricFieldEnergy(idx, c.E, omega, type),
              WithinRel(posts[0]->GetDomainElectricFieldEnergy(idx, c.E, omega, type),
                        1.0e-10));
        }
      }
    }
  }

  // Two Debye terms 1/(s + 1) + 1/(s + 2) store 1/2 + 2/5 = 0.9 at ω = 1 as typed terms.
  // The same χ has realizations with other stored energies, so as RAW terms it has none.
  const std::vector<config::PermittivityTermData> debyes = {config::DebyeData{1.0, 1.0},
                                                            config::DebyeData{0.5, 0.5}};
  MaterialOperator typed_op({MakeMaterial({1, 2}, debyes)}, periodic, ProblemType::DRIVEN,
                            mesh);
  MaterialOperator raw_op({MakeMaterial({1, 2}, RawTerms(debyes))}, periodic,
                          ProblemType::DRIVEN, mesh);
  CHECK_THAT(typed_op.GetFrequencyDependentStoredEnergyPermittivity(0, 1.0) - eps_inf,
             WithinRel(0.9, 1.0e-12));
  CHECK_FALSE(raw_op.HasStoredElectricEnergy(0));
  CHECK_THROWS(raw_op.GetFrequencyDependentStoredEnergyPermittivity(0, 1.0));
}

TEST_CASE("Frequency-dependent material without a stored electric energy",
          "[frequencydependent][Serial][Parallel]")
{
  // A RAW −0.5/(s + 1) compensated by a typed Debye 1/(1 + s) or a native
  // Djordjevic-Sarkar 2 log((s + 2)/(s + 1)) is passive, as are a RAW Debye 1/(s + 1) and
  // nearly coincident RAW poles. None has gain, and none has a specified stored electric
  // energy: the material states the reason, the stored energy is not evaluated (an error,
  // not a value), and so neither are the energies of the regions with it. Its response, the
  // field energy at Re ε, and the regions without it are unaffected. The status is that of
  // the configured model, the same on every process, also where the material (attribute 1,
  // one element) has no elements.
  test::IntegrationSettingsGuard guard;
  fem::DefaultIntegrationOrder::p_trial = 1;
  auto mesh = MakeBoxMesh(3, 2, mfem::Element::HEXAHEDRON);
  config::PeriodicBoundaryData periodic;
  EnergyConsumers c(mesh);
  using ElectricEnergy = DomainPostOperator::ElectricEnergy;
  const double delta = 1.0e-9;
  std::vector<config::PermittivityTermData> close;
  for (int k = 0; k < 4; k++)
  {
    close.push_back(config::PermittivityPoleData{-1.0 - k * delta, 1.0});
  }
  const std::vector<std::vector<config::PermittivityTermData>> cases = {
      {config::DebyeData{1.0, 1.0}, config::PermittivityPoleData{-1.0, -0.5}},
      {config::DjordjevicSarkarData{2.0, 1.0, 2.0},
       config::PermittivityPoleData{-1.0, -0.5}},
      {config::PermittivityPoleData{-1.0, 1.0}},
      close};
  config::MaterialData constant;
  constant.attributes = {2};
  constant.epsilon_r.s = {3.0, 3.0, 3.0};
  // The compensated Debye has the response of a typed Debye 0.5/(1 + s), and the regions
  // without the material, those of any other material on attribute 1.
  MaterialOperator reference_op(
      {MakeMaterial({1}, {config::DebyeData{0.5, 1.0}}), constant}, periodic,
      ProblemType::DRIVEN, mesh);
  DomainPostOperator reference_post(c.postpro, reference_op, c.nd_fespace, c.rt_fespace);
  for (std::size_t k = 0; k < cases.size(); k++)
  {
    CAPTURE(k);
    const auto &terms = cases[k];
    CHECK(permittivity::FindDefects(eps_inf, terms).gain.empty());
    MaterialOperator mat_op({MakeMaterial({1}, terms), constant}, periodic,
                            ProblemType::DRIVEN, mesh);
    CHECK_THAT(mat_op.GetStoredElectricEnergyReason(0),
               ContainsSubstring("no stored-energy realization"));
    CHECK_FALSE(mat_op.HasStoredElectricEnergy(0));
    CHECK(mat_op.HasStoredElectricEnergy(1));
    CHECK_FALSE(mat_op.HasStoredElectricEnergy());
    CHECK(mat_op.HasStoredElectricEnergy(std::vector<int>{2}));
    CHECK_FALSE(mat_op.HasStoredElectricEnergy(std::vector<int>{1, 2}));
    DomainPostOperator post(c.postpro, mat_op, c.nd_fespace, c.rt_fespace);
    CHECK_FALSE(post.HasStoredElectricFieldEnergy());
    CHECK(post.HasStoredDomainElectricFieldEnergy(1));
    CHECK_FALSE(post.HasStoredDomainElectricFieldEnergy(2));
    for (double omega : {0.3, 1.1})
    {
      CAPTURE(omega);
      const std::complex<double> s(0.0, omega);
      const auto eps = mat_op.EvaluateFrequencyDependentPermittivity(0, s);
      CHECK(std::abs(eps - (eps_inf + permittivity::Susceptibility(terms, s))) <=
            1.0e-12 * std::abs(eps));
      CHECK_THROWS(mat_op.GetFrequencyDependentStoredEnergyPermittivity(0, omega));
      CHECK_THROWS(mat_op.EvaluateStoredEnergyPermittivity(omega));
      CHECK_THROWS(post.GetElectricFieldEnergy(c.E, omega));
      CHECK_THROWS(post.GetDomainElectricFieldEnergy(2, c.E, omega));
      CHECK_THAT(
          post.GetDomainElectricFieldEnergy(1, c.E, omega),
          WithinRel(reference_post.GetDomainElectricFieldEnergy(1, c.E, omega), 1.0e-12));
      const double material =
          post.GetElectricFieldEnergy(c.E, omega, ElectricEnergy::MATERIAL);
      CHECK(std::isfinite(material));
      // Both energies from one set of products: the same values, without the stored one.
      const auto energies = post.GetElectricFieldEnergies(c.E, omega);
      CHECK(!energies.stored);
      CHECK(energies.material == material);
      const auto reference = reference_post.GetElectricFieldEnergies(c.E, omega);
      CHECK(reference.stored.has_value() == reference_post.HasStoredElectricFieldEnergy());
      if (reference.stored)
      {
        CHECK(*reference.stored == reference_post.GetElectricFieldEnergy(c.E, omega));
      }
      CHECK(reference.material ==
            reference_post.GetElectricFieldEnergy(c.E, omega, ElectricEnergy::MATERIAL));
      if (k == 0)
      {
        CHECK_THAT(material, WithinRel(reference_post.GetElectricFieldEnergy(
                                           c.E, omega, ElectricEnergy::MATERIAL),
                                       1.0e-12));
      }
    }
  }
}

TEST_CASE("Electric energies at an exact lossless resonance",
          "[frequencydependent][Serial][Parallel]")
{
  // At the resonance ω₀ = 1 of a lossless Lorentz, typed or as a RAW conjugate pair, Re χ
  // and the typed stored energy are not finite: the material weights and the full and
  // regional energies (and the weights of the reduced energies) are an error, not inf or
  // NaN. Near the resonance they are finite, the typed stored energy positive, and the
  // energies at Re ε of typed and RAW agree. The RAW pair specifies no stored energy.
  test::IntegrationSettingsGuard guard;
  fem::DefaultIntegrationOrder::p_trial = 1;
  auto mesh = MakeBoxMesh(3, 2, mfem::Element::HEXAHEDRON);
  config::PeriodicBoundaryData periodic;
  EnergyConsumers c(mesh);
  using ElectricEnergy = DomainPostOperator::ElectricEnergy;
  const std::vector<config::PermittivityTermData> typed = {
      config::LorentzData{1.5, 1.0, 0.0}};
  std::vector<std::unique_ptr<MaterialOperator>> ops;
  std::vector<std::unique_ptr<DomainPostOperator>> posts;
  for (const auto &terms : {typed, RawTerms(typed)})
  {
    ops.push_back(std::make_unique<MaterialOperator>(
        std::vector<config::MaterialData>{MakeMaterial({1, 2}, terms)}, periodic,
        ProblemType::DRIVEN, mesh));
    REQUIRE(ops.back()->HasStoredElectricEnergy(0) == (ops.size() == 1));
    posts.push_back(std::make_unique<DomainPostOperator>(c.postpro, *ops.back(),
                                                         c.nd_fespace, c.rt_fespace));
  }
  for (std::size_t i = 0; i < ops.size(); i++)
  {
    CAPTURE(i);
    const auto &mat_op = *ops[i];
    const auto &post = *posts[i];
    const double omega = 1.0;
    const auto stored = (i == 0) ? ContainsSubstring("singular")
                                 : ContainsSubstring("no stored electric energy");
    CHECK_THROWS_WITH(mat_op.GetFrequencyDependentStoredEnergySusceptibility(0, omega),
                      stored);
    CHECK_THROWS_WITH(mat_op.GetFrequencyDependentStoredEnergyPermittivity(0, omega),
                      stored);
    CHECK_THROWS_WITH(mat_op.GetFrequencyDependentSusceptibilityReal(0, omega),
                      ContainsSubstring("singular"));
    const auto types = (i == 0)
                           ? std::vector{ElectricEnergy::STORED, ElectricEnergy::MATERIAL}
                           : std::vector{ElectricEnergy::MATERIAL};
    for (auto type : types)
    {
      CAPTURE(static_cast<int>(type));
      CHECK_THROWS_WITH(post.GetFrequencyDependentEnergyWeight(0, omega, type),
                        ContainsSubstring("singular"));
      CHECK_THROWS_WITH(post.GetElectricFieldEnergy(c.E, omega, type),
                        ContainsSubstring("singular"));
      for (int idx : {1, 2})
      {
        CHECK_THROWS_WITH(post.GetDomainElectricFieldEnergy(idx, c.E, omega, type),
                          ContainsSubstring("singular"));
      }
    }
  }
  const double near = 1.01;
  for (auto type : {ElectricEnergy::STORED, ElectricEnergy::MATERIAL})
  {
    CAPTURE(static_cast<int>(type));
    const double energy = posts[0]->GetElectricFieldEnergy(c.E, near, type);
    REQUIRE(std::isfinite(energy));
    if (type == ElectricEnergy::STORED)
    {
      CHECK(energy > 0.0);
      continue;
    }
    CHECK_THAT(posts[1]->GetElectricFieldEnergy(c.E, near, type),
               WithinRel(energy, 1.0e-10));
    for (int idx : {1, 2})
    {
      CHECK_THAT(
          posts[1]->GetDomainElectricFieldEnergy(idx, c.E, near, type),
          WithinRel(posts[0]->GetDomainElectricFieldEnergy(idx, c.E, near, type), 1.0e-10));
    }
  }
}

}  // namespace palace
