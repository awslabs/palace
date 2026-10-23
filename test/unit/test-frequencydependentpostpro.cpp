// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Consumers of the complete permittivity ε(iω) and the stored-energy permittivity W(ω) of
// frequency-dependent materials: material evaluation, postprocessing, and the flux error
// estimators. Expected values come from the model formulas, independent of the
// implementation (internal units, ε₀ = 1).

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
#include "fem/point_field_evaluator.hpp"
#include "fixtures.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/farfieldboundaryoperator.hpp"
#include "models/floquetportoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/spaceoperator.hpp"
#include "models/surfacepostoperator.hpp"
#include "test-helpers.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
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
  // The absorbing boundary scale r(ω) = √(ε(iω)/ε∞) of the neighboring material uses the
  // complete permittivity, including the zero pole folded into the conductivity: its
  // imaginary part is the loss which A2(s)/s² alone omits (and whose omission reversed the
  // sign of the boundary loss). A configured model with only a zero pole has no nonlinear
  // term but still a frequency-dependent boundary. At complex ω, ε is continued
  // analytically.
  auto mesh = MakeBoxMesh(3, 2, mfem::Element::HEXAHEDRON);
  config::PeriodicBoundaryData periodic;
  for (const auto &model : {Models()[0], Models()[3]})
  {
    CAPTURE(model.name);
    MaterialOperator mat_op({MakeMaterial({1, 2}, model.terms)}, periodic,
                            ProblemType::DRIVEN, mesh);
    config::FarfieldBoundaryData farfield;
    farfield.attributes = {1, 2, 3, 4, 5, 6};
    FarfieldBoundaryOperator farfield_op(farfield, ProblemType::DRIVEN, mat_op, mesh.Get());
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
      const auto r = std::sqrt((eps_inf + chi) / eps_inf);
      const auto value = farfield_op.GetDispersiveScale(0, omega);
      CHECK_THAT(value.real(), WithinRel(r.real(), 1.0e-12));
      CHECK_THAT(value.imag(), WithinRel(r.imag(), 1.0e-12));
      // Passive at real ω: the boundary dissipates, Re{iω r} ≤ 0 is the loss sign of the
      // first-order term i ω r (here r has Im r < 0).
      if (omega.imag() == 0.0)
      {
        CHECK(value.imag() < 0.0);
      }
    }
  }
}

}  // namespace palace
