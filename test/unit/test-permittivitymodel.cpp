// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <complex>
#include <limits>
#include <numbers>
#include <random>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "utils/configfile.hpp"
#include "utils/permittivitymodel.hpp"

using namespace palace;
using Catch::Matchers::WithinRel;

namespace
{

// A section of the model in its auxiliary differential equation form, driven by E(t) =
// cos(ωt): a relaxation ẋ + a₀x = E with P = b₀x, or an oscillator ẍ + a₁ẋ + a₀x = E with
// P = b₀x + b₁ẋ (ε₀ = 1).
struct Section
{
  bool oscillator;
  double a0, a1, b0, b1;
};

// Stored energy and dissipated power of the section's state, the forms of
// utils/permittivitymodel.hpp, and the power E Ṗ delivered to the polarization, along the
// steady state at time t.
struct Powers
{
  double W, dW, D, EdP;
};
Powers SteadyState(const Section &sec, double omega, double t)
{
  const std::complex<double> s(0.0, omega), phase = std::exp(s * t);
  const std::complex<double> xh =
      sec.oscillator ? 1.0 / (s * s + sec.a1 * s + sec.a0) : 1.0 / (s + sec.a0);
  const double x = (xh * phase).real(), xd = (s * xh * phase).real(),
               xdd = (s * s * xh * phase).real(), xddd = (s * s * s * xh * phase).real();
  const double E = std::cos(omega * t);
  if (!sec.oscillator)
  {
    return {0.5 * sec.a0 * sec.b0 * x * x, sec.a0 * sec.b0 * x * xd, sec.b0 * xd * xd,
            E * sec.b0 * xd};
  }
  const double c = sec.b0 + sec.a1 * sec.b1;
  return {0.5 * c * xd * xd + sec.a0 * sec.b1 * x * xd + 0.5 * sec.a0 * sec.b0 * x * x,
          c * xd * xdd + sec.a0 * sec.b1 * (xd * xd + x * xdd) + sec.a0 * sec.b0 * x * xd,
          sec.b1 * xdd * xdd + (sec.a1 * sec.b0 - sec.a0 * sec.b1) * xd * xd,
          E * (sec.b0 * xd + sec.b1 * xdd)};
}

// Time average of the stored energy over a period, by the trapezoidal rule (exact for a
// trigonometric polynomial of this degree).
double AverageStoredEnergy(const Section &sec, double omega)
{
  constexpr int n = 64;
  const double period = 2.0 * std::numbers::pi / omega;
  double sum = 0.0;
  for (int k = 0; k < n; k++)
  {
    const auto p = SteadyState(sec, omega, period * k / n);
    sum += p.W;
    // The exact power balance E Ṗ = dW/dt + D at every time.
    CHECK_THAT(p.EdP, Catch::Matchers::WithinAbs(
                          p.dW + p.D, 1.0e-12 * (std::abs(p.EdP) + std::abs(p.D))));
  }
  return sum / n;
}

}  // namespace

TEST_CASE("Permittivity model stored energy is that of each term's model",
          "[permittivitymodel][Serial]")
{
  // The time-averaged stored energy of E(t) = cos(ωt) is ¼ε₀ w(ω), for a permittivity
  // ε∞ + w(ω), so the time average of each section's W is w/4.
  const double omega = 1.7;
  auto Check = [omega](const config::PermittivityTermData &term, const Section &sec)
  {
    const std::vector<config::PermittivityTermData> terms = {term};
    CHECK_THAT(4.0 * AverageStoredEnergy(sec, omega),
               WithinRel(permittivity::StoredEnergySusceptibility(terms, omega), 1.0e-12));
  };
  SECTION("Typed terms")
  {
    // Debye Δε/(1 + sτ) = (Δε/τ)/(s + 1/τ).
    Check(config::DebyeData{1.5, 0.4}, {false, 1.0 / 0.4, 0.0, 1.5 / 0.4, 0.0});
    // Lorentz Δε ω₀²/(s² + γs + ω₀²), under- and overdamped.
    for (const double gamma : {0.3, 7.0})
    {
      Check(config::LorentzData{0.8, 2.5, gamma},
            {true, 2.5 * 2.5, gamma, 0.8 * 2.5 * 2.5, 0.0});
    }
    // Drude ω_p²/(s² + γs), with the kinetic energy of the carriers.
    Check(config::DrudeData{3.0, 0.2}, {true, 0.0, 0.2, 9.0, 0.0});
  }
  SECTION("PoleResidue terms")
  {
    // A real pole r/(s - p) is a relaxation, and a conjugate pair r/(s - p) + conj an
    // oscillator with a₁ = -2 Re p, a₀ = |p|², b₁ = 2 Re r, and b₀ = -2 Re(r p̄).
    Check(config::PermittivityPoleData{{-0.7, 0.0}, {2.0, 0.0}},
          {false, 0.7, 0.0, 2.0, 0.0});
    const std::complex<double> p(-0.4, 2.2), r(0.3, -1.1);
    Check(config::PermittivityPoleData{p, r},
          {true, std::norm(p), -2.0 * p.real(), -2.0 * (r * std::conj(p)).real(),
           2.0 * r.real()});
  }
  SECTION("Djordjevic-Sarkar terms are a continuum of relaxations")
  {
    // Strength log((s + ω_hi)/(s + ω_lo)) = ∫ Strength/(s + a) da over [ω_lo, ω_hi], with
    // the energies of the relaxations b₀/(s + a₀) = Strength/(s + a) summed by
    // Gauss-Legendre quadrature in log a.
    const double strength = 0.2, lo = 0.01, hi = 100.0;
    const std::vector<config::PermittivityTermData> terms = {
        config::DjordjevicSarkarData{strength, lo, hi}};
    constexpr int n = 400;
    double w = 0.0;
    for (int k = 0; k < n; k++)
    {
      // Midpoint rule in u = log a, with da = a du.
      const double a = lo * std::exp((k + 0.5) * std::log(hi / lo) / n);
      const double da = a * std::log(hi / lo) / n;
      w += 4.0 * AverageStoredEnergy({false, a, 0.0, strength * da, 0.0}, omega);
    }
    CHECK_THAT(w,
               WithinRel(permittivity::StoredEnergySusceptibility(terms, omega), 1.0e-4));
  }
}

TEST_CASE("Permittivity model pole form and susceptibility agree",
          "[permittivitymodel][Serial]")
{
  const std::vector<config::PermittivityTermData> terms = {
      config::DrudeData{3.0, 0.2},
      config::DebyeData{1.5, 0.4},
      config::LorentzData{0.8, 2.5, 0.3},
      config::LorentzData{0.5, 1.0, 3.0},
      config::PermittivityPoleData{{-0.4, 2.2}, {0.3, -1.1}},
      config::PermittivityPoleData{{-0.7, 0.0}, {-2.0, 0.0}},
      config::DjordjevicSarkarData{0.2, 0.01, 100.0}};
  for (const std::complex<double> s : {std::complex<double>(0.0, 1.7), {0.3, 0.9}})
  {
    std::complex<double> chi = 0.0;
    for (const auto &term : permittivity::PoleTerms(terms))
    {
      chi += term.residue / (s - term.pole);
      if (term.pole.imag() > 0.0)
      {
        chi += std::conj(term.residue) / (s - std::conj(term.pole));
      }
    }
    for (const auto &term : permittivity::DjordjevicSarkarTerms(terms))
    {
      chi += term.strength * std::log((s + term.omega_upper) / (s + term.omega_lower));
    }
    CHECK(std::abs(chi - permittivity::Susceptibility(terms, s)) <=
          1.0e-12 * std::abs(chi));
  }
  // χ(0) is the static susceptibility, infinite with free carriers.
  const std::vector<config::PermittivityTermData> bound(terms.begin() + 1, terms.end());
  CHECK_THAT(permittivity::StaticSusceptibility(bound),
             WithinRel(permittivity::Susceptibility(bound, 0.0).real(), 1.0e-12));
  CHECK(std::isinf(permittivity::StaticSusceptibility(terms)));
}

TEST_CASE("Permittivity model defects", "[permittivitymodel][Serial]")
{
  SECTION("Typed models are passive")
  {
    const std::vector<config::PermittivityTermData> terms = {
        config::DrudeData{3.0, 0.2}, config::DebyeData{1.5, 0.4},
        config::LorentzData{0.8, 2.5, 0.3}, config::LorentzData{0.8, 2.5, 0.0},
        config::DjordjevicSarkarData{0.2, 0.01, 100.0}};
    CHECK(permittivity::IsPassiveByConstruction(terms));
    const auto defects = permittivity::FindDefects(1.0, terms);
    CHECK(defects.gain.empty());
    CHECK(defects.negative_energy.empty());
  }

  SECTION("A negative relaxation has gain at every frequency")
  {
    const std::vector<config::PermittivityTermData> terms = {
        config::PermittivityPoleData{{-1.0, 0.0}, {-0.5, 0.0}}};
    const auto defects = permittivity::FindDefects(1.0, terms);
    REQUIRE(defects.gain.size() == 1);
    CHECK(defects.gain[0].lower == 0.0);
    CHECK(std::isinf(defects.gain[0].upper));
    CHECK(defects.negative_energy.empty());  // ε∞ = 1 exceeds -Δε = 0.5
  }

  SECTION("Intervals match dense sampling")
  {
    // Passive Lorentz terms with a narrow negative pair, whose gain is not compensated near
    // its resonance, and a negative relaxation, which drives the stored energy negative at
    // low frequencies for a small ε∞.
    const std::vector<config::PermittivityTermData> terms = {
        config::LorentzData{0.8, 2.5, 0.3},
        config::PermittivityPoleData{{-1.0e-3, 4.0}, {0.0, 1.0e-3}},
        config::PermittivityPoleData{{-0.05, 0.0}, {-0.1, 0.0}},
        config::DjordjevicSarkarData{0.05, 0.01, 100.0}};
    const double eps_inf = 0.5;
    const auto defects = permittivity::FindDefects(eps_inf, terms);
    REQUIRE_FALSE(defects.gain.empty());
    REQUIRE_FALSE(defects.negative_energy.empty());

    auto Inside = [](const std::vector<permittivity::FrequencyInterval> &intervals,
                     double omega, double margin)
    {
      for (const auto &i : intervals)
      {
        if (omega > i.lower * (1.0 + margin) && omega < i.upper * (1.0 - margin))
        {
          return 1;
        }
        if (omega >= i.lower * (1.0 - margin) && omega <= i.upper * (1.0 + margin))
        {
          return 0;  // At an endpoint, within the margin
        }
      }
      return -1;
    };
    constexpr int n = 200000;
    int checked = 0;
    for (int k = 0; k < n; k++)
    {
      const double omega = 1.0e-6 * std::pow(10.0, 12.0 * k / n);
      const double d = permittivity::DissipationSusceptibility(terms, omega);
      const double w = eps_inf + permittivity::StoredEnergySusceptibility(terms, omega);
      const int in_gain = Inside(defects.gain, omega, 1.0e-6);
      const int in_energy = Inside(defects.negative_energy, omega, 1.0e-6);
      if (in_gain != 0)
      {
        CAPTURE(omega, d);
        CHECK((d < 0.0) == (in_gain == 1));
        checked++;
      }
      if (in_energy != 0)
      {
        CAPTURE(omega, w);
        CHECK((w < 0.0) == (in_energy == 1));
      }
    }
    CHECK(checked > n / 2);
  }
}
