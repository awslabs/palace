// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <cmath>
#include <complex>
#include <numbers>
#include <utility>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "utils/conductivitymodel.hpp"
#include "utils/configfile.hpp"

namespace palace
{

using namespace Catch::Matchers;

namespace
{

// Complete elliptic integrals K(k) and E(k) by the arithmetic-geometric mean.
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

}  // namespace

TEST_CASE("Mattis-Bardeen conductivity", "[conductivitymodel][Serial]")
{
  SECTION("Zero temperature closed forms")
  {
    // Mattis and Bardeen (1958), Eqs. (3.11)-(3.13), with x = ħω/Δ.
    for (double x : {0.2, 1.0, 1.8, 2.4, 4.0, 10.0})
    {
      const auto sigma = conductivity::MattisBardeenRatio(x, 0.0);
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
      const auto sigma = conductivity::MattisBardeenRatio(x, 0.25);
      CHECK_THAT(sigma.real(), WithinRel(sigma_1, 1.0e-9));
      CHECK_THAT(-sigma.imag(), WithinRel(sigma_2, 1.0e-9));
    }
  }

  SECTION("London and normal-state limits")
  {
    // ħω σ₂/σₙ → π Δ tanh(Δ/2kB T) as ω → 0, and σ → σₙ for ħω >> Δ.
    const double theta = 0.3;
    const auto sigma_low = conductivity::MattisBardeenRatio(1.0e-5, theta);
    CHECK_THAT(-1.0e-5 * sigma_low.imag(),
               WithinRel(std::numbers::pi * std::tanh(0.5 / theta), 1.0e-5));
    const auto sigma_high = conductivity::MattisBardeenRatio(2.0e3, theta);
    CHECK_THAT(sigma_high.real(), WithinAbs(1.0, 1.0e-4));
    CHECK(-sigma_high.imag() > 0.0);
    CHECK(-sigma_high.imag() < 1.0e-3);
  }

  SECTION("Conductivity and frequency derivative")
  {
    // σ(ω) = σₙ r(ω/ω_gap, ω_T/ω_gap), and the derivative against a Richardson-extrapolated
    // difference with larger steps. Well below the gap, -ω dσ₂/dω = σ₂ (London).
    const config::MattisBardeenData mb{5.0, 2.0, 0.5};
    for (double omega : {0.3, 1.5, 3.0, 6.0})
    {
      CHECK(conductivity::Conductivity(mb, omega) ==
            5.0 * conductivity::MattisBardeenRatio(omega / 2.0, 0.25));
      auto difference = [&](double h)
      {
        return (conductivity::Conductivity(mb, omega * (1.0 + h)) -
                conductivity::Conductivity(mb, omega * (1.0 - h))) /
               (2.0 * h * omega);
      };
      const auto ref = (4.0 * difference(5.0e-3) - difference(1.0e-2)) / 3.0;
      const auto dsigma = conductivity::ConductivityDerivative(mb, omega);
      CHECK_THAT(dsigma.real(), WithinRel(ref.real(), 1.0e-6));
      CHECK_THAT(dsigma.imag(), WithinRel(ref.imag(), 1.0e-6));
    }
    const config::MattisBardeenData cold{5.0, 2.0, 0.1};
    const double omega = 0.02;
    CHECK_THAT(omega * conductivity::ConductivityDerivative(cold, omega).imag(),
               WithinRel(-conductivity::Conductivity(cold, omega).imag(), 1.0e-3));
  }
}

TEST_CASE("BCS gap", "[conductivitymodel][Serial]")
{
  CHECK(conductivity::BCSGapRatio(0.0) == 1.0);
  CHECK(conductivity::BCSGapRatio(1.0) == 0.0);
  CHECK(conductivity::BCSGapRatio(1.5) == 0.0);
  // Reference roots of the weak-coupling gap equation (independent SciPy evaluation).
  CHECK_THAT(conductivity::BCSGapRatio(0.3), WithinRel(0.9971237624901917, 1.0e-10));
  CHECK_THAT(conductivity::BCSGapRatio(0.5), WithinRel(0.9568847158748907, 1.0e-10));
  CHECK_THAT(conductivity::BCSGapRatio(0.9), WithinRel(0.5263419345511171, 1.0e-10));
  CHECK_THAT(conductivity::BCSGapRatio(0.99), WithinRel(0.17295385407970149, 1.0e-9));
  // Ginzburg-Landau limit Δ(T)/Δ(0) ≈ 1.74 √(1 - t) near Tc.
  CHECK_THAT(conductivity::BCSGapRatio(0.9999), WithinRel(1.74 * 0.01, 1.0e-2));
}

TEST_CASE("Tabulated conductivity", "[conductivitymodel][Serial]")
{
  SECTION("Interpolation")
  {
    // Two-fluid data: σ₁ constant and ωσ₂ constant are reproduced exactly between the
    // nodes, and held constant outside the table.
    config::TabulatedConductivityData data{
        {1.0, 2.0, 5.0, 10.0}, {3.0, 3.0, 3.0, 3.0}, {4.0, 2.0, 0.8, 0.4}};
    for (double omega : {1.0, 1.3, 2.0, 3.7, 9.9, 10.0})
    {
      const auto sigma = conductivity::Conductivity(data, omega);
      CHECK_THAT(sigma.real(), WithinRel(3.0, 1.0e-14));
      CHECK_THAT(-omega * sigma.imag(), WithinRel(4.0, 1.0e-14));
    }
    CHECK_THAT(-0.5 * conductivity::Conductivity(data, 0.5).imag(),
               WithinRel(4.0, 1.0e-14));
    CHECK_THAT(-20.0 * conductivity::Conductivity(data, 20.0).imag(),
               WithinRel(4.0, 1.0e-14));

    // Nodes are interpolated, and the interpolant is monotone between monotone data.
    data.sigma1 = {1.0, 2.0, 2.5, 6.0};
    for (std::size_t i = 0; i < data.omega.size(); i++)
    {
      CHECK_THAT(conductivity::Conductivity(data, data.omega[i]).real(),
                 WithinRel(data.sigma1[i], 1.0e-14));
    }
    double prev = conductivity::Conductivity(data, 1.0).real();
    for (double omega = 1.05; omega <= 10.0; omega += 0.05)
    {
      const double s = conductivity::Conductivity(data, omega).real();
      CHECK(s >= prev - 1.0e-14);
      prev = s;
    }
  }

  SECTION("Frequency derivative")
  {
    // σ₁ linear in log(ω) and ωσ₂ constant are reproduced exactly, so dσ/dω =
    // 2/ω + 4i/ω² (σ = σ₁ - iσ₂), and dσ₁/dω = 0 outside the table.
    config::TabulatedConductivityData data;
    data.omega = {1.0, 2.0, 5.0, 10.0};
    for (double omega : data.omega)
    {
      data.sigma1.push_back(1.0 + 2.0 * std::log(omega));
      data.sigma2.push_back(4.0 / omega);
    }
    for (double omega : {1.0, 1.3, 2.0, 3.7, 9.9, 10.0})
    {
      const auto dsigma = conductivity::ConductivityDerivative(data, omega);
      CHECK_THAT(dsigma.real(), WithinRel(2.0 / omega, 1.0e-12));
      CHECK_THAT(dsigma.imag(), WithinRel(4.0 / (omega * omega), 1.0e-12));
    }
    const auto dsigma = conductivity::ConductivityDerivative(data, 20.0);
    CHECK(dsigma.real() == 0.0);
    CHECK_THAT(dsigma.imag(), WithinRel(4.0 / 400.0, 1.0e-14));
  }
}

}  // namespace palace
