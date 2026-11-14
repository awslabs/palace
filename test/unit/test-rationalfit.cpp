// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <complex>
#include <numbers>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include "utils/communication.hpp"
#include "utils/rationalfit.hpp"

namespace palace
{

namespace
{

using Complex = std::complex<double>;
constexpr double two_pi = 2.0 * std::numbers::pi;
const Complex i(0.0, 1.0);

// The Lorentz medium of the absorbing boundary box of the eigenmode driver tests, in
// λ = i 2π f with f in GHz: ε(λ) = ε∞ + Δε ω₀² / (ω₀² + γ λ + λ²).
constexpr double eps_inf = 2.0, delta_eps = 0.1, w0 = two_pi * 14.0, gamma = two_pi * 0.5;
Complex LorentzScale(Complex l)
{
  return std::sqrt(1.0 + (delta_eps / eps_inf) * w0 * w0 / (w0 * w0 + gamma * l + l * l));
}

// The poles and residues of ε / ε∞ - 1 of the Lorentz medium.
void LorentzPoles(std::vector<Complex> &p, std::vector<Complex> &r)
{
  const Complex pole(-0.5 * gamma, std::sqrt(w0 * w0 - 0.25 * gamma * gamma));
  p = {pole, std::conj(pole)};
  const double a = delta_eps / eps_inf * w0 * w0;
  r = {a / (p[0] - p[1]), a / (p[1] - p[0])};
}

// The largest error of the fit, relative to the largest |f|, on a dense polar grid of the
// window where it is valid, independent of the points of the fit and of its check.
double ValidError(const utils::ScalarFit &fit, const utils::FitWindow &window,
                  const std::function<Complex(Complex)> &f, int &num_valid,
                  int &num_invalid)
{
  double err = 0.0, scale = 0.0;
  num_valid = num_invalid = 0;
  const double phi0 = -window.margin, phi1 = 0.5 * std::numbers::pi + window.margin;
  for (int a = 0; a <= 400; a++)
  {
    for (int b = 0; b <= 400; b++)
    {
      const auto w = window.omega_t + std::polar(window.radius * (a + 0.37) / 400.37,
                                                 phi0 + (phi1 - phi0) * b / 400.0);
      if (!fit.IsValid(window, w))
      {
        num_invalid++;
        continue;
      }
      num_valid++;
      const auto l = i * w;
      err = std::max(err, std::abs(fit.Eval(l) - f(l)));
      scale = std::max(scale, std::abs(f(l)));
    }
  }
  return err / scale;
}

}  // namespace

TEST_CASE("Rational fit of the branch cut of the square root of a Lorentz permittivity",
          "[rationalfit][Serial]")
{
  // The window of the Lorentz box (target 8 GHz, TargetUpper 24 GHz) holds the branch
  // cut of √(ε / ε∞) between the material pole and the permittivity zero, at Im f = 0.25
  // GHz. No rational fit follows the cut: the fit is valid outside a strip around it,
  // which the check must find, and the independent dense check must confirm.
  const utils::FitWindow window = {two_pi * 8.0, two_pi * 16.0};
  std::vector<Complex> p, r;
  LorentzPoles(p, r);
  const auto cut = utils::SquareRootCut(p, r, window.radius / 4096.0);
  // The cut runs from the zeros of ε to its poles, and ε is real and not positive on it.
  const Complex zero = i * two_pi * Complex(14.343554, 0.25);
  CHECK(std::ranges::any_of(cut, [&](Complex w)
                            { return std::abs(w - p[0] / i) < 1e-12 * w0; }));
  CHECK(std::ranges::any_of(cut,
                            [&](Complex w) { return std::abs(i * w - zero) < 1e-4 * w0; }));
  for (const auto w : cut)
  {
    if (std::abs(i * w - p[0]) < 1.0e-12 * w0 || std::abs(i * w - p[1]) < 1.0e-12 * w0)
    {
      continue;  // A pole, where e is not defined
    }
    const auto e = LorentzScale(i * w) * LorentzScale(i * w);
    CAPTURE(w / two_pi, e);
    CHECK(std::abs(e.imag()) <= 1.0e-8 * std::max(1.0, std::abs(e)));
    CHECK(e.real() <= 1.0e-8);
  }

  // The damping term ω r(ω) of the first-order absorbing boundary next to the medium.
  constexpr double tol = 1.0e-10;
  auto g = [](Complex l) { return LorentzScale(l); };
  const auto fit = utils::FitScalar(g, -i, 1, cut, window, tol);
  REQUIRE(fit.success);
  CHECK(fit.error <= tol);
  CHECK(fit.strip > 0.0);
  CHECK(fit.strip <= window.radius / 8.0);
  CHECK(fit.poles.size() <= 24);
  CHECK(fit.num_samples < fit.num_checks);
  int num_valid, num_invalid;
  auto f = [](Complex l) { return -i * l * LorentzScale(l); };
  const double err = ValidError(fit, window, f, num_valid, num_invalid);
  CAPTURE(fit.strip / two_pi, fit.poles.size(), fit.num_samples, fit.num_checks, fit.error,
          err, num_valid, num_invalid);
  CHECK(err <= 10.0 * tol);
  CHECK(num_invalid > 0);
  CHECK(num_valid > 100 * num_invalid);
  // And on the boundary of the strip, just outside it, where the error is largest: the cut
  // is the segment from the pole to the zero (Re λ = -γ / 2).
  double err_strip = 0.0, scale = 0.0;
  const Complex a = p[0] / i, b = zero / i;
  const double d = 1.0001 * fit.strip, L = std::abs(b - a);
  const double perimeter = 2.0 * L + 2.0 * std::numbers::pi * d;
  for (int k = 0; k < 4000; k++)
  {
    // The stadium of half-width d around the segment: two sides and two half circles.
    const double u = perimeter * k / 4000.0;
    Complex w;
    if (u < L)
    {
      w = a + u + Complex(0.0, d);
    }
    else if (u < L + std::numbers::pi * d)
    {
      w = b + std::polar(d, 0.5 * std::numbers::pi - (u - L) / d);
    }
    else if (u < 2.0 * L + std::numbers::pi * d)
    {
      w = b - (u - L - std::numbers::pi * d) - Complex(0.0, d);
    }
    else
    {
      w = a +
          std::polar(d, -0.5 * std::numbers::pi - (u - 2.0 * L - std::numbers::pi * d) / d);
    }
    REQUIRE(fit.IsValid(window, w));
    err_strip = std::max(err_strip, std::abs(fit.Eval(i * w) - f(i * w)));
    scale = std::max(scale, std::abs(f(i * w)));
  }
  CHECK(err_strip <= 10.0 * tol * scale);
  Mpi::Print("Lorentz fit: strip {:.4e} GHz, {:d} poles, {:d} samples, {:d} check points, "
             "error {:.2e}; independent errors {:.2e} on a polar grid ({:d} valid, {:d} "
             "in the strip), {:.2e} next to the strip\n",
             fit.strip / two_pi, fit.poles.size(), fit.num_samples, fit.num_checks,
             fit.error, err, num_valid, num_invalid, err_strip / scale);
  // A grid point in the strip, next to the cut, is not valid.
  CHECK_FALSE(fit.IsValid(window, (p[0] + zero) / (2.0 * i) + Complex(0.0, 1.0e-3)));
}

TEST_CASE("Rational fit of functions analytic in the window", "[rationalfit][Serial]")
{
  const utils::FitWindow window = {two_pi * 18.0, two_pi * 36.0};
  constexpr double tol = 1.0e-10;
  int num_valid, num_invalid;

  // The Debye medium of the absorbing boundary box: its cut is on the imaginary axis of f,
  // outside the window, so the fit is valid in all of it with few poles. The curl-curl
  // term 0.5 / (ω r(ω)) = 0.5 i / (λ r(λ)).
  constexpr double delta = 0.5, tau = 0.01;  // ns, with λ in rad / ns
  auto debye = [=](Complex l)
  { return std::sqrt(1.0 + (delta / eps_inf) / (1.0 + tau * l)); };
  const std::vector<Complex> p = {-1.0 / tau}, r = {(delta / eps_inf) / tau};
  const auto cut = utils::SquareRootCut(p, r, window.radius / 4096.0);
  const auto fit = utils::FitScalar([&](Complex l) { return 1.0 / debye(l); }, 0.5 * i, -1,
                                    cut, window, tol);
  REQUIRE(fit.success);
  CHECK(fit.strip == 0.0);
  CHECK(fit.poles.size() <= 8);
  const double err_debye = ValidError(
      fit, window, [&](Complex l) { return 0.5 * i / (l * debye(l)); }, num_valid,
      num_invalid);
  CHECK(err_debye <= 10.0 * tol);
  CHECK(num_invalid == 0);
  Mpi::Print("Debye fit: {:d} poles, {:d} samples, {:d} check points, error {:.2e}; "
             "independent error {:.2e}\n",
             fit.poles.size(), fit.num_samples, fit.num_checks, fit.error, err_debye);

  // The Djordjevic-Sarkar term λ² S log((λ + ω_u) / (λ + ω_l)), fitted as λ² g(λ): the cut
  // of the logarithm is on the imaginary axis of ω, from i ω_l to i ω_u.
  constexpr double strength = 0.05, w_lo = two_pi * 0.1, w_hi = two_pi * 100.0;
  auto ds = [=](Complex l) { return strength * std::log((l + w_hi) / (l + w_lo)); };
  std::vector<Complex> ds_cut = {i * w_lo, i * w_hi};
  for (int k = 1; k < 4096; k++)
  {
    ds_cut.push_back(i * (w_lo + (w_hi - w_lo) * k / 4096.0));
  }
  const auto ds_fit = utils::FitScalar(ds, 1.0, 2, ds_cut, window, tol);
  REQUIRE(ds_fit.success);
  CHECK(ds_fit.strip == 0.0);
  const double err_ds = ValidError(
      ds_fit, window, [&](Complex l) { return l * l * ds(l); }, num_valid, num_invalid);
  CHECK(err_ds <= 10.0 * tol);
  CHECK(num_invalid == 0);
  Mpi::Print("Djordjevic-Sarkar fit: {:d} poles, {:d} samples, {:d} check points, error "
             "{:.2e}; independent error {:.2e}\n",
             ds_fit.poles.size(), ds_fit.num_samples, ds_fit.num_checks, ds_fit.error,
             err_ds);

  // The skin effect term √λ of a finite conductivity boundary, fitted as λ · λ^(-1/2): its
  // branch point is at λ = 0, and its cut on the negative real axis of f.
  std::vector<Complex> ray;
  for (int k = 0; k <= 64; k++)
  {
    ray.push_back(-window.radius * k / 16.0);
  }
  auto skin = [](Complex l) { return std::sqrt(-i * l); };
  const auto conductivity =
      utils::FitScalar([&](Complex l) { return skin(l) / l; }, 1.0, 1, ray, window, tol);
  REQUIRE(conductivity.success);
  CHECK(conductivity.strip == 0.0);
  const double err_skin = ValidError(conductivity, window, skin, num_valid, num_invalid);
  CHECK(err_skin <= 10.0 * tol);
  CHECK(num_invalid == 0);
  Mpi::Print("Skin-effect fit: {:d} poles, {:d} samples, {:d} check points, error {:.2e}; "
             "independent error {:.2e}\n",
             conductivity.poles.size(), conductivity.num_samples, conductivity.num_checks,
             conductivity.error, err_skin);
}

}  // namespace palace
