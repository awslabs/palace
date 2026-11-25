// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "conductivitymodel.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <queue>
#include <mfem.hpp>
#include "utils/configfile.hpp"

namespace palace::conductivity
{

namespace
{

// Globally adaptive 7-point Gauss / 15-point Kronrod quadrature of a bounded integrand on
// [a, b], to a relative tolerance on the total error estimate.
template <typename Function>
double IntegrateGK(const Function &f, double a, double b, double rel_tol = 1.0e-13)
{
  static constexpr double xgk[8] = {
      0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
      0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
      0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
      0.207784955007898467600689403773245, 0.000000000000000000000000000000000};
  static constexpr double wgk[8] = {
      0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
      0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
      0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
      0.204432940075298892414161999234649, 0.209482141084727828012999174891714};
  static constexpr double wg[4] = {
      0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
      0.381830050505118944950369775488975, 0.417959183673469387755102040816327};
  constexpr std::size_t max_intervals = 4096;
  struct Segment
  {
    double a, b, integral, error;
    bool operator<(const Segment &other) const { return error < other.error; }
  };
  auto rule = [&f](double a, double b)
  {
    const double c = 0.5 * (a + b), h = 0.5 * (b - a);
    const double fc = f(c);
    double res_k = fc * wgk[7], res_g = fc * wg[3];
    for (int j = 0; j < 7; j++)
    {
      const double fsum = f(c - h * xgk[j]) + f(c + h * xgk[j]);
      res_k += wgk[j] * fsum;
      if (j % 2 == 1)
      {
        res_g += wg[j / 2] * fsum;
      }
    }
    return Segment{a, b, res_k * h, std::abs((res_k - res_g) * h)};
  };
  std::priority_queue<Segment> segments;
  segments.push(rule(a, b));
  double integral = segments.top().integral, error = segments.top().error;
  while (error > rel_tol * std::abs(integral) && segments.size() < max_intervals)
  {
    const Segment s = segments.top();
    segments.pop();
    const double c = 0.5 * (s.a + s.b);
    if (!(c > s.a && c < s.b))
    {
      // Not divisible in floating point: keep it, without its error estimate.
      segments.push({s.a, s.b, s.integral, 0.0});
      error -= s.error;
      continue;
    }
    const Segment l = rule(s.a, c), r = rule(c, s.b);
    integral += l.integral + r.integral - s.integral;
    error += l.error + r.error - s.error;
    segments.push(l);
    segments.push(r);
  }
  return integral;
}

// 1 - 2f(e) = tanh(e/2θ) for the Fermi function f, and sign(e) at θ = 0.
double OneMinusTwoFermi(double e, double theta)
{
  return (theta > 0.0) ? std::tanh(0.5 * e / theta) : ((e > 0.0) - (e < 0.0));
}

// Monotone cubic Hermite (Fritsch-Carlson) slope at node i of the data y(k) at strictly
// increasing x(k), k < n, from its neighbors.
template <typename X, typename Y>
double PCHIPSlope(std::size_t n, std::size_t i, const X &x, const Y &y)
{
  auto delta = [&](std::size_t k) { return (y(k + 1) - y(k)) / (x(k + 1) - x(k)); };
  if (n == 2)
  {
    return delta(0);
  }
  if (i == 0 || i == n - 1)
  {
    // Three-point end slope, limited to preserve the shape.
    const std::size_t k = (i == 0) ? 0 : n - 2, l = (i == 0) ? 1 : n - 3;
    const double h0 = x(k + 1) - x(k), h1 = x(l + 1) - x(l);
    const double del0 = delta(k), del1 = delta(l);
    const double s = ((2.0 * h0 + h1) * del0 - h0 * del1) / (h0 + h1);
    if (s * del0 <= 0.0)
    {
      return 0.0;
    }
    return (del0 * del1 <= 0.0 && std::abs(s) > std::abs(3.0 * del0)) ? 3.0 * del0 : s;
  }
  const double del0 = delta(i - 1), del1 = delta(i);
  if (del0 * del1 <= 0.0)
  {
    return 0.0;
  }
  const double w1 = 2.0 * (x(i + 1) - x(i)) + (x(i) - x(i - 1)),
               w2 = (x(i + 1) - x(i)) + 2.0 * (x(i) - x(i - 1));
  return (w1 + w2) / (w1 / del0 + w2 / del1);
}

// Value, or derivative in x, of the PCHIP interpolant of y(k) at xq, held constant outside
// the data.
template <typename X, typename Y>
double PCHIP(std::size_t n, const X &x, const Y &y, double xq, bool derivative)
{
  if (xq < x(0) || xq > x(n - 1))
  {
    return derivative ? 0.0 : y((xq < x(0)) ? 0 : n - 1);
  }
  std::size_t i = 0, j = n - 1;
  while (j - i > 1)
  {
    const std::size_t m = (i + j) / 2;
    (x(m) <= xq ? i : j) = m;
  }
  const double h = x(i + 1) - x(i), s = (xq - x(i)) / h;
  const double d0 = PCHIPSlope(n, i, x, y), d1 = PCHIPSlope(n, i + 1, x, y);
  if (derivative)
  {
    return 6.0 * s * (s - 1.0) / h * (y(i) - y(i + 1)) + (1.0 - s) * (1.0 - 3.0 * s) * d0 +
           s * (3.0 * s - 2.0) * d1;
  }
  return (1.0 + 2.0 * s) * (1.0 - s) * (1.0 - s) * y(i) +
         s * (1.0 - s) * (1.0 - s) * h * d0 + s * s * (3.0 - 2.0 * s) * y(i + 1) +
         s * s * (s - 1.0) * h * d1;
}

}  // namespace

double BCSGapRatio(double t)
{
  MFEM_VERIFY(t >= 0.0, "Reduced temperature must be non-negative!");
  if (t >= 1.0)
  {
    return 0.0;
  }
  if (t < 0.02)
  {
    return 1.0;  // Thermal corrections ~ exp(-BCS_GAP_RATIO / t) < 1e-38
  }
  // Root δ = Δ(T)/Δ(0) of G(δ) = -ln δ - 2 ∫₀^∞ f(Δ cosh u) du, decreasing on (0, 1].
  const double a = BCS_GAP_RATIO / t;
  auto G = [a](double delta)
  {
    const double b = a * delta;
    const double I =
        IntegrateGK([b](double u) { return 1.0 / (std::exp(b * std::cosh(u)) + 1.0); }, 0.0,
                    std::acosh(std::max(1.0, 50.0 / b)) + 1.0);
    return -std::log(delta) - 2.0 * I;
  };
  double lo = 1.0e-12, hi = 1.0;
  MFEM_VERIFY(G(lo) > 0.0 && G(hi) < 0.0, "Failed to bracket the BCS gap equation root!");
  while (hi - lo > 1.0e-15)
  {
    const double mid = 0.5 * (lo + hi);
    (G(mid) > 0.0 ? lo : hi) = mid;
  }
  return 0.5 * (lo + hi);
}

std::complex<double> MattisBardeenRatio(double x, double theta)
{
  MFEM_VERIFY(x > 0.0 && theta >= 0.0,
              "Mattis-Bardeen conductivity requires ħω/Δ > 0 and kB T/Δ >= 0!");
  // Energies e = E/Δ. The inverse square root endpoint singularities are removed with
  // e = cosh u on [1, ∞) and e = c + h cos φ on finite intervals.

  // Thermal term of σ1: (2/x) ∫₁^∞ [f(e) - f(e+x)] g(e) de, g(e) = (e² + 1 + xe) /
  // (√(e²-1) √((e+x)²-1)), truncated where f(e) < e^-50 f(1). The Fermi factors are scaled
  // by exp(1/θ) to avoid underflow.
  double A = 0.0;
  if (theta > 0.0)
  {
    const double c = std::exp(-1.0 / theta), em1 = -std::expm1(-x / theta);
    auto f = [x, theta, c, em1](double u)
    {
      const double e = std::cosh(u), a = std::exp(-(e - 1.0) / theta), b = a * (1.0 - em1);
      return a * em1 / ((1.0 + c * a) * (1.0 + c * b)) * (e * e + 1.0 + x * e) /
             std::sqrt((e + x - 1.0) * (e + x + 1.0));
    };
    A = c * IntegrateGK(f, 0.0, std::acosh(1.0 + 50.0 * theta));
  }

  // Pair-breaking term of σ1 for x > 2: (1/x) ∫_{1-x}^{-1} [1 - 2f(e+x)] |g(e)| de.
  double B = 0.0;
  if (x > 2.0)
  {
    const double c = -0.5 * x, h = 0.5 * x - 1.0;
    auto f = [x, theta, c, h](double phi)
    {
      const double e = c + h * std::cos(phi);
      return OneMinusTwoFermi(e + x, theta) * std::abs(e * e + 1.0 + x * e) /
             std::sqrt((1.0 - e) * (e + x + 1.0));
    };
    B = IntegrateGK(f, 0.0, std::numbers::pi);
  }

  // σ2: (1/x) ∫_a^1 [1 - 2f(e+x)] (e² + 1 + xe) / (√(1-e²) √((e+x)²-1)) de, with
  // a = max(1-x, -1).
  const bool below_gap = (x < 2.0);
  const double a = below_gap ? 1.0 - x : -1.0, c = 0.5 * (a + 1.0), h = 0.5 * (1.0 - a);
  auto f = [x, theta, c, h, below_gap](double phi)
  {
    const double e = c + h * std::cos(phi);
    const double r = below_gap ? std::sqrt((1.0 + e) * (e + x + 1.0))
                               : std::sqrt((e + x - 1.0) * (e + x + 1.0));
    return OneMinusTwoFermi(e + x, theta) * (e * e + 1.0 + x * e) / r;
  };
  const double C = IntegrateGK(f, 0.0, std::numbers::pi);

  return {(2.0 * A + B) / x, -C / x};
}

std::complex<double> Conductivity(const config::MattisBardeenData &data, double omega)
{
  MFEM_VERIFY(omega > 0.0, "Mattis-Bardeen conductivity requires a positive frequency!");
  return data.sigma_n *
         MattisBardeenRatio(omega / data.omega_gap, data.omega_T / data.omega_gap);
}

std::complex<double> ConductivityDerivative(const config::MattisBardeenData &data,
                                            double omega)
{
  // The step balances truncation and the 1e-13 quadrature tolerance.
  constexpr double h = 1.0e-4;
  return (Conductivity(data, omega * (1.0 + h)) - Conductivity(data, omega * (1.0 - h))) /
         (2.0 * h * omega);
}

namespace
{

// The interpolants Y₁ + iY₂ of σ₁ = Y₁(log ω) and σ₂ = Y₂(log ω)/ω, or their derivative
// in log ω.
std::complex<double> TableInterpolant(const config::TabulatedConductivityData &data,
                                      double omega, bool derivative)
{
  MFEM_VERIFY(omega > 0.0, "Tabulated conductivity requires a positive frequency!");
  const std::size_t n = data.omega.size();
  auto x = [&data](std::size_t k) { return std::log(data.omega[k]); };
  auto y1 = [&data](std::size_t k) { return data.sigma1[k]; };
  auto y2 = [&data](std::size_t k) { return data.omega[k] * data.sigma2[k]; };
  const double xq = std::log(omega);
  return {PCHIP(n, x, y1, xq, derivative), PCHIP(n, x, y2, xq, derivative)};
}

}  // namespace

std::complex<double> Conductivity(const config::TabulatedConductivityData &data,
                                  double omega)
{
  const auto Y = TableInterpolant(data, omega, false);
  return {Y.real(), -Y.imag() / omega};
}

std::complex<double> ConductivityDerivative(const config::TabulatedConductivityData &data,
                                            double omega)
{
  const auto Y = TableInterpolant(data, omega, false),
             dY = TableInterpolant(data, omega, true);
  return {dY.real() / omega, -(dY.imag() - Y.imag()) / (omega * omega)};
}

}  // namespace palace::conductivity
