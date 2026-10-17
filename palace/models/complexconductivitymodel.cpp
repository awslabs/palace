// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "complexconductivitymodel.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <queue>
#include <mfem.hpp>
#include "utils/configfile.hpp"

namespace palace
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

// Monotone cubic Hermite (Fritsch-Carlson) slopes for data y at strictly increasing x.
std::vector<double> PCHIPSlopes(const std::vector<double> &x, const std::vector<double> &y)
{
  const std::size_t n = x.size();
  std::vector<double> d(n, 0.0), delta(n - 1);
  for (std::size_t i = 0; i + 1 < n; i++)
  {
    delta[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
  }
  if (n == 2)
  {
    d[0] = d[1] = delta[0];
    return d;
  }
  for (std::size_t i = 1; i + 1 < n; i++)
  {
    if (delta[i - 1] * delta[i] > 0.0)
    {
      const double w1 = 2.0 * (x[i + 1] - x[i]) + (x[i] - x[i - 1]),
                   w2 = (x[i + 1] - x[i]) + 2.0 * (x[i] - x[i - 1]);
      d[i] = (w1 + w2) / (w1 / delta[i - 1] + w2 / delta[i]);
    }
  }
  auto end_slope = [](double h0, double h1, double del0, double del1)
  {
    const double s = ((2.0 * h0 + h1) * del0 - h0 * del1) / (h0 + h1);
    if (s * del0 <= 0.0)
    {
      return 0.0;
    }
    return (del0 * del1 <= 0.0 && std::abs(s) > std::abs(3.0 * del0)) ? 3.0 * del0 : s;
  };
  d[0] = end_slope(x[1] - x[0], x[2] - x[1], delta[0], delta[1]);
  d[n - 1] =
      end_slope(x[n - 1] - x[n - 2], x[n - 2] - x[n - 3], delta[n - 2], delta[n - 3]);
  return d;
}

}  // namespace

namespace superconductor
{

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

}  // namespace superconductor

ComplexConductivityModel::ComplexConductivityModel(
    const config::ComplexConductivityData &data)
  : type(data.type)
{
  switch (type)
  {
    case ComplexConductivityType::MATTIS_BARDEEN:
      sigma_n = data.sigma_n;
      omega_T = data.temperature;
      omega_gap = (data.gap > 0.0)
                      ? data.gap
                      : superconductor::BCS_GAP_RATIO * data.critical_temperature *
                            superconductor::BCSGapRatio(data.temperature /
                                                        data.critical_temperature);
      MFEM_VERIFY(sigma_n > 0.0 && omega_gap > 0.0,
                  "Mattis-Bardeen conductivity requires a positive normal conductivity and "
                  "gap!");
      break;
    case ComplexConductivityType::TABLE:
      {
        const std::size_t n = data.table_f.size();
        MFEM_VERIFY(n >= 2 && data.table_sigma1.size() == n &&
                        data.table_sigma2.size() == n,
                    "Invalid complex conductivity table!");
        log_omega.resize(n);
        y1 = data.table_sigma1;
        y2.resize(n);
        for (std::size_t i = 0; i < n; i++)
        {
          log_omega[i] = std::log(data.table_f[i]);
          y2[i] = data.table_f[i] * data.table_sigma2[i];
        }
        d1 = PCHIPSlopes(log_omega, y1);
        d2 = PCHIPSlopes(log_omega, y2);
      }
      break;
    case ComplexConductivityType::TWO_FLUID:
      MFEM_ABORT("Two-fluid conductivity is resolved into the material conductivity and "
                 "London depth!");
      break;
  }
}

double ComplexConductivityModel::EvalTable(double log_w, const std::vector<double> &y,
                                           const std::vector<double> &d) const
{
  const std::size_t n = log_omega.size();
  if (log_w <= log_omega[0])
  {
    return y[0];
  }
  if (log_w >= log_omega[n - 1])
  {
    return y[n - 1];
  }
  const std::size_t i =
      std::upper_bound(log_omega.begin(), log_omega.end(), log_w) - log_omega.begin() - 1;
  const double h = log_omega[i + 1] - log_omega[i], s = (log_w - log_omega[i]) / h;
  return (1.0 + 2.0 * s) * (1.0 - s) * (1.0 - s) * y[i] +
         s * (1.0 - s) * (1.0 - s) * h * d[i] + s * s * (3.0 - 2.0 * s) * y[i + 1] +
         s * s * (s - 1.0) * h * d[i + 1];
}

std::complex<double> ComplexConductivityModel::operator()(double omega) const
{
  MFEM_VERIFY(omega > 0.0, "Complex conductivity requires a positive frequency!");
  if (type == ComplexConductivityType::MATTIS_BARDEEN)
  {
    return sigma_n *
           superconductor::MattisBardeenRatio(omega / omega_gap, omega_T / omega_gap);
  }
  const double log_w = std::log(omega);
  return {EvalTable(log_w, y1, d1), -EvalTable(log_w, y2, d2) / omega};
}

}  // namespace palace
