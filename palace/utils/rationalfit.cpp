// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "rationalfit.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <Eigen/Dense>
#include "utils/aaa.hpp"

namespace palace::utils
{

namespace
{

using Complex = std::complex<double>;

double Distance(Complex w, const std::vector<Complex> &points)
{
  double d = std::numeric_limits<double>::infinity();
  for (const auto p : points)
  {
    d = std::min(d, std::abs(w - p));
  }
  return d;
}

// The product of two polynomials (lowest degree first).
std::vector<Complex> Product(const std::vector<Complex> &a, const std::vector<Complex> &b)
{
  std::vector<Complex> c(a.size() + b.size() - 1, 0.0);
  for (std::size_t i = 0; i < a.size(); i++)
  {
    for (std::size_t j = 0; j < b.size(); j++)
    {
      c[i + j] += a[i] * b[j];
    }
  }
  return c;
}

// The roots of a monic polynomial (lowest degree first), from its companion matrix.
std::vector<Complex> Roots(const std::vector<Complex> &a)
{
  const int n = static_cast<int>(a.size()) - 1;
  Eigen::MatrixXcd companion = Eigen::MatrixXcd::Zero(n, n);
  for (int j = 0; j < n; j++)
  {
    companion(0, j) = -a[n - 1 - j] / a[n];
  }
  for (int i = 1; i < n; i++)
  {
    companion(i, i - 1) = 1.0;
  }
  const Eigen::ComplexEigenSolver<Eigen::MatrixXcd> roots(companion, false);
  return {roots.eigenvalues().begin(), roots.eigenvalues().end()};
}

}  // namespace

bool FitWindow::Contains(std::complex<double> omega) const
{
  constexpr double eps = 1.0e-12;
  const auto d = omega - omega_t;
  if (std::abs(d) > radius * (1.0 + eps))
  {
    return false;
  }
  const double phi = std::arg(d);
  return std::abs(d) <= eps * radius ||
         (phi >= -margin - eps && phi <= 0.5 * std::numbers::pi + margin + eps);
}

std::complex<double> ScalarFit::Eval(std::complex<double> l) const
{
  std::complex<double> f = poly[0] + l * (poly[1] + l * poly[2]);
  for (std::size_t j = 0; j < poles.size(); j++)
  {
    f += residues[j] / (l - poles[j]);
  }
  return f;
}

bool ScalarFit::IsValid(const FitWindow &window, std::complex<double> omega) const
{
  return window.Contains(omega) && (strip == 0.0 || Distance(omega, singular) >= strip);
}

ScalarFit FitScalar(const std::function<std::complex<double>(std::complex<double>)> &g,
                    std::complex<double> alpha, int m,
                    const std::vector<std::complex<double>> &singular,
                    const FitWindow &window, double tol, int max_poles)
{
  const Complex i(0.0, 1.0);
  const double R = window.radius;
  const double phi0 = -window.margin, phi1 = 0.5 * std::numbers::pi + window.margin;
  std::vector<double> strips;
  if (std::ranges::none_of(singular, [&](auto s) { return window.Contains(s); }))
  {
    strips.push_back(0.0);
  }
  for (int j = 10; j >= 3; j--)
  {
    strips.push_back(R / (1 << j));
  }
  ScalarFit fit;
  fit.singular = singular;
  for (const double d : strips)
  {
    // The points of the boundary of the window minus the strip, at most a quarter of their
    // distance to the singular points apart, and those of the boundary of the strip, from
    // the singular points (whose first are the branch points). With density 4 they are
    // the samples, and with density 32, offset by half a step and with a grid inside, the
    // check points, which the samples cannot alias.
    auto keep = [&](Complex w)
    { return window.Contains(w) && (d == 0.0 || Distance(w, singular) >= 0.999 * d); };
    auto Points = [&](int density, double offset, bool interior)
    {
      std::vector<Complex> points;
      auto add = [&](Complex w)
      {
        if (keep(w) && Distance(w, points) > 1.0e-6 * std::max(d, R / 4096.0))
        {
          points.push_back(w);
        }
      };
      auto march = [&](auto &&point, double length)
      {
        for (double s = 0.0; s < length;)
        {
          const double h = std::clamp(Distance(point(s), singular) / density,
                                      std::max(d, R / 4096.0) / density, R / (4 * density));
          add(point(s + offset * h));
          s += h;
        }
      };
      march([&](double s) { return window.omega_t + std::polar(s, phi0); }, R);
      march([&](double s) { return window.omega_t + std::polar(R, phi0 + s / R); },
            R * (phi1 - phi0));
      march([&](double s) { return window.omega_t + std::polar(R - s, phi1); }, R);
      if (d > 0.0)
      {
        std::vector<Complex> centers;
        for (const auto s : singular)
        {
          if (Distance(s, centers) >= d / density)
          {
            centers.push_back(s);
          }
        }
        for (const auto s : centers)
        {
          for (int k = 0; k < 3 * density; k++)
          {
            add(s + std::polar(d, 2.0 * std::numbers::pi * (k + offset) / (3 * density)));
          }
        }
      }
      for (int a = 0; interior && a < 8; a++)
      {
        for (int b = 1; b < 8; b++)
        {
          add(window.omega_t +
              std::polar(R * (a + 0.5) / 8.0, phi0 + (phi1 - phi0) * b / 8.0));
        }
      }
      return points;
    };
    const auto samples = Points(4, 0.0, false), checks = Points(32, 0.5, true);
    const auto n = static_cast<long>(samples.size());
    Eigen::VectorXcd z(n), F(n);
    for (long k = 0; k < n; k++)
    {
      z(k) = i * samples[k];
      F(k) = g(z(k));
    }

    // AAA on the samples, to a tenth of the tolerance, or with more poles to a hundredth
    // when the check fails.
    for (const double factor : {0.1, 0.01})
    {
      const auto aaa = RunAAA(z, F, factor * tol, static_cast<std::size_t>(max_poles + 1));
      if (!aaa.converged)
      {
        break;
      }

      // The poles of the fit, from the arrowhead pencil of its barycentric form (an
      // expansion into monomials is not stable for many support points), and its residues
      // and linear part g(λ) ≈ b₀ + b₁ λ + Σⱼ cⱼ / (λ - pⱼ) by least squares on the samples
      // (with the SVD: the columns of poles close to each other are close).
      std::vector<Complex> pg;
      {
        int na = static_cast<int>(aaa.zj.size()) + 1;
        Eigen::MatrixXcd E = Eigen::MatrixXcd::Zero(na, na),
                         Bm = Eigen::MatrixXcd::Identity(na, na);
        Bm(0, 0) = 0.0;
        for (int k = 1; k < na; k++)
        {
          E(0, k) = aaa.wj(k - 1);
          E(k, 0) = 1.0;
          E(k, k) = aaa.zj(k - 1);
        }
        char jobvl = 'N', jobvr = 'N';
        int lwork = 4 * na, ldv = 1, info = 0;
        std::vector<Complex> alpha_(na), beta_(na), work(lwork);
        std::vector<double> rwork(8 * na);
        Complex vl = 0.0, vr = 0.0;
        zggev_(&jobvl, &jobvr, &na, E.data(), &na, Bm.data(), &na, alpha_.data(),
               beta_.data(), &vl, &ldv, &vr, &ldv, work.data(), &lwork, rwork.data(),
               &info);
        if (info != 0)
        {
          continue;
        }
        for (int k = 0; k < na; k++)
        {
          if (std::abs(beta_[k]) > 1.0e-13 * std::abs(alpha_[k]))
          {
            pg.push_back(alpha_[k] / beta_[k]);
          }
        }
      }
      // Without the linear part for m = 2, where λ² b₁ λ would be cubic.
      const int np = (m == 2) ? 1 : 2;
      Eigen::MatrixXcd A(n, pg.size() + np);
      for (long k = 0; k < n; k++)
      {
        A(k, 0) = 1.0;
        if (np == 2)
        {
          A(k, 1) = z(k) / R;
        }
        for (std::size_t j = 0; j < pg.size(); j++)
        {
          A(k, j + np) = 1.0 / (z(k) - pg[j]);
        }
      }
      const Eigen::VectorXcd coef =
          A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(F);

      // The term α λᵐ g(λ): for m ≥ 0, λᵐ / (λ - p) = Σⱼ₌₀ᵐ⁻¹ pᵐ⁻¹⁻ʲ λʲ + pᵐ / (λ - p), and
      // 1 / (λ (λ - p)) = (1 / p) (1 / (λ - p) - 1 / λ).
      const Complex b[2] = {coef(0), (np == 2) ? coef(1) / R : 0.0};
      std::vector<Complex> poly(4, 0.0), poles, residues;
      bool valid = true;
      if (m == -1)
      {
        poles.push_back(0.0);
        residues.push_back(alpha * b[0]);
        poly[0] = alpha * b[1];
      }
      else
      {
        poly[m] = alpha * b[0];
        poly[m + 1] = alpha * b[1];
      }
      for (std::size_t j = 0; j < pg.size(); j++)
      {
        const Complex p = pg[j], c = coef(j + np);
        valid = valid && (m >= 0 || p != 0.0);
        for (int k = 0; k < m; k++)
        {
          poly[k] += alpha * c * std::pow(p, m - 1 - k);
        }
        poles.push_back(p);
        residues.push_back(alpha * c * std::pow(p, m));
        if (m == -1)
        {
          residues[0] -= alpha * c / p;
        }
        // A pole of the fit where it is valid makes its error there unbounded.
        valid = valid && !keep(p / i);
      }
      if (!valid)
      {
        continue;
      }
      ScalarFit candidate;
      candidate.poly = {poly[0], poly[1], poly[2]};
      candidate.poles = std::move(poles);
      candidate.residues = std::move(residues);
      double err = 0.0, scale = 0.0;
      for (const auto w : checks)
      {
        const Complex l = i * w, f = alpha * std::pow(l, m) * g(l);
        err = std::max(err, std::abs(candidate.Eval(l) - f));
        scale = std::max(scale, std::abs(f));
      }
      if (err > tol * scale)
      {
        continue;
      }
      candidate.success = true;
      candidate.num_samples = static_cast<int>(n);
      candidate.num_checks = static_cast<int>(checks.size());
      candidate.error = err / scale;
      candidate.strip = d;
      candidate.singular = singular;
      return candidate;
    }
  }
  return fit;
}

std::vector<std::complex<double>> SquareRootCut(const std::vector<std::complex<double>> &p,
                                                const std::vector<std::complex<double>> &r,
                                                double spacing)
{
  const Complex i(0.0, 1.0);
  std::vector<Complex> D = {1.0}, P = {0.0};
  for (std::size_t k = 0; k < p.size(); k++)
  {
    std::vector<Complex> others = {r[k]};
    for (std::size_t j = 0; j < p.size(); j++)
    {
      if (j != k)
      {
        others = Product(others, {-p[j], 1.0});
      }
    }
    P.resize(std::max(P.size(), others.size()), 0.0);
    for (std::size_t j = 0; j < others.size(); j++)
    {
      P[j] += others[j];
    }
    D = Product(D, {-p[k], 1.0});
  }
  auto RootsAt = [&](double u)
  {
    if (u == 0.0)
    {
      return p;
    }
    auto a = D;
    for (std::size_t j = 0; j < P.size(); j++)
    {
      a[j] += u * P[j];
    }
    return Roots(a);
  };
  // The largest distance from a root of one set to the nearest of the other.
  auto Gap = [](const std::vector<Complex> &a, const std::vector<Complex> &b)
  {
    double gap = 0.0;
    for (const auto x : a)
    {
      gap = std::max(gap, Distance(x, b));
    }
    for (const auto x : b)
    {
      gap = std::max(gap, Distance(x, a));
    }
    return gap;
  };
  std::vector<Complex> cut;
  auto Add = [&](const std::vector<Complex> &roots)
  {
    for (const auto x : roots)
    {
      cut.push_back(x / i);
    }
  };
  // Bisect in u until the roots at the ends of each interval are close (depth limited).
  auto Trace = [&](auto &&self, double u0, const std::vector<Complex> &r0, double u1,
                   const std::vector<Complex> &r1, int depth) -> void
  {
    if (depth == 0 || Gap(r0, r1) <= spacing)
    {
      return;
    }
    const double u = 0.5 * (u0 + u1);
    const auto ru = RootsAt(u);
    Add(ru);
    self(self, u0, r0, u, ru, depth - 1);
    self(self, u, ru, u1, r1, depth - 1);
  };
  // The branch points first: the poles (u = 0) and the zeros (u = 1).
  const auto r_end = RootsAt(1.0);
  std::vector<Complex> r_prev = RootsAt(0.0);
  Add(r_prev);
  Add(r_end);
  for (int k = 1; k <= 16; k++)
  {
    const auto r_next = (k < 16) ? RootsAt(k / 16.0) : r_end;
    if (k < 16)
    {
      Add(r_next);
    }
    Trace(Trace, (k - 1) / 16.0, r_prev, k / 16.0, r_next, 24);
    r_prev = r_next;
  }
  return cut;
}

}  // namespace palace::utils
