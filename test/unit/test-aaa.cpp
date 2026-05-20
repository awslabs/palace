// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>
#include <Eigen/Dense>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "utils/aaa.hpp"

using namespace palace::utils;
using namespace Catch::Matchers;

namespace
{

// Sample F at omega ∈ [a, b] uniformly with N points.
Eigen::VectorXcd LinSpace(double a, double b, std::size_t N)
{
  Eigen::VectorXcd z(N);
  for (std::size_t i = 0; i < N; i++)
  {
    z(static_cast<long>(i)) =
        a + (b - a) * static_cast<double>(i) / static_cast<double>(N - 1);
  }
  return z;
}

}  // namespace

TEST_CASE("AAA: trivial constant function", "[aaa][Serial]")
{
  // f(z) = 7.0 — should converge with m=1 to machine precision.
  auto z = LinSpace(0.0, 10.0, 21);
  Eigen::VectorXcd F = Eigen::VectorXcd::Constant(z.size(), 7.0);
  auto r = RunAAA(z, F, 1e-12, 5);
  REQUIRE(r.converged);
  REQUIRE(r.zj.size() == 1);
  for (long i = 0; i < z.size(); i++)
  {
    auto v = EvaluateAAA(r, z(i));
    REQUIRE_THAT(std::abs(v - F(i)), WithinAbs(0.0, 1e-13));
  }
}

TEST_CASE("AAA: rational function recovered exactly", "[aaa][Serial]")
{
  // f(z) = 1/(z − 2) + 3/(z − 5). AAA should converge with m = 2 supports
  // (m+1 = 3 weights → degree-2 numerator / degree-2 denominator) up to roundoff.
  auto z = LinSpace(7.0, 30.0, 30);
  Eigen::VectorXcd F(z.size());
  for (long i = 0; i < z.size(); i++)
  {
    F(i) = 1.0 / (z(i) - 2.0) + 3.0 / (z(i) - 5.0);
  }
  auto r = RunAAA(z, F, 1e-13, 8);
  REQUIRE(r.converged);
  REQUIRE(r.zj.size() <= 4);  // typically 2-3 enough
  // Verify barycentric form on a fine grid (interior points).
  auto z_test = LinSpace(8.0, 28.0, 100);
  for (long i = 0; i < z_test.size(); i++)
  {
    auto truth = 1.0 / (z_test(i) - 2.0) + 3.0 / (z_test(i) - 5.0);
    auto v = EvaluateAAA(r, z_test(i));
    REQUIRE_THAT(std::abs(v - truth) / std::abs(truth), WithinAbs(0.0, 1e-9));
  }
  // Verify pole/residue extraction recovers the same function. With f(∞)=0 the
  // asymptote d should be very small; both poles and residues should be real.
  auto pr = AAAToPoleResidue(r);
  REQUIRE_THAT(std::abs(pr.d), WithinAbs(0.0, 1e-12));
  // Match poles to known set {2, 5} (order-independent).
  std::vector<std::complex<double>> sorted_poles(pr.poles.data(),
                                                 pr.poles.data() + pr.poles.size());
  std::sort(sorted_poles.begin(), sorted_poles.end(),
            [](auto a, auto b) { return a.real() < b.real(); });
  REQUIRE(sorted_poles.size() >= 2);
  REQUIRE_THAT(std::abs(sorted_poles[0] - 2.0) / 2.0, WithinAbs(0.0, 1e-9));
  REQUIRE_THAT(std::abs(sorted_poles[1] - 5.0) / 5.0, WithinAbs(0.0, 1e-9));
  for (long i = 0; i < z_test.size(); i++)
  {
    auto truth = 1.0 / (z_test(i) - 2.0) + 3.0 / (z_test(i) - 5.0);
    auto v = EvaluatePoleResidue(pr, z_test(i));
    REQUIRE_THAT(std::abs(v - truth) / std::abs(truth), WithinAbs(0.0, 1e-9));
  }
}

TEST_CASE("AAA: textbook waveguide dispersion residual converges", "[aaa][Serial]")
{
  // The polynomial-residual case from prom-waveport-validation: kₙ(ω) = √(γω² + α)
  // (textbook waveguide, γ = 1/c², α = -ω_c²/c²), polynomial-fitted by α₀+α₁ω+α₂ω²
  // on a dense grid in [ω_min, ω_max]. AAA fits the residual.
  //
  // Numbers approximating the adapter's rectangular port: ω_c/(2π) = 6.52 GHz,
  // c_eff = c₀ ≈ 2.998e8, sweep band 7.5 - 16 GHz.
  const double c_eff = 2.998e8;
  const double omega_c = 2.0 * M_PI * 6.5172e9;
  const double gamma_kn2 = 1.0 / (c_eff * c_eff);
  const double alpha_kn2 = -omega_c * omega_c * gamma_kn2;
  const double omega_lo = 2.0 * M_PI * 7.5e9;
  const double omega_hi = 2.0 * M_PI * 16.0e9;
  constexpr int n_fit = 30;
  Eigen::VectorXcd z_full(n_fit);
  Eigen::VectorXd kn_full(n_fit);
  for (int i = 0; i < n_fit; i++)
  {
    double w = omega_lo + (omega_hi - omega_lo) * i / (n_fit - 1.0);
    z_full(i) = w;
    kn_full(i) = std::sqrt(gamma_kn2 * w * w + alpha_kn2);
  }
  // Polynomial fit (degree 2 in ω).
  Eigen::MatrixXd vand(n_fit, 3);
  for (int i = 0; i < n_fit; i++)
  {
    double w = z_full(i).real();
    vand(i, 0) = 1.0;
    vand(i, 1) = w;
    vand(i, 2) = w * w;
  }
  Eigen::Vector3d coeffs = vand.colPivHouseholderQr().solve(kn_full);
  Eigen::VectorXcd delta(n_fit);
  for (int i = 0; i < n_fit; i++)
  {
    double w = z_full(i).real();
    delta(i) = kn_full(i) - (coeffs(0) + coeffs(1) * w + coeffs(2) * w * w);
  }
  const double delta_max = delta.cwiseAbs().maxCoeff();
  REQUIRE(delta_max >
          1e-3 * std::abs(coeffs(2) * omega_hi * omega_hi));  // nontrivial residual
  // AAA should hit 1e-3 relative on δkₙ with m ≤ 6.
  auto r = RunAAA(z_full, delta, 1e-6, 10);
  REQUIRE(r.converged);
  REQUIRE(r.zj.size() <= 8);
  // Spot-check on a finer grid (in band).
  auto z_test = LinSpace(omega_lo + 1e8, omega_hi - 1e8, 100);
  double max_rel = 0.0;
  for (long i = 0; i < z_test.size(); i++)
  {
    double w = z_test(i).real();
    double truth = std::sqrt(gamma_kn2 * w * w + alpha_kn2) -
                   (coeffs(0) + coeffs(1) * w + coeffs(2) * w * w);
    auto v = EvaluateAAA(r, z_test(i));
    max_rel = std::max(max_rel, std::abs(v.real() - truth) / std::abs(truth));
  }
  REQUIRE(max_rel < 1e-3);  // looser since this is OUT-of-sample evaluation
}
