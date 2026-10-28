// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <random>
#include <vector>
#include <Eigen/Dense>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "utils/configfile.hpp"
#include "utils/permittivitymodel.hpp"

using namespace palace;
using permittivity::PreparedStorage;
using permittivity::RawStorageStatus;
using Terms = std::vector<config::PermittivityTermData>;
using C = std::complex<double>;

namespace
{

const std::vector<double> omegas = {0.0, 0.02, 0.2, 0.5, 0.9, 1.3, 3.0, 17.0};

// The PoleResidue (RAW) form of the terms.
Terms Raw(const Terms &terms)
{
  Terms raw;
  for (const auto &pole : permittivity::PoleTerms(terms))
  {
    raw.push_back(pole);
  }
  return raw;
}

Terms Poles(const std::vector<std::pair<std::complex<double>, std::complex<double>>> &poles)
{
  Terms raw;
  for (const auto &[p, r] : poles)
  {
    raw.push_back(config::PermittivityPoleData{p, r});
  }
  return raw;
}

double Energy(const PreparedStorage &storage, double omega)
{
  const auto w = storage.StoredEnergySusceptibility(omega);
  REQUIRE(w);
  return *w;
}

double RelativeError(double a, double b)
{
  return std::abs(a - b) / std::max(std::abs(b), 1.0);
}

// Independent certificate from the realization: it reproduces χ_raw, X is positive
// definite, the positive-real (KYP) matrix is [L W]ᵀ[L W] ⪰ 0, and the energy is κ z*Xz.
void CheckCertificate(const Terms &raw, const PreparedStorage &storage, bool strict = true)
{
  REQUIRE(storage.Status() == RawStorageStatus::CheckedAvailable);
  const auto r = storage.Realization();
  const auto &cert = storage.Certificate();
  const int n = r.n;
  CHECK(n == cert.states);
  const double rate = cert.rate_scale, amp = cert.amplitude_scale;
  Eigen::MatrixXd A = Eigen::Map<const Eigen::MatrixXd>(r.A.data(), n, n).transpose();
  Eigen::MatrixXd X = Eigen::Map<const Eigen::MatrixXd>(r.X.data(), n, n).transpose();
  const Eigen::VectorXd B = Eigen::Map<const Eigen::VectorXd>(r.B.data(), n),
                        C = Eigen::Map<const Eigen::VectorXd>(r.C.data(), n),
                        L = Eigen::Map<const Eigen::VectorXd>(r.L.data(), n);
  for (const std::complex<double> s : {std::complex<double>(0.0, 0.37 * rate),
                                       {0.0, 2.9 * rate},
                                       {0.4 * rate, 1.1 * rate}})
  {
    const std::complex<double> u = s / rate;
    std::complex<double> y = r.D;
    if (n > 0)
    {
      const Eigen::VectorXcd z =
          (u * Eigen::MatrixXcd::Identity(n, n) - A.cast<std::complex<double>>())
              .partialPivLu()
              .solve(B.cast<std::complex<double>>());
      y += C.cast<std::complex<double>>().dot(z);
    }
    for (const auto &[beta, strength] : r.lossless)
    {
      y += strength * u / (u * u + beta * beta);
    }
    const auto chi = permittivity::Susceptibility(raw, s);
    CAPTURE(s, chi, amp * y / u);
    CHECK(std::abs(amp * y / u - chi) <= 1.0e-12 * std::abs(chi));
  }
  if (n == 0)
  {
    return;
  }
  // A strictly positive factor, unless nearly coincident poles make X nearly singular.
  if (strict)
  {
    CHECK(cert.zeroed_pivots == 0);
    CHECK(cert.min_pivot > 0.0);
    CHECK(Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(X).eigenvalues().minCoeff() > 0.0);
  }
  CHECK(Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(X).eigenvalues().minCoeff() >=
        -1.0e-14 * X.norm());
  for (double residual :
       {cert.zero_backward_error, cert.reality_error, cert.lyapunov_residual,
        cert.cross_residual, cert.feedthrough_residual, cert.factor_residual})
  {
    CHECK(residual <= PreparedStorage::tolerance);
  }
  Eigen::MatrixXd M(n + 1, n + 1);
  M.topLeftCorner(n, n) = -A.transpose() * X - X * A;
  M.topRightCorner(n, 1) = C - X * B;
  M.bottomLeftCorner(1, n) = (C - X * B).transpose();
  M(n, n) = 2.0 * r.D;
  Eigen::VectorXd LW(n + 1);
  LW << L, r.W;
  // Normwise, within the blockwise backward error the kernel accepts.
  const double scale = (n + 1) * PreparedStorage::tolerance *
                       (2.0 * A.norm() * X.norm() + LW.squaredNorm() + C.norm());
  CHECK((M - LW * LW.transpose()).norm() <= scale);
  CHECK(Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(M).eigenvalues().minCoeff() >=
        -scale);
  for (double omega : {0.1 * rate, 1.0 * rate, 10.0 * rate})
  {
    const std::complex<double> u(0.0, omega / rate);
    const Eigen::VectorXcd z =
        (u * Eigen::MatrixXcd::Identity(n, n) - A.cast<std::complex<double>>())
            .partialPivLu()
            .solve(B.cast<std::complex<double>>());
    double w = (z.adjoint() * X.cast<std::complex<double>>() * z)(0, 0).real();
    for (const auto &[beta, strength] : r.lossless)
    {
      const double v = omega / rate;
      w += strength * (v * v + beta * beta) / std::pow(beta * beta - v * v, 2);
    }
    CHECK_THAT(amp * w, Catch::Matchers::WithinRel(
                            *storage.RawStoredEnergySusceptibility(omega), 1.0e-12));
  }
}

// Scale every rate by a: χ_a(s) = χ(s/a), so w_a(aω) = w(ω).
Terms ScaleRates(const Terms &raw, double a)
{
  Terms scaled;
  for (const auto &term : raw)
  {
    auto t = std::get<config::PermittivityPoleData>(term);
    scaled.push_back(config::PermittivityPoleData{a * t.pole, a * t.residue});
  }
  return scaled;
}

// Reverse the terms and split each residue r at the same pole into 2r, −r/2, −r/2, and
// into c, a, −c, r − a with a = 0.7r, c = 0.123r (r − a is exact, Sterbenz): the exact sum
// is r, but summing in this order rounds.
Terms SplitReversed(const Terms &raw)
{
  Terms split;
  for (auto it = raw.rbegin(); it != raw.rend(); ++it)
  {
    const auto t = std::get<config::PermittivityPoleData>(*it);
    const std::complex<double> a = 0.7 * t.residue, c = 0.123 * t.residue;
    for (const auto r : {-0.5 * t.residue, 2.0 * t.residue, -0.5 * t.residue})
    {
      split.push_back(config::PermittivityPoleData{t.pole, 0.5 * r});
    }
    for (const auto r : {c, a, -c, t.residue - a})
    {
      split.push_back(config::PermittivityPoleData{t.pole, 0.5 * r});
    }
  }
  return split;
}

// A RAW block is certified, independently checked, invariant to order and split of its
// residues (exactly), and to a common rate scale, and matches an oracle if given.
template <typename Oracle>
void CheckRaw(const Terms &raw, Oracle &&oracle)
{
  const PreparedStorage storage(raw);
  CAPTURE(storage.Reason());
  CheckCertificate(raw, storage);
  const PreparedStorage split(SplitReversed(raw));
  REQUIRE(split.Status() == RawStorageStatus::CheckedAvailable);
  for (double omega : omegas)
  {
    const double w = Energy(storage, omega);
    CAPTURE(omega, w);
    CHECK(w >= 0.0);
    CHECK(Energy(split, omega) == w);
    if (const auto ref = oracle(omega))
    {
      CHECK(RelativeError(w, *ref) <= 1.0e-12);
    }
  }
  for (double a : {1.0e-12, 1.0e12})
  {
    const PreparedStorage scaled(ScaleRates(raw, a));
    CAPTURE(a, scaled.Reason());
    REQUIRE(scaled.Status() == RawStorageStatus::CheckedAvailable);
    for (double omega : omegas)
    {
      CHECK(RelativeError(Energy(scaled, a * omega), Energy(storage, omega)) <= 1.0e-12);
    }
  }
}

auto NoOracle = [](double) -> std::optional<double> { return std::nullopt; };

}  // namespace

TEST_CASE("RAW storage of a single typed model is its physical energy",
          "[permittivitystorage][Serial]")
{
  // A single Drude, Debye, or Lorentz section has one storage, so the RAW available storage
  // of its pole-residue form is the energy of the typed model.
  for (const Terms &typed :
       {Terms{config::DrudeData{1.0, 0.1}}, Terms{config::DrudeData{3.0, 0.2}},
        Terms{config::DebyeData{1.5, 0.4}}, Terms{config::LorentzData{0.8, 2.5, 0.3}},
        Terms{config::LorentzData{0.8, 2.5, 7.0}},
        Terms{config::LorentzData{0.8, 2.5, 0.0}},
        Terms{config::LorentzData{1.0, 1.0, 1.0e-10}}})
  {
    const PreparedStorage physical(typed);
    CHECK_FALSE(physical.HasRaw());
    CHECK(physical.Status() == RawStorageStatus::CheckedAvailable);
    CheckRaw(Raw(typed),
             [&](double omega) -> std::optional<double>
             {
               // Away from an exact lossless resonance.
               const double w = permittivity::StoredEnergySusceptibility(typed, omega);
               CHECK(Energy(physical, omega) == w);
               return w;
             });
  }

  // Raw Drude ω_p = 1, γ = 0.1, where the sum of the termwise energies of the two poles
  // is −3.846 at ω = 0.5.
  const PreparedStorage drude(Poles({{0.0, 10.0}, {-0.1, -10.0}}));
  CHECK_THAT(Energy(drude, 0.5), Catch::Matchers::WithinRel(1.0 / 0.26, 1.0e-14));
  CHECK(drude.Certificate().states == 1);
  CHECK(drude.Certificate().spectral_zeros == 0);
}

TEST_CASE("RAW storage is the available storage of the whole block",
          "[permittivitystorage][Serial]")
{
  SECTION("Two Debye relaxations")
  {
    // G(s) = 2s(s + √(5/2))/((s + 1)(s + 2)): the RAW storage is not the 0.9 of two typed
    // Debye terms with the same χ.
    const auto raw = Poles({{-1.0, 1.0}, {-2.0, 1.0}});
    CheckRaw(raw, NoOracle);
    const PreparedStorage storage(raw);
    CHECK_THAT(Energy(storage, 1.0),
               Catch::Matchers::WithinRel(0.8675444679663243, 1.0e-13));
    const Terms typed = {config::DebyeData{1.0, 1.0}, config::DebyeData{0.5, 0.5}};
    CHECK_THAT(Energy(PreparedStorage(typed), 1.0),
               Catch::Matchers::WithinRel(0.9, 1.0e-14));
  }

  SECTION("Four real poles with negative residues")
  {
    // No pairing of poles: D = 0 and D > 0 aggregates.
    const auto d0 = Poles({{-1.0, 1.0}, {-2.0, 1.0}, {-3.0, -1.0}, {-4.0, -1.0}});
    CheckRaw(d0, NoOracle);
    CHECK_THAT(Energy(PreparedStorage(d0), 0.5),
               Catch::Matchers::WithinRel(0.7814457115712845, 1.0e-12));
    const auto dp = Poles({{-1.0, 2.0}, {-2.0, -0.4}, {-3.0, 0.7}, {-4.0, 0.3}});
    CheckRaw(dp, NoOracle);
    CHECK_THAT(Energy(PreparedStorage(dp), 0.5),
               Catch::Matchers::WithinRel(1.6954993541048196, 1.0e-12));
  }

  SECTION("Real and conjugate poles")
  {
    auto raw = Raw({config::LorentzData{1.0, 2.0, 0.3}});
    raw.push_back(config::PermittivityPoleData{-1.0, 2.0});
    raw.push_back(config::PermittivityPoleData{-3.0, -0.3});
    CheckRaw(raw, NoOracle);
  }

  SECTION("Exact lossless resonance with a stable pole")
  {
    auto raw = Raw({config::LorentzData{1.0, 1.0, 0.0}});
    raw.push_back(config::PermittivityPoleData{-2.0, 3.0});
    CheckRaw(raw, [](double w) -> std::optional<double>
             { return (w * w + 1.0) / std::pow(1.0 - w * w, 2) + 6.0 / (w * w + 4.0); });
    CHECK(PreparedStorage(raw).Certificate().lossless == 1);
  }

  SECTION("Zero pole feedthrough")
  {
    // The zero pole is a conductance in the current, which removes its integrator but not
    // its contribution D.
    CheckRaw(Poles({{0.0, 3.0}, {-2.0, 1.0}}), [](double w) -> std::optional<double>
             { return 2.0 * std::pow(2.0 - std::sqrt(3.0), 2) / (w * w + 4.0); });
    const PreparedStorage conductance(Poles({{0.0, 3.0}}));
    CHECK(conductance.Status() == RawStorageStatus::CheckedAvailable);
    CHECK(conductance.Certificate().states == 0);
    CHECK(Energy(conductance, 0.7) == 0.0);
  }

  SECTION("Spectral zero on the imaginary axis")
  {
    // G = (s² + 1)/((s + 1)(s + 2)) has a double zero of F at s = ±i.
    CheckRaw(Poles({{0.0, 1.0 / 8.0}, {-1.0, -2.0 / 3.0}, {-2.0, 25.0 / 24.0}}),
             [](double w) -> std::optional<double>
             { return (19.0 * w * w + 11.0) / (12.0 * (w * w + 1.0) * (w * w + 4.0)); });
  }
}

TEST_CASE("RAW storage of passive composites", "[permittivitystorage][Serial]")
{
  SECTION("Drude terms sharing the zero pole")
  {
    // The zero-pole residues 45 + 1/0.7 round when aggregated, but D = 0 is exact for the
    // input: D is summed from the input residues.
    const auto raw = Raw({config::DrudeData{3.0, 0.2}, config::DrudeData{1.0, 0.7}});
    CheckRaw(raw, NoOracle);
    CHECK(PreparedStorage(raw).Certificate().feedthrough == 0.0);
  }

  SECTION("Weakly damped resonances")
  {
    // The two complex terms of each pair cancel to a value proportional to its damping,
    // and the real coordinates of a pair are one mode for the residual scale.
    CheckRaw(Poles({{0.0, 307.86387837760401},
                    {-0.015137046177438112, -307.86387837760401},
                    {C(-1.3952386676258726e-06, 0.032587360311030683),
                     C(0.0, -0.0021542149678251761)}}),
             NoOracle);
    CheckRaw(
        Poles({{C(-5.718340278900861e-07, 54.21075909130127), C(0.0, -54.634319983629453)},
               {C(-1.7759425872441289e-06, 0.020296880498810295),
                C(0.0, -0.011087688304958456)},
               {C(-6.5283770158256266e-05, 5.3563390737354206),
                C(0.0, -3.6186345068191708)}}),
        NoOracle);
  }

  SECTION("Random typed composites")
  {
    // The typed physical energy is a storage of the same χ, so the available storage of its
    // RAW form is no larger. Rates span six decades, with Lorentz damping down to 1e-6 of
    // the resonance.
    std::mt19937 gen(12345);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    auto log_uniform = [&](double lo, double hi)
    { return lo * std::pow(hi / lo, uniform(gen)); };
    int checked = 0;
    constexpr int trials = 200;
    for (int trial = 0; trial < trials; trial++)
    {
      Terms typed;
      for (int k = 0; k <= trial % 10; k++)
      {
        const double kind = uniform(gen), delta = log_uniform(0.1, 3.0);
        if (kind < 0.35)
        {
          typed.push_back(config::DebyeData{delta, 1.0 / log_uniform(1.0e-3, 1.0e3)});
        }
        else if (kind < 0.8)
        {
          const double w0 = log_uniform(1.0e-2, 1.0e2),
                       gamma = (kind < 0.6) ? log_uniform(1.0e-6, 1.9) * w0
                                            : log_uniform(2.1, 50.0) * w0;
          typed.push_back(config::LorentzData{delta, w0, gamma});
        }
        else
        {
          typed.push_back(
              config::DrudeData{log_uniform(0.1, 10.0), log_uniform(1.0e-2, 1.0e2)});
        }
      }
      const PreparedStorage storage(Raw(typed));
      CAPTURE(trial, storage.Reason());
      CHECK(storage.Status() != RawStorageStatus::NonpassiveRaw);
      if (storage.Status() != RawStorageStatus::CheckedAvailable)
      {
        continue;
      }
      checked++;
      for (int j = 0; j < 100; j++)
      {
        const double omega = 1.0e-5 * std::pow(10.0, 10.0 * j / 99.0);
        const double w = Energy(storage, omega),
                     physical = permittivity::StoredEnergySusceptibility(typed, omega);
        CAPTURE(omega, w, physical);
        CHECK(w >= 0.0);
        CHECK(w <= physical * (1.0 + 1.0e-9));
      }
    }
    CHECK(checked >= 0.95 * trials);
  }
}

TEST_CASE("RAW storage is minimal: KYP alone does not select it",
          "[permittivitystorage][Serial]")
{
  // χ = 1/(s + 1) + 1/(s + 2): Y = s χ = D + Σ qⱼ/(s − pⱼ), D = 2, q = (−1, −2), A =
  // diag(−1, −2), B = 1. The outer factor G = W + Σ lⱼ/(s − pⱼ) has W = 2, zero −√(5/2)
  // (and 0); reflecting that zero into the right half-plane gives another factor of the
  // same F. Both storages X = lᵢlⱼ/(aᵢ + aⱼ) satisfy every positive-real equation.
  const double a[2] = {1.0, 2.0}, cq[2] = {-1.0, -2.0}, W = 2.0;
  std::vector<double> frequencies = omegas;
  frequencies.push_back(1.0);
  std::vector<std::vector<double>> energies;
  for (const double k : {std::sqrt(2.5), -std::sqrt(2.5)})
  {
    // G = 2s(s + k)/((s + 1)(s + 2)), l₁ = G residue at −1, l₂ at −2.
    const double l[2] = {2.0 - 2.0 * k, 4.0 * k - 8.0};
    Eigen::Matrix2d X;
    for (int i = 0; i < 2; i++)
    {
      for (int j = 0; j < 2; j++)
      {
        X(i, j) = l[i] * l[j] / (a[i] + a[j]);
      }
    }
    for (int i = 0; i < 2; i++)
    {
      CHECK(std::abs(X.row(i).sum() - cq[i] + l[i] * W) <= 1.0e-14);  // XB = Cᵀ − LᵀW
    }
    CHECK(Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d>(X).eigenvalues().minCoeff() > 0.0);
    std::vector<double> w;
    for (double omega : frequencies)
    {
      const Eigen::Vector2cd z(1.0 / std::complex<double>(a[0], omega),
                               1.0 / std::complex<double>(a[1], omega));
      w.push_back((z.adjoint() * X.cast<std::complex<double>>() * z)(0, 0).real());
    }
    energies.push_back(w);
  }
  const PreparedStorage storage(Poles({{-1.0, 1.0}, {-2.0, 1.0}}));
  for (std::size_t i = 0; i < frequencies.size(); i++)
  {
    CAPTURE(frequencies[i]);
    CHECK(RelativeError(Energy(storage, frequencies[i]), energies[0][i]) <= 1.0e-13);
    CHECK(energies[1][i] >= energies[0][i]);  // Equal at ω = 0
  }
  CHECK_THAT(energies[0].back(), Catch::Matchers::WithinRel(0.8675444679663243, 1.0e-13));
  CHECK_THAT(energies[1].back(), Catch::Matchers::WithinRel(2.132455532033677, 1.0e-13));
}

TEST_CASE("RAW storage status of nonpassive and unresolved blocks",
          "[permittivitystorage][Serial]")
{
  auto Status = [](const Terms &raw)
  {
    const PreparedStorage storage(raw);
    UNSCOPED_INFO(storage.Reason());
    if (storage.Status() != RawStorageStatus::CheckedAvailable)
    {
      CHECK_FALSE(storage.StoredEnergySusceptibility(1.0));
      CHECK_FALSE(storage.RawStoredEnergySusceptibility(1.0));
      CHECK(storage.Realization().n == 0);
    }
    return storage.Status();
  };

  // Negative relaxation (D < 0), negative conductance, and invalid lossless residues.
  CHECK(Status(Poles({{-1.0, -0.5}})) == RawStorageStatus::NonpassiveRaw);
  CHECK(Status(Poles({{0.0, -1.0}, {-1.0, 2.0}})) == RawStorageStatus::NonpassiveRaw);
  CHECK(Status(Poles({{C(0.0, 1.0), C(0.0, 0.5)}})) == RawStorageStatus::NonpassiveRaw);
  CHECK(Status(Poles({{C(0.0, 1.0), C(0.1, -0.5)}})) == RawStorageStatus::NonpassiveRaw);

  // A narrow negative resonance on a Debye: D > 0 but gain near ω = 4, found by sampling.
  CHECK(Status(Poles({{-1.0, 1.0}, {C(-1.0e-3, 4.0), C(0.0, 1.0e-3)}})) ==
        RawStorageStatus::NonpassiveRaw);

  // A passive Debye plus a passive weakly damped Lorentz pair whose spectral density has no
  // resolved leading coefficient: unresolved, and the search for gain must not mistake the
  // rounding of an exactly zero coefficient (n₀ of a Lorentz pair) at ω = 1e-23 for gain.
  CHECK(Status(Poles({{-1.0e-20, 1.0e-60}, {C(-1.0e-18, 1.1), C(0.0, -0.7)}})) ==
        RawStorageStatus::NumericalUnresolved);

  // Outside the range of exact scaling: the damping of −2⁻¹⁰⁷⁴ + 2i would underflow to an
  // exact lossless pole when scaled, and the residues −1e308, −1e308, 1e308 of one lossless
  // pole have a finite sum but an overflowing intermediate one. Neither is checked, nor
  // reclassified as lossless, nor evaluated to NaN.
  for (const auto &raw :
       {Poles({{C(-std::numeric_limits<double>::denorm_min(), 2.0), C(0.0, -1.0)}}),
        Poles({{C(0.0, 1.0), C(0.0, -1.0e308)},
               {C(0.0, 1.0), C(0.0, -1.0e308)},
               {C(0.0, 1.0), C(0.0, 1.0e308)}})})
  {
    CHECK(Status(raw) == RawStorageStatus::NumericalUnresolved);
    CHECK(PreparedStorage(raw).Certificate().lossless == 0);
  }

  // Nearly coincident passive poles −(1 + kδ), residues 1, k = 0..3: passive, so never
  // NonpassiveRaw. Their storage tends to that of the merged Debye 4/(s + a) at the mean
  // a = 1 + 3δ/2, which stores 4a/(ω² + a²), as O(δ²). The
  // storage matrix is then nearly singular, and its factor drops trailing pivots within
  // rounding. At δ = 1e-9 the certificate fails, and the block is unresolved.
  for (const double delta : {1.0e-4, 1.0e-6, 1.0e-9})
  {
    const auto close = Poles({{-1.0, 1.0},
                              {-1.0 - delta, 1.0},
                              {-1.0 - 2.0 * delta, 1.0},
                              {-1.0 - 3.0 * delta, 1.0}});
    const PreparedStorage storage(close);
    CAPTURE(delta, storage.Reason());
    CHECK(Status(close) != RawStorageStatus::NonpassiveRaw);
    if (delta == 1.0e-9)
    {
      CHECK(storage.Status() == RawStorageStatus::NumericalUnresolved);
    }
    // At δ = 1e-6 the largest residual is 9.1e-11, just below the tolerance, so either
    // status is correct there; a certified result must still be accurate.
    if (delta == 1.0e-4 || storage.Status() == RawStorageStatus::CheckedAvailable)
    {
      CheckCertificate(close, storage, false);
      for (double omega : omegas)
      {
        const double a = 1.0 + 1.5 * delta;
        CHECK(RelativeError(Energy(storage, omega), 4.0 * a / (omega * omega + a * a)) <=
              10.0 * delta * delta);
      }
    }
  }
}

TEST_CASE("RAW storage is finite at extreme frequencies", "[permittivitystorage][Serial]")
{
  // A checked block evaluates to a finite value at any finite ω other than an exact
  // lossless resonance, underflowing to zero rather than to NaN through ∞/∞.
  auto Finite = [](const Terms &raw, double omega)
  {
    const PreparedStorage storage(raw);
    REQUIRE(storage.Status() == RawStorageStatus::CheckedAvailable);
    const double w = Energy(storage, omega);
    CAPTURE(omega, w);
    CHECK(std::isfinite(w));
    return w;
  };

  // Lossless (ω² + 1)/(1 − ω²)² of p = i, r = −i/2: about 1/ω².
  const auto lossless = Poles({{C(0.0, 1.0), C(0.0, -0.5)}});
  CHECK_THAT(Finite(lossless, 1.0e100), Catch::Matchers::WithinRel(1.0e-200, 1.0e-12));
  CHECK(Finite(lossless, 1.0e200) == 0.0);

  // Debye 6/(ω² + 4) of p = −2, r = 3.
  const auto debye = Poles({{-2.0, 3.0}});
  CHECK_THAT(Finite(debye, 1.0e150), Catch::Matchers::WithinRel(6.0e-300, 1.0e-12));
  CHECK(Finite(debye, 1.0e200) == 0.0);

  // Lorentz Δε ω₀²(ω² + ω₀²)/((ω₀² − ω²)² + γ²ω²), as ω₀²Δε(1/ω² + ω₀²/ω⁴)/((1 − ω₀²/ω²)²
  // + (γ/ω)²) for large ω (the typed formula itself overflows there).
  const double delta = 0.8, w0 = 2.5, gamma = 0.3;
  const auto lorentz = Raw({config::LorentzData{delta, w0, gamma}});
  for (const double omega : {1.0e100, 1.0e150})
  {
    const double x = w0 / omega;
    const double ref = delta * w0 * w0 *
                       (1.0 / (omega * omega) + x * x * x * x / (w0 * w0)) /
                       ((1.0 - x * x) * (1.0 - x * x) + (gamma / omega) * (gamma / omega));
    CHECK_THAT(Finite(lorentz, omega), Catch::Matchers::WithinRel(ref, 1.0e-12));
  }
  CHECK(Finite(lorentz, 1.0e200) == 0.0);

  // p = t(−1 + i), r = −itκ', χ = 2κ't²/(s² + 2ts + 2t²): the normalized model is that of
  // t = κ' = 1, w = 2κ'(v² + 2)/((2 − v²)² + 4v²) at v = ω/t, but v overflows at ω = 2⁷⁰⁰
  // for t = 2⁻⁷⁰⁰, where the energy rounds to zero. With r = −i2³²³ (κ' = 2¹⁰²³) at ω =
  // 2³³⁰, v = 2¹⁰³⁰ still overflows while the energy 2¹⁰²⁴/v²(1 + O(1/v²)) = 2⁻¹⁰³⁶ is a
  // representable (subnormal) number, which a normalized energy times κ would underflow to
  // zero.
  const double t = std::ldexp(1.0, -700), big = std::ldexp(1.0, 1023);
  auto Pair = [t](double kappa) { return Poles({{C(-t, t), C(0.0, -t * kappa)}}); };
  auto Exact = [](double kappa, double v)
  { return kappa * (2.0 * (v * v + 2.0) / ((2.0 - v * v) * (2.0 - v * v) + 4.0 * v * v)); };
  for (const double v : {0.0, 0.3, 1.0, 7.0, 1024.0})
  {
    for (const double kappa : {1.0, big})
    {
      CAPTURE(v, kappa);
      CHECK_THAT(Finite(Pair(kappa), v * t),
                 Catch::Matchers::WithinRel(Exact(kappa, v), 1.0e-12));
    }
  }
  CHECK(Finite(Pair(1.0), std::ldexp(1.0, 700)) == 0.0);
  // Subnormal: 2⁻¹⁰³⁶ has 38 significant bits.
  CHECK_THAT(Finite(Pair(big), std::ldexp(1.0, 330)),
             Catch::Matchers::WithinRel(std::ldexp(1.0, -1036), 1.0e-10));
}

TEST_CASE("A RAW block compensated by other terms is not material gain",
          "[permittivitystorage][Serial]")
{
  // RAW −0.5/(s + 1) has gain alone; with a typed Debye 1/(s + 1) or a native
  // Djordjevic-Sarkar 2 log((s + 2)/(s + 1)) the material is passive. The RAW block alone
  // has no available storage, and the whole-material gain search finds no gain.
  for (const auto &other :
       {config::PermittivityTermData{config::DebyeData{1.0, 1.0}},
        config::PermittivityTermData{config::DjordjevicSarkarData{2.0, 1.0, 2.0}}})
  {
    const Terms terms = {config::PermittivityPoleData{-1.0, -0.5}, other};
    const PreparedStorage storage(terms);
    CHECK(storage.HasRaw());
    CHECK(storage.Status() == RawStorageStatus::NonpassiveRaw);
    CHECK_FALSE(storage.StoredEnergySusceptibility(1.0));
    CHECK(permittivity::FindDefects(1.0, terms).gain.empty());
  }

  // A passive RAW block with typed and native terms adds their physical energies.
  const Terms typed = {config::DrudeData{3.0, 0.2}, config::LorentzData{0.8, 2.5, 0.3},
                       config::DjordjevicSarkarData{0.2, 0.01, 100.0}};
  Terms terms = typed;
  for (const auto &t : Poles({{-1.0, 1.0}, {-2.0, 1.0}}))
  {
    terms.push_back(t);
  }
  const PreparedStorage storage(terms);
  const PreparedStorage raw(Poles({{-1.0, 1.0}, {-2.0, 1.0}}));
  for (double omega : omegas)
  {
    CHECK_THAT(
        Energy(storage, omega),
        Catch::Matchers::WithinRel(permittivity::StoredEnergySusceptibility(typed, omega) +
                                       Energy(raw, omega),
                                   1.0e-15));
  }
}
