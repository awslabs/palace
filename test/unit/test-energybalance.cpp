// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <complex>
#include <functional>
#include <limits>
#include <memory>
#include <Eigen/Dense>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "linalg/energybalance.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

namespace palace
{

using namespace std::complex_literals;

namespace
{

// Dense complex matrix as a ComplexOperator, for small serial tests.
class DenseTestOperator : public ComplexOperator
{
  Eigen::MatrixXcd A;

public:
  DenseTestOperator(const Eigen::MatrixXcd &A) : ComplexOperator(A.rows(), A.cols()), A(A)
  {
  }

  void Mult(const ComplexVector &x, ComplexVector &y) const override
  {
    y.SetSize(height);
    const auto *XR = x.Real().HostRead(), *XI = x.Imag().HostRead();
    Eigen::VectorXcd v(width);
    for (int i = 0; i < width; i++)
    {
      v(i) = {XR[i], XI[i]};
    }
    const Eigen::VectorXcd w = A * v;
    auto *YR = y.Real().HostWrite(), *YI = y.Imag().HostWrite();
    for (int i = 0; i < height; i++)
    {
      YR[i] = w(i).real();
      YI[i] = w(i).imag();
    }
  }

  void AddMult(const ComplexVector &x, ComplexVector &y,
               const std::complex<double> a = 1.0) const override
  {
    ComplexVector t;
    Mult(x, t);
    y.Add(a, t);
  }
};

ComplexVector ToComplexVector(const Eigen::VectorXcd &v)
{
  ComplexVector x(v.size());
  auto *XR = x.Real().HostWrite(), *XI = x.Imag().HostWrite();
  for (int i = 0; i < v.size(); i++)
  {
    XR[i] = v(i).real();
    XI[i] = v(i).imag();
  }
  return x;
}

// Loss rate of the eigenvalue λ = iω, Re{λ} = -Im{ω}.
double LossRateError(std::complex<double> l, std::complex<double> l_ref)
{
  return std::abs(l.real() - l_ref.real()) / std::abs(l_ref.real());
}

// Weakly damped test matrices: K and M Hermitian positive definite, with a small complex
// symmetric (lossy) part of M, and a low-rank damping matrix C, so that the modes have
// Q ~ 1e6.
struct TestMatrices
{
  Eigen::MatrixXcd K, C, M;
};

TestMatrices MakeTestMatrices(int n)
{
  std::srand(1234);
  const Eigen::MatrixXcd R = Eigen::MatrixXcd::Random(n, n);
  const Eigen::MatrixXcd S = Eigen::MatrixXcd::Random(n, n);
  TestMatrices T;
  const Eigen::MatrixXd Sr = S.real();
  const Eigen::VectorXd u = Eigen::VectorXcd::Random(n).real();
  T.K = R * R.adjoint() + static_cast<double>(n) * Eigen::MatrixXcd::Identity(n, n);
  T.M = (Sr * Sr.transpose() + static_cast<double>(n) * Eigen::MatrixXd::Identity(n, n))
            .cast<std::complex<double>>();
  T.M -= std::complex<double>(0.0, 1.0e-7) * Eigen::MatrixXcd::Identity(n, n);
  T.C = (1.0e-6 * u * u.transpose()).cast<std::complex<double>>();
  return T;
}

// All eigenpairs of the quadratic eigenvalue problem (K + λ C + λ² M) x = 0, from its
// companion linearization.
Eigen::ComplexEigenSolver<Eigen::MatrixXcd> SolveQuadratic(const TestMatrices &T)
{
  const int n = T.K.rows();
  Eigen::MatrixXcd L = Eigen::MatrixXcd::Zero(2 * n, 2 * n);
  const Eigen::MatrixXcd Minv = T.M.inverse();
  L.topRightCorner(n, n) = Eigen::MatrixXcd::Identity(n, n);
  L.bottomLeftCorner(n, n) = -Minv * T.K;
  L.bottomRightCorner(n, n) = -Minv * T.C;
  return Eigen::ComplexEigenSolver<Eigen::MatrixXcd>(L);
}

}  // namespace

TEST_CASE("Energy balance roots", "[energybalance][Serial]")
{
  // Quadratic: m λ² + c λ + k = 0 with roots -γ ± iω₀, for a mode with Q = ω₀ / (2γ) = 5e9.
  // The loss rate γ is computed to round-off relative to it.
  {
    const double w0 = 2.0, g = 2.0e-10;
    const auto l = QuadraticRootNearest(1.0, 2.0 * g, w0 * w0 + g * g, 1i * w0);
    CHECK_THAT(l.imag(), Catch::Matchers::WithinRel(w0, 1.0e-15));
    CHECK_THAT(-l.real(), Catch::Matchers::WithinRel(g, 1.0e-12));
    const auto l_minus = QuadraticRootNearest(1.0, 2.0 * g, w0 * w0 + g * g, -1i * w0);
    CHECK_THAT(l_minus.imag(), Catch::Matchers::WithinRel(-w0, 1.0e-15));
    CHECK(QuadraticRootNearest(0.0, 2.0, 4.0, 1.0) == -2.0);
  }

  // Complex coefficients (a lossy mass term) and a root which solves the quadratic to
  // round-off.
  {
    const std::complex<double> m = 1.0 - 1.0e-6i, c = 3.0e-7, k = 4.0 + 1.0e-8i;
    const auto l = QuadraticRootNearest(m, c, k, 2.0i);
    CHECK(std::abs(m * l * l + c * l + k) <= 1.0e-15 * std::abs(k));
    CHECK(std::abs(l - 2.0i) < 1.0e-3);
  }

  // Frequency-dependent terms: a(λ) = β √λ (finite conductivity) or β / λ (pole of the
  // second-order absorbing boundary condition). The root is checked against Newton's method
  // with the exact derivative.
  const int type = GENERATE(0, 1, 2);
  const double beta = 1.0e-3;
  std::function<std::complex<double>(std::complex<double>)> a, da;
  switch (type)
  {
    case 0:
      a = [](std::complex<double>) { return 0.0; };
      da = [](std::complex<double>) { return 0.0; };
      break;
    case 1:
      a = [=](std::complex<double> l) { return beta * std::sqrt(l); };
      da = [=](std::complex<double> l) { return 0.5 * beta / std::sqrt(l); };
      break;
    case 2:
      a = [=](std::complex<double> l) { return beta / l; };
      da = [=](std::complex<double> l) { return -beta / (l * l); };
      break;
  }
  const std::complex<double> k = 4.0, c = 2.0e-9, m = 1.0;
  std::complex<double> l_ref = 2.0i;
  for (int it = 0; it < 50; it++)
  {
    l_ref -=
        (k + l_ref * c + l_ref * l_ref * m + a(l_ref)) / (c + 2.0 * l_ref * m + da(l_ref));
  }
  const auto [l, num_it] = EnergyBalanceRoot(k, c, m, a, 2.0i * (1.0 + 1.0e-6));
  CHECK(num_it > 0);
  CHECK(num_it <= 6);
  CHECK(std::abs(l - l_ref) <= 1.0e-14 * std::abs(l_ref));
  CHECK(LossRateError(l, l_ref) <= 1.0e-8);
}

TEST_CASE("Energy balance of eigenvectors", "[energybalance][Serial]")
{
  // Weakly damped quadratic eigenvalue problem with a lossy mass term: the energy balance
  // of an exact eigenvector gives the eigenvalue whatever the eigenvalue it is started
  // from, and the relative error of its quality factor is of the order of the error of the
  // eigenvector, while an error of the same relative size of the eigenvalue is amplified by
  // 2Q in its quality factor.
  constexpr int n = 12;
  const auto T = MakeTestMatrices(n);
  const auto es = SolveQuadratic(T);
  DenseTestOperator K(T.K), C(T.C), M(T.M);
  MPI_Comm comm = Mpi::World();
  ComplexVector y(n);
  int num_checked = 0;
  for (int j = 0; j < 2 * n; j++)
  {
    const std::complex<double> l = es.eigenvalues()(j);
    if (l.imag() <= 0.0)
    {
      continue;
    }
    const ComplexVector x = ToComplexVector(es.eigenvectors().col(j).head(n));
    const double Q = GetQualityFactor(l / 1i);
    REQUIRE(Q > 1.0e5);

    // Exact eigenpair, started from a perturbed eigenvalue.
    const auto l_eb =
        GetEnergyBalanceEigenvalue(comm, x, l * (1.0 + 1.0e-6), K, &C, M, nullptr, y);
    CHECK(std::abs(l_eb - l) <= 1.0e-13 * std::abs(l));
    CHECK(GetQualityFactorError(l / 1i, l_eb / 1i) <= 1.0e-7);

    // Perturbed eigenvector, x + ε r: a first-order error of the energy balance quality
    // factor, against an error amplified by 2Q for a much smaller relative perturbation of
    // the eigenvalue (the loss rate of the energy balance is only accurate to round-off
    // relative to the frequency, about 1e-9 relative to the loss rate here).
    constexpr double eps = 1.0e-6;
    ComplexVector r(n), xp(n);
    linalg::SetRandom(comm, r, j);
    xp = x;
    xp.Add(eps * linalg::Norml2(comm, x) / linalg::Norml2(comm, r), r);
    const auto lp_eb = GetEnergyBalanceEigenvalue(comm, xp, l, K, &C, M, nullptr, y);
    CHECK(GetQualityFactorError(l_eb / 1i, lp_eb / 1i) <= 1.0e2 * eps);
    const double eps_l = 1.0e-3 / Q;
    CHECK(GetQualityFactorError(l * (1.0 + eps_l * (1.0 + 1.0i)) / 1i, l_eb / 1i) >=
          Q * eps_l);
    num_checked++;
  }
  CHECK(num_checked == n);

  // Linear problem (K - μ M) x = 0, μ = ω², with a lossy mass term.
  {
    const Eigen::ComplexEigenSolver<Eigen::MatrixXcd> les(T.M.inverse() * T.K);
    for (int j = 0; j < n; j++)
    {
      const std::complex<double> mu = les.eigenvalues()(j);
      const ComplexVector x = ToComplexVector(les.eigenvectors().col(j));
      const auto mu_eb = GetEnergyBalanceEigenvalue(comm, x, K, M, y);
      CHECK(std::abs(mu_eb - mu) <= 1.0e-13 * std::abs(mu));
      CHECK(GetQualityFactorError(std::sqrt(mu), std::sqrt(mu_eb)) <= 1.0e-6);
    }
  }

  // Nonlinear problem (K + λ C + λ² M + β √λ B) x = 0 with diagonal matrices, whose
  // eigenvectors are the unit vectors: the energy balance with the frequency-dependent term
  // gives the eigenvalue from a perturbed start.
  {
    Eigen::MatrixXcd Kd = Eigen::MatrixXcd::Zero(n, n), Cd = Kd, Md = Kd, Bd = Kd;
    for (int i = 0; i < n; i++)
    {
      Kd(i, i) = 1.0 + i;
      Cd(i, i) = 1.0e-8 * (1.0 + 0.5 * i);
      Md(i, i) = 1.0 - 1.0e-9i;
      Bd(i, i) = 1.0e-4 * (2.0 + i);
    }
    DenseTestOperator Kop(Kd), Cop(Cd), Mop(Md);
    const ExtraSystemMatrixFunction funcA2 =
        [&Bd](std::complex<double> l) -> std::unique_ptr<ComplexOperator>
    { return std::make_unique<DenseTestOperator>(std::sqrt(l) * Bd); };
    for (int i = 0; i < n; i++)
    {
      Eigen::VectorXcd e = Eigen::VectorXcd::Zero(n);
      e(i) = 1.0;
      const ComplexVector x = ToComplexVector(e);
      std::complex<double> l_ref = 1i * std::sqrt(Kd(i, i));
      for (int it = 0; it < 50; it++)
      {
        const auto f = Kd(i, i) + l_ref * Cd(i, i) + l_ref * l_ref * Md(i, i) +
                       std::sqrt(l_ref) * Bd(i, i);
        const auto df =
            Cd(i, i) + 2.0 * l_ref * Md(i, i) + 0.5 * Bd(i, i) / std::sqrt(l_ref);
        l_ref -= f / df;
      }
      int num_it = 0;
      const auto l_eb = GetEnergyBalanceEigenvalue(comm, x, l_ref * (1.0 + 1.0e-5), Kop,
                                                   &Cop, Mop, &funcA2, y, &num_it);
      CHECK(num_it > 0);
      CHECK(std::abs(l_eb - l_ref) <= 1.0e-13 * std::abs(l_ref));
      CHECK(LossRateError(l_eb, l_ref) <= 1.0e-8);
    }
  }
}

TEST_CASE("Quality factor error", "[energybalance][Serial]")
{
  // Relative difference of the quality factors, |Q - Q_EB| / Q_EB.
  CHECK_THAT(GetQualityFactorError(1.0 + 1.0e-6i, 1.0 + 1.2e-6i),
             Catch::Matchers::WithinRel(0.2, 1.0e-9));
  CHECK_THAT(GetQualityFactor(1.0 + 1.0e-6i), Catch::Matchers::WithinRel(5.0e5, 1.0e-12));
  // Lossless modes to round-off: no difference, but an eigenvalue without loss for an
  // eigenvector which dissipates differs by an infinite factor.
  CHECK(GetQualityFactorError(1.0 + 1.0e-17i, 1.0) == 0.0);
  CHECK(std::isinf(GetQualityFactorError(1.0, 1.0 + 1.0e-6i)));
  CHECK(std::isinf(GetQualityFactor(2.0)));
}

}  // namespace palace
