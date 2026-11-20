// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "energybalance.hpp"

#include <algorithm>
#include <limits>
#include <mfem.hpp>

namespace palace
{

std::complex<double> QuadraticRootNearest(std::complex<double> m, std::complex<double> c,
                                          std::complex<double> k,
                                          std::complex<double> lambda0)
{
  if (m == 0.0)
  {
    return (c == 0.0) ? lambda0 : -k / c;
  }
  // The roots are q / m and k / q with q = -(c ± √(c² - 4 m k)) / 2, where the sign avoids
  // cancellation. For a weakly damped mode, both keep the loss rate -c / (2 m) accurate
  // to round-off relative to it, although it is much smaller than the frequency.
  const auto d = std::sqrt(c * c - 4.0 * m * k);
  const auto q = -0.5 * ((std::real(std::conj(c) * d) >= 0.0) ? c + d : c - d);
  if (q == 0.0)
  {
    return 0.0;  // Double root at zero (k = c = 0)
  }
  const auto l1 = q / m, l2 = k / q;
  return (std::abs(l1 - lambda0) <= std::abs(l2 - lambda0)) ? l1 : l2;
}

std::pair<std::complex<double>, int>
EnergyBalanceRoot(std::complex<double> k, std::complex<double> c, std::complex<double> m,
                  const std::function<std::complex<double>(std::complex<double>)> &a,
                  std::complex<double> lambda0, int max_it)
{
  // Fixed-point map: the root nearest to λ of the energy balance with the
  // frequency-dependent term frozen at λ. Its fixed points are roots of the energy balance,
  // and since it depends on λ only through a(λ), the secant method applied to
  // g(λ) = λ - F(λ) converges in a few iterations.
  int num_it = 0;
  auto F = [&](std::complex<double> l)
  {
    num_it++;
    return QuadraticRootNearest(m, c, k + a(l), l);
  };

  // The iteration has converged when the step changes the loss rate Re{λ} by a small
  // fraction of it, or λ by round-off (for a mode without losses).
  auto Converged = [](std::complex<double> l, std::complex<double> dl)
  {
    constexpr double eps = std::numeric_limits<double>::epsilon();
    return std::abs(dl) <= 1.0e-10 * std::abs(l.real()) + 1.0e2 * eps * std::abs(l);
  };
  std::complex<double> l0 = lambda0, g0 = l0 - F(l0), l1 = l0 - g0;
  if (Converged(l1, g0))
  {
    return {l1, num_it};
  }
  while (num_it < max_it)
  {
    const auto f1 = F(l1), g1 = l1 - f1;
    if (Converged(f1, g1))
    {
      return {f1, num_it};
    }
    const auto dg = g1 - g0;
    const auto l2 = (dg != 0.0) ? l1 - g1 * (l1 - l0) / dg : f1;
    l0 = l1;
    g0 = g1;
    l1 = l2;
  }
  return {l1, -num_it};
}

namespace
{

// Energy form xᴴ A x, using y as workspace.
std::complex<double> EnergyForm(MPI_Comm comm, const ComplexVector &x,
                                const ComplexOperator &A, ComplexVector &y)
{
  A.Mult(x, y);
  return linalg::Dot(comm, y, x);
}

}  // namespace

std::complex<double> GetEnergyBalanceEigenvalue(
    MPI_Comm comm, const ComplexVector &x, std::complex<double> lambda,
    const ComplexOperator &K, const ComplexOperator *C, const ComplexOperator &M,
    const ExtraSystemMatrixFunction *funcA2, ComplexVector &y, int *num_it)
{
  const auto k = EnergyForm(comm, x, K, y), m = EnergyForm(comm, x, M, y);
  const std::complex<double> c = C ? EnergyForm(comm, x, *C, y) : 0.0;
  if (!funcA2 || !*funcA2)
  {
    if (num_it)
    {
      *num_it = 0;
    }
    return QuadraticRootNearest(m, c, k, lambda);
  }
  auto a = [&](std::complex<double> l) -> std::complex<double>
  {
    auto A2 = (*funcA2)(l);
    return A2 ? EnergyForm(comm, x, *A2, y) : 0.0;
  };
  const auto [lambda_eb, it] = EnergyBalanceRoot(k, c, m, a, lambda);
  if (num_it)
  {
    *num_it = it;
  }
  return lambda_eb;
}

std::complex<double> GetEnergyBalanceEigenvalue(MPI_Comm comm, const ComplexVector &x,
                                                const ComplexOperator &K,
                                                const ComplexOperator &M, ComplexVector &y)
{
  const auto k = EnergyForm(comm, x, K, y), m = EnergyForm(comm, x, M, y);
  MFEM_VERIFY(m != 0.0, "Zero mass energy form for the energy balance of an eigenvector!");
  return k / m;
}

double GetQualityFactor(std::complex<double> omega)
{
  return (omega == 0.0) ? mfem::infinity() : 0.5 * std::abs(omega) / std::abs(omega.imag());
}

double GetQualityFactorError(std::complex<double> omega, std::complex<double> omega_eb)
{
  // Loss rates 1 / Q = 2 |Im{ω}| / |ω|, below r_min (Q above about 4.5e12) at the level of
  // round-off.
  auto LossRate = [](std::complex<double> w)
  { return (w == 0.0) ? 0.0 : 2.0 * std::abs(w.imag()) / std::abs(w); };
  constexpr double r_min = 1.0e3 * std::numeric_limits<double>::epsilon();
  const double r = LossRate(omega), r_eb = LossRate(omega_eb);
  if (std::max(r, r_eb) <= r_min)
  {
    return 0.0;
  }
  return (r > 0.0) ? std::abs(r_eb - r) / r : mfem::infinity();  // |Q - Q_EB| / Q_EB
}

}  // namespace palace
