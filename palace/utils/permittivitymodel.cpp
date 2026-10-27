// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "permittivitymodel.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <optional>
#include <utility>
#include <variant>
#include <mfem.hpp>

namespace palace::permittivity
{

namespace
{

template <class... Ts>
struct Overloaded : Ts...
{
  using Ts::operator()...;
};

std::complex<double> TermSusceptibility(const config::PermittivityTermData &term,
                                        std::complex<double> s)
{
  return std::visit(
      Overloaded{
          [s](const config::DrudeData &t) -> std::complex<double>
          { return t.omega_p * t.omega_p / (s * (s + t.gamma)); },
          [s](const config::DebyeData &t) -> std::complex<double>
          { return t.delta / (1.0 + s * t.tau); },
          [s](const config::LorentzData &t) -> std::complex<double>
          {
            const double w0_2 = t.omega_0 * t.omega_0;
            return t.delta * w0_2 / (s * s + t.gamma * s + w0_2);
          },
          [s](const config::PermittivityPoleData &t) -> std::complex<double>
          {
            std::complex<double> chi = t.residue / (s - t.pole);
            if (t.pole.imag() > 0.0)
            {
              chi += std::conj(t.residue) / (s - std::conj(t.pole));
            }
            return chi;
          },
          [s](const config::DjordjevicSarkarData &t) -> std::complex<double>
          { return t.strength * std::log((s + t.omega_upper) / (s + t.omega_lower)); }},
      term);
}

// Time-averaged stored energy of a term as a susceptibility (see the header): a relaxation
// b₀/(s + a₀) stores a₀b₀/(ω² + a₀²) = Re χ, and an oscillator (b₀ + b₁s)/(s² + a₁s + a₀)
// stores ((b₀ + a₁b₁)ω² + a₀b₀)/((a₀ − ω²)² + a₁²ω²), the averages of W with ⟨x ẋ⟩ = 0.
double TermStoredEnergy(const config::PermittivityTermData &term, double omega)
{
  const double w2 = omega * omega;
  auto oscillator = [w2](double a0, double a1, double b0, double b1)
  { return ((b0 + a1 * b1) * w2 + a0 * b0) / ((a0 - w2) * (a0 - w2) + a1 * a1 * w2); };
  return std::visit(Overloaded{[w2](const config::DrudeData &t)
                               { return t.omega_p * t.omega_p / (w2 + t.gamma * t.gamma); },
                               [w2](const config::DebyeData &t)
                               { return t.delta / (1.0 + w2 * t.tau * t.tau); },
                               [&](const config::LorentzData &t)
                               {
                                 const double w0_2 = t.omega_0 * t.omega_0;
                                 return oscillator(w0_2, t.gamma, t.delta * w0_2, 0.0);
                               },
                               [&](const config::PermittivityPoleData &t)
                               {
                                 if (t.pole.imag() > 0.0)
                                 {
                                   return oscillator(
                                       std::norm(t.pole), -2.0 * t.pole.real(),
                                       -2.0 * (t.residue * std::conj(t.pole)).real(),
                                       2.0 * t.residue.real());
                                 }
                                 const double a0 = -t.pole.real(), b0 = t.residue.real();
                                 return (a0 == 0.0) ? 0.0 : a0 * b0 / (w2 + a0 * a0);
                               },
                               [w2](const config::DjordjevicSarkarData &t)
                               {
                                 return 0.5 * t.strength *
                                        std::log((w2 + t.omega_upper * t.omega_upper) /
                                                 (w2 + t.omega_lower * t.omega_lower));
                               }},
                    term);
}

void AppendPole(std::complex<double> pole, std::complex<double> residue,
                std::vector<config::PermittivityPoleData> &poles)
{
  MFEM_VERIFY(std::isfinite(pole.real()) && std::isfinite(pole.imag()) &&
                  std::isfinite(residue.real()) && std::isfinite(residue.imag()),
              "Material permittivity poles and residues must be finite!");
  poles.push_back({pole, residue});
}

// The rates of a term, which set the frequency range of the defect search, and its
// resonances with their half-widths.
void TermRates(const config::PermittivityTermData &term, std::vector<double> &rates,
               std::vector<std::pair<double, double>> &resonances)
{
  std::visit(Overloaded{[&](const config::DrudeData &t)
                        {
                          rates.push_back(t.omega_p);
                          rates.push_back(t.gamma);
                        },
                        [&](const config::DebyeData &t) { rates.push_back(1.0 / t.tau); },
                        [&](const config::LorentzData &t)
                        {
                          rates.push_back(t.omega_0);
                          if (t.gamma > 0.0)
                          {
                            rates.push_back(t.gamma);
                          }
                          if (t.gamma < 2.0 * t.omega_0)
                          {
                            resonances.emplace_back(t.omega_0, 0.5 * t.gamma);
                          }
                        },
                        [&](const config::PermittivityPoleData &t)
                        {
                          for (double r : {std::abs(t.pole), -t.pole.real(), t.pole.imag(),
                                           std::abs(t.residue)})
                          {
                            if (r > 0.0)
                            {
                              rates.push_back(r);
                            }
                          }
                          if (t.pole.imag() > 0.0)
                          {
                            resonances.emplace_back(t.pole.imag(), -t.pole.real());
                          }
                        },
                        [&](const config::DjordjevicSarkarData &t)
                        {
                          rates.push_back(t.omega_lower);
                          rates.push_back(t.omega_upper);
                        }},
             term);
}

// Intervals of ω where bad(ω) holds, from samples on a sorted grid, with the endpoints
// bisected between neighboring samples. Samples where bad(ω) is undefined (a lossless
// resonance) are skipped.
std::vector<FrequencyInterval>
FindIntervals(const std::vector<double> &grid,
              const std::function<std::optional<bool>(double)> &bad)
{
  auto bisect = [&](double good, double bad_omega)
  {
    // Logarithmic bisection to a relative 1e-12.
    for (int it = 0; it < 200 && std::abs(bad_omega - good) > 1.0e-12 * bad_omega; it++)
    {
      const double mid = std::sqrt(good * bad_omega);
      const auto b = bad(mid);
      if (!b)
      {
        break;
      }
      (*b ? bad_omega : good) = mid;
    }
    return bad_omega;
  };
  std::vector<FrequencyInterval> intervals;
  std::optional<double> prev;
  bool prev_bad = false;
  for (double omega : grid)
  {
    const auto b = bad(omega);
    if (!b)
    {
      continue;
    }
    if (*b && !prev_bad)
    {
      intervals.push_back(
          {prev ? bisect(*prev, omega) : 0.0, std::numeric_limits<double>::infinity()});
    }
    else if (!*b && prev_bad)
    {
      intervals.back().upper = bisect(omega, *prev);
    }
    prev = omega;
    prev_bad = *b;
  }
  return intervals;
}

}  // namespace

std::vector<config::PermittivityPoleData>
PoleTerms(const std::vector<config::PermittivityTermData> &terms)
{
  std::vector<config::PermittivityPoleData> poles;
  for (const auto &term : terms)
  {
    std::visit(Overloaded{[&](const config::DrudeData &t)
                          {
                            const double residue = t.omega_p * t.omega_p / t.gamma;
                            AppendPole(0.0, residue, poles);
                            AppendPole(-t.gamma, -residue, poles);
                          },
                          [&](const config::DebyeData &t)
                          { AppendPole(-1.0 / t.tau, t.delta / t.tau, poles); },
                          [&](const config::LorentzData &t)
                          {
                            const double w0 = t.omega_0, gamma = t.gamma;
                            if (gamma < 2.0 * w0)
                            {
                              const double wd =
                                  std::sqrt((w0 - 0.5 * gamma) * (w0 + 0.5 * gamma));
                              AppendPole({-0.5 * gamma, wd},
                                         {0.0, -0.5 * t.delta * w0 * w0 / wd}, poles);
                            }
                            else
                            {
                              const double root =
                                  std::sqrt((gamma - 2.0 * w0) * (gamma + 2.0 * w0));
                              const double p_fast = -0.5 * (gamma + root);
                              const double p_slow = -2.0 * w0 * w0 / (gamma + root);
                              const double residue = t.delta * w0 * w0 / (p_slow - p_fast);
                              AppendPole(p_slow, residue, poles);
                              AppendPole(p_fast, -residue, poles);
                            }
                          },
                          [&](const config::PermittivityPoleData &t)
                          { AppendPole(t.pole, t.residue, poles); },
                          [](const config::DjordjevicSarkarData &) {}},
               term);
  }
  return poles;
}

std::vector<config::DjordjevicSarkarData>
DjordjevicSarkarTerms(const std::vector<config::PermittivityTermData> &terms)
{
  std::vector<config::DjordjevicSarkarData> ds;
  for (const auto &term : terms)
  {
    if (const auto *t = std::get_if<config::DjordjevicSarkarData>(&term))
    {
      ds.push_back(*t);
    }
  }
  return ds;
}

std::complex<double> Susceptibility(const std::vector<config::PermittivityTermData> &terms,
                                    std::complex<double> s)
{
  std::complex<double> chi = 0.0;
  for (const auto &term : terms)
  {
    chi += TermSusceptibility(term, s);
  }
  return chi;
}

double StaticSusceptibility(const std::vector<config::PermittivityTermData> &terms)
{
  constexpr double inf = std::numeric_limits<double>::infinity();
  double chi = 0.0;
  for (const auto &term : terms)
  {
    chi += std::visit(
        Overloaded{[](const config::DrudeData &) { return inf; },
                   [](const config::DebyeData &t) { return t.delta; },
                   [](const config::LorentzData &t) { return t.delta; },
                   [](const config::PermittivityPoleData &t)
                   {
                     if (t.pole == std::complex<double>(0.0, 0.0))
                     {
                       return inf;
                     }
                     const auto c = -t.residue / t.pole;
                     return (t.pole.imag() > 0.0) ? 2.0 * c.real() : c.real();
                   },
                   [](const config::DjordjevicSarkarData &t)
                   { return t.strength * std::log(t.omega_upper / t.omega_lower); }},
        term);
  }
  return chi;
}

double StoredEnergySusceptibility(const std::vector<config::PermittivityTermData> &terms,
                                  double omega)
{
  double w = 0.0;
  for (const auto &term : terms)
  {
    w += TermStoredEnergy(term, omega);
  }
  return w;
}

double DissipationSusceptibility(const std::vector<config::PermittivityTermData> &terms,
                                 double omega)
{
  return -Susceptibility(terms, std::complex<double>(0.0, omega)).imag();
}

void CheckTerms(double eps_inf, const std::vector<config::PermittivityTermData> &terms)
{
  MFEM_VERIFY(std::isfinite(eps_inf) && eps_inf > 0.0,
              "Material permittivity HighFrequency must be positive for a positive field "
              "energy!");
  for (const auto &term : terms)
  {
    std::visit(
        Overloaded{
            [](const config::DrudeData &t)
            {
              MFEM_VERIFY(std::isfinite(t.omega_p) && std::isfinite(t.gamma) &&
                              t.omega_p > 0.0 && t.gamma > 0.0,
                          "Drude PlasmaFrequency and CollisionFrequency must be positive!");
            },
            [](const config::DebyeData &t)
            {
              MFEM_VERIFY(std::isfinite(t.delta) && t.delta > 0.0,
                          "Debye DeltaPermittivity must be positive: a negative increment "
                          "has gain at every frequency!");
              MFEM_VERIFY(std::isfinite(t.tau) && t.tau > 0.0,
                          "Debye RelaxationTime must be positive!");
            },
            [](const config::LorentzData &t)
            {
              MFEM_VERIFY(
                  std::isfinite(t.delta) && t.delta > 0.0,
                  "Lorentz DeltaPermittivity must be positive: a negative increment "
                  "has gain, or negative stored energy if undamped, at every "
                  "frequency!");
              MFEM_VERIFY(std::isfinite(t.omega_0) && std::isfinite(t.gamma) &&
                              t.omega_0 > 0.0 && t.gamma >= 0.0,
                          "Lorentz ResonanceFrequency must be positive and "
                          "DampingFrequency must be nonnegative!");
              // Near critical damping the two poles of the solver's pole-residue form merge
              // and the residues grow without bound with opposite signs, so that form loses
              // about log10(1/|γ - 2ω₀|) digits.
              MFEM_VERIFY(
                  std::abs(t.gamma - 2.0 * t.omega_0) > 1.0e-6 * 2.0 * t.omega_0,
                  "A critically damped Lorentz term (DampingFrequency within a "
                  "relative 1e-6 of 2 * ResonanceFrequency) has a repeated pole and "
                  "cannot be represented by simple pole-residue terms!");
            },
            [](const config::PermittivityPoleData &t)
            {
              MFEM_VERIFY(std::isfinite(t.pole.real()) && std::isfinite(t.pole.imag()) &&
                              std::isfinite(t.residue.real()) &&
                              std::isfinite(t.residue.imag()),
                          "Material permittivity poles and residues must be finite!");
              MFEM_VERIFY(t.pole.real() <= 0.0,
                          "Material permittivity poles must satisfy Re(Pole) <= 0!");
              MFEM_VERIFY(t.residue != std::complex<double>(0.0, 0.0),
                          "Material permittivity pole residues must be nonzero!");
              MFEM_VERIFY(t.pole.imag() >= 0.0,
                          "Complex material permittivity poles must be in the upper "
                          "half-plane!");
              MFEM_VERIFY(t.pole.imag() > 0.0 || t.residue.imag() == 0.0,
                          "A real material permittivity pole requires a real residue!");
            },
            [](const config::DjordjevicSarkarData &t)
            {
              MFEM_VERIFY(std::isfinite(t.strength) && t.strength > 0.0,
                          "DjordjevicSarkar Strength must be positive: a negative strength "
                          "has gain at every frequency!");
              MFEM_VERIFY(std::isfinite(t.omega_lower) && std::isfinite(t.omega_upper) &&
                              t.omega_lower > 0.0 && t.omega_upper > t.omega_lower,
                          "DjordjevicSarkar frequencies must satisfy UpperFrequency > "
                          "LowerFrequency > 0!");
            }},
        term);
  }
  const auto poles = PoleTerms(terms);
  for (const auto &term : poles)
  {
    std::complex<double> residue = 0.0;
    for (const auto &other : poles)
    {
      residue += (other.pole == term.pole) ? other.residue : 0.0;
    }
    MFEM_VERIFY(residue != 0.0, "Material permittivity terms with a common pole cancel "
                                "exactly; remove them!");
  }
}

bool Contains(const std::vector<FrequencyInterval> &intervals, double omega)
{
  return std::any_of(intervals.begin(), intervals.end(), [omega](const auto &i)
                     { return omega >= i.lower && omega <= i.upper; });
}

bool IsPassiveByConstruction(const std::vector<config::PermittivityTermData> &terms)
{
  return std::none_of(
      terms.begin(), terms.end(), [](const auto &term)
      { return std::holds_alternative<config::PermittivityPoleData>(term); });
}

ModelDefects FindDefects(double eps_inf,
                         const std::vector<config::PermittivityTermData> &terms)
{
  std::vector<double> rates;
  std::vector<std::pair<double, double>> resonances;
  for (const auto &term : terms)
  {
    TermRates(term, rates, resonances);
  }
  if (rates.empty())
  {
    return {};
  }

  // A logarithmic grid three decades beyond the model's rates, at 64 points per decade,
  // with every resonance sampled at a quarter of its half-width over 16 half-widths. A
  // lossless resonance has no width, so its sign is resolved at 1e-6 of its frequency.
  const auto [rate_min, rate_max] = std::minmax_element(rates.begin(), rates.end());
  const double omega_min = 1.0e-3 * *rate_min, omega_max = 1.0e3 * *rate_max;
  constexpr int per_decade = 64;
  const int n = static_cast<int>(std::ceil(per_decade * std::log10(omega_max / omega_min)));
  std::vector<double> grid;
  for (int k = 0; k <= n; k++)
  {
    grid.push_back(omega_min * std::pow(10.0, static_cast<double>(k) / per_decade));
  }
  for (const auto &[center, half_width] : resonances)
  {
    const double h = std::max(half_width, 1.0e-6 * center);
    for (int k = -64; k <= 64; k++)
    {
      const double omega = center + 0.25 * k * h;
      if (omega > 0.0)
      {
        grid.push_back(omega);
      }
    }
  }
  std::sort(grid.begin(), grid.end());

  // A value is negative beyond rounding when it is below a relative 1e-10 of the sum of
  // the magnitudes of its contributions.
  constexpr double tol = 1.0e-10;
  auto gain = [&](double omega) -> std::optional<bool>
  {
    const std::complex<double> s(0.0, omega);
    double d = 0.0, scale = 0.0;
    for (const auto &term : terms)
    {
      const auto chi = TermSusceptibility(term, s);
      d -= chi.imag();
      scale += std::abs(chi.imag());
    }
    if (!std::isfinite(d))
    {
      return std::nullopt;
    }
    return d < -tol * scale;
  };
  auto negative_energy = [&](double omega) -> std::optional<bool>
  {
    double w = eps_inf, scale = eps_inf;
    for (const auto &term : terms)
    {
      const double wk = TermStoredEnergy(term, omega);
      w += wk;
      scale += std::abs(wk);
    }
    if (!std::isfinite(w))
    {
      return std::nullopt;
    }
    return w < -tol * scale;
  };
  return {FindIntervals(grid, gain), FindIntervals(grid, negative_energy)};
}

}  // namespace palace::permittivity
