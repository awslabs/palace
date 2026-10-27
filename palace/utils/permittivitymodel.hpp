// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_PERMITTIVITY_MODEL_HPP
#define PALACE_UTILS_PERMITTIVITY_MODEL_HPP

#include <complex>
#include <vector>
#include "utils/configfile.hpp"

namespace palace::permittivity
{

//
// Frequency-dependent relative permittivity ε(s) = ε∞ + Σₖ χₖ(s), s = iω, of a material
// with additive terms (config::PermittivityTermData). Rates may be in any consistent unit.
//
// Each term is the polarization P = ε₀χE of an auxiliary differential equation driven by E,
// for which the power delivered to the polarization splits exactly as E·∂P/∂t = dW/dt + D,
// with W the stored energy and D ≥ 0 the dissipated power:
//
//   Relaxation  ẋ + a₀x = E, P = ε₀b₀x:  W = ½ε₀a₀b₀x², D = ε₀b₀ẋ².
//   Oscillator  ẍ + a₁ẋ + a₀x = E, P = ε₀(b₀x + b₁ẋ):
//               W = ε₀(½(b₀ + a₁b₁)ẋ² + a₀b₁xẋ + ½a₀b₀x²), D = ε₀(b₁ẍ² + (a₁b₀ − a₀b₁)ẋ²).
//
// Debye is a relaxation and Djordjevic-Sarkar a positive continuum of them; Lorentz and
// Drude (a₀ = 0) are oscillators with b₁ = 0. For these, W is the only quadratic form of
// the state for which the identity holds for every E, so the stored energy is that of the
// physical model, and positive parameters make W and D nonnegative. A PoleResidue term is a
// relaxation for a real pole and an oscillator for a conjugate pair. With b₁ ≠ 0 the stored
// energy of a passive term is not unique: W + ½δẋ², with |δ| ≤ 2√(b₁(a₁b₀ − a₀b₁)), keeps D
// nonnegative for the same χ. The W above is the midpoint of this interval.
//

// Pole-residue form of the terms other than Djordjevic-Sarkar, for the solver: real poles
// and upper-half-plane poles, whose conjugates are implied, including exact zero poles.
std::vector<config::PermittivityPoleData>
PoleTerms(const std::vector<config::PermittivityTermData> &terms);

// The Djordjevic-Sarkar terms.
std::vector<config::DjordjevicSarkarData>
DjordjevicSarkarTerms(const std::vector<config::PermittivityTermData> &terms);

// Susceptibility χ(s) = Σₖ χₖ(s) at complex s.
std::complex<double> Susceptibility(const std::vector<config::PermittivityTermData> &terms,
                                    std::complex<double> s);

// χ(0), infinite for free carriers (Drude, or a PoleResidue pole at zero).
double StaticSusceptibility(const std::vector<config::PermittivityTermData> &terms);

// The stored electric energy of the polarization at real ω > 0, as a susceptibility: the
// time-averaged stored energy of a field E e^{iωt} in a material with these terms is that
// of a constant permittivity ε∞ + StoredEnergySusceptibility(terms, ω). For a relaxation
// this is Re χ; an oscillator adds the kinetic energy of its motion.
double StoredEnergySusceptibility(const std::vector<config::PermittivityTermData> &terms,
                                  double omega);

// The dissipated power at real ω > 0, as a susceptibility: −Im χ(iω), so that the
// time-averaged dissipated power of the field is that of a loss ω ε₀ (−Im χ) |E|²/2.
double DissipationSusceptibility(const std::vector<config::PermittivityTermData> &terms,
                                 double omega);

// Verify the parameters of the terms. A typed term must be passive on its own: a negative
// strength has gain at every frequency. ε∞ must be positive for a positive field energy.
void CheckTerms(double eps_inf, const std::vector<config::PermittivityTermData> &terms);

// Whether every term is a typed, passive model. Otherwise, use FindDefects.
bool IsPassiveByConstruction(const std::vector<config::PermittivityTermData> &terms);

// Frequency ranges ω ∈ [lower, upper] (upper possibly infinite, lower possibly zero) where
// the model has gain, a negative dissipated power, or a negative stored electric energy
// ε∞ + StoredEnergySusceptibility. Sampled on a logarithmic grid spanning the model's
// rates, refined at every resonance on the scale of its linewidth, with bisected endpoints.
struct FrequencyInterval
{
  double lower, upper;
};
struct ModelDefects
{
  std::vector<FrequencyInterval> gain, negative_energy;
};
ModelDefects FindDefects(double eps_inf,
                         const std::vector<config::PermittivityTermData> &terms);

// Whether ω lies in one of the intervals.
bool Contains(const std::vector<FrequencyInterval> &intervals, double omega);

}  // namespace palace::permittivity

#endif  // PALACE_UTILS_PERMITTIVITY_MODEL_HPP
