// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_CONDUCTIVITY_MODEL_HPP
#define PALACE_UTILS_CONDUCTIVITY_MODEL_HPP

#include <complex>

namespace palace
{

namespace config
{

struct MattisBardeenData;
struct TabulatedConductivityData;

}  // namespace config

namespace conductivity
{

//
// Complex conductivity σ(ω) = σ₁(ω) - iσ₂(ω) of the conductor terms of a
// frequency-dependent permittivity (config::MattisBardeenData,
// config::TabulatedConductivityData), divided by ε₀: an angular rate, with χ(s) =
// σ(ω)/(ε₀s) for s = iω. Defined at real ω > 0 only.
//

// Weak-coupling BCS ratio Δ(0)/(kB Tc) = π exp(-γ), with γ the Euler-Mascheroni constant.
inline constexpr double BCS_GAP_RATIO = 1.7638769888620456;

// Weak-coupling BCS gap Δ(T)/Δ(0) at reduced temperature t = T/Tc (zero for t >= 1).
double BCSGapRatio(double t);

// Mattis-Bardeen conductivity σ/σₙ = (σ₁ - iσ₂)/σₙ in the local (dirty) limit, at reduced
// frequency x = ħω/Δ > 0 and temperature θ = kB T/Δ >= 0 (Mattis and Bardeen, Phys. Rev.
// 111, 412 (1958), Eqs. (3.9)-(3.10)).
std::complex<double> MattisBardeenRatio(double x, double theta);

// σ(ω)/ε₀ at ω > 0, and its derivative with respect to ω (by central difference for
// Mattis-Bardeen). A tabulated conductivity interpolates σ₁ and ωσ₂ with monotone cubic
// (PCHIP) interpolants in log(ω), constant outside the table.
std::complex<double> Conductivity(const config::MattisBardeenData &data, double omega);
std::complex<double> ConductivityDerivative(const config::MattisBardeenData &data,
                                            double omega);
std::complex<double> Conductivity(const config::TabulatedConductivityData &data,
                                  double omega);
std::complex<double> ConductivityDerivative(const config::TabulatedConductivityData &data,
                                            double omega);

}  // namespace conductivity

}  // namespace palace

#endif  // PALACE_UTILS_CONDUCTIVITY_MODEL_HPP
