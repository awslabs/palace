// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_COMPLEX_CONDUCTIVITY_MODEL_HPP
#define PALACE_MODELS_COMPLEX_CONDUCTIVITY_MODEL_HPP

#include <complex>
#include <vector>
#include "utils/labels.hpp"

namespace palace
{

namespace config
{

struct ComplexConductivityData;

}  // namespace config

namespace superconductor
{

// Weak-coupling BCS ratio Δ(0)/(kB Tc) = π exp(-γ), with γ the Euler-Mascheroni constant.
inline constexpr double BCS_GAP_RATIO = 1.7638769888620456;

// Weak-coupling BCS gap Δ(T)/Δ(0) at reduced temperature t = T/Tc (zero for t >= 1).
double BCSGapRatio(double t);

// Mattis-Bardeen conductivity σ/σn = (σ1 - iσ2)/σn in the local (dirty) limit, at reduced
// frequency x = ħω/Δ > 0 and temperature θ = kB T/Δ >= 0 (Mattis and Bardeen, Phys. Rev.
// 111, 412 (1958), Eqs. (3.9)-(3.10)).
std::complex<double> MattisBardeenRatio(double x, double theta);

}  // namespace superconductor

//
// Frequency-dependent complex conductivity σ(ω) = σ1(ω) - iσ2(ω) of a superconducting
// domain in nondimensional units, from the Mattis-Bardeen model or tabulated data.
//
class ComplexConductivityModel
{
private:
  ComplexConductivityType type;

  // Mattis-Bardeen: normal-state conductivity, and gap and temperature as angular
  // frequencies Δ/ħ and kB T/ħ.
  double sigma_n = 0.0, omega_gap = 0.0, omega_T = 0.0;

  // Table: monotone cubic (PCHIP) interpolants of σ1 and ωσ2 in log(ω), with values y and
  // slopes d, held constant outside the table.
  std::vector<double> log_omega, y1, d1, y2, d2;

  // Interpolant value, or its derivative in log(ω).
  double EvalTable(double log_w, const std::vector<double> &y, const std::vector<double> &d,
                   bool derivative = false) const;

public:
  ComplexConductivityModel(const config::ComplexConductivityData &data);

  // Complex conductivity at angular frequency ω > 0.
  std::complex<double> operator()(double omega) const;

  // Frequency derivative dσ/dω at angular frequency ω > 0 (central difference for
  // Mattis-Bardeen).
  std::complex<double> Derivative(double omega) const;

  ComplexConductivityType GetType() const { return type; }

  // Gap angular frequency Δ/ħ (Mattis-Bardeen).
  double GetGapFrequency() const { return omega_gap; }
};

}  // namespace palace

#endif  // PALACE_MODELS_COMPLEX_CONDUCTIVITY_MODEL_HPP
