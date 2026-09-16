// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_CONSTANTS_HPP
#define PALACE_UTILS_CONSTANTS_HPP

#include <cmath>

namespace palace::electromagnetics
{

//
// Define physical constants for nondimensionalization.
//
// SI units switched definition of the Ampere in 2019 so that mu0 is now defined in terms of
// the fine-structure constant alpha. Updates values: Mohr, P. J.; Newell, D. B.; Taylor, B.
// N.; Tiesinga, E. CODATA Recommended Values of the Fundamental Physical Constants: 2022.
// Rev. Mod. Phys. 2025, 97 (2), 025002. https://doi.org/10.1103/RevModPhys.97.025002.

// Speed of light in free space [m/s].
// Exact.
static constexpr double c0_ = 299'792'458;

// Permeability of free space [H/m].
// CODATA value with error: 1.256 637 061 27(20) * 10^(−6)
static constexpr double mu0_ = 1.256'637'061'27e-6;

// Permittivity of free space [F/m].
// CODATA value with error: 8.854 187 8188(14) * 10^(−12)
static constexpr double epsilon0_ = 1.0 / (mu0_ * c0_ * c0_);

// Impedance of free space [Ω].
// CODATA value with error: 376.730 313 412(59)
static constexpr double Z0_ = mu0_ * c0_;  // = sqrt(mu0 / epsilon0)

//
// Time-harmonic (phasor) conventions.
//
// Frequency domain solutions (driven, eigenmode, boundary mode) are peak phasors: the
// physical field is Re{E e^{iωt}}. Quadratic quantities (stored energies, powers, Poynting
// vectors) are reported as time averages over one period, which carry a factor of 1/2
// relative to the same product of the complex phasors:
//                    <Re{a e^{iωt}} Re{b e^{iωt}}> = 1/2 Re{a b⋆}.
// Real-valued fields from time domain or static solutions are instantaneous and carry no
// such factor. These weights are applied by the owners of the physical output quantities,
// not by generic norms, mass matrices, or mode-overlap functionals, which keep their
// mathematical definitions.
//

// Weight relating the product of two fields to its physical (time-averaged for peak
// phasors, instantaneous for real fields) value.
constexpr double TimeAverageWeight(bool complex_phasor)
{
  return complex_phasor ? 0.5 : 1.0;
}

// Amplitude of a source normalized to unit physical incident power, relative to a source
// whose raw amplitude product is unity: for peak phasors |a|² = 2 P_avg, so a unit time-
// averaged power requires a peak amplitude of sqrt(2).
inline double UnitPowerAmplitude(bool complex_phasor)
{
  return complex_phasor ? std::sqrt(2.0) : 1.0;
}

}  // namespace palace::electromagnetics

#endif  // PALACE_UTILS_CONSTANTS_HPP
