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

// Permittivity of free space [F/m].
static constexpr double epsilon0_ = 8.8541878176e-12;

// Permeability of free space [H/m].
static constexpr double mu0_ = 4.0e-7 * M_PI;

// Speed of light in free space [m/s].
static
#if defined(PALACE_WITH_CONSTEXPR_SQRT)
    constexpr
#else
    const
#endif
    double c0_ = 1.0 / std::sqrt(epsilon0_ * mu0_);

// Impedance of free space [Ω].
static
#if defined(PALACE_WITH_CONSTEXPR_SQRT)
    constexpr
#else
    const
#endif
    double Z0_ = std::sqrt(mu0_ / epsilon0_);

}  // namespace palace::electromagnetics

#endif  // PALACE_UTILS_CONSTANTS_HPP
