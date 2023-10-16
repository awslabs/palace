// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_DORFLER_HPP
#define PALACE_UTILS_DORFLER_HPP

#include <mfem.hpp>

namespace palace::utils
{

// Given a vector of estimates, e, and a fraction, compute a partition value E, such that
// that the set of all estimates with value greater than E, K_E, is the smallest set to
// achieve sum_{K_E} e² >= fraction * sum e². Namely the smallest set of elements that
// will mark the top fraction of the sum of the squared error.
// Reference: Willy Dörfler. A convergent adaptive algorithm for Poisson’s equation. SIAM J.
//            Numer. Anal. (1996).
double ComputeDorflerThreshold(double fraction, const mfem::Vector &e);

}  // namespace palace::utils

#endif  // PALACE_UTILS_DORFLER_HPP
