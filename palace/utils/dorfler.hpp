// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DORFLER_HPP
#define PALACE_DORFLER_HPP

#include <mfem.hpp>

namespace palace::utils
{

// Given a vector of estimates, e, and a fraction, compute a partition value E,
// such that that the set of all estimates with value greater than e, K_E is the
// smallest number to achieve sum_{K_E} e >= fraction * sum e. Namely the
// smallest set of elements that will mark the top fraction of the sum of the error.
double ComputeDorflerThreshold(double fraction, const mfem::Vector &e);

}  // namespace palace::utils

#endif  // PALACE_DORFLER_HPP
