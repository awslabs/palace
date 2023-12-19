// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_COEFFICIENT_HPP
#define PALACE_LIBCEED_COEFFICIENT_HPP

#include <vector>

union CeedIntScalar;

namespace palace
{

class MaterialPropertyCoefficient;

namespace ceed
{

template <int DIM>
std::vector<CeedIntScalar> PopulateCoefficientContext(const MaterialPropertyCoefficient *Q,
                                                      double a = 1.0);

template <int DIM, int DIM_MASS>
std::vector<CeedIntScalar>
PopulateCoefficientContext(const MaterialPropertyCoefficient *Q,
                           const MaterialPropertyCoefficient *Q_mass, double a = 1.0,
                           double a_mass = 1.0);

}  // namespace ceed

}  // namespace palace

#endif  // PALACE_LIBCEED_COEFFICIENT_HPP
