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

std::vector<CeedIntScalar> PopulateCoefficientContext(int dim,
                                                      const MaterialPropertyCoefficient *Q,
                                                      bool transpose = false,
                                                      double a = 1.0);

std::vector<CeedIntScalar>
PopulateCoefficientContext(int dim_mass, const MaterialPropertyCoefficient *Q_mass, int dim,
                           const MaterialPropertyCoefficient *Q, bool tranpose_mass = false,
                           bool transpose = false, double a_mass = 1.0, double a = 1.0);

}  // namespace ceed

}  // namespace palace

#endif  // PALACE_LIBCEED_COEFFICIENT_HPP
