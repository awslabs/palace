// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_COEFFICIENT_HPP
#define PALACE_LIBCEED_COEFFICIENT_HPP

#include <vector>
#include <mfem.hpp>

union CeedIntScalar;

namespace palace
{

class MaterialPropertyCoefficient;

namespace ceed
{

std::vector<CeedIntScalar> PopulateCoefficientContext(int dim,
                                                      const MaterialPropertyCoefficient *Q,
                                                      mfem::real_t a = 1.0);

std::vector<CeedIntScalar>
PopulateCoefficientContext(int dim_mass, const MaterialPropertyCoefficient *Q_mass, int dim,
                           const MaterialPropertyCoefficient *Q, mfem::real_t a_mass = 1.0,
                           mfem::real_t a = 1.0);

}  // namespace ceed

}  // namespace palace

#endif  // PALACE_LIBCEED_COEFFICIENT_HPP
