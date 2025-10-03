// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_STRATTONCHU_HPP
#define PALACE_MODELS_STRATTONCHU_HPP

#include <vector>
#include <mfem.hpp>
#include <linalg/vector.hpp>
#include "fem/gridfunction.hpp"
#include "models/materialoperator.hpp"

namespace palace
{
void AddStrattonChuIntegrandAtElement(const GridFunction &E, const GridFunction &B,
                                      const MaterialOperator &mat_op, double omega_re,
                                      double omega_im,
                                      std::vector<std::array<double, 3>> &r_naughts,
                                      mfem::ElementTransformation &T,
                                      const mfem::IntegrationRule &ir,
                                      std::vector<std::array<double, 3>> &integrand_r,
                                      std::vector<std::array<double, 3>> &integrand_i);

}  // namespace palace

#endif  // PALACE_MODELS_STRATTONCHU_HPP
