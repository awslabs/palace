// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_BASIS_HPP
#define PALACE_LIBCEED_BASIS_HPP

#include "fem/libceed/ceed.hpp"

namespace mfem
{

class FiniteElement;
class IntegrationRule;

}  // namespace mfem

namespace palace::ceed
{

void InitBasis(const mfem::FiniteElement &fe, const mfem::IntegrationRule &ir, int num_comp,
               Ceed ceed, CeedBasis *basis);

void InitInterpolatorBasis(const mfem::FiniteElement &trial_fe,
                           const mfem::FiniteElement &test_fe, int trial_num_comp,
                           int test_num_comp, Ceed ceed, CeedBasis *basis);

}  // namespace palace::ceed

#endif  // PALACE_LIBCEED_BASIS_HPP
