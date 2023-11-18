// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_BASIS_HPP
#define PALACE_LIBCEED_BASIS_HPP

#include <ceed.h>
#include <mfem.hpp>

namespace palace::ceed
{

void InitBasis(const mfem::FiniteElement &fe, const mfem::IntegrationRule &ir,
               CeedInt ncomp, Ceed ceed, CeedBasis *basis);

void InitInterpolatorBasis(const mfem::FiniteElement &trial_fe,
                           const mfem::FiniteElement &test_fe, CeedInt trial_ncomp,
                           CeedInt test_ncomp, Ceed ceed, CeedBasis *basis);

}  // namespace palace::ceed

#endif  // PALACE_LIBCEED_BASIS_HPP
