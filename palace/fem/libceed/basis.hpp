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

// Initialize an ordinary libCEED basis from an arbitrary fixed integration rule (which
// need not be tensor-product, for example a rule mapped from a face into an element).
// This always uses MFEM FULL tabulation in native dof ordering, so the corresponding
// element restriction must also use native ordering (is_interp = true for scalar tensor
// product spaces). The rule is fixed when the operator is assembled.
void InitBasisFromRule(const mfem::FiniteElement &fe, const mfem::IntegrationRule &ir,
                       int num_comp, Ceed ceed, CeedBasis *basis);

// Reuse a process-lifetime ordinary basis for finite mapped face/domain rule families.
// Do not use for arbitrary probe point clouds.
void InitCachedBasisFromRule(const mfem::FiniteElement &fe, const mfem::IntegrationRule &ir,
                             int num_comp, Ceed ceed, CeedBasis *basis);

void InitInterpolatorBasis(const mfem::FiniteElement &trial_fe,
                           const mfem::FiniteElement &test_fe, int trial_num_comp,
                           int test_num_comp, Ceed ceed, CeedBasis *basis);

}  // namespace palace::ceed

#endif  // PALACE_LIBCEED_BASIS_HPP
