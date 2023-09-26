// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_BASIS_HPP
#define PALACE_LIBCEED_BASIS_HPP

#include <unordered_map>
#include <vector>
#include <ceed.h>
#include <mfem.hpp>

namespace palace::ceed
{

void InitBasis(const mfem::ParFiniteElementSpace &fespace, const mfem::FiniteElement &fe,
               const mfem::IntegrationRule &ir, Ceed ceed, CeedBasis *basis);

inline void InitBasis(const mfem::ParFiniteElementSpace &fespace,
                      const mfem::IntegrationRule &ir, const std::vector<int> &indices,
                      bool use_bdr, Ceed ceed, CeedBasis *basis)
{
  const mfem::FiniteElement &fe =
      use_bdr ? *fespace.GetBE(indices[0]) : *fespace.GetFE(indices[0]);
  InitBasis(fespace, fe, ir, ceed, basis);
}

void InitInterpolatorBasis(const mfem::ParFiniteElementSpace &trial_fes,
                           const mfem::ParFiniteElementSpace &test_fes,
                           const mfem::FiniteElement &trial_fe,
                           const mfem::FiniteElement &test_fe, Ceed ceed, CeedBasis *basis);

inline void InitInterpolatorBasis(const mfem::ParFiniteElementSpace &trial_fespace,
                                  const mfem::ParFiniteElementSpace &test_fespace,
                                  const std::vector<int> &indices, Ceed ceed,
                                  CeedBasis *basis)
{
  const mfem::FiniteElement &trial_fe = *trial_fespace.GetFE(indices[0]);
  const mfem::FiniteElement &test_fe = *test_fespace.GetFE(indices[0]);
  InitInterpolatorBasis(trial_fespace, test_fespace, trial_fe, test_fe, ceed, basis);
}

namespace internal
{

// Destroy the cached CeedBasis objects.
void ClearBasisCache();

}  // namespace internal

}  // namespace palace::ceed

#endif  // PALACE_LIBCEED_BASIS_HPP
