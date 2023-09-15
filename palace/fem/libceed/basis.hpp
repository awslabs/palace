// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_BASIS_HPP
#define PALACE_LIBCEED_BASIS_HPP

#include <tuple>
#include <unordered_map>
#include <vector>
#include <ceed.h>
#include <mfem.hpp>
#include "utils.hpp"

namespace palace::ceed
{

namespace internal
{

using BasisKey = std::tuple<Ceed, void *, void *, int>;

struct BasisHash
{
  std::size_t operator()(const BasisKey &k) const
  {
    return CeedHashCombine(
        CeedHashCombine(CeedHash(std::get<0>(k)), CeedHash(std::get<1>(k))),
        CeedHashCombine(CeedHash(std::get<2>(k)), CeedHash(std::get<3>(k))));
  }
};

extern std::unordered_map<BasisKey, CeedBasis, BasisHash> basis_map;

}  // namespace internal

void InitBasis(const mfem::FiniteElementSpace &fespace, const mfem::FiniteElement &fe,
               const mfem::IntegrationRule &ir, Ceed ceed, CeedBasis *basis);

inline void InitBasis(const mfem::FiniteElementSpace &fespace,
                      const mfem::IntegrationRule &ir, const std::vector<int> &indices,
                      bool use_bdr, Ceed ceed, CeedBasis *basis)
{
  const mfem::FiniteElement &fe =
      use_bdr ? *fespace.GetBE(indices[0]) : *fespace.GetFE(indices[0]);
  InitBasis(fespace, fe, ir, ceed, basis);
}

void InitInterpolatorBasis(const mfem::FiniteElementSpace &trial_fes,
                           const mfem::FiniteElementSpace &test_fes,
                           const mfem::FiniteElement &trial_fe,
                           const mfem::FiniteElement &test_fe, Ceed ceed, CeedBasis *basis);

inline void InitInterpolatorBasis(const mfem::FiniteElementSpace &trial_fespace,
                                  const mfem::FiniteElementSpace &test_fespace,
                                  const std::vector<int> &indices, Ceed ceed,
                                  CeedBasis *basis)
{
  const mfem::FiniteElement &trial_fe = *trial_fespace.GetFE(indices[0]);
  const mfem::FiniteElement &test_fe = *test_fespace.GetFE(indices[0]);
  InitInterpolatorBasis(trial_fespace, test_fespace, trial_fe, test_fe, ceed, basis);
}

}  // namespace palace::ceed

#endif  // PALACE_LIBCEED_BASIS_HPP
