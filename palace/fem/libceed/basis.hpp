// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_BASIS_HPP
#define PALACE_LIBCEED_BASIS_HPP

#include <unordered_map>
#include <vector>
#include <ceed.h>
#include <mfem.hpp>
#include "utils.hpp"

namespace palace::ceed
{

namespace internal
{

struct BasisKey
{
  Ceed ceed;
  const mfem::FiniteElement *trial_fe;
  const mfem::FiniteElement *test_fe;
  const mfem::IntegrationRule *ir;
  int ncomp;
  bool operator==(const BasisKey &k) const
  {
    return (ceed == k.ceed && trial_fe == k.trial_fe && test_fe == k.test_fe && ir == k.ir && ncomp == k.ncomp);
  }
};

struct BasisHash
{
  std::size_t operator()(const BasisKey &k) const
  {
    std::size_t hash = 0;
    CeedHashCombine(hash, k.ceed, k.trial_fe, k.test_fe, k.ir, k.ncomp);
    return hash;
  }
};

// struct BasisEqual : public std::binary_function<BasisKey, BasisKey, bool>
// {
//   bool operator()(const BasisKey &k0, const BasisKey &k1) const
//   {
//     return (k0.ceed == k1.ceed && k0.trial_fe == k1.trial_fe && k0.test_fe == k1.test_fe && k0.ir == k1.ir && k0.ncomp == k1.ncomp);
//   }
// };

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
