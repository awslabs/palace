// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_RESTRICTION_HPP
#define PALACE_LIBCEED_RESTRICTION_HPP

#include <unordered_map>
#include <vector>
#include <ceed.h>
#include <mfem.hpp>
#include "utils.hpp"

namespace palace::ceed
{

namespace internal
{


///XX TODO WIP

struct RestrKey
{
  Ceed ceed;
  const mfem::FiniteElementSpace *fespace;
  const mfem::FiniteElement *fe;
  int ncomp;
  int unique_range_restr;
  bool operator==(const RestrKey &k) const
  {
    return (ceed == k.ceed && fespace == k.fespace && fe == k.fe && ncomp == k.ncomp && unique_range_restr == k.unique_range_restr);
  }
};

struct RestrHash
{
  std::size_t operator()(const RestrKey &k) const
  {
    std::size_t hash = 0;
    CeedHashCombine(hash, k.ceed, k.fespace, k.fe, k.ncomp, k.unique_range_restr);
    return hash;
  }
};

// struct RestrEqual : public std::binary_function<RestrKey, RestrKey, bool>
// {
//   bool operator()(const RestrKey &k0, const RestrKey &k1) const
//   {
//     return (k0.ceed == k1.ceed && k0.fespace == k1.fespace && k0.fe == k1.fe && k0.ncomp == k1.ncomp && k0.unique_range_restr == k1.unique_range_restr);
//   }
// };

extern std::unordered_map<RestrKey, CeedElemRestriction, RestrHash> restr_map;

}  // namespace internal

void InitRestriction(const mfem::FiniteElementSpace &fespace,
                     const std::vector<int> &indices, bool use_bdr, bool is_interp,
                     bool is_range, Ceed ceed, CeedElemRestriction *restr);

inline void InitRestriction(const mfem::FiniteElementSpace &fespace,
                            const std::vector<int> &indices, bool use_bdr, Ceed ceed,
                            CeedElemRestriction *restr)
{
  InitRestriction(fespace, indices, use_bdr, false, false, ceed, restr);
}

}  // namespace palace::ceed

#endif  // PALACE_LIBCEED_RESTRICTION_HPP
