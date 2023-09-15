// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_RESTRICTION_HPP
#define PALACE_LIBCEED_RESTRICTION_HPP

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

using RestrKey = std::tuple<Ceed, void *, int, int, int>;

struct RestrHash
{
  std::size_t operator()(const RestrKey &k) const
  {
    return CeedHashCombine(
        CeedHashCombine(CeedHash(std::get<0>(k)), CeedHash(std::get<1>(k))),
        CeedHashCombine(
            CeedHash(std::get<2>(k)),
            CeedHashCombine(CeedHash(std::get<3>(k)), CeedHash(std::get<4>(k)))));
  }
};

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
