// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_RESTRICTION_HPP
#define PALACE_LIBCEED_RESTRICTION_HPP

#include <unordered_map>
#include <vector>
#include <ceed.h>
#include <mfem.hpp>
#include "hash.hpp"
#include "utils.hpp"

namespace palace::ceed
{

namespace internal
{

extern std::unordered_map<RestrKey, CeedElemRestriction, RestrHash> restr_map;

void ClearRestrictionCache();

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
