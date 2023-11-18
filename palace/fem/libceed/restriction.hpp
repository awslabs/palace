// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_RESTRICTION_HPP
#define PALACE_LIBCEED_RESTRICTION_HPP

#include <vector>
#include <ceed.h>
#include "fem/fespace.hpp"

namespace palace::ceed
{

void InitRestriction(const FiniteElementSpace &fespace,
                     const std::vector<std::size_t> &indices, bool use_bdr, bool is_interp,
                     bool is_interp_range, Ceed ceed, CeedElemRestriction *restr);

}  // namespace palace::ceed

#endif  // PALACE_LIBCEED_RESTRICTION_HPP
