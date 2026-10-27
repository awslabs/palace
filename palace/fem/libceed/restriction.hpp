// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_RESTRICTION_HPP
#define PALACE_LIBCEED_RESTRICTION_HPP

#include <functional>
#include <vector>
#include <mfem.hpp>
#include "fem/libceed/ceed.hpp"

namespace mfem
{

class FiniteElementSpace;

}  // namespace mfem

namespace palace::ceed
{

// Setup-only mapping from MFEM's signed coordinate DOFs to a space's actual L-vector.
// The space owns incident-side resolution, including traces lacking an MFEM boundary FE.
struct ElementDofRemap
{
  std::function<void(int, bool, mfem::Array<int> &)> apply;
  CeedSize l_size = 0;
};

void InitRestriction(const mfem::FiniteElementSpace &fespace,
                     const std::vector<int> &indices, bool use_bdr, bool is_interp,
                     bool is_interp_range, Ceed ceed, CeedElemRestriction *restr,
                     const ElementDofRemap *remap = nullptr);

}  // namespace palace::ceed

#endif  // PALACE_LIBCEED_RESTRICTION_HPP
