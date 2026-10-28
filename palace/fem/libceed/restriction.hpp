// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_RESTRICTION_HPP
#define PALACE_LIBCEED_RESTRICTION_HPP

#include <vector>
#include "fem/libceed/ceed.hpp"

namespace mfem
{

class FiniteElementSpace;

}  // namespace mfem

namespace palace::ceed
{

// Optional remapping of element DOFs to L-vector indices, used by broken finite element
// spaces: for local element e, entries [offsets[e], offsets[e + 1]) give native local DOF
// indices (local) and the new, unsigned L-vector indices (ldof) replacing the ones of the
// underlying finite element space, whose orientations are kept. The L-vector size of the
// restriction is l_size.
struct ElementDofRemap
{
  const int *offsets = nullptr, *local = nullptr, *ldof = nullptr;
  CeedSize l_size = 0;
};

void InitRestriction(const mfem::FiniteElementSpace &fespace,
                     const std::vector<int> &indices, bool use_bdr, bool is_interp,
                     bool is_interp_range, Ceed ceed, CeedElemRestriction *restr,
                     const ElementDofRemap *remap = nullptr);

}  // namespace palace::ceed

#endif  // PALACE_LIBCEED_RESTRICTION_HPP
