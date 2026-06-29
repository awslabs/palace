// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_CEED_GROUP_OPERATOR_HPP
#define PALACE_FEM_CEED_GROUP_OPERATOR_HPP

#include <array>
#include <string>
#include <utility>
#include <vector>
#include "fem/libceed/ceed.hpp"
#include "linalg/vector.hpp"

namespace palace
{

namespace fem
{

// An assembled libCEED operator over one group of elements, with the named passive field
// inputs re-pointed at caller data (by source vector index) on each evaluation.
struct CeedGroupOperator
{
  Ceed ceed;
  CeedOperator op;
  std::vector<std::pair<std::string, int>> field_sources;
  // Optional retained QFunction context handle for in-place runtime updates (e.g.
  // far-field frequency) without reassembly; nullptr if the operator has no context or
  // the context is not updated. Owned by the group (destroyed with it).
  CeedQFunctionContext ctx = nullptr;
  // Cached passive field vectors for field_sources, populated on first apply to avoid
  // repeated string lookups in libCEED during ParaView point-field output.
  mutable std::vector<std::pair<CeedVector, int>> field_vec_sources;
  // Reusable output vector wrapper. The pointed-to MFEM Vector data is supplied at
  // apply time, but the libCEED vector object itself can be retained across repeated
  // postprocessing evaluations instead of being created/destroyed for every group apply.
  mutable CeedVector out_vec = nullptr;
  mutable CeedSize out_size = 0;
};

// Re-point the passive field inputs of each group operator at the given source vectors
// and accumulate into the output vector with CeedOperatorApplyAdd. A field source index
// of 4 (out of the srcs range) selects the optional imported vector instead, used to
// feed face-neighbor (ghost) field values exchanged for two-sided interior boundaries on
// parallel interfaces.
void ApplyAddGroupOperators(const std::vector<CeedGroupOperator> &groups,
                            const std::array<const Vector *, 4> &srcs, const Vector &out,
                            const Vector *imported = nullptr);

}  // namespace fem

}  // namespace palace

#endif  // PALACE_FEM_CEED_GROUP_OPERATOR_HPP
