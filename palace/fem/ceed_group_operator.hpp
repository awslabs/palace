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
  Ceed ceed = nullptr;
  CeedOperator op = nullptr;
  std::vector<std::pair<std::string, int>> field_sources;
  // Optional mesh-node source and operator input names. Unlike fixed coefficient data,
  // mesh coordinates may be scaled in place between applies for dimensional output, so
  // every listed libCEED vector is re-pointed at current MFEM memory on each evaluation.
  const Vector *mesh_nodes = nullptr;
  std::vector<std::string> mesh_node_fields;
  mutable std::vector<CeedVector> mesh_node_vecs;
  // Cached passive field vectors for field_sources, populated on first apply to avoid
  // repeated string lookups during pointwise libCEED operator application.
  mutable std::vector<std::pair<CeedVector, int>> field_vec_sources;
  // Reusable output vector wrapper. The pointed-to MFEM Vector data is supplied at
  // apply time, but the libCEED vector object itself can be retained across repeated
  // postprocessing evaluations instead of being created/destroyed for every group apply.
  mutable CeedVector out_vec = nullptr;
  mutable CeedSize out_size = 0;
};

// Populate cached libCEED vector handles for the named passive fields of a group
// operator. Calling this right after assembly moves field-name lookup overhead out of
// the first timed postprocessing apply; ApplyAddGroupOperators also calls it lazily for
// older call sites or any groups assembled without explicit caching.
void CacheGroupOperatorFieldVectors(const CeedGroupOperator &group);

// Re-point the passive field inputs of each group operator at the given source vectors
// and accumulate into the output vector with CeedOperatorApplyAdd.
void ApplyAddGroupOperators(const std::vector<CeedGroupOperator> &groups,
                            const std::array<const Vector *, 4> &srcs, const Vector &out);

// Destroy the libCEED references owned by a group-operator vector and clear it. The Ceed
// context itself is borrowed and is not destroyed here.
void DestroyGroupOperators(std::vector<CeedGroupOperator> &groups);

}  // namespace fem

}  // namespace palace

#endif  // PALACE_FEM_CEED_GROUP_OPERATOR_HPP
