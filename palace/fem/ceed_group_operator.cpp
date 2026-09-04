// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/ceed_group_operator.hpp"

#include <mfem.hpp>

namespace palace
{

namespace fem
{

void DestroyGroupOperators(std::vector<CeedGroupOperator> &groups)
{
  for (auto &group : groups)
  {
    for (auto &[field_vec, source] : group.field_vec_sources)
    {
      (void)source;
      PalaceCeedCall(group.ceed, CeedVectorDestroy(&field_vec));
    }
    group.field_vec_sources.clear();
    group.field_vectors_detached = false;
    for (auto &mesh_nodes_vec : group.mesh_node_vecs)
    {
      if (mesh_nodes_vec)
      {
        PalaceCeedCall(group.ceed, CeedVectorDestroy(&mesh_nodes_vec));
      }
    }
    group.mesh_node_vecs.clear();
    group.mesh_node_fields.clear();
    group.mesh_nodes = nullptr;
    if (group.out_vec)
    {
      PalaceCeedCall(group.ceed, CeedVectorDestroy(&group.out_vec));
      group.out_size = 0;
    }
    if (group.ctx)
    {
      PalaceCeedCall(group.ceed, CeedQFunctionContextDestroy(&group.ctx));
    }
    if (group.op)
    {
      PalaceCeedCall(group.ceed, CeedOperatorDestroy(&group.op));
    }
    group.field_sources.clear();
  }
  groups.clear();
}

void CacheGroupOperatorFieldVectors(const CeedGroupOperator &group)
{
  if (group.mesh_node_vecs.size() != group.mesh_node_fields.size())
  {
    group.mesh_node_vecs.clear();
    group.mesh_node_vecs.reserve(group.mesh_node_fields.size());
    for (const auto &name : group.mesh_node_fields)
    {
      CeedOperatorField field;
      CeedVector mesh_nodes_vec = nullptr;
      PalaceCeedCall(group.ceed,
                     CeedOperatorGetFieldByName(group.op, name.c_str(), &field));
      // libCEED may omit a passive input which the selected QFunction does not use.
      if (field)
      {
        PalaceCeedCall(group.ceed, CeedOperatorFieldGetVector(field, &mesh_nodes_vec));
      }
      group.mesh_node_vecs.push_back(mesh_nodes_vec);
    }
  }
  if (group.field_vec_sources.size() == group.field_sources.size())
  {
    return;
  }
  group.field_vec_sources.clear();
  group.field_vec_sources.reserve(group.field_sources.size());
  for (const auto &[name, source] : group.field_sources)
  {
    CeedOperatorField field;
    CeedVector field_vec;
    PalaceCeedCall(group.ceed, CeedOperatorGetFieldByName(group.op, name.c_str(), &field));
    PalaceCeedCall(group.ceed, CeedOperatorFieldGetVector(field, &field_vec));
    group.field_vec_sources.emplace_back(field_vec, source);
  }
}

void DetachGroupOperatorFieldVectors(const std::vector<CeedGroupOperator> &groups)
{
  for (const auto &group : groups)
  {
    CacheGroupOperatorFieldVectors(group);
    if (group.field_vec_sources.empty() || group.field_vectors_detached)
    {
      continue;
    }
    CeedMemType mem;
    PalaceCeedCall(group.ceed, CeedGetPreferredMemType(group.ceed, &mem));
    if (!mfem::Device::Allows(mfem::Backend::DEVICE_MASK) && mem == CEED_MEM_DEVICE)
    {
      mem = CEED_MEM_HOST;
    }
    for (auto &[field_vec, source] : group.field_vec_sources)
    {
      (void)source;
      PalaceCeedCall(group.ceed, CeedVectorTakeArray(field_vec, mem, nullptr));
    }
    group.field_vectors_detached = true;
  }
}

void ApplyAddGroupOperators(const std::vector<CeedGroupOperator> &groups,
                            const std::array<const Vector *, 4> &srcs, const Vector &out,
                            const Vector *imported)
{
  if (groups.empty())
  {
    return;
  }

  CeedMemType out_mem;
  PalaceCeedCall(groups.front().ceed,
                 CeedGetPreferredMemType(groups.front().ceed, &out_mem));
  if (!mfem::Device::Allows(mfem::Backend::DEVICE_MASK) && out_mem == CEED_MEM_DEVICE)
  {
    out_mem = CEED_MEM_HOST;
  }
  auto *out_data = const_cast<Vector &>(out).ReadWrite(out_mem == CEED_MEM_DEVICE);
  const CeedSize out_size = out.Size();

  for (const auto &group : groups)
  {
    CacheGroupOperatorFieldVectors(group);
    if (!group.mesh_node_vecs.empty())
    {
      MFEM_ASSERT(group.mesh_nodes, "Missing mesh-node source for libCEED geometry input!");
      for (auto &mesh_nodes_vec : group.mesh_node_vecs)
      {
        if (mesh_nodes_vec)
        {
          ceed::InitCeedVector(*group.mesh_nodes, group.ceed, &mesh_nodes_vec, false);
        }
      }
    }
    for (auto &[field_vec, source] : group.field_vec_sources)
    {
      // Source index 4 selects an optional imported vector used by surface
      // postprocessing operators for face-neighbor field values. The restriction slices
      // and transposes the shared vector to the per-element layout.
      const Vector *sv = (source < 4) ? srcs[source] : imported;
      MFEM_ASSERT(sv, "Missing source vector for libCEED field input!");
      ceed::InitCeedVector(*sv, group.ceed, &field_vec, false,
                           !group.field_vectors_detached);
    }
    group.field_vectors_detached = false;
    if (!group.out_vec || group.out_size != out_size)
    {
      if (group.out_vec)
      {
        PalaceCeedCall(group.ceed, CeedVectorDestroy(&group.out_vec));
      }
      PalaceCeedCall(group.ceed, CeedVectorCreate(group.ceed, out_size, &group.out_vec));
      group.out_size = out_size;
    }
    PalaceCeedCall(group.ceed,
                   CeedVectorSetArray(group.out_vec, out_mem, CEED_USE_POINTER, out_data));
    PalaceCeedCall(group.ceed, CeedOperatorApplyAdd(group.op, CEED_VECTOR_NONE,
                                                    group.out_vec, CEED_REQUEST_IMMEDIATE));
    PalaceCeedCall(group.ceed, CeedVectorTakeArray(group.out_vec, out_mem, nullptr));
  }
}

}  // namespace fem

}  // namespace palace
