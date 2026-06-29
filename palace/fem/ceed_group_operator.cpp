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

void ApplyAddGroupOperators(const std::vector<CeedGroupOperator> &groups,
                            const std::array<const Vector *, 4> &srcs, const Vector &out,
                            const Vector *imported)
{
  for (const auto &group : groups)
  {
    for (const auto &[name, source] : group.field_sources)
    {
      // Source index 4 selects an optional imported vector, used by surface reductions
      // and boundary point fields for face-neighbor field values. The operator's
      // restriction slices and transposes the shared vector to the per-element layout.
      const Vector *sv = (source < 4) ? srcs[source] : imported;
      MFEM_ASSERT(sv, "Missing source vector for libCEED field input!");
      CeedOperatorField field;
      CeedVector field_vec;
      PalaceCeedCall(group.ceed, CeedOperatorGetFieldByName(group.op, name.c_str(), &field));
      PalaceCeedCall(group.ceed, CeedOperatorFieldGetVector(field, &field_vec));
      ceed::InitCeedVector(*sv, group.ceed, &field_vec, false);
      PalaceCeedCall(group.ceed, CeedVectorDestroy(&field_vec));
    }
    CeedMemType out_mem;
    PalaceCeedCall(group.ceed, CeedGetPreferredMemType(group.ceed, &out_mem));
    if (!mfem::Device::Allows(mfem::Backend::DEVICE_MASK) && out_mem == CEED_MEM_DEVICE)
    {
      out_mem = CEED_MEM_HOST;
    }
    auto *out_data = const_cast<Vector &>(out).ReadWrite(out_mem == CEED_MEM_DEVICE);
    const CeedSize out_size = out.Size();
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
