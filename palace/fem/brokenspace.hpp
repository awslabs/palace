// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_BROKEN_SPACE_HPP
#define PALACE_FEM_BROKEN_SPACE_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>

namespace palace
{
class FiniteElementSpace;

namespace fem
{

// Private representation of an owned broken space. Element maps retain MFEM's signed
// native ordering; neighbour maps use ExchangeFaceNbrData's element ordering. Only local
// element reads vote for true-DOF ownership. Trace-only rows add storage, not unknowns.
struct BrokenDofLayout
{
  std::vector<int> offsets, dofs, nbr_offsets, nbr_dofs, row_origin;
  mfem::Array<HYPRE_BigInt> row_offsets, col_offsets;
  std::unique_ptr<mfem::HypreParMatrix> P;
  std::unique_ptr<mfem::SparseMatrix> R;
};

// Build on the current mesh, collectively. Each entity-side is a component of its
// incident elements connected through unmarked faces. MFEM's P supplies the coefficient
// chart of a master entity, including its hanging trace restrictions; its coefficients
// are unchanged. No root mesh, refinement history, or transported side map is used.
std::unique_ptr<BrokenDofLayout>
BuildBrokenDofLayout(mfem::ParFiniteElementSpace &space,
                     const std::vector<int> &break_attributes);

// Exact, rank-local left inverse. Direct rows are selected; remaining coefficients are
// recovered from complete injective trace blocks. Remote columns must cancel as well.
// Construction is collective; application never communicates or solves a system.
std::unique_ptr<mfem::SparseMatrix> BuildBrokenRestriction(const FiniteElementSpace &space);

}  // namespace fem
}  // namespace palace

#endif  // PALACE_FEM_BROKEN_SPACE_HPP
