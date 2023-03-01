// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FE_UTILS_HPP
#define PALACE_FE_UTILS_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>

namespace palace
{

// Construct a Finite Element Collection. If pc_pmg is true, construct a
// p-multigrid collection, in increasing polynomial order. If pc_lor is true,
// use mfem::BasisTypeIntegratedGLL for Vector Finite Element Collections. p is
// the polynomial order and dim is the spatial dimension.
template <typename FECollection>
std::vector<std::unique_ptr<FECollection>> ConstructFECollections(bool pc_pmg, bool pc_lor,
                                                                  int p, int dim);

// Construct a Finite Element Space Hierarchy given a sequence of meshes and
// finite element collections. Dirichlet boundary conditions are additionally
// marked.
template <typename FECollection>
mfem::ParFiniteElementSpaceHierarchy ConstructFiniteElementSpaceHierarchy(
    const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
    const std::vector<std::unique_ptr<FECollection>> &fecs,
    const mfem::Array<int> &dbc_marker);

// Construct a degenerate Finite Element Space Hierarchy given a single mesh and
// finite element collection. Unnecessary to pass the dirichlet boundary
// conditions as they need not be incorporated in any inter-space projectors.
template <typename FECollection>
mfem::ParFiniteElementSpaceHierarchy
ConstructFiniteElementSpaceHierarchy(mfem::ParMesh &mesh, const FECollection &fec);

}  // namespace palace

#endif  // PALACE_FE_UTILS_HPP