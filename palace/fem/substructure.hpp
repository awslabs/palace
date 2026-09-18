// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_SUBSTRUCTURE_HPP
#define PALACE_FEM_SUBSTRUCTURE_HPP

#include <vector>
#include <mfem.hpp>

namespace palace
{

// For a conforming submesh finite element space cut from a parent space, build the signed
// map from each submesh (local) DOF to its parent DOF: par_dof[i] is the parent DOF and
// sign[i] the ±1 orientation (always +1 for H1; ±1 for H(curl)/H(div)). Read from the
// element DOF tables via the submesh parent-element map (exact for order-1 without DOF
// transformations). Returns false for any DOF not mapping to a single ±1 parent DOF (e.g.
// nonconforming/hanging DOFs), which substructuring does not support.
bool BuildSubMeshDofMap(const mfem::ParFiniteElementSpace &sub_fespace,
                        const mfem::ParFiniteElementSpace &parent_fespace,
                        const mfem::Array<int> &parent_element_ids,
                        std::vector<int> &par_dof, std::vector<double> &sign);

}  // namespace palace

#endif  // PALACE_FEM_SUBSTRUCTURE_HPP
