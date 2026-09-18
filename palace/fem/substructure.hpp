// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_SUBSTRUCTURE_HPP
#define PALACE_FEM_SUBSTRUCTURE_HPP

#include <memory>
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

// A substructure is a region of the parent mesh (selected by domain attributes) treated as
// its own finite element problem: it owns a ParSubMesh and matching FE space, and the
// signed map from its local DOFs to the parent DOFs (see BuildSubMeshDofMap). Substructures
// sharing parent DOFs are coupled across the interface those shared DOFs form.
class Substructure
{
public:
  Substructure(mfem::ParFiniteElementSpace &parent_fespace,
               const mfem::Array<int> &domain_attrs, mfem::FiniteElementCollection &fec);

  mfem::ParSubMesh &GetSubMesh() { return *submesh; }
  mfem::ParFiniteElementSpace &GetFESpace() { return *fespace; }

  // Signed submesh-DOF -> parent-DOF map.
  const std::vector<int> &GetParentDof() const { return par_dof; }
  const std::vector<double> &GetSign() const { return sign; }

  // False if the submesh->parent DOF map is not a clean signed permutation (nonconforming).
  bool ConformingMap() const { return map_ok; }

private:
  std::unique_ptr<mfem::ParSubMesh> submesh;
  std::unique_ptr<mfem::ParFiniteElementSpace> fespace;
  std::vector<int> par_dof;
  std::vector<double> sign;
  bool map_ok;
};

// Fill owner (size = number of parent DOFs) with a bitmask: bit (1 << k) is set when
// subs[k] covers that parent DOF. A DOF covered by more than one substructure lies on a
// shared interface between them.
void MarkParentDofOwnership(const std::vector<const Substructure *> &subs,
                            int n_parent_dofs, std::vector<int> &owner);

}  // namespace palace

#endif  // PALACE_FEM_SUBSTRUCTURE_HPP
