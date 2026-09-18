// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_SUBSTRUCTURE_HPP
#define PALACE_FEM_SUBSTRUCTURE_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

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
  const mfem::ParFiniteElementSpace &GetFESpace() const { return *fespace; }

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

// Compact interface indexing from the ownership bitmask: gamma_index[p] in [0, nG) for
// parent DOFs on an interface (more than one owner bit set), else -1. Returns nG.
int BuildInterfaceIndex(const std::vector<int> &owner, std::vector<int> &gamma_index);

// Environment Dirichlet-to-Neumann operator: the dense Schur complement S_E and load g_E of
// an environment substructure, condensed onto the interface in parent orientation and the
// compact interface indexing. Reusable across region redesigns that preserve the interface.
//
// A_env / f_env are the environment substructure's assembled operator and load (in its own
// submesh DOFs); the problem-specific bilinear/linear forms are the caller's
// responsibility.
class DtNBoundaryOperator
{
public:
  DtNBoundaryOperator(const Substructure &environment, const std::vector<int> &gamma_index,
                      const mfem::SparseMatrix &A_env, const mfem::Vector &f_env);

  int Size() const { return S_E.Height(); }
  const mfem::DenseMatrix &Schur() const { return S_E; }  // nG x nG, parent orientation
  const mfem::Vector &Load() const { return g_E; }        // nG

private:
  mfem::DenseMatrix S_E;
  mfem::Vector g_E;
};

// Applies the environment DtN Schur complement S_E as a Palace Operator on a region finite
// element space's DOFs: y = (S_E on the shared interface) x, zero on region-interior DOFs.
// The region's local<->parent orientation signs are applied so the action is in the
// region's local orientation (matching a region operator assembled on its submesh). Serial
// for now; the parallel interface gather is a later step.
class RegionDtNOperator : public Operator
{
public:
  RegionDtNOperator(const Substructure &region, const std::vector<int> &gamma_index,
                    const DtNBoundaryOperator &dtn);
  void Mult(const Vector &x, Vector &y) const override;
  void MultTranspose(const Vector &x, Vector &y) const override { Mult(x, y); }

private:
  std::vector<int> loc_gamma;    // region local DOF -> gamma index, or -1
  std::vector<double> loc_sign;  // region local DOF -> +/-1 (interface DOFs)
  const mfem::DenseMatrix &S;    // dtn.Schur(), parent orientation
  mutable Vector xg, rg;         // interface work vectors
};

}  // namespace palace

#endif  // PALACE_FEM_SUBSTRUCTURE_HPP
