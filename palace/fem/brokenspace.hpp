// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_BROKEN_SPACE_HPP
#define PALACE_FEM_BROKEN_SPACE_HPP

#include <cstdint>
#include <vector>
#include <mfem.hpp>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

//
// Support for finite element spaces which are discontinuous ("broken") across selected
// interior boundaries without modifying the mesh.
//
// An interior boundary ("crack") is the set of interior faces carrying one of a given list
// of boundary attributes. Every mesh entity (vertex, edge, face) lying on such a face is
// split into one copy per group of neighboring elements which are connected to each other
// without crossing the interior boundary, unless there is only one such group (for example
// at the free edge of an interior boundary). The groups on the same side of a connected
// interior boundary form a region. The regions, adjacent when they are on different sides
// of a split entity, are two-colored: those of the first color keep the original degrees
// of freedom of their split entities, while those of the second color read a single copy
// (junctions of three or more regions therefore share one copy). Each element then reads
// either the original degrees of freedom of all of its split entities or the copies of all
// of them, as required by the constraints of hanging entities on a nonconforming mesh,
// which can involve the split entities of several interior boundaries. This is not
// possible where the regions cannot be two-colored: at junctions of three pairwise adjacent
// regions, and along an interior boundary whose two sides are joined into a single region
// through the split entities of another one (an air bridge standing on a ground plane, for
// example).
//
// The sides are stored as bitmasks per local element, with bit b for the local entity b of
// the element, where local entities are numbered as vertices [0, nv), then edges
// [nv, nv + ne), then (3D only) faces [nv + ne, nv + ne + nf), in the local ordering of
// mfem::Mesh::GetElementVertices, GetElementEdges, and GetElementFaces.
//
struct CrackSides
{
  // Bit b of copy[e] is set when element e reads the copy of the DOFs of its local entity
  // b. On a nonconforming mesh, this also includes entities which are not split, but whose
  // DOFs are constrained to those of split entities (hanging entities next to an interior
  // boundary): these read the copies of the split DOFs in their constraints.
  std::vector<std::uint32_t> copy;

  // Bit b of split[e] is set when the local entity b of element e is split.
  std::vector<std::uint32_t> split;

  // Return whether any local element reads copies.
  bool Any() const;
};

namespace mesh
{

// Compute the interior boundary sides for a conforming mesh (or a nonconforming mesh
// without hanging entities). The marker is indexed by boundary attribute - 1. Collective.
CrackSides ComputeCrackSides(const mfem::ParMesh &mesh,
                             const mfem::Array<int> &bdr_attr_marker);

// Return whether the mesh has hanging (nonconforming) vertices, edges, or faces, in which
// case the interior boundary sides cannot be discovered from the mesh topology and need to
// be inherited through refinement instead. Collective.
bool HasHangingEntities(const mfem::ParMesh &mesh);

// Interior boundary sides after a refinement of the mesh, given the sides of the coarse
// elements before refinement: an entity of a fine element lying on a split entity of its
// parent element inherits its split and copy flags, and one lying on the boundary of the
// parent element otherwise reads copies if the parent does for any entity in the closure
// of the smallest parent entity containing it (refinement moves no element across an
// interior boundary). The coarse-to-fine transformations are those returned by
// mfem::Mesh::GetRefinementTransforms. Not collective.
CrackSides InheritCrackSides(const mfem::Mesh &fine_mesh,
                             const mfem::CoarseFineTransformations &cf,
                             const CrackSides &coarse_sides);

// Compute the interior boundary sides for a nonconforming mesh with hanging entities, for
// example one refined in a previous simulation, from those of the coarsest mesh without
// hanging entities in its refinement hierarchy: a copy of the mesh is gathered on the root
// process and derefined until it has no hanging entities, and the sides discovered there
// are inherited back through the refinements. Returns false, leaving the sides unchanged,
// when this is not possible: for 3D meshes with anisotropic refinements or with pyramids,
// or when the coarsest mesh of the hierarchy has hanging entities. Collective.
bool ReconstructCrackSides(const mfem::ParMesh &mesh,
                           const mfem::Array<int> &bdr_attr_marker, CrackSides &sides);

}  // namespace mesh

namespace fem
{

// For local element e, return the local entity (bit index, see above) owning each element
// DOF in the native element DOF ordering of mfem::FiniteElementSpace::GetElementDofs, or -1
// for interior (bubble) DOFs.
void GetElementDofEntities(const mfem::FiniteElementSpace &fespace, int e,
                           mfem::Array<int> &entities);

}  // namespace fem

//
// Prolongation operator from the true DOFs of a broken space to its L-vector. The true DOF
// vector of the broken space is [x; x_s], with x the true DOFs of the underlying conforming
// space and x_s the copies of the split true DOFs owned by this process. The L-vector is
// [y; y_c], with y the L-vector of the conforming space and y_c the copied L-DOFs, so that
// y = P x and y_c = (P z)[copy], where z = x with the split true DOFs replaced by their
// copies x_s (constrained L-DOFs can depend on both split and unsplit true DOFs). Parallel
// communication is that of the conforming prolongation P.
//
class BrokenProlongation : public Operator
{
private:
  const Operator &P;
  const mfem::HypreParMatrix *hP;
  mfem::Array<int> copy_ldofs, split_tdofs;
  const int vsize, tsize;
  mutable Vector tx, lx;

  // Transpose of the copy rows: with w = Pᵀ x_c (or |P|ᵀ x_c), x_c scattered into an
  // L-vector, add a w to y_t for the unsplit true DOFs, and set y_s = a w[split] + b y_s.
  void AddCopyTranspose(const Vector &x, double a, double b, Vector &y, bool abs) const;

public:
  BrokenProlongation(const Operator &P, const mfem::Array<int> &copy_ldofs,
                     const mfem::Array<int> &split_tdofs);

  void Mult(const Vector &x, Vector &y) const override;

  void MultTranspose(const Vector &x, Vector &y) const override;

  // y = a |P_broken|ᵀ x + b y, used for diagonal assembly on nonconforming meshes.
  void AbsMultTranspose(double a, const Vector &x, double b, Vector &y) const;
};

}  // namespace palace

#endif  // PALACE_FEM_BROKEN_SPACE_HPP
