// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_BROKEN_SPACE_HPP
#define PALACE_FEM_BROKEN_SPACE_HPP

#include <array>
#include <cstdint>
#include <memory>
#include <vector>
#include <mfem.hpp>

namespace palace
{

//
// Support for finite element spaces which are discontinuous ("broken") across selected
// interior boundaries without modifying the mesh.
//
// An interior boundary ("crack") is the set of interior faces carrying one of a given list
// of boundary attributes. Every mesh entity (vertex, edge, face) lying on such a face is
// split into one version per group of neighboring elements which are connected to each
// other without crossing the interior boundary, unless there is only one such group (for
// example at the free edge of an interior boundary). The group with the smallest label
// reads the original degrees of freedom (version 0), and each other group its own copy
// (versions 1, 2, ...), also at junctions of three or more groups, as on a cracked mesh.
//
// On a nonconforming mesh, the degrees of freedom of a hanging entity are constrained by
// those of the entities of the closure of its master entity, possibly split, which a broken
// space has to read in the version of the side of the element reading the hanging entity.
// Each element therefore also carries a side label, the global number of its ancestor on
// the mesh on which the sides were computed: all elements with the same label are on the
// same side of every split entity of the closure of the ancestor, which contains the
// master entities of their hanging entities.
//
// Unlike on a cracked mesh, an interior boundary face refined on one side only remains a
// master face for the hanging entities on the refined side, which read their own version of
// its DOFs but are still constrained by them: the recovery on the refined side is limited
// along the interior boundary to the trace of the coarse face, which makes the error
// estimate larger (more conservative) there.
//
// The sides are stored per local element, with bit (or 4-bit version) b for the local
// entity b of the element, where local entities are numbered as vertices [0, nv), then
// edges [nv, nv + ne), then (3D only) faces [nv + ne, nv + ne + nf), in the local ordering
// of mfem::Mesh::GetElementVertices, GetElementEdges, and GetElementFaces.
//
struct CrackSides
{
  // Bit b of split[e] is set when the local entity b of element e is split.
  std::vector<std::uint32_t> split;

  // Version of the DOFs of each split local entity read by the element (0 for the original
  // DOFs, k > 0 for the k-th copy), 4 bits per local entity.
  std::vector<std::array<std::uint32_t, 4>> version;

  // Dimension of the entity of the mesh on which the sides were computed which carries each
  // split local entity of the element (the entity itself, or the one it was refined from),
  // 2 bits per local entity. The DOFs of a hanging entity are constrained by those of the
  // entities of the closure of its master entity, whose versions are those of its carrier
  // for its own DOFs, but can differ for the entities of lower dimension of its closure.
  std::vector<std::uint64_t> carrier;

  // Side label of each element (the global number of its ancestor on the mesh on which the
  // sides were computed).
  std::vector<std::int64_t> side;

  static constexpr int max_version = 15;

  int GetCarrier(int e, int b) const
  {
    return static_cast<int>((carrier[e] >> (2 * b)) & 0x3u);
  }
  void SetCarrier(int e, int b, int d)
  {
    carrier[e] = (carrier[e] & ~(std::uint64_t(0x3) << (2 * b))) |
                 (static_cast<std::uint64_t>(d) << (2 * b));
  }

  int GetVersion(int e, int b) const
  {
    return static_cast<int>((version[e][b / 8] >> (4 * (b % 8))) & 0xFu);
  }
  void SetVersion(int e, int b, int v)
  {
    auto &word = version[e][b / 8];
    word = (word & ~(std::uint32_t(0xF) << (4 * (b % 8)))) |
           (static_cast<std::uint32_t>(v) << (4 * (b % 8)));
  }

  // Bitmask of the local entities of element e read as a copy (version > 0).
  std::uint32_t CopyBits(int e) const;

  // Resize for the given number of local elements, with no split entities and zero labels.
  void Reset(std::size_t ne);

  // Return whether any local entity is split.
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

// Prolongation of a broken space as a HypreParMatrix, with the global offsets of its rows
// and columns (in the format of mfem::ParFiniteElementSpace::GetDofOffsets) and the column
// map of its off-diagonal block, which the matrix references.
struct BrokenProlongationMatrix
{
  mfem::Array<HYPRE_BigInt> row_offsets, col_offsets;
  std::vector<HYPRE_BigInt> col_map;
  std::unique_ptr<mfem::HypreParMatrix> P;
};

// Build the prolongation of a broken space from the prolongation P of the underlying space
// (from its true DOFs to its L-vector). The true DOFs of each process are its true DOFs
// followed by the copies of its true DOFs with several versions (num_versions[t] - 1 copies
// of the local true DOF t, contiguous), and the L-vector is the L-vector of the underlying
// space followed by the copied L-DOFs. The row of the copied L-DOF c is the row of its
// original L-DOF copy_ldofs[c], with the version given for each of its entries (in the
// order of the diagonal then off-diagonal entries of the row of P) in
// copy_versions[copy_offsets[c], copy_offsets[c + 1]), the L-DOFs of constrained entities
// depending on several true DOFs, possibly in different versions. Collective.
void BuildBrokenProlongation(const mfem::ParMesh &mesh, mfem::HypreParMatrix &P,
                             const std::vector<int> &copy_ldofs,
                             const std::vector<int> &copy_offsets,
                             const std::vector<std::uint8_t> &copy_versions,
                             const std::vector<int> &num_versions,
                             BrokenProlongationMatrix &out);

}  // namespace fem

}  // namespace palace

#endif  // PALACE_FEM_BROKEN_SPACE_HPP
