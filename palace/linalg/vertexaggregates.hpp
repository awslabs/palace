// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_VERTEX_AGGREGATES_HPP
#define PALACE_LINALG_VERTEX_AGGREGATES_HPP

#include <cstdint>
#include <vector>

namespace palace
{

class Mesh;

//
// Aggregates of the vertices of a mesh along the directions in which the elements are
// flat, for the blocks of the block Jacobi smoother.
//
// On an element which is much thinner in one direction than in the others, a second-order
// operator splits into a strong part, which penalizes the variation of a function across
// the element, and a weak part for the variation along it. A relaxation is robust in the
// aspect ratio of the elements only if it solves exactly for the functions which the
// strong part does not see, and these are not local to a vertex: across a stack of flat
// elements they are constant on the lines of vertices through the stack.
//
// The strong part is read from the element matrices of the Laplacian on the lowest order
// nodal space of the mesh, scaled by the assembled diagonal. An element is flat when the
// nonzero eigenvalues of its scaled matrix have a gap: the eigenvectors above the gap span
// the strong part, and the entries of the orthogonal projector onto them tie the vertices
// which the strong part couples. The aggregates are the connected components of the ties
// over all elements: lines of vertices through stacks of flat elements, sheets where flat
// elements are unstructured, and single vertices elsewhere. No length, angle or aspect
// ratio of the mesh is used.
//
struct VertexAggregates
{
  // Global identifier of the aggregate of each local mesh vertex.
  std::vector<std::uint64_t> aggregate;

  // Whether that aggregate has more than one vertex.
  std::vector<char> nontrivial;
  // Whether any aggregate of the mesh, on any process, has more than one vertex: without,
  // the mesh has no flat elements and the smoothers have no blocks.
  bool any = false;

  // The ties which the elements of this process found: pairs of local vertices with the
  // entry of the projector onto the strong part of the element. The strong part penalizes
  // the difference of a function between two vertices on opposite sides of a flat element,
  // whose entry is negative; two vertices on the same side move together in it, and their
  // entry is positive. An aggregate which is too large for one block is split along the
  // pairs which are not tied across the flat direction (see linalg/aggregateblocks.hpp).
  std::vector<int> tie_a, tie_b;
  std::vector<double> tie_w;
};

// The vertex aggregates of a mesh, computed on first use and cached with the mesh.
// Collective.
const VertexAggregates &GetVertexAggregates(const Mesh &mesh);

}  // namespace palace

#endif  // PALACE_LINALG_VERTEX_AGGREGATES_HPP
