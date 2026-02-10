// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_MESH_TOPOLOGY_HPP
#define PALACE_FEM_MESH_TOPOLOGY_HPP

#include <array>
#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>
#include <mpi.h>

namespace mfem
{
class Mesh;
}  // namespace mfem

namespace palace
{

//
// Prototype mesh topology for simplex elements with conformal bisection refinement.
// Stores elements, vertices, and boundary elements with single-source attributes.
// Maintains an append-only refinement history enabling conformal coarsening.
//
// Current scope: serial, tet/tri only (simplex elements).
//
class MeshTopology
{
public:
  struct Element
  {
    std::vector<int> vertices;
    int attribute;
    int parent;   // -1 if original (from initial mesh)
    int level;    // refinement level (0 = original)
    bool active;  // false if refined into children
  };

  struct BdrElement
  {
    std::vector<int> vertices;
    int attribute;
  };

  struct RefinementRecord
  {
    int parent;
    std::array<int, 2> children;
    // The bisected edge (global vertex indices) and the midpoint vertex.
    int edge_v0, edge_v1;
    int midpoint;
  };

private:
  int dim;   // topological dimension (2 = tri, 3 = tet)
  int sdim;  // space dimension

  std::vector<std::array<double, 3>> vertices;
  std::vector<Element> elements;
  std::vector<BdrElement> bdr_elements;
  std::vector<RefinementRecord> history;

  // Edge midpoint cache: maps sorted (v0, v1) pair to midpoint vertex index.
  // Ensures shared edges produce a single midpoint.
  std::unordered_map<int64_t, int> edge_midpoints;

  static int64_t EdgeKey(int v0, int v1)
  {
    if (v0 > v1)
    {
      std::swap(v0, v1);
    }
    return (static_cast<int64_t>(v0) << 32) | static_cast<int64_t>(v1);
  }

  // Get or create the midpoint of an edge.
  int GetMidpoint(int v0, int v1);

  // Find the longest edge of an element, returned as (local_idx_a, local_idx_b).
  std::pair<int, int> FindLongestEdge(int elem_idx) const;

  // Bisect a single tet/tri element. When force_va/force_vb are non-negative,
  // bisect along that edge instead of the longest edge (used for closure).
  std::array<int, 2> BisectElement(int elem_idx, int force_va = -1, int force_vb = -1);

  // Update boundary elements after refinement of a set of elements.
  void UpdateBdrElements(const std::vector<int> &refined_elements);

public:
  // Import from an MFEM serial mesh (simplex elements only).
  static MeshTopology FromMFEM(const mfem::Mesh &mesh);

  // Export to an MFEM serial mesh.
  std::unique_ptr<mfem::Mesh> ToMFEM() const;

  // Uniform refinement: bisect all active elements.
  void UniformRefinement();

  // Adaptive refinement: bisect marked elements (global indices) with closure to
  // maintain conformity. The result depends only on the mesh topology and the marked
  // element set, not on how the mesh is partitioned across ranks.
  void Refine(const std::vector<int> &marked_elements);

  // Distributed adaptive refinement: each rank provides local element marks using a
  // local-to-global map. Marks are gathered globally, then refinement + closure runs
  // on the global topology. comm is used for the Allgather.
  void RefineDistributed(MPI_Comm comm, const std::vector<int> &local_marks,
                         const std::vector<int> &local_to_global_elem);

  // Coarsen: merge sibling pairs whose combined error is below threshold.
  void Coarsen(const std::vector<double> &elem_error, double threshold);

  // Queries.
  int Dimension() const { return dim; }
  int SpaceDimension() const { return sdim; }
  int GetNE() const;  // count of active elements
  int GetNV() const { return static_cast<int>(vertices.size()); }
  int GetNBE() const { return static_cast<int>(bdr_elements.size()); }
  int GetTotalElements() const { return static_cast<int>(elements.size()); }
  const std::vector<RefinementRecord> &GetHistory() const { return history; }
  const Element &GetElement(int i) const { return elements[i]; }
  int GetActiveElementCount() const;
};

}  // namespace palace

#endif  // PALACE_FEM_MESH_TOPOLOGY_HPP
