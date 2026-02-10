// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_MESH_TOPOLOGY_HPP
#define PALACE_FEM_MESH_TOPOLOGY_HPP

#include <array>
#include <cstdint>
#include <map>
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
// Mesh topology for simplex and hexahedral elements with refinement support.
// Simplex elements (tet/tri) use conformal bisection refinement.
// Hexahedral elements (hex/quad) use isotropic refinement (1 hex → 8 children).
// Stores elements, vertices, and boundary elements with single-source attributes.
// Maintains an append-only refinement history enabling coarsening.
//
// Current scope: serial, tet/tri and hex/quad elements.
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
    std::vector<int> children;  // 2 for bisection, 8 for hex refinement
    // The bisected edge (global vertex indices) and the midpoint vertex.
    // Set to -1 for hex refinement (which uses face/body centers instead).
    int edge_v0, edge_v1;
    int midpoint;
  };

private:
  int dim;   // topological dimension (2 = tri/quad, 3 = tet/hex)
  int sdim;  // space dimension

  std::vector<std::array<double, 3>> vertices;
  std::vector<Element> elements;
  std::vector<BdrElement> bdr_elements;
  std::vector<RefinementRecord> history;

  // Edge midpoint cache: maps sorted (v0, v1) pair to midpoint vertex index.
  // Ensures shared edges produce a single midpoint.
  std::unordered_map<int64_t, int> edge_midpoints;

  // Face center cache: maps sorted (v0, v1, v2, v3) to center vertex index.
  // Used for hex/quad refinement where quad faces need center vertices.
  std::map<std::array<int, 4>, int> face_centers;

  // Body center cache: maps sorted 8-vertex key to center vertex index.
  // Used for hex refinement where the hex body needs a center vertex.
  std::map<std::array<int, 8>, int> body_centers;

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

  // Get or create the center vertex of a quad face (4 vertices).
  int GetFaceCenter(int v0, int v1, int v2, int v3);

  // Get or create the body center vertex of a hex element (8 vertices).
  int GetBodyCenter(const std::vector<int> &hex_vertices);

  // Refine a hex element isotropically into 8 children. Returns child indices.
  std::array<int, 8> RefineHexElement(int elem_idx);

  // Refine a quad element isotropically into 4 children. Returns child indices.
  std::array<int, 4> RefineQuadElement(int elem_idx);

  // Update boundary elements after refinement of a set of elements.
  void UpdateBdrElements(const std::vector<int> &refined_elements);

public:
  // Import from an MFEM serial mesh (simplex and hexahedral elements).
  static MeshTopology FromMFEM(const mfem::Mesh &mesh);

  // Export to an MFEM serial mesh.
  std::unique_ptr<mfem::Mesh> ToMFEM() const;

  // Uniform refinement: bisect all active simplex elements, isotropically refine hex/quad.
  void UniformRefinement();

  // Adaptive refinement: bisect marked simplex elements with closure for conformity,
  // or isotropically refine marked hex/quad elements (no closure needed).
  // The result depends only on the mesh topology and the marked element set,
  // not on how the mesh is partitioned across ranks.
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
