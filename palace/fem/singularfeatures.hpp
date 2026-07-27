// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_SINGULARFEATURES_HPP
#define PALACE_FEM_SINGULARFEATURES_HPP

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>
#include <mfem.hpp>

namespace palace
{

namespace fem
{

namespace singular
{

using GlobalVertexId = std::int64_t;

enum class FeatureVertexType
{
  REGULAR,
  CORNER,
  ENDPOINT
};

struct FeatureSegment
{
  int mesh_edge;
  std::array<int, 2> mesh_vertices;
  std::size_t feature;
  std::vector<int> boundary_attributes;
};

struct FeatureVertex
{
  std::size_t id;
  int mesh_vertex;
  FeatureVertexType type;
  std::vector<std::size_t> segments;
  std::vector<std::size_t> features;
  double nu;
};

struct StraightFeature
{
  std::size_t id;
  std::vector<std::size_t> segments;
  std::vector<int> mesh_vertices;
  double nu;
  bool closed;
};

struct SheetFace
{
  // Source-serial mesh face ID and sorted source-serial vertex IDs.
  int mesh_face;
  std::array<GlobalVertexId, 3> mesh_vertices;
  int boundary_attribute;
};

struct ElementNodeFeature
{
  // Index of the unique singular-node record in FeatureTopology::vertices.
  // A corner shared by several straight features still has only one node set,
  // as in Fig. 2(c) of Elkin et al.
  std::size_t vertex;
  int mesh_vertex;

  // Element-local barycentric node indices. Entry zero is the singular node;
  // the remaining entries are ordered by mesh vertex ID.
  std::array<int, 4> canonical_nodes;
};

struct ElementEdgeFeature
{
  std::size_t feature;
  std::size_t segment;
  int mesh_edge;

  // Element-local barycentric node indices. Entries zero and one are the
  // singular edge ordered by mesh vertex ID; the remaining entries are ordered
  // by mesh vertex ID.
  std::array<int, 4> canonical_nodes;
};

struct ElementFeatureIncidence
{
  std::vector<ElementNodeFeature> nodes;
  std::vector<ElementEdgeFeature> edges;
};

struct FeatureTopology
{
  // The records below form a decomposition-independent blueprint. Their mesh
  // vertex and edge IDs refer to the source serial mesh. Element incidence
  // instead refers to the mesh on which it is used, which is the source mesh
  // for serial extraction and the rank-local mesh after distribution.
  std::vector<FeatureVertex> vertices;
  std::vector<FeatureSegment> segments;
  std::vector<StraightFeature> features;
  std::vector<SheetFace> sheet_faces;
  std::vector<ElementFeatureIncidence> elements;

  bool Empty() const { return features.empty(); }
};

struct TriangleSelectedSegment
{
  int boundary_element;
  int mesh_edge;
  std::array<int, 2> mesh_vertices;
  int boundary_attribute;
};

struct TriangleTipVertex
{
  std::size_t id;
  int mesh_vertex;
  std::vector<std::size_t> selected_segments;
  double nu;
};

struct TriangleElementNodeFeature
{
  std::size_t vertex;
  int mesh_vertex;

  // Element-local barycentric node indices. Entry zero is the singular node;
  // the other two entries are ordered by source mesh vertex ID.
  std::array<int, 3> canonical_nodes;
};

struct TriangleElementFeatureIncidence
{
  std::vector<TriangleElementNodeFeature> nodes;
};

struct TriangleFeatureTopology
{
  std::vector<TriangleTipVertex> vertices;
  std::vector<TriangleSelectedSegment> selected_segments;
  std::vector<TriangleElementFeatureIncidence> elements;

  bool Empty() const { return vertices.empty(); }
};

// Extract the topological perimeter of selected zero-thickness internal PEC
// boundary attributes in a serial, conforming tetrahedral mesh. Curved
// elements away from and adjoining the sheet are accepted, but every selected
// perimeter edge must be geometrically straight.
//
// This phase deliberately uses serial mesh vertex and edge IDs as stable entity
// IDs. A ParMesh with more than one rank is rejected. Parallel extraction must
// first provide decomposition-independent global entity IDs and is part of the
// enrichment DOF-manager phase.
FeatureTopology ExtractSerialSheetFeatures(const mfem::Mesh &mesh,
                                           const std::vector<int> &boundary_attributes,
                                           double nu = 0.5);

// Extract endpoint singularities of selected straight, zero-thickness internal
// PEC line chains in a serial, conforming triangular mesh. Curved elements are
// accepted, but selected PEC segments must be geometrically straight.
// Endpoints on the exterior mesh boundary are excluded. Bends, branches,
// selected exterior edges, and non-triangular meshes are rejected until their
// wedge exponents and continuity contracts are implemented.
TriangleFeatureTopology
ExtractSerialLineTipFeatures(const mfem::Mesh &mesh,
                             const std::vector<int> &boundary_attributes, double nu = 0.5);

// Reconstruct rank-local feature incidence from a serial feature blueprint.
// serial_vertex_ids maps every rank-local mesh vertex to its source serial
// vertex ID. Global feature, segment, and singular-node IDs are preserved.
//
// Canonical barycentric tuples are ordered by source serial vertex ID, never
// by rank-local mesh numbering. The map is invalid after refinement or
// rebalancing unless it is explicitly propagated through those operations.
FeatureTopology
DistributeSerialSheetFeatures(const mfem::Mesh &serial_mesh,
                              const FeatureTopology &serial_features,
                              const mfem::ParMesh &parallel_mesh,
                              const std::vector<GlobalVertexId> &serial_vertex_ids);

// Broadcast an immutable serial feature blueprint from rank zero. Element
// incidence is serialized sparsely, so elements not touching a singular
// feature incur no per-element payload beyond the global element count.
void BroadcastSerialSheetFeatures(FeatureTopology &serial_features, MPI_Comm comm);

// Broadcast and distribute the immutable 2D line-tip blueprint. Selected
// segments retain source-serial IDs; only element incidence is reconstructed
// with rank-local mesh vertex numbers.
void BroadcastSerialLineTipFeatures(TriangleFeatureTopology &serial_features,
                                    MPI_Comm comm);

// Reconstruct rank-local feature incidence from exact source-serial entity
// maps transported by mesh::Partition. source_element_ids maps local element
// ordering to the element ordering used by serial_features.elements.
FeatureTopology
DistributeSerialSheetFeatures(const FeatureTopology &serial_features,
                              const mfem::ParMesh &parallel_mesh,
                              const std::vector<GlobalVertexId> &serial_vertex_ids,
                              const std::vector<GlobalVertexId> &source_element_ids);

TriangleFeatureTopology
DistributeSerialLineTipFeatures(const TriangleFeatureTopology &serial_features,
                                const mfem::ParMesh &parallel_mesh,
                                const std::vector<GlobalVertexId> &serial_vertex_ids,
                                const std::vector<GlobalVertexId> &source_element_ids);

}  // namespace singular

}  // namespace fem

}  // namespace palace

#endif  // PALACE_FEM_SINGULARFEATURES_HPP
