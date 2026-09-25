// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_METAL_EDGE_HPP
#define PALACE_UTILS_METAL_EDGE_HPP

#include <array>
#include <cstddef>
#include <functional>
#include <optional>
#include <vector>
#include <mfem.hpp>

namespace palace
{

namespace config
{
struct BoundaryData;
}  // namespace config

enum class MetalBoundaryConditionType : char
{
  PEC,
  CONDUCTIVITY,
  IMPEDANCE,
  RATIONAL_IMPEDANCE
};

struct MetalBoundaryCondition
{
  MetalBoundaryConditionType type;

  // Zero-based index into the corresponding BoundaryData vector. PEC-like boundaries
  // (including electrostatic terminals and prescribed potentials) use index zero.
  int index;
};

enum class MetalEdgeVertexType : char
{
  REGULAR,
  CORNER,
  ENDPOINT,
  JUNCTION
};

// Classification of a geometric metal perimeter edge by the distinct metal faces which
// support it (coincident crack copies and duplicate boundary elements count once):
// PHYSICAL / TRUNCATION edges are one-sided (metal on exactly one side within the face
// plane), a FOLD edge joins two non-coplanar metal faces (the metal turns around it: a
// staple, a box edge), a NONMANIFOLD edge is shared by three or more face directions (a
// wall standing on a sheet). Edges with metal continuing on both sides are not perimeter.
enum class MetalEdgeSegmentType : char
{
  PHYSICAL,
  TRUNCATION,
  FOLD,
  NONMANIFOLD,
  // A one-sided edge bordering a port boundary (LumpedPort / WavePort attribute of the
  // configuration): the port is not metal, the metal edge along it is a cut like a
  // truncation (decision 82(5)).
  PORT
};

struct MetalEdgeVertex
{
  std::array<double, 3> coordinate{};
  std::vector<std::size_t> segments;

  // Topology of the complete metal perimeter and of the physical-edge graph after
  // simulation-boundary truncation segments have been removed. The latter is empty for a
  // vertex which only belongs to truncation segments.
  MetalEdgeVertexType type = MetalEdgeVertexType::REGULAR;
  std::optional<MetalEdgeVertexType> physical_type;
  bool on_truncation_boundary = false;
  bool on_port_boundary = false;
};

struct MetalEdgeSegment
{
  std::array<std::size_t, 2> vertices{};
  int component = -1;
  int physical_component = -1;
  int physical_chain = -1;
  int metal_component = -1;
  MetalEdgeSegmentType type = MetalEdgeSegmentType::PHYSICAL;

  // Metal face attributes and boundary conditions which geometrically support this
  // perimeter segment.
  std::vector<int> metal_attributes;
  std::vector<MetalBoundaryCondition> conditions;

  // Existing dielectric postprocessing indices whose geometric perimeters coincide with
  // this metal edge. Empty or multiple classifications are retained for diagnostics.
  std::vector<int> sa_interfaces;
  std::vector<int> ms_interfaces;
  std::vector<int> ma_interfaces;

  // Nonmetal, non-interface boundary attributes which geometrically support this
  // segment. Such a segment is an artificial termination at a simulation cut surface,
  // rather than a fabricated metal edge.
  std::vector<int> truncation_attributes;
  // Port boundary attributes (LumpedPort / WavePort) whose faces border this segment.
  std::vector<int> port_attributes;

  // Distinct metal faces supporting the segment, their unit normals (sign-canonical), and
  // the domain element attributes adjacent to those faces over every coincident copy. A
  // sheet with the same material on both sides (an airbridge span, metal embedded in one
  // dielectric) has a single side attribute.
  int face_count = 0;
  std::vector<std::array<double, 3>> face_normals;
  std::vector<int> side_attributes;

  // Every supporting face lies on the bounding box of the mesh (a PEC simulation box).
  bool on_bounding_box = false;
};

struct MetalSurfaceFace
{
  // Ordered polygonal boundary of this surface facet. High-order boundary edges are
  // sampled consistently with the extracted metal perimeter.
  std::vector<std::array<double, 3>> vertices;
  int component = -1;

  // Global (deduplicated) faces only: metal attribute, sign-canonical unit normal, and
  // whether the face lies on the bounding box of the mesh (a PEC simulation box).
  int attribute = -1;
  std::array<double, 3> normal{};
  bool on_bounding_box = false;
};

struct MetalEdgeGeometry
{
  std::vector<MetalEdgeVertex> vertices;
  std::vector<MetalEdgeSegment> segments;

  // Rank-local facets (MetalSurfaceExtraction::retain_faces) and the global deduplicated
  // metal faces replicated on every rank (MetalSurfaceExtraction::retain_global_faces).
  std::vector<MetalSurfaceFace> surface_faces;
  std::vector<MetalSurfaceFace> global_faces;

  // Area-weighted principal direction of the metal face normals (sign-canonical): the
  // normal of the plane which carries most of the metal, i.e. the process plane.
  std::array<double, 3> layer_normal{};
  int components = 0;
  int physical_components = 0;
  int physical_chains = 0;

  // Connected components of the metal surface through shared face edges (coincident
  // crack copies are one face): the geometric conductor identity, independent of the
  // boundary attribute numbering.
  int metal_components = 0;

  bool Empty() const { return segments.empty(); }
};

struct MetalSurfaceExtraction
{
  // Metal components are always classified; this flag is retained for callers.
  bool classify_components = false;

  // Retain rank-local surface facets for exact spatial-neighborhood matching; callers
  // should exchange only clipped facets.
  bool retain_faces = false;

  // Retain the global deduplicated metal faces (replicated on every rank).
  bool retain_global_faces = false;
};

// Automatically extract the geometric perimeter of all PEC-like, conductivity, and
// impedance metal surfaces in a 3D mesh. The perimeter is derived from the distinct
// geometric metal faces (coincident crack copies and duplicate boundary elements of one
// attribute count once), so it does not depend on CrackInternalBoundaryElements or on the
// materials adjacent to a sheet. The result is replicated on every rank and classified
// using only existing boundary conditions and dielectric postprocessing surfaces;
// InterfaceDielectricData::edge_attributes is intentionally not used.
MetalEdgeGeometry ExtractMetalEdgeGeometry(const mfem::ParMesh &mesh,
                                           const config::BoundaryData &boundaries,
                                           MetalSurfaceExtraction surface = {});

// Infer a process normal for selected physical metal-edge segments from their supporting
// metal faces. The material score must increase from the process/substrate side toward
// the air/vacuum side. The optional fallback orients materially ambiguous surfaces (the
// same material on both sides); when given, ambiguous_side reports per segment whether the
// side had to be taken from the fallback or, without one, from the dominant component of
// the face normal (a guess the caller should not rely on).
std::vector<std::array<double, 3>> BuildMetalEdgeProcessNormals(
    const mfem::ParMesh &mesh, const MetalEdgeGeometry &geometry,
    const std::vector<std::size_t> &segment_indices,
    const std::function<double(int)> &material_score,
    const std::optional<std::array<double, 3>> &fallback = std::nullopt,
    std::vector<bool> *ambiguous_side = nullptr);

// Infer the in-plane direction from metal toward the adjacent gap for selected physical
// edge segments. The process normals must correspond one-to-one with segment_indices.
std::vector<std::array<double, 3>>
BuildMetalEdgeGapDirections(const mfem::ParMesh &mesh, const MetalEdgeGeometry &geometry,
                            const std::vector<std::size_t> &segment_indices,
                            const std::vector<std::array<double, 3>> &process_normals);

}  // namespace palace

#endif  // PALACE_UTILS_METAL_EDGE_HPP
