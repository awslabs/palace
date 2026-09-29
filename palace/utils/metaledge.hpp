// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_METAL_EDGE_HPP
#define PALACE_UTILS_METAL_EDGE_HPP

#include <array>
#include <cmath>
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

// Geometric joint noise threshold (USER decision 121 (B), 2026-09-28; replaces the 1 deg
// angular threshold of decision 117(4) and the 30 deg corner class of decision 73). A
// perimeter vertex with two segments, turning by t between two straight PIECES (collinear
// mesh segments merged: refinement cannot change the pieces), is a straight continuation —
// a REGULAR joint of its chain — when the deviation from straight it implies at the
// resolution of the correction is below kJointNoiseSagittaOverRadius x R: the implied
// sagitta (c / 2) tan(t / 4) of a chord of length c = the SHORTER adjacent piece read as one
// chord of a circle turning t per chord (sagitta = rho (1 - cos(t / 2)) with c = 2 rho
// sin(t / 2)). Every other joint is a CORNER unless the identification's arc rule absorbs
// it (>= 3 concyclic joints, each turning less than ArcMaxJointTurnDegrees). Why this
// quantity: (i) it is the very quantity the arc rule records per chord
// (Arcs[].MaxChordSagittaOverR, warning above the same 0.05 R), so a joint below it is a
// chord joint of a curve the correction cannot resolve — one constant for the noise
// threshold and the mesh-coarseness diagnostic; (ii) sub-nm mesh slivers (1-2 nm segments
// on DS-SCT-002's flux lines whose roundoff directions turned by > 1 deg and fragmented
// the stacks) are bounded by c / 2 whatever their turn, a U-turn included (tan(pi / 4) =
// 1); (iii) spline steps of 1-6 deg on 1-5 um chords (0.5-2.6 R) imply 2-65 nm =
// 0.001-0.034 R and read as the smooth curves they discretise (the alternative, the
// vertex distance from the line through its neighbours c sin(t), is 4x larger at small t
// and would keep 6 deg spline joints on 1 um chords as 174 deg corners); (iv) a real
// corner between long arms is never noise: 5 um arms turning 20 deg imply 0.23 R. Taking
// the SHORTER piece makes a kink next to a short mesh segment on a long straight edge a
// property of the piece pair, not of the mesh segment (a refinement midpoint is collinear
// and merges). The knife edge is on the implied sagitta at 0.05 R (Conventions
// JointNoiseSagittaOverR; knife-edge census JointNoiseSagittaOverR.0.05). Recorded as
// Identification.Conventions.JointNoiseSagittaOverR; the extraction receives it in mesh
// units (MetalSurfaceExtraction::joint_noise_sagitta = this constant x R).
constexpr double kJointNoiseSagittaOverRadius = 0.05;

// The implied sagitta of a joint turning by turn_radians between two straight pieces whose
// shorter one has length shorter_piece (mesh units): (c / 2) tan(t / 4).
inline double ImpliedJointSagitta(double turn_radians, double shorter_piece)
{
  return 0.5 * shorter_piece * std::tan(0.25 * turn_radians);
}

// The joint noise rule on a fixed relative grid (1e-9 of the threshold), so that a
// roundoff-level perturbation of the vertex coordinates cannot flip a joint between REGULAR
// and CORNER: true when the implied sagitta is below the threshold.
inline bool JointIsNoise(double turn_radians, double shorter_piece, double noise_sagitta)
{
  const double grid = 1.0e-9 * noise_sagitta;
  return std::round(ImpliedJointSagitta(turn_radians, shorter_piece) / grid) <
         std::round(noise_sagitta / grid);
}

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

  // Joint rule for two-segment vertices (CORNER and chain break, or REGULAR). Geometric
  // (USER decision 121 (B)): REGULAR iff JointIsNoise(turn, shorter adjacent straight
  // piece, joint_noise_sagitta), joint_noise_sagitta = kJointNoiseSagittaOverRadius x R in
  // mesh units (the identification passes the library's matching radius, the
  // post-processing edge tools their largest edge distance). Angular (the legacy per-group
  // classifier, comparison only): REGULAR iff turn <= corner_turn_tolerance_degrees.
  // Exactly one of the two must be set.
  double joint_noise_sagitta = 0.0;
  std::optional<double> corner_turn_tolerance_degrees;
};

// The extraction options with the geometric joint noise rule at matching radius R (mesh
// units).
inline MetalSurfaceExtraction JointNoiseExtraction(double radius)
{
  MetalSurfaceExtraction surface;
  surface.joint_noise_sagitta = kJointNoiseSagittaOverRadius * radius;
  return surface;
}

// The same at the length scale of the post-processing edge tools: the largest edge distance
// of the automatic-edge interface dielectrics (the correction path requires the largest
// EdgeDistances value to equal the library's matching radius). Fails when no automatic-edge
// interface configures edge distances.
MetalSurfaceExtraction JointNoiseExtractionFor(const config::BoundaryData &boundaries);

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
