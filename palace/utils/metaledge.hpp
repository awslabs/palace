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

enum class MetalEdgeSegmentType : char
{
  PHYSICAL,
  TRUNCATION
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
};

struct MetalEdgeSegment
{
  std::array<std::size_t, 2> vertices{};
  int component = -1;
  int physical_component = -1;
  int physical_chain = -1;
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
};

struct MetalEdgeGeometry
{
  std::vector<MetalEdgeVertex> vertices;
  std::vector<MetalEdgeSegment> segments;
  int components = 0;
  int physical_components = 0;
  int physical_chains = 0;

  bool Empty() const { return segments.empty(); }
};

// Automatically extract the physical perimeter of all PEC-like, conductivity, and
// impedance metal surfaces in a 3D mesh. The result is replicated on every rank and
// classified using only existing boundary conditions and dielectric postprocessing
// surfaces; InterfaceDielectricData::edge_attributes is intentionally not used.
MetalEdgeGeometry ExtractMetalEdgeGeometry(const mfem::ParMesh &mesh,
                                           const config::BoundaryData &boundaries);

// Infer a process normal for selected physical metal-edge segments from their supporting
// metal faces. The material score must increase from the process/substrate side toward
// the air/vacuum side. The optional fallback orients materially ambiguous surfaces.
std::vector<std::array<double, 3>> BuildMetalEdgeProcessNormals(
    const mfem::ParMesh &mesh, const MetalEdgeGeometry &geometry,
    const std::vector<std::size_t> &segment_indices,
    const std::function<double(int)> &material_score,
    const std::optional<std::array<double, 3>> &fallback = std::nullopt);

// Infer the in-plane direction from metal toward the adjacent gap for selected physical
// edge segments. The process normals must correspond one-to-one with segment_indices.
std::vector<std::array<double, 3>> BuildMetalEdgeGapDirections(
    const mfem::ParMesh &mesh, const MetalEdgeGeometry &geometry,
    const std::vector<std::size_t> &segment_indices,
    const std::vector<std::array<double, 3>> &process_normals);

}  // namespace palace

#endif  // PALACE_UTILS_METAL_EDGE_HPP
