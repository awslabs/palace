// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_EDGE_DISTANCE_HPP
#define PALACE_UTILS_EDGE_DISTANCE_HPP

#include <array>
#include <limits>
#include <memory>
#include <optional>
#include <vector>
#include <mfem.hpp>
#include "utils/geodata.hpp"
#include "utils/metaledge.hpp"

namespace palace
{

namespace config
{
struct BoundaryData;
}  // namespace config

enum class InterfaceDielectric : char;

class EdgeDistanceTree
{
public:
  struct NearestResult
  {
    double distance_squared;
    std::size_t segment;
  };

  struct SegmentMetadata
  {
    int component = -1;
    int chain = -1;
    std::array<int, 2> vertex_types{-1, -1};
    std::array<double, 2> vertex_distances{std::numeric_limits<double>::infinity(),
                                           std::numeric_limits<double>::infinity()};
  };

private:
  static constexpr std::size_t leaf_size = 8;

  struct Node
  {
    std::array<double, 3> min;
    std::array<double, 3> max;
    std::size_t begin = 0;
    std::size_t end = 0;
    int left = -1;
    int right = -1;

    bool IsLeaf() const { return left < 0; }
  };

  std::vector<mesh::BoundaryEdgeSegment> segments;
  std::vector<std::array<double, 3>> process_normals;
  std::vector<SegmentMetadata> metadata;
  std::vector<std::size_t> indices;
  std::vector<Node> nodes;

  static double DistanceSquaredToSegment(const mfem::Vector &point,
                                         const mesh::BoundaryEdgeSegment &segment);
  static double DistanceSquaredToBox(const mfem::Vector &point, const Node &node);
  int Build(std::size_t begin, std::size_t end);
  NearestResult Search(int node_index, const mfem::Vector &point, NearestResult best) const;

public:
  explicit EdgeDistanceTree(std::vector<mesh::BoundaryEdgeSegment> segments,
                            std::vector<std::array<double, 3>> process_normals = {},
                            std::vector<SegmentMetadata> metadata = {});

  double DistanceSquared(const mfem::Vector &point) const;
  NearestResult Nearest(const mfem::Vector &point) const;
  // Distance to ONE segment of the tree (the per-segment test of the sagitta-aware
  // exclusion, ExcludeCoincidentEdgeSegments).
  double DistanceSquaredToSegment(const mfem::Vector &point, std::size_t index) const
  {
    return DistanceSquaredToSegment(point, segments[index]);
  }
  double DistanceAlongEdgeToNonregularVertex(const mfem::Vector &point,
                                             std::size_t segment) const;
  std::size_t Size() const { return segments.size(); }
  const mesh::BoundaryEdgeSegment &GetSegment(std::size_t index) const
  {
    return segments[index];
  }
  bool HasProcessNormals() const { return process_normals.size() == segments.size(); }
  const std::array<double, 3> &GetProcessNormal(std::size_t index) const
  {
    MFEM_ASSERT(HasProcessNormals(), "Edge-distance tree has no process normals!");
    return process_normals[index];
  }
  bool HasMetadata() const { return metadata.size() == segments.size(); }
  const SegmentMetadata &GetMetadata(std::size_t index) const
  {
    MFEM_ASSERT(HasMetadata(), "Edge-distance tree has no metal-edge metadata!");
    return metadata[index];
  }
};

struct EdgeFrame
{
  std::array<double, 3> normal;
  std::array<double, 3> transverse;
  std::array<double, 3> tangent;
};

EdgeFrame BuildEdgeFrame(const mesh::BoundaryEdgeSegment &segment,
                         const std::array<double, 3> &configured_normal, int dimension);

std::array<double, 6> GetPolarizedEdgeEnergyDensity(
    const mfem::Vector &point, const mesh::BoundaryEdgeSegment &segment,
    const EdgeFrame &frame, const mfem::Vector &field, const mfem::Vector &displacement);

// Explicit segment exclusion (EdgeExcludeSegments, decision 282 rule B5): a perimeter
// segment whose two ends and midpoint each lie within `tolerance` + the listed segment's
// chord sagitta (`sagittas`, parallel to `exclude_segments`; empty = every sagitta 0) of
// one listed exclusion segment [x0, y0, z0, x1, y1, z1] is dropped from the edge-distance
// perimeter (the decision-611 rule: a mesh perimeter following the TRUE arc of a listed
// chord lies up to the sagitta off it; decision 657 (a)). An empty list leaves the
// perimeter unchanged; removing the whole perimeter is an error. `record`, when given,
// receives the exclusion's coverage record (the listed length, the perimeter length
// actually removed, and the KEPT perimeter length sharing a vertex with a removed segment:
// the segments straddling the ends of the removed runs, at most one per end of a listed
// piece when the piece is removed in one run; decision 665's end allowance, since a listed
// end falling inside a perimeter segment keeps that segment whole).
struct EdgeExclusionRecord
{
  std::size_t exclude_segments_count = 0;
  double exclude_segments_length = 0.0;
  std::size_t perimeter_count = 0;
  double perimeter_length = 0.0;
  std::size_t excluded_perimeter_count = 0;
  double excluded_perimeter_length = 0.0;
  double straddling_perimeter_length = 0.0;
  std::size_t removed_runs = 0;
  std::size_t listed_chains = 0;
  double tolerance = 0.0;
  double max_sagitta = 0.0;
};

// The end allowance of decision 667 (1) over a classified perimeter (`excluded` parallels
// `segments`): straddling_perimeter_length = the length of the kept perimeter segments that
// contain a listed chain's end point within the tolerance (at most one kept segment per
// listed end, each kept segment credited once), removed_runs = the maximal runs of removed
// perimeter segments chained by shared vertices, listed_chains = the chains of listed
// segments joined end to end within the tolerance. Public so the unit tests read it on
// exclusion patterns the coincidence test cannot produce.
void RecordEdgeExclusionEnds(const std::vector<mesh::BoundaryEdgeSegment> &segments,
                             const std::vector<bool> &excluded,
                             const std::vector<std::array<double, 6>> &exclude_segments,
                             double tolerance, EdgeExclusionRecord &record);

std::shared_ptr<const EdgeDistanceTree> BuildEdgeDistanceTree(
    const mfem::ParMesh &mesh, const std::vector<int> &edge_attributes,
    const std::vector<int> &edge_exclude_attributes,
    const std::optional<std::array<double, 3>> &process_normal = std::nullopt,
    const std::vector<std::array<double, 6>> &edge_exclude_segments = {},
    double edge_exclude_segment_tolerance = 0.0,
    const std::vector<double> &edge_exclude_segment_sagittas = {},
    EdgeExclusionRecord *record = nullptr);

void ExcludeCoincidentEdgeSegments(
    const std::vector<std::array<double, 6>> &exclude_segments,
    const std::vector<double> &sagittas, double tolerance, int space_dimension,
    std::vector<mesh::BoundaryEdgeSegment> &edge_segments,
    EdgeExclusionRecord *record = nullptr);

void ExcludeCoincidentMetalEdgeSegmentIndices(
    const MetalEdgeGeometry &geometry,
    const std::vector<std::array<double, 6>> &exclude_segments,
    const std::vector<double> &sagittas, double tolerance, int space_dimension,
    std::vector<std::size_t> &segment_indices, EdgeExclusionRecord *record = nullptr);

std::vector<std::size_t>
GetInterfaceMetalEdgeSegmentIndices(const MetalEdgeGeometry &geometry, int interface_index,
                                    InterfaceDielectric interface_type);

void ExcludeMetalEdgeSegmentIndices(const mfem::ParMesh &mesh,
                                    const MetalEdgeGeometry &geometry,
                                    const std::vector<int> &exclude_attributes,
                                    std::vector<std::size_t> &segment_indices);

std::shared_ptr<const EdgeDistanceTree>
BuildEdgeDistanceTree(const MetalEdgeGeometry &geometry,
                      const std::vector<std::size_t> &segment_indices,
                      std::vector<std::array<double, 3>> process_normals = {});

struct EdgeRefinementContext
{
  std::shared_ptr<const EdgeDistanceTree> distance_tree;
  double radius;
  double target_size;
  double outer_radius;
  double core_indicator_weight;
};

std::vector<EdgeRefinementContext>
BuildEdgeRefinementContexts(const mfem::ParMesh &mesh,
                            const config::BoundaryData &boundaries);

mfem::Array<int>
MarkEdgeRefinementElements(const mfem::ParMesh &mesh,
                           const std::vector<EdgeRefinementContext> &contexts);

int WeightEdgeCoreIndicators(const mfem::ParMesh &mesh,
                             const std::vector<EdgeRefinementContext> &contexts,
                             mfem::Vector &indicators);

}  // namespace palace

#endif  // PALACE_UTILS_EDGE_DISTANCE_HPP
