// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_EDGE_DISTANCE_HPP
#define PALACE_UTILS_EDGE_DISTANCE_HPP

#include <array>
#include <limits>
#include <memory>
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

std::shared_ptr<const EdgeDistanceTree>
BuildEdgeDistanceTree(const mfem::ParMesh &mesh, const std::vector<int> &edge_attributes,
                      const std::vector<int> &edge_exclude_attributes);

std::vector<std::size_t>
GetInterfaceMetalEdgeSegmentIndices(const MetalEdgeGeometry &geometry, int interface_index,
                                    InterfaceDielectric interface_type);

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
