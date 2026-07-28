// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "edgedistance.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <queue>
#include <tuple>
#include "utils/configfile.hpp"
#include "utils/diagnostic.hpp"
#include "utils/metaledge.hpp"

namespace palace
{

double EdgeDistanceTree::DistanceSquaredToSegment(const mfem::Vector &point,
                                                  const mesh::BoundaryEdgeSegment &segment)
{
  double direction[3], offset[3];
  double length_squared = 0.0;
  double projection = 0.0;
  for (int d = 0; d < point.Size(); d++)
  {
    direction[d] = segment.p1[d] - segment.p0[d];
    offset[d] = point[d] - segment.p0[d];
    length_squared += direction[d] * direction[d];
    projection += offset[d] * direction[d];
  }
  const double t =
      (length_squared > 0.0) ? std::clamp(projection / length_squared, 0.0, 1.0) : 0.0;
  double distance_squared = 0.0;
  for (int d = 0; d < point.Size(); d++)
  {
    const double delta = offset[d] - t * direction[d];
    distance_squared += delta * delta;
  }
  return distance_squared;
}

double EdgeDistanceTree::DistanceSquaredToBox(const mfem::Vector &point, const Node &node)
{
  double distance_squared = 0.0;
  for (int d = 0; d < point.Size(); d++)
  {
    const double delta = (point[d] < node.min[d])
                             ? node.min[d] - point[d]
                             : ((point[d] > node.max[d]) ? point[d] - node.max[d] : 0.0);
    distance_squared += delta * delta;
  }
  return distance_squared;
}

int EdgeDistanceTree::Build(std::size_t begin, std::size_t end)
{
  Node node;
  node.begin = begin;
  node.end = end;
  node.min.fill(mfem::infinity());
  node.max.fill(-mfem::infinity());
  for (std::size_t i = begin; i < end; i++)
  {
    const auto &segment = segments[indices[i]];
    for (int d = 0; d < 3; d++)
    {
      node.min[d] = std::min({node.min[d], segment.p0[d], segment.p1[d]});
      node.max[d] = std::max({node.max[d], segment.p0[d], segment.p1[d]});
    }
  }

  const int node_index = static_cast<int>(nodes.size());
  nodes.push_back(node);
  if (end - begin <= leaf_size)
  {
    return node_index;
  }

  int axis = 0;
  for (int d = 1; d < 3; d++)
  {
    if (node.max[d] - node.min[d] > node.max[axis] - node.min[axis])
    {
      axis = d;
    }
  }
  const std::size_t mid = begin + (end - begin) / 2;
  std::nth_element(indices.begin() + begin, indices.begin() + mid, indices.begin() + end,
                   [this, axis](std::size_t a, std::size_t b)
                   {
                     const auto &sa = segments[a];
                     const auto &sb = segments[b];
                     return sa.p0[axis] + sa.p1[axis] < sb.p0[axis] + sb.p1[axis];
                   });
  const int left = Build(begin, mid);
  const int right = Build(mid, end);
  nodes[node_index].left = left;
  nodes[node_index].right = right;
  return node_index;
}

EdgeDistanceTree::NearestResult EdgeDistanceTree::Search(int node_index,
                                                         const mfem::Vector &point,
                                                         NearestResult best) const
{
  const Node &node = nodes[node_index];
  if (DistanceSquaredToBox(point, node) > best.distance_squared)
  {
    return best;
  }
  if (node.IsLeaf())
  {
    for (std::size_t i = node.begin; i < node.end; i++)
    {
      const std::size_t segment = indices[i];
      const double distance_squared = DistanceSquaredToSegment(point, segments[segment]);
      if (distance_squared < best.distance_squared ||
          (distance_squared == best.distance_squared && segment < best.segment))
      {
        best = {distance_squared, segment};
      }
    }
    return best;
  }

  const double left_distance = DistanceSquaredToBox(point, nodes[node.left]);
  const double right_distance = DistanceSquaredToBox(point, nodes[node.right]);
  if (left_distance < right_distance)
  {
    best = Search(node.left, point, best);
    return Search(node.right, point, best);
  }
  best = Search(node.right, point, best);
  return Search(node.left, point, best);
}

EdgeDistanceTree::EdgeDistanceTree(std::vector<mesh::BoundaryEdgeSegment> segments_,
                                   std::vector<std::array<double, 3>> process_normals_,
                                   std::vector<SegmentMetadata> metadata_)
  : segments(std::move(segments_)), process_normals(std::move(process_normals_)),
    metadata(std::move(metadata_)), indices(segments.size())
{
  MFEM_VERIFY(!segments.empty(), "Cannot build an empty edge-distance tree!");
  MFEM_VERIFY(process_normals.empty() || process_normals.size() == segments.size(),
              "Edge-distance tree process normals must match its segments!");
  MFEM_VERIFY(metadata.empty() || metadata.size() == segments.size(),
              "Edge-distance tree metadata must match its segments!");
  std::iota(indices.begin(), indices.end(), 0);
  nodes.reserve(2 * segments.size());
  Build(0, segments.size());
}

double EdgeDistanceTree::DistanceSquared(const mfem::Vector &point) const
{
  return Nearest(point).distance_squared;
}

EdgeDistanceTree::NearestResult EdgeDistanceTree::Nearest(const mfem::Vector &point) const
{
  return Search(0, point, {mfem::infinity(), segments.size()});
}

double
EdgeDistanceTree::DistanceAlongEdgeToNonregularVertex(const mfem::Vector &point,
                                                      std::size_t segment_index) const
{
  if (!HasMetadata())
  {
    return mfem::infinity();
  }
  MFEM_ASSERT(segment_index < segments.size(), "Invalid edge segment index!");
  const auto &segment = segments[segment_index];
  double length_squared = 0.0;
  double projection = 0.0;
  for (int d = 0; d < point.Size(); d++)
  {
    const double direction = segment.p1[d] - segment.p0[d];
    length_squared += direction * direction;
    projection += (point[d] - segment.p0[d]) * direction;
  }
  MFEM_VERIFY(length_squared > 0.0,
              "Automatic 3D edge topology contains a zero-length segment!");
  const double length = std::sqrt(length_squared);
  const double coordinate = std::clamp(projection / length, 0.0, length);
  const auto &distances = metadata[segment_index].vertex_distances;
  return std::min(distances[0] + coordinate, distances[1] + length - coordinate);
}

EdgeFrame BuildEdgeFrame(const mesh::BoundaryEdgeSegment &segment,
                         const std::array<double, 3> &configured_normal, int dimension)
{
  EdgeFrame frame{};
  frame.normal = configured_normal;
  if (dimension == 2)
  {
    frame.normal[2] = 0.0;
    frame.tangent[2] = 1.0;
  }
  else
  {
    double tangent_norm_squared = 0.0;
    for (int d = 0; d < 3; d++)
    {
      frame.tangent[d] = segment.p1[d] - segment.p0[d];
      tangent_norm_squared += frame.tangent[d] * frame.tangent[d];
    }
    MFEM_VERIFY(tangent_norm_squared > 0.0,
                "A 3D polarized edge diagnostic encountered a zero-length segment!");
    const double inverse_tangent_norm = 1.0 / std::sqrt(tangent_norm_squared);
    for (double &value : frame.tangent)
    {
      value *= inverse_tangent_norm;
    }
    double normal_tangent = 0.0;
    for (int d = 0; d < 3; d++)
    {
      normal_tangent += frame.normal[d] * frame.tangent[d];
    }
    for (int d = 0; d < 3; d++)
    {
      frame.normal[d] -= normal_tangent * frame.tangent[d];
    }
  }

  double normal_norm_squared = 0.0;
  for (int d = 0; d < 3; d++)
  {
    normal_norm_squared += frame.normal[d] * frame.normal[d];
  }
  MFEM_VERIFY(normal_norm_squared > 1.0e-12,
              "Interface dielectric \"EdgeFrameNormal\" is parallel to an edge segment!");
  const double inverse_normal_norm = 1.0 / std::sqrt(normal_norm_squared);
  for (double &value : frame.normal)
  {
    value *= inverse_normal_norm;
  }

  frame.transverse = {
      frame.normal[1] * frame.tangent[2] - frame.normal[2] * frame.tangent[1],
      frame.normal[2] * frame.tangent[0] - frame.normal[0] * frame.tangent[2],
      frame.normal[0] * frame.tangent[1] - frame.normal[1] * frame.tangent[0]};
  return frame;
}

std::array<double, 6> GetPolarizedEdgeEnergyDensity(
    const mfem::Vector &point, const mesh::BoundaryEdgeSegment &segment,
    const EdgeFrame &frame, const mfem::Vector &field, const mfem::Vector &displacement)
{
  MFEM_ASSERT(point.Size() == field.Size() && field.Size() == displacement.Size(),
              "Edge energy vectors must have matching dimensions!");
  double side_coordinate = 0.0;
  for (int d = 0; d < point.Size(); d++)
  {
    side_coordinate += (point[d] - segment.p0[d]) * frame.normal[d];
  }
  const std::size_t side_offset = side_coordinate >= 0.0 ? 0 : 3;

  std::array<double, 6> energy{};
  for (std::size_t component = 0; component < 3; component++)
  {
    const auto &axis = component == 0   ? frame.normal
                       : component == 1 ? frame.transverse
                                        : frame.tangent;
    double field_component = 0.0;
    double displacement_component = 0.0;
    for (int d = 0; d < field.Size(); d++)
    {
      field_component += field[d] * axis[d];
      displacement_component += displacement[d] * axis[d];
    }
    energy[side_offset + component] = 0.5 * field_component * displacement_component;
  }
  return energy;
}

namespace
{

bool IsCoincidentWithExcludedBoundary(const mesh::BoundaryEdgeSegment &segment,
                                      const EdgeDistanceTree &excluded_tree,
                                      int space_dimension,
                                      double distance_tolerance_squared)
{
  mfem::Vector point(space_dimension);
  for (const double t : {0.0, 0.5, 1.0})
  {
    for (int d = 0; d < space_dimension; d++)
    {
      point[d] = (1.0 - t) * segment.p0[d] + t * segment.p1[d];
    }
    if (excluded_tree.DistanceSquared(point) > distance_tolerance_squared)
    {
      return false;
    }
  }
  return true;
}

struct ElementBounds
{
  mfem::Vector center;
  double radius;
  double diameter;
};

ElementBounds GetElementBounds(const mfem::ParMesh &mesh, int element)
{
  const int dim = mesh.SpaceDimension();
  mfem::DenseMatrix points;
  if (mesh.GetNodes())
  {
    auto &mutable_mesh = const_cast<mfem::ParMesh &>(mesh);
    auto &T = *mutable_mesh.GetElementTransformation(element);
    const int order = mesh.GetNodes()->FESpace()->GetMaxElementOrder();
    const auto *refined =
        mfem::GlobGeometryRefiner.Refine(T.GetGeometryType(), std::max(1, order));
    T.Transform(refined->RefPts, points);
  }
  else
  {
    const int *vertices = mesh.GetElement(element)->GetVertices();
    const int count = mesh.GetElement(element)->GetNVertices();
    points.SetSize(dim, count);
    for (int j = 0; j < count; j++)
    {
      const double *coordinate = mesh.GetVertex(vertices[j]);
      for (int d = 0; d < dim; d++)
      {
        points(d, j) = coordinate[d];
      }
    }
  }

  mfem::Vector center(dim);
  center = 0.0;
  for (int j = 0; j < points.Width(); j++)
  {
    for (int d = 0; d < dim; d++)
    {
      center[d] += points(d, j);
    }
  }
  center /= points.Width();

  double radius_squared = 0.0;
  for (int j = 0; j < points.Width(); j++)
  {
    double distance_squared = 0.0;
    for (int d = 0; d < dim; d++)
    {
      const double delta = points(d, j) - center[d];
      distance_squared += delta * delta;
    }
    radius_squared = std::max(radius_squared, distance_squared);
  }
  const double radius = std::sqrt(radius_squared);
  return {std::move(center), radius, 2.0 * radius};
}

}  // namespace

std::shared_ptr<const EdgeDistanceTree>
BuildEdgeDistanceTree(const mfem::ParMesh &mesh, const std::vector<int> &edge_attributes,
                      const std::vector<int> &edge_exclude_attributes,
                      const std::optional<std::array<double, 3>> &process_normal)
{
  auto edge_marker = mesh::BdrAttrToMarker(mesh, edge_attributes, true);
  auto edge_segments = mesh::GetBoundaryEdgeSegments(mesh, edge_marker);
  MFEM_VERIFY(!edge_segments.empty(),
              "No perimeter was found for interface dielectric edge attributes!");
  if (!edge_exclude_attributes.empty())
  {
    auto exclude_marker = mesh::BdrAttrToMarker(mesh, edge_exclude_attributes, true);
    auto excluded_segments = mesh::GetBoundaryElementEdgeSegments(mesh, exclude_marker);
    MFEM_VERIFY(!excluded_segments.empty(),
                "No boundary geometry was found for interface dielectric edge exclusion "
                "attributes!");
    const EdgeDistanceTree excluded_tree(std::move(excluded_segments));
    mfem::Vector bbmin, bbmax;
    mesh::GetAxisAlignedBoundingBox(mesh, bbmin, bbmax);
    double extent = 0.0;
    for (int d = 0; d < mesh.SpaceDimension(); d++)
    {
      extent = std::max(extent, bbmax[d] - bbmin[d]);
    }
    MFEM_VERIFY(extent > 0.0,
                "Degenerate mesh geometry for interface dielectric edge exclusion!");
    const double tolerance_squared = 1.0e-20 * extent * extent;
    edge_segments.erase(std::remove_if(edge_segments.begin(), edge_segments.end(),
                                       [&](const auto &segment)
                                       {
                                         return IsCoincidentWithExcludedBoundary(
                                             segment, excluded_tree, mesh.SpaceDimension(),
                                             tolerance_squared);
                                       }),
                        edge_segments.end());
    MFEM_VERIFY(!edge_segments.empty(),
                "Interface dielectric edge exclusion removed the entire perimeter!");
  }
  if (process_normal && mesh.SpaceDimension() == 3)
  {
    const double normal_squared = std::inner_product(
        process_normal->begin(), process_normal->end(), process_normal->begin(), 0.0);
    MFEM_VERIFY(normal_squared > 0.0,
                "Interface dielectric \"EdgeFrameNormal\" must be nonzero!");
    edge_segments.erase(std::remove_if(edge_segments.begin(), edge_segments.end(),
                                       [&](const auto &segment)
                                       {
                                         std::array<double, 3> tangent{};
                                         double tangent_squared = 0.0;
                                         double normal_tangent = 0.0;
                                         for (int d = 0; d < 3; d++)
                                         {
                                           tangent[d] = segment.p1[d] - segment.p0[d];
                                           tangent_squared += tangent[d] * tangent[d];
                                           normal_tangent +=
                                               (*process_normal)[d] * tangent[d];
                                         }
                                         return normal_squared * tangent_squared -
                                                    normal_tangent * normal_tangent <=
                                                1.0e-12 * normal_squared * tangent_squared;
                                       }),
                        edge_segments.end());
    MFEM_VERIFY(!edge_segments.empty(),
                "Interface dielectric edge extraction retained no segments transverse "
                "to \"EdgeFrameNormal\"!");
  }
  return std::make_shared<EdgeDistanceTree>(std::move(edge_segments));
}

std::vector<std::size_t>
GetInterfaceMetalEdgeSegmentIndices(const MetalEdgeGeometry &geometry, int interface_index,
                                    InterfaceDielectric interface_type)
{
  MFEM_VERIFY(interface_type != InterfaceDielectric::DEFAULT,
              "Automatic metal edge selection requires an explicit MA, MS, or SA "
              "interface type!");
  std::vector<std::size_t> indices;
  for (std::size_t i = 0; i < geometry.segments.size(); i++)
  {
    const auto &segment = geometry.segments[i];
    if (segment.type != MetalEdgeSegmentType::PHYSICAL)
    {
      continue;
    }
    const auto &interfaces =
        interface_type == InterfaceDielectric::SA   ? segment.sa_interfaces
        : interface_type == InterfaceDielectric::MS ? segment.ms_interfaces
                                                    : segment.ma_interfaces;
    if (std::find(interfaces.begin(), interfaces.end(), interface_index) !=
        interfaces.end())
    {
      indices.push_back(i);
    }
  }
  MFEM_VERIFY(!indices.empty(),
              "No physical metal perimeter was found for automatically extracted "
              "interface "
                  << interface_index << "!");
  return indices;
}

void ExcludeMetalEdgeSegmentIndices(const mfem::ParMesh &mesh,
                                    const MetalEdgeGeometry &geometry,
                                    const std::vector<int> &exclude_attributes,
                                    std::vector<std::size_t> &segment_indices)
{
  if (exclude_attributes.empty())
  {
    return;
  }
  auto exclude_marker = mesh::BdrAttrToMarker(mesh, exclude_attributes, true);
  auto excluded_segments = mesh::GetBoundaryElementEdgeSegments(mesh, exclude_marker);
  MFEM_VERIFY(!excluded_segments.empty(),
              "No boundary geometry was found for interface dielectric edge exclusion "
              "attributes!");
  const EdgeDistanceTree excluded_tree(std::move(excluded_segments));
  mfem::Vector bbmin, bbmax;
  mesh::GetAxisAlignedBoundingBox(mesh, bbmin, bbmax);
  double extent = 0.0;
  for (int d = 0; d < mesh.SpaceDimension(); d++)
  {
    extent = std::max(extent, bbmax[d] - bbmin[d]);
  }
  MFEM_VERIFY(extent > 0.0,
              "Degenerate mesh geometry for interface dielectric edge exclusion!");
  const double tolerance_squared = 1.0e-20 * extent * extent;
  segment_indices.erase(
      std::remove_if(segment_indices.begin(), segment_indices.end(),
                     [&](std::size_t index)
                     {
                       MFEM_VERIFY(index < geometry.segments.size(),
                                   "Invalid automatically extracted metal edge segment "
                                   "index!");
                       const auto &source = geometry.segments[index];
                       const mesh::BoundaryEdgeSegment segment{
                           geometry.vertices[source.vertices[0]].coordinate,
                           geometry.vertices[source.vertices[1]].coordinate};
                       return IsCoincidentWithExcludedBoundary(segment, excluded_tree,
                                                               mesh.SpaceDimension(),
                                                               tolerance_squared);
                     }),
      segment_indices.end());
  MFEM_VERIFY(!segment_indices.empty(),
              "Interface dielectric edge exclusion removed the entire automatically "
              "extracted perimeter!");
}

std::shared_ptr<const EdgeDistanceTree>
BuildEdgeDistanceTree(const MetalEdgeGeometry &geometry,
                      const std::vector<std::size_t> &segment_indices,
                      std::vector<std::array<double, 3>> process_normals)
{
  using VertexDistance = std::pair<double, std::size_t>;
  std::vector<double> vertex_distances(geometry.vertices.size(), mfem::infinity());
  std::priority_queue<VertexDistance, std::vector<VertexDistance>,
                      std::greater<VertexDistance>>
      queue;
  for (std::size_t vertex = 0; vertex < geometry.vertices.size(); vertex++)
  {
    const auto type = geometry.vertices[vertex].physical_type;
    if (type && *type != MetalEdgeVertexType::REGULAR)
    {
      vertex_distances[vertex] = 0.0;
      queue.emplace(0.0, vertex);
    }
  }
  while (!queue.empty())
  {
    const auto [distance, vertex] = queue.top();
    queue.pop();
    if (distance != vertex_distances[vertex])
    {
      continue;
    }
    for (const std::size_t segment_index : geometry.vertices[vertex].segments)
    {
      const auto &segment = geometry.segments[segment_index];
      if (segment.type != MetalEdgeSegmentType::PHYSICAL)
      {
        continue;
      }
      const std::size_t other =
          segment.vertices[0] == vertex ? segment.vertices[1] : segment.vertices[0];
      double length_squared = 0.0;
      for (int d = 0; d < 3; d++)
      {
        const double delta = geometry.vertices[other].coordinate[d] -
                             geometry.vertices[vertex].coordinate[d];
        length_squared += delta * delta;
      }
      MFEM_VERIFY(length_squared > 0.0,
                  "Automatic 3D edge topology contains a zero-length segment!");
      const double candidate = distance + std::sqrt(length_squared);
      if (candidate < vertex_distances[other])
      {
        vertex_distances[other] = candidate;
        queue.emplace(candidate, other);
      }
    }
  }

  std::vector<mesh::BoundaryEdgeSegment> segments;
  std::vector<EdgeDistanceTree::SegmentMetadata> metadata;
  segments.reserve(segment_indices.size());
  metadata.reserve(segment_indices.size());
  for (const std::size_t index : segment_indices)
  {
    MFEM_VERIFY(index < geometry.segments.size(),
                "Invalid automatically extracted metal edge segment index!");
    const auto &segment = geometry.segments[index];
    MFEM_VERIFY(segment.type == MetalEdgeSegmentType::PHYSICAL,
                "Automatic edge-distance trees cannot contain truncation segments!");
    segments.push_back({geometry.vertices[segment.vertices[0]].coordinate,
                        geometry.vertices[segment.vertices[1]].coordinate});
    EdgeDistanceTree::SegmentMetadata data;
    data.component = segment.physical_component;
    data.chain = segment.physical_chain;
    for (int endpoint = 0; endpoint < 2; endpoint++)
    {
      const auto type = geometry.vertices[segment.vertices[endpoint]].physical_type;
      MFEM_VERIFY(type, "Physical metal edge segment has an unclassified endpoint!");
      data.vertex_types[endpoint] = static_cast<int>(*type);
      data.vertex_distances[endpoint] = vertex_distances[segment.vertices[endpoint]];
    }
    metadata.push_back(data);
  }
  MFEM_VERIFY(!segments.empty(),
              "Cannot build an empty automatic metal edge-distance tree!");
  return std::make_shared<EdgeDistanceTree>(std::move(segments), std::move(process_normals),
                                            std::move(metadata));
}

std::vector<EdgeRefinementContext>
BuildEdgeRefinementContexts(const mfem::ParMesh &mesh,
                            const config::BoundaryData &boundaries)
{
  using TreeKey =
      std::tuple<bool, std::vector<int>, std::vector<int>, std::vector<std::size_t>>;
  using ContextKey = std::tuple<TreeKey, double, int, double, double>;
  std::map<TreeKey, std::shared_ptr<const EdgeDistanceTree>> trees;
  std::map<ContextKey, EdgeRefinementContext> contexts;
  MetalEdgeGeometry metal_edges;
  bool have_metal_edges = false;
  for (const auto &[index, dielectric] : boundaries.postpro.dielectric)
  {
    if (!dielectric.edge_refinement)
    {
      continue;
    }
    const auto &refinement = *dielectric.edge_refinement;
    std::vector<std::size_t> segment_indices;
    if (dielectric.automatic_edges)
    {
      if (!have_metal_edges)
      {
        metal_edges = ExtractMetalEdgeGeometry(mesh, boundaries);
        have_metal_edges = true;
      }
      segment_indices =
          GetInterfaceMetalEdgeSegmentIndices(metal_edges, index, dielectric.type);
      ExcludeMetalEdgeSegmentIndices(mesh, metal_edges, dielectric.edge_exclude_attributes,
                                     segment_indices);
    }
    const TreeKey tree_key{dielectric.automatic_edges, dielectric.edge_attributes,
                           dielectric.edge_exclude_attributes, segment_indices};
    auto tree = trees.find(tree_key);
    if (tree == trees.end())
    {
      auto distance_tree = dielectric.automatic_edges
                               ? BuildEdgeDistanceTree(metal_edges, segment_indices)
                               : BuildEdgeDistanceTree(mesh, dielectric.edge_attributes,
                                                       dielectric.edge_exclude_attributes);
      tree = trees.try_emplace(tree_key, std::move(distance_tree)).first;
    }
    const ContextKey context_key{
        tree_key, refinement.radius, refinement.elements_per_radius,
        refinement.outer_radius_factor, refinement.core_indicator_weight};
    contexts.try_emplace(
        context_key,
        EdgeRefinementContext{tree->second, refinement.radius,
                              refinement.radius / refinement.elements_per_radius,
                              refinement.outer_radius_factor * refinement.radius,
                              refinement.core_indicator_weight});
  }

  std::vector<EdgeRefinementContext> result;
  result.reserve(contexts.size());
  for (auto &[key, context] : contexts)
  {
    (void)key;
    result.push_back(std::move(context));
  }
  return result;
}

mfem::Array<int>
MarkEdgeRefinementElements(const mfem::ParMesh &mesh,
                           const std::vector<EdgeRefinementContext> &contexts)
{
  mfem::Array<int> marked;
  marked.Reserve(mesh.GetNE());
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    const auto bounds = GetElementBounds(mesh, element);
    for (const auto &context : contexts)
    {
      const double distance =
          std::sqrt(context.distance_tree->DistanceSquared(bounds.center));
      if (distance <= context.outer_radius + bounds.radius &&
          bounds.diameter > context.target_size * (1.0 + 1.0e-12))
      {
        marked.Append(element);
        break;
      }
    }
  }
  return marked;
}

int WeightEdgeCoreIndicators(const mfem::ParMesh &mesh,
                             const std::vector<EdgeRefinementContext> &contexts,
                             mfem::Vector &indicators)
{
  MFEM_VERIFY(indicators.Size() == mesh.GetNE(),
              "Edge-core indicator weighting requires one value per mesh element!");
  int weighted = 0;
  auto *values = indicators.HostReadWrite();
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    const auto bounds = GetElementBounds(mesh, element);
    double weight = 1.0;
    for (const auto &context : contexts)
    {
      const double distance =
          std::sqrt(context.distance_tree->DistanceSquared(bounds.center));
      if (distance + bounds.radius < context.radius)
      {
        weight = std::min(weight, context.core_indicator_weight);
      }
    }
    if (weight < 1.0)
    {
      values[element] *= weight;
      weighted++;
    }
  }
  return weighted;
}

}  // namespace palace
