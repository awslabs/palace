// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "singularfeatures.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <stdexcept>
#include <tuple>
#include <nlohmann/json.hpp>

#include "fem/singulargeometry.hpp"
#include "utils/communication.hpp"

namespace palace
{

namespace fem
{

namespace singular
{

namespace
{

using EdgeKey = std::array<int, 2>;
using StableEdgeKey = std::array<GlobalVertexId, 2>;
using StableElementKey = std::array<GlobalVertexId, 4>;
using json = nlohmann::json;

struct SelectedEdge
{
  int mesh_edge = -1;
  std::vector<int> boundary_elements;
  std::vector<int> boundary_attributes;
};

EdgeKey GetEdgeKey(const mfem::Mesh &mesh, int edge)
{
  mfem::Array<int> vertices;
  mesh.GetEdgeVertices(edge, vertices);
  if (vertices.Size() != 2 || vertices[0] == vertices[1])
  {
    throw std::runtime_error("Singular-feature extraction found an invalid mesh edge!");
  }
  EdgeKey key{vertices[0], vertices[1]};
  std::sort(key.begin(), key.end());
  return key;
}

std::array<double, 3> DirectionFromVertex(const mfem::Mesh &mesh, int vertex,
                                          const std::array<int, 2> &segment_vertices)
{
  if (mesh.SpaceDimension() < 1 || mesh.SpaceDimension() > 3 ||
      (segment_vertices[0] != vertex && segment_vertices[1] != vertex))
  {
    throw std::runtime_error(
        "Singular-feature segment does not contain its incident vertex!");
  }
  const int other =
      segment_vertices[0] == vertex ? segment_vertices[1] : segment_vertices[0];
  std::array<double, 3> x{}, y{};
  mesh.GetNode(vertex, x.data());
  mesh.GetNode(other, y.data());

  std::array<double, 3> direction{};
  double norm_squared = 0.0;
  for (int d = 0; d < mesh.SpaceDimension(); d++)
  {
    direction[d] = y[d] - x[d];
    norm_squared += direction[d] * direction[d];
  }
  if (!(norm_squared > 0.0))
  {
    throw std::runtime_error(
        "Singular-feature extraction found a zero-length perimeter edge!");
  }
  const double inverse_norm = 1.0 / std::sqrt(norm_squared);
  for (double &value : direction)
  {
    value *= inverse_norm;
  }
  return direction;
}

template <typename GetVertexId>
std::array<int, 4> CanonicalNodeTuple(const mfem::Element &element, int singular_vertex,
                                      GetVertexId &&get_vertex_id)
{
  const int *vertices = element.GetVertices();
  int singular_local = -1;
  std::array<std::pair<GlobalVertexId, int>, 3> remaining;
  int next = 0;
  for (int local = 0; local < 4; local++)
  {
    if (vertices[local] == singular_vertex)
    {
      singular_local = local;
    }
    else
    {
      remaining[next++] = {get_vertex_id(vertices[local]), local};
    }
  }
  if (singular_local < 0 || next != 3)
  {
    throw std::runtime_error(
        "Singular-feature node is not a vertex of its incident tetrahedron!");
  }
  std::sort(remaining.begin(), remaining.end());
  return {singular_local, remaining[0].second, remaining[1].second, remaining[2].second};
}

std::array<int, 4> CanonicalNodeTuple(const mfem::Element &element, int singular_vertex)
{
  return CanonicalNodeTuple(element, singular_vertex,
                            [](int vertex) { return static_cast<GlobalVertexId>(vertex); });
}

template <typename GetVertexId>
std::array<int, 3> CanonicalTriangleNodeTuple(const mfem::Element &element,
                                              int singular_vertex,
                                              GetVertexId &&get_vertex_id)
{
  const int *vertices = element.GetVertices();
  int singular_local = -1;
  std::array<std::pair<GlobalVertexId, int>, 2> remaining;
  int next = 0;
  for (int local = 0; local < 3; local++)
  {
    if (vertices[local] == singular_vertex)
    {
      singular_local = local;
    }
    else
    {
      remaining[next++] = {get_vertex_id(vertices[local]), local};
    }
  }
  if (singular_local < 0 || next != 2)
  {
    throw std::runtime_error(
        "Singular-feature tip is not a vertex of its incident triangle!");
  }
  std::sort(remaining.begin(), remaining.end());
  return {singular_local, remaining[0].second, remaining[1].second};
}

std::array<int, 3> CanonicalTriangleNodeTuple(const mfem::Element &element,
                                              int singular_vertex)
{
  return CanonicalTriangleNodeTuple(element, singular_vertex, [](int vertex)
                                    { return static_cast<GlobalVertexId>(vertex); });
}

template <typename GetVertexId>
std::array<int, 4> CanonicalEdgeTuple(const mfem::Element &element,
                                      const std::array<int, 2> &singular_edge,
                                      GetVertexId &&get_vertex_id)
{
  const int *vertices = element.GetVertices();
  std::array<std::pair<GlobalVertexId, int>, 2> singular_local{
      std::pair<GlobalVertexId, int>{-1, -1}, std::pair<GlobalVertexId, int>{-1, -1}};
  std::array<std::pair<GlobalVertexId, int>, 2> remaining;
  int remaining_count = 0;
  for (int local = 0; local < 4; local++)
  {
    if (vertices[local] == singular_edge[0])
    {
      singular_local[0] = {get_vertex_id(vertices[local]), local};
    }
    else if (vertices[local] == singular_edge[1])
    {
      singular_local[1] = {get_vertex_id(vertices[local]), local};
    }
    else
    {
      remaining[remaining_count++] = {get_vertex_id(vertices[local]), local};
    }
  }
  if (singular_local[0].second < 0 || singular_local[1].second < 0 || remaining_count != 2)
  {
    throw std::runtime_error(
        "Singular-feature edge is not an edge of its incident tetrahedron!");
  }
  std::sort(singular_local.begin(), singular_local.end());
  std::sort(remaining.begin(), remaining.end());
  return {singular_local[0].second, singular_local[1].second, remaining[0].second,
          remaining[1].second};
}

std::array<int, 4> CanonicalEdgeTuple(const mfem::Element &element,
                                      const EdgeKey &singular_edge)
{
  return CanonicalEdgeTuple(element, singular_edge,
                            [](int vertex) { return static_cast<GlobalVertexId>(vertex); });
}

template <std::size_t N, typename GetVertexId>
std::array<GlobalVertexId, N> StableVertexKey(const int *vertices,
                                              GetVertexId &&get_vertex_id)
{
  std::array<GlobalVertexId, N> result;
  for (std::size_t i = 0; i < N; i++)
  {
    result[i] = get_vertex_id(vertices[i]);
  }
  std::sort(result.begin(), result.end());
  if (std::adjacent_find(result.begin(), result.end()) != result.end())
  {
    throw std::invalid_argument(
        "A singular-feature mesh entity contains duplicate stable vertex IDs!");
  }
  return result;
}

StableEdgeKey StableMeshEdgeKey(const mfem::Mesh &mesh, int edge,
                                const std::vector<GlobalVertexId> &vertex_ids)
{
  mfem::Array<int> vertices;
  mesh.GetEdgeVertices(edge, vertices);
  if (vertices.Size() != 2)
  {
    throw std::invalid_argument("A rank-local singular-feature edge has invalid topology!");
  }
  return StableVertexKey<2>(vertices.GetData(),
                            [&vertex_ids](int vertex) { return vertex_ids[vertex]; });
}

void ValidateMesh(const mfem::Mesh &mesh, bool allow_parallel = false)
{
  const auto *parallel_mesh = dynamic_cast<const mfem::ParMesh *>(&mesh);
  if (!allow_parallel && parallel_mesh && parallel_mesh->GetNRanks() > 1)
  {
    throw std::invalid_argument(
        "Serial singular-feature extraction cannot use a multi-rank ParMesh!");
  }
  if (mesh.Dimension() != 3 || mesh.SpaceDimension() != 3)
  {
    throw std::invalid_argument(
        "Singular-feature extraction requires a three-dimensional mesh!");
  }
  if (mesh.Nonconforming())
  {
    throw std::invalid_argument(
        "Singular-feature extraction does not support nonconforming meshes!");
  }
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    if (mesh.GetElementGeometry(element) != mfem::Geometry::TETRAHEDRON)
    {
      throw std::invalid_argument(
          "Singular-feature extraction initially supports only tetrahedral meshes!");
    }
  }
}

void ValidateSerialFeatureBlueprint(const mfem::Mesh &mesh, const FeatureTopology &features)
{
  if (features.elements.size() != static_cast<std::size_t>(mesh.GetNE()))
  {
    throw std::invalid_argument(
        "Singular-feature blueprint does not match its source serial mesh!");
  }
  for (std::size_t i = 0; i < features.vertices.size(); i++)
  {
    const auto &vertex = features.vertices[i];
    if (vertex.id != i || vertex.mesh_vertex < 0 || vertex.mesh_vertex >= mesh.GetNV() ||
        !std::isfinite(vertex.nu) || vertex.nu <= 0.0 || vertex.nu >= 1.0)
    {
      throw std::invalid_argument(
          "Singular-feature blueprint contains an invalid node record!");
    }
  }
  for (std::size_t i = 0; i < features.features.size(); i++)
  {
    const auto &feature = features.features[i];
    if (feature.id != i || !std::isfinite(feature.nu) || feature.nu <= 0.0 ||
        feature.nu >= 1.0)
    {
      throw std::invalid_argument(
          "Singular-feature blueprint contains an invalid straight feature!");
    }
  }
  for (std::size_t i = 0; i < features.segments.size(); i++)
  {
    const auto &segment = features.segments[i];
    if (segment.mesh_edge < 0 || segment.mesh_edge >= mesh.GetNEdges() ||
        segment.feature >= features.features.size() ||
        GetEdgeKey(mesh, segment.mesh_edge) != segment.mesh_vertices)
    {
      throw std::invalid_argument(
          "Singular-feature blueprint contains an invalid segment record!");
    }
    const auto &feature_segments = features.features[segment.feature].segments;
    if (std::find(feature_segments.begin(), feature_segments.end(), i) ==
        feature_segments.end())
    {
      throw std::invalid_argument(
          "Singular-feature segment is absent from its straight feature!");
    }
  }
  mfem::Array<int> face_vertices;
  std::set<std::array<GlobalVertexId, 3>> sheet_face_keys;
  for (const auto &face : features.sheet_faces)
  {
    if (face.mesh_face < 0 || face.mesh_face >= mesh.GetNumFaces() ||
        face.boundary_attribute <= 0)
    {
      throw std::invalid_argument(
          "Singular-feature blueprint contains an invalid sheet face!");
    }
    mesh.GetFaceVertices(face.mesh_face, face_vertices);
    if (face_vertices.Size() != 3)
    {
      throw std::invalid_argument(
          "Singular-feature blueprint sheet face is not triangular!");
    }
    const auto key = StableVertexKey<3>(face_vertices.GetData(), [](int vertex)
                                        { return static_cast<GlobalVertexId>(vertex); });
    if (key != face.mesh_vertices || !sheet_face_keys.insert(key).second)
    {
      throw std::invalid_argument(
          "Singular-feature blueprint contains an inconsistent sheet face!");
    }
  }

  for (int element = 0; element < mesh.GetNE(); element++)
  {
    const auto &tetrahedron = *mesh.GetElement(element);
    for (const auto &node : features.elements[element].nodes)
    {
      if (node.vertex >= features.vertices.size() ||
          node.mesh_vertex != features.vertices[node.vertex].mesh_vertex ||
          node.canonical_nodes != CanonicalNodeTuple(tetrahedron, node.mesh_vertex))
      {
        throw std::invalid_argument(
            "Singular-feature blueprint has inconsistent node incidence!");
      }
    }
    for (const auto &edge : features.elements[element].edges)
    {
      if (edge.segment >= features.segments.size())
      {
        throw std::invalid_argument(
            "Singular-feature blueprint has invalid edge incidence!");
      }
      const auto &segment = features.segments[edge.segment];
      if (edge.feature != segment.feature || edge.mesh_edge != segment.mesh_edge ||
          edge.canonical_nodes != CanonicalEdgeTuple(tetrahedron, segment.mesh_vertices))
      {
        throw std::invalid_argument(
            "Singular-feature blueprint has inconsistent edge incidence!");
      }
    }
  }
}

bool IsPermutation(const std::array<int, 4> &indices)
{
  auto sorted = indices;
  std::sort(sorted.begin(), sorted.end());
  return sorted == std::array<int, 4>{0, 1, 2, 3};
}

void ValidateFeatureBlueprintStructure(const FeatureTopology &features)
{
  std::set<int> feature_vertices;
  for (std::size_t i = 0; i < features.vertices.size(); i++)
  {
    const auto &vertex = features.vertices[i];
    if (vertex.id != i || vertex.mesh_vertex < 0 ||
        !feature_vertices.insert(vertex.mesh_vertex).second || !std::isfinite(vertex.nu) ||
        vertex.nu <= 0.0 || vertex.nu >= 1.0)
    {
      throw std::invalid_argument(
          "Singular-feature blueprint contains an invalid node record!");
    }
    for (std::size_t segment : vertex.segments)
    {
      if (segment >= features.segments.size())
      {
        throw std::invalid_argument("Singular-feature node references an invalid segment!");
      }
    }
    for (std::size_t feature : vertex.features)
    {
      if (feature >= features.features.size())
      {
        throw std::invalid_argument(
            "Singular-feature node references an invalid straight feature!");
      }
    }
  }
  for (std::size_t i = 0; i < features.features.size(); i++)
  {
    const auto &feature = features.features[i];
    if (feature.id != i || !std::isfinite(feature.nu) || feature.nu <= 0.0 ||
        feature.nu >= 1.0)
    {
      throw std::invalid_argument(
          "Singular-feature blueprint contains an invalid straight feature!");
    }
    for (std::size_t segment : feature.segments)
    {
      if (segment >= features.segments.size())
      {
        throw std::invalid_argument(
            "Straight singular feature references an invalid segment!");
      }
    }
  }
  for (std::size_t i = 0; i < features.segments.size(); i++)
  {
    const auto &segment = features.segments[i];
    if (segment.mesh_edge < 0 || segment.mesh_vertices[0] < 0 ||
        segment.mesh_vertices[0] >= segment.mesh_vertices[1] ||
        segment.feature >= features.features.size() ||
        std::find(features.features[segment.feature].segments.begin(),
                  features.features[segment.feature].segments.end(),
                  i) == features.features[segment.feature].segments.end())
    {
      throw std::invalid_argument(
          "Singular-feature blueprint contains an invalid segment record!");
    }
  }
  std::set<std::array<GlobalVertexId, 3>> sheet_faces;
  for (const auto &face : features.sheet_faces)
  {
    if (face.mesh_face < 0 || face.boundary_attribute <= 0 || face.mesh_vertices[0] < 0 ||
        !std::is_sorted(face.mesh_vertices.begin(), face.mesh_vertices.end()) ||
        std::adjacent_find(face.mesh_vertices.begin(), face.mesh_vertices.end()) !=
            face.mesh_vertices.end() ||
        !sheet_faces.insert(face.mesh_vertices).second)
    {
      throw std::invalid_argument(
          "Singular-feature blueprint contains an invalid sheet face!");
    }
  }
  for (const auto &element : features.elements)
  {
    for (const auto &node : element.nodes)
    {
      if (node.vertex >= features.vertices.size() ||
          node.mesh_vertex != features.vertices[node.vertex].mesh_vertex ||
          !IsPermutation(node.canonical_nodes))
      {
        throw std::invalid_argument(
            "Singular-feature blueprint has inconsistent node incidence!");
      }
    }
    for (const auto &edge : element.edges)
    {
      if (edge.segment >= features.segments.size() ||
          edge.feature != features.segments[edge.segment].feature ||
          edge.mesh_edge != features.segments[edge.segment].mesh_edge ||
          !IsPermutation(edge.canonical_nodes))
      {
        throw std::invalid_argument(
            "Singular-feature blueprint has inconsistent edge incidence!");
      }
    }
  }
}

json PackFeatureTopology(const FeatureTopology &features)
{
  json packed{{"num_elements", features.elements.size()},
              {"vertices", json::array()},
              {"segments", json::array()},
              {"features", json::array()},
              {"sheet_faces", json::array()},
              {"elements", json::array()}};
  for (const auto &vertex : features.vertices)
  {
    packed["vertices"].push_back({vertex.id, vertex.mesh_vertex,
                                  static_cast<int>(vertex.type), vertex.segments,
                                  vertex.features, vertex.nu});
  }
  for (const auto &segment : features.segments)
  {
    packed["segments"].push_back({segment.mesh_edge, segment.mesh_vertices, segment.feature,
                                  segment.boundary_attributes});
  }
  for (const auto &feature : features.features)
  {
    packed["features"].push_back(
        {feature.id, feature.segments, feature.mesh_vertices, feature.nu, feature.closed});
  }
  for (const auto &face : features.sheet_faces)
  {
    packed["sheet_faces"].push_back(
        {face.mesh_face, face.mesh_vertices, face.boundary_attribute});
  }
  for (std::size_t element = 0; element < features.elements.size(); element++)
  {
    const auto &incidence = features.elements[element];
    if (incidence.nodes.empty() && incidence.edges.empty())
    {
      continue;
    }
    json nodes = json::array(), edges = json::array();
    for (const auto &node : incidence.nodes)
    {
      nodes.push_back({node.vertex, node.mesh_vertex, node.canonical_nodes});
    }
    for (const auto &edge : incidence.edges)
    {
      edges.push_back({edge.feature, edge.segment, edge.mesh_edge, edge.canonical_nodes});
    }
    packed["elements"].push_back({element, std::move(nodes), std::move(edges)});
  }
  return packed;
}

FeatureTopology UnpackFeatureTopology(const json &packed)
{
  FeatureTopology features;
  features.elements.resize(packed.at("num_elements").get<std::size_t>());
  for (const auto &entry : packed.at("vertices"))
  {
    features.vertices.push_back({entry.at(0).get<std::size_t>(), entry.at(1).get<int>(),
                                 static_cast<FeatureVertexType>(entry.at(2).get<int>()),
                                 entry.at(3).get<std::vector<std::size_t>>(),
                                 entry.at(4).get<std::vector<std::size_t>>(),
                                 entry.at(5).get<double>()});
  }
  for (const auto &entry : packed.at("segments"))
  {
    features.segments.push_back(
        {entry.at(0).get<int>(), entry.at(1).get<std::array<int, 2>>(),
         entry.at(2).get<std::size_t>(), entry.at(3).get<std::vector<int>>()});
  }
  for (const auto &entry : packed.at("features"))
  {
    features.features.push_back({entry.at(0).get<std::size_t>(),
                                 entry.at(1).get<std::vector<std::size_t>>(),
                                 entry.at(2).get<std::vector<int>>(),
                                 entry.at(3).get<double>(), entry.at(4).get<bool>()});
  }
  for (const auto &entry : packed.at("sheet_faces"))
  {
    features.sheet_faces.push_back({entry.at(0).get<int>(),
                                    entry.at(1).get<std::array<GlobalVertexId, 3>>(),
                                    entry.at(2).get<int>()});
  }
  std::set<std::size_t> unpacked_elements;
  for (const auto &entry : packed.at("elements"))
  {
    const std::size_t element = entry.at(0).get<std::size_t>();
    if (element >= features.elements.size() || !unpacked_elements.insert(element).second)
    {
      throw std::invalid_argument(
          "Serialized singular-feature blueprint contains an invalid element!");
    }
    auto &incidence = features.elements[element];
    for (const auto &node : entry.at(1))
    {
      incidence.nodes.push_back({node.at(0).get<std::size_t>(), node.at(1).get<int>(),
                                 node.at(2).get<std::array<int, 4>>()});
    }
    for (const auto &edge : entry.at(2))
    {
      incidence.edges.push_back({edge.at(0).get<std::size_t>(),
                                 edge.at(1).get<std::size_t>(), edge.at(2).get<int>(),
                                 edge.at(3).get<std::array<int, 4>>()});
    }
  }
  ValidateFeatureBlueprintStructure(features);
  return features;
}

bool IsTrianglePermutation(const std::array<int, 3> &indices)
{
  auto sorted = indices;
  std::sort(sorted.begin(), sorted.end());
  return sorted == std::array<int, 3>{0, 1, 2};
}

void ValidateTriangleFeatureBlueprintStructure(const TriangleFeatureTopology &features)
{
  std::set<int> tip_vertices;
  for (std::size_t i = 0; i < features.vertices.size(); i++)
  {
    const auto &vertex = features.vertices[i];
    if (vertex.id != i || vertex.mesh_vertex < 0 ||
        !tip_vertices.insert(vertex.mesh_vertex).second ||
        vertex.selected_segments.empty() || !std::isfinite(vertex.nu) || vertex.nu <= 0.0 ||
        vertex.nu >= 1.0)
    {
      throw std::invalid_argument(
          "Singular line-tip blueprint contains an invalid tip record!");
    }
    for (std::size_t segment : vertex.selected_segments)
    {
      if (segment >= features.selected_segments.size())
      {
        throw std::invalid_argument(
            "Singular line tip references an invalid selected segment!");
      }
    }
  }
  std::set<StableEdgeKey> selected_edges;
  for (const auto &segment : features.selected_segments)
  {
    const StableEdgeKey edge{segment.mesh_vertices[0], segment.mesh_vertices[1]};
    if (segment.boundary_element < 0 || segment.mesh_edge < 0 ||
        segment.boundary_attribute <= 0 || edge[0] < 0 || edge[0] >= edge[1] ||
        !selected_edges.insert(edge).second)
    {
      throw std::invalid_argument(
          "Singular line-tip blueprint contains an invalid selected segment!");
    }
  }
  for (const auto &element : features.elements)
  {
    std::set<std::size_t> element_tips;
    for (const auto &node : element.nodes)
    {
      if (node.vertex >= features.vertices.size() ||
          node.mesh_vertex != features.vertices[node.vertex].mesh_vertex ||
          !IsTrianglePermutation(node.canonical_nodes) ||
          !element_tips.insert(node.vertex).second)
      {
        throw std::invalid_argument(
            "Singular line-tip blueprint has inconsistent element incidence!");
      }
    }
  }
}

json PackTriangleFeatureTopology(const TriangleFeatureTopology &features)
{
  json packed{{"num_elements", features.elements.size()},
              {"vertices", json::array()},
              {"selected_segments", json::array()},
              {"elements", json::array()}};
  for (const auto &vertex : features.vertices)
  {
    packed["vertices"].push_back(
        {vertex.id, vertex.mesh_vertex, vertex.selected_segments, vertex.nu});
  }
  for (const auto &segment : features.selected_segments)
  {
    packed["selected_segments"].push_back({segment.boundary_element, segment.mesh_edge,
                                           segment.mesh_vertices,
                                           segment.boundary_attribute});
  }
  for (std::size_t element = 0; element < features.elements.size(); element++)
  {
    const auto &incidence = features.elements[element];
    if (incidence.nodes.empty())
    {
      continue;
    }
    json nodes = json::array();
    for (const auto &node : incidence.nodes)
    {
      nodes.push_back({node.vertex, node.mesh_vertex, node.canonical_nodes});
    }
    packed["elements"].push_back({element, std::move(nodes)});
  }
  return packed;
}

TriangleFeatureTopology UnpackTriangleFeatureTopology(const json &packed)
{
  TriangleFeatureTopology features;
  features.elements.resize(packed.at("num_elements").get<std::size_t>());
  for (const auto &entry : packed.at("vertices"))
  {
    features.vertices.push_back({entry.at(0).get<std::size_t>(), entry.at(1).get<int>(),
                                 entry.at(2).get<std::vector<std::size_t>>(),
                                 entry.at(3).get<double>()});
  }
  for (const auto &entry : packed.at("selected_segments"))
  {
    features.selected_segments.push_back({entry.at(0).get<int>(), entry.at(1).get<int>(),
                                          entry.at(2).get<std::array<int, 2>>(),
                                          entry.at(3).get<int>()});
  }
  std::set<std::size_t> unpacked_elements;
  for (const auto &entry : packed.at("elements"))
  {
    const std::size_t element = entry.at(0).get<std::size_t>();
    if (element >= features.elements.size() || !unpacked_elements.insert(element).second)
    {
      throw std::invalid_argument(
          "Serialized singular line-tip blueprint contains an invalid element!");
    }
    for (const auto &node : entry.at(1))
    {
      features.elements[element].nodes.push_back({node.at(0).get<std::size_t>(),
                                                  node.at(1).get<int>(),
                                                  node.at(2).get<std::array<int, 3>>()});
    }
  }
  ValidateTriangleFeatureBlueprintStructure(features);
  return features;
}

}  // namespace

TriangleFeatureTopology
ExtractSerialLineTipFeatures(const mfem::Mesh &mesh,
                             const std::vector<int> &boundary_attributes, double nu)
{
  const auto *parallel_mesh = dynamic_cast<const mfem::ParMesh *>(&mesh);
  if (parallel_mesh && parallel_mesh->GetNRanks() > 1)
  {
    throw std::invalid_argument(
        "Serial singular-line extraction cannot use a multi-rank ParMesh!");
  }
  if (mesh.Dimension() != 2 || mesh.SpaceDimension() != 2)
  {
    throw std::invalid_argument(
        "Singular-line extraction requires a two-dimensional mesh!");
  }
  if (mesh.Nonconforming())
  {
    throw std::invalid_argument(
        "Singular-line extraction does not support nonconforming meshes!");
  }
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    if (mesh.GetElementGeometry(element) != mfem::Geometry::TRIANGLE)
    {
      throw std::invalid_argument(
          "Singular-line extraction initially supports only triangular meshes!");
    }
  }
  if (!std::isfinite(nu) || nu <= 0.0 || nu >= 1.0)
  {
    throw std::invalid_argument(
        "Singular-line exponent must be finite and satisfy 0 < nu < 1!");
  }

  std::set<int> selected_attributes;
  for (int attribute : boundary_attributes)
  {
    if (attribute <= 0)
    {
      throw std::invalid_argument("Singular-line boundary attributes must be positive!");
    }
    selected_attributes.insert(attribute);
  }

  TriangleFeatureTopology result;
  result.elements.resize(mesh.GetNE());
  if (selected_attributes.empty())
  {
    return result;
  }

  std::set<int> present_attributes;
  std::map<EdgeKey, TriangleSelectedSegment> selected_edges;
  for (int boundary_element = 0; boundary_element < mesh.GetNBE(); boundary_element++)
  {
    const int attribute = mesh.GetBdrAttribute(boundary_element);
    if (selected_attributes.find(attribute) == selected_attributes.end())
    {
      continue;
    }
    present_attributes.insert(attribute);
    if (mesh.GetBdrElementGeometry(boundary_element) != mfem::Geometry::SEGMENT)
    {
      throw std::invalid_argument(
          "Selected singular-line boundary elements must be segments!");
    }

    int edge, orientation, element_1, element_2;
    mesh.GetBdrElementFace(boundary_element, &edge, &orientation);
    mesh.GetFaceElements(edge, &element_1, &element_2);
    if (element_1 < 0 || element_2 < 0)
    {
      throw std::invalid_argument(
          "Selected zero-thickness PEC lines must be internal mesh boundaries!");
    }
    mfem::IsoparametricTransformation edge_transformation;
    mesh.GetEdgeTransformation(edge, &edge_transformation);
    if (!IsGeometricallyStraightSegmentTransformation(edge_transformation))
    {
      throw std::invalid_argument(
          "Selected zero-thickness PEC line segments must be geometrically straight!");
    }
    const auto key = GetEdgeKey(mesh, edge);
    if (!selected_edges
             .emplace(key, TriangleSelectedSegment{boundary_element, edge, key, attribute})
             .second)
    {
      throw std::invalid_argument(
          "A selected zero-thickness PEC segment is represented more than once!");
    }
  }
  if (present_attributes != selected_attributes)
  {
    throw std::invalid_argument(
        "At least one requested singular-line boundary attribute is absent!");
  }

  std::set<int> exterior_vertices;
  mfem::Array<int> face_vertices;
  for (int face = 0; face < mesh.GetNumFaces(); face++)
  {
    int element_1, element_2;
    mesh.GetFaceElements(face, &element_1, &element_2);
    if (element_1 >= 0 && element_2 >= 0)
    {
      continue;
    }
    mesh.GetFaceVertices(face, face_vertices);
    if (face_vertices.Size() != 2)
    {
      throw std::runtime_error("Triangular mesh edge has invalid topology!");
    }
    exterior_vertices.insert(face_vertices[0]);
    exterior_vertices.insert(face_vertices[1]);
  }

  std::map<int, std::vector<std::size_t>> vertex_segments;
  for (const auto &[key, segment] : selected_edges)
  {
    (void)key;
    const std::size_t index = result.selected_segments.size();
    result.selected_segments.push_back(segment);
    for (int vertex : segment.mesh_vertices)
    {
      vertex_segments[vertex].push_back(index);
    }
  }

  constexpr double straight_tolerance = 1.0e-10;
  for (const auto &[vertex, segments] : vertex_segments)
  {
    if (segments.size() > 2)
    {
      throw std::invalid_argument(
          "Selected zero-thickness PEC lines contain a branch or junction!");
    }
    if (segments.size() == 2)
    {
      std::array<std::array<double, 2>, 2> direction;
      for (int side = 0; side < 2; side++)
      {
        const auto &segment = result.selected_segments[segments[side]];
        const auto physical_direction =
            DirectionFromVertex(mesh, vertex, segment.mesh_vertices);
        direction[side] = {physical_direction[0], physical_direction[1]};
      }
      const double dot =
          direction[0][0] * direction[1][0] + direction[0][1] * direction[1][1];
      if (dot > -1.0 + straight_tolerance)
      {
        throw std::invalid_argument(
            "Selected zero-thickness PEC lines must be straight; corner exponents "
            "are not implemented yet!");
      }
      continue;
    }

    if (exterior_vertices.find(vertex) == exterior_vertices.end())
    {
      result.vertices.push_back({result.vertices.size(), vertex, segments, nu});
    }
  }

  std::map<int, std::size_t> tip_by_vertex;
  for (const auto &vertex : result.vertices)
  {
    tip_by_vertex.emplace(vertex.mesh_vertex, vertex.id);
  }
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    const auto &triangle = *mesh.GetElement(element);
    const int *vertices = triangle.GetVertices();
    for (int local = 0; local < 3; local++)
    {
      const auto tip = tip_by_vertex.find(vertices[local]);
      if (tip == tip_by_vertex.end())
      {
        continue;
      }
      result.elements[element].nodes.push_back(
          {tip->second, vertices[local],
           CanonicalTriangleNodeTuple(triangle, vertices[local])});
    }
  }
  ValidateTriangleFeatureBlueprintStructure(result);
  return result;
}

FeatureTopology ExtractSerialSheetFeatures(const mfem::Mesh &mesh,
                                           const std::vector<int> &boundary_attributes,
                                           double nu)
{
  ValidateMesh(mesh);
  if (!std::isfinite(nu) || nu <= 0.0 || nu >= 1.0)
  {
    throw std::invalid_argument(
        "Singular-feature exponent must be finite and satisfy 0 < nu < 1!");
  }

  std::set<int> selected_attributes;
  for (int attribute : boundary_attributes)
  {
    if (attribute <= 0)
    {
      throw std::invalid_argument("Singular-feature boundary attributes must be positive!");
    }
    selected_attributes.insert(attribute);
  }

  FeatureTopology result;
  result.elements.resize(mesh.GetNE());
  if (selected_attributes.empty())
  {
    return result;
  }

  std::set<int> present_attributes;
  std::map<int, int> selected_faces;
  std::map<EdgeKey, SelectedEdge> selected_edges;
  mfem::Array<int> edges, orientations;
  for (int boundary_element = 0; boundary_element < mesh.GetNBE(); boundary_element++)
  {
    const int attribute = mesh.GetBdrAttribute(boundary_element);
    if (selected_attributes.find(attribute) == selected_attributes.end())
    {
      continue;
    }
    present_attributes.insert(attribute);
    if (mesh.GetBdrElementGeometry(boundary_element) != mfem::Geometry::TRIANGLE)
    {
      throw std::invalid_argument(
          "Selected singular-feature boundary elements must be triangles!");
    }

    int face, orientation, element_1, element_2;
    mesh.GetBdrElementFace(boundary_element, &face, &orientation);
    mesh.GetFaceElements(face, &element_1, &element_2);
    if (element_1 < 0 || element_2 < 0)
    {
      throw std::invalid_argument(
          "Selected zero-thickness PEC sheets must be internal mesh boundaries!");
    }
    if (!selected_faces.emplace(face, boundary_element).second)
    {
      throw std::invalid_argument(
          "A selected zero-thickness PEC face is represented more than once!");
    }

    mesh.GetBdrElementEdges(boundary_element, edges, orientations);
    if (edges.Size() != 3)
    {
      throw std::runtime_error(
          "Selected triangular PEC face has an invalid edge topology!");
    }
    for (int edge : edges)
    {
      const auto key = GetEdgeKey(mesh, edge);
      auto &[mesh_edge, boundary_elements, attributes] = selected_edges[key];
      if (mesh_edge >= 0 && mesh_edge != edge)
      {
        throw std::runtime_error(
            "One topological PEC edge maps to multiple mesh edge IDs!");
      }
      mesh_edge = edge;
      boundary_elements.push_back(boundary_element);
      attributes.push_back(attribute);
    }
  }

  if (present_attributes != selected_attributes)
  {
    throw std::invalid_argument(
        "At least one requested singular-feature boundary attribute is absent!");
  }

  mfem::Array<int> face_vertices;
  for (const auto &[face, boundary_element] : selected_faces)
  {
    mesh.GetFaceVertices(face, face_vertices);
    if (face_vertices.Size() != 3)
    {
      throw std::runtime_error("Selected PEC sheet face has invalid topology!");
    }
    result.sheet_faces.push_back(
        {face,
         StableVertexKey<3>(face_vertices.GetData(),
                            [](int vertex) { return static_cast<GlobalVertexId>(vertex); }),
         mesh.GetBdrAttribute(boundary_element)});
  }

  for (const auto &[key, edge] : selected_edges)
  {
    (void)key;
    if (edge.boundary_elements.size() > 2)
    {
      throw std::invalid_argument("Selected PEC sheets are nonmanifold along a mesh edge!");
    }
  }

  for (auto &[key, edge] : selected_edges)
  {
    if (edge.boundary_elements.size() != 1)
    {
      continue;
    }
    std::sort(edge.boundary_attributes.begin(), edge.boundary_attributes.end());
    edge.boundary_attributes.erase(
        std::unique(edge.boundary_attributes.begin(), edge.boundary_attributes.end()),
        edge.boundary_attributes.end());
    mfem::IsoparametricTransformation edge_transformation;
    mesh.GetEdgeTransformation(edge.mesh_edge, &edge_transformation);
    if (!IsGeometricallyStraightSegmentTransformation(edge_transformation))
    {
      throw std::invalid_argument(
          "Selected PEC sheet perimeter edges must be geometrically straight!");
    }
    result.segments.push_back(
        {edge.mesh_edge, key, 0, std::move(edge.boundary_attributes)});
  }

  std::map<int, std::size_t> vertex_index;
  for (std::size_t segment = 0; segment < result.segments.size(); segment++)
  {
    for (int vertex : result.segments[segment].mesh_vertices)
    {
      auto [it, inserted] = vertex_index.emplace(vertex, result.vertices.size());
      if (inserted)
      {
        result.vertices.push_back(
            {result.vertices.size(), vertex, FeatureVertexType::REGULAR, {}, {}, nu});
      }
      result.vertices[it->second].segments.push_back(segment);
    }
  }

  for (auto &vertex : result.vertices)
  {
    std::sort(vertex.segments.begin(), vertex.segments.end());
    if (vertex.segments.size() == 1)
    {
      vertex.type = FeatureVertexType::ENDPOINT;
    }
    else if (vertex.segments.size() == 2)
    {
      const auto direction_0 = DirectionFromVertex(
          mesh, vertex.mesh_vertex, result.segments[vertex.segments[0]].mesh_vertices);
      const auto direction_1 = DirectionFromVertex(
          mesh, vertex.mesh_vertex, result.segments[vertex.segments[1]].mesh_vertices);
      double dot = 0.0;
      for (int d = 0; d < 3; d++)
      {
        dot += direction_0[d] * direction_1[d];
      }
      constexpr double straight_tolerance = 1.0e-10;
      vertex.type = dot <= -1.0 + straight_tolerance ? FeatureVertexType::REGULAR
                                                     : FeatureVertexType::CORNER;
    }
    else
    {
      throw std::invalid_argument(
          "Selected PEC sheet perimeter contains an unsupported junction!");
    }
  }

  std::vector<int> segment_feature(result.segments.size(), -1);
  for (std::size_t seed = 0; seed < result.segments.size(); seed++)
  {
    if (segment_feature[seed] >= 0)
    {
      continue;
    }
    const std::size_t feature = result.features.size();
    StraightFeature record{feature, {}, {}, nu, true};
    std::queue<std::size_t> queue;
    queue.push(seed);
    segment_feature[seed] = static_cast<int>(feature);
    while (!queue.empty())
    {
      const std::size_t segment = queue.front();
      queue.pop();
      record.segments.push_back(segment);
      result.segments[segment].feature = feature;
      for (int vertex : result.segments[segment].mesh_vertices)
      {
        record.mesh_vertices.push_back(vertex);
        const auto &vertex_record = result.vertices[vertex_index.at(vertex)];
        if (vertex_record.type != FeatureVertexType::REGULAR)
        {
          record.closed = false;
          continue;
        }
        for (std::size_t neighbor : vertex_record.segments)
        {
          if (segment_feature[neighbor] < 0)
          {
            segment_feature[neighbor] = static_cast<int>(feature);
            queue.push(neighbor);
          }
        }
      }
    }
    std::sort(record.segments.begin(), record.segments.end());
    std::sort(record.mesh_vertices.begin(), record.mesh_vertices.end());
    record.mesh_vertices.erase(
        std::unique(record.mesh_vertices.begin(), record.mesh_vertices.end()),
        record.mesh_vertices.end());
    result.features.push_back(std::move(record));
  }

  for (auto &vertex : result.vertices)
  {
    for (std::size_t segment : vertex.segments)
    {
      vertex.features.push_back(result.segments[segment].feature);
    }
    std::sort(vertex.features.begin(), vertex.features.end());
    vertex.features.erase(std::unique(vertex.features.begin(), vertex.features.end()),
                          vertex.features.end());
  }

  std::map<int, std::size_t> edge_segment;
  std::map<int, std::size_t> singular_vertices;
  for (std::size_t segment = 0; segment < result.segments.size(); segment++)
  {
    edge_segment.emplace(result.segments[segment].mesh_edge, segment);
  }
  for (std::size_t vertex = 0; vertex < result.vertices.size(); vertex++)
  {
    singular_vertices.emplace(result.vertices[vertex].mesh_vertex, vertex);
  }

  mfem::Array<int> element_edges, edge_orientations;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    const auto &tetrahedron = *mesh.GetElement(element);
    auto &incidence = result.elements[element];
    mesh.GetElementEdges(element, element_edges, edge_orientations);
    for (int edge : element_edges)
    {
      const auto segment_it = edge_segment.find(edge);
      if (segment_it == edge_segment.end())
      {
        continue;
      }
      const std::size_t segment = segment_it->second;
      const auto &segment_record = result.segments[segment];
      incidence.edges.push_back(
          {segment_record.feature, segment, edge,
           CanonicalEdgeTuple(tetrahedron, segment_record.mesh_vertices)});
    }

    const int *element_vertices = tetrahedron.GetVertices();
    for (int local = 0; local < 4; local++)
    {
      const int mesh_vertex = element_vertices[local];
      const auto vertex_it = singular_vertices.find(mesh_vertex);
      if (vertex_it == singular_vertices.end())
      {
        continue;
      }
      incidence.nodes.push_back(
          {vertex_it->second, mesh_vertex, CanonicalNodeTuple(tetrahedron, mesh_vertex)});
    }
    std::sort(incidence.edges.begin(), incidence.edges.end(),
              [](const auto &a, const auto &b)
              { return std::tie(a.feature, a.segment) < std::tie(b.feature, b.segment); });
    std::sort(
        incidence.nodes.begin(), incidence.nodes.end(), [](const auto &a, const auto &b)
        { return std::tie(a.vertex, a.mesh_vertex) < std::tie(b.vertex, b.mesh_vertex); });
  }

  return result;
}

FeatureTopology
DistributeSerialSheetFeatures(const mfem::Mesh &serial_mesh,
                              const FeatureTopology &serial_features,
                              const mfem::ParMesh &parallel_mesh,
                              const std::vector<GlobalVertexId> &serial_vertex_ids)
{
  ValidateMesh(serial_mesh);
  ValidateMesh(parallel_mesh, true);
  ValidateSerialFeatureBlueprint(serial_mesh, serial_features);
  if (serial_mesh.Dimension() != parallel_mesh.Dimension() ||
      serial_mesh.SpaceDimension() != parallel_mesh.SpaceDimension() ||
      serial_vertex_ids.size() != static_cast<std::size_t>(parallel_mesh.GetNV()))
  {
    throw std::invalid_argument(
        "Serial feature blueprint and rank-local mesh are incompatible!");
  }

  std::map<GlobalVertexId, int> local_vertices;
  for (int local_vertex = 0; local_vertex < parallel_mesh.GetNV(); local_vertex++)
  {
    const GlobalVertexId serial_vertex = serial_vertex_ids[local_vertex];
    if (serial_vertex < 0 || serial_vertex >= serial_mesh.GetNV() ||
        !local_vertices.emplace(serial_vertex, local_vertex).second)
    {
      throw std::invalid_argument(
          "Rank-local mesh has invalid or duplicate source serial vertex IDs!");
    }
    std::array<double, 3> serial_coordinate{}, parallel_coordinate{};
    serial_mesh.GetNode(static_cast<int>(serial_vertex), serial_coordinate.data());
    parallel_mesh.GetNode(local_vertex, parallel_coordinate.data());
    for (int d = 0; d < serial_mesh.SpaceDimension(); d++)
    {
      const double scale =
          1.0 + std::max(std::abs(serial_coordinate[d]), std::abs(parallel_coordinate[d]));
      if (std::abs(serial_coordinate[d] - parallel_coordinate[d]) >
          32.0 * std::numeric_limits<double>::epsilon() * scale)
      {
        throw std::invalid_argument(
            "Rank-local vertex and source serial vertex have different coordinates!");
      }
    }
  }

  std::map<StableElementKey, int> serial_elements;
  for (int element = 0; element < serial_mesh.GetNE(); element++)
  {
    const auto &tetrahedron = *serial_mesh.GetElement(element);
    const auto key = StableVertexKey<4>(tetrahedron.GetVertices(), [](int vertex)
                                        { return static_cast<GlobalVertexId>(vertex); });
    if (!serial_elements.emplace(key, element).second)
    {
      throw std::invalid_argument(
          "Source serial mesh contains duplicate tetrahedron vertex sets!");
    }
  }

  std::map<StableEdgeKey, int> local_edges;
  for (int edge = 0; edge < parallel_mesh.GetNEdges(); edge++)
  {
    if (!local_edges
             .emplace(StableMeshEdgeKey(parallel_mesh, edge, serial_vertex_ids), edge)
             .second)
    {
      throw std::invalid_argument("Rank-local mesh contains duplicate edge vertex sets!");
    }
  }

  FeatureTopology result = serial_features;
  result.elements.assign(parallel_mesh.GetNE(), {});
  std::set<int> matched_serial_elements;
  const auto local_vertex_id = [&serial_vertex_ids](int vertex)
  { return serial_vertex_ids[vertex]; };
  for (int local_element = 0; local_element < parallel_mesh.GetNE(); local_element++)
  {
    const auto &tetrahedron = *parallel_mesh.GetElement(local_element);
    const auto key = StableVertexKey<4>(tetrahedron.GetVertices(), local_vertex_id);
    const auto serial_element = serial_elements.find(key);
    if (serial_element == serial_elements.end() ||
        !matched_serial_elements.insert(serial_element->second).second)
    {
      throw std::invalid_argument(
          "Rank-local tetrahedron does not map uniquely to the source serial mesh!");
    }

    const auto &source = serial_features.elements[serial_element->second];
    auto &incidence = result.elements[local_element];
    for (const auto &node : source.nodes)
    {
      const GlobalVertexId source_vertex = node.mesh_vertex;
      const auto local_vertex = local_vertices.find(source_vertex);
      if (local_vertex == local_vertices.end())
      {
        throw std::invalid_argument(
            "Rank-local tetrahedron is missing an incident singular node!");
      }
      incidence.nodes.push_back(
          {node.vertex, local_vertex->second,
           CanonicalNodeTuple(tetrahedron, local_vertex->second, local_vertex_id)});
    }
    for (const auto &edge : source.edges)
    {
      if (edge.segment >= serial_features.segments.size())
      {
        throw std::invalid_argument("Source singular-feature edge incidence is invalid!");
      }
      const auto &segment = serial_features.segments[edge.segment];
      const StableEdgeKey source_edge{
          static_cast<GlobalVertexId>(segment.mesh_vertices[0]),
          static_cast<GlobalVertexId>(segment.mesh_vertices[1])};
      const auto local_edge = local_edges.find(source_edge);
      const auto first = local_vertices.find(source_edge[0]);
      const auto second = local_vertices.find(source_edge[1]);
      if (local_edge == local_edges.end() || first == local_vertices.end() ||
          second == local_vertices.end())
      {
        throw std::invalid_argument(
            "Rank-local tetrahedron is missing an incident singular edge!");
      }
      const std::array<int, 2> local_edge_vertices{first->second, second->second};
      incidence.edges.push_back(
          {edge.feature, edge.segment, local_edge->second,
           CanonicalEdgeTuple(tetrahedron, local_edge_vertices, local_vertex_id)});
    }
  }
  return result;
}

void BroadcastSerialSheetFeatures(FeatureTopology &serial_features, MPI_Comm comm)
{
  std::vector<std::uint8_t> bytes;
  if (Mpi::Root(comm))
  {
    ValidateFeatureBlueprintStructure(serial_features);
    bytes = json::to_cbor(PackFeatureTopology(serial_features));
  }
  std::int64_t size = static_cast<std::int64_t>(bytes.size());
  if (Mpi::Root(comm) && static_cast<std::size_t>(size) != bytes.size())
  {
    throw std::overflow_error("Serialized singular-feature blueprint is too large!");
  }
  Mpi::Broadcast(1, &size, 0, comm);
  if (size < 0)
  {
    throw std::runtime_error("Serialized singular-feature blueprint has invalid size!");
  }
  if (!Mpi::Root(comm))
  {
    bytes.resize(static_cast<std::size_t>(size));
  }
  Mpi::BroadcastLarge(size, bytes.data(), 0, comm);
  if (!Mpi::Root(comm))
  {
    serial_features = UnpackFeatureTopology(json::from_cbor(bytes));
  }
}

void BroadcastSerialLineTipFeatures(TriangleFeatureTopology &serial_features, MPI_Comm comm)
{
  std::vector<std::uint8_t> bytes;
  if (Mpi::Root(comm))
  {
    ValidateTriangleFeatureBlueprintStructure(serial_features);
    bytes = json::to_cbor(PackTriangleFeatureTopology(serial_features));
  }
  std::int64_t size = static_cast<std::int64_t>(bytes.size());
  if (Mpi::Root(comm) && static_cast<std::size_t>(size) != bytes.size())
  {
    throw std::overflow_error("Serialized singular line-tip blueprint is too large!");
  }
  Mpi::Broadcast(1, &size, 0, comm);
  if (size < 0)
  {
    throw std::runtime_error("Serialized singular line-tip blueprint has invalid size!");
  }
  if (!Mpi::Root(comm))
  {
    bytes.resize(static_cast<std::size_t>(size));
  }
  Mpi::BroadcastLarge(size, bytes.data(), 0, comm);
  if (!Mpi::Root(comm))
  {
    serial_features = UnpackTriangleFeatureTopology(json::from_cbor(bytes));
  }
}

FeatureTopology
DistributeSerialSheetFeatures(const FeatureTopology &serial_features,
                              const mfem::ParMesh &parallel_mesh,
                              const std::vector<GlobalVertexId> &serial_vertex_ids,
                              const std::vector<GlobalVertexId> &source_element_ids)
{
  ValidateMesh(parallel_mesh, true);
  ValidateFeatureBlueprintStructure(serial_features);
  if (serial_vertex_ids.size() != static_cast<std::size_t>(parallel_mesh.GetNV()) ||
      source_element_ids.size() != static_cast<std::size_t>(parallel_mesh.GetNE()))
  {
    throw std::invalid_argument(
        "Serial feature blueprint and rank-local source maps are incompatible!");
  }

  std::map<GlobalVertexId, int> local_vertices;
  for (int local_vertex = 0; local_vertex < parallel_mesh.GetNV(); local_vertex++)
  {
    const GlobalVertexId serial_vertex = serial_vertex_ids[local_vertex];
    if (serial_vertex < 0 || !local_vertices.emplace(serial_vertex, local_vertex).second)
    {
      throw std::invalid_argument(
          "Rank-local mesh has invalid or duplicate source serial vertex IDs!");
    }
  }
  std::map<StableEdgeKey, int> local_edges;
  for (int edge = 0; edge < parallel_mesh.GetNEdges(); edge++)
  {
    if (!local_edges
             .emplace(StableMeshEdgeKey(parallel_mesh, edge, serial_vertex_ids), edge)
             .second)
    {
      throw std::invalid_argument("Rank-local mesh contains duplicate edge vertex sets!");
    }
  }

  FeatureTopology result = serial_features;
  result.elements.assign(parallel_mesh.GetNE(), {});
  std::set<GlobalVertexId> matched_serial_elements;
  const auto local_vertex_id = [&serial_vertex_ids](int vertex)
  { return serial_vertex_ids[vertex]; };
  for (int local_element = 0; local_element < parallel_mesh.GetNE(); local_element++)
  {
    const GlobalVertexId serial_element = source_element_ids[local_element];
    if (serial_element < 0 ||
        serial_element >= static_cast<GlobalVertexId>(serial_features.elements.size()) ||
        !matched_serial_elements.insert(serial_element).second)
    {
      throw std::invalid_argument(
          "Rank-local tetrahedron has an invalid source serial element ID!");
    }

    const auto &tetrahedron = *parallel_mesh.GetElement(local_element);
    const auto &source = serial_features.elements[serial_element];
    auto &incidence = result.elements[local_element];
    for (const auto &node : source.nodes)
    {
      const GlobalVertexId source_vertex = node.mesh_vertex;
      const auto local_vertex = local_vertices.find(source_vertex);
      if (local_vertex == local_vertices.end())
      {
        throw std::invalid_argument(
            "Rank-local tetrahedron is missing an incident singular node!");
      }
      incidence.nodes.push_back(
          {node.vertex, local_vertex->second,
           CanonicalNodeTuple(tetrahedron, local_vertex->second, local_vertex_id)});
    }
    for (const auto &edge : source.edges)
    {
      const auto &segment = serial_features.segments[edge.segment];
      const StableEdgeKey source_edge{
          static_cast<GlobalVertexId>(segment.mesh_vertices[0]),
          static_cast<GlobalVertexId>(segment.mesh_vertices[1])};
      const auto local_edge = local_edges.find(source_edge);
      const auto first = local_vertices.find(source_edge[0]);
      const auto second = local_vertices.find(source_edge[1]);
      if (local_edge == local_edges.end() || first == local_vertices.end() ||
          second == local_vertices.end())
      {
        throw std::invalid_argument(
            "Rank-local tetrahedron is missing an incident singular edge!");
      }
      const std::array<int, 2> local_edge_vertices{first->second, second->second};
      incidence.edges.push_back(
          {edge.feature, edge.segment, local_edge->second,
           CanonicalEdgeTuple(tetrahedron, local_edge_vertices, local_vertex_id)});
    }
  }
  return result;
}

TriangleFeatureTopology
DistributeSerialLineTipFeatures(const TriangleFeatureTopology &serial_features,
                                const mfem::ParMesh &parallel_mesh,
                                const std::vector<GlobalVertexId> &serial_vertex_ids,
                                const std::vector<GlobalVertexId> &source_element_ids)
{
  ValidateTriangleFeatureBlueprintStructure(serial_features);
  if (parallel_mesh.Dimension() != 2 || parallel_mesh.SpaceDimension() != 2 ||
      parallel_mesh.Nonconforming() ||
      serial_vertex_ids.size() != static_cast<std::size_t>(parallel_mesh.GetNV()) ||
      source_element_ids.size() != static_cast<std::size_t>(parallel_mesh.GetNE()))
  {
    throw std::invalid_argument(
        "Singular line-tip blueprint and rank-local source maps are incompatible!");
  }
  std::map<GlobalVertexId, int> local_vertices;
  for (int local_vertex = 0; local_vertex < parallel_mesh.GetNV(); local_vertex++)
  {
    const GlobalVertexId serial_vertex = serial_vertex_ids[local_vertex];
    if (serial_vertex < 0 || !local_vertices.emplace(serial_vertex, local_vertex).second)
    {
      throw std::invalid_argument(
          "Rank-local triangular mesh has invalid or duplicate source vertex IDs!");
    }
  }

  TriangleFeatureTopology result = serial_features;
  result.elements.assign(parallel_mesh.GetNE(), {});
  std::set<GlobalVertexId> matched_serial_elements;
  const auto local_vertex_id = [&serial_vertex_ids](int vertex)
  { return serial_vertex_ids[vertex]; };
  for (int local_element = 0; local_element < parallel_mesh.GetNE(); local_element++)
  {
    if (parallel_mesh.GetElementGeometry(local_element) != mfem::Geometry::TRIANGLE)
    {
      throw std::invalid_argument(
          "Singular line-tip distribution requires triangular elements!");
    }
    const GlobalVertexId serial_element = source_element_ids[local_element];
    if (serial_element < 0 ||
        serial_element >= static_cast<GlobalVertexId>(serial_features.elements.size()) ||
        !matched_serial_elements.insert(serial_element).second)
    {
      throw std::invalid_argument(
          "Rank-local triangle has an invalid source serial element ID!");
    }

    const auto &triangle = *parallel_mesh.GetElement(local_element);
    const auto &source = serial_features.elements[serial_element];
    auto &incidence = result.elements[local_element];
    for (const auto &node : source.nodes)
    {
      if (node.vertex >= serial_features.vertices.size())
      {
        throw std::invalid_argument(
            "Source singular line-tip incidence references an invalid tip!");
      }
      const GlobalVertexId source_vertex = node.mesh_vertex;
      const auto local_vertex = local_vertices.find(source_vertex);
      if (local_vertex == local_vertices.end())
      {
        throw std::invalid_argument(
            "Rank-local triangle is missing an incident singular tip!");
      }
      incidence.nodes.push_back(
          {node.vertex, local_vertex->second,
           CanonicalTriangleNodeTuple(triangle, local_vertex->second, local_vertex_id)});
    }
  }
  return result;
}

}  // namespace singular

}  // namespace fem

}  // namespace palace
