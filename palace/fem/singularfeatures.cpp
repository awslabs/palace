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
using FaceKey = std::array<int, 3>;
using StableEdgeKey = std::array<GlobalVertexId, 2>;
using StableElementKey = std::array<GlobalVertexId, 4>;
using json = nlohmann::json;

struct SelectedEdge
{
  int mesh_edge = -1;
  std::vector<int> boundary_elements;
  std::vector<int> boundary_attributes;
};

FaceKey GetFaceKey(const mfem::Mesh &mesh, int face)
{
  mfem::Array<int> vertices;
  mesh.GetFaceVertices(face, vertices);
  if (vertices.Size() != 3)
  {
    throw std::runtime_error("Singular-feature extraction found an invalid mesh face!");
  }
  FaceKey key{vertices[0], vertices[1], vertices[2]};
  std::sort(key.begin(), key.end());
  if (std::adjacent_find(key.begin(), key.end()) != key.end())
  {
    throw std::runtime_error("Singular-feature extraction found a degenerate mesh face!");
  }
  return key;
}

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

bool IsGeometricallyStraightElementEdge(const mfem::Mesh &mesh, int element, int edge)
{
  if (element < 0 || element >= mesh.GetNE())
  {
    throw std::invalid_argument(
        "Singular-feature edge requires a valid incident mesh element!");
  }

  mfem::Array<int> element_edges, orientations;
  mesh.GetElementEdges(element, element_edges, orientations);
  const int *local_edge = std::find(element_edges.begin(), element_edges.end(), edge);
  if (local_edge == element_edges.end())
  {
    throw std::invalid_argument(
        "Singular-feature edge is not incident to its selected mesh element!");
  }
  const int local_edge_index = static_cast<int>(local_edge - element_edges.begin());

  const auto *mesh_element = mesh.GetElement(element);
  const int geometry = mesh_element->GetGeometryType();
  const auto *reference_vertices = mfem::Geometries.GetVertices(geometry);
  const int *local_vertices = mesh_element->GetEdgeVertices(local_edge_index);
  if (!reference_vertices || !local_vertices || local_vertices[0] < 0 ||
      local_vertices[1] < 0 || local_vertices[0] >= reference_vertices->GetNPoints() ||
      local_vertices[1] >= reference_vertices->GetNPoints())
  {
    throw std::runtime_error(
        "Singular-feature edge has invalid reference-element topology!");
  }

  mfem::IsoparametricTransformation element_transformation;
  mesh.GetElementTransformation(element, &element_transformation);
  const int map_order = std::max(1, element_transformation.Order());
  mfem::H1_SegmentElement segment_element(map_order);
  mfem::IsoparametricTransformation edge_transformation;
  edge_transformation.SetFE(&segment_element);
  mfem::DenseMatrix physical_points(mesh.SpaceDimension(), segment_element.GetDof());
  mfem::Vector physical_point(mesh.SpaceDimension());
  const auto &start = reference_vertices->IntPoint(local_vertices[0]);
  const auto &end = reference_vertices->IntPoint(local_vertices[1]);
  const auto &segment_nodes = segment_element.GetNodes();
  for (int i = 0; i < segment_nodes.GetNPoints(); i++)
  {
    const double t = segment_nodes.IntPoint(i).x;
    mfem::IntegrationPoint point;
    point.x = (1.0 - t) * start.x + t * end.x;
    point.y = (1.0 - t) * start.y + t * end.y;
    point.z = (1.0 - t) * start.z + t * end.z;
    element_transformation.Transform(point, physical_point);
    physical_points.SetCol(i, physical_point);
  }
  edge_transformation.SetPointMat(physical_points);
  return IsGeometricallyStraightSegmentTransformation(edge_transformation);
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
        (segment.type != FeatureSegmentType::SHEET_EDGE &&
         segment.type != FeatureSegmentType::TRANSMISSION_WEDGE) ||
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
        (segment.type != FeatureSegmentType::SHEET_EDGE &&
         segment.type != FeatureSegmentType::TRANSMISSION_WEDGE) ||
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
                                  segment.boundary_attributes,
                                  static_cast<int>(segment.type)});
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
         entry.at(2).get<std::size_t>(), entry.at(3).get<std::vector<int>>(),
         static_cast<FeatureSegmentType>(entry.at(4).get<int>())});
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
        vertex.nu >= 1.0 ||
        (vertex.type != FeatureVertexType::ENDPOINT &&
         vertex.type != FeatureVertexType::CORNER))
    {
      throw std::invalid_argument(
          "Singular line-tip blueprint contains an invalid tip record!");
    }
    if ((vertex.type == FeatureVertexType::ENDPOINT &&
         (vertex.selected_segments.size() != 1 || !vertex.sectors.empty())) ||
        (vertex.type == FeatureVertexType::CORNER &&
         (vertex.selected_segments.size() != 2 || vertex.sectors.empty())))
    {
      throw std::invalid_argument(
          "Singular line-tip blueprint contains inconsistent endpoint/corner metadata!");
    }
    for (std::size_t segment : vertex.selected_segments)
    {
      if (segment >= features.selected_segments.size())
      {
        throw std::invalid_argument(
            "Singular line tip references an invalid selected segment!");
      }
    }
    double opening_angle = 0.0;
    for (const auto &sector : vertex.sectors)
    {
      if (sector.domain_attribute <= 0 || !std::isfinite(sector.angle) ||
          !(sector.angle > 0.0) || !std::isfinite(sector.permittivity) ||
          !(sector.permittivity > 0.0))
      {
        throw std::invalid_argument(
            "Singular PEC corner contains an invalid material sector!");
      }
      opening_angle += sector.angle;
    }
    if (vertex.type == FeatureVertexType::CORNER &&
        (!std::isfinite(opening_angle) || !(opening_angle > 0.0) ||
         !(opening_angle < 6.28318530717958647693)))
    {
      throw std::invalid_argument(
          "Singular PEC corner contains an invalid dielectric opening angle!");
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
    json sectors = json::array();
    for (const auto &sector : vertex.sectors)
    {
      sectors.push_back({sector.domain_attribute, sector.angle, sector.permittivity});
    }
    packed["vertices"].push_back({vertex.id, vertex.mesh_vertex, vertex.selected_segments,
                                  vertex.nu, static_cast<int>(vertex.type),
                                  std::move(sectors)});
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
    TriangleTipVertex vertex{entry.at(0).get<std::size_t>(),
                             entry.at(1).get<int>(),
                             entry.at(2).get<std::vector<std::size_t>>(),
                             entry.at(3).get<double>(),
                             static_cast<FeatureVertexType>(entry.at(4).get<int>()),
                             {}};
    for (const auto &sector : entry.at(5))
    {
      vertex.sectors.push_back({sector.at(0).get<int>(), sector.at(1).get<double>(),
                                sector.at(2).get<double>()});
    }
    features.vertices.push_back(std::move(vertex));
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

long double
EvaluateDirichletWedgeCharacteristic(const std::vector<TriangleWedgeSector> &sectors,
                                     long double nu)
{
  if (!(nu > 0.0L))
  {
    long double limit = 0.0L;
    for (const auto &sector : sectors)
    {
      limit += static_cast<long double>(sector.angle) / sector.permittivity;
    }
    return limit;
  }

  // Propagate [f, epsilon f'/nu] through each constant-material sector. Starting
  // from [0, 1] enforces the first PEC face; the first component at the second
  // face is the characteristic function. Division by nu removes its trivial
  // zero at the origin.
  long double f = 0.0L;
  long double flux = 1.0L;
  for (const auto &sector : sectors)
  {
    const long double phase = nu * sector.angle;
    const long double cosine = std::cos(phase);
    const long double sine = std::sin(phase);
    const long double epsilon = sector.permittivity;
    const long double next_f = cosine * f + sine * flux / epsilon;
    const long double next_flux = -epsilon * sine * f + cosine * flux;
    f = next_f;
    flux = next_flux;
  }
  return f / nu;
}

}  // namespace

double ComputeDirichletWedgeExponent(const std::vector<TriangleWedgeSector> &sectors)
{
  if (sectors.empty())
  {
    throw std::invalid_argument(
        "A Dirichlet material wedge requires at least one angular sector!");
  }
  double opening_angle = 0.0;
  for (const auto &sector : sectors)
  {
    if (sector.domain_attribute <= 0 || !std::isfinite(sector.angle) ||
        !(sector.angle > 0.0) || !std::isfinite(sector.permittivity) ||
        !(sector.permittivity > 0.0))
    {
      throw std::invalid_argument("A Dirichlet material wedge contains an invalid sector!");
    }
    opening_angle += sector.angle;
  }
  if (!std::isfinite(opening_angle) || !(opening_angle > 0.0) ||
      opening_angle > 6.28318530717958647693 + 1.0e-10)
  {
    throw std::invalid_argument(
        "A Dirichlet material wedge has an invalid total opening angle!");
  }

  constexpr int scan_intervals = 16384;
  long double lower = 0.0L;
  long double lower_value = EvaluateDirichletWedgeCharacteristic(sectors, lower);
  for (int interval = 1; interval <= scan_intervals; interval++)
  {
    const long double upper = static_cast<long double>(interval) / scan_intervals;
    const long double upper_value = EvaluateDirichletWedgeCharacteristic(sectors, upper);
    if (!std::isfinite(lower_value) || !std::isfinite(upper_value))
    {
      throw std::runtime_error(
          "Dirichlet material-wedge eigenvalue evaluation produced nonfinite data!");
    }
    if (upper_value == 0.0L || std::signbit(lower_value) != std::signbit(upper_value))
    {
      long double left = lower;
      long double right = upper;
      long double left_value = lower_value;
      for (int iteration = 0; iteration < 100; iteration++)
      {
        const long double middle = 0.5L * (left + right);
        const long double middle_value =
            EvaluateDirichletWedgeCharacteristic(sectors, middle);
        if (middle_value == 0.0L)
        {
          left = right = middle;
          break;
        }
        if (std::signbit(left_value) != std::signbit(middle_value))
        {
          right = middle;
        }
        else
        {
          left = middle;
          left_value = middle_value;
        }
      }
      const double root = static_cast<double>(0.5L * (left + right));
      return root < 1.0 - 1.0e-12 ? root : 1.0;
    }
    lower = upper;
    lower_value = upper_value;
  }
  return 1.0;
}

namespace
{

struct TetrahedronEdgeAngleData
{
  double dot;
  double oriented_area;
  double angle;
};

TetrahedronEdgeAngleData GetTetrahedronEdgeAngleDataAt(const mfem::Mesh &mesh, int element,
                                                       const EdgeKey &edge_vertices,
                                                       double coordinate)
{
  const auto *tetrahedron = mesh.GetElement(element);
  if (!tetrahedron || tetrahedron->GetGeometryType() != mfem::Geometry::TETRAHEDRON ||
      !std::isfinite(coordinate) || coordinate < 0.0 || coordinate > 1.0)
  {
    throw std::invalid_argument(
        "A PEC edge-wedge fan contains an invalid tetrahedron or edge coordinate!");
  }
  std::array<int, 2> edge_local{-1, -1};
  std::array<int, 2> opposite_local{-1, -1};
  int opposite_count = 0;
  for (int local = 0; local < 4; local++)
  {
    const int vertex = tetrahedron->GetVertices()[local];
    if (vertex == edge_vertices[0])
    {
      edge_local[0] = local;
    }
    else if (vertex == edge_vertices[1])
    {
      edge_local[1] = local;
    }
    else
    {
      opposite_local[opposite_count++] = local;
    }
  }
  if (edge_local[0] < 0 || edge_local[1] < 0 || opposite_count != 2)
  {
    throw std::invalid_argument(
        "A PEC edge-wedge tetrahedron does not contain the selected mesh edge!");
  }

  const auto *reference_vertices =
      mfem::Geometries.GetVertices(mfem::Geometry::TETRAHEDRON);
  if (!reference_vertices || reference_vertices->GetNPoints() != 4)
  {
    throw std::runtime_error("The reference tetrahedron has invalid vertex data!");
  }
  const auto &start = reference_vertices->IntPoint(edge_local[0]);
  const auto &end = reference_vertices->IntPoint(edge_local[1]);
  mfem::IntegrationPoint point;
  point.x = (1.0 - coordinate) * start.x + coordinate * end.x;
  point.y = (1.0 - coordinate) * start.y + coordinate * end.y;
  point.z = (1.0 - coordinate) * start.z + coordinate * end.z;

  mfem::IsoparametricTransformation transformation;
  mesh.GetElementTransformation(element, &transformation);
  mfem::Vector physical_start(3), physical_end(3);
  transformation.Transform(start, physical_start);
  transformation.Transform(end, physical_end);
  mfem::Vector tangent(physical_end);
  tangent -= physical_start;
  const double tangent_norm = tangent.Norml2();
  if (!std::isfinite(tangent_norm) || !(tangent_norm > 0.0))
  {
    throw std::invalid_argument(
        "A PEC edge-wedge tetrahedron has a degenerate edge tangent!");
  }
  tangent /= tangent_norm;

  transformation.SetIntPoint(&point);
  const auto &jacobian = transformation.Jacobian();
  if (jacobian.Height() != 3 || jacobian.Width() != 3)
  {
    throw std::runtime_error(
        "A PEC edge-wedge tetrahedron has an invalid physical Jacobian!");
  }
  const auto physical_direction = [&](const mfem::IntegrationPoint &target)
  {
    mfem::Vector reference(3), physical(3);
    reference[0] = target.x - point.x;
    reference[1] = target.y - point.y;
    reference[2] = target.z - point.z;
    jacobian.Mult(reference, physical);
    return physical;
  };

  std::array<mfem::Vector, 2> radial{
      physical_direction(reference_vertices->IntPoint(opposite_local[0])),
      physical_direction(reference_vertices->IntPoint(opposite_local[1]))};
  for (auto &direction : radial)
  {
    direction.Add(-(direction * tangent), tangent);
    const double norm = direction.Norml2();
    if (!std::isfinite(norm) || !(norm > 0.0))
    {
      throw std::invalid_argument(
          "A PEC edge-wedge tetrahedron is degenerate normal to its edge!");
    }
  }
  mfem::Vector cross(3);
  cross[0] = radial[0][1] * radial[1][2] - radial[0][2] * radial[1][1];
  cross[1] = radial[0][2] * radial[1][0] - radial[0][0] * radial[1][2];
  cross[2] = radial[0][0] * radial[1][1] - radial[0][1] * radial[1][0];
  const double dot = radial[0] * radial[1];
  const double oriented_area = cross * tangent;
  const double angle = std::atan2(std::abs(oriented_area), dot);
  if (!std::isfinite(angle) || !(angle > 0.0) || !(angle < 3.14159265358979323846))
  {
    throw std::invalid_argument(
        "A PEC edge-wedge tetrahedron has an invalid cross-sectional angle!");
  }
  return {dot, oriented_area, angle};
}

long double BinomialCoefficient(int n, int k)
{
  if (n < 0 || k < 0 || k > n)
  {
    throw std::invalid_argument("Invalid Bernstein polynomial degree!");
  }
  k = std::min(k, n - k);
  long double value = 1.0L;
  for (int i = 1; i <= k; i++)
  {
    value *= static_cast<long double>(n - k + i) / i;
  }
  return value;
}

mfem::Vector MultiplyBernsteinPolynomials(const mfem::Vector &left,
                                          const mfem::Vector &right)
{
  if (left.Size() < 1 || right.Size() < 1)
  {
    throw std::invalid_argument("Cannot multiply empty Bernstein polynomials!");
  }
  const int left_order = left.Size() - 1;
  const int right_order = right.Size() - 1;
  mfem::Vector product(left_order + right_order + 1);
  product = 0.0;
  for (int i = 0; i <= left_order; i++)
  {
    for (int j = 0; j <= right_order; j++)
    {
      const int k = i + j;
      const long double weight = BinomialCoefficient(left_order, i) *
                                 BinomialCoefficient(right_order, j) /
                                 BinomialCoefficient(left_order + right_order, k);
      product[k] += static_cast<double>(weight * left[i] * right[j]);
    }
  }
  return product;
}

mfem::Vector DifferentiateBernsteinPolynomial(const mfem::Vector &polynomial)
{
  if (polynomial.Size() < 2)
  {
    throw std::invalid_argument("Cannot differentiate a constant Bernstein polynomial!");
  }
  const int order = polynomial.Size() - 1;
  mfem::Vector derivative(order);
  for (int i = 0; i < order; i++)
  {
    derivative[i] = order * (polynomial[i + 1] - polynomial[i]);
  }
  return derivative;
}

bool HasStrictlyPositiveBernsteinPolynomial(const mfem::Vector &polynomial,
                                            double tolerance, int depth = 0)
{
  constexpr int maximum_depth = 30;
  if (polynomial.Size() < 1 || !std::isfinite(tolerance) || tolerance < 0.0)
  {
    return false;
  }
  double minimum = std::numeric_limits<double>::infinity();
  double maximum = -std::numeric_limits<double>::infinity();
  for (double coefficient : polynomial)
  {
    if (!std::isfinite(coefficient))
    {
      return false;
    }
    minimum = std::min(minimum, coefficient);
    maximum = std::max(maximum, coefficient);
  }
  if (minimum > tolerance)
  {
    return true;
  }
  if (maximum <= tolerance || depth == maximum_depth)
  {
    return false;
  }
  mfem::Vector left, right;
  detail::SplitBernsteinCoefficients(polynomial, left, right);
  return HasStrictlyPositiveBernsteinPolynomial(left, tolerance, depth + 1) &&
         HasStrictlyPositiveBernsteinPolynomial(right, tolerance, depth + 1);
}

double BernsteinCoefficientScale(const mfem::Vector &polynomial)
{
  double scale = 0.0;
  for (double coefficient : polynomial)
  {
    scale = std::max(scale, std::abs(coefficient));
  }
  return scale;
}

double GetTetrahedronEdgeAngle(const mfem::Mesh &mesh, int element,
                               const EdgeKey &edge_vertices)
{
  mfem::IsoparametricTransformation transformation;
  mesh.GetElementTransformation(element, &transformation);
  const auto midpoint = GetTetrahedronEdgeAngleDataAt(mesh, element, edge_vertices, 0.5);
  if (transformation.OrderJ() == 0)
  {
    return midpoint.angle;
  }

  // With a globally straight physical edge, projection onto its fixed normal
  // plane makes the two radial vectors polynomial in the edge coordinate. Their
  // dot product a(t) and oriented area b(t) have degree at most twice the map
  // order. The wedge angle satisfies
  //
  //   theta'(t) = (a b' - b a') / (a^2 + b^2).
  //
  // Bernstein bounds therefore certify a constant cross-sectional angle over
  // the complete edge, rather than only at a finite set of sample points.
  const int map_order = std::max(1, transformation.Order());
  const int invariant_order = 2 * map_order;
  const int coefficient_count = invariant_order + 1;
  mfem::DenseMatrix bernstein_values(coefficient_count);
  mfem::DenseMatrix invariant_values(coefficient_count, 2);
  mfem::Vector bernstein_shape(coefficient_count);
  const double orientation = std::signbit(midpoint.oriented_area) ? -1.0 : 1.0;
  for (int q = 0; q < coefficient_count; q++)
  {
    const double coordinate = static_cast<double>(q) / invariant_order;
    mfem::Poly_1D::CalcBernstein(invariant_order, coordinate, bernstein_shape);
    bernstein_values.SetRow(q, bernstein_shape);
    const auto data =
        GetTetrahedronEdgeAngleDataAt(mesh, element, edge_vertices, coordinate);
    invariant_values(q, 0) = data.dot;
    invariant_values(q, 1) = orientation * data.oriented_area;
  }
  mfem::DenseMatrix invariant_control;
  mfem::DenseMatrixInverse(bernstein_values).Mult(invariant_values, invariant_control);
  mfem::Vector dot_control(coefficient_count), area_control(coefficient_count);
  invariant_control.GetColumn(0, dot_control);
  invariant_control.GetColumn(1, area_control);

  // Independently verify the polynomial recovery. An ill-conditioned nodal to
  // Bernstein conversion is rejected instead of weakening the global bound.
  const auto &rule = mfem::IntRules.Get(mfem::Geometry::SEGMENT, 2 * invariant_order + 2);
  double invariant_scale = std::max(BernsteinCoefficientScale(dot_control),
                                    BernsteinCoefficientScale(area_control));
  for (int q = -2; q < rule.GetNPoints(); q++)
  {
    const double coordinate = q < 0 ? (q == -2 ? 0.0 : 1.0) : rule.IntPoint(q).x;
    const auto data =
        GetTetrahedronEdgeAngleDataAt(mesh, element, edge_vertices, coordinate);
    mfem::Poly_1D::CalcBernstein(invariant_order, coordinate, bernstein_shape);
    double reconstructed_dot = 0.0;
    double reconstructed_area = 0.0;
    for (int i = 0; i < coefficient_count; i++)
    {
      reconstructed_dot += bernstein_shape[i] * dot_control[i];
      reconstructed_area += bernstein_shape[i] * area_control[i];
    }
    const double actual_area = orientation * data.oriented_area;
    invariant_scale =
        std::max({invariant_scale, std::abs(data.dot), std::abs(actual_area)});
    const double representation_tolerance =
        1.0e-10 * std::max(invariant_scale, std::numeric_limits<double>::min());
    if (std::abs(reconstructed_dot - data.dot) > representation_tolerance ||
        std::abs(reconstructed_area - actual_area) > representation_tolerance)
    {
      throw std::invalid_argument(
          "A curved PEC edge wedge could not certify its cross-sectional angle "
          "polynomial!");
    }
  }

  const double area_tolerance =
      8192.0 * std::numeric_limits<double>::epsilon() *
      std::max(BernsteinCoefficientScale(area_control), std::numeric_limits<double>::min());
  if (!HasStrictlyPositiveBernsteinPolynomial(area_control, area_tolerance))
  {
    throw std::invalid_argument(
        "A curved PEC edge wedge has a degenerate or reversing cross section!");
  }

  const auto dot_derivative = DifferentiateBernsteinPolynomial(dot_control);
  const auto area_derivative = DifferentiateBernsteinPolynomial(area_control);
  auto angle_numerator = MultiplyBernsteinPolynomials(dot_control, area_derivative);
  const auto second_numerator = MultiplyBernsteinPolynomials(area_control, dot_derivative);
  for (int i = 0; i < angle_numerator.Size(); i++)
  {
    angle_numerator[i] -= second_numerator[i];
  }
  auto angle_denominator = MultiplyBernsteinPolynomials(dot_control, dot_control);
  const auto area_squared = MultiplyBernsteinPolynomials(area_control, area_control);
  for (int i = 0; i < angle_denominator.Size(); i++)
  {
    angle_denominator[i] += area_squared[i];
  }
  const mfem::Vector constant_one({1.0, 1.0});
  const auto elevated_numerator =
      MultiplyBernsteinPolynomials(angle_numerator, constant_one);
  if (elevated_numerator.Size() != angle_denominator.Size())
  {
    throw std::logic_error(
        "PEC edge-wedge angle certification produced inconsistent polynomial degrees!");
  }

  const double derivative_tolerance = 1.0e-9 * std::max(1.0, midpoint.angle);
  mfem::Vector positive_margin(angle_denominator.Size()),
      negative_margin(angle_denominator.Size());
  for (int i = 0; i < angle_denominator.Size(); i++)
  {
    positive_margin[i] =
        derivative_tolerance * angle_denominator[i] + elevated_numerator[i];
    negative_margin[i] =
        derivative_tolerance * angle_denominator[i] - elevated_numerator[i];
  }
  const double margin_scale = std::max(BernsteinCoefficientScale(positive_margin),
                                       BernsteinCoefficientScale(negative_margin));
  const double margin_tolerance =
      8192.0 * std::numeric_limits<double>::epsilon() *
      std::max(margin_scale, std::numeric_limits<double>::min());
  if (!HasStrictlyPositiveBernsteinPolynomial(positive_margin, margin_tolerance) ||
      !HasStrictlyPositiveBernsteinPolynomial(negative_margin, margin_tolerance))
  {
    throw std::invalid_argument(
        "A curved PEC edge wedge has a cross-sectional angle which varies along the "
        "feature!");
  }
  return midpoint.angle;
}

std::vector<TriangleWedgeSector>
BuildTetrahedronEdgeSectors(const mfem::Mesh &mesh, int mesh_edge,
                            const EdgeKey &edge_vertices,
                            const std::vector<int> &selected_boundary_elements,
                            const std::map<int, double> &permittivity)
{
  if (selected_boundary_elements.size() != 2)
  {
    throw std::invalid_argument(
        "A one-sided PEC wedge edge requires exactly two selected boundary faces!");
  }

  struct FanTetrahedron
  {
    int element;
    std::array<FaceKey, 2> radial_faces;
    double angle;
    int domain_attribute;
    double permittivity;
  };
  std::vector<FanTetrahedron> fan;
  std::map<FaceKey, std::vector<std::size_t>> face_tetrahedra;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    const auto *tetrahedron = mesh.GetElement(element);
    const int *vertices = tetrahedron->GetVertices();
    if (std::find(vertices, vertices + 4, edge_vertices[0]) == vertices + 4 ||
        std::find(vertices, vertices + 4, edge_vertices[1]) == vertices + 4)
    {
      continue;
    }
    if (!IsGeometricallyStraightElementEdge(mesh, element, mesh_edge))
    {
      throw std::invalid_argument(
          "Selected finite-metal wedge edges must be geometrically straight!");
    }
    std::array<int, 2> opposite_vertices;
    int next = 0;
    for (int local = 0; local < 4; local++)
    {
      if (vertices[local] != edge_vertices[0] && vertices[local] != edge_vertices[1])
      {
        opposite_vertices[next++] = vertices[local];
      }
    }
    if (next != 2 || opposite_vertices[0] == opposite_vertices[1])
    {
      throw std::invalid_argument(
          "A PEC edge-wedge tetrahedron has invalid radial-face topology!");
    }
    std::array<FaceKey, 2> radial_faces{
        FaceKey{edge_vertices[0], edge_vertices[1], opposite_vertices[0]},
        FaceKey{edge_vertices[0], edge_vertices[1], opposite_vertices[1]}};
    std::sort(radial_faces[0].begin(), radial_faces[0].end());
    std::sort(radial_faces[1].begin(), radial_faces[1].end());
    const int attribute = mesh.GetAttribute(element);
    const auto material = permittivity.find(attribute);
    if (material == permittivity.end())
    {
      throw std::invalid_argument(
          "A one-sided PEC wedge edge has no isotropic permittivity for domain "
          "attribute " +
          std::to_string(attribute) + "!");
    }
    const std::size_t index = fan.size();
    fan.push_back({element, radial_faces,
                   GetTetrahedronEdgeAngle(mesh, element, edge_vertices), attribute,
                   material->second});
    face_tetrahedra[radial_faces[0]].push_back(index);
    face_tetrahedra[radial_faces[1]].push_back(index);
  }
  if (fan.empty())
  {
    throw std::invalid_argument(
        "A one-sided PEC wedge edge has no incident tetrahedral fan!");
  }

  std::array<FaceKey, 2> boundary_faces;
  for (int side = 0; side < 2; side++)
  {
    int face, orientation;
    mesh.GetBdrElementFace(selected_boundary_elements[side], &face, &orientation);
    boundary_faces[side] = GetFaceKey(mesh, face);
    const auto incidence = face_tetrahedra.find(boundary_faces[side]);
    if (incidence == face_tetrahedra.end() || incidence->second.size() != 1)
    {
      throw std::invalid_argument(
          "A one-sided PEC wedge boundary face does not bound exactly one "
          "tetrahedron!");
    }
  }
  if (boundary_faces[0] == boundary_faces[1])
  {
    throw std::invalid_argument("A one-sided PEC wedge edge has duplicate boundary faces!");
  }

  std::vector<TriangleWedgeSector> sectors;
  std::set<std::size_t> visited;
  FaceKey current_face = boundary_faces[0];
  std::size_t current_tetrahedron = face_tetrahedra.at(current_face)[0];
  while (true)
  {
    if (!visited.insert(current_tetrahedron).second)
    {
      throw std::invalid_argument(
          "A one-sided PEC edge-wedge fan contains a cycle or branch!");
    }
    const auto &tetrahedron = fan[current_tetrahedron];
    if (!sectors.empty() &&
        sectors.back().domain_attribute == tetrahedron.domain_attribute &&
        sectors.back().permittivity == tetrahedron.permittivity)
    {
      sectors.back().angle += tetrahedron.angle;
    }
    else
    {
      sectors.push_back(
          {tetrahedron.domain_attribute, tetrahedron.angle, tetrahedron.permittivity});
    }

    FaceKey next_face;
    if (tetrahedron.radial_faces[0] == current_face)
    {
      next_face = tetrahedron.radial_faces[1];
    }
    else if (tetrahedron.radial_faces[1] == current_face)
    {
      next_face = tetrahedron.radial_faces[0];
    }
    else
    {
      throw std::logic_error("A one-sided PEC edge-wedge fan walk lost its incoming face!");
    }
    if (next_face == boundary_faces[1])
    {
      break;
    }
    const auto incidence = face_tetrahedra.find(next_face);
    if (incidence == face_tetrahedra.end() || incidence->second.size() != 2)
    {
      throw std::invalid_argument(
          "A one-sided PEC edge-wedge fan is interrupted by another exterior face!");
    }
    current_tetrahedron = incidence->second[0] == current_tetrahedron
                              ? incidence->second[1]
                              : incidence->second[0];
    current_face = next_face;
  }
  if (visited.size() != fan.size())
  {
    throw std::invalid_argument(
        "A one-sided PEC wedge edge has disconnected or nonmanifold incident "
        "tetrahedra!");
  }
  double opening_angle = 0.0;
  for (const auto &sector : sectors)
  {
    opening_angle += sector.angle;
  }
  if (!std::isfinite(opening_angle) || !(opening_angle > 0.0) ||
      !(opening_angle < 6.28318530717958647693 - 1.0e-10))
  {
    throw std::invalid_argument(
        "A one-sided PEC wedge edge has an invalid dielectric opening angle!");
  }
  return sectors;
}

double GetTriangleCornerAngle(const mfem::Mesh &mesh, int element, int vertex)
{
  const auto *triangle = mesh.GetElement(element);
  if (!triangle || triangle->GetGeometryType() != mfem::Geometry::TRIANGLE)
  {
    throw std::invalid_argument(
        "A PEC material-wedge fan contains a nontriangular element!");
  }
  int local_vertex = -1;
  for (int local = 0; local < 3; local++)
  {
    if (triangle->GetVertices()[local] == vertex)
    {
      local_vertex = local;
      break;
    }
  }
  if (local_vertex < 0)
  {
    throw std::invalid_argument(
        "A PEC material-wedge fan does not contain its corner vertex!");
  }

  constexpr std::array<std::array<double, 2>, 3> reference_vertices{
      std::array<double, 2>{0.0, 0.0}, std::array<double, 2>{1.0, 0.0},
      std::array<double, 2>{0.0, 1.0}};
  mfem::IntegrationPoint point;
  point.Set2(reference_vertices[local_vertex][0], reference_vertices[local_vertex][1]);
  mfem::IsoparametricTransformation transformation;
  mesh.GetElementTransformation(element, &transformation);
  transformation.SetIntPoint(&point);
  const auto &jacobian = transformation.Jacobian();
  if (jacobian.Height() != 2 || jacobian.Width() != 2)
  {
    throw std::runtime_error(
        "A PEC material-wedge triangle has an invalid physical Jacobian!");
  }

  std::array<std::array<double, 2>, 2> tangent;
  int next = 0;
  for (int local = 0; local < 3; local++)
  {
    if (local == local_vertex)
    {
      continue;
    }
    mfem::Vector reference_direction(2), physical_direction(2);
    for (int d = 0; d < 2; d++)
    {
      reference_direction[d] =
          reference_vertices[local][d] - reference_vertices[local_vertex][d];
    }
    jacobian.Mult(reference_direction, physical_direction);
    const double norm = physical_direction.Norml2();
    if (!std::isfinite(norm) || !(norm > 0.0))
    {
      throw std::invalid_argument(
          "A PEC material-wedge triangle is degenerate at its corner!");
    }
    tangent[next++] = {physical_direction[0] / norm, physical_direction[1] / norm};
  }
  const double dot = tangent[0][0] * tangent[1][0] + tangent[0][1] * tangent[1][1];
  const double cross = tangent[0][0] * tangent[1][1] - tangent[0][1] * tangent[1][0];
  const double angle = std::atan2(std::abs(cross), std::clamp(dot, -1.0, 1.0));
  if (!std::isfinite(angle) || !(angle > 0.0) || !(angle < 3.14159265358979323846))
  {
    throw std::invalid_argument(
        "A PEC material-wedge triangle has an invalid corner angle!");
  }
  return angle;
}

std::vector<TriangleWedgeSector>
BuildTriangleCornerSectors(const mfem::Mesh &mesh, int vertex,
                           const std::vector<std::size_t> &selected_segment_indices,
                           const std::vector<TriangleSelectedSegment> &selected_segments,
                           const std::vector<int> &incident_elements,
                           const std::map<int, double> &permittivity)
{
  if (selected_segment_indices.size() != 2 || incident_elements.empty())
  {
    throw std::invalid_argument(
        "A one-sided PEC corner requires two boundary segments and an element fan!");
  }

  struct FanTriangle
  {
    int element;
    std::array<EdgeKey, 2> rays;
    double angle;
    int domain_attribute;
    double permittivity;
  };
  std::vector<FanTriangle> fan;
  fan.reserve(incident_elements.size());
  std::map<EdgeKey, std::vector<std::size_t>> ray_triangles;
  for (int element : incident_elements)
  {
    const auto *triangle = mesh.GetElement(element);
    if (!triangle || triangle->GetGeometryType() != mfem::Geometry::TRIANGLE)
    {
      throw std::invalid_argument(
          "A one-sided PEC corner contains a nontriangular element!");
    }
    std::array<EdgeKey, 2> rays;
    int next = 0;
    for (int local = 0; local < 3; local++)
    {
      const int other = triangle->GetVertices()[local];
      if (other == vertex)
      {
        continue;
      }
      rays[next] = {vertex, other};
      std::sort(rays[next].begin(), rays[next].end());
      next++;
    }
    if (next != 2 || rays[0] == rays[1])
    {
      throw std::invalid_argument(
          "A one-sided PEC corner has invalid radial-edge topology!");
    }
    const int attribute = mesh.GetAttribute(element);
    const auto material = permittivity.find(attribute);
    if (material == permittivity.end())
    {
      throw std::invalid_argument(
          "A one-sided PEC corner has no isotropic permittivity for domain attribute " +
          std::to_string(attribute) + "!");
    }
    const std::size_t index = fan.size();
    fan.push_back({element, rays, GetTriangleCornerAngle(mesh, element, vertex), attribute,
                   material->second});
    ray_triangles[rays[0]].push_back(index);
    ray_triangles[rays[1]].push_back(index);
  }

  std::array<EdgeKey, 2> boundary_rays;
  for (int side = 0; side < 2; side++)
  {
    const std::size_t segment = selected_segment_indices[side];
    if (segment >= selected_segments.size())
    {
      throw std::invalid_argument(
          "A one-sided PEC corner references an invalid boundary segment!");
    }
    boundary_rays[side] = selected_segments[segment].mesh_vertices;
    const auto incidence = ray_triangles.find(boundary_rays[side]);
    if (incidence == ray_triangles.end() || incidence->second.size() != 1)
    {
      throw std::invalid_argument(
          "A one-sided PEC corner boundary ray does not bound exactly one triangle!");
    }
  }
  if (boundary_rays[0] == boundary_rays[1])
  {
    throw std::invalid_argument("A one-sided PEC corner has duplicate boundary rays!");
  }

  std::vector<TriangleWedgeSector> sectors;
  std::set<std::size_t> visited;
  EdgeKey current_ray = boundary_rays[0];
  std::size_t current_triangle = ray_triangles.at(current_ray)[0];
  while (true)
  {
    if (!visited.insert(current_triangle).second)
    {
      throw std::invalid_argument(
          "A one-sided PEC corner element fan contains a cycle or branch!");
    }
    const auto &triangle = fan[current_triangle];
    if (!sectors.empty() && sectors.back().domain_attribute == triangle.domain_attribute &&
        sectors.back().permittivity == triangle.permittivity)
    {
      sectors.back().angle += triangle.angle;
    }
    else
    {
      sectors.push_back({triangle.domain_attribute, triangle.angle, triangle.permittivity});
    }

    EdgeKey next_ray;
    if (triangle.rays[0] == current_ray)
    {
      next_ray = triangle.rays[1];
    }
    else if (triangle.rays[1] == current_ray)
    {
      next_ray = triangle.rays[0];
    }
    else
    {
      throw std::logic_error(
          "A one-sided PEC corner fan walk lost its incoming radial edge!");
    }
    if (next_ray == boundary_rays[1])
    {
      break;
    }
    const auto incidence = ray_triangles.find(next_ray);
    if (incidence == ray_triangles.end() || incidence->second.size() != 2)
    {
      throw std::invalid_argument(
          "A one-sided PEC corner fan is interrupted by another exterior boundary!");
    }
    current_triangle = incidence->second[0] == current_triangle ? incidence->second[1]
                                                                : incidence->second[0];
    current_ray = next_ray;
  }
  if (visited.size() != fan.size())
  {
    throw std::invalid_argument(
        "A one-sided PEC corner has disconnected or nonmanifold incident triangles!");
  }

  double opening_angle = 0.0;
  for (const auto &sector : sectors)
  {
    opening_angle += sector.angle;
  }
  if (!std::isfinite(opening_angle) || !(opening_angle > 0.0) ||
      !(opening_angle < 6.28318530717958647693 - 1.0e-10))
  {
    throw std::invalid_argument(
        "A one-sided PEC corner has an invalid dielectric opening angle!");
  }
  return sectors;
}

}  // namespace

TriangleFeatureTopology ExtractSerialLineFeatures(
    const mfem::Mesh &mesh, const std::vector<int> &boundary_attributes,
    const std::vector<TriangleMaterial> &materials, double line_tip_nu)
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
  if (!std::isfinite(line_tip_nu) || line_tip_nu <= 0.0 || line_tip_nu >= 1.0)
  {
    throw std::invalid_argument(
        "Singular-line exponent must be finite and satisfy 0 < nu < 1!");
  }
  std::map<int, double> material_permittivity;
  for (const auto &material : materials)
  {
    if (material.domain_attribute <= 0 || !std::isfinite(material.permittivity) ||
        !(material.permittivity > 0.0) ||
        !material_permittivity.emplace(material.domain_attribute, material.permittivity)
             .second)
    {
      throw std::invalid_argument(
          "Singular-line extraction received invalid or duplicate isotropic material "
          "data!");
    }
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
  std::map<EdgeKey, bool> selected_edge_internal;
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
    if (element_1 < 0)
    {
      throw std::invalid_argument(
          "A selected PEC line segment has no incident mesh element!");
    }
    const bool internal = element_2 >= 0;
    if (!internal && material_permittivity.empty())
    {
      throw std::invalid_argument(
          "Selected zero-thickness PEC lines must be internal mesh boundaries!");
    }
    if (!IsGeometricallyStraightElementEdge(mesh, element_1, edge))
    {
      throw std::invalid_argument(
          "Selected singular PEC line segments must be geometrically straight!");
    }
    const auto key = GetEdgeKey(mesh, edge);
    if (!selected_edges
             .emplace(key, TriangleSelectedSegment{boundary_element, edge, key, attribute})
             .second)
    {
      throw std::invalid_argument(
          "A selected singular PEC segment is represented more than once!");
    }
    selected_edge_internal.emplace(key, internal);
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

  std::vector<std::vector<int>> vertex_elements(mesh.GetNV());
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    const auto *triangle = mesh.GetElement(element);
    for (int local = 0; local < 3; local++)
    {
      vertex_elements[triangle->GetVertices()[local]].push_back(element);
    }
  }

  constexpr double straight_tolerance = 1.0e-10;
  for (const auto &[vertex, segments] : vertex_segments)
  {
    if (segments.size() > 2)
    {
      throw std::invalid_argument(
          "Selected singular PEC lines contain a branch or junction!");
    }
    bool has_internal = false;
    bool has_one_sided = false;
    for (std::size_t segment : segments)
    {
      const auto &selected = result.selected_segments[segment];
      const bool internal = selected_edge_internal.at(selected.mesh_vertices);
      has_internal = has_internal || internal;
      has_one_sided = has_one_sided || !internal;
    }
    if (has_internal && has_one_sided)
    {
      throw std::invalid_argument(
          "A singular PEC vertex mixes internal-sheet and one-sided boundary segments!");
    }

    if (has_internal && segments.size() == 2)
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
            "Selected internal zero-thickness PEC lines must be straight!");
      }
      continue;
    }

    if (has_internal)
    {
      if (exterior_vertices.find(vertex) == exterior_vertices.end())
      {
        result.vertices.push_back({result.vertices.size(),
                                   vertex,
                                   segments,
                                   line_tip_nu,
                                   FeatureVertexType::ENDPOINT,
                                   {}});
      }
      continue;
    }

    // A single selected one-sided segment ends where another, unselected boundary
    // condition begins. That mixed-boundary corner is intentionally outside the
    // selected finite conductor and is not enriched.
    if (segments.size() == 1)
    {
      continue;
    }

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
    if (dot <= -1.0 + straight_tolerance)
    {
      continue;
    }

    auto sectors =
        BuildTriangleCornerSectors(mesh, vertex, segments, result.selected_segments,
                                   vertex_elements[vertex], material_permittivity);
    const double corner_nu = ComputeDirichletWedgeExponent(sectors);
    if (corner_nu < 1.0)
    {
      result.vertices.push_back({result.vertices.size(), vertex, segments, corner_nu,
                                 FeatureVertexType::CORNER, std::move(sectors)});
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

TriangleFeatureTopology
ExtractSerialLineTipFeatures(const mfem::Mesh &mesh,
                             const std::vector<int> &boundary_attributes, double nu)
{
  return ExtractSerialLineFeatures(mesh, boundary_attributes, {}, nu);
}

FeatureTopology ExtractSerialSheetFeatures(const mfem::Mesh &mesh,
                                           const std::vector<int> &boundary_attributes,
                                           double nu)
{
  return ExtractSerialSheetFeatures(mesh, boundary_attributes, {}, nu);
}

FeatureTopology ExtractSerialSheetFeatures(const mfem::Mesh &mesh,
                                           const std::vector<int> &boundary_attributes,
                                           const std::vector<TriangleMaterial> &materials,
                                           double nu)
{
  ValidateMesh(mesh);
  if (!std::isfinite(nu) || nu <= 0.0 || nu >= 1.0)
  {
    throw std::invalid_argument(
        "Singular-feature exponent must be finite and satisfy 0 < nu < 1!");
  }
  std::map<int, double> material_permittivity;
  for (const auto &material : materials)
  {
    if (material.domain_attribute <= 0 || !std::isfinite(material.permittivity) ||
        !(material.permittivity > 0.0) ||
        !material_permittivity.emplace(material.domain_attribute, material.permittivity)
             .second)
    {
      throw std::invalid_argument(
          "Singular sheet/wedge extraction received invalid or duplicate isotropic "
          "material data!");
    }
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
  std::map<int, bool> selected_face_internal;
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
    if (element_1 < 0)
    {
      throw std::invalid_argument("A selected PEC face has no incident mesh element!");
    }
    const bool internal = element_2 >= 0;
    if (!internal && material_permittivity.empty())
    {
      throw std::invalid_argument(
          "Selected one-sided finite-metal PEC faces require isotropic material data!");
    }
    if (!selected_faces.emplace(face, boundary_element).second)
    {
      throw std::invalid_argument(
          "A selected zero-thickness PEC face is represented more than once!");
    }
    selected_face_internal.emplace(face, internal);

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

  std::vector<double> segment_exponents;
  for (auto &[key, edge] : selected_edges)
  {
    bool has_internal_face = false;
    bool has_external_face = false;
    for (int boundary_element : edge.boundary_elements)
    {
      int face, orientation;
      mesh.GetBdrElementFace(boundary_element, &face, &orientation);
      const bool internal = selected_face_internal.at(face);
      has_internal_face = has_internal_face || internal;
      has_external_face = has_external_face || !internal;
    }
    if (has_internal_face && has_external_face)
    {
      throw std::invalid_argument(
          "One selected PEC mesh edge mixes internal-sheet and one-sided boundary "
          "faces!");
    }

    double edge_nu = 1.0;
    if (has_internal_face)
    {
      if (edge.boundary_elements.size() != 1)
      {
        continue;
      }
      edge_nu = nu;
    }
    else
    {
      if (edge.boundary_elements.size() != 2)
      {
        // A single selected one-sided face ends at an unselected boundary
        // condition. That mixed-condition edge is outside the selected finite
        // conductor wedge.
        continue;
      }
      const auto sectors = BuildTetrahedronEdgeSectors(
          mesh, edge.mesh_edge, key, edge.boundary_elements, material_permittivity);
      edge_nu = ComputeDirichletWedgeExponent(sectors);
      if (!(edge_nu < 1.0))
      {
        continue;
      }
    }

    if (!std::isfinite(edge_nu) || !(edge_nu > 0.0) || !(edge_nu < 1.0))
    {
      throw std::runtime_error("Selected PEC edge produced an invalid singular exponent!");
    }
    std::sort(edge.boundary_attributes.begin(), edge.boundary_attributes.end());
    edge.boundary_attributes.erase(
        std::unique(edge.boundary_attributes.begin(), edge.boundary_attributes.end()),
        edge.boundary_attributes.end());
    int face, orientation, element_1, element_2;
    mesh.GetBdrElementFace(edge.boundary_elements.front(), &face, &orientation);
    mesh.GetFaceElements(face, &element_1, &element_2);
    if (element_1 < 0 ||
        !IsGeometricallyStraightElementEdge(mesh, element_1, edge.mesh_edge) ||
        (element_2 >= 0 &&
         !IsGeometricallyStraightElementEdge(mesh, element_2, edge.mesh_edge)))
    {
      throw std::invalid_argument(
          "Selected PEC sheet or wedge edges must be geometrically straight!");
    }
    result.segments.push_back({edge.mesh_edge, key, 0, std::move(edge.boundary_attributes),
                               has_internal_face ? FeatureSegmentType::SHEET_EDGE
                                                 : FeatureSegmentType::TRANSMISSION_WEDGE});
    segment_exponents.push_back(edge_nu);
  }

  if (segment_exponents.size() != result.segments.size())
  {
    throw std::logic_error("Singular-feature extraction lost a segment exponent!");
  }

  std::map<int, std::size_t> vertex_index;
  for (std::size_t segment = 0; segment < result.segments.size(); segment++)
  {
    for (int vertex : result.segments[segment].mesh_vertices)
    {
      auto [it, inserted] = vertex_index.emplace(vertex, result.vertices.size());
      if (inserted)
      {
        result.vertices.push_back({result.vertices.size(),
                                   vertex,
                                   FeatureVertexType::REGULAR,
                                   {},
                                   {},
                                   segment_exponents[segment]});
      }
      else
      {
        result.vertices[it->second].nu =
            std::min(result.vertices[it->second].nu, segment_exponents[segment]);
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
      continue;
    }
    if (vertex.segments.size() == 2)
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
      const double nu_0 = segment_exponents[vertex.segments[0]];
      const double nu_1 = segment_exponents[vertex.segments[1]];
      constexpr double straight_tolerance = 1.0e-10;
      constexpr double exponent_tolerance = 1.0e-12;
      const bool same_exponent =
          std::abs(nu_0 - nu_1) <= exponent_tolerance * std::max({1.0, nu_0, nu_1});
      if (!same_exponent)
      {
        throw std::invalid_argument(
            "A selected PEC wedge edge changes singular exponent at a mesh vertex; "
            "three-dimensional point-transition enrichment is not implemented!");
      }
      vertex.type = dot <= -1.0 + straight_tolerance ? FeatureVertexType::REGULAR
                                                     : FeatureVertexType::CORNER;
      continue;
    }
    throw std::invalid_argument(
        "Selected PEC sheet or wedge edges contain an unsupported junction!");
  }

  std::vector<int> segment_feature(result.segments.size(), -1);
  for (std::size_t seed = 0; seed < result.segments.size(); seed++)
  {
    if (segment_feature[seed] >= 0)
    {
      continue;
    }
    const std::size_t feature = result.features.size();
    StraightFeature record{feature, {}, {}, segment_exponents[seed], true};
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
