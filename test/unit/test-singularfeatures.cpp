// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <tuple>
#include <utility>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>

#include "fem/singularfeatures.hpp"
#include "test-helpers.hpp"
#include "utils/communication.hpp"

namespace palace
{

namespace
{

mfem::Mesh RectangleSheetMesh(bool reverse_elements = false, bool split_attributes = true,
                              bool duplicate_face = false)
{
  mfem::Mesh mesh(3, 6, 4, duplicate_face ? 3 : 2, 3);
  mesh.AddVertex(0.0, 0.0, 0.0);
  mesh.AddVertex(1.0, 0.0, 0.0);
  mesh.AddVertex(0.0, 1.0, 0.0);
  mesh.AddVertex(1.0, 1.0, 0.0);
  mesh.AddVertex(0.5, 0.5, 1.0);
  mesh.AddVertex(0.5, 0.5, -1.0);

  const std::array<std::array<int, 4>, 4> elements{
      std::array<int, 4>{0, 1, 2, 4}, std::array<int, 4>{0, 2, 1, 5},
      std::array<int, 4>{1, 3, 2, 4}, std::array<int, 4>{1, 2, 3, 5}};
  if (reverse_elements)
  {
    for (auto element = elements.rbegin(); element != elements.rend(); ++element)
    {
      mesh.AddTet(element->data(), 1);
    }
  }
  else
  {
    for (const auto &element : elements)
    {
      mesh.AddTet(element.data(), 1);
    }
  }
  mesh.AddBdrTriangle(0, 1, 2, 7);
  mesh.AddBdrTriangle(1, 3, 2, split_attributes ? 8 : 7);
  if (duplicate_face)
  {
    mesh.AddBdrTriangle(0, 1, 2, 7);
  }
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

mfem::Mesh SubdividedTriangleSheetMesh()
{
  mfem::Mesh mesh(3, 6, 4, 2, 3);
  mesh.AddVertex(0.0, 0.0, 0.0);
  mesh.AddVertex(0.5, 0.0, 0.0);
  mesh.AddVertex(1.0, 0.0, 0.0);
  mesh.AddVertex(0.0, 1.0, 0.0);
  mesh.AddVertex(0.4, 0.3, 1.0);
  mesh.AddVertex(0.4, 0.3, -1.0);
  mesh.AddTet(0, 1, 3, 4, 1);
  mesh.AddTet(0, 3, 1, 5, 1);
  mesh.AddTet(1, 2, 3, 4, 1);
  mesh.AddTet(1, 3, 2, 5, 1);
  mesh.AddBdrTriangle(0, 1, 3, 7);
  mesh.AddBdrTriangle(1, 2, 3, 7);
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

mfem::Mesh NonmanifoldSheetEdgeMesh()
{
  mfem::Mesh mesh(3, 11, 6, 3, 3);
  mesh.AddVertex(0.0, 0.0, -1.0);
  mesh.AddVertex(0.0, 0.0, 1.0);
  constexpr double pi = 3.14159265358979323846;
  for (int branch = 0; branch < 3; branch++)
  {
    const double angle = 2.0 * pi * branch / 3.0;
    const double x = std::cos(angle);
    const double y = std::sin(angle);
    const int sheet_vertex = 2 + 3 * branch;
    const int side_0 = sheet_vertex + 1;
    const int side_1 = sheet_vertex + 2;
    mesh.AddVertex(x, y, 0.0);
    mesh.AddVertex(0.5 * x - y, 0.5 * y + x, 0.0);
    mesh.AddVertex(0.5 * x + y, 0.5 * y - x, 0.0);
    mesh.AddTet(0, 1, sheet_vertex, side_0, 1);
    mesh.AddTet(0, sheet_vertex, 1, side_1, 1);
    mesh.AddBdrTriangle(0, 1, sheet_vertex, 7);
  }
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

mfem::Mesh PerimeterJunctionMesh()
{
  mfem::Mesh mesh(3, 13, 6, 3, 3);
  mesh.AddVertex(0.0, 0.0, 0.0);
  constexpr double pi = 3.14159265358979323846;
  for (int branch = 0; branch < 3; branch++)
  {
    const double angle_0 = 2.0 * pi * branch / 3.0;
    const double angle_1 = angle_0 + pi / 3.0;
    const double x_0 = std::cos(angle_0);
    const double y_0 = std::sin(angle_0);
    const double x_1 = std::cos(angle_1);
    const double y_1 = std::sin(angle_1);
    const int face_0 = 1 + 4 * branch;
    const int face_1 = face_0 + 1;
    const int side_0 = face_0 + 2;
    const int side_1 = face_0 + 3;
    mesh.AddVertex(x_0, y_0, 0.0);
    mesh.AddVertex(x_1, y_1, 0.0);
    mesh.AddVertex(0.25 * (x_0 + x_1), 0.25 * (y_0 + y_1), 1.0);
    mesh.AddVertex(0.25 * (x_0 + x_1), 0.25 * (y_0 + y_1), -1.0);
    mesh.AddTet(0, face_0, face_1, side_0, 1);
    mesh.AddTet(0, face_1, face_0, side_1, 1);
    mesh.AddBdrTriangle(0, face_0, face_1, 7);
  }
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

std::map<std::array<int, 4>,
         std::pair<std::vector<std::pair<std::size_t, std::array<int, 4>>>,
                   std::vector<std::pair<std::size_t, std::array<int, 4>>>>>
GlobalElementIncidence(const mfem::Mesh &mesh,
                       const fem::singular::FeatureTopology &topology)
{
  using Incidence = std::pair<std::vector<std::pair<std::size_t, std::array<int, 4>>>,
                              std::vector<std::pair<std::size_t, std::array<int, 4>>>>;
  std::map<std::array<int, 4>, Incidence> result;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    const int *vertices = mesh.GetElement(element)->GetVertices();
    std::array<int, 4> element_key{vertices[0], vertices[1], vertices[2], vertices[3]};
    std::sort(element_key.begin(), element_key.end());
    auto &[nodes, edges] = result[element_key];
    for (const auto &node : topology.elements[element].nodes)
    {
      std::array<int, 4> global;
      for (int i = 0; i < 4; i++)
      {
        global[i] = vertices[node.canonical_nodes[i]];
      }
      nodes.push_back({node.vertex, global});
    }
    for (const auto &edge : topology.elements[element].edges)
    {
      std::array<int, 4> global;
      for (int i = 0; i < 4; i++)
      {
        global[i] = vertices[edge.canonical_nodes[i]];
      }
      edges.push_back({edge.feature, global});
    }
    std::sort(nodes.begin(), nodes.end());
    std::sort(edges.begin(), edges.end());
  }
  return result;
}

void CheckCanonicalElementTuples(const mfem::Mesh &mesh,
                                 const fem::singular::FeatureTopology &topology)
{
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    const int *vertices = mesh.GetElement(element)->GetVertices();
    for (const auto &node : topology.elements[element].nodes)
    {
      CAPTURE(element, node.vertex, node.mesh_vertex, node.canonical_nodes);
      CHECK(vertices[node.canonical_nodes[0]] == node.mesh_vertex);
      CHECK(vertices[node.canonical_nodes[1]] < vertices[node.canonical_nodes[2]]);
      CHECK(vertices[node.canonical_nodes[2]] < vertices[node.canonical_nodes[3]]);
    }
    for (const auto &edge : topology.elements[element].edges)
    {
      CAPTURE(element, edge.feature, edge.mesh_edge, edge.canonical_nodes);
      CHECK(vertices[edge.canonical_nodes[0]] < vertices[edge.canonical_nodes[1]]);
      CHECK(vertices[edge.canonical_nodes[2]] < vertices[edge.canonical_nodes[3]]);
      const auto &segment = topology.segments[edge.segment];
      CHECK(vertices[edge.canonical_nodes[0]] == segment.mesh_vertices[0]);
      CHECK(vertices[edge.canonical_nodes[1]] == segment.mesh_vertices[1]);
    }
  }
}

}  // namespace

TEST_CASE("Singular sheet features extract topological perimeter and element incidence",
          "[singularfeatures][Serial]")
{
  auto mesh = RectangleSheetMesh();
  const auto topology = fem::singular::ExtractSerialSheetFeatures(mesh, {8, 7}, 0.5);
  REQUIRE(topology.segments.size() == 4);
  REQUIRE(topology.vertices.size() == 4);
  REQUIRE(topology.features.size() == 4);
  REQUIRE(topology.sheet_faces.size() == 2);
  REQUIRE(topology.elements.size() == 4);
  CHECK(topology.sheet_faces[0].mesh_vertices ==
        std::array<fem::singular::GlobalVertexId, 3>{0, 1, 2});
  CHECK(topology.sheet_faces[0].boundary_attribute == 7);
  CHECK(topology.sheet_faces[1].mesh_vertices ==
        std::array<fem::singular::GlobalVertexId, 3>{1, 2, 3});
  CHECK(topology.sheet_faces[1].boundary_attribute == 8);

  const std::array<std::array<int, 2>, 4> expected_edges{
      std::array<int, 2>{0, 1}, std::array<int, 2>{0, 2}, std::array<int, 2>{1, 3},
      std::array<int, 2>{2, 3}};
  for (std::size_t i = 0; i < expected_edges.size(); i++)
  {
    CAPTURE(i);
    CHECK(topology.segments[i].mesh_vertices == expected_edges[i]);
    CHECK(topology.segments[i].feature == i);
    CHECK(topology.features[i].id == i);
    CHECK(topology.features[i].segments == std::vector<std::size_t>{i});
    CHECK(topology.features[i].nu == 0.5);
    CHECK_FALSE(topology.features[i].closed);
  }
  for (const auto &vertex : topology.vertices)
  {
    CHECK(vertex.type == fem::singular::FeatureVertexType::CORNER);
    CHECK(vertex.segments.size() == 2);
    CHECK(vertex.features.size() == 2);
    CHECK(vertex.id < topology.vertices.size());
    CHECK(vertex.nu == 0.5);
  }

  std::size_t node_incidence = 0;
  std::size_t edge_incidence = 0;
  for (const auto &element : topology.elements)
  {
    node_incidence += element.nodes.size();
    edge_incidence += element.edges.size();
    std::set<int> singular_vertices;
    for (const auto &node : element.nodes)
    {
      CHECK(singular_vertices.insert(node.mesh_vertex).second);
      CHECK(topology.vertices[node.vertex].mesh_vertex == node.mesh_vertex);
    }
    for (const auto &edge : element.edges)
    {
      const auto &segment = topology.segments[edge.segment];
      CHECK(singular_vertices.count(segment.mesh_vertices[0]) == 1);
      CHECK(singular_vertices.count(segment.mesh_vertices[1]) == 1);
    }
  }
  CHECK(node_incidence == 12);
  CHECK(edge_incidence == 8);
  CheckCanonicalElementTuples(mesh, topology);

  const auto one_triangle = fem::singular::ExtractSerialSheetFeatures(mesh, {7});
  REQUIRE(one_triangle.segments.size() == 3);
  CHECK(std::find_if(one_triangle.segments.begin(), one_triangle.segments.end(),
                     [](const auto &segment)
                     { return segment.mesh_vertices == std::array<int, 2>{1, 2}; }) !=
        one_triangle.segments.end());
}

TEST_CASE("Singular sheet features group collinear mesh edges into straight features",
          "[singularfeatures][Serial]")
{
  auto mesh = SubdividedTriangleSheetMesh();
  const auto topology = fem::singular::ExtractSerialSheetFeatures(mesh, {7}, 2.0 / 3.0);
  REQUIRE(topology.segments.size() == 4);
  REQUIRE(topology.features.size() == 3);

  const auto vertex =
      std::find_if(topology.vertices.begin(), topology.vertices.end(),
                   [](const auto &entry) { return entry.mesh_vertex == 1; });
  REQUIRE(vertex != topology.vertices.end());
  CHECK(vertex->type == fem::singular::FeatureVertexType::REGULAR);
  REQUIRE(vertex->features.size() == 1);
  const auto &straight = topology.features[vertex->features[0]];
  CHECK(straight.segments.size() == 2);
  CHECK(straight.mesh_vertices == std::vector<int>{0, 1, 2});
  CHECK(straight.nu == 2.0 / 3.0);
  CHECK_FALSE(straight.closed);
}

TEST_CASE("Singular sheet feature IDs and basis tuples ignore element ordering",
          "[singularfeatures][Serial]")
{
  auto forward_mesh = RectangleSheetMesh(false);
  auto reverse_mesh = RectangleSheetMesh(true);
  const auto forward = fem::singular::ExtractSerialSheetFeatures(forward_mesh, {7, 8});
  const auto reverse = fem::singular::ExtractSerialSheetFeatures(reverse_mesh, {7, 8});

  REQUIRE(reverse.segments.size() == forward.segments.size());
  REQUIRE(reverse.features.size() == forward.features.size());
  for (std::size_t i = 0; i < forward.segments.size(); i++)
  {
    CHECK(reverse.segments[i].mesh_vertices == forward.segments[i].mesh_vertices);
    CHECK(reverse.segments[i].feature == forward.segments[i].feature);
    CHECK(reverse.segments[i].boundary_attributes ==
          forward.segments[i].boundary_attributes);
  }
  for (std::size_t i = 0; i < forward.features.size(); i++)
  {
    CHECK(reverse.features[i].segments == forward.features[i].segments);
    CHECK(reverse.features[i].mesh_vertices == forward.features[i].mesh_vertices);
  }
  CHECK(GlobalElementIncidence(reverse_mesh, reverse) ==
        GlobalElementIncidence(forward_mesh, forward));
}

TEST_CASE("Singular sheet feature extraction rejects unsupported inputs",
          "[singularfeatures][Serial]")
{
  auto mesh = RectangleSheetMesh();
  const auto empty = fem::singular::ExtractSerialSheetFeatures(mesh, {});
  CHECK(empty.Empty());
  CHECK(empty.vertices.empty());
  CHECK(empty.segments.empty());
  CHECK(empty.sheet_faces.empty());
  CHECK(empty.elements.size() == static_cast<std::size_t>(mesh.GetNE()));
  for (const auto &element : empty.elements)
  {
    CHECK(element.nodes.empty());
    CHECK(element.edges.empty());
  }

  CHECK_THROWS_AS(fem::singular::ExtractSerialSheetFeatures(mesh, {0}),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::ExtractSerialSheetFeatures(mesh, {99}),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::ExtractSerialSheetFeatures(mesh, {7}, 0.0),
                  std::invalid_argument);

  auto external = SingleTetMesh();
  CHECK_THROWS_AS(fem::singular::ExtractSerialSheetFeatures(external, {1}),
                  std::invalid_argument);

  auto curved = RectangleSheetMesh();
  curved.SetCurvature(2);
  CHECK_THROWS_AS(fem::singular::ExtractSerialSheetFeatures(curved, {7, 8}),
                  std::invalid_argument);

  auto hexahedron =
      mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0);
  CHECK_THROWS_AS(fem::singular::ExtractSerialSheetFeatures(hexahedron, {1}),
                  std::invalid_argument);

  auto duplicate = RectangleSheetMesh(false, true, true);
  CHECK_THROWS_AS(fem::singular::ExtractSerialSheetFeatures(duplicate, {7}),
                  std::invalid_argument);

  auto nonmanifold = NonmanifoldSheetEdgeMesh();
  CHECK_THROWS_AS(fem::singular::ExtractSerialSheetFeatures(nonmanifold, {7}),
                  std::invalid_argument);

  auto junction = PerimeterJunctionMesh();
  CHECK_THROWS_AS(fem::singular::ExtractSerialSheetFeatures(junction, {7}),
                  std::invalid_argument);
}

TEST_CASE("Serial singular sheet extraction rejects a multi-rank ParMesh",
          "[singularfeatures][Parallel]")
{
  auto serial_mesh = RectangleSheetMesh();
  mfem::ParMesh parallel_mesh(Mpi::World(), serial_mesh);
  if (parallel_mesh.GetNRanks() == 1)
  {
    SUCCEED("Multi-rank rejection is exercised by the [Parallel] test run.");
    return;
  }
  CHECK_THROWS_AS(fem::singular::ExtractSerialSheetFeatures(parallel_mesh, {7, 8}),
                  std::invalid_argument);
}

TEST_CASE("Serial singular feature blueprints broadcast sparsely",
          "[singularfeatures][Parallel]")
{
  auto mesh = RectangleSheetMesh();
  const auto expected = fem::singular::ExtractSerialSheetFeatures(mesh, {7, 8});
  auto broadcast = Mpi::Root(Mpi::World()) ? expected : fem::singular::FeatureTopology{};
  fem::singular::BroadcastSerialSheetFeatures(broadcast, Mpi::World());

  CHECK(broadcast.vertices.size() == expected.vertices.size());
  CHECK(broadcast.segments.size() == expected.segments.size());
  CHECK(broadcast.features.size() == expected.features.size());
  CHECK(broadcast.sheet_faces.size() == expected.sheet_faces.size());
  CHECK(broadcast.elements.size() == expected.elements.size());
  CHECK(GlobalElementIncidence(mesh, broadcast) == GlobalElementIncidence(mesh, expected));
  for (std::size_t i = 0; i < expected.vertices.size(); i++)
  {
    CHECK(broadcast.vertices[i].id == expected.vertices[i].id);
    CHECK(broadcast.vertices[i].mesh_vertex == expected.vertices[i].mesh_vertex);
    CHECK(broadcast.vertices[i].type == expected.vertices[i].type);
    CHECK(broadcast.vertices[i].segments == expected.vertices[i].segments);
    CHECK(broadcast.vertices[i].features == expected.vertices[i].features);
    CHECK(broadcast.vertices[i].nu == expected.vertices[i].nu);
  }
  for (std::size_t i = 0; i < expected.sheet_faces.size(); i++)
  {
    CHECK(broadcast.sheet_faces[i].mesh_face == expected.sheet_faces[i].mesh_face);
    CHECK(broadcast.sheet_faces[i].mesh_vertices == expected.sheet_faces[i].mesh_vertices);
    CHECK(broadcast.sheet_faces[i].boundary_attribute ==
          expected.sheet_faces[i].boundary_attribute);
  }
}

}  // namespace palace
