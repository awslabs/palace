// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <numeric>
#include <set>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "drivers/singularsolver.hpp"
#include "fem/singulardofs.hpp"
#include "fem/singularfeatures.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

namespace palace
{

namespace
{

IoData MakeLossData(double first, double second)
{
  IoData iodata(Units(1.0, 1.0));
  auto &material_1 = iodata.domains.materials.emplace_back();
  material_1.attributes = {1};
  material_1.tandelta.s = {first, first, first};
  auto &material_2 = iodata.domains.materials.emplace_back();
  material_2.attributes = {2};
  material_2.tandelta.s = {second, second, second};
  return iodata;
}

mfem::Mesh TwoMaterialTetrahedra()
{
  mfem::Mesh mesh(3, 8, 2, 0, 3);
  mesh.AddVertex(0.0, 0.0, 0.0);
  mesh.AddVertex(1.0, 0.0, 0.0);
  mesh.AddVertex(0.0, 1.0, 0.0);
  mesh.AddVertex(0.0, 0.0, 1.0);
  mesh.AddVertex(2.0, 0.0, 0.0);
  mesh.AddVertex(3.0, 0.0, 0.0);
  mesh.AddVertex(2.0, 1.0, 0.0);
  mesh.AddVertex(2.0, 0.0, 1.0);
  mesh.AddTet(0, 1, 2, 3, 1);
  mesh.AddTet(4, 5, 6, 7, 2);
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

mfem::Mesh LongInternalLineMesh()
{
  constexpr int nx = 4;
  constexpr int ny = 2;
  mfem::Mesh mesh(2, (nx + 1) * (ny + 1), 2 * nx * ny, 1, 2);
  const auto vertex = [](int x, int y) { return y * (nx + 1) + x; };
  for (int y = 0; y <= ny; y++)
  {
    for (int x = 0; x <= nx; x++)
    {
      mesh.AddVertex(static_cast<double>(x), static_cast<double>(y));
    }
  }
  for (int y = 0; y < ny; y++)
  {
    for (int x = 0; x < nx; x++)
    {
      mesh.AddTriangle(vertex(x, y), vertex(x + 1, y), vertex(x + 1, y + 1), 1);
      mesh.AddTriangle(vertex(x, y), vertex(x + 1, y + 1), vertex(x, y + 1), 1);
    }
  }
  mesh.AddBdrSegment(vertex(1, 1), vertex(2, 1), 7);
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

mfem::Mesh InternalSheetMesh()
{
  mfem::Mesh mesh(3, 6, 4, 1, 3);
  mesh.AddVertex(0.0, 0.0, 0.0);
  mesh.AddVertex(1.0, 0.0, 0.0);
  mesh.AddVertex(0.0, 1.0, 0.0);
  mesh.AddVertex(1.0, 1.0, 0.0);
  mesh.AddVertex(0.5, 0.5, 1.0);
  mesh.AddVertex(0.5, 0.5, -1.0);
  mesh.AddTet(0, 1, 2, 4, 1);
  mesh.AddTet(0, 2, 1, 5, 1);
  mesh.AddTet(1, 3, 2, 4, 1);
  mesh.AddTet(1, 2, 3, 5, 1);
  mesh.AddBdrTriangle(0, 1, 2, 7);
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

mfem::Mesh TwoTrianglePartitionMesh()
{
  mfem::Mesh mesh(2, 4, 2, 0, 2);
  mesh.AddVertex(0.0, 0.0);
  mesh.AddVertex(1.0, 0.0);
  mesh.AddVertex(0.0, 1.0);
  mesh.AddVertex(1.0, 1.0);
  mesh.AddTriangle(0, 1, 2, 1);
  mesh.AddTriangle(1, 3, 2, 1);
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

mfem::Mesh CoincidentDisconnectedTriangles()
{
  mfem::Mesh mesh(2, 6, 2, 0, 2);
  for (int copy = 0; copy < 2; copy++)
  {
    mesh.AddVertex(0.0, 0.0);
    mesh.AddVertex(1.0, 0.0);
    mesh.AddVertex(0.0, 1.0);
    mesh.AddTriangle(3 * copy, 3 * copy + 1, 3 * copy + 2, 1);
  }
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

// Two disconnected triangles at identical coordinates but with distinct domain
// attributes. Refinement creates pairs of distinct new vertices at identical positions,
// and the differing attributes make a wrong identity assignment observable: swapping the
// identities of two coincident vertices changes the identity-to-attribute incidence,
// whereas a fully symmetric fixture would make the swap undetectable.
mfem::Mesh CoincidentDistinctAttributeTriangles()
{
  mfem::Mesh mesh(2, 6, 2, 0, 2);
  for (int copy = 0; copy < 2; copy++)
  {
    mesh.AddVertex(0.0, 0.0);
    mesh.AddVertex(1.0, 0.0);
    mesh.AddVertex(0.0, 1.0);
    mesh.AddTriangle(3 * copy, 3 * copy + 1, 3 * copy + 2, copy + 1);
  }
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

// Map each persistent vertex identity to the sorted set of domain attributes of the
// elements incident on it, gathered globally. This is invariant under partitioning for a
// fixed refinement history, and it distinguishes coincident vertices belonging to
// different attribute regions.
std::map<fem::singular::GlobalVertexId, std::set<int>>
GatherIdentityAttributes(const mfem::ParMesh &mesh,
                         const std::vector<fem::singular::GlobalVertexId> &vertex_ids)
{
  MFEM_VERIFY(vertex_ids.size() == static_cast<std::size_t>(mesh.GetNV()),
              "Identity attribute test data has an invalid local size!");
  std::vector<fem::singular::GlobalVertexId> flat;
  mfem::Array<int> element_vertices;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    mesh.GetElementVertices(element, element_vertices);
    for (int i = 0; i < element_vertices.Size(); i++)
    {
      flat.push_back(vertex_ids[element_vertices[i]]);
      flat.push_back(
          static_cast<fem::singular::GlobalVertexId>(mesh.GetAttribute(element)));
    }
  }
  const int local_count = static_cast<int>(flat.size());
  std::vector<int> counts(Mpi::Size(mesh.GetComm()));
  Mpi::Allgather(1, &local_count, counts.data(), mesh.GetComm());
  std::vector<int> offsets(counts.size(), 0);
  std::partial_sum(counts.begin(), counts.end() - 1, offsets.begin() + 1);
  const std::size_t total = std::accumulate(counts.begin(), counts.end(), std::size_t{0});
  std::vector<fem::singular::GlobalVertexId> gathered(total);
  Mpi::Allgatherv(local_count, flat.data(), gathered.data(), counts.data(), offsets.data(),
                  mesh.GetComm());

  std::map<fem::singular::GlobalVertexId, std::set<int>> result;
  for (std::size_t i = 0; i + 1 < gathered.size(); i += 2)
  {
    result[gathered[i]].insert(static_cast<int>(gathered[i + 1]));
  }
  return result;
}

fem::singular::TriangleFeatureTopology
TriangleFeature(fem::singular::FeatureVertexType type)
{
  fem::singular::TriangleFeatureTopology features;
  features.vertices.push_back(
      {0,
       0,
       {},
       0.5,
       type,
       type == fem::singular::FeatureVertexType::CORNER
           ? std::vector<fem::singular::TriangleWedgeSector>{{1, 1.0, 2.3}, {2, 2.0, 4.7}}
           : std::vector<fem::singular::TriangleWedgeSector>{}});
  return features;
}

fem::singular::FeatureTopology TetrahedralFeature(fem::singular::FeatureSegmentType type)
{
  fem::singular::FeatureTopology features;
  features.segments.push_back({0, {0, 1}, 0, {7}, type});
  features.elements.resize(2);
  for (int element = 0; element < 2; element++)
  {
    features.elements[element].edges.push_back({0, 0, 0, {0, 1, 2, 3}});
  }
  return features;
}

std::map<fem::singular::GlobalVertexId, std::array<double, 3>>
GatherVertexIdentities(const mfem::ParMesh &mesh,
                       const std::vector<fem::singular::GlobalVertexId> &vertex_ids)
{
  MFEM_VERIFY(vertex_ids.size() == static_cast<std::size_t>(mesh.GetNV()),
              "Vertex identity test data has an invalid local size!");
  const int local_vertices = mesh.GetNV();
  std::vector<int> counts(Mpi::Size(mesh.GetComm()));
  Mpi::Allgather(1, &local_vertices, counts.data(), mesh.GetComm());
  std::vector<int> offsets(counts.size());
  std::partial_sum(counts.begin(), counts.end() - 1, offsets.begin() + 1);
  const int global_vertices = std::accumulate(counts.begin(), counts.end(), 0);

  std::vector<double> local_coordinates(3 * local_vertices, 0.0);
  for (int vertex = 0; vertex < local_vertices; vertex++)
  {
    const double *coordinate = mesh.GetVertex(vertex);
    for (int d = 0; d < mesh.SpaceDimension(); d++)
    {
      local_coordinates[3 * vertex + d] = coordinate[d];
    }
  }
  std::vector<int> coordinate_counts(counts), coordinate_offsets(offsets);
  for (int &count : coordinate_counts)
  {
    count *= 3;
  }
  for (int &offset : coordinate_offsets)
  {
    offset *= 3;
  }
  std::vector<double> coordinates(3 * global_vertices);
  std::vector<fem::singular::GlobalVertexId> identities(global_vertices);
  Mpi::Allgatherv(3 * local_vertices, local_coordinates.data(), coordinates.data(),
                  coordinate_counts.data(), coordinate_offsets.data(), mesh.GetComm());
  Mpi::Allgatherv(local_vertices, vertex_ids.data(), identities.data(), counts.data(),
                  offsets.data(), mesh.GetComm());

  std::map<fem::singular::GlobalVertexId, std::array<double, 3>> result;
  for (int vertex = 0; vertex < global_vertices; vertex++)
  {
    const std::array<double, 3> coordinate{
        coordinates[3 * vertex], coordinates[3 * vertex + 1], coordinates[3 * vertex + 2]};
    const auto [record, inserted] = result.emplace(identities[vertex], coordinate);
    MFEM_VERIFY(inserted || record->second == coordinate,
                "One persistent vertex identity has inconsistent coordinates!");
  }
  return result;
}

}  // namespace

TEST_CASE("Singular loss-tangent validation distinguishes sheets from transmission wedges",
          "[singularsolver][singularelements][Serial]")
{
  constexpr double common_loss_tangent = 1.0e-3;
  auto common = MakeLossData(common_loss_tangent, common_loss_tangent);
  auto unequal = MakeLossData(1.0e-3, 2.0e-3);

  SECTION("Two-dimensional thin-sheet endpoint")
  {
    const auto features = TriangleFeature(fem::singular::FeatureVertexType::ENDPOINT);
    CHECK_NOTHROW(ValidateSingularLossTangents(unequal, features));
  }
  SECTION("Two-dimensional finite transmission wedge")
  {
    const auto features = TriangleFeature(fem::singular::FeatureVertexType::CORNER);
    CHECK_NOTHROW(ValidateSingularLossTangents(common, features));
    CHECK_THROWS_WITH(
        ValidateSingularLossTangents(unequal, features),
        Catch::Matchers::ContainsSubstring("exact transmission-wedge exponent is then "
                                           "complex"));
  }
  SECTION("Three-dimensional thin-sheet edge")
  {
    const auto mesh = TwoMaterialTetrahedra();
    const auto features = TetrahedralFeature(fem::singular::FeatureSegmentType::SHEET_EDGE);
    CHECK_NOTHROW(ValidateSingularLossTangents(unequal, mesh, features));
  }
  SECTION("Three-dimensional finite transmission wedge")
  {
    const auto mesh = TwoMaterialTetrahedra();
    const auto features =
        TetrahedralFeature(fem::singular::FeatureSegmentType::TRANSMISSION_WEDGE);
    CHECK_NOTHROW(ValidateSingularLossTangents(common, mesh, features));
    CHECK_THROWS_WITH(
        ValidateSingularLossTangents(unequal, mesh, features),
        Catch::Matchers::ContainsSubstring("exact transmission-wedge exponent is then "
                                           "complex"));
  }
}

TEST_CASE("Singular impedance validation constrains nonintegrable homogeneous traces",
          "[singularsolver][singularelements][impedance][Serial]")
{
  IoData iodata(Units(1.0, 1.0));
  auto &impedance = iodata.boundaries.impedance.emplace_back();
  impedance.attributes = {7};
  impedance.Ls = 1.0;

  SECTION("Two-dimensional finite and thin corners")
  {
    fem::singular::TriangleFeatureTopology features;
    features.selected_segments.push_back({0, 0, {0, 1}, 7});
    features.vertices.push_back(
        {0, 0, {0}, 0.6, fem::singular::FeatureVertexType::CORNER, {}});
    CHECK_NOTHROW(ValidateSingularImpedanceFeatures(iodata, features));
    CHECK(GetConstrainedSingularImpedanceAttributes(iodata, features).empty());

    features.vertices[0].nu = 0.5;
    CHECK_NOTHROW(ValidateSingularImpedanceFeatures(iodata, features));
    CHECK(GetConstrainedSingularImpedanceAttributes(iodata, features) == std::set<int>{7});
  }

  SECTION("Two-dimensional mixed PEC and impedance corner")
  {
    iodata.boundaries.pec.attributes = {8};
    fem::singular::TriangleFeatureTopology features;
    features.selected_segments.push_back({0, 0, {0, 1}, 7});
    features.selected_segments.push_back({1, 1, {0, 2}, 8});
    features.vertices.push_back(
        {0, 0, {0, 1}, 2.0 / 3.0, fem::singular::FeatureVertexType::CORNER, {}});
    CHECK_THROWS_WITH(
        ValidateSingularImpedanceFeatures(iodata, features),
        Catch::Matchers::ContainsSubstring("joins PEC and surface-impedance"));
  }

  SECTION("Three-dimensional finite and thin edges")
  {
    fem::singular::FeatureTopology features;
    features.features.push_back({0, {0}, {0, 1}, 0.6, false});
    features.segments.push_back(
        {0, {0, 1}, 0, {7}, fem::singular::FeatureSegmentType::TRANSMISSION_WEDGE});
    CHECK_NOTHROW(ValidateSingularImpedanceFeatures(iodata, features));
    CHECK(GetConstrainedSingularImpedanceAttributes(iodata, features).empty());

    features.features[0].nu = 0.5;
    CHECK_NOTHROW(ValidateSingularImpedanceFeatures(iodata, features));
    CHECK(GetConstrainedSingularImpedanceAttributes(iodata, features) == std::set<int>{7});
  }

  SECTION("Three-dimensional mixed PEC and impedance edge")
  {
    iodata.boundaries.pec.attributes = {8};
    fem::singular::FeatureTopology features;
    features.features.push_back({0, {0}, {0, 1}, 2.0 / 3.0, false});
    features.segments.push_back(
        {0, {0, 1}, 0, {7, 8}, fem::singular::FeatureSegmentType::TRANSMISSION_WEDGE});
    CHECK_THROWS_WITH(
        ValidateSingularImpedanceFeatures(iodata, features),
        Catch::Matchers::ContainsSubstring("joins PEC and surface-impedance"));
  }
}

TEST_CASE("Conforming AMR preserves singular source identities and line-tip incidence",
          "[singularsolver][singularelements][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  auto serial_mesh = LongInternalLineMesh();
  auto serial_features = fem::singular::ExtractSerialLineTipFeatures(serial_mesh, {7});
  mfem::ParMesh mesh(Mpi::World(), serial_mesh);

  std::vector<fem::singular::GlobalVertexId> source_vertex_ids(mesh.GetNV());
  std::iota(source_vertex_ids.begin(), source_vertex_ids.end(), 0);
  std::vector<fem::singular::GlobalVertexId> source_element_ids(mesh.GetNE());
  std::iota(source_element_ids.begin(), source_element_ids.end(), 0);
  NonconformingVertexIdentity vertex_identity;
  auto local_features = fem::singular::DistributeSerialLineTipFeatures(
      serial_features, mesh, source_vertex_ids, source_element_ids);
  mfem::Array<int> marked;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    if (!local_features.elements[element].nodes.empty())
    {
      marked.Append(element);
      break;
    }
  }
  REQUIRE(marked.Size() == 1);
  const int old_vertices = mesh.GetNV();
  mesh.GeneralRefinement(marked, -1, 1);
  REQUIRE_FALSE(mesh.Nonconforming());

  UpdateSingularSourceEntityIds(mesh, source_vertex_ids, source_element_ids,
                                vertex_identity);
  REQUIRE(source_vertex_ids.size() == static_cast<std::size_t>(mesh.GetNV()));
  for (int vertex = 0; vertex < old_vertices; vertex++)
  {
    CHECK(source_vertex_ids[vertex] == vertex);
  }
  CHECK(std::set(source_vertex_ids.begin(), source_vertex_ids.end()).size() ==
        source_vertex_ids.size());

  RebuildRefinedSingularFeatures(mesh, {7}, source_vertex_ids, serial_features);
  local_features = fem::singular::DistributeSerialLineTipFeatures(
      serial_features, mesh, source_vertex_ids, source_element_ids);
  CHECK_NOTHROW(BuildSingularRefinementProtection(mesh, local_features, source_vertex_ids));
  CHECK_NOTHROW(fem::singular::BuildLocalTriangleDofTopology(mesh, local_features,
                                                             source_vertex_ids, 1));
}

TEST_CASE("Conforming AMR distinguishes coincident topological vertices",
          "[singularsolver][singularelements][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  auto serial_mesh = CoincidentDisconnectedTriangles();
  mfem::ParMesh mesh(Mpi::World(), serial_mesh);
  std::vector<fem::singular::GlobalVertexId> source_vertex_ids(mesh.GetNV());
  std::iota(source_vertex_ids.begin(), source_vertex_ids.end(), 0);
  std::vector<fem::singular::GlobalVertexId> source_element_ids(mesh.GetNE());
  std::iota(source_element_ids.begin(), source_element_ids.end(), 0);
  NonconformingVertexIdentity vertex_identity;

  mfem::Array<int> marked(mesh.GetNE());
  std::iota(marked.begin(), marked.end(), 0);
  const int old_vertices = mesh.GetNV();
  mesh.GeneralRefinement(marked, -1, 1);
  REQUIRE_FALSE(mesh.Nonconforming());
  REQUIRE(mesh.GetNV() > old_vertices);

  CHECK_NOTHROW(UpdateSingularSourceEntityIds(mesh, source_vertex_ids, source_element_ids,
                                              vertex_identity));
  CHECK(std::set(source_vertex_ids.begin(), source_vertex_ids.end()).size() ==
        source_vertex_ids.size());

  std::map<std::array<double, 2>, int> coordinate_counts;
  for (int vertex = old_vertices; vertex < mesh.GetNV(); vertex++)
  {
    const double *coordinate = mesh.GetVertex(vertex);
    coordinate_counts[{coordinate[0], coordinate[1]}]++;
  }
  CHECK(std::any_of(coordinate_counts.begin(), coordinate_counts.end(),
                    [](const auto &entry) { return entry.second > 1; }));
}

TEST_CASE("Nonconforming AMR refines and rebuilds a conforming singular patch",
          "[singularsolver][singularelements][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  const int max_nc_levels = GENERATE(0, 1);
  CAPTURE(max_nc_levels);
  auto serial_mesh = LongInternalLineMesh();
  const int original_vertices = serial_mesh.GetNV();
  auto serial_features = fem::singular::ExtractSerialLineTipFeatures(serial_mesh, {7});
  serial_mesh.EnsureNCMesh(true);
  mfem::ParMesh mesh(Mpi::World(), serial_mesh);

  std::vector<fem::singular::GlobalVertexId> source_vertex_ids;
  std::vector<fem::singular::GlobalVertexId> source_element_ids;
  NonconformingVertexIdentity vertex_identity;
  UpdateSingularSourceEntityIds(mesh, source_vertex_ids, source_element_ids,
                                vertex_identity);
  for (int vertex = 0; vertex < original_vertices; vertex++)
  {
    CHECK(source_vertex_ids[vertex] == vertex);
  }
  auto local_features = fem::singular::DistributeSerialLineTipFeatures(
      serial_features, mesh, source_vertex_ids, source_element_ids);
  for (int iteration = 0; iteration < 2; iteration++)
  {
    mfem::Array<int> protection;
    bool conforming = false;
    INFO("Checking the singular patch before nonconforming refinement " << iteration + 1);
    REQUIRE_NOTHROW(protection = BuildSingularRefinementProtection(
                        mesh, local_features, source_vertex_ids, &conforming));
    REQUIRE(conforming);

    mfem::Array<int> marked;
    for (int element = 0; element < protection.Size(); element++)
    {
      if (protection[element])
      {
        marked.Append(element);
      }
    }
    REQUIRE(marked.Size() > 0);
    mesh.GeneralRefinement(marked, -1, max_nc_levels);
    REQUIRE(mesh.Nonconforming());

    UpdateSingularSourceEntityIds(mesh, source_vertex_ids, source_element_ids,
                                  vertex_identity);
    CHECK(std::set(source_vertex_ids.begin(), source_vertex_ids.end()).size() ==
          source_vertex_ids.size());
    RebuildRefinedSingularFeatures(mesh, {7}, source_vertex_ids, serial_features);
    local_features = fem::singular::DistributeSerialLineTipFeatures(
        serial_features, mesh, source_vertex_ids, source_element_ids);
    INFO("Checking the singular patch after nonconforming refinement " << iteration + 1);
    conforming = false;
    CHECK_NOTHROW(BuildSingularRefinementProtection(mesh, local_features, source_vertex_ids,
                                                    &conforming));
    CHECK(conforming);
    CHECK_NOTHROW(fem::singular::BuildLocalTriangleDofTopology(mesh, local_features,
                                                               source_vertex_ids, 1));
  }
}

TEST_CASE("Protected nonconforming AMR keeps singular traces conforming",
          "[singularsolver][singularelements][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  auto serial_mesh = LongInternalLineMesh();
  auto serial_features = fem::singular::ExtractSerialLineTipFeatures(serial_mesh, {7});
  serial_mesh.EnsureNCMesh(true);
  mfem::ParMesh mesh(Mpi::World(), serial_mesh);

  std::vector<fem::singular::GlobalVertexId> source_vertex_ids;
  std::vector<fem::singular::GlobalVertexId> source_element_ids;
  NonconformingVertexIdentity vertex_identity;
  UpdateSingularSourceEntityIds(mesh, source_vertex_ids, source_element_ids,
                                vertex_identity);
  auto local_features = fem::singular::DistributeSerialLineTipFeatures(
      serial_features, mesh, source_vertex_ids, source_element_ids);
  const auto CountEnrichedElements = [](const auto &features)
  {
    return std::count_if(features.elements.begin(), features.elements.end(),
                         [](const auto &element) { return !element.nodes.empty(); });
  };
  const int initial_enriched_elements = CountEnrichedElements(local_features);
  REQUIRE(initial_enriched_elements > 0);

  for (int iteration = 0; iteration < 3; iteration++)
  {
    bool conforming = false;
    const auto protection = BuildSingularRefinementProtection(
        mesh, local_features, source_vertex_ids, &conforming);
    REQUIRE(conforming);

    mfem::Array<int> marked;
    for (int element = 0; element < protection.Size(); element++)
    {
      if (!protection[element])
      {
        marked.Append(element);
      }
    }
    REQUIRE(marked.Size() > 0);
    mesh.GeneralRefinement(marked, -1, 0);

    UpdateSingularSourceEntityIds(mesh, source_vertex_ids, source_element_ids,
                                  vertex_identity);
    RebuildRefinedSingularFeatures(mesh, {7}, source_vertex_ids, serial_features);
    local_features = fem::singular::DistributeSerialLineTipFeatures(
        serial_features, mesh, source_vertex_ids, source_element_ids);
    CHECK(CountEnrichedElements(local_features) == initial_enriched_elements);
    CHECK_NOTHROW(BuildSingularRefinementProtection(mesh, local_features, source_vertex_ids,
                                                    &conforming));
    CHECK(conforming);
    CHECK_NOTHROW(fem::singular::BuildLocalTriangleDofTopology(mesh, local_features,
                                                               source_vertex_ids, 1));
  }
}

TEST_CASE("Nonconforming AMR detects and repairs a hanging singular interface",
          "[singularsolver][singularelements][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  const int dimension = GENERATE(2, 3);
  const int max_nc_levels = GENERATE(0, 1);
  CAPTURE(dimension, max_nc_levels);

  auto serial_mesh = dimension == 2 ? LongInternalLineMesh() : InternalSheetMesh();
  fem::singular::TriangleFeatureTopology serial_line_features;
  fem::singular::FeatureTopology serial_sheet_features;
  if (dimension == 2)
  {
    serial_line_features = fem::singular::ExtractSerialLineTipFeatures(serial_mesh, {7});
  }
  else
  {
    serial_sheet_features = fem::singular::ExtractSerialSheetFeatures(serial_mesh, {7});
  }
  serial_mesh.EnsureNCMesh(true);
  mfem::ParMesh mesh(Mpi::World(), serial_mesh);
  std::vector<fem::singular::GlobalVertexId> source_vertex_ids;
  std::vector<fem::singular::GlobalVertexId> source_element_ids;
  NonconformingVertexIdentity vertex_identity;
  UpdateSingularSourceEntityIds(mesh, source_vertex_ids, source_element_ids,
                                vertex_identity);

  if (dimension == 2)
  {
    auto local_features = fem::singular::DistributeSerialLineTipFeatures(
        serial_line_features, mesh, source_vertex_ids, source_element_ids);
    const auto initial_closure =
        BuildSingularRefinementProtection(mesh, local_features, source_vertex_ids);
    const int initial_closure_size =
        std::count(initial_closure.begin(), initial_closure.end(), 1);
    mfem::Array<int> marked;
    for (int element = 0; element < mesh.GetNE(); element++)
    {
      if (!local_features.elements[element].nodes.empty())
      {
        marked.Append(element);
        break;
      }
    }
    REQUIRE(marked.Size() == 1);
    REQUIRE(marked.Size() < initial_closure_size);
    mesh.GeneralRefinement(marked, -1, max_nc_levels);
    UpdateSingularSourceEntityIds(mesh, source_vertex_ids, source_element_ids,
                                  vertex_identity);
    RebuildRefinedSingularFeatures(mesh, {7}, source_vertex_ids, serial_line_features);
    local_features = fem::singular::DistributeSerialLineTipFeatures(
        serial_line_features, mesh, source_vertex_ids, source_element_ids);

    bool conforming = true;
    mfem::Array<int> repair_elements;
    auto closure = BuildSingularRefinementProtection(
        mesh, local_features, source_vertex_ids, &conforming, &repair_elements);
    REQUIRE_FALSE(conforming);
    for (int repair = 0; repair < 8 && !conforming; repair++)
    {
      marked.SetSize(0);
      for (int element = 0; element < repair_elements.Size(); element++)
      {
        if (repair_elements[element])
        {
          marked.Append(element);
        }
      }
      REQUIRE(marked.Size() > 0);
      mesh.GeneralRefinement(marked, -1, max_nc_levels);
      UpdateSingularSourceEntityIds(mesh, source_vertex_ids, source_element_ids,
                                    vertex_identity);
      RebuildRefinedSingularFeatures(mesh, {7}, source_vertex_ids, serial_line_features);
      local_features = fem::singular::DistributeSerialLineTipFeatures(
          serial_line_features, mesh, source_vertex_ids, source_element_ids);
      closure = BuildSingularRefinementProtection(mesh, local_features, source_vertex_ids,
                                                  &conforming, &repair_elements);
    }
    CHECK(conforming);
    CHECK_NOTHROW(fem::singular::BuildLocalTriangleDofTopology(mesh, local_features,
                                                               source_vertex_ids, 1));
  }
  else
  {
    auto local_features = fem::singular::DistributeSerialSheetFeatures(
        serial_sheet_features, mesh, source_vertex_ids, source_element_ids);
    const auto initial_closure =
        BuildSingularRefinementProtection(mesh, local_features, source_vertex_ids);
    const int initial_closure_size =
        std::count(initial_closure.begin(), initial_closure.end(), 1);
    mfem::Array<int> marked;
    for (int element = 0; element < mesh.GetNE(); element++)
    {
      if (!local_features.elements[element].nodes.empty() ||
          !local_features.elements[element].edges.empty())
      {
        marked.Append(element);
        break;
      }
    }
    REQUIRE(marked.Size() == 1);
    REQUIRE(marked.Size() < initial_closure_size);
    mesh.GeneralRefinement(marked, -1, max_nc_levels);
    UpdateSingularSourceEntityIds(mesh, source_vertex_ids, source_element_ids,
                                  vertex_identity);
    RebuildRefinedSingularFeatures(mesh, {7}, source_vertex_ids, serial_sheet_features);
    local_features = fem::singular::DistributeSerialSheetFeatures(
        serial_sheet_features, mesh, source_vertex_ids, source_element_ids);

    bool conforming = true;
    mfem::Array<int> repair_elements;
    auto closure = BuildSingularRefinementProtection(
        mesh, local_features, source_vertex_ids, &conforming, &repair_elements);
    REQUIRE_FALSE(conforming);
    for (int repair = 0; repair < 8 && !conforming; repair++)
    {
      marked.SetSize(0);
      for (int element = 0; element < repair_elements.Size(); element++)
      {
        if (repair_elements[element])
        {
          marked.Append(element);
        }
      }
      REQUIRE(marked.Size() > 0);
      mesh.GeneralRefinement(marked, -1, max_nc_levels);
      UpdateSingularSourceEntityIds(mesh, source_vertex_ids, source_element_ids,
                                    vertex_identity);
      RebuildRefinedSingularFeatures(mesh, {7}, source_vertex_ids, serial_sheet_features);
      local_features = fem::singular::DistributeSerialSheetFeatures(
          serial_sheet_features, mesh, source_vertex_ids, source_element_ids);
      closure = BuildSingularRefinementProtection(mesh, local_features, source_vertex_ids,
                                                  &conforming, &repair_elements);
    }
    CHECK(conforming);
    CHECK_NOTHROW(
        fem::singular::BuildLocalDofTopology(mesh, local_features, source_vertex_ids, 1));
  }
}

TEST_CASE("Conforming AMR rebuilds refined three-dimensional singular edges",
          "[singularsolver][singularelements][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  auto serial_mesh = InternalSheetMesh();
  auto serial_features = fem::singular::ExtractSerialSheetFeatures(serial_mesh, {7});
  mfem::ParMesh mesh(Mpi::World(), serial_mesh);
  std::vector<fem::singular::GlobalVertexId> source_vertex_ids(mesh.GetNV());
  std::iota(source_vertex_ids.begin(), source_vertex_ids.end(), 0);
  std::vector<fem::singular::GlobalVertexId> source_element_ids(mesh.GetNE());
  std::iota(source_element_ids.begin(), source_element_ids.end(), 0);
  NonconformingVertexIdentity vertex_identity;

  mfem::Array<int> marked(1);
  marked[0] = 0;
  mesh.GeneralRefinement(marked, -1, 1);
  UpdateSingularSourceEntityIds(mesh, source_vertex_ids, source_element_ids,
                                vertex_identity);
  RebuildRefinedSingularFeatures(mesh, {7}, source_vertex_ids, serial_features);
  const auto local_features = fem::singular::DistributeSerialSheetFeatures(
      serial_features, mesh, source_vertex_ids, source_element_ids);
  CHECK(serial_features.segments.size() >= 4);
  CHECK_NOTHROW(BuildSingularRefinementProtection(mesh, local_features, source_vertex_ids));
  CHECK_NOTHROW(
      fem::singular::BuildLocalDofTopology(mesh, local_features, source_vertex_ids, 1));
}

TEST_CASE("Singular refinement protection crosses MPI partition faces",
          "[singularsolver][singularelements][Parallel]")
{
  if (Mpi::Size(Mpi::World()) == 1)
  {
    SUCCEED("Cross-rank singular-patch protection is exercised by the parallel suite.");
    return;
  }
  REQUIRE(Mpi::Size(Mpi::World()) == 2);
  auto serial_mesh = TwoTrianglePartitionMesh();
  const std::array<int, 2> partition{0, 1};
  mfem::ParMesh mesh(Mpi::World(), serial_mesh, partition.data());
  const auto source_vertex_ids =
      fem::singular::MapPartitionedSerialVertexIds(serial_mesh, mesh, partition.data());
  REQUIRE(mesh.GetNE() == 1);

  fem::singular::TriangleFeatureTopology local_features;
  local_features.elements.resize(mesh.GetNE());
  if (Mpi::Rank(Mpi::World()) == 0)
  {
    local_features.elements[0].nodes.push_back({});
  }
  const auto protection =
      BuildSingularRefinementProtection(mesh, local_features, source_vertex_ids);
  REQUIRE(protection.Size() == 1);
  CHECK(protection[0] == 1);
}

TEST_CASE("Conforming refinement gives decomposition-independent vertex identities",
          "[singularsolver][singularelements][Parallel]")
{
  if (Mpi::Size(Mpi::World()) == 1)
  {
    SUCCEED("Conforming refined identities are exercised by the parallel suite.");
    return;
  }
  REQUIRE(Mpi::Size(Mpi::World()) == 2);

  auto reference_serial_mesh = TwoTrianglePartitionMesh();
  mfem::ParMesh reference_mesh(MPI_COMM_SELF, reference_serial_mesh);
  std::vector<fem::singular::GlobalVertexId> reference_vertex_ids(reference_mesh.GetNV());
  std::iota(reference_vertex_ids.begin(), reference_vertex_ids.end(), 0);
  std::vector<fem::singular::GlobalVertexId> reference_element_ids(reference_mesh.GetNE());
  std::iota(reference_element_ids.begin(), reference_element_ids.end(), 0);
  NonconformingVertexIdentity reference_vertex_identity;
  mfem::Array<int> reference_marked(reference_mesh.GetNE());
  std::iota(reference_marked.begin(), reference_marked.end(), 0);
  reference_mesh.GeneralRefinement(reference_marked, -1, 1);
  UpdateSingularSourceEntityIds(reference_mesh, reference_vertex_ids, reference_element_ids,
                                reference_vertex_identity);
  std::map<std::array<double, 2>, fem::singular::GlobalVertexId> reference_by_coordinate;
  for (int vertex = 0; vertex < reference_mesh.GetNV(); vertex++)
  {
    const double *coordinate = reference_mesh.GetVertex(vertex);
    REQUIRE(reference_by_coordinate
                .emplace(std::array<double, 2>{coordinate[0], coordinate[1]},
                         reference_vertex_ids[vertex])
                .second);
  }

  auto serial_mesh = TwoTrianglePartitionMesh();
  const std::array<int, 2> partition{0, 1};
  mfem::ParMesh mesh(Mpi::World(), serial_mesh, partition.data());
  auto source_vertex_ids =
      fem::singular::MapPartitionedSerialVertexIds(serial_mesh, mesh, partition.data());
  std::vector<fem::singular::GlobalVertexId> source_element_ids(mesh.GetNE());
  NonconformingVertexIdentity vertex_identity;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    source_element_ids[element] = mesh.GetGlobalElementNum(element);
  }
  mfem::Array<int> marked(mesh.GetNE());
  std::iota(marked.begin(), marked.end(), 0);
  mesh.GeneralRefinement(marked, -1, 1);
  UpdateSingularSourceEntityIds(mesh, source_vertex_ids, source_element_ids,
                                vertex_identity);

  for (int vertex = 0; vertex < mesh.GetNV(); vertex++)
  {
    const double *coordinate = mesh.GetVertex(vertex);
    const auto reference = reference_by_coordinate.find({coordinate[0], coordinate[1]});
    REQUIRE(reference != reference_by_coordinate.end());
    CHECK(source_vertex_ids[vertex] == reference->second);
  }
}

TEST_CASE("Nonconforming refinement gives shared vertices rank-consistent identities",
          "[singularsolver][singularelements][Parallel]")
{
  if (Mpi::Size(Mpi::World()) == 1)
  {
    SUCCEED("Shared nonconforming identities are exercised by the parallel suite.");
    return;
  }
  REQUIRE(Mpi::Size(Mpi::World()) == 2);
  auto serial_mesh = TwoTrianglePartitionMesh();
  serial_mesh.EnsureNCMesh(true);
  const std::array<int, 2> partition{0, 1};
  mfem::ParMesh mesh(Mpi::World(), serial_mesh, partition.data());
  REQUIRE(mesh.Nonconforming());
  REQUIRE(mesh.GetNE() == 1);

  std::vector<fem::singular::GlobalVertexId> source_vertex_ids;
  std::vector<fem::singular::GlobalVertexId> source_element_ids;
  NonconformingVertexIdentity vertex_identity;
  UpdateSingularSourceEntityIds(mesh, source_vertex_ids, source_element_ids,
                                vertex_identity);
  for (int iteration = 0; iteration < 4; iteration++)
  {
    mfem::Array<int> marked(1);
    marked[0] = 0;
    mesh.GeneralRefinement(marked, -1, 0);
    UpdateSingularSourceEntityIds(mesh, source_vertex_ids, source_element_ids,
                                  vertex_identity);
  }

  const int local_vertices = mesh.GetNV();
  std::vector<int> counts(Mpi::Size(Mpi::World()));
  Mpi::Allgather(1, &local_vertices, counts.data(), Mpi::World());
  std::vector<int> offsets(counts.size());
  std::partial_sum(counts.begin(), counts.end() - 1, offsets.begin() + 1);
  const int global_vertices = std::accumulate(counts.begin(), counts.end(), 0);

  std::vector<double> local_coordinates(2 * local_vertices);
  for (int vertex = 0; vertex < local_vertices; vertex++)
  {
    const auto *coordinate = mesh.GetVertex(vertex);
    local_coordinates[2 * vertex] = coordinate[0];
    local_coordinates[2 * vertex + 1] = coordinate[1];
  }
  std::vector<int> coordinate_counts(counts), coordinate_offsets(offsets);
  for (int &count : coordinate_counts)
  {
    count *= 2;
  }
  for (int &offset : coordinate_offsets)
  {
    offset *= 2;
  }
  std::vector<double> coordinates(2 * global_vertices);
  std::vector<fem::singular::GlobalVertexId> identities(global_vertices);
  Mpi::Allgatherv(2 * local_vertices, local_coordinates.data(), coordinates.data(),
                  coordinate_counts.data(), coordinate_offsets.data(), Mpi::World());
  Mpi::Allgatherv(local_vertices, source_vertex_ids.data(), identities.data(),
                  counts.data(), offsets.data(), Mpi::World());

  std::map<std::array<double, 2>, fem::singular::GlobalVertexId> identity_by_coordinate;
  std::map<fem::singular::GlobalVertexId, std::array<double, 2>> coordinate_by_identity;
  int repeated_coordinates = 0;
  for (int vertex = 0; vertex < global_vertices; vertex++)
  {
    const std::array<double, 2> coordinate{coordinates[2 * vertex],
                                           coordinates[2 * vertex + 1]};
    const auto [record, inserted] =
        identity_by_coordinate.emplace(coordinate, identities[vertex]);
    const auto [reverse, inserted_reverse] =
        coordinate_by_identity.emplace(identities[vertex], coordinate);
    CHECK((inserted_reverse || reverse->second == coordinate));
    if (!inserted)
    {
      repeated_coordinates++;
      CHECK(record->second == identities[vertex]);
    }
  }
  CHECK(repeated_coordinates >= 3);
}

TEST_CASE("Rebalancing preserves singular feature incidence and global DOFs",
          "[singularsolver][singularelements][Parallel]")
{
  if (Mpi::Size(Mpi::World()) == 1)
  {
    SUCCEED("Singular feature migration is exercised by the parallel suite.");
    return;
  }
  REQUIRE(Mpi::Size(Mpi::World()) == 2);

  const bool nonconforming = GENERATE(false, true);
  CAPTURE(nonconforming);
  auto serial_mesh = LongInternalLineMesh();
  auto serial_features = fem::singular::ExtractSerialLineTipFeatures(serial_mesh, {7});
  if (nonconforming)
  {
    serial_mesh.EnsureNCMesh(true);
  }

  std::vector<int> partition(serial_mesh.GetNE(), 0);
  partition.back() = 1;
  auto parallel_mesh =
      std::make_unique<mfem::ParMesh>(Mpi::World(), serial_mesh, partition.data());
  mesh::PartitionMetadata metadata;
  NonconformingVertexIdentity vertex_identity;
  if (nonconforming)
  {
    UpdateSingularSourceEntityIds(*parallel_mesh, metadata.source_vertex_ids,
                                  metadata.source_element_ids, vertex_identity);
  }
  else
  {
    metadata.source_vertex_ids = fem::singular::MapPartitionedSerialVertexIds(
        serial_mesh, *parallel_mesh, partition.data());
    metadata.source_element_ids.resize(parallel_mesh->GetNE());
    std::iota(metadata.source_element_ids.begin(), metadata.source_element_ids.end(), 0);
  }

  const auto CheckTopology =
      [&](const mfem::ParMesh &mesh, const mesh::PartitionMetadata &source)
  {
    const auto local_features = fem::singular::DistributeSerialLineTipFeatures(
        serial_features, mesh, source.source_vertex_ids, source.source_element_ids);
    const auto local_dofs = fem::singular::BuildLocalTriangleDofTopology(
        mesh, local_features, source.source_vertex_ids, 1);
    const auto numbering =
        fem::singular::BuildParallelDofNumbering(Mpi::World(), local_dofs);
    bool conforming = false;
    CHECK_NOTHROW(BuildSingularRefinementProtection(mesh, local_features,
                                                    source.source_vertex_ids, &conforming));
    CHECK(conforming);
    int enriched_elements = 0;
    for (const auto &element : local_features.elements)
    {
      enriched_elements += !element.nodes.empty();
    }
    Mpi::GlobalSum(1, &enriched_elements, Mpi::World());
    return std::array<HYPRE_BigInt, 3>{numbering.h1.global_size, numbering.nd.global_size,
                                       enriched_elements};
  };

  auto local_features = fem::singular::DistributeSerialLineTipFeatures(
      serial_features, *parallel_mesh, metadata.source_vertex_ids,
      metadata.source_element_ids);
  const auto closure = BuildSingularRefinementProtection(*parallel_mesh, local_features,
                                                         metadata.source_vertex_ids);
  mfem::Array<int> marked;
  for (int element = 0; element < closure.Size(); element++)
  {
    if (closure[element])
    {
      marked.Append(element);
    }
  }
  int global_marked = marked.Size();
  Mpi::GlobalSum(1, &global_marked, Mpi::World());
  REQUIRE(global_marked > 0);
  parallel_mesh->GeneralRefinement(marked, -1, 1);
  UpdateSingularSourceEntityIds(*parallel_mesh, metadata.source_vertex_ids,
                                metadata.source_element_ids, vertex_identity);
  RebuildRefinedSingularFeatures(*parallel_mesh, {7}, metadata.source_vertex_ids,
                                 serial_features);
  local_features = fem::singular::DistributeSerialLineTipFeatures(
      serial_features, *parallel_mesh, metadata.source_vertex_ids,
      metadata.source_element_ids);
  bool refined_closure_conforming = false;
  CHECK_NOTHROW(BuildSingularRefinementProtection(*parallel_mesh, local_features,
                                                  metadata.source_vertex_ids,
                                                  &refined_closure_conforming));
  CHECK(refined_closure_conforming);
  int refined_enriched_elements = 0;
  for (const auto &element : local_features.elements)
  {
    refined_enriched_elements += !element.nodes.empty();
  }
  Mpi::GlobalSum(1, &refined_enriched_elements, Mpi::World());
  const auto refined_vertex_identities =
      GatherVertexIdentities(*parallel_mesh, metadata.source_vertex_ids);

  IoData iodata(Units(1.0, 1.0));
  iodata.model.refinement.maximum_imbalance = 1.01;
  const double initial_imbalance = mesh::RebalanceMesh(iodata, parallel_mesh, &metadata);
  CHECK(initial_imbalance > iodata.model.refinement.maximum_imbalance);
  CHECK(parallel_mesh->GetNE() > 0);
  if (nonconforming)
  {
    CHECK_NOTHROW(vertex_identity.Observe(*parallel_mesh, metadata.source_vertex_ids));
  }
  CHECK(GatherVertexIdentities(*parallel_mesh, metadata.source_vertex_ids) ==
        refined_vertex_identities);
  const auto rebalanced_topology = CheckTopology(*parallel_mesh, metadata);
  CHECK(rebalanced_topology[0] > 0);
  CHECK(rebalanced_topology[1] >= rebalanced_topology[0]);
  CHECK(rebalanced_topology[2] == refined_enriched_elements);

  if (nonconforming)
  {
    local_features = fem::singular::DistributeSerialLineTipFeatures(
        serial_features, *parallel_mesh, metadata.source_vertex_ids,
        metadata.source_element_ids);
    const auto rebalanced_closure = BuildSingularRefinementProtection(
        *parallel_mesh, local_features, metadata.source_vertex_ids);
    marked.DeleteAll();
    for (int element = 0; element < rebalanced_closure.Size(); element++)
    {
      if (rebalanced_closure[element])
      {
        marked.Append(element);
      }
    }
    global_marked = marked.Size();
    Mpi::GlobalSum(1, &global_marked, Mpi::World());
    REQUIRE(global_marked > 0);

    parallel_mesh->GeneralRefinement(marked, -1, 1);
    CHECK_NOTHROW(UpdateSingularSourceEntityIds(*parallel_mesh, metadata.source_vertex_ids,
                                                metadata.source_element_ids,
                                                vertex_identity));
    RebuildRefinedSingularFeatures(*parallel_mesh, {7}, metadata.source_vertex_ids,
                                   serial_features);
    const auto twice_refined_topology = CheckTopology(*parallel_mesh, metadata);
    CHECK(twice_refined_topology[0] > 0);
    CHECK(twice_refined_topology[1] >= twice_refined_topology[0]);
    CHECK(twice_refined_topology[2] >= refined_enriched_elements);
  }
}

TEST_CASE("Three-dimensional singular refinement survives MPI rebalancing",
          "[singularsolver][singularelements][Parallel]")
{
  if (Mpi::Size(Mpi::World()) == 1)
  {
    SUCCEED("Three-dimensional singular AMR migration is exercised in parallel.");
    return;
  }
  REQUIRE(Mpi::Size(Mpi::World()) == 2);

  const bool nonconforming = GENERATE(false, true);
  CAPTURE(nonconforming);
  auto serial_mesh = InternalSheetMesh();
  auto serial_features = fem::singular::ExtractSerialSheetFeatures(serial_mesh, {7});
  if (nonconforming)
  {
    serial_mesh.EnsureNCMesh(true);
  }

  std::vector<int> partition(serial_mesh.GetNE(), 0);
  partition.back() = 1;
  auto parallel_mesh =
      std::make_unique<mfem::ParMesh>(Mpi::World(), serial_mesh, partition.data());
  mesh::PartitionMetadata metadata;
  NonconformingVertexIdentity vertex_identity;
  if (nonconforming)
  {
    UpdateSingularSourceEntityIds(*parallel_mesh, metadata.source_vertex_ids,
                                  metadata.source_element_ids, vertex_identity);
  }
  else
  {
    metadata.source_vertex_ids = fem::singular::MapPartitionedSerialVertexIds(
        serial_mesh, *parallel_mesh, partition.data());
    metadata.source_element_ids.resize(parallel_mesh->GetNE());
    std::iota(metadata.source_element_ids.begin(), metadata.source_element_ids.end(), 0);
  }

  auto local_features = fem::singular::DistributeSerialSheetFeatures(
      serial_features, *parallel_mesh, metadata.source_vertex_ids,
      metadata.source_element_ids);
  const auto closure = BuildSingularRefinementProtection(*parallel_mesh, local_features,
                                                         metadata.source_vertex_ids);
  mfem::Array<int> marked;
  for (int element = 0; element < closure.Size(); element++)
  {
    if (closure[element])
    {
      marked.Append(element);
    }
  }
  int global_marked = marked.Size();
  Mpi::GlobalSum(1, &global_marked, Mpi::World());
  REQUIRE(global_marked > 0);
  parallel_mesh->GeneralRefinement(marked, -1, 1);
  UpdateSingularSourceEntityIds(*parallel_mesh, metadata.source_vertex_ids,
                                metadata.source_element_ids, vertex_identity);
  RebuildRefinedSingularFeatures(*parallel_mesh, {7}, metadata.source_vertex_ids,
                                 serial_features);
  local_features = fem::singular::DistributeSerialSheetFeatures(
      serial_features, *parallel_mesh, metadata.source_vertex_ids,
      metadata.source_element_ids);
  bool refined_closure_conforming = false;
  CHECK_NOTHROW(BuildSingularRefinementProtection(*parallel_mesh, local_features,
                                                  metadata.source_vertex_ids,
                                                  &refined_closure_conforming));
  CHECK(refined_closure_conforming);
  int refined_enriched_elements = 0;
  for (const auto &element : local_features.elements)
  {
    refined_enriched_elements += !element.nodes.empty() || !element.edges.empty();
  }
  Mpi::GlobalSum(1, &refined_enriched_elements, Mpi::World());
  const auto refined_vertex_identities =
      GatherVertexIdentities(*parallel_mesh, metadata.source_vertex_ids);

  const auto CheckTopology =
      [&](const mfem::ParMesh &mesh, const mesh::PartitionMetadata &source)
  {
    const auto features = fem::singular::DistributeSerialSheetFeatures(
        serial_features, mesh, source.source_vertex_ids, source.source_element_ids);
    const auto dofs =
        fem::singular::BuildLocalDofTopology(mesh, features, source.source_vertex_ids, 1);
    const auto numbering = fem::singular::BuildParallelDofNumbering(Mpi::World(), dofs);
    bool conforming = false;
    CHECK_NOTHROW(BuildSingularRefinementProtection(mesh, features,
                                                    source.source_vertex_ids, &conforming));
    CHECK(conforming);
    int enriched_elements = 0;
    for (const auto &element : features.elements)
    {
      enriched_elements += !element.nodes.empty() || !element.edges.empty();
    }
    Mpi::GlobalSum(1, &enriched_elements, Mpi::World());
    return std::array<HYPRE_BigInt, 3>{numbering.h1.global_size, numbering.nd.global_size,
                                       enriched_elements};
  };
  IoData iodata(Units(1.0, 1.0));
  iodata.model.refinement.maximum_imbalance = 1.01;
  const double initial_imbalance = mesh::RebalanceMesh(iodata, parallel_mesh, &metadata);
  CHECK(initial_imbalance > iodata.model.refinement.maximum_imbalance);
  CHECK(parallel_mesh->GetNE() > 0);
  if (nonconforming)
  {
    CHECK_NOTHROW(vertex_identity.Observe(*parallel_mesh, metadata.source_vertex_ids));
  }
  CHECK(GatherVertexIdentities(*parallel_mesh, metadata.source_vertex_ids) ==
        refined_vertex_identities);
  const auto rebalanced_topology = CheckTopology(*parallel_mesh, metadata);
  CHECK(rebalanced_topology[0] > 0);
  CHECK(rebalanced_topology[1] >= rebalanced_topology[0]);
  CHECK(rebalanced_topology[2] == refined_enriched_elements);

  if (nonconforming)
  {
    local_features = fem::singular::DistributeSerialSheetFeatures(
        serial_features, *parallel_mesh, metadata.source_vertex_ids,
        metadata.source_element_ids);
    const auto rebalanced_closure = BuildSingularRefinementProtection(
        *parallel_mesh, local_features, metadata.source_vertex_ids);
    marked.DeleteAll();
    for (int element = 0; element < rebalanced_closure.Size(); element++)
    {
      if (rebalanced_closure[element])
      {
        marked.Append(element);
      }
    }
    global_marked = marked.Size();
    Mpi::GlobalSum(1, &global_marked, Mpi::World());
    REQUIRE(global_marked > 0);

    parallel_mesh->GeneralRefinement(marked, -1, 1);
    CHECK_NOTHROW(UpdateSingularSourceEntityIds(*parallel_mesh, metadata.source_vertex_ids,
                                                metadata.source_element_ids,
                                                vertex_identity));
    RebuildRefinedSingularFeatures(*parallel_mesh, {7}, metadata.source_vertex_ids,
                                   serial_features);
    const auto twice_refined_topology = CheckTopology(*parallel_mesh, metadata);
    CHECK(twice_refined_topology[0] > 0);
    CHECK(twice_refined_topology[1] >= twice_refined_topology[0]);
    CHECK(twice_refined_topology[2] >= refined_enriched_elements);
  }
}

TEST_CASE("Conforming refinement ancestry keys distinct coincident vertices",
          "[singularsolver][singularelements][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  auto serial_mesh = CoincidentDisconnectedTriangles();
  mfem::ParMesh mesh(Mpi::World(), serial_mesh);
  std::vector<fem::singular::GlobalVertexId> source_vertex_ids(mesh.GetNV());
  std::iota(source_vertex_ids.begin(), source_vertex_ids.end(), 0);
  std::vector<fem::singular::GlobalVertexId> source_element_ids(mesh.GetNE());
  std::iota(source_element_ids.begin(), source_element_ids.end(), 0);
  NonconformingVertexIdentity vertex_identity;
  ConformingVertexAncestry conforming_ancestry;

  // Snapshot the coarse edge ancestry before refining, exactly as the AMR loop does.
  conforming_ancestry.Observe(mesh, source_vertex_ids);
  mfem::Array<int> marked(mesh.GetNE());
  std::iota(marked.begin(), marked.end(), 0);
  const int old_vertices = mesh.GetNV();
  mesh.GeneralRefinement(marked, -1, 1);
  REQUIRE_FALSE(mesh.Nonconforming());
  REQUIRE(mesh.GetNV() > old_vertices);

  REQUIRE_NOTHROW(UpdateSingularSourceEntityIds(mesh, source_vertex_ids, source_element_ids,
                                                vertex_identity, &conforming_ancestry));
  REQUIRE(source_vertex_ids.size() == static_cast<std::size_t>(mesh.GetNV()));

  // Retained vertices keep their identities and every identity stays unique even though
  // the two disconnected triangles are geometrically coincident.
  for (int vertex = 0; vertex < old_vertices; vertex++)
  {
    CHECK(source_vertex_ids[vertex] == vertex);
  }
  CHECK(std::set(source_vertex_ids.begin(), source_vertex_ids.end()).size() ==
        source_vertex_ids.size());

  // The fixture really does produce coincident new vertices, so the uniqueness check
  // above is not vacuous.
  std::map<std::array<double, 2>, int> coordinate_counts;
  for (int vertex = old_vertices; vertex < mesh.GetNV(); vertex++)
  {
    const double *coordinate = mesh.GetVertex(vertex);
    coordinate_counts[{coordinate[0], coordinate[1]}]++;
  }
  CHECK(std::any_of(coordinate_counts.begin(), coordinate_counts.end(),
                    [](const auto &entry) { return entry.second > 1; }));
}

TEST_CASE("Conforming refinement ancestry is decomposition independent",
          "[singularsolver][singularelements][Parallel]")
{
  // The same refinement history must produce the same identities regardless of how the
  // mesh is partitioned. Compare a one-rank reference against the current rank count.
  // The coincident case is the discriminating one: distinct new vertices share
  // coordinates, so a coordinate-ordered assignment must break ties using transient
  // global numbering, which changes with the partition.
  const int fixture = GENERATE(0, 1, 2);
  CAPTURE(fixture);
  auto serial_mesh = fixture == 0   ? LongInternalLineMesh()
                     : fixture == 1 ? InternalSheetMesh()
                                    : CoincidentDistinctAttributeTriangles();
  const int dimension = serial_mesh.Dimension();
  CAPTURE(dimension);

  // Seed identities from the serial mesh numbering so the starting point is itself
  // partition independent; GetGlobalVertexIndices is a transient parallel numbering and
  // would differ between the two decompositions being compared.
  std::vector<int> partition(serial_mesh.GetNE(), 0);
  for (int element = 0; element < serial_mesh.GetNE(); element++)
  {
    partition[element] = element % Mpi::Size(Mpi::World());
  }
  mfem::ParMesh mesh(Mpi::World(), serial_mesh, partition.data());
  std::vector<fem::singular::GlobalVertexId> source_vertex_ids =
      fem::singular::MapPartitionedSerialVertexIds(serial_mesh, mesh, partition.data());
  std::vector<fem::singular::GlobalVertexId> source_element_ids(mesh.GetNE());
  std::iota(source_element_ids.begin(), source_element_ids.end(), 0);
  NonconformingVertexIdentity vertex_identity;
  ConformingVertexAncestry conforming_ancestry;

  // Refine every element so the refinement history is partition independent.
  conforming_ancestry.Observe(mesh, source_vertex_ids);
  mfem::Array<int> marked(mesh.GetNE());
  std::iota(marked.begin(), marked.end(), 0);
  mesh.GeneralRefinement(marked, -1, 1);
  REQUIRE_FALSE(mesh.Nonconforming());
  REQUIRE_NOTHROW(UpdateSingularSourceEntityIds(mesh, source_vertex_ids, source_element_ids,
                                                vertex_identity, &conforming_ancestry));

  // A vertex identity must be a function of the edge ancestry alone, so the set of
  // identities carried by each element's vertices is a partition-independent fingerprint.
  // Collect the globally sorted multiset of per-element identity tuples.
  std::vector<std::vector<fem::singular::GlobalVertexId>> element_keys;
  mfem::Array<int> element_vertices;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    mesh.GetElementVertices(element, element_vertices);
    std::vector<fem::singular::GlobalVertexId> key;
    key.reserve(element_vertices.Size());
    for (int i = 0; i < element_vertices.Size(); i++)
    {
      key.push_back(source_vertex_ids[element_vertices[i]]);
    }
    std::sort(key.begin(), key.end());
    element_keys.push_back(std::move(key));
  }

  // Flatten, gather, and compare against the same quantity computed serially.
  const int nodes_per_element = dimension + 1;
  std::vector<fem::singular::GlobalVertexId> flat;
  for (const auto &key : element_keys)
  {
    REQUIRE(static_cast<int>(key.size()) == nodes_per_element);
    flat.insert(flat.end(), key.begin(), key.end());
  }
  const int local_count = static_cast<int>(flat.size());
  std::vector<int> counts(Mpi::Size(Mpi::World()));
  Mpi::Allgather(1, &local_count, counts.data(), Mpi::World());
  std::vector<int> offsets(counts.size(), 0);
  std::partial_sum(counts.begin(), counts.end() - 1, offsets.begin() + 1);
  const std::size_t total = std::accumulate(counts.begin(), counts.end(), std::size_t{0});
  std::vector<fem::singular::GlobalVertexId> gathered(total);
  Mpi::Allgatherv(local_count, flat.data(), gathered.data(), counts.data(), offsets.data(),
                  Mpi::World());
  std::vector<std::vector<fem::singular::GlobalVertexId>> parallel_keys;
  for (std::size_t i = 0; i + nodes_per_element <= gathered.size(); i += nodes_per_element)
  {
    parallel_keys.emplace_back(gathered.begin() + i,
                               gathered.begin() + i + nodes_per_element);
  }
  std::sort(parallel_keys.begin(), parallel_keys.end());

  // Serial reference on a private communicator so the comparison is independent.
  MPI_Comm self = MPI_COMM_SELF;
  std::vector<int> reference_partition(serial_mesh.GetNE(), 0);
  mfem::ParMesh reference_mesh(self, serial_mesh, reference_partition.data());
  std::vector<fem::singular::GlobalVertexId> reference_vertex_ids =
      fem::singular::MapPartitionedSerialVertexIds(serial_mesh, reference_mesh,
                                                   reference_partition.data());
  std::vector<fem::singular::GlobalVertexId> reference_element_ids(reference_mesh.GetNE());
  std::iota(reference_element_ids.begin(), reference_element_ids.end(), 0);
  NonconformingVertexIdentity reference_identity;
  ConformingVertexAncestry reference_ancestry;
  reference_ancestry.Observe(reference_mesh, reference_vertex_ids);
  mfem::Array<int> reference_marked(reference_mesh.GetNE());
  std::iota(reference_marked.begin(), reference_marked.end(), 0);
  reference_mesh.GeneralRefinement(reference_marked, -1, 1);
  REQUIRE_NOTHROW(UpdateSingularSourceEntityIds(reference_mesh, reference_vertex_ids,
                                                reference_element_ids, reference_identity,
                                                &reference_ancestry));
  std::vector<std::vector<fem::singular::GlobalVertexId>> reference_keys;
  for (int element = 0; element < reference_mesh.GetNE(); element++)
  {
    reference_mesh.GetElementVertices(element, element_vertices);
    std::vector<fem::singular::GlobalVertexId> key;
    for (int i = 0; i < element_vertices.Size(); i++)
    {
      key.push_back(reference_vertex_ids[element_vertices[i]]);
    }
    std::sort(key.begin(), key.end());
    reference_keys.push_back(std::move(key));
  }
  std::sort(reference_keys.begin(), reference_keys.end());

  // Guard against a vacuous comparison: both sequences must be nonempty, must cover the
  // whole refined mesh, and must actually contain identities created by refinement.
  REQUIRE(!reference_keys.empty());
  REQUIRE(parallel_keys.size() == reference_keys.size());
  REQUIRE(reference_keys.size() == static_cast<std::size_t>(reference_mesh.GetNE()));
  const auto coarse_vertices = static_cast<fem::singular::GlobalVertexId>(
      mfem::ParMesh(MPI_COMM_SELF, serial_mesh).GetNV());
  CHECK(std::any_of(reference_keys.begin(), reference_keys.end(),
                    [coarse_vertices](const auto &key)
                    {
                      return std::any_of(key.begin(), key.end(), [coarse_vertices](auto id)
                                         { return id >= coarse_vertices; });
                    }));

  CHECK(parallel_keys == reference_keys);
}

TEST_CASE("Conforming refinement ancestry separates coincident vertices across ranks",
          "[singularsolver][singularelements][Parallel]")
{
  // The decisive case for ancestry keys. Two disconnected but geometrically coincident
  // triangles are split across two ranks, so refinement creates pairs of distinct new
  // vertices at identical coordinates. A coordinate-ordered scheme cannot separate them
  // without falling back on transient global numbering, which is decomposition
  // dependent; parent-edge ancestry separates them exactly.
  if (Mpi::Size(Mpi::World()) != 2)
  {
    return;
  }
  auto serial_mesh = CoincidentDistinctAttributeTriangles();
  REQUIRE(serial_mesh.GetNE() == 2);
  std::vector<int> partition{0, 1};
  mfem::ParMesh mesh(Mpi::World(), serial_mesh, partition.data());
  auto source_vertex_ids =
      fem::singular::MapPartitionedSerialVertexIds(serial_mesh, mesh, partition.data());
  std::vector<fem::singular::GlobalVertexId> source_element_ids(mesh.GetNE());
  std::iota(source_element_ids.begin(), source_element_ids.end(), 0);
  NonconformingVertexIdentity vertex_identity;
  ConformingVertexAncestry conforming_ancestry;

  conforming_ancestry.Observe(mesh, source_vertex_ids);
  mfem::Array<int> marked(mesh.GetNE());
  std::iota(marked.begin(), marked.end(), 0);
  mesh.GeneralRefinement(marked, -1, 1);
  REQUIRE_FALSE(mesh.Nonconforming());
  REQUIRE_NOTHROW(UpdateSingularSourceEntityIds(mesh, source_vertex_ids, source_element_ids,
                                                vertex_identity, &conforming_ancestry));

  // Every identity must map to exactly one coordinate. GatherVertexIdentities throws if
  // one identity carries two different coordinates.
  std::map<fem::singular::GlobalVertexId, std::array<double, 3>> identities;
  REQUIRE_NOTHROW(identities = GatherVertexIdentities(mesh, source_vertex_ids));

  // The fixture really does produce coincident vertices globally: fewer distinct
  // coordinates than distinct identities.
  std::set<std::array<double, 3>> distinct_coordinates;
  for (const auto &[id, coordinate] : identities)
  {
    MFEM_CONTRACT_VAR(id);
    distinct_coordinates.insert(coordinate);
  }
  REQUIRE(distinct_coordinates.size() < identities.size());

  // The discriminating check. Each triangle sits in its own domain attribute, so every
  // vertex identity must belong to exactly one attribute. Confusing two coincident
  // vertices merges their incidence and produces an identity spanning both attributes.
  const auto attributes = GatherIdentityAttributes(mesh, source_vertex_ids);
  REQUIRE(attributes.size() == identities.size());
  for (const auto &[id, incident] : attributes)
  {
    CAPTURE(id);
    CHECK(incident.size() == 1);
  }
}

TEST_CASE("Conforming refinement ancestry composes across repeated refinement",
          "[singularsolver][singularelements][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  const int dimension = GENERATE(2, 3);
  CAPTURE(dimension);
  auto serial_mesh = dimension == 2 ? LongInternalLineMesh() : InternalSheetMesh();
  mfem::ParMesh mesh(Mpi::World(), serial_mesh);
  std::vector<fem::singular::GlobalVertexId> source_vertex_ids(mesh.GetNV());
  std::iota(source_vertex_ids.begin(), source_vertex_ids.end(), 0);
  std::vector<fem::singular::GlobalVertexId> source_element_ids(mesh.GetNE());
  std::iota(source_element_ids.begin(), source_element_ids.end(), 0);
  NonconformingVertexIdentity vertex_identity;
  ConformingVertexAncestry conforming_ancestry;

  // Identities already assigned must never change under later unrelated refinement.
  std::vector<fem::singular::GlobalVertexId> previous;
  for (int pass = 0; pass < 3; pass++)
  {
    previous = source_vertex_ids;
    conforming_ancestry.Observe(mesh, source_vertex_ids);
    mfem::Array<int> marked(mesh.GetNE());
    std::iota(marked.begin(), marked.end(), 0);
    mesh.GeneralRefinement(marked, -1, 1);
    REQUIRE_FALSE(mesh.Nonconforming());
    REQUIRE_NOTHROW(UpdateSingularSourceEntityIds(mesh, source_vertex_ids,
                                                  source_element_ids, vertex_identity,
                                                  &conforming_ancestry));
    REQUIRE(source_vertex_ids.size() == static_cast<std::size_t>(mesh.GetNV()));
    for (std::size_t vertex = 0; vertex < previous.size(); vertex++)
    {
      CHECK(source_vertex_ids[vertex] == previous[vertex]);
    }
    CHECK(std::set(source_vertex_ids.begin(), source_vertex_ids.end()).size() ==
          source_vertex_ids.size());
  }
}

TEST_CASE("Localized singular repair converges across an MPI partition face",
          "[singularsolver][singularelements][Parallel]")
{
  // Phase 5 of the singular AMR plan. Refine one enriched element on one rank so the
  // resulting hanging singular face straddles a partition boundary, then drive the
  // detect/refine/rebuild loop to convergence. The coarse master is deliberately placed
  // on each rank in turn, because MFEM is not guaranteed to visit a master face on the
  // rank owning the coarse element.
  if (Mpi::Size(Mpi::World()) != 2)
  {
    return;
  }
  const int dimension = GENERATE(2, 3);
  const int max_nc_levels = GENERATE(0, 1);
  // Which rank owns the element that gets refined first.
  const int refining_rank = GENERATE(0, 1);
  CAPTURE(dimension, max_nc_levels, refining_rank);

  auto serial_mesh = dimension == 2 ? LongInternalLineMesh() : InternalSheetMesh();
  fem::singular::TriangleFeatureTopology serial_line_features;
  fem::singular::FeatureTopology serial_sheet_features;
  if (dimension == 2)
  {
    serial_line_features = fem::singular::ExtractSerialLineTipFeatures(serial_mesh, {7});
  }
  else
  {
    serial_sheet_features = fem::singular::ExtractSerialSheetFeatures(serial_mesh, {7});
  }
  serial_mesh.EnsureNCMesh(true);

  // Split the mesh so the enriched region spans both ranks.
  std::vector<int> partition(serial_mesh.GetNE(), 0);
  for (int element = 0; element < serial_mesh.GetNE(); element++)
  {
    partition[element] = element % 2;
  }
  mfem::ParMesh mesh(Mpi::World(), serial_mesh, partition.data());
  std::vector<fem::singular::GlobalVertexId> source_vertex_ids, source_element_ids;
  NonconformingVertexIdentity vertex_identity;
  UpdateSingularSourceEntityIds(mesh, source_vertex_ids, source_element_ids,
                                vertex_identity);

  // Rebuild the serial feature blueprint from the current mesh. Collective.
  const auto rebuild = [&]()
  {
    if (dimension == 2)
    {
      RebuildRefinedSingularFeatures(mesh, {7}, source_vertex_ids, serial_line_features);
    }
    else
    {
      RebuildRefinedSingularFeatures(mesh, {7}, source_vertex_ids, serial_sheet_features);
    }
  };

  // Detect conformity and collect owner-local repair markers. Collective.
  const auto detect = [&](bool &conforming, mfem::Array<int> &repair)
  {
    if (dimension == 2)
    {
      const auto features = fem::singular::DistributeSerialLineTipFeatures(
          serial_line_features, mesh, source_vertex_ids, source_element_ids);
      return BuildSingularRefinementProtection(mesh, features, source_vertex_ids,
                                               &conforming, &repair);
    }
    const auto features = fem::singular::DistributeSerialSheetFeatures(
        serial_sheet_features, mesh, source_vertex_ids, source_element_ids);
    return BuildSingularRefinementProtection(mesh, features, source_vertex_ids, &conforming,
                                             &repair);
  };

  // Build the enriched marker collectively on every rank: the Distribute* calls are
  // collective, so they must not be made inside a rank-conditional branch.
  const auto enriched_marker = [&]()
  {
    std::vector<bool> enriched(mesh.GetNE(), false);
    if (dimension == 2)
    {
      const auto features = fem::singular::DistributeSerialLineTipFeatures(
          serial_line_features, mesh, source_vertex_ids, source_element_ids);
      for (int element = 0; element < mesh.GetNE(); element++)
      {
        enriched[element] = !features.elements[element].nodes.empty();
      }
      return enriched;
    }
    const auto features = fem::singular::DistributeSerialSheetFeatures(
        serial_sheet_features, mesh, source_vertex_ids, source_element_ids);
    for (int element = 0; element < mesh.GetNE(); element++)
    {
      enriched[element] = !features.elements[element].nodes.empty() ||
                          !features.elements[element].edges.empty();
    }
    return enriched;
  }();

  // Refine exactly one enriched element, owned by refining_rank, to create the hanging
  // singular interface. Only the mark selection is rank-local; no collective runs here.
  mfem::Array<int> marked;
  if (Mpi::Rank(Mpi::World()) == refining_rank)
  {
    for (int element = 0; element < mesh.GetNE(); element++)
    {
      if (enriched_marker[element])
      {
        marked.Append(element);
        break;
      }
    }
  }
  int global_marked = marked.Size();
  Mpi::GlobalSum(1, &global_marked, Mpi::World());
  REQUIRE(global_marked == 1);

  mesh.GeneralRefinement(marked, -1, max_nc_levels);
  UpdateSingularSourceEntityIds(mesh, source_vertex_ids, source_element_ids,
                                vertex_identity);
  rebuild();

  bool conforming = true;
  mfem::Array<int> repair;
  detect(conforming, repair);

  // Drive the repair loop. Each pass must make progress: either the subcomplex becomes
  // conforming or the global repair-candidate count is positive.
  constexpr int repair_pass_limit = 8;
  int passes = 0;
  long long total_repair_marked = 0;
  while (!conforming && passes < repair_pass_limit)
  {
    marked.SetSize(0);
    for (int element = 0; element < repair.Size(); element++)
    {
      if (repair[element])
      {
        marked.Append(element);
      }
    }
    int global_repair = marked.Size();
    Mpi::GlobalSum(1, &global_repair, Mpi::World());
    // Fail closed: violations remain but nothing is repairable anywhere.
    REQUIRE(global_repair > 0);
    total_repair_marked += global_repair;

    mesh.GeneralRefinement(marked, -1, max_nc_levels);
    UpdateSingularSourceEntityIds(mesh, source_vertex_ids, source_element_ids,
                                  vertex_identity);
    rebuild();
    detect(conforming, repair);
    passes++;
  }

  CHECK(conforming);
  CAPTURE(passes, total_repair_marked);
  CHECK(passes < repair_pass_limit);

  // The repaired mesh must still support a valid enriched DOF topology on both ranks.
  if (dimension == 2)
  {
    const auto features = fem::singular::DistributeSerialLineTipFeatures(
        serial_line_features, mesh, source_vertex_ids, source_element_ids);
    CHECK_NOTHROW(
        fem::singular::BuildLocalTriangleDofTopology(mesh, features, source_vertex_ids, 1));
  }
  else
  {
    const auto features = fem::singular::DistributeSerialSheetFeatures(
        serial_sheet_features, mesh, source_vertex_ids, source_element_ids);
    CHECK_NOTHROW(
        fem::singular::BuildLocalDofTopology(mesh, features, source_vertex_ids, 1, {}));
  }
}

TEST_CASE("Tetrahedral singular face mask exempts only zero-trace faces",
          "[singularsolver][singularelements][Serial]")
{
  // Phase 7 classifier. The companion test in test-singularelements.cpp proves, by direct
  // basis evaluation, WHICH faces carry zero trace; this pins the mask the refinement code
  // derives from feature incidence to that same rule, over all vertex permutations.
  std::array<int, 4> permutation{0, 1, 2, 3};
  int permutations = 0;
  do
  {
    CAPTURE(permutation);
    permutations++;

    // Node feature: only the face opposite the singular node is inactive.
    {
      fem::singular::ElementFeatureIncidence incidence;
      incidence.nodes.push_back({0, 0, permutation});
      const int mask = GetTetrahedronSingularFaceMask(incidence);
      const int singular = permutation[0];
      for (int face = 0; face < 4; face++)
      {
        CAPTURE(face, singular, mask);
        CHECK(static_cast<bool>(mask & (1 << face)) == (face != singular));
      }
    }

    // Edge feature: only the two faces containing both endpoints are active.
    {
      fem::singular::ElementFeatureIncidence incidence;
      incidence.edges.push_back({0, 0, 0, permutation});
      const int mask = GetTetrahedronSingularFaceMask(incidence);
      const int first = permutation[0], second = permutation[1];
      int active = 0;
      for (int face = 0; face < 4; face++)
      {
        const bool expected = (face != first) && (face != second);
        CAPTURE(face, first, second, mask);
        CHECK(static_cast<bool>(mask & (1 << face)) == expected);
        active += expected ? 1 : 0;
      }
      CHECK(active == 2);
    }

    // Multiple features reduce with a logical OR: a face is exempt only if every feature
    // leaves it inactive. A node feature plus an edge feature sharing no exempt face must
    // therefore activate everything those two would individually.
    {
      fem::singular::ElementFeatureIncidence combined;
      combined.nodes.push_back({0, 0, permutation});
      combined.edges.push_back({0, 0, 0, permutation});
      fem::singular::ElementFeatureIncidence node_only, edge_only;
      node_only.nodes = combined.nodes;
      edge_only.edges = combined.edges;
      const int mask = GetTetrahedronSingularFaceMask(combined);
      const int expected = GetTetrahedronSingularFaceMask(node_only) |
                           GetTetrahedronSingularFaceMask(edge_only);
      CAPTURE(mask, expected);
      CHECK(mask == expected);
    }
  } while (std::next_permutation(permutation.begin(), permutation.end()));
  CHECK(permutations == 24);

  // An element with no singular features has no active faces at all.
  CHECK(GetTetrahedronSingularFaceMask({}) == 0);
}

TEST_CASE("Singular face mask saturates when several features share an element",
          "[singularsolver][singularelements][Serial]")
{
  // Explains the Phase 7 measurement. On the transmon the trace-aware exemption reduced
  // repair marks by only ~0.5%, far short of the 25-50% the single-feature face count
  // suggests. The reason is the multi-feature OR reduction: the transmon carries ~5.6 H1
  // and ~10.9 ND basis incidences per enriched element, i.e. many features per element,
  // and features with different singular entities exempt different faces. Their union
  // saturates.
  //
  // Two node features at distinct nodes already activate all four faces: feature at node a
  // exempts only the face opposite a, which feature at node b activates.
  for (int a = 0; a < 4; a++)
  {
    for (int b = 0; b < 4; b++)
    {
      if (a == b)
      {
        continue;
      }
      std::array<int, 4> first{a, 0, 0, 0}, second{b, 0, 0, 0};
      int next = 0;
      for (int node = 0; node < 4; node++)
      {
        if (node != a)
        {
          first[++next] = node;
        }
      }
      next = 0;
      for (int node = 0; node < 4; node++)
      {
        if (node != b)
        {
          second[++next] = node;
        }
      }
      fem::singular::ElementFeatureIncidence incidence;
      incidence.nodes.push_back({0, 0, first});
      incidence.nodes.push_back({1, 1, second});
      CAPTURE(a, b);
      CHECK(GetTetrahedronSingularFaceMask(incidence) == SingularFaceMaskAllActive);
    }
  }

  // Three edge features covering all four vertices likewise saturate, which is the regime a
  // sheet edge threading a tetrahedron produces.
  {
    fem::singular::ElementFeatureIncidence incidence;
    incidence.edges.push_back({0, 0, 0, {0, 1, 2, 3}});
    incidence.edges.push_back({1, 1, 1, {1, 2, 3, 0}});
    incidence.edges.push_back({2, 2, 2, {2, 3, 0, 1}});
    CHECK(GetTetrahedronSingularFaceMask(incidence) == SingularFaceMaskAllActive);
  }

  // A single feature is the only case that actually exempts anything, so the benefit is
  // confined to elements touched by exactly one feature.
  {
    fem::singular::ElementFeatureIncidence node_only;
    node_only.nodes.push_back({0, 0, {0, 1, 2, 3}});
    CHECK(GetTetrahedronSingularFaceMask(node_only) != SingularFaceMaskAllActive);
    fem::singular::ElementFeatureIncidence edge_only;
    edge_only.edges.push_back({0, 0, 0, {0, 1, 2, 3}});
    CHECK(GetTetrahedronSingularFaceMask(edge_only) != SingularFaceMaskAllActive);
  }
}

TEST_CASE("Conforming closure refines a superset of the seed with amplified error",
          "[singularsolver][singularelements][Serial]")
{
  // Pins the two quantities the AMR loop reports about conforming closure, which are easy
  // to conflate: the closure RATIO (refined parents / seed elements, a topology quantity)
  // and the ERROR AMPLIFICATION (theta_closure / theta_seed, an indicator quantity). Both
  // are computed here the same way basesolver.cpp does, from the refinement transforms.
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  auto serial_mesh = InternalSheetMesh();
  mfem::ParMesh mesh(Mpi::World(), serial_mesh);
  const int coarse_elements = mesh.GetNE();
  REQUIRE(coarse_elements > 1);

  // A synthetic indicator concentrated on one element, so the seed is a single element and
  // any closure beyond it is unambiguous.
  Vector indicator(coarse_elements);
  indicator = 1.0;
  indicator[0] = 10.0;
  double total_squared = 0.0;
  for (int element = 0; element < coarse_elements; element++)
  {
    total_squared += indicator[element] * indicator[element];
  }
  mfem::Array<int> marked;
  marked.Append(0);
  const double seed_squared = indicator[0] * indicator[0];

  mesh.GeneralRefinement(marked, -1, 0);
  const auto &transforms = mesh.GetRefinementTransforms();
  REQUIRE(transforms.embeddings.Size() == mesh.GetNE());

  // Count children per coarse parent; a parent with more than one child was refined.
  std::vector<int> children(coarse_elements, 0);
  for (int child = 0; child < transforms.embeddings.Size(); child++)
  {
    const int parent = transforms.embeddings[child].parent;
    REQUIRE(parent >= 0);
    REQUIRE(parent < coarse_elements);
    children[parent]++;
  }
  long long refined_parents = 0;
  double closure_squared = 0.0;
  for (int element = 0; element < coarse_elements; element++)
  {
    if (children[element] > 1)
    {
      refined_parents++;
      closure_squared += indicator[element] * indicator[element];
    }
  }

  // Every child maps to some parent, and the total child count exceeds the parent count.
  CHECK(mesh.GetNE() > coarse_elements);
  // The seed element must itself have been refined: closure is a SUPERSET of the seed.
  CHECK(children[0] > 1);
  // Conforming bisection propagates, so at least the seed is refined and typically more.
  CHECK(refined_parents >= marked.Size());
  // Error amplification is closure/seed and can never be below one, because the closure
  // contains the seed.
  const double amplification = closure_squared / seed_squared;
  CAPTURE(refined_parents, closure_squared, seed_squared, amplification);
  CHECK(amplification >= 1.0);
  // theta_closure is a fraction of the total and must be a valid fraction.
  const double theta_closure = closure_squared / total_squared;
  CAPTURE(theta_closure);
  CHECK(theta_closure > 0.0);
  CHECK(theta_closure <= 1.0);
}

TEST_CASE("Trace-active face marks survive the cross-rank orientation reduction",
          "[singularsolver][singularelements][Parallel]")
{
  // Regression test for signed cancellation. RT DOFs are oriented, so a naive signed
  // prolongation round trip can sum (+1)*1 + (-1)*1 = 0 across a shared face and report an
  // active face as inactive. Exercise a sheet feature spanning a partition so shared faces
  // are marked from both sides, and require the closure to remain nonempty and conforming.
  if (Mpi::Size(Mpi::World()) != 2)
  {
    return;
  }
  auto serial_mesh = InternalSheetMesh();
  auto serial_features = fem::singular::ExtractSerialSheetFeatures(serial_mesh, {7});
  std::vector<int> partition(serial_mesh.GetNE(), 0);
  for (int element = 0; element < serial_mesh.GetNE(); element++)
  {
    partition[element] = element % 2;
  }
  mfem::ParMesh mesh(Mpi::World(), serial_mesh, partition.data());
  auto source_vertex_ids =
      fem::singular::MapPartitionedSerialVertexIds(serial_mesh, mesh, partition.data());
  std::vector<fem::singular::GlobalVertexId> source_element_ids(mesh.GetNE());
  std::iota(source_element_ids.begin(), source_element_ids.end(), 0);
  const auto features = fem::singular::DistributeSerialSheetFeatures(
      serial_features, mesh, source_vertex_ids, source_element_ids);

  bool conforming = true;
  mfem::Array<int> repair;
  const auto closure = BuildSingularRefinementProtection(mesh, features, source_vertex_ids,
                                                         &conforming, &repair);

  // The unrefined partitioned mesh is conforming, and the closure must still contain the
  // enriched elements: a cancelled reduction would silently empty it.
  CHECK(conforming);
  int local_closure = 0, local_enriched = 0;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    local_closure += (closure.Size() > element && closure[element]) ? 1 : 0;
    const auto &incidence = features.elements[element];
    local_enriched += (!incidence.nodes.empty() || !incidence.edges.empty()) ? 1 : 0;
  }
  int global_closure = local_closure, global_enriched = local_enriched;
  Mpi::GlobalSum(1, &global_closure, Mpi::World());
  Mpi::GlobalSum(1, &global_enriched, Mpi::World());
  CAPTURE(global_closure, global_enriched);
  REQUIRE(global_enriched > 0);
  CHECK(global_closure >= global_enriched);
}

}  // namespace palace
