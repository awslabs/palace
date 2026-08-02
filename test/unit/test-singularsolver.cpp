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
}

}  // namespace palace
