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

TEST_CASE("Conforming AMR preserves singular source identities and line-tip incidence",
          "[singularsolver][singularelements][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  auto serial_mesh = LongInternalLineMesh();
  const auto serial_features =
      fem::singular::ExtractSerialLineTipFeatures(serial_mesh, {7});
  mfem::ParMesh mesh(Mpi::World(), serial_mesh);

  std::vector<fem::singular::GlobalVertexId> source_vertex_ids(mesh.GetNV());
  std::iota(source_vertex_ids.begin(), source_vertex_ids.end(), 0);
  std::vector<fem::singular::GlobalVertexId> source_element_ids(mesh.GetNE());
  std::iota(source_element_ids.begin(), source_element_ids.end(), 0);
  auto local_features = fem::singular::DistributeSerialLineTipFeatures(
      serial_features, mesh, source_vertex_ids, source_element_ids);
  const auto protection =
      BuildSingularRefinementProtection(mesh, local_features, source_vertex_ids);

  mfem::Array<int> marked;
  for (int element = protection.Size() - 1; element >= 0; element--)
  {
    if (!protection[element])
    {
      marked.Append(element);
      break;
    }
  }
  REQUIRE(marked.Size() == 1);
  const int old_vertices = mesh.GetNV();
  mesh.GeneralRefinement(marked, -1, 1);
  REQUIRE_FALSE(mesh.Nonconforming());

  UpdateSingularSourceEntityIds(mesh, source_vertex_ids, source_element_ids);
  REQUIRE(source_vertex_ids.size() == static_cast<std::size_t>(mesh.GetNV()));
  for (int vertex = 0; vertex < old_vertices; vertex++)
  {
    CHECK(source_vertex_ids[vertex] == vertex);
  }
  CHECK(std::set(source_vertex_ids.begin(), source_vertex_ids.end()).size() ==
        source_vertex_ids.size());

  local_features = fem::singular::DistributeSerialLineTipFeatures(
      serial_features, mesh, source_vertex_ids, source_element_ids);
  CHECK_NOTHROW(BuildSingularRefinementProtection(mesh, local_features, source_vertex_ids));
}

TEST_CASE("Nonconforming AMR preserves singular source identities and a conforming patch",
          "[singularsolver][singularelements][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  auto serial_mesh = LongInternalLineMesh();
  const int original_vertices = serial_mesh.GetNV();
  const auto serial_features =
      fem::singular::ExtractSerialLineTipFeatures(serial_mesh, {7});
  serial_mesh.EnsureNCMesh(true);
  mfem::ParMesh mesh(Mpi::World(), serial_mesh);

  std::vector<fem::singular::GlobalVertexId> source_vertex_ids;
  std::vector<fem::singular::GlobalVertexId> source_element_ids;
  UpdateSingularSourceEntityIds(mesh, source_vertex_ids, source_element_ids);
  for (int vertex = 0; vertex < original_vertices; vertex++)
  {
    CHECK(source_vertex_ids[vertex] == vertex);
  }
  auto local_features = fem::singular::DistributeSerialLineTipFeatures(
      serial_features, mesh, source_vertex_ids, source_element_ids);
  for (int iteration = 0; iteration < 2; iteration++)
  {
    mfem::Array<int> protection;
    INFO("Checking the singular patch before nonconforming refinement " << iteration + 1);
    REQUIRE_NOTHROW(protection = BuildSingularRefinementProtection(mesh, local_features,
                                                                   source_vertex_ids));

    mfem::Array<int> marked;
    for (int element = protection.Size() - 1; element >= 0; element--)
    {
      if (!protection[element])
      {
        marked.Append(element);
        break;
      }
    }
    REQUIRE(marked.Size() == 1);
    mesh.GeneralRefinement(marked, -1, 0);
    REQUIRE(mesh.Nonconforming());

    UpdateSingularSourceEntityIds(mesh, source_vertex_ids, source_element_ids);
    CHECK(std::set(source_vertex_ids.begin(), source_vertex_ids.end()).size() ==
          source_vertex_ids.size());
    local_features = fem::singular::DistributeSerialLineTipFeatures(
        serial_features, mesh, source_vertex_ids, source_element_ids);
    INFO("Checking the singular patch after nonconforming refinement " << iteration + 1);
    CHECK_NOTHROW(
        BuildSingularRefinementProtection(mesh, local_features, source_vertex_ids));
  }
}

TEST_CASE("Three-dimensional refinement through a singular edge fails closed",
          "[singularsolver][singularelements][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  auto serial_mesh = InternalSheetMesh();
  const auto serial_features = fem::singular::ExtractSerialSheetFeatures(serial_mesh, {7});
  mfem::ParMesh mesh(Mpi::World(), serial_mesh);
  std::vector<fem::singular::GlobalVertexId> source_vertex_ids(mesh.GetNV());
  std::iota(source_vertex_ids.begin(), source_vertex_ids.end(), 0);
  std::vector<fem::singular::GlobalVertexId> source_element_ids(mesh.GetNE());
  std::iota(source_element_ids.begin(), source_element_ids.end(), 0);

  mfem::Array<int> marked(1);
  marked[0] = 0;
  mesh.GeneralRefinement(marked, -1, 1);
  UpdateSingularSourceEntityIds(mesh, source_vertex_ids, source_element_ids);
  CHECK_THROWS_WITH(
      fem::singular::DistributeSerialSheetFeatures(serial_features, mesh, source_vertex_ids,
                                                   source_element_ids),
      Catch::Matchers::ContainsSubstring(
          "refinement changed a protected three-dimensional singular feature"));
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

  mfem::Array<int> marked(1);
  marked[0] = 0;
  mesh.GeneralRefinement(marked, -1, 0);

  std::vector<fem::singular::GlobalVertexId> source_vertex_ids;
  std::vector<fem::singular::GlobalVertexId> source_element_ids;
  UpdateSingularSourceEntityIds(mesh, source_vertex_ids, source_element_ids);

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
  int repeated_coordinates = 0;
  for (int vertex = 0; vertex < global_vertices; vertex++)
  {
    const std::array<double, 2> coordinate{coordinates[2 * vertex],
                                           coordinates[2 * vertex + 1]};
    const auto [record, inserted] =
        identity_by_coordinate.emplace(coordinate, identities[vertex]);
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
  const auto serial_features =
      fem::singular::ExtractSerialLineTipFeatures(serial_mesh, {7});
  const auto serial_dofs =
      fem::singular::BuildSerialTriangleDofTopology(serial_mesh, serial_features, 1);
  if (nonconforming)
  {
    serial_mesh.EnsureNCMesh(true);
  }

  std::vector<int> partition(serial_mesh.GetNE(), 0);
  auto parallel_mesh =
      std::make_unique<mfem::ParMesh>(Mpi::World(), serial_mesh, partition.data());
  mesh::PartitionMetadata metadata;
  if (nonconforming)
  {
    UpdateSingularSourceEntityIds(*parallel_mesh, metadata.source_vertex_ids,
                                  metadata.source_element_ids);
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
    CHECK(numbering.h1.global_size ==
          static_cast<HYPRE_BigInt>(serial_dofs.h1_dofs.size()));
    CHECK(numbering.nd.global_size ==
          static_cast<HYPRE_BigInt>(serial_dofs.nd_dofs.size()));
    CHECK_NOTHROW(
        BuildSingularRefinementProtection(mesh, local_features, source.source_vertex_ids));
    int enriched_elements = 0;
    for (const auto &element : local_features.elements)
    {
      enriched_elements += !element.nodes.empty();
    }
    Mpi::GlobalSum(1, &enriched_elements, Mpi::World());
    return enriched_elements;
  };

  const auto initial_features = fem::singular::DistributeSerialLineTipFeatures(
      serial_features, *parallel_mesh, metadata.source_vertex_ids,
      metadata.source_element_ids);
  int initial_enriched_elements = 0;
  for (const auto &element : initial_features.elements)
  {
    initial_enriched_elements += !element.nodes.empty();
  }
  Mpi::GlobalSum(1, &initial_enriched_elements, Mpi::World());
  IoData iodata(Units(1.0, 1.0));
  iodata.model.refinement.maximum_imbalance = 1.01;
  const double initial_imbalance = mesh::RebalanceMesh(iodata, parallel_mesh, &metadata);
  CHECK(std::isinf(initial_imbalance));
  CHECK(parallel_mesh->GetNE() > 0);
  CHECK(CheckTopology(*parallel_mesh, metadata) == initial_enriched_elements);
}

}  // namespace palace
