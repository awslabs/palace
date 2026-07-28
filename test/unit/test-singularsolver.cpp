// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "drivers/singularsolver.hpp"
#include "fem/singularfeatures.hpp"
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

}  // namespace palace
