// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <cmath>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "models/materialoperator.hpp"
#include "models/surfacepostoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/interfaceregion.hpp"
#include "utils/units.hpp"

namespace palace
{

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace
{

constexpr std::array<double, 3> kZ = {0.0, 0.0, 1.0};

bool Inside(const InterfaceRegion &region, double x, double y, double z)
{
  const double point[3] = {x, y, z};
  return region.Contains(point);
}

}  // namespace

// Decision 352 follow-up (2): the region's membership rule — the axis-aligned box, the
// translational cell of a segment (along-range between the perpendicular end cuts,
// in-plane transverse distance with the normal projected out, closed intervals), the
// union over segments intersected with the box, the plane normal, and the refused inputs.
TEST_CASE("InterfaceRegion membership", "[interfaceregion][Serial][Parallel]")
{
  const std::optional<std::array<double, 3>> none;
  // The box alone.
  {
    InterfaceRegion box(std::array<double, 3>{0.0, -1.0, 2.0},
                        std::array<double, 3>{4.0, 1.0, 3.0}, {}, 0.0, kZ);
    CHECK(box.HasBox());
    CHECK(box.SegmentCount() == 0);
    CHECK(Inside(box, 2.0, 0.0, 2.5));
    CHECK(Inside(box, 0.0, -1.0, 2.0));  // closed
    CHECK(Inside(box, 4.0, 1.0, 3.0));
    CHECK_FALSE(Inside(box, 4.1, 0.0, 2.5));
    CHECK_FALSE(Inside(box, 2.0, 0.0, 1.9));
  }
  // One segment along x from (0, 0, 0) to (10, 0, 0), transverse distance 2: the strip
  // |y| <= 2 for 0 <= x <= 10 at ANY z (in-plane geometry, the z normal projected out).
  {
    InterfaceRegion cell(none, none, {{0.0, 0.0, 0.0, 10.0, 0.0, 0.0}}, 2.0, kZ);
    CHECK_FALSE(cell.HasBox());
    CHECK(cell.SegmentCount() == 1);
    CHECK(Inside(cell, 5.0, 1.5, 0.0));
    CHECK(Inside(cell, 5.0, -2.0, 7.0));  // transverse closed, z free
    CHECK(Inside(cell, 0.0, 0.0, 0.0));   // along closed at the start cut
    CHECK(Inside(cell, 10.0, 2.0, 0.0));  // ... and at the end cut
    CHECK_FALSE(Inside(cell, 5.0, 2.01, 0.0));
    CHECK_FALSE(Inside(cell, -0.01, 0.0, 0.0));  // beyond the start cut, however close
    CHECK_FALSE(Inside(cell, 10.01, 0.0, 0.0));
    CHECK_FALSE(Inside(cell, 11.0, 0.5, 0.0));  // the end cut is flat, not a ball
  }
  // A diagonal segment: the along-range and the transverse distance follow its direction;
  // the orientation of the segment does not matter.
  {
    const double s = std::sqrt(0.5);
    InterfaceRegion diagonal(none, none, {{0.0, 0.0, 0.0, 10.0, 10.0, 0.0}}, 1.0, kZ);
    InterfaceRegion reversed(none, none, {{10.0, 10.0, 0.0, 0.0, 0.0, 0.0}}, 1.0, kZ);
    for (const auto *region : {&diagonal, &reversed})
    {
      CHECK(Inside(*region, 5.0, 5.0, 0.0));
      CHECK(Inside(*region, 5.0 - 0.99 * s, 5.0 + 0.99 * s, 0.0));
      CHECK_FALSE(Inside(*region, 5.0 - 1.01 * s, 5.0 + 1.01 * s, 0.0));
      CHECK(Inside(*region, 10.0, 10.0, 0.0));        // the far end cut, closed
      CHECK_FALSE(Inside(*region, 10.5, 10.5, 0.0));  // beyond the end cut
    }
  }
  // The union over segments, intersected with the box.
  {
    InterfaceRegion region(
        std::array<double, 3>{0.0, -5.0, -0.5}, std::array<double, 3>{6.0, 5.0, 0.5},
        {{0.0, 0.0, 0.0, 10.0, 0.0, 0.0}, {0.0, 4.0, 0.0, 10.0, 4.0, 0.0}}, 1.0, kZ);
    CHECK(Inside(region, 3.0, 0.5, 0.0));
    CHECK(Inside(region, 3.0, 4.5, 0.0));
    CHECK_FALSE(Inside(region, 3.0, 2.0, 0.0));   // between the two cells
    CHECK_FALSE(Inside(region, 7.0, 0.0, 0.0));   // in a cell, outside the box
    CHECK_FALSE(Inside(region, 3.0, 0.0, 0.75));  // the box's z range
  }
  // A non-z normal: the in-plane geometry of a vertical interface (x normal projected
  // out): the segment along z at (y, z) = (0, 0)..(0, 10), distance 1 -> the strip |y| <= 1
  // for 0 <= z <= 10 at any x; a segment whose in-plane length vanishes is refused.
  {
    InterfaceRegion vertical(none, none, {{0.0, 0.0, 0.0, 0.0, 0.0, 10.0}}, 1.0,
                             std::array<double, 3>{2.0, 0.0, 0.0});
    CHECK(Inside(vertical, 3.0, 0.5, 5.0));
    CHECK_FALSE(Inside(vertical, 3.0, 1.5, 5.0));
    CHECK_FALSE(Inside(vertical, 3.0, 0.0, 10.5));
    CHECK_THROWS(InterfaceRegion(none, none, {{0.0, 0.0, 0.0, 0.0, 0.0, 10.0}}, 1.0, kZ));
  }
  // Refused inputs.
  CHECK_THROWS(InterfaceRegion(none, none, {}, 0.0, kZ));  // neither box nor segments
  CHECK_THROWS(InterfaceRegion(std::array<double, 3>{0.0, 0.0, 0.0}, none, {}, 0.0, kZ));
  CHECK_THROWS(InterfaceRegion(std::array<double, 3>{1.0, 0.0, 0.0},
                               std::array<double, 3>{0.0, 1.0, 1.0}, {}, 0.0, kZ));
  CHECK_THROWS(InterfaceRegion(none, none, {{0.0, 0.0, 0.0, 1.0, 0.0, 0.0}}, 0.0, kZ));
  CHECK_THROWS(InterfaceRegion(none, none, {{0.0, 0.0, 0.0, 1.0, 0.0, 0.0}}, 1.0,
                               std::array<double, 3>{0.0, 0.0, 0.0}));
  // A two-dimensional point lies in z = 0.
  {
    InterfaceRegion region(std::array<double, 3>{0.0, 0.0, -0.5},
                           std::array<double, 3>{1.0, 1.0, 0.5}, {}, 0.0, kZ);
    mfem::Vector point(2);
    point(0) = 0.5;
    point(1) = 0.5;
    CHECK(region.Contains(point));
    point(0) = 1.5;
    CHECK_FALSE(region.Contains(point));
  }
}

// The region as the quadrature-level filter of an interface dielectric entry: on the top
// face z = 1 of the unit cube (an MA interface, t = 0.2, eps = 2) under the exactly
// represented field E = (1 + x) z-hat, whose MA energy density is 0.05 (1 + x)^2 per unit
// area, (a) regions cut along mesh lines reproduce the analytic energies, (b) regions cut
// anywhere partition the total, the EdgeDistances outside / annulus energies, and the
// response matrices to roundoff (every quadrature point in exactly one cell), (c) an
// unfiltered entry is unchanged, (d) a segment at another z selects the plane (in-plane
// geometry) while the box's z range excludes it.
TEST_CASE("InterfaceRegion exhaustive surface quadrature",
          "[interfaceregion][Serial][Parallel]")
{
  auto serial = mfem::Mesh::MakeCartesian3D(2, 2, 1, mfem::Element::TETRAHEDRON);
  int top = -1;
  mfem::Vector center(3);
  for (int be = 0; be < serial.GetNBE(); be++)
  {
    auto *T = serial.GetBdrElementTransformation(be);
    T->Transform(mfem::Geometries.GetCenter(T->GetGeometryType()), center);
    if (std::abs(center(2) - 1.0) < 1e-12)
    {
      top = serial.GetBdrAttribute(be);
      break;
    }
  }
  REQUIRE(top > 0);
  Mesh mesh(std::make_unique<mfem::ParMesh>(Mpi::World(), serial));
  mfem::H1_FECollection h1(2, 3);
  mfem::ND_FECollection nd(2, 3);
  FiniteElementSpace h1_space(mesh, &h1), nd_space(mesh, &nd);
  fem::DefaultIntegrationOrder::p_trial = 2;
  config::MaterialData material;
  material.attributes = {1};
  config::PeriodicBoundaryData periodic;
  MaterialOperator materials({material}, periodic, ProblemType::ELECTROSTATIC, mesh);
  GridFunction field(nd_space, false);
  mfem::VectorFunctionCoefficient coefficient(3,
                                              [](const mfem::Vector &x, mfem::Vector &E)
                                              {
                                                E.SetSize(3);
                                                E = 0.0;
                                                E(2) = 1.0 + x(0);
                                              });
  field.Real().ProjectCoefficient(coefficient);

  config::InterfaceDielectricData data;
  data.attributes = {top};
  data.type = InterfaceDielectric::MA;
  data.t = 0.2;
  data.epsilon_r = 2.0;
  data.edge_attributes = {top};
  data.edge_distances = {0.25};
  // The analytic energy of the strip a <= x <= b over the full y range.
  auto Strip = [](double a, double b)
  { return 0.05 * (std::pow(1.0 + b, 3) - std::pow(1.0 + a, 3)) / 3.0; };

  config::BoundaryPostData plain;
  plain.dielectric.emplace(1, data);
  SurfacePostOperator unmasked(plain, ProblemType::ELECTROSTATIC, materials, h1_space,
                               nd_space);
  const double reference = unmasked.GetInterfaceElectricFieldEnergy(1, field);
  CHECK_THAT(reference, WithinRel(Strip(0.0, 1.0), 1e-12));
  const auto reference_edges = unmasked.GetInterfaceEdgeElectricFieldEnergies(1, field);
  REQUIRE(reference_edges.size() == 1);
  CHECK(reference_edges[0].energy_outside > 0.0);
  CHECK(reference_edges[0].energy_annulus > 0.0);
  CHECK(reference_edges[0].energy_outside < reference);

  // (a) Cuts along mesh lines: the box x <= 0.5, the cell of a segment along y at x = 0.25
  // with distance 0.25 (the same strip), the cell of the segment (0.5, 0, z) - (0.5, 0.5,
  // z) at z = 0 with distance 0.5 (the quadrant x in [0, 1], y in [0, 0.5]).
  config::BoundaryPostData aligned;
  {
    auto entry = data;
    entry.region.emplace();
    entry.region->box_min = std::array<double, 3>{0.0, 0.0, 0.5};
    entry.region->box_max = std::array<double, 3>{0.5, 1.0, 1.5};
    aligned.dielectric.emplace(1, entry);
    entry.region.reset();
    entry.region.emplace();
    entry.region->segments = {{0.25, 0.0, 1.0, 0.25, 1.0, 1.0}};
    entry.region->distance = 0.25;
    aligned.dielectric.emplace(2, entry);
    entry.region->segments = {{0.5, 0.0, 0.0, 0.5, 0.5, 0.0}};
    entry.region->distance = 0.5;
    aligned.dielectric.emplace(3, entry);
    // The same cell, but the box's z range [-0.5, 0.5] excludes the plane z = 1.
    entry.region->box_min = std::array<double, 3>{-1.0, -1.0, -0.5};
    entry.region->box_max = std::array<double, 3>{2.0, 2.0, 0.5};
    aligned.dielectric.emplace(4, entry);
    // (c) The unfiltered entry next to them.
    aligned.dielectric.emplace(5, data);
  }
  SurfacePostOperator aligned_post(aligned, ProblemType::ELECTROSTATIC, materials, h1_space,
                                   nd_space);
  CHECK_THAT(aligned_post.GetInterfaceElectricFieldEnergy(1, field),
             WithinRel(Strip(0.0, 0.5), 1e-12));
  CHECK_THAT(aligned_post.GetInterfaceElectricFieldEnergy(2, field),
             WithinRel(Strip(0.0, 0.5), 1e-12));
  CHECK_THAT(aligned_post.GetInterfaceElectricFieldEnergy(3, field),
             WithinRel(0.5 * Strip(0.0, 1.0), 1e-12));
  CHECK(aligned_post.GetInterfaceElectricFieldEnergy(4, field) == 0.0);
  CHECK(aligned_post.GetInterfaceElectricFieldEnergy(5, field) == reference);
  {
    const auto edges = aligned_post.GetInterfaceEdgeElectricFieldEnergies(5, field);
    REQUIRE(edges.size() == 1);
    CHECK(edges[0].energy_outside == reference_edges[0].energy_outside);
    CHECK(edges[0].energy_annulus == reference_edges[0].energy_annulus);
  }

  // (b) A partition by cuts through the faces: the cells of the segment along x at y = 0.37
  // (distance 0.37: y in [0, 0.74]) with along-ranges [0, 0.61] and [0.61, 1], and the box
  // y >= 0.74; every quadrature point of the plane lies in exactly one of them.
  config::BoundaryPostData partition;
  {
    auto entry = data;
    entry.region.emplace();
    entry.region->segments = {{0.0, 0.37, 1.0, 0.61, 0.37, 1.0}};
    entry.region->distance = 0.37;
    partition.dielectric.emplace(1, entry);
    entry.region->segments = {{0.61, 0.37, 1.0, 1.0, 0.37, 1.0}};
    partition.dielectric.emplace(2, entry);
    entry.region.reset();
    entry.region.emplace();
    entry.region->box_min = std::array<double, 3>{-1.0, 0.74, 0.5};
    entry.region->box_max = std::array<double, 3>{2.0, 2.0, 1.5};
    partition.dielectric.emplace(3, entry);
  }
  SurfacePostOperator partition_post(partition, ProblemType::ELECTROSTATIC, materials,
                                     h1_space, nd_space);
  double total = 0.0, outside = 0.0, annulus = 0.0;
  for (int idx : {1, 2, 3})
  {
    const double energy = partition_post.GetInterfaceElectricFieldEnergy(idx, field);
    CHECK(energy > 0.0);
    CHECK(energy < reference);
    total += energy;
    const auto edges = partition_post.GetInterfaceEdgeElectricFieldEnergies(idx, field);
    REQUIRE(edges.size() == 1);
    outside += edges[0].energy_outside;
    annulus += edges[0].energy_annulus;
  }
  CHECK_THAT(total, WithinRel(reference, 1e-12));
  CHECK_THAT(outside, WithinRel(reference_edges[0].energy_outside, 1e-12));
  CHECK_THAT(annulus, WithinRel(reference_edges[0].energy_annulus, 1e-12));

  // The response matrices (localized edge energies) of the partition sum to the unfiltered
  // ones; a region filters surfaces only (no localized volume diagnostics).
  auto matrix_data = partition;
  matrix_data.dielectric.emplace(4, data);
  for (auto &[index, entry] : matrix_data.dielectric)
  {
    entry.localize_edge_energy = true;
    entry.save_local_edge_energy = false;
    entry.edge_frame_normal = std::array<double, 3>{0.0, 0.0, 1.0};
  }
  SurfacePostOperator matrix_post(matrix_data, ProblemType::ELECTROSTATIC, materials,
                                  h1_space, nd_space);
  GridFunction second(nd_space, false);
  mfem::Vector value(3);
  value = 0.0;
  value(2) = 2.0;
  mfem::VectorConstantCoefficient second_coefficient(value);
  second.Real().ProjectCoefficient(second_coefficient);
  const auto matrices =
      matrix_post.GetInterfaceElectricFieldEnergyMatrices({&field, &second}, {});
  REQUIRE(matrices.size() == 4);
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      double sum_total = 0.0, sum_inside = 0.0;
      for (int idx : {1, 2, 3})
      {
        sum_total += matrices.at(idx).at(0).energy_total(i, j);
        sum_inside += matrices.at(idx).at(0).energy_inside(i, j);
      }
      CHECK_THAT(sum_total, WithinRel(matrices.at(4).at(0).energy_total(i, j), 1e-12));
      CHECK_THAT(sum_inside, WithinRel(matrices.at(4).at(0).energy_inside(i, j), 1e-12));
    }
  }
  CHECK_THROWS(
      matrix_post.GetInterfaceLocalEdgeElectricFieldEnergies(1, field, nullptr,
                                                             /*include_volume=*/true));

  // The configuration's region is scaled with the mesh.
  {
    auto scaled = data;
    scaled.region.emplace();
    scaled.region->box_min = std::array<double, 3>{0.0, 0.0, 0.0};
    scaled.region->box_max = std::array<double, 3>{4.0, 8.0, 12.0};
    scaled.region->segments = {{0.0, 0.0, 0.0, 4.0, 0.0, 0.0}};
    scaled.region->distance = 2.0;
    scaled.region->normal = std::array<double, 3>{0.0, 0.0, 3.0};
    Units units(1e-6, 4e-6);
    config::Nondimensionalize(units, scaled);
    CHECK((*scaled.region->box_max)[0] == 1.0);
    CHECK((*scaled.region->box_max)[2] == 3.0);
    CHECK(scaled.region->segments[0][3] == 1.0);
    CHECK(scaled.region->distance == 0.5);
    CHECK(scaled.region->normal[2] == 3.0);  // a direction, not a length
  }
}

}  // namespace palace
