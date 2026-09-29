// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// The corner family's angle-interpolation rule (corner-qualification block 2026-09-29;
// design doc SURFACE-RESPONSE-IDENTIFICATION.md, Conventions CornerTraceBasisRule): the
// geometric events of the trace basis (a knot of a metal ring passing a fixed-layout
// vertex: the hats jump there), the segment connectivity that fixes the band triangulation
// over a segment, and the stencil rule that never straddles a knot-corner passage (legacy
// coupons exact only, per-side coupons at the passages, fail closed on a segment across a
// passage). Pinned to the same numbers as the Python mirror corner_family_interpolation.py
// (test_corner_family_interpolation.py).

#include "fixtures.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/mesh.hpp"
#include "models/cornertracebasis.hpp"
#include "models/surfaceresponseoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

namespace fs = std::filesystem;

using json = nlohmann::json;
using namespace Catch::Matchers;

namespace
{

constexpr double kR = 1.9, kT = 0.1, kOE = 0.05;
constexpr double kDeg = M_PI / 180.0;

double FractionDistanceForTest(double a, double b)
{
  const double d = std::abs(std::fmod(a - b + 2.0, 1.0));
  return std::min(d, 1.0 - d);
}

std::set<std::array<int, 3>> TriangleSet(const ConstructedCornerTraceBasis &basis)
{
  std::set<std::array<int, 3>> set;
  for (auto t : basis.triangles)
  {
    std::sort(t.begin(), t.end());
    set.insert(t);
  }
  return set;
}

ConstructedCornerTraceBasis Build(double angle, bool convex,
                                  std::optional<double> connectivity = std::nullopt)
{
  const CornerTraceBasisRule rule;
  const auto seed = MakeCornerBoxSeed(kR, kT, kOE, convex, rule);
  std::optional<double> connectivity_radians;
  if (connectivity)
  {
    connectivity_radians = *connectivity * kDeg;
  }
  return BuildCornerTraceBasis(seed.points, seed.contour_groups, seed.zero_trace_indices,
                               angle * kDeg, convex, rule, connectivity_radians);
}

std::vector<double> CornerEvents(bool convex, double min_angle)
{
  std::vector<double> angles;
  for (const auto &event : CornerBasisEvents(convex, CornerTraceBasisRule{}))
  {
    if (event.corner && event.angle_degrees >= min_angle &&
        (angles.empty() || event.angle_degrees - angles.back() > 1.0e-6))
    {
      angles.push_back(event.angle_degrees);
    }
  }
  return angles;
}

// The qualified family's node set of one convexity: (angle, connectivity angle) per
// segment.
std::vector<CornerFamilyNode> QualifiedFamily(bool convex)
{
  const double passage = convex ? 153.434948822922 : 158.198590513648;
  const std::vector<std::pair<double, double>> nodes =
      convex ? std::vector<std::pair<double, double>>{{75.0, 82.5},     {80.0, 82.5},
                                                      {85.0, 82.5},     {90.0, 82.5},
                                                      {90.0, 112.5},    {105.0, 112.5},
                                                      {120.0, 112.5},   {135.0, 112.5},
                                                      {135.0, 144.2},   {144.0, 144.2},
                                                      {150.0, 144.2},   {passage, 144.2},
                                                      {passage, 166.7}, {159.0, 166.7},
                                                      {165.0, 166.7},   {180.0, 166.7}}
             : std::vector<std::pair<double, double>>{
                   {75.0, 82.5},     {80.0, 82.5},   {85.0, 82.5},   {90.0, 82.5},
                   {90.0, 112.5},    {105.0, 112.5}, {120.0, 112.5}, {135.0, 112.5},
                   {135.0, 146.6},   {143.0, 146.6}, {150.0, 146.6}, {passage, 146.6},
                   {passage, 169.1}, {165.0, 169.1}, {172.0, 169.1}, {180.0, 169.1}};
  std::vector<CornerFamilyNode> family;
  for (std::size_t i = 0; i < nodes.size(); i++)
  {
    family.push_back({nodes[i].first, nodes[i].second, i});
  }
  return family;
}

}  // namespace

TEST_CASE("CornerFamilyBasisEvents", "[cornerfamily][Serial][Parallel]")
{
  const CornerTraceBasisRule rule;
  // Knot-corner passages in the family range: convex 90 (the fixed layout: free 1 / 3 / 5
  // and metal 1 at the four corners), 135 (the second crossing at (-R, R)), 153.435 =
  // 180 - atan(1 / 2) (free 2 at (-R, -R)); concave 90 (free 3 at (R, R), metal 1 at
  // (-R, -R)), 135 (free 2 at (R, R) and the crossing), 158.199 = 180 - atan(2 / 5) (free 5
  // at (-R, R)).
  const auto convex = CornerEvents(true, 75.0);
  const auto concave = CornerEvents(false, 75.0);
  REQUIRE(convex.size() == 3);
  REQUIRE(concave.size() == 3);
  CHECK_THAT(convex[0], WithinAbs(90.0, 1.0e-9));
  CHECK_THAT(convex[1], WithinAbs(135.0, 1.0e-9));
  CHECK_THAT(convex[2], WithinAbs(153.434948822922, 1.0e-9));
  CHECK_THAT(convex[2], WithinAbs(180.0 - std::atan(0.5) / kDeg, 1.0e-9));
  CHECK_THAT(concave[0], WithinAbs(90.0, 1.0e-9));
  CHECK_THAT(concave[1], WithinAbs(135.0, 1.0e-9));
  CHECK_THAT(concave[2], WithinAbs(158.198590513648, 1.0e-9));
  CHECK_THAT(concave[2], WithinAbs(180.0 - std::atan(0.4) / kDeg, 1.0e-9));
  // Side-midpoint passages inside the range (removed by the segment connectivity): convex
  // free 1 through (-R, 0) at 141.340 = 180 - atan(4 / 5), concave free 5 through (0, R) at
  // 111.801 = 180 - atan(5 / 2); the anchor's crossing reaches (-R, 0) at 180.
  for (const bool is_convex : {true, false})
  {
    // (Midpoint passages at a knot-corner passage angle, e.g. the fixed layout at 90, are
    // absorbed by that event.)
    std::vector<CornerBasisEvent> midpoints;
    const auto corners = CornerEvents(is_convex, 0.0);
    for (const auto &event : CornerBasisEvents(is_convex, rule))
    {
      const bool at_corner_event =
          std::any_of(corners.begin(), corners.end(), [&](double corner)
                      { return std::abs(corner - event.angle_degrees) <= 1.0e-6; });
      if (!event.corner && !at_corner_event && event.angle_degrees > 75.0 &&
          event.angle_degrees < 180.0)
      {
        midpoints.push_back(event);
      }
    }
    REQUIRE(midpoints.size() == 1);
    CHECK(midpoints[0].role == (is_convex ? "free1" : "free5"));
    CHECK_THAT(midpoints[0].fraction, WithinAbs(is_convex ? 0.0 : 0.75, 1.0e-12));
    CHECK_THAT(midpoints[0].angle_degrees,
               WithinAbs(is_convex ? 141.340191745910 : 111.801409486352, 1.0e-9));
    // Every event is a real passage: the knot of that role sits on the fixed vertex.
    for (const auto &event : CornerBasisEvents(is_convex, rule))
    {
      if (event.angle_degrees < 75.0 || event.angle_degrees > 180.0)
      {
        continue;
      }
      const auto layout =
          CornerMetalRingLayout(kR, event.angle_degrees * kDeg, is_convex, rule);
      bool found = false;
      for (const auto &vertex : layout)
      {
        if (vertex.kind != CornerRingVertex::Kind::SLAVE &&
            FractionDistanceForTest(vertex.fraction, event.fraction) <= 1.0e-9)
        {
          found = true;
        }
      }
      CHECK(found);
    }
  }
}

TEST_CASE("CornerFamilySegmentConnectivity", "[cornerfamily][Serial][Parallel]")
{
  // With a connectivity angle the bands next to the metal rings are merged in the order of
  // the rule's layout at that angle: one triangulation over the segment [135, 153.435]
  // (keys 144.2) although free 1 passes the side midpoint (-R, 0) at 141.34 inside it; the
  // knots stay at the angle's positions and the slave table is the perimeter-ordered one;
  // without a connectivity angle the construction is the perimeter-ordered merge (equal to
  // the keyed one once the knot has passed the midpoint, different before).
  std::optional<std::set<std::array<int, 3>>> reference;
  for (const double angle : {137.0, 141.0, 142.0, 150.0, 153.0})
  {
    const auto keyed = Build(angle, true, 144.2);
    const auto plain = Build(angle, true);
    REQUIRE(keyed.knots.size() == plain.knots.size());
    for (std::size_t k = 0; k < keyed.knots.size(); k++)
    {
      CHECK(keyed.knots[k] == plain.knots[k]);
    }
    CHECK(keyed.zero == plain.zero);
    REQUIRE(keyed.vertices.size() == plain.vertices.size());
    for (std::size_t v = 0; v < keyed.vertices.size(); v++)
    {
      CHECK(keyed.vertices[v].point == plain.vertices[v].point);
      CHECK(keyed.vertices[v].parent_a == plain.vertices[v].parent_a);
      CHECK(keyed.vertices[v].parent_b == plain.vertices[v].parent_b);
    }
    CHECK(keyed.triangles.size() == plain.triangles.size());
    const auto set = TriangleSet(keyed);
    CHECK((TriangleSet(plain) == set) == (angle > 141.3402));
    if (!reference)
    {
      reference = set;
    }
    CHECK(set == *reference);
  }
  // Across a knot-corner passage the keyed merge would fold: refused.
  CHECK_THROWS(Build(130.0, true, 144.2));
  CHECK_THROWS(Build(160.0, true, 144.2));
  // The recorded coupons (perimeter order at their own angle) at an event angle equal
  // neither neighbouring segment's triangulation; away from events they equal their
  // segment's.
  CHECK(TriangleSet(Build(90.0, true)) != TriangleSet(Build(90.0, true, 82.5)));
  CHECK(TriangleSet(Build(90.0, true)) != TriangleSet(Build(90.0, true, 112.5)));
  CHECK(TriangleSet(Build(90.0, true, 82.5)) != TriangleSet(Build(90.0, true, 112.5)));
  CHECK(TriangleSet(Build(150.0, true)) == TriangleSet(Build(150.0, true, 144.2)));
  CHECK(TriangleSet(Build(105.0, true)) == TriangleSet(Build(105.0, true, 112.5)));
  CHECK(TriangleSet(Build(120.0, false)) == TriangleSet(Build(120.0, false, 112.5)));
  CHECK(TriangleSet(Build(105.0, false)) != TriangleSet(Build(105.0, false, 112.5)));
}

TEST_CASE("CornerFamilyStencil", "[cornerfamily][Serial][Parallel]")
{
  const CornerTraceBasisRule rule;
  constexpr double tol = 1.0e-2;  // kSignatureAngleToleranceDegrees
  for (const bool convex : {true, false})
  {
    const auto family = QualifiedFamily(convex);
    const auto passages = CornerEvents(convex, 75.0);
    // Cubic in every segment; the stencil's nodes share the segment's connectivity angle
    // and lie on the angle's side of every knot-corner passage; the weights sum to one.
    for (const double angle : {78.0, 82.5, 100.0, 112.5, 127.5, 142.5, 147.0, 172.5, 176.0})
    {
      const auto stencil = SelectCornerFamilyStencil(family, angle, convex, rule, tol);
      REQUIRE(stencil.reason.empty());
      CHECK(stencil.rule == "cubic");
      REQUIRE(stencil.connectivity_angle_degrees.has_value());
      double sum = 0.0;
      for (const auto &[index, weight] : stencil.nodes)
      {
        sum += weight;
        CHECK(family[index].connectivity_angle_degrees ==
              stencil.connectivity_angle_degrees);
        for (const double passage : passages)
        {
          const bool straddles =
              std::min(angle, family[index].angle_degrees) + tol < passage &&
              passage < std::max(angle, family[index].angle_degrees) - tol;
          CHECK_FALSE(straddles);
        }
      }
      CHECK_THAT(sum, WithinAbs(1.0, 1.0e-12));
      CHECK(std::abs(family[stencil.base].angle_degrees - angle) <= 7.5 + tol);
    }
    // The [75, 90] segment: cubic on 75 / 80 / 85 / 90- (was linear on 75 / 90).
    const auto held_out = SelectCornerFamilyStencil(family, 82.5, convex, rule, tol);
    REQUIRE(held_out.nodes.size() == 4);
    CHECK(family[held_out.nodes[0].first].angle_degrees == 75.0);
    CHECK(family[held_out.nodes[3].first].angle_degrees == 90.0);
    CHECK(*held_out.connectivity_angle_degrees == 82.5);
    // Exact node at an event angle: the lower-angle segment's coupon; with a legacy tie
    // coupon at that angle the legacy one.
    auto exact = SelectCornerFamilyStencil(family, 90.0, convex, rule, tol);
    CHECK(exact.rule == "exact");
    REQUIRE(exact.nodes.size() == 1);
    CHECK(family[exact.nodes[0].first].connectivity_angle_degrees == 82.5);
    exact = SelectCornerFamilyStencil(family, 135.0, convex, rule, tol);
    CHECK(family[exact.nodes[0].first].connectivity_angle_degrees == 112.5);
    auto with_legacy = family;
    with_legacy.push_back({90.0, std::nullopt, 99});
    exact = SelectCornerFamilyStencil(with_legacy, 90.0, convex, rule, tol);
    CHECK(exact.rule == "exact");
    CHECK(exact.nodes.at(0).first == 99);
    CHECK_FALSE(exact.connectivity_angle_degrees.has_value());
    // Refusals: outside the range (no extrapolation).
    CHECK_THAT(SelectCornerFamilyStencil(family, 70.0, convex, rule, tol).reason,
               ContainsSubstring("sharper than"));
    std::vector<CornerFamilyNode> no_anchor;
    for (const auto &node : family)
    {
      if (node.angle_degrees < 180.0)
      {
        no_anchor.push_back(node);
      }
    }
    CHECK_THAT(SelectCornerFamilyStencil(no_anchor, 175.0, convex, rule, tol).reason,
               ContainsSubstring("wider than"));
    // A gap between segments (no coupon on the upper side of 135).
    std::vector<CornerFamilyNode> gap;
    for (const auto &node : family)
    {
      if (*node.connectivity_angle_degrees != (convex ? 144.2 : 146.6))
      {
        gap.push_back(node);
      }
    }
    CHECK_THAT(SelectCornerFamilyStencil(gap, 142.5, convex, rule, tol).reason,
               ContainsSubstring("no segment"));
  }
  // Legacy family (the recorded coupons, no connectivity records): exact matches only, an
  // interpolation refused with the reason; a legacy coupon beside keyed segments does not
  // extend them.
  {
    std::vector<CornerFamilyNode> legacy;
    std::size_t index = 0;
    for (const double angle : {75.0, 90.0, 105.0, 120.0, 135.0, 150.0, 165.0, 180.0})
    {
      legacy.push_back({angle, std::nullopt, index++});
    }
    CHECK(SelectCornerFamilyStencil(legacy, 120.0, true, rule, tol).rule == "exact");
    CHECK_THAT(SelectCornerFamilyStencil(legacy, 112.5, true, rule, tol).reason,
               ContainsSubstring("without segment connectivity records"));
    CHECK_THAT(SelectCornerFamilyStencil(legacy, 172.5, true, rule, tol).reason,
               ContainsSubstring("without segment connectivity records"));
    auto mixed = QualifiedFamily(true);
    mixed.push_back({60.0, std::nullopt, 99});
    CHECK_THAT(SelectCornerFamilyStencil(mixed, 70.0, true, rule, tol).reason,
               ContainsSubstring("without segment connectivity"));
  }
  // Fail closed: a segment across a knot-corner passage, a connectivity angle on one, two
  // coupons at one angle in a segment, overlapping segments.
  CHECK_THROWS(SelectCornerFamilyStencil({{120.0, 112.5, 0}, {150.0, 112.5, 1}}, 130.0,
                                         true, rule, tol));
  CHECK_THROWS(SelectCornerFamilyStencil({{120.0, 135.0, 0}, {130.0, 135.0, 1}}, 125.0,
                                         true, rule, tol));
  CHECK_THROWS(SelectCornerFamilyStencil(
      {{105.0, 112.5, 0}, {105.0, 112.5, 1}, {120.0, 112.5, 2}}, 110.0, true, rule, tol));
  CHECK_THROWS(SelectCornerFamilyStencil(
      {{95.0, 100.0, 0}, {125.0, 100.0, 1}, {110.0, 115.0, 2}, {130.0, 115.0, 3}}, 120.0,
      true, rule, tol));
  // Lower orders where a segment has fewer nodes.
  const std::vector<CornerFamilyNode> shortened = {
      {135.0, 144.2, 0}, {150.0, 144.2, 1}, {153.434948822922, 144.2, 2}};
  CHECK(SelectCornerFamilyStencil(shortened, 142.5, true, rule, tol).rule == "quadratic");
  CHECK(SelectCornerFamilyStencil({shortened[0], shortened[1]}, 142.5, true, rule, tol)
            .rule == "linear");
}

TEST_CASE("CornerFamilyNodeToleranceBoundary", "[cornerfamily][Serial][Parallel]")
{
  // Qualification review 2026-09-29 m1: an angle at |node - angle| = tolerance (in floating
  // point, however it was formed: node +- tol, the decimal literal, the review's 74.5 +
  // 0.005 i grid) is EXACT when the difference is <= tol and interior (cubic) when > tol,
  // by the same difference in both tests, never "lies in no segment"; beyond the family
  // range by more than the tolerance it is refused, within the tolerance it is the end
  // node.
  const CornerTraceBasisRule rule;
  constexpr double tol = 1.0e-2;
  for (const bool convex : {true, false})
  {
    const auto family = QualifiedFamily(convex);
    std::set<double> angles;
    for (const auto &node : family)
    {
      const double n = node.angle_degrees;
      for (const double angle : {n + tol, n - tol, std::round((n + tol) * 1.0e6) / 1.0e6,
                                 std::round((n - tol) * 1.0e6) / 1.0e6,
                                 n + tol * (1.0 - 1.0e-12), n - tol * (1.0 - 1.0e-12),
                                 n + tol * (1.0 + 1.0e-12), n - tol * (1.0 + 1.0e-12)})
      {
        angles.insert(angle);
      }
    }
    for (int i = 0; i <= 21100; i++)
    {
      const double angle = 74.5 + 0.005 * i;
      for (const auto &node : family)
      {
        if (std::abs(angle - node.angle_degrees) <= 2.0 * tol)
        {
          angles.insert(angle);
        }
      }
    }
    int exact_count = 0, interior_count = 0, refused_count = 0;
    for (const double angle : angles)
    {
      bool exact = false;
      for (const auto &node : family)
      {
        exact = exact || !(std::abs(node.angle_degrees - angle) > tol);
      }
      const auto stencil = SelectCornerFamilyStencil(family, angle, convex, rule, tol);
      if (!exact && (angle < 75.0 || angle > 180.0))
      {
        CHECK_THAT(stencil.reason, ContainsSubstring("than"));
        refused_count++;
        continue;
      }
      CHECK(stencil.reason.empty());
      CHECK((stencil.rule == "exact") == exact);
      if (!exact)
      {
        CHECK(stencil.rule == "cubic");
        interior_count++;
      }
      else
      {
        exact_count++;
      }
    }
    CHECK(exact_count > 0);
    CHECK(interior_count > 0);
    CHECK(refused_count > 0);
  }
}

TEST_CASE("CornerFamilySegmentStructure", "[cornerfamily][Serial][Parallel]")
{
  // The segment structure check (decision 137 (1); verified at library load by
  // ReadProcessLibrary, CheckCornerFamilySegments): the qualified family with or without
  // the legacy tie coupons at 90 / 135 / 180 passes; a segment across a knot-corner
  // passage, a connectivity angle on one, two coupons at one angle in a segment and
  // overlapping segments are refused with the reason; a shared node angle between two
  // segments and a legacy coupon inside a segment's range are not overlaps.
  const CornerTraceBasisRule rule;
  constexpr double tol = 1.0e-2;
  for (const bool convex : {true, false})
  {
    auto family = QualifiedFamily(convex);
    CHECK(CheckCornerFamilySegments(family, convex, rule, tol).empty());
    const std::vector<CornerFamilyNode> ties = {
        {90.0, std::nullopt, 90}, {135.0, std::nullopt, 91}, {180.0, std::nullopt, 92}};
    family.insert(family.end(), ties.begin(), ties.end());
    CHECK(CheckCornerFamilySegments(family, convex, rule, tol).empty());
    CHECK(CheckCornerFamilySegments(ties, convex, rule, tol).empty());
  }
  CHECK_THAT(
      CheckCornerFamilySegments({{120.0, 112.5, 0}, {150.0, 112.5, 1}}, true, rule, tol),
      ContainsSubstring("node at 150") &&
          ContainsSubstring("across a knot-corner passage"));
  CHECK_THAT(
      CheckCornerFamilySegments({{120.0, 135.0, 0}, {130.0, 135.0, 1}}, true, rule, tol),
      ContainsSubstring("lies on a knot-corner passage"));
  CHECK_THAT(
      CheckCornerFamilySegments({{105.0, 112.5, 0}, {105.0, 112.5, 1}, {120.0, 112.5, 2}},
                                true, rule, tol),
      ContainsSubstring("two coupons at 105"));
  CHECK_THAT(
      CheckCornerFamilySegments(
          {{95.0, 100.0, 0}, {125.0, 100.0, 1}, {110.0, 115.0, 2}, {130.0, 115.0, 3}}, true,
          rule, tol),
      ContainsSubstring("overlap beyond a shared node angle"));
  CHECK(
      CheckCornerFamilySegments(
          {{90.0, 82.5, 0}, {90.0, 112.5, 1}, {105.0, 112.5, 2}, {105.0, std::nullopt, 3}},
          true, rule, tol)
          .empty());
}

namespace
{

// A unit-test corner coupon on the rule's layout at `angle` (with the segment connectivity
// when given): basis points, trace mesh and synthetic matrices written to `directory`
// (calling rank), the library model entry returned.
json WriteCornerCoupon(const fs::path &directory, const std::string &tag, double angle,
                       std::optional<double> connectivity)
{
  const CornerTraceBasisRule rule;
  const auto seed = MakeCornerBoxSeed(kR, kT, kOE, true, rule);
  std::optional<double> connectivity_radians;
  if (connectivity)
  {
    connectivity_radians = *connectivity * kDeg;
  }
  const auto basis =
      BuildCornerTraceBasis(seed.points, seed.contour_groups, seed.zero_trace_indices,
                            angle * kDeg, true, rule, connectivity_radians);
  const fs::path points = directory / ("coupon-" + tag + "-points.csv");
  const fs::path vertices = directory / ("coupon-" + tag + "-vertices.csv");
  const fs::path triangles = directory / ("coupon-" + tag + "-triangles.csv");
  const fs::path domain = directory / ("coupon-" + tag + "-domain.csv");
  const fs::path surface = directory / ("coupon-" + tag + "-surface.csv");
  std::vector<int> zero;
  {
    std::ofstream output(points);
    output << std::setprecision(17) << "x,y,z\n";
    for (std::size_t k = 0; k < basis.knots.size(); k++)
    {
      output << basis.knots[k][0] << "," << basis.knots[k][1] << "," << basis.knots[k][2]
             << "\n";
      if (basis.zero[k])
      {
        zero.push_back(static_cast<int>(k) + 1);
      }
    }
  }
  {
    std::ofstream output(vertices);
    output << std::setprecision(17)
           << "vertex,x,y,z,basis,conductor,parent_a,parent_b,weight_a\n";
    for (std::size_t v = 0; v < basis.vertices.size(); v++)
    {
      const auto &vertex = basis.vertices[v];
      output << v + 1 << "," << vertex.point[0] << "," << vertex.point[1] << ","
             << vertex.point[2] << ",";
      if (vertex.basis >= 0)
      {
        output << vertex.basis + 1 << "," << (basis.zero[vertex.basis] ? 1 : 0)
               << ",0,0,0\n";
      }
      else
      {
        output << "0,0," << vertex.parent_a + 1 << "," << vertex.parent_b + 1 << ","
               << vertex.weight_a << "\n";
      }
    }
  }
  {
    std::ofstream output(triangles);
    output << "triangle,vertex_i,vertex_j,vertex_k\n";
    for (std::size_t t = 0; t < basis.triangles.size(); t++)
    {
      output << t + 1 << "," << basis.triangles[t][0] + 1 << ","
             << basis.triangles[t][1] + 1 << "," << basis.triangles[t][2] + 1 << "\n";
    }
  }
  {
    std::ofstream domain_output(domain), surface_output(surface);
    domain_output << "basis_i,basis_j,Q_ij (J)\n";
    surface_output << "interface,edge,R (m),basis_i,basis_j,Q_ij (J),Q_total_ij (J)\n";
    for (std::size_t i = 0; i < basis.knots.size(); i++)
    {
      for (std::size_t j = i; j < basis.knots.size(); j++)
      {
        const double value = (i == j ? 3.0 : 0.05 / (1.0 + double(j - i))) * 1.0e-12;
        domain_output << i + 1 << "," << j + 1 << "," << value << "\n";
        surface_output << "1,1," << kR << "," << i + 1 << "," << j + 1 << "," << value
                       << "," << value << "\n";
      }
    }
  }
  json trace_basis = {{"RingSize", rule.ring_size},
                      {"MetalInteriorKnots", rule.metal_interior_knots},
                      {"FreeKnots", rule.free_knots},
                      {"Fractions", rule.fractions}};
  if (connectivity)
  {
    trace_basis["ConnectivityAngleDegrees"] = *connectivity;
  }
  return {
      {"Name", "convex-corner-" + tag},
      {"Topology", "ConvexCorner"},
      {"Angle", angle},
      {"AngleDegrees", angle},
      {"Convexity", "Convex"},
      {"AngleTolerance", 1.0e-6},
      {"CornerRadius", 0.0},
      {"CornerRadiusTolerance", 0.0},
      {"FabricatedMatrix", domain.string()},
      {"ThinMatrix", domain.string()},
      {"FabricatedSurfaceMatrix", surface.string()},
      {"ThinSurfaceMatrix", surface.string()},
      {"BasisPoints", points.string()},
      {"TraceMesh", {{"Vertices", vertices.string()}, {"Triangles", triangles.string()}}},
      {"ContourGroups", seed.contour_groups},
      {"ZeroTraceIndices", zero},
      {"Interfaces", {{{"Type", "SA"}, {"Coupon", 1}}}},
      {"TraceBasis", trace_basis}};
}

// A hexahedral box with a square metal island (cracked boundary attribute 9 on the plane
// y = 0.5, [2, 6] x [2, 6] of the 8 x 8 (x, z) extent): four exact 90-degree convex
// corners.
std::unique_ptr<mfem::ParMesh> MakeSquareIslandMesh()
{
  constexpr double extent = 8.0;
  mfem::Mesh serial = mfem::Mesh::MakeCartesian3D(16, 4, 16, mfem::Element::HEXAHEDRON,
                                                  extent, 1.0, extent);
  for (int face = 0; face < serial.GetNumFaces(); face++)
  {
    int element1, element2;
    serial.GetFaceElements(face, &element1, &element2);
    if (element1 < 0 || element2 < 0)
    {
      continue;
    }
    mfem::Array<int> vertices;
    serial.GetFaceVertices(face, vertices);
    bool inside = true;
    for (const int vertex : vertices)
    {
      const double *point = serial.GetVertex(vertex);
      inside = inside && std::abs(point[1] - 0.5) < 1.0e-12 && point[0] >= 2.0 - 1.0e-12 &&
               point[0] <= 6.0 + 1.0e-12 && point[2] >= 2.0 - 1.0e-12 &&
               point[2] <= 6.0 + 1.0e-12;
    }
    if (inside)
    {
      serial.AddBdrElement(serial.GetFace(face)->Duplicate(&serial));
      serial.SetBdrAttribute(serial.GetNBE() - 1, 9);
    }
  }
  serial.FinalizeTopology();
  serial.Finalize();
  return std::make_unique<mfem::ParMesh>(Mpi::World(), serial);
}

}  // namespace

TEST_CASE("CornerFamilyLibraryLoad", "[cornerfamily][Serial][Parallel]")
{
  // The segment structure fails closed at LIBRARY LOAD (decision 137 (1);
  // ReadProcessLibrary through the preflight): a family whose segment spans the 135-degree
  // passage is refused before any corner is matched; a consistent family loads and, on a
  // square island, the exact 90-degree corners match the LEGACY tie coupon (its own
  // triangulation, the verified model) ahead of the per-side 90- / 90+ coupons (decision
  // 140 (1)).
  test::SharedTempDir temp;
  const fs::path isolated_points = temp.temp_dir / "isolated-points.csv";
  const fs::path isolated_domain = temp.temp_dir / "isolated-domain.csv";
  const fs::path isolated_surface = temp.temp_dir / "isolated-surface.csv";
  std::map<std::string, fs::path> libraries = {
      {"consistent", temp.temp_dir / "library-consistent.json"},
      {"across-passage", temp.temp_dir / "library-across-passage.json"},
      {"two-coupons", temp.temp_dir / "library-two-coupons.json"}};
  if (Mpi::Root(Mpi::World()))
  {
    {
      std::ofstream output(isolated_points);
      output << "x,y,z\n-0.16,-0.12,0.0\n0.16,-0.12,0.0\n0.16,0.12,0.0\n-0.16,0.12,0.0\n";
      std::ofstream domain(isolated_domain);
      domain << "basis_i,basis_j,Q_ij (J)\n";
      std::ofstream surface(isolated_surface);
      surface << "interface,edge,basis_i,basis_j,Q_total_ij (J)\n";
      for (int i = 1; i <= 4; i++)
      {
        for (int j = i; j <= 4; j++)
        {
          const double value = (i == j ? 2.0 : 0.2) * 1.0e-12;
          domain << i << "," << j << "," << value << "\n";
          surface << "1,1," << i << "," << j << "," << value << "\n";
        }
      }
    }
    const json base = {
        {"Version", 3},
        {"TraceLiftVersion", 2},
        {"MatchingRadius", kR},
        {"Fabrication",
         {{"InterfaceLayers", {{"SA", {{"Thickness", 0.002}, {"Permittivity", 4.0}}}}}}},
        {"Models",
         {{{"Name", "isolated"},
           {"Topology", "IsolatedEdge"},
           {"CouponDepth", kR},
           {"FabricatedMatrix", isolated_domain.string()},
           {"ThinMatrix", isolated_domain.string()},
           {"FabricatedSurfaceMatrix", isolated_surface.string()},
           {"ThinSurfaceMatrix", isolated_surface.string()},
           {"BasisPoints", isolated_points.string()},
           {"Interfaces", {{{"Type", "SA"}, {"Coupon", 1}}}}}}}};
    // The segment [90, 135] keyed 112.5 with the per-side 90- (82.5) and the legacy tie 90.
    const json legacy_90 = WriteCornerCoupon(temp.temp_dir, "90", 90.0, std::nullopt);
    const json lower_90 = WriteCornerCoupon(temp.temp_dir, "90-c82.5", 90.0, 82.5);
    const json upper_90 = WriteCornerCoupon(temp.temp_dir, "90-c112.5", 90.0, 112.5);
    const json node_105 = WriteCornerCoupon(temp.temp_dir, "105-c112.5", 105.0, 112.5);
    const json node_120 = WriteCornerCoupon(temp.temp_dir, "120-c112.5", 120.0, 112.5);
    const json node_135 = WriteCornerCoupon(temp.temp_dir, "135-c112.5", 135.0, 112.5);
    auto Write = [&](const std::string &name, const std::vector<json> &corners)
    {
      json library = base;
      library["Name"] = "unit-test-corner-family-" + name;
      for (const auto &corner : corners)
      {
        library["Models"].push_back(corner);
      }
      std::ofstream output(libraries.at(name));
      output << library.dump(2) << "\n";
    };
    Write("consistent", {upper_90, node_105, legacy_90, node_120, lower_90, node_135});
    // A 150-degree node stamped 112.5 (built at its own segment's key, the builder refuses
    // a keyed merge across the passage): the segment would span the 135-degree passage.
    json stamped_150 = WriteCornerCoupon(temp.temp_dir, "150-c144.2", 150.0, 144.2);
    stamped_150["TraceBasis"]["ConnectivityAngleDegrees"] = 112.5;
    Write("across-passage", {upper_90, node_105, node_120, node_135, stamped_150});
    // Two coupons at 105 degrees in one segment.
    json second_105 = node_105;
    second_105["Name"] = "convex-corner-105-c112.5-duplicate";
    Write("two-coupons", {upper_90, node_105, second_105, node_120, node_135});
  }
  Mpi::Barrier(Mpi::World());

  json config = {
      {"Problem", {{"Type", "Electrostatic"}, {"Output", temp.temp_dir.string()}}},
      {"Model", {{"Mesh", "unused.msh"}}},
      {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
      {"Boundaries",
       {{"Ground", {{"Attributes", {1, 2, 3, 4, 5, 6}}}},
        {"Terminal", {{{"Index", 1}, {"Attributes", {9}}}}},
        {"Postprocessing",
         {{"Dielectric",
           {{{"Index", 4},
             {"Attributes", {9}},
             {"Type", "SA"},
             {"Thickness", 0.002},
             {"Permittivity", 4.0},
             {"AutomaticEdges", true},
             {"EdgeDistances", {kR}},
             {"EdgeFrameNormal", {0.0, 1.0, 0.0}}}}}}}}},
      {"Solver",
       {{"Order", 1},
        {"Electrostatic",
         {{"ResponseCorrection",
           {{"Library", ""}, {"TargetInterfaces", {4}}, {"UnmatchedPolicy", "Warn"}}}}}}}};
  auto mesh = MakeSquareIslandMesh();
  auto Preflight = [&](const std::string &name)
  {
    auto library_config = config;
    library_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
        libraries.at(name).string();
    IoData iodata(library_config, false);
    iodata.boundaries.cracked_attributes.insert(9);
    const auto manifest_path = temp.temp_dir / ("requirements-" + name + ".json");
    WriteSurfaceResponseRequirements(iodata, *mesh, manifest_path.string());
    Mpi::Barrier(Mpi::World());
    std::ifstream input(manifest_path);
    REQUIRE(input);
    return json::parse(input);
  };
  // "Fabrication-process response library ... corner family" pins the LOAD-time check of
  // ReadProcessLibrary (the match-time precondition of SelectCornerFamilyStencil says
  // "Corner family segment structure: ..." instead).
  CHECK_THROWS_WITH(Preflight("across-passage"),
                    ContainsSubstring("Fabrication-process response library") &&
                        ContainsSubstring("corner family") &&
                        ContainsSubstring("node at 150") &&
                        ContainsSubstring("across a knot-corner passage"));
  CHECK_THROWS_WITH(Preflight("two-coupons"),
                    ContainsSubstring("Fabrication-process response library") &&
                        ContainsSubstring("corner family") &&
                        ContainsSubstring("two coupons at 105"));
  const json manifest = Preflight("consistent");
  int corners = 0;
  for (const auto &feature : manifest["Identification"]["Features"])
  {
    if (feature["Type"] != "ConvexCorner")
    {
      continue;
    }
    corners++;
    CHECK(feature["Match"]["Status"] == "Matched");
    CHECK(feature["Match"]["Model"] == "convex-corner-90");
  }
  CHECK(corners == 4);
}

}  // namespace palace
