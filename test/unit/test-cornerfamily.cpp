// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// The corner family's angle-interpolation rule (corner-qualification block 2026-09-29;
// design doc SURFACE-RESPONSE-IDENTIFICATION.md, Conventions CornerTraceBasisRule): the
// geometric events of the trace basis (a knot of a metal ring passing a fixed-layout vertex:
// the hats jump there), the segment connectivity that fixes the band triangulation over a
// segment, and the stencil rule that never straddles a knot-corner passage (legacy coupons
// exact only, per-side coupons at the passages, fail closed on a segment across a passage).
// Pinned to the same numbers as the Python mirror corner_family_interpolation.py
// (test_corner_family_interpolation.py).

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "models/cornertracebasis.hpp"

namespace palace
{

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

// The qualified family's node set of one convexity: (angle, connectivity angle) per segment.
std::vector<CornerFamilyNode> QualifiedFamily(bool convex)
{
  const double passage = convex ? 153.434948822922 : 158.198590513648;
  const std::vector<std::pair<double, double>> nodes =
      convex ? std::vector<std::pair<double, double>>{
                   {75.0, 82.5},     {80.0, 82.5},   {85.0, 82.5},   {90.0, 82.5},
                   {90.0, 112.5},    {105.0, 112.5}, {120.0, 112.5}, {135.0, 112.5},
                   {135.0, 144.2},   {144.0, 144.2}, {150.0, 144.2}, {passage, 144.2},
                   {passage, 166.7}, {159.0, 166.7}, {165.0, 166.7}, {180.0, 166.7}}
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
        CHECK(family[index].connectivity_angle_degrees == stencil.connectivity_angle_degrees);
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
  CHECK_THROWS(SelectCornerFamilyStencil({{120.0, 112.5, 0}, {150.0, 112.5, 1}}, 130.0, true,
                                         rule, tol));
  CHECK_THROWS(SelectCornerFamilyStencil({{120.0, 135.0, 0}, {130.0, 135.0, 1}}, 125.0, true,
                                         rule, tol));
  CHECK_THROWS(SelectCornerFamilyStencil({{105.0, 112.5, 0}, {105.0, 112.5, 1}, {120.0, 112.5, 2}},
                                         110.0, true, rule, tol));
  CHECK_THROWS(SelectCornerFamilyStencil(
      {{95.0, 100.0, 0}, {125.0, 100.0, 1}, {110.0, 115.0, 2}, {130.0, 115.0, 3}}, 120.0,
      true, rule, tol));
  // Lower orders where a segment has fewer nodes.
  const std::vector<CornerFamilyNode> shortened = {
      {135.0, 144.2, 0}, {150.0, 144.2, 1}, {153.434948822922, 144.2, 2}};
  CHECK(SelectCornerFamilyStencil(shortened, 142.5, true, rule, tol).rule == "quadratic");
  CHECK(SelectCornerFamilyStencil({shortened[0], shortened[1]}, 142.5, true, rule, tol).rule ==
        "linear");
}

}  // namespace palace
