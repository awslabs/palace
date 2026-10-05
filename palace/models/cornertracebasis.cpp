// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "cornertracebasis.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <mfem.hpp>

namespace palace
{

namespace
{

// Knot coincidence in perimeter fraction (corner-basis-fix review m3, 2026-09-29; the same
// value as KNOT_COINCIDENCE_FRACTION in generate_corner_response.py, where the rationale is
// recorded): a free or metal-interior knot within this fraction of a fixed-layout fraction
// k / RingSize snaps to it (SnapFraction) and a box corner within it of any knot gets no
// slave vertex, so the smallest triangle of a constructed basis has an edge of 8e-6 R, six
// orders above the surface mortar's degenerate-triangle threshold, while a snap moves a
// knot by at most 8e-6 R (15 pm at R = 1.9 um). Crossing knots are never moved.
constexpr double kKnotCoincidenceFraction = 1.0e-6;

double WrapFraction(double fraction)
{
  fraction = std::fmod(fraction, 1.0);
  if (fraction < 0.0)
  {
    fraction += 1.0;
  }
  return fraction == 1.0 ? 0.0 : fraction;
}

// The generator's snap_fraction: a fraction within kKnotCoincidenceFraction of k /
// ring_size is that fraction exactly.
double SnapFraction(double fraction, int ring_size)
{
  const double nearest = std::round(fraction * ring_size);
  if (std::abs(fraction - nearest / ring_size) <= kKnotCoincidenceFraction)
  {
    return WrapFraction(nearest / ring_size);
  }
  return fraction;
}

double FractionDistance(double a, double b)
{
  const double d = std::abs(WrapFraction(a) - WrapFraction(b));
  return std::min(d, 1.0 - d);
}

// Points of the fixed layout (the generator's square_ring): knots at the fractions
// k / size from (-half_width, 0) counterclockwise. SquarePerimeterPoint is the generator's
// square_perimeter_point formula; at the fixed fractions it reproduces square_ring's
// coordinates to floating-point rounding (~1e-15 R residues of either arithmetic; the
// generator's ring_points reuses square_ring's coordinates there for byte identity of the
// recorded files, MatchCornerFamily compares at 1e-9 R; the unit test
// SurfaceResponseOperatorCornerTraceBasis and the generator's test pin the eight points).
std::vector<std::array<double, 3>> FixedRingPoints(double half_width, double z, int size)
{
  std::vector<std::array<double, 3>> points;
  points.reserve(size);
  for (int k = 0; k < size; k++)
  {
    points.push_back(SquarePerimeterPoint(half_width, z, static_cast<double>(k) / size));
  }
  return points;
}

std::vector<std::string> RingRoleOrder(bool convex, const CornerTraceBasisRule &rule)
{
  std::vector<std::string> free, metal;
  for (int k = 1; k <= rule.free_knots; k++)
  {
    free.push_back("free" + std::to_string(k));
  }
  for (int m = 1; m <= rule.metal_interior_knots; m++)
  {
    metal.push_back("metal" + std::to_string(m));
  }
  std::vector<std::string> order;
  if (convex)
  {
    order.insert(order.end(), free.begin() + 1, free.end());
    order.push_back("crossing1");
    order.insert(order.end(), metal.begin(), metal.end());
    order.push_back("crossing2");
    order.push_back(free.front());
  }
  else
  {
    order.push_back("crossing1");
    order.insert(order.end(), free.begin(), free.end());
    order.push_back("crossing2");
    order.insert(order.end(), metal.begin(), metal.end());
  }
  return order;
}

// The generator's connect_rings_by_fraction: the band between two closed rings sharing the
// perimeter parameter (both start at fraction 0, counterclockwise), possibly with different
// vertices. Equal vertex sets give the lane-2 connect_rings triangles in the same order.
void ConnectRingsByFraction(std::vector<std::array<int, 3>> &triangles,
                            const std::vector<double> &first_fractions,
                            const std::vector<int> &first_vertices,
                            const std::vector<double> &second_fractions,
                            const std::vector<int> &second_vertices)
{
  const int first_count = static_cast<int>(first_fractions.size());
  const int second_count = static_cast<int>(second_fractions.size());
  constexpr double tolerance = 1.0e-12;
  constexpr double infinity = std::numeric_limits<double>::infinity();
  int i = 0, j = 0;
  while (i < first_count || j < second_count)
  {
    const double next_first =
        i < first_count ? (i + 1 < first_count ? first_fractions[i + 1] : 1.0) : infinity;
    const double next_second = j < second_count
                                   ? (j + 1 < second_count ? second_fractions[j + 1] : 1.0)
                                   : infinity;
    const int first_vertex = first_vertices[i % first_count];
    const int second_vertex = second_vertices[j % second_count];
    if (next_first < next_second - tolerance)
    {
      triangles.push_back(
          {first_vertex, first_vertices[(i + 1) % first_count], second_vertex});
      i++;
    }
    else if (next_second < next_first - tolerance)
    {
      triangles.push_back(
          {first_vertex, second_vertices[(j + 1) % second_count], second_vertex});
      j++;
    }
    else
    {
      const int next_first_vertex = first_vertices[(i + 1) % first_count];
      const int next_second_vertex = second_vertices[(j + 1) % second_count];
      triangles.push_back({first_vertex, next_first_vertex, next_second_vertex});
      triangles.push_back({first_vertex, next_second_vertex, second_vertex});
      i++;
      j++;
    }
  }
}

// Closed metal footprint of a sharp corner in the coupon frame (the generator's
// metal_footprint_mask for CornerRadius 0): the sector between the arms for a convex
// corner (its complement for a concave one), arms included.
bool OnMetalFootprint(const std::array<double, 3> &point, double angle_radians, bool convex,
                      double tolerance)
{
  const double first_distance = point[1];
  const double second_distance =
      point[0] * std::sin(angle_radians) - point[1] * std::cos(angle_radians);
  const bool in_wedge = first_distance >= -tolerance && second_distance >= -tolerance;
  if (convex)
  {
    return in_wedge;
  }
  const bool on_arm =
      (std::abs(first_distance) <= tolerance && second_distance >= -tolerance) ||
      (std::abs(second_distance) <= tolerance && first_distance >= -tolerance);
  return !in_wedge || on_arm;
}

// The generator's free_knot_fractions: the graded knots at free_knot_grading x R from both
// ends of the free arc (start, end) — scaled by FreeKnotGradingScale on a free arc shorter
// than the rule's reference — the remaining at equal fractions between the innermost
// graded ones (equal fractions of the whole arc without a grading).
std::vector<double> FreeKnotFractions(const std::pair<double, double> &free,
                                      const CornerTraceBasisRule &rule)
{
  std::vector<double> fractions;
  const auto [start, end] = free;
  if (rule.free_knot_grading.empty())
  {
    for (int k = 1; k <= rule.free_knots; k++)
    {
      fractions.push_back(start + (end - start) * k / (rule.free_knots + 1));
    }
    return fractions;
  }
  // A perimeter distance of g R is the fraction g / 8 (the perimeter is 8 R). The scale is
  // applied only below the reference so that every unscaled layout keeps its exact
  // floating-point expressions.
  const double scale = FreeKnotGradingScale(8.0 * (end - start), rule);
  auto Offset = [&](double g) { return scale < 1.0 ? (g / 8.0) * scale : g / 8.0; };
  const double outermost = Offset(rule.free_knot_grading.back());
  const double inner_start = start + outermost, inner_end = end - outermost;
  MFEM_VERIFY(inner_end - inner_start > kKnotCoincidenceFraction,
              "The corner trace basis free arc is too short for the free knot grading!");
  for (const double g : rule.free_knot_grading)
  {
    fractions.push_back(start + Offset(g));
    fractions.push_back(end - Offset(g));
  }
  const int remaining =
      rule.free_knots - 2 * static_cast<int>(rule.free_knot_grading.size());
  for (int k = 1; k <= remaining; k++)
  {
    fractions.push_back(inner_start + (inner_end - inner_start) * k / (remaining + 1));
  }
  std::sort(fractions.begin(), fractions.end());
  return fractions;
}

}  // namespace

double FreeKnotGradingScale(double free_arc_over_r, const CornerTraceBasisRule &rule)
{
  const double reference = rule.free_knot_grading_reference_free_arc_over_r;
  if (rule.free_knot_grading.empty() ||
      free_arc_over_r >= reference - kFreeKnotGradingReferenceToleranceOverR)
  {
    return 1.0;
  }
  return free_arc_over_r / reference;
}

CornerTraceBasisRule RefinedCornerTraceBasisRule()
{
  CornerTraceBasisRule rule;
  rule.ring_layout = CornerRingLayout::ALL_RINGS_FOLLOW_METAL;
  rule.metal_interior_knots = 5;
  rule.free_knots = 9;
  rule.ring_size = 16;
  rule.free_knot_grading = {1.0 / 3.0, 2.0 / 3.0};
  rule.extra_levels_above_over_overetch = {1.0, 4.0};
  rule.free_knot_grading_reference_free_arc_over_r = kFreeKnotGradingReferenceFreeArcOverR;
  return rule;
}

std::vector<double> CornerRuleLevels(const CornerTraceBasisRule &rule, double radius,
                                     double metal_thickness, double overetch_depth)
{
  std::set<double> level_set = {-radius,         -radius / 3.0, -overetch_depth, 0.0,
                                metal_thickness, radius / 3.0,  radius};
  for (const double k : rule.extra_levels_above_over_overetch)
  {
    level_set.insert(metal_thickness + k * overetch_depth);
  }
  return std::vector<double>(level_set.begin(), level_set.end());
}

std::string CheckCornerTraceBasisRule(const CornerTraceBasisRule &rule)
{
  if (rule.ring_size < 3 || rule.metal_interior_knots < 0 || rule.free_knots < 1 ||
      rule.ring_size != 2 + rule.metal_interior_knots + rule.free_knots)
  {
    return "RingSize must equal 2 + MetalInteriorKnots + FreeKnots";
  }
  if (rule.fractions != "PerimeterArcLength")
  {
    return "Fractions must be PerimeterArcLength";
  }
  if (!rule.extra_levels_above_over_overetch.empty())
  {
    if (!rule.AllRings())
    {
      return "ExtraLevelsAboveOverOveretch is an AllRingsFollowMetal layout option";
    }
    for (std::size_t k = 0; k < rule.extra_levels_above_over_overetch.size(); k++)
    {
      if (!(rule.extra_levels_above_over_overetch[k] > 0.0) ||
          (k > 0 && !(rule.extra_levels_above_over_overetch[k] >
                      rule.extra_levels_above_over_overetch[k - 1])))
      {
        return "ExtraLevelsAboveOverOveretch must be increasing positive multiples of "
               "OveretchDepth";
      }
    }
  }
  if (!rule.free_knot_grading.empty())
  {
    if (!rule.AllRings())
    {
      return "FreeKnotGrading is an AllRingsFollowMetal layout option";
    }
    if (2 * static_cast<int>(rule.free_knot_grading.size()) > rule.free_knots)
    {
      return "FreeKnotGrading places more knots than FreeKnots";
    }
    for (std::size_t g = 0; g < rule.free_knot_grading.size(); g++)
    {
      if (!(rule.free_knot_grading[g] > 0.0) ||
          (g > 0 && !(rule.free_knot_grading[g] > rule.free_knot_grading[g - 1])))
      {
        return "FreeKnotGrading must be increasing positive distances over R";
      }
    }
    // The unscaled grading needs a free arc longer than twice its outermost distance; the
    // reference (the arc at and above which the grading is unscaled) must lie above that.
    if (!(rule.free_knot_grading_reference_free_arc_over_r >
          2.0 * rule.free_knot_grading.back()))
    {
      return "FreeKnotGradingReferenceFreeArcOverR must exceed twice the outermost "
             "FreeKnotGrading distance";
    }
  }
  return "";
}

std::array<double, 3> SquarePerimeterPoint(double half_width, double z, double fraction)
{
  const double coordinate = WrapFraction(fraction) * 8.0 * half_width;
  if (coordinate < half_width)
  {
    return {-half_width, -coordinate, z};
  }
  if (coordinate < 3.0 * half_width)
  {
    return {-half_width + coordinate - half_width, -half_width, z};
  }
  if (coordinate < 5.0 * half_width)
  {
    return {half_width, -half_width + coordinate - 3.0 * half_width, z};
  }
  if (coordinate < 7.0 * half_width)
  {
    return {half_width - (coordinate - 5.0 * half_width), half_width, z};
  }
  return {-half_width, half_width - (coordinate - 7.0 * half_width), z};
}

double SquarePerimeterFraction(double half_width, const std::array<double, 3> &point)
{
  const double x = point[0], y = point[1];
  const double tolerance = 1.0e-9 * half_width;
  double coordinate;
  if (std::abs(x + half_width) <= tolerance && y <= 0.0)
  {
    coordinate = -y;
  }
  else if (std::abs(y + half_width) <= tolerance)
  {
    coordinate = half_width + (x + half_width);
  }
  else if (std::abs(x - half_width) <= tolerance)
  {
    coordinate = 3.0 * half_width + (y + half_width);
  }
  else if (std::abs(y - half_width) <= tolerance)
  {
    coordinate = 5.0 * half_width + (half_width - x);
  }
  else if (std::abs(x + half_width) <= tolerance)
  {
    coordinate = 7.0 * half_width + (half_width - y);
  }
  else
  {
    throw std::runtime_error("Corner trace basis point does not lie on the box perimeter");
  }
  return WrapFraction(coordinate / (8.0 * half_width));
}

std::pair<double, double> ArmCrossingFractions(double radius, double angle_radians)
{
  std::array<double, 2> fractions{};
  const std::array<std::array<double, 2>, 2> directions = {
      {{1.0, 0.0}, {std::cos(angle_radians), std::sin(angle_radians)}}};
  for (int arm = 0; arm < 2; arm++)
  {
    const auto &direction = directions[arm];
    const double scale = radius / std::max(std::abs(direction[0]), std::abs(direction[1]));
    std::array<double, 3> point = {scale * direction[0], scale * direction[1], 0.0};
    for (int d = 0; d < 2; d++)
    {
      if (std::abs(std::abs(point[d]) - radius) <= 1.0e-12 * radius)
      {
        point[d] = point[d] < 0.0 ? -radius : radius;
      }
    }
    fractions[arm] = SquarePerimeterFraction(radius, point);
  }
  double first = fractions[0], second = fractions[1];
  if (second <= first + kKnotCoincidenceFraction)
  {
    second += 1.0;
  }
  return {first, second};
}

std::vector<CornerRingVertex> CornerMetalRingLayout(double radius, double angle_radians,
                                                    bool convex,
                                                    const CornerTraceBasisRule &rule)
{
  return CornerRuleRingLayout(radius, angle_radians, convex, rule, true);
}

std::vector<CornerRingVertex> CornerRuleRingLayout(double radius, double angle_radians,
                                                   bool convex,
                                                   const CornerTraceBasisRule &rule,
                                                   bool pec)
{
  {
    const std::string reason = CheckCornerTraceBasisRule(rule);
    MFEM_VERIFY(reason.empty(), "Invalid corner trace basis rule: " << reason << "!");
  }
  const auto [first, second] = ArmCrossingFractions(radius, angle_radians);
  std::pair<double, double> metal, free;
  if (convex)
  {
    metal = {first, second};
    free = {second, first + 1.0};
  }
  else
  {
    metal = {second, first + 1.0};
    free = {first, second};
  }
  std::map<std::string, double> roles = {{"crossing1", first},
                                         {"crossing2", WrapFraction(second)}};
  for (int m = 1; m <= rule.metal_interior_knots; m++)
  {
    roles["metal" + std::to_string(m)] =
        SnapFraction(WrapFraction(metal.first + (metal.second - metal.first) * m /
                                                    (rule.metal_interior_knots + 1)),
                     rule.ring_size);
  }
  {
    const auto fractions = FreeKnotFractions(free, rule);
    for (int k = 1; k <= rule.free_knots; k++)
    {
      roles["free" + std::to_string(k)] =
          SnapFraction(WrapFraction(fractions[k - 1]), rule.ring_size);
    }
  }
  std::vector<CornerRingVertex> knots;
  const auto order = RingRoleOrder(convex, rule);
  for (int slot = 0; slot < static_cast<int>(order.size()); slot++)
  {
    const auto &role = order[slot];
    knots.push_back({roles.at(role),
                     (role.rfind("free", 0) == 0 || !pec) ? CornerRingVertex::Kind::FREE
                                                          : CornerRingVertex::Kind::ZERO,
                     slot});
  }
  std::sort(knots.begin(), knots.end(),
            [](const auto &a, const auto &b) { return a.fraction < b.fraction; });
  for (std::size_t k = 1; k < knots.size(); k++)
  {
    MFEM_VERIFY(knots[k].fraction - knots[k - 1].fraction > kKnotCoincidenceFraction,
                "Two knots of the corner trace basis rule coincide!");
  }
  std::vector<CornerRingVertex> vertices = knots;
  for (const double corner : {0.125, 0.375, 0.625, 0.875})
  {
    double separation = 1.0;
    for (const auto &knot : knots)
    {
      separation = std::min(separation, FractionDistance(corner, knot.fraction));
    }
    if (separation > kKnotCoincidenceFraction)
    {
      vertices.push_back({corner, CornerRingVertex::Kind::SLAVE, -1});
    }
  }
  std::sort(vertices.begin(), vertices.end(),
            [](const auto &a, const auto &b) { return a.fraction < b.fraction; });
  return vertices;
}

std::vector<int> CornerZeroSlots(bool convex, const CornerTraceBasisRule &rule)
{
  std::vector<int> slots;
  const auto order = RingRoleOrder(convex, rule);
  for (int slot = 0; slot < static_cast<int>(order.size()); slot++)
  {
    if (order[slot].rfind("free", 0) != 0)
    {
      slots.push_back(slot);
    }
  }
  return slots;
}

CornerBoxRings DescribeCornerBoxRings(const std::vector<std::array<double, 3>> &points,
                                      const std::vector<int> &contour_groups,
                                      const std::vector<int> &zero_trace_indices)
{
  CornerBoxRings box;
  MFEM_VERIFY(!contour_groups.empty(),
              "A corner coupon's trace basis needs ContourGroups (closed box rings)!");
  const std::set<int> zero(zero_trace_indices.begin(), zero_trace_indices.end());
  int offset = 0;
  for (const int size : contour_groups)
  {
    MFEM_VERIFY(offset + size <= static_cast<int>(points.size()),
                "Corner coupon ContourGroups exceed its BasisPoints!");
    CornerBoxRings::Ring ring;
    ring.offset = offset;
    ring.size = size;
    ring.z = points[offset][2];
    for (int i = 0; i < size; i++)
    {
      const auto &point = points[offset + i];
      MFEM_VERIFY(std::abs(point[2] - ring.z) <= 1.0e-9 * std::max(1.0, std::abs(ring.z)),
                  "A corner coupon box ring is not at one height!");
      ring.half_width =
          std::max(ring.half_width, std::max(std::abs(point[0]), std::abs(point[1])));
      ring.metal = ring.metal || zero.count(offset + i) > 0;
    }
    box.rings.push_back(ring);
    offset += size;
  }
  MFEM_VERIFY(offset == static_cast<int>(points.size()),
              "Corner coupon ContourGroups do not partition its BasisPoints!");
  box.radius = 0.0;
  for (const auto &ring : box.rings)
  {
    box.radius = std::max(box.radius, ring.half_width);
  }
  // The generator's ring order: the outer rings (half width R) by ascending height, then
  // the top inner cap ring (z = +R) and the bottom inner cap ring (z = -R).
  const double tolerance = 1.0e-9 * box.radius;
  box.outer_count = 0;
  while (box.outer_count < static_cast<int>(box.rings.size()) &&
         std::abs(box.rings[box.outer_count].half_width - box.radius) <= tolerance)
  {
    box.outer_count++;
  }
  MFEM_VERIFY(
      box.outer_count >= 2 && box.outer_count + 2 == static_cast<int>(box.rings.size()) &&
          std::abs(box.rings[box.outer_count].z - box.radius) <= tolerance &&
          std::abs(box.rings[box.outer_count + 1].z + box.radius) <= tolerance &&
          !box.rings[box.outer_count].metal && !box.rings[box.outer_count + 1].metal,
      "A corner coupon's box rings are not the generator's (outer rings, then the "
      "top and bottom inner cap rings)!");
  for (int r = 1; r < box.outer_count; r++)
  {
    MFEM_VERIFY(box.rings[r].z > box.rings[r - 1].z,
                "A corner coupon's outer box rings are not ordered by height!");
  }
  return box;
}

CornerBoxSeed MakeCornerBoxSeed(double radius, double metal_thickness,
                                double overetch_depth, bool convex,
                                const CornerTraceBasisRule &rule)
{
  MFEM_VERIFY(radius > 0.0 && metal_thickness > 0.0 && metal_thickness < radius / 3.0 &&
                  overetch_depth >= 0.0 && overetch_depth < radius / 3.0 &&
                  overetch_depth != metal_thickness,
              "Invalid corner box seed dimensions!");
  const std::vector<double> levels =
      CornerRuleLevels(rule, radius, metal_thickness, overetch_depth);
  CornerBoxSeed seed;
  const auto zero_slots = CornerZeroSlots(convex, rule);
  for (const double level : levels)
  {
    const int offset = static_cast<int>(seed.points.size());
    const auto ring = FixedRingPoints(radius, level, rule.ring_size);
    seed.points.insert(seed.points.end(), ring.begin(), ring.end());
    seed.contour_groups.push_back(rule.ring_size);
    if (level == 0.0 || level == metal_thickness)
    {
      for (const int slot : zero_slots)
      {
        seed.zero_trace_indices.push_back(offset + slot);
      }
    }
  }
  for (const double level : {radius, -radius})
  {
    const auto ring = FixedRingPoints(radius / 3.0, level, rule.ring_size);
    seed.points.insert(seed.points.end(), ring.begin(), ring.end());
    seed.contour_groups.push_back(rule.ring_size);
  }
  std::sort(seed.zero_trace_indices.begin(), seed.zero_trace_indices.end());
  return seed;
}

ConstructedCornerTraceBasis
BuildCornerTraceBasis(const std::vector<std::array<double, 3>> &node_points,
                      const std::vector<int> &contour_groups,
                      const std::vector<int> &zero_trace_indices, double angle_radians,
                      bool convex, const CornerTraceBasisRule &rule,
                      std::optional<double> connectivity_angle_radians)
{
  const CornerBoxRings box =
      DescribeCornerBoxRings(node_points, contour_groups, zero_trace_indices);
  const double radius = box.radius;
  const double tolerance = 1.0e-9 * radius;
  const bool all_rings = rule.AllRings();
  MFEM_VERIFY(!(all_rings && connectivity_angle_radians),
              "The AllRingsFollowMetal corner trace basis has no events and takes no "
              "connectivity angle!");
  ConstructedCornerTraceBasis basis;
  basis.knots.assign(node_points.size(), {});
  basis.zero.assign(node_points.size(), false);
  basis.contour_groups = contour_groups;
  const int basis_size = static_cast<int>(node_points.size());
  // Per ring: the vertices in MERGE order (key, vertex index); the key is the perimeter
  // fraction, or with a connectivity angle the role's fraction at that angle (slaves keyed
  // at their corner).
  std::vector<std::vector<double>> ring_fractions(box.rings.size());
  std::vector<std::vector<int>> ring_vertices(box.rings.size());
  for (std::size_t r = 0; r < box.rings.size(); r++)
  {
    const auto &ring = box.rings[r];
    MFEM_VERIFY(ring.size == rule.ring_size, "A corner coupon box ring has "
                                                 << ring.size
                                                 << " knots; the trace basis rule needs "
                                                 << rule.ring_size << "!");
    if (!ring.metal && !all_rings)
    {
      // The fixed layout: the node's own points, at the fractions k / RingSize.
      for (int k = 0; k < ring.size; k++)
      {
        const auto &point = node_points[ring.offset + k];
        const double fraction = SquarePerimeterFraction(ring.half_width, point);
        MFEM_VERIFY(FractionDistance(fraction, static_cast<double>(k) / ring.size) <=
                        kKnotCoincidenceFraction,
                    "A corner coupon ring that does not meet the metal is not the fixed "
                    "layout (knots at the fractions k / RingSize)!");
        basis.knots[ring.offset + k] = point;
        ring_fractions[r].push_back(static_cast<double>(k) / ring.size);
        ring_vertices[r].push_back(ring.offset + k);
      }
      continue;
    }
    MFEM_VERIFY(!ring.metal || std::abs(ring.half_width - radius) <= tolerance,
                "A corner coupon ring that meets the metal is not an outer box ring!");
    // The rule's layout (perimeter fractions are scale-free: the cap rings of half width
    // R / 3 share them); PEC knots on the rings the node's zero set marks.
    const auto layout =
        CornerRuleRingLayout(radius, angle_radians, convex, rule, ring.metal);
    // Merge keys (the generator's connectivity_keys).
    std::vector<double> keys(layout.size());
    for (std::size_t v = 0; v < layout.size(); v++)
    {
      keys[v] = layout[v].fraction;
    }
    if (connectivity_angle_radians)
    {
      const auto key_layout =
          CornerMetalRingLayout(radius, *connectivity_angle_radians, convex, rule);
      for (std::size_t v = 0; v < layout.size(); v++)
      {
        if (layout[v].kind == CornerRingVertex::Kind::SLAVE)
        {
          continue;
        }
        bool found = false;
        for (const auto &key_vertex : key_layout)
        {
          if (key_vertex.kind != CornerRingVertex::Kind::SLAVE &&
              key_vertex.slot == layout[v].slot)
          {
            keys[v] = key_vertex.fraction;
            found = true;
          }
        }
        MFEM_VERIFY(found, "Corner trace basis rule: a knot role has no key!");
      }
      std::vector<std::size_t> by_key(layout.size());
      for (std::size_t v = 0; v < by_key.size(); v++)
      {
        by_key[v] = v;
      }
      std::sort(by_key.begin(), by_key.end(),
                [&](std::size_t a, std::size_t b) { return keys[a] < keys[b]; });
      int descents = 0;
      for (std::size_t v = 0; v < by_key.size(); v++)
      {
        if (layout[by_key[(v + 1) % by_key.size()]].fraction < layout[by_key[v]].fraction)
        {
          descents++;
        }
      }
      MFEM_VERIFY(descents <= 1,
                  "Corner trace basis: the connectivity angle "
                      << *connectivity_angle_radians * 180.0 / std::acos(-1.0)
                      << " deg lies across a knot-corner passage from the angle "
                      << angle_radians * 180.0 / std::acos(-1.0)
                      << " deg (the band triangles would fold)!");
    }
    std::vector<int> knot_positions;
    std::vector<int> pending_slaves;
    std::vector<int> perimeter_vertices(layout.size(), -1);
    for (std::size_t v = 0; v < layout.size(); v++)
    {
      const auto &vertex = layout[v];
      if (vertex.kind == CornerRingVertex::Kind::SLAVE)
      {
        pending_slaves.push_back(static_cast<int>(v));
        continue;
      }
      const int index = ring.offset + vertex.slot;
      basis.knots[index] = SquarePerimeterPoint(ring.half_width, ring.z, vertex.fraction);
      basis.zero[index] = vertex.kind == CornerRingVertex::Kind::ZERO;
      perimeter_vertices[v] = index;
      knot_positions.push_back(static_cast<int>(v));
    }
    for (const int position : pending_slaves)
    {
      int previous = -1, following = -1;
      for (const int p : knot_positions)
      {
        if (p < position)
        {
          previous = p;
        }
        if (p > position && following < 0)
        {
          following = p;
        }
      }
      if (previous < 0)
      {
        previous = knot_positions.back();
      }
      if (following < 0)
      {
        following = knot_positions.front();
      }
      const double f_previous =
          layout[previous].fraction - (previous > position ? 1.0 : 0.0);
      const double f_following =
          layout[following].fraction + (following < position ? 1.0 : 0.0);
      const double weight_a =
          (f_following - layout[position].fraction) / (f_following - f_previous);
      ConstructedCornerTraceBasis::Vertex slave;
      slave.point =
          SquarePerimeterPoint(ring.half_width, ring.z, layout[position].fraction);
      slave.basis = -1;
      slave.parent_a = perimeter_vertices[previous];
      slave.parent_b = perimeter_vertices[following];
      slave.weight_a = weight_a;
      perimeter_vertices[position] = basis_size + static_cast<int>(basis.vertices.size());
      basis.vertices.push_back(slave);
    }
    // The ring in merge order (= the perimeter order without a connectivity angle).
    std::vector<std::size_t> order(layout.size());
    for (std::size_t v = 0; v < order.size(); v++)
    {
      order[v] = v;
    }
    std::stable_sort(order.begin(), order.end(),
                     [&](std::size_t a, std::size_t b) { return keys[a] < keys[b]; });
    for (const std::size_t v : order)
    {
      ring_fractions[r].push_back(keys[v]);
      ring_vertices[r].push_back(perimeter_vertices[v]);
    }
  }
  // Vertex list: the knots (basis order), then the slaves.
  std::vector<ConstructedCornerTraceBasis::Vertex> vertices;
  vertices.reserve(basis_size + basis.vertices.size());
  for (int k = 0; k < basis_size; k++)
  {
    ConstructedCornerTraceBasis::Vertex vertex;
    vertex.point = basis.knots[k];
    vertex.basis = k;
    vertices.push_back(vertex);
  }
  vertices.insert(vertices.end(), basis.vertices.begin(), basis.vertices.end());
  basis.vertices = std::move(vertices);

  // The zero set of the rule must be the node's zero set (the same semantics).
  std::vector<int> rule_zero;
  for (int k = 0; k < basis_size; k++)
  {
    if (basis.zero[k])
    {
      rule_zero.push_back(k);
    }
  }
  std::vector<int> node_zero(zero_trace_indices);
  std::sort(node_zero.begin(), node_zero.end());
  MFEM_VERIFY(rule_zero == node_zero,
              "The corner trace basis rule's zero set differs from the coupon's "
              "ZeroTraceIndices!");

  // Triangulation: consecutive outer rings, the top outer ring to the top inner cap ring,
  // the bottom inner cap ring to the bottom outer ring, and the cap fans (the generator's
  // build_surface).
  const int outer = box.outer_count;
  auto Connect = [&](int first, int second)
  {
    ConnectRingsByFraction(basis.triangles, ring_fractions[first], ring_vertices[first],
                           ring_fractions[second], ring_vertices[second]);
  };
  for (int r = 0; r + 1 < outer; r++)
  {
    Connect(r, r + 1);
  }
  Connect(outer - 1, outer);
  Connect(outer + 1, 0);
  const auto &top = ring_vertices[outer];
  const auto &bottom = ring_vertices[outer + 1];
  if (all_rings)
  {
    // Cap fans from centre slaves at the mean of the cap ring's two crossing-slot knots
    // (the generator's build_surface; a fan from a ring vertex is degenerate as soon as two
    // consecutive vertices share the apex's side).
    const auto order = RingRoleOrder(convex, rule);
    int crossing1 = -1, crossing2 = -1;
    for (int slot = 0; slot < static_cast<int>(order.size()); slot++)
    {
      if (order[slot] == "crossing1")
      {
        crossing1 = slot;
      }
      if (order[slot] == "crossing2")
      {
        crossing2 = slot;
      }
    }
    MFEM_VERIFY(crossing1 >= 0 && crossing2 >= 0,
                "Corner trace basis rule without crossings!");
    for (int cap = 0; cap < 2; cap++)
    {
      const auto &ring = box.rings[outer + cap];
      const auto &vertices = cap == 0 ? top : bottom;
      ConstructedCornerTraceBasis::Vertex centre;
      centre.point = {0.0, 0.0, ring.z};
      centre.basis = -1;
      centre.parent_a = ring.offset + crossing1;
      centre.parent_b = ring.offset + crossing2;
      centre.weight_a = 0.5;
      const int centre_index = static_cast<int>(basis.vertices.size());
      basis.vertices.push_back(centre);
      const int count = static_cast<int>(vertices.size());
      for (int i = 0; i < count; i++)
      {
        const int a = vertices[i], b = vertices[(i + 1) % count];
        if (cap == 0)
        {
          basis.triangles.push_back({a, b, centre_index});
        }
        else
        {
          basis.triangles.push_back({b, a, centre_index});
        }
      }
    }
    return basis;
  }
  for (int i = 1; i + 1 < rule.ring_size; i++)
  {
    basis.triangles.push_back({top[0], top[i], top[i + 1]});
    basis.triangles.push_back({bottom[i + 1], bottom[i], bottom[0]});
  }
  return basis;
}

std::string CheckCornerBasisCrossings(const std::vector<std::array<double, 3>> &points,
                                      const std::vector<int> &contour_groups,
                                      const std::vector<int> &zero_trace_indices,
                                      double angle_radians, bool convex, double tolerance)
{
  const std::set<int> zero(zero_trace_indices.begin(), zero_trace_indices.end());
  int offset = 0;
  auto Describe = [](const std::array<double, 3> &point)
  {
    return "(" + std::to_string(point[0]) + ", " + std::to_string(point[1]) + ", " +
           std::to_string(point[2]) + ")";
  };
  for (const int size : contour_groups)
  {
    bool metal = false;
    double half_width = 0.0;
    for (int i = 0; i < size; i++)
    {
      metal = metal || zero.count(offset + i) > 0;
      half_width = std::max(half_width, std::max(std::abs(points[offset + i][0]),
                                                 std::abs(points[offset + i][1])));
    }
    // The gate is about the box contour of the coupon frame: a square ring centred on the
    // apex (every knot on the perimeter of |x|, |y| <= R, the generator's matching box). A
    // ring with PEC knots that is not such a ring is not a corner-coupon box contour and is
    // not judged here.
    bool centred_box = true;
    for (int i = 0; i < size; i++)
    {
      centred_box = centred_box && std::abs(std::max(std::abs(points[offset + i][0]),
                                                     std::abs(points[offset + i][1])) -
                                            half_width) <= tolerance;
    }
    if (metal && centred_box)
    {
      const double z = points[offset][2];
      const auto [first, second] = ArmCrossingFractions(half_width, angle_radians);
      for (const double fraction : {first, second})
      {
        const auto crossing = SquarePerimeterPoint(half_width, z, fraction);
        bool found = false;
        for (int i = 0; i < size && !found; i++)
        {
          const auto &point = points[offset + i];
          const double distance =
              std::hypot(point[0] - crossing[0], point[1] - crossing[1]);
          found = distance <= tolerance && zero.count(offset + i) > 0;
        }
        if (!found)
        {
          return "the metal arm crosses the box ring at z = " + std::to_string(z) + " at " +
                 Describe(crossing) +
                 " where no PEC (ZeroTraceIndices) knot lies: a free trace hat has support "
                 "on the PEC part of the box contour (corner-family review 2026-09-29)";
        }
      }
      for (int i = 0; i < size; i++)
      {
        if (zero.count(offset + i) == 0 &&
            OnMetalFootprint(points[offset + i], angle_radians, convex, tolerance))
        {
          return "free trace knot " + std::to_string(offset + i + 1) + " at " +
                 Describe(points[offset + i]) + " lies on the metal footprint";
        }
      }
    }
    offset += size;
  }
  return "";
}

std::vector<CornerBasisEvent> CornerBasisEvents(bool convex,
                                                const CornerTraceBasisRule &rule)
{
  if (rule.AllRings())
  {
    // Every ring carries the same fractions: no fixed vertex is ever passed, and a knot
    // passing a box-corner slave on every ring at once leaves the interpolant continuous
    // (the slave's value is the linear interpolant of its neighbours on both rings and the
    // column diagonals keep one orientation). Measured on the option-(c) interpolant at the
    // MetalRingsOnly event angles: 1-3e-5, the smooth-angle level (corner-basis-refinement
    // 2026-09-30, hat_continuity_probe).
    return {};
  }
  MFEM_VERIFY(rule.free_knot_grading.empty(),
              "Corner basis events are defined for the equal-fraction free knots!");
  // Every knot's perimeter fraction is affine in the second crossing's unwrapped fraction
  // s2 in (first, first + 1] (first = 0.5, the +x arm): fraction = a s2 + b (see
  // CornerMetalRingLayout). A knot passes the fixed-layout vertex j / RingSize (+ n) at
  // s2 = (j / RingSize + n - b) / a; the second arm through that perimeter point gives the
  // corner angle (in (0, 180] for s2 in (first, first + 0.5]).
  const double first = 0.5;
  std::vector<std::pair<std::string, std::pair<double, double>>> roles;  // role, (a, b)
  roles.push_back({"crossing2", {1.0, 0.0}});
  for (int m = 1; m <= rule.metal_interior_knots; m++)
  {
    const double t = static_cast<double>(m) / (rule.metal_interior_knots + 1);
    if (convex)
    {
      roles.push_back({"metal" + std::to_string(m), {t, first * (1.0 - t)}});
    }
    else
    {
      roles.push_back({"metal" + std::to_string(m), {1.0 - t, (first + 1.0) * t}});
    }
  }
  for (int k = 1; k <= rule.free_knots; k++)
  {
    const double t = static_cast<double>(k) / (rule.free_knots + 1);
    if (convex)
    {
      roles.push_back({"free" + std::to_string(k), {1.0 - t, (first + 1.0) * t}});
    }
    else
    {
      roles.push_back({"free" + std::to_string(k), {t, first * (1.0 - t)}});
    }
  }
  std::vector<CornerBasisEvent> events;
  for (const auto &[role, coefficients] : roles)
  {
    const auto [a, b] = coefficients;
    for (int j = 0; j < rule.ring_size; j++)
    {
      const double c = static_cast<double>(j) / rule.ring_size;
      for (int n = 0; n <= 2; n++)
      {
        const double s2 = (c + n - b) / a;
        if (s2 <= first + kKnotCoincidenceFraction || s2 > first + 0.5 + 1.0e-12)
        {
          continue;
        }
        const auto point = SquarePerimeterPoint(1.0, 0.0, s2);
        const double angle = std::atan2(point[1], point[0]) * 180.0 / std::acos(-1.0);
        CornerBasisEvent event;
        event.angle_degrees = angle <= 0.0 ? angle + 360.0 : angle;
        event.role = role;
        event.fraction = c;
        // The box corners of the fixed layout lie at the odd multiples of 1 / (2 x 4).
        event.corner = std::abs(std::fmod(c * 4.0 + 0.5, 1.0) - 0.0) <= 1.0e-12 ||
                       std::abs(std::fmod(c * 4.0 + 0.5, 1.0) - 1.0) <= 1.0e-12;
        events.push_back(event);
      }
    }
  }
  std::sort(events.begin(), events.end(),
            [](const auto &x, const auto &y)
            {
              return x.angle_degrees != y.angle_degrees ? x.angle_degrees < y.angle_degrees
                                                        : x.role < y.role;
            });
  return events;
}

namespace
{

// The segment structure of a corner family: the knot-corner passages (segment boundaries)
// and the coupons grouped by connectivity angle (sorted by angle), legacy coupons apart.
struct CornerFamilySegments
{
  std::vector<double> boundaries;
  std::map<double, std::vector<const CornerFamilyNode *>> segments;
  bool legacy = false;
};

CornerFamilySegments GroupCornerFamilySegments(const std::vector<CornerFamilyNode> &nodes,
                                               bool convex,
                                               const CornerTraceBasisRule &rule,
                                               double tolerance)
{
  CornerFamilySegments grouped;
  for (const auto &event : CornerBasisEvents(convex, rule))
  {
    if (event.corner && (grouped.boundaries.empty() ||
                         event.angle_degrees - grouped.boundaries.back() > tolerance))
    {
      grouped.boundaries.push_back(event.angle_degrees);
    }
  }
  if (rule.AllRings())
  {
    // No events: the whole family is ONE segment (no connectivity angle; a coupon carrying
    // one is refused at library load).
    std::vector<const CornerFamilyNode *> members;
    for (const auto &node : nodes)
    {
      MFEM_VERIFY(!node.connectivity_angle_degrees,
                  "An AllRingsFollowMetal corner coupon carries a segment connectivity "
                  "angle (the layout has no events)!");
      members.push_back(&node);
    }
    std::sort(members.begin(), members.end(), [](const auto *a, const auto *b)
              { return a->angle_degrees < b->angle_degrees; });
    grouped.segments[0.0] = std::move(members);
    return grouped;
  }
  for (const auto &node : nodes)
  {
    if (!node.connectivity_angle_degrees)
    {
      grouped.legacy = true;
      continue;
    }
    bool merged = false;
    for (auto &[key, members] : grouped.segments)
    {
      if (std::abs(key - *node.connectivity_angle_degrees) <= tolerance)
      {
        members.push_back(&node);
        merged = true;
        break;
      }
    }
    if (!merged)
    {
      grouped.segments[*node.connectivity_angle_degrees] = {&node};
    }
  }
  for (auto &[key, members] : grouped.segments)
  {
    std::sort(members.begin(), members.end(), [](const auto *a, const auto *b)
              { return a->angle_degrees < b->angle_degrees; });
  }
  return grouped;
}

// The structural check of the segments (the reason, empty when consistent): every segment
// in one event-free interval with its connectivity angle, one coupon per angle in a
// segment, segments overlapping at shared node angles only.
std::string CornerFamilySegmentsReason(const CornerFamilySegments &grouped,
                                       double tolerance)
{
  const auto &boundaries = grouped.boundaries;
  auto Interval = [&](double angle)
  {
    // The event-free open interval containing `angle` (an angle within the tolerance of an
    // event belongs to the interval below and above: index of the first boundary above).
    std::size_t i = 0;
    while (i < boundaries.size() && boundaries[i] < angle - tolerance)
    {
      i++;
    }
    return i;
  };
  auto OnBoundary = [&](double angle)
  {
    for (const double boundary : boundaries)
    {
      if (std::abs(angle - boundary) <= tolerance)
      {
        return true;
      }
    }
    return false;
  };
  std::ostringstream text;
  std::vector<std::pair<double, double>> ranges;
  for (const auto &[key, members] : grouped.segments)
  {
    if (OnBoundary(key))
    {
      text << "corner family segment connectivity angle " << key
           << " deg lies on a knot-corner passage of the trace basis";
      return text.str();
    }
    const std::size_t interval = Interval(key);
    for (const auto *member : members)
    {
      if (!(std::abs(member->angle_degrees - key) <= tolerance ||
            Interval(member->angle_degrees) == interval ||
            (OnBoundary(member->angle_degrees) &&
             (Interval(member->angle_degrees) == interval ||
              Interval(member->angle_degrees) + 1 == interval))))
      {
        text << "corner family segment with connectivity angle " << key
             << " deg has a node at " << member->angle_degrees
             << " deg across a knot-corner passage of the trace basis (its band "
                "triangulation jumps there): the segment must be split";
        return text.str();
      }
    }
    for (std::size_t i = 1; i < members.size(); i++)
    {
      if (!(members[i]->angle_degrees - members[i - 1]->angle_degrees > tolerance))
      {
        text << "corner family segment with connectivity angle " << key
             << " deg has two coupons at " << members[i]->angle_degrees << " deg";
        return text.str();
      }
    }
    ranges.push_back({members.front()->angle_degrees, members.back()->angle_degrees});
  }
  std::sort(ranges.begin(), ranges.end());
  for (std::size_t i = 1; i < ranges.size(); i++)
  {
    if (!(ranges[i].first >= ranges[i - 1].second - tolerance))
    {
      text << "corner family segments [" << ranges[i - 1].first << ", "
           << ranges[i - 1].second << "] and [" << ranges[i].first << ", "
           << ranges[i].second << "] deg overlap beyond a shared node angle";
      return text.str();
    }
  }
  return "";
}

}  // namespace

std::string CheckCornerFamilySegments(const std::vector<CornerFamilyNode> &nodes,
                                      bool convex, const CornerTraceBasisRule &rule,
                                      double angle_tolerance_degrees)
{
  return CornerFamilySegmentsReason(
      GroupCornerFamilySegments(nodes, convex, rule, angle_tolerance_degrees),
      angle_tolerance_degrees);
}

CornerFamilyStencil SelectCornerFamilyStencil(const std::vector<CornerFamilyNode> &nodes,
                                              double angle_degrees, bool convex,
                                              const CornerTraceBasisRule &rule,
                                              double angle_tolerance_degrees)
{
  CornerFamilyStencil stencil;
  MFEM_VERIFY(!nodes.empty(), "A corner family stencil needs nodes!");
  const double tolerance = angle_tolerance_degrees;
  stencil.min_angle_degrees = std::numeric_limits<double>::infinity();
  stencil.max_angle_degrees = -std::numeric_limits<double>::infinity();
  for (const auto &node : nodes)
  {
    stencil.min_angle_degrees = std::min(stencil.min_angle_degrees, node.angle_degrees);
    stencil.max_angle_degrees = std::max(stencil.max_angle_degrees, node.angle_degrees);
  }
  // Node tolerance (qualification review m1): a node is EXACT when |node - angle| <= tol
  // and the angle is strictly inside a segment when both end distances are > tol, computed
  // as the same floating-point differences, so an angle at the tolerance boundary of a
  // node is exact or interior, never in no segment.
  auto Beyond = [tolerance](double first, double second)
  { return std::abs(first - second) > tolerance; };
  auto Weights = [&](const std::vector<const CornerFamilyNode *> &window)
  {
    stencil.nodes.clear();
    double best = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < window.size(); i++)
    {
      double weight = 1.0;
      for (std::size_t j = 0; j < window.size(); j++)
      {
        if (j != i)
        {
          weight *= (angle_degrees - window[j]->angle_degrees) /
                    (window[i]->angle_degrees - window[j]->angle_degrees);
        }
      }
      stencil.nodes.push_back({window[i]->index, weight});
      const double distance = std::abs(window[i]->angle_degrees - angle_degrees);
      if (distance < best)
      {
        best = distance;
        stencil.base = window[i]->index;
      }
    }
    stencil.rule =
        window.size() == 4 ? "cubic" : (window.size() == 3 ? "quadratic" : "linear");
  };

  // Exact node: the legacy coupon at that angle (its own tie triangulation) if the family
  // has one, else the coupon of the lower-angle segment (deterministic).
  {
    const CornerFamilyNode *exact = nullptr;
    for (const auto &node : nodes)
    {
      if (Beyond(node.angle_degrees, angle_degrees))
      {
        continue;
      }
      if (!exact ||
          (!node.connectivity_angle_degrees && exact->connectivity_angle_degrees) ||
          (node.connectivity_angle_degrees && exact->connectivity_angle_degrees &&
           *node.connectivity_angle_degrees < *exact->connectivity_angle_degrees))
      {
        exact = &node;
      }
    }
    if (exact)
    {
      stencil.nodes = {{exact->index, 1.0}};
      stencil.rule = "exact";
      stencil.base = exact->index;
      stencil.connectivity_angle_degrees = exact->connectivity_angle_degrees;
      return stencil;
    }
  }
  if (angle_degrees < stencil.min_angle_degrees)
  {
    std::ostringstream text;
    text << "corner angle " << angle_degrees << " deg is sharper than the smallest "
         << (convex ? "convex" : "concave") << " coupon angle " << stencil.min_angle_degrees
         << " deg (no extrapolation)";
    stencil.reason = text.str();
    return stencil;
  }
  if (angle_degrees > stencil.max_angle_degrees)
  {
    std::ostringstream text;
    text << "corner angle " << angle_degrees << " deg is wider than the widest "
         << (convex ? "convex" : "concave") << " coupon angle " << stencil.max_angle_degrees
         << " deg (no extrapolation; the first-order regime needs the straight anchor)";
    stencil.reason = text.str();
    return stencil;
  }
  // Segments: the coupons sharing a connectivity angle, each in one event-free interval
  // with its connectivity angle, ordered and overlapping at shared node angles only. The
  // structure is verified at library load (ReadProcessLibrary, CheckCornerFamilySegments;
  // decision 137 (1)); the precondition is asserted here.
  const auto grouped = GroupCornerFamilySegments(nodes, convex, rule, tolerance);
  {
    const std::string reason = CornerFamilySegmentsReason(grouped, tolerance);
    MFEM_VERIFY(reason.empty(), "Corner family segment structure: " << reason << "!");
  }
  const std::vector<const CornerFamilyNode *> *segment = nullptr;
  double segment_key = 0.0;
  for (const auto &[key, members] : grouped.segments)
  {
    if (angle_degrees > members.front()->angle_degrees &&
        Beyond(members.front()->angle_degrees, angle_degrees) &&
        angle_degrees < members.back()->angle_degrees &&
        Beyond(members.back()->angle_degrees, angle_degrees))
    {
      segment = &members;
      segment_key = key;
    }
  }
  if (!segment)
  {
    std::ostringstream text;
    if (grouped.legacy)
    {
      text << "corner angle " << angle_degrees
           << " deg needs interpolation but the corner family is built without segment "
              "connectivity records (TraceBasis ConnectivityAngleDegrees): its trace bases "
              "jump at every knot passage of a box vertex; rebuild the family on segment "
              "connectivity (corner-qualification block 2026-09-29)";
    }
    else
    {
      text
          << "corner angle " << angle_degrees
          << " deg lies in no segment of the corner family (segments end at the "
             "knot-corner "
             "passages of the trace basis: a node is needed on each side of every passage)";
    }
    stencil.reason = text.str();
    return stencil;
  }
  std::size_t interval = 0;
  while (interval + 1 < segment->size() &&
         (*segment)[interval + 1]->angle_degrees < angle_degrees)
  {
    interval++;
  }
  std::size_t begin = 0, end = segment->size();
  if (segment->size() > 4)
  {
    begin = std::min(interval > 0 ? interval - 1 : 0, segment->size() - 4);
    end = begin + 4;
  }
  Weights(std::vector<const CornerFamilyNode *>(segment->begin() + begin,
                                                segment->begin() + end));
  if (!rule.AllRings())
  {
    stencil.connectivity_angle_degrees = segment_key;
  }
  return stencil;
}

}  // namespace palace
