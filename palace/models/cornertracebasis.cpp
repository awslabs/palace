// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "cornertracebasis.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <mfem.hpp>

namespace palace
{

namespace
{

constexpr double kKnotCoincidenceFraction = 1.0e-9;

double WrapFraction(double fraction)
{
  fraction = std::fmod(fraction, 1.0);
  if (fraction < 0.0)
  {
    fraction += 1.0;
  }
  return fraction == 1.0 ? 0.0 : fraction;
}

double FractionDistance(double a, double b)
{
  const double d = std::abs(WrapFraction(a) - WrapFraction(b));
  return std::min(d, 1.0 - d);
}

// Points of the fixed layout (the generator's square_ring): knots at the fractions
// k / size from (-half_width, 0) counterclockwise.
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
    const double next_second =
        j < second_count ? (j + 1 < second_count ? second_fractions[j + 1] : 1.0)
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
  const bool in_wedge =
      first_distance >= -tolerance && second_distance >= -tolerance;
  if (convex)
  {
    return in_wedge;
  }
  const bool on_arm = (std::abs(first_distance) <= tolerance && second_distance >= 0.0) ||
                      (std::abs(second_distance) <= tolerance && first_distance >= 0.0);
  return !in_wedge || on_arm;
}

}  // namespace

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
    const double scale =
        radius / std::max(std::abs(direction[0]), std::abs(direction[1]));
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
  MFEM_VERIFY(rule.ring_size == 2 + rule.metal_interior_knots + rule.free_knots,
              "The corner trace basis rule needs RingSize = 2 crossings + "
              "MetalInteriorKnots + FreeKnots!");
  MFEM_VERIFY(rule.fractions == "PerimeterArcLength",
              "Unsupported corner trace basis fraction parametrisation \"" << rule.fractions
                                                                          << "\"!");
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
    roles["metal" + std::to_string(m)] = WrapFraction(
        metal.first + (metal.second - metal.first) * m / (rule.metal_interior_knots + 1));
  }
  for (int k = 1; k <= rule.free_knots; k++)
  {
    roles["free" + std::to_string(k)] =
        WrapFraction(free.first + (free.second - free.first) * k / (rule.free_knots + 1));
  }
  std::vector<CornerRingVertex> knots;
  const auto order = RingRoleOrder(convex, rule);
  for (int slot = 0; slot < static_cast<int>(order.size()); slot++)
  {
    const auto &role = order[slot];
    knots.push_back({roles.at(role), role.rfind("free", 0) == 0 ? CornerRingVertex::Kind::FREE
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
  // The generator's ring order: the outer rings (half width R) by ascending height, then the
  // top inner cap ring (z = +R) and the bottom inner cap ring (z = -R).
  const double tolerance = 1.0e-9 * box.radius;
  box.outer_count = 0;
  while (box.outer_count < static_cast<int>(box.rings.size()) &&
         std::abs(box.rings[box.outer_count].half_width - box.radius) <= tolerance)
  {
    box.outer_count++;
  }
  MFEM_VERIFY(box.outer_count >= 2 &&
                  box.outer_count + 2 == static_cast<int>(box.rings.size()) &&
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

CornerBoxSeed MakeCornerBoxSeed(double radius, double metal_thickness, double overetch_depth,
                                bool convex, const CornerTraceBasisRule &rule)
{
  MFEM_VERIFY(radius > 0.0 && metal_thickness > 0.0 && metal_thickness < radius / 3.0 &&
                  overetch_depth >= 0.0 && overetch_depth < radius / 3.0 &&
                  overetch_depth != metal_thickness,
              "Invalid corner box seed dimensions!");
  std::set<double> level_set = {-radius,         -radius / 3.0,   -overetch_depth, 0.0,
                                metal_thickness, radius / 3.0,    radius};
  std::vector<double> levels(level_set.begin(), level_set.end());
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

ConstructedCornerTraceBasis BuildCornerTraceBasis(
    const std::vector<std::array<double, 3>> &node_points,
    const std::vector<int> &contour_groups, const std::vector<int> &zero_trace_indices,
    double angle_radians, bool convex, const CornerTraceBasisRule &rule)
{
  const CornerBoxRings box =
      DescribeCornerBoxRings(node_points, contour_groups, zero_trace_indices);
  const double radius = box.radius;
  const double tolerance = 1.0e-9 * radius;
  ConstructedCornerTraceBasis basis;
  basis.knots.assign(node_points.size(), {});
  basis.zero.assign(node_points.size(), false);
  basis.contour_groups = contour_groups;
  const int basis_size = static_cast<int>(node_points.size());
  // Per ring: the vertices in perimeter order (fraction, vertex index).
  std::vector<std::vector<double>> ring_fractions(box.rings.size());
  std::vector<std::vector<int>> ring_vertices(box.rings.size());
  for (std::size_t r = 0; r < box.rings.size(); r++)
  {
    const auto &ring = box.rings[r];
    MFEM_VERIFY(ring.size == rule.ring_size,
                "A corner coupon box ring has " << ring.size
                                                << " knots; the trace basis rule needs "
                                                << rule.ring_size << "!");
    if (!ring.metal)
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
    MFEM_VERIFY(std::abs(ring.half_width - radius) <= tolerance,
                "A corner coupon ring that meets the metal is not an outer box ring!");
    const auto layout = CornerMetalRingLayout(radius, angle_radians, convex, rule);
    std::vector<int> knot_positions;
    std::vector<int> pending_slaves;
    for (std::size_t v = 0; v < layout.size(); v++)
    {
      const auto &vertex = layout[v];
      ring_fractions[r].push_back(vertex.fraction);
      if (vertex.kind == CornerRingVertex::Kind::SLAVE)
      {
        ring_vertices[r].push_back(-1);  // resolved below
        pending_slaves.push_back(static_cast<int>(v));
        continue;
      }
      const int index = ring.offset + vertex.slot;
      basis.knots[index] = SquarePerimeterPoint(radius, ring.z, vertex.fraction);
      basis.zero[index] = vertex.kind == CornerRingVertex::Kind::ZERO;
      ring_vertices[r].push_back(index);
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
      slave.point = SquarePerimeterPoint(radius, ring.z, layout[position].fraction);
      slave.basis = -1;
      slave.parent_a = ring_vertices[r][previous];
      slave.parent_b = ring_vertices[r][following];
      slave.weight_a = weight_a;
      ring_vertices[r][position] = basis_size + static_cast<int>(basis.vertices.size());
      basis.vertices.push_back(slave);
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
  double radius = 0.0;
  for (const auto &point : points)
  {
    radius = std::max(radius, std::max(std::abs(point[0]), std::abs(point[1])));
  }
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
      centred_box = centred_box &&
                    std::abs(std::max(std::abs(points[offset + i][0]),
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
          const double distance = std::hypot(point[0] - crossing[0], point[1] - crossing[1]);
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
  (void)radius;
  return "";
}

}  // namespace palace
