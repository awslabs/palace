// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_CORNER_TRACE_BASIS_HPP
#define PALACE_MODELS_CORNER_TRACE_BASIS_HPP

#include <array>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace palace
{

// Trace basis rule of the angle-interpolated corner family (supervisor decision on the
// corner-family review, 2026-09-29; the Python reference is
// examples/cpw3d_surface/corner_coupon/generate_corner_response.py, record "TraceBasis" of
// every corner coupon). The matching box |x|, |y| <= R of a corner coupon carries closed
// trace rings at several heights; every ring that meets the metal (the sheet / slab foot at
// z = 0 and the slab top at z = MetalThickness) has knots whose SEMANTICS do not depend on
// the corner angle: the two crossings of the metal arms with the ring (PEC, in the zero
// set), MetalInteriorKnots knots at equal perimeter-arc-length fractions of the metal arc
// between them (PEC) and FreeKnots knots at equal fractions of the free arc (free). Every
// node of the family therefore has the same knot count, the same zero set and like-to-like
// free knots at the same basis indices whose POSITIONS vary smoothly with the angle, so the
// entrywise Lagrange blend of the nodes' matrices is well posed. The knots are ordered by
// role within the ring (convex: free 2 .. free 5, first crossing, metal interior, second
// crossing, free 1 — the lane-2 order of the 90-degree node; concave: first crossing, free
// 1 .. free 5, second crossing, metal interior). Box corners that are no knot are SLAVE
// vertices of the trace triangulation (the linear interpolation in the perimeter fraction
// between the neighbouring knots) so the trace surface lies on the box. Rings that do not
// meet the metal keep the fixed layout (knots at the fractions k / RingSize, the corners
// among them). Fractions are perimeter arc length from (-R, 0) counterclockwise (down the
// left side first), the traversal of the generator's square_ring.
struct CornerTraceBasisRule
{
  int ring_size = 8;
  int metal_interior_knots = 1;
  int free_knots = 5;
  std::string fractions = "PerimeterArcLength";

  bool operator==(const CornerTraceBasisRule &other) const
  {
    return ring_size == other.ring_size &&
           metal_interior_knots == other.metal_interior_knots &&
           free_knots == other.free_knots && fractions == other.fractions;
  }
  bool operator!=(const CornerTraceBasisRule &other) const { return !(*this == other); }
};

// A vertex of a ring that meets the metal: its perimeter fraction in [0, 1), whether it is a
// PEC knot, a free knot or a slave (box corner), and for a knot its basis slot in the ring.
struct CornerRingVertex
{
  enum class Kind : char
  {
    ZERO,
    FREE,
    SLAVE
  };
  double fraction = 0.0;
  Kind kind = Kind::FREE;
  int slot = -1;
};

// Point of the square |x|, |y| <= half_width at the perimeter fraction (see above).
std::array<double, 3> SquarePerimeterPoint(double half_width, double z, double fraction);

// Perimeter fraction of a point on the square's perimeter (throws when it is not).
double SquarePerimeterFraction(double half_width, const std::array<double, 3> &point);

// Perimeter fractions where the two metal arms of a corner (the first along +x, the second
// at the angle counterclockwise; the straight anchor's second arm is -x) meet the box:
// (first, second) with first < second <= first + 1.
std::pair<double, double> ArmCrossingFractions(double radius, double angle_radians);

// The knots and slave vertices of a ring that meets the metal, in perimeter order.
std::vector<CornerRingVertex> CornerMetalRingLayout(double radius, double angle_radians,
                                                    bool convex,
                                                    const CornerTraceBasisRule &rule);

// The zero-set slots (basis positions within a ring that meets the metal).
std::vector<int> CornerZeroSlots(bool convex, const CornerTraceBasisRule &rule);

// A trace basis constructed by the rule: the knots (basis points), their PEC flags, the
// trace triangulation with the slave vertices appended after the knots (vertex indices
// 0-based; slave parents 0-based) and the per-ring knot counts.
struct ConstructedCornerTraceBasis
{
  struct Vertex
  {
    std::array<double, 3> point{};
    int basis = -1;  // 0-based knot index, -1 for a slave
    int parent_a = -1;
    int parent_b = -1;
    double weight_a = 0.0;
  };
  std::vector<std::array<double, 3>> knots;
  std::vector<bool> zero;
  std::vector<Vertex> vertices;
  std::vector<std::array<int, 3>> triangles;
  std::vector<int> contour_groups;
};

// The description of a family node's box rings taken from its BasisPoints / ContourGroups /
// ZeroTraceIndices: the rings' half widths and heights, which rings meet the metal, and the
// ring adjacency of the generator (outer rings by ascending height, then the top inner cap
// ring at z = +R, then the bottom inner cap ring at z = -R). Throws when the node's basis
// does not have that structure.
struct CornerBoxRings
{
  struct Ring
  {
    int offset = 0;
    int size = 0;
    double half_width = 0.0;
    double z = 0.0;
    bool metal = false;
  };
  double radius = 0.0;
  std::vector<Ring> rings;  // the generator's order
  int outer_count = 0;      // rings 0 .. outer_count - 1 are the outer rings
};

CornerBoxRings DescribeCornerBoxRings(const std::vector<std::array<double, 3>> &points,
                                      const std::vector<int> &contour_groups,
                                      const std::vector<int> &zero_trace_indices);

// The generator's box rings before the rule is applied: the fixed layout on every ring
// (outer rings at z = -R, -R / 3, -OveretchDepth, 0, MetalThickness, R / 3, R, then the top
// and bottom inner cap rings of half width R / 3) with the rule's zero set on the two rings
// that meet the metal. BuildCornerTraceBasis on this seed gives the rule's basis at any
// angle (unit tests; the Python generator is the reference for coupon files).
struct CornerBoxSeed
{
  std::vector<std::array<double, 3>> points;
  std::vector<int> contour_groups;
  std::vector<int> zero_trace_indices;  // 0-based
};

CornerBoxSeed MakeCornerBoxSeed(double radius, double metal_thickness, double overetch_depth,
                                bool convex, const CornerTraceBasisRule &rule);

// The trace basis of the family at any corner angle: the node's fixed rings (its own
// points) and the metal rings laid out by the rule at `angle_radians`, triangulated by
// perimeter fraction (equal vertex sets reproduce the generator's connect_rings). Used for
// the runtime model of an interpolated corner (the nodes' knot semantics are the same; the
// positions at the device angle are the rule's) and, at a node's own angle, to check a
// coupon's files against the rule.
ConstructedCornerTraceBasis BuildCornerTraceBasis(
    const std::vector<std::array<double, 3>> &node_points,
    const std::vector<int> &contour_groups, const std::vector<int> &zero_trace_indices,
    double angle_radians, bool convex, const CornerTraceBasisRule &rule);

// The basis gate of a corner coupon (corner-family review 2026-09-29, root cause of the
// non-90-degree MS / MA over-correction): on every ring that meets the metal (a ring with
// PEC knots), the two crossings of the metal arms with the ring must be PEC knots (within
// `tolerance` of a zero-set knot) and no free knot may lie on the metal footprint (the
// closed sector between the arms for a convex corner, its complement for a concave one; a
// rounded corner's arms are straight where they cross the box). Only rings that are square
// box rings centred on the apex (the generator's matching box) are judged. Returns the
// empty string when the basis passes, else the reason.
std::string CheckCornerBasisCrossings(const std::vector<std::array<double, 3>> &points,
                                      const std::vector<int> &contour_groups,
                                      const std::vector<int> &zero_trace_indices,
                                      double angle_radians, bool convex, double tolerance);

}  // namespace palace

#endif  // PALACE_MODELS_CORNER_TRACE_BASIS_HPP
