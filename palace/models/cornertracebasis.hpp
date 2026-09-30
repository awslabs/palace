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
//
// Two ring LAYOUTS (record TraceBasis.RingLayout). MetalRingsOnly (the recorded rule above,
// the default when the record has no RingLayout): the two metal rings follow the rule, every
// other ring is the fixed layout; the family has EVENTS (a knot passing a fixed vertex flips
// the band triangulation) and interpolates within segments of one connectivity.
// AllRingsFollowMetal (corner-basis refinement, USER decision 161 (1), 2026-09-30): EVERY
// ring — the outer rings (the standard levels plus the extra ring at MetalThickness +
// OveretchDepth mirroring the trench ring) and the two inner cap rings — carries the same
// angle-dependent fractions, so every band is a regular column grid with one diagonal
// orientation; the knots are PEC on the two metal rings only; the box corners are slaves on
// every ring; each cap is a fan from a centre slave at the mean of the cap ring's two
// crossing-slot knots. No knot ever passes a fixed vertex and a knot passing a corner slave
// on every ring at once leaves the hats continuous: NO events, one segment per convexity,
// no connectivity angle. The free knots are GRADED: free_knot_grading lists perimeter
// distances (over R) from each end of the free arc (both crossings) carrying a knot, the
// remaining free knots at equal fractions between the innermost graded ones (the option-(c)
// held-out trace ramps over R / 3 from the metal arc).
enum class CornerRingLayout : char
{
  METAL_RINGS_ONLY,
  ALL_RINGS_FOLLOW_METAL
};

struct CornerTraceBasisRule
{
  int ring_size = 8;
  int metal_interior_knots = 1;
  int free_knots = 5;
  std::string fractions = "PerimeterArcLength";
  CornerRingLayout ring_layout = CornerRingLayout::METAL_RINGS_ONLY;
  std::vector<double> free_knot_grading;
  // Extra outer rings above the metal top at MetalThickness + k OveretchDepth
  // (AllRingsFollowMetal; k = 1 mirrors the trench ring, k = 4 resolves the trace right above
  // the metal top over the metal arc: the concave family's fabricated MA read +7 % without it).
  std::vector<double> extra_levels_above_over_overetch;

  bool AllRings() const { return ring_layout == CornerRingLayout::ALL_RINGS_FOLLOW_METAL; }

  bool operator==(const CornerTraceBasisRule &other) const
  {
    return ring_size == other.ring_size &&
           metal_interior_knots == other.metal_interior_knots &&
           free_knots == other.free_knots && fractions == other.fractions &&
           ring_layout == other.ring_layout && free_knot_grading == other.free_knot_grading &&
           extra_levels_above_over_overetch == other.extra_levels_above_over_overetch;
  }
  bool operator!=(const CornerTraceBasisRule &other) const { return !(*this == other); }
};

// The refined rule of the rebuilt corner family (generate_corner_response.REFINED_RULE;
// corner-basis-refinement-20260930 Phase 1 candidate C14): 5 metal-interior + 9 graded free
// knots (R/3, 2R/3 from each crossing), rings at t + d and t + 4d; 11 rings x 16 knots = 176.
CornerTraceBasisRule RefinedCornerTraceBasisRule();

// The outer ring levels of an AllRingsFollowMetal coupon: -R, -R/3, -d, 0, t, t + k d for the
// rule's k, R/3, R (sorted).
std::vector<double> CornerRuleLevels(const CornerTraceBasisRule &rule, double radius,
                                     double metal_thickness, double overetch_depth);

// The reason a rule is invalid (RingSize = 2 + MetalInteriorKnots + FreeKnots, the
// fraction parametrisation, an increasing positive grading fitting the free knots, no
// grading on MetalRingsOnly), empty when valid.
std::string CheckCornerTraceBasisRule(const CornerTraceBasisRule &rule);

// A vertex of a ring that meets the metal: its perimeter fraction in [0, 1), whether it is
// a PEC knot, a free knot or a slave (box corner), and for a knot its basis slot in the
// ring.
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

// The knots and slave vertices of a ring laid out by the rule, in perimeter order: `pec`
// true for a ring that meets the metal (crossings and metal-interior knots in the zero set),
// false for a ring that follows the metal fractions off the metal (every knot free).
std::vector<CornerRingVertex> CornerRuleRingLayout(double radius, double angle_radians,
                                                   bool convex,
                                                   const CornerTraceBasisRule &rule, bool pec);

// The knots and slave vertices of a ring that meets the metal, in perimeter order
// (CornerRuleRingLayout with pec).
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
// (outer rings at z = -R, -R / 3, -OveretchDepth, 0, MetalThickness, R / 3, R — plus the
// rule's extra levels above the metal top for the AllRingsFollowMetal layout — then the top and
// bottom inner cap rings of half width R / 3) with the rule's zero set on the two rings
// that meet the metal. BuildCornerTraceBasis on this seed gives the rule's basis at any
// angle (unit tests; the Python generator is the reference for coupon files).
struct CornerBoxSeed
{
  std::vector<std::array<double, 3>> points;
  std::vector<int> contour_groups;
  std::vector<int> zero_trace_indices;  // 0-based
};

CornerBoxSeed MakeCornerBoxSeed(double radius, double metal_thickness,
                                double overetch_depth, bool convex,
                                const CornerTraceBasisRule &rule);

// The trace basis of the family at any corner angle: the node's fixed rings (its own
// points) and the metal rings laid out by the rule at `angle_radians`, triangulated by
// perimeter fraction (equal vertex sets reproduce the generator's connect_rings). With a
// connectivity angle the bands next to the metal rings are merged in the order of the
// rule's layout at THAT angle (the generator's connectivity_keys: a knot's key is its
// role's fraction at the connectivity angle, a slave's its corner), so the triangulation is
// the same for every angle of a segment (throws when a knot-corner passage lies between the
// two angles: the triangles would fold). Under the AllRingsFollowMetal layout every ring
// (the node's ring heights and half widths, the rule's fractions; PEC on the rings the
// node's zero set marks) is laid out by the rule, the caps are fans from centre slaves and
// no connectivity angle is accepted. Used for the runtime model of an interpolated corner
// (the nodes' knot semantics are the same; the positions at the device angle are the
// rule's, the connectivity the segment's) and, at a node's own angle, to check a coupon's
// files against the rule.
ConstructedCornerTraceBasis
BuildCornerTraceBasis(const std::vector<std::array<double, 3>> &node_points,
                      const std::vector<int> &contour_groups,
                      const std::vector<int> &zero_trace_indices, double angle_radians,
                      bool convex, const CornerTraceBasisRule &rule,
                      std::optional<double> connectivity_angle_radians = std::nullopt);

// Geometric events of the rule's trace basis in the corner angle (corner-qualification
// block 2026-09-29): a knot of a ring that meets the metal passes a vertex of the fixed
// layout — a box corner (`corner` true) or a side midpoint (`corner` false). At every such
// angle the band triangulation of the perimeter-ordered merge flips a quad diagonal and the
// hats JUMP (measured O(1)); a knot passing a box corner also swaps its order with the
// corner's slave vertex, so no connectivity is continuous across it. Every event in
// (0, 180) degrees of every knot (crossings, metal-interior and free), sorted by angle. The
// `fraction` is the fixed-layout vertex passed. The AllRingsFollowMetal layout has NO
// events (no fixed vertex; a knot passing a corner slave on every ring at once is
// continuous, measured 1e-5 at the MetalRingsOnly event angles): the empty list.
struct CornerBasisEvent
{
  double angle_degrees = 0.0;
  std::string role;
  double fraction = 0.0;
  bool corner = false;
};

std::vector<CornerBasisEvent> CornerBasisEvents(bool convex,
                                                const CornerTraceBasisRule &rule);

// Angle-interpolation stencil of the corner family (MatchCornerFamily's window rule,
// mirrored by corner_family_interpolation.py). Nodes are the family's sharp coupons of one
// convexity: angle and, for coupons built with a segment connectivity (TraceBasis
// ConnectivityAngleDegrees), that angle; a coupon without one is a LEGACY node (the
// perimeter-ordered merge at its own angle). Rule: coupons sharing a connectivity angle
// form a SEGMENT whose node angles must lie in one corner-event-free interval together
// with the connectivity angle (fail closed otherwise); a device angle strictly inside a
// segment's node range is interpolated by Lagrange on the segment's nodes nearest to it
// (cubic on four, else quadratic / linear: never across a corner event); an angle equal to
// a node (within the tolerance) is exact — with several coupons at that angle the legacy
// one (the tie triangulation of the recorded 90 / 135 / 180 coupons) is preferred, else the
// coupon of the lower-angle segment; legacy nodes are never interpolated (their bases jump
// at every event: refused with the reason); outside the node range refused (no
// extrapolation). `reason` is empty when a stencil was found. The node tolerance is
// closed on the exact side and open on the interior side with the same floating-point
// difference (|node - angle| <= tol exact, > tol interior): an angle at the boundary is
// never in no segment (qualification review 2026-09-29 m1).
struct CornerFamilyNode
{
  double angle_degrees = 0.0;
  std::optional<double> connectivity_angle_degrees;
  std::size_t index = 0;
};

struct CornerFamilyStencil
{
  std::vector<std::pair<std::size_t, double>> nodes;  // (node index, Lagrange weight)
  std::string rule;                                   // exact / linear / quadratic / cubic
  std::optional<double> connectivity_angle_degrees;   // the segment's (interpolated)
  std::size_t base = 0;                               // the nearest node of the stencil
  double min_angle_degrees = 0.0, max_angle_degrees = 0.0;  // the family's node range
  std::string reason;
};

CornerFamilyStencil SelectCornerFamilyStencil(const std::vector<CornerFamilyNode> &nodes,
                                              double angle_degrees, bool convex,
                                              const CornerTraceBasisRule &rule,
                                              double angle_tolerance_degrees);

// The structural check of a corner family's segments (decision 137 (1), verified at library
// load by ReadProcessLibrary, fail closed): every segment's connectivity angle off the
// knot-corner passages and its nodes in one event-free interval with it, one coupon per
// angle within a segment, segments overlapping at shared node angles only. Returns the
// empty string when consistent, else the reason.
std::string CheckCornerFamilySegments(const std::vector<CornerFamilyNode> &nodes,
                                      bool convex, const CornerTraceBasisRule &rule,
                                      double angle_tolerance_degrees);

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
