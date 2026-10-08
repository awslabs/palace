// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfaceresponseidentification.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <iomanip>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <mfem.hpp>
#include "utils/enum_string.hpp"

namespace palace
{

namespace
{

using Point3D = std::array<double, 3>;
using Interval = std::pair<double, double>;

// The decision-69 quantizer: lengths on the 1e-8 R grid, direction cosines on the 1e-12
// grid, every "within" test a strict less-than on the quantized values (design (c)).
constexpr double kLengthQuantumOverRadius = 1.0e-8;
constexpr double kDirectionQuantum = 1.0e-12;
// The corner threshold is the geometric joint noise rule kJointNoiseSagittaOverRadius of
// metaledge.hpp (USER decision 121 (B)): every joint whose implied sagitta on its shorter
// adjacent piece reaches 0.05 R is a corner unless a fitted arc absorbs it; the constant
// lives with the perimeter extraction that classifies the vertices, so that the chains, the
// identification (the arc rule's end-joint test) and the audit mirror read one value.
constexpr double kInteractionDistanceOverRadius = 2.0;
constexpr double kThroughVertexZoneOverRadius = 2.0;
constexpr double kClusterBallOverRadius = 1.0;
// One interaction distance (decision 82(1)): every join / candidate / cluster decision is
// the 3D distance strictly below 2R on the quantized grid — events, through-vertex zones,
// core merging, site-site cores and the vertex-site join (a site within 2R of a core: its
// radius-R window overlaps the core's radius-R ball; was 3R, which joined the L2 junction
// regions of DS-SCT-002 to the L1 loop ends 2.4R below). The cluster ball (R) is the claim
// radius of the region (union of radius-R balls), not an interaction decision.
constexpr double kVertexJoinsClusterOverRadius = 2.0;
constexpr double kVertexWindowOverRadius = 1.0;
constexpr double kParallelCosineTolerance = 1.0e-8;
// Arc rule (decision 82(3), design (b) 4 / 7; the CONCYCLICITY form of USER decisions 121 /
// 122, 2026-09-28, amending the sagitta form of 117(4)): a run of at least three
// consecutive joints of the perimeter path (through corner vertices) turning the same way
// by at most 180 deg in total is ONE arc iff its joint vertices lie on one circle within
// the signature parameter tolerance (kArcFitToleranceOverRadius x R = 1e-3 R: the circle
// tangent to both arms when the joints lie on it, else — for a bend of radius >= R over at
// least FOUR joints — the least-squares circle of the joints, the arms meeting it within
// the joint noise rule of its tangents or as chords of it) AND every joint turns less than
// kArcMaxJointTurnDegrees (50: regular polygons such as squares and hexagons stay corners).
// The chord sagitta rho (1 - cos(central angle / 2)) no longer decides membership: it is
// recorded per arc (max_sagitta -> Arcs[].MaxChordSagittaOverR) and arcs whose largest
// chord sagitta reaches kArcSagittaOverRadius x R (0.05 R, the resolution of the
// correction) are listed in the manifest's MeshCoarsenessWarning — a coarsely meshed design
// curve is the curve (the transmon's 19.4 R CPW bends at 11-16 deg per chord are bends,
// not 387 corners of 164-175 deg), and its coarseness is reported. An arc of radius below R
// is ONE rounded corner (total turn, radius / R in the signature; its arms are tangent to
// the circle by construction), an arc of radius >= R is a bend of exactly that radius for
// the curved-edge chain rule (its joints, corner vertices included, are no features and its
// density is 1 / radius over the arc). Every joint that no arc absorbs and that is not
// noise under the geometric rule is a corner feature: a 4-chord semicircle of radius 20 R
// (joints of 45 deg on one circle) is a bend, a square hole is four corners, a chamfer (two
// joints) stays two corners. The pieces between joints carry no length constraint of their
// own; the former 2R piece rule with its sub-corner exception is replaced.
// kArcFitToleranceOverRadius / kArcSagittaOverRadius: surfaceresponseidentification.hpp.
// Cap of the wedge a claimed piece end extends by in the cluster extension's across rule
// (2R tan(cap) = 1.15 R): formerly the corner threshold itself; kept at 30 deg so that the
// extension is unchanged by the noise threshold (a piece end at a sub-noise joint or at an
// arc joint turns by less anyway; a corner site's window end reaches at most this wedge).
// Recorded as Conventions.ClusterExtensionWedgeCapDegrees.
constexpr double kClusterExtensionWedgeCapDegrees = 30.0;
// Curve-aware cluster geometry (option A, decision 91(1)): every run on a fitted arc
// (rounded corner or bend) is evaluated on the arc — event cores, core merging, the vertex
// join, the radius-R claims, the extension and the signature portions — so that a cluster's
// extent and signature are functions of the design curve, not of the chord count. The
// distance along an arc is not convex: its sublevel sets are bracketed by samples at most
// this spacing apart (in units of R; at least kArcSampleMinimumIntervals per piece) and
// their crossings located by bisection; recorded as Conventions.ArcSampleSpacingOverR. The
// library builder chords a signature arc at the canonical angular step
// kClusterArcChordStepDegrees (or finer so that no chord exceeds
// kClusterArcChordMaxLengthOverR x R), so that the coupon geometry is mesh independent.
constexpr double kArcSampleSpacingOverRadius = 0.05;
constexpr double kClusterArcChordStepDegrees = 5.0;
constexpr double kClusterArcChordMaxLengthOverRadius = 0.25;

// Curved-edge chain rule (decision 73(1), design (b) 7). The turn at every sub-corner joint
// of a chain is spread over the two adjacent half-chords (the polyline's discrete curvature
// density, exact for a polygon inscribed in a circle at any chord length); the windowed
// curvature at a chain point is the mean density over a window of length R centred on it
// (the response at a point integrates the geometry within ~R). A point whose windowed bend
// radius is below kStraightBendRadiusOverRadius x R is curved: the first-order curvature
// correction to an edge response scales as R / radius (the response integrates the field
// over distances <= R from the edge, and an in-plane bend perturbs that geometry at
// relative order R / radius), so at 10 R it is <= 10 % of the local edge correction, i.e.
// ~1e-3 of the corrected edge participation for edge corrections of a few per cent
// (decision 75: the transmon's 38.9 um = 19.4 R CPW bends are straight-like; curved coupons
// for 1 < radius / R < 10 are deferred to the library regeneration). Two chains are a pair
// along a bend where their closest-point
// separation is LOCALLY constant: a sample of one chain's facing region (within the
// candidate reach 2R (1 + kPairSeparationTolerance) of the other chain, not beyond its
// ends, outside the shared-vertex zones; samples at most R / 2 apart) is constant when the
// sampled distances within R of it along its own chain vary by at most
// kPairSeparationTolerance x their minimum (a polyline of sub-noise turns at constant
// width varies by less than 1 / cos(0.5 deg) - 1; a bend arc's chords of sagitta <= 0.05 R
// at constant width vary by up to 0.05 R / separation, 2.5 % at 2R; the pair response
// sensitivity d dR/dd is O(1)). The constant portions are the pair; the portions that are
// not (divergence at tees and port ends, fast tapers, acute corner arms) keep the event
// rule, so a slow taper is a pair and a tee is a cluster. Whether a constant portion
// interacts is decided on the separation of the underlying curves, so that a discretisation
// never changes the classification: the chords of a polyline inscribed in a curve lie
// inside it (two concentric inscribed polylines are w cos(turn / 2) apart mid-chord and
// exactly w apart at their vertices), while the exact offset polyline of a bent path keeps
// corresponding chords at the design separation and the outer side's samples near the
// joints project onto the inner vertices at up to w / cos(turn / 2). In both constructions
// the sampled closest-point distance from one chain to the other reaches the curve
// separation w as its maximum on the side whose maximum is smaller: the sample's chord
// reading C = min over the two chains of the maximum sampled distance within a window of
// half-width max(R, local chord) where the chain bends and R on straight runs, about the
// sample on its own chain and about its foot on the other chain (a straight taper is read
// locally). C is exact for an offset polyline; two polylines inscribed in the curves at
// aligned angles are C = w cos(turn / 2) apart everywhere (chords and vertex-to-polyline
// alike) for a curve separation w, and the polyline pair alone cannot tell the two
// constructions apart (they differ at order turn^2): the inscribed reading is C / cos(turn
// / 2) with turn = the larger local joint turn of the two chains. A portion interacts iff
// BOTH readings are below 2R on the quantized grid — the same strict-less decision as a
// straight parallel pair at that separation, taken on the non-interacting side of the
// recorded ambiguity w (1 / cos(turn / 2) - 1) (below 1e-4 w for joints under 1.6 deg; a
// CPW gap of exactly 2R along a bend is isolated edges like a straight one; DS-SCT-001's 4
// um gaps at R = 2 um read 3.9998 mid-chord and became 3 mm clusters). The pair feature's
// separation is the mean chord reading over its samples; the cross-chord interactions of a
// locally constant portion are never event cores, whether or not it interacts.
// (kStraightBendRadiusOverRadius = 10 is declared in the header: the curvature families'
// first-order rule is keyed to it.)
constexpr double kCurvatureWindowOverRadius = 1.0;
constexpr double kPairSeparationTolerance = 0.05;
constexpr int kPairSeparationSamplesPerInterval = 16;
// Locality of the bent-pair readings (USER decision 184 (3), implementation corrected
// 2026-10-01): a chain "bends within R" of a sample when a joint with a turn or a fitted
// bend arc lies within kPairBendProximityOverRadius x R of it ALONG THE CHAIN; there and
// only there the sample has no exact straight reading and its chord reading is the window
// maximum over a half-width of max(R, min(local chord, kPairChordWindowCapOverRadius x R))
// instead of R. The windowed curvature is NOT this test: the curvature rule spreads every
// joint's turn over its two half-runs, so a 0.1 deg taper kink between a 200 um lead and a
// 400 um taper run made kappa > 0 over 100 um of the lead, and the former window of a full
// run length (200 um) read the taper from the lead: the DS-CTX-003 532 um 1 / 2 / 1 um
// stacks read 2.2 / 2.7 / 2.0 um, the straight leads of the re-mesh gate's 10 R curved
// stacks mesh-dependent. The cap keeps the window of a coarse chord (beyond 2R: chord
// readings anyway) from reaching the neighbouring piece.
constexpr double kPairBendProximityOverRadius = 1.0;
constexpr double kPairChordWindowCapOverRadius = 2.0;
// Exact stretches (USER decision 203, 2026-10-02): within a locally constant sub-piece the
// samples with an exact reading agreeing within the signature parameter tolerance form
// exact stretches; one splits off as its own exact sub-piece only when its exact samples
// span at least kExactStretchMinLengthOverRadius x R (a shorter one merges into the
// adjacent non-exact stretch; a sub-piece that is one exact stretch stays exact at any
// length). The near-bend samples (not judged) join an adjacent stretch within
// kPairBendProximityOverRadius x R of its last exact / judged sample.
constexpr double kExactStretchMinLengthOverRadius = 1.0;
// Self-pairing (decision 82(2) addition, 2026-09-25): a chain folding back onto itself
// within 2R through a bend of radius >= R (a hairpin, a meander with smooth bends, the
// U-turn of a narrow strip whose fold is not a rounded corner) pairs with itself. The
// local neighbourhood along the chain takes no part: two points of one chain closer than
// pi R along the chain are within 2R of each other along any bend of radius >= R (on the
// tightest bend, radius R, the chord reaches 2R exactly after half a turn = pi R of arc;
// on a wider bend earlier), so only points more than pi R apart along the chain can face
// each other across a fold. Runs closer than pi R of arc length are never partners.
constexpr double kSelfPairNeighbourhoodOverRadius = 3.14159265358979323846;
// Stack composition cap (decision 82(2) assembly): the breadth-first traversal of the
// links from a chain position stops after this many members. A cross-section of k edges
// spans at least (k - 1) x (the smallest gap) laterally; the cap is far above any physical
// stack (the chip's widest stack has 12 edges) and only bounds a runaway traversal through
// inconsistent links; every hit is counted under Diagnostics.StackCompositionCapHits.
constexpr std::size_t kStackCompositionCap = 64;
// Stack-end images (decision 93, DS-CTX-003 defect 1): the foot of a breakpoint on a
// partner chain within this distance (units of R) of an existing breakpoint of that chain
// is the same cut and creates no new breakpoint; the signature parameter tolerance, below
// which the contract resolves no parameter.
constexpr double kStackImageToleranceOverRadius = kSignatureParameterToleranceOverRadius;
// Cluster extension closure (decision 93, DS-CTX-003 defect 2): the extension / stack
// recomposition loop stops after this many passes, or when a pass repeats the previous one
// (the same absorbed length and portion count: the recomposed stacks re-cut the same
// slivers at the moved cut images, DS-CTX-003 80 nm / 16 portions in passes 2-5); the last
// pass is applied without a further recomposition, as for a sub-tolerance pass.
constexpr std::size_t kClusterExtensionMaxPasses = 12;
// Knife-edge census (decision 82(4), 2026-09-25): every rule of the identification is a
// strict comparison with a threshold (the interaction distance 2R, the cluster ball /
// vertex window R, the straight-bend radius 10R, the corner (joint noise) turn, the arc
// sagitta 0.05 R); a design whose dimensions sit on a threshold is decided by roundoff. The
// manifest reports the perimeter length whose distance to the nearest other perimeter point
// (3D; the same chain beyond the self-pair neighbourhood) lies within the band of R and of
// 2R, the chain length whose windowed bend radius lies within the band of 10R, the vertex
// count whose turn lies within the band of the corner threshold and the vertex count whose
// implied chord sagitta lies within the band of 0.05 R, each split into the below / above
// sides, sampled every kKnifeEdgeSampleSpacingOverR; and (USER decision 184 (4)) the
// cluster-composition band: the claimed length and count of the clusters whose edge count
// or member vertices differ when the same perimeter is identified at R (1 -/+ band)
// (Identifier::ClusterCompositionBand).
constexpr double kKnifeEdgeBandRelative = 0.01;
constexpr double kKnifeEdgeSampleSpacingOverRadius = 0.5;

// ---------------------------------------------------------------------------------------
// Vector helpers
// ---------------------------------------------------------------------------------------

double Dot(const Point3D &a, const Point3D &b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Point3D Sub(const Point3D &a, const Point3D &b)
{
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}
Point3D Add(const Point3D &a, const Point3D &b)
{
  return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}
Point3D Scale(double s, const Point3D &a)
{
  return {s * a[0], s * a[1], s * a[2]};
}
Point3D Cross(const Point3D &a, const Point3D &b)
{
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
double Norm(const Point3D &a)
{
  return std::sqrt(Dot(a, a));
}
Point3D Normalize(const Point3D &a)
{
  const double n = Norm(a);
  MFEM_VERIFY(n > 0.0, "Cannot normalize a zero vector in the geometry identification!");
  return Scale(1.0 / n, a);
}
double Distance(const Point3D &a, const Point3D &b)
{
  return Norm(Sub(a, b));
}

double PointSegmentDistance(const Point3D &p, const Point3D &a, const Point3D &b)
{
  const Point3D d = Sub(b, a);
  const double dd = Dot(d, d);
  double t = dd > 0.0 ? Dot(Sub(p, a), d) / dd : 0.0;
  t = std::clamp(t, 0.0, 1.0);
  return Distance(p, Add(a, Scale(t, d)));
}

// Closest approach between two segments (Ericson, Real-Time Collision Detection 5.1.9).
double SegmentSegmentDistance(const Point3D &p0, const Point3D &p1, const Point3D &q0,
                              const Point3D &q1)
{
  const Point3D d1 = Sub(p1, p0), d2 = Sub(q1, q0), r = Sub(p0, q0);
  const double a = Dot(d1, d1), e = Dot(d2, d2), f = Dot(d2, r);
  double s = 0.0, t = 0.0;
  if (a <= 0.0 && e <= 0.0)
  {
    return Distance(p0, q0);
  }
  if (a <= 0.0)
  {
    t = std::clamp(f / e, 0.0, 1.0);
  }
  else
  {
    const double c = Dot(d1, r);
    if (e <= 0.0)
    {
      s = std::clamp(-c / a, 0.0, 1.0);
    }
    else
    {
      const double b = Dot(d1, d2);
      const double denom = a * e - b * b;
      s = denom != 0.0 ? std::clamp((b * f - c * e) / denom, 0.0, 1.0) : 0.0;
      t = (b * s + f) / e;
      if (t < 0.0)
      {
        t = 0.0;
        s = std::clamp(-c / a, 0.0, 1.0);
      }
      else if (t > 1.0)
      {
        t = 1.0;
        s = std::clamp((b - c) / a, 0.0, 1.0);
      }
    }
  }
  return Distance(Add(p0, Scale(s, d1)), Add(q0, Scale(t, d2)));
}

// Distance from a point to a planar polygon (convex, or the sampled loop of a curved face
// treated as convex): the normal distance when the projection falls inside, else the
// distance to the nearest polygon edge.
double PointPolygonDistance(const Point3D &p, const std::vector<Point3D> &polygon,
                            const Point3D &normal)
{
  const std::size_t n = polygon.size();
  bool inside = n >= 3;
  for (std::size_t i = 0; i < n && inside; i++)
  {
    const Point3D &a = polygon[i];
    const Point3D &b = polygon[(i + 1) % n];
    const Point3D &c = polygon[(i + 2) % n];
    const Point3D edge = Sub(b, a);
    const Point3D interior = Cross(normal, edge);  // towards the polygon interior (or away)
    const double orientation = Dot(Sub(c, b), interior);
    const double side = Dot(Sub(p, a), interior);
    inside = orientation * side >= 0.0;
  }
  if (inside)
  {
    return std::abs(Dot(Sub(p, polygon[0]), normal));
  }
  double best = std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < n; i++)
  {
    best = std::min(best, PointSegmentDistance(p, polygon[i], polygon[(i + 1) % n]));
  }
  return best;
}

double QuantizeDirection(double cosine)
{
  return std::round(cosine / kDirectionQuantum);
}
bool DirectionLess(double first, double second)
{
  return QuantizeDirection(first) < QuantizeDirection(second);
}

class Quantizer
{
public:
  explicit Quantizer(double radius) : quantum(kLengthQuantumOverRadius * radius)
  {
    MFEM_VERIFY(std::isfinite(quantum) && quantum > 0.0,
                "Invalid matching radius for the geometry identification!");
  }
  double Q(double length) const { return std::round(length / quantum); }
  bool Less(double a, double b) const { return Q(a) < Q(b); }
  bool Equal(double a, double b) const { return Q(a) == Q(b); }
  double Snap(double length) const
  {
    const double snapped = Q(length) * quantum;
    return snapped == 0.0 ? 0.0 : snapped;
  }
  double Quantum() const { return quantum; }

private:
  double quantum;
};

// ---------------------------------------------------------------------------------------
// Interval utilities (closed intervals on a run parameter, merged canonically)
// ---------------------------------------------------------------------------------------

std::vector<Interval> MergeIntervals(std::vector<Interval> intervals, double tolerance)
{
  std::sort(intervals.begin(), intervals.end());
  std::vector<Interval> merged;
  for (const auto &interval : intervals)
  {
    if (interval.second - interval.first <= tolerance)
    {
      continue;
    }
    if (!merged.empty() && interval.first <= merged.back().second + tolerance)
    {
      merged.back().second = std::max(merged.back().second, interval.second);
    }
    else
    {
      merged.push_back(interval);
    }
  }
  return merged;
}

std::vector<Interval> SubtractIntervals(const std::vector<Interval> &from,
                                        const std::vector<Interval> &minus,
                                        double tolerance)
{
  std::vector<Interval> result = from;
  for (const auto &cut : minus)
  {
    std::vector<Interval> next;
    for (const auto &interval : result)
    {
      if (cut.second <= interval.first + tolerance ||
          cut.first >= interval.second - tolerance)
      {
        next.push_back(interval);
        continue;
      }
      if (cut.first > interval.first + tolerance)
      {
        next.emplace_back(interval.first, cut.first);
      }
      if (cut.second < interval.second - tolerance)
      {
        next.emplace_back(cut.second, interval.second);
      }
    }
    result = std::move(next);
  }
  return MergeIntervals(result, tolerance);
}

std::vector<Interval> IntersectIntervals(const std::vector<Interval> &a,
                                         const std::vector<Interval> &b, double tolerance)
{
  std::vector<Interval> result;
  for (const auto &x : a)
  {
    for (const auto &y : b)
    {
      const double lo = std::max(x.first, y.first), hi = std::min(x.second, y.second);
      if (hi - lo > tolerance)
      {
        result.emplace_back(lo, hi);
      }
    }
  }
  return MergeIntervals(result, tolerance);
}

// Sublevel set {s in [lo, hi] : f(s) < level} of a convex function, as one interval. The
// minimizer is located by golden-section search and the crossings by bisection; the
// emptiness decision is the quantized comparison of the minimum with the level.
std::optional<Interval> ConvexSublevelInterval(const std::function<double(double)> &f,
                                               double lo, double hi, double level,
                                               const Quantizer &quantizer)
{
  if (hi - lo <= 0.0)
  {
    return std::nullopt;
  }
  constexpr double golden = 0.6180339887498949;
  double a = lo, b = hi;
  double c = b - golden * (b - a), d = a + golden * (b - a);
  double fc = f(c), fd = f(d);
  for (int it = 0; it < 90 && (b - a) > 1.0e-14 * (hi - lo); it++)
  {
    if (fc < fd)
    {
      b = d;
      d = c;
      fd = fc;
      c = b - golden * (b - a);
      fc = f(c);
    }
    else
    {
      a = c;
      c = d;
      fc = fd;
      d = a + golden * (b - a);
      fd = f(d);
    }
  }
  double s_min = 0.5 * (a + b);
  double f_min = f(s_min);
  for (const double candidate : {lo, hi})
  {
    const double value = f(candidate);
    if (value < f_min)
    {
      f_min = value;
      s_min = candidate;
    }
  }
  if (!quantizer.Less(f_min, level))
  {
    return std::nullopt;
  }
  auto Bisect = [&](double inside, double outside)
  {
    // f(inside) < level <= f(outside) (or outside is the domain end with f < level).
    if (f(outside) < level)
    {
      return outside;
    }
    for (int it = 0; it < 100 && std::abs(outside - inside) > 1.0e-15 * (hi - lo); it++)
    {
      const double mid = 0.5 * (inside + outside);
      (f(mid) < level ? inside : outside) = mid;
    }
    return inside;
  };
  return Interval{Bisect(s_min, lo), Bisect(s_min, hi)};
}

// Minimizer of a continuous function on [a, b] by golden-section search (the bracket comes
// from a sampled local minimum: the function is unimodal there).
std::pair<double, double> GoldenMinimum(const std::function<double(double)> &f, double a,
                                        double b)
{
  constexpr double golden = 0.6180339887498949;
  const double span = b - a;
  double c = b - golden * (b - a), d = a + golden * (b - a);
  double fc = f(c), fd = f(d);
  for (int it = 0; it < 90 && (b - a) > 1.0e-14 * span; it++)
  {
    if (fc < fd)
    {
      b = d;
      d = c;
      fd = fc;
      c = b - golden * (b - a);
      fc = f(c);
    }
    else
    {
      a = c;
      c = d;
      fc = fd;
      d = a + golden * (b - a);
      fd = f(d);
    }
  }
  const double s = 0.5 * (a + b);
  return {s, f(s)};
}

// Sublevel set {s in [lo, hi] : f(s) < level} of a continuous function that need not be
// convex (the distance along a circular arc, or a distance restricted to a projection
// domain, +inf outside it), as merged intervals: f is sampled every `step` (at least
// kArcSampleMinimumIntervals intervals), every local minimum of the samples is refined by
// golden section, the emptiness is decided by the quantized comparison of the smallest
// refined minimum with the level (as for ConvexSublevelInterval) and every crossing between
// consecutive points of different status is located by bisection on the predicate. The
// result depends on the geometry alone up to the bisection precision (the sample positions
// only bracket the crossings), so that a re-meshed arc gives the same intervals.
constexpr std::size_t kArcSampleMinimumIntervals = 8;
std::vector<Interval> SampledSublevelIntervals(const std::function<double(double)> &f,
                                               double lo, double hi, double level,
                                               const Quantizer &quantizer, double step)
{
  std::vector<Interval> result;
  if (hi - lo <= 0.0 || !(step > 0.0))
  {
    return result;
  }
  const std::size_t n = std::max<std::size_t>(
      kArcSampleMinimumIntervals, static_cast<std::size_t>(std::ceil((hi - lo) / step)));
  std::vector<std::pair<double, double>> points;  // (s, f(s)), sorted by s
  points.reserve(2 * n + 2);
  for (std::size_t k = 0; k <= n; k++)
  {
    const double s =
        k == n ? hi : lo + (hi - lo) * static_cast<double>(k) / static_cast<double>(n);
    points.emplace_back(s, f(s));
  }
  const std::size_t samples = points.size();
  std::vector<std::pair<double, double>> refined;
  for (std::size_t k = 0; k < samples; k++)
  {
    const double value = points[k].second;
    if (!std::isfinite(value))
    {
      continue;
    }
    const bool left_ok = k == 0 || !(points[k - 1].second < value);
    const bool right_ok = k + 1 == samples || !(points[k + 1].second < value);
    if (left_ok && right_ok)
    {
      const double a = points[k > 0 ? k - 1 : 0].first;
      const double b = points[std::min(k + 1, samples - 1)].first;
      auto minimum = GoldenMinimum(f, a, b);
      if (!(minimum.second < value))
      {
        minimum = points[k];
      }
      refined.push_back(minimum);
    }
  }
  double f_min = std::numeric_limits<double>::infinity();
  for (const auto &point : refined)
  {
    f_min = std::min(f_min, point.second);
  }
  for (const auto &point : points)
  {
    f_min = std::min(f_min, point.second);
  }
  if (!std::isfinite(f_min) || !quantizer.Less(f_min, level))
  {
    return result;
  }
  points.insert(points.end(), refined.begin(), refined.end());
  std::sort(points.begin(), points.end());
  auto Bisect = [&](double inside, double outside)
  {
    for (int it = 0; it < 100 && std::abs(outside - inside) > 1.0e-15 * (hi - lo); it++)
    {
      const double mid = 0.5 * (inside + outside);
      (f(mid) < level ? inside : outside) = mid;
    }
    return inside;
  };
  std::size_t i = 0;
  while (i < points.size())
  {
    if (!(points[i].second < level))
    {
      i++;
      continue;
    }
    std::size_t j = i;
    while (j + 1 < points.size() && points[j + 1].second < level)
    {
      j++;
    }
    const double start = i == 0 ? lo : Bisect(points[i].first, points[i - 1].first);
    const double end =
        j + 1 == points.size() ? hi : Bisect(points[j].first, points[j + 1].first);
    result.emplace_back(start, end);
    i = j + 1;
  }
  return MergeIntervals(std::move(result), 0.0);
}

// ---------------------------------------------------------------------------------------
// Circular arcs (option A, decision 91(1)): the geometry of the perimeter on a fitted arc
// ---------------------------------------------------------------------------------------

// A circular arc piece: centre, radius, orthonormal in-plane basis (u, v) and the signed
// angular range from theta0 to theta1 (point(theta) = c + r (cos theta u + sin theta v)).
// The runs of a fitted arc (rounded corner or bend) are chords of one circle; their claims
// are evaluated on the arc so that a cluster's extent is a function of the design curve and
// not of the chord count (fix-4 residual 1: the hairpin cluster moved 12.08-13.40 um and
// the u-ring cluster 73.12-73.36 um across re-meshes on chord geometry).
struct ArcPiece
{
  Point3D center{};
  double radius = 0.0;
  Point3D u{}, v{};
  double theta0 = 0.0, theta1 = 0.0;

  Point3D At(double theta) const
  {
    return Add(center,
               Add(Scale(radius * std::cos(theta), u), Scale(radius * std::sin(theta), v)));
  }
  // Point at the fraction t in [0, 1] of the angular range.
  Point3D AtFraction(double t) const { return At(theta0 + t * (theta1 - theta0)); }
  double Sweep() const { return std::abs(theta1 - theta0); }
  double Length() const { return radius * Sweep(); }
  Point3D Normal() const { return Cross(u, v); }
  // Unit tangent in the direction of travel (theta0 -> theta1) and the outward radial
  // direction at theta.
  Point3D Tangent(double theta) const
  {
    const double sign = theta1 >= theta0 ? 1.0 : -1.0;
    return Scale(sign, Add(Scale(-std::sin(theta), u), Scale(std::cos(theta), v)));
  }
  Point3D Radial(double theta) const
  {
    return Add(Scale(std::cos(theta), u), Scale(std::sin(theta), v));
  }
  // The sub-arc between the fractions t0 and t1 of the angular range.
  ArcPiece SubArc(double t0, double t1) const
  {
    ArcPiece piece = *this;
    piece.theta0 = theta0 + t0 * (theta1 - theta0);
    piece.theta1 = theta0 + t1 * (theta1 - theta0);
    return piece;
  }
  // Angle of a point's in-plane direction from the centre, reduced into the range
  // [theta_min - pi, theta_min + pi) about the middle of the angular range.
  double AngleOf(const Point3D &p) const
  {
    const Point3D d = Sub(p, center);
    double theta = std::atan2(Dot(d, v), Dot(d, u));
    const double middle = 0.5 * (theta0 + theta1);
    const double two_pi = 2.0 * std::acos(-1.0);
    while (theta < middle - std::acos(-1.0))
    {
      theta += two_pi;
    }
    while (theta >= middle + std::acos(-1.0))
    {
      theta -= two_pi;
    }
    return theta;
  }
  bool Contains(double theta) const
  {
    return theta >= std::min(theta0, theta1) && theta <= std::max(theta0, theta1);
  }
};

// Exact distance from a point to an arc piece: the radial distance in the arc plane where
// the point's angle falls inside the range, else the distance to the nearer end.
double PointArcDistance(const Point3D &p, const ArcPiece &arc)
{
  const Point3D d = Sub(p, arc.center);
  const double x = Dot(d, arc.u), y = Dot(d, arc.v);
  const double in_plane = std::hypot(x, y);
  const Point3D n = arc.Normal();
  const double off_plane = Dot(d, n);
  if (in_plane > 0.0)
  {
    const double theta = arc.AngleOf(p);
    if (arc.Contains(theta))
    {
      return std::hypot(in_plane - arc.radius, off_plane);
    }
  }
  return std::min(Distance(p, arc.At(arc.theta0)), Distance(p, arc.At(arc.theta1)));
}

// A straight segment or an arc piece (the geometry of a run interval or of a claimed core).
struct CurvePiece
{
  Point3D a{}, b{};
  std::optional<ArcPiece> arc;

  Point3D At(double t) const  // t in [0, 1]
  {
    return arc ? arc->AtFraction(t) : Add(a, Scale(t, Sub(b, a)));
  }
  double Length() const { return arc ? arc->Length() : Distance(a, b); }
  // Every point of an arc piece lies within its sagitta of its chord: the chord distance
  // minus the sagitta is a lower bound of the distance to the arc.
  double Sagitta() const
  {
    return arc ? arc->radius *
                     (1.0 - std::cos(0.5 * std::min(arc->Sweep(), std::acos(-1.0))))
               : 0.0;
  }
};

double PointPieceDistance(const Point3D &p, const CurvePiece &piece)
{
  return piece.arc ? PointArcDistance(p, *piece.arc)
                   : PointSegmentDistance(p, piece.a, piece.b);
}

// Minimum of a continuous function of one parameter sampled every `step` over [0, 1] in
// parameter length `length` (at least kArcSampleMinimumIntervals intervals), every local
// minimum of the samples refined by golden section.
double SampledMinimum(const std::function<double(double)> &f, double length, double step)
{
  const std::size_t n = std::max<std::size_t>(
      kArcSampleMinimumIntervals,
      static_cast<std::size_t>(std::ceil(length / std::max(step, 1.0e-300))));
  std::vector<double> values(n + 1);
  for (std::size_t k = 0; k <= n; k++)
  {
    values[k] = f(static_cast<double>(k) / static_cast<double>(n));
  }
  double best = std::numeric_limits<double>::infinity();
  for (std::size_t k = 0; k <= n; k++)
  {
    const bool left_ok = k == 0 || !(values[k - 1] < values[k]);
    const bool right_ok = k == n || !(values[k + 1] < values[k]);
    if (left_ok && right_ok)
    {
      const double a = static_cast<double>(k > 0 ? k - 1 : 0) / static_cast<double>(n);
      const double b = static_cast<double>(std::min(k + 1, n)) / static_cast<double>(n);
      best = std::min({best, values[k], GoldenMinimum(f, a, b).second});
    }
  }
  return best;
}

// Minimal distance between two pieces: the exact segment formula for two segments; with an
// arc the arc side is sampled (the distance from a point of one piece to the other piece is
// exact) at `step` and every local minimum refined, so that the result is the geometric
// minimum to the golden-section precision.
double PiecePieceDistance(const CurvePiece &A, const CurvePiece &B, double step)
{
  if (!A.arc && !B.arc)
  {
    return SegmentSegmentDistance(A.a, A.b, B.a, B.b);
  }
  const CurvePiece &sampled = A.arc ? A : B;
  const CurvePiece &other = A.arc ? B : A;
  return SampledMinimum([&](double t) { return PointPieceDistance(sampled.At(t), other); },
                        sampled.Length(), step);
}

// ---------------------------------------------------------------------------------------
// SHA-256 (FIPS 180-4), for the stable feature hashes and the geometry digest
// ---------------------------------------------------------------------------------------

std::string Sha256HexImpl(const std::string &text)
{
  static constexpr std::uint32_t k[64] = {
      0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4,
      0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe,
      0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f,
      0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
      0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc,
      0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
      0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116,
      0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
      0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7,
      0xc67178f2};
  std::uint32_t h[8] = {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                        0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
  std::string message = text;
  const std::uint64_t bit_length = static_cast<std::uint64_t>(text.size()) * 8;
  message.push_back(static_cast<char>(0x80));
  while (message.size() % 64 != 56)
  {
    message.push_back(static_cast<char>(0));
  }
  for (int i = 7; i >= 0; i--)
  {
    message.push_back(static_cast<char>((bit_length >> (8 * i)) & 0xff));
  }
  auto rotr = [](std::uint32_t x, int n) { return (x >> n) | (x << (32 - n)); };
  for (std::size_t offset = 0; offset < message.size(); offset += 64)
  {
    std::uint32_t w[64];
    for (int i = 0; i < 16; i++)
    {
      w[i] = 0;
      for (int j = 0; j < 4; j++)
      {
        w[i] = (w[i] << 8) | static_cast<unsigned char>(message[offset + 4 * i + j]);
      }
    }
    for (int i = 16; i < 64; i++)
    {
      const std::uint32_t s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
      const std::uint32_t s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
      w[i] = w[i - 16] + s0 + w[i - 7] + s1;
    }
    std::uint32_t a = h[0], b = h[1], c = h[2], d = h[3], e = h[4], f = h[5], g = h[6],
                  hh = h[7];
    for (int i = 0; i < 64; i++)
    {
      const std::uint32_t S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
      const std::uint32_t ch = (e & f) ^ (~e & g);
      const std::uint32_t t1 = hh + S1 + ch + k[i] + w[i];
      const std::uint32_t S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
      const std::uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
      const std::uint32_t t2 = S0 + maj;
      hh = g;
      g = f;
      f = e;
      e = d + t1;
      d = c;
      c = b;
      b = a;
      a = t1 + t2;
    }
    h[0] += a;
    h[1] += b;
    h[2] += c;
    h[3] += d;
    h[4] += e;
    h[5] += f;
    h[6] += g;
    h[7] += hh;
  }
  std::ostringstream out;
  out << std::hex << std::setfill('0');
  for (const std::uint32_t value : h)
  {
    out << std::setw(8) << value;
  }
  return out.str();
}

// ---------------------------------------------------------------------------------------
// Straight runs of the perimeter chains
// ---------------------------------------------------------------------------------------

struct RunSegment
{
  std::size_t segment = 0;
  bool forward = true;  // p0 -> p1 of the mesh segment follows the run tangent
  double t0 = 0.0;
  double t1 = 0.0;
};

struct Run
{
  int chain = -1;
  std::size_t index_in_chain = 0;
  Point3D start{};
  Point3D end{};
  Point3D tangent{};
  double length = 0.0;
  int conductor = 0;
  std::map<InterfaceDielectric, int> targets;
  Point3D gap_direction{};
  Point3D process_normal{};
  std::string boundary_law;
  std::vector<RunSegment> segments;
  std::size_t start_vertex = 0;
  std::size_t end_vertex = 0;
  bool excluded = false;

  Point3D At(double s) const { return Add(start, Scale(s, tangent)); }
};

struct Chain
{
  int id = -1;
  std::vector<std::size_t> runs;  // ordered along the chain
  bool closed = false;
  std::set<std::size_t> vertices;

  // Arc-length parameterisation and the curved-edge chain rule (filled by
  // Identifier::ComputeCurvature): run_offset[k] is the chain position of run k's start;
  // joint_turn[k] the turn (radians) at the joint before run k (0 at an open chain's start
  // and where a rounded corner accounts for the joint); kappa_nodes the piecewise-linear
  // windowed curvature; curved the chain intervals whose windowed bend radius is below the
  // straight threshold, with the maximal windowed curvature of each.
  std::vector<double> run_offset;
  double length = 0.0;
  std::vector<double> joint_turn;
  std::vector<bool> joint_excluded;
  std::vector<std::pair<double, double>> kappa_nodes;
  std::vector<Interval> curved;
  std::vector<double> curved_max_kappa;
  // The signed windowed curvature at the same nodes: positive where the chain turns toward
  // its metal (a convex bend: metal inside the bend, the edge of a disk), negative toward
  // the gap (concave, the edge of a hole). Its integral over a portion is the portion's
  // signed turn toward the metal (radians), the weight of the first-order curvature term.
  std::vector<double> signed_kappa_nodes;
  // Fitted bend arcs on this chain: (arc, x0, x1) of the arc runs (Identifier::arcs).
  std::vector<std::tuple<int, double, double>> arc_spans;

  bool Rigid() const { return runs.size() == 1; }
};

std::vector<std::string> InterfaceNames(const std::map<InterfaceDielectric, int> &targets)
{
  std::vector<std::string> names;
  for (const auto &[type, index] : targets)
  {
    (void)index;
    names.push_back(ToString(type));
  }
  return names;
}

bool SameFrame(const IdentificationSegment &a, const IdentificationSegment &b)
{
  return a.conductor == b.conductor && a.targets == b.targets &&
         a.boundary_law == b.boundary_law &&
         !DirectionLess(Dot(a.gap_direction, b.gap_direction),
                        1.0 - kParallelCosineTolerance);
}

// Order the segments of every chain along the chain (canonical start: the lexicographically
// smallest end vertex of an open chain, the smallest vertex of a loop, heading towards the
// smaller neighbour) and split the chain into straight runs where consecutive segments are
// not collinear on the direction grid or change conductor / interface / law / gap side.
void BuildRuns(const IdentificationInput &input, std::vector<Run> &runs,
               std::vector<Chain> &chains, std::map<int, std::size_t> &chain_index)
{
  std::map<int, std::vector<std::size_t>> segments_by_chain;
  for (std::size_t i = 0; i < input.segments.size(); i++)
  {
    const auto &segment = input.segments[i];
    if (segment.truncation || segment.exclusion || segment.chain < 0)
    {
      continue;
    }
    segments_by_chain[segment.chain].push_back(i);
  }
  // A chain interrupted by excluded segments (untargeted, undetermined process side, ...)
  // continues as separate chains on either side of the exclusion; the cut vertices are
  // ExclusionCut entries of the vertex table, not endpoints.
  int next_chain_id = segments_by_chain.empty() ? 0 : segments_by_chain.rbegin()->first + 1;
  std::map<int, std::vector<std::size_t>> split_chains;
  for (const auto &[chain_id, members] : segments_by_chain)
  {
    std::map<std::size_t, std::size_t> member_index;
    for (std::size_t k = 0; k < members.size(); k++)
    {
      member_index.emplace(members[k], k);
    }
    std::vector<std::size_t> parent(members.size());
    std::iota(parent.begin(), parent.end(), 0);
    auto Find = [&](std::size_t i)
    {
      while (parent[i] != i)
      {
        parent[i] = parent[parent[i]];
        i = parent[i];
      }
      return i;
    };
    std::map<std::size_t, std::size_t> first_member_at_vertex;
    for (std::size_t k = 0; k < members.size(); k++)
    {
      for (const std::size_t v : input.segments[members[k]].vertices)
      {
        auto [entry, inserted] = first_member_at_vertex.try_emplace(v, k);
        if (!inserted)
        {
          const std::size_t a = Find(entry->second), b = Find(k);
          if (a != b)
          {
            parent[std::max(a, b)] = std::min(a, b);
          }
        }
      }
    }
    std::map<std::size_t, int> id_by_root;
    for (std::size_t k = 0; k < members.size(); k++)
    {
      const std::size_t root = Find(k);
      auto [entry, inserted] =
          id_by_root.try_emplace(root, id_by_root.empty() ? chain_id : next_chain_id);
      if (inserted && entry->second == next_chain_id)
      {
        next_chain_id++;
      }
      split_chains[entry->second].push_back(members[k]);
    }
  }
  auto Less = [](const Point3D &a, const Point3D &b) { return a < b; };
  for (const auto &[chain_id, members] : split_chains)
  {
    std::map<std::size_t, std::vector<std::size_t>> incident;
    for (const std::size_t s : members)
    {
      for (const std::size_t v : input.segments[s].vertices)
      {
        incident[v].push_back(s);
      }
    }
    std::vector<std::size_t> ends;
    for (const auto &[v, list] : incident)
    {
      if (list.size() == 1)
      {
        ends.push_back(v);
      }
      MFEM_VERIFY(list.size() <= 2,
                  "A perimeter chain has a vertex with more than two chain segments!");
    }
    Chain chain;
    chain.id = chain_id;
    chain.closed = ends.empty();
    std::size_t start_vertex;
    if (!chain.closed)
    {
      MFEM_VERIFY(ends.size() == 2, "An open perimeter chain must have two ends!");
      start_vertex =
          Less(input.vertices[ends[0]].coordinate, input.vertices[ends[1]].coordinate)
              ? ends[0]
              : ends[1];
    }
    else
    {
      start_vertex = incident.begin()->first;
      for (const auto &[v, list] : incident)
      {
        if (Less(input.vertices[v].coordinate, input.vertices[start_vertex].coordinate))
        {
          start_vertex = v;
        }
      }
    }
    // Walk.
    std::vector<std::pair<std::size_t, bool>> ordered;  // (segment, forward)
    std::set<std::size_t> used;
    std::size_t current = start_vertex;
    std::optional<std::size_t> first_segment;
    if (chain.closed)
    {
      const auto &list = incident[start_vertex];
      auto Other = [&](std::size_t s)
      {
        const auto &v = input.segments[s].vertices;
        return v[0] == start_vertex ? v[1] : v[0];
      };
      first_segment = Less(input.vertices[Other(list[0])].coordinate,
                           input.vertices[Other(list[1])].coordinate)
                          ? list[0]
                          : list[1];
    }
    while (true)
    {
      std::optional<std::size_t> next;
      for (const std::size_t s : incident[current])
      {
        if (used.find(s) == used.end() &&
            (!first_segment || ordered.size() > 0 || s == *first_segment))
        {
          next = s;
          break;
        }
      }
      if (!next)
      {
        break;
      }
      used.insert(*next);
      const auto &v = input.segments[*next].vertices;
      const bool forward = v[0] == current;
      ordered.emplace_back(*next, forward);
      current = forward ? v[1] : v[0];
      chain.vertices.insert(current);
    }
    chain.vertices.insert(start_vertex);
    MFEM_VERIFY(ordered.size() == members.size(),
                "A perimeter chain is not a simple path or loop!");

    // Split into runs.
    std::size_t run_begin = 0;
    auto Direction = [&](std::size_t k)
    {
      const auto &segment = input.segments[ordered[k].first];
      const Point3D d =
          ordered[k].second ? Sub(segment.p1, segment.p0) : Sub(segment.p0, segment.p1);
      return Normalize(d);
    };
    auto Flush = [&](std::size_t begin, std::size_t end)
    {
      Run run;
      run.chain = chain_id;
      run.index_in_chain = chain.runs.size();
      const auto &first = input.segments[ordered[begin].first];
      const auto &last = input.segments[ordered[end - 1].first];
      run.start = ordered[begin].second ? first.p0 : first.p1;
      run.end = ordered[end - 1].second ? last.p1 : last.p0;
      run.start_vertex = ordered[begin].second ? first.vertices[0] : first.vertices[1];
      run.end_vertex = ordered[end - 1].second ? last.vertices[1] : last.vertices[0];
      run.length = Distance(run.start, run.end);
      run.tangent = Normalize(Sub(run.end, run.start));
      run.conductor = first.conductor;
      run.targets = first.targets;
      run.boundary_law = first.boundary_law;
      Point3D gap{}, normal{};
      double t = 0.0;
      for (std::size_t k = begin; k < end; k++)
      {
        const auto &segment = input.segments[ordered[k].first];
        const double length = Distance(segment.p0, segment.p1);
        run.segments.push_back({ordered[k].first, ordered[k].second, t, t + length});
        t += length;
        gap = Add(gap, Scale(length, segment.gap_direction));
        normal = Add(normal, Scale(length, segment.process_normal));
      }
      // Straight runs: the summed segment lengths equal the chord length to roundoff;
      // record the parameterisation on the chord so that s is a true arc length.
      const double ratio = run.length / t;
      for (auto &rs : run.segments)
      {
        rs.t0 *= ratio;
        rs.t1 *= ratio;
      }
      run.gap_direction = Normalize(gap);
      run.process_normal = Normalize(normal);
      chain.runs.push_back(runs.size());
      runs.push_back(std::move(run));
    };
    for (std::size_t k = 1; k < ordered.size(); k++)
    {
      const bool collinear =
          !DirectionLess(Dot(Direction(k - 1), Direction(k)), 1.0 - kDirectionQuantum);
      if (!collinear || !SameFrame(input.segments[ordered[k - 1].first],
                                   input.segments[ordered[k].first]))
      {
        Flush(run_begin, k);
        run_begin = k;
      }
    }
    Flush(run_begin, ordered.size());
    chain_index[chain_id] = chains.size();
    chains.push_back(std::move(chain));
  }
}

// ---------------------------------------------------------------------------------------
// Frames, planes and exclusions
// ---------------------------------------------------------------------------------------

Point3D SignCanonical(Point3D v)
{
  for (int d = 0; d < 3; d++)
  {
    if (std::abs(v[d]) > kDirectionQuantum)
    {
      return v[d] < 0.0 ? Scale(-1.0, v) : v;
    }
  }
  return v;
}

std::array<long long int, 3> DirectionKey(const Point3D &v, double quantum)
{
  return {static_cast<long long int>(std::llround(v[0] / quantum)),
          static_cast<long long int>(std::llround(v[1] / quantum)),
          static_cast<long long int>(std::llround(v[2] / quantum))};
}

Point3D ReferenceProcessNormal(const std::vector<Run> &runs)
{
  std::map<std::array<long long int, 3>, std::pair<double, Point3D>> buckets;
  for (const auto &run : runs)
  {
    const Point3D n = SignCanonical(run.process_normal);
    auto &bucket = buckets[DirectionKey(n, 1.0e-6)];
    bucket.first += run.length;
    bucket.second = Add(bucket.second, Scale(run.length, n));
  }
  MFEM_VERIFY(!buckets.empty(), "Geometry identification found no perimeter runs!");
  const auto best =
      std::max_element(buckets.begin(), buckets.end(),
                       [](const auto &a, const auto &b)
                       {
                         return a.second.first < b.second.first ||
                                (a.second.first == b.second.first && a.first > b.first);
                       });
  return Normalize(best->second.second);
}

// ---------------------------------------------------------------------------------------
// Canonical cluster signature
// ---------------------------------------------------------------------------------------

double RoundTo(double value, double quantum)
{
  const double r = std::round(value / quantum) * quantum;
  return r == 0.0 ? 0.0 : r;
}

std::array<double, 2> LocalCoordinates(const Point3D &p, const Point3D &origin,
                                       const Point3D &x, const Point3D &y, double radius)
{
  const Point3D r = Sub(p, origin);
  return {RoundTo(Dot(r, x) / radius, kSignatureLengthQuantumOverRadius),
          RoundTo(Dot(r, y) / radius, kSignatureLengthQuantumOverRadius)};
}

// The geometry entry of one portion in the frame (its sort key in the serialisation is the
// dump of this object; the conductor label is added after sorting).
// Point of a signature arc at the angle theta.
Point3D SignatureArcPoint(const SignatureArc &arc, double theta)
{
  return Add(arc.center, Add(Scale(arc.radius * std::cos(theta), arc.u),
                             Scale(arc.radius * std::sin(theta), arc.v)));
}

// Unit direction at the angle theta: the tangent in the direction of travel (theta0 ->
// theta1) or the outward radial direction.
Point3D SignatureArcTangent(const SignatureArc &arc, double theta)
{
  const double sign = arc.theta1 >= arc.theta0 ? 1.0 : -1.0;
  return Scale(sign, Add(Scale(-std::sin(theta), arc.u), Scale(std::cos(theta), arc.v)));
}

Point3D SignatureArcRadial(const SignatureArc &arc, double theta)
{
  return Add(Scale(std::cos(theta), arc.u), Scale(std::sin(theta), arc.v));
}

double SignatureArcLength(const SignatureArc &arc)
{
  return arc.radius * std::abs(arc.theta1 - arc.theta0);
}

// Centroid of the arc (length weighted): centre + radius sin(sweep / 2) / (sweep / 2) along
// the middle radial direction.
Point3D SignatureArcCentroid(const SignatureArc &arc)
{
  const double half = 0.5 * std::abs(arc.theta1 - arc.theta0);
  const double factor = half > 1.0e-12 ? std::sin(half) / half : 1.0;
  return Add(arc.center, Scale(arc.radius * factor,
                               SignatureArcRadial(arc, 0.5 * (arc.theta0 + arc.theta1))));
}

// The geometry entry of one portion in the frame (its sort key in the serialisation is the
// dump of this object; the conductor label is added after sorting). An arc portion is the
// sorted pair of its end points, its centre and its midpoint (the point midway along the
// arc: with the ends it fixes the circle and the side; a closed circle has equal ends and
// the midpoint at the antipode) and its radial gap sign, so that a re-meshed arc with other
// chords serialises identically.
nlohmann::json PortionGeometryInFrame(const SignaturePortion &portion,
                                      const Point3D &origin, const Point3D &x,
                                      const Point3D &y, double radius)
{
  auto a = LocalCoordinates(portion.p0, origin, x, y, radius),
       b = LocalCoordinates(portion.p1, origin, x, y, radius);
  if (b < a)
  {
    std::swap(a, b);
  }
  if (portion.arc)
  {
    const SignatureArc &arc = *portion.arc;
    const auto c = LocalCoordinates(arc.center, origin, x, y, radius);
    auto m = LocalCoordinates(SignatureArcPoint(arc, 0.5 * (arc.theta0 + arc.theta1)),
                              origin, x, y, radius);
    if (std::abs(arc.theta1 - arc.theta0) >= 2.0 * std::acos(-1.0) - 1.0e-9)
    {
      // A closed circle has no ends of its own (the mesh's first joint is not geometry):
      // serialised as centre + radius, with the ends at the frame's +x point of the circle
      // and the midpoint at its antipode, so every chording / start vertex agrees.
      const double r = RoundTo(arc.radius / radius, kSignatureLengthQuantumOverRadius);
      a = {c[0] + r, c[1]};
      b = a;
      m = {c[0] - r, c[1]};
    }
    return nlohmann::json{{"P", {a[0], a[1], b[0], b[1]}},
                          {"Arc", {c[0], c[1], m[0], m[1]}},
                          {"GapRadial", arc.gap_radial},
                          {"Interfaces", portion.interfaces},
                          {"Law", portion.boundary_law}};
  }
  return nlohmann::json{
      {"P", {a[0], a[1], b[0], b[1]}},
      {"Gap",
       std::array<double, 2>{
           RoundTo(Dot(portion.gap_direction, x), kSignatureLengthQuantumOverRadius),
           RoundTo(Dot(portion.gap_direction, y), kSignatureLengthQuantumOverRadius)}},
      {"Interfaces", portion.interfaces},
      {"Law", portion.boundary_law}};
}

// The smallest portion entry key of a frame: the serialisation in that frame begins with
// this entry (labelled conductor 1), so two frames whose smallest keys differ compare as
// their smallest keys do (two dumps of complete JSON objects differ at a position inside
// both) — an exact filter on the candidate frames of CanonicalClusterSignature.
std::string SmallestPortionKeyInFrame(const std::vector<SignaturePortion> &portions,
                                      const Point3D &origin, const Point3D &x,
                                      const Point3D &y, double radius)
{
  std::string smallest;
  bool have = false;
  for (const auto &portion : portions)
  {
    std::string key = PortionGeometryInFrame(portion, origin, x, y, radius).dump();
    if (!have || key < smallest)
    {
      smallest = std::move(key);
      have = true;
    }
  }
  return smallest;
}

nlohmann::json SerializeInFrame(const std::vector<SignaturePortion> &portions,
                                const std::vector<SignatureVertex> &vertices,
                                const Point3D &origin, const Point3D &x, const Point3D &y,
                                double radius)
{
  auto Local = [&](const Point3D &p) { return LocalCoordinates(p, origin, x, y, radius); };
  struct Entry
  {
    nlohmann::json geometry;
    std::string key;  // geometry.dump(), the sort key (computed once per entry)
    int conductor;
  };
  std::vector<Entry> entries;
  for (const auto &portion : portions)
  {
    nlohmann::json geometry = PortionGeometryInFrame(portion, origin, x, y, radius);
    std::string key = geometry.dump();
    entries.push_back({std::move(geometry), std::move(key), portion.conductor});
  }
  std::sort(entries.begin(), entries.end(),
            [](const Entry &a, const Entry &b) { return a.key < b.key; });
  std::map<int, int> labels;
  nlohmann::json portion_list = nlohmann::json::array();
  for (auto &entry : entries)
  {
    const auto [it, inserted] =
        labels.emplace(entry.conductor, static_cast<int>(labels.size()) + 1);
    (void)inserted;
    entry.geometry["Conductor"] = it->second;
    portion_list.push_back(entry.geometry);
  }
  std::vector<std::string> vertex_keys;
  nlohmann::json vertex_list = nlohmann::json::array();
  std::vector<nlohmann::json> vertex_entries;
  for (const auto &vertex : vertices)
  {
    const auto p = Local(vertex.point);
    vertex_entries.push_back(
        {{"P", {p[0], p[1]}},
         {"Type", vertex.type},
         {"TurnDegrees", RoundTo(vertex.turn_degrees, kSignatureAngleQuantumDegrees)}});
  }
  std::sort(vertex_entries.begin(), vertex_entries.end(),
            [](const auto &a, const auto &b) { return a.dump() < b.dump(); });
  for (auto &entry : vertex_entries)
  {
    vertex_list.push_back(std::move(entry));
  }
  return {{"Portions", portion_list}, {"Vertices", vertex_list}};
}

}  // namespace

std::string Sha256Hex(const std::string &text)
{
  return Sha256HexImpl(text);
}

std::pair<std::string, std::string> SignatureKeyAndHash(nlohmann::json signature,
                                                        const std::string &type)
{
  signature["Type"] = type;
  const std::string key = signature.dump();
  return {key, Sha256HexImpl(key)};
}

double DistanceToSerializedPortion(const nlohmann::json &portion,
                                   const std::array<double, 2> &q)
{
  const auto P = portion.at("P").get<std::array<double, 4>>();
  const std::array<double, 2> a = {P[0], P[1]}, b = {P[2], P[3]};
  if (!portion.contains("Arc"))
  {
    const double dx = b[0] - a[0], dy = b[1] - a[1];
    const double length2 = dx * dx + dy * dy;
    double s = length2 > 0.0 ? ((q[0] - a[0]) * dx + (q[1] - a[1]) * dy) / length2 : 0.0;
    s = std::clamp(s, 0.0, 1.0);
    return std::hypot(q[0] - (a[0] + s * dx), q[1] - (a[1] + s * dy));
  }
  const auto arc = portion.at("Arc").get<std::array<double, 4>>();
  const std::array<double, 2> c = {arc[0], arc[1]}, m = {arc[2], arc[3]};
  const double r = std::hypot(a[0] - c[0], a[1] - c[1]);
  const double radial = std::abs(std::hypot(q[0] - c[0], q[1] - c[1]) - r);
  const bool closed = std::hypot(a[0] - b[0], a[1] - b[1]) <= 1.0e-9 * std::max(r, 1.0);
  if (closed)
  {
    return radial;
  }
  const double two_pi = 2.0 * std::acos(-1.0);
  auto Angle = [&](const std::array<double, 2> &p)
  { return std::atan2(p[1] - c[1], p[0] - c[0]); };
  const double ta = Angle(a);
  // The arc runs from a to b through m: counterclockwise when m lies on the
  // counterclockwise sweep from a to b, clockwise otherwise.
  const double ccw = std::fmod(Angle(b) - ta + two_pi, two_pi);
  const bool counterclockwise = std::fmod(Angle(m) - ta + two_pi, two_pi) <= ccw + 1.0e-12;
  const double sweep = counterclockwise ? ccw : ccw - two_pi;
  double tq = std::fmod(Angle(q) - ta + two_pi, two_pi);
  if (!counterclockwise)
  {
    tq = tq > 0.0 ? tq - two_pi : tq;
  }
  const bool within = counterclockwise ? tq <= sweep + 1.0e-12 : tq >= sweep - 1.0e-12;
  if (within)
  {
    return radial;
  }
  return std::min(std::hypot(q[0] - a[0], q[1] - a[1]),
                  std::hypot(q[0] - b[0], q[1] - b[1]));
}

nlohmann::json FrameToJson(const std::array<double, 3> &origin,
                           const std::array<std::array<double, 3>, 3> &axes, double radius,
                           double length_scale)
{
  const double scaled_radius = radius * length_scale;
  const double tolerance =
      std::max(1.0e-10 * scaled_radius, 64.0 * std::numeric_limits<double>::epsilon());
  const double step = std::pow(10.0, std::floor(std::log10(tolerance)));
  auto L = [&](double value)
  {
    const double snapped = std::round(value * length_scale / step) * step;
    return snapped == 0.0 ? 0.0 : snapped;
  };
  auto D = [](const std::array<double, 3> &v)
  {
    auto S = [](double x)
    {
      const double s = std::round(x * 1.0e12) * 1.0e-12;
      return s == 0.0 ? 0.0 : s;
    };
    return nlohmann::json{S(v[0]), S(v[1]), S(v[2])};
  };
  return {{"Origin", {L(origin[0]), L(origin[1]), L(origin[2])}},
          {"Axes", {D(axes[0]), D(axes[1]), D(axes[2])}}};
}

namespace
{

// Continuous entries of the signatures (design (b)): lengths in units of R and angles in
// degrees; every other entry is topology.
bool IsLengthParameter(const std::string &key)
{
  return key == "OffsetOverR" || key == "SeparationOverR" || key == "RadiusOverR" ||
         key == "CornerRadiusOverR";
}

bool IsAngleParameter(const std::string &key)
{
  return key == "AngleDegrees" || key == "ArmAnglesDegrees";
}

void SplitParameters(nlohmann::json &node, SignatureParameters &out)
{
  if (node.is_object())
  {
    for (auto it = node.begin(); it != node.end(); ++it)
    {
      const bool length = IsLengthParameter(it.key());
      const bool angle = IsAngleParameter(it.key());
      if (length || angle)
      {
        auto &list = length ? out.lengths_over_R : out.angles_degrees;
        if (it.value().is_array())
        {
          for (const auto &v : it.value())
          {
            list.push_back(v.get<double>());
          }
        }
        else
        {
          list.push_back(it.value().get<double>());
        }
        it.value() = nullptr;
      }
      else
      {
        SplitParameters(it.value(), out);
      }
    }
  }
  else if (node.is_array())
  {
    for (auto &v : node)
    {
      SplitParameters(v, out);
    }
  }
}

}  // namespace

namespace
{

// The cluster branch of SplitSignatureParameters: every number of P / Arc / Gap / Box is a
// length on the 1e-6 R grid, every TurnDegrees an angle on the 1e-6 deg grid; everything
// else (Type, EdgeCount, Conductor, Interfaces, Law, GapRadial, Chain, Unboxable) is
// topology. The entry ORDER of the Portions / Context / Vertices arrays is part of the
// topology key (the canonical serialisation sorts them; a permutation is another key).
bool IsClusterLengthKey(const std::string &key)
{
  return key == "P" || key == "Arc" || key == "Gap" || key == "Box";
}

bool IsClusterAngleKey(const std::string &key)
{
  return key == "TurnDegrees";
}

void SplitClusterParameters(nlohmann::json &node, const std::string &path,
                            SignatureParameters &out)
{
  if (node.is_object())
  {
    for (auto it = node.begin(); it != node.end(); ++it)
    {
      const bool length = IsClusterLengthKey(it.key());
      const bool angle = IsClusterAngleKey(it.key());
      const std::string entry = path.empty() ? it.key() : path + "." + it.key();
      if (length || angle)
      {
        auto &list = length ? out.lengths_over_R : out.angles_degrees;
        auto &paths = length ? out.length_paths : out.angle_paths;
        if (it.value().is_array())
        {
          std::size_t k = 0;
          for (const auto &v : it.value())
          {
            MFEM_VERIFY(v.is_number(),
                        "Cluster signature entry " << entry << " is not numeric!");
            list.push_back(v.get<double>());
            paths.push_back(entry + "[" + std::to_string(k++) + "]");
          }
        }
        else
        {
          MFEM_VERIFY(it.value().is_number(),
                      "Cluster signature entry " << entry << " is not numeric!");
          list.push_back(it.value().get<double>());
          paths.push_back(entry);
        }
        it.value() = nullptr;
      }
      else
      {
        SplitClusterParameters(it.value(), entry, out);
      }
    }
  }
  else if (node.is_array())
  {
    std::size_t k = 0;
    for (auto &v : node)
    {
      SplitClusterParameters(v, path + "[" + std::to_string(k++) + "]", out);
    }
  }
}

bool IsClusterSignature(const nlohmann::json &signature)
{
  return signature.is_object() && signature.contains("Type") &&
         signature["Type"] == "SpatialEdgeCluster";
}

}  // namespace

SignatureParameters SplitSignatureParameters(const nlohmann::json &signature)
{
  SignatureParameters out;
  if (IsClusterSignature(signature))
  {
    nlohmann::json topology = signature;
    SplitClusterParameters(topology, "", out);
    out.topology_key = topology.dump();
    return out;
  }
  nlohmann::json topology = signature;
  SplitParameters(topology, out);
  out.topology_key = topology.dump();
  return out;
}

std::optional<ClusterSignatureDifference>
ClusterSignatureQuantumDifference(const nlohmann::json &a, const nlohmann::json &b)
{
  const SignatureParameters pa = SplitSignatureParameters(a),
                            pb = SplitSignatureParameters(b);
  if (pa.topology_key != pb.topology_key ||
      pa.lengths_over_R.size() != pb.lengths_over_R.size() ||
      pa.angles_degrees.size() != pb.angles_degrees.size())
  {
    return std::nullopt;
  }
  ClusterSignatureDifference out;
  for (std::size_t i = 0; i < pa.lengths_over_R.size(); i++)
  {
    const double quanta = std::abs(pa.lengths_over_R[i] - pb.lengths_over_R[i]) /
                          kSignatureLengthQuantumOverRadius;
    if (quanta > 0.5)  // a difference below half a quantum is the same grid value
    {
      out.differing_paths.push_back(pa.length_paths[i]);
    }
    out.max_delta_quanta = std::max(out.max_delta_quanta, quanta);
  }
  for (std::size_t i = 0; i < pa.angles_degrees.size(); i++)
  {
    const double quanta = std::abs(pa.angles_degrees[i] - pb.angles_degrees[i]) /
                          kSignatureAngleQuantumDegrees;
    if (quanta > 0.5)
    {
      out.differing_paths.push_back(pa.angle_paths[i]);
    }
    out.max_delta_quanta = std::max(out.max_delta_quanta, quanta);
  }
  return out;
}

void ValidateSpanCapAllowances(std::vector<SpanCapAllowance> &allowances)
{
  for (auto &allowance : allowances)
  {
    MFEM_VERIFY(
        allowance.claims_signature.is_object() &&
            allowance.claims_signature.contains("Portions") &&
            allowance.claims_signature["Portions"].is_array() &&
            !allowance.claims_signature["Portions"].empty(),
        "SpanCapAllowances entry \""
            << allowance.label
            << "\": ClaimsSignature must be a SpatialEdgeCluster claims-only "
               "signature object with Portions (copy SpatialSupport.ClaimsSignature "
               "of the refused cluster from the preflight inventory)!");
    const std::string type =
        allowance.claims_signature.value("Type", std::string("SpatialEdgeCluster"));
    MFEM_VERIFY(type == "SpatialEdgeCluster", "SpanCapAllowances entry \""
                                                  << allowance.label
                                                  << "\": ClaimsSignature Type must be "
                                                     "SpatialEdgeCluster, not "
                                                  << type << "!");
    allowance.claims_signature["Type"] = "SpatialEdgeCluster";
    // The Missing placeholder key of a refused cluster is its claims-only signature plus
    // "Unboxable": true (a manifest written before the ClaimsSignature export): accepted.
    allowance.claims_signature.erase("Unboxable");
    MFEM_VERIFY(!allowance.claims_signature.contains("Box") &&
                    !allowance.claims_signature.contains("Context"),
                "SpanCapAllowances entry \""
                    << allowance.label
                    << "\": ClaimsSignature must be the CLAIMS-ONLY signature (no Box / "
                       "Context: the allowance is resolved before any box)!");
    MFEM_VERIFY(allowance.span_cap_over_R >= kSupportSpanCapOverRadius,
                "SpanCapAllowances entry \""
                    << allowance.label << "\": SpanCapOverR " << allowance.span_cap_over_R
                    << " is below the plan span cap " << kSupportSpanCapOverRadius
                    << " (an allowance never lowers the cap)!");
  }
  for (std::size_t i = 0; i < allowances.size(); i++)
  {
    for (std::size_t j = i + 1; j < allowances.size(); j++)
    {
      const auto difference = ClusterSignatureQuantumDifference(
          allowances[i].claims_signature, allowances[j].claims_signature);
      MFEM_VERIFY(!difference || !ClusterQuantumDuplicate(difference->max_delta_quanta),
                  "SpanCapAllowances entries \""
                      << allowances[i].label << "\" and \"" << allowances[j].label
                      << "\" embed claims-only signatures of one topology within "
                      << (difference ? difference->max_delta_quanta : 0.0)
                      << " signature quanta of each other (<= "
                      << 2 * kClusterQuantumNearMatchMaxQuanta << " + "
                      << kClusterQuantumInclusiveMargin
                      << "): two allowances, one geometry (block (b) DESIGN A4 (4))!");
    }
  }
}

std::optional<ResolvedSpanCapAllowance>
ResolveSpanCapAllowance(const nlohmann::json &claims_signature,
                        const std::vector<SpanCapAllowance> &allowances)
{
  std::optional<ResolvedSpanCapAllowance> best;
  for (std::size_t a = 0; a < allowances.size(); a++)
  {
    const auto difference =
        ClusterSignatureQuantumDifference(claims_signature, allowances[a].claims_signature);
    if (!difference || !WithinClusterQuantumNearMatch(difference->max_delta_quanta))
    {
      continue;
    }
    if (!best || difference->max_delta_quanta < best->difference.max_delta_quanta ||
        (difference->max_delta_quanta == best->difference.max_delta_quanta &&
         allowances[a].label < allowances[best->index].label))
    {
      best = ResolvedSpanCapAllowance{a, *difference};
    }
  }
  return best;
}

nlohmann::json QuantumNearMatchRecord(const std::string &model_key,
                                      const std::string &feature_key,
                                      const ClusterSignatureDifference &difference)
{
  return nlohmann::json{
      {"ModelKey", model_key},
      {"FeatureKey", feature_key},
      {"MaxDeltaQuanta", difference.max_delta_quanta},
      {"DifferingNumbers",
       {{"Count", difference.differing_paths.size()},
        {"Paths", difference.differing_paths}}},
      {"MaxQuanta", kClusterQuantumNearMatchMaxQuanta},
      {"InclusiveMargin", kClusterQuantumInclusiveMargin},
      {"Rule", "block (b) DESIGN section 4 (decision 303): a SpatialEdgeCluster key whose "
               "topology (every entry with its numbers nulled, order preserved) equals the "
               "model's and whose numbers lie within MaxQuanta signature quanta (1e-6 R / "
               "1e-6 deg; read half-quantum inclusive: <= MaxQuanta + InclusiveMargin, "
               "decisions 287 / 288 / 317) of the model's is the same geometry at the "
               "grid: matched Exact "
               "to the model (placed in its own canonical frame with M = identity, the A10 "
               "checks unchanged); the FeatureKey is recorded beside the ModelKey"}};
}

nlohmann::json MirrorTranslationalSignature(const nlohmann::json &signature)
{
  if (!signature.is_object() || !signature.contains("Edges") ||
      !signature["Edges"].is_array() || signature["Edges"].empty())
  {
    return signature;
  }
  nlohmann::json mirror = signature;
  const auto &edges = signature["Edges"];
  const double span = edges.back()["OffsetOverR"].get<double>();
  nlohmann::json list = nlohmann::json::array();
  std::map<int, int> labels;
  for (auto it = edges.rbegin(); it != edges.rend(); ++it)
  {
    nlohmann::json edge = *it;
    edge["OffsetOverR"] = RoundTo(span - (*it)["OffsetOverR"].get<double>(),
                                  kSignatureLengthQuantumOverRadius);
    edge["GapSide"] = -(*it)["GapSide"].get<int>();
    const auto [label, inserted] =
        labels.emplace((*it)["Conductor"].get<int>(), static_cast<int>(labels.size()) + 1);
    (void)inserted;
    edge["Conductor"] = label->second;
    list.push_back(std::move(edge));
  }
  mirror["Edges"] = std::move(list);
  return mirror;
}

namespace
{

void SubstituteParameters(nlohmann::json &node, const std::vector<double> &lengths,
                          const std::vector<double> &angles, std::size_t &next_length,
                          std::size_t &next_angle)
{
  if (node.is_object())
  {
    for (auto it = node.begin(); it != node.end(); ++it)
    {
      const bool length = IsLengthParameter(it.key());
      const bool angle = IsAngleParameter(it.key());
      if (length || angle)
      {
        const auto &list = length ? lengths : angles;
        auto &next = length ? next_length : next_angle;
        const double quantum =
            length ? kSignatureLengthQuantumOverRadius : kSignatureAngleQuantumDegrees;
        if (it.value().is_array())
        {
          for (auto &v : it.value())
          {
            MFEM_VERIFY(next < list.size(),
                        "Signature parameter substitution out of range!");
            v = RoundTo(list[next++], quantum);
          }
        }
        else
        {
          MFEM_VERIFY(next < list.size(), "Signature parameter substitution out of range!");
          it.value() = RoundTo(list[next++], quantum);
        }
      }
      else
      {
        SubstituteParameters(it.value(), lengths, angles, next_length, next_angle);
      }
    }
  }
  else if (node.is_array())
  {
    for (auto &v : node)
    {
      SubstituteParameters(v, lengths, angles, next_length, next_angle);
    }
  }
}

}  // namespace

nlohmann::json SubstituteSignatureParameters(const nlohmann::json &signature,
                                             const std::vector<double> &lengths_over_R,
                                             const std::vector<double> &angles_degrees)
{
  nlohmann::json out = signature;
  if (IsClusterSignature(signature))
  {
    return out;
  }
  std::size_t next_length = 0, next_angle = 0;
  SubstituteParameters(out, lengths_over_R, angles_degrees, next_length, next_angle);
  MFEM_VERIFY(next_length == lengths_over_R.size() && next_angle == angles_degrees.size(),
              "Signature parameter substitution left parameters unused!");
  return out;
}

nlohmann::json RepresentativeSignature(const std::vector<nlohmann::json> &signatures)
{
  MFEM_VERIFY(!signatures.empty(),
              "A representative signature needs at least one instance!");
  // The lead: the lexicographically smallest serialisation (a set function).
  const nlohmann::json *lead = &signatures.front();
  std::string lead_key = lead->dump();
  for (const auto &signature : signatures)
  {
    const std::string key = signature.dump();
    if (key < lead_key)
    {
      lead_key = key;
      lead = &signature;
    }
  }
  const SignatureParameters lead_parameters = SplitSignatureParameters(*lead);
  std::vector<double> lo_length = lead_parameters.lengths_over_R,
                      hi_length = lead_parameters.lengths_over_R;
  std::vector<double> lo_angle = lead_parameters.angles_degrees,
                      hi_angle = lead_parameters.angles_degrees;
  for (const auto &signature : signatures)
  {
    // The orientation of this instance nearest to the lead.
    std::optional<double> best;
    SignatureParameters aligned;
    for (const nlohmann::json &candidate :
         {signature, MirrorTranslationalSignature(signature)})
    {
      const SignatureParameters parameters = SplitSignatureParameters(candidate);
      if (parameters.topology_key != lead_parameters.topology_key ||
          parameters.lengths_over_R.size() != lo_length.size() ||
          parameters.angles_degrees.size() != lo_angle.size())
      {
        continue;
      }
      double deviation = 0.0;
      for (std::size_t i = 0; i < lo_length.size(); i++)
      {
        deviation = std::max(deviation, std::abs(parameters.lengths_over_R[i] -
                                                 lead_parameters.lengths_over_R[i]));
      }
      for (std::size_t i = 0; i < lo_angle.size(); i++)
      {
        deviation = std::max(deviation, std::abs(parameters.angles_degrees[i] -
                                                 lead_parameters.angles_degrees[i]));
      }
      if (!best || deviation < *best)
      {
        best = deviation;
        aligned = parameters;
      }
    }
    MFEM_VERIFY(best.has_value(),
                "A representative signature over instances of different topologies!");
    for (std::size_t i = 0; i < lo_length.size(); i++)
    {
      lo_length[i] = std::min(lo_length[i], aligned.lengths_over_R[i]);
      hi_length[i] = std::max(hi_length[i], aligned.lengths_over_R[i]);
    }
    for (std::size_t i = 0; i < lo_angle.size(); i++)
    {
      lo_angle[i] = std::min(lo_angle[i], aligned.angles_degrees[i]);
      hi_angle[i] = std::max(hi_angle[i], aligned.angles_degrees[i]);
    }
  }
  std::vector<double> lengths(lo_length.size()), angles(lo_angle.size());
  for (std::size_t i = 0; i < lengths.size(); i++)
  {
    lengths[i] = 0.5 * (lo_length[i] + hi_length[i]);
  }
  for (std::size_t i = 0; i < angles.size(); i++)
  {
    angles[i] = 0.5 * (lo_angle[i] + hi_angle[i]);
  }
  return SubstituteSignatureParameters(*lead, lengths, angles);
}

std::optional<double> SignatureDeviation(const nlohmann::json &a, const nlohmann::json &b)
{
  if (IsClusterSignature(a) || IsClusterSignature(b))
  {
    // The quantum near-match (DESIGN section 4): no mirror orientation (the chirality is
    // folded into the canonical key), the tolerance k quanta read half-quantum inclusive
    // (decision 317 MINOR-1): <= 1 iff max |delta| <= k + 1/2 quanta.
    const auto difference = ClusterSignatureQuantumDifference(a, b);
    if (!difference)
    {
      return std::nullopt;
    }
    return difference->max_delta_quanta /
           (kClusterQuantumNearMatchMaxQuanta + kClusterQuantumInclusiveMargin);
  }
  const SignatureParameters pa = SplitSignatureParameters(a);
  std::optional<double> best;
  for (const nlohmann::json &candidate : {b, MirrorTranslationalSignature(b)})
  {
    const SignatureParameters pb = SplitSignatureParameters(candidate);
    if (pa.topology_key != pb.topology_key ||
        pa.lengths_over_R.size() != pb.lengths_over_R.size() ||
        pa.angles_degrees.size() != pb.angles_degrees.size())
    {
      continue;
    }
    double deviation = 0.0;
    for (std::size_t i = 0; i < pa.lengths_over_R.size(); i++)
    {
      deviation =
          std::max(deviation, std::abs(pa.lengths_over_R[i] - pb.lengths_over_R[i]) /
                                  kSignatureParameterToleranceOverRadius);
    }
    for (std::size_t i = 0; i < pa.angles_degrees.size(); i++)
    {
      deviation =
          std::max(deviation, std::abs(pa.angles_degrees[i] - pb.angles_degrees[i]) /
                                  kSignatureAngleToleranceDegrees);
    }
    if (!best || deviation < *best)
    {
      best = deviation;
    }
  }
  return best;
}

nlohmann::json CanonicalCornerSignature(const std::vector<std::string> &interfaces,
                                        const std::string &boundary_law,
                                        double angle_degrees, double corner_radius_over_R)
{
  return {{"Interfaces", interfaces},
          {"Law", boundary_law},
          {"AngleDegrees", RoundTo(angle_degrees, kSignatureAngleQuantumDegrees)},
          {"CornerRadiusOverR",
           RoundTo(corner_radius_over_R, kSignatureLengthQuantumOverRadius)}};
}

nlohmann::json CanonicalJunctionSignature(const std::vector<std::string> &interfaces,
                                          const std::string &boundary_law,
                                          std::vector<double> arm_angles_degrees,
                                          std::vector<int> arm_conductors,
                                          JunctionCanonicalOrder *order)
{
  for (double &angle : arm_angles_degrees)
  {
    angle = RoundTo(angle, kSignatureAngleQuantumDegrees);
  }
  if (arm_conductors.empty())
  {
    arm_conductors.assign(arm_angles_degrees.size(), 0);
  }
  MFEM_VERIFY(arm_conductors.size() == arm_angles_degrees.size(),
              "Junction arm conductors must correspond to the arm angles!");
  // Difference k is the gap after arm k (from arm k to arm k + 1). In the reversed
  // (mirror) orientation the gap after arm k is difference k - 1, so the reversed
  // difference sequence pairs with the reversed conductors rotated right by one.
  auto Relabel = [](const std::vector<int> &conductors)
  {
    std::vector<int> labels;
    std::map<int, int> first_appearance;
    for (const int conductor : conductors)
    {
      labels.push_back(first_appearance.try_emplace(conductor, first_appearance.size() + 1)
                           .first->second);
    }
    return labels;
  };
  // Canonical cyclic order: minimal (angles, conductor labels), both orientations (mirror).
  // Difference k is the gap after arm k in the angular order; in the reversed sequence the
  // element at position k is arm (n - k) mod n with the gap after it clockwise.
  const std::size_t n = arm_angles_degrees.size();
  std::pair<std::vector<double>, std::vector<int>> best{arm_angles_degrees,
                                                        Relabel(arm_conductors)};
  JunctionCanonicalOrder best_order;
  for (const bool reverse : {false, true})
  {
    std::vector<double> sequence = arm_angles_degrees;
    std::vector<int> conductors = arm_conductors;
    if (reverse)
    {
      std::reverse(sequence.begin(), sequence.end());
      std::reverse(conductors.begin(), conductors.end());
      std::rotate(conductors.begin(), conductors.end() - 1, conductors.end());
    }
    for (std::size_t shift = 0; shift < sequence.size(); shift++)
    {
      std::rotate(sequence.begin(), sequence.begin() + 1, sequence.end());
      std::rotate(conductors.begin(), conductors.begin() + 1, conductors.end());
      const std::pair<std::vector<double>, std::vector<int>> candidate{sequence,
                                                                       Relabel(conductors)};
      if (candidate < best)
      {
        best = candidate;
        const std::size_t position = (shift + 1) % n;
        best_order = {reverse ? (n - position) % n : position, reverse};
      }
    }
  }
  if (order)
  {
    *order = best_order;
  }
  return {{"Interfaces", interfaces},
          {"Law", boundary_law},
          {"ArmAnglesDegrees", best.first},
          {"ArmConductors", best.second}};
}

TranslationalSignature CanonicalTranslationalSignature(std::vector<TranslationalEdge> edges,
                                                       double radius)
{
  MFEM_VERIFY(edges.size() >= 2, "A translational signature needs at least two edges!");
  std::sort(edges.begin(), edges.end(),
            [](const TranslationalEdge &a, const TranslationalEdge &b)
            { return a.offset < b.offset; });
  TranslationalSignature best;
  std::string best_key;
  for (const int orientation : {1, -1})
  {
    std::vector<TranslationalEdge> ordered = edges;
    if (orientation < 0)
    {
      std::reverse(ordered.begin(), ordered.end());
    }
    const double w0 = ordered.front().offset;
    nlohmann::json list = nlohmann::json::array();
    std::map<int, int> labels;
    for (const auto &edge : ordered)
    {
      const auto [it, inserted] =
          labels.emplace(edge.conductor, static_cast<int>(labels.size()) + 1);
      (void)inserted;
      list.push_back({{"OffsetOverR", RoundTo(orientation * (edge.offset - w0) / radius,
                                              kSignatureLengthQuantumOverRadius)},
                      {"GapSide", orientation * edge.gap_sign},
                      {"Conductor", it->second},
                      {"Interfaces", edge.interfaces},
                      {"Law", edge.boundary_law}});
    }
    nlohmann::json candidate = {{"Edges", list}};
    if (edges.size() == 2)
    {
      candidate["SeparationOverR"] =
          RoundTo(std::abs(edges.back().offset - edges.front().offset) / radius,
                  kSignatureLengthQuantumOverRadius);
    }
    const std::string key = candidate.dump();
    if (best_key.empty() || key < best_key)
    {
      best_key = key;
      best.signature = std::move(candidate);
      best.chirality = orientation;
    }
    else if (key == best_key)
    {
      best.chirality = 0;  // symmetric under the lateral reflection
    }
  }
  return best;
}

namespace
{

// Origin (the length-weighted centroid of the portions, arcs by their analytic centroids)
// and the candidate in-plane x axes of a cluster frame: portion tangents and their
// perpendiculars, both signs; an arc portion contributes its end tangents and end radial
// directions (a closed circle none: its ends are a property of the mesh); closed circles
// only: the directions from the centroid to the circle centres, else any fixed in-plane
// axis. Shared by the claims-only (v2) and the Box + Context (v3) canonicalisations.
struct ClusterFrameCandidates
{
  Point3D origin{};
  std::vector<Point3D> axes;
};

ClusterFrameCandidates
ClusterFrameCandidatesOf(const std::vector<SignaturePortion> &portions, const Point3D &n)
{
  MFEM_VERIFY(!portions.empty(), "A cluster signature needs at least one edge portion!");
  ClusterFrameCandidates result;
  Point3D &origin = result.origin;
  double total = 0.0;
  for (const auto &portion : portions)
  {
    if (portion.arc)
    {
      const double length = SignatureArcLength(*portion.arc);
      origin = Add(origin, Scale(length, SignatureArcCentroid(*portion.arc)));
      total += length;
      continue;
    }
    const double length = Distance(portion.p0, portion.p1);
    origin = Add(origin, Scale(0.5 * length, Add(portion.p0, portion.p1)));
    total += length;
  }
  MFEM_VERIFY(total > 0.0, "A cluster signature needs portions of positive length!");
  origin = Scale(1.0 / total, origin);
  std::vector<Point3D> &candidates = result.axes;
  std::set<std::array<long long int, 3>> seen;
  auto AddCandidate = [&](Point3D v)
  {
    v = Sub(v, Scale(Dot(v, n), n));
    if (Norm(v) <= 1.0e-12)
    {
      return;
    }
    v = Normalize(v);
    if (seen.insert(DirectionKey(v, 1.0e-9)).second)
    {
      candidates.push_back(v);
    }
  };
  for (const auto &portion : portions)
  {
    std::vector<Point3D> directions;
    if (portion.arc)
    {
      const SignatureArc &arc = *portion.arc;
      if (std::abs(arc.theta1 - arc.theta0) < 2.0 * std::acos(-1.0) - 1.0e-9)
      {
        for (const double theta : {arc.theta0, arc.theta1})
        {
          directions.push_back(SignatureArcTangent(arc, theta));
          directions.push_back(SignatureArcRadial(arc, theta));
        }
      }
    }
    else
    {
      directions.push_back(Normalize(Sub(portion.p1, portion.p0)));
    }
    for (const Point3D &t : directions)
    {
      for (const double sign : {1.0, -1.0})
      {
        AddCandidate(Scale(sign, t));
        AddCandidate(Scale(sign, Cross(n, t)));
      }
    }
  }
  if (candidates.empty())
  {
    for (const auto &portion : portions)
    {
      if (portion.arc)
      {
        AddCandidate(Sub(portion.arc->center, origin));
      }
    }
    if (candidates.empty())
    {
      int axis = 0;
      for (int d = 1; d < 3; d++)
      {
        axis = std::abs(n[d]) < std::abs(n[axis]) ? d : axis;
      }
      Point3D e{};
      e[axis] = 1.0;
      AddCandidate(e);
    }
  }
  return result;
}

// The context entries of a v3 serialisation in the frame: the portion encoding of every
// piece plus its ownership class, sorted by the geometry dump, conductor labels continuing
// the claims' label map (first appearance over Portions THEN Context).
nlohmann::json SerializeContextInFrame(const std::vector<SupportContextPiece> &context,
                                       const Point3D &origin, const Point3D &x,
                                       const Point3D &y, double radius,
                                       std::map<int, int> &labels)
{
  struct Entry
  {
    nlohmann::json geometry;
    std::string key;
    int conductor;
  };
  std::vector<Entry> entries;
  for (const auto &piece : context)
  {
    nlohmann::json geometry = PortionGeometryInFrame(piece.portion, origin, x, y, radius);
    geometry["Chain"] = piece.chain;
    std::string key = geometry.dump();
    entries.push_back({std::move(geometry), std::move(key), piece.portion.conductor});
  }
  std::sort(entries.begin(), entries.end(),
            [](const Entry &a, const Entry &b) { return a.key < b.key; });
  nlohmann::json list = nlohmann::json::array();
  for (auto &entry : entries)
  {
    const auto [it, inserted] =
        labels.emplace(entry.conductor, static_cast<int>(labels.size()) + 1);
    (void)inserted;
    entry.geometry["Conductor"] = it->second;
    list.push_back(std::move(entry.geometry));
  }
  return list;
}

// The conductor label map of SerializeInFrame (first appearance over the sorted portion
// entries), rebuilt from the serialised portions so that the context labels continue it.
std::map<int, int> PortionLabelMap(const std::vector<SignaturePortion> &portions,
                                   const Point3D &origin, const Point3D &x,
                                   const Point3D &y, double radius)
{
  std::vector<std::pair<std::string, int>> entries;
  for (const auto &portion : portions)
  {
    entries.emplace_back(PortionGeometryInFrame(portion, origin, x, y, radius).dump(),
                         portion.conductor);
  }
  std::sort(entries.begin(), entries.end());
  std::map<int, int> labels;
  for (const auto &[key, conductor] : entries)
  {
    (void)key;
    labels.emplace(conductor, static_cast<int>(labels.size()) + 1);
  }
  return labels;
}

nlohmann::json BoxJson(const std::array<double, 4> &box)
{
  return nlohmann::json{RoundTo(box[0], kSignatureLengthQuantumOverRadius),
                        RoundTo(box[1], kSignatureLengthQuantumOverRadius),
                        RoundTo(box[2], kSignatureLengthQuantumOverRadius),
                        RoundTo(box[3], kSignatureLengthQuantumOverRadius)};
}

}  // namespace

CanonicalSignature
CanonicalClusterSignature(const std::vector<SignaturePortion> &portions,
                          const std::vector<SignatureVertex> &vertices,
                          const Point3D &process_normal, double radius,
                          const std::function<void(std::size_t, std::size_t)> &progress)
{
  const Point3D n = Normalize(process_normal);
  const auto frame = ClusterFrameCandidatesOf(portions, n);
  const Point3D &origin = frame.origin;
  const std::vector<Point3D> &candidates = frame.axes;
  CanonicalSignature best;
  std::string best_key, best_smallest_portion;
  bool have = false;
  std::set<int> minimal_handedness;
  for (std::size_t c = 0; c < candidates.size(); c++)
  {
    const Point3D &x = candidates[c];
    if (progress)
    {
      progress(c, candidates.size());
    }
    for (const int handedness : {1, -1})
    {
      const Point3D y = Scale(static_cast<double>(handedness), Cross(n, x));
      // Exact filter: a frame whose smallest portion entry exceeds the best frame's cannot
      // serialise below (or equal to) the best key.
      std::string smallest_portion =
          SmallestPortionKeyInFrame(portions, origin, x, y, radius);
      if (have && smallest_portion > best_smallest_portion)
      {
        continue;
      }
      auto serialized = SerializeInFrame(portions, vertices, origin, x, y, radius);
      const std::string key = serialized.dump();
      if (!have || key < best_key)
      {
        have = true;
        best_key = key;
        best_smallest_portion = std::move(smallest_portion);
        best.signature = std::move(serialized);
        best.chirality = handedness;
        best.origin = origin;
        best.axes = {x, y, n};
        minimal_handedness = {handedness};
      }
      else if (key == best_key)
      {
        minimal_handedness.insert(handedness);
      }
    }
  }
  // A mirror-symmetric cluster reaches the minimal serialisation with both handedness
  // values: chirality 0 (its mirror image is itself).
  if (minimal_handedness.size() == 2)
  {
    best.chirality = 0;
  }
  best.key = best_key;
  best.hash = Sha256HexImpl(best_key);
  return best;
}

std::optional<CanonicalSignature> CanonicalClusterSignatureWithSupport(
    const std::vector<SignaturePortion> &portions,
    const std::vector<SignatureVertex> &vertices, const Point3D &process_normal,
    double radius,
    const std::function<std::optional<FrameSupport>(const Point3D &, const Point3D &,
                                                    const Point3D &)> &support,
    const std::function<void(std::size_t, std::size_t)> &progress)
{
  const Point3D n = Normalize(process_normal);
  const auto frame = ClusterFrameCandidatesOf(portions, n);
  const Point3D &origin = frame.origin;
  const std::vector<Point3D> &candidates = frame.axes;
  CanonicalSignature best;
  std::string best_key;
  bool have = false;
  std::set<int> minimal_handedness;
  for (std::size_t c = 0; c < candidates.size(); c++)
  {
    const Point3D &x = candidates[c];
    if (progress)
    {
      progress(c, candidates.size());
    }
    for (const int handedness : {1, -1})
    {
      const Point3D y = Scale(static_cast<double>(handedness), Cross(n, x));
      const auto frame_support = support(origin, x, y);
      if (!frame_support)
      {
        continue;  // unboxable in this frame
      }
      // The dump of a complete JSON object starts with its "Box" entry (keys sorted), so a
      // frame whose box dump exceeds the best frame's cannot serialise below it.
      nlohmann::json serialized =
          SerializeInFrame(portions, vertices, origin, x, y, radius);
      serialized["Box"] = BoxJson(frame_support->box);
      std::map<int, int> labels = PortionLabelMap(portions, origin, x, y, radius);
      serialized["Context"] =
          SerializeContextInFrame(frame_support->context, origin, x, y, radius, labels);
      const std::string key = serialized.dump();
      if (!have || key < best_key)
      {
        have = true;
        best_key = key;
        best.signature = std::move(serialized);
        best.chirality = handedness;
        best.origin = origin;
        best.axes = {x, y, n};
        minimal_handedness = {handedness};
      }
      else if (key == best_key)
      {
        minimal_handedness.insert(handedness);
      }
    }
  }
  if (!have)
  {
    return std::nullopt;
  }
  if (minimal_handedness.size() == 2)
  {
    best.chirality = 0;
  }
  best.key = best_key;
  best.hash = Sha256HexImpl(best_key);
  return best;
}

std::string SpatialSupportContextDigest(const nlohmann::json &signature)
{
  if (!signature.is_object() || !signature.contains("Box"))
  {
    return {};
  }
  const nlohmann::json context = {
      {"Box", signature.at("Box")},
      {"Context", signature.value("Context", nlohmann::json::array())}};
  return Sha256HexImpl(context.dump());
}

// The chords of one serialised portion (the builder's chording,
// signature_library.cluster_plan_view_edges): a straight portion is its own chord; an arc
// portion is cut into n equal chords from a to b through the midpoint m, n = max(ceil(sweep
// / step), ceil(arc length / max chord), 1), a closed circle (equal ends) over 2 pi, each
// chord's gap direction the arc's radial direction at the chord's middle.
std::vector<SerializedPortionChord> ChordSerializedPortion(const nlohmann::json &portion)
{
  using Point2 = std::array<double, 2>;
  std::vector<SerializedPortionChord> chords;
  const double pi = std::acos(-1.0);
  const auto &P = portion.at("P");
  const Point2 a = {P[0].get<double>(), P[1].get<double>()},
               b = {P[2].get<double>(), P[3].get<double>()};
  if (!portion.contains("Arc"))
  {
    const auto &G = portion.at("Gap");
    chords.push_back({{a[0], a[1], b[0], b[1]}, {G[0].get<double>(), G[1].get<double>()}});
    return chords;
  }
  const auto &A = portion.at("Arc");
  const Point2 c = {A[0].get<double>(), A[1].get<double>()},
               m = {A[2].get<double>(), A[3].get<double>()};
  const double r = std::hypot(a[0] - c[0], a[1] - c[1]);
  auto Angle = [&](const Point2 &q) { return std::atan2(q[1] - c[1], q[0] - c[0]); };
  const double ta = Angle(a), tb = Angle(b), tm = Angle(m);
  const bool closed = std::hypot(a[0] - b[0], a[1] - b[1]) <= 1.0e-9 * std::max(r, 1.0);
  double sweep;
  if (closed)
  {
    sweep = 2.0 * pi;
  }
  else
  {
    auto Mod = [&](double v)
    {
      double w = std::fmod(v, 2.0 * pi);
      return w < 0.0 ? w + 2.0 * pi : w;
    };
    const double ccw = Mod(tb - ta);
    sweep = Mod(tm - ta) <= ccw + 1.0e-12 ? ccw : ccw - 2.0 * pi;
  }
  const int n =
      std::max({static_cast<int>(std::ceil(
                    std::abs(sweep) / (kClusterArcChordStepDegrees * pi / 180.0) - 1.0e-9)),
                static_cast<int>(std::ceil(
                    r * std::abs(sweep) / kClusterArcChordMaxLengthOverRadius - 1.0e-9)),
                1});
  const double sign = portion.at("GapRadial").get<int>();
  for (int k = 0; k < n; k++)
  {
    const double t0 = ta + sweep * k / n, t1 = ta + sweep * (k + 1) / n,
                 tmid = 0.5 * (t0 + t1);
    chords.push_back({{c[0] + r * std::cos(t0), c[1] + r * std::sin(t0),
                       c[0] + r * std::cos(t1), c[1] + r * std::sin(t1)},
                      {sign * std::cos(tmid), sign * std::sin(tmid)}});
  }
  return chords;
}

std::vector<ContextPieceChord> ContextPieceChords(const nlohmann::json &signature)
{
  std::vector<ContextPieceChord> pieces;
  if (!signature.contains("Context"))
  {
    return pieces;
  }
  for (const auto &entry : signature.at("Context"))
  {
    const bool chain = entry.value("Chain", false);
    const int conductor = entry.value("Conductor", 0);
    for (const auto &chord : ChordSerializedPortion(entry))
    {
      pieces.push_back({chord.P, chain, conductor});
    }
  }
  return pieces;
}

// Rule B2 (the coupon generator's `coupon_bounds` / `edge_rows` / `extended_interval` in
// the canonical frame, units of R): see SupportBoxFromSignature in the header.
std::array<double, 4> SupportBoxFromSignature(const nlohmann::json &signature,
                                              std::size_t *band_hits)
{
  using Point2 = std::array<double, 2>;
  struct Edge
  {
    Point2 p0, p1, gap;
  };
  std::vector<Edge> edges;
  std::size_t hits = 0;
  // The box rule's two thresholds read within the knife-edge band (review MINOR-4): the
  // end coincidence and the "row end at or beyond R continues by 2R" rule.
  auto Band = [&](double value, double threshold)
  {
    if (std::abs(value / threshold - 1.0) <= kKnifeEdgeBandRelative)
    {
      hits++;
    }
  };
  const auto &portion_list = signature.at("Portions");
  MFEM_VERIFY(portion_list.is_array() && !portion_list.empty(),
              "A support box needs the claimed portions!");
  for (const auto &portion : portion_list)
  {
    for (const auto &chord : ChordSerializedPortion(portion))
    {
      edges.push_back({{chord.P[0], chord.P[1]}, {chord.P[2], chord.P[3]}, chord.gap});
    }
  }
  std::vector<Point2> vertex_points;
  if (signature.contains("Vertices"))
  {
    for (const auto &vertex : signature["Vertices"])
    {
      vertex_points.push_back({vertex["P"][0].get<double>(), vertex["P"][1].get<double>()});
    }
  }
  // End states (cluster_signature_geometry.end_states): an end touching a vertex or another
  // edge's end within the coincidence tolerance is connected; every other end is a claim
  // cut (free).
  auto Near = [](const Point2 &a, const Point2 &b, double tolerance)
  { return std::hypot(a[0] - b[0], a[1] - b[1]) <= tolerance; };
  const double coincidence = kSupportEndCoincidenceOverRadius;
  auto Connected = [&](std::size_t i, const Point2 &point)
  {
    for (const auto &v : vertex_points)
    {
      Band(std::hypot(point[0] - v[0], point[1] - v[1]), coincidence);
      if (Near(point, v, coincidence))
      {
        return true;
      }
    }
    for (std::size_t j = 0; j < edges.size(); j++)
    {
      if (j == i)
      {
        continue;
      }
      for (const Point2 &other : {edges[j].p0, edges[j].p1})
      {
        Band(std::hypot(point[0] - other[0], point[1] - other[1]), coincidence);
      }
      if (Near(point, edges[j].p0, coincidence) || Near(point, edges[j].p1, coincidence))
      {
        return true;
      }
    }
    return false;
  };
  double x0 = std::numeric_limits<double>::infinity(), y0 = x0, x1 = -x0, y1 = -x0;
  for (std::size_t i = 0; i < edges.size(); i++)
  {
    const Edge &edge = edges[i];
    const double gap_norm = std::hypot(edge.gap[0], edge.gap[1]);
    MFEM_VERIFY(gap_norm > 0.0, "A support box portion has no gap direction!");
    const Point2 gap = {edge.gap[0] / gap_norm, edge.gap[1] / gap_norm};
    const Point2 tangent = {gap[1], -gap[0]};  // gap x normal for the frame's +z normal
    const double length = std::hypot(edge.p1[0] - edge.p0[0], edge.p1[1] - edge.p0[1]);
    MFEM_VERIFY(length > 0.0, "A support box portion has zero length!");
    const Point2 midpoint = {0.5 * (edge.p0[0] + edge.p1[0]),
                             0.5 * (edge.p0[1] + edge.p1[1])};
    const bool begin_free = !Connected(i, edge.p0), end_free = !Connected(i, edge.p1);
    const bool forward_is_p1 =
        (edge.p1[0] - midpoint[0]) * tangent[0] + (edge.p1[1] - midpoint[1]) * tangent[1] >
        0.0;
    const bool end_is_free = forward_is_p1 ? end_free : begin_free;
    const bool begin_is_free = forward_is_p1 ? begin_free : end_free;
    const double half = 0.5 * length;
    double begin = begin_is_free ? -std::max(half, 1.0) : -half;
    double end = end_is_free ? std::max(half, 1.0) : half;
    // extended_interval: a row end at or beyond R from the row point continues by 2R.
    const double tolerance = 1.0e-10;
    // (A free end is lengthened to >= R and always continues: no threshold there.)
    if (!begin_is_free)
    {
      Band(half, 1.0);
    }
    if (!end_is_free)
    {
      Band(half, 1.0);
    }
    if (begin <= -1.0 + tolerance)
    {
      begin -= kSupportContinuationOverRadius;
    }
    if (end >= 1.0 - tolerance)
    {
      end += kSupportContinuationOverRadius;
    }
    for (const double coordinate : {begin, end})
    {
      for (const double side : {-1.0, 1.0})
      {
        const double px = midpoint[0] + coordinate * tangent[0] + side * gap[0];
        const double py = midpoint[1] + coordinate * tangent[1] + side * gap[1];
        x0 = std::min(x0, px);
        y0 = std::min(y0, py);
        x1 = std::max(x1, px);
        y1 = std::max(y1, py);
      }
    }
  }
  const double padding = kSupportPaddingOverRadius;
  if (band_hits)
  {
    *band_hits = hits;
  }
  return {RoundTo(x0 - padding, kSignatureLengthQuantumOverRadius),
          RoundTo(y0 - padding, kSignatureLengthQuantumOverRadius),
          RoundTo(x1 + padding, kSignatureLengthQuantumOverRadius),
          RoundTo(y1 + padding, kSignatureLengthQuantumOverRadius)};
}

// ---------------------------------------------------------------------------------------
// Identification
// ---------------------------------------------------------------------------------------

namespace
{

struct Claim
{
  int feature = -1;
  int priority = 0;  // 0 cluster, 1 vertex window, 2 translational, 3 isolated
  Interval interval;
  // Side of a pair / parallel cluster in its canonical lateral order (0 = the lowest offset
  // along the feature's lateral axis); 0 for every other feature.
  int side = 0;
};

struct EventCore
{
  std::size_t run = 0;
  Interval interval;
  Point3D p0{}, p1{};
  int plane = 0;  // cores of two metal planes never merge (decision 82(1))
  // The core's geometry: the run interval on its arc where the run lies on one (option A),
  // else the chord [p0, p1] (a site core is the point p0 = p1).
  CurvePiece piece;
};

struct UnionFind
{
  std::vector<std::size_t> parent;
  explicit UnionFind(std::size_t n) : parent(n)
  {
    std::iota(parent.begin(), parent.end(), 0);
  }
  std::size_t Find(std::size_t i)
  {
    while (parent[i] != i)
    {
      parent[i] = parent[parent[i]];
      i = parent[i];
    }
    return i;
  }
  void Union(std::size_t a, std::size_t b)
  {
    a = Find(a);
    b = Find(b);
    if (a != b)
    {
      parent[std::max(a, b)] = std::min(a, b);
    }
  }
};

// Uniform grid over axis-aligned boxes: every candidate search of the identification (runs
// against runs, faces, sites, cores) asks the grid for the items whose boxes may lie within
// a margin of a query box (an exact superset: an item is stored in every cell its box
// overlaps, the query visits every cell the enlarged box overlaps) and then applies the
// original geometric test to the candidates IN THE ORIGINAL ORDER (sorted indices), so that
// the results are identical to the former all-pairs loops. Pure acceleration; no rule lives
// here.
class UniformGrid
{
public:
  UniformGrid(double cell_size, const Point3D &lower) : cell(cell_size), origin(lower) {}

  void Insert(std::size_t id, const Point3D &lo, const Point3D &hi)
  {
    const auto c0 = Cell(lo), c1 = Cell(hi);
    for (int d = 0; d < 3; d++)
    {
      min_cell[d] = std::min(min_cell[d], c0[d]);
      max_cell[d] = std::max(max_cell[d], c1[d]);
    }
    for (long long int ix = c0[0]; ix <= c1[0]; ix++)
    {
      for (long long int iy = c0[1]; iy <= c1[1]; iy++)
      {
        for (long long int iz = c0[2]; iz <= c1[2]; iz++)
        {
          cells[Key(ix, iy, iz)].push_back(id);
        }
      }
    }
  }

  // Sorted, unique ids of the items stored in the cells overlapping [lo - margin, hi +
  // margin].
  std::vector<std::size_t> Query(const Point3D &lo, const Point3D &hi, double margin) const
  {
    std::vector<std::size_t> out;
    Point3D qlo = lo, qhi = hi;
    for (int d = 0; d < 3; d++)
    {
      qlo[d] -= margin;
      qhi[d] += margin;
    }
    const auto c0 = Cell(qlo), c1 = Cell(qhi);
    for (long long int ix = std::max(c0[0], min_cell[0]);
         ix <= std::min(c1[0], max_cell[0]); ix++)
    {
      for (long long int iy = std::max(c0[1], min_cell[1]);
           iy <= std::min(c1[1], max_cell[1]); iy++)
      {
        for (long long int iz = std::max(c0[2], min_cell[2]);
             iz <= std::min(c1[2], max_cell[2]); iz++)
        {
          if (auto it = cells.find(Key(ix, iy, iz)); it != cells.end())
          {
            out.insert(out.end(), it->second.begin(), it->second.end());
          }
        }
      }
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
  }

  // Items in the cells at Chebyshev cell distance exactly k from the cell of p (ring 0 =
  // the cell itself); every point of such a cell is at least (k - 1) cells away from p.
  void Ring(const Point3D &p, long long int k, std::vector<std::size_t> &out) const
  {
    const auto c = Cell(p);
    for (long long int ix = c[0] - k; ix <= c[0] + k; ix++)
    {
      if (ix < min_cell[0] || ix > max_cell[0])
      {
        continue;
      }
      for (long long int iy = c[1] - k; iy <= c[1] + k; iy++)
      {
        if (iy < min_cell[1] || iy > max_cell[1])
        {
          continue;
        }
        const bool face_xy = std::abs(ix - c[0]) == k || std::abs(iy - c[1]) == k;
        for (long long int iz = c[2] - k; iz <= c[2] + k; iz++)
        {
          if (iz < min_cell[2] || iz > max_cell[2] ||
              (!face_xy && std::abs(iz - c[2]) != k))
          {
            continue;
          }
          if (auto it = cells.find(Key(ix, iy, iz)); it != cells.end())
          {
            out.insert(out.end(), it->second.begin(), it->second.end());
          }
        }
      }
    }
  }

  // Rings beyond this index hold no cell.
  long long int MaxRing(const Point3D &p) const
  {
    const auto c = Cell(p);
    long long int k = 0;
    for (int d = 0; d < 3; d++)
    {
      k = std::max({k, std::abs(c[d] - min_cell[d]), std::abs(max_cell[d] - c[d])});
    }
    return k;
  }

  double CellSize() const { return cell; }
  std::size_t Cells() const { return cells.size(); }

private:
  double cell;
  Point3D origin;
  std::array<long long int, 3> min_cell{std::numeric_limits<long long int>::max(),
                                        std::numeric_limits<long long int>::max(),
                                        std::numeric_limits<long long int>::max()};
  std::array<long long int, 3> max_cell{std::numeric_limits<long long int>::min(),
                                        std::numeric_limits<long long int>::min(),
                                        std::numeric_limits<long long int>::min()};
  std::unordered_map<std::uint64_t, std::vector<std::size_t>> cells;

  std::array<long long int, 3> Cell(const Point3D &p) const
  {
    std::array<long long int, 3> c{};
    for (int d = 0; d < 3; d++)
    {
      c[d] = static_cast<long long int>(std::floor((p[d] - origin[d]) / cell));
    }
    return c;
  }
  static std::uint64_t Key(long long int ix, long long int iy, long long int iz)
  {
    constexpr std::uint64_t offset = std::uint64_t{1} << 20;
    return ((static_cast<std::uint64_t>(ix) + offset) << 42) |
           ((static_cast<std::uint64_t>(iy) + offset) << 21) |
           (static_cast<std::uint64_t>(iz) + offset);
  }
};

void BoundingBox(const Point3D &a, const Point3D &b, Point3D &lo, Point3D &hi)
{
  for (int d = 0; d < 3; d++)
  {
    lo[d] = std::min(a[d], b[d]);
    hi[d] = std::max(a[d], b[d]);
  }
}

// Wall-clock stage timer and throttled progress lines through IdentificationInput::log.
class StageLog
{
public:
  explicit StageLog(const std::function<void(const std::string &)> &log_) : log(log_) {}

  void Begin(const std::string &stage)
  {
    name = stage;
    started = last = std::chrono::steady_clock::now();
  }
  // A progress line at most every 10 s: fraction of the stage's loop done (optionally of a
  // named inner loop, e.g. the candidate frames of one cluster signature).
  void Progress(std::size_t done, std::size_t total, const std::string &inner = "")
  {
    if (!log)
    {
      return;
    }
    const auto now = std::chrono::steady_clock::now();
    if (std::chrono::duration<double>(now - last).count() < 10.0)
    {
      return;
    }
    last = now;
    std::ostringstream text;
    text << "  Identification " << name << (inner.empty() ? "" : " (" + inner + ")") << ": "
         << done << " / " << total << " (" << std::fixed << std::setprecision(1)
         << (total ? 100.0 * static_cast<double>(done) / static_cast<double>(total) : 100.0)
         << " %), " << std::setprecision(1) << Elapsed() << " s\n";
    log(text.str());
  }
  void End(const std::string &counts)
  {
    if (!log)
    {
      return;
    }
    std::ostringstream text;
    text << "  Identification " << name << ": " << counts << " (" << std::fixed
         << std::setprecision(2) << Elapsed() << " s)\n";
    log(text.str());
  }
  double Elapsed() const
  {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - started)
        .count();
  }

private:
  const std::function<void(const std::string &)> &log;
  std::string name;
  std::chrono::steady_clock::time_point started, last;
};

struct VertexFeatureSite
{
  // A mesh vertex (corner / endpoint / junction) or the virtual corner of a fillet.
  std::optional<std::size_t> vertex;
  Point3D point{};
  std::string type;  // ConvexCorner | ConcaveCorner | Endpoint | Junction
  double turn_degrees = 0.0;
  double angle_degrees = 0.0;
  double corner_radius = 0.0;
  std::vector<double> arm_angles;
  std::vector<int> arm_conductors;
  // Unit directions away from the site along its arms (corner: two, endpoint: one,
  // junction: in the angular order of arm_angles) and the endpoint's gap direction: the
  // feature frame of the patch construction is built from them.
  std::vector<Point3D> arm_directions;
  Point3D gap_direction{};
  std::vector<std::pair<std::size_t, Interval>> window;  // run, interval claimed
  std::vector<std::size_t> runs_at_site;                 // runs incident to the site
  std::string boundary_law;
  std::vector<std::string> interfaces;
  int cluster = -1;
};

class Identifier
{
public:
  explicit Identifier(const IdentificationInput &input_)
    : original_input(input_), faces(input_.faces), quantizer(input_.radius),
      R(input_.radius)
  {
    // Working copy of the segments and vertices: the arc rule (design (b) 4 / 7, decision
    // 82(3)) merges the chains through the corner vertices it absorbs and demotes those
    // vertices to regular; the faces are read only.
    input.radius = input_.radius;
    input.segments = input_.segments;
    input.vertices = input_.vertices;
    input.log = input_.log;
    input.span_cap_allowances = input_.span_cap_allowances;
    ValidateSpanCapAllowances(input.span_cap_allowances);
  }
  // The same perimeter identified at another matching radius, silently and without a
  // census of its own (the knife-edge census's cluster-composition band).
  Identifier(const IdentificationInput &input_, double radius)
    : original_input(input_), faces(input_.faces), quantizer(radius), R(radius),
      census_enabled(false)
  {
    input.radius = radius;
    input.segments = input_.segments;
    input.vertices = input_.vertices;
  }

  IdentificationResult Identify();

private:
  const IdentificationInput &original_input;
  IdentificationInput input;
  bool census_enabled = true;
  const std::vector<IdentificationFace> &faces;
  Quantizer quantizer;
  double R;
  std::vector<Run> runs;
  std::vector<Chain> chains;
  std::map<int, std::size_t> chain_index;
  Point3D n_ref{};
  std::vector<std::optional<std::pair<std::string, std::string>>> segment_exclusion;
  // Metal plane of every run (distinct offsets along the reference process normal on the
  // decision grid): features never span two planes (decision 82(1)); metal of another
  // plane within 2R is the CrossLayer exclusion.
  std::vector<int> run_plane;
  // Parallel class of every rigid (single-run, non-excluded) run, -1 otherwise
  // (BuildDirectionClasses): the ONE parallel relation of the identification (USER decision
  // 184 (1)). Two rigid runs of one plane are parallel iff they share a class; parallel
  // rigid runs pair by the translational rule and by no other (the bent-pair and the event
  // rules skip them), so that no pair of rigid runs is handled by two rules or by none.
  std::vector<int> run_direction_class;
  std::size_t direction_class_count = 0;
  std::vector<std::vector<Claim>> claims;  // per run
  // Per run: (other chain, interval) of the pieces of a constant-separation pair along a
  // bend with that chain, interacting (a pair feature) or not; their cross-chord
  // interactions are not event cores.
  std::vector<std::vector<std::pair<int, Interval>>> bent_claims;
  // Decision 73(3) CrossLayer zones: per run, the intervals within 2R of metal off the
  // run's plane (facing layers, walls, staples); vertices within 2R of such metal.
  std::vector<std::vector<Interval>> cross_layer;
  std::set<std::size_t> cross_layer_vertices;
  std::vector<IdentifiedFeature> features;
  std::vector<double> feature_max_kappa;  // per feature, over its assigned portions
  // Number of sides a pair / parallel cluster must end up with after claim resolution
  // (0 for every other feature).
  std::vector<int> feature_sides;
  std::vector<VertexFeatureSite> sites;
  // Fitted arcs of the perimeter path (DetectArcs): joints (mesh vertices, path order),
  // the arc segments between the tangent points, the arms, the circle.
  struct Arc
  {
    std::vector<std::size_t> joints;
    std::vector<std::size_t> segments;
    std::size_t segment_before = 0, segment_after = 0;  // the arm segments
    Point3D tangent_a{}, tangent_b{};  // unit arm tangents (path direction)
    Point3D center{};
    Point3D origin{};  // virtual corner (turn < 180 deg) or arc midpoint
    double radius = 0.0;
    double turn = 0.0;         // radians, total
    double max_sagitta = 0.0;  // the largest chord sagitta on the fitted circle
    bool corner = false;       // radius < R (tangent arms): a rounded corner
    std::vector<std::size_t> absorbed_corners;  // input corner vertices demoted to regular
    int chain = -1;  // the arc's own chain (its segments are split off the arms' chains)
    int feature = -1;
  };
  // Zones of the chains connected through an arc (the arc plays the shared vertex): per
  // chain pair (ascending ids), the arcs joining them.
  std::map<std::pair<int, int>, std::vector<int>> through_arc;
  // Chains of the rounded-corner arcs: the arc plays the shared vertex of its arms, is
  // claimed whole by its corner feature and never pairs (a vertex has no pair).
  std::set<int> corner_arc_chains;
  // Through-vertex zones of chain A's run a with respect to chain B: within 2R of every
  // shared vertex and of every run of an arc joining the two chains. With on_arcs (the
  // cluster machinery, option A) the run and the arc's runs are evaluated on their fitted
  // circles; the pair rules keep the chord reading (candidate-gathering margins).
  std::vector<std::vector<Interval>>
  ThroughZones(std::size_t a, const Chain &A, const Chain &B, bool on_arcs = false) const
  {
    std::vector<std::vector<Interval>> zones;
    for (const std::size_t v : SharedVertices(A, B))
    {
      const Point3D &p = input.vertices[v].coordinate;
      zones.push_back(on_arcs
                          ? RunIntervalWithinPiece(a, CurvePiece{p, p, std::nullopt},
                                                   kThroughVertexZoneOverRadius * R)
                          : RunIntervalWithin(a, p, p, kThroughVertexZoneOverRadius * R));
    }
    const auto it =
        through_arc.find(std::make_pair(std::min(A.id, B.id), std::max(A.id, B.id)));
    if (it != through_arc.end())
    {
      for (const int arc : it->second)
      {
        const auto ci = chain_index.find(arcs[static_cast<std::size_t>(arc)].chain);
        if (ci == chain_index.end())
        {
          continue;
        }
        std::vector<Interval> zone;
        for (const std::size_t r : chains[ci->second].runs)
        {
          const auto within = on_arcs
                                  ? RunIntervalWithinPiece(a, WholeRunPiece(r),
                                                           kThroughVertexZoneOverRadius * R)
                                  : RunIntervalWithin(a, runs[r].start, runs[r].end,
                                                      kThroughVertexZoneOverRadius * R);
          zone.insert(zone.end(), within.begin(), within.end());
        }
        zones.push_back(MergeIntervals(std::move(zone), Tol()));
      }
    }
    return zones;
  }
  std::vector<Arc> arcs;
  std::vector<int> vertex_arc;   // per mesh vertex: the arc holding it as a joint, or -1
  std::vector<int> segment_arc;  // per segment: the BEND arc it lies on, or -1
  // Arc geometry of every run lying on a fitted arc (rounded corner or bend; option A): the
  // run's chord as the piece of the fitted circle between its end angles (theta0 at s = 0,
  // theta1 at s = length; the run parameter maps linearly onto the angle). Filled by
  // BuildRunArcGeometry after BuildRunIndex.
  struct RunArcGeometry
  {
    int arc = -1;
    ArcPiece piece;
  };
  std::vector<std::optional<RunArcGeometry>> run_arcs;
  // The largest sagitta of an arc run (0 without arcs): a run's arc lies within it of the
  // run's chord box, so a box query for the geometry within a margin adds it.
  double max_run_sagitta = 0.0;
  void BuildRunArcGeometry();
  double ArcSampleStep() const { return kArcSampleSpacingOverRadius * R; }
  // The point of run r at parameter s (on its arc when it lies on one, else on the chord).
  Point3D RunPoint(std::size_t r, double s) const
  {
    const auto &geometry = run_arcs[r];
    if (geometry)
    {
      const double t = runs[r].length > 0.0 ? s / runs[r].length : 0.0;
      return geometry->piece.AtFraction(t);
    }
    return runs[r].At(s);
  }
  // The geometry of the run interval [s0, s1] (an arc piece or a segment).
  CurvePiece RunPiece(std::size_t r, double s0, double s1) const
  {
    CurvePiece piece;
    piece.a = runs[r].At(s0);
    piece.b = runs[r].At(s1);
    const auto &geometry = run_arcs[r];
    if (geometry)
    {
      const double length = runs[r].length;
      piece.arc = geometry->piece.SubArc(length > 0.0 ? s0 / length : 0.0,
                                         length > 0.0 ? s1 / length : 1.0);
      piece.a = piece.arc->At(piece.arc->theta0);
      piece.b = piece.arc->At(piece.arc->theta1);
    }
    return piece;
  }
  CurvePiece WholeRunPiece(std::size_t r) const { return RunPiece(r, 0.0, runs[r].length); }
  // Unit tangent of run r at parameter s in the run direction (on its arc when on one).
  Point3D RunTangent(std::size_t r, double s) const
  {
    const auto &geometry = run_arcs[r];
    if (geometry)
    {
      const double t = runs[r].length > 0.0 ? s / runs[r].length : 0.0;
      const ArcPiece &arc = geometry->piece;
      return arc.Tangent(arc.theta0 + t * (arc.theta1 - arc.theta0));
    }
    return runs[r].tangent;
  }
  // Turn at the joint before run k of the chain (Chain::joint_turn[k]) read on the fitted
  // geometry: between the end tangent of the previous run and the start tangent of run k on
  // their arcs where they lie on one (two chords of one circle, or a chord and its tangent
  // arm, turn by nothing); the chord reading where neither run is on an arc.
  double GeometricJointTurn(const Chain &C, std::size_t k) const
  {
    const std::size_t m = C.runs.size();
    if (k >= m || C.joint_turn.size() != m)
    {
      return 0.0;
    }
    const std::size_t prev = (k + m - 1) % m;
    if (!run_arcs[C.runs[k]] && !run_arcs[C.runs[prev]])
    {
      return C.joint_turn[k];
    }
    if ((k == 0 && !C.closed) || C.joint_turn[k] <= 0.0)
    {
      return 0.0;
    }
    const Point3D t_in = RunTangent(C.runs[prev], runs[C.runs[prev]].length);
    const Point3D t_out = RunTangent(C.runs[k], 0.0);
    return std::acos(std::clamp(Dot(t_in, t_out), -1.0, 1.0));
  }
  // Minimal distance between two run pieces (exact for two chords; the arc side sampled).
  double PieceDistance(const CurvePiece &A, const CurvePiece &B) const
  {
    return PiecePieceDistance(A, B, ArcSampleStep());
  }
  // Interval of the run parameter where the run's point (on its arc or chord) lies within
  // `distance` of the piece: the convex solver on two chords (the former
  // RunIntervalWithin(run, a, b, distance) exactly), the sampled solver otherwise,
  // bracketed by the convex sublevel set of the chord distance below distance + sagitta.
  std::vector<Interval> RunIntervalWithinPiece(std::size_t run, const CurvePiece &piece,
                                               double distance) const;
  // The candidate runs whose boxes lie within margin of a piece's box (the sagitta of an
  // arc piece is inside the box of its chord ends only when the box is enlarged by it).
  std::vector<std::size_t> RunsNearPiece(const CurvePiece &piece, double margin) const
  {
    Point3D lo, hi;
    BoundingBox(piece.a, piece.b, lo, hi);
    return RunsNear(lo, hi, margin + piece.Sagitta() + max_run_sagitta);
  }
  static void PieceBoundingBox(const CurvePiece &piece, Point3D &lo, Point3D &hi)
  {
    BoundingBox(piece.a, piece.b, lo, hi);
    const double sagitta = piece.Sagitta();
    for (int d = 0; d < 3; d++)
    {
      lo[d] -= sagitta;
      hi[d] += sagitta;
    }
  }
  // The bend arc a run lies on (its segments all lie on one arc or on none), or -1.
  int RunArc(std::size_t run) const
  {
    const auto &segments = runs[run].segments;
    if (segments.empty() || segment_arc.empty())
    {
      return -1;
    }
    const int arc = segment_arc[segments.front().segment];
    for (const auto &rs : segments)
    {
      if (segment_arc[rs.segment] != arc)
      {
        return -1;
      }
    }
    return arc;
  }
  // Chains joined through the corner vertices absorbed by a bend are one chain (union-find
  // over the input chain ids, applied to the segments before the runs are built), so that a
  // chord count that puts a corner vertex inside a bend does not split the chain.
  std::map<int, int> chain_group;
  int ChainGroup(int chain) const
  {
    auto it = chain_group.find(chain);
    while (it != chain_group.end() && it->second != chain)
    {
      chain = it->second;
      it = chain_group.find(chain);
    }
    return chain;
  }
  std::vector<std::vector<EventCore>> cluster_cores;
  std::vector<std::vector<std::size_t>> cluster_sites;
  // Claimed portions per cluster (run, interval), sorted; extended by ExtendClusters.
  std::vector<std::vector<std::pair<std::size_t, Interval>>> cluster_claimed;
  // Cluster extension diagnostics (decision 85(2)).
  std::size_t extension_passes = 0, extension_portions = 0, extension_sites = 0;
  double extension_length = 0.0;
  double stack_end_third_body_length = 0.0;
  // Translational (pair / stack) stretches absorbed by the cluster bounding them on both
  // sides (decision 224): count, length, longest.
  std::size_t extension_translational_pieces = 0, extension_translational_two_sided = 0;
  double extension_translational_length = 0.0, extension_translational_max_length = 0.0,
         extension_translational_two_sided_length = 0.0;
  std::map<std::size_t, int> vertex_feature;  // mesh vertex -> feature id
  // Raw material of the pair / stack assembly (decision 82(2), BuildPairsAndStacks): the
  // locally constant, interacting facing relations of the bent-pair rule (one link per
  // chain pair and separation group; chain_a == chain_b for a chain facing itself across a
  // fold) and the translational spans of the rigid parallel runs. Links and spans sharing
  // a run over a common interval are one cross-section: a multi-edge stack.
  struct LinkPiece
  {
    std::size_t run;
    Interval interval;
    bool curved;
    double max_kappa;
    double separation;  // the piece's mean separation (exact or chord, see `exact`)
    bool exact;
  };
  // Separation of a link per curvature class (index 0 straight-like, 1 curved), decision
  // 85(1): the length-weighted mean of the class's exact piece separations (two exactly
  // parallel straight runs, two concentric fitted arcs) when it has any, else of its chord
  // readings (`exact` false: a polyline bend with chords beyond 2R, a taper, a transition).
  struct ClassSeparation
  {
    double value = 0.0;
    bool exact = false;
    bool present = false;
  };
  using LinkSeparation = std::array<ClassSeparation, 2>;
  static LinkSeparation LinkSeparations(const std::array<std::vector<LinkPiece>, 2> &pieces)
  {
    LinkSeparation result;
    for (int curved_class = 0; curved_class < 2; curved_class++)
    {
      double exact_sum = 0.0, exact_length = 0.0, chord_sum = 0.0, chord_length = 0.0;
      for (const auto &side : pieces)
      {
        for (const auto &piece : side)
        {
          if (piece.curved != (curved_class == 1))
          {
            continue;
          }
          const double length = piece.interval.second - piece.interval.first;
          if (piece.exact)
          {
            exact_sum += piece.separation * length;
            exact_length += length;
          }
          else
          {
            chord_sum += piece.separation * length;
            chord_length += length;
          }
        }
      }
      auto &entry = result[static_cast<std::size_t>(curved_class)];
      if (exact_length > 0.0)
      {
        entry = {exact_sum / exact_length, true, true};
      }
      else if (chord_length > 0.0)
      {
        entry = {chord_sum / chord_length, false, true};
      }
    }
    for (int curved_class = 0; curved_class < 2; curved_class++)
    {
      auto &entry = result[static_cast<std::size_t>(curved_class)];
      const auto &other = result[static_cast<std::size_t>(1 - curved_class)];
      if (!entry.present && other.present)
      {
        entry = {other.value, other.exact, false};
      }
    }
    return result;
  }
  struct PairLink
  {
    int chain_a = -1, chain_b = -1;
    LinkSeparation separation;  // per curvature class
    Point3D lead_point{}, lateral_ab{};
    std::size_t run_a = 0, run_b = 0;              // the runs at the lead pair of points
    std::array<std::vector<LinkPiece>, 2> pieces;  // side 0 on chain a, side 1 on chain b
    bool Self() const { return chain_a == chain_b; }
  };
  std::vector<PairLink> pair_links;
  struct TranslationalSpan
  {
    struct Member
    {
      std::size_t run;
      double w;
      int gap_sign;
      Interval interval;  // run parameter
    };
    std::vector<Member> members;  // sorted by lateral offset w
    double lo = 0.0, hi = 0.0;
    Point3D axis{}, lateral{};
  };
  std::vector<TranslationalSpan> translational_spans;
  // Per chain id: the chain intervals (sorted, disjoint) claimed by clusters and by vertex
  // windows before the pair / stack assembly; a member of a cross-section taken there is
  // not part of it (the stack-end rule).
  std::map<int, std::vector<Interval>> taken_on_chain;
  // Claims of the same priority overlapping on a run would be resolved by feature id — a
  // tie-break the rules must never need (decision 82(2)); counted and reported.
  std::size_t same_priority_claim_overlaps = 0;
  // Stack assembly diagnostics (reported under Diagnostics): elementary intervals whose
  // offsets took a geometric lateral distance for want of a consecutive link, and those
  // whose composition reached kStackCompositionCap members.
  std::size_t stack_geometric_offsets = 0;
  std::size_t stack_composition_cap_hits = 0;
  std::size_t stack_images_merged = 0;
  // Closure of the extension loop (decision 93): stopped by the pass cap or by a repeated
  // pass (reported).
  bool extension_cap_reached = false, extension_repeat_detected = false;
  // Acceleration only (every result is the former all-pairs answer): the runs by their
  // boxes and the runs incident to every vertex.
  std::optional<UniformGrid> run_grid;
  std::vector<std::vector<std::size_t>> runs_at_vertex;
  std::vector<long long int> run_of_segment;  // -1 for segments outside every run
  // Shared vertices of two chains (sorted), computed once per chain pair from the smaller
  // chain's vertices and the runs incident to them (the former set intersection of the two
  // vertex sets, which cost the size of both chains for every run pair).
  mutable std::map<std::pair<int, int>, std::vector<std::size_t>> shared_vertices_cache;
  StageLog stage{input.log};

  double Tol() const { return quantizer.Quantum(); }

  const std::vector<std::size_t> &SharedVertices(const Chain &A, const Chain &B) const
  {
    const auto key = std::make_pair(std::min(A.id, B.id), std::max(A.id, B.id));
    auto it = shared_vertices_cache.find(key);
    if (it != shared_vertices_cache.end())
    {
      return it->second;
    }
    const Chain &small = A.vertices.size() <= B.vertices.size() ? A : B;
    const Chain &large = &small == &A ? B : A;
    std::vector<std::size_t> shared;
    for (const std::size_t v : small.vertices)  // std::set: ascending
    {
      // A chain's vertex set is the set of endpoints of its segments.
      const auto &segments = input.vertices[v].segments;
      if (std::any_of(segments.begin(), segments.end(),
                      [&](std::size_t s)
                      {
                        return run_of_segment[s] >= 0 &&
                               runs[static_cast<std::size_t>(run_of_segment[s])].chain ==
                                   large.id;
                      }))
      {
        shared.push_back(v);
      }
    }
    return shared_vertices_cache.emplace(key, std::move(shared)).first->second;
  }

  // Candidate runs (sorted) whose boxes lie within margin of [lo, hi].
  std::vector<std::size_t> RunsNear(const Point3D &lo, const Point3D &hi,
                                    double margin) const
  {
    return run_grid->Query(lo, hi, margin);
  }
  std::vector<std::size_t> RunsNearRun(std::size_t r, double margin) const
  {
    Point3D lo, hi;
    BoundingBox(runs[r].start, runs[r].end, lo, hi);
    return RunsNear(lo, hi, margin);
  }
  void BuildRunIndex()
  {
    Point3D lo{}, hi{};
    bool first = true;
    for (const auto &run : runs)
    {
      Point3D rlo, rhi;
      BoundingBox(run.start, run.end, rlo, rhi);
      for (int d = 0; d < 3; d++)
      {
        lo[d] = first ? rlo[d] : std::min(lo[d], rlo[d]);
        hi[d] = first ? rhi[d] : std::max(hi[d], rhi[d]);
      }
      first = false;
    }
    run_grid.emplace(4.0 * R, lo);
    for (std::size_t r = 0; r < runs.size(); r++)
    {
      Point3D rlo, rhi;
      BoundingBox(runs[r].start, runs[r].end, rlo, rhi);
      run_grid->Insert(r, rlo, rhi);
    }
    run_of_segment.assign(input.segments.size(), -1);
    for (std::size_t r = 0; r < runs.size(); r++)
    {
      for (const auto &rs : runs[r].segments)
      {
        run_of_segment[rs.segment] = static_cast<long long int>(r);
      }
    }
    runs_at_vertex.assign(input.vertices.size(), {});
    for (std::size_t r = 0; r < runs.size(); r++)
    {
      runs_at_vertex[runs[r].start_vertex].push_back(r);
      if (runs[r].end_vertex != runs[r].start_vertex)
      {
        runs_at_vertex[runs[r].end_vertex].push_back(r);
      }
    }
    for (auto &incident : runs_at_vertex)
    {
      std::sort(incident.begin(), incident.end());
    }
  }

  int NewFeature(const std::string &type, nlohmann::json signature, int chirality = 1)
  {
    IdentifiedFeature feature;
    feature.id = static_cast<int>(features.size());
    feature.type = type;
    signature["Type"] = type;
    feature.signature = std::move(signature);
    std::tie(feature.signature_key, feature.hash) =
        SignatureKeyAndHash(feature.signature, type);
    feature.chirality = chirality;
    features.push_back(std::move(feature));
    feature_max_kappa.push_back(0.0);
    feature_sides.push_back(0);
    return features.back().id;
  }

  void ExcludeRun(std::size_t run, const std::string &cls, const std::string &reason)
  {
    runs[run].excluded = true;
    for (const auto &rs : runs[run].segments)
    {
      segment_exclusion[rs.segment] = std::make_pair(cls, reason);
    }
  }

  bool IsFeatureVertex(std::size_t v) const
  {
    const auto &vertex = input.vertices[v];
    if (!vertex.physical_type || *vertex.physical_type == MetalEdgeVertexType::REGULAR)
    {
      return false;
    }
    if (*vertex.physical_type == MetalEdgeVertexType::ENDPOINT &&
        (vertex.on_truncation_boundary || vertex.on_port_boundary))
    {
      return false;  // simulation cut / port cut
    }
    return cross_layer_vertices.find(v) == cross_layer_vertices.end();
  }

  // Metal of different edge-connected components meets at this vertex (a point contact).
  bool IsPointContact(std::size_t v) const
  {
    std::set<int> conductors;
    for (const std::size_t s : input.vertices[v].segments)
    {
      const auto &segment = input.segments[s];
      if (!segment.truncation && !segment.exclusion)
      {
        conductors.insert(segment.conductor);
      }
    }
    return conductors.size() > 1;
  }

  int SitePlane(const VertexFeatureSite &site) const
  {
    MFEM_VERIFY(!site.runs_at_site.empty(), "A vertex site without runs!");
    return run_plane[site.runs_at_site.front()];
  }

  // Signed process normal (substrate -> vacuum) of a spatial feature from its own runs:
  // length weighted like ReferenceProcessNormal but WITHOUT SignCanonical, so that a
  // feature of a flipped plane (a flip-chip top chip, normal -z) is framed with w pointing
  // into ITS vacuum and the coupon's substrate half-space (w < 0) lands in its substrate
  // (decision 266). The device-global n_ref is a planarity reference only (|dot| tests,
  // plane numbering): it is sign-canonical and puts every flipped plane's coupon upside
  // down. A feature whose runs disagree in sign has no process side: fail closed naming it.
  // The returned direction is +-n_ref (every framed run is parallel to n_ref within
  // kParallelCosineTolerance by the NonPlanar rule), oriented by the feature's own runs:
  // single-plane devices keep their frames bit-identical.
  Point3D FeatureProcessNormal(const std::vector<std::pair<std::size_t, double>> &weighted,
                               const std::string &what) const
  {
    MFEM_VERIFY(!weighted.empty(), "Process normal of " << what << " without runs!");
    Point3D sum{};
    const Point3D first = runs[weighted.front().first].process_normal;
    for (const auto &[r, weight] : weighted)
    {
      const Point3D &n = runs[r].process_normal;
      MFEM_VERIFY(Dot(n, first) > 0.0,
                  "Runs of " << what << " disagree on the sign of their process normal (["
                             << first[0] << ", " << first[1] << ", " << first[2] << "] vs ["
                             << n[0] << ", " << n[1] << ", " << n[2] << "] on run " << r
                             << "): no common substrate -> vacuum side!");
      sum = Add(sum, Scale(std::max(weight, 0.0), n));
    }
    if (Norm(sum) <= 0.0)
    {
      sum = first;  // zero-length weights: the (agreeing) direction itself
    }
    return SignedReferenceNormal(Normalize(sum), what);
  }

  // The reference normal with the sign of the given (parallel) process normal.
  Point3D SignedReferenceNormal(const Point3D &n, const std::string &what) const
  {
    const double cosine = Dot(n, n_ref);
    MFEM_VERIFY(!DirectionLess(std::abs(cosine), 1.0 - kParallelCosineTolerance),
                "Process normal of "
                    << what << " ([" << n[0] << ", " << n[1] << ", " << n[2]
                    << "]) is not parallel to the reference process normal!");
    return cosine < 0.0 ? Scale(-1.0, n_ref) : n_ref;
  }

  Point3D SiteProcessNormal(const VertexFeatureSite &site) const
  {
    std::vector<std::pair<std::size_t, double>> weighted;
    for (const std::size_t r : site.runs_at_site)
    {
      weighted.emplace_back(r, runs[r].length);
    }
    std::ostringstream what;
    what << site.type << " site at (" << site.point[0] << ", " << site.point[1] << ", "
         << site.point[2] << ")";
    return FeatureProcessNormal(weighted, what.str());
  }

  std::vector<std::size_t> RunsAtVertex(std::size_t v) const
  {
    std::vector<std::size_t> result;
    for (const std::size_t r : runs_at_vertex[v])
    {
      if (!runs[r].excluded)
      {
        result.push_back(r);
      }
    }
    return result;
  }

  // Direction away from the vertex along the run.
  Point3D ArmDirection(std::size_t run, std::size_t v) const
  {
    return runs[run].start_vertex == v ? runs[run].tangent : Scale(-1.0, runs[run].tangent);
  }

  void DetectArcs();
  void ClassifyPlanes();
  void BuildDirectionClasses();
  bool ParallelRigidRuns(std::size_t a, std::size_t b) const
  {
    return run_direction_class[a] >= 0 && run_direction_class[a] == run_direction_class[b];
  }
  void ClassifyVertices();
  void BuildArcSites();
  void ComputeCurvature();
  void BuildTranslationalFeatures();
  void BuildBentPairs();
  void BuildPairsAndStacks();
  void AssembleStack(const std::vector<std::size_t> &link_items,
                     const std::vector<std::size_t> &span_items);
  void BuildClusters();
  double ExtendClusters(bool measure_joint_claims = false);
  void EmitClusters();
  void RecordContextFeatureIds(int feature);
  std::vector<SignaturePortion> ClusterSignaturePortions(
      const std::vector<std::pair<std::size_t, Interval>> &claimed) const;
  // Spatial-support contract v3 (decision 282): the support of cluster c in the candidate
  // frame (origin, x, y) — the claims-derived box, the face rules T1-T3 with growth, the
  // device-plan context — with its manifest record (units of R).
  struct ClusterSupportResult
  {
    std::optional<FrameSupport> support;  // nullopt: unboxable in this frame
    nlohmann::json record;
    // Every context piece is implied by the decision-236 contract (a straight chain piece
    // collinear with the claim it continues from a claim-cut end to a box face): the
    // device plan clipped to the box IS the legacy coupon geometry.
    bool legacy_equivalent = false;
    bool grown = false;
    bool exceeds_span_cap = false;
    double chain_length = 0.0;  // mesh units
    double foreign_length = 0.0;
    double fictitious_continuation_length = 0.0;
    std::size_t band_hits = 0;
  };
  // `span_cap_over_R`: the plan span cap of this cluster (kSupportSpanCapOverRadius, or the
  // resolved per-case allowance of block (b) DESIGN section 3 (a) / A4).
  ClusterSupportResult ClusterSupport(std::size_t c,
                                      const std::vector<SignaturePortion> &portions,
                                      const std::vector<SignatureVertex> &vertices,
                                      const Point3D &origin, const Point3D &x,
                                      const Point3D &y, double span_cap_over_R) const;
  // Sites (vertex features) incident to every run, for the context's vertex census.
  std::multimap<std::size_t, std::size_t> sites_by_run;
  // Every cluster's claimed intervals per run (interval, cluster): the context pieces of a
  // cluster's box lying on ANOTHER cluster's claims are split at the claim cuts and never
  // chain (decision 285 (1): owned by the other coupon).
  std::map<std::size_t, std::vector<std::pair<Interval, std::size_t>>> claims_by_run;
  // Feature index of every cluster (EmitClusters) and of every free site (Assign), for the
  // non-hashed ownership records of the context census.
  std::vector<int> cluster_feature, site_feature;
  // Spatial-support summary over the clusters (IdentificationResult::spatial_support).
  IdentificationResult::SpatialSupportSummary spatial_support_summary;
  // The span-cap allowances (indices into input.span_cap_allowances) some cluster resolved.
  std::set<std::size_t> used_span_cap_allowances;
  void BuildVertexWindows();
  void Assign(IdentificationResult &result);
  nlohmann::json KnifeEdgeCensus() const;
  nlohmann::json ClusterCompositionBand(const IdentificationResult &result) const;
  std::vector<Interval> RunIntervalWithin(std::size_t run, const Point3D &a,
                                          const Point3D &b, double distance) const;
  std::vector<Interval>
  ChainWindow(std::size_t run, std::size_t from_vertex, double length,
              std::vector<std::pair<std::size_t, Interval>> &out) const;

  // Curved-edge chain rule helpers on the chain arc length x.
  double WindowedCurvature(const Chain &chain, double x) const;
  // The signed windowed curvature at x (positive where the chain turns toward its metal).
  double SignedWindowedCurvature(const Chain &chain, double x) const;
  double MaxCurvature(const Chain &chain, double x0, double x1) const;
  // Integral of the signed windowed curvature over the chain interval [x0, x1]: the signed
  // turn toward the metal (radians; exact on the piecewise-linear nodes).
  double SignedTurn(const Chain &chain, double x0, double x1) const;
  // Extremes of the signed windowed curvature over the chain interval [x0, x1]: the
  // largest value toward the metal (positive) and toward the gap (negative), accumulated
  // into `extremes` (so that several intervals of one feature side combine).
  struct SignedCurvatureExtremes
  {
    double toward_metal = 0.0, toward_gap = 0.0;
  };
  void AccumulateSignedCurvature(const Chain &chain, double x0, double x1,
                                 SignedCurvatureExtremes &extremes) const;
  // Convexity read on the extremes: "Convex" when the signed windowed curvature bends
  // around the metal (the edge of a disk), "Concave" when around the gap (the edge of a
  // hole), "Mixed" when both senses reach the curved regime (windowed bend radius below
  // kStraightBendRadiusOverRadius R: an S-bend inside one curved section — reported in
  // the signature, never read as one convexity; a straight-like wobble of the minority
  // sense is not a bend of the class), empty when the interval carries no curvature.
  std::string ConvexityName(const SignedCurvatureExtremes &extremes) const;
  bool IsCurvedAt(const Chain &chain, double x, std::size_t *section = nullptr) const;
  struct ChainPoint
  {
    std::size_t run = 0;
    double s = 0.0;
    double x = 0.0;
    double distance = 0.0;
  };
  // Closest point of the chain to p; with exclude_x, the runs within the self-pair
  // neighbourhood (kSelfPairNeighbourhoodOverRadius R of arc length, wrapping on a closed
  // chain) of the chain position exclude_x are not candidates (a chain facing itself); with
  // max_distance, the search stops beyond it and returns an infinite distance when no
  // candidate lies within it (the callers use the result only within their reach).
  ChainPoint ClosestPointOnChain(const Chain &chain, const Point3D &p,
                                 std::optional<double> exclude_x = std::nullopt,
                                 std::optional<double> max_distance = std::nullopt) const;
  // Image of the chain position x of `source` on `other` read on the fitted arcs (option
  // A): the point on the source's arc where it lies on one, and the foot on the other
  // chain's arc run holding the radial projection of that point (two concentric bends image
  // radially, whatever their chords); the chord reading of ClosestPointOnChain otherwise.
  ChainPoint ArcAwareImage(const Chain &source, double x, const Chain &other,
                           std::optional<double> exclude_x,
                           std::optional<double> max_distance) const;
  // Arc-length distance between the chain interval [x0, x1] and the chain position x
  // (wrapping on a closed chain).
  double ChainArcDistance(const Chain &chain, double x0, double x1, double x) const;
  Point3D ChainAt(const Chain &chain, double x, std::size_t *run_out = nullptr,
                  double *s_out = nullptr) const;
  std::size_t RunIndexInChain(const Chain &chain, std::size_t run) const;
};

// Runs whose process normal is not parallel to the reference normal are walls / staples
// (NonPlanar). Metal off a run's own plane within 2R of it (a facing layer across a gap, a
// wall or staple standing nearby) invalidates the planar coupon there: those parts of the
// run, solved analytically on the distance to every such face, are the CrossLayer zones
// (decision 73(3)); feature vertices within 2R of such metal are excluded the same way.
// Every other plane carrying metal is identified in its own right.
void Identifier::ClassifyPlanes()
{
  n_ref = ReferenceProcessNormal(runs);
  cross_layer.assign(runs.size(), {});
  for (std::size_t r = 0; r < runs.size(); r++)
  {
    if (DirectionLess(std::abs(Dot(runs[r].process_normal, n_ref)),
                      1.0 - kParallelCosineTolerance))
    {
      ExcludeRun(r, "NonPlanar",
                 "metal edge whose process normal is not parallel to the reference process "
                 "normal (wall, staple, via)");
    }
  }
  // Plane of every run: the distinct offsets of the run midpoints along the reference
  // normal, merged on the decision grid, numbered in increasing offset.
  {
    std::vector<double> offsets;
    for (const auto &run : runs)
    {
      offsets.push_back(Dot(Scale(0.5, Add(run.start, run.end)), n_ref));
    }
    std::vector<double> planes;
    for (const double offset : std::set<double>(offsets.begin(), offsets.end()))
    {
      if (planes.empty() || !quantizer.Equal(offset, planes.back()))
      {
        planes.push_back(offset);
      }
    }
    run_plane.assign(runs.size(), 0);
    for (std::size_t r = 0; r < runs.size(); r++)
    {
      auto it = std::lower_bound(planes.begin(), planes.end(), offsets[r] - Tol());
      MFEM_VERIFY(it != planes.end() && quantizer.Equal(*it, offsets[r]),
                  "Run plane not found!");
      run_plane[r] = static_cast<int>(it - planes.begin());
    }
  }
  if (faces.empty())
  {
    return;
  }

  struct Face
  {
    std::size_t index;
    bool planar;
    double offset;  // along n_ref (planar faces)
    Point3D lower, upper;
  };
  std::vector<Face> face_boxes;
  face_boxes.reserve(faces.size());
  for (std::size_t f = 0; f < faces.size(); f++)
  {
    const auto &face = faces[f];
    if (face.vertices.size() < 3)
    {
      continue;
    }
    Face entry{
        f,
        !DirectionLess(std::abs(Dot(face.normal, n_ref)), 1.0 - kParallelCosineTolerance),
        0.0, face.vertices.front(), face.vertices.front()};
    for (const auto &p : face.vertices)
    {
      entry.offset += Dot(p, n_ref) / static_cast<double>(face.vertices.size());
      for (int d = 0; d < 3; d++)
      {
        entry.lower[d] = std::min(entry.lower[d], p[d]);
        entry.upper[d] = std::max(entry.upper[d], p[d]);
      }
    }
    face_boxes.push_back(entry);
  }
  const double reach = 2.0 * R;
  // Faces by their boxes: the candidates of a run / vertex are the faces whose boxes lie
  // within the reach of it, visited in face order (the former loop over every face).
  Point3D faces_lo{}, faces_hi{};
  for (std::size_t k = 0; k < face_boxes.size(); k++)
  {
    for (int d = 0; d < 3; d++)
    {
      faces_lo[d] =
          k == 0 ? face_boxes[k].lower[d] : std::min(faces_lo[d], face_boxes[k].lower[d]);
      faces_hi[d] =
          k == 0 ? face_boxes[k].upper[d] : std::max(faces_hi[d], face_boxes[k].upper[d]);
    }
  }
  UniformGrid face_grid(8.0 * R, faces_lo);
  for (std::size_t k = 0; k < face_boxes.size(); k++)
  {
    face_grid.Insert(k, face_boxes[k].lower, face_boxes[k].upper);
  }
  const double margin = reach + 2.0 * Tol();
  auto OffPlane = [&](const Face &face, double offset)
  {
    // Metal in the run's own plane is the identified perimeter itself; a parallel layer
    // is a candidate only when it lies within reach along the normal.
    if (!face.planar)
    {
      return true;
    }
    return !quantizer.Equal(face.offset, offset) &&
           quantizer.Less(std::abs(face.offset - offset), reach);
  };
  auto NearBox = [&](const Face &face, const Point3D &lower, const Point3D &upper)
  {
    for (int d = 0; d < 3; d++)
    {
      if (face.lower[d] > upper[d] + reach || face.upper[d] < lower[d] - reach)
      {
        return false;
      }
    }
    return true;
  };
  for (std::size_t r = 0; r < runs.size(); r++)
  {
    if (runs[r].excluded)
    {
      continue;
    }
    const Run &run = runs[r];
    const double offset = Dot(Scale(0.5, Add(run.start, run.end)), n_ref);
    Point3D lower = run.start, upper = run.start;
    for (int d = 0; d < 3; d++)
    {
      lower[d] = std::min(run.start[d], run.end[d]);
      upper[d] = std::max(run.start[d], run.end[d]);
    }
    stage.Progress(r, runs.size());
    std::vector<Interval> zones;
    for (const std::size_t k : face_grid.Query(lower, upper, margin))
    {
      const Face &face = face_boxes[k];
      if (!OffPlane(face, offset) || !NearBox(face, lower, upper))
      {
        continue;
      }
      const auto &polygon = faces[face.index];
      const auto zone = ConvexSublevelInterval(
          [&](double t)
          { return PointPolygonDistance(run.At(t), polygon.vertices, polygon.normal); },
          0.0, run.length, reach, quantizer);
      if (zone)
      {
        zones.push_back(*zone);
      }
    }
    cross_layer[r] = MergeIntervals(std::move(zones), Tol());
  }
  for (std::size_t v = 0; v < input.vertices.size(); v++)
  {
    const auto &vertex = input.vertices[v];
    if (!vertex.physical_type || *vertex.physical_type == MetalEdgeVertexType::REGULAR)
    {
      continue;
    }
    std::optional<double> offset;
    for (const std::size_t s : vertex.segments)
    {
      const auto &segment = input.segments[s];
      if (!segment.truncation && !segment.exclusion && segment.chain >= 0)
      {
        offset = Dot(vertex.coordinate, n_ref);
        break;
      }
    }
    if (!offset)
    {
      continue;
    }
    for (const std::size_t k :
         face_grid.Query(vertex.coordinate, vertex.coordinate, margin))
    {
      const Face &face = face_boxes[k];
      if (!OffPlane(face, *offset) || !NearBox(face, vertex.coordinate, vertex.coordinate))
      {
        continue;
      }
      const auto &polygon = faces[face.index];
      if (quantizer.Less(
              PointPolygonDistance(vertex.coordinate, polygon.vertices, polygon.normal),
              reach))
      {
        cross_layer_vertices.insert(v);
        break;
      }
    }
  }
}

// Corner / endpoint / junction sites at the feature vertices of the non-excluded runs.
// Parallel classes of the rigid runs (USER decision 184 (1), 2026-10-01). Two straight runs
// are parallel when their tangents agree within the parallel cosine tolerance
// (|cos| > 1 - kParallelCosineTolerance: 1.4e-4 rad, the tolerance every other parallel
// test of the identification uses; dimensionless, so independent of the length unit and
// of R). The classes are the connected components of that relation per metal plane
// (single linkage on the run tangents: a mesh-noise tilt of a straight edge (DS-SCT-002:
// 1.8e-4 um over 186 um, 9e-7 rad) leaves it in the class of its exactly axis-aligned
// partner). Before this rule the translational stage grouped the runs on the 1e-9
// DirectionKey grid while the bent-pair and event stages skipped rigid pairs within the
// cosine tolerance: runs tilted by between 1e-9 and 1.4e-4 rad were paired by NEITHER rule
// (the E8-1 flux lines: a 3-edge stack + an isolated fourth edge). The components are
// built from the DirectionKey buckets (a canonical, order-independent partition) sorted by
// their in-plane angle: consecutive buckets within the tolerance are one class (the wrap
// at 0 / pi included), so the class of a run does not depend on the input numbering.
void Identifier::BuildDirectionClasses()
{
  run_direction_class.assign(runs.size(), -1);
  // In-plane orthonormal basis about the reference process normal.
  Point3D e1{};
  {
    int least = 0;
    for (int d = 1; d < 3; d++)
    {
      least = std::abs(n_ref[d]) < std::abs(n_ref[least]) ? d : least;
    }
    e1[least] = 1.0;
    e1 = Normalize(Sub(e1, Scale(Dot(e1, n_ref), n_ref)));
  }
  const Point3D e2 = Normalize(Cross(n_ref, e1));
  struct Bucket
  {
    int plane;
    std::array<long long int, 3> key;
    Point3D sum{};
    std::vector<std::size_t> members;
    Point3D representative{};
    double angle = 0.0;  // in-plane direction angle in [0, pi)
  };
  std::map<std::pair<int, std::array<long long int, 3>>, std::size_t> bucket_index;
  std::vector<Bucket> buckets;
  for (std::size_t r = 0; r < runs.size(); r++)
  {
    if (runs[r].excluded || !chains[chain_index.at(runs[r].chain)].Rigid())
    {
      continue;  // chains with joints pair through the curved-edge chain rule
    }
    const Point3D t = SignCanonical(runs[r].tangent);
    const auto key = std::make_pair(run_plane[r], DirectionKey(t, 1.0e-9));
    auto it = bucket_index.find(key);
    if (it == bucket_index.end())
    {
      it = bucket_index.emplace(key, buckets.size()).first;
      buckets.push_back({run_plane[r], key.second, {}, {}, {}, 0.0});
    }
    Bucket &bucket = buckets[it->second];
    bucket.sum = Add(bucket.sum, Scale(runs[r].length, t));
    bucket.members.push_back(r);
  }
  const double pi = std::acos(-1.0);
  for (auto &bucket : buckets)
  {
    bucket.representative = Normalize(bucket.sum);
    double angle =
        std::atan2(Dot(bucket.representative, e2), Dot(bucket.representative, e1));
    if (angle < 0.0)
    {
      angle += pi;
    }
    bucket.angle = angle >= pi ? angle - pi : angle;
  }
  auto Parallel = [](const Point3D &u, const Point3D &v)
  { return !DirectionLess(std::abs(Dot(u, v)), 1.0 - kParallelCosineTolerance); };
  // Per plane: the buckets sorted by angle, consecutive ones within the tolerance linked.
  std::vector<std::size_t> order(buckets.size());
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(),
            [&](std::size_t a, std::size_t b)
            {
              return std::make_tuple(buckets[a].plane, buckets[a].angle, buckets[a].key) <
                     std::make_tuple(buckets[b].plane, buckets[b].angle, buckets[b].key);
            });
  UnionFind uf(buckets.size());
  for (std::size_t i = 0; i < order.size();)
  {
    std::size_t j = i;
    while (j < order.size() && buckets[order[j]].plane == buckets[order[i]].plane)
    {
      j++;
    }
    for (std::size_t k = i; k + 1 < j; k++)
    {
      if (Parallel(buckets[order[k]].representative, buckets[order[k + 1]].representative))
      {
        uf.Union(order[k], order[k + 1]);
      }
    }
    if (j - i > 2 &&
        Parallel(buckets[order[i]].representative, buckets[order[j - 1]].representative))
    {
      uf.Union(order[i], order[j - 1]);  // the wrap of the direction angle at 0 / pi
    }
    i = j;
  }
  std::map<std::size_t, int> class_of_root;
  for (const std::size_t b : order)
  {
    const std::size_t root = uf.Find(b);
    auto it = class_of_root.find(root);
    if (it == class_of_root.end())
    {
      it = class_of_root.emplace(root, static_cast<int>(class_of_root.size())).first;
    }
    for (const std::size_t r : buckets[b].members)
    {
      run_direction_class[r] = it->second;
    }
  }
  direction_class_count = class_of_root.size();
}

void Identifier::ClassifyVertices()
{
  for (std::size_t v = 0; v < input.vertices.size(); v++)
  {
    if (!IsFeatureVertex(v))
    {
      continue;
    }
    const auto incident = RunsAtVertex(v);
    if (incident.empty())
    {
      continue;  // every incident run is excluded; recorded with the runs
    }
    VertexFeatureSite site;
    site.vertex = v;
    site.point = input.vertices[v].coordinate;
    site.runs_at_site = incident;
    site.boundary_law = runs[incident.front()].boundary_law;
    std::set<std::string> interfaces;
    for (const std::size_t r : incident)
    {
      for (const auto &name : InterfaceNames(runs[r].targets))
      {
        interfaces.insert(name);
      }
    }
    site.interfaces.assign(interfaces.begin(), interfaces.end());
    const auto type = *input.vertices[v].physical_type;
    if (incident.size() == 1 || type == MetalEdgeVertexType::ENDPOINT)
    {
      site.type = "Endpoint";
      site.arm_directions = {ArmDirection(incident.front(), v)};
      site.gap_direction = runs[incident.front()].gap_direction;
    }
    else if (incident.size() == 2)
    {
      const Point3D da = ArmDirection(incident[0], v);
      const Point3D db = ArmDirection(incident[1], v);
      site.arm_directions = {da, db};
      const double wedge =
          std::acos(std::clamp(Dot(da, db), -1.0, 1.0)) * 180.0 / std::acos(-1.0);
      site.turn_degrees = 180.0 - wedge;
      site.angle_degrees = wedge;
      // The gap lies inside the wedge -> concave (metal angle > 180), else convex.
      const Point3D bisector = Add(da, db);
      const Point3D gap =
          Add(runs[incident[0]].gap_direction, runs[incident[1]].gap_direction);
      site.type = Dot(bisector, gap) > 0.0 ? "ConcaveCorner" : "ConvexCorner";
    }
    else
    {
      site.type = "Junction";
      // Arm angles around the site's own signed process normal (the frame normal of the
      // feature, decision 266), as consecutive differences.
      const Point3D n_site = SiteProcessNormal(site);
      const Point3D x =
          Normalize(Sub(ArmDirection(incident[0], v),
                        Scale(Dot(ArmDirection(incident[0], v), n_site), n_site)));
      const Point3D y = Cross(n_site, x);
      std::vector<std::tuple<double, int, std::size_t>> arms;  // angle, conductor, run
      for (const std::size_t r : incident)
      {
        const Point3D d = ArmDirection(r, v);
        double angle = std::atan2(Dot(d, y), Dot(d, x)) * 180.0 / std::acos(-1.0);
        if (angle < 0.0)
        {
          angle += 360.0;
        }
        arms.emplace_back(angle, runs[r].conductor, r);
      }
      std::sort(arms.begin(), arms.end());
      for (std::size_t i = 0; i < arms.size(); i++)
      {
        const double next =
            i + 1 < arms.size() ? std::get<0>(arms[i + 1]) : std::get<0>(arms[0]) + 360.0;
        site.arm_angles.push_back(next - std::get<0>(arms[i]));
        site.arm_conductors.push_back(std::get<1>(arms[i]));
        site.arm_directions.push_back(ArmDirection(std::get<2>(arms[i]), v));
      }
    }
    sites.push_back(std::move(site));
  }
}

// Arc rule (decision 82(3); concyclicity form of USER decisions 121 / 122): fitted arcs
// along the perimeter path through corner vertices. A path continues through every vertex
// with exactly two path segments (regular or corner); it stops at endpoints, junctions and
// cuts. Its joints are the non-collinear vertices; from every unconsumed joint the largest
// range of at least three following joints that (i) turn the same way, each by less than
// kArcMaxJointTurnDegrees, (ii) total at most 180 deg and (iii) lie on one circle within
// kArcFitToleranceOverRadius x R is an arc, whatever its chord sagitta (recorded). The
// corner vertices of an arc make no vertex feature (a 90 deg fillet meshed with two 45 deg
// chords has one in its middle), the arc's curvature is its exact radius and its curved
// features are one across those vertices, so that the description does not depend on the
// number of chords; a polygon with joints of 50 deg or more is corners.
void Identifier::DetectArcs()
{
  const double fit_tolerance = kArcFitToleranceOverRadius * R;
  const double joint_turn_cap = kArcMaxJointTurnDegrees * std::acos(-1.0) / 180.0;
  const double noise_sagitta = kJointNoiseSagittaOverRadius * R;
  arcs.clear();
  vertex_arc.assign(input.vertices.size(), -1);
  segment_arc.assign(input.segments.size(), -1);
  // Path segments and their vertex incidence.
  std::vector<std::vector<std::size_t>> incident(input.vertices.size());
  std::vector<bool> on_path(input.segments.size(), false);
  for (std::size_t i = 0; i < input.segments.size(); i++)
  {
    const auto &segment = input.segments[i];
    if (segment.truncation || segment.exclusion || segment.chain < 0)
    {
      continue;
    }
    on_path[i] = true;
    incident[segment.vertices[0]].push_back(i);
    incident[segment.vertices[1]].push_back(i);
  }
  auto Continues = [&](std::size_t v)
  {
    const auto &vertex = input.vertices[v];
    return incident[v].size() == 2 && vertex.physical_type &&
           (*vertex.physical_type == MetalEdgeVertexType::REGULAR ||
            *vertex.physical_type == MetalEdgeVertexType::CORNER) &&
           !vertex.on_truncation_boundary && !vertex.on_port_boundary;
  };
  auto Other = [&](std::size_t s, std::size_t v)
  {
    const auto &ends = input.segments[s].vertices;
    return ends[0] == v ? ends[1] : ends[0];
  };
  auto Direction = [&](std::size_t s, std::size_t from)
  {
    const auto &segment = input.segments[s];
    return Normalize(segment.vertices[0] == from ? Sub(segment.p1, segment.p0)
                                                 : Sub(segment.p0, segment.p1));
  };
  int next_chain_id = 0;
  for (const auto &segment : input.segments)
  {
    next_chain_id = std::max(next_chain_id, segment.chain + 1);
  }
  std::vector<bool> visited(input.segments.size(), false);
  std::size_t paths = 0;
  for (std::size_t seed = 0; seed < input.segments.size(); seed++)
  {
    if (!on_path[seed] || visited[seed])
    {
      continue;
    }
    stage.Progress(seed, input.segments.size());
    // Walk backwards to the path start (or around a loop), then forwards.
    std::size_t start_segment = seed;
    std::size_t start_vertex = input.segments[seed].vertices[0];
    {
      std::size_t s = seed, v = input.segments[seed].vertices[0];
      std::set<std::size_t> seen = {seed};
      while (Continues(v))
      {
        const std::size_t previous = incident[v][0] == s ? incident[v][1] : incident[v][0];
        if (seen.count(previous))
        {
          break;  // a closed loop: start anywhere
        }
        seen.insert(previous);
        s = previous;
        v = Other(previous, v);
      }
      start_segment = s;
      start_vertex = v;
    }
    std::vector<std::size_t> path_segments;  // in order
    std::vector<std::size_t> path_vertices;  // path_vertices[k] precedes path_segments[k]
    std::size_t s = start_segment, v = start_vertex;
    while (true)
    {
      visited[s] = true;
      path_segments.push_back(s);
      path_vertices.push_back(v);
      v = Other(s, v);
      if (!Continues(v))
      {
        path_vertices.push_back(v);
        break;
      }
      const std::size_t next = incident[v][0] == s ? incident[v][1] : incident[v][0];
      if (visited[next])
      {
        path_vertices.push_back(v);  // loop closed
        break;
      }
      s = next;
    }
    paths++;
    const bool closed = path_vertices.front() == path_vertices.back() &&
                        path_segments.size() > 2 && Continues(path_vertices.front());
    // The arcs of the path in one traversal direction, without side effects (the greedy
    // scan from every unconsumed joint depends on the direction where pseudo-arcs of a
    // spline compete with a true arc; both directions are scanned below and the one
    // absorbing more joints is applied, so that the result does not depend on the input
    // orientation and a mirrored mesh gives the mirrored arcs).
    auto ScanPath = [&](const std::vector<std::size_t> &path_segments,
                        const std::vector<std::size_t> &path_vertices) -> std::vector<Arc>
    {
      std::vector<Arc> found;
      // Pieces (maximal collinear runs of segments) and the joints between them.
      struct Joint
      {
        std::size_t vertex;
        std::size_t index;  // path index of the vertex (segments index-1 | index)
        Point3D in, out;    // unit directions before / after
        double turn;        // |turn| in radians
        int sign;           // about the local process normal
      };
      std::vector<Joint> joints;
      const std::size_t n = path_segments.size();
      for (std::size_t k = closed ? 0 : 1; k < n; k++)
      {
        const std::size_t before = path_segments[(k + n - 1) % n], after = path_segments[k];
        const std::size_t vertex = path_vertices[k];
        const Point3D in = Direction(before, Other(before, vertex));
        const Point3D out = Direction(after, vertex);
        const double dot = std::clamp(Dot(in, out), -1.0, 1.0);
        if (!DirectionLess(dot, 1.0 - kDirectionQuantum))
        {
          continue;  // collinear: no joint
        }
        const Point3D normal = Normalize(Add(input.segments[before].process_normal,
                                             input.segments[after].process_normal));
        const double sign = Dot(Cross(in, out), normal);
        joints.push_back({vertex, k, in, out, std::acos(dot), sign >= 0.0 ? 1 : -1});
      }
      if (joints.size() < 2)
      {
        return found;
      }
      // Path position of every vertex.
      std::vector<double> position(path_vertices.size(), 0.0);
      for (std::size_t k = 0; k < n; k++)
      {
        const auto &segment = input.segments[path_segments[k]];
        position[k + 1] = position[k] + Distance(segment.p0, segment.p1);
      }
      const double path_length = position[n];
      if (closed)
      {
        // Start the cyclic scan at the joint following the longest piece (a straight arm),
        // so that no arc is split by the arbitrary loop start.
        std::size_t best = 0;
        double longest = -1.0;
        for (std::size_t j = 0; j < joints.size(); j++)
        {
          const std::size_t prev = (j + joints.size() - 1) % joints.size();
          double gap = position[joints[j].index] - position[joints[prev].index];
          if (gap <= 0.0)
          {
            gap += path_length;
          }
          if (gap > longest)
          {
            longest = gap;
            best = j;
          }
        }
        std::rotate(joints.begin(), joints.begin() + static_cast<std::ptrdiff_t>(best),
                    joints.end());
      }
      const std::size_t m = joints.size();
      // Algebraic least-squares circle (Kasa fit) through the joint vertices in the process
      // plane: exact for vertices on one circle, a set function of the joints.
      auto LeastSquaresCircle = [&](const std::vector<std::size_t> &joint_vertices,
                                    const Point3D &normal, const Point3D &tangent,
                                    Point3D &center, double &radius)
      {
        const Point3D o = input.vertices[joint_vertices.front()].coordinate;
        const Point3D u = Normalize(Sub(tangent, Scale(Dot(tangent, normal), normal)));
        const Point3D v = Normalize(Cross(normal, u));
        double sxx = 0.0, sxy = 0.0, syy = 0.0, sx = 0.0, sy = 0.0, s1 = 0.0;
        double sxz = 0.0, syz = 0.0, sz = 0.0;
        for (const std::size_t jv : joint_vertices)
        {
          const Point3D d = Sub(input.vertices[jv].coordinate, o);
          const double x = Dot(d, u), y = Dot(d, v), z = x * x + y * y;
          sxx += x * x, sxy += x * y, syy += y * y, sx += x, sy += y, s1 += 1.0;
          sxz += x * z, syz += y * z, sz += z;
        }
        // Normal equations for (D, E, F) of x^2 + y^2 + D x + E y + F = 0 (Cramer's rule).
        const double a11 = sxx, a12 = sxy, a13 = sx, a22 = syy, a23 = sy, a33 = s1;
        const double b1 = -sxz, b2 = -syz, b3 = -sz;
        const double det = a11 * (a22 * a33 - a23 * a23) - a12 * (a12 * a33 - a23 * a13) +
                           a13 * (a12 * a23 - a22 * a13);
        if (std::abs(det) <= 1.0e-30)
        {
          return false;
        }
        const double D = (b1 * (a22 * a33 - a23 * a23) - a12 * (b2 * a33 - a23 * b3) +
                          a13 * (b2 * a23 - a22 * b3)) /
                         det;
        const double E = (a11 * (b2 * a33 - a23 * b3) - b1 * (a12 * a33 - a23 * a13) +
                          a13 * (a12 * b3 - b2 * a13)) /
                         det;
        const double F = (a11 * (a22 * b3 - b2 * a23) - a12 * (a12 * b3 - b2 * a13) +
                          b1 * (a12 * a23 - a22 * a13)) /
                         det;
        const double cx = -0.5 * D, cy = -0.5 * E;
        const double radius2 = cx * cx + cy * cy - F;
        if (!(radius2 > 0.0))
        {
          return false;
        }
        radius = std::sqrt(radius2);
        center = Add(o, Add(Scale(cx, u), Scale(cy, v)));
        return true;
      };
      // Every joint of the range within the fit tolerance of the circle (strict on the
      // grid).
      auto JointsOnCircle = [&](const std::vector<std::size_t> &joint_vertices,
                                const Point3D &center, double radius)
      {
        return std::all_of(
            joint_vertices.begin(), joint_vertices.end(),
            [&](std::size_t jv)
            {
              return quantizer.Less(
                  std::abs(Distance(input.vertices[jv].coordinate, center) - radius),
                  fit_tolerance);
            });
      };
      // The largest sagitta of the chords between consecutive joints of the range on the
      // fitted circle (rho - sqrt(rho^2 - (c / 2)^2); a chord longer than the diameter has
      // none: the range is no arc of that circle) — recorded per arc as the mesh-coarseness
      // diagnostic (Arcs[].MaxChordSagittaOverR), no longer a membership test.
      auto MaxChordSagitta =
          [&](const std::vector<std::size_t> &joint_vertices, double radius, bool cyclic)
      {
        double worst = 0.0;
        const std::size_t count = joint_vertices.size();
        for (std::size_t q = 0; q + (cyclic ? 0 : 1) < count; q++)
        {
          const double chord =
              Distance(input.vertices[joint_vertices[q]].coordinate,
                       input.vertices[joint_vertices[(q + 1) % count]].coordinate);
          if (chord >= 2.0 * radius)
          {
            return std::numeric_limits<double>::infinity();
          }
          worst =
              std::max(worst, radius - std::sqrt(radius * radius - 0.25 * chord * chord));
        }
        return worst;
      };
      // A joint turns less than the arc rule's joint-turn cap (strict on the direction
      // grid).
      auto JointTurnBelowCap = [&](const Joint &joint)
      { return DirectionLess(std::cos(joint_turn_cap), std::cos(joint.turn)); };
      // The straight pieces on either side of joint j in the scan direction (between it and
      // the neighbouring joints, or the path ends of an open path).
      // On a closed path the joint list is rotated (below) so that the scan starts after
      // the longest piece: the pair of consecutive joints straddling the path's original
      // start has a negative position difference, wrapped by the path length (formerly
      // only the pair of the list's first and last joints was wrapped: the pieces of the
      // straddling pair read negative, which made the joint noise rule of the arc tests
      // accept any turn there).
      auto CyclicGap = [&](double gap)
      { return closed && gap <= 0.0 ? gap + path_length : gap; };
      auto PieceBefore = [&](std::size_t j)
      {
        const std::size_t idx = joints[j].index;
        if (j == 0)
        {
          return closed ? CyclicGap(position[idx] - position[joints[m - 1].index])
                        : position[idx];
        }
        return CyclicGap(position[idx] - position[joints[j - 1].index]);
      };
      auto PieceAfter = [&](std::size_t j)
      {
        const std::size_t idx = joints[j].index;
        if (j + 1 == m)
        {
          return closed ? CyclicGap(position[joints[0].index] - position[idx])
                        : path_length - position[idx];
        }
        return CyclicGap(position[joints[j + 1].index] - position[idx]);
      };
      // The concyclicity test on EVERY point of the range (USER decision 184 (2),
      // 2026-10-01): the interior points of a piece between two consecutive joints lie on
      // the chord, off the circle by up to the chord sagitta. That deviation is admissible
      // where it is below the joint noise resolution (kJointNoiseSagittaOverRadius x R, the
      // deviation from straight the noise rule cannot resolve: a finely chorded curve) or
      // where the chord is a design chord of a coarsely discretised curve — one of its end
      // joints is a real turn (not noise under the geometric joint rule on its shorter
      // piece; USER decision 122: a coarse polyline of real joints IS the curve whatever
      // the sagitta, recorded as mesh coarseness). A chord deviating from the circle by
      // more than the resolution between two NOISE joints is a straight edge: a 5,500 um
      // straight trace edge whose ends carry two 1.2 / 2.4 deg joints on 12 um pieces
      // (exactly concyclic by mirror symmetry) read as a 132 mm bend bowing 15 R off the
      // metal (DS-OSC-003 E8-3: false 2-edge clusters over 3.9 mm; DS-SCT-002 E8-4). For an
      // inscribed polyline of equal chords the joint's implied sagitta and the chord
      // sagitta are the same quantity, so a uniformly chorded arc passes iff its joints
      // pass the noise rule as real joints or its chords are below the resolution — no
      // new threshold enters. Chords in the range: consecutive joint pairs (cyclic on a
      // closed circle).
      auto RealJoint = [&](std::size_t j)
      {
        return !JointIsNoise(joints[j].turn, std::min(PieceBefore(j), PieceAfter(j)),
                             noise_sagitta);
      };
      auto ChordsResolved =
          [&](std::size_t i, std::size_t count, double radius, bool cyclic)
      {
        for (std::size_t q = 0; q + (cyclic ? 0 : 1) < count; q++)
        {
          const std::size_t j0 = (i + q) % m, j1 = (i + q + 1) % m;
          const double chord = Distance(input.vertices[joints[j0].vertex].coordinate,
                                        input.vertices[joints[j1].vertex].coordinate);
          if (chord >= 2.0 * radius)
          {
            return false;
          }
          const double sagitta = radius - std::sqrt(radius * radius - 0.25 * chord * chord);
          if (!quantizer.Less(sagitta, noise_sagitta) && !RealJoint(j0) && !RealJoint(j1))
          {
            if (std::getenv("PALACE_IDENTIFICATION_DEBUG_ARCS") && input.log)
            {
              std::ostringstream dbg;
              dbg << std::setprecision(10) << "  DEBUG unresolved chord: range " << i << "+"
                  << count << " joints " << j0 << " (turn "
                  << joints[j0].turn * 180.0 / std::acos(-1.0) << " deg, pieces "
                  << PieceBefore(j0) << " / " << PieceAfter(j0) << ") -> " << j1
                  << " (turn " << joints[j1].turn * 180.0 / std::acos(-1.0)
                  << " deg, pieces " << PieceBefore(j1) << " / " << PieceAfter(j1)
                  << ") chord " << chord << " sagitta " << sagitta / R << " R radius "
                  << radius / R << " R at ("
                  << input.vertices[joints[j0].vertex].coordinate[0] << ", "
                  << input.vertices[joints[j0].vertex].coordinate[1] << ") m " << m
                  << (closed ? " closed" : " open") << " positions "
                  << position[joints[j0].index] << " / " << position[joints[j1].index]
                  << " path length " << path_length << " vertices " << joints[j0].vertex
                  << " / " << joints[j1].vertex << " j1 at ("
                  << input.vertices[joints[j1].vertex].coordinate[0] << ", "
                  << input.vertices[joints[j1].vertex].coordinate[1] << ")\n";
              input.log(dbg.str());
            }
            return false;
          }
        }
        return true;
      };
      // Fit the arc over joints [i, i + count) (cyclic indices on a loop); returns the
      // circle.
      struct Fit
      {
        bool ok = false;
        bool tangent =
            false;  // the arms are tangent to the circle (a rounded corner may be)
        Point3D center{}, origin{};
        double radius = 0.0;
        double turn = 0.0;
        double sagitta = 0.0;  // the largest chord sagitta
      };
      std::vector<bool> consumed(m, false);
      auto TryFit = [&](std::size_t i, std::size_t count) -> Fit
      {
        Fit fit;
        const Joint &first = joints[i];
        const Joint &last = joints[(i + count - 1) % m];
        const Point3D ta = first.in, tb = last.out;
        const Point3D Ta = input.vertices[first.vertex].coordinate;
        const Point3D Tb = input.vertices[last.vertex].coordinate;
        double turn = 0.0;
        std::vector<std::size_t> range_vertices;
        for (std::size_t j = 0; j < count; j++)
        {
          turn += joints[(i + j) % m].turn;
          range_vertices.push_back(joints[(i + j) % m].vertex);
        }
        const double angle = std::acos(std::clamp(Dot(ta, tb), -1.0, 1.0));
        // The arm tangents must enclose the accumulated turn (a monotone arc below 180
        // deg).
        if (std::abs(angle - turn) > 1.0e-6 &&
            std::abs((2.0 * std::acos(-1.0) - angle) - turn) > 1.0e-6)
        {
          return fit;
        }
        const Point3D normal = Normalize(
            Add(input.segments[path_segments[first.index % n]].process_normal,
                input.segments[path_segments[(first.index + n - 1) % n]].process_normal));
        const Point3D na =
            Scale(static_cast<double>(first.sign), Normalize(Cross(normal, ta)));
        // The circle tangent to both arms at the end joints (a fillet, a rounded corner, a
        // bend between tangent arms): centre Ta + radius na, radius from the tangent
        // lengths. Every joint, the last one included, must lie on it within the fit
        // tolerance (unequal tangent lengths put Tb off the circle).
        double radius = 0.0;
        Point3D center{}, origin{};
        bool tangent_circle = false;
        const double sin_turn = std::sin(turn);
        if (std::abs(sin_turn) > 1.0e-9 && turn < std::acos(-1.0) - 1.0e-9)
        {
          // Virtual corner X = Ta + a ta = Tb - b tb, solved in the (ta, na) frame.
          const Point3D w = Sub(Tb, Ta);
          const double wx = Dot(w, ta), wy = Dot(w, na);
          const double tbx = Dot(tb, ta), tby = Dot(tb, na);
          if (std::abs(tby) > 1.0e-12)
          {
            const double b = wy / tby;
            const double a = wx - b * tbx;
            if (a > 0.0 && b > 0.0)
            {
              radius = 0.5 * (a + b) / std::tan(0.5 * turn);
              origin = Add(Ta, Scale(a, ta));
              center = Add(Ta, Scale(radius, na));
              tangent_circle = radius > 0.0;
            }
          }
        }
        else
        {
          // Antiparallel arms (a U-turn): the radius is half the arm separation and the
          // tangent points face each other across it.
          const Point3D w = Sub(Tb, Ta);
          const double across = Dot(w, na);
          if (across > 0.0)
          {
            radius = 0.5 * across;
            center = Add(Ta, Scale(radius, na));
            origin = Add(center, Scale(radius, ta));  // the arc's midpoint
            tangent_circle = true;
          }
        }
        if (tangent_circle && JointsOnCircle(range_vertices, center, radius) &&
            ChordsResolved(i, count, radius, false))
        {
          const double sagitta = MaxChordSagitta(range_vertices, radius, false);
          if (std::isfinite(sagitta))
          {
            fit.ok = true;
            fit.tangent = true;
            fit.center = center;
            fit.origin = origin;
            fit.radius = radius;
            fit.turn = turn;
            fit.sagitta = sagitta;
          }
          return fit;
        }
        // Arms not tangent to the circle (a route bend preceded by a spline piece,
        // DS-SCT-001: 0.5 % radius error and a 0.8 um centre offset of the tangent
        // construction on a 370 um bend): the least-squares circle of the joints, exact for
        // an inscribed polyline and a set function of the joints, whose radius is at least
        // R — a bend inside its chain. A rounded corner (radius below R) keeps the tangent
        // construction: its arms are tangent by construction, and its site (virtual corner,
        // arm directions) has no meaning otherwise. Without the arms the circle is
        // constrained by the joints alone: at least FOUR (three points are concyclic
        // whatever they are, and the inscribed-angle relation tests interior joints only —
        // a 90 deg lead-end corner, the first joint of an arc and its neighbour passed as a
        // 3-joint "bend" of radius 350 um whose 6 um lead chord had a sagitta of 0.006 R).
        if (count < 4)
        {
          return fit;
        }
        // The end joints of a least-squares bend: the kink between the arm and the circle's
        // tangent at the end joint (read from the fitted circle: a 1e-3 R position error
        // moves the tangent by 1e-3 R / rho, far less than a chord-based estimate on short
        // chords) must be noise under the geometric joint rule on the shorter of the arm
        // piece and the first chord (a spline piece joining the bend tangent-continuously;
        // a mitred offset polyline) — else the arm is a chord of the same circle (its far
        // vertex on the circle within the fit tolerance: an arc starting at a corner that
        // lies on its circle, the sharp end of a rounded slot; the corner stays a corner,
        // the arc starts at its next joint). A lead meeting a circular arc at a 90 deg
        // corner is neither and the corner is never absorbed as a bend vertex. The arm's
        // far vertex is the far end of its straight PIECE (the neighbouring joint, or the
        // path end): the rigid-run joint, not the far end of the adjacent mesh segment —
        // inserting collinear vertices on the chord is a geometric no-op and must not
        // change the reading (decision 212; VALIDATION-PLAN (h)-8: a coarse round pad whose
        // lead attaches through joints above the turn cap read a bend over its interior
        // joints on the chip mesh, one mesh edge per design chord, and sharp corners on a
        // thin window mesh subdividing the same chords at 4 um, whose first sub-vertex lies
        // off the circle by most of the chord sagitta).
        auto ArmFarVertex = [&](std::size_t j, bool before)
        {
          if (before)
          {
            return j == 0 && !closed ? path_vertices.front()
                                     : joints[(j + m - 1) % m].vertex;
          }
          return j + 1 == m && !closed ? path_vertices.back() : joints[(j + 1) % m].vertex;
        };
        // The kink between an arm and the circle's tangent at an end joint is noise under
        // the geometric joint rule on the shorter of the arm piece and the first chord (a
        // tangent arm); at_start orients the circle tangent along the path (toward the
        // arc's interior at the first joint, away from it at the last; the arm direction is
        // the path direction there in both cases).
        auto ArmKinkIsNoise = [&](const Joint &end, const Point3D &arm_direction,
                                  std::size_t neighbour_vertex, bool at_start,
                                  const Point3D &center, double arm_piece,
                                  double first_chord)
        {
          const Point3D at = input.vertices[end.vertex].coordinate;
          Point3D circle_tangent = Normalize(Cross(normal, Sub(at, center)));
          const Point3D toward_neighbour =
              Sub(input.vertices[neighbour_vertex].coordinate, at);
          const double along = Dot(circle_tangent, toward_neighbour);
          if (at_start != (along > 0.0))
          {
            circle_tangent = Scale(-1.0, circle_tangent);
          }
          const double kink =
              std::acos(std::clamp(Dot(arm_direction, circle_tangent), -1.0, 1.0));
          return JointIsNoise(kink, std::min(arm_piece, first_chord), noise_sagitta);
        };
        // The arm is a chord of the circle: its far vertex lies on it within the fit
        // tolerance.
        auto ArmFarOnCircle =
            [&](std::size_t arm_far_vertex, const Point3D &center, double radius)
        {
          const Point3D far = input.vertices[arm_far_vertex].coordinate;
          return quantizer.Less(std::abs(Distance(far, center) - radius), fit_tolerance);
        };
        // A chord arm at the FIRST joint whose far joint p is itself absorbable —
        // unconsumed, the same sign, below the turn cap, and its own arm consistent with
        // the circle (tangent or a chord; p lies on the circle: that is the chord-arm
        // clause itself) — is no arm: the range is not maximal at its start. Without this a
        // scan starting INSIDE an arc (a closed loop whose longest piece is a chord of the
        // arc, the recorded exposure of the closed-loop start rule) accepted a sub-arc
        // anchored at the loop start and chopped the arc there (joint-only meshes; on a
        // subdivided mesh the sub-vertex failed the former mesh-segment test and the whole
        // arc was found from its real start: the two discretisations disagreed). The scan
        // reaches the arc's real start later (cyclically on a loop). A corner on the circle
        // whose own arm is neither tangent nor a chord (a lead meeting a round pad at a
        // sub-cap angle) is not absorbable: the arc still starts at its next joint and the
        // corner stays a corner. At the LAST joint the next joint is a chord arm as before
        // (a polyline turning more than 180 deg is cut there under the current rule).
        auto FirstArmJointAbsorbable =
            [&](std::size_t j, const Point3D &center, double radius)
        {
          if (j == 0 && !closed)
          {
            return false;
          }
          const std::size_t p = (j + m - 1) % m;
          return !consumed[p] && joints[p].sign == joints[j].sign &&
                 JointTurnBelowCap(joints[p]) &&
                 (ArmKinkIsNoise(joints[p], joints[p].in, joints[j].vertex, true, center,
                                 PieceBefore(p), PieceAfter(p)) ||
                  ArmFarOnCircle(ArmFarVertex(p, true), center, radius));
        };
        auto FirstJointConsistent = [&](std::size_t j, const Point3D &center, double radius)
        {
          return ArmKinkIsNoise(joints[j], joints[j].in, range_vertices[1], true, center,
                                PieceBefore(j), PieceAfter(j)) ||
                 (!FirstArmJointAbsorbable(j, center, radius) &&
                  ArmFarOnCircle(ArmFarVertex(j, true), center, radius));
        };
        auto LastJointConsistent = [&](std::size_t j, const Point3D &center, double radius)
        {
          return ArmKinkIsNoise(joints[j], joints[j].out, range_vertices[count - 2], false,
                                center, PieceAfter(j), PieceBefore(j)) ||
                 ArmFarOnCircle(ArmFarVertex(j, false), center, radius);
        };
        if (LeastSquaresCircle(range_vertices, normal, first.in, center, radius) &&
            !quantizer.Less(radius, R) && JointsOnCircle(range_vertices, center, radius) &&
            ChordsResolved(i, count, radius, false) &&
            FirstJointConsistent(i, center, radius) &&
            LastJointConsistent((i + count - 1) % m, center, radius))
        {
          const double sagitta = MaxChordSagitta(range_vertices, radius, false);
          if (std::isfinite(sagitta))
          {
            fit.ok = true;
            fit.tangent = false;
            fit.center = center;
            fit.origin = center;
            fit.radius = radius;
            fit.turn = turn;
            fit.sagitta = sagitta;
          }
        }
        return fit;
      };
      // A closed path whose joints all turn one way through 360 deg and lie on one circle
      // (a round pad, hole or via) is ONE arc of total turn 2 pi: a bend of exact radius
      // whatever its radius (the rounded-corner semantics of arms meeting through a fillet
      // do not apply to a closed circle; a circle of radius < R is one CurvedEdge with
      // RadiusOverR < 1 instead of two 180 deg "rounded corners" split at a
      // numbering-dependent joint). The joint-turn cap tells a circle from a polygon: a
      // square or hexagonal hole is corners whatever its size, an octagon (45 deg per
      // joint) is a circle.
      if (closed && m >= 3)
      {
        double total = 0.0;
        bool same_sign = true;
        bool below_cap = true;
        for (std::size_t j = 0; j < m; j++)
        {
          total += joints[j].turn;
          same_sign = same_sign && joints[j].sign == joints.front().sign;
          below_cap = below_cap && JointTurnBelowCap(joints[j]);
        }
        if (same_sign && below_cap && std::abs(total - 2.0 * std::acos(-1.0)) < 1.0e-6)
        {
          std::vector<std::size_t> all_joints;
          for (const auto &joint : joints)
          {
            all_joints.push_back(joint.vertex);
          }
          const Point3D normal = Normalize(
              Add(input.segments[path_segments[joints.front().index % n]].process_normal,
                  input.segments[path_segments[(joints.front().index + n - 1) % n]]
                      .process_normal));
          Point3D center{};
          double radius = 0.0;
          double sagitta = 0.0;
          if (LeastSquaresCircle(all_joints, normal, joints.front().in, center, radius) &&
              JointsOnCircle(all_joints, center, radius) &&
              ChordsResolved(0, m, radius, true) &&
              std::isfinite(sagitta = MaxChordSagitta(all_joints, radius, true)))
          {
            Arc arc;
            arc.joints = all_joints;
            arc.segments = path_segments;
            arc.segment_before = path_segments.front();
            arc.segment_after = path_segments.front();
            arc.tangent_a = joints.front().in;
            arc.tangent_b = joints.front().in;
            arc.center = center;
            arc.origin = center;
            arc.radius = radius;
            arc.turn = 2.0 * std::acos(-1.0);
            arc.max_sagitta = sagitta;
            arc.corner = false;
            found.push_back(std::move(arc));
            std::fill(consumed.begin(), consumed.end(), true);
          }
        }
      }
      const std::size_t first_start = 0;
      for (std::size_t i = first_start; i < m; i++)
      {
        if (consumed[i] || !JointTurnBelowCap(joints[i]))
        {
          continue;
        }
        // Extend while the turn is monotone and at most 180 deg; the largest range that
        // fits wins (a perturbed joint inside a fillet is absorbed by a larger range whose
        // circle it lies on; a spline scanned past its exact-fit range simply fails the
        // later fits).
        std::size_t best_count = 0;
        Fit best;
        double turn = joints[i].turn;
        // At least three joints (two chords): a polyline with two joints is a chamfer or a
        // square strip end, which no test can tell from a one-chord arc (the diameter chord
        // of a semicircle is the square end): those stay corners.
        for (std::size_t count = 2; count <= m; count++)
        {
          const std::size_t k = (i + count - 1) % m;
          if ((!closed && i + count - 1 >= m) || k == i || consumed[k] ||
              joints[k].sign != joints[i].sign || !JointTurnBelowCap(joints[k]))
          {
            break;
          }
          turn += joints[k].turn;
          if (turn > std::acos(-1.0) + 1.0e-9)
          {
            break;
          }
          if (count < 3)
          {
            continue;
          }
          const Fit fit = TryFit(i, count);
          if (fit.ok)
          {
            best = fit;
            best_count = count;
          }
        }
        if (best_count == 0)
        {
          continue;
        }
        Arc arc;
        for (std::size_t j = 0; j < best_count; j++)
        {
          const std::size_t k = (i + j) % m;
          consumed[k] = true;
          arc.joints.push_back(joints[k].vertex);
        }
        const Joint &first = joints[i];
        const Joint &last = joints[(i + best_count - 1) % m];
        for (std::size_t k = first.index; k != last.index; k = (k + 1) % n)
        {
          arc.segments.push_back(path_segments[k]);
        }
        arc.segment_before = path_segments[(first.index + n - 1) % n];
        arc.segment_after = path_segments[last.index % n];
        arc.tangent_a = first.in;
        arc.tangent_b = last.out;
        arc.center = best.center;
        arc.origin = best.origin;
        arc.radius = best.radius;
        arc.turn = best.turn;
        arc.max_sagitta = best.sagitta;
        // A rounded corner is an arc of radius below R (tangent arms by construction);
        // every arc of radius >= R is a bend, whatever its turn.
        arc.corner = best.tangent && quantizer.Less(arc.radius, R);
        if (std::getenv("PALACE_IDENTIFICATION_DEBUG_ARCS") && input.log)
        {
          std::ostringstream line;
          line << "    arc joints " << arc.joints.size() << " radius " << arc.radius << " ("
               << arc.radius / R << " R) turn " << arc.turn * 180.0 / std::acos(-1.0)
               << " deg sagitta " << arc.max_sagitta / R << " R "
               << (best.tangent ? "tangent" : "least-squares")
               << (arc.corner ? " rounded corner" : " bend") << "\n";
          input.log(line.str());
        }
        found.push_back(std::move(arc));
      }
      return found;
    };
    // Tie-break serialisation of an arc set: per arc the radius, the total turn, the joint
    // count, the distance of its centre from the path's joint centroid and the arc-length
    // position of its FIRST joint in the scan direction measured from the nearer path end,
    // on the signature grid. Radius, turn, joint count and centre distance are invariant
    // under translation, rotation, mirroring and path reversal; the first-joint position is
    // not a set function of the geometric arc (the reversed scan reads the arc from its
    // other end), but the PAIR of serialisations {forward, backward} is invariant: a
    // mirrored or reversed path scanned the other way gives key(forward') = key(backward),
    // so the comparison below picks the congruent arc set whatever the input orientation
    // (verified by the mirror gate and mesh-level mirror / rotation checks; review fix-4
    // m-C). The former serialisation of absolute centre coordinates made the choice between
    // two equally absorbing arc sets depend on the input orientation: DS-SCT-001 mirrored
    // read another chopping of a spline bend into exact-fit arcs, moved a curved boundary
    // 0.036 um and a cluster boundary 0.03 um (review fix-3 m2). On a CLOSED path both
    // scans start at the joint after the longest piece (ties: the first such joint from the
    // seed vertex, which the canonical numbering makes coordinate-dependent) and from_end
    // is 0 for every arc; a different start point of a loop is not covered by the
    // two-direction scan: the exposure (a loop whose longest piece is a chord inside a
    // bend) is closed for bends of four or more joints by the first-joint absorbability
    // clause of the end-joint test (decision 213; rotation variants of the mirror and
    // subdivision gates).
    Point3D path_centroid{};
    std::vector<double> path_position(path_vertices.size(), 0.0);
    {
      const std::size_t distinct = closed ? path_vertices.size() - 1 : path_vertices.size();
      for (std::size_t k = 0; k < distinct; k++)
      {
        path_centroid = Add(path_centroid, input.vertices[path_vertices[k]].coordinate);
      }
      path_centroid = Scale(1.0 / static_cast<double>(distinct), path_centroid);
      for (std::size_t k = 0; k < path_segments.size(); k++)
      {
        const auto &segment = input.segments[path_segments[k]];
        path_position[k + 1] = path_position[k] + Distance(segment.p0, segment.p1);
      }
    }
    const double path_total = path_position.back();
    auto Score = [&](const std::vector<Arc> &list)
    {
      std::size_t absorbed = 0;
      std::string serial;
      for (const auto &arc : list)
      {
        absorbed += arc.joints.size();
      }
      std::vector<std::string> keys;
      for (const auto &arc : list)
      {
        double x_first = 0.0;
        for (std::size_t k = 0; k < path_vertices.size(); k++)
        {
          if (path_vertices[k] == arc.joints.front())
          {
            x_first = path_position[k];
            break;
          }
        }
        const double from_end = closed ? 0.0 : std::min(x_first, path_total - x_first);
        std::ostringstream key;
        key << RoundTo(arc.radius / R, kSignatureLengthQuantumOverRadius) << ","
            << RoundTo(arc.turn * 180.0 / std::acos(-1.0), kSignatureAngleQuantumDegrees)
            << "," << arc.joints.size() << ","
            << RoundTo(Distance(arc.center, path_centroid) / R,
                       kSignatureLengthQuantumOverRadius)
            << "," << RoundTo(from_end / R, kSignatureLengthQuantumOverRadius);
        keys.push_back(key.str());
      }
      std::sort(keys.begin(), keys.end());
      for (const auto &key : keys)
      {
        serial += key + ";";
      }
      return std::make_tuple(absorbed, list.size(), serial);
    };
    std::vector<Arc> forward = ScanPath(path_segments, path_vertices);
    std::vector<std::size_t> reversed_segments(path_segments.rbegin(),
                                               path_segments.rend());
    std::vector<std::size_t> reversed_vertices(path_vertices.rbegin(),
                                               path_vertices.rend());
    std::vector<Arc> backward = ScanPath(reversed_segments, reversed_vertices);
    const auto score_forward = Score(forward), score_backward = Score(backward);
    // More joints absorbed, then FEWER arcs (one arc over a perturbed fillet rather than
    // two sub-ranges), then the smaller serialisation (a set function).
    std::vector<Arc> &chosen =
        std::make_tuple(std::get<0>(score_backward),
                        -static_cast<long long>(std::get<1>(score_backward)),
                        std::get<2>(score_forward)) >
                std::make_tuple(std::get<0>(score_forward),
                                -static_cast<long long>(std::get<1>(score_forward)),
                                std::get<2>(score_backward))
            ? backward
            : forward;
    for (Arc &arc : chosen)
    {
      if (arc.corner)
      {
        // A rounded corner is a corner: the arc is a chain of its own between its tangent
        // points and the arms on either side are distinct chains that meet through it
        // (ThroughZones), so that the pair and event rules read a filleted corner like a
        // sharp one whatever the chords (a U-turned narrow strip keeps its strip pair).
        arc.chain = next_chain_id++;
        for (const std::size_t s : arc.segments)
        {
          input.segments[s].chain = arc.chain;
        }
      }
      else
      {
        // A bend continues its edge: it stays inside the chain (the pairs along a route are
        // read on whole chains), and the chains a corner vertex inside it separated are
        // merged below, so that a coarse bend with a super-threshold joint is the same
        // chain as a fine one.
        for (const std::size_t s : arc.segments)
        {
          segment_arc[s] = static_cast<int>(arcs.size());
        }
      }
      for (const std::size_t v : arc.joints)
      {
        vertex_arc[v] = static_cast<int>(arcs.size());
        auto &vertex = input.vertices[v];
        if (vertex.physical_type && *vertex.physical_type == MetalEdgeVertexType::CORNER)
        {
          // Absorbed corner: no vertex site (the arc accounts for its turn). The chains on
          // either side stay separate (a chain never folds back onto itself: the pair rules
          // pair distinct chains); the arc's curved features are merged across it below.
          vertex.physical_type = MetalEdgeVertexType::REGULAR;
          arc.absorbed_corners.push_back(v);
          if (!arc.corner)
          {
            const int ca = ChainGroup(input.segments[incident[v][0]].chain);
            const int cb = ChainGroup(input.segments[incident[v][1]].chain);
            chain_group.try_emplace(ca, ca);
            chain_group.try_emplace(cb, cb);
            if (ca != cb)
            {
              chain_group[std::max(ca, cb)] = std::min(ca, cb);
            }
          }
        }
      }
      arcs.push_back(std::move(arc));
    }
  }
  if (!chain_group.empty())
  {
    for (auto &segment : input.segments)
    {
      if (segment.chain >= 0 && chain_group.count(segment.chain))
      {
        segment.chain = ChainGroup(segment.chain);
      }
    }
  }
}

// Rounded-corner sites of the arcs of radius below R with a turn above the corner
// threshold: one site per arc (total turn, fillet radius) claiming the arc runs and R along
// each arm from the tangent points; the arc joints take no part in the curvature
// (BuildArcSites runs after BuildRuns / BuildRunIndex).
void Identifier::BuildArcSites()
{
  through_arc.clear();
  for (std::size_t a = 0; a < arcs.size(); a++)
  {
    Arc &arc = arcs[a];
    const long long int before = run_of_segment[arc.segment_before];
    const long long int after = run_of_segment[arc.segment_after];
    if (before < 0 || after < 0)
    {
      continue;  // an arm excluded after the fact (non-planar): no site
    }
    if (!arc.corner)
    {
      continue;
    }
    corner_arc_chains.insert(arc.chain);
    {
      const int ca = runs[static_cast<std::size_t>(before)].chain;
      const int cb = runs[static_cast<std::size_t>(after)].chain;
      if (ca != cb)
      {
        through_arc[std::make_pair(std::min(ca, cb), std::max(ca, cb))].push_back(
            static_cast<int>(a));
      }
    }
    const Run &arm_a = runs[static_cast<std::size_t>(before)];
    const Run &arm_b = runs[static_cast<std::size_t>(after)];
    const std::size_t Ta = arc.joints.front(), Tb = arc.joints.back();
    VertexFeatureSite site;
    site.point = arc.origin;
    site.turn_degrees = arc.turn * 180.0 / std::acos(-1.0);
    site.angle_degrees = 180.0 - site.turn_degrees;
    site.corner_radius = arc.radius;
    const Point3D ta = arc.tangent_a, tb = arc.tangent_b;
    // Convex when the arc turns toward the metal: its centre lies on the metal side of arm
    // A (opposite to the gap direction); well defined for a U-turn, where the bisector test
    // of a sharp corner degenerates.
    const Point3D toward_center = Normalize(Sub(arc.center, input.vertices[Ta].coordinate));
    site.type =
        Dot(arm_a.gap_direction, toward_center) < 0.0 ? "ConvexCorner" : "ConcaveCorner";
    site.arm_directions = {Scale(-1.0, ta), tb};  // away from the virtual corner
    site.boundary_law = arm_a.boundary_law;
    std::set<std::string> interfaces;
    for (const auto &name : InterfaceNames(arm_a.targets))
    {
      interfaces.insert(name);
    }
    for (const auto &name : InterfaceNames(arm_b.targets))
    {
      interfaces.insert(name);
    }
    site.interfaces.assign(interfaces.begin(), interfaces.end());
    // Claims: the arc runs completely (the runs of the arc segments), R along each arm from
    // its tangent point.
    std::set<std::size_t> arc_runs;
    for (const std::size_t s : arc.segments)
    {
      if (run_of_segment[s] >= 0)
      {
        arc_runs.insert(static_cast<std::size_t>(run_of_segment[s]));
      }
    }
    for (const std::size_t r : arc_runs)
    {
      site.window.emplace_back(r, Interval{0.0, runs[r].length});
      site.runs_at_site.push_back(r);
    }
    auto ArmWindow = [&](const Run &arm, std::size_t r, std::size_t tangent_vertex)
    {
      if (arm.end_vertex == tangent_vertex)
      {
        site.window.emplace_back(r, Interval{std::max(0.0, arm.length - R), arm.length});
      }
      else
      {
        site.window.emplace_back(r, Interval{0.0, std::min(R, arm.length)});
      }
      site.runs_at_site.push_back(r);
    };
    ArmWindow(arm_a, static_cast<std::size_t>(before), Ta);
    ArmWindow(arm_b, static_cast<std::size_t>(after), Tb);
    arc.feature = static_cast<int>(sites.size());  // site index; the feature id follows
    sites.push_back(std::move(site));
  }
}

// ---------------------------------------------------------------------------------------
// Curved-edge chain rule: windowed curvature along every chain
// ---------------------------------------------------------------------------------------

std::size_t Identifier::RunIndexInChain(const Chain &chain, std::size_t run) const
{
  MFEM_VERIFY(runs[run].chain == chain.id && runs[run].index_in_chain < chain.runs.size() &&
                  chain.runs[runs[run].index_in_chain] == run,
              "Run missing from its chain!");
  return runs[run].index_in_chain;
}

// Windowed curvature at chain position x from the piecewise-linear nodes.
double Identifier::WindowedCurvature(const Chain &chain, double x) const
{
  const auto &nodes = chain.kappa_nodes;
  if (nodes.empty())
  {
    return 0.0;
  }
  if (x <= nodes.front().first)
  {
    return nodes.front().second;
  }
  if (x >= nodes.back().first)
  {
    return nodes.back().second;
  }
  const auto upper =
      std::upper_bound(nodes.begin(), nodes.end(), std::make_pair(x, 0.0),
                       [](const auto &a, const auto &b) { return a.first < b.first; });
  const auto lower = std::prev(upper);
  const double span = upper->first - lower->first;
  if (span <= 0.0)
  {
    return std::max(lower->second, upper->second);
  }
  const double w = (x - lower->first) / span;
  return (1.0 - w) * lower->second + w * upper->second;
}

// Maximum of the piecewise-linear windowed curvature over [x0, x1] (nodes inside and the
// two ends).
double Identifier::MaxCurvature(const Chain &chain, double x0, double x1) const
{
  if (x1 < x0)
  {
    std::swap(x0, x1);
  }
  double best = std::max(WindowedCurvature(chain, x0), WindowedCurvature(chain, x1));
  // The nodes are sorted by position: only those strictly inside (x0, x1).
  const auto &nodes = chain.kappa_nodes;
  auto first =
      std::upper_bound(nodes.begin(), nodes.end(), std::make_pair(x0, 0.0),
                       [](const auto &a, const auto &b) { return a.first < b.first; });
  for (auto it = first; it != nodes.end() && it->first < x1; ++it)
  {
    best = std::max(best, it->second);
  }
  return best;
}

double Identifier::SignedTurn(const Chain &chain, double x0, double x1) const
{
  if (x1 < x0)
  {
    std::swap(x0, x1);
  }
  const auto &nodes = chain.kappa_nodes;
  if (nodes.size() < 2 || x1 <= x0)
  {
    return 0.0;
  }
  // Piecewise-linear signed curvature on the nodes; exact trapezoids on [x0, x1].
  auto At = [&](std::size_t i, double x)
  {
    const double span = nodes[i + 1].first - nodes[i].first;
    const double w = span > 0.0 ? std::clamp((x - nodes[i].first) / span, 0.0, 1.0) : 0.0;
    return (1.0 - w) * chain.signed_kappa_nodes[i] + w * chain.signed_kappa_nodes[i + 1];
  };
  double turn = 0.0;
  for (std::size_t i = 0; i + 1 < nodes.size(); i++)
  {
    const double a = std::max(x0, nodes[i].first), b = std::min(x1, nodes[i + 1].first);
    if (b > a)
    {
      turn += 0.5 * (At(i, a) + At(i, b)) * (b - a);
    }
  }
  return turn;
}

void Identifier::AccumulateSignedCurvature(const Chain &chain, double x0, double x1,
                                           SignedCurvatureExtremes &extremes) const
{
  if (x1 < x0)
  {
    std::swap(x0, x1);
  }
  const auto &nodes = chain.kappa_nodes;
  if (nodes.empty())
  {
    return;
  }
  auto Read = [&](double x)
  {
    const double value = SignedWindowedCurvature(chain, x);
    extremes.toward_metal = std::max(extremes.toward_metal, value);
    extremes.toward_gap = std::max(extremes.toward_gap, -value);
  };
  Read(x0);
  Read(x1);
  auto first =
      std::upper_bound(nodes.begin(), nodes.end(), std::make_pair(x0, 0.0),
                       [](const auto &a, const auto &b) { return a.first < b.first; });
  for (auto it = first; it != nodes.end() && it->first < x1; ++it)
  {
    Read(it->first);
  }
}

double Identifier::SignedWindowedCurvature(const Chain &chain, double x) const
{
  const auto &nodes = chain.kappa_nodes;
  if (nodes.empty())
  {
    return 0.0;
  }
  const auto upper =
      std::upper_bound(nodes.begin(), nodes.end(), std::make_pair(x, 0.0),
                       [](const auto &a, const auto &b) { return a.first < b.first; });
  if (upper == nodes.begin())
  {
    return chain.signed_kappa_nodes.front();
  }
  if (upper == nodes.end())
  {
    return chain.signed_kappa_nodes.back();
  }
  const auto i = static_cast<std::size_t>(upper - nodes.begin()) - 1;
  const double span = nodes[i + 1].first - nodes[i].first;
  const double w = span > 0.0 ? (x - nodes[i].first) / span : 0.0;
  return (1.0 - w) * chain.signed_kappa_nodes[i] + w * chain.signed_kappa_nodes[i + 1];
}

std::string Identifier::ConvexityName(const SignedCurvatureExtremes &extremes) const
{
  if (std::max(extremes.toward_metal, extremes.toward_gap) <= 0.0)
  {
    return "";
  }
  // The same reading as the curved-section rule: a sense is a bend of the class when its
  // windowed bend radius is (quantized strictly) below the straight threshold.
  const double straight_radius = kStraightBendRadiusOverRadius * R;
  auto Curved = [&](double kappa)
  { return kappa > 0.0 && quantizer.Less(1.0 / kappa, straight_radius); };
  if (Curved(extremes.toward_metal) && Curved(extremes.toward_gap))
  {
    return "Mixed";
  }
  return extremes.toward_gap > extremes.toward_metal ? "Concave" : "Convex";
}

bool Identifier::IsCurvedAt(const Chain &chain, double x, std::size_t *section) const
{
  // The curved intervals are disjoint and ascending: the first whose end reaches x is the
  // only candidate (an earlier one ends before x, a later one starts after that end).
  const auto j =
      static_cast<std::size_t>(std::lower_bound(chain.curved.begin(), chain.curved.end(), x,
                                                [&](const Interval &c, double value)
                                                { return c.second + Tol() < value; }) -
                               chain.curved.begin());
  if (j < chain.curved.size() && x >= chain.curved[j].first - Tol())
  {
    if (section)
    {
      *section = j;
    }
    return true;
  }
  return false;
}

double Identifier::ChainArcDistance(const Chain &chain, double x0, double x1,
                                    double x) const
{
  double distance = x < x0 ? x0 - x : (x > x1 ? x - x1 : 0.0);
  if (chain.closed && chain.length > 0.0)
  {
    // Wrapping: the position x shifted by one period on either side.
    for (const double shift : {-chain.length, chain.length})
    {
      const double xs = x + shift;
      distance = std::min(distance, xs < x0 ? x0 - xs : (xs > x1 ? xs - x1 : 0.0));
    }
  }
  return distance;
}

Point3D Identifier::ChainAt(const Chain &chain, double x, std::size_t *run_out,
                            double *s_out) const
{
  const std::size_t k =
      std::min(static_cast<std::size_t>(
                   std::upper_bound(chain.run_offset.begin(), chain.run_offset.end(), x) -
                   chain.run_offset.begin()) -
                   1,
               chain.runs.size() - 1);
  const std::size_t r = chain.runs[k];
  const double s = std::clamp(x - chain.run_offset[k], 0.0, runs[r].length);
  if (run_out)
  {
    *run_out = r;
  }
  if (s_out)
  {
    *s_out = s;
  }
  return runs[r].At(s);
}

Identifier::ChainPoint
Identifier::ClosestPointOnChain(const Chain &chain, const Point3D &p,
                                std::optional<double> exclude_x,
                                std::optional<double> max_distance) const
{
  const double neighbourhood = kSelfPairNeighbourhoodOverRadius * R;
  // The part of run k at least the self-pair neighbourhood of arc length from exclude_x
  // (point-wise: the former run-level exclusion dropped a whole 28 um leg of a hairpin for
  // every point of its fold, which then had no partner at all), as run-parameter intervals.
  auto Allowed = [&](std::size_t k)
  {
    const double x0 = chain.run_offset[k], x1 = x0 + runs[chain.runs[k]].length;
    std::vector<Interval> allowed = {Interval{0.0, x1 - x0}};
    if (!exclude_x)
    {
      return allowed;
    }
    std::vector<Interval> near;
    for (const double shift : chain.closed
                                  ? std::vector<double>{-chain.length, 0.0, chain.length}
                                  : std::vector<double>{0.0})
    {
      const double lo = std::max(x0, *exclude_x + shift - neighbourhood) - x0;
      const double hi = std::min(x1, *exclude_x + shift + neighbourhood) - x0;
      if (hi - lo > Tol())
      {
        near.emplace_back(lo, hi);
      }
    }
    return SubtractIntervals(allowed, MergeIntervals(near, Tol()), Tol());
  };
  ChainPoint best;
  best.distance = std::numeric_limits<double>::infinity();
  auto Evaluate = [&](std::size_t k)
  {
    const Run &r = runs[chain.runs[k]];
    for (const auto &interval : Allowed(k))
    {
      const double s =
          std::clamp(Dot(Sub(p, r.start), r.tangent), interval.first, interval.second);
      const double distance = Distance(p, r.At(s));
      if (distance < best.distance)
      {
        best = {chain.runs[k], s, chain.run_offset[k] + s, distance};
      }
    }
  };
  // Grid rings around p until no farther cell can hold a run as close as the best so far
  // (a cell at ring k lies at least (k - 1) cells away): the candidates then contain every
  // run of the chain at the minimal distance, and the minimum is taken in chain order (the
  // former loop over the whole chain, with its first-in-chain tie rule).
  std::vector<std::size_t> candidates, ring;
  const long long int max_ring = run_grid->MaxRing(p);
  for (long long int k = 0; k <= max_ring; k++)
  {
    if (k > 0 && static_cast<double>(k - 1) * run_grid->CellSize() > best.distance)
    {
      break;
    }
    if (max_distance && k > 0 &&
        static_cast<double>(k - 1) * run_grid->CellSize() > *max_distance + Tol())
    {
      break;  // nothing closer than max_distance beyond this ring
    }
    ring.clear();
    run_grid->Ring(p, k, ring);
    for (const std::size_t r : ring)
    {
      if (runs[r].chain == chain.id && !runs[r].excluded)
      {
        candidates.push_back(runs[r].index_in_chain);
        Evaluate(runs[r].index_in_chain);
      }
    }
  }
  best = ChainPoint{};
  best.distance = std::numeric_limits<double>::infinity();
  std::sort(candidates.begin(), candidates.end());
  candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
  for (const std::size_t k : candidates)
  {
    Evaluate(k);
  }
  if (max_distance && best.distance > *max_distance + Tol())
  {
    best = ChainPoint{};
    best.distance = std::numeric_limits<double>::infinity();
  }
  return best;
}

Identifier::ChainPoint Identifier::ArcAwareImage(const Chain &source, double x,
                                                 const Chain &other,
                                                 std::optional<double> exclude_x,
                                                 std::optional<double> max_distance) const
{
  std::size_t source_run = 0;
  double source_s = 0.0;
  ChainAt(source, x, &source_run, &source_s);
  const Point3D p = RunPoint(source_run, source_s);
  ChainPoint foot = ClosestPointOnChain(other, p, exclude_x, max_distance);
  if (!std::isfinite(foot.distance) || !run_arcs[foot.run])
  {
    return foot;
  }
  // The foot on the arc: the run of the same arc (in chain order from the chord foot's run,
  // both ways) whose angular range holds the radial projection of p; else the chord foot's
  // run clamped at the nearer end.
  const int arc = run_arcs[foot.run]->arc;
  const std::size_t m = other.runs.size();
  const std::size_t k0 = runs[foot.run].index_in_chain;
  auto OnRun = [&](std::size_t k, bool clamp) -> std::optional<ChainPoint>
  {
    const std::size_t r = other.runs[k];
    if (!run_arcs[r] || run_arcs[r]->arc != arc || runs[r].excluded)
    {
      return std::nullopt;
    }
    const ArcPiece &piece = run_arcs[r]->piece;
    double theta = piece.AngleOf(p);
    if (!piece.Contains(theta))
    {
      if (!clamp)
      {
        return std::nullopt;
      }
      theta = std::clamp(theta, std::min(piece.theta0, piece.theta1),
                         std::max(piece.theta0, piece.theta1));
    }
    const double sweep = piece.theta1 - piece.theta0;
    const double t = std::abs(sweep) > 0.0 ? (theta - piece.theta0) / sweep : 0.0;
    const double s = std::clamp(t, 0.0, 1.0) * runs[r].length;
    ChainPoint image{r, s, other.run_offset[k] + s, Distance(p, piece.At(theta))};
    if (exclude_x)
    {
      const double neighbourhood = kSelfPairNeighbourhoodOverRadius * R;
      if (ChainArcDistance(other, image.x, image.x, *exclude_x) < neighbourhood)
      {
        return std::nullopt;
      }
    }
    return image;
  };
  if (auto image = OnRun(k0, false))
  {
    return *image;
  }
  for (std::size_t step = 1; step < m; step++)
  {
    for (const int direction : {1, -1})
    {
      const long long int index =
          static_cast<long long int>(k0) + direction * static_cast<long long int>(step);
      if (!other.closed && (index < 0 || index >= static_cast<long long int>(m)))
      {
        continue;
      }
      const std::size_t k = static_cast<std::size_t>(
          (index % static_cast<long long int>(m) + static_cast<long long int>(m)) %
          static_cast<long long int>(m));
      if (!run_arcs[other.runs[k]] || run_arcs[other.runs[k]]->arc != arc)
      {
        continue;
      }
      if (auto image = OnRun(k, false))
      {
        return *image;
      }
    }
    // Both directions left the arc: stop at the first step where neither neighbour is on
    // it.
    bool any = false;
    for (const int direction : {1, -1})
    {
      const long long int index =
          static_cast<long long int>(k0) + direction * static_cast<long long int>(step);
      if (!other.closed && (index < 0 || index >= static_cast<long long int>(m)))
      {
        continue;
      }
      const std::size_t k = static_cast<std::size_t>(
          (index % static_cast<long long int>(m) + static_cast<long long int>(m)) %
          static_cast<long long int>(m));
      any = any || (run_arcs[other.runs[k]] && run_arcs[other.runs[k]]->arc == arc);
    }
    if (!any)
    {
      break;
    }
  }
  if (auto image = OnRun(k0, true))
  {
    return *image;
  }
  return foot;
}

// The turn at every joint of a chain (a sub-corner vertex between two runs) is spread over
// the two adjacent half-chords; the windowed curvature at x is the integral of that density
// over [x - W/2, x + W/2] (W = kCurvatureWindowOverRadius R, clipped at the ends of an open
// chain, periodic on a closed one) divided by W. Curved intervals are the sublevel set of
// the windowed bend radius below kStraightBendRadiusOverRadius R, with the crossings solved
// on the piecewise-linear function so that they do not depend on the mesh.
void Identifier::ComputeCurvature()
{
  const double W = kCurvatureWindowOverRadius * R;
  const double straight_radius = kStraightBendRadiusOverRadius * R;
  for (std::size_t c = 0; c < chains.size(); c++)
  {
    stage.Progress(c, chains.size());
    auto &chain = chains[c];
    const std::size_t m = chain.runs.size();
    chain.run_offset.assign(m, 0.0);
    chain.length = 0.0;
    for (std::size_t k = 0; k < m; k++)
    {
      chain.run_offset[k] = chain.length;
      chain.length += runs[chain.runs[k]].length;
    }
    chain.joint_turn.assign(m, 0.0);
    if (chain.joint_excluded.size() != m)
    {
      chain.joint_excluded.assign(m, false);
    }
    chain.kappa_nodes.clear();
    chain.signed_kappa_nodes.clear();
    chain.curved.clear();
    chain.curved_max_kappa.clear();
    if (m < 2 || chain.length <= 0.0)
    {
      continue;
    }
    // Side of the gap along the chain (constant: the chain bounds one metal face): +1 when
    // the gap lies to the left of the traversal (process normal x tangent), -1 to the
    // right. A turn toward the metal is a turn away from the gap.
    auto GapLeft = [&](const Run &run)
    { return Dot(run.gap_direction, Cross(run.process_normal, run.tangent)) > 0.0; };
    // Joints of a fitted arc (DetectArcs) take no part as joints: a rounded corner accounts
    // for them; a bend contributes its exact density 1 / radius over the arc below.
    std::set<int> chain_arcs;
    for (std::size_t k = 0; k < m; k++)
    {
      if (k == 0 && !chain.closed)
      {
        continue;
      }
      const Run &a = runs[chain.runs[(k + m - 1) % m]];
      const Run &b = runs[chain.runs[k]];
      if (a.excluded || b.excluded || a.end_vertex != b.start_vertex)
      {
        continue;
      }
      // The geometric turn of every joint (the pair rule reads the local joint turn of the
      // polyline whether or not an arc accounts for it); an arc joint contributes no
      // density of its own.
      chain.joint_turn[k] = std::acos(std::clamp(Dot(a.tangent, b.tangent), -1.0, 1.0));
      if (const int arc = vertex_arc[b.start_vertex]; arc >= 0)
      {
        chain.joint_excluded[k] = true;
        if (!arcs[static_cast<std::size_t>(arc)].corner)
        {
          chain_arcs.insert(arc);
        }
      }
    }
    // Density contributions: every joint's turn over its two half-chords, every bend arc's
    // 1 / radius over [x(Ta), x(Tb)]; summed on the elementary intervals between their
    // breakpoints (cumulative turn F(x) = int_0^x density).
    struct Piece
    {
      double x0, x1, density;
      double signed_density;  // toward the metal positive
    };
    std::vector<Piece> contributions;
    for (std::size_t k = 0; k < m; k++)
    {
      if (chain.joint_turn[k] <= 0.0 || chain.joint_excluded[k])
      {
        continue;
      }
      const std::size_t prev = (k + m - 1) % m;
      const Run &a = runs[chain.runs[prev]];
      const Run &b = runs[chain.runs[k]];
      const double half_prev = 0.5 * a.length;
      const double half = 0.5 * b.length;
      const double window = half_prev + half;
      if (window <= 0.0)
      {
        continue;
      }
      const double density = chain.joint_turn[k] / window;
      const bool turn_left = Dot(Cross(a.tangent, b.tangent), b.process_normal) > 0.0;
      const double toward_metal = turn_left != GapLeft(b) ? 1.0 : -1.0;
      const double x = chain.run_offset[k];
      if (x - half_prev < 0.0)  // the joint before run 0 of a closed chain wraps
      {
        contributions.push_back(
            {x - half_prev + chain.length, chain.length, density, toward_metal * density});
        contributions.push_back({0.0, x + half, density, toward_metal * density});
      }
      else
      {
        contributions.push_back({x - half_prev, x + half, density, toward_metal * density});
      }
    }
    // Bend arcs touching this chain (through a joint or an arc segment of one of its runs):
    // the chain's part of every such arc carries the density 1 / radius.
    chain.arc_spans.clear();
    for (std::size_t k = 0; k < m; k++)
    {
      for (const auto &rs : runs[chain.runs[k]].segments)
      {
        if (const int a = segment_arc[rs.segment]; a >= 0)
        {
          chain_arcs.insert(a);
        }
      }
    }
    for (const int a : chain_arcs)
    {
      const Arc &arc = arcs[static_cast<std::size_t>(a)];
      std::vector<Interval> spans;
      std::optional<double> toward_metal;
      for (const std::size_t s : arc.segments)
      {
        const long long int r = run_of_segment[s];
        if (r < 0 || runs[static_cast<std::size_t>(r)].chain != chain.id)
        {
          continue;
        }
        const Run &run = runs[static_cast<std::size_t>(r)];
        const std::size_t k = run.index_in_chain;
        for (const auto &rs : run.segments)
        {
          if (rs.segment == s)
          {
            spans.emplace_back(chain.run_offset[k] + rs.t0, chain.run_offset[k] + rs.t1);
          }
        }
        if (!toward_metal)
        {
          // The arc turns toward its centre: toward the metal when the centre and the gap
          // lie on opposite sides of the run.
          const bool centre_left = Dot(Sub(arc.center, run.At(0.5 * run.length)),
                                       Cross(run.process_normal, run.tangent)) > 0.0;
          toward_metal = centre_left != GapLeft(run) ? 1.0 : -1.0;
        }
      }
      for (const auto &span : MergeIntervals(std::move(spans), Tol()))
      {
        contributions.push_back(
            {span.first, span.second, 1.0 / arc.radius, *toward_metal / arc.radius});
        chain.arc_spans.emplace_back(a, span.first, span.second);
      }
    }
    // Sweep over the contribution ends (+density at x0, -density at x1); the signed
    // density is carried alongside (the same breakpoints).
    struct Event
    {
      double x, density, signed_density;
      bool operator<(const Event &other) const
      {
        return std::tie(x, density, signed_density) <
               std::tie(other.x, other.density, other.signed_density);
      }
    };
    std::vector<Event> events = {{0.0, 0.0, 0.0}, {chain.length, 0.0, 0.0}};
    for (const auto &c : contributions)
    {
      events.push_back({std::clamp(c.x0, 0.0, chain.length), c.density, c.signed_density});
      events.push_back(
          {std::clamp(c.x1, 0.0, chain.length), -c.density, -c.signed_density});
    }
    std::sort(events.begin(), events.end());
    std::vector<Piece> pieces;
    double density = 0.0, signed_density = 0.0;
    for (std::size_t i = 0; i < events.size(); i++)
    {
      density += events[i].density;
      signed_density += events[i].signed_density;
      if (i + 1 < events.size() && events[i + 1].x > events[i].x)
      {
        pieces.push_back({events[i].x, events[i + 1].x, std::max(0.0, density),
                          std::abs(signed_density) <= 1.0e-12 * std::abs(density)
                              ? 0.0
                              : signed_density});
      }
    }
    std::vector<double> cumulative(pieces.size() + 1, 0.0),
        signed_cumulative(pieces.size() + 1, 0.0);
    for (std::size_t i = 0; i < pieces.size(); i++)
    {
      cumulative[i + 1] = cumulative[i] + pieces[i].density * (pieces[i].x1 - pieces[i].x0);
      signed_cumulative[i + 1] =
          signed_cumulative[i] + pieces[i].signed_density * (pieces[i].x1 - pieces[i].x0);
    }
    const double total_turn = cumulative.back();
    const double total_signed_turn = signed_cumulative.back();
    auto Cumulative = [&](double x, bool with_sign)
    {
      const auto &sums = with_sign ? signed_cumulative : cumulative;
      double shift = 0.0;
      if (chain.closed)
      {
        const double periods = std::floor(x / chain.length);
        x -= periods * chain.length;
        shift = periods * (with_sign ? total_signed_turn : total_turn);
      }
      else
      {
        x = std::clamp(x, 0.0, chain.length);
      }
      // The last piece starting before x (the pieces are consecutive, x0 non-decreasing).
      const auto after = std::lower_bound(pieces.begin(), pieces.end(), x,
                                          [](const Piece &piece, double value)
                                          { return piece.x0 < value; });
      if (after == pieces.begin())
      {
        return shift;
      }
      const auto i = static_cast<std::size_t>(after - pieces.begin()) - 1;
      const double rho = with_sign ? pieces[i].signed_density : pieces[i].density;
      return sums[i] + rho * (std::min(x, pieces[i].x1) - pieces[i].x0) + shift;
    };
    auto F = [&](double x) { return Cumulative(x, false); };
    auto Kappa = [&](double x) { return (F(x + 0.5 * W) - F(x - 0.5 * W)) / W; };
    auto SignedKappa = [&](double x)
    { return (Cumulative(x + 0.5 * W, true) - Cumulative(x - 0.5 * W, true)) / W; };
    // Nodes: the ends, every density breakpoint shifted by +-W/2 (periodic images for a
    // closed chain).
    std::vector<double> xs = {0.0, chain.length};
    for (const auto &piece : pieces)
    {
      for (const double b : {piece.x0, piece.x1})
      {
        for (const double shift : {-0.5 * W, 0.5 * W})
        {
          double x = b + shift;
          if (chain.closed)
          {
            x -= std::floor(x / chain.length) * chain.length;
          }
          if (x > 0.0 && x < chain.length)
          {
            xs.push_back(x);
          }
        }
      }
    }
    std::sort(xs.begin(), xs.end());
    xs.erase(std::unique(xs.begin(), xs.end(), [&](double a, double b)
                         { return std::abs(a - b) <= 1.0e-12 * R; }),
             xs.end());
    for (const double x : xs)
    {
      chain.kappa_nodes.emplace_back(x, Kappa(x));
      chain.signed_kappa_nodes.push_back(SignedKappa(x));
    }
    // Curved: windowed bend radius 1 / kappa below the straight threshold (quantized).
    auto Curved = [&](double kappa)
    { return kappa > 0.0 && quantizer.Less(1.0 / kappa, straight_radius); };
    const double level = 1.0 / straight_radius;
    std::optional<double> open_start;
    auto Close = [&](double x)
    {
      if (open_start && x - *open_start > Tol())
      {
        chain.curved.emplace_back(*open_start, x);
        // A section on a fitted bend arc reads the arc's exact curvature (the windowed
        // value is at most 1 / radius, less where the window is clipped or the arc shorter
        // than the window), so that the class value does not depend on the chords.
        double kappa_max = MaxCurvature(chain, *open_start, x);
        for (const auto &[a, x0, x1] : chain.arc_spans)
        {
          if (std::min(x, x1) - std::max(*open_start, x0) > Tol())
          {
            kappa_max = std::max(kappa_max, 1.0 / arcs[static_cast<std::size_t>(a)].radius);
          }
        }
        chain.curved_max_kappa.push_back(kappa_max);
      }
      open_start.reset();
    };
    for (std::size_t i = 0; i < chain.kappa_nodes.size(); i++)
    {
      const auto [x, kappa] = chain.kappa_nodes[i];
      const bool curved = Curved(kappa);
      if (i == 0)
      {
        if (curved)
        {
          open_start = x;
        }
        continue;
      }
      const auto [xp, kp] = chain.kappa_nodes[i - 1];
      const bool previous = Curved(kp);
      if (curved == previous)
      {
        continue;
      }
      // Linear crossing of the level between the two nodes.
      const double denominator = kappa - kp;
      const double xc =
          std::abs(denominator) > 0.0 ? xp + (level - kp) / denominator * (x - xp) : x;
      const double crossing = std::clamp(xc, xp, x);
      if (curved)
      {
        open_start = crossing;
      }
      else
      {
        Close(crossing);
      }
    }
    Close(chain.length);
  }
}

// Pairs along bends: two chains that are not both single straight runs and whose
// closest-point separation over the mutually paired intervals (within 2R, not beyond either
// chain's ends, outside the through-vertex zones of shared vertices) is constant within
// kPairSeparationTolerance form one pair feature per curvature class (straight-like /
// curved) claiming both chains; the curved class carries the tightest windowed bend radius.
void Identifier::BuildBentPairs()
{
  const double interaction = kInteractionDistanceOverRadius * R;
  // Candidate facing region: within 2R (1 + tolerance), so that a pair at the threshold
  // is sampled up to its vertices (see the rule at kPairSeparationTolerance).
  const double reach = interaction * (1.0 + kPairSeparationTolerance);
  const double zone = kThroughVertexZoneOverRadius * R;
  struct Bounds
  {
    Point3D lo{}, hi{};
  };
  std::vector<Bounds> bounds(chains.size());
  std::vector<bool> usable(chains.size(), false);
  for (std::size_t c = 0; c < chains.size(); c++)
  {
    bool first = true;
    for (const std::size_t r : chains[c].runs)
    {
      if (runs[r].excluded)
      {
        continue;
      }
      for (const Point3D *p : {&runs[r].start, &runs[r].end})
      {
        for (int d = 0; d < 3; d++)
        {
          bounds[c].lo[d] = first ? (*p)[d] : std::min(bounds[c].lo[d], (*p)[d]);
          bounds[c].hi[d] = first ? (*p)[d] : std::max(bounds[c].hi[d], (*p)[d]);
        }
        first = false;
        usable[c] = true;
      }
    }
  }
  // A sample of chain A's facing region: run, run parameter, chain position, closest point
  // on B (chain position and run), distance.
  struct Sample
  {
    std::size_t run;
    double s, x, d, qx;
    // The partner run the foot lies on: a non-exact stretch is cut where it changes (one
    // non-exact sub-piece faces one partner run, so both sides are discretised alike).
    std::size_t foot_run;
    // Window half-widths for the curve separation: a full chord (capped at
    // kPairChordWindowCapOverRadius R) where the chain bends within R of the sample / the
    // foot (a joint or an arc along the chain: the chords of an inscribed polyline dip
    // mid-chord), R on a straight run (no dip; a taper must be read locally).
    double half_own, half_other;
    // The larger joint turn (radians) at the ends of the sample's run and of the foot's
    // run: the chord reading C and the inscribed-vertex reading C / cos(turn / 2) of the
    // pair separation differ by this discretisation ambiguity.
    double turn;
    // Exact geometric separation of the sample (decision 85(1)): the perpendicular distance
    // of two exactly parallel straight runs, or the radius difference of two concentric
    // fitted bend arcs; absent where only the chord reading exists (a polyline bend with
    // chords beyond 2R, a taper, the transition between a lead and an arc).
    std::optional<double> exact;
    // Either chain bends within kPairBendProximityOverRadius R of the sample / its foot
    // (BendsWithinR): the sample is read off the chords and, without an exact reading, is
    // the junction's geometry, not its run's — the sub-piece's exactness is not judged on
    // it (the end sample of a lead at the joint reads the partner's first chord).
    bool near_bend;
  };
  // The larger joint turn at the two ends of run k of a chain (0 at an open chain's ends).
  auto LocalTurn = [&](const Chain &chain, std::size_t k)
  {
    const std::size_t m = chain.runs.size();
    double turn = chain.joint_turn.empty() ? 0.0 : chain.joint_turn[k];
    if (k + 1 < m)
    {
      turn = std::max(turn, chain.joint_turn[k + 1]);
    }
    else if (chain.closed && m > 0)
    {
      turn = std::max(turn, chain.joint_turn[0]);
    }
    return turn;
  };
  // The chain bends within R of chain position x (kPairBendProximityOverRadius): a joint
  // with a turn (not excluded; the wrap joint of a closed chain included) or a fitted bend
  // arc span within that distance along the chain. Binary searches on the run offsets and
  // on the chain's arc spans sorted along the chain (built once per chain here: Chain::
  // arc_spans is in arc order).
  const double proximity = kPairBendProximityOverRadius * R;
  std::vector<std::vector<Interval>> sorted_arc_spans(chains.size());
  for (std::size_t c = 0; c < chains.size(); c++)
  {
    for (const auto &[a, x0, x1] : chains[c].arc_spans)
    {
      (void)a;
      sorted_arc_spans[c].emplace_back(x0, x1);
    }
    std::sort(sorted_arc_spans[c].begin(), sorted_arc_spans[c].end());
  }
  auto BendsWithinR = [&](const Chain &chain, double x)
  {
    const std::size_t m = chain.runs.size();
    if (chain.closed && chain.length > 0.0)
    {
      x -= std::floor(x / chain.length) * chain.length;
    }
    auto JointNear = [&](std::size_t k, double joint_x)
    {
      if (k >= m || chain.joint_turn.empty() || chain.joint_turn[k] <= 0.0 ||
          chain.joint_excluded[k])
      {
        return false;
      }
      double dx = std::abs(x - joint_x);
      if (chain.closed)
      {
        dx = std::min(dx, chain.length - dx);
      }
      return !quantizer.Less(proximity, dx);
    };
    // The joints (before run k, at run_offset[k]) from x - R on, until beyond x + R.
    const auto first = std::lower_bound(chain.run_offset.begin(), chain.run_offset.end(),
                                        x - proximity - Tol());
    for (auto it = first; it != chain.run_offset.end() && *it <= x + proximity + Tol();
         ++it)
    {
      if (JointNear(static_cast<std::size_t>(it - chain.run_offset.begin()), *it))
      {
        return true;
      }
    }
    if (chain.closed && m > 0 && JointNear(0, chain.length))
    {
      return true;  // the wrap joint seen from the chain's end
    }
    // Arc spans [x0, x1] (disjoint, ascending): the last span starting at or before x + R
    // is the only one that can reach x - R; on a closed chain the spans at either end wrap.
    const auto &spans = sorted_arc_spans[chain_index.at(chain.id)];
    const auto after =
        std::upper_bound(spans.begin(), spans.end(), x + proximity,
                         [](double value, const Interval &s) { return value < s.first; });
    if (after != spans.begin() && !quantizer.Less(std::prev(after)->second, x - proximity))
    {
      return true;
    }
    if (chain.closed && !spans.empty() &&
        (!quantizer.Less(spans.back().second - chain.length, x - proximity) ||
         !quantizer.Less(x + proximity - chain.length, spans.front().first)))
    {
      return true;
    }
    return false;
  };
  // Window half-width of the chord reading where the chain bends: a full chord, at least R
  // and at most kPairChordWindowCapOverRadius R.
  auto ChordWindow = [&](double chord)
  { return std::max(R, std::min(chord, kPairChordWindowCapOverRadius * R)); };
  struct Piece
  {
    std::size_t run;
    Interval interval;
    bool curved;
    double max_kappa;
    std::vector<Sample> samples;  // ordered along the run
  };
  // Exact separation of a sample on run a whose foot lies on run b (decision 85(1)): two
  // concentric fitted bend arcs (the radius difference; the arc rule's centres and radii
  // are exact for polylines inscribed in the design circles at any chord count) or two
  // exactly parallel straight runs where neither chain bends within R of the sample / the
  // foot (the leads of a route: the perpendicular distance of their lines; a chord of a
  // coarse polyline bend has a joint within R and is read off the chords like the rest).
  // The straight reading is the sample's own separation only where its foot is the
  // perpendicular projection onto run b (the closest-point distance equals the line
  // distance within the signature tolerance); a foot clamped at an end of run b (the
  // sample lies past the parallel run, opposite a diverging piece) is no exact reading
  // (USER decision 184 (3): a piece's separation reflects its own geometry). "Bends within
  // R" is BendsWithinR (a joint or arc along the chain), not the windowed curvature.
  auto ExactSeparation = [&](std::size_t a, std::size_t b, bool bends_a, bool bends_b,
                             double chord_distance) -> std::optional<double>
  {
    if (a == b)
    {
      return std::nullopt;
    }
    const Run &ra = runs[a];
    const Run &rb = runs[b];
    const int arc_a = RunArc(a), arc_b = RunArc(b);
    if (arc_a >= 0 && arc_b >= 0)
    {
      if (arc_a == arc_b)
      {
        return std::nullopt;  // a chain facing itself along one arc
      }
      // Concentric within the signature parameter tolerance: the radius difference is then
      // exact to that tolerance (the centre offset bounds its error); two local circles of
      // a spline (a spiral fitted piecewise) are not concentric and read off the chords.
      const Arc &A = arcs[static_cast<std::size_t>(arc_a)];
      const Arc &B = arcs[static_cast<std::size_t>(arc_b)];
      if (Distance(A.center, B.center) <= kSignatureParameterToleranceOverRadius * R)
      {
        return std::abs(A.radius - B.radius);
      }
      return std::nullopt;
    }
    if (arc_a < 0 && arc_b < 0 && !bends_a && !bends_b &&
        !DirectionLess(std::abs(Dot(ra.tangent, rb.tangent)),
                       1.0 - kParallelCosineTolerance))
    {
      const Point3D delta = Sub(rb.start, ra.start);
      const double line_distance =
          Norm(Sub(delta, Scale(Dot(delta, ra.tangent), ra.tangent)));
      if (std::abs(chord_distance - line_distance) <
          kSignatureParameterToleranceOverRadius * R)
      {
        return line_distance;
      }
    }
    return std::nullopt;
  };
  // Candidate pieces of chain A's runs with respect to chain B (within the reach, not
  // beyond either end of B, outside the shared-vertex zones, cut at the curved boundaries
  // of both chains) with their samples: at least kPairSeparationSamplesPerInterval + 1 per
  // piece and at most R / 2 apart, so that a 2R window always holds several samples.
  const double margin = reach + 2.0 * Tol();
  // The points of every chain's curved boundaries (chain position x mapped onto the run
  // containing it), indexed by their positions: a run's cuts at the partner's curved
  // boundaries come from the points within its reach instead of every boundary of the
  // partner.
  struct CurvedBoundary
  {
    int chain;
    Point3D point;
  };
  std::vector<CurvedBoundary> curved_boundaries;
  std::optional<UniformGrid> curved_boundary_grid;
  {
    Point3D lo{}, hi{};
    for (const auto &B : chains)
    {
      for (const auto &ci : B.curved)
      {
        for (const double x : {ci.first, ci.second})
        {
          const std::size_t kb =
              std::min(static_cast<std::size_t>(
                           std::upper_bound(B.run_offset.begin(), B.run_offset.end(), x) -
                           B.run_offset.begin()) -
                           1,
                       B.runs.size() - 1);
          const std::size_t rb = B.runs[kb];
          const Point3D q =
              runs[rb].At(std::clamp(x - B.run_offset[kb], 0.0, runs[rb].length));
          for (int d = 0; d < 3; d++)
          {
            lo[d] = curved_boundaries.empty() ? q[d] : std::min(lo[d], q[d]);
            hi[d] = curved_boundaries.empty() ? q[d] : std::max(hi[d], q[d]);
          }
          curved_boundaries.push_back({B.id, q});
        }
      }
    }
    if (!curved_boundaries.empty())
    {
      curved_boundary_grid.emplace(4.0 * R, lo);
      for (std::size_t id = 0; id < curved_boundaries.size(); id++)
      {
        curved_boundary_grid->Insert(id, curved_boundaries[id].point,
                                     curved_boundaries[id].point);
      }
    }
  }
  // Runs of chain B (in chain order) whose boxes lie within the reach of the given box.
  auto RunsOfChainNear = [&](const Chain &B, const Bounds &box)
  {
    std::vector<std::size_t> indices;
    for (const std::size_t r : RunsNear(box.lo, box.hi, margin))
    {
      if (runs[r].chain == B.id && !runs[r].excluded)
      {
        indices.push_back(runs[r].index_in_chain);
      }
    }
    std::sort(indices.begin(), indices.end());
    return indices;
  };
  // Runs of chain A (in chain order) whose boxes lie within the reach of a run of chain B:
  // the only runs of A whose facing region can be non-empty. Found from the smaller chain
  // (a chain's box is not its extent: the ground chain around a route spans the chip).
  auto RunsOfChainFacing = [&](const Chain &A, const Chain &B)
  {
    std::vector<std::size_t> indices;
    if (A.runs.size() <= B.runs.size())
    {
      for (const std::size_t a : A.runs)
      {
        if (runs[a].excluded)
        {
          continue;
        }
        const auto near = RunsNearRun(a, margin);
        if (std::any_of(near.begin(), near.end(), [&](std::size_t r)
                        { return runs[r].chain == B.id && !runs[r].excluded; }))
        {
          indices.push_back(runs[a].index_in_chain);
        }
      }
    }
    else
    {
      for (const std::size_t b : B.runs)
      {
        if (runs[b].excluded)
        {
          continue;
        }
        for (const std::size_t r : RunsNearRun(b, margin))
        {
          if (runs[r].chain == A.id && !runs[r].excluded)
          {
            indices.push_back(runs[r].index_in_chain);
          }
        }
      }
      std::sort(indices.begin(), indices.end());
      indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
    }
    return indices;
  };
  auto Pieces = [&](const Chain &A, const Chain &B, std::vector<Piece> &pieces)
  {
    // A chain facing itself across a fold (A == B): partner runs at least pi R of arc
    // length away along the chain (kSelfPairNeighbourhoodOverRadius), the foot taken
    // outside that neighbourhood, no shared-vertex zones (every vertex is shared).
    const bool self = A.id == B.id;
    const double neighbourhood = kSelfPairNeighbourhoodOverRadius * R;
    // A chain facing itself: a run pair whose points are all within the neighbourhood of
    // each other along the chain (the largest arc distance between the two runs below pi R)
    // never faces; otherwise every point of run ka is a candidate and its partner is the
    // closest point of the chain OUTSIDE the neighbourhood of that point (the per-sample
    // foot below), so that the pair starts exactly where the facing point is pi R of arc
    // away.
    auto SelfAllowed = [&](std::size_t ka, std::size_t kb)
    {
      std::vector<Interval> allowed;
      const double a0 = A.run_offset[ka], a1 = a0 + runs[A.runs[ka]].length;
      const double b0 = A.run_offset[kb], b1 = b0 + runs[A.runs[kb]].length;
      if (!self)
      {
        allowed.emplace_back(0.0, a1 - a0);
        return allowed;
      }
      if (ka == kb)
      {
        return allowed;
      }
      double farthest = std::max(std::abs(b1 - a0), std::abs(a1 - b0));
      if (A.closed)
      {
        // The shorter way round for the two far ends.
        farthest = std::max(std::min(std::abs(b1 - a0), A.length - std::abs(b1 - a0)),
                            std::min(std::abs(a1 - b0), A.length - std::abs(a1 - b0)));
      }
      if (quantizer.Less(farthest, neighbourhood))
      {
        return allowed;
      }
      allowed.emplace_back(0.0, a1 - a0);
      return allowed;
    };
    for (const std::size_t ka : RunsOfChainFacing(A, B))
    {
      const std::size_t a = A.runs[ka];
      const Run &ra = runs[a];
      if (ra.excluded)
      {
        continue;
      }
      std::vector<Interval> within;
      Point3D a_lo, a_hi;
      BoundingBox(ra.start, ra.end, a_lo, a_hi);
      for (const std::size_t kb : RunsOfChainNear(B, Bounds{a_lo, a_hi}))
      {
        const std::size_t b = B.runs[kb];
        const Run &rb = runs[b];
        if (rb.excluded ||
            !quantizer.Less(SegmentSegmentDistance(ra.start, ra.end, rb.start, rb.end),
                            reach))
        {
          continue;
        }
        auto found = RunIntervalWithin(a, rb.start, rb.end, reach);
        if (self)
        {
          found = IntersectIntervals(found, SelfAllowed(ka, kb), Tol());
        }
        within.insert(within.end(), found.begin(), found.end());
      }
      within = MergeIntervals(within, Tol());
      if (within.empty())
      {
        continue;
      }
      std::vector<Interval> excluded;
      if (!B.closed)
      {
        // Beyond either end of B: past the end along its outward tangent AND within the
        // reach of that end (the half-plane alone is right for a straight B, but past the
        // end plane of a bent route it covers points that face B's interior: DS-SCT-001's
        // ground-side edges lost their whole facing region to the trace's end planes and
        // the pairs were described on one side only).
        for (const bool at_end : {false, true})
        {
          const Run &terminal = runs[B.runs[at_end ? B.runs.size() - 1 : 0]];
          const Point3D e = at_end ? terminal.end : terminal.start;
          const Point3D t_out = at_end ? terminal.tangent : Scale(-1.0, terminal.tangent);
          const double g0 = Dot(Sub(ra.start, e), t_out), slope = Dot(ra.tangent, t_out);
          // g(s) = g0 + slope s > 0.
          std::vector<Interval> past;
          if (std::abs(slope) <= kDirectionQuantum)
          {
            if (g0 > 0.0)
            {
              past.emplace_back(0.0, ra.length);
            }
          }
          else
          {
            const double root = -g0 / slope;
            if (slope > 0.0)
            {
              past.emplace_back(std::max(root, 0.0), ra.length);
            }
            else
            {
              past.emplace_back(0.0, std::min(root, ra.length));
            }
          }
          const auto near_end =
              IntersectIntervals(past, RunIntervalWithin(a, e, e, reach), Tol());
          excluded.insert(excluded.end(), near_end.begin(), near_end.end());
        }
      }
      if (!self)
      {
        for (const auto &zone_intervals : ThroughZones(a, A, B))
        {
          excluded.insert(excluded.end(), zone_intervals.begin(), zone_intervals.end());
        }
      }
      const auto paired = SubtractIntervals(within, MergeIntervals(excluded, Tol()), Tol());
      if (std::getenv("PALACE_IDENTIFICATION_DEBUG_PIECES") && input.log)
      {
        std::ostringstream dbg;
        dbg << std::setprecision(10) << "  DEBUG facing run " << a << " of chain " << A.id
            << " vs chain " << B.id << " start (" << ra.start[0] << ", " << ra.start[1]
            << ") length " << ra.length << " within";
        for (const auto &w : within)
        {
          dbg << " [" << w.first << ", " << w.second << "]";
        }
        dbg << " excluded";
        for (const auto &w : excluded)
        {
          dbg << " [" << w.first << ", " << w.second << "]";
        }
        input.log(dbg.str() + "\n");
      }
      for (const auto &interval : paired)
      {
        // Cuts at A's curved boundaries and at B's curved boundaries mapped onto the run.
        // (The cuts are sorted below and zero-length pieces skipped: the candidate order
        // and duplicates do not matter, only the set of cut positions.)
        std::vector<double> cuts = {interval.first, interval.second};
        {
          // A's curved boundaries (disjoint, ascending) near the interval (a superset by
          // R on either side; the exact test below is the former one).
          const double x_lo = A.run_offset[ka] + interval.first - R;
          const double x_hi = A.run_offset[ka] + interval.second + R;
          auto ci = std::lower_bound(A.curved.begin(), A.curved.end(), x_lo,
                                     [](const Interval &c, double value)
                                     { return c.second <= value; });
          for (; ci != A.curved.end() && ci->first < x_hi; ++ci)
          {
            for (const double x : {ci->first, ci->second})
            {
              const double s = x - A.run_offset[ka];
              if (s > interval.first + Tol() && s < interval.second - Tol())
              {
                cuts.push_back(s);
              }
            }
          }
        }
        if (curved_boundary_grid)
        {
          // B's curved boundaries within the reach of run a (the closest point on the run
          // is the projection): candidates from the grid of every chain's boundary points.
          Point3D a_lo, a_hi;
          BoundingBox(ra.start, ra.end, a_lo, a_hi);
          for (const std::size_t id : curved_boundary_grid->Query(a_lo, a_hi, margin))
          {
            if (curved_boundaries[id].chain != B.id)
            {
              continue;
            }
            const Point3D &q = curved_boundaries[id].point;
            const double s = std::clamp(Dot(Sub(q, ra.start), ra.tangent), 0.0, ra.length);
            if (quantizer.Less(Distance(ra.At(s), q), reach) &&
                s > interval.first + Tol() && s < interval.second - Tol())
            {
              cuts.push_back(s);
            }
          }
        }
        std::sort(cuts.begin(), cuts.end());
        for (std::size_t i = 0; i + 1 < cuts.size(); i++)
        {
          const Interval piece{cuts[i], cuts[i + 1]};
          if (piece.second - piece.first <= Tol())
          {
            continue;
          }
          const double x0 = A.run_offset[ka] + piece.first;
          const double x1 = A.run_offset[ka] + piece.second;
          const double mid = 0.5 * (x0 + x1);
          const auto q_mid =
              ClosestPointOnChain(B, ra.At(0.5 * (piece.first + piece.second)),
                                  self ? std::optional<double>(mid) : std::nullopt,
                                  self ? std::optional<double>(reach) : std::nullopt);
          if (!std::isfinite(q_mid.distance))
          {
            continue;  // self: no partner outside the neighbourhood
          }
          const bool curved = IsCurvedAt(A, mid) || IsCurvedAt(B, q_mid.x);
          const double max_kappa =
              std::max(MaxCurvature(A, x0, x1), WindowedCurvature(B, q_mid.x));
          const int n_samples =
              std::max(kPairSeparationSamplesPerInterval,
                       static_cast<int>(std::ceil(2.0 * (piece.second - piece.first) / R)));
          Piece result{a, piece, curved, max_kappa, {}};
          for (int i_s = 0; i_s <= n_samples; i_s++)
          {
            const double s = piece.first + (piece.second - piece.first) * i_s / n_samples;
            const double x = A.run_offset[ka] + s;
            const auto q = ClosestPointOnChain(
                B, ra.At(s), self ? std::optional<double>(x) : std::nullopt,
                self ? std::optional<double>(reach) : std::nullopt);
            if (!std::isfinite(q.distance))
            {
              continue;
            }
            const bool bends_own = BendsWithinR(A, x);
            const bool bends_other = BendsWithinR(B, q.x);
            const double half_own = bends_own ? ChordWindow(ra.length) : R;
            const double half_other = bends_other ? ChordWindow(runs[q.run].length) : R;
            const double turn =
                std::max(LocalTurn(A, ka), LocalTurn(B, RunIndexInChain(B, q.run)));
            result.samples.push_back(
                {a, s, x, q.distance, q.x, q.run, half_own, half_other, turn,
                 ExactSeparation(a, q.run, bends_own, bends_other, q.distance),
                 bends_own || bends_other});
          }
          pieces.push_back(std::move(result));
        }
      }
    }
  };
  // Local constancy and the curve separation of every sample (rule at
  // kPairSeparationTolerance): a sample is constant when the distances within R of it
  // along its own chain vary by at most the tolerance; its curve separation is the smaller
  // of the two directional maxima over the windows of half-width max(R, local chord) where
  // the chain bends (a window that always contains a vertex of an inscribed polyline and a
  // full chord of an offset polyline) and R on straight runs — on its own chain about x and
  // on the other chain about the foot qx.
  auto Classify = [&](std::vector<Piece> &own, const std::vector<Piece> &other,
                      std::vector<char> &constant, std::vector<double> &separation,
                      std::vector<double> &upper)
  {
    std::vector<const Sample *> own_samples, other_samples;
    for (const auto &piece : own)
    {
      for (const auto &sample : piece.samples)
      {
        own_samples.push_back(&sample);
      }
    }
    for (const auto &piece : other)
    {
      for (const auto &sample : piece.samples)
      {
        other_samples.push_back(&sample);
      }
    }
    auto ByX = [](const Sample *u, const Sample *v) { return u->x < v->x; };
    std::sort(own_samples.begin(), own_samples.end(), ByX);
    std::sort(other_samples.begin(), other_samples.end(), ByX);
    // Window maximum / minimum of the sampled distances within half of x along one chain;
    // with a constancy mask only the locally constant samples count (the chord reading of
    // a constant sample must not reach into a diverging neighbour: the last chord before a
    // taper kink, whose window of a full chord crosses the kink, read 2.15 - 2.7 um for a
    // 2 um gap on DS-SCT-001 and formed a separation group, a sliver pair and a hole in the
    // stack).
    auto WindowMax = [&](const std::vector<const Sample *> &list,
                         const std::vector<char> *mask, double x, double half,
                         double *min_out)
    {
      Sample probe{};
      probe.x = x - half;
      auto lo = std::lower_bound(list.begin(), list.end(), &probe, ByX);
      probe.x = x + half;
      auto hi = std::upper_bound(list.begin(), list.end(), &probe, ByX);
      double max_d = 0.0, min_d = std::numeric_limits<double>::infinity();
      for (auto it = lo; it != hi; ++it)
      {
        if (mask && !(*mask)[static_cast<std::size_t>(it - list.begin())])
        {
          continue;
        }
        max_d = std::max(max_d, (*it)->d);
        min_d = std::min(min_d, (*it)->d);
      }
      if (min_out)
      {
        *min_out = min_d;
      }
      return max_d;
    };
    // Local constancy of every sample of both chains (in the sorted order).
    auto Constancy = [&](const std::vector<const Sample *> &list)
    {
      std::vector<char> flags(list.size());
      for (std::size_t i = 0; i < list.size(); i++)
      {
        double min_d = 0.0;
        const double max_d = WindowMax(list, nullptr, list[i]->x, R, &min_d);
        flags[i] = !quantizer.Less(kPairSeparationTolerance * min_d, max_d - min_d);
      }
      return flags;
    };
    const std::vector<char> own_constant = Constancy(own_samples);
    const std::vector<char> other_constant = Constancy(other_samples);
    std::map<const Sample *, std::size_t> own_index;
    for (std::size_t i = 0; i < own_samples.size(); i++)
    {
      own_index.emplace(own_samples[i], i);
    }
    for (auto &piece : own)
    {
      for (const auto &sample : piece.samples)
      {
        const bool is_constant = own_constant[own_index.at(&sample)] != 0;
        constant.push_back(is_constant);
        const double w_own =
            WindowMax(own_samples, &own_constant, sample.x, sample.half_own, nullptr);
        const double w_other = other_samples.empty()
                                   ? w_own
                                   : WindowMax(other_samples, &other_constant, sample.qx,
                                               sample.half_other, nullptr);
        // Chord reading C (exact for an offset polyline) and the inscribed-vertex reading
        // C / cos(turn / 2) (exact for two polylines inscribed in the curves at aligned
        // angles): the pair interacts only when both readings are below 2R.
        const double c = w_other > 0.0 ? std::min(w_own, w_other) : w_own;
        separation.push_back(c);
        upper.push_back(c / std::cos(0.5 * sample.turn));
      }
    }
  };
  // Sub-pieces of one status along a piece: consecutive samples with the same (constant,
  // interacting) flags; the cut between two samples of different status is halfway. Such a
  // block is then cut into STRETCHES (USER decision 203, 2026-10-02, completing decision
  // 184 (3): a sub-piece's separation reflects its own geometry, and the exactness is
  // LOCAL along the run like the locality of the readings):
  //  - an EXACT stretch is a maximal contiguous stretch of samples that HAVE an exact
  //    reading (ExactSeparation with its perpendicular-foot guard) whose readings all lie
  //    within the signature parameter tolerance (1e-3 R) of the stretch's exact mean (a
  //    reading that would take the mean more than the tolerance from any member starts a
  //    new stretch: an extremely slow taper of chords parallel within the cosine tolerance
  //    keys a staircase of exact stretches, the parameter-tolerance semantics);
  //  - the samples WITHOUT an exact reading that are judged (away from any bend within
  //    kPairBendProximityOverRadius R) form NON-EXACT stretches: their status decides (a
  //    judged sample lacks an exact reading only where its partner is not parallel within
  //    the tolerance or its foot is not perpendicular), they are no longer compared with
  //    the exact mean;
  //  - the near-bend samples (not judged: the junction's geometry within R of a bend) never
  //    form, split or decide a stretch on their own: they join an adjacent stretch within
  //    kPairBendProximityOverRadius R along the run of that stretch's last exact / judged
  //    sample (between an exact and a non-exact stretch the first R goes to the exact one);
  //    a short gap of them between two exact stretches whose readings agree is transparent
  //    (one exact stretch across the partner's noise joint); contiguous near-bend samples
  //    farther than R from every exact / judged sample form a NON-EXACT stretch (the bend
  //    exemption is a 1 R zone, not a licence for an unbounded unjudged region: a straight
  //    run facing a finely chorded one-sided taper, every foot within R of a noise joint,
  //    read exact at the lead's value over the whole taper before);
  //  - an exact stretch splits off as its own EXACT sub-piece only when it is at least
  //    kExactStretchMinLengthOverRadius R long (the span of its exact samples); a shorter
  //    one merges into the adjacent non-exact stretch(es). A block that is ONE exact
  //    stretch stays exact whatever its length, as before (the sub-R exact stack ends are
  //    not re-keyed by this clause); consecutive non-exact stretches are one sub-piece.
  // The exact sub-piece's separation is the mean of its exact readings; a non-exact
  // sub-piece reads the length-weighted mean over its samples of the chord readings
  // (equally spaced samples; the exact reading where a sample has one), not exact. A
  // straight 450 um run facing a partner that is parallel over 150 um and then tapers
  // one-sidedly keys the facing part an exact strip and the rest a non-exact strip at its
  // mean; before (exactness all-or-nothing per block) one judged sample at the far end made
  // the whole run non-exact at its 450 um chord mean, the lead's exact group had no partner
  // piece and the 2.0 um strip (0.53 R) read as two IsolatedEdges. Bounded consequence: an
  // exact stretch absorbs up to R of a slowly changing partner (a mis-keying of at most
  // slope x R; 0.006 R in that reproducer). Considered alternative, not taken (it changes
  // the judging globally): judging the near-bend samples too at the tolerance widened by
  // their discretisation ambiguity d (1 - cos(turn)).
  struct SubPiece
  {
    std::size_t run;
    Interval interval;
    bool curved;
    double max_kappa;
    bool constant;
    bool interacting;
    double separation;
    bool exact;
    // The one partner run every sample's foot lies on (a non-exact sub-piece after the
    // partner-run cut; an exact one when its feet happen to lie on one run), else absent.
    std::optional<std::size_t> foot_run;
  };
  // Stretch statistics for the stage log: blocks cut into several stretches, exact
  // stretches, the fragmentation census asked for with the rule (a non-exact stretch
  // shorter than R between two exact stretches of the same mean) and the pieces moved by
  // the mutual-facing group invariant.
  std::size_t stretch_split_blocks = 0, exact_stretch_count = 0, sliver_count = 0,
              moved_count = 0;
  double sliver_length = 0.0, moved_length = 0.0;
  std::ostringstream sliver_sites, moved_sites;
  auto Split = [&](const std::vector<Piece> &pieces, const std::vector<char> &constant,
                   const std::vector<double> &separation, const std::vector<double> &upper,
                   std::vector<SubPiece> &out)
  {
    const double parameter_tolerance = kSignatureParameterToleranceOverRadius * R;
    const double min_exact_span = kExactStretchMinLengthOverRadius * R;
    enum class Kind
    {
      Exact,
      NonExact,
      NearBend
    };
    struct Stretch
    {
      Kind kind;
      std::size_t first, last;  // sample indices (inclusive) within the piece
      // Exact stretches: the exact readings' statistics and the span of the exact samples.
      double sum = 0.0, lo = 0.0, hi = 0.0;
      int count = 0;
      std::size_t first_exact = 0, last_exact = 0;
    };
    // The reading joins the stretch when every member stays within the tolerance of the
    // new mean.
    auto Admits = [&](const Stretch &st, double sum, int count, double lo, double hi)
    {
      const double mean = (st.sum + sum) / (st.count + count);
      return quantizer.Less(std::max(st.hi, hi) - mean, parameter_tolerance) &&
             quantizer.Less(mean - std::min(st.lo, lo), parameter_tolerance);
    };
    auto Absorb = [](Stretch &st, const Stretch &other)
    {
      if (other.count > 0)
      {
        if (st.count == 0)
        {
          st.lo = other.lo;
          st.hi = other.hi;
          st.first_exact = other.first_exact;
        }
        else
        {
          st.lo = std::min(st.lo, other.lo);
          st.hi = std::max(st.hi, other.hi);
        }
        st.sum += other.sum;
        st.count += other.count;
        st.last_exact = other.last_exact;
      }
      st.last = other.last;
    };
    std::size_t k = 0;
    for (const auto &piece : pieces)
    {
      const std::size_t n = piece.samples.size();
      std::size_t i = 0;
      while (i < n)
      {
        const bool c = constant[k + i];
        const bool inter = quantizer.Less(upper[k + i], interaction);
        std::size_t j = i;
        while (j < n && constant[k + j] == c &&
               quantizer.Less(upper[k + j], interaction) == inter)
        {
          j++;
        }
        // Raw stretches of the block: consecutive samples of one kind, the exact ones cut
        // by the tolerance rule.
        std::vector<Stretch> raw;
        for (std::size_t q = i; q < j; q++)
        {
          const auto &sample = piece.samples[q];
          const Kind kind = sample.exact       ? Kind::Exact
                            : sample.near_bend ? Kind::NearBend
                                               : Kind::NonExact;
          bool joins = !raw.empty() && raw.back().kind == kind;
          if (joins && kind == Kind::Exact)
          {
            joins = Admits(raw.back(), *sample.exact, 1, *sample.exact, *sample.exact);
          }
          if (!joins)
          {
            raw.push_back({kind, q, q});
          }
          Stretch &st = raw.back();
          st.last = q;
          if (kind == Kind::Exact)
          {
            Stretch one{kind, q, q, *sample.exact, *sample.exact, *sample.exact, 1, q, q};
            Absorb(st, one);
          }
        }
        // The near-bend samples within R along the run of a stretch's boundary sample,
        // measured between the samples' cells (one sample spacing, at most R / 2, added):
        // the exact samples end just before (joint - R), so the near-bend sample AT the
        // joint lies R + up to one spacing from the last exact sample.
        const double spacing =
            n > 1 ? (piece.interval.second - piece.interval.first) / (n - 1) : 0.0;
        auto WithinProximity = [&](std::size_t q, std::size_t anchor)
        {
          return !quantizer.Less(proximity + spacing,
                                 std::abs(piece.samples[q].s - piece.samples[anchor].s));
        };
        // A short gap of near-bend samples (each within R of an adjacent exact sample)
        // between two exact stretches whose readings agree is transparent.
        std::vector<Stretch> merged;
        for (std::size_t m = 0; m < raw.size(); m++)
        {
          if (raw[m].kind == Kind::NearBend && m + 1 < raw.size() && !merged.empty() &&
              merged.back().kind == Kind::Exact && raw[m + 1].kind == Kind::Exact)
          {
            bool transparent = true;
            for (std::size_t q = raw[m].first; transparent && q <= raw[m].last; q++)
            {
              transparent = WithinProximity(q, merged.back().last_exact) ||
                            WithinProximity(q, raw[m + 1].first_exact);
            }
            const Stretch &next = raw[m + 1];
            if (transparent &&
                Admits(merged.back(), next.sum, next.count, next.lo, next.hi))
            {
              Absorb(merged.back(), raw[m]);
              Absorb(merged.back(), next);
              m++;
              continue;
            }
          }
          merged.push_back(raw[m]);
        }
        // The remaining near-bend stretches: a prefix joins the stretch before it and a
        // suffix the one after it (within R of their boundary samples; the exact side
        // first, else the nearer); the samples farther than R from both are a non-exact
        // stretch of their own.
        std::vector<Stretch> resolved;
        for (std::size_t m = 0; m < merged.size(); m++)
        {
          if (merged[m].kind != Kind::NearBend)
          {
            resolved.push_back(merged[m]);
            continue;
          }
          Stretch *left = resolved.empty() ? nullptr : &resolved.back();
          const Stretch *right = m + 1 < merged.size() ? &merged[m + 1] : nullptr;
          const std::size_t left_anchor =
              left ? (left->kind == Kind::Exact ? left->last_exact : left->last) : 0;
          const std::size_t right_anchor =
              right ? (right->kind == Kind::Exact ? right->first_exact : right->first) : 0;
          std::size_t q_lo = merged[m].first, q_hi = merged[m].last + 1;  // the orphans
          auto ToLeft = [&](std::size_t q)
          {
            if (!left || !WithinProximity(q, left_anchor))
            {
              return false;
            }
            if (!right || !WithinProximity(q, right_anchor))
            {
              return true;
            }
            if ((left->kind == Kind::Exact) != (right->kind == Kind::Exact))
            {
              return left->kind == Kind::Exact;
            }
            return piece.samples[q].s - piece.samples[left_anchor].s <=
                   piece.samples[right_anchor].s - piece.samples[q].s;
          };
          while (q_lo < q_hi && ToLeft(q_lo))
          {
            left->last = q_lo;
            q_lo++;
          }
          std::size_t right_first = right ? right->first : 0;
          while (q_lo < q_hi && right && WithinProximity(q_hi - 1, right_anchor))
          {
            q_hi--;
            right_first = q_hi;
          }
          if (q_lo < q_hi)
          {
            resolved.push_back({Kind::NonExact, q_lo, q_hi - 1});
          }
          if (right && right_first < right->first)
          {
            merged[m + 1].first = right_first;
          }
        }
        // Exact stretches split off only when their exact samples span at least the
        // minimum; a block that is one exact stretch stays exact. Consecutive non-exact
        // stretches are one.
        std::vector<Stretch> final_stretches;
        for (const auto &st : resolved)
        {
          Stretch item = st;
          if (item.kind == Kind::Exact && resolved.size() > 1 &&
              quantizer.Less(piece.samples[item.last_exact].s -
                                 piece.samples[item.first_exact].s,
                             min_exact_span))
          {
            item.kind = Kind::NonExact;
          }
          if (!final_stretches.empty() && final_stretches.back().kind == Kind::NonExact &&
              item.kind == Kind::NonExact)
          {
            final_stretches.back().last = item.last;
            continue;
          }
          final_stretches.push_back(item);
        }
        if (final_stretches.size() > 1)
        {
          stretch_split_blocks++;
        }
        // Emission: an exact stretch is one sub-piece; a non-exact stretch is cut where
        // its samples' foot crosses a run boundary of the partner chain (one non-exact
        // sub-piece faces one partner run and reads that run's local chord mean, the value
        // the partner's own piece reads: the two sides are discretised alike before the
        // link grouping — one long non-exact stretch over a changing partner, mean 2.58 um
        // over a 2.0 -> 3.7 um taper, could not pair with the partner's per-run pieces,
        // whose first 73 um joined the exact 2.0 um group within the 5 % pair tolerance).
        struct Emitted
        {
          std::size_t stretch, first, last;
        };
        std::vector<Emitted> emitted;
        for (std::size_t m = 0; m < final_stretches.size(); m++)
        {
          const Stretch &st = final_stretches[m];
          std::size_t first = st.first;
          for (std::size_t q = st.first + 1; st.kind == Kind::NonExact && q <= st.last; q++)
          {
            if (piece.samples[q].foot_run != piece.samples[q - 1].foot_run)
            {
              emitted.push_back({m, first, q - 1});
              first = q;
            }
          }
          emitted.push_back({m, first, st.last});
        }
        for (const Emitted &item : emitted)
        {
          const std::size_t m = item.stretch;
          const Stretch &st = final_stretches[m];
          const std::size_t lo_q = item.first, hi_q = item.last + 1;
          const double s_lo =
              lo_q == 0 ? piece.interval.first
                        : 0.5 * (piece.samples[lo_q - 1].s + piece.samples[lo_q].s);
          const double s_hi =
              hi_q == n ? piece.interval.second
                        : 0.5 * (piece.samples[hi_q - 1].s + piece.samples[hi_q].s);
          if (s_hi - s_lo <= Tol())
          {
            continue;
          }
          const bool exact_piece = st.kind == Kind::Exact;
          double value = 0.0;
          if (exact_piece)
          {
            value = st.sum / st.count;
            exact_stretch_count++;  // one sub-piece per exact stretch
          }
          else
          {
            // Mean of the chord readings over the stretch's samples (equally spaced along
            // the run: the length-weighted mean), the exact reading where there is one.
            for (std::size_t q = lo_q; q < hi_q; q++)
            {
              const auto &exact = piece.samples[q].exact;
              value += exact ? *exact : separation[k + q];
            }
            value /= static_cast<double>(hi_q - lo_q);
            // Fragmentation census: a non-exact stretch shorter than R between two exact
            // stretches of the same mean (counted once per stretch, on its first
            // sub-piece, with the stretch's span).
            const double stretch_span =
                piece.samples[st.last].s - piece.samples[st.first].s + spacing;
            if (lo_q == st.first && m > 0 && m + 1 < final_stretches.size() &&
                quantizer.Less(stretch_span, min_exact_span))
            {
              const Stretch &before = final_stretches[m - 1],
                            &after = final_stretches[m + 1];
              if (before.kind == Kind::Exact && after.kind == Kind::Exact &&
                  quantizer.Less(
                      std::abs(before.sum / before.count - after.sum / after.count),
                      parameter_tolerance))
              {
                sliver_count++;
                sliver_length += stretch_span;
                if (sliver_count <= 50)
                {
                  const Point3D p = runs[piece.run].At(s_lo);
                  sliver_sites << std::setprecision(6) << " run " << piece.run << " s "
                               << s_lo << " at (" << p[0] << ", " << p[1] << ", " << p[2]
                               << ") length " << stretch_span;
                }
              }
            }
          }
          if (std::getenv("PALACE_IDENTIFICATION_DEBUG_PIECES") && input.log &&
              final_stretches.size() > 1)
          {
            int n_exact = 0, n_judged = 0, n_near_bend = 0;
            for (std::size_t q = lo_q; q < hi_q; q++)
            {
              (piece.samples[q].exact       ? n_exact
               : piece.samples[q].near_bend ? n_near_bend
                                            : n_judged)++;
            }
            std::ostringstream dbg;
            dbg << std::setprecision(10) << "  DEBUG stretch " << m + 1 << " / "
                << final_stretches.size() << " of run " << piece.run << " s [" << s_lo
                << ", " << s_hi << "] samples " << lo_q << ".." << hi_q - 1 << " (exact "
                << n_exact << ", judged " << n_judged << ", near-bend " << n_near_bend
                << ") " << (exact_piece ? "exact " : "non-exact ") << value << "\n";
            input.log(dbg.str());
          }
          std::optional<std::size_t> foot_run = piece.samples[lo_q].foot_run;
          for (std::size_t q = lo_q + 1; foot_run && q < hi_q; q++)
          {
            if (piece.samples[q].foot_run != *foot_run)
            {
              foot_run.reset();
            }
          }
          out.push_back({piece.run,
                         {s_lo, s_hi},
                         piece.curved,
                         piece.max_kappa,
                         c,
                         inter,
                         value,
                         exact_piece,
                         foot_run});
        }
        i = j;
      }
      k += n;
    }
  };
  std::size_t chain_pairs_examined = 0, chain_pairs_paired = 0, self_pairs = 0;
  std::vector<std::size_t> partners;
  for (std::size_t ca = 0; ca < chains.size(); ca++)
  {
    stage.Progress(ca, chains.size());
    // Candidate partners cb > ca (ascending: the former loop over every later chain): the
    // chains with a run within the reach of a run of this chain. A pair without such runs
    // has no facing pieces and contributed nothing. The chain itself (cb == ca) is a
    // candidate when it has joints: a chain folding back within 2R through a bend >= R
    // faces itself (kSelfPairNeighbourhoodOverRadius).
    partners.clear();
    for (const std::size_t a : chains[ca].runs)
    {
      if (runs[a].excluded)
      {
        continue;
      }
      for (const std::size_t r : RunsNearRun(a, margin))
      {
        const std::size_t cb = chain_index.at(runs[r].chain);
        if ((cb > ca || (cb == ca && !chains[ca].Rigid())) && !runs[r].excluded)
        {
          partners.push_back(cb);
        }
      }
    }
    std::sort(partners.begin(), partners.end());
    partners.erase(std::unique(partners.begin(), partners.end()), partners.end());
    for (const std::size_t cb : partners)
    {
      const Chain &A = chains[ca];
      const Chain &B = chains[cb];
      const bool self = ca == cb;
      if (!usable[ca] || !usable[cb] ||
          run_plane[A.runs.front()] != run_plane[B.runs.front()])
      {
        continue;  // a pair never spans two metal planes (decision 82(1))
      }
      if (ParallelRigidRuns(A.runs.front(), B.runs.front()))
      {
        continue;  // two parallel straight runs (one class): the translational rule
      }
      bool near = true;
      for (int d = 0; d < 3; d++)
      {
        near = near && bounds[ca].lo[d] <= bounds[cb].hi[d] + interaction &&
               bounds[cb].lo[d] <= bounds[ca].hi[d] + interaction;
      }
      if (!near)
      {
        continue;
      }
      chain_pairs_examined++;
      std::vector<Piece> pieces_a, pieces_b;
      Pieces(A, B, pieces_a);
      if (!self)
      {
        Pieces(B, A, pieces_b);
      }
      if (pieces_a.empty() || (!self && pieces_b.empty()))
      {
        continue;
      }
      chain_pairs_paired++;
      std::vector<char> constant_a, constant_b;
      std::vector<double> separation_a, separation_b, upper_a, upper_b;
      Classify(pieces_a, self ? pieces_a : pieces_b, constant_a, separation_a, upper_a);
      if (!self)
      {
        Classify(pieces_b, pieces_a, constant_b, separation_b, upper_b);
      }
      std::vector<SubPiece> sub_a, sub_b;
      Split(pieces_a, constant_a, separation_a, upper_a, sub_a);
      if (!self)
      {
        Split(pieces_b, constant_b, separation_b, upper_b, sub_b);
      }
      if (std::getenv("PALACE_IDENTIFICATION_DEBUG_PIECES") && input.log)
      {
        std::ostringstream dbg;
        dbg << "  DEBUG pieces of chains " << A.id << " - " << B.id << ":\n";
        for (const auto *list : {&sub_a, &sub_b})
        {
          for (const auto &piece : *list)
          {
            const Run &rr = runs[piece.run];
            dbg << "    side " << (list == &sub_a ? 0 : 1) << " run " << piece.run << " s ["
                << piece.interval.first << ", " << piece.interval.second << "] at ("
                << rr.At(piece.interval.first)[0] << ", " << rr.At(piece.interval.first)[1]
                << ") constant " << piece.constant << " interacting " << piece.interacting
                << " separation " << piece.separation << (piece.exact ? " exact" : " chord")
                << (piece.curved ? " curved" : " straight") << "\n";
          }
        }
        input.log(dbg.str());
      }
      // Locally constant portions that do not interact (separation at or beyond 2R): no
      // feature, but their cross-chord interactions are not events. Portions that are not
      // locally constant (divergence at tees and port ends, fast tapers, acute arms) keep
      // the event rule.
      bool any_interacting = false;
      for (const auto *list : {&sub_a, &sub_b})
      {
        const int other = list == &sub_a ? B.id : A.id;
        for (const auto &piece : *list)
        {
          if (!piece.constant)
          {
            continue;
          }
          bent_claims[piece.run].emplace_back(other, piece.interval);
          any_interacting = any_interacting || piece.interacting;
        }
      }
      if (!any_interacting)
      {
        continue;
      }
      // A rounded corner is a vertex: its arc never pairs, and a foreign chain concentric
      // with it (a constant facing piece) is no event either — like the arms of a sharp
      // corner nested at the same separations, whose arms face the foreign arms at exactly
      // the separation and never closer. (Tried and rejected, 2026-09-25: keeping the
      // arc's foreign facing pieces event-eligible made the DS-SCT-001 flux-loop ends one
      // 160 um cluster each: the event set of a chain facing an arc reaches sqrt(3) R past
      // the tangent points, the R-balls another R, and the cores merged with the trace
      // junction clusters 1.5 R away.) The member the corner turns away is recomposed out
      // of the stack by the taken rule (AssembleStack); the pair / stack side facing the
      // corner's arc or window is a recorded neighbour of a vertex feature.
      if (corner_arc_chains.count(A.id) > 0 || corner_arc_chains.count(B.id) > 0)
      {
        continue;
      }
      // The constant, interacting pieces of both sides grouped by separation: a closed
      // chain (a trace loop) can face the partner at two separations (its near side at the
      // gap, its far side across the strip), and one pair is one separation (the design's
      // "one separation per pair" is exact for a constant pair and the mean of a slow
      // taper). Grouping (decision 85(1) as amended by USER decision 184 (3),
      // 2026-10-01): the EXACT pieces agreeing within the signature parameter tolerance
      // (1e-3 R) form one exact group each (in ascending order of separation); a chord
      // piece joins the exact group whose separation is within the pair tolerance (5 %,
      // the same as the local constancy) of its own — the nearest when several — as a
      // chord reading of that design separation; the remaining chord pieces form chord
      // groups where consecutive separations differ by at most the pair tolerance. Before
      // the amendment the pieces were chained by the 5 % step alone, so a slow taper
      // linked an exact 2 um lead to a 3.7 um strip through its intermediate readings and
      // the whole link took the exact mean of the lead (the E8-7 key); now a piece's group
      // holds only pieces within the tolerance of ITS separation.
      struct GroupedPiece
      {
        const SubPiece *piece;
        int side;
        double mean_separation;
      };
      std::vector<GroupedPiece> grouped;
      for (const auto *list : {&sub_a, &sub_b})
      {
        const int side = list == &sub_a ? 0 : 1;  // A at offset 0, B at the separation
        for (const auto &piece : *list)
        {
          if (piece.constant && piece.interacting)
          {
            grouped.push_back({&piece, side, piece.separation});
          }
        }
      }
      std::sort(grouped.begin(), grouped.end(),
                [](const GroupedPiece &u, const GroupedPiece &v)
                { return u.mean_separation < v.mean_separation; });
      // Groups as index lists into `grouped` (ascending separation within a group).
      std::vector<std::vector<std::size_t>> groups;
      {
        struct ExactGroup
        {
          std::vector<std::size_t> members;
          double value = 0.0;  // length-weighted mean of the exact separations
          double length = 0.0;
        };
        std::vector<ExactGroup> exact_groups;
        const double parameter_tolerance = kSignatureParameterToleranceOverRadius * R;
        for (std::size_t i = 0; i < grouped.size(); i++)
        {
          if (!grouped[i].piece->exact)
          {
            continue;
          }
          if (exact_groups.empty() ||
              !quantizer.Less(
                  grouped[i].mean_separation -
                      grouped[exact_groups.back().members.back()].mean_separation,
                  parameter_tolerance))
          {
            exact_groups.push_back({});
          }
          auto &group = exact_groups.back();
          const double length =
              grouped[i].piece->interval.second - grouped[i].piece->interval.first;
          group.members.push_back(i);
          group.value = (group.value * group.length + grouped[i].mean_separation * length) /
                        (group.length + length);
          group.length += length;
        }
        std::vector<std::size_t> chord_pieces;
        for (std::size_t i = 0; i < grouped.size(); i++)
        {
          if (grouped[i].piece->exact)
          {
            continue;
          }
          // The nearest exact group within the pair tolerance of the piece's separation.
          std::size_t nearest = exact_groups.size();
          for (std::size_t g = 0; g < exact_groups.size(); g++)
          {
            const double gap = std::abs(grouped[i].mean_separation - exact_groups[g].value);
            if (!quantizer.Less(kPairSeparationTolerance * exact_groups[g].value, gap) &&
                (nearest == exact_groups.size() ||
                 gap < std::abs(grouped[i].mean_separation - exact_groups[nearest].value)))
            {
              nearest = g;
            }
          }
          if (nearest < exact_groups.size())
          {
            exact_groups[nearest].members.push_back(i);
          }
          else
          {
            chord_pieces.push_back(i);
          }
        }
        for (auto &group : exact_groups)
        {
          std::sort(group.members.begin(), group.members.end());
          groups.push_back(std::move(group.members));
        }
        for (std::size_t p = 0; p < chord_pieces.size();)
        {
          std::size_t q = p + 1;
          while (q < chord_pieces.size() &&
                 !quantizer.Less(kPairSeparationTolerance *
                                     grouped[chord_pieces[q - 1]].mean_separation,
                                 grouped[chord_pieces[q]].mean_separation -
                                     grouped[chord_pieces[q - 1]].mean_separation))
          {
            q++;
          }
          groups.emplace_back(chord_pieces.begin() + static_cast<std::ptrdiff_t>(p),
                              chord_pieces.begin() + static_cast<std::ptrdiff_t>(q));
          p = q;
        }
        // Ascending order of the groups' smallest separation (the former group order).
        std::sort(groups.begin(), groups.end(),
                  [&](const std::vector<std::size_t> &u, const std::vector<std::size_t> &v)
                  {
                    return std::make_pair(grouped[u.front()].mean_separation, u.front()) <
                           std::make_pair(grouped[v.front()].mean_separation, v.front());
                  });
      }
      // Invariant of the grouping (USER decision 203 follow-up, 2026-10-02): mutually
      // facing pieces belong to the same group. Two non-exact pieces P and Q on the two
      // sides face each other mutually when every foot of P lies on Q's run, every foot of
      // Q on P's run, and each is the only non-exact piece of its run facing the other's
      // run. Where they straddle a group threshold (the two sides sample one geometry with
      // different sample sets and agree to ~1e-3 only: at the 5 % boundary of a taper one
      // side's piece joined the exact group and its facing piece the chord group, leaving
      // both unpaired, a strip edge read as an IsolatedEdge sliver), the piece sampled over
      // exactly its own run (interval = the whole run) decides and the other joins its
      // group; when both or neither span their run, the piece on the lower run index
      // decides. Each piece is in at most one mutual pair and every move is decided on the
      // original grouping: deterministic and order-independent. The 5 % grouping itself is
      // unchanged.
      if (groups.size() > 1)
      {
        std::vector<std::size_t> group_of(grouped.size());
        for (std::size_t g = 0; g < groups.size(); g++)
        {
          for (const std::size_t i : groups[g])
          {
            group_of[i] = g;
          }
        }
        // Non-exact pieces per (side, run).
        std::map<std::pair<int, std::size_t>, std::vector<std::size_t>> by_run;
        for (std::size_t i = 0; i < grouped.size(); i++)
        {
          if (!grouped[i].piece->exact)
          {
            by_run[{grouped[i].side, grouped[i].piece->run}].push_back(i);
          }
        }
        // The single non-exact piece of the other side on run `run` whose feet all lie on
        // `facing`.
        auto PartnerPiece = [&](int side, std::size_t run,
                                std::size_t facing) -> std::optional<std::size_t>
        {
          const auto it = by_run.find({1 - side, run});
          if (it == by_run.end())
          {
            return std::nullopt;
          }
          std::optional<std::size_t> found;
          for (const std::size_t j : it->second)
          {
            if (grouped[j].piece->foot_run && *grouped[j].piece->foot_run == facing)
            {
              if (found)
              {
                return std::nullopt;  // not the only one
              }
              found = j;
            }
          }
          return found;
        };
        auto SpansRun = [&](const SubPiece &piece)
        {
          return piece.interval.first <= Tol() &&
                 piece.interval.second >= runs[piece.run].length - Tol();
        };
        std::vector<std::pair<std::size_t, std::size_t>> moves;  // (piece, target group)
        for (std::size_t i = 0; i < grouped.size(); i++)
        {
          const SubPiece &P = *grouped[i].piece;
          if (P.exact || !P.foot_run || self)
          {
            continue;
          }
          const auto q = PartnerPiece(grouped[i].side, *P.foot_run, P.run);
          if (!q || group_of[*q] == group_of[i] ||
              PartnerPiece(grouped[*q].side, P.run, *P.foot_run) != i)
          {
            continue;
          }
          const SubPiece &Q = *grouped[*q].piece;
          const bool p_spans = SpansRun(P), q_spans = SpansRun(Q);
          const bool q_decides = p_spans == q_spans ? Q.run < P.run : q_spans;
          if (q_decides)
          {
            moves.emplace_back(i, group_of[*q]);
          }
        }
        for (const auto &[i, target] : moves)
        {
          auto &from = groups[group_of[i]];
          from.erase(std::find(from.begin(), from.end(), i));
          groups[target].insert(
              std::upper_bound(groups[target].begin(), groups[target].end(), i), i);
          const SubPiece &P = *grouped[i].piece;
          moved_count++;
          moved_length += P.interval.second - P.interval.first;
          if (moved_count <= 50)
          {
            const Point3D p = runs[P.run].At(P.interval.first);
            moved_sites << std::setprecision(6) << " run " << P.run << " s ["
                        << P.interval.first << ", " << P.interval.second << "] at (" << p[0]
                        << ", " << p[1] << ", " << p[2] << ") " << P.separation;
          }
        }
        if (!moves.empty())
        {
          groups.erase(std::remove_if(groups.begin(), groups.end(),
                                      [](const std::vector<std::size_t> &g)
                                      { return g.empty(); }),
                       groups.end());
        }
      }
      if (std::getenv("PALACE_IDENTIFICATION_DEBUG_GROUPS") && input.log &&
          groups.size() > 1)
      {
        std::ostringstream dbg;
        dbg << "  DEBUG groups of chains " << A.id << " - " << B.id << ":\n";
        for (const auto &group : groups)
        {
          double length = 0.0;
          for (const std::size_t i : group)
          {
            length += grouped[i].piece->interval.second - grouped[i].piece->interval.first;
          }
          dbg << "    group mean " << grouped[group.front()].mean_separation << " .. "
              << grouped[group.back()].mean_separation << " pieces " << group.size()
              << " length " << length << "\n";
          if (group.size() <= 4)
          {
            for (const std::size_t i : group)
            {
              const auto &pc = *grouped[i].piece;
              const Run &rr = runs[pc.run];
              dbg << "      side " << grouped[i].side << " run " << pc.run << " s ["
                  << pc.interval.first << ", " << pc.interval.second << "] of " << rr.length
                  << " at (" << rr.At(pc.interval.first)[0] << ", "
                  << rr.At(pc.interval.first)[1] << ") " << (pc.exact ? "exact" : "chord")
                  << " mean " << grouped[i].mean_separation << "\n";
            }
          }
        }
        input.log(dbg.str());
      }
      for (const auto &group : groups)
      {
        // Lead piece (on A when the group has one, else on B) and the lateral A -> B there:
        // the frame of a two-edge feature (side 0 = chain A); a chain facing itself takes
        // its foot outside the self-pair neighbourhood.
        const GroupedPiece *lead = nullptr;
        for (const std::size_t i : group)
        {
          if (!lead || (lead->side == 1 && grouped[i].side == 0))
          {
            lead = &grouped[i];
          }
        }
        const bool lead_on_a = lead->side == 0;
        const Chain &lead_other = lead_on_a ? B : A;
        const double lead_s =
            0.5 * (lead->piece->interval.first + lead->piece->interval.second);
        const Point3D pa = runs[lead->piece->run].At(lead_s);
        const Chain &lead_chain = lead_on_a ? A : B;
        const auto qb = ClosestPointOnChain(
            lead_other, pa,
            self
                ? std::optional<double>(
                      lead_chain.run_offset[RunIndexInChain(lead_chain, lead->piece->run)] +
                      lead_s)
                : std::nullopt,
            self ? std::optional<double>(reach) : std::nullopt);
        if (!std::isfinite(qb.distance))
        {
          continue;
        }
        const Point3D lateral = Normalize(Sub(runs[qb.run].At(qb.s), pa));
        PairLink link;
        link.chain_a = A.id;
        link.chain_b = B.id;
        link.run_a = lead_on_a ? lead->piece->run : qb.run;
        link.run_b = lead_on_a ? qb.run : lead->piece->run;
        link.lateral_ab = lead_on_a ? lateral : Scale(-1.0, lateral);
        link.lead_point = lead_on_a ? pa : runs[qb.run].At(qb.s);
        for (const std::size_t i : group)
        {
          const SubPiece &piece = *grouped[i].piece;
          link.pieces[static_cast<std::size_t>(grouped[i].side)].push_back(
              {piece.run, piece.interval, piece.curved, piece.max_kappa, piece.separation,
               piece.exact});
        }
        // One separation per link and curvature class (decision 85(1)): the length-weighted
        // mean of the class's EXACT piece separations when it has any, else of its chord
        // readings; a class without pieces takes the other class's value (the composition
        // reads the class of the whole cross-section).
        link.separation = LinkSeparations(link.pieces);
        if (std::getenv("PALACE_IDENTIFICATION_DEBUG_EXACT") && input.log)
        {
          std::ostringstream dbg;
          dbg << "  DEBUG link " << A.id << " - " << B.id << " separation straight "
              << link.separation[0].value
              << (link.separation[0].exact ? " exact" : " chord") << " curved "
              << link.separation[1].value
              << (link.separation[1].exact ? " exact" : " chord") << "\n";
          for (int side = 0; side < 2; side++)
          {
            for (const auto &piece : link.pieces[static_cast<std::size_t>(side)])
            {
              const Run &rr = runs[piece.run];
              dbg << "    side " << side << " run " << piece.run << " len " << rr.length
                  << " s [" << piece.interval.first << ", " << piece.interval.second
                  << "] arc " << RunArc(piece.run);
              if (RunArc(piece.run) >= 0)
              {
                const Arc &arc = arcs[static_cast<std::size_t>(RunArc(piece.run))];
                dbg << " (r " << arc.radius << " c " << arc.center[0] << ","
                    << arc.center[1] << " joints " << arc.joints.size() << " turn "
                    << arc.turn << ")";
              }
              dbg << (piece.curved ? " curved" : " straight")
                  << (piece.exact ? " EXACT " : " chord ") << piece.separation << " at ("
                  << rr.At(piece.interval.first)[0] << ", "
                  << rr.At(piece.interval.first)[1] << ")\n";
            }
          }
          input.log(dbg.str());
        }
        self_pairs += self ? 1 : 0;
        pair_links.push_back(std::move(link));
      }
    }
  }
  stage.End(
      std::to_string(chains.size()) + " chains, " + std::to_string(chain_pairs_examined) +
      " chain pairs within reach, " + std::to_string(chain_pairs_paired) +
      " with facing pieces, " + std::to_string(pair_links.size()) + " pair links (" +
      std::to_string(self_pairs) + " of chains facing themselves); exact stretches: " +
      std::to_string(stretch_split_blocks) + " sub-pieces cut into stretches, " +
      std::to_string(exact_stretch_count) + " exact stretches, " +
      std::to_string(sliver_count) +
      " sub-R non-exact slivers between agreeing exact stretches" +
      (sliver_count > 0
           ? " (" + std::to_string(sliver_length) + " length):" + sliver_sites.str()
           : std::string("")) +
      "; group invariant: " + std::to_string(moved_count) +
      " pieces joined their facing partner's group" +
      (moved_count > 0
           ? " (" + std::to_string(moved_length) + " length):" + moved_sites.str()
           : std::string("")));
}

// Interval of run parameter s where the distance from run(s) to the segment [a, b] is below
// the given distance (convex sublevel set).
std::vector<Interval> Identifier::RunIntervalWithin(std::size_t run, const Point3D &a,
                                                    const Point3D &b, double distance) const
{
  const Run &r = runs[run];
  auto f = [&](double s) { return PointSegmentDistance(r.At(s), a, b); };
  std::vector<Interval> result;
  if (auto interval = ConvexSublevelInterval(f, 0.0, r.length, distance, quantizer))
  {
    result.push_back(*interval);
  }
  return result;
}

// The cluster machinery's form of RunIntervalWithin (option A): the run's points on its arc
// where it lies on one, the piece an arc where it is one; two chords take the former path.
std::vector<Interval> Identifier::RunIntervalWithinPiece(std::size_t run,
                                                         const CurvePiece &piece,
                                                         double distance) const
{
  const Run &r = runs[run];
  if (!run_arcs[run] && !piece.arc)
  {
    return RunIntervalWithin(run, piece.a, piece.b, distance);
  }
  auto f = [&](double s) { return PointPieceDistance(RunPoint(run, s), piece); };
  double lo = 0.0, hi = r.length;
  if (!run_arcs[run])
  {
    // A straight run against an arc piece: the sublevel set lies inside the convex sublevel
    // set of the distance to the piece's chord below distance + sagitta.
    auto chord = [&](double s) { return PointSegmentDistance(r.At(s), piece.a, piece.b); };
    const auto bracket = ConvexSublevelInterval(
        chord, 0.0, r.length, distance + piece.Sagitta() + 2.0 * Tol(), quantizer);
    if (!bracket)
    {
      return {};
    }
    lo = bracket->first;
    hi = bracket->second;
  }
  // The sample spacing in run parameter: the arc length of an arc run exceeds its chord.
  const double scale =
      run_arcs[run] && r.length > 0.0 ? run_arcs[run]->piece.Length() / r.length : 1.0;
  return SampledSublevelIntervals(f, lo, hi, distance, quantizer,
                                  ArcSampleStep() / std::max(scale, 1.0e-300));
}

// The runs of every fitted arc (rounded corners and bends alike) as pieces of the fitted
// circle: a run belongs to an arc when all its segments do (a run is a maximal collinear
// sequence, so a chord between two joints is one run); its end angles are those of its
// end points about the centre in the plane of the run's process normal.
void Identifier::BuildRunArcGeometry()
{
  run_arcs.assign(runs.size(), std::nullopt);
  std::vector<int> any_arc(input.segments.size(), -1);
  for (std::size_t a = 0; a < arcs.size(); a++)
  {
    for (const std::size_t s : arcs[a].segments)
    {
      any_arc[s] = static_cast<int>(a);
    }
  }
  for (std::size_t r = 0; r < runs.size(); r++)
  {
    const Run &run = runs[r];
    if (run.segments.empty() || run.length <= 0.0)
    {
      continue;
    }
    const int arc = any_arc[run.segments.front().segment];
    if (arc < 0 || std::any_of(
                       run.segments.begin(), run.segments.end(),
                       [&](const RunSegment &rs) { return any_arc[rs.segment] != arc; }))
    {
      continue;
    }
    const Arc &fitted = arcs[static_cast<std::size_t>(arc)];
    if (!(fitted.radius > 0.0))
    {
      continue;
    }
    const Point3D n = Normalize(run.process_normal);
    Point3D d0 = Sub(run.start, fitted.center), d1 = Sub(run.end, fitted.center);
    d0 = Sub(d0, Scale(Dot(d0, n), n));
    d1 = Sub(d1, Scale(Dot(d1, n), n));
    if (Norm(d0) <= 1.0e-12 * fitted.radius || Norm(d1) <= 1.0e-12 * fitted.radius)
    {
      continue;
    }
    RunArcGeometry geometry;
    geometry.arc = arc;
    geometry.piece.center = fitted.center;
    geometry.piece.radius = fitted.radius;
    geometry.piece.u = Normalize(d0);
    geometry.piece.v = Normalize(Cross(n, geometry.piece.u));
    geometry.piece.theta0 = 0.0;
    geometry.piece.theta1 =
        std::atan2(Dot(d1, geometry.piece.v), Dot(d1, geometry.piece.u));
    if (std::abs(geometry.piece.theta1) <= 1.0e-12)
    {
      continue;  // degenerate chord
    }
    run_arcs[r] = geometry;
    max_run_sagitta = std::max(max_run_sagitta, WholeRunPiece(r).Sagitta());
  }
}

// Parallel classes -> elementary longitudinal intervals -> connected active runs.
void Identifier::BuildTranslationalFeatures()
{
  struct ClassMember
  {
    std::size_t run;
    double u0, u1, w;
    int gap_sign;
  };
  // Direction classes (BuildDirectionClasses: the parallel relation within the cosine
  // tolerance, per metal plane): parallel runs of two planes never pair (decision 82(1));
  // metal of another plane within 2R is the CrossLayer exclusion.
  std::map<int, std::vector<std::size_t>> classes;
  for (std::size_t r = 0; r < runs.size(); r++)
  {
    if (runs[r].excluded || run_direction_class[r] < 0)
    {
      continue;  // chains with joints pair through the curved-edge chain rule
    }
    classes[run_direction_class[r]].push_back(r);
  }
  const double interaction = kInteractionDistanceOverRadius * R;
  std::size_t total_members = 0, total_spans = 0;
  for (const auto &[key, members_in_class] : classes)
  {
    (void)key;
    total_members += members_in_class.size();
    if (members_in_class.size() < 2)
    {
      continue;
    }
    // Class axis: tangent of the run with the smallest start point.
    std::size_t reference = members_in_class.front();
    for (const std::size_t r : members_in_class)
    {
      if (std::make_pair(runs[r].start, runs[r].end) <
          std::make_pair(runs[reference].start, runs[reference].end))
      {
        reference = r;
      }
    }
    const Point3D axis = SignCanonical(runs[reference].tangent);
    const Point3D lateral = Normalize(Cross(n_ref, axis));
    std::vector<ClassMember> members;
    std::vector<double> splits;
    for (const std::size_t r : members_in_class)
    {
      double u0 = Dot(runs[r].start, axis), u1 = Dot(runs[r].end, axis);
      if (u0 > u1)
      {
        std::swap(u0, u1);
      }
      const double w = Dot(Scale(0.5, Add(runs[r].start, runs[r].end)), lateral);
      const int gap_sign = Dot(runs[r].gap_direction, lateral) > 0.0 ? 1 : -1;
      members.push_back({r, u0, u1, w, gap_sign});
      splits.push_back(u0);
      splits.push_back(u1);
    }
    std::sort(splits.begin(), splits.end());
    std::vector<double> unique_splits;
    for (const double u : splits)
    {
      if (unique_splits.empty() || !quantizer.Equal(u, unique_splits.back()))
      {
        unique_splits.push_back(u);
      }
    }
    // Components per elementary interval, then merge consecutive equal components. The
    // members active on an interval (u0 <= lo + Tol and u1 >= hi - Tol) are maintained by a
    // sweep (both bounds grow with k: a member enters once and leaves once), in member
    // order; the interacting pairs are found through the members sorted by lateral offset
    // and united in the former (i, j) order.
    struct Span
    {
      std::vector<std::size_t> component;  // member indices sorted by w
      double lo, hi;
    };
    std::vector<Span> spans;
    std::vector<std::size_t> by_u0(members.size());
    std::iota(by_u0.begin(), by_u0.end(), 0);
    std::sort(
        by_u0.begin(), by_u0.end(), [&](std::size_t a, std::size_t b)
        { return std::make_pair(members[a].u0, a) < std::make_pair(members[b].u0, b); });
    std::set<std::size_t> active_set;
    std::size_t next_entering = 0;
    std::vector<std::pair<std::size_t, std::size_t>> pairs;
    for (std::size_t k = 0; k + 1 < unique_splits.size(); k++)
    {
      stage.Progress(k, unique_splits.size());
      const double lo = unique_splits[k], hi = unique_splits[k + 1];
      while (next_entering < by_u0.size() && members[by_u0[next_entering]].u0 <= lo + Tol())
      {
        active_set.insert(by_u0[next_entering]);
        next_entering++;
      }
      std::vector<std::size_t> active;
      for (auto it = active_set.begin(); it != active_set.end();)
      {
        if (members[*it].u1 < hi - Tol())
        {
          it = active_set.erase(it);  // never active again: hi only grows
        }
        else
        {
          active.push_back(*it);
          ++it;
        }
      }
      if (active.size() < 2)
      {
        continue;
      }
      UnionFind uf(active.size());
      std::vector<std::size_t> by_w(active.size());
      std::iota(by_w.begin(), by_w.end(), 0);
      std::sort(by_w.begin(), by_w.end(), [&](std::size_t a, std::size_t b)
                { return members[active[a]].w < members[active[b]].w; });
      pairs.clear();
      for (std::size_t p = 0; p < by_w.size(); p++)
      {
        for (std::size_t q = p + 1;
             q < by_w.size() && members[active[by_w[q]]].w - members[active[by_w[p]]].w <
                                    interaction + 2.0 * Tol();
             q++)
        {
          const std::size_t i = std::min(by_w[p], by_w[q]), j = std::max(by_w[p], by_w[q]);
          const auto &mi = members[active[i]], &mj = members[active[j]];
          if (runs[mi.run].chain != runs[mj.run].chain &&
              quantizer.Less(std::abs(mi.w - mj.w), interaction))
          {
            pairs.emplace_back(i, j);
          }
        }
      }
      std::sort(pairs.begin(), pairs.end());
      for (const auto &[i, j] : pairs)
      {
        uf.Union(i, j);
      }
      std::map<std::size_t, std::vector<std::size_t>> components;
      for (std::size_t i = 0; i < active.size(); i++)
      {
        components[uf.Find(i)].push_back(active[i]);
      }
      for (auto &[root, component] : components)
      {
        (void)root;
        if (component.size() < 2)
        {
          continue;
        }
        std::sort(component.begin(), component.end(), [&](std::size_t a, std::size_t b)
                  { return members[a].w < members[b].w; });
        if (!spans.empty() && spans.back().component == component &&
            quantizer.Equal(spans.back().hi, lo))
        {
          spans.back().hi = hi;
        }
        else
        {
          spans.push_back({component, lo, hi});
        }
      }
    }
    for (const auto &span : spans)
    {
      TranslationalSpan record;
      record.lo = span.lo;
      record.hi = span.hi;
      record.axis = axis;
      record.lateral = lateral;
      for (const std::size_t m : span.component)
      {
        const auto &member = members[m];
        const Run &run = runs[member.run];
        const bool forward = Dot(run.tangent, axis) > 0.0;
        const double us = Dot(run.start, axis);
        double s0 = forward ? span.lo - us : us - span.hi;
        double s1 = forward ? span.hi - us : us - span.lo;
        s0 = std::clamp(s0, 0.0, run.length);
        s1 = std::clamp(s1, 0.0, run.length);
        record.members.push_back({member.run, member.w, member.gap_sign, {s0, s1}});
      }
      translational_spans.push_back(std::move(record));
    }
    total_spans += spans.size();
  }
  stage.End(std::to_string(classes.size()) + " direction classes, " +
            std::to_string(total_members) + " rigid runs, " + std::to_string(total_spans) +
            " translational spans");
}

// Pair / stack assembly (decision 82(2)). The pair links of the bent-pair rule and the
// translational spans of the rigid parallel runs are the pairwise facing relations of the
// perimeter; where several of them share a run over a common interval, the edges are one
// translation-invariant cross-section of k >= 3 edges with consecutive separations below
// 2R — a stack — and the stack is ONE feature (ParallelEdgeCluster, or
// CurvedParallelEdgeCluster along a bend, with the ordered offsets / R, the gap pattern,
// the conductor identities and the bend radius in its signature), on straight runs and
// along curved routes alike. The pairwise candidates inside a stack are superseded by the
// stack: no two claims of the pair priority ever overlap, so claim resolution never keeps
// or drops a feature by its id (decision 78). A link or span alone is the pair / span
// feature of the former rules, unchanged.
void Identifier::BuildPairsAndStacks()
{
  // Chain intervals taken before the assembly: cluster portions (priority 0 claims) and the
  // windows of the vertex features outside clusters (priority 1, claimed in Assign).
  taken_on_chain.clear();
  auto Take = [&](std::size_t r, const Interval &interval)
  {
    const Chain &C = chains[chain_index.at(runs[r].chain)];
    const double offset = C.run_offset[runs[r].index_in_chain];
    taken_on_chain[C.id].emplace_back(offset + interval.first, offset + interval.second);
  };
  for (const auto &claimed : cluster_claimed)
  {
    for (const auto &[r, interval] : claimed)
    {
      Take(r, interval);
    }
  }
  for (const auto &site : sites)
  {
    if (site.cluster < 0)
    {
      for (const auto &[r, interval] : site.window)
      {
        Take(r, interval);
      }
    }
  }
  for (auto &[chain, list] : taken_on_chain)
  {
    (void)chain;
    list = MergeIntervals(std::move(list), Tol());
  }
  const std::size_t n_links = pair_links.size(), n_spans = translational_spans.size();
  const std::size_t n_items = n_links + n_spans;
  struct RunInterval
  {
    double lo, hi;
    std::size_t item;
  };
  std::map<std::size_t, std::vector<RunInterval>> by_run;
  for (std::size_t i = 0; i < n_links; i++)
  {
    for (const auto &side : pair_links[i].pieces)
    {
      for (const auto &piece : side)
      {
        by_run[piece.run].push_back({piece.interval.first, piece.interval.second, i});
      }
    }
  }
  for (std::size_t j = 0; j < n_spans; j++)
  {
    for (const auto &member : translational_spans[j].members)
    {
      by_run[member.run].push_back(
          {member.interval.first, member.interval.second, n_links + j});
    }
  }
  // Items sharing a run over an interval longer than the decision quantum are one
  // cross-section (union-find over the items; the result does not depend on the order).
  UnionFind uf(n_items);
  for (auto &[run, list] : by_run)
  {
    (void)run;
    std::sort(list.begin(), list.end(), [](const RunInterval &a, const RunInterval &b)
              { return std::tie(a.lo, a.hi, a.item) < std::tie(b.lo, b.hi, b.item); });
    for (std::size_t i = 0; i < list.size(); i++)
    {
      for (std::size_t j = i + 1; j < list.size() && list[j].lo < list[i].hi - Tol(); j++)
      {
        if (list[i].item != list[j].item)
        {
          uf.Union(list[i].item, list[j].item);
        }
      }
    }
  }
  std::map<std::size_t, std::pair<std::vector<std::size_t>, std::vector<std::size_t>>>
      components;  // root -> (links, spans) in ascending item order
  std::vector<std::size_t> component_order;
  for (std::size_t item = 0; item < n_items; item++)
  {
    const std::size_t root = uf.Find(item);
    auto [it, inserted] = components.emplace(
        root, std::make_pair(std::vector<std::size_t>{}, std::vector<std::size_t>{}));
    if (inserted)
    {
      component_order.push_back(root);
    }
    if (item < n_links)
    {
      it->second.first.push_back(item);
    }
    else
    {
      it->second.second.push_back(item - n_links);
    }
  }
  std::size_t pairs = 0, spans = 0, stacks = 0;
  const std::size_t features_before = features.size();
  for (std::size_t k = 0; k < component_order.size(); k++)
  {
    stage.Progress(k, component_order.size(), "pair / stack components");
    const auto &[links, span_items] = components.at(component_order[k]);
    // Every component — a lone two-edge link or span included — is assembled per
    // cross-section with the taken rule (decision 85(2)): a pair side whose partner is
    // cluster / window metal is no pair there and returns to the single-edge remainder,
    // where the cluster extension tests it (the former lone-link path left such a side to
    // the mutual-sides trimming of the claim resolution, after the extension).
    if (links.size() == 1 && span_items.empty())
    {
      pairs++;
    }
    else if (links.empty() && span_items.size() == 1 &&
             translational_spans[span_items.front()].members.size() == 2)
    {
      spans++;
    }
    else
    {
      stacks++;
    }
    AssembleStack(links, span_items);
  }
  stage.End(std::to_string(component_order.size()) +
            " cross-section components: " + std::to_string(pairs) + " single links, " +
            std::to_string(spans) + " single spans, " + std::to_string(stacks) +
            " multi-link (stacks, chains facing themselves, mixed), all assembled per "
            "cross-section, " +
            std::to_string(features.size() - features_before) + " features, " +
            std::to_string(features.size()) + " features so far");
}

// One cross-section component of several links / spans (or a chain facing itself). Along
// every member chain the link pieces are cut at every breakpoint — the piece ends of every
// link on the chain and the images of the other chains' breakpoints through the links
// (the closest point on the partner chain, so that a stack end on one edge cuts every
// edge of the stack at the same cross-section) — and the composition of the cross-section
// at the middle of every elementary interval is read off the links active there: from the
// chain, the partner chains reached through active links (breadth first), each at the foot
// of the previous point, ordered by their lateral position along the in-plane normal. The
// signature offsets are the sums of the CONSECUTIVE links' separations (path independent:
// a non-consecutive link within 2R, e.g. the outer edges of a 4-edge stack at 4 um and R =
// 2.1 um, does not enter), gap sides and conductors read at the members; one feature per
// (type, signature without the bend radius, curvature class, member chains), its
// RadiusOverR the tightest windowed radius over its pieces as for the two-edge pairs. A
// member taken by a cluster portion or a vertex window (taken_on_chain) is not part of the
// cross-section and is not traversed: the stack ends at the claim boundary (a breakpoint on
// every member through the links) and the remaining members are recomposed there (a
// smaller stack, a pair, or nothing). Sides are ranked along the lateral axis in the
// canonical order of the signature (chirality +1; a symmetric cross-section, chirality 0,
// puts the member with the smallest (chain, position) on side 0).
void Identifier::AssembleStack(const std::vector<std::size_t> &link_items,
                               const std::vector<std::size_t> &span_items)
{
  const double reach =
      kInteractionDistanceOverRadius * R * (1.0 + kPairSeparationTolerance);
  const double neighbourhood = kSelfPairNeighbourhoodOverRadius * R;
  const bool debug_log = std::getenv("PALACE_IDENTIFICATION_DEBUG") && input.log;
  // Phase timing of a component (PALACE_IDENTIFICATION_DEBUG_TIMING), for the chip-scale
  // routes.
  const bool debug_timing = std::getenv("PALACE_IDENTIFICATION_DEBUG_TIMING") && input.log;
  auto timer_started = std::chrono::steady_clock::now();
  std::string timer_phase = "setup";
  auto StackTimer = [&](const std::string &phase)
  {
    if (debug_timing)
    {
      const auto now = std::chrono::steady_clock::now();
      const double seconds = std::chrono::duration<double>(now - timer_started).count();
      if (seconds > 0.05)
      {
        std::ostringstream line;
        line << "    stack component (" << link_items.size() << " links, "
             << span_items.size() << " spans): " << timer_phase << " " << std::fixed
             << std::setprecision(2) << seconds << " s\n";
        input.log(line.str());
      }
      timer_started = now;
      timer_phase = phase;
    }
  };
  struct ELink
  {
    int chain_a, chain_b;
    LinkSeparation separation;  // per curvature class
    std::array<std::vector<LinkPiece>, 2> pieces;
    bool Self() const { return chain_a == chain_b; }
  };
  std::vector<ELink> elinks;
  for (const std::size_t i : link_items)
  {
    const PairLink &link = pair_links[i];
    elinks.push_back({link.chain_a, link.chain_b, link.separation, link.pieces});
  }
  for (const std::size_t j : span_items)
  {
    const TranslationalSpan &span = translational_spans[j];
    for (std::size_t m = 0; m + 1 < span.members.size(); m++)
    {
      const auto &lower = span.members[m], &upper = span.members[m + 1];
      // Two rigid parallel runs: the lateral offset difference is exact for both classes
      // when it holds over the span within the signature parameter tolerance. The members
      // of a parallel class may be tilted by up to the cosine tolerance (1.4e-4 rad) and
      // their lateral offsets are read at the runs' midpoints, so the lines' separation at
      // the span's ends deviates by the tilt times the half-length (above 1e-3 R beyond
      // ~7 R of span at the maximum tilt): such a pair reads the midpoint value, not exact
      // (USER decision 203 follow-up, 2026-10-02).
      auto LateralAt = [&](const TranslationalSpan::Member &member, double u)
      {
        const Run &run = runs[member.run];
        const double along = Dot(run.tangent, span.axis);
        const Point3D p =
            Add(run.start, Scale((u - Dot(run.start, span.axis)) / along, run.tangent));
        return Dot(p, span.lateral);
      };
      const double midpoint_value = upper.w - lower.w;
      double deviation = 0.0;
      for (const double u : {span.lo, span.hi})
      {
        deviation = std::max(deviation, std::abs(LateralAt(upper, u) - LateralAt(lower, u) -
                                                 midpoint_value));
      }
      const bool tilt_exact =
          quantizer.Less(deviation, kSignatureParameterToleranceOverRadius * R);
      const ClassSeparation exact{midpoint_value, tilt_exact, true};
      ELink e{runs[lower.run].chain, runs[upper.run].chain, {exact, exact}, {}};
      e.pieces[0].push_back(
          {lower.run, lower.interval, false, 0.0, exact.value, tilt_exact});
      e.pieces[1].push_back(
          {upper.run, upper.interval, false, 0.0, exact.value, tilt_exact});
      elinks.push_back(std::move(e));
    }
  }
  auto ChainOf = [&](int id) -> const Chain & { return chains[chain_index.at(id)]; };
  auto ChainInterval = [&](const LinkPiece &piece)
  {
    const Chain &c = ChainOf(runs[piece.run].chain);
    const double offset = c.run_offset[runs[piece.run].index_in_chain];
    return Interval{offset + piece.interval.first, offset + piece.interval.second};
  };
  // Chain intervals of every link side, sorted (the pieces of one side are disjoint): the
  // membership tests below are binary searches instead of scans over every piece.
  std::vector<std::array<std::vector<Interval>, 2>> side_intervals(elinks.size());
  for (std::size_t e = 0; e < elinks.size(); e++)
  {
    for (int side = 0; side < 2; side++)
    {
      auto &list = side_intervals[e][static_cast<std::size_t>(side)];
      for (const auto &piece : elinks[e].pieces[static_cast<std::size_t>(side)])
      {
        list.push_back(ChainInterval(piece));
      }
      // Touching pieces (the sub-pieces of consecutive runs, whose analytic ends differ by
      // roundoff below the signature grid) are one interval: the composition cannot change
      // inside it, so its interior run joints are no breakpoints.
      list = MergeIntervals(std::move(list), kSignatureLengthQuantumOverRadius * R);
    }
  }
  // The interval of the link side holding x within the tolerance (closed), or null.
  auto IntervalAt = [&](std::size_t e, int side, double x,
                        double tolerance) -> const Interval *
  {
    const auto &list = side_intervals[e][static_cast<std::size_t>(side)];
    auto it = std::upper_bound(list.begin(), list.end(), x,
                               [](double value, const Interval &interval)
                               { return value < interval.first; });
    // Candidates: the interval starting at or before x (it - 1) and, within the tolerance,
    // the one starting just after.
    for (auto candidate : {it == list.begin() ? list.end() : it - 1, it})
    {
      if (candidate != list.end() && x >= candidate->first - tolerance &&
          x <= candidate->second + tolerance)
      {
        return &*candidate;
      }
    }
    return nullptr;
  };
  auto ActiveAt = [&](std::size_t e, int side, double x)
  { return IntervalAt(e, side, x, Tol()) != nullptr; };
  auto StrictlyActiveAt = [&](std::size_t e, int side, double x)
  {
    const Interval *interval = IntervalAt(e, side, x, Tol());
    return interval && x > interval->first + Tol() && x < interval->second - Tol();
  };
  // Links on every chain: (link, side), in link order. Breakpoints per chain (sorted, no
  // two within the decision quantum): the piece ends of every link on the chain and the
  // ends of the higher-priority claims (cluster portions, vertex windows) on it — where a
  // member of the cross-section is taken by a cluster or a vertex window the stack ends and
  // the remaining members are recomposed (the stack-end rule of decision 82(2)).
  std::map<int, std::vector<std::pair<std::size_t, int>>> links_on_chain;
  std::map<int, std::set<double>> breakpoints;
  // Breakpoints closer than the signature grid are one cut (an elementary interval below
  // the grid is roundoff between two claim boundaries and joins its neighbour anyway).
  const double breakpoint_grid = kSignatureLengthQuantumOverRadius * R;
  auto AddBreakpoint = [&](int chain, double x)
  {
    auto &list = breakpoints[chain];
    auto it = list.lower_bound(x - breakpoint_grid);
    if (it != list.end() && std::abs(*it - x) <= breakpoint_grid)
    {
      return false;
    }
    list.insert(x);
    return true;
  };
  // An IMAGE (the foot of a breakpoint on a partner chain) within the signature parameter
  // tolerance of an existing breakpoint is the same cut (decision 93, DS-CTX-003 defect 1):
  // on bent partner chains that are not concentric the foot of a foot drifts away from its
  // source by less than the tolerance per hop, and under the 1e-6 R grid every hop
  // multiplied the breakpoints (10 seeds -> 10^6 per chain, 744 s per pass) into pm-long
  // elementary intervals of one composition. An elementary interval shorter than the
  // tolerance carries no parameter the contract resolves.
  const double image_grid = kStackImageToleranceOverRadius * R;
  auto AddImage = [&](int chain, double x)
  {
    auto &list = breakpoints[chain];
    auto it = list.lower_bound(x - image_grid);
    if (it != list.end() && *it <= x + image_grid)
    {
      stack_images_merged++;
      return false;
    }
    list.insert(x);
    return true;
  };
  for (std::size_t e = 0; e < elinks.size(); e++)
  {
    for (int side = 0; side < 2; side++)
    {
      const int chain = side == 0 ? elinks[e].chain_a : elinks[e].chain_b;
      if (elinks[e].pieces[static_cast<std::size_t>(side)].empty())
      {
        continue;
      }
      links_on_chain[chain].emplace_back(e, side);
      // Breakpoints: the ends of the merged pieces of each curvature class (a run joint
      // inside a class is no breakpoint; the curved boundary between the classes is one).
      for (const bool curved_class : {false, true})
      {
        std::vector<Interval> of_class;
        for (const auto &piece : elinks[e].pieces[static_cast<std::size_t>(side)])
        {
          if (piece.curved == curved_class)
          {
            of_class.push_back(ChainInterval(piece));
          }
        }
        for (const auto &interval :
             MergeIntervals(std::move(of_class), kSignatureLengthQuantumOverRadius * R))
        {
          AddBreakpoint(chain, interval.first);
          AddBreakpoint(chain, interval.second);
        }
      }
    }
  }
  // The taken (cluster portion / vertex window) ends are imaged on the fitted arcs (option
  // A: a cluster boundary on a bend cuts the concentric partner radially, whatever the
  // chords), their images likewise; the link piece ends keep the chord reading.
  std::map<int, std::set<double>> taken_breakpoints;
  for (const auto &[chain, links] : links_on_chain)
  {
    (void)links;
    const auto it = taken_on_chain.find(chain);
    if (it == taken_on_chain.end())
    {
      continue;
    }
    for (const auto &interval : it->second)
    {
      for (const double x : {interval.first, interval.second})
      {
        if (AddBreakpoint(chain, x))
        {
          taken_breakpoints[chain].insert(x);
        }
      }
    }
  }
  // A chain position inside a cluster portion or a vertex window (taken_on_chain: sorted,
  // disjoint) takes no part in a cross-section.
  auto Taken = [&](int chain, double x)
  {
    const auto it = taken_on_chain.find(chain);
    if (it == taken_on_chain.end())
    {
      return false;
    }
    const auto &list = it->second;
    auto iv = std::upper_bound(list.begin(), list.end(), x,
                               [](double value, const Interval &interval)
                               { return value < interval.first; });
    return iv != list.begin() && x < (iv - 1)->second - Tol() &&
           x > (iv - 1)->first + Tol();
  };
  StackTimer("images");
  // Images of the breakpoints through the links, breadth first over the chains (bounded by
  // the number of links: an image travels at most once over every link).
  {
    std::map<int, std::vector<std::pair<double, bool>>> frontier;  // (x, from a taken end)
    for (const auto &[chain, xs] : breakpoints)
    {
      const auto taken_it = taken_breakpoints.find(chain);
      for (const double x : xs)
      {
        frontier[chain].emplace_back(x, taken_it != taken_breakpoints.end() &&
                                            taken_it->second.count(x) > 0);
      }
    }
    for (std::size_t depth = 0; depth <= elinks.size() && !frontier.empty(); depth++)
    {
      std::map<int, std::vector<std::pair<double, bool>>> next;
      for (const auto &[chain, xs] : frontier)
      {
        const auto it = links_on_chain.find(chain);
        if (it == links_on_chain.end())
        {
          continue;
        }
        for (const auto &[x, from_taken] : xs)
        {
          for (const auto &[e, side] : it->second)
          {
            const ELink &link = elinks[e];
            if (!StrictlyActiveAt(e, side, x))
            {
              continue;
            }
            const int other = side == 0 ? link.chain_b : link.chain_a;
            const auto exclude = link.Self() ? std::optional<double>(x) : std::nullopt;
            const auto foot =
                from_taken
                    ? ArcAwareImage(ChainOf(chain), x, ChainOf(other), exclude, reach)
                    : ClosestPointOnChain(ChainOf(other), ChainAt(ChainOf(chain), x),
                                          exclude, reach);
            if (!std::isfinite(foot.distance) || !quantizer.Less(foot.distance, reach))
            {
              continue;
            }
            if (AddImage(other, foot.x))
            {
              next[other].emplace_back(foot.x, from_taken);
            }
          }
        }
      }
      frontier = std::move(next);
    }
  }
  StackTimer("compose-def");
  // Composition of the cross-section through (chain, x).
  struct Node
  {
    int chain;
    double x;
    Point3D point;
    std::size_t run;
    double position;  // along the lateral axis
  };
  struct Composition
  {
    std::vector<Node> nodes;  // in increasing position along `lateral`
    Point3D lateral;
    std::vector<double> offsets;
    bool curved = false;  // a member is curved at the cross-section
    bool exact = true;    // every consecutive separation entering the offsets is exact
    bool geometric_fallback = false;
    bool cap_reached = false;
    // Lateral positions of the partner feet that are cluster / window metal (taken): a
    // taken member between two members of the cross-section interrupts it.
    std::vector<double> taken_positions;
  };
  std::size_t fallbacks = 0, cap_hits = 0;
  auto Compose = [&](int chain, double x)
  {
    Composition composition;
    std::size_t run = 0;
    double s = 0.0;
    const Point3D p = ChainAt(ChainOf(chain), x, &run, &s);
    composition.lateral = Normalize(Cross(n_ref, runs[run].tangent));
    if (Taken(chain, x))
    {
      return composition;  // a cluster portion / vertex window: no cross-section of its own
    }
    composition.nodes.push_back({chain, x, p, run, 0.0});
    for (std::size_t i = 0; i < composition.nodes.size(); i++)
    {
      if (composition.nodes.size() >= kStackCompositionCap)
      {
        composition.cap_reached = true;
        break;
      }
      const Node node = composition.nodes[i];
      const auto it = links_on_chain.find(node.chain);
      if (it == links_on_chain.end())
      {
        continue;
      }
      for (const auto &[e, side] : it->second)
      {
        const ELink &link = elinks[e];
        if (!ActiveAt(e, side, node.x))
        {
          continue;
        }
        const int other = side == 0 ? link.chain_b : link.chain_a;
        const auto foot = ClosestPointOnChain(
            ChainOf(other), node.point,
            link.Self() ? std::optional<double>(node.x) : std::nullopt, reach);
        if (!std::isfinite(foot.distance) || !quantizer.Less(foot.distance, reach))
        {
          continue;  // no partner there
        }
        if (Taken(other, foot.x))
        {
          // The partner is cluster / window metal: not a member, but it stands between the
          // members on either side of it (a 2 um strip inside a cluster between two edges 4
          // um apart must not leave those two as a "pair" across it: DS-SCT-002's 60
          // UnclassifiedParallelPair of 144 um at R = 2.1 um).
          const Point3D q = runs[foot.run].At(foot.s);
          composition.taken_positions.push_back(Dot(Sub(q, p), composition.lateral));
          continue;
        }
        const bool present = std::any_of(
            composition.nodes.begin(), composition.nodes.end(),
            [&](const Node &n)
            {
              return n.chain == other &&
                     quantizer.Less(ChainArcDistance(ChainOf(other), n.x, n.x, foot.x),
                                    neighbourhood);
            });
        if (present)
        {
          continue;
        }
        const Point3D q = runs[foot.run].At(foot.s);
        composition.nodes.push_back(
            {other, foot.x, q, foot.run, Dot(Sub(q, p), composition.lateral)});
      }
    }
    auto Less = [](const Node &a, const Node &b)
    { return std::tie(a.position, a.chain, a.x) < std::tie(b.position, b.chain, b.x); };
    std::sort(composition.nodes.begin(), composition.nodes.end(), Less);
    // A taken member laterally between two members cuts the cross-section: the part holding
    // the seed (position 0) is the cross-section through (chain, x); the members beyond the
    // taken one are composed from their own chains.
    if (!composition.taken_positions.empty())
    {
      double lo = -std::numeric_limits<double>::infinity();
      double hi = std::numeric_limits<double>::infinity();
      for (const double t : composition.taken_positions)
      {
        if (t < -Tol() && t > lo)
        {
          lo = t;
        }
        if (t > Tol() && t < hi)
        {
          hi = t;
        }
      }
      composition.nodes.erase(
          std::remove_if(composition.nodes.begin(), composition.nodes.end(),
                         [&](const Node &n) { return n.position < lo || n.position > hi; }),
          composition.nodes.end());
    }
    // Class of the cross-section: curved where a member is curved (the class decides which
    // separation of every consecutive link enters the offsets, decision 85(1)).
    for (const Node &node : composition.nodes)
    {
      composition.curved = composition.curved || IsCurvedAt(ChainOf(node.chain), node.x);
    }
    // Provisional orientation (final: the canonical signature's, below), which decides only
    // for a cross-section that is its own mirror image (chirality 0). A straight one has no
    // observable orientation: the member with the smallest (chain, position along its
    // chain) is side 0 (the two-edge convention: side 0 = chain A, the lower chain index).
    // A curved one has: the bend tells its inner from its outer edge and the signature's
    // Convexity is read on side 0, while the chain ids follow the coordinate-sorted
    // canonical numbering, so that the key would make the Convexity of a symmetric curved
    // stack a function of the frame (decision 214 (i): the subdivision gate's rotation
    // variant flipped the k4 / k6 curved stacks Convex <-> Concave). The bend sense orients
    // it: the lateral axis points toward the centre of curvature (the sum over the members
    // of the signed windowed curvature toward the metal times the metal side along the
    // lateral axis), so that side 0 is the outermost edge of the bend; the key decides only
    // where that sum vanishes.
    auto Key = [](const Node &n) { return std::make_pair(n.chain, n.x); };
    bool reverse = Key(composition.nodes.back()) < Key(composition.nodes.front());
    if (composition.curved)
    {
      double toward_lateral = 0.0;
      for (const Node &node : composition.nodes)
      {
        toward_lateral -= SignedWindowedCurvature(ChainOf(node.chain), node.x) *
                          Dot(runs[node.run].gap_direction, composition.lateral);
      }
      if (toward_lateral != 0.0)
      {
        reverse = toward_lateral < 0.0;
      }
    }
    if (reverse)
    {
      std::reverse(composition.nodes.begin(), composition.nodes.end());
      composition.lateral = Scale(-1.0, composition.lateral);
      for (auto &n : composition.nodes)
      {
        n.position = -n.position;
      }
    }
    // Offsets from the consecutive links' separations of the cross-section's class.
    composition.offsets.assign(composition.nodes.size(), 0.0);
    for (std::size_t i = 0; i + 1 < composition.nodes.size(); i++)
    {
      const Node &lower = composition.nodes[i], &upper = composition.nodes[i + 1];
      std::optional<double> separation;
      const auto it = links_on_chain.find(lower.chain);
      if (it != links_on_chain.end())
      {
        for (const auto &[e, side] : it->second)
        {
          const ELink &link = elinks[e];
          const int other = side == 0 ? link.chain_b : link.chain_a;
          if (other != upper.chain || separation)
          {
            continue;
          }
          if (ActiveAt(e, side, lower.x))
          {
            const ClassSeparation &of_class = link.separation[composition.curved ? 1 : 0];
            separation = of_class.value;
            composition.exact = composition.exact && of_class.exact;
          }
        }
      }
      if (!separation)
      {
        separation = upper.position - lower.position;  // geometric (no consecutive link)
        composition.geometric_fallback = true;
        composition.exact = false;
      }
      composition.offsets[i + 1] = composition.offsets[i] + *separation;
    }
    return composition;
  };
  StackTimer("intervals");
  if (debug_timing)
  {
    std::ostringstream line;
    line << "    stack component breakpoints per chain:";
    for (const auto &[chain, list] : breakpoints)
    {
      line << " " << chain << ":" << list.size();
    }
    line << "; link side intervals:";
    for (std::size_t e = 0; e < elinks.size(); e++)
    {
      line << " " << side_intervals[e][0].size() << "/" << side_intervals[e][1].size();
      const auto &list = side_intervals[e][0];
      for (std::size_t i = 0; i < std::min<std::size_t>(list.size(), 6); i++)
      {
        line << " [" << std::setprecision(8) << list[i].first << ", " << list[i].second
             << "]";
      }
    }
    input.log(line.str() + "\n");
  }
  // Elementary intervals of every chain, their compositions, grouped into features.
  struct Assigned
  {
    int chain;
    double x0, x1;
    std::string key;
    int side;
    bool curved;
    double max_kappa;
    std::string type;
    nlohmann::json signature;
    int chirality;
    std::optional<std::string> reason;
    Point3D origin;
    std::array<Point3D, 3> axes;
    int sides;
    bool exact;
  };
  std::vector<Assigned> assigned;
  for (const auto &[chain, sorted_breakpoints] : breakpoints)
  {
    const std::vector<double> list(sorted_breakpoints.begin(), sorted_breakpoints.end());
    const auto it = links_on_chain.find(chain);
    if (it == links_on_chain.end())
    {
      continue;
    }
    const Chain &C = ChainOf(chain);
    std::vector<Assigned> on_chain;
    for (std::size_t i = 0; i + 1 < list.size(); i++)
    {
      const double x0 = list[i], x1 = list[i + 1];
      if (x1 - x0 <= breakpoint_grid)
      {
        continue;  // below the signature grid: roundoff between two cuts
      }
      const double mid = 0.5 * (x0 + x1);
      bool covered = false;
      for (const auto &[e, side] : it->second)
      {
        covered = covered || ActiveAt(e, side, mid);
      }
      if (!covered)
      {
        continue;
      }
      auto composition = Compose(chain, mid);
      const std::size_t k = composition.nodes.size();
      if (debug_log && (composition.nodes.size() < 2 || x1 - x0 > R))
      {
        std::ostringstream line;
        line << "    interval chain " << chain << " [" << x0 << ", " << x1 << "] nodes "
             << k;
        for (const auto &node : composition.nodes)
        {
          line << " (" << node.chain << " @ " << node.x << ")";
        }
        input.log(line.str() + "\n");
      }
      if (k < 2)
      {
        continue;
      }
      fallbacks += composition.geometric_fallback ? 1 : 0;
      cap_hits += composition.cap_reached ? 1 : 0;
      std::vector<TranslationalEdge> edges;
      auto Edges = [&]()
      {
        edges.clear();
        for (std::size_t n = 0; n < k; n++)
        {
          const Run &run = runs[composition.nodes[n].run];
          edges.push_back({composition.offsets[n],
                           Dot(run.gap_direction, composition.lateral) > 0.0 ? 1 : -1,
                           run.conductor, InterfaceNames(run.targets), run.boundary_law});
        }
      };
      Edges();
      auto translational = CanonicalTranslationalSignature(edges, R);
      // Sides in the canonical order: for a cross-section that is not its own mirror image
      // the lateral axis is oriented so that the signature's first edge is side 0
      // (chirality +1 always; the same sides at every cross-section of the feature whatever
      // the chain numbering), a symmetric one keeps the provisional orientation (chirality
      // 0).
      if (translational.chirality < 0)
      {
        std::reverse(composition.nodes.begin(), composition.nodes.end());
        composition.lateral = Scale(-1.0, composition.lateral);
        const double span = composition.offsets.back();
        for (std::size_t n = 0; n < k; n++)
        {
          composition.nodes[n].position = -composition.nodes[n].position;
        }
        std::reverse(composition.offsets.begin(), composition.offsets.end());
        for (auto &offset : composition.offsets)
        {
          offset = span - offset;
        }
        Edges();
        translational = CanonicalTranslationalSignature(edges, R);
        MFEM_VERIFY(translational.chirality > 0,
                    "The reversed cross-section is not the canonical orientation!");
      }
      const bool curved = composition.curved;
      double max_kappa = 0.0;
      std::vector<int> member_chains;
      int own_side = -1;
      for (std::size_t n = 0; n < k; n++)
      {
        const Node &node = composition.nodes[n];
        max_kappa = std::max(max_kappa, WindowedCurvature(ChainOf(node.chain), node.x));
        member_chains.push_back(node.chain);
        if (node.chain == chain && own_side < 0 &&
            quantizer.Less(ChainArcDistance(C, node.x, node.x, mid), neighbourhood))
        {
          own_side = static_cast<int>(n);
        }
      }
      if (own_side < 0)
      {
        if (debug_log)
        {
          input.log("      no own side within the neighbourhood: skipped\n");
        }
        continue;
      }
      max_kappa = std::max(max_kappa, MaxCurvature(C, x0, x1));
      std::string type;
      std::optional<std::string> reason;
      if (k == 2)
      {
        const bool same = runs[composition.nodes[0].run].conductor ==
                          runs[composition.nodes[1].run].conductor;
        const int gap_lower = edges[0].gap_sign, gap_upper = edges[1].gap_sign;
        if (gap_lower > 0 && gap_upper < 0)
        {
          type = same ? "SameConductorGap" : "DifferentConductorGap";
        }
        else if (gap_lower < 0 && gap_upper > 0)
        {
          type = "SameConductorStrip";
        }
        else
        {
          type = "UnclassifiedParallelPair";
          reason = "parallel metal edges within 2R with the gap on the same side "
                   "(overlapping metal in one process plane)";
        }
      }
      else
      {
        type = "ParallelEdgeCluster";
      }
      std::sort(member_chains.begin(), member_chains.end());
      std::string key = type + "|" + translational.signature.dump() + "|" +
                        (curved ? "curved" : "straight") + "|";
      for (const int m : member_chains)
      {
        key += std::to_string(m) + ",";
      }
      const Node &first = composition.nodes.front();
      on_chain.push_back({chain,
                          x0,
                          x1,
                          key,
                          own_side,
                          curved,
                          max_kappa,
                          type,
                          translational.signature,
                          translational.chirality,
                          reason,
                          first.point,
                          {runs[first.run].tangent, composition.lateral, n_ref},
                          static_cast<int>(k),
                          composition.exact});
    }
    // Adjacent elementary intervals of one feature and side merge.
    for (auto &entry : on_chain)
    {
      if (!assigned.empty() && assigned.back().chain == entry.chain &&
          assigned.back().key == entry.key && assigned.back().side == entry.side &&
          std::abs(assigned.back().x1 - entry.x0) <= Tol())
      {
        assigned.back().x1 = entry.x1;
        assigned.back().max_kappa = std::max(assigned.back().max_kappa, entry.max_kappa);
        assigned.back().exact = assigned.back().exact && entry.exact;
      }
      else
      {
        assigned.push_back(std::move(entry));
      }
    }
  }
  if (std::getenv("PALACE_IDENTIFICATION_DEBUG") && input.log)
  {
    std::ostringstream dbg;
    dbg << "  DEBUG component: " << elinks.size() << " elinks\n";
    for (const auto &e : elinks)
    {
      dbg << "    link " << e.chain_a << " - " << e.chain_b << " sep "
          << e.separation[0].value << (e.separation[0].exact ? " (exact)" : " (chord)")
          << " / " << e.separation[1].value
          << (e.separation[1].exact ? " (exact)" : " (chord)") << " pieces "
          << e.pieces[0].size() << " / " << e.pieces[1].size() << "\n";
      for (const int c : {e.chain_a, e.chain_b})
      {
        const Chain &C = ChainOf(c);
        const Run &r0 = runs[C.runs.front()];
        dbg << "      chain " << c << " runs " << C.runs.size() << " length " << C.length
            << " start (" << r0.start[0] << ", " << r0.start[1] << ") curved intervals "
            << C.curved.size() << "\n";
      }
      for (int side = 0; side < 2; side++)
        for (const auto &piece : e.pieces[static_cast<std::size_t>(side)])
        {
          const auto iv = ChainInterval(piece);
          dbg << "      side " << side << " run " << piece.run << " x [" << iv.first << ", "
              << iv.second << "] curved " << piece.curved << "\n";
        }
    }
    for (const auto &a : assigned)
    {
      dbg << "    assigned chain " << a.chain << " [" << a.x0 << ", " << a.x1 << "] side "
          << a.side << " k " << a.sides << " curved " << a.curved << " " << a.type
          << " key " << std::hash<std::string>{}(a.key) % 100000 << "\n";
    }
    input.log(dbg.str());
  }
  StackTimer("features");
  // Features: one per key (in first-appearance order over the chains), RadiusOverR from
  // the tightest windowed radius over the feature's pieces.
  std::map<std::string, std::vector<std::size_t>> groups;
  std::vector<std::string> group_order;
  for (std::size_t i = 0; i < assigned.size(); i++)
  {
    auto [it, inserted] = groups.emplace(assigned[i].key, std::vector<std::size_t>{});
    if (inserted)
    {
      group_order.push_back(assigned[i].key);
    }
    it->second.push_back(i);
  }
  for (const auto &key : group_order)
  {
    const auto &members = groups.at(key);
    const Assigned &lead = assigned[members.front()];
    double length = 0.0, max_kappa = 0.0;
    for (const std::size_t i : members)
    {
      length += assigned[i].x1 - assigned[i].x0;
      max_kappa = std::max(max_kappa, assigned[i].max_kappa);
    }
    if (length <= kSignatureLengthQuantumOverRadius * R)
    {
      continue;  // slivers between cuts, below the signature grid
    }
    nlohmann::json signature = lead.signature;
    std::string type = lead.type;
    if (lead.reason)
    {
      signature["Reason"] = *lead.reason;
    }
    if (lead.curved)
    {
      type = "Curved" + type;
      MFEM_VERIFY(max_kappa > 0.0, "A curved pair / stack without curvature!");
      signature["RadiusOverR"] =
          RoundTo(1.0 / (max_kappa * R), kSignatureLengthQuantumOverRadius);
      // Convexity of the signature's first edge (side 0 for chirality >= 0, the last side
      // for -1: the edge the library model's e1 is placed on): metal inside / outside the
      // bend relative to its gap direction, read on that side's pieces — or as the
      // opposite of the far side's when it carries no curvature of its own (a straight
      // edge facing a bent one within the pair tolerance; concentric edges bend in
      // opposite senses relative to their gaps). The curvature family of that convexity
      // models the pair.
      const int first_side = lead.chirality < 0 ? lead.sides - 1 : 0;
      const int last_side = lead.chirality < 0 ? 0 : lead.sides - 1;
      SignedCurvatureExtremes first, last;
      for (const std::size_t i : members)
      {
        const Assigned &entry = assigned[i];
        if (entry.side == first_side || entry.side == last_side)
        {
          AccumulateSignedCurvature(ChainOf(entry.chain), entry.x0, entry.x1,
                                    entry.side == first_side ? first : last);
        }
      }
      std::string convexity = ConvexityName(first);
      if (convexity.empty())
      {
        const std::string other = ConvexityName(last);
        convexity = other == "Convex" ? "Concave" : other == "Concave" ? "Convex" : other;
      }
      // A curved k > 2 stack whose two outer sides are both exactly straight (the bend on
      // an interior side only) cannot be composed by IdentifyMetalPerimeter: a bent side
      // facing a straight partner at separation d < 2R leaves the pair-constancy band (5 %
      // of d over +-R) exactly where its windowed bend radius r drops below 10R: "curved"
      // needs x > x_b + 0.1 r - R/2 and "constant" needs x < x_b + sqrt(0.1 r d) - R, which
      // together require sqrt(0.1 r d) > 0.1 r + R/2, impossible for d <= 2R (equality at
      // r = 5R, d = 2R); concentric bends make an outer side at least as curved as an
      // interior one. The interior-side fallback that stood here (review J/C4 m5) was
      // therefore unreachable and is replaced by this fail-closed check (drift-block review
      // m9); the GapSide relation it encoded is pinned by the curved 3-edge stack test.
      MFEM_VERIFY(!convexity.empty(),
                  "A curved pair / stack whose outer sides carry no signed curvature "
                  "cannot be composed (an interior-only bend is unreachable through "
                  "IdentifyMetalPerimeter)!");
      signature["Convexity"] = convexity;
    }
    const int feature = NewFeature(type, signature, lead.chirality);
    features[feature].origin = lead.origin;
    features[feature].axes = lead.axes;
    feature_sides[feature] = lead.sides;
    features[feature].exact_parameters = std::all_of(
        members.begin(), members.end(), [&](std::size_t i) { return assigned[i].exact; });
    for (const std::size_t i : members)
    {
      const Assigned &entry = assigned[i];
      const Chain &C = ChainOf(entry.chain);
      // Runs overlapping [x0, x1] (run_offset ascending).
      std::size_t kr =
          std::min(static_cast<std::size_t>(std::upper_bound(C.run_offset.begin(),
                                                             C.run_offset.end(), entry.x0) -
                                            C.run_offset.begin()) -
                       1,
                   C.runs.size() - 1);
      for (; kr < C.runs.size() && C.run_offset[kr] < entry.x1 - Tol(); kr++)
      {
        const double offset = C.run_offset[kr];
        const Run &run = runs[C.runs[kr]];
        const double s0 = std::max(entry.x0, offset) - offset;
        const double s1 = std::min(entry.x1, offset + run.length) - offset;
        if (s1 - s0 > Tol())
        {
          claims[C.runs[kr]].push_back({feature, 2, {s0, s1}, entry.side});
        }
      }
    }
  }
  StackTimer("end");
  stack_geometric_offsets += fallbacks;
  stack_composition_cap_hits += cap_hits;
  if (fallbacks > 0 && input.log)
  {
    input.log("  Identification stacks: " + std::to_string(fallbacks) +
              " elementary intervals with a geometric offset (no consecutive link)\n");
  }
  if (cap_hits > 0 && input.log)
  {
    input.log("  Identification stacks: " + std::to_string(cap_hits) +
              " elementary intervals reached the composition cap of " +
              std::to_string(kStackCompositionCap) + " members\n");
  }
}

// Window of a vertex site along one incident run and, for curved chains, the following
// runs: R of arc length, shortened to half the chain distance to the next feature vertex
// when that distance is below 2R (two corners R apart on a strip end share the end edge).
std::vector<Interval>
Identifier::ChainWindow(std::size_t run, std::size_t from_vertex, double length,
                        std::vector<std::pair<std::size_t, Interval>> &out) const
{
  const Chain &chain = chains.at(chain_index.at(runs[run].chain));
  const std::size_t start_index = RunIndexInChain(chain, run);
  const bool forward = runs[run].start_vertex == from_vertex;
  const std::size_t m = chain.runs.size();
  // Distance along the chain to the next feature vertex (or chain end).
  double distance = 0.0;
  bool found_feature = false;
  {
    std::size_t k = start_index;
    for (std::size_t steps = 0; steps < m; steps++)
    {
      const Run &r = runs[chain.runs[k]];
      if (r.excluded)
      {
        break;
      }
      distance += r.length;
      const std::size_t far = forward ? r.end_vertex : r.start_vertex;
      if (IsFeatureVertex(far))
      {
        found_feature = true;
        break;
      }
      if (forward ? k + 1 >= m : k == 0)
      {
        if (!chain.closed)
        {
          break;
        }
        k = forward ? 0 : m - 1;
      }
      else
      {
        k = forward ? k + 1 : k - 1;
      }
    }
  }
  double limit = length;
  if (found_feature && quantizer.Less(distance, 2.0 * length))
  {
    limit = 0.5 * distance;
  }
  std::vector<Interval> result;
  double remaining = limit;
  std::size_t k = start_index;
  for (std::size_t steps = 0; steps < m && remaining > Tol(); steps++)
  {
    const Run &r = runs[chain.runs[k]];
    if (r.excluded)
    {
      break;
    }
    const double take = std::min(remaining, r.length);
    const Interval interval =
        forward ? Interval{0.0, take} : Interval{r.length - take, r.length};
    out.emplace_back(chain.runs[k], interval);
    result.push_back(interval);
    remaining -= take;
    if (forward ? k + 1 >= m : k == 0)
    {
      if (!chain.closed)
      {
        break;
      }
      k = forward ? 0 : m - 1;
    }
    else
    {
      k = forward ? k + 1 : k - 1;
    }
  }
  return result;
}

void Identifier::BuildVertexWindows()
{
  for (auto &site : sites)
  {
    if (!site.vertex)
    {
      continue;  // rounded corners carry their window from detection
    }
    for (const std::size_t r : site.runs_at_site)
    {
      ChainWindow(r, *site.vertex, kVertexWindowOverRadius * R, site.window);
    }
  }
}

// Events, cores, clusters and their claimed portions.
void Identifier::BuildClusters()
{
  const double interaction = kInteractionDistanceOverRadius * R;
  const double zone = kThroughVertexZoneOverRadius * R;
  const double ball = kClusterBallOverRadius * R;
  const double join = kVertexJoinsClusterOverRadius * R;

  std::vector<EventCore> cores;
  // self_partner (a chain facing itself): the part of run b at least pi R of arc length
  // from the chain position x of a point of run a; the piece of a is then subdivided into
  // sub-pieces of at most R / 2 and each faces the partner region of its midpoint (the
  // former run-level zone excluded a whole run within pi R of the other run's nearest end,
  // which hid the fold end of a hairpin facing the far part of a long leg).
  auto CoresOnRun =
      [&](std::size_t a, std::size_t b, const std::vector<std::vector<Interval>> &zones_a,
          const std::vector<std::vector<Interval>> &zones_b,
          const std::function<std::vector<Interval>(double)> *self_partner = nullptr)
  {
    const Run &ra = runs[a];
    const Run &rb = runs[b];
    std::vector<double> points = {0.0, ra.length};
    for (const auto &zone_intervals : zones_a)
    {
      for (const auto &interval : zone_intervals)
      {
        points.push_back(std::clamp(interval.first, 0.0, ra.length));
        points.push_back(std::clamp(interval.second, 0.0, ra.length));
      }
    }
    std::sort(points.begin(), points.end());
    std::vector<Interval> found;
    for (std::size_t k = 0; k + 1 < points.size(); k++)
    {
      const double lo = points[k], hi = points[k + 1];
      if (hi - lo <= Tol())
      {
        continue;
      }
      // Through-vertex pairs have at least one point inside a shared vertex's 2R zone: a
      // piece of run a inside a zone has no events at all, and a piece outside pairs only
      // with the part of run b outside every zone.
      const double mid = 0.5 * (lo + hi);
      const bool inside_zone =
          std::any_of(zones_a.begin(), zones_a.end(),
                      [&](const std::vector<Interval> &zone)
                      {
                        return std::any_of(zone.begin(), zone.end(), [&](const Interval &i)
                                           { return i.first <= mid && mid <= i.second; });
                      });
      if (inside_zone)
      {
        continue;
      }
      std::vector<Interval> excluded;
      for (const auto &zone : zones_b)
      {
        excluded.insert(excluded.end(), zone.begin(), zone.end());
      }
      const auto complement = SubtractIntervals({Interval{0.0, rb.length}},
                                                MergeIntervals(excluded, Tol()), Tol());
      const std::size_t sub_pieces =
          self_partner ? std::max<std::size_t>(
                             1, static_cast<std::size_t>(std::ceil((hi - lo) / (0.5 * R))))
                       : 1;
      for (std::size_t q = 0; q < sub_pieces; q++)
      {
        const double sub_lo = lo + (hi - lo) * q / sub_pieces;
        const double sub_hi = lo + (hi - lo) * (q + 1) / sub_pieces;
        std::vector<Interval> partner = complement;
        if (self_partner)
        {
          partner = IntersectIntervals(complement, (*self_partner)(0.5 * (sub_lo + sub_hi)),
                                       Tol());
        }
        for (const auto &piece : partner)
        {
          if (run_arcs[a] || run_arcs[b])
          {
            // Option A: the event set on the fitted arcs (the sampled sublevel solver; the
            // distance along an arc is not convex).
            const CurvePiece target = RunPiece(b, piece.first, piece.second);
            auto f = [&](double s) { return PointPieceDistance(RunPoint(a, s), target); };
            const double scale = run_arcs[a] && ra.length > 0.0
                                     ? run_arcs[a]->piece.Length() / ra.length
                                     : 1.0;
            double lo_a = sub_lo, hi_a = sub_hi;
            if (!run_arcs[a])
            {
              // A straight run against an arc piece: bracket by the convex sublevel set of
              // the distance to the piece's chord below 2R + sagitta (a long ground edge is
              // not sampled along its whole length for every fillet chord near it).
              auto chord = [&](double s)
              { return PointSegmentDistance(ra.At(s), target.a, target.b); };
              const auto bracket = ConvexSublevelInterval(
                  chord, sub_lo, sub_hi, interaction + target.Sagitta() + 2.0 * Tol(),
                  quantizer);
              if (!bracket)
              {
                continue;
              }
              lo_a = bracket->first;
              hi_a = bracket->second;
            }
            for (const auto &interval :
                 SampledSublevelIntervals(f, lo_a, hi_a, interaction, quantizer,
                                          ArcSampleStep() / std::max(scale, 1.0e-300)))
            {
              found.push_back(interval);
            }
            continue;
          }
          const Point3D q0 = rb.At(piece.first), q1 = rb.At(piece.second);
          auto f = [&](double s) { return PointSegmentDistance(ra.At(s), q0, q1); };
          if (auto interval =
                  ConvexSublevelInterval(f, sub_lo, sub_hi, interaction, quantizer))
          {
            found.push_back(*interval);
          }
        }
      }
    }
    for (const auto &interval : MergeIntervals(found, Tol()))
    {
      const CurvePiece piece = RunPiece(a, interval.first, interval.second);
      cores.push_back({a, interval, piece.a, piece.b, run_plane[a], piece});
    }
  };

  std::size_t run_pairs_examined = 0;
  for (std::size_t a = 0; a < runs.size(); a++)
  {
    if (runs[a].excluded)
    {
      continue;
    }
    stage.Progress(a, runs.size());
    // Candidates b > a whose boxes lie within the interaction distance (sorted: the former
    // order of the loop over every b).
    for (const std::size_t b : RunsNearRun(
             a, interaction + 2.0 * Tol() + WholeRunPiece(a).Sagitta() + max_run_sagitta))
    {
      if (b <= a || runs[b].excluded || run_plane[a] != run_plane[b])
      {
        continue;  // no events between two metal planes (decision 82(1))
      }
      const Chain &ca = chains[chain_index.at(runs[a].chain)];
      const Chain &cb = chains[chain_index.at(runs[b].chain)];
      // Two runs of ONE chain (a chain facing itself, decision 82(2) addition): no events
      // within the self-pair neighbourhood of pi R along the chain (the local neighbourhood
      // of a bend); beyond it the runs interact like two chains — a strip flaring into a
      // pad diverges (not constant) and is a cluster, as it would be for two chains. The
      // zone of run a is the part within pi R of arc length of run b (its ends), in place
      // of the shared-vertex zones (every vertex of a chain is shared with itself).
      const bool self = runs[a].chain == runs[b].chain;
      std::function<std::vector<Interval>(double)> self_partner_on_b, self_partner_on_a;
      if (self)
      {
        if (ca.Rigid())
        {
          continue;  // one run: no other run of the chain
        }
        const double neighbourhood = kSelfPairNeighbourhoodOverRadius * R;
        const double a0 = ca.run_offset[runs[a].index_in_chain], a1 = a0 + runs[a].length;
        const double b0 = ca.run_offset[runs[b].index_in_chain], b1 = b0 + runs[b].length;
        double farthest = std::max(std::abs(b1 - a0), std::abs(a1 - b0));
        if (ca.closed)
        {
          farthest = std::max(std::min(std::abs(b1 - a0), ca.length - std::abs(b1 - a0)),
                              std::min(std::abs(a1 - b0), ca.length - std::abs(a1 - b0)));
        }
        if (quantizer.Less(farthest, neighbourhood))
        {
          continue;
        }
        // The partner region on run [y0, y1] of a point at run parameter s of run [x0, x1]:
        // the run's points at least pi R of arc length away (the shorter way round on a
        // closed chain), as run parameters.
        auto Partner = [&](double x0, double y0, double y1)
        {
          return [=, &ca](double s)
          {
            const double x = x0 + s;
            std::vector<Interval> near;
            for (const double shift : ca.closed
                                          ? std::vector<double>{-ca.length, 0.0, ca.length}
                                          : std::vector<double>{0.0})
            {
              const double lo = std::max(y0, x + shift - neighbourhood) - y0;
              const double hi = std::min(y1, x + shift + neighbourhood) - y0;
              if (hi - lo > Tol())
              {
                near.emplace_back(lo, hi);
              }
            }
            return SubtractIntervals({Interval{0.0, y1 - y0}}, MergeIntervals(near, Tol()),
                                     Tol());
          };
        };
        self_partner_on_b = Partner(a0, b0, b1);
        self_partner_on_a = Partner(b0, a0, a1);
      }
      run_pairs_examined++;
      if (ParallelRigidRuns(a, b))
      {
        continue;  // parallel straight runs (one class): a translational interaction
      }
      if (!quantizer.Less(PieceDistance(WholeRunPiece(a), WholeRunPiece(b)), interaction))
      {
        continue;
      }
      std::vector<std::vector<Interval>> zones_a =
                                             self ? std::vector<std::vector<Interval>>{}
                                                  : ThroughZones(a, ca, cb, true),
                                         zones_b =
                                             self ? std::vector<std::vector<Interval>>{}
                                                  : ThroughZones(b, cb, ca, true);
      // Portions of a constant-separation pair along a bend with the other chain (a pair
      // feature, or a non-interacting pair at or beyond 2R) are not events.
      std::vector<Interval> bent_a, bent_b;
      for (const auto &[other, interval] : bent_claims[a])
      {
        if (other == runs[b].chain)
        {
          bent_a.push_back(interval);
        }
      }
      for (const auto &[other, interval] : bent_claims[b])
      {
        if (other == runs[a].chain)
        {
          bent_b.push_back(interval);
        }
      }
      if (!bent_a.empty())
      {
        zones_a.push_back(MergeIntervals(bent_a, Tol()));
      }
      if (!bent_b.empty())
      {
        zones_b.push_back(MergeIntervals(bent_b, Tol()));
      }
      CoresOnRun(a, b, zones_a, zones_b, self ? &self_partner_on_b : nullptr);
      CoresOnRun(b, a, zones_b, zones_a, self ? &self_partner_on_a : nullptr);
    }
  }

  const std::size_t run_cores = cores.size();
  // Two vertex sites within 2R have overlapping radius-R windows (invariant A2): they are
  // an interaction event of their own, with the two points as degenerate cores.
  Point3D sites_lo{}, sites_hi{};
  for (std::size_t i = 0; i < sites.size(); i++)
  {
    for (int d = 0; d < 3; d++)
    {
      sites_lo[d] = i == 0 ? sites[i].point[d] : std::min(sites_lo[d], sites[i].point[d]);
      sites_hi[d] = i == 0 ? sites[i].point[d] : std::max(sites_hi[d], sites[i].point[d]);
    }
  }
  UniformGrid site_grid(4.0 * R, sites_lo);
  for (std::size_t i = 0; i < sites.size(); i++)
  {
    site_grid.Insert(i, sites[i].point, sites[i].point);
  }
  for (std::size_t i = 0; i < sites.size(); i++)
  {
    for (const std::size_t j :
         site_grid.Query(sites[i].point, sites[i].point, interaction + 2.0 * Tol()))
    {
      if (j <= i || SitePlane(sites[i]) != SitePlane(sites[j]))
      {
        continue;
      }
      if (quantizer.Less(Distance(sites[i].point, sites[j].point), interaction))
      {
        cores.push_back({std::numeric_limits<std::size_t>::max(),
                         {0.0, 0.0},
                         sites[i].point,
                         sites[i].point,
                         SitePlane(sites[i]),
                         CurvePiece{sites[i].point, sites[i].point, std::nullopt}});
        cores.push_back({std::numeric_limits<std::size_t>::max(),
                         {0.0, 0.0},
                         sites[j].point,
                         sites[j].point,
                         SitePlane(sites[j]),
                         CurvePiece{sites[j].point, sites[j].point, std::nullopt}});
      }
    }
  }

  // Connected union of the radius-R balls: cores whose distance is below 2R (candidates
  // from the cores' boxes, united in the former (i, j) order).
  Point3D cores_lo{}, cores_hi{};
  for (std::size_t i = 0; i < cores.size(); i++)
  {
    Point3D lo, hi;
    PieceBoundingBox(cores[i].piece, lo, hi);
    for (int d = 0; d < 3; d++)
    {
      cores_lo[d] = i == 0 ? lo[d] : std::min(cores_lo[d], lo[d]);
      cores_hi[d] = i == 0 ? hi[d] : std::max(cores_hi[d], hi[d]);
    }
  }
  UniformGrid core_grid(4.0 * R, cores_lo);
  for (std::size_t i = 0; i < cores.size(); i++)
  {
    Point3D lo, hi;
    PieceBoundingBox(cores[i].piece, lo, hi);
    core_grid.Insert(i, lo, hi);
  }
  UnionFind uf(cores.size());
  for (std::size_t i = 0; i < cores.size(); i++)
  {
    Point3D lo, hi;
    PieceBoundingBox(cores[i].piece, lo, hi);
    for (const std::size_t j : core_grid.Query(lo, hi, interaction + 2.0 * Tol()))
    {
      if (j <= i || cores[i].plane != cores[j].plane)
      {
        continue;
      }
      // Cores already in one region need no distance: the union-find partition is a set
      // function of the "< 2R" relation, and a pair already connected cannot change it
      // (unit-test profile 2026-09-28: the sampled / golden arc-arc PieceDistance of every
      // grid pair was 98.8 % of the identification's time on the high-order island meshes;
      // stage counts and manifests are identical with the guard).
      if (uf.Find(i) == uf.Find(j))
      {
        continue;
      }
      if (quantizer.Less(PieceDistance(cores[i].piece, cores[j].piece), interaction))
      {
        uf.Union(i, j);
      }
    }
  }
  // A vertex site within 2R of an event core of its own plane (its radius-R window
  // overlaps the core's radius-R ball) joins that region (and merges the regions it
  // reaches).
  std::vector<std::optional<std::size_t>> site_root(sites.size());
  for (std::size_t s = 0; s < sites.size(); s++)
  {
    std::optional<std::size_t> root;
    for (const std::size_t i :
         core_grid.Query(sites[s].point, sites[s].point, join + 2.0 * Tol()))
    {
      if (cores[i].plane == SitePlane(sites[s]) &&
          quantizer.Less(PointPieceDistance(sites[s].point, cores[i].piece), join))
      {
        if (root)
        {
          uf.Union(*root, i);
        }
        root = uf.Find(i);
      }
    }
    site_root[s] = root;
  }
  std::map<std::size_t, std::size_t> cluster_of_root;
  for (std::size_t i = 0; i < cores.size(); i++)
  {
    const std::size_t root = uf.Find(i);
    const auto [it, inserted] = cluster_of_root.emplace(root, cluster_of_root.size());
    (void)inserted;
    if (cluster_cores.size() <= it->second)
    {
      cluster_cores.resize(it->second + 1);
      cluster_sites.resize(it->second + 1);
    }
    cluster_cores[it->second].push_back(cores[i]);
  }
  for (std::size_t s = 0; s < sites.size(); s++)
  {
    if (site_root[s])
    {
      const std::size_t c = cluster_of_root.at(uf.Find(*site_root[s]));
      sites[s].cluster = static_cast<int>(c);
      cluster_sites[c].push_back(s);
    }
  }

  // Claimed portions: run intervals within R of a core, plus the member sites' windows. The
  // candidate runs of a cluster are the runs whose boxes lie within the ball of one of its
  // cores and the runs of its sites' windows, in run order (the former loop over every
  // run). The features are emitted by EmitClusters after the extension (decision 85(2)).
  std::size_t largest_cluster_edges = 0;
  cluster_claimed.assign(cluster_cores.size(), {});
  for (std::size_t c = 0; c < cluster_cores.size(); c++)
  {
    stage.Progress(c, cluster_cores.size());
    std::vector<std::size_t> candidate_runs;
    for (const auto &core : cluster_cores[c])
    {
      const auto near = RunsNearPiece(core.piece, ball + 2.0 * Tol());
      candidate_runs.insert(candidate_runs.end(), near.begin(), near.end());
    }
    for (const std::size_t s : cluster_sites[c])
    {
      for (const auto &[wr, interval] : sites[s].window)
      {
        (void)interval;
        candidate_runs.push_back(wr);
      }
    }
    std::sort(candidate_runs.begin(), candidate_runs.end());
    candidate_runs.erase(std::unique(candidate_runs.begin(), candidate_runs.end()),
                         candidate_runs.end());
    for (const std::size_t r : candidate_runs)
    {
      if (runs[r].excluded)
      {
        continue;
      }
      std::vector<Interval> intervals;
      const CurvePiece run_piece = WholeRunPiece(r);
      for (const auto &core : cluster_cores[c])
      {
        if (!quantizer.Less(PieceDistance(run_piece, core.piece), ball))
        {
          continue;
        }
        const auto within = RunIntervalWithinPiece(r, core.piece, ball);
        intervals.insert(intervals.end(), within.begin(), within.end());
      }
      for (const std::size_t s : cluster_sites[c])
      {
        for (const auto &[wr, interval] : sites[s].window)
        {
          if (wr == r)
          {
            intervals.push_back(interval);
          }
        }
      }
      for (const auto &interval : MergeIntervals(intervals, Tol()))
      {
        cluster_claimed[c].emplace_back(r, interval);
      }
    }
    MFEM_VERIFY(!cluster_claimed[c].empty(), "A spatial cluster claims no perimeter!");
    largest_cluster_edges = std::max(largest_cluster_edges, cluster_claimed[c].size());
  }
  {
    std::ostringstream counts;
    counts << run_pairs_examined << " run pairs within reach, " << run_cores
           << " event cores on runs + " << (cores.size() - run_cores) << " site cores, "
           << cluster_cores.size() << " clusters (largest " << largest_cluster_edges
           << " edges before the extension)";
    stage.End(counts.str());
  }
  if (std::getenv("PALACE_IDENTIFICATION_DEBUG_EXTENSION") && input.log)
  {
    std::ostringstream dbg;
    dbg << std::setprecision(10);
    for (std::size_t c = 0; c < cluster_claimed.size(); c++)
    {
      for (const auto &core : cluster_cores[c])
      {
        dbg << "    base core cluster " << c << " run " << core.run << " from ("
            << core.p0[0] << ", " << core.p0[1] << ") to (" << core.p1[0] << ", "
            << core.p1[1] << ")\n";
      }
      for (const auto &[r, interval] : cluster_claimed[c])
      {
        const Run &run = runs[r];
        dbg << "    base claim cluster " << c << " run " << r << " chain " << run.chain
            << " s [" << interval.first << ", " << interval.second << "] of " << run.length
            << " from (" << RunPoint(r, interval.first)[0] << ", "
            << RunPoint(r, interval.first)[1] << ") to (" << RunPoint(r, interval.second)[0]
            << ", " << RunPoint(r, interval.second)[1] << ")\n";
      }
    }
    input.log(dbg.str());
  }
}

// Cluster extension (decision 85(2)). Every single-edge portion — a run interval outside
// every cluster portion, vertex window, pair / stack claim and CrossLayer zone, i.e. what
// would become an IsolatedEdge / CurvedEdge — whose 3D distance to a cluster's claimed
// perimeter on another chain (or on its own chain beyond the self-pair neighbourhood) is
// strictly below 2R, outside the through-vertex zones of the two chains, joins that cluster
// over the sub-interval within 2R; a single-edge portion within 2R of the window of a
// vertex feature outside every cluster makes that vertex a cluster together with the
// portion. A portion within 2R of several clusters (or vertex features) merges them
// (union-find: the result does not depend on the order). Pairs and stacks are joint
// descriptions and are never absorbed; the stack-end recomposition then runs again on the
// enlarged claims, and the extension iterates to closure over the single-edge portions
// only. Returns the absorbed length of this pass (the caller stops when it is at most the
// signature parameter tolerance 1e-3 R: that last pass IS applied without a further
// recomposition — its pieces come from the unclaimed remainder, so no pair / stack claim
// overlaps them and the partition holds; the stack cut images on the other members are then
// consistent with the final claims to within the tolerance, below the resolution of every
// parameter; NOT applying it would leave sub-tolerance unclaimed slivers as isolated
// features facing the cluster metal, e.g. a 0.24 nm IsolatedEdge on DS-SCT-001 — review
// fix-3 m8, stated). With measure_joint_claims the SAME across rule is evaluated on the
// pair / stack claims (priority 2) instead of the single-edge remainder and nothing is
// absorbed: the returned length is Diagnostics.StackEndThirdBodyLength (decision 85(2)),
// the one definition the facing gate's StackEndRecomposition / StackEndThirdBody classes
// sample (review m5).
double Identifier::ExtendClusters(bool measure_joint_claims)
{
  const double interaction = kInteractionDistanceOverRadius * R;
  const double neighbourhood = kSelfPairNeighbourhoodOverRadius * R;
  const double sliver = kSignatureLengthQuantumOverRadius * R;
  // Taken intervals per run: cluster portions, non-cluster windows, pair / stack claims,
  // CrossLayer zones. The remainder of every run is the single-edge candidate set.
  std::vector<std::vector<Interval>> taken(runs.size());
  for (std::size_t c = 0; c < cluster_claimed.size(); c++)
  {
    for (const auto &[r, interval] : cluster_claimed[c])
    {
      taken[r].push_back(interval);
    }
  }
  std::vector<std::size_t> free_sites;
  for (std::size_t s = 0; s < sites.size(); s++)
  {
    if (sites[s].cluster < 0)
    {
      free_sites.push_back(s);
      for (const auto &[r, interval] : sites[s].window)
      {
        taken[r].push_back(interval);
      }
    }
  }
  for (std::size_t r = 0; r < runs.size(); r++)
  {
    for (const auto &claim : claims[r])
    {
      taken[r].push_back(claim.interval);
    }
    for (const auto &zone : cross_layer[r])
    {
      taken[r].push_back(zone);
    }
  }
  // Claimed pieces of the clusters and the windows of the free sites, indexed by their
  // boxes: piece -> (cluster index, or cluster_claimed.size() + free site index), with the
  // piece's chain interval and the projection extension at its ends (below).
  struct Piece
  {
    std::size_t run;
    Interval interval;
    std::size_t owner;            // cluster c, or cluster_claimed.size() + free site index
    double extend_lo, extend_hi;  // projection-domain extension at the piece ends
  };
  std::vector<Piece> pieces;
  // A single-edge point joins the owner where its PERPENDICULAR projection onto a claimed
  // piece falls inside the piece and its distance to the piece is below 2R: the claimed
  // perimeter is faced across, never reached diagonally past its end (a diagonal reach
  // would creep along a sub-2R strip: the partner beyond the claim end pairs with the
  // free continuation instead). Inside a claimed chain interval the perpendicular domains
  // of consecutive pieces leave a wedge at every joint; each piece's domain is extended
  // there by 2R tan(turn) (the wedge's width at 2R; never at the end of a claimed
  // interval).
  auto AddPieces =
      [&](std::size_t owner, const std::vector<std::pair<std::size_t, Interval>> &claimed)
  {
    std::map<int, std::vector<Interval>> by_chain;
    for (const auto &[r, interval] : claimed)
    {
      const Chain &C = chains[chain_index.at(runs[r].chain)];
      const double offset = C.run_offset[runs[r].index_in_chain];
      by_chain[C.id].emplace_back(offset + interval.first, offset + interval.second);
    }
    for (auto &[chain, list] : by_chain)
    {
      (void)chain;
      list = MergeIntervals(std::move(list), Tol());
    }
    for (const auto &[r, interval] : claimed)
    {
      const Chain &C = chains[chain_index.at(runs[r].chain)];
      const std::size_t k = runs[r].index_in_chain;
      const double offset = C.run_offset[k];
      const double x0 = offset + interval.first, x1 = offset + interval.second;
      const auto &list = by_chain.at(C.id);
      const Interval *whole = nullptr;
      for (const auto &candidate : list)
      {
        if (x0 >= candidate.first - Tol() && x1 <= candidate.second + Tol())
        {
          whole = &candidate;
          break;
        }
      }
      auto Extension = [&](bool at_start)
      {
        const double x = at_start ? x0 : x1;
        if (!whole || x <= whole->first + Tol() || x >= whole->second - Tol())
        {
          return 0.0;  // the end of the claimed chain interval: no reach past it
        }
        // The joint at this piece end: the run boundary (joint_turn[k] before run k,
        // joint_turn[k + 1] after; wrapping on a closed chain), 0 inside a run. At a joint
        // of an arc run the turn is read between the arc tangents (option A): consecutive
        // chords of one circle and a chord meeting its tangent arm turn by nothing on the
        // fitted geometry, so the wedge does not depend on the chord count.
        const double s_end = at_start ? interval.first : interval.second;
        double turn = 0.0;
        if (at_start && s_end <= Tol())
        {
          turn = C.joint_turn.empty() ? 0.0 : GeometricJointTurn(C, k);
        }
        else if (!at_start && s_end >= runs[r].length - Tol())
        {
          if (k + 1 < C.runs.size())
          {
            turn = GeometricJointTurn(C, k + 1);
          }
          else if (C.closed && !C.joint_turn.empty())
          {
            turn = GeometricJointTurn(C, 0);
          }
        }
        turn = std::min(turn, kClusterExtensionWedgeCapDegrees * std::acos(-1.0) / 180.0);
        return interaction * std::tan(turn);
      };
      pieces.push_back({r, interval, owner, Extension(true), Extension(false)});
    }
  };
  for (std::size_t c = 0; c < cluster_claimed.size(); c++)
  {
    AddPieces(c, cluster_claimed[c]);
  }
  for (std::size_t k = 0; k < free_sites.size(); k++)
  {
    AddPieces(cluster_claimed.size() + k, sites[free_sites[k]].window);
  }
  if (pieces.empty())
  {
    return 0.0;
  }
  Point3D lo{}, hi{};
  for (std::size_t i = 0; i < pieces.size(); i++)
  {
    Point3D plo, phi;
    PieceBoundingBox(
        RunPiece(pieces[i].run, pieces[i].interval.first, pieces[i].interval.second), plo,
        phi);
    for (int d = 0; d < 3; d++)
    {
      lo[d] = i == 0 ? plo[d] : std::min(lo[d], plo[d]);
      hi[d] = i == 0 ? phi[d] : std::max(hi[d], phi[d]);
    }
  }
  UniformGrid piece_grid(4.0 * R, lo);
  for (std::size_t i = 0; i < pieces.size(); i++)
  {
    Point3D plo, phi;
    PieceBoundingBox(
        RunPiece(pieces[i].run, pieces[i].interval.first, pieces[i].interval.second), plo,
        phi);
    piece_grid.Insert(i, plo, phi);
  }
  // Site of every feature vertex (the through-vertex zones of a cluster's own member
  // vertices do not exclude: the cluster describes them).
  std::map<std::size_t, std::size_t> site_of_vertex;
  for (std::size_t i = 0; i < sites.size(); i++)
  {
    if (sites[i].vertex)
    {
      site_of_vertex.try_emplace(*sites[i].vertex, i);
    }
  }
  auto ThroughZonesExcept = [&](std::size_t a, const Chain &A, const Chain &B, int cluster)
  {
    std::vector<std::vector<Interval>> zones;
    for (const std::size_t v : SharedVertices(A, B))
    {
      const auto sv = site_of_vertex.find(v);
      if (cluster >= 0 && sv != site_of_vertex.end() &&
          sites[sv->second].cluster == cluster)
      {
        continue;
      }
      const Point3D &p = input.vertices[v].coordinate;
      zones.push_back(RunIntervalWithinPiece(a, CurvePiece{p, p, std::nullopt},
                                             kThroughVertexZoneOverRadius * R));
    }
    const auto it =
        through_arc.find(std::make_pair(std::min(A.id, B.id), std::max(A.id, B.id)));
    if (it != through_arc.end())
    {
      for (const int arc : it->second)
      {
        const Arc &fitted = arcs[static_cast<std::size_t>(arc)];
        if (cluster >= 0 && fitted.feature >= 0 &&
            sites[static_cast<std::size_t>(fitted.feature)].cluster == cluster)
        {
          continue;
        }
        const auto ci = chain_index.find(fitted.chain);
        if (ci == chain_index.end())
        {
          continue;
        }
        std::vector<Interval> zone;
        for (const std::size_t r : chains[ci->second].runs)
        {
          const auto within =
              RunIntervalWithinPiece(a, WholeRunPiece(r), kThroughVertexZoneOverRadius * R);
          zone.insert(zone.end(), within.begin(), within.end());
        }
        zones.push_back(MergeIntervals(std::move(zone), Tol()));
      }
    }
    return zones;
  };
  // Absorptions of this pass: (run, interval, owners) computed from the state before the
  // pass; applied afterwards.
  struct Absorption
  {
    std::size_t run;
    Interval interval;
    std::vector<std::size_t> owners;
    bool translational = false;  // a pair / stack stretch (decision 224), else single-edge
  };
  std::vector<Absorption> absorptions;
  for (std::size_t r = 0; r < runs.size(); r++)
  {
    if (runs[r].excluded)
    {
      continue;
    }
    stage.Progress(r, runs.size(),
                   measure_joint_claims ? "stack-end third body" : "cluster extension");
    std::vector<Interval> remainder;
    if (measure_joint_claims)
    {
      for (const auto &claim : claims[r])
      {
        if (claim.priority == 2)
        {
          remainder.push_back(claim.interval);
        }
      }
      remainder = MergeIntervals(std::move(remainder), Tol());
    }
    else
    {
      remainder = SubtractIntervals({Interval{0.0, runs[r].length}},
                                    MergeIntervals(taken[r], Tol()), Tol());
    }
    if (remainder.empty())
    {
      continue;
    }
    const Chain &A = chains[chain_index.at(runs[r].chain)];
    const double x_a0 = A.run_offset[runs[r].index_in_chain];
    Point3D rlo, rhi;
    BoundingBox(runs[r].start, runs[r].end, rlo, rhi);
    std::map<std::size_t, std::vector<Interval>> within_by_owner;
    std::map<std::pair<int, int>, std::vector<std::vector<Interval>>> zones_by_chain;
    for (const std::size_t i : piece_grid.Query(rlo, rhi, interaction + 2.0 * Tol()))
    {
      const Piece &piece = pieces[i];
      const Run &rp = runs[piece.run];
      if (rp.excluded || piece.run == r || run_plane[piece.run] != run_plane[r])
      {
        continue;
      }
      const Chain &B = chains[chain_index.at(rp.chain)];
      const int owner_cluster =
          piece.owner < cluster_claimed.size() ? static_cast<int>(piece.owner) : -1;
      // The part of the claimed piece outside the through-vertex zones of the shared
      // vertices that are not the owner's members (an interaction has BOTH points outside
      // the zones, as for the events: the window of a free corner lies inside its own zone
      // and never absorbs its arms; a cluster's member vertices are the cluster's).
      std::vector<Interval> piece_parts = {piece.interval};
      if (rp.chain != runs[r].chain)
      {
        std::vector<Interval> excluded;
        for (const auto &zone : ThroughZonesExcept(piece.run, B, A, owner_cluster))
        {
          excluded.insert(excluded.end(), zone.begin(), zone.end());
        }
        piece_parts =
            SubtractIntervals(piece_parts, MergeIntervals(excluded, Tol()), Tol());
      }
      std::vector<Interval> found;
      for (const auto &part : piece_parts)
      {
        const Point3D q0 = rp.At(part.first), q1 = rp.At(part.second);
        const double lo_p =
            part.first -
            (std::abs(part.first - piece.interval.first) <= Tol() ? piece.extend_lo : 0.0);
        const double hi_p =
            part.second + (std::abs(part.second - piece.interval.second) <= Tol()
                               ? piece.extend_hi
                               : 0.0);
        if (run_arcs[r] || run_arcs[piece.run])
        {
          // Option A: the claimed part and / or the candidate run on their fitted arcs. The
          // across domain of an arc part is the radial projection inside its angular range
          // (extended by the wedge as an angle); a straight part keeps the perpendicular
          // projection of the candidate's arc points onto its line. The sublevel set of the
          // distance restricted to the domain is solved on the samples.
          const CurvePiece target = RunPiece(piece.run, part.first, part.second);
          if (!quantizer.Less(PieceDistance(WholeRunPiece(r), target), interaction))
          {
            continue;
          }
          std::function<bool(const Point3D &)> in_domain;
          if (target.arc)
          {
            const ArcPiece &arc = *target.arc;
            const double angle_lo =
                std::min(arc.theta0, arc.theta1) - (part.first - lo_p) / arc.radius;
            const double angle_hi =
                std::max(arc.theta0, arc.theta1) + (hi_p - part.second) / arc.radius;
            in_domain = [&arc, angle_lo, angle_hi](const Point3D &p)
            {
              const Point3D d = Sub(p, arc.center);
              if (std::hypot(Dot(d, arc.u), Dot(d, arc.v)) <= 0.0)
              {
                return false;
              }
              const double theta = arc.AngleOf(p);
              return theta >= angle_lo && theta <= angle_hi;
            };
          }
          else
          {
            in_domain = [&rp, lo_p, hi_p](const Point3D &p)
            {
              const double g = Dot(Sub(p, rp.start), rp.tangent);
              return g >= lo_p && g <= hi_p;
            };
          }
          auto f = [&](double s)
          {
            const Point3D p = RunPoint(r, s);
            return in_domain(p) ? PointPieceDistance(p, target)
                                : std::numeric_limits<double>::infinity();
          };
          const double scale = run_arcs[r] && runs[r].length > 0.0
                                   ? run_arcs[r]->piece.Length() / runs[r].length
                                   : 1.0;
          for (const auto &bracket : RunIntervalWithinPiece(r, target, interaction))
          {
            const auto within = SampledSublevelIntervals(
                f, bracket.first, bracket.second, interaction, quantizer,
                ArcSampleStep() / std::max(scale, 1.0e-300));
            found.insert(found.end(), within.begin(), within.end());
          }
        }
        else
        {
          if (!quantizer.Less(SegmentSegmentDistance(runs[r].start, runs[r].end, q0, q1),
                              interaction))
          {
            continue;
          }
          auto within = RunIntervalWithin(r, q0, q1, interaction);
          // Perpendicular projection onto the piece's line inside the (extended) part.
          const double g0 = Dot(Sub(runs[r].start, rp.start), rp.tangent);
          const double slope = Dot(runs[r].tangent, rp.tangent);
          std::vector<Interval> domain;
          if (std::abs(slope) <= kDirectionQuantum)
          {
            if (g0 >= lo_p - Tol() && g0 <= hi_p + Tol())
            {
              domain.emplace_back(0.0, runs[r].length);
            }
          }
          else
          {
            double s0 = (lo_p - g0) / slope, s1 = (hi_p - g0) / slope;
            if (s0 > s1)
            {
              std::swap(s0, s1);
            }
            s0 = std::clamp(s0, 0.0, runs[r].length);
            s1 = std::clamp(s1, 0.0, runs[r].length);
            if (s1 - s0 > Tol())
            {
              domain.emplace_back(s0, s1);
            }
          }
          within = IntersectIntervals(within, domain, Tol());
          found.insert(found.end(), within.begin(), within.end());
        }
        // A claimed piece ending where its CHAIN ends (a strip end, the tangent point of a
        // rounded corner's arc chain) faces the metal around that end within the full 2R
        // ball: nothing continues beyond it that could pair with the neighbour (the
        // perpendicular rule only guards a claim end inside a continuing chain).
        for (const bool at_start : {true, false})
        {
          const double s_end = at_start ? part.first : part.second;
          const std::size_t k = rp.index_in_chain;
          const bool run_end = at_start ? s_end <= Tol() : s_end >= rp.length - Tol();
          const bool chain_end = !B.closed && (at_start ? k == 0 : k + 1 == B.runs.size());
          if (run_end && chain_end)
          {
            const Point3D e = RunPoint(piece.run, s_end);
            const auto ball =
                RunIntervalWithinPiece(r, CurvePiece{e, e, std::nullopt}, interaction);
            found.insert(found.end(), ball.begin(), ball.end());
          }
        }
      }
      found = MergeIntervals(std::move(found), Tol());
      if (found.empty())
      {
        continue;
      }
      if (std::getenv("PALACE_IDENTIFICATION_DEBUG_EXTENSION_PIECES") && input.log)
      {
        std::ostringstream dbg;
        dbg << std::setprecision(10) << "    extension run " << r << " piece run "
            << piece.run << " s [" << piece.interval.first << ", " << piece.interval.second
            << "] ext " << piece.extend_lo << " / " << piece.extend_hi << " owner "
            << piece.owner << " found";
        for (const auto &interval : found)
        {
          dbg << " [" << interval.first << ", " << interval.second << "]";
        }
        input.log(dbg.str() + "\n");
      }
      if (rp.chain == runs[r].chain)
      {
        // The chain's own neighbourhood: only points more than pi R of arc length from the
        // piece can face it across a fold.
        const double b0 = B.run_offset[rp.index_in_chain] + piece.interval.first;
        const double b1 = B.run_offset[rp.index_in_chain] + piece.interval.second;
        std::vector<Interval> allowed;
        for (const auto &interval : found)
        {
          std::vector<Interval> near;
          for (const double shift : A.closed ? std::vector<double>{-A.length, 0.0, A.length}
                                             : std::vector<double>{0.0})
          {
            const double n0 = b0 + shift - neighbourhood - x_a0;
            const double n1 = b1 + shift + neighbourhood - x_a0;
            if (n1 > n0)
            {
              near.emplace_back(n0, n1);
            }
          }
          const auto far =
              SubtractIntervals({interval}, MergeIntervals(near, Tol()), Tol());
          allowed.insert(allowed.end(), far.begin(), far.end());
        }
        found = allowed;
      }
      else
      {
        const auto zone_key = std::make_pair(rp.chain, owner_cluster);
        auto zit = zones_by_chain.find(zone_key);
        if (zit == zones_by_chain.end())
        {
          zit = zones_by_chain.emplace(zone_key, ThroughZonesExcept(r, A, B, owner_cluster))
                    .first;
        }
        std::vector<Interval> excluded;
        for (const auto &zone : zit->second)
        {
          excluded.insert(excluded.end(), zone.begin(), zone.end());
        }
        found = SubtractIntervals(found, MergeIntervals(excluded, Tol()), Tol());
      }
      auto &list = within_by_owner[piece.owner];
      list.insert(list.end(), found.begin(), found.end());
    }
    if (within_by_owner.empty())
    {
      continue;
    }
    // Per remainder piece: the sub-intervals within 2R of every owner; overlapping
    // sub-intervals of different owners are one absorption with several owners (merged).
    for (const auto &free : remainder)
    {
      std::vector<std::pair<Interval, std::size_t>> hits;
      for (auto &[owner, list] : within_by_owner)
      {
        for (const auto &interval :
             IntersectIntervals({free}, MergeIntervals(list, Tol()), Tol()))
        {
          if (interval.second - interval.first > sliver)
          {
            hits.emplace_back(interval, owner);
          }
        }
      }
      if (hits.empty())
      {
        continue;
      }
      std::sort(hits.begin(), hits.end());
      // Merge overlapping / touching hits into absorptions carrying every owner.
      std::vector<Absorption> merged;
      for (const auto &[interval, owner] : hits)
      {
        if (!merged.empty() && interval.first <= merged.back().interval.second + Tol())
        {
          merged.back().interval.second =
              std::max(merged.back().interval.second, interval.second);
          merged.back().owners.push_back(owner);
        }
        else
        {
          merged.push_back({r, interval, {owner}});
        }
      }
      absorptions.insert(absorptions.end(), merged.begin(), merged.end());
    }
  }
  // Decision 224 (the S1p 41-edge loop end, 2026-10-02): the ONE exception to "pairs /
  // stacks are never absorbed" — pair / stack stretches that exist ONLY because a cluster's
  // claim boundary cut them are absorbed by that cluster. A stretch is a maximal contiguous
  // run of one cross-section's claims (one signature key) along a chain: the stack-end
  // recomposition piece follows the longer stack's claim on the same chain and is judged
  // on its own, while consecutive claims of one cross-section are one stretch whatever
  // their feature ids (the stack assembly can emit one feature per mesh segment). It must
  // lie ENTIRELY within the cluster ball radius ClusterBallOverR x R of the cluster's
  // claimed pieces (3D distance, the same sublevel sets as the cluster's own claims, arcs
  // included) — shorter than the interaction distance from the cluster's own claims, so the
  // cluster's coupon describes it while the stack's translational patches would correct the
  // same surface a second time — and be one of two classes, counted separately:
  // TWO-SIDED, bounded at both ends (within the decision quantum) by claimed intervals of
  // the SAME cluster: a piece shorter than 2R between two of its claims (on S1p the three
  // leads of the loop end were each split by a 1.738 um piece of the 3-edge stack between
  // the cluster's 5.131 um claims); STACK-END RECOMPOSITION, adjacent to the cluster at one
  // end and continuing, at the other, a pair / stack claim of a cross-section with strictly
  // more members containing its own — the members the cluster takes later continue past
  // the member it takes first (the S1p 3-edge stack end, 3 x 1.0 um at y = -136.3, between
  // the 4-edge stack and the loop end's claims on the leads). Never between two DIFFERENT
  // clusters (the placement's spatial-vs-spatial overlap check covers overlapping coupon
  // boxes), never a stretch whose far end is free, a bend or a smaller cross-section:
  // strips and gaps alongside a cluster, bent-strip halves and the gap between two pad-end
  // clusters are genuine translational features the 2D models describe (the broader "ball
  // only" form tried the same day ate them: 10 identity cells, DS-SCT-001 +43.6 um of
  // cluster). Applied per pass like the single-edge absorptions (the next pass recomposes
  // the stacks around the enlarged claims); the stretch is cut into its runs' intervals.
  // Diagnostics.ClusterExtension.TranslationalPiecesAbsorbed. A translational stretch the
  // rule leaves inside a spatial coupon's volume is recorded by the placement's ownership
  // record (decision 236: Continuation / Foreign, never an abort).
  if (!measure_joint_claims)
  {
    const double ball = kClusterBallOverRadius * R;
    std::map<int, std::vector<std::tuple<double, double, std::size_t>>> cluster_on_chain;
    for (std::size_t c = 0; c < cluster_claimed.size(); c++)
    {
      for (const auto &[r, interval] : cluster_claimed[c])
      {
        const Chain &C = chains[chain_index.at(runs[r].chain)];
        const double offset = C.run_offset[runs[r].index_in_chain];
        cluster_on_chain[C.id].emplace_back(offset + interval.first,
                                            offset + interval.second, c);
      }
    }
    // The pair / stack claims per (chain, cross-section = the feature's signature key): a
    // stretch is one cross-section's claim (the stack-end recomposition piece follows the
    // longer stack's claim on the same chain and is judged on its own; consecutive claims
    // of one cross-section are one stretch whatever their feature ids).
    std::map<std::pair<int, std::string>, std::vector<Interval>> translational_on_chain;
    // Every pair / stack claim per chain with its feature, and every feature's member
    // chains (the chains its claims lie on): the stack-end recomposition test.
    std::map<int, std::vector<std::tuple<double, double, int>>> translational_claims;
    std::map<int, std::set<int>> feature_members;
    for (std::size_t r = 0; r < runs.size(); r++)
    {
      if (runs[r].excluded)
      {
        continue;
      }
      const Chain &C = chains[chain_index.at(runs[r].chain)];
      const double offset = C.run_offset[runs[r].index_in_chain];
      for (const auto &claim : claims[r])
      {
        if (claim.priority == 2)
        {
          translational_on_chain[std::make_pair(
                                     C.id, features[static_cast<std::size_t>(claim.feature)]
                                               .signature_key)]
              .emplace_back(offset + claim.interval.first, offset + claim.interval.second);
          translational_claims[C.id].emplace_back(
              offset + claim.interval.first, offset + claim.interval.second, claim.feature);
          feature_members[claim.feature].insert(C.id);
        }
      }
    }
    // The stretch's intervals on its runs lie entirely within the ball radius of the
    // cluster's claimed pieces (the pieces indexed in piece_grid with owner c; same plane).
    auto WithinBall =
        [&](std::size_t c, const std::vector<std::pair<std::size_t, Interval>> &on_runs)
    {
      for (const auto &[r, interval] : on_runs)
      {
        Point3D lo, hi;
        PieceBoundingBox(RunPiece(r, interval.first, interval.second), lo, hi);
        std::vector<Interval> covered;
        for (const std::size_t i : piece_grid.Query(lo, hi, ball + 2.0 * Tol()))
        {
          const Piece &piece = pieces[i];
          if (piece.owner != c || runs[piece.run].excluded ||
              run_plane[piece.run] != run_plane[r])
          {
            continue;
          }
          const auto within = RunIntervalWithinPiece(
              r, RunPiece(piece.run, piece.interval.first, piece.interval.second), ball);
          covered.insert(covered.end(), within.begin(), within.end());
        }
        if (!SubtractIntervals({interval}, MergeIntervals(std::move(covered), Tol()), Tol())
                 .empty())
        {
          return false;
        }
      }
      return true;
    };
    for (auto &[chain_feature, list] : translational_on_chain)
    {
      const int chain_id = chain_feature.first;
      const std::string &cross_section = chain_feature.second;
      const auto bounding = cluster_on_chain.find(chain_id);
      if (bounding == cluster_on_chain.end())
      {
        continue;
      }
      const Chain &C = chains[chain_index.at(chain_id)];
      const auto &chain_claims = translational_claims.at(chain_id);
      // The stretch's members: the union of the member chains of the features whose claims
      // form it; the stretch at x_far continues a larger stack when a claim of another
      // cross-section with strictly more members, containing them all, ends (or starts)
      // there.
      auto ContinuesLargerStack = [&](const Interval &stretch, double x_far)
      {
        std::set<int> members;
        for (const auto &[x0, x1, feature] : chain_claims)
        {
          if (x0 >= stretch.first - Tol() && x1 <= stretch.second + Tol() &&
              features[static_cast<std::size_t>(feature)].signature_key == cross_section)
          {
            const auto &chains_of = feature_members[feature];
            members.insert(chains_of.begin(), chains_of.end());
          }
        }
        for (const auto &[x0, x1, feature] : chain_claims)
        {
          const bool adjacent =
              std::abs(x1 - x_far) <= Tol() || std::abs(x0 - x_far) <= Tol();
          const bool inside = x0 >= stretch.first - Tol() && x1 <= stretch.second + Tol();
          if (!adjacent || inside ||
              features[static_cast<std::size_t>(feature)].signature_key == cross_section)
          {
            continue;  // not at the far end, or one of the stretch's own claims
          }
          const auto &larger = feature_members[feature];
          if (larger.size() > members.size() &&
              std::includes(larger.begin(), larger.end(), members.begin(), members.end()))
          {
            return true;
          }
        }
        return false;
      };
      // The cluster whose claim ends at x (a claim ending at the chain's end bounds a
      // stretch starting at 0 on a closed chain, and vice versa); -1 when none.
      auto ClusterEndingAt = [&](double x)
      {
        for (const auto &[x0, x1, c] : bounding->second)
        {
          if (std::abs(x1 - x) <= Tol() ||
              (C.closed && x <= Tol() && x1 >= C.length - Tol()))
          {
            return static_cast<int>(c);
          }
        }
        return -1;
      };
      auto ClusterStartingAt = [&](double x)
      {
        for (const auto &[x0, x1, c] : bounding->second)
        {
          if (std::abs(x0 - x) <= Tol() ||
              (C.closed && x >= C.length - Tol() && x0 <= Tol()))
          {
            return static_cast<int>(c);
          }
        }
        return -1;
      };
      for (const auto &stretch : MergeIntervals(std::move(list), Tol()))
      {
        const int before = ClusterEndingAt(stretch.first);
        const int after = ClusterStartingAt(stretch.second);
        if (before < 0 && after < 0)
        {
          continue;
        }
        std::vector<std::pair<std::size_t, Interval>> on_runs;
        for (std::size_t k = 0; k < C.runs.size(); k++)
        {
          const std::size_t r = C.runs[k];
          const double offset = C.run_offset[k];
          const double s0 = std::max(stretch.first, offset) - offset;
          const double s1 = std::min(stretch.second, offset + runs[r].length) - offset;
          if (!runs[r].excluded && s1 - s0 > Tol())
          {
            on_runs.emplace_back(r, Interval{s0, s1});
          }
        }
        // The owner: the same cluster at both ends (two-sided), or the cluster at one end
        // when the other end continues the larger stack the stretch's members belong to
        // (stack-end recomposition). Two different clusters at the two ends: no absorption.
        int owner = -1;
        bool two_sided = false;
        if (before >= 0 && before == after)
        {
          owner = before;
          two_sided = true;
        }
        else if (before >= 0 && after < 0 && ContinuesLargerStack(stretch, stretch.second))
        {
          owner = before;
        }
        else if (after >= 0 && before < 0 && ContinuesLargerStack(stretch, stretch.first))
        {
          owner = after;
        }
        if (owner < 0 || !WithinBall(static_cast<std::size_t>(owner), on_runs))
        {
          continue;
        }
        for (const auto &[r, interval] : on_runs)
        {
          absorptions.push_back({r, interval, {static_cast<std::size_t>(owner)}, true});
        }
        const double length = stretch.second - stretch.first;
        extension_translational_pieces++;
        extension_translational_length += length;
        extension_translational_max_length =
            std::max(extension_translational_max_length, length);
        if (two_sided)
        {
          extension_translational_two_sided++;
          extension_translational_two_sided_length += length;
        }
        if (std::getenv("PALACE_IDENTIFICATION_DEBUG_EXTENSION") && input.log)
        {
          std::ostringstream line;
          line << std::setprecision(10) << "    translational stretch absorbed ("
               << (two_sided ? "two-sided" : "stack-end recomposition") << "): chain "
               << chain_id << " x [" << stretch.first << ", " << stretch.second
               << "] (length " << length << ") by cluster " << owner << "\n";
          input.log(line.str());
        }
      }
    }
  }
  double candidate_length = 0.0;
  for (const auto &absorption : absorptions)
  {
    candidate_length += absorption.interval.second - absorption.interval.first;
  }
  if (measure_joint_claims)
  {
    if (std::getenv("PALACE_IDENTIFICATION_DEBUG_EXTENSION") && input.log)
    {
      // The measured pair / stack claim intervals (the analytic reading of
      // Diagnostics.StackEndThirdBodyLength), for the comparison with the facing gate's
      // sampled StackEndThirdBody / ThroughVertex classes (review fix-4).
      std::ostringstream dbg;
      dbg << std::setprecision(10);
      for (const auto &absorption : absorptions)
      {
        const Run &run = runs[absorption.run];
        dbg << "    third body run " << absorption.run << " chain " << run.chain << " s ["
            << absorption.interval.first << ", " << absorption.interval.second << "] of "
            << run.length << " from ("
            << RunPoint(absorption.run, absorption.interval.first)[0] << ", "
            << RunPoint(absorption.run, absorption.interval.first)[1] << ") to ("
            << RunPoint(absorption.run, absorption.interval.second)[0] << ", "
            << RunPoint(absorption.run, absorption.interval.second)[1] << ") owners";
        for (const std::size_t owner : absorption.owners)
        {
          dbg << " " << owner;
        }
        dbg << "\n";
      }
      input.log(dbg.str());
    }
    return candidate_length;  // measurement only: nothing absorbed
  }
  // Owners of one absorption are one cluster (union-find over clusters and free sites).
  const std::size_t n_owner = cluster_claimed.size() + free_sites.size();
  UnionFind uf(n_owner);
  for (const auto &absorption : absorptions)
  {
    for (std::size_t k = 1; k < absorption.owners.size(); k++)
    {
      uf.Union(absorption.owners[0], absorption.owners[k]);
    }
  }
  // Every free site that took part becomes a cluster of its own (its window), then the
  // roots are merged: new cluster indices in the order of the smallest member owner.
  std::vector<std::size_t> owner_cluster(n_owner, std::numeric_limits<std::size_t>::max());
  std::vector<std::vector<std::pair<std::size_t, Interval>>> new_claimed;
  std::vector<std::vector<EventCore>> new_cores;
  std::vector<std::vector<std::size_t>> new_sites;
  std::set<std::size_t> involved_free_sites;
  for (const auto &absorption : absorptions)
  {
    for (const std::size_t owner : absorption.owners)
    {
      if (owner >= cluster_claimed.size())
      {
        involved_free_sites.insert(owner);
      }
    }
  }
  std::map<std::size_t, std::size_t> root_cluster;
  for (std::size_t owner = 0; owner < n_owner; owner++)
  {
    const bool is_site = owner >= cluster_claimed.size();
    if (is_site && involved_free_sites.count(owner) == 0)
    {
      continue;  // a vertex feature nothing joined stays a vertex feature
    }
    const std::size_t root = uf.Find(owner);
    auto [it, inserted] = root_cluster.emplace(root, new_claimed.size());
    if (inserted)
    {
      new_claimed.emplace_back();
      new_cores.emplace_back();
      new_sites.emplace_back();
    }
    const std::size_t c = it->second;
    owner_cluster[owner] = c;
    if (is_site)
    {
      const std::size_t s = free_sites[owner - cluster_claimed.size()];
      new_sites[c].push_back(s);
      for (const auto &[r, interval] : sites[s].window)
      {
        new_claimed[c].emplace_back(r, interval);
      }
    }
    else
    {
      new_claimed[c].insert(new_claimed[c].end(), cluster_claimed[owner].begin(),
                            cluster_claimed[owner].end());
      new_cores[c].insert(new_cores[c].end(), cluster_cores[owner].begin(),
                          cluster_cores[owner].end());
      new_sites[c].insert(new_sites[c].end(), cluster_sites[owner].begin(),
                          cluster_sites[owner].end());
    }
  }
  double absorbed = 0.0, absorbed_translational = 0.0;
  const bool debug_extension =
      std::getenv("PALACE_IDENTIFICATION_DEBUG_EXTENSION") && input.log;
  for (const auto &absorption : absorptions)
  {
    new_claimed[owner_cluster[absorption.owners[0]]].emplace_back(absorption.run,
                                                                  absorption.interval);
    absorbed += absorption.interval.second - absorption.interval.first;
    if (absorption.translational)
    {
      absorbed_translational += absorption.interval.second - absorption.interval.first;
    }
    else
    {
      extension_portions++;
    }
    if (debug_extension)
    {
      const Run &run = runs[absorption.run];
      std::ostringstream line;
      line << "    absorbed run " << absorption.run << " chain " << run.chain << " s ["
           << absorption.interval.first << ", " << absorption.interval.second << "] of "
           << run.length << " from ("
           << RunPoint(absorption.run, absorption.interval.first)[0] << ", "
           << RunPoint(absorption.run, absorption.interval.first)[1] << ") to ("
           << RunPoint(absorption.run, absorption.interval.second)[0] << ", "
           << RunPoint(absorption.run, absorption.interval.second)[1] << ") owners";
      for (const std::size_t owner : absorption.owners)
      {
        line << " " << owner;
      }
      input.log(line.str() + "\n");
    }
  }
  extension_length += absorbed - absorbed_translational;  // single-edge portions
  extension_sites += involved_free_sites.size();
  // Merge the intervals per run inside every cluster.
  for (auto &claimed : new_claimed)
  {
    std::map<std::size_t, std::vector<Interval>> by_run;
    for (const auto &[r, interval] : claimed)
    {
      by_run[r].push_back(interval);
    }
    claimed.clear();
    for (auto &[r, list] : by_run)
    {
      for (const auto &interval : MergeIntervals(std::move(list), Tol()))
      {
        claimed.emplace_back(r, interval);
      }
    }
    std::sort(claimed.begin(), claimed.end());
  }
  cluster_claimed = std::move(new_claimed);
  cluster_cores = std::move(new_cores);
  cluster_sites = std::move(new_sites);
  for (std::size_t c = 0; c < cluster_sites.size(); c++)
  {
    std::sort(cluster_sites[c].begin(), cluster_sites[c].end());
    for (const std::size_t s : cluster_sites[c])
    {
      sites[s].cluster = static_cast<int>(c);
    }
  }
  return absorbed;
}

// The cluster features from the final claimed portions and member sites (after the
// extension): canonical signature, priority-0 claims, vertex membership.
// The signature portions of a cluster's claimed run intervals: a chord portion per
// straight run interval; the intervals on the runs of one fitted arc merged, in angular
// order, wherever they abut (within the decision tolerance along the arc) with the same
// conductor, interfaces, law and radial gap sign, into ONE arc portion (option A: the
// serialisation of a cluster holding a rounded corner or a bend does not depend on the
// chord count). A closed circle claimed whole is one portion of sweep 2 pi.
std::vector<SignaturePortion> Identifier::ClusterSignaturePortions(
    const std::vector<std::pair<std::size_t, Interval>> &claimed) const
{
  std::vector<SignaturePortion> portions;
  struct ArcPart
  {
    double lo, hi;  // angles in the arc's frame
    std::size_t run;
    int gap_radial;
  };
  std::map<int, std::vector<ArcPart>> parts_by_arc;
  std::map<int, SignatureArc> frame_by_arc;
  const double two_pi = 2.0 * std::acos(-1.0);
  for (const auto &[r, interval] : claimed)
  {
    if (!run_arcs[r])
    {
      portions.push_back({runs[r].At(interval.first), runs[r].At(interval.second),
                          runs[r].gap_direction, runs[r].conductor,
                          InterfaceNames(runs[r].targets), runs[r].boundary_law});
      continue;
    }
    const int a = run_arcs[r]->arc;
    const Arc &fitted = arcs[static_cast<std::size_t>(a)];
    auto fit = frame_by_arc.find(a);
    if (fit == frame_by_arc.end())
    {
      // The arc's angular frame: u along the radial direction of the arc's MIDPOINT (the
      // arc then spans [-turn / 2, turn / 2]: no wrap-around for a turn below 2 pi), v = n
      // x u with n the plane normal of the arc's runs.
      const Point3D n = Normalize(runs[r].process_normal);
      Point3D d_a = Sub(input.vertices[fitted.joints.front()].coordinate, fitted.center);
      d_a = Normalize(Sub(d_a, Scale(Dot(d_a, n), n)));
      const Point3D v_a = Normalize(Cross(n, d_a));
      double sign = 1.0;
      if (fitted.joints.size() > 1)
      {
        const Point3D d_1 = Sub(input.vertices[fitted.joints[1]].coordinate, fitted.center);
        sign = Dot(d_1, v_a) >= 0.0 ? 1.0 : -1.0;
      }
      const double mid = 0.5 * sign * fitted.turn;
      SignatureArc frame;
      frame.center = fitted.center;
      frame.radius = fitted.radius;
      frame.u = Add(Scale(std::cos(mid), d_a), Scale(std::sin(mid), v_a));
      frame.v = Normalize(Cross(n, frame.u));
      fit = frame_by_arc.emplace(a, frame).first;
    }
    const SignatureArc &frame = fit->second;
    auto AngleOf = [&](const Point3D &p)
    {
      const Point3D d = Sub(p, frame.center);
      return std::atan2(Dot(d, frame.v), Dot(d, frame.u));
    };
    const double s_mid = 0.5 * (interval.first + interval.second);
    const Point3D radial = SignatureArcRadial(frame, AngleOf(RunPoint(r, s_mid)));
    const int gap_radial = Dot(runs[r].gap_direction, radial) >= 0.0 ? 1 : -1;
    double t0 = AngleOf(RunPoint(r, interval.first));
    double t1 = AngleOf(RunPoint(r, interval.second));
    if (t0 > t1)
    {
      std::swap(t0, t1);
    }
    if (t1 - t0 > std::acos(-1.0))
    {
      // The chord crosses the seam at +-pi (a closed circle's antipode of the first joint):
      // split it there.
      parts_by_arc[a].push_back({t1, std::acos(-1.0), r, gap_radial});
      t1 = t0;
      t0 = -std::acos(-1.0);
    }
    parts_by_arc[a].push_back({t0, t1, r, gap_radial});
  }
  for (auto &[a, parts] : parts_by_arc)
  {
    const SignatureArc &frame = frame_by_arc.at(a);
    const double angle_tolerance = Tol() / frame.radius;
    std::sort(parts.begin(), parts.end(),
              [](const ArcPart &x, const ArcPart &y) { return x.lo < y.lo; });
    auto SameKind = [&](const ArcPart &x, const ArcPart &y)
    {
      const Run &rx = runs[x.run], &ry = runs[y.run];
      return x.gap_radial == y.gap_radial && rx.conductor == ry.conductor &&
             rx.targets == ry.targets && rx.boundary_law == ry.boundary_law;
    };
    std::vector<ArcPart> merged;
    for (const auto &part : parts)
    {
      if (!merged.empty() && part.lo <= merged.back().hi + angle_tolerance &&
          SameKind(merged.back(), part))
      {
        merged.back().hi = std::max(merged.back().hi, part.hi);
      }
      else
      {
        merged.push_back(part);
      }
    }
    // Wrap-around at the seam (+-pi): the last part reaching pi joins the first starting at
    // -pi.
    if (merged.size() >= 2 && merged.back().hi >= std::acos(-1.0) - angle_tolerance &&
        merged.front().lo <= -std::acos(-1.0) + angle_tolerance &&
        SameKind(merged.back(), merged.front()))
    {
      merged.back().hi += merged.front().hi + std::acos(-1.0);
      merged.erase(merged.begin());
    }
    for (const auto &part : merged)
    {
      SignatureArc arc = frame;
      arc.theta0 = part.lo;
      arc.theta1 = part.hi;
      arc.gap_radial = part.gap_radial;
      if (arc.theta1 - arc.theta0 >= two_pi - angle_tolerance)
      {
        arc.theta0 = 0.0;
        arc.theta1 = two_pi;
      }
      const Run &run = runs[part.run];
      SignaturePortion portion;
      portion.p0 = SignatureArcPoint(arc, arc.theta0);
      portion.p1 = SignatureArcPoint(arc, arc.theta1);
      portion.gap_direction = run.gap_direction;
      portion.conductor = run.conductor;
      portion.interfaces = InterfaceNames(run.targets);
      portion.boundary_law = run.boundary_law;
      portion.arc = arc;
      portions.push_back(std::move(portion));
    }
  }
  return portions;
}

// The support of cluster c in the frame (origin, x, y); every length below is in units of
// R in the frame unless stated. Rule B2 gives the box from the claims serialised in the
// frame; the device plan of the cluster's plane (every non-excluded run, straight or on its
// fitted arc) is clipped to the box, the cluster's own claims removed, pieces shorter than
// the snap quantum dropped (T1); T2 is evaluated per face and every failing face grows by
// one step (T3), all failing faces per step, until every face passes or the step cap / the
// span cap refuses the growth (unboxable).
Identifier::ClusterSupportResult
Identifier::ClusterSupport(std::size_t c, const std::vector<SignaturePortion> &portions,
                           const std::vector<SignatureVertex> &vertices,
                           const Point3D &origin, const Point3D &x, const Point3D &y,
                           double span_cap_over_R) const
{
  using Point2 = std::array<double, 2>;
  ClusterSupportResult out;
  const double snap = kSupportFaceSnapOverRadius;
  const double clearance = kSupportFaceClearanceOverRadius;
  const double step = kSupportFaceGrowthStepOverRadius;
  const double span_cap = span_cap_over_R;
  const double band = kKnifeEdgeBandRelative;
  const double pi = std::acos(-1.0);
  // Quantised face-rule comparisons (decision 287 (b)): the box faces sit on the
  // kSignatureLengthQuantumOverRadius grid and a step moves them by a grid multiple, so a
  // device edge lying ON a claims-box face reads exactly the clearance from the moved face
  // up to the quantisation residual (|r| <= half a quantum) and mesh / float noise. Every
  // length threshold is therefore read on the grid: a value within half a quantum of the
  // threshold is AT the threshold and decided by the rule's own inclusive side (>= passes
  // the clearance, <= the snap is on the face), never by sub-quantum noise.
  const double half_quantum = 0.5 * kSignatureLengthQuantumOverRadius;
  auto AtMostSnap = [&](double distance) { return distance < snap + half_quantum; };
  auto BelowClearance = [&](double distance)
  { return distance < clearance - half_quantum; };
  auto IsSliverLength = [&](double length) { return length / R < snap - half_quantum; };
  auto Local = [&](const Point3D &p) -> Point2
  {
    const Point3D r = Sub(p, origin);
    return {Dot(r, x) / R, Dot(r, y) / R};
  };
  auto LocalDirection = [&](const Point3D &d) -> Point2
  {
    const Point2 v = {Dot(d, x), Dot(d, y)};
    const double norm = std::hypot(v[0], v[1]);
    return norm > 0.0 ? Point2{v[0] / norm, v[1] / norm} : Point2{0.0, 0.0};
  };
  auto Q = [](double v) { return RoundTo(v, kSignatureLengthQuantumOverRadius); };

  // Rule B2: the claims-derived box from the claims serialised in this frame.
  const nlohmann::json claims = SerializeInFrame(portions, vertices, origin, x, y, R);
  std::size_t box_band_hits = 0;
  const std::array<double, 4> claims_box = SupportBoxFromSignature(claims, &box_band_hits);
  std::array<double, 4> box = claims_box;
  const int plane = run_plane[cluster_claimed[c].front().first];
  std::map<std::size_t, std::vector<Interval>> claimed_by_run;
  for (const auto &[r, interval] : cluster_claimed[c])
  {
    claimed_by_run[r].push_back(interval);
  }
  for (auto &[r, intervals] : claimed_by_run)
  {
    intervals = MergeIntervals(std::move(intervals), Tol());
  }

  // Face f: 0 = x0, 1 = y0, 2 = x1, 3 = y1; its line coordinate, axis (0 = x, 1 = y) and
  // the extent of the face along the other axis.
  auto FaceLine = [&](int f) { return box[static_cast<std::size_t>(f)]; };
  auto FaceAxis = [](int f) { return f % 2; };
  auto Inside = [&](const Point2 &p, double tolerance)
  {
    return p[0] >= box[0] - tolerance && p[0] <= box[2] + tolerance &&
           p[1] >= box[1] - tolerance && p[1] <= box[3] + tolerance;
  };

  // A plan piece: a run-parameter interval of a run, or (run == kExcludedRun) a parameter
  // interval of a perimeter SEGMENT excluded before the identification (a port cut, an
  // undetermined process side, a non-manifold or non-planar edge: rule B1, nothing inside
  // the box is omitted — such a segment is device geometry the coupon must carry, but it
  // is no feature edge: foreign context, never a chain piece, excluded from the within-R
  // accounting like the device's own perimeter excludes it).
  constexpr std::size_t kExcludedRun = std::numeric_limits<std::size_t>::max();
  struct Piece
  {
    std::size_t run;
    Interval s;
    // The piece lies on another cluster's claimed interval (decision 285 (1)): owned by
    // that coupon, never part of this cluster's chain.
    std::optional<std::size_t> claimed_by;
    std::size_t segment = kExcludedRun;  // the excluded segment of a kExcludedRun piece
    // The end at s.first / s.second exists only because another cluster's claim cut split
    // the run there (nothing geometric: the device edge continues): exempt from the T2
    // end tests (R1 final review MINOR-1), still a crossing when it lies on a face.
    bool cut_a = false, cut_b = false;
    bool Excluded() const { return run == kExcludedRun; }
    std::pair<std::size_t, std::size_t> Key() const { return {run, segment}; }
  };
  auto SegmentLength = [&](std::size_t i)
  { return Distance(input.segments[i].p0, input.segments[i].p1); };
  auto SegmentAt = [&](std::size_t i, double s)
  {
    const auto &segment = input.segments[i];
    const double length = SegmentLength(i);
    return length > 0.0 ? Add(segment.p0, Scale(s / length, Sub(segment.p1, segment.p0)))
                        : segment.p0;
  };
  // The geometry of a piece between the parameters s0 and s1 (a run piece on its fitted
  // arc, an excluded segment straight).
  auto PieceCurve = [&](const Piece &piece, double s0, double s1)
  {
    if (!piece.Excluded())
    {
      return RunPiece(piece.run, s0, s1);
    }
    CurvePiece curve;
    curve.a = SegmentAt(piece.segment, s0);
    curve.b = SegmentAt(piece.segment, s1);
    return curve;
  };
  // The in-plane gap direction of an excluded segment (unset on the segment: a port cut has
  // no target interface): the normal pointing away from the metal, read at an end vertex
  // from an adjacent feature run — a point just inside the run's metal (one step along the
  // run away from the vertex and one step away from its gap) lies on the metal side of the
  // segment. Zero when no adjacent run exists (the piece is then recorded, not oriented).
  auto ExcludedGapDirection = [&](std::size_t i) -> Point3D
  {
    const auto &segment = input.segments[i];
    if (Norm(segment.gap_direction) > 0.0)
    {
      return Normalize(segment.gap_direction);
    }
    const Point3D normal = Normalize(Cross(x, y));
    const Point3D t = Normalize(Sub(segment.p1, segment.p0));
    const Point3D n0 = Normalize(Cross(t, normal));
    for (const std::size_t v : segment.vertices)
    {
      if (v >= input.vertices.size())
      {
        continue;
      }
      const Point3D &vertex = input.vertices[v].coordinate;
      for (const std::size_t s : input.vertices[v].segments)
      {
        if (s == i || s >= run_of_segment.size() || run_of_segment[s] < 0)
        {
          continue;
        }
        const Run &run = runs[static_cast<std::size_t>(run_of_segment[s])];
        if (run.excluded || Norm(run.gap_direction) <= 0.0)
        {
          continue;
        }
        const auto &adjacent = input.segments[s];
        const Point3D far = Distance(adjacent.p0, vertex) > Distance(adjacent.p1, vertex)
                                ? adjacent.p0
                                : adjacent.p1;
        const Point3D along = Normalize(Sub(far, vertex));
        const Point3D inside = Sub(along, Normalize(run.gap_direction));
        const double side = Dot(inside, n0);
        if (std::abs(side) > 1.0e-9)
        {
          return side < 0.0 ? n0 : Scale(-1.0, n0);
        }
      }
    }
    return Point3D{};
  };
  auto PieceGapDirection = [&](const Piece &piece) -> Point3D
  {
    return piece.Excluded() ? ExcludedGapDirection(piece.segment)
                            : runs[piece.run].gap_direction;
  };
  struct ArcLocal
  {
    Point2 center;
    double radius;
    double phi_x, amplitude_x;  // x(theta) = cx + r A_x cos(theta - phi_x)
    double phi_y, amplitude_y;
    Point2 u, v;
  };
  auto ArcInFrame = [&](const ArcPiece &arc)
  {
    ArcLocal a;
    a.center = Local(arc.center);
    a.radius = arc.radius / R;
    a.u = {Dot(arc.u, x), Dot(arc.u, y)};
    a.v = {Dot(arc.v, x), Dot(arc.v, y)};
    a.amplitude_x = std::hypot(a.u[0], a.v[0]);
    a.phi_x = std::atan2(a.v[0], a.u[0]);
    a.amplitude_y = std::hypot(a.u[1], a.v[1]);
    a.phi_y = std::atan2(a.v[1], a.u[1]);
    return a;
  };
  auto ArcPointLocal = [&](const ArcLocal &a, double theta) -> Point2
  {
    return {a.center[0] + a.radius * (std::cos(theta) * a.u[0] + std::sin(theta) * a.v[0]),
            a.center[1] + a.radius * (std::cos(theta) * a.u[1] + std::sin(theta) * a.v[1])};
  };
  // Length of a run piece (mesh units): the chord, or the arc length on the fitted arc.
  auto PieceLength = [&](const Piece &piece)
  {
    if (piece.Excluded())
    {
      return piece.s.second - piece.s.first;
    }
    const auto &geometry = run_arcs[piece.run];
    if (geometry)
    {
      const double length = runs[piece.run].length;
      const double sweep = std::abs(geometry->piece.theta1 - geometry->piece.theta0);
      return length > 0.0 ? geometry->piece.radius * sweep *
                                (piece.s.second - piece.s.first) / length
                          : 0.0;
    }
    return piece.s.second - piece.s.first;
  };

  // The device plan clipped to a box (the current box, or the box dilated by the clearance
  // for the two-sided face rule): run-parameter intervals inside it minus the cluster's
  // claims, split at the other clusters' claim cuts when `split_at_claim_cuts` (the box
  // pieces: the Chain / ClaimedByOtherFeature classification; the shell pieces of the
  // two-sided rule are never classed and stay whole so that a cut inside the shell is no
  // piece end), pieces shorter than the snap quantum dropped (T1).
  std::size_t slivers_dropped = 0;
  auto ClipPlanTo = [&](const std::array<double, 4> &clip, bool split_at_claim_cuts)
  {
    std::vector<Piece> pieces;
    slivers_dropped = 0;
    auto InsideClip = [&](const Point2 &p, double tolerance)
    {
      return p[0] >= clip[0] - tolerance && p[0] <= clip[2] + tolerance &&
             p[1] >= clip[1] - tolerance && p[1] <= clip[3] + tolerance;
    };
    Point3D lo{}, hi{};
    bool first = true;
    for (const double bx : {clip[0], clip[2]})
    {
      for (const double by : {clip[1], clip[3]})
      {
        const Point3D corner = Add(origin, Add(Scale(bx * R, x), Scale(by * R, y)));
        if (first)
        {
          lo = hi = corner;
          first = false;
          continue;
        }
        for (int d = 0; d < 3; d++)
        {
          lo[d] = std::min(lo[d], corner[d]);
          hi[d] = std::max(hi[d], corner[d]);
        }
      }
    }
    const double margin = max_run_sagitta + 2.0 * Tol() + snap * R;
    for (const std::size_t r : RunsNear(lo, hi, margin))
    {
      if (runs[r].excluded || run_plane[r] != plane || runs[r].length <= 0.0)
      {
        continue;
      }
      std::vector<Interval> inside;
      const double length = runs[r].length;
      if (!run_arcs[r])
      {
        // Liang-Barsky on the chord in the frame.
        const Point2 p = Local(runs[r].start), q = Local(runs[r].end);
        double t0 = 0.0, t1 = 1.0;
        bool empty = false;
        for (int axis = 0; axis < 2 && !empty; axis++)
        {
          const double lo_a = clip[static_cast<std::size_t>(axis)],
                       hi_a = clip[static_cast<std::size_t>(axis + 2)];
          const double d = q[axis] - p[axis];
          if (std::abs(d) <= 1.0e-15)
          {
            if (p[axis] < lo_a - snap || p[axis] > hi_a + snap)
            {
              empty = true;
            }
            continue;
          }
          double ta = (lo_a - p[axis]) / d, tb = (hi_a - p[axis]) / d;
          if (ta > tb)
          {
            std::swap(ta, tb);
          }
          t0 = std::max(t0, ta);
          t1 = std::min(t1, tb);
          empty = t0 >= t1;
        }
        if (!empty)
        {
          inside.emplace_back(t0 * length, t1 * length);
        }
      }
      else
      {
        const ArcPiece &arc = run_arcs[r]->piece;
        const ArcLocal a = ArcInFrame(arc);
        const double ta = std::min(arc.theta0, arc.theta1),
                     tb = std::max(arc.theta0, arc.theta1);
        std::vector<double> breaks = {ta, tb};
        for (int f = 0; f < 4; f++)
        {
          const int axis = FaceAxis(f);
          const double amplitude = axis == 0 ? a.amplitude_x : a.amplitude_y;
          const double phi = axis == 0 ? a.phi_x : a.phi_y;
          if (a.radius * amplitude <= 0.0)
          {
            continue;
          }
          const double cosine =
              (clip[static_cast<std::size_t>(f)] - a.center[axis]) / (a.radius * amplitude);
          if (std::abs(cosine) > 1.0)
          {
            continue;
          }
          const double delta = std::acos(std::clamp(cosine, -1.0, 1.0));
          for (const double root : {phi + delta, phi - delta})
          {
            // Every representative of the root inside [ta, tb].
            const double base = root - 2.0 * pi * std::floor((root - ta) / (2.0 * pi));
            for (double theta = base; theta <= tb + 1.0e-12; theta += 2.0 * pi)
            {
              if (theta >= ta - 1.0e-12)
              {
                breaks.push_back(std::clamp(theta, ta, tb));
              }
            }
          }
        }
        std::sort(breaks.begin(), breaks.end());
        auto ToS = [&](double theta)
        {
          const double t = (theta - arc.theta0) / (arc.theta1 - arc.theta0);
          return std::clamp(t, 0.0, 1.0) * length;
        };
        for (std::size_t k = 0; k + 1 < breaks.size(); k++)
        {
          if (breaks[k + 1] - breaks[k] <= 1.0e-12)
          {
            continue;
          }
          const double mid = 0.5 * (breaks[k] + breaks[k + 1]);
          if (InsideClip(ArcPointLocal(a, mid), 0.0))
          {
            double s0 = ToS(breaks[k]), s1 = ToS(breaks[k + 1]);
            if (s0 > s1)
            {
              std::swap(s0, s1);
            }
            inside.emplace_back(s0, s1);
          }
        }
        inside = MergeIntervals(std::move(inside), Tol());
      }
      if (inside.empty())
      {
        continue;
      }
      if (auto it = claimed_by_run.find(r); it != claimed_by_run.end())
      {
        inside = SubtractIntervals(inside, it->second, Tol());
      }
      // Split at the other clusters' claim cuts on this run and class every piece lying on
      // such a claim by its owner (decision 285 (1)).
      std::vector<std::pair<Interval, std::size_t>> other_claims;
      if (auto it = claims_by_run.find(r); it != claims_by_run.end())
      {
        for (const auto &[interval, owner] : it->second)
        {
          if (owner != c)
          {
            other_claims.emplace_back(interval, owner);
          }
        }
      }
      // (interval, end a is a cut, end b is a cut)
      std::vector<std::tuple<Interval, bool, bool>> split;
      if (!other_claims.empty() && split_at_claim_cuts)
      {
        for (const auto &interval : inside)
        {
          std::vector<double> cuts = {interval.first, interval.second};
          for (const auto &[claim, owner] : other_claims)
          {
            (void)owner;
            for (const double cut : {claim.first, claim.second})
            {
              if (cut > interval.first + Tol() && cut < interval.second - Tol())
              {
                cuts.push_back(cut);
              }
            }
          }
          std::sort(cuts.begin(), cuts.end());
          for (std::size_t k = 0; k + 1 < cuts.size(); k++)
          {
            if (cuts[k + 1] > cuts[k] + Tol())
            {
              split.emplace_back(Interval{cuts[k], cuts[k + 1]}, k > 0,
                                 k + 2 < cuts.size());
            }
          }
        }
      }
      else
      {
        for (const auto &interval : inside)
        {
          split.emplace_back(interval, false, false);
        }
      }
      for (const auto &[interval, cut_a, cut_b] : split)
      {
        Piece piece{r, interval, std::nullopt};
        piece.cut_a = cut_a;
        piece.cut_b = cut_b;
        const double mid = 0.5 * (interval.first + interval.second);
        if (split_at_claim_cuts)
        {
          for (const auto &[claim, owner] : other_claims)
          {
            if (mid >= claim.first - Tol() && mid <= claim.second + Tol())
            {
              piece.claimed_by = owner;
              break;
            }
          }
        }
        if (IsSliverLength(PieceLength(piece)))
        {
          slivers_dropped++;
          continue;
        }
        pieces.push_back(piece);
      }
    }
    // The excluded perimeter segments on the plane inside the clip (rule B1): straight,
    // clipped by Liang-Barsky, never on the cluster's claims (claims are feature runs).
    const Point3D plane_normal = Normalize(Cross(x, y));
    for (std::size_t i = 0; i < input.segments.size(); i++)
    {
      const auto &segment = input.segments[i];
      const bool excluded =
          !segment.truncation &&
          (segment.exclusion || (i < segment_exclusion.size() && segment_exclusion[i]));
      if (!excluded)
      {
        continue;
      }
      if (std::abs(Dot(Sub(segment.p0, origin), plane_normal)) > 2.0 * Tol() ||
          std::abs(Dot(Sub(segment.p1, origin), plane_normal)) > 2.0 * Tol())
      {
        continue;
      }
      const double length = SegmentLength(i);
      if (length <= 0.0)
      {
        continue;
      }
      const Point2 p = Local(segment.p0), q = Local(segment.p1);
      double t0 = 0.0, t1 = 1.0;
      bool empty = false;
      for (int axis = 0; axis < 2 && !empty; axis++)
      {
        const double lo_a = clip[static_cast<std::size_t>(axis)],
                     hi_a = clip[static_cast<std::size_t>(axis + 2)];
        const double d = q[axis] - p[axis];
        if (std::abs(d) <= 1.0e-15)
        {
          empty = p[axis] < lo_a - snap || p[axis] > hi_a + snap;
          continue;
        }
        double ta = (lo_a - p[axis]) / d, tb = (hi_a - p[axis]) / d;
        if (ta > tb)
        {
          std::swap(ta, tb);
        }
        t0 = std::max(t0, ta);
        t1 = std::min(t1, tb);
        empty = t0 >= t1;
      }
      if (empty)
      {
        continue;
      }
      Piece piece{kExcludedRun, Interval{t0 * length, t1 * length}, std::nullopt, i};
      if (IsSliverLength(PieceLength(piece)))
      {
        slivers_dropped++;
        continue;
      }
      pieces.push_back(piece);
    }
    std::sort(pieces.begin(), pieces.end(), [](const Piece &a, const Piece &b)
              { return a.Key() != b.Key() ? a.Key() < b.Key() : a.s < b.s; });
    return pieces;
  };
  auto ClipPlan = [&]() { return ClipPlanTo(box, true); };

  // T2 on the current box and pieces: the failing faces and the record of the pass.
  struct FaceCheck
  {
    std::array<bool, 4> failing{};
    std::array<std::size_t, 4> crossings{};
    std::optional<double> min_clearance, min_crossing_sine, min_cross_section;
    std::size_t narrow_cross_sections = 0;
    std::size_t band_hits = 0;
    // Two-sided T2 (decision 285 (2)): device vertices / edges inside the clearance shell
    // OUTSIDE a face that do not cross it.
    std::size_t exterior_vertices = 0, exterior_edges = 0;
    std::string first_failure;
    std::vector<std::string> failures;  // every failure of the pass, in reading order
  };
  auto CheckFaces = [&](const std::vector<Piece> &pieces, const std::vector<Piece> &shell)
  {
    FaceCheck check;
    auto Fail = [&](int f, const std::string &why)
    {
      check.failing[static_cast<std::size_t>(f)] = true;
      if (check.first_failure.empty())
      {
        check.first_failure = why;
      }
      check.failures.push_back(why);
    };
    auto Band = [&](double value, double threshold)
    {
      if (std::abs(value / threshold - 1.0) <= band)
      {
        check.band_hits++;
      }
    };
    auto Clearance = [&](double value)
    {
      check.min_clearance = std::min(check.min_clearance.value_or(value), value);
      Band(value, clearance);
    };
    struct Crossing
    {
      double along;
      double gap_along;  // the gap direction's component along the face
    };
    std::array<std::vector<Crossing>, 4> crossings;
    for (const Piece &piece : pieces)
    {
      const CurvePiece geometry = PieceCurve(piece, piece.s.first, piece.s.second);
      const Point2 a = Local(geometry.a), b = Local(geometry.b);
      std::optional<ArcLocal> arc;
      if (geometry.arc)
      {
        arc = ArcInFrame(*geometry.arc);
      }
      // Tangent (direction of travel a -> b) and gap direction at an end.
      auto EndTangent = [&](bool at_a) -> Point2
      {
        if (geometry.arc)
        {
          const double theta = at_a ? geometry.arc->theta0 : geometry.arc->theta1;
          return LocalDirection(geometry.arc->Tangent(theta));
        }
        const Point2 d = {b[0] - a[0], b[1] - a[1]};
        const double norm = std::hypot(d[0], d[1]);
        return norm > 0.0 ? Point2{d[0] / norm, d[1] / norm} : Point2{0.0, 0.0};
      };
      auto EndGap = [&](bool at_a) -> Point2
      {
        if (geometry.arc)
        {
          const double theta = at_a ? geometry.arc->theta0 : geometry.arc->theta1;
          const Point3D radial = geometry.arc->Radial(theta);
          const double sign = Dot(PieceGapDirection(piece),
                                  geometry.arc->Radial(0.5 * (geometry.arc->theta0 +
                                                              geometry.arc->theta1))) >= 0.0
                                  ? 1.0
                                  : -1.0;
          return LocalDirection(Scale(sign, radial));
        }
        return LocalDirection(PieceGapDirection(piece));
      };
      for (const bool at_a : {true, false})
      {
        const Point2 &e = at_a ? a : b;
        const bool cut_end = at_a ? piece.cut_a : piece.cut_b;
        std::array<bool, 4> on_face{};
        for (int f = 0; f < 4; f++)
        {
          const double distance = std::abs(e[FaceAxis(f)] - FaceLine(f));
          on_face[static_cast<std::size_t>(f)] = AtMostSnap(distance);
          Band(distance, snap);
        }
        for (int f = 0; f < 4; f++)
        {
          const double distance = std::abs(e[FaceAxis(f)] - FaceLine(f));
          if (on_face[static_cast<std::size_t>(f)])
          {
            // A crossing: its angle to the face and its position along the face.
            const Point2 tangent = EndTangent(at_a);
            const double sine = std::abs(tangent[FaceAxis(f)]);
            check.min_crossing_sine =
                std::min(check.min_crossing_sine.value_or(sine), sine);
            Band(sine, clearance);
            if (sine + 1.0e-12 < clearance)
            {
              Fail(f, "a device edge crosses face " + std::to_string(f) +
                          " at sin(theta) = " + std::to_string(sine));
            }
            const int other = 1 - FaceAxis(f);
            const Point2 gap = EndGap(at_a);
            crossings[static_cast<std::size_t>(f)].push_back({e[other], gap[other]});
            check.crossings[static_cast<std::size_t>(f)]++;
            continue;
          }
          if (cut_end)
          {
            continue;  // another cluster's claim cut: no device vertex here
          }
          // An interior device vertex, a claim end or a crossing of another face: its
          // clearance from this face (a crossing within the clearance of an adjacent face
          // is a crossing near a box corner).
          Clearance(distance);
          if (BelowClearance(distance))
          {
            Fail(f, "a device vertex or face crossing lies " + std::to_string(distance) +
                        " R from face " + std::to_string(f));
          }
        }
      }
      // The piece's clearance from the faces neither end touches. A straight piece is
      // closest to a face at an end; an end made only by another cluster's claim cut reads
      // nothing (the edge continues into the neighbouring piece, which reads its own ends).
      for (int f = 0; f < 4; f++)
      {
        const int axis = FaceAxis(f);
        const bool touches = AtMostSnap(std::abs(a[axis] - FaceLine(f))) ||
                             AtMostSnap(std::abs(b[axis] - FaceLine(f)));
        if (touches)
        {
          continue;
        }
        double distance = std::numeric_limits<double>::infinity();
        if (!piece.cut_a)
        {
          distance = std::min(distance, std::abs(a[axis] - FaceLine(f)));
        }
        if (!piece.cut_b)
        {
          distance = std::min(distance, std::abs(b[axis] - FaceLine(f)));
        }
        if (arc)
        {
          // The coordinate's extrema inside the angular range.
          const double amplitude = axis == 0 ? arc->amplitude_x : arc->amplitude_y;
          const double phi = axis == 0 ? arc->phi_x : arc->phi_y;
          const double ta = std::min(geometry.arc->theta0, geometry.arc->theta1),
                       tb = std::max(geometry.arc->theta0, geometry.arc->theta1);
          for (const double extremum : {phi, phi + pi})
          {
            const double base =
                extremum - 2.0 * pi * std::floor((extremum - ta) / (2.0 * pi));
            for (double theta = base; theta <= tb; theta += 2.0 * pi)
            {
              if (theta >= ta)
              {
                const double value =
                    arc->center[axis] + arc->radius * amplitude * std::cos(theta - phi);
                distance = std::min(distance, std::abs(value - FaceLine(f)));
              }
            }
          }
        }
        if (!std::isfinite(distance))
        {
          continue;  // a straight piece between two claim cuts: nothing to read
        }
        Clearance(distance);
        if (BelowClearance(distance))
        {
          Fail(f, "a device edge runs " + std::to_string(distance) + " R from face " +
                      std::to_string(f) + " without crossing it");
        }
      }
    }
    // Two-sided T2 (decision 285 (2)): the plan inside the shell between a face and the
    // face moved outward by the clearance. A run end there (a device vertex or joint) that
    // lies on neither boundary, or a piece that enters and leaves the shell without
    // crossing the face, sits within the clearance of the face trace without being
    // representable on it: the face fails and grows over it (T3). The exterior tail of a
    // face crossing (from the face to the shell boundary) is the crossing itself.
    {
      std::map<std::pair<std::size_t, std::size_t>, std::vector<Interval>> inside_by_run;
      for (const Piece &piece : pieces)
      {
        inside_by_run[piece.Key()].push_back(piece.s);
      }
      std::array<double, 4> dilated = box;
      for (int f = 0; f < 4; f++)
      {
        dilated[static_cast<std::size_t>(f)] += (f < 2 ? -1.0 : 1.0) * clearance;
      }
      auto OnBoundary = [&](const Point2 &e, const std::array<double, 4> &b)
      {
        for (int f = 0; f < 4; f++)
        {
          const int axis = FaceAxis(f), other = 1 - axis;
          if (AtMostSnap(std::abs(e[axis] - b[static_cast<std::size_t>(f)])) &&
              AtMostSnap(b[static_cast<std::size_t>(other)] - e[other]) &&
              AtMostSnap(e[other] - b[static_cast<std::size_t>(other + 2)]))
          {
            return true;
          }
        }
        return false;
      };
      for (const Piece &piece : shell)
      {
        std::vector<Interval> exterior = {piece.s};
        if (auto it = inside_by_run.find(piece.Key()); it != inside_by_run.end())
        {
          exterior = SubtractIntervals(exterior, it->second, Tol());
        }
        for (const auto &interval : exterior)
        {
          Piece part = piece;
          part.s = interval;
          part.claimed_by.reset();
          if (IsSliverLength(PieceLength(part)))
          {
            continue;
          }
          const CurvePiece geometry = PieceCurve(part, interval.first, interval.second);
          const Point2 a = Local(geometry.a), b = Local(geometry.b),
                       m = Local(geometry.At(0.5));
          if (Inside(m, snap + half_quantum))
          {
            continue;  // an arc's exterior part read inside: nothing outside this face
          }
          std::vector<int> beyond;
          for (int f = 0; f < 4; f++)
          {
            const double side = (f < 2 ? -1.0 : 1.0) * (m[FaceAxis(f)] - FaceLine(f));
            if (!AtMostSnap(side))
            {
              beyond.push_back(f);
            }
          }
          if (beyond.empty())
          {
            continue;
          }
          bool interior_end = false, box_end = false;
          double vertex_distance = std::numeric_limits<double>::infinity();
          for (const Point2 &e : {a, b})
          {
            if (OnBoundary(e, box))
            {
              box_end = true;
            }
            else if (!OnBoundary(e, dilated))
            {
              interior_end = true;
              for (const int f : beyond)
              {
                vertex_distance =
                    std::min(vertex_distance, std::abs(e[FaceAxis(f)] - FaceLine(f)));
              }
            }
          }
          if (interior_end)
          {
            Clearance(vertex_distance);
            check.exterior_vertices++;
            for (const int f : beyond)
            {
              Fail(f, "a device vertex lies " + std::to_string(vertex_distance) +
                          " R outside face " + std::to_string(f));
            }
          }
          else if (!box_end)
          {
            double edge_distance = std::numeric_limits<double>::infinity();
            for (const int f : beyond)
            {
              for (const Point2 &e : {a, b, m})
              {
                edge_distance =
                    std::min(edge_distance, std::abs(e[FaceAxis(f)] - FaceLine(f)));
              }
            }
            Clearance(edge_distance);
            check.exterior_edges++;
            for (const int f : beyond)
            {
              Fail(f, "a device edge runs " + std::to_string(edge_distance) +
                          " R outside face " + std::to_string(f) + " without crossing it");
            }
          }
        }
      }
    }
    // Two crossings of one face closer than the clearance bound metal (a narrow lead:
    // allowed, recorded) or gap (a channel the trace basis cannot resolve: the face fails).
    for (int f = 0; f < 4; f++)
    {
      auto &list = crossings[static_cast<std::size_t>(f)];
      std::sort(list.begin(), list.end(),
                [](const Crossing &p, const Crossing &q) { return p.along < q.along; });
      for (std::size_t k = 0; k + 1 < list.size(); k++)
      {
        const double separation = list[k + 1].along - list[k].along;
        Band(separation, clearance);
        if (!BelowClearance(separation))
        {
          continue;
        }
        const bool metal_between = list[k].gap_along < 0.0 && list[k + 1].gap_along > 0.0;
        if (metal_between)
        {
          check.min_cross_section =
              std::min(check.min_cross_section.value_or(separation), separation);
          check.narrow_cross_sections++;
          continue;
        }
        Fail(f, "two crossings of face " + std::to_string(f) + " bound a gap of " +
                    std::to_string(separation) + " R");
      }
    }
    return check;
  };

  // T3: grow every failing face by one step until every face passes or a cap refuses.
  // `steps` counts the applied steps per face, `attempted` every step a failing face asked
  // for (a refused step is attempted, not applied; review MINOR-2).
  std::array<int, 4> steps{}, attempted{};
  std::vector<std::string> step_reasons;  // the first failure behind every applied step
  // Every T2 pass (decision 287 (a)): its failing faces, every failure, and its band /
  // exterior readings — the final pass is the last entry; the counters below sum them.
  nlohmann::json passes = nlohmann::json::array();
  std::size_t total_band_hits = 0, total_exterior_vertices = 0, total_exterior_edges = 0;
  std::vector<Piece> pieces;
  FaceCheck check;
  std::optional<std::string> unboxable, unboxable_reason;
  for (;;)
  {
    pieces = ClipPlan();
    std::array<double, 4> dilated = box;
    for (int f = 0; f < 4; f++)
    {
      dilated[static_cast<std::size_t>(f)] += (f < 2 ? -1.0 : 1.0) * clearance;
    }
    const std::vector<Piece> shell = ClipPlanTo(dilated, false);
    const std::size_t box_slivers = slivers_dropped;
    check = CheckFaces(pieces, shell);
    slivers_dropped = box_slivers;
    passes.push_back({{"Box", BoxJson(box)},
                      {"Failing", check.failing},
                      {"Failures", check.failures},
                      {"ThresholdBandHits", check.band_hits},
                      {"ExteriorVertices", check.exterior_vertices},
                      {"ExteriorEdges", check.exterior_edges}});
    total_band_hits += check.band_hits;
    total_exterior_vertices += check.exterior_vertices;
    total_exterior_edges += check.exterior_edges;
    if (std::none_of(check.failing.begin(), check.failing.end(), [](bool b) { return b; }))
    {
      break;
    }
    std::array<double, 4> grown = box;
    for (int f = 0; f < 4; f++)
    {
      if (!check.failing[static_cast<std::size_t>(f)])
      {
        continue;
      }
      attempted[static_cast<std::size_t>(f)]++;
      if (steps[static_cast<std::size_t>(f)] >= kSupportFaceGrowthMaxSteps)
      {
        unboxable = "face " + std::to_string(f) + " failed the clearance rule after " +
                    std::to_string(kSupportFaceGrowthMaxSteps) +
                    " growth steps: " + check.first_failure;
        unboxable_reason = "GrowthStepsExhausted";
        break;
      }
      grown[static_cast<std::size_t>(f)] += (f < 2 ? -1.0 : 1.0) * step;
    }
    if (unboxable)
    {
      break;
    }
    const double span = std::max(grown[2] - grown[0], grown[3] - grown[1]);
    if (span > span_cap * (1.0 + 1.0e-12))
    {
      unboxable = "the growth needed by the clearance rule takes the plan span to " +
                  std::to_string(span) + " R beyond the span cap " +
                  std::to_string(span_cap) + " R: " + check.first_failure;
      unboxable_reason = "SpanCapRefusedGrowth";
      break;
    }
    for (int f = 0; f < 4; f++)
    {
      if (check.failing[static_cast<std::size_t>(f)])
      {
        steps[static_cast<std::size_t>(f)]++;
      }
    }
    step_reasons.push_back(check.first_failure);
    for (double &v : grown)
    {
      v = Q(v);
    }
    box = grown;
    out.grown = true;
  }
  const double span = std::max(box[2] - box[0], box[3] - box[1]);
  out.exceeds_span_cap = span > span_cap * (1.0 + 1.0e-12);

  // The context pieces in the portion encoding (arc pieces of one arc merged as the
  // claims are, never across another cluster's claim cut), classed chain (connected to
  // the claims inside the box through run ends and device vertices; a piece on another
  // cluster's claims stops the chain and never joins it, decision 285 (1)) or foreign.
  std::vector<SignaturePortion> context_portions;
  std::vector<std::optional<std::size_t>> context_claimed_by;
  std::vector<bool> context_excluded;
  std::size_t excluded_pieces = 0;
  double excluded_length = 0.0;
  nlohmann::json excluded_entries = nlohmann::json::array();
  {
    std::map<std::optional<std::size_t>, std::vector<std::pair<std::size_t, Interval>>>
        groups;
    for (const Piece &piece : pieces)
    {
      if (piece.Excluded())
      {
        const auto &segment = input.segments[piece.segment];
        SignaturePortion portion;
        portion.p0 = SegmentAt(piece.segment, piece.s.first);
        portion.p1 = SegmentAt(piece.segment, piece.s.second);
        portion.gap_direction = ExcludedGapDirection(piece.segment);
        portion.conductor = segment.conductor;
        portion.interfaces = InterfaceNames(segment.targets);
        portion.boundary_law = segment.boundary_law;
        const Point2 c0 = Local(portion.p0), c1 = Local(portion.p1);
        const auto &exclusion =
            segment.exclusion ? *segment.exclusion : *segment_exclusion[piece.segment];
        excluded_entries.push_back({{"P", {Q(c0[0]), Q(c0[1]), Q(c1[0]), Q(c1[1])}},
                                    {"LengthOverR", Q(PieceLength(piece) / R)},
                                    {"Class", exclusion.first},
                                    {"Reason", exclusion.second}});
        excluded_pieces++;
        excluded_length += PieceLength(piece);
        context_portions.push_back(std::move(portion));
        context_claimed_by.push_back(std::nullopt);
        context_excluded.push_back(true);
        continue;
      }
      groups[piece.claimed_by].emplace_back(piece.run, piece.s);
    }
    for (const auto &[owner, intervals] : groups)
    {
      for (auto &portion : ClusterSignaturePortions(intervals))
      {
        context_portions.push_back(std::move(portion));
        context_claimed_by.push_back(owner);
        context_excluded.push_back(false);
      }
    }
  }
  const double join = kSignatureParameterToleranceOverRadius * R;
  std::vector<Point3D> chain_ends;
  for (const auto &portion : portions)
  {
    chain_ends.push_back(portion.p0);
    chain_ends.push_back(portion.p1);
  }
  std::vector<bool> is_chain(context_portions.size(), false);
  for (bool changed = true; changed;)
  {
    changed = false;
    for (std::size_t k = 0; k < context_portions.size(); k++)
    {
      if (is_chain[k] || context_claimed_by[k] || context_excluded[k])
      {
        continue;
      }
      const auto &portion = context_portions[k];
      const bool touches = std::any_of(
          chain_ends.begin(), chain_ends.end(), [&](const Point3D &e)
          { return Distance(e, portion.p0) <= join || Distance(e, portion.p1) <= join; });
      if (touches)
      {
        is_chain[k] = true;
        chain_ends.push_back(portion.p0);
        chain_ends.push_back(portion.p1);
        changed = true;
      }
    }
  }
  FrameSupport support;
  support.box = box;
  std::set<int> claim_conductors, foreign_conductors;
  for (const auto &portion : portions)
  {
    claim_conductors.insert(portion.conductor);
  }
  std::size_t chain_pieces = 0, foreign_pieces = 0, other_claimed_pieces = 0;
  double other_claimed_length = 0.0;
  nlohmann::json claimed_by_other = nlohmann::json::array();
  for (std::size_t k = 0; k < context_portions.size(); k++)
  {
    const auto &portion = context_portions[k];
    const double length =
        portion.arc ? SignatureArcLength(*portion.arc) : Distance(portion.p0, portion.p1);
    if (is_chain[k])
    {
      chain_pieces++;
      out.chain_length += length;
    }
    else
    {
      foreign_pieces++;
      out.foreign_length += length;
      if (!claim_conductors.count(portion.conductor))
      {
        foreign_conductors.insert(portion.conductor);
      }
      if (context_claimed_by[k])
      {
        other_claimed_pieces++;
        other_claimed_length += length;
        const Point2 c0 = Local(portion.p0), c1 = Local(portion.p1);
        claimed_by_other.push_back({{"P", {Q(c0[0]), Q(c0[1]), Q(c1[0]), Q(c1[1])}},
                                    {"LengthOverR", Q(length / R)},
                                    {"Cluster", *context_claimed_by[k]}});
      }
    }
    support.context.push_back({portion, is_chain[k]});
  }

  // Vertex features inside the box that are not members of the cluster: on a chain piece
  // end (the device corner the chain turns at: owned by the coupon at placement, rule B4)
  // or foreign.
  nlohmann::json chain_vertices = nlohmann::json::array(),
                 foreign_vertices = nlohmann::json::array();
  {
    std::set<std::size_t> candidate_sites;
    for (const Piece &piece : pieces)
    {
      if (piece.Excluded())
      {
        continue;
      }
      for (auto it = sites_by_run.lower_bound(piece.run);
           it != sites_by_run.end() && it->first == piece.run; ++it)
      {
        candidate_sites.insert(it->second);
      }
    }
    for (const std::size_t s : candidate_sites)
    {
      const auto &site = sites[s];
      if (site.cluster == static_cast<int>(c))
      {
        continue;
      }
      const Point2 p = Local(site.point);
      if (!Inside(p, snap))
      {
        continue;
      }
      const bool on_chain =
          std::any_of(chain_ends.begin(), chain_ends.end(),
                      [&](const Point3D &e) { return Distance(e, site.point) <= join; });
      // With the vertex's distance from the nearest face (units of R): a corner closer
      // than R to a face has an arm partly outside the box (review MINOR-1 (a); the
      // placement's vertex-ownership rule reads it).
      double face_distance = std::numeric_limits<double>::infinity();
      for (int f = 0; f < 4; f++)
      {
        face_distance = std::min(face_distance, std::abs(p[FaceAxis(f)] - FaceLine(f)));
      }
      (on_chain ? chain_vertices : foreign_vertices)
          .push_back({{"P", {Q(p[0]), Q(p[1])}},
                      {"Type", site.type},
                      {"FaceDistanceOverR", Q(face_distance)},
                      {"Site", s},
                      {"Cluster", site.cluster >= 0 ? nlohmann::json(site.cluster)
                                                    : nlohmann::json(nullptr)},
                      {"Feature", nullptr}});
    }
  }

  // The legacy contract's census of the same feature (decision 236: every claim-cut end
  // continued straight to the claims-derived box face): the straight continuation length
  // and the part of it lying on no device edge (fictitious metal boundary).
  double continuation_length = 0.0;
  {
    std::vector<Point3D> ends;
    for (const auto &portion : portions)
    {
      ends.push_back(portion.p0);
      ends.push_back(portion.p1);
    }
    for (const auto &vertex : vertices)
    {
      ends.push_back(vertex.point);
    }
    const double coincidence = kSupportEndCoincidenceOverRadius * R;
    for (std::size_t k = 0; k < portions.size(); k++)
    {
      const auto &portion = portions[k];
      if (portion.arc)
      {
        continue;
      }
      for (const bool at_p0 : {true, false})
      {
        const Point3D &end = at_p0 ? portion.p0 : portion.p1;
        bool connected = false;
        for (std::size_t j = 0; j < ends.size() && !connected; j++)
        {
          const bool own = j / 2 == k && j < 2 * portions.size();
          connected = !own && Distance(ends[j], end) <= coincidence;
        }
        if (connected)
        {
          continue;
        }
        const Point2 e = Local(end);
        const Point2 d = LocalDirection(at_p0 ? Sub(portion.p0, portion.p1)
                                              : Sub(portion.p1, portion.p0));
        // Exit of the ray e + t d from the claims-derived box.
        double exit = std::numeric_limits<double>::infinity();
        for (int axis = 0; axis < 2; axis++)
        {
          if (std::abs(d[axis]) > 1.0e-15)
          {
            exit = std::min(
                exit,
                std::max((claims_box[static_cast<std::size_t>(axis)] - e[axis]) / d[axis],
                         (claims_box[static_cast<std::size_t>(axis + 2)] - e[axis]) /
                             d[axis]));
          }
        }
        if (!std::isfinite(exit) || exit <= 0.0)
        {
          continue;
        }
        continuation_length += exit * R;
        // Covered by straight device pieces collinear with the ray.
        std::vector<Interval> covered;
        for (const auto &context : context_portions)
        {
          if (context.arc)
          {
            continue;
          }
          const Point2 c0 = Local(context.p0), c1 = Local(context.p1);
          const Point2 dir = {c1[0] - c0[0], c1[1] - c0[1]};
          const double norm = std::hypot(dir[0], dir[1]);
          if (norm <= 0.0 || std::abs(dir[0] * d[1] - dir[1] * d[0]) / norm > 1.0e-6)
          {
            continue;
          }
          const double off = (c0[0] - e[0]) * d[1] - (c0[1] - e[1]) * d[0];
          if (std::abs(off) > snap)
          {
            continue;
          }
          double t0 = (c0[0] - e[0]) * d[0] + (c0[1] - e[1]) * d[1];
          double t1 = (c1[0] - e[0]) * d[0] + (c1[1] - e[1]) * d[1];
          if (t0 > t1)
          {
            std::swap(t0, t1);
          }
          t0 = std::max(t0, 0.0);
          t1 = std::min(t1, exit);
          if (t1 > t0)
          {
            covered.emplace_back(t0, t1);
          }
        }
        double covered_length = 0.0;
        for (const auto &interval : MergeIntervals(std::move(covered), 0.0))
        {
          covered_length += interval.second - interval.first;
        }
        out.fictitious_continuation_length += std::max(exit - covered_length, 0.0) * R;
      }
    }
  }

  // Legacy equivalence (ruling on R0 review MAJOR-1, option (i)): the context is "empty"
  // in the sense of the key rule when every context piece is what the decision-236
  // contract already draws — a straight chain piece abutting a claim-cut end, collinear
  // with that claim and reaching a face of the box — so that the device plan clipped to
  // the box is the legacy coupon geometry and the claims-only key stands.
  out.legacy_equivalent = !out.grown;
  {
    std::vector<Point3D> ends;
    for (const auto &portion : portions)
    {
      ends.push_back(portion.p0);
      ends.push_back(portion.p1);
    }
    for (const auto &vertex : vertices)
    {
      ends.push_back(vertex.point);
    }
    const double coincidence = kSupportEndCoincidenceOverRadius * R;
    auto FreeEnd = [&](std::size_t k, const Point3D &end)
    {
      for (std::size_t j = 0; j < ends.size(); j++)
      {
        const bool own = j < 2 * portions.size() && j / 2 == k;
        if (!own && Distance(ends[j], end) <= coincidence)
        {
          return false;
        }
      }
      return true;
    };
    auto OnFace = [&](const Point2 &p)
    {
      for (int f = 0; f < 4; f++)
      {
        if (std::abs(p[FaceAxis(f)] - FaceLine(f)) <= snap)
        {
          return true;
        }
      }
      return false;
    };
    for (std::size_t k = 0; k < context_portions.size() && out.legacy_equivalent; k++)
    {
      const auto &context = context_portions[k];
      bool implied = false;
      if (!context.arc && is_chain[k])
      {
        const Point2 c0 = Local(context.p0), c1 = Local(context.p1);
        const Point2 dir = LocalDirection(Sub(context.p1, context.p0));
        for (std::size_t j = 0; j < portions.size() && !implied; j++)
        {
          const auto &claim = portions[j];
          if (claim.arc)
          {
            continue;
          }
          const Point2 t = LocalDirection(Sub(claim.p1, claim.p0));
          if (std::abs(dir[0] * t[1] - dir[1] * t[0]) > 1.0e-6)
          {
            continue;  // not parallel to this claim
          }
          for (const bool at_p0 : {true, false})
          {
            const Point3D &end = at_p0 ? claim.p0 : claim.p1;
            if (!FreeEnd(j, end))
            {
              continue;
            }
            const Point2 e = Local(end);
            const bool abuts_0 = Distance(context.p0, end) <= join,
                       abuts_1 = Distance(context.p1, end) <= join;
            if (!abuts_0 && !abuts_1)
            {
              continue;
            }
            const Point2 &far = abuts_0 ? c1 : c0;
            const double off = (far[0] - e[0]) * t[1] - (far[1] - e[1]) * t[0];
            if (std::abs(off) <= snap && OnFace(far))
            {
              implied = true;
              break;
            }
          }
        }
      }
      out.legacy_equivalent = implied;
    }
  }

  // Truncation segments (window cuts) inside the box: the plan is cut there (the excluded
  // perimeter segments inside the box are context pieces, Context.ExcludedSegments: review
  // MINOR-5, rule B1).
  std::size_t truncation_segments = 0;
  double truncation_length = 0.0;
  const Point3D plane_normal = Normalize(Cross(x, y));
  for (std::size_t i = 0; i < input.segments.size(); i++)
  {
    const auto &segment = input.segments[i];
    if (!segment.truncation)
    {
      continue;
    }
    const Point2 p = Local(segment.p0), q = Local(segment.p1);
    if (std::max(p[0], q[0]) < box[0] || std::min(p[0], q[0]) > box[2] ||
        std::max(p[1], q[1]) < box[1] || std::min(p[1], q[1]) > box[3])
    {
      continue;
    }
    // Plane test: the segment's distance along the cluster normal.
    if (std::abs(Dot(Sub(segment.p0, origin), plane_normal)) > 2.0 * Tol() ||
        std::abs(Dot(Sub(segment.p1, origin), plane_normal)) > 2.0 * Tol())
    {
      continue;
    }
    double t0 = 0.0, t1 = 1.0;
    bool empty = false;
    for (int axis = 0; axis < 2 && !empty; axis++)
    {
      const double d = q[axis] - p[axis];
      const double lo_a = box[static_cast<std::size_t>(axis)],
                   hi_a = box[static_cast<std::size_t>(axis + 2)];
      if (std::abs(d) <= 1.0e-15)
      {
        empty = p[axis] < lo_a || p[axis] > hi_a;
        continue;
      }
      double ta = (lo_a - p[axis]) / d, tb = (hi_a - p[axis]) / d;
      if (ta > tb)
      {
        std::swap(ta, tb);
      }
      t0 = std::max(t0, ta);
      t1 = std::min(t1, tb);
      empty = t0 >= t1;
    }
    if (!empty)
    {
      truncation_segments++;
      truncation_length += (t1 - t0) * Distance(segment.p0, segment.p1);
    }
  }

  out.band_hits = total_band_hits + box_band_hits;
  nlohmann::json growth = {{"StepOverR", step},
                           {"MaxSteps", kSupportFaceGrowthMaxSteps},
                           {"Steps", steps},
                           {"AttemptedSteps", attempted},
                           {"StepReasons", step_reasons},
                           {"Passes", passes},
                           {"Grown", out.grown}};
  nlohmann::json face_rules = {
      {"SnapOverR", snap},
      {"ClearanceOverR", clearance},
      {"ComparisonQuantumOverR", kSignatureLengthQuantumOverRadius},
      {"SliversDropped", slivers_dropped},
      {"Crossings", check.crossings},
      {"MinClearanceOverR", check.min_clearance ? nlohmann::json(Q(*check.min_clearance))
                                                : nlohmann::json(nullptr)},
      {"MinCrossingSine", check.min_crossing_sine
                              ? nlohmann::json(Q(*check.min_crossing_sine))
                              : nlohmann::json(nullptr)},
      {"MinCrossSectionOverR", check.min_cross_section
                                   ? nlohmann::json(Q(*check.min_cross_section))
                                   : nlohmann::json(nullptr)},
      {"NarrowCrossSections", check.narrow_cross_sections},
      {"ExteriorVertices", total_exterior_vertices},
      {"ExteriorEdges", total_exterior_edges},
      {"ThresholdBandRelative", band},
      {"ThresholdBandHits", total_band_hits},
      {"FinalPassThresholdBandHits", check.band_hits},
      {"Passes", passes.size()},
      {"BoxRuleThresholdBandHits", box_band_hits}};
  out.record = {
      {"ClaimsBox", BoxJson(claims_box)},
      {"Box", BoxJson(box)},
      {"SpanOverR", Q(span)},
      {"SpanCapOverR", span_cap},
      {"ExceedsSpanCap", out.exceeds_span_cap},
      {"LegacyEquivalent", out.legacy_equivalent},
      {"Growth", growth},
      {"FaceRules", face_rules},
      {"Context",
       {{"Pieces", context_portions.size()},
        {"ChainPieces", chain_pieces},
        {"ForeignPieces", foreign_pieces},
        {"ChainLengthOverR", Q(out.chain_length / R)},
        {"ForeignLengthOverR", Q(out.foreign_length / R)},
        {"ForeignConductors", foreign_conductors.size()},
        {"ClaimedByOtherFeature",
         {{"Pieces", other_claimed_pieces},
          {"LengthOverR", Q(other_claimed_length / R)},
          {"Entries", claimed_by_other}}},
        {"ExcludedSegments",
         {{"Pieces", excluded_pieces},
          {"LengthOverR", Q(excluded_length / R)},
          {"Entries", excluded_entries}}},
        {"ChainVertices", chain_vertices},
        {"ForeignVertices", foreign_vertices}}},
      {"LegacyContinuation",
       {{"StraightContinuationLengthOverR", Q(continuation_length / R)},
        {"FictitiousContinuationLengthOverR", Q(out.fictitious_continuation_length / R)}}},
      {"Truncation",
       {{"Segments", truncation_segments}, {"LengthOverR", Q(truncation_length / R)}}},
      {"Unboxable", unboxable ? nlohmann::json(*unboxable) : nlohmann::json(nullptr)},
      {"UnboxableReason",
       unboxable_reason ? nlohmann::json(*unboxable_reason) : nlohmann::json(nullptr)}};
  if (!unboxable)
  {
    out.support = std::move(support);
  }
  return out;
}

// The non-hashed ownership ids of a cluster feature's context census (review MINOR-1):
// every ChainVertices / ForeignVertices entry names its vertex feature (a free site) or the
// cluster whose claims hold it; every ClaimedByOtherFeature entry names the other cluster's
// feature. Called after EmitClusters (cluster ids) and after Assign (site ids).
void Identifier::RecordContextFeatureIds(int feature)
{
  if (feature < 0 || static_cast<std::size_t>(feature) >= features.size())
  {
    return;
  }
  nlohmann::json &record = features[static_cast<std::size_t>(feature)].spatial_support;
  if (!record.is_object() || !record.contains("Context"))
  {
    return;
  }
  auto ClusterFeature = [&](const nlohmann::json &cluster) -> nlohmann::json
  {
    if (!cluster.is_number_integer())
    {
      return nullptr;
    }
    const auto c = cluster.get<std::size_t>();
    return c < cluster_feature.size() && cluster_feature[c] >= 0
               ? nlohmann::json(cluster_feature[c])
               : nlohmann::json(nullptr);
  };
  nlohmann::json &context = record["Context"];
  for (const char *list : {"ChainVertices", "ForeignVertices"})
  {
    for (auto &entry : context[list])
    {
      if (entry["Cluster"].is_number_integer())
      {
        entry["Feature"] = ClusterFeature(entry["Cluster"]);
      }
      else if (entry.contains("Site") && entry["Site"].is_number_integer())
      {
        const auto s = entry["Site"].get<std::size_t>();
        entry["Feature"] = s < site_feature.size() && site_feature[s] >= 0
                               ? nlohmann::json(site_feature[s])
                               : nlohmann::json(nullptr);
      }
    }
  }
  for (auto &entry : context["ClaimedByOtherFeature"]["Entries"])
  {
    entry["Feature"] = ClusterFeature(entry["Cluster"]);
  }
}

void Identifier::EmitClusters()
{
  double signature_seconds = 0.0;
  std::size_t largest_cluster_edges = 0;
  sites_by_run.clear();
  for (std::size_t s = 0; s < sites.size(); s++)
  {
    for (const std::size_t r : sites[s].runs_at_site)
    {
      sites_by_run.emplace(r, s);
    }
  }
  claims_by_run.clear();
  for (std::size_t c = 0; c < cluster_claimed.size(); c++)
  {
    for (const auto &[r, interval] : cluster_claimed[c])
    {
      claims_by_run[r].emplace_back(interval, c);
    }
  }
  cluster_feature.assign(cluster_claimed.size(), -1);
  spatial_support_summary = {};
  for (std::size_t c = 0; c < cluster_claimed.size(); c++)
  {
    stage.Progress(c, cluster_claimed.size());
    std::vector<SignaturePortion> portions = ClusterSignaturePortions(cluster_claimed[c]);
    std::vector<SignatureVertex> vertices;
    for (const std::size_t s : cluster_sites[c])
    {
      vertices.push_back({sites[s].point, sites[s].type, sites[s].turn_degrees});
    }
    MFEM_VERIFY(!portions.empty(), "A spatial cluster claims no perimeter!");
    largest_cluster_edges = std::max(largest_cluster_edges, portions.size());
    // The cluster's own signed process normal (claimed-length weighted; decision 266), not
    // the sign-canonical n_ref: on a flipped plane the canonical frame's w must point into
    // that plane's vacuum.
    Point3D n_cluster;
    {
      std::vector<std::pair<std::size_t, double>> weighted;
      for (const auto &[r, interval] : cluster_claimed[c])
      {
        weighted.emplace_back(r, interval.second - interval.first);
      }
      n_cluster =
          FeatureProcessNormal(weighted, "spatial cluster " + std::to_string(c) + " (" +
                                             std::to_string(portions.size()) + " edges)");
    }
    const auto signature_started = std::chrono::steady_clock::now();
    const auto canonical = CanonicalClusterSignature(
        portions, vertices, n_cluster, R,
        [&](std::size_t done, std::size_t total)
        {
          stage.Progress(done, total,
                         "candidate frames of cluster " + std::to_string(c) + " / " +
                             std::to_string(cluster_claimed.size()) + ", " +
                             std::to_string(portions.size()) + " edges");
        });
    // Spatial-support contract v3 (decision 282): the support in the claims-only canonical
    // frame decides the contract — an empty context without growth keeps the claims-only
    // key (today's coupons keep their keys); otherwise Box + Context enter the hashed
    // signature, minimised over the candidate frames with the support of each frame; a
    // cluster no frame can box is an UnboxableFeature (a Missing placeholder key).
    // The claims-only signature (contract 2 key; the key a legacy model carries; the key
    // the span-cap allowances are resolved by, BEFORE any box: block (b) DESIGN A4).
    nlohmann::json claims_signature = canonical.signature;
    claims_signature["EdgeCount"] = portions.size();
    claims_signature["Type"] = "SpatialEdgeCluster";
    double span_cap_over_R = kSupportSpanCapOverRadius;
    nlohmann::json span_cap_allowance;
    if (const auto resolved =
            ResolveSpanCapAllowance(claims_signature, input.span_cap_allowances))
    {
      const SpanCapAllowance &allowance = input.span_cap_allowances[resolved->index];
      span_cap_over_R = allowance.span_cap_over_R;
      used_span_cap_allowances.insert(resolved->index);
      span_cap_allowance = {
          {"Label", allowance.label},
          {"SpanCapOverR", allowance.span_cap_over_R},
          {"Reason", allowance.reason},
          {"Approval", allowance.approval},
          {"MatchedQuanta", resolved->difference.max_delta_quanta},
          {"DifferingNumbers",
           {{"Count", resolved->difference.differing_paths.size()},
            {"Paths", resolved->difference.differing_paths}}},
          {"Rule", "block (b) DESIGN section 3 (a) / A4 (decision 303): the allowance's "
                   "embedded claims-only signature near-matches this cluster's claims-only "
                   "signature within ClusterQuantumNearMatchMaxQuanta quanta (half-quantum "
                   "inclusive); its "
                   "SpanCapOverR replaces the plan span cap for this cluster's box growth "
                   "and ExceedsSpanCap reading (never below the default cap)"}};
    }
    ClusterSupportResult support =
        ClusterSupport(c, portions, vertices, canonical.origin, canonical.axes[0],
                       canonical.axes[1], span_cap_over_R);
    nlohmann::json signature;
    int chirality = canonical.chirality;
    Point3D origin = canonical.origin;
    std::array<Point3D, 3> axes = canonical.axes;
    int contract = 2;
    if (support.support && support.legacy_equivalent)
    {
      signature = canonical.signature;
    }
    else
    {
      const auto v3 = CanonicalClusterSignatureWithSupport(
          portions, vertices, n_cluster, R,
          [&](const Point3D &o, const Point3D &x,
              const Point3D &y) -> std::optional<FrameSupport>
          {
            ClusterSupportResult frame_support =
                ClusterSupport(c, portions, vertices, o, x, y, span_cap_over_R);
            return frame_support.support;
          },
          [&](std::size_t done, std::size_t total)
          {
            stage.Progress(done, total,
                           "support frames of cluster " + std::to_string(c) + " / " +
                               std::to_string(cluster_claimed.size()));
          });
      if (v3)
      {
        contract = 3;
        signature = v3->signature;
        chirality = v3->chirality;
        origin = v3->origin;
        axes = v3->axes;
        // The record of the frame the key was formed in.
        support = ClusterSupport(c, portions, vertices, origin, axes[0], axes[1],
                                 span_cap_over_R);
      }
      else
      {
        contract = 0;
        signature = canonical.signature;
        signature["Unboxable"] = true;
      }
    }
    signature_seconds +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - signature_started)
            .count();
    signature["EdgeCount"] = portions.size();
    const int feature = NewFeature("SpatialEdgeCluster", signature, chirality);
    cluster_feature[c] = feature;
    features[feature].origin = origin;
    features[feature].axes = axes;
    features[feature].claims_origin = canonical.origin;
    features[feature].claims_axes = canonical.axes;
    features[feature].claims_chirality = canonical.chirality;
    support.record["Contract"] = contract;
    // The claims-only (contract-2) key of the cluster: the key a legacy model built under
    // the decision-236 contract carries; equal to Hash for contract 2. A legacy-contract
    // alias (USER decision 283) is verified against it.
    support.record["ClaimsKey"] =
        SignatureKeyAndHash(claims_signature, "SpatialEdgeCluster").second;
    if (!span_cap_allowance.is_null())
    {
      support.record["SpanCapAllowance"] = span_cap_allowance;
    }
    if (contract == 0)
    {
      // A refused cluster (SpanCapRefusedGrowth / UnboxableFeature) exports its claims-only
      // signature verbatim: the object an operator copies into a SpanCapAllowances entry
      // (block (b) DESIGN A4 (2)).
      support.record["ClaimsSignature"] = claims_signature;
    }
    support.record["ContextDigest"] =
        contract == 3
            ? nlohmann::json(SpatialSupportContextDigest(features[feature].signature))
            : nlohmann::json(nullptr);
    features[feature].spatial_support = support.record;
    spatial_support_summary.clusters++;
    spatial_support_summary.claims_keyed += contract == 2 ? 1 : 0;
    spatial_support_summary.context_keyed += contract == 3 ? 1 : 0;
    spatial_support_summary.unboxable += contract == 0 ? 1 : 0;
    spatial_support_summary.grown += support.grown ? 1 : 0;
    spatial_support_summary.exceeding_span_cap += support.exceeds_span_cap ? 1 : 0;
    spatial_support_summary.chain_length += support.chain_length;
    spatial_support_summary.foreign_length += support.foreign_length;
    spatial_support_summary.fictitious_continuation_length +=
        support.fictitious_continuation_length;
    spatial_support_summary.threshold_band_hits += support.band_hits;
    if (input.log && contract != 2)
    {
      std::ostringstream text;
      text << "  Identification spatial support: cluster " << c << " (" << portions.size()
           << " edges) " << (contract == 3 ? "keyed with Box + Context" : "UNBOXABLE")
           << ": context pieces "
           << (support.support ? support.support->context.size() : std::size_t(0))
           << ", chain " << support.chain_length / R << " R, foreign "
           << support.foreign_length / R << " R, grown " << support.grown
           << ", fictitious straight continuation "
           << support.fictitious_continuation_length / R << " R, threshold-band readings "
           << support.band_hits << " (every T2 pass + the box rule)\n";
      if (support.band_hits > 0)
      {
        text << "  WARNING: cluster " << c << "'s key rests on " << support.band_hits
             << " face-rule reading(s) within " << kKnifeEdgeBandRelative
             << " of a threshold (decision 287); see "
                "Features[].SpatialSupport.FaceRules.ThresholdBandHits and "
                "Growth.Passes\n";
      }
      input.log(text.str());
    }
    for (const auto &[r, interval] : cluster_claimed[c])
    {
      claims[r].push_back({feature, 0, interval});
    }
    for (const std::size_t s : cluster_sites[c])
    {
      if (sites[s].vertex)
      {
        vertex_feature[*sites[s].vertex] = feature;
        features[feature].vertices.push_back(*sites[s].vertex);
      }
    }
  }
  // The feature ids of the other clusters in the context census (known only now that
  // every cluster has its feature; the free sites' ids follow in Assign).
  for (std::size_t c = 0; c < cluster_claimed.size(); c++)
  {
    RecordContextFeatureIds(cluster_feature[c]);
  }
  {
    std::ostringstream counts;
    counts << cluster_claimed.size() << " clusters (largest " << largest_cluster_edges
           << " edges), canonical signatures " << std::fixed << std::setprecision(2)
           << signature_seconds << " s, " << features.size() << " features so far";
    stage.End(counts.str());
  }
}

// Knife-edge census (rule at kKnifeEdgeBandRelative): distances sampled along every
// non-excluded run to the nearest other run (3D, any plane; the same chain only beyond the
// self-pair neighbourhood of arc length; adjacent runs sharing a vertex excluded), the
// windowed bend radius along every chain, the geometric turn of every vertex of two runs.
nlohmann::json Identifier::KnifeEdgeCensus() const
{
  const double band = kKnifeEdgeBandRelative * R;
  const double spacing = kKnifeEdgeSampleSpacingOverRadius * R;
  const double neighbourhood = kSelfPairNeighbourhoodOverRadius * R;
  struct Band
  {
    double below = 0.0, above = 0.0;
    void Add(double value, double threshold, double width, double weight)
    {
      if (std::abs(value - threshold) < width)
      {
        (value < threshold ? below : above) += weight;
      }
    }
    nlohmann::json ToJson() const
    {
      return {{"Below", below}, {"Above", above}, {"Total", below + above}};
    }
  };
  Band at_r, at_2r, at_bend, at_noise, at_arc_turn, at_sagitta;
  double sampled = 0.0;
  const double reach = 2.0 * R + band + 2.0 * Tol();
  for (std::size_t a = 0; a < runs.size(); a++)
  {
    const Run &ra = runs[a];
    if (ra.excluded || ra.length <= 0.0)
    {
      continue;
    }
    const Chain &A = chains[chain_index.at(ra.chain)];
    const double xa0 = A.run_offset[ra.index_in_chain];
    const int n = std::max(1, static_cast<int>(std::ceil(ra.length / spacing)));
    const double weight = ra.length / n;
    const auto candidates = RunsNearRun(a, reach);
    for (int i = 0; i < n; i++)
    {
      const double s = ra.length * (i + 0.5) / n;
      const Point3D p = ra.At(s);
      // Every other run within the reach is an interaction candidate of the sample: the
      // sample counts once per threshold when any of them lies within the band (the far
      // ground of a 2 / 2 / 2 um stack is 2R from the trace edge whose nearest neighbour is
      // R away).
      std::optional<double> near_r, near_2r;
      for (const std::size_t b : candidates)
      {
        const Run &rb = runs[b];
        if (b == a || rb.excluded || rb.start_vertex == ra.start_vertex ||
            rb.start_vertex == ra.end_vertex || rb.end_vertex == ra.start_vertex ||
            rb.end_vertex == ra.end_vertex)
        {
          continue;
        }
        if (rb.chain == ra.chain)
        {
          const double xb0 = A.run_offset[rb.index_in_chain];
          if (ChainArcDistance(A, xb0, xb0 + rb.length, xa0 + s) < neighbourhood)
          {
            continue;
          }
        }
        const double t = std::clamp(Dot(Sub(p, rb.start), rb.tangent), 0.0, rb.length);
        const double d = Distance(p, rb.At(t));
        if (std::abs(d - R) < band && (!near_r || std::abs(d - R) < std::abs(*near_r - R)))
        {
          near_r = d;
        }
        if (std::abs(d - 2.0 * R) < band &&
            (!near_2r || std::abs(d - 2.0 * R) < std::abs(*near_2r - 2.0 * R)))
        {
          near_2r = d;
        }
      }
      sampled += weight;
      if (near_r)
      {
        at_r.Add(*near_r, R, band, weight);
      }
      if (near_2r)
      {
        at_2r.Add(*near_2r, 2.0 * R, band, weight);
      }
      const double kappa = WindowedCurvature(A, xa0 + s);
      if (kappa > 0.0)
      {
        at_bend.Add(1.0 / kappa, kStraightBendRadiusOverRadius * R,
                    kKnifeEdgeBandRelative * kStraightBendRadiusOverRadius * R, weight);
      }
    }
  }
  const double noise_sagitta = kJointNoiseSagittaOverRadius * R;
  const double arc_turn_cap = kArcMaxJointTurnDegrees;
  const double sagitta_cap = kArcSagittaOverRadius * R;
  // Threshold keys of the census ("50", "0.05"): the shortest decimal of the constant.
  auto ThresholdKey = [](double value)
  {
    std::ostringstream key;
    key << std::setprecision(12) << value;
    return key.str();
  };
  for (std::size_t v = 0; v < input.vertices.size(); v++)
  {
    const auto incident = RunsAtVertex(v);
    if (incident.size() != 2)
    {
      continue;
    }
    const Point3D ta = ArmDirection(incident[0], v), tb = ArmDirection(incident[1], v);
    const double turn =
        180.0 - std::acos(std::clamp(Dot(ta, tb), -1.0, 1.0)) * 180.0 / std::acos(-1.0);
    if (turn > 0.0)
    {
      // The joint noise rule: the sagitta the turn implies on the SHORTER run (the rule's
      // quantity) against JointNoiseSagittaOverR x R; the arc rule's joint-turn cap; and
      // the mesh-coarseness diagnostic: the LONGER run read as a chord of one circle at the
      // joint's turn (an inscribed polyline's chord of central angle t has sagitta
      // (c / 2) tan(t / 4); a fitted arc's chords are its runs between joints).
      const double shorter = std::min(runs[incident[0]].length, runs[incident[1]].length);
      const double longer = std::max(runs[incident[0]].length, runs[incident[1]].length);
      const double tan_quarter = std::tan(0.25 * turn * std::acos(-1.0) / 180.0);
      at_noise.Add(0.5 * shorter * tan_quarter, noise_sagitta,
                   kKnifeEdgeBandRelative * noise_sagitta, 1.0);
      at_arc_turn.Add(turn, arc_turn_cap, kKnifeEdgeBandRelative * arc_turn_cap, 1.0);
      at_sagitta.Add(0.5 * longer * tan_quarter, sagitta_cap,
                     kKnifeEdgeBandRelative * sagitta_cap, 1.0);
    }
  }
  return {
      {"BandRelative", kKnifeEdgeBandRelative},
      {"SampleSpacingOverR", kKnifeEdgeSampleSpacingOverRadius},
      {"Rule", "perimeter length with another perimeter point (3D; the same chain beyond "
               "the self-pair neighbourhood) at a distance within the band of the "
               "threshold, split into the below / above sides; the chain length "
               "whose windowed bend radius lies within the band of the straight-bend "
               "radius; the two-run vertices whose implied sagitta (c / 2) tan(turn / 4) "
               "on the SHORTER run lies within the band of JointNoiseSagittaOverR x R "
               "(the joint noise rule), whose turn lies within the band of "
               "ArcMaxJointTurnDegrees (the arc rule's cap), and whose LONGER run read "
               "as a chord at the joint's turn has a sagitta within the band of "
               "SagittaOverR x R (the mesh-coarseness diagnostic) (counts, not lengths)"},
      {"SampledLength", sampled},
      {"Distance", {{"R", at_r.ToJson()}, {"2R", at_2r.ToJson()}}},
      {"BendRadius", {{"10R", at_bend.ToJson()}}},
      {"JointNoiseSagittaOverR",
       {{ThresholdKey(kJointNoiseSagittaOverRadius), at_noise.ToJson()}}},
      {"ArcMaxJointTurnDegrees", {{ThresholdKey(arc_turn_cap), at_arc_turn.ToJson()}}},
      {"ArcSagittaOverR", {{ThresholdKey(kArcSagittaOverRadius), at_sagitta.ToJson()}}}};
}

// Cluster-composition band of the knife-edge census (USER decision 184 (4), 2026-10-01;
// the E6 finding of stage 0: the DS-CTX-003 loop-end clusters went from 42 to 48 edges at
// R 1.85 while every distance band stayed quiet). A cluster's membership is the outcome of
// the whole cluster machinery (events, cores, joins, claims, the extension), not of one
// distance, so the band is measured the way the stage-0 study measured it: the same
// perimeter (the input as received: segments, vertices, faces; the vertex classification
// of the perimeter extraction at R is kept) is identified again at R (1 - band) and at
// R (1 + band), silently, and every cluster at R is matched to the cluster at the other
// radius sharing the most claimed perimeter with it. Its composition is unchanged iff the
// match has the same edge count (Signature.EdgeCount) and the same set of member vertices
// (corner / endpoint / junction mesh vertices); otherwise the cluster's claimed length
// counts on the Below (R (1 - band)) or Above (R (1 + band)) side, with the number of
// clusters affected. Reported like the other bands (Below / Above / Total in length units)
// and in the identification log; costs two more identifications on the root rank.
nlohmann::json Identifier::ClusterCompositionBand(const IdentificationResult &result) const
{
  struct ClusterRecord
  {
    int edge_count = 0;
    std::vector<std::size_t> vertices;
    double length = 0.0;
    std::map<std::size_t, std::vector<Interval>> portions;  // per segment
  };
  auto Clusters = [](const IdentificationResult &r)
  {
    std::vector<ClusterRecord> records;
    for (const auto &feature : r.features)
    {
      if (feature.type != "SpatialEdgeCluster")
      {
        continue;
      }
      ClusterRecord record;
      record.edge_count = feature.signature.value("EdgeCount", 0);
      record.vertices = feature.vertices;
      std::sort(record.vertices.begin(), record.vertices.end());
      record.length = feature.length;
      for (const auto &portion : feature.portions)
      {
        record.portions[portion.segment].emplace_back(std::min(portion.s0, portion.s1),
                                                      std::max(portion.s0, portion.s1));
      }
      records.push_back(std::move(record));
    }
    return records;
  };
  const std::vector<ClusterRecord> at_R = Clusters(result);
  nlohmann::json band = {{"BandRelative", kKnifeEdgeBandRelative},
                         {"Clusters", at_R.size()},
                         {"Rule",
                          "the perimeter identified again at R (1 - BandRelative) and at "
                          "R (1 + BandRelative); every SpatialEdgeCluster at R is matched "
                          "to the cluster at the other radius sharing the most claimed "
                          "perimeter with it and is unchanged iff the match has the same "
                          "EdgeCount and the same member vertices; the claimed length (at "
                          "R) of the changed clusters on the Below / Above side, and their "
                          "counts (ClustersBelow / ClustersAbove of Clusters)"}};
  double below = 0.0, above = 0.0;
  std::size_t clusters_below = 0, clusters_above = 0;
  for (const double factor : {1.0 - kKnifeEdgeBandRelative, 1.0 + kKnifeEdgeBandRelative})
  {
    Identifier other(original_input, R * factor);
    const IdentificationResult other_result = other.Identify();
    const std::vector<ClusterRecord> at_other = Clusters(other_result);
    std::map<std::size_t, std::vector<std::size_t>> by_segment;  // segment -> clusters
    for (std::size_t c = 0; c < at_other.size(); c++)
    {
      for (const auto &[segment, intervals] : at_other[c].portions)
      {
        (void)intervals;
        by_segment[segment].push_back(c);
      }
    }
    for (const ClusterRecord &cluster : at_R)
    {
      std::map<std::size_t, double> shared;
      for (const auto &[segment, intervals] : cluster.portions)
      {
        const auto it = by_segment.find(segment);
        if (it == by_segment.end())
        {
          continue;
        }
        for (const std::size_t c : it->second)
        {
          for (const auto &[a0, a1] : intervals)
          {
            for (const auto &[b0, b1] : at_other[c].portions.at(segment))
            {
              shared[c] += std::max(0.0, std::min(a1, b1) - std::max(a0, b0));
            }
          }
        }
      }
      const ClusterRecord *match = nullptr;
      double best = 0.0;
      for (const auto &[c, length] : shared)
      {
        if (length > best)
        {
          best = length;
          match = &at_other[c];
        }
      }
      const bool unchanged = match && match->edge_count == cluster.edge_count &&
                             match->vertices == cluster.vertices;
      if (!unchanged)
      {
        (factor < 1.0 ? below : above) += cluster.length;
        (factor < 1.0 ? clusters_below : clusters_above)++;
      }
    }
  }
  band["Below"] = below;
  band["Above"] = above;
  band["Total"] = below + above;
  band["ClustersBelow"] = clusters_below;
  band["ClustersAbove"] = clusters_above;
  return band;
}

void Identifier::Assign(IdentificationResult &result)
{
  site_feature.assign(sites.size(), -1);
  // Vertex features outside clusters.
  for (auto &site : sites)
  {
    if (site.cluster >= 0)
    {
      continue;
    }
    nlohmann::json signature = {{"Interfaces", site.interfaces},
                                {"Law", site.boundary_law}};
    // Frame of the vertex feature (the patch frame of the library model, design (b) 4):
    // x = the first arm away from the site, y in the process plane, z = the site's OWN
    // signed process normal n (substrate -> vacuum; decision 266 — not the sign-canonical
    // n_ref, which turns a flipped plane's coupon upside down). Corner: the arms ordered so
    // that the second is counterclockwise about n from the first (a corner is its own
    // mirror image), y = n x x (right-handed); endpoint: y toward the gap; junction: x =
    // the canonical first arm, y = +-(n x x) so that the canonical arm order proceeds
    // counterclockwise in (x, y) (the arm angles are measured about the same n).
    const Point3D n_site = SiteProcessNormal(site);
    std::array<Point3D, 3> axes = {Point3D{}, Point3D{}, n_site};
    if (site.type == "ConvexCorner" || site.type == "ConcaveCorner")
    {
      signature = CanonicalCornerSignature(site.interfaces, site.boundary_law,
                                           site.angle_degrees, site.corner_radius / R);
      Point3D first = site.arm_directions[0], second = site.arm_directions[1];
      if (Dot(Cross(n_site, first), second) < 0.0)
      {
        std::swap(first, second);
      }
      axes[0] = Normalize(Sub(first, Scale(Dot(first, n_site), n_site)));
      axes[1] = Normalize(Cross(n_site, axes[0]));
    }
    else if (site.type == "Junction")
    {
      JunctionCanonicalOrder order;
      signature = CanonicalJunctionSignature(site.interfaces, site.boundary_law,
                                             site.arm_angles, site.arm_conductors, &order);
      const Point3D first = site.arm_directions[order.first_arm];
      axes[0] = Normalize(Sub(first, Scale(Dot(first, n_site), n_site)));
      axes[1] = Scale(order.reversed ? -1.0 : 1.0, Normalize(Cross(n_site, axes[0])));
    }
    else if (!site.arm_directions.empty())
    {
      const Point3D first = site.arm_directions.front();
      axes[0] = Normalize(Sub(first, Scale(Dot(first, n_site), n_site)));
      axes[1] = Normalize(Cross(n_site, axes[0]));
      if (Dot(axes[1], site.gap_direction) < 0.0)
      {
        axes[1] = Scale(-1.0, axes[1]);
      }
    }
    const int feature = NewFeature(site.type, signature);
    site_feature[static_cast<std::size_t>(&site - sites.data())] = feature;
    features[feature].origin = site.point;
    features[feature].axes = axes;
    for (const auto &[r, interval] : site.window)
    {
      claims[r].push_back({feature, 1, interval});
    }
    if (site.vertex)
    {
      vertex_feature[*site.vertex] = feature;
      features[feature].vertices.push_back(*site.vertex);
    }
  }
  for (const int feature : cluster_feature)
  {
    RecordContextFeatureIds(feature);
  }

  // Resolve the claims by priority and hand the remainder to the chain's isolated edge.
  std::map<std::pair<int, std::string>, int> isolated_features;
  std::vector<std::vector<std::tuple<double, double, int, int>>> assigned(runs.size());
  for (std::size_t r = 0; r < runs.size(); r++)
  {
    if (runs[r].excluded)
    {
      continue;
    }
    stage.Progress(r, runs.size(), "claims");
    auto &run_claims = claims[r];
    // Claims by priority, then by position along the run (never by feature id, review m1):
    // claims of one priority never overlap (clusters disjoint, windows shortened to abut,
    // pairs / stacks assembled per cross-section; Diagnostics.SamePriorityClaimOverlaps is
    // gated at 0 by the audit), so their order within a priority is immaterial to the
    // assignment; were one to overlap, the earlier interval along the run would win.
    std::sort(run_claims.begin(), run_claims.end(),
              [](const Claim &a, const Claim &b)
              {
                return std::tie(a.priority, a.interval, a.feature) <
                       std::tie(b.priority, b.interval, b.feature);
              });
    // Two claims of one priority by different features overlapping on the run: counted.
    for (std::size_t i = 0; i < run_claims.size(); i++)
    {
      for (std::size_t j = i + 1;
           j < run_claims.size() && run_claims[j].priority == run_claims[i].priority; j++)
      {
        if (run_claims[i].feature != run_claims[j].feature &&
            std::min(run_claims[i].interval.second, run_claims[j].interval.second) -
                    std::max(run_claims[i].interval.first, run_claims[j].interval.first) >
                Tol())
        {
          same_priority_claim_overlaps++;
        }
      }
    }
    // CrossLayer zones are excluded before any feature claims the run.
    std::vector<Interval> taken = cross_layer[r];
    for (const auto &claim : run_claims)
    {
      for (const auto &piece : SubtractIntervals({claim.interval}, taken, Tol()))
      {
        assigned[r].emplace_back(piece.first, piece.second, claim.feature, claim.side);
        taken.push_back(piece);
      }
      taken = MergeIntervals(taken, Tol());
    }
    if (const char *debug_segment = std::getenv("PALACE_IDENTIFICATION_DEBUG_SEGMENT");
        debug_segment && input.log)
    {
      const auto wanted = static_cast<std::size_t>(std::atol(debug_segment));
      // A negative value: every run whose sorted claims of one feature leave a gap longer
      // than 0.1 R between them.
      bool gap = false;
      if (std::atol(debug_segment) < 0)
      {
        std::vector<std::tuple<double, double, int>> pieces;
        for (const auto &claim : run_claims)
        {
          pieces.emplace_back(claim.interval.first, claim.interval.second, claim.feature);
        }
        std::sort(pieces.begin(), pieces.end());
        for (std::size_t i = 1; i < pieces.size(); i++)
        {
          if (std::get<2>(pieces[i]) == std::get<2>(pieces[i - 1]) &&
              std::get<0>(pieces[i]) - std::get<1>(pieces[i - 1]) > 0.1 * R)
          {
            gap = true;
          }
        }
      }
      if (gap || std::any_of(runs[r].segments.begin(), runs[r].segments.end(),
                             [&](const RunSegment &rs) { return rs.segment == wanted; }))
      {
        std::ostringstream dbg;
        dbg << "  DEBUG run " << r << " (segment " << wanted << ") length "
            << runs[r].length << " from (" << runs[r].start[0] << ", " << runs[r].start[1]
            << ") to (" << runs[r].end[0] << ", " << runs[r].end[1] << ") claims:\n";
        for (const auto &claim : run_claims)
        {
          dbg << "    feature " << claim.feature << " priority " << claim.priority << " ["
              << claim.interval.first << ", " << claim.interval.second << "] side "
              << claim.side << " " << features[static_cast<std::size_t>(claim.feature)].type
              << "\n";
        }
        for (const auto &rs : runs[r].segments)
        {
          dbg << "    segment " << rs.segment << " t [" << rs.t0 << ", " << rs.t1 << "]\n";
        }
        input.log(dbg.str());
      }
    }
    // Pieces shorter than the signature parameter tolerance between the boundaries of two
    // claims (a cluster ball cutting a pair piece, a claim ending next to a run end, the
    // foot of a neighbour's cut end) are joined to their adjacent portion by the sliver
    // rule below, once every portion of the chain is known (decision 222).
    std::sort(assigned[r].begin(), assigned[r].end());
  }

  // A pair / parallel cluster whose claim on one of its sides was lost entirely to an
  // earlier claim of the same priority (a chain facing two partners within 2R along a bend:
  // DS-SCT-001's 2 um gaps next to 2 um strips, a case the bent-pair rule describes per
  // chain pair) is not a pair: its surviving pieces return to their runs and are handed to
  // the chain's isolated / curved edge below (recorded as a claim-resolution rule; the
  // multi-partner bent neighbourhood itself is a PENDING rule item).
  {
    // A pair's two sides face each other: the part of one side whose facing point lies
    // beyond the other side's surviving pieces (a cluster ball or a vertex window took the
    // partner piece but not this one) is not paired there and returns to the run; the
    // facing test is the pair tolerance on the separation — the LOCAL separation of the
    // piece (its ends and midpoint to the partner chain), never below the feature's mean:
    // a slow taper's wider end is still facing (the pair rule admitted it on its local
    // reading below 2R), while the mean alone cut every pair whose separation grows by more
    // than the tolerance along its length (DS-CTX-003's flux-line launcher, 1 -> 6 um over
    // 560 um: the inner gap's pair vanished where the outer links ended and the four edges
    // read isolated at 3.7 um; review m11). A constant pair reads its mean everywhere.
    // Recorded as claim resolution: both sides of a pair are then mutual.
    struct SidePiece
    {
      std::size_t run;
      double lo, hi;
    };
    std::map<int, std::array<std::vector<SidePiece>, 2>> pair_pieces;
    for (std::size_t r = 0; r < runs.size(); r++)
    {
      for (const auto &[lo, hi, feature, side] : assigned[r])
      {
        if (feature >= 0 && feature_sides[feature] == 2 && (side == 0 || side == 1))
        {
          pair_pieces[feature][side].push_back({r, lo, hi});
        }
      }
    }
    std::size_t pair_index = 0;
    for (const auto &[feature, sides] : pair_pieces)
    {
      stage.Progress(pair_index++, pair_pieces.size(), "mutual pair sides");
      if (!features[feature].signature.contains("SeparationOverR"))
      {
        continue;
      }
      const double mean_separation =
          features[feature].signature["SeparationOverR"].get<double>() * R;
      // The other side's pieces by run: a piece faces only the pieces on runs whose boxes
      // lie within the reach of its own run (the facing intervals of every other piece are
      // empty; the union is sorted, so the candidate order does not matter).
      std::array<std::map<std::size_t, std::vector<const SidePiece *>>, 2> pieces_by_run;
      std::array<std::set<int>, 2> side_chains;
      for (int k = 0; k < 2; k++)
      {
        for (const auto &piece : sides[k])
        {
          pieces_by_run[k][piece.run].push_back(&piece);
          side_chains[k].insert(runs[piece.run].chain);
        }
      }
      // The local separation of a side piece: the largest distance from its ends and
      // midpoint to the partner side's chains (within the pair rule's candidate reach).
      const double candidate_reach =
          kInteractionDistanceOverRadius * R * (1.0 + kPairSeparationTolerance);
      auto LocalSeparation = [&](const SidePiece &piece, int k)
      {
        double local = 0.0;
        for (const double s : {piece.lo, 0.5 * (piece.lo + piece.hi), piece.hi})
        {
          const Point3D p = runs[piece.run].At(s);
          double nearest = std::numeric_limits<double>::infinity();
          for (const int other_chain : side_chains[1 - k])
          {
            const auto foot = ClosestPointOnChain(chains[chain_index.at(other_chain)], p,
                                                  std::nullopt, candidate_reach);
            if (std::isfinite(foot.distance))
            {
              nearest = std::min(nearest, foot.distance);
            }
          }
          if (std::isfinite(nearest))
          {
            local = std::max(local, nearest);
          }
        }
        return local;
      };
      for (int k = 0; k < 2; k++)
      {
        for (const auto &piece : sides[k])
        {
          const double reach = std::max(mean_separation, LocalSeparation(piece, k)) *
                               (1.0 + kPairSeparationTolerance);
          std::vector<Interval> facing;
          for (const std::size_t other_run : RunsNearRun(piece.run, reach + 2.0 * Tol()))
          {
            const auto it = pieces_by_run[1 - k].find(other_run);
            if (it == pieces_by_run[1 - k].end())
            {
              continue;
            }
            for (const SidePiece *other : it->second)
            {
              const auto found =
                  RunIntervalWithin(piece.run, runs[other->run].At(other->lo),
                                    runs[other->run].At(other->hi), reach);
              facing.insert(facing.end(), found.begin(), found.end());
            }
          }
          const auto keep = IntersectIntervals({Interval{piece.lo, piece.hi}},
                                               MergeIntervals(facing, Tol()), Tol());
          auto &pieces = assigned[piece.run];
          const auto it = std::find_if(pieces.begin(), pieces.end(),
                                       [&](const auto &entry)
                                       {
                                         return std::get<2>(entry) == feature &&
                                                std::get<3>(entry) == k &&
                                                std::get<0>(entry) == piece.lo &&
                                                std::get<1>(entry) == piece.hi;
                                       });
          if (it == pieces.end())
          {
            continue;
          }
          pieces.erase(it);
          for (const auto &interval : keep)
          {
            if (interval.second - interval.first > kSignatureLengthQuantumOverRadius * R)
            {
              pieces.emplace_back(interval.first, interval.second, feature, k);
            }
          }
        }
      }
    }
    for (auto &pieces : assigned)
    {
      std::sort(pieces.begin(), pieces.end());
    }

    std::map<int, std::map<int, double>> side_lengths;
    for (std::size_t r = 0; r < runs.size(); r++)
    {
      for (const auto &[lo, hi, feature, side] : assigned[r])
      {
        if (feature >= 0 && feature_sides[feature] > 0)
        {
          side_lengths[feature][side] += hi - lo;
        }
      }
    }
    std::set<int> degenerate;
    for (int feature = 0; feature < static_cast<int>(features.size()); feature++)
    {
      if (feature_sides[feature] == 0)
      {
        continue;
      }
      const auto &lengths = side_lengths[feature];
      for (int side = 0; side < feature_sides[feature]; side++)
      {
        const auto it = lengths.find(side);
        // A side whose surviving length is below the signature parameter tolerance is no
        // side (its pieces would join their neighbours under the sliver rule, decision
        // 222).
        if (it == lengths.end() || it->second < kSignatureParameterToleranceOverRadius * R)
        {
          degenerate.insert(feature);
          break;
        }
      }
    }
    if (!degenerate.empty())
    {
      for (auto &pieces : assigned)
      {
        pieces.erase(std::remove_if(pieces.begin(), pieces.end(), [&](const auto &piece)
                                    { return degenerate.count(std::get<2>(piece)) > 0; }),
                     pieces.end());
      }
    }
  }

  for (std::size_t r = 0; r < runs.size(); r++)
  {
    if (runs[r].excluded)
    {
      continue;
    }
    stage.Progress(r, runs.size(), "isolated / curved edges");
    std::vector<Interval> taken = cross_layer[r];
    for (const auto &[lo, hi, feature, side] : assigned[r])
    {
      (void)feature;
      (void)side;
      taken.push_back({lo, hi});
    }
    taken = MergeIntervals(taken, Tol());
    const auto remainder = SubtractIntervals({Interval{0.0, runs[r].length}}, taken, Tol());
    if (!remainder.empty())
    {
      // Curved sections of the chain (windowed bend radius below the straight threshold)
      // are CurvedEdge features, one per section; the rest is the chain's isolated edge.
      const Chain &chain = chains[chain_index.at(runs[r].chain)];
      const double offset = chain.run_offset[RunIndexInChain(chain, r)];
      std::vector<Interval> curved_on_run;
      std::vector<std::size_t> curved_section;
      // The chain's curved intervals (disjoint, ascending) that can overlap the run: from
      // the first ending after the run's start (less R) to the last starting before its end
      // (plus R); the exact intersection below is the former one.
      const auto first_curved = std::lower_bound(
          chain.curved.begin(), chain.curved.end(), offset - R,
          [](const Interval &c, double value) { return c.second <= value; });
      for (auto ci = first_curved;
           ci != chain.curved.end() && ci->first < offset + runs[r].length + R; ++ci)
      {
        const auto j = static_cast<std::size_t>(ci - chain.curved.begin());
        const Interval local{chain.curved[j].first - offset,
                             chain.curved[j].second - offset};
        const auto overlap = IntersectIntervals({local}, remainder, Tol());
        for (const auto &piece : overlap)
        {
          curved_on_run.push_back(piece);
          curved_section.push_back(j);
        }
      }
      for (std::size_t i = 0; i < curved_on_run.size(); i++)
      {
        nlohmann::json signature = {{"Interfaces", InterfaceNames(runs[r].targets)},
                                    {"Law", runs[r].boundary_law}};
        signature["RadiusOverR"] =
            RoundTo(1.0 / (chain.curved_max_kappa[curved_section[i]] * R),
                    kSignatureLengthQuantumOverRadius);
        {
          // Convexity of the section (metal inside / outside the bend relative to the gap
          // direction): the curvature family of that convexity models it.
          const Interval &section = chain.curved[curved_section[i]];
          SignedCurvatureExtremes extremes;
          AccumulateSignedCurvature(chain, section.first, section.second, extremes);
          const std::string convexity = ConvexityName(extremes);
          MFEM_VERIFY(!convexity.empty(), "A curved section without signed curvature!");
          signature["Convexity"] = convexity;
        }
        signature["Type"] = "CurvedEdge";
        // One feature per curved section of the chain — or per fitted bend arc when the
        // section lies on one (the arc's parts on either side of an absorbed corner vertex
        // are one feature; key by the smallest overlapping arc).
        int section_arc = -1;
        for (const auto &[a, x0, x1] : chain.arc_spans)
        {
          const Interval &section = chain.curved[curved_section[i]];
          if (std::min(section.second, x1) - std::max(section.first, x0) > Tol() &&
              (section_arc < 0 || a < section_arc))
          {
            section_arc = a;
          }
        }
        const auto key =
            section_arc >= 0
                ? std::make_pair(-1,
                                 signature.dump() + "#arc" + std::to_string(section_arc))
                : std::make_pair(ChainGroup(runs[r].chain),
                                 signature.dump() + "#" + std::to_string(runs[r].chain) +
                                     "/" + std::to_string(curved_section[i]));
        auto it = isolated_features.find(key);
        if (it == isolated_features.end())
        {
          signature.erase("Type");
          const int feature = NewFeature("CurvedEdge", signature);
          features[feature].origin = runs[r].At(curved_on_run[i].first);
          // Representative frame (informative: the patches use the segment frames): the
          // run's own signed process normal, right-handed (decision 266).
          const Point3D n_run = SignedReferenceNormal(runs[r].process_normal,
                                                      "isolated run " + std::to_string(r));
          features[feature].axes = {runs[r].tangent, Cross(n_run, runs[r].tangent), n_run};
          it = isolated_features.emplace(key, feature).first;
        }
        assigned[r].emplace_back(curved_on_run[i].first, curved_on_run[i].second,
                                 it->second, 0);
      }
      const auto straight = SubtractIntervals(remainder, curved_on_run, Tol());
      if (!straight.empty())
      {
        nlohmann::json signature = {{"Interfaces", InterfaceNames(runs[r].targets)},
                                    {"Law", runs[r].boundary_law}};
        signature["Type"] = "IsolatedEdge";
        const auto key = std::make_pair(ChainGroup(runs[r].chain), signature.dump());
        auto it = isolated_features.find(key);
        if (it == isolated_features.end())
        {
          signature.erase("Type");
          const int feature = NewFeature("IsolatedEdge", signature);
          features[feature].origin = runs[r].start;
          // Representative frame (informative: the patches use the segment frames): the
          // run's own signed process normal, right-handed (decision 266).
          const Point3D n_run = SignedReferenceNormal(runs[r].process_normal,
                                                      "isolated run " + std::to_string(r));
          features[feature].axes = {runs[r].tangent, Cross(n_run, runs[r].tangent), n_run};
          it = isolated_features.emplace(key, feature).first;
        }
        for (const auto &piece : straight)
        {
          assigned[r].emplace_back(piece.first, piece.second, it->second, 0);
        }
      }
    }
    std::sort(assigned[r].begin(), assigned[r].end());
  }

  // Sliver rule (supervisor decision 222, 2026-10-02): NO portion shorter than the
  // signature parameter tolerance exists. A portion here is a maximal contiguous stretch of
  // one feature side along a chain, across the runs and mesh segments of the chain (a
  // sub-tolerance MESH SEGMENT inside a long edge continues that edge's portion and is not
  // one). A stretch shorter than kSignatureParameterToleranceOverRadius x R is roundoff
  // between the boundaries of two claims — a cluster ball cutting a pair piece, a claim
  // ending next to a run end, the perpendicular foot of a neighbour's cut end on a
  // near-parallel member (the S1p stage-1 window: the 3-edge stack owned a 1.94e-6 um piece
  // at the y = -200 truncation cut 64 um from its stack end, 2 um x the 9.7e-7 rad tilt of
  // the window set, and its twin on the next member read as a 2e-6 um IsolatedEdge; a
  // sample placed on that piece pointed its lateral axis along the segment and the
  // placement aborted) — which no signature parameter resolves (every length is matched
  // within the tolerance) and no coupon models: it joins its adjacent portion on the chain,
  // the LONGER of its two neighbours (ties: the one before it along the chain), taking that
  // neighbour's feature and side. The former rule joined only pieces at or below the
  // signature grid (1e-6 R), a knife-edge the S1p piece missed by 2.3 % (DS-SCT-001's
  // three 1.6e-7 um CurvedSameConductorStrip pieces are the same class). A stretch with no
  // adjacent portion on its chain (a whole chain, or a piece between two CrossLayer zones,
  // shorter than the tolerance) has nothing to join and stays, counted. Deterministic and
  // numbering-independent: the chain order and the neighbours' lengths decide, never a
  // feature id. Reported as Diagnostics.SubTolerancePortionsJoined (count, length, longest,
  // isolated).
  // The stretch index of every piece (per feature, in chain order) is recorded on its
  // portions for the placement's ownership record (a translational stretch wholly inside a
  // spatial support is recorded, decisions 224 / 236).
  std::vector<std::vector<int>> piece_stretch(runs.size());
  {
    const double tolerance = kSignatureParameterToleranceOverRadius * R;
    struct ChainPiece
    {
      std::size_t run;
      std::size_t index;  // into assigned[run]
    };
    struct Stretch
    {
      std::size_t first, last;  // chain piece indices (inclusive; last < first wraps)
      int feature, side;
      double length;
    };
    auto Entry = [&](const ChainPiece &p) -> std::tuple<double, double, int, int> &
    { return assigned[p.run][p.index]; };
    auto ChainPieces = [&](const Chain &chain)
    {
      std::vector<ChainPiece> pieces;
      for (const std::size_t r : chain.runs)
      {
        if (runs[r].excluded)
        {
          continue;
        }
        for (std::size_t i = 0; i < assigned[r].size(); i++)
        {
          pieces.push_back({r, i});
        }
      }
      return pieces;
    };
    // Piece b follows piece a contiguously along the chain: on one run with touching
    // intervals, or at the joint of two consecutive non-excluded runs of the chain (a's
    // interval reaching its run's end, b's starting at its run's start); on a closed chain
    // the last run's end meets the first run's start.
    auto Adjacent = [&](const Chain &chain, const ChainPiece &a, const ChainPiece &b)
    {
      const double ahi = std::get<1>(Entry(a)), blo = std::get<0>(Entry(b));
      if (a.run == b.run)
      {
        return std::abs(ahi - blo) <= Tol();
      }
      const std::size_t ka = runs[a.run].index_in_chain, kb = runs[b.run].index_in_chain;
      const bool consecutive =
          kb == ka + 1 || (chain.closed && ka + 1 == chain.runs.size() && kb == 0);
      return consecutive && ahi >= runs[a.run].length - Tol() && blo <= Tol();
    };
    auto SameFeature = [&](const ChainPiece &a, const ChainPiece &b)
    {
      return std::get<2>(Entry(a)) == std::get<2>(Entry(b)) &&
             std::get<3>(Entry(a)) == std::get<3>(Entry(b));
    };
    auto Length = [&](const ChainPiece &p)
    { return std::get<1>(Entry(p)) - std::get<0>(Entry(p)); };
    // Maximal contiguous stretches of one feature side, in chain order; on a closed chain
    // the first and the last stretch may be one (wraps: the last piece is adjacent to the
    // first).
    auto Stretches =
        [&](const Chain &chain, const std::vector<ChainPiece> &pieces, bool &wraps)
    {
      std::vector<Stretch> stretches;
      for (std::size_t i = 0; i < pieces.size(); i++)
      {
        if (!stretches.empty() && Adjacent(chain, pieces[i - 1], pieces[i]) &&
            SameFeature(pieces[i - 1], pieces[i]))
        {
          stretches.back().last = i;
          stretches.back().length += Length(pieces[i]);
        }
        else
        {
          stretches.push_back({i, i, std::get<2>(Entry(pieces[i])),
                               std::get<3>(Entry(pieces[i])), Length(pieces[i])});
        }
      }
      wraps = chain.closed && pieces.size() > 1 &&
              Adjacent(chain, pieces.back(), pieces.front());
      if (wraps && stretches.size() > 1 && SameFeature(pieces.back(), pieces.front()))
      {
        stretches.back().last = stretches.front().last;
        stretches.back().length += stretches.front().length;
        stretches.erase(stretches.begin());
      }
      return stretches;
    };
    std::set<std::size_t> touched_runs;
    for (const Chain &chain : chains)
    {
      const std::vector<ChainPiece> pieces = ChainPieces(chain);
      if (pieces.empty())
      {
        continue;
      }
      for (bool changed = true; changed;)
      {
        changed = false;
        bool wraps = false;
        const std::vector<Stretch> stretches = Stretches(chain, pieces, wraps);
        const std::size_t n = stretches.size();
        for (std::size_t g = 0; g < n && !changed; g++)
        {
          const Stretch &stretch = stretches[g];
          if (stretch.length >= tolerance)
          {
            continue;
          }
          // The neighbours: the stretch before (its last piece adjacent to this stretch's
          // first) and after (this stretch's last piece adjacent to its first); on a closed
          // chain of two stretches both are the other one.
          const Stretch *before = nullptr, *after = nullptr;
          if (g > 0 || wraps)
          {
            const std::size_t gb = g > 0 ? g - 1 : n - 1;
            if (gb != g &&
                Adjacent(chain, pieces[stretches[gb].last], pieces[stretch.first]))
            {
              before = &stretches[gb];
            }
          }
          if (g + 1 < n || wraps)
          {
            const std::size_t ga = g + 1 < n ? g + 1 : 0;
            if (ga != g &&
                Adjacent(chain, pieces[stretch.last], pieces[stretches[ga].first]))
            {
              after = &stretches[ga];
            }
          }
          const Stretch *target = before && after
                                      ? (after->length > before->length ? after : before)
                                      : (before ? before : after);
          if (!target)
          {
            continue;  // nothing to join on this chain
          }
          // Relabel every piece of the stretch (wrapping on a closed chain).
          for (std::size_t i = stretch.first;; i = (i + 1) % pieces.size())
          {
            auto &entry = Entry(pieces[i]);
            std::get<2>(entry) = target->feature;
            std::get<3>(entry) = target->side;
            touched_runs.insert(pieces[i].run);
            if (i == stretch.last)
            {
              break;
            }
          }
          result.sub_tolerance_portions.count++;
          result.sub_tolerance_portions.length += stretch.length;
          result.sub_tolerance_portions.max_length =
              std::max(result.sub_tolerance_portions.max_length, stretch.length);
          changed = true;
        }
        if (!changed)
        {
          for (const Stretch &stretch : stretches)
          {
            if (stretch.length < tolerance)
            {
              result.sub_tolerance_portions.isolated++;
            }
          }
        }
      }
    }
    // Pieces of one feature side made contiguous by a join are one piece of the run.
    for (const std::size_t r : touched_runs)
    {
      auto &pieces = assigned[r];
      std::vector<std::tuple<double, double, int, int>> merged;
      for (const auto &piece : pieces)
      {
        if (!merged.empty() && std::get<2>(merged.back()) == std::get<2>(piece) &&
            std::get<3>(merged.back()) == std::get<3>(piece) &&
            std::abs(std::get<1>(merged.back()) - std::get<0>(piece)) <= Tol())
        {
          std::get<1>(merged.back()) = std::get<1>(piece);
        }
        else
        {
          merged.push_back(piece);
        }
      }
      pieces = std::move(merged);
    }
    // Stretch numbering: per feature, in the order of the chains and along each chain.
    std::map<int, int> stretches_of_feature;
    for (std::size_t r = 0; r < runs.size(); r++)
    {
      piece_stretch[r].assign(assigned[r].size(), -1);
    }
    for (const Chain &chain : chains)
    {
      const std::vector<ChainPiece> pieces = ChainPieces(chain);
      bool wraps = false;
      for (const Stretch &stretch : Stretches(chain, pieces, wraps))
      {
        const int index = stretches_of_feature[stretch.feature]++;
        for (std::size_t i = stretch.first;; i = (i + 1) % pieces.size())
        {
          piece_stretch[pieces[i].run][pieces[i].index] = index;
          if (i == stretch.last)
          {
            break;
          }
        }
      }
    }
  }

  // Curvature annotation of every feature over its assigned portions.
  for (std::size_t r = 0; r < runs.size(); r++)
  {
    if (runs[r].excluded)
    {
      continue;
    }
    const Chain &chain = chains[chain_index.at(runs[r].chain)];
    const double offset = chain.run_offset[RunIndexInChain(chain, r)];
    for (const auto &[lo, hi, feature, side] : assigned[r])
    {
      (void)side;
      feature_max_kappa[feature] = std::max(feature_max_kappa[feature],
                                            MaxCurvature(chain, offset + lo, offset + hi));
    }
  }

  // Fitted arcs (rounded corners and bends).
  std::vector<int> segment_any_arc(input.segments.size(), -1);
  for (std::size_t a = 0; a < arcs.size(); a++)
  {
    const Arc &arc = arcs[a];
    result.arcs.push_back({arc.center, arc.radius, arc.turn * 180.0 / std::acos(-1.0),
                           arc.max_sagitta / R, arc.corner ? "RoundedCorner" : "Bend",
                           arc.joints.size(), arc.segments.size()});
    for (const std::size_t s : arc.segments)
    {
      segment_any_arc[s] = static_cast<int>(a);
    }
  }
  // Per-segment table.
  result.segments.resize(input.segments.size());
  for (std::size_t i = 0; i < input.segments.size(); i++)
  {
    const auto &segment = input.segments[i];
    const double length = Distance(segment.p0, segment.p1);
    result.perimeter_length += length;
    result.segments[i].key =
        segment.p0 < segment.p1
            ? std::array<std::array<double, 3>, 2>{segment.p0, segment.p1}
            : std::array<std::array<double, 3>, 2>{segment.p1, segment.p0};
    result.segments[i].length = length;
    result.segments[i].chain = segment.chain;
    result.segments[i].arc = i < segment_any_arc.size() ? segment_any_arc[i] : -1;
    if (segment.truncation)
    {
      result.segments[i].exclusion =
          std::make_pair("TruncationCut", "metal perimeter on the simulation boundary");
    }
    else if (segment.exclusion)
    {
      result.segments[i].exclusion = segment.exclusion;
    }
    else if (segment_exclusion[i])
    {
      result.segments[i].exclusion = segment_exclusion[i];
    }
    else if (segment.chain < 0)
    {
      result.segments[i].exclusion = std::make_pair(
          "Unchained", "physical metal segment outside every perimeter chain");
    }
    if (result.segments[i].exclusion)
    {
      result.excluded_length += length;
    }
  }
  std::map<std::pair<std::string, std::string>, std::pair<int, double>> exclusions;
  for (std::size_t i = 0; i < input.segments.size(); i++)
  {
    if (result.segments[i].exclusion)
    {
      auto &entry = exclusions[*result.segments[i].exclusion];
      entry.first++;
      entry.second += Distance(input.segments[i].p0, input.segments[i].p1);
    }
  }
  for (const auto &[key, value] : exclusions)
  {
    result.exclusions.push_back({key.first, key.second, value.first, value.second});
  }
  // The CrossLayer record collects the analytic zones (count = zones, length = their sum).
  const int cross_layer_record = static_cast<int>(result.exclusions.size());
  result.exclusions.push_back({"CrossLayer",
                               "planar metal edge within 2R of metal off its own plane "
                               "(facing layer, wall, staple)",
                               0, 0.0});
  for (std::size_t r = 0; r < runs.size(); r++)
  {
    if (runs[r].excluded)
    {
      continue;
    }
    stage.Progress(r, runs.size(), "segment table");
    for (const auto &rs : runs[r].segments)
    {
      auto &table = result.segments[rs.segment];
      MFEM_VERIFY(!table.exclusion, "An excluded segment received a feature portion!");
      const double segment_length =
          Distance(input.segments[rs.segment].p0, input.segments[rs.segment].p1);
      const double scale = segment_length / (rs.t1 - rs.t0);
      auto SegmentPortion = [&](double lo,
                                double hi) -> std::optional<std::array<double, 2>>
      {
        const double a = std::max(lo, rs.t0), b = std::min(hi, rs.t1);
        if (b - a <= Tol())
        {
          return std::nullopt;
        }
        double s0 = (a - rs.t0) * scale, s1 = (b - rs.t0) * scale;
        if (!rs.forward)
        {
          const double t = s0;
          s0 = segment_length - s1;
          s1 = segment_length - t;
        }
        s0 = std::clamp(s0, 0.0, segment_length);
        s1 = std::clamp(s1, 0.0, segment_length);
        if (!(input.segments[rs.segment].p0 < input.segments[rs.segment].p1))
        {
          const double t = s0;
          s0 = segment_length - s1;
          s1 = segment_length - t;
        }
        return std::array<double, 2>{s0, s1};
      };
      const Chain &chain = chains[chain_index.at(runs[r].chain)];
      const double offset = chain.run_offset[RunIndexInChain(chain, r)];
      for (std::size_t i = 0; i < assigned[r].size(); i++)
      {
        const auto &[lo, hi, feature, side] = assigned[r][i];
        const auto portion = SegmentPortion(lo, hi);
        if (!portion)
        {
          continue;
        }
        const auto [s0, s1] = *portion;
        table.portions.push_back({s0, s1, static_cast<double>(feature)});
        const double turn =
            SignedTurn(chain, offset + std::max(lo, rs.t0), offset + std::min(hi, rs.t1));
        features[feature].portions.push_back(
            {rs.segment, s0, s1, side, turn, piece_stretch[r][i]});
        features[feature].length += s1 - s0;
        result.assigned_length += s1 - s0;
      }
      std::sort(table.portions.begin(), table.portions.end());
      for (const auto &zone : cross_layer[r])
      {
        const auto portion = SegmentPortion(zone.first, zone.second);
        if (!portion)
        {
          continue;
        }
        const auto [s0, s1] = *portion;
        table.excluded_portions.push_back(
            {s0, s1, static_cast<double>(cross_layer_record)});
        result.exclusions[cross_layer_record].count++;
        result.exclusions[cross_layer_record].length += s1 - s0;
        result.excluded_length += s1 - s0;
      }
      std::sort(table.excluded_portions.begin(), table.excluded_portions.end());
    }
  }
  if (result.exclusions[cross_layer_record].count == 0)
  {
    result.exclusions.erase(result.exclusions.begin() + cross_layer_record);
  }

  // Drop features that ended up claiming nothing (translational spans fully inside clusters
  // or vertex windows) and renumber the survivors so that ids are dense.
  {
    std::vector<int> renumber(features.size(), -1);
    std::vector<IdentifiedFeature> kept;
    for (auto &feature : features)
    {
      if (feature.length > Tol() || !feature.vertices.empty())
      {
        if (feature_max_kappa[feature.id] > 0.0)
        {
          feature.bend_radius_over_R = 1.0 / (feature_max_kappa[feature.id] * R);
        }
        renumber[feature.id] = static_cast<int>(kept.size());
        feature.id = static_cast<int>(kept.size());
        kept.push_back(std::move(feature));
      }
    }
    features = std::move(kept);
    for (auto &segment : result.segments)
    {
      for (auto &portion : segment.portions)
      {
        portion[2] = renumber[static_cast<int>(portion[2])];
        MFEM_VERIFY(portion[2] >= 0, "A segment portion refers to a dropped feature!");
      }
    }
    for (auto &[vertex, feature] : vertex_feature)
    {
      (void)vertex;
      feature = renumber[feature];
    }
  }

  // Vertex table. A regular vertex where a chain was cut by an excluded segment is an
  // ExclusionCut (the metal edge continues into the exclusion; no endpoint feature).
  std::map<std::size_t, int> runs_at_vertex;
  for (const auto &run : runs)
  {
    if (!run.excluded)
    {
      runs_at_vertex[run.start_vertex]++;
      runs_at_vertex[run.end_vertex]++;
    }
  }
  // The first site at every mesh vertex (the former linear search over the sites).
  std::map<std::size_t, std::size_t> site_of_vertex;
  for (std::size_t i = 0; i < sites.size(); i++)
  {
    if (sites[i].vertex)
    {
      site_of_vertex.try_emplace(*sites[i].vertex, i);
    }
  }
  // Corner features outside clusters by their origin (the last feature at a point wins, as
  // in the former loop over every feature) and, per segment, the spatial cluster features
  // with a portion on it (ascending: the former first match over the features).
  std::map<Point3D, int> corner_feature_at;
  std::map<std::size_t, std::vector<int>> cluster_features_on_segment;
  for (const auto &feature : features)
  {
    if (feature.type == "ConvexCorner" || feature.type == "ConcaveCorner")
    {
      corner_feature_at[feature.origin] = feature.id;
    }
    else if (feature.type == "SpatialEdgeCluster")
    {
      for (const auto &portion : feature.portions)
      {
        auto &list = cluster_features_on_segment[portion.segment];
        if (list.empty() || list.back() != feature.id)
        {
          list.push_back(feature.id);
        }
      }
    }
  }
  for (std::size_t v = 0; v < input.vertices.size(); v++)
  {
    const auto &vertex = input.vertices[v];
    if (!vertex.physical_type || *vertex.physical_type == MetalEdgeVertexType::REGULAR)
    {
      const auto ends = runs_at_vertex.find(v);
      if (vertex.physical_type && ends != runs_at_vertex.end() && ends->second == 1)
      {
        IdentifiedVertex entry;
        entry.vertex = v;
        entry.type = "ExclusionCut";
        entry.point_contact = IsPointContact(v);
        result.vertices.push_back(entry);
      }
      else if (vertex.mirror_joint)
      {
        // A real truncation vertex joined straight to its image chain on a mirror plane
        // (boundary-cut DESIGN 2.2.2): an interior joint, listed so that the record shows
        // where the chain continues through the plane; never a feature.
        IdentifiedVertex entry;
        entry.vertex = v;
        entry.type = "MirrorJoint";
        entry.point_contact = IsPointContact(v);
        result.vertices.push_back(entry);
      }
      continue;
    }
    IdentifiedVertex entry;
    entry.vertex = v;
    entry.point_contact = IsPointContact(v);
    if (!IsFeatureVertex(v))
    {
      entry.type =
          cross_layer_vertices.find(v) != cross_layer_vertices.end()
              ? "Excluded"
              : (vertex.image_band_cut
                     ? "ImageBandCut"
                     : (vertex.on_truncation_boundary ? "TruncationCut" : "PortCut"));
      result.vertices.push_back(entry);
      continue;
    }
    const auto site_at = site_of_vertex.find(v);
    if (site_at == site_of_vertex.end())
    {
      entry.type = "Excluded";
      result.vertices.push_back(entry);
      continue;
    }
    const VertexFeatureSite *site = &sites[site_at->second];
    entry.type = site->type;
    entry.turn_degrees = site->turn_degrees;
    const auto feature = vertex_feature.find(v);
    entry.feature = feature == vertex_feature.end() ? -1 : feature->second;
    result.vertices.push_back(entry);
  }
  std::map<std::size_t, int> rounded_site_feature;  // site index -> feature
  for (std::size_t site_index = 0; site_index < sites.size(); site_index++)
  {
    const auto &site = sites[site_index];
    if (site.vertex)
    {
      continue;
    }
    IdentifiedVertex entry;
    entry.vertex = std::numeric_limits<std::size_t>::max();
    entry.type = "RoundedCorner";
    entry.turn_degrees = site.turn_degrees;
    entry.feature = -1;
    if (site.cluster < 0)
    {
      if (const auto at = corner_feature_at.find(site.point); at != corner_feature_at.end())
      {
        entry.feature = at->second;
      }
    }
    else
    {
      // The first spatial cluster feature (by id) with a portion on a segment of the
      // site's window runs.
      for (const auto &w : site.window)
      {
        for (const RunSegment &rs : runs[w.first].segments)
        {
          const auto on = cluster_features_on_segment.find(rs.segment);
          if (on != cluster_features_on_segment.end() &&
              (entry.feature < 0 || on->second.front() < entry.feature))
          {
            entry.feature = on->second.front();
          }
        }
      }
    }
    rounded_site_feature[site_index] = entry.feature;
    result.vertices.push_back(entry);
  }
  // Corner vertices absorbed by a fitted arc (design (b) 4 / 7, decision 82(3)): members of
  // the rounded corner (RoundedCornerVertex, its feature) or of a bend (BendVertex); never
  // a corner feature of their own.
  for (const auto &arc : arcs)
  {
    for (const std::size_t v : arc.absorbed_corners)
    {
      IdentifiedVertex entry;
      entry.vertex = v;
      entry.type = arc.corner ? "RoundedCornerVertex" : "BendVertex";
      entry.turn_degrees = arc.turn * 180.0 / std::acos(-1.0);
      entry.feature = -1;
      if (arc.corner && arc.feature >= 0)
      {
        const auto it = rounded_site_feature.find(static_cast<std::size_t>(arc.feature));
        entry.feature = it == rounded_site_feature.end() ? -1 : it->second;
      }
      entry.point_contact = IsPointContact(v);
      result.vertices.push_back(entry);
    }
  }

  // Digest over the sorted feature signatures (with multiplicity) and the exclusion
  // classes. Lengths are continuous quantities (their sums differ at roundoff between
  // meshes of the same layout) and are compared with a tolerance by the audit instead of
  // being hashed.
  std::vector<std::string> lines;
  for (const auto &feature : features)
  {
    lines.push_back(feature.signature_key);
  }
  std::sort(lines.begin(), lines.end());
  for (const auto &exclusion : result.exclusions)
  {
    lines.push_back("X|" + exclusion.cls + "|" + exclusion.reason);
  }
  std::string digest_input;
  for (const auto &line : lines)
  {
    digest_input += line;
    digest_input += '\n';
  }
  result.geometry_digest = Sha256HexImpl(digest_input);
  result.same_priority_claim_overlaps = same_priority_claim_overlaps;
  result.stack_geometric_offsets = stack_geometric_offsets;
  result.stack_composition_cap_hits = stack_composition_cap_hits;
  result.stack_images_merged = stack_images_merged;
  result.spatial_support = spatial_support_summary;
  for (std::size_t a = 0; a < input.span_cap_allowances.size(); a++)
  {
    if (!used_span_cap_allowances.count(a))
    {
      result.unused_span_cap_allowances.push_back(input.span_cap_allowances[a].label);
    }
  }
  if (input.log && !result.unused_span_cap_allowances.empty())
  {
    std::ostringstream text;
    text << "  Identification spatial support: " << result.unused_span_cap_allowances.size()
         << " span-cap allowance(s) matched no cluster of this geometry (recorded "
            "UnusedSpanCapAllowances):";
    for (const auto &label : result.unused_span_cap_allowances)
    {
      text << " " << label;
    }
    text << "\n";
    input.log(text.str());
  }
  result.extension = {extension_passes,
                      extension_portions,
                      extension_sites,
                      extension_length,
                      stack_end_third_body_length,
                      extension_cap_reached,
                      extension_repeat_detected,
                      extension_translational_pieces,
                      extension_translational_length,
                      extension_translational_max_length,
                      extension_translational_two_sided,
                      extension_translational_two_sided_length};
  result.features = features;
  result.radius = R;
  result.reference_process_normal = n_ref;
}

IdentificationResult Identifier::Identify()
{
  IdentificationResult result;
  const auto started = std::chrono::steady_clock::now();
  segment_exclusion.assign(input.segments.size(), std::nullopt);
  stage.Begin("arcs");
  DetectArcs();
  {
    std::size_t corners = 0, bends = 0, absorbed = 0, coarse = 0;
    double worst_sagitta = 0.0;
    for (const auto &arc : arcs)
    {
      (corners += arc.corner ? 1 : 0), (bends += arc.corner ? 0 : 1);
      if (arc.max_sagitta / R >= kArcSagittaOverRadius)
      {
        coarse++;
        worst_sagitta = std::max(worst_sagitta, arc.max_sagitta / R);
      }
    }
    for (const int a : vertex_arc)
    {
      absorbed += a >= 0 ? 1 : 0;
    }
    std::ostringstream text;
    text << corners << " rounded corners, " << bends << " bends of exact radius, "
         << absorbed << " joints absorbed";
    if (coarse > 0)
    {
      text << "; MESH COARSENESS WARNING: " << coarse
           << " arcs with chord sagitta >= " << kArcSagittaOverRadius << " R (worst "
           << std::setprecision(3) << worst_sagitta << " R)";
    }
    stage.End(text.str());
  }
  stage.Begin("runs / chains");
  BuildRuns(input, runs, chains, chain_index);
  BuildRunIndex();
  BuildRunArcGeometry();
  claims.assign(runs.size(), {});
  bent_claims.assign(runs.size(), {});
  stage.End(std::to_string(input.segments.size()) + " segments, " +
            std::to_string(input.vertices.size()) + " vertices, " +
            std::to_string(faces.size()) + " faces, " + std::to_string(runs.size()) +
            " runs, " + std::to_string(chains.size()) + " chains, run grid " +
            std::to_string(run_grid->Cells()) + " cells");
  if (!runs.empty())
  {
    stage.Begin("planes / exclusions");
    ClassifyPlanes();
    {
      std::size_t excluded_runs = 0, cross_layer_runs = 0;
      for (std::size_t r = 0; r < runs.size(); r++)
      {
        excluded_runs += runs[r].excluded ? 1 : 0;
        cross_layer_runs += cross_layer[r].empty() ? 0 : 1;
      }
      stage.End(std::to_string(excluded_runs) + " non-planar runs, " +
                std::to_string(cross_layer_runs) + " runs with cross-layer zones, " +
                std::to_string(cross_layer_vertices.size()) + " cross-layer vertices");
    }
    stage.Begin("parallel classes");
    BuildDirectionClasses();
    stage.End(std::to_string(direction_class_count) + " parallel classes of rigid runs");
    stage.Begin("vertex sites");
    ClassifyVertices();
    stage.End(std::to_string(sites.size()) + " corner / endpoint / junction sites");
    stage.Begin("fillets");
    BuildArcSites();
    stage.End(std::to_string(sites.size()) + " sites incl. rounded corners");
    stage.Begin("curvature");
    ComputeCurvature();
    {
      std::size_t curved = 0;
      for (const auto &chain : chains)
      {
        curved += chain.curved.size();
      }
      stage.End(std::to_string(curved) + " curved chain intervals");
      if (std::getenv("PALACE_IDENTIFICATION_DEBUG_CHAINS") && input.log)
      {
        std::ostringstream dbg;
        dbg << std::setprecision(10);
        for (const auto &chain : chains)
        {
          const Run &first = runs[chain.runs.front()], &last = runs[chain.runs.back()];
          dbg << "  DEBUG chain " << chain.id << (chain.closed ? " closed" : " open")
              << " runs " << chain.runs.size() << " length " << chain.length << " from ("
              << first.start[0] << ", " << first.start[1] << ") to (" << last.end[0] << ", "
              << last.end[1] << ") curved";
          for (const auto &c : chain.curved)
          {
            dbg << " [" << c.first << ", " << c.second << "]";
          }
          dbg << " joints";
          for (std::size_t k = 0; k < chain.joint_turn.size(); k++)
          {
            if (chain.joint_turn[k] > 1.0e-6)
            {
              const Run &run = runs[chain.runs[k]];
              dbg << " (" << run.start[0] << ", " << run.start[1] << ": "
                  << chain.joint_turn[k] * 180.0 / std::acos(-1.0) << " deg)";
            }
          }
          dbg << " arcs";
          for (const auto &[arc, x0, x1] : chain.arc_spans)
          {
            dbg << " (" << arc << " [" << x0 << ", " << x1 << "] r " << arcs[arc].radius
                << ")";
          }
          input.log(dbg.str() + "\n");
          dbg.str("");
        }
      }
    }
    stage.Begin("vertex windows");
    BuildVertexWindows();
    stage.End(std::to_string(sites.size()) + " sites");
    stage.Begin("bent pairs");
    BuildBentPairs();
    stage.Begin("translational spans");
    BuildTranslationalFeatures();
    stage.Begin("events / cores / clusters");
    BuildClusters();
    // Pairs / stacks, then the cluster extension over the single-edge remainder (decision
    // 85(2)), iterated to closure: an enlarged cluster claim moves the stack ends, and the
    // recomposed stacks leave new single-edge portions to test.
    double previous_absorbed = -1.0;
    std::size_t previous_portions = 0;
    for (std::size_t pass = 0;; pass++)
    {
      const std::size_t n_features = features.size();
      stage.Begin("pairs / stacks (pass " + std::to_string(pass + 1) + ")");
      stack_geometric_offsets = 0;
      stack_composition_cap_hits = 0;
      stack_images_merged = 0;
      BuildPairsAndStacks();
      stage.Begin("cluster extension (pass " + std::to_string(pass + 1) + ")");
      const std::size_t portions_before = extension_portions;
      const double absorbed = ExtendClusters();
      extension_passes++;
      const std::size_t pass_portions = extension_portions - portions_before;
      stage.End("absorbed " + std::to_string(absorbed) + " (mesh units) in " +
                std::to_string(extension_portions) + " portions so far, " +
                std::to_string(cluster_claimed.size()) + " clusters");
      // Closure: a pass absorbing at most the signature parameter tolerance (1e-3 R) in
      // total moves no parameter and no gate (DS-SCT-002 at R 2.1: pass 1 absorbed the
      // neighbours, passes 3-15 re-cut 1 nm slivers at the new breakpoints for 8 s each).
      // Its absorptions are applied (they come from the unclaimed remainder: no pair /
      // stack claim overlaps them) and the pairs / stacks are not recomposed again: their
      // cut images differ from the final claims by less than the tolerance (review m8,
      // stated in the design doc; not applying them left sub-tolerance isolated slivers).
      if (absorbed <= kSignatureParameterToleranceOverRadius * R)
      {
        break;
      }
      // Decision 93 (DS-CTX-003 defect 2): a pass that repeats the previous one (the same
      // absorbed length within the signature parameter tolerance and the same portion
      // count: the recomposed stacks re-cut the same slivers at the moved cut images;
      // DS-CTX-003 absorbed 80 nm in 9 portions on every pass from the 3rd, each pass's
      // total differing from the last by picometres) or the pass cap ends the loop the same
      // way (applied, not recomposed; reported under Diagnostics).
      if (previous_absorbed >= 0.0 &&
          std::abs(absorbed - previous_absorbed) <=
              kSignatureParameterToleranceOverRadius * R &&
          pass_portions == previous_portions)
      {
        extension_repeat_detected = true;
        if (input.log)
        {
          input.log("  Identification cluster extension: pass " + std::to_string(pass + 1) +
                    " repeats pass " + std::to_string(pass) + " (absorbed " +
                    std::to_string(absorbed) + " in " + std::to_string(pass_portions) +
                    " portions): closure by repetition\n");
        }
        break;
      }
      if (pass + 1 >= kClusterExtensionMaxPasses)
      {
        extension_cap_reached = true;
        if (input.log)
        {
          input.log("  Identification cluster extension: pass cap " +
                    std::to_string(kClusterExtensionMaxPasses) + " reached (absorbed " +
                    std::to_string(absorbed) + " in the last pass)\n");
        }
        break;
      }
      previous_absorbed = absorbed;
      previous_portions = pass_portions;
      // The pair / stack features and claims are rebuilt on the enlarged claims.
      features.resize(n_features);
      feature_max_kappa.resize(n_features);
      feature_sides.resize(n_features);
      for (auto &run_claims : claims)
      {
        run_claims.erase(std::remove_if(run_claims.begin(), run_claims.end(),
                                        [](const Claim &c) { return c.priority == 2; }),
                         run_claims.end());
      }
    }
    stage.Begin("stack-end third body");
    stack_end_third_body_length = ExtendClusters(true);
    stage.End(std::to_string(stack_end_third_body_length) + " (mesh units)");
    stage.Begin("cluster signatures");
    EmitClusters();
  }
  stage.Begin("claim resolution / assignment / tables");
  Assign(result);
  stage.End(std::to_string(result.features.size()) + " features, " +
            std::to_string(result.exclusions.size()) + " exclusion classes");
  if (!runs.empty() && census_enabled)
  {
    stage.Begin("knife-edge census");
    nlohmann::json census = KnifeEdgeCensus();
    census["ClusterComposition"] = ClusterCompositionBand(result);
    result.knife_edge_census = census.dump();
    {
      // The census is in mesh units (the manifest scales it to its length unit); the log
      // also gives the lengths in units of R so that the two readings cannot be confused.
      const double at_R = census["Distance"]["R"]["Total"].get<double>();
      const double at_2R = census["Distance"]["2R"]["Total"].get<double>();
      std::ostringstream text;
      const auto &composition = census["ClusterComposition"];
      text << "distance within " << kKnifeEdgeBandRelative << " R of R: " << at_R
           << " mesh units = " << std::fixed << std::setprecision(1) << at_R / R
           << " R, of 2R: " << std::defaultfloat << at_2R << " mesh units = " << std::fixed
           << std::setprecision(1) << at_2R / R << " R; cluster composition changing at R "
           << std::defaultfloat << "(1 -/+ " << kKnifeEdgeBandRelative
           << "): " << composition["ClustersBelow"].get<std::size_t>() << " / "
           << composition["ClustersAbove"].get<std::size_t>() << " of "
           << composition["Clusters"].get<std::size_t>() << " clusters, "
           << composition["Below"].get<double>() << " / "
           << composition["Above"].get<double>() << " mesh units";
      stage.End(text.str());
    }
  }
  if (input.log)
  {
    std::ostringstream text;
    text
        << "  Identification total: " << std::fixed << std::setprecision(2)
        << std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count()
        << " s\n";
    input.log(text.str());
  }
  return result;
}

}  // namespace

IdentificationResult IdentifyMetalPerimeter(const IdentificationInput &input)
{
  Identifier identifier(input);
  return identifier.Identify();
}

namespace
{

// Byte writer / reader of the result: PODs raw, strings and arrays length-prefixed, read
// back in the order written.
class ByteWriter
{
public:
  template <typename T>
  void Pod(const T &value)
  {
    static_assert(std::is_trivially_copyable_v<T>);
    const auto *bytes = reinterpret_cast<const char *>(&value);
    buffer.append(bytes, sizeof(T));
  }
  void Size(std::size_t n) { Pod(static_cast<std::uint64_t>(n)); }
  void String(const std::string &text)
  {
    Size(text.size());
    buffer.append(text);
  }
  void Point(const std::array<double, 3> &p)
  {
    for (const double v : p)
    {
      Pod(v);
    }
  }
  std::string buffer;
};

class ByteReader
{
public:
  explicit ByteReader(const std::string &buffer_) : buffer(buffer_) {}
  template <typename T>
  T Pod()
  {
    static_assert(std::is_trivially_copyable_v<T>);
    MFEM_VERIFY(pos + sizeof(T) <= buffer.size(),
                "Identification result buffer ends early!");
    T value;
    std::memcpy(&value, buffer.data() + pos, sizeof(T));
    pos += sizeof(T);
    return value;
  }
  std::size_t Size() { return static_cast<std::size_t>(Pod<std::uint64_t>()); }
  std::string String()
  {
    const std::size_t n = Size();
    MFEM_VERIFY(pos + n <= buffer.size(), "Identification result buffer ends early!");
    std::string text(buffer.data() + pos, n);
    pos += n;
    return text;
  }
  std::array<double, 3> Point()
  {
    std::array<double, 3> p;
    for (double &v : p)
    {
      v = Pod<double>();
    }
    return p;
  }
  bool Done() const { return pos == buffer.size(); }

private:
  const std::string &buffer;
  std::size_t pos = 0;
};

}  // namespace

std::string SerializeIdentificationResult(const IdentificationResult &result)
{
  ByteWriter w;
  w.Pod(result.radius);
  w.Point(result.reference_process_normal);
  w.Pod(result.perimeter_length);
  w.Size(result.real_segments);
  w.Pod(result.assigned_length);
  w.Pod(result.excluded_length);
  w.String(result.geometry_digest);
  w.Size(result.same_priority_claim_overlaps);
  w.Size(result.stack_geometric_offsets);
  w.Size(result.stack_composition_cap_hits);
  w.Size(result.extension.passes);
  w.Size(result.extension.portions);
  w.Size(result.extension.sites);
  w.Pod(result.extension.length);
  w.Pod(result.extension.stack_end_third_body_length);
  w.Pod(result.extension.cap_reached);
  w.Pod(result.extension.repeat_detected);
  w.Size(result.extension.translational_pieces);
  w.Pod(result.extension.translational_length);
  w.Pod(result.extension.translational_max_length);
  w.Size(result.extension.translational_two_sided);
  w.Pod(result.extension.translational_two_sided_length);
  w.Size(result.stack_images_merged);
  w.Size(result.sub_tolerance_portions.count);
  w.Pod(result.sub_tolerance_portions.length);
  w.Pod(result.sub_tolerance_portions.max_length);
  w.Size(result.sub_tolerance_portions.isolated);
  w.String(result.knife_edge_census);
  w.Size(result.unused_span_cap_allowances.size());
  for (const auto &label : result.unused_span_cap_allowances)
  {
    w.String(label);
  }
  w.Size(result.features.size());
  for (const auto &f : result.features)
  {
    w.Pod(f.id);
    w.String(f.type);
    w.String(f.signature.dump());
    w.String(f.signature_key);
    w.String(f.hash);
    w.Pod(f.chirality);
    w.Pod(f.length);
    w.Size(f.portions.size());
    for (const auto &p : f.portions)
    {
      w.Size(p.segment);
      w.Pod(p.s0);
      w.Pod(p.s1);
      w.Pod(p.side);
      w.Pod(p.turn);
      w.Pod(p.stretch);
    }
    w.Size(f.vertices.size());
    for (const std::size_t v : f.vertices)
    {
      w.Size(v);
    }
    w.Point(f.origin);
    for (const auto &axis : f.axes)
    {
      w.Point(axis);
    }
    w.Pod(f.bend_radius_over_R.has_value());
    w.Pod(f.bend_radius_over_R.value_or(0.0));
    w.Pod(f.exact_parameters);
    w.Pod(f.matched_model.has_value());
    w.String(f.matched_model.value_or(""));
    w.Pod(f.match_deviation.value_or(-1.0));
    w.Pod(f.match_note.has_value());
    w.String(f.match_note.value_or(""));
    w.String(f.quantum_near_match.is_null() ? std::string() : f.quantum_near_match.dump());
    w.String(f.spatial_support.is_null() ? std::string() : f.spatial_support.dump());
    w.String(f.legacy_contract.is_null() ? std::string() : f.legacy_contract.dump());
    w.String(f.blend_eigenvalues.is_null() ? std::string() : f.blend_eigenvalues.dump());
    w.String(f.interpolation_fallback.is_null() ? std::string()
                                                : f.interpolation_fallback.dump());
    w.String(f.mirror.is_null() ? std::string() : f.mirror.dump());
    w.Point(f.claims_origin);
    for (const auto &axis : f.claims_axes)
    {
      w.Point(axis);
    }
    w.Pod(f.claims_chirality);
  }
  w.Size(result.spatial_support.clusters);
  w.Size(result.spatial_support.claims_keyed);
  w.Size(result.spatial_support.context_keyed);
  w.Size(result.spatial_support.grown);
  w.Size(result.spatial_support.unboxable);
  w.Size(result.spatial_support.exceeding_span_cap);
  w.Pod(result.spatial_support.chain_length);
  w.Pod(result.spatial_support.foreign_length);
  w.Pod(result.spatial_support.fictitious_continuation_length);
  w.Size(result.spatial_support.threshold_band_hits);
  w.Size(result.segments.size());
  for (const auto &s : result.segments)
  {
    w.Point(s.key[0]);
    w.Point(s.key[1]);
    w.Pod(s.length);
    w.Pod(s.chain);
    w.Pod(s.arc);
    w.Size(s.portions.size());
    for (const auto &p : s.portions)
    {
      w.Point(p);
    }
    w.Size(s.excluded_portions.size());
    for (const auto &p : s.excluded_portions)
    {
      w.Point(p);
    }
    w.Pod(s.exclusion.has_value());
    w.String(s.exclusion ? s.exclusion->first : "");
    w.String(s.exclusion ? s.exclusion->second : "");
  }
  w.Size(result.vertices.size());
  for (const auto &v : result.vertices)
  {
    w.Size(v.vertex);
    w.String(v.type);
    w.Pod(v.turn_degrees);
    w.Pod(v.feature);
    w.Pod(v.point_contact);
  }
  w.Size(result.exclusions.size());
  for (const auto &x : result.exclusions)
  {
    w.String(x.cls);
    w.String(x.reason);
    w.Pod(x.count);
    w.Pod(x.length);
  }
  w.Size(result.arcs.size());
  for (const auto &a : result.arcs)
  {
    w.Point(a.center);
    w.Pod(a.radius);
    w.Pod(a.turn_degrees);
    w.Pod(a.max_sagitta_over_R);
    w.String(a.kind);
    w.Size(a.joints);
    w.Size(a.segments);
  }
  return std::move(w.buffer);
}

IdentificationResult DeserializeIdentificationResult(const std::string &buffer)
{
  ByteReader r(buffer);
  IdentificationResult result;
  result.radius = r.Pod<double>();
  result.reference_process_normal = r.Point();
  result.perimeter_length = r.Pod<double>();
  result.real_segments = r.Size();
  result.assigned_length = r.Pod<double>();
  result.excluded_length = r.Pod<double>();
  result.geometry_digest = r.String();
  result.same_priority_claim_overlaps = r.Size();
  result.stack_geometric_offsets = r.Size();
  result.stack_composition_cap_hits = r.Size();
  result.extension.passes = r.Size();
  result.extension.portions = r.Size();
  result.extension.sites = r.Size();
  result.extension.length = r.Pod<double>();
  result.extension.stack_end_third_body_length = r.Pod<double>();
  result.extension.cap_reached = r.Pod<bool>();
  result.extension.repeat_detected = r.Pod<bool>();
  result.extension.translational_pieces = r.Size();
  result.extension.translational_length = r.Pod<double>();
  result.extension.translational_max_length = r.Pod<double>();
  result.extension.translational_two_sided = r.Size();
  result.extension.translational_two_sided_length = r.Pod<double>();
  result.stack_images_merged = r.Size();
  result.sub_tolerance_portions.count = r.Size();
  result.sub_tolerance_portions.length = r.Pod<double>();
  result.sub_tolerance_portions.max_length = r.Pod<double>();
  result.sub_tolerance_portions.isolated = r.Size();
  result.knife_edge_census = r.String();
  result.unused_span_cap_allowances.resize(r.Size());
  for (auto &label : result.unused_span_cap_allowances)
  {
    label = r.String();
  }
  result.features.resize(r.Size());
  for (auto &f : result.features)
  {
    f.id = r.Pod<int>();
    f.type = r.String();
    f.signature = nlohmann::json::parse(r.String());
    f.signature_key = r.String();
    f.hash = r.String();
    f.chirality = r.Pod<int>();
    f.length = r.Pod<double>();
    f.portions.resize(r.Size());
    for (auto &p : f.portions)
    {
      p.segment = r.Size();
      p.s0 = r.Pod<double>();
      p.s1 = r.Pod<double>();
      p.side = r.Pod<int>();
      p.turn = r.Pod<double>();
      p.stretch = r.Pod<int>();
    }
    f.vertices.resize(r.Size());
    for (auto &v : f.vertices)
    {
      v = r.Size();
    }
    f.origin = r.Point();
    for (auto &axis : f.axes)
    {
      axis = r.Point();
    }
    const bool has_bend = r.Pod<bool>();
    const double bend = r.Pod<double>();
    if (has_bend)
    {
      f.bend_radius_over_R = bend;
    }
    f.exact_parameters = r.Pod<bool>();
    const bool has_model = r.Pod<bool>();
    const std::string model = r.String();
    const double deviation = r.Pod<double>();
    if (deviation >= 0.0)
    {
      f.match_deviation = deviation;
    }
    if (has_model)
    {
      f.matched_model = model;
    }
    const bool has_note = r.Pod<bool>();
    const std::string note = r.String();
    if (has_note)
    {
      f.match_note = note;
    }
    const std::string near_match = r.String();
    f.quantum_near_match =
        near_match.empty() ? nlohmann::json(nullptr) : nlohmann::json::parse(near_match);
    const std::string support = r.String();
    f.spatial_support =
        support.empty() ? nlohmann::json(nullptr) : nlohmann::json::parse(support);
    const std::string legacy = r.String();
    f.legacy_contract =
        legacy.empty() ? nlohmann::json(nullptr) : nlohmann::json::parse(legacy);
    const std::string eigenvalues = r.String();
    f.blend_eigenvalues =
        eigenvalues.empty() ? nlohmann::json(nullptr) : nlohmann::json::parse(eigenvalues);
    const std::string fallback = r.String();
    f.interpolation_fallback =
        fallback.empty() ? nlohmann::json(nullptr) : nlohmann::json::parse(fallback);
    const std::string mirror = r.String();
    f.mirror = mirror.empty() ? nlohmann::json(nullptr) : nlohmann::json::parse(mirror);
    f.claims_origin = r.Point();
    for (auto &axis : f.claims_axes)
    {
      axis = r.Point();
    }
    f.claims_chirality = r.Pod<int>();
  }
  result.spatial_support.clusters = r.Size();
  result.spatial_support.claims_keyed = r.Size();
  result.spatial_support.context_keyed = r.Size();
  result.spatial_support.grown = r.Size();
  result.spatial_support.unboxable = r.Size();
  result.spatial_support.exceeding_span_cap = r.Size();
  result.spatial_support.chain_length = r.Pod<double>();
  result.spatial_support.foreign_length = r.Pod<double>();
  result.spatial_support.fictitious_continuation_length = r.Pod<double>();
  result.spatial_support.threshold_band_hits = r.Size();
  result.segments.resize(r.Size());
  for (auto &s : result.segments)
  {
    s.key[0] = r.Point();
    s.key[1] = r.Point();
    s.length = r.Pod<double>();
    s.chain = r.Pod<int>();
    s.arc = r.Pod<int>();
    s.portions.resize(r.Size());
    for (auto &p : s.portions)
    {
      p = r.Point();
    }
    s.excluded_portions.resize(r.Size());
    for (auto &p : s.excluded_portions)
    {
      p = r.Point();
    }
    const bool has_exclusion = r.Pod<bool>();
    const std::string cls = r.String(), reason = r.String();
    if (has_exclusion)
    {
      s.exclusion = std::make_pair(cls, reason);
    }
  }
  result.vertices.resize(r.Size());
  for (auto &v : result.vertices)
  {
    v.vertex = r.Size();
    v.type = r.String();
    v.turn_degrees = r.Pod<double>();
    v.feature = r.Pod<int>();
    v.point_contact = r.Pod<bool>();
  }
  result.exclusions.resize(r.Size());
  for (auto &x : result.exclusions)
  {
    x.cls = r.String();
    x.reason = r.String();
    x.count = r.Pod<int>();
    x.length = r.Pod<double>();
  }
  result.arcs.resize(r.Size());
  for (auto &a : result.arcs)
  {
    a.center = r.Point();
    a.radius = r.Pod<double>();
    a.turn_degrees = r.Pod<double>();
    a.max_sagitta_over_R = r.Pod<double>();
    a.kind = r.String();
    a.joints = r.Size();
    a.segments = r.Size();
  }
  MFEM_VERIFY(r.Done(), "Identification result buffer was not consumed exactly!");
  return result;
}

// The census lengths (mesh units) converted to the manifest units; counts unchanged.
nlohmann::json ScaledCensus(const std::string &text, double length_scale)
{
  if (text.empty())
  {
    return nullptr;
  }
  nlohmann::json census = nlohmann::json::parse(text);
  for (const char *section : {"Distance", "BendRadius"})
  {
    if (!census.contains(section))
    {
      continue;
    }
    for (auto &[threshold, entry] : census[section].items())
    {
      (void)threshold;
      for (const char *key : {"Below", "Above", "Total"})
      {
        if (entry.contains(key))
        {
          entry[key] = entry[key].get<double>() * length_scale;
        }
      }
    }
  }
  if (census.contains("SampledLength"))
  {
    census["SampledLength"] = census["SampledLength"].get<double>() * length_scale;
  }
  if (census.contains("ClusterComposition"))
  {
    auto &entry = census["ClusterComposition"];
    for (const char *key : {"Below", "Above", "Total"})
    {
      if (entry.contains(key))
      {
        entry[key] = entry[key].get<double>() * length_scale;
      }
    }
  }
  return census;
}

nlohmann::json IdentificationResult::ToJson(double length_scale) const
{
  const double scaled_radius = radius * length_scale;
  const double tolerance =
      std::max(1.0e-10 * scaled_radius, 64.0 * std::numeric_limits<double>::epsilon());
  const double step = std::pow(10.0, std::floor(std::log10(tolerance)));
  auto L = [&](double value)
  {
    const double snapped = std::round(value * length_scale / step) * step;
    return snapped == 0.0 ? 0.0 : snapped;
  };
  auto P = [&](const std::array<double, 3> &p)
  { return nlohmann::json{L(p[0]), L(p[1]), L(p[2])}; };
  auto D = [&](const std::array<double, 3> &v)
  {
    auto S = [](double x)
    {
      const double s = std::round(x * 1.0e12) * 1.0e-12;
      return s == 0.0 ? 0.0 : s;
    };
    return nlohmann::json{S(v[0]), S(v[1]), S(v[2])};
  };
  // Mirror band: image segments [real_segments, size) are listed apart (ImagePortions /
  // ImageSegments), never as perimeter.
  const std::size_t real_segment_count =
      real_segments > 0 && real_segments < segments.size() ? real_segments
                                                           : segments.size();
  nlohmann::json feature_list = nlohmann::json::array();
  for (const auto &feature : features)
  {
    nlohmann::json portions = nlohmann::json::array();
    nlohmann::json image_portions = nlohmann::json::array();
    nlohmann::json sides = nlohmann::json::array();
    nlohmann::json turns = nlohmann::json::array();
    bool multi_sided = false;
    double total_turn = 0.0;
    for (const auto &portion : feature.portions)
    {
      if (portion.segment >= real_segment_count)
      {
        image_portions.push_back(
            {portion.segment, L(portion.s0), L(portion.s1), portion.side});
        continue;
      }
      portions.push_back({portion.segment, L(portion.s0), L(portion.s1)});
      sides.push_back(portion.side);
      multi_sided = multi_sided || portion.side != 0;
      const double turn = std::round(portion.turn * 1.0e9) * 1.0e-9;
      turns.push_back(turn == 0.0 ? 0.0 : turn);
      total_turn += portion.turn;
    }
    nlohmann::json entry = {
        {"Id", feature.id},
        {"Type", feature.type},
        {"Signature", feature.signature},
        {"Hash", feature.hash},
        {"Chirality", feature.chirality},
        {"Length", L(feature.length)},
        {"Portions", portions},
        {"Vertices", feature.vertices},
        {"Frame", FrameToJson(feature.origin, feature.axes, radius, length_scale)},
        {"BendRadiusOverR",
         feature.bend_radius_over_R
             ? nlohmann::json(std::round(*feature.bend_radius_over_R /
                                         kSignatureLengthQuantumOverRadius) *
                              kSignatureLengthQuantumOverRadius)
             : nlohmann::json(nullptr)},
        {"ExactParameters", feature.exact_parameters},
        {"Match", {{"Status", feature.matched_model ? "Matched" : "Missing"}}}};
    if (feature.matched_model)
    {
      entry["Match"]["Model"] = *feature.matched_model;
      entry["Match"]["Deviation"] = feature.match_deviation.value_or(0.0);
    }
    if (feature.match_note)
    {
      entry["Match"]["Note"] = *feature.match_note;
    }
    if (!feature.quantum_near_match.is_null())
    {
      // Matched within the quantum near-match (block (b) DESIGN section 4): the model's key
      // and the feature's own key differ by <= kClusterQuantumNearMatchMaxQuanta quanta.
      entry["Match"]["QuantumNearMatch"] = feature.quantum_near_match;
    }
    if (!feature.legacy_contract.is_null())
    {
      // Matched through an explicit legacy-contract alias of the library (USER decision
      // 283): the legacy model, the aliased v3 key, the verified context digest and the
      // recorded context.
      entry["Match"]["LegacyContract"] = feature.legacy_contract;
    }
    if (!feature.mirror.is_null())
    {
      // The feature touches a mirror image across a natural truncation plane (boundary-cut
      // DESIGN 2.2.2; decisions 442 / 454): {Planes, RealLength, ImageLength, Status}; its
      // portions on image segments (never placed, never counted) listed apart. An unmerged
      // mirror-formed configuration's Contract (decision 557) carries its Frame in the
      // manifest's units like Features[].Frame (the record keeps the identification's).
      entry["Mirror"] = feature.mirror;
      if (feature.mirror.contains("Contract") && feature.mirror["Contract"].is_object())
      {
        entry["Mirror"]["Contract"]["Frame"] =
            FrameToJson(feature.origin, feature.axes, radius, length_scale);
        entry["Mirror"]["Contract"]["Frame"]["Chirality"] = feature.chirality;
      }
      if (!image_portions.empty())
      {
        entry["ImagePortions"] = image_portions;
      }
    }
    if (!feature.blend_eigenvalues.is_null())
    {
      // An angle-interpolated corner (decision 376): the applied blend's min eigenvalues
      // on the free knots relative to max |eigenvalue| per blended matrix.
      entry["Match"]["BlendEigenvalues"] = feature.blend_eigenvalues;
    }
    if (!feature.interpolation_fallback.is_null())
    {
      // The stencil's blended fabricated or thin domain matrix was not PSD beyond roundoff:
      // the convex linear blend of the two bracketing nodes is applied (decision 374 (B)).
      entry["Match"]["InterpolationFallback"] = feature.interpolation_fallback;
    }
    if (!feature.spatial_support.is_null())
    {
      // Spatial-support record (contract v3, decision 282; units of R in the feature's
      // frame): the claims-derived and the grown box, the face-rule outcome, the context
      // census (chain / foreign pieces and vertices), the legacy contract's fictitious
      // straight continuation, the window truncation inside the box.
      entry["SpatialSupport"] = feature.spatial_support;
      // The claims-only canonical frame (the placement frame of a legacy-contract alias,
      // USER decision 283); equal to Frame unless the key is contract 3.
      entry["ClaimsFrame"] = {{"Origin", P(feature.claims_origin)},
                              {"Axes",
                               {D(feature.claims_axes[0]), D(feature.claims_axes[1]),
                                D(feature.claims_axes[2])}},
                              {"Chirality", feature.claims_chirality}};
    }
    if (multi_sided)
    {
      // Side of every portion of a pair / parallel cluster (parallel to Portions): the
      // signature's edge order for chirality +1, reversed for -1.
      entry["Sides"] = sides;
    }
    if (feature.bend_radius_over_R)
    {
      // Signed turn toward the metal of every portion (radians, parallel to Portions; not
      // hashed): the weight of the first-order curvature term on a straight-like feature
      // (design (b)7); absent on features of straight chains. Their sum for a one-sided
      // feature (the sides of a pair turn in opposite senses).
      entry["PortionTurns"] = turns;
      if (!multi_sided)
      {
        const double sum = std::round(total_turn * 1.0e9) * 1.0e-9;
        entry["TurnTowardMetal"] = sum == 0.0 ? 0.0 : sum;
      }
    }
    feature_list.push_back(std::move(entry));
  }
  nlohmann::json segment_list = nlohmann::json::array();
  nlohmann::json image_segment_list = nlohmann::json::array();
  for (std::size_t s = 0; s < segments.size(); s++)
  {
    const auto &segment = segments[s];
    if (s >= real_segment_count)
    {
      image_segment_list.push_back({{"Key", {P(segment.key[0]), P(segment.key[1])}},
                                    {"Length", L(segment.length)},
                                    {"Chain", segment.chain}});
      continue;
    }
    nlohmann::json entry = {{"Key", {P(segment.key[0]), P(segment.key[1])}},
                            {"Length", L(segment.length)},
                            {"Chain", segment.chain}};
    if (segment.arc >= 0)
    {
      entry["Arc"] = segment.arc;
    }
    if (segment.exclusion)
    {
      entry["Exclusion"] = {{"Class", segment.exclusion->first},
                            {"Reason", segment.exclusion->second}};
    }
    else
    {
      nlohmann::json portions = nlohmann::json::array();
      for (const auto &portion : segment.portions)
      {
        portions.push_back({L(portion[0]), L(portion[1]), static_cast<int>(portion[2])});
      }
      entry["Portions"] = portions;
      if (!segment.excluded_portions.empty())
      {
        nlohmann::json excluded = nlohmann::json::array();
        for (const auto &portion : segment.excluded_portions)
        {
          excluded.push_back({L(portion[0]), L(portion[1]), static_cast<int>(portion[2])});
        }
        entry["ExcludedPortions"] = excluded;
      }
    }
    segment_list.push_back(std::move(entry));
  }
  nlohmann::json arc_list = nlohmann::json::array();
  for (const auto &arc : arcs)
  {
    arc_list.push_back(
        {{"Center", P(arc.center)},
         {"Radius", L(arc.radius)},
         {"RadiusOverR", RoundTo(arc.radius / radius, kSignatureLengthQuantumOverRadius)},
         {"TurnDegrees", RoundTo(arc.turn_degrees, kSignatureAngleQuantumDegrees)},
         {"MaxChordSagittaOverR",
          RoundTo(arc.max_sagitta_over_R, kSignatureLengthQuantumOverRadius)},
         {"Kind", arc.kind},
         {"Joints", arc.joints},
         {"Segments", arc.segments}});
  }
  // Mesh-coarseness WARNING (USER decision 122): the arcs whose largest chord sagitta
  // reaches SagittaOverR x R are arcs (the concyclicity rule decides membership) whose
  // polyline the correction cannot resolve; count, arc length, worst sagitta and the arc
  // indices are listed so that a coarse mesh is never silent.
  nlohmann::json coarse_arcs = nlohmann::json::array();
  double coarse_length = 0.0, worst_sagitta = 0.0;
  for (std::size_t a = 0; a < arcs.size(); a++)
  {
    if (arcs[a].max_sagitta_over_R >= kArcSagittaOverRadius)
    {
      coarse_arcs.push_back(a);
      coarse_length += arcs[a].radius * arcs[a].turn_degrees * std::acos(-1.0) / 180.0;
      worst_sagitta = std::max(worst_sagitta, arcs[a].max_sagitta_over_R);
    }
  }
  const nlohmann::json mesh_coarseness_warning = {
      {"SagittaOverR", kArcSagittaOverRadius},
      {"Rule", "fitted arcs (table Arcs) whose largest chord sagitta rho (1 - cos(central "
               "angle / 2)) is at least SagittaOverR x R: the polyline chords these arcs "
               "coarser than the resolution of the correction (a diagnostic of the mesh; "
               "membership is by concyclicity and the joint-turn cap, ArcRule)"},
      {"Count", coarse_arcs.size()},
      {"Length", L(coarse_length)},
      {"WorstMaxChordSagittaOverR",
       RoundTo(worst_sagitta, kSignatureLengthQuantumOverRadius)},
      {"Arcs", coarse_arcs}};
  nlohmann::json vertex_list = nlohmann::json::array();
  for (const auto &vertex : vertices)
  {
    nlohmann::json entry = {{"Type", vertex.type}};
    if (vertex.vertex != std::numeric_limits<std::size_t>::max())
    {
      entry["Vertex"] = vertex.vertex;
    }
    if (vertex.type != "TruncationCut" && vertex.type != "PortCut" &&
        vertex.type != "Excluded" && vertex.type != "ExclusionCut")
    {
      entry["TurnDegrees"] = std::round(vertex.turn_degrees * 1.0e6) * 1.0e-6;
      if (vertex.type != "BendVertex")
      {
        entry["Feature"] = vertex.feature;
      }
    }
    if (vertex.point_contact)
    {
      entry["PointContact"] = true;
    }
    vertex_list.push_back(std::move(entry));
  }
  nlohmann::json exclusion_list = nlohmann::json::array();
  for (const auto &exclusion : exclusions)
  {
    exclusion_list.push_back({{"Class", exclusion.cls},
                              {"Reason", exclusion.reason},
                              {"Count", exclusion.count},
                              {"Length", L(exclusion.length)}});
  }
  nlohmann::json manifest = {
      {"Version", 2},
      {"MatchingRadius", scaled_radius},
      {"Conventions",
       {{"JointNoiseSagittaOverR", kJointNoiseSagittaOverRadius},
        {"CornerRule",
         "geometric joint NOISE threshold (USER decision 121 (B); was the 1 deg angular "
         "threshold of decision 117(4) and the 30 deg corner class of decision 73): a "
         "vertex "
         "with two path segments turning by t between two straight pieces (collinear mesh "
         "segments merged) is a straight continuation of its chain (its turn feeds the "
         "windowed curvature) when the implied sagitta (c / 2) tan(t / 4) of the SHORTER "
         "adjacent piece c is below JointNoiseSagittaOverR x R — the deviation from "
         "straight "
         "the joint implies at the resolution of the correction, the same quantity the arc "
         "rule records per chord (sub-nm mesh slivers are noise whatever their turn, "
         "spline "
         "steps of 1-6 deg on 0.5-2.6 R chords imply 0.001-0.034 R, 5 um arms turning 20 "
         "deg imply 0.23 R); every other joint is a corner feature unless a fitted arc "
         "absorbs it (ArcRule). The chains of the perimeter extraction (metaledge.cpp) "
         "break "
         "at the same rule (metaledge.hpp kJointNoiseSagittaOverRadius, JointIsNoise, "
         "compared on a 1e-9 relative grid); a bend arc merges the chains its absorbed "
         "corners separated"},
        {"InteractionDistanceOverR", kInteractionDistanceOverRadius},
        {"ThroughVertexZoneOverR", kThroughVertexZoneOverRadius},
        {"ClusterBallOverR", kClusterBallOverRadius},
        {"VertexJoinsClusterOverR", kVertexJoinsClusterOverRadius},
        {"VertexWindowOverR", kVertexWindowOverRadius},
        {"ParallelCosineTolerance", kParallelCosineTolerance},
        {"ParallelClassRule",
         "one parallel relation for every rule (USER decision 184 (1), 2026-10-01): two "
         "straight (single-run) runs of one metal plane are parallel iff |cos| of their "
         "tangents exceeds 1 - ParallelCosineTolerance (1.4e-4 rad; dimensionless), the "
         "parallel classes are the connected components of that relation per plane (built "
         "from the runs' 1e-9 direction buckets sorted by their in-plane angle, "
         "consecutive "
         "buckets within the tolerance linked, the 0 / pi wrap included: a canonical, "
         "numbering-independent partition). Parallel rigid runs pair and stack by the "
         "translational rule ONLY (their class), every other rigid pair by the bent-pair / "
         "event rules; before this rule the translational stage grouped runs on the 1e-9 "
         "direction grid while the other two stages skipped rigid pairs within the cosine "
         "tolerance, so runs tilted by between 1e-9 and 1.4e-4 rad (a mesh-noise wobble of "
         "1.8e-4 um over a 186 um chip edge) were paired by neither. The linkage is "
         "transitive (single linkage over the sorted buckets): a dense continuum of "
         "directions would chain a class beyond the tolerance (not a chip configuration: "
         "rigid runs are isolated design edges). A class pair's separation (its lateral "
         "offset difference at the runs' midpoints) is exact only when the lines' "
         "separation at both ends of the span stays within "
         "SignatureParameterToleranceOverR "
         "R of it (tilt x half-length), else the stack reads it with ExactParameters "
         "false"},
        {"ArcFitToleranceOverR", kArcFitToleranceOverRadius},
        {"ArcMaxJointTurnDegrees", kArcMaxJointTurnDegrees},
        {"CornerTraceBasisRule",
         "corner coupons of an angle-interpolated family are built on the trace basis rule "
         "(corner-family review 2026-09-29): on every box ring that meets the metal the "
         "knots are the two metal-arm crossings (PEC), MetalInteriorKnots = 1 knot at "
         "equal "
         "perimeter-arc-length fractions of the metal arc (PEC) and FreeKnots = 5 knots at "
         "equal fractions of the free arc, ordered by role (one knot semantics, one zero "
         "set, like-to-like free knots for every node); box corners that are no knot are "
         "slave trace vertices; the library load refuses a spatial corner coupon whose "
         "metal arm crosses a box ring at no PEC knot; the runtime basis of an "
         "interpolated "
         "corner is constructed by the rule at the feature's angle (TraceBasis record: "
         "RingSize 8, MetalInteriorKnots 1, FreeKnots 5, Fractions PerimeterArcLength; "
         "knot coincidence 1e-6 of the perimeter: a free or metal-interior knot within it "
         "of a fixed fraction k / RingSize takes that fraction, a box corner within it of "
         "any knot gets no slave vertex). Kink-aware interpolation (corner-qualification "
         "block 2026-09-29): a stencil never straddles a geometric EVENT of the trace "
         "basis "
         "(a knot of a metal ring passing a fixed-layout vertex: the perimeter-ordered "
         "band "
         "triangulation flips a quad diagonal and the hats jump; convex knot-corner "
         "passages 90 / 135 / 153.434948822922 deg, concave 90 / 135 / 158.198590513648 "
         "deg; side-midpoint passages convex 141.340191745910, concave 111.801409486352 "
         "deg); coupons sharing TraceBasis ConnectivityAngleDegrees (the band merge keyed "
         "by the rule's layout at that angle: one triangulation per SEGMENT, which removes "
         "the midpoint flips) form a segment with one coupon per side at every knot-corner "
         "passage; the stencil is Lagrange on the segment's nodes nearest to the angle "
         "(cubic on four, else the highest order the segment supports), the runtime basis "
         "constructed with the segment's connectivity; an exact node (|node - angle| <= "
         "the signature angle tolerance, > it interior: the boundary is never in no "
         "segment) prefers the legacy coupon at that angle (production libraries carry the "
         "legacy tie coupons at 90 / 135 / 180, decision 140 (1)), else the lower-angle "
         "segment's; coupons without the record (legacy) are exact matches only; a segment "
         "across a passage, two coupons at one angle in a segment or overlapping segments "
         "fail closed at library load (CheckCornerFamilySegments), a node's trace "
         "triangulation not the rule's for its connectivity at match time. REFINED rule "
         "(corner-basis refinement 2026-09-30, USER decision 161, TraceBasis RingLayout "
         "AllRingsFollowMetal): EVERY ring of the box (outer rings at -R, -R/3, "
         "-OveretchDepth, 0, MetalThickness, MetalThickness + k OveretchDepth for k in "
         "ExtraLevelsAboveOverOveretch = [1, 4], R/3, R and the two inner cap rings) "
         "carries "
         "the same fractions — the two crossings, MetalInteriorKnots = 5 at equal "
         "fractions "
         "of the metal arc, FreeKnots = 9 with FreeKnotGrading [1/3, 2/3] (knots at R/3 "
         "and "
         "2R/3 along the perimeter from each crossing on the free side, 5 at equal "
         "fractions "
         "between; on a free arc shorter than FreeKnotGradingReferenceFreeArcOverR = "
         "2 - cot 60 deg, the concave 60-degree node's, the graded distances scale by "
         "FreeArc / Reference so acute concave nodes keep that node's proportions, slots "
         "and zero set, block (b) family 4, decision 318; the key is recorded only where "
         "the scaling is active, absent = the default) — PEC on the two metal rings only "
         "(14 of 176 knots), box corners slaves "
         "on "
         "every ring, cap centres slaves at the mean of the two crossing knots, every band "
         "a "
         "regular column grid: NO events (CornerBasisEvents empty; a knot passing a corner "
         "slave on every ring at once leaves the hats continuous, measured 1e-5), one "
         "segment per convexity (no connectivity angle; a coupon carrying one is refused), "
         "the stencil the cubic sliding window on the four nearest nodes; the library load "
         "checks every coupon's ring levels, basis points and trace mesh against the rule "
         "at "
         "its own angle (CheckCornerRuleCouponFiles). The coupon self-check judges the "
         "basis "
         "against a held-out reference DECOUPLED from the basis levels (64 knots per ring, "
         "rings every OveretchDepth across both R/3 ramps + MetalThickness + R/3: the "
         "vertical representation is gated too); the option-(c) traces gate the family's "
         "interpolation check (USER decision 161 (2))"},
        {"SagittaOverR", kArcSagittaOverRadius},
        {"ArcSampleSpacingOverR", kArcSampleSpacingOverRadius},
        {"ClusterArcChordStepDegrees", kClusterArcChordStepDegrees},
        {"ClusterArcChordMaxLengthOverR", kClusterArcChordMaxLengthOverRadius},
        {"ClusterGeometryRule",
         "option A (decision 91(1)): every run on a fitted arc (rounded corner or bend, "
         "table Arcs; Segments[].Arc) is evaluated on the arc by the cluster machinery — "
         "event cores, core merging, the vertex join, the radius-R claims, the through "
         "zones, the extension's across rule (radial projection inside an arc piece; the "
         "wedge turn between arc tangents), the stack-end images of cluster / window ends "
         "on concentric partner arcs — so that a cluster's extent is a function of the "
         "design curves; the sublevel sets along an arc are bracketed by samples at most "
         "ArcSampleSpacingOverR x R apart and their crossings bisected; a cluster's claims "
         "on "
         "one arc are ONE signature portion serialised as {P: sorted ends, Arc: centre + "
         "midpoint, GapRadial: +1 metal inside the circle / -1 outside} in the frame "
         "(frame "
         "candidates: the arc's end tangents and end radial directions; origin: the "
         "length-weighted centroid with the arcs' analytic centroids); the library builder "
         "chords a signature arc at ClusterArcChordStepDegrees, finer so that no chord "
         "exceeds ClusterArcChordMaxLengthOverR x R (canonical: the coupon geometry does "
         "not "
         "depend on the device mesh); two chords keep the former exact formulas"},
        {"ArcRule",
         "concyclicity rule (USER decisions 121 / 122, 2026-09-28; amends the sagitta form "
         "of 117(4)): a run of >= 3 consecutive joints of a perimeter path (through corner "
         "vertices), same turn sign, <= 180 deg in total, is ONE arc iff its joint "
         "vertices "
         "lie on one circle within ArcFitToleranceOverR x R AND every joint turns less "
         "than "
         "ArcMaxJointTurnDegrees (regular polygons such as squares / hexagons stay "
         "corners); "
         "bends and rounded corners alike, whatever the chord sagitta and with no "
         "piece-length rule, PROVIDED the concyclicity holds at every point of the range "
         "(USER decision 184 (2), 2026-10-01): the interior points of a chord between two "
         "consecutive joints deviate from the circle by the chord sagitta, admissible only "
         "where it is below the joint noise resolution (JointNoiseSagittaOverR x R) or "
         "where one of the chord's end joints is a real turn (not noise under CornerRule "
         "on its shorter piece: a design chord of a coarsely discretised curve); a chord "
         "beyond the resolution between two noise joints is a straight edge and its range "
         "is no arc (a long straight edge with tiny same-sign end joints, exactly "
         "concyclic by mirror symmetry, is never a bend). The chord sagitta rho (1 - "
         "cos(central angle / 2)) is recorded "
         "per arc (Arcs[].MaxChordSagittaOverR) and arcs at or above SagittaOverR x R are "
         "listed in MeshCoarsenessWarning (count, length, worst): a mesh-coarseness "
         "diagnostic, not a membership test. The circle is the one tangent to both arms at "
         "the end joints when every joint lies on it (tangent-length radius); otherwise, "
         "for a radius >= R over at least FOUR joints only, the least-squares circle of "
         "the "
         "joints (exact for an inscribed polyline; three points are always concyclic, so a "
         "3-joint least-squares bend would absorb any three same-sign joints — refused), "
         "each arm meeting the circle's tangent at its end joint within the joint noise "
         "rule (on the shorter of the arm piece and the first chord) or lying on the "
         "circle "
         "as a chord (a bend between non-tangent arms after a spline piece; an arc "
         "starting "
         "at a corner on its circle). Radius < R = one rounded corner (own chain, arms "
         "meet "
         "through it, AngleDegrees = 180 - total turn, CornerRadiusOverR from the tangent "
         "lengths; a rounded corner of small turn is a corner of large angle); radius >= R "
         "= a bend inside its chain. A closed path of joints turning one way through 360 "
         "deg on one circle within the fit tolerance with every joint below the cap is one "
         "arc (a round pad or hole; a square or hexagonal hole is corners whatever its "
         "size, an octagon is a circle). Every joint no arc absorbs and that is not noise "
         "is a corner: a two-joint chamfer, a square strip end, a 90 deg lead-end corner. "
         "Both traversal directions of a path are scanned and the set absorbing more "
         "joints "
         "wins, then fewer arcs, then the smaller serialisation of (radius, turn, joints, "
         "centre distance from the path centroid, first joint's distance from the nearer "
         "path end) on the signature grid (the pair of serialisations of the two scans is "
         "orientation invariant): a translated, rotated or mirrored mesh gives the "
         "congruent arcs; a closed path is scanned from the joint after its longest piece "
         "(ties: the first from the seed vertex): no arc is split by the loop start"},
        {"LengthQuantumOverR", kLengthQuantumOverRadius},
        {"DirectionQuantum", kDirectionQuantum},
        {"SignatureLengthQuantumOverR", kSignatureLengthQuantumOverRadius},
        {"SignatureAngleQuantumDegrees", kSignatureAngleQuantumDegrees},
        {"ClusterQuantumNearMatchMaxQuanta", kClusterQuantumNearMatchMaxQuanta},
        {"ClusterQuantumInclusiveMargin", kClusterQuantumInclusiveMargin},
        {"ClusterQuantumNearMatch",
         "block (b) DESIGN section 4 (decision 303): a SpatialEdgeCluster feature matches "
         "a library model whose topology key (the signature with every number of "
         "Portions / Context [].P / Arc / Gap, Box and Vertices[].P / TurnDegrees nulled, "
         "entry order preserved) equals its own and whose numbers lie within "
         "ClusterQuantumNearMatchMaxQuanta signature quanta of its own (read half-quantum "
         "inclusive, <= MaxQuanta + ClusterQuantumInclusiveMargin, decisions 287 / 288 / "
         "317; the nearest such "
         "model; ties by name); recorded Features[].Match.QuantumNearMatch and in Summary; "
         "a permuted entry order is another key (Missing); two library models within "
         "twice that many quanta (+ the margin) are refused at load"},
        {"StraightBendRadiusOverR", kStraightBendRadiusOverRadius},
        {"CurvatureWindowOverR", kCurvatureWindowOverRadius},
        {"PairSeparationToleranceRelative", kPairSeparationTolerance},
        {"PairSeparationSamplesPerInterval", kPairSeparationSamplesPerInterval},
        {"PairBendProximityOverR", kPairBendProximityOverRadius},
        {"PairChordWindowCapOverR", kPairChordWindowCapOverRadius},
        {"ExactStretchMinLengthOverR", kExactStretchMinLengthOverRadius},
        {"PairSeparationEstimate",
         "per sample: chord reading C = min over the two chains of the maximum sampled "
         "closest-point distance within max(R, min(local chord, PairChordWindowCapOverR "
         "R)) of the sample / its foot where the chain bends within PairBendProximityOverR "
         "R of it along the chain (a joint with a turn or a fitted bend arc; not the "
         "windowed curvature, which spreads a joint's turn over its two half-runs), R "
         "elsewhere; inscribed reading C / cos(turn / 2) with the larger local joint turn; "
         "interacting iff both < 2R (the decision). The feature separation (decision "
         "85(1), amended by USER decision 184 (3), 2026-10-01) per link and curvature "
         "class is EXACT where the geometry allows: the perpendicular distance of two "
         "exactly parallel straight runs neither of which bends within "
         "PairBendProximityOverR R of the sample / foot and whose foot is the "
         "perpendicular "
         "projection (closest-point distance = line distance within the parameter "
         "tolerance), or the radius difference of two fitted bend arcs whose centres "
         "coincide within the parameter tolerance. A sub-piece (consecutive samples of one "
         "constancy / interaction status along a run) is cut into STRETCHES (USER decision "
         "203, 2026-10-02): an exact stretch = a maximal contiguous stretch of samples "
         "with "
         "an exact reading whose readings lie within the parameter tolerance of the "
         "stretch's exact mean; it splits off as its own exact sub-piece (reading that "
         "mean) only when its exact samples span at least ExactStretchMinLengthOverR R, "
         "else it merges into the adjacent non-exact stretch (a sub-piece that is one "
         "exact "
         "stretch stays exact at any length); the judged samples without an exact reading "
         "(away from any bend within PairBendProximityOverR R) form non-exact stretches by "
         "status; the near-bend samples (not judged) join an adjacent stretch within "
         "PairBendProximityOverR R along the run of its last exact / judged sample (the "
         "exact side first; a short agreeing gap between two exact stretches is "
         "transparent) and farther than that form a non-exact stretch; a non-exact stretch "
         "is cut where its feet cross a partner run boundary (one non-exact sub-piece "
         "faces "
         "one partner run); a non-exact "
         "sub-piece reads the mean of the chord readings "
         "(equally spaced samples: the length-weighted mean over the sub-piece) and is not "
         "exact. The pieces of a chain "
         "pair are grouped into links so that a link holds ONE separation: exact pieces "
         "within the "
         "parameter tolerance of each other form an exact group (value = their "
         "length-weighted mean), a chord piece joins the exact group within the pair "
         "tolerance (PairSeparationToleranceRelative) of its own separation (nearest), the "
         "remaining chord pieces form chord groups by consecutive steps within the pair "
         "tolerance (value = the length-weighted mean chord reading, ExactParameters "
         "false); a piece is never keyed by the exact readings of pieces beyond its own "
         "tolerance (a 3.7 um strip next to a 2 um lead reads 3.7); never a "
         "sample-count-weighted mean over pieces. Invariant: mutually facing non-exact "
         "pieces (each the only non-exact piece of its run whose feet all lie on the "
         "other's run) belong to one group; where they straddle a group threshold the "
         "piece "
         "spanning exactly its own run decides (ties: the lower run index) and the other "
         "joins its group"},
        {"SignatureParameterToleranceOverR", kSignatureParameterToleranceOverRadius},
        {"SignatureAngleToleranceDegrees", kSignatureAngleToleranceDegrees},
        {"PortionMinimumLengthOverR", kSignatureParameterToleranceOverRadius},
        {"SliverRule",
         "no portion shorter than PortionMinimumLengthOverR x R (= the signature parameter "
         "tolerance) exists (supervisor decision 222, 2026-10-02): a portion is a maximal "
         "contiguous stretch of one feature side along a chain (across its runs and mesh "
         "segments); a shorter stretch is roundoff between two claim boundaries that no "
         "signature parameter resolves and joins its adjacent portion on the chain, the "
         "longer of its two neighbours (ties: the one before it along the chain), taking "
         "its feature and side; a pair / stack side surviving below it is no side (the "
         "feature dissolves); a stretch with no adjacent portion on its chain stays "
         "(Diagnostics.SubTolerancePortionsJoined.Isolated). Dimensionless; replaces the "
         "former join of pieces at or below the SignatureLengthQuantumOverR grid, the "
         "knife-edge the S1p window's 1.0229e-6 R cut-end sliver missed"},
        {"SignatureMatching",
         "topology key = the signature with every continuous parameter (OffsetOverR, "
         "SeparationOverR, RadiusOverR, CornerRadiusOverR; AngleDegrees, ArmAnglesDegrees) "
         "removed; the library groups instances of one topology whose parameters agree "
         "within the tolerances (single linkage over the distinct signatures, both "
         "orientations of a translational signature) into one coupon at the representative "
         "(midpoint of every parameter's range, in the orientation nearest to the "
         "lexicographically smallest instance); the matcher accepts a model of the same "
         "topology within the tolerance and takes the nearest (ties by model name); a "
         "SpatialEdgeCluster signature matches exactly"},
        {"PairConstancyWindowOverR", 1.0},
        {"PairSampleSpacingOverR", 0.5},
        {"PairCandidateReachOverR",
         kInteractionDistanceOverRadius * (1.0 + kPairSeparationTolerance)},
        {"SamplingMargins",
         "the bent-pair candidate reach 2R (1 + 0.05) and the mutual-sides "
         "facing test at the pair separation x 1.05 gather candidates / "
         "samples only; every interaction decision is the 3D distance "
         "strictly below 2R on the quantized grid"},
        {"SelfPairNeighbourhoodOverR", kSelfPairNeighbourhoodOverRadius},
        {"StackRule",
         "links (bent pairs) and translational spans sharing a run over a common "
         "interval are one cross-section: k >= 3 edges with consecutive "
         "separations < 2R are one ParallelEdgeCluster / "
         "CurvedParallelEdgeCluster (offsets from the consecutive links' "
         "separations of the cross-section's curvature class, gap pattern, "
         "conductors, bend radius), straight and along bends; the pairwise "
         "candidates inside it are superseded; every component (a lone pair "
         "included) is assembled per cross-section: a member taken by a cluster "
         "portion or a vertex window is recomposed out of the cross-section "
         "(the stack ends there); a chain folding back within 2R beyond pi R of "
         "arc length pairs with itself"},
        {"StackCompositionCap", kStackCompositionCap},
        {"ClusterExtensionWedgeCapDegrees", kClusterExtensionWedgeCapDegrees},
        {"ClusterExtensionClosureOverR", kSignatureParameterToleranceOverRadius},
        {"ClusterExtensionMaxPasses", kClusterExtensionMaxPasses},
        {"StackImageToleranceOverR", kStackImageToleranceOverRadius},
        {"ClusterExtensionRule",
         "every single-edge portion (the isolated / curved remainder after clusters, "
         "windows and pairs / stacks) within 2R (3D, strict) of a cluster's claimed "
         "perimeter or of a free vertex feature's window on another chain (own chain "
         "beyond "
         "pi R), faced ACROSS (perpendicular projection inside the claimed piece, extended "
         "by 2R tan(turn) at interior joints of the claimed chain interval with the turn "
         "capped at ClusterExtensionWedgeCapDegrees, never past the end of the claimed "
         "interval; the full 2R ball around a claimed piece ending where its chain ends), "
         "outside the through-vertex zones of vertices not in that cluster, joins the "
         "cluster (a vertex feature so joined becomes a cluster; several owners merge); "
         "iterated to closure with the stack recomposition, stopping after a pass that "
         "absorbs at most ClusterExtensionClosureOverR x R in total (that pass is applied, "
         "the stacks are not recomposed again: sub-tolerance) (ratified as the 'across' "
         "rule, decision 88(2)); pairs / stacks are never absorbed and their claimed "
         "length satisfying the same across rule is Diagnostics.StackEndThirdBodyLength "
         "(decision 85(2)) — with ONE exception (decision 224): a stretch of one "
         "cross-section's pair / stack claims along a chain that exists only because a "
         "cluster's claim boundary cut it — bounded at both ends by claims of the same "
         "cluster, or adjacent to the cluster at one end and continuing a larger stack "
         "containing its members at the other — and lying entirely within "
         "ClusterBallOverR x R of that cluster's claimed pieces is absorbed by that "
         "cluster (Diagnostics.ClusterExtension.TranslationalPiecesAbsorbed, by class "
         "TwoSided / StackEndRecomposition)"},
        {"MutualSidesOverhangOverSeparation",
         std::sqrt((1.0 + kPairSeparationTolerance) * (1.0 + kPairSeparationTolerance) -
                   1.0)},
        {"KnifeEdgeBandRelative", kKnifeEdgeBandRelative},
        {"KnifeEdgeSampleSpacingOverR", kKnifeEdgeSampleSpacingOverRadius},
        {"FacingGateExclusions",
         "audit A8 (facing_check.py): isolated / curved edges facing "
         "another edge and pair / stack sides facing a third edge "
         "within 2R (every facing segment tested; own sides = the "
         "portions + member chains within the feature's reach) must "
         "be zero apart from AtExactly2R, ThroughVertex, "
         "SelfNeighbourhood and, for pair / stack sides only, "
         "StackEndThirdBody (facing cluster / vertex metal) and "
         "StackEndRecomposition (facing a pair / stack sharing a "
         "member chain; bounded to 1 R per site, the operative "
         "bound, and 1 R per feature in total) and, for isolated / "
         "curved portions shorter than the signature parameter "
         "tolerance that are an unclaimed remainder bounded on "
         "both sides along their run by other features' claims, "
         "or a whole feature shorter than the tolerance, "
         "SubToleranceFeature (decision 89), each recorded with "
         "its length (decision 85(2): ClusterNeighbour / "
         "VertexNeighbour are gone)"},
        {"CrossLayerReachOverR", kInteractionDistanceOverRadius},
        {"PlaneRule", "features (pairs, clusters, vertex joins, translational classes) "
                      "never span two metal planes; metal of another plane within the "
                      "interaction distance is the CrossLayer exclusion"},
        {"PortRule", "metal perimeter bordering a LumpedPort / WavePort boundary face is "
                     "the Port exclusion; its vertices are PortCut, never features"},
        {"SupportContinuationOverR", kSupportContinuationOverRadius},
        {"SupportPaddingOverR", kSupportPaddingOverRadius},
        {"SupportEndCoincidenceOverR", kSupportEndCoincidenceOverRadius},
        {"SupportFaceSnapOverR", kSupportFaceSnapOverRadius},
        {"SupportFaceClearanceOverR", kSupportFaceClearanceOverRadius},
        {"SupportFaceGrowthStepOverR", kSupportFaceGrowthStepOverRadius},
        {"SupportFaceGrowthMaxSteps", kSupportFaceGrowthMaxSteps},
        {"SupportSpanCapOverR", kSupportSpanCapOverRadius},
        {"SpanCapAllowances",
         "block (b) DESIGN section 3 (a) / A4 (decision 303): a per-case allowance of the "
         "response correction config (SpatialSupport.SpanCapAllowances[] {ClaimsSignature, "
         "SpanCapOverR, Label, Reason, Approval}) replaces SupportSpanCapOverR for the "
         "cluster whose CLAIMS-ONLY signature near-matches the embedded ClaimsSignature "
         "within ClusterQuantumNearMatchMaxQuanta quanta (half-quantum inclusive; resolved "
         "before any box, so one "
         "entry serves every window instance); recorded Features[].SpatialSupport."
         "SpanCapAllowance {Label, SpanCapOverR, Reason, Approval, MatchedQuanta, "
         "DifferingNumbers}; a refused cluster exports SpatialSupport.ClaimsSignature; an "
         "allowance below the cap or two allowances within twice the quanta are refused at "
         "load; an allowance matching no cluster is listed under UnusedSpanCapAllowances"},
        {"SupportComparisonQuantumOverR", kSignatureLengthQuantumOverRadius},
        {"SpatialSupportContract",
         "contract v3 (USER decision 281, supervisor decision 282, 2026-10-03; replaces "
         "the "
         "decision-236 straight-continuation geometry): the coupon metal of a "
         "SpatialEdgeCluster is the DEVICE PLAN of the cluster's plane (every perimeter "
         "run, straight or on its fitted arc, with its conductor) clipped to the claims-"
         "derived support box in the cluster's canonical frame "
         "(Features[].SpatialSupport). "
         "Box (rule B2 = the coupon generator's coupon_bounds in the canonical frame, "
         "units "
         "of R): every claimed portion is a row along its own tangent; a claim-cut end "
         "(touching no vertex and no other portion end within SupportEndCoincidenceOverR) "
         "is lengthened to at least R from the row midpoint, every row end at or beyond R "
         "continues by SupportContinuationOverR R, the rows are widened by R on both "
         "sides and the bounding box padded by SupportPaddingOverR R; arc portions are "
         "chorded at ClusterArcChordStepDegrees / ClusterArcChordMaxLengthOverR as the "
         "builder chords them. Context: every run piece of the plane inside the box that "
         "is not a claim, in the portion encoding with Chain true on the pieces connected "
         "to the claims inside the box through run ends and device vertices (the own "
         "continuation chains, rule B3: the within-R accounting and the continuation "
         "ownership follow them) and false on foreign edges (present in both coupon twins, "
         "excluded from the within-R accounting; their patches untouched); conductor "
         "labels by first appearance over the sorted Portions THEN the sorted Context. "
         "Face rules, dimensionless: T1 a piece end within SupportFaceSnapOverR R of a "
         "face "
         "is ON the face and a context piece shorter than that is dropped (the sliver "
         "quantum); T2 every device edge inside the box keeps SupportFaceClearanceOverR R "
         "from every face it does not cross, every piece end not on a face (device vertex, "
         "claim end; an end made only by another cluster's claim cut is no vertex and is "
         "exempt) and every face crossing keep that clearance from every other face "
         "(a crossing near a box corner), every crossing meets its face at sin(theta) >= "
         "SupportFaceClearanceOverR, and two crossings of one face closer than that may "
         "bound metal (a narrow lead: allowed, FaceRules.MinCrossSectionOverR) but not "
         "gap; every length threshold (snap, clearance, separation) is read on the "
         "SupportComparisonQuantumOverR grid (decision 287 (b): a value within half a "
         "quantum of the threshold is AT it and takes the rule's inclusive side, so an "
         "edge on a claims-box face reads exactly the clearance from the moved face and "
         "passes deterministically; keys are geometry-precise to the quantum and "
         "sub-quantum noise cannot move them); "
         "T3 every face failing T2 moves outward by SupportFaceGrowthStepOverR R per step "
         "(all failing faces of a step together), T1 / T2 re-evaluated on the grown box "
         "(new edges enter), at most SupportFaceGrowthMaxSteps steps per face and never "
         "past the plan span cap SupportSpanCapOverR R (a claims-derived box already "
         "beyond the cap is keyed and recorded ExceedsSpanCap — the builder's cap and its "
         "per-case override decide as before — but may not grow); a cluster no box "
         "satisfies is an UnboxableFeature: its signature carries Unboxable true (a "
         "Missing placeholder no builder makes, never a knife-edge coupon) and the record "
         "names the face and the reason. Key (ruling on R0 review MAJOR-1, option (i)): "
         "Box + Context enter the hashed signature ONLY when the context is non-empty or "
         "the box grew (Contract 3; the serialisation {Box, Context, Portions, Vertices} "
         "minimised over the candidate frames of the CLAIMS, both handedness values, with "
         "the box and context of each frame — the box frame IS the canonical frame, ruling "
         "MAJOR-2 — so mirror images keep one key with opposite Chirality and rotated / "
         "translated copies the same key); otherwise the key is the claims-only key of "
         "contract 2 byte for byte (unchanged coupons keep their keys; no library "
         "migration); the record is written in every case. The straight continuation of "
         "the "
         "decision-236 contract and the part of it lying on no device edge are recorded "
         "per feature (LegacyContinuation: the fictitious island metal of D2 / D3-C) and "
         "summed "
         "under Diagnostics.SpatialSupport; vertex features (corners) inside a box on a "
         "chain are listed (Context.ChainVertices: [x, y, Type, distance from the nearest "
         "face / R]) for the placement's vertex ownership (rule B4; a corner closer than R "
         "to a face has an arm partly outside the box); their signatures are unchanged"},
        {"Comparison", "strict less on the quantized grid"},
        {"MirrorRule",
         "boundary-cut DESIGN 2.2 (decisions 442 / 454 / 481 / 512): the perimeter within "
         "BandOverR x R of every planar NATURAL vertical truncation plane is reflected "
         "into "
         "the identification input and the result merged onto the unextended run; a "
         "straight joint (collinear within the direction quantum or a sub-noise turn) "
         "continues its chain (Mirror Status Continued: the real feature verbatim, "
         "whenever every real portion of the extended chain is identified identically - "
         "type, key, side, turn within the joint noise rule - whatever the split of the "
         "portions and whether the real chain extends further), an "
         "oblique meeting is a corner of 2 theta, a parallel edge at d < R a strip / gap "
         "of "
         "2 d (Modelled: placed on the real half - a vertex coupon ON the plane at weight "
         "1 "
         "/ 2 with the real arm's cells from s_half = (R + s) / 2, a pair's real side with "
         "its side factor - with the trace by even extension; Missing: the real portions "
         "kept raw as <Type>:MirrorFormed uncovered portions); a mirror-formed stack, "
         "cluster or curved pair has no mirror placement (Unmerged) and is read exactly "
         "as a real Missing feature (decision 481): its OWN cells - the patches whose "
         "own-edge footprint overlaps one of its REAL portions whose identification "
         "DIFFERS between the unextended and the extended run (decision 512 (c)) - are "
         "DomainBoundary with their raw claims (Reason UnmergedTopology), identically "
         "identified cells and neighbouring cells keep their "
         "Applied / Mirrored classification, the coupon supports' reach into the "
         "configuration is recorded as information (UnmergedSupportReach); image "
         "perimeter is never placed or counted"}}},
      {"ReferenceProcessNormal", D(reference_process_normal)},
      {"Features", feature_list},
      {"Segments", segment_list},
      {"Vertices", vertex_list},
      {"Arcs", arc_list},
      {"MeshCoarsenessWarning", mesh_coarseness_warning},
      {"Exclusions", exclusion_list},
      {"Totals",
       {{"PerimeterLength", L(perimeter_length)},
        {"AssignedLength", L(assigned_length)},
        {"ExcludedLength", L(excluded_length)}}},
      {"Diagnostics",
       {{"SamePriorityClaimOverlaps", same_priority_claim_overlaps},
        {"Rule", "claims of one priority by different features never overlap on a run "
                 "(clusters disjoint, windows abut, pairs / stacks assembled per "
                 "cross-section): claim resolution never decides by feature id"},
        {"StackGeometricOffsetIntervals", stack_geometric_offsets},
        {"StackCompositionCapHits", stack_composition_cap_hits},
        {"StackCompositionCap", kStackCompositionCap},
        {"StackImagesMerged", stack_images_merged},
        {"StackImageToleranceOverR", kStackImageToleranceOverRadius},
        {"SubTolerancePortionsJoined",
         {{"Count", sub_tolerance_portions.count},
          {"Length", L(sub_tolerance_portions.length)},
          {"MaxLength", L(sub_tolerance_portions.max_length)},
          {"Isolated", sub_tolerance_portions.isolated},
          {"ToleranceOverR", kSignatureParameterToleranceOverRadius},
          {"Rule",
           "sliver rule (supervisor decision 222, 2026-10-02): no portion shorter than "
           "SignatureParameterToleranceOverR x R exists — a portion being a maximal "
           "contiguous stretch of one feature side along a chain, across its runs and mesh "
           "segments (a sub-tolerance mesh segment inside a long edge continues that "
           "edge); "
           "a shorter stretch (roundoff between two claim boundaries: a cluster ball "
           "cutting "
           "a pair piece, a claim ending next to a run end, the foot of a neighbour's cut "
           "end on a near-parallel member) joins its adjacent portion on the chain, the "
           "longer of its two neighbours (ties: the one before it along the chain), taking "
           "its feature and side; Count / Length / MaxLength are the stretches joined, "
           "Isolated the stretches with no adjacent portion on their chain (kept). "
           "Replaces "
           "the former join of pieces at or below the 1e-6 R signature grid"}}},
        {"ClusterExtension",
         {{"Passes", extension.passes},
          {"MaxPasses", kClusterExtensionMaxPasses},
          {"CapReached", extension.cap_reached},
          {"RepeatDetected", extension.repeat_detected},
          {"AbsorbedPortions", extension.portions},
          {"AbsorbedLength", L(extension.length)},
          {"VertexFeaturesJoined", extension.sites},
          {"TranslationalPiecesAbsorbed",
           {{"Count", extension.translational_pieces},
            {"Length", L(extension.translational_length)},
            {"MaxLength", L(extension.translational_max_length)},
            {"TwoSided",
             {{"Count", extension.translational_two_sided},
              {"Length", L(extension.translational_two_sided_length)}}},
            {"StackEndRecomposition",
             {{"Count", extension.translational_pieces - extension.translational_two_sided},
              {"Length", L(extension.translational_length -
                           extension.translational_two_sided_length)}}},
            {"Rule",
             "decision 224 (2026-10-02), the ONE exception to 'pairs / stacks are never "
             "absorbed': pair / stack stretches that exist ONLY because a cluster's claim "
             "boundary cut them are absorbed by that cluster. A stretch is a maximal "
             "contiguous run of one cross-section's claims (one signature key) along a "
             "chain; it must lie entirely within ClusterBallOverR x R (3D) of the "
             "cluster's "
             "claimed pieces (shorter than the interaction distance from the cluster's own "
             "claims: the cluster's coupon describes it, the stack's translational patches "
             "would correct the same surface twice) and be either TwoSided = bounded at "
             "both ends by claims of the SAME cluster (a piece shorter than 2R between two "
             "of its claims) or a StackEndRecomposition piece = adjacent to the cluster at "
             "one end and continuing, at the other, a pair / stack claim of a "
             "cross-section "
             "with strictly more members containing its own (the members the cluster takes "
             "later, continuing past the member it takes first). Never between two "
             "different clusters (the placement's spatial-vs-spatial overlap check covers "
             "overlapping coupon boxes), never a stretch whose far end is free, a bend or "
             "a smaller cross-section: strips and gaps alongside a cluster are genuine "
             "translational features. Applied per pass like the single-edge absorptions, "
             "the stacks recomposed around the enlarged claims; a stretch the rule leaves "
             "inside a spatial coupon's volume is recorded by the placement's ownership "
             "record (decision 236: Continuation / Foreign, never an abort)"}}},
          {"Rule",
           "every single-edge portion (isolated / curved edge remainder) within 2R "
           "(3D, strict) of a cluster's claimed perimeter or of a vertex feature's "
           "window on another chain (or its own chain beyond pi R), faced across "
           "(Conventions.ClusterExtensionRule), outside the through-vertex zones of "
           "non-member vertices, joins that cluster (a vertex feature so joined "
           "becomes a cluster); iterated to closure over single-edge portions (the "
           "loop stops after a pass absorbing at most ClusterExtensionClosureOverR x "
           "R, after a pass repeating the previous one (same absorbed length and "
           "portion count) or after MaxPasses passes (decision 93), the last pass "
           "applied without a further recomposition); pairs / stacks are never "
           "absorbed (decision 85(2), across rule ratified 88(2)) — with ONE exception "
           "(decision 224): the claim-boundary artefacts of TranslationalPiecesAbsorbed "
           "(TwoSided / StackEndRecomposition within ClusterBallOverR x R of the "
           "cluster's claims)"}}},
        {"SpatialSupport",
         {{"Clusters", spatial_support.clusters},
          {"ClaimsKeyed", spatial_support.claims_keyed},
          {"ContextKeyed", spatial_support.context_keyed},
          {"Grown", spatial_support.grown},
          {"Unboxable", spatial_support.unboxable},
          {"ExceedingSpanCap", spatial_support.exceeding_span_cap},
          {"ChainLength", L(spatial_support.chain_length)},
          {"ForeignLength", L(spatial_support.foreign_length)},
          {"FictitiousContinuationLength",
           L(spatial_support.fictitious_continuation_length)},
          {"ThresholdBandHits", spatial_support.threshold_band_hits},
          {"Rule", "Conventions.SpatialSupportContract: clusters keyed by their claims "
                   "alone (contract 2: empty context, no growth), with Box + Context "
                   "(contract 3), grown by T3, unboxable (Missing placeholder keys) and "
                   "with a claims-derived box beyond the span cap; the context lengths "
                   "(own chains / foreign edges inside the boxes) and the legacy "
                   "contract's fictitious straight continuation, in the manifest's length "
                   "unit; ThresholdBandHits counts the face-rule readings within "
                   "KnifeEdgeBandRelative of a threshold (snap, clearance, crossing "
                   "sine, crossing separation) over EVERY T2 pass of every cluster plus "
                   "the box-rule readings (decision 287 (a): a key whose growth sequence "
                   "was decided at a threshold is flagged; the per-pass readings are in "
                   "Features[].SpatialSupport.Growth.Passes)"}}},
        {"StackEndThirdBodyLength", L(extension.stack_end_third_body_length)},
        {"StackEndThirdBodyRule",
         "pair / stack claimed length within 2R (3D, strict) of a cluster's claimed "
         "perimeter or a free vertex feature's window, faced across by the cluster "
         "extension's rule, outside the through-vertex zones of non-member vertices "
         "(the same predicate as the extension, evaluated on the pair / stack claims and "
         "never absorbed); the facing gate's StackEndThirdBody class is the sampled "
         "reading of the same definition (0.5 R samples, across = the facing direction "
         "within 60 deg of the sample's normal, through-vertex zones excluded first)"}}},
      {"KnifeEdgeCensus", ScaledCensus(knife_edge_census, length_scale)},
      {"UnusedSpanCapAllowances", unused_span_cap_allowances},
      {"GeometryDigest", geometry_digest}};
  if (!image_segment_list.empty())
  {
    // Mirror band: the image segments (never perimeter) listed apart.
    manifest["ImageSegments"] = image_segment_list;
  }
  return manifest;
}

}  // namespace palace
