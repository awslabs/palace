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
constexpr double kCornerTurnToleranceDegrees = 30.0;
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
// Arc rule (decision 82(3), design (b) 4 / 7): consecutive joints of the perimeter path
// (through corner vertices) connected by pieces shorter than the interaction distance 2R
// and turning the same way are one arc when the arm tangents and every joint vertex fit one
// circle within kArcFitToleranceRelative (tangent lengths equal, radial deviation of every
// vertex below the tolerance x radius). An arc of radius below R whose total turn exceeds
// the corner threshold is ONE rounded corner (total turn, radius / R in the signature); an
// arc of radius >= R is a bend of exactly that radius for the curved-edge chain rule (its
// joints, corner vertices included, are no features and its density is 1 / radius over the
// arc) — the description does not depend on the number of chords.
constexpr double kArcFitToleranceRelative = 0.05;
constexpr double kRoundedCornerTangentTolerance = kArcFitToleranceRelative;

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
// kPairSeparationTolerance x their minimum (a polyline of sub-corner turns <= 30 deg at
// constant width varies by 1 / cos(15 deg) - 1 = 3.5 %; the pair response sensitivity d
// dR/dd is O(1)). The constant portions are the pair; the portions that are not (divergence
// at tees and port ends, fast tapers, acute corner arms) keep the event rule, so a slow
// taper is a pair and a tee is a cluster. Whether a constant portion interacts is decided
// on the separation of the underlying curves, so that a discretisation never changes the
// classification: the chords of a polyline inscribed in a curve lie inside it (two
// concentric inscribed polylines are w cos(turn / 2) apart mid-chord and exactly w apart at
// their vertices), while the exact offset polyline of a bent path keeps corresponding
// chords at the design separation and the outer side's samples near the joints project onto
// the inner vertices at up to w / cos(turn / 2). In both constructions the sampled
// closest-point distance from one chain to the other reaches the curve separation w as its
// maximum on the side whose maximum is smaller: the sample's chord reading C = min over the
// two chains of the maximum sampled distance within a window of half-width max(R, local
// chord) where the chain bends and R on straight runs, about the sample on its own chain
// and about its foot on the other chain (a straight taper is read locally). C is exact for
// an offset polyline; two polylines inscribed in the curves at aligned angles are C = w
// cos(turn / 2) apart everywhere (chords and vertex-to-polyline alike) for a curve
// separation w, and the polyline pair alone cannot tell the two constructions apart (they
// differ at order turn^2): the inscribed reading is C / cos(turn / 2) with turn = the
// larger local joint turn of the two chains. A portion interacts iff BOTH readings are
// below 2R on the quantized grid — the same strict-less decision as a straight parallel
// pair at that separation, taken on the non-interacting side of the recorded ambiguity w (1
// / cos(turn / 2) - 1) (below 1e-4 w for joints under 1.6 deg; a CPW gap of exactly 2R
// along a bend is isolated edges like a straight one; DS-SCT-001's 4 um gaps at R = 2 um
// read 3.9998 mid-chord and became 3 mm clusters). The pair feature's separation is the
// mean chord reading over its samples; the cross-chord interactions of a locally constant
// portion are never event cores, whether or not it interacts.
constexpr double kStraightBendRadiusOverRadius = 10.0;
constexpr double kCurvatureWindowOverRadius = 1.0;
constexpr double kPairSeparationTolerance = 0.05;
constexpr int kPairSeparationSamplesPerInterval = 16;
// Self-pairing (decision 82(2) addition, 2026-09-25): a chain folding back onto itself
// within 2R through a bend of radius >= R (a hairpin, a meander with smooth bends, the
// U-turn of a narrow strip whose fold is not a rounded corner) pairs with itself. The
// local neighbourhood along the chain takes no part: two points of one chain closer than
// pi R along the chain are within 2R of each other along any bend of radius >= R (on the
// tightest bend, radius R, the chord reaches 2R exactly after half a turn = pi R of arc;
// on a wider bend earlier), so only points more than pi R apart along the chain can face
// each other across a fold. Runs closer than pi R of arc length are never partners.
constexpr double kSelfPairNeighbourhoodOverRadius = 3.14159265358979323846;

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

CanonicalSignature
CanonicalClusterSignature(const std::vector<SignaturePortion> &portions,
                          const std::vector<SignatureVertex> &vertices,
                          const Point3D &process_normal, double radius,
                          const std::function<void(std::size_t, std::size_t)> &progress)
{
  MFEM_VERIFY(!portions.empty(), "A cluster signature needs at least one edge portion!");
  const Point3D n = Normalize(process_normal);
  // Origin: length-weighted centroid of the portions.
  Point3D origin{};
  double total = 0.0;
  for (const auto &portion : portions)
  {
    const double length = Distance(portion.p0, portion.p1);
    origin = Add(origin, Scale(0.5 * length, Add(portion.p0, portion.p1)));
    total += length;
  }
  MFEM_VERIFY(total > 0.0, "A cluster signature needs portions of positive length!");
  origin = Scale(1.0 / total, origin);
  // Candidate in-plane axes: portion tangents and their perpendiculars, both signs.
  std::vector<Point3D> candidates;
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
    const Point3D t = Normalize(Sub(portion.p1, portion.p0));
    for (const double sign : {1.0, -1.0})
    {
      AddCandidate(Scale(sign, t));
      AddCandidate(Scale(sign, Cross(n, t)));
    }
  }
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
    : faces(input_.faces), quantizer(input_.radius), R(input_.radius)
  {
    // Working copy of the segments and vertices: the arc rule (design (b) 4 / 7, decision
    // 82(3)) merges the chains through the corner vertices it absorbs and demotes those
    // vertices to regular; the faces are read only.
    input.radius = input_.radius;
    input.segments = input_.segments;
    input.vertices = input_.vertices;
    input.log = input_.log;
  }

  IdentificationResult Identify();

private:
  IdentificationInput input;
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
    double turn = 0.0;    // radians, total
    bool corner = false;  // radius < R and turn > the corner threshold: a rounded corner
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
  // shared vertex and of every run of an arc joining the two chains.
  std::vector<std::vector<Interval>> ThroughZones(std::size_t a, const Chain &A,
                                                  const Chain &B) const
  {
    std::vector<std::vector<Interval>> zones;
    for (const std::size_t v : SharedVertices(A, B))
    {
      const Point3D &p = input.vertices[v].coordinate;
      zones.push_back(RunIntervalWithin(a, p, p, kThroughVertexZoneOverRadius * R));
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
          const auto within = RunIntervalWithin(a, runs[r].start, runs[r].end,
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
    double weighted_separation;
    int sample_count;
  };
  struct PairLink
  {
    int chain_a = -1, chain_b = -1;
    double separation = 0.0;  // mean chord reading over the samples
    Point3D lead_point{}, lateral_ab{};
    std::size_t run_a = 0, run_b = 0;  // the runs at the lead pair of points
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
  void ClassifyVertices();
  void BuildArcSites();
  void ComputeCurvature();
  void BuildTranslationalFeatures();
  void BuildBentPairs();
  void BuildPairsAndStacks();
  void EmitPairLink(const PairLink &link);
  void EmitTranslationalSpan(const TranslationalSpan &span);
  void AssembleStack(const std::vector<std::size_t> &link_items,
                     const std::vector<std::size_t> &span_items);
  void BuildClusters();
  void BuildVertexWindows();
  void Assign(IdentificationResult &result);
  std::vector<Interval> RunIntervalWithin(std::size_t run, const Point3D &a,
                                          const Point3D &b, double distance) const;
  std::vector<Interval>
  ChainWindow(std::size_t run, std::size_t from_vertex, double length,
              std::vector<std::pair<std::size_t, Interval>> &out) const;

  // Curved-edge chain rule helpers on the chain arc length x.
  double WindowedCurvature(const Chain &chain, double x) const;
  double MaxCurvature(const Chain &chain, double x0, double x1) const;
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
      // Arm angles around the process normal, as consecutive differences.
      const Point3D x =
          Normalize(Sub(ArmDirection(incident[0], v),
                        Scale(Dot(ArmDirection(incident[0], v), n_ref), n_ref)));
      const Point3D y = Cross(n_ref, x);
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

// Arc rule (decision 82(3)): fitted arcs along the perimeter path through corner vertices.
// A path continues through every vertex with exactly two path segments (regular or corner);
// it stops at endpoints, junctions and cuts. Its joints are the non-collinear vertices;
// from every unconsumed joint the largest range of following joints that (i) are connected
// by pieces shorter than 2R, (ii) turn the same way, (iii) total at most 180 deg and (iv)
// fit one circle with the two arm tangents is an arc. The corner vertices of an arc make no
// vertex feature (a 90 deg fillet meshed with two 45 deg chords has one in its middle), the
// arc's curvature is its exact radius and its curved features are one across those
// vertices, so that the description does not depend on the number of chords.
void Identifier::DetectArcs()
{
  const double interaction = kInteractionDistanceOverRadius * R;
  const double corner_turn = kCornerTurnToleranceDegrees * std::acos(-1.0) / 180.0;
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
      const Point3D normal = Normalize(
          Add(input.segments[before].process_normal, input.segments[after].process_normal));
      const double sign = Dot(Cross(in, out), normal);
      joints.push_back({vertex, k, in, out, std::acos(dot), sign >= 0.0 ? 1 : -1});
    }
    if (joints.size() < 2)
    {
      continue;
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
      // Start the cyclic scan at the joint following the longest piece (a straight arm), so
      // that no arc is split by the arbitrary loop start.
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
    auto Gap = [&](std::size_t ja, std::size_t jb)
    {
      // Path length from joint ja to the following joint jb (cyclic on a loop).
      double d = position[joints[jb].index] - position[joints[ja].index];
      if (d < 0.0)
      {
        d += path_length;
      }
      return d;
    };
    const std::size_t m = joints.size();
    // Fit the arc over joints [i, k] (cyclic indices on a loop); returns the circle.
    struct Fit
    {
      bool ok = false;
      Point3D center{}, origin{};
      double radius = 0.0;
      double turn = 0.0;
    };
    auto TryFit = [&](std::size_t i, std::size_t count) -> Fit
    {
      Fit fit;
      const Joint &first = joints[i];
      const Joint &last = joints[(i + count - 1) % m];
      const Point3D ta = first.in, tb = last.out;
      const Point3D Ta = input.vertices[first.vertex].coordinate;
      const Point3D Tb = input.vertices[last.vertex].coordinate;
      double turn = 0.0;
      for (std::size_t j = 0; j < count; j++)
      {
        turn += joints[(i + j) % m].turn;
      }
      const double angle = std::acos(std::clamp(Dot(ta, tb), -1.0, 1.0));
      // The arm tangents must enclose the accumulated turn (a monotone arc below 180 deg).
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
      double radius = 0.0;
      Point3D center{}, origin{};
      const double sin_turn = std::sin(turn);
      if (std::abs(sin_turn) > 1.0e-9 && turn < std::acos(-1.0) - 1.0e-9)
      {
        // Virtual corner X = Ta + a ta = Tb - b tb, solved in the (ta, na) frame.
        const Point3D w = Sub(Tb, Ta);
        const double wx = Dot(w, ta), wy = Dot(w, na);
        const double tbx = Dot(tb, ta), tby = Dot(tb, na);
        if (std::abs(tby) <= 1.0e-12)
        {
          return fit;
        }
        const double b = wy / tby;
        const double a = wx - b * tbx;
        if (!(a > 0.0 && b > 0.0) ||
            std::abs(a - b) > kArcFitToleranceRelative * std::max(a, b))
        {
          return fit;
        }
        radius = 0.5 * (a + b) / std::tan(0.5 * turn);
        origin = Add(Ta, Scale(a, ta));
        center = Add(Ta, Scale(radius, na));
      }
      else
      {
        // Antiparallel arms (a U-turn): the radius is half the arm separation and the
        // tangent points face each other across it.
        const Point3D w = Sub(Tb, Ta);
        const double along = Dot(w, ta), across = Dot(w, na);
        if (across <= 0.0 || std::abs(along) > kArcFitToleranceRelative * across)
        {
          return fit;
        }
        radius = 0.5 * across;
        center = Add(Ta, Scale(radius, na));
        origin = Add(center, Scale(radius, ta));  // the arc's midpoint
      }
      if (!(radius > 0.0))
      {
        return fit;
      }
      for (std::size_t j = 0; j < count; j++)
      {
        const Point3D p = input.vertices[joints[(i + j) % m].vertex].coordinate;
        if (std::abs(Distance(p, center) - radius) > kArcFitToleranceRelative * radius)
        {
          return fit;
        }
      }
      fit.ok = true;
      fit.center = center;
      fit.origin = origin;
      fit.radius = radius;
      fit.turn = turn;
      return fit;
    };
    std::vector<bool> consumed(m, false);
    const std::size_t first_start = 0;
    for (std::size_t i = first_start; i < m; i++)
    {
      if (consumed[i])
      {
        continue;
      }
      // Extend while the pieces are short, the turn monotone and at most 180 deg.
      std::size_t best_count = 0;
      Fit best;
      double turn = joints[i].turn;
      // At least three joints (two chords): a polyline with two joints is a chamfer or a
      // square strip end, which no test can tell from a one-chord arc (the diameter chord
      // of a semicircle is the square end): those stay corners.
      for (std::size_t count = 2; count <= m; count++)
      {
        const std::size_t prev = (i + count - 2) % m, k = (i + count - 1) % m;
        if ((!closed && i + count - 1 >= m) || k == i || consumed[k] ||
            joints[k].sign != joints[i].sign || !quantizer.Less(Gap(prev, k), interaction))
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
      arc.corner = quantizer.Less(arc.radius, R) &&
                   arc.turn * 180.0 / std::acos(-1.0) > kCornerTurnToleranceDegrees;
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

double Identifier::ChainArcDistance(const Chain &chain, double x0, double x1, double x) const
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
  const std::size_t k = std::min(
      static_cast<std::size_t>(
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

Identifier::ChainPoint Identifier::ClosestPointOnChain(const Chain &chain,
                                                       const Point3D &p,
                                                       std::optional<double> exclude_x,
                                                       std::optional<double> max_distance) const
{
  const double neighbourhood = kSelfPairNeighbourhoodOverRadius * R;
  auto Excluded = [&](std::size_t k)
  {
    return exclude_x &&
           quantizer.Less(ChainArcDistance(chain, chain.run_offset[k],
                                           chain.run_offset[k] + runs[chain.runs[k]].length,
                                           *exclude_x),
                          neighbourhood);
  };
  ChainPoint best;
  best.distance = std::numeric_limits<double>::infinity();
  auto Evaluate = [&](std::size_t k)
  {
    if (Excluded(k))
    {
      return;
    }
    const Run &r = runs[chain.runs[k]];
    const double s = std::clamp(Dot(Sub(p, r.start), r.tangent), 0.0, r.length);
    const double distance = Distance(p, r.At(s));
    if (distance < best.distance)
    {
      best = {chain.runs[k], s, chain.run_offset[k] + s, distance};
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
    chain.curved.clear();
    chain.curved_max_kappa.clear();
    if (m < 2 || chain.length <= 0.0)
    {
      continue;
    }
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
    };
    std::vector<Piece> contributions;
    for (std::size_t k = 0; k < m; k++)
    {
      if (chain.joint_turn[k] <= 0.0 || chain.joint_excluded[k])
      {
        continue;
      }
      const std::size_t prev = (k + m - 1) % m;
      const double half_prev = 0.5 * runs[chain.runs[prev]].length;
      const double half = 0.5 * runs[chain.runs[k]].length;
      const double window = half_prev + half;
      if (window <= 0.0)
      {
        continue;
      }
      const double density = chain.joint_turn[k] / window;
      const double x = chain.run_offset[k];
      if (x - half_prev < 0.0)  // the joint before run 0 of a closed chain wraps
      {
        contributions.push_back({x - half_prev + chain.length, chain.length, density});
        contributions.push_back({0.0, x + half, density});
      }
      else
      {
        contributions.push_back({x - half_prev, x + half, density});
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
      for (const std::size_t s : arc.segments)
      {
        const long long int r = run_of_segment[s];
        if (r < 0 || runs[static_cast<std::size_t>(r)].chain != chain.id)
        {
          continue;
        }
        const std::size_t k = runs[static_cast<std::size_t>(r)].index_in_chain;
        for (const auto &rs : runs[static_cast<std::size_t>(r)].segments)
        {
          if (rs.segment == s)
          {
            spans.emplace_back(chain.run_offset[k] + rs.t0, chain.run_offset[k] + rs.t1);
          }
        }
      }
      for (const auto &span : MergeIntervals(std::move(spans), Tol()))
      {
        contributions.push_back({span.first, span.second, 1.0 / arc.radius});
        chain.arc_spans.emplace_back(a, span.first, span.second);
      }
    }
    // Sweep over the contribution ends (+density at x0, -density at x1).
    std::vector<std::pair<double, double>> events = {{0.0, 0.0}, {chain.length, 0.0}};
    for (const auto &c : contributions)
    {
      events.emplace_back(std::clamp(c.x0, 0.0, chain.length), c.density);
      events.emplace_back(std::clamp(c.x1, 0.0, chain.length), -c.density);
    }
    std::sort(events.begin(), events.end());
    std::vector<Piece> pieces;
    double density = 0.0;
    for (std::size_t i = 0; i < events.size(); i++)
    {
      density += events[i].second;
      if (i + 1 < events.size() && events[i + 1].first > events[i].first)
      {
        pieces.push_back({events[i].first, events[i + 1].first, std::max(0.0, density)});
      }
    }
    std::vector<double> cumulative(pieces.size() + 1, 0.0);
    for (std::size_t i = 0; i < pieces.size(); i++)
    {
      cumulative[i + 1] = cumulative[i] + pieces[i].density * (pieces[i].x1 - pieces[i].x0);
    }
    const double total_turn = cumulative.back();
    auto F = [&](double x)
    {
      double shift = 0.0;
      if (chain.closed)
      {
        const double periods = std::floor(x / chain.length);
        x -= periods * chain.length;
        shift = periods * total_turn;
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
      return cumulative[i] +
             pieces[i].density * (std::min(x, pieces[i].x1) - pieces[i].x0) + shift;
    };
    auto Kappa = [&](double x) { return (F(x + 0.5 * W) - F(x - 0.5 * W)) / W; };
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
    // Window half-widths for the curve separation: a full chord where the chain bends
    // (windowed curvature > 0: the chords of an inscribed polyline dip mid-chord), R on a
    // straight run (no dip; a taper must be read locally).
    double half_own, half_other;
    // The larger joint turn (radians) at the ends of the sample's run and of the foot's
    // run: the chord reading C and the inscribed-vertex reading C / cos(turn / 2) of the
    // pair separation differ by this discretisation ambiguity.
    double turn;
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
  struct Piece
  {
    std::size_t run;
    Interval interval;
    bool curved;
    double max_kappa;
    std::vector<Sample> samples;  // ordered along the run
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
    auto SelfExcluded = [&](std::size_t ka, std::size_t kb)
    {
      if (!self)
      {
        return false;
      }
      const double a0 = A.run_offset[ka], a1 = a0 + runs[A.runs[ka]].length;
      const double b0 = A.run_offset[kb], b1 = b0 + runs[A.runs[kb]].length;
      const double gap = std::min(ChainArcDistance(A, a0, a1, b0), ChainArcDistance(A, a0, a1, b1));
      return quantizer.Less(gap, neighbourhood);
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
        if (rb.excluded || SelfExcluded(ka, kb) ||
            !quantizer.Less(SegmentSegmentDistance(ra.start, ra.end, rb.start, rb.end),
                            reach))
        {
          continue;
        }
        const auto found = RunIntervalWithin(a, rb.start, rb.end, reach);
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
          const auto q_mid = ClosestPointOnChain(
              B, ra.At(0.5 * (piece.first + piece.second)),
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
            const auto q = ClosestPointOnChain(B, ra.At(s),
                                               self ? std::optional<double>(x) : std::nullopt,
                                               self ? std::optional<double>(reach) : std::nullopt);
            if (!std::isfinite(q.distance))
            {
              continue;
            }
            const double half_own =
                WindowedCurvature(A, x) > 0.0 ? std::max(R, ra.length) : R;
            const double half_other =
                WindowedCurvature(B, q.x) > 0.0 ? std::max(R, runs[q.run].length) : R;
            const double turn =
                std::max(LocalTurn(A, ka), LocalTurn(B, RunIndexInChain(B, q.run)));
            result.samples.push_back(
                {a, s, x, q.distance, q.x, half_own, half_other, turn});
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
    auto WindowMax =
        [&](const std::vector<const Sample *> &list, double x, double half, double *min_out)
    {
      Sample probe{};
      probe.x = x - half;
      auto lo = std::lower_bound(list.begin(), list.end(), &probe, ByX);
      probe.x = x + half;
      auto hi = std::upper_bound(list.begin(), list.end(), &probe, ByX);
      double max_d = 0.0, min_d = std::numeric_limits<double>::infinity();
      for (auto it = lo; it != hi; ++it)
      {
        max_d = std::max(max_d, (*it)->d);
        min_d = std::min(min_d, (*it)->d);
      }
      if (min_out)
      {
        *min_out = min_d;
      }
      return max_d;
    };
    for (auto &piece : own)
    {
      for (const auto &sample : piece.samples)
      {
        double min_d = 0.0;
        const double max_own = WindowMax(own_samples, sample.x, R, &min_d);
        constant.push_back(
            !quantizer.Less(kPairSeparationTolerance * min_d, max_own - min_d));
        const double w_own = WindowMax(own_samples, sample.x, sample.half_own, nullptr);
        const double w_other =
            other_samples.empty()
                ? w_own
                : WindowMax(other_samples, sample.qx, sample.half_other, nullptr);
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
  // interacting) flags; the cut between two samples of different status is halfway.
  struct SubPiece
  {
    std::size_t run;
    Interval interval;
    bool curved;
    double max_kappa;
    bool constant;
    bool interacting;
    double weighted_separation;  // sum over samples of the separation, for the mean
    int sample_count;
  };
  auto Split = [&](const std::vector<Piece> &pieces, const std::vector<char> &constant,
                   const std::vector<double> &separation, const std::vector<double> &upper,
                   std::vector<SubPiece> &out)
  {
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
        double weighted = 0.0;
        while (j < n && constant[k + j] == c &&
               quantizer.Less(upper[k + j], interaction) == inter)
        {
          weighted += separation[k + j];
          j++;
        }
        const double s_lo = i == 0 ? piece.interval.first
                                   : 0.5 * (piece.samples[i - 1].s + piece.samples[i].s);
        const double s_hi = j == n ? piece.interval.second
                                   : 0.5 * (piece.samples[j - 1].s + piece.samples[j].s);
        if (s_hi - s_lo > Tol())
        {
          out.push_back({piece.run,
                         {s_lo, s_hi},
                         piece.curved,
                         piece.max_kappa,
                         c,
                         inter,
                         weighted,
                         static_cast<int>(j - i)});
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
      if (A.Rigid() && B.Rigid() &&
          !DirectionLess(
              std::abs(Dot(runs[A.runs.front()].tangent, runs[B.runs.front()].tangent)),
              1.0 - kParallelCosineTolerance))
      {
        continue;  // two exactly parallel straight runs: the translational rule
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
      // taper). The groups are split where consecutive mean separations differ by more than
      // the pair tolerance (the same 5 % as the local constancy).
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
          if (piece.constant && piece.interacting && piece.sample_count > 0)
          {
            grouped.push_back(
                {&piece, side, piece.weighted_separation / piece.sample_count});
          }
        }
      }
      std::sort(grouped.begin(), grouped.end(),
                [](const GroupedPiece &u, const GroupedPiece &v)
                { return u.mean_separation < v.mean_separation; });
      std::vector<std::pair<std::size_t, std::size_t>> groups;  // [begin, end) in grouped
      for (std::size_t i = 0; i < grouped.size();)
      {
        std::size_t j = i + 1;
        while (j < grouped.size() &&
               !quantizer.Less(kPairSeparationTolerance * grouped[j - 1].mean_separation,
                               grouped[j].mean_separation - grouped[j - 1].mean_separation))
        {
          j++;
        }
        groups.emplace_back(i, j);
        i = j;
      }
      for (const auto &[g0, g1] : groups)
      {
        // Lead piece (on A when the group has one, else on B) and the lateral A -> B there:
        // the frame of a two-edge feature (side 0 = chain A); a chain facing itself takes
        // its foot outside the self-pair neighbourhood.
        const GroupedPiece *lead = nullptr;
        for (std::size_t i = g0; i < g1; i++)
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
            self ? std::optional<double>(
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
        double weighted = 0.0;
        int samples = 0;
        for (std::size_t i = g0; i < g1; i++)
        {
          const SubPiece &piece = *grouped[i].piece;
          link.pieces[static_cast<std::size_t>(grouped[i].side)].push_back(
              {piece.run, piece.interval, piece.curved, piece.max_kappa,
               piece.weighted_separation, piece.sample_count});
          weighted += piece.weighted_separation;
          samples += piece.sample_count;
        }
        // One separation per link (and per curvature class of its feature, taken again
        // from the class's own samples when emitted): the mean chord reading.
        link.separation = weighted / samples;
        self_pairs += self ? 1 : 0;
        pair_links.push_back(std::move(link));
      }
    }
  }
  stage.End(std::to_string(chains.size()) + " chains, " +
            std::to_string(chain_pairs_examined) + " chain pairs within reach, " +
            std::to_string(chain_pairs_paired) + " with facing pieces, " +
            std::to_string(pair_links.size()) + " pair links (" +
            std::to_string(self_pairs) + " of chains facing themselves)");
}

// A single link is the two-edge feature of the bent-pair rule (side 0 = chain A, the
// lateral from A to B at the lead piece), split by curvature class.
void Identifier::EmitPairLink(const PairLink &link)
{
  const Run &ra = runs[link.run_a];
  const Run &rb = runs[link.run_b];
  const int gap_a = Dot(ra.gap_direction, link.lateral_ab) > 0.0 ? 1 : -1;
  const int gap_b = Dot(rb.gap_direction, link.lateral_ab) > 0.0 ? 1 : -1;
  const bool same = ra.conductor == rb.conductor;
  std::string base;
  std::optional<std::string> reason;
  if (gap_a > 0 && gap_b < 0)
  {
    base = same ? "SameConductorGap" : "DifferentConductorGap";
  }
  else if (gap_a < 0 && gap_b > 0)
  {
    base = "SameConductorStrip";
  }
  else
  {
    base = "UnclassifiedParallelPair";
    reason = "parallel metal edges within 2R with the gap on the same side (overlapping "
             "metal in one process plane)";
  }
  for (const bool curved_class : {false, true})
  {
    double length = 0.0, max_kappa = 0.0, weighted = 0.0;
    int samples = 0;
    for (int side = 0; side < 2; side++)
    {
      for (const auto &piece : link.pieces[static_cast<std::size_t>(side)])
      {
        if (piece.curved != curved_class)
        {
          continue;
        }
        length += piece.interval.second - piece.interval.first;
        max_kappa = std::max(max_kappa, piece.max_kappa);
        weighted += piece.weighted_separation;
        samples += piece.sample_count;
      }
    }
    if (length <= kSignatureLengthQuantumOverRadius * R)
    {
      continue;  // slivers between cuts, below the signature grid
    }
    // One separation per pair and curvature class: the mean chord reading of its samples
    // (exact for a constant pair; a slow taper is described by its mean).
    const double separation = weighted / samples;
    std::vector<TranslationalEdge> edges = {
        {0.0, gap_a, ra.conductor, InterfaceNames(ra.targets), ra.boundary_law},
        {separation, gap_b, rb.conductor, InterfaceNames(rb.targets), rb.boundary_law}};
    auto translational = CanonicalTranslationalSignature(edges, R);
    nlohmann::json signature = std::move(translational.signature);
    std::string type = base;
    if (reason)
    {
      signature["Reason"] = *reason;
    }
    if (curved_class)
    {
      type = "Curved" + base;
      MFEM_VERIFY(max_kappa > 0.0, "A curved pair without curvature!");
      signature["RadiusOverR"] =
          RoundTo(1.0 / (max_kappa * R), kSignatureLengthQuantumOverRadius);
    }
    const int feature = NewFeature(type, signature, translational.chirality);
    features[feature].origin = link.lead_point;
    features[feature].axes = {ra.tangent, link.lateral_ab, n_ref};
    feature_sides[feature] = 2;
    for (int side = 0; side < 2; side++)
    {
      for (const auto &piece : link.pieces[static_cast<std::size_t>(side)])
      {
        if (piece.curved == curved_class)
        {
          claims[piece.run].push_back({feature, 2, piece.interval, side});
        }
      }
    }
  }
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

// Parallel classes -> elementary longitudinal intervals -> connected active runs.
void Identifier::BuildTranslationalFeatures()
{
  struct ClassMember
  {
    std::size_t run;
    double u0, u1, w;
    int gap_sign;
  };
  // Direction classes per metal plane: parallel runs of two planes never pair (decision
  // 82(1)); metal of another plane within 2R is the CrossLayer exclusion.
  std::map<std::pair<int, std::array<long long int, 3>>, std::vector<std::size_t>> classes;
  for (std::size_t r = 0; r < runs.size(); r++)
  {
    if (runs[r].excluded || !chains[chain_index.at(runs[r].chain)].Rigid())
    {
      continue;  // chains with joints pair through the curved-edge chain rule
    }
    classes[std::make_pair(run_plane[r],
                           DirectionKey(SignCanonical(runs[r].tangent), 1.0e-9))]
        .push_back(r);
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

// A span of rigid parallel runs alone is the translational feature of design (b) 2: a pair
// or a ParallelEdgeCluster with the members' lateral offsets, side k = the k-th edge in
// increasing lateral offset.
void Identifier::EmitTranslationalSpan(const TranslationalSpan &span)
{
  std::vector<TranslationalEdge> edges;
  for (const auto &member : span.members)
  {
    edges.push_back({member.w, member.gap_sign, runs[member.run].conductor,
                     InterfaceNames(runs[member.run].targets),
                     runs[member.run].boundary_law});
  }
  auto translational = CanonicalTranslationalSignature(edges, R);
  nlohmann::json best = std::move(translational.signature);
  const int chirality = translational.chirality;
  std::string type;
  std::optional<std::pair<std::string, std::string>> exclusion;
  if (span.members.size() == 2)
  {
    const auto &lower = span.members[0], &upper = span.members[1];
    const bool same = runs[lower.run].conductor == runs[upper.run].conductor;
    if (lower.gap_sign > 0 && upper.gap_sign < 0)
    {
      type = same ? "SameConductorGap" : "DifferentConductorGap";
    }
    else if (lower.gap_sign < 0 && upper.gap_sign > 0)
    {
      type = "SameConductorStrip";
    }
    else
    {
      exclusion = std::make_pair(
          "IncompatibleParallelPair",
          "parallel metal edges within 2R with the gap on the same side (overlapping "
          "metal in one process plane)");
    }
  }
  else
  {
    type = "ParallelEdgeCluster";
  }
  // A same-side pair cannot occur in one process plane; it is described (and never
  // matched) as an UnclassifiedParallelPair feature so that the partition stays exact.
  if (exclusion)
  {
    type = "UnclassifiedParallelPair";
    best["Reason"] = exclusion->second;
  }
  const int feature = NewFeature(type, best, chirality);
  features[feature].origin = Scale(0.5 * (span.lo + span.hi), span.axis);
  features[feature].axes = {span.axis, span.lateral, n_ref};
  feature_sides[feature] = static_cast<int>(span.members.size());
  for (std::size_t k = 0; k < span.members.size(); k++)
  {
    const auto &member = span.members[k];
    if (member.interval.second - member.interval.first <= Tol())
    {
      continue;
    }
    // Side k = the k-th edge in increasing lateral offset (the signature's order for
    // chirality +1, reversed for -1).
    claims[member.run].push_back({feature, 2, member.interval, static_cast<int>(k)});
  }
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
  for (std::size_t r = 0; r < runs.size(); r++)
  {
    for (const auto &claim : claims[r])
    {
      if (claim.priority == 0)
      {
        Take(r, claim.interval);
      }
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
    std::sort(list.begin(), list.end(),
              [](const RunInterval &a, const RunInterval &b)
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
    auto [it, inserted] = components.emplace(root, std::make_pair(std::vector<std::size_t>{},
                                                                  std::vector<std::size_t>{}));
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
    if (links.size() == 1 && span_items.empty() && !pair_links[links.front()].Self())
    {
      EmitPairLink(pair_links[links.front()]);
      pairs++;
    }
    else if (links.empty() && span_items.size() == 1)
    {
      EmitTranslationalSpan(translational_spans[span_items.front()]);
      spans++;
    }
    else
    {
      AssembleStack(links, span_items);
      stacks++;
    }
  }
  stage.End(std::to_string(component_order.size()) + " cross-section components: " +
            std::to_string(pairs) + " single links, " + std::to_string(spans) +
            " single spans, " + std::to_string(stacks) +
            " assembled (stacks, chains facing themselves, mixed), " +
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
  const double reach = kInteractionDistanceOverRadius * R * (1.0 + kPairSeparationTolerance);
  const double neighbourhood = kSelfPairNeighbourhoodOverRadius * R;
  struct ELink
  {
    int chain_a, chain_b;
    double separation;
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
      ELink e{runs[lower.run].chain, runs[upper.run].chain, upper.w - lower.w, {}};
      e.pieces[0].push_back({lower.run, lower.interval, false, 0.0, 0.0, 0});
      e.pieces[1].push_back({upper.run, upper.interval, false, 0.0, 0.0, 0});
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
  // Links on every chain: (link, side), in link order. Breakpoints per chain (sorted, no
  // two within the decision quantum): the piece ends of every link on the chain and the
  // ends of the higher-priority claims (cluster portions, vertex windows) on it — where a
  // member of the cross-section is taken by a cluster or a vertex window the stack ends and
  // the remaining members are recomposed (the stack-end rule of decision 82(2)).
  std::map<int, std::vector<std::pair<std::size_t, int>>> links_on_chain;
  std::map<int, std::set<double>> breakpoints;
  auto AddBreakpoint = [&](int chain, double x)
  {
    auto &list = breakpoints[chain];
    auto it = list.lower_bound(x - Tol());
    if (it != list.end() && std::abs(*it - x) <= Tol())
    {
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
      for (const auto &piece : elinks[e].pieces[static_cast<std::size_t>(side)])
      {
        const auto interval = ChainInterval(piece);
        AddBreakpoint(chain, interval.first);
        AddBreakpoint(chain, interval.second);
      }
    }
  }
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
      AddBreakpoint(chain, interval.first);
      AddBreakpoint(chain, interval.second);
    }
  }
  auto Contains = [&](const LinkPiece &piece, double x)
  {
    const auto interval = ChainInterval(piece);
    return x >= interval.first - Tol() && x <= interval.second + Tol();
  };
  auto StrictlyInside = [&](const LinkPiece &piece, double x)
  {
    const auto interval = ChainInterval(piece);
    return x > interval.first + Tol() && x < interval.second - Tol();
  };
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
    return iv != list.begin() && x < (iv - 1)->second - Tol() && x > (iv - 1)->first + Tol();
  };
  // Images of the breakpoints through the links, breadth first over the chains (bounded by
  // the number of links: an image travels at most once over every link).
  {
    std::map<int, std::vector<double>> frontier;
    for (const auto &[chain, xs] : breakpoints)
    {
      frontier[chain].assign(xs.begin(), xs.end());
    }
    for (std::size_t depth = 0; depth <= elinks.size() && !frontier.empty(); depth++)
    {
      std::map<int, std::vector<double>> next;
      for (const auto &[chain, xs] : frontier)
      {
        const auto it = links_on_chain.find(chain);
        if (it == links_on_chain.end())
        {
          continue;
        }
        for (const double x : xs)
        {
          for (const auto &[e, side] : it->second)
          {
            const ELink &link = elinks[e];
            bool active = false;
            for (const auto &piece : link.pieces[static_cast<std::size_t>(side)])
            {
              active = active || StrictlyInside(piece, x);
            }
            if (!active)
            {
              continue;
            }
            const int other = side == 0 ? link.chain_b : link.chain_a;
            const Point3D p = ChainAt(ChainOf(chain), x);
            const auto foot = ClosestPointOnChain(
                ChainOf(other), p, link.Self() ? std::optional<double>(x) : std::nullopt,
                reach);
            if (!std::isfinite(foot.distance) || !quantizer.Less(foot.distance, reach))
            {
              continue;
            }
            if (AddBreakpoint(other, foot.x))
            {
              next[other].push_back(foot.x);
            }
          }
        }
      }
      frontier = std::move(next);
    }
  }
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
    bool geometric_fallback = false;
  };
  std::size_t fallbacks = 0;
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
    for (std::size_t i = 0; i < composition.nodes.size() && composition.nodes.size() < 64; i++)
    {
      const Node node = composition.nodes[i];
      const auto it = links_on_chain.find(node.chain);
      if (it == links_on_chain.end())
      {
        continue;
      }
      for (const auto &[e, side] : it->second)
      {
        const ELink &link = elinks[e];
        bool active = false;
        for (const auto &piece : link.pieces[static_cast<std::size_t>(side)])
        {
          active = active || Contains(piece, node.x);
        }
        if (!active)
        {
          continue;
        }
        const int other = side == 0 ? link.chain_b : link.chain_a;
        const auto foot = ClosestPointOnChain(
            ChainOf(other), node.point,
            link.Self() ? std::optional<double>(node.x) : std::nullopt, reach);
        if (!std::isfinite(foot.distance) || !quantizer.Less(foot.distance, reach) ||
            Taken(other, foot.x))
        {
          continue;  // no partner there, or the partner is cluster / window metal
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
    // Provisional orientation (final: the canonical signature's, below): the member with
    // the smallest (chain, position along its chain) on side 0 (the two-edge convention:
    // side 0 = chain A, the lower chain index), which decides only for a cross-section that
    // is its own mirror image.
    auto Key = [](const Node &n) { return std::make_pair(n.chain, n.x); };
    if (Key(composition.nodes.back()) < Key(composition.nodes.front()))
    {
      std::reverse(composition.nodes.begin(), composition.nodes.end());
      composition.lateral = Scale(-1.0, composition.lateral);
      for (auto &n : composition.nodes)
      {
        n.position = -n.position;
      }
    }
    // Offsets from the consecutive links' separations.
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
          for (const auto &piece : link.pieces[static_cast<std::size_t>(side)])
          {
            if (Contains(piece, lower.x))
            {
              separation = link.separation;
              break;
            }
          }
        }
      }
      if (!separation)
      {
        separation = upper.position - lower.position;  // geometric (no consecutive link)
        composition.geometric_fallback = true;
      }
      composition.offsets[i + 1] = composition.offsets[i] + *separation;
    }
    return composition;
  };
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
      if (x1 - x0 <= Tol())
      {
        continue;
      }
      const double mid = 0.5 * (x0 + x1);
      bool covered = false;
      for (const auto &[e, side] : it->second)
      {
        for (const auto &piece : elinks[e].pieces[static_cast<std::size_t>(side)])
        {
          covered = covered || Contains(piece, mid);
        }
      }
      if (!covered)
      {
        continue;
      }
      auto composition = Compose(chain, mid);
      const std::size_t k = composition.nodes.size();
      if (k < 2)
      {
        continue;
      }
      fallbacks += composition.geometric_fallback ? 1 : 0;
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
      // the lateral axis is oriented so that the signature's first edge is side 0 (chirality
      // +1 always; the same sides at every cross-section of the feature whatever the chain
      // numbering), a symmetric one keeps the provisional orientation (chirality 0).
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
      bool curved = false;
      double max_kappa = 0.0;
      std::vector<int> member_chains;
      int own_side = -1;
      for (std::size_t n = 0; n < k; n++)
      {
        const Node &node = composition.nodes[n];
        curved = curved || IsCurvedAt(ChainOf(node.chain), node.x);
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
      on_chain.push_back({chain, x0, x1, key, own_side, curved, max_kappa, type,
                          translational.signature, translational.chirality, reason,
                          first.point,
                          {runs[first.run].tangent, composition.lateral, n_ref},
                          static_cast<int>(k)});
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
      dbg << "    link " << e.chain_a << " - " << e.chain_b << " sep " << e.separation
          << " pieces " << e.pieces[0].size() << " / " << e.pieces[1].size() << "\n";
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
          dbg << "      side " << side << " run " << piece.run << " x [" << iv.first << ", " << iv.second << "] curved " << piece.curved << "\n";
        }
    }
    for (const auto &a : assigned)
    {
      dbg << "    assigned chain " << a.chain << " [" << a.x0 << ", " << a.x1 << "] side " << a.side
          << " k " << a.sides << " curved " << a.curved << " " << a.type << "\n";
    }
    input.log(dbg.str());
  }
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
    }
    const int feature = NewFeature(type, signature, lead.chirality);
    features[feature].origin = lead.origin;
    features[feature].axes = lead.axes;
    feature_sides[feature] = lead.sides;
    for (const std::size_t i : members)
    {
      const Assigned &entry = assigned[i];
      const Chain &C = ChainOf(entry.chain);
      // Runs overlapping [x0, x1] (run_offset ascending).
      std::size_t kr = std::min(
          static_cast<std::size_t>(
              std::upper_bound(C.run_offset.begin(), C.run_offset.end(), entry.x0) -
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
  if (fallbacks > 0 && input.log)
  {
    input.log("  Identification stacks: " + std::to_string(fallbacks) +
              " elementary intervals with a geometric offset (no consecutive link)\n");
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
  auto CoresOnRun = [&](std::size_t a, std::size_t b,
                        const std::vector<std::vector<Interval>> &zones_a,
                        const std::vector<std::vector<Interval>> &zones_b)
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
      for (const auto &piece : complement)
      {
        const Point3D q0 = rb.At(piece.first), q1 = rb.At(piece.second);
        auto f = [&](double s) { return PointSegmentDistance(ra.At(s), q0, q1); };
        if (auto interval = ConvexSublevelInterval(f, lo, hi, interaction, quantizer))
        {
          found.push_back(*interval);
        }
      }
    }
    for (const auto &interval : MergeIntervals(found, Tol()))
    {
      cores.push_back(
          {a, interval, ra.At(interval.first), ra.At(interval.second), run_plane[a]});
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
    for (const std::size_t b : RunsNearRun(a, interaction + 2.0 * Tol()))
    {
      if (b <= a || runs[b].excluded || runs[a].chain == runs[b].chain ||
          run_plane[a] != run_plane[b])
      {
        continue;  // no events between two metal planes (decision 82(1))
      }
      run_pairs_examined++;
      const Chain &ca = chains[chain_index.at(runs[a].chain)];
      const Chain &cb = chains[chain_index.at(runs[b].chain)];
      if (ca.Rigid() && cb.Rigid() &&
          !DirectionLess(std::abs(Dot(runs[a].tangent, runs[b].tangent)),
                         1.0 - kParallelCosineTolerance))
      {
        continue;  // parallel straight runs: a translational interaction
      }
      if (!quantizer.Less(SegmentSegmentDistance(runs[a].start, runs[a].end, runs[b].start,
                                                 runs[b].end),
                          interaction))
      {
        continue;
      }
      std::vector<std::vector<Interval>> zones_a = ThroughZones(a, ca, cb),
                                         zones_b = ThroughZones(b, cb, ca);
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
      CoresOnRun(a, b, zones_a, zones_b);
      CoresOnRun(b, a, zones_b, zones_a);
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
                         SitePlane(sites[i])});
        cores.push_back({std::numeric_limits<std::size_t>::max(),
                         {0.0, 0.0},
                         sites[j].point,
                         sites[j].point,
                         SitePlane(sites[j])});
      }
    }
  }

  // Connected union of the radius-R balls: cores whose distance is below 2R (candidates
  // from the cores' boxes, united in the former (i, j) order).
  Point3D cores_lo{}, cores_hi{};
  for (std::size_t i = 0; i < cores.size(); i++)
  {
    Point3D lo, hi;
    BoundingBox(cores[i].p0, cores[i].p1, lo, hi);
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
    BoundingBox(cores[i].p0, cores[i].p1, lo, hi);
    core_grid.Insert(i, lo, hi);
  }
  UnionFind uf(cores.size());
  for (std::size_t i = 0; i < cores.size(); i++)
  {
    Point3D lo, hi;
    BoundingBox(cores[i].p0, cores[i].p1, lo, hi);
    for (const std::size_t j : core_grid.Query(lo, hi, interaction + 2.0 * Tol()))
    {
      if (j <= i || cores[i].plane != cores[j].plane)
      {
        continue;
      }
      if (quantizer.Less(
              SegmentSegmentDistance(cores[i].p0, cores[i].p1, cores[j].p0, cores[j].p1),
              interaction))
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
          quantizer.Less(PointSegmentDistance(sites[s].point, cores[i].p0, cores[i].p1),
                         join))
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
  // run).
  std::size_t largest_cluster_edges = 0;
  double signature_seconds = 0.0;
  for (std::size_t c = 0; c < cluster_cores.size(); c++)
  {
    stage.Progress(c, cluster_cores.size());
    std::vector<SignaturePortion> portions;
    std::vector<SignatureVertex> vertices;
    std::vector<std::pair<std::size_t, Interval>> claimed;
    std::vector<std::size_t> candidate_runs;
    for (const auto &core : cluster_cores[c])
    {
      Point3D lo, hi;
      BoundingBox(core.p0, core.p1, lo, hi);
      const auto near = RunsNear(lo, hi, ball + 2.0 * Tol());
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
      for (const auto &core : cluster_cores[c])
      {
        if (!quantizer.Less(
                SegmentSegmentDistance(runs[r].start, runs[r].end, core.p0, core.p1), ball))
        {
          continue;
        }
        const auto within = RunIntervalWithin(r, core.p0, core.p1, ball);
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
        claimed.emplace_back(r, interval);
        portions.push_back({runs[r].At(interval.first), runs[r].At(interval.second),
                            runs[r].gap_direction, runs[r].conductor,
                            InterfaceNames(runs[r].targets), runs[r].boundary_law});
      }
    }
    for (const std::size_t s : cluster_sites[c])
    {
      vertices.push_back({sites[s].point, sites[s].type, sites[s].turn_degrees});
    }
    MFEM_VERIFY(!portions.empty(), "A spatial cluster claims no perimeter!");
    largest_cluster_edges = std::max(largest_cluster_edges, portions.size());
    const auto signature_started = std::chrono::steady_clock::now();
    const auto canonical = CanonicalClusterSignature(
        portions, vertices, n_ref, R,
        [&](std::size_t done, std::size_t total)
        {
          stage.Progress(done, total,
                         "candidate frames of cluster " + std::to_string(c) + " / " +
                             std::to_string(cluster_cores.size()) + ", " +
                             std::to_string(portions.size()) + " edges");
        });
    signature_seconds +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - signature_started)
            .count();
    nlohmann::json signature = canonical.signature;
    signature["EdgeCount"] = portions.size();
    const int feature = NewFeature("SpatialEdgeCluster", signature, canonical.chirality);
    features[feature].origin = canonical.origin;
    features[feature].axes = canonical.axes;
    for (const auto &[r, interval] : claimed)
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
  {
    std::ostringstream counts;
    counts << run_pairs_examined << " run pairs within reach, " << run_cores
           << " event cores on runs + " << (cores.size() - run_cores) << " site cores, "
           << cluster_cores.size() << " clusters (largest " << largest_cluster_edges
           << " edges), canonical signatures " << std::fixed << std::setprecision(2)
           << signature_seconds << " s, " << features.size() << " features so far";
    stage.End(counts.str());
  }
}

void Identifier::Assign(IdentificationResult &result)
{
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
    // x = the first arm away from the site, y in the process plane, z = n_ref. Corner: the
    // arms ordered so that the second is counterclockwise from the first (a corner is its
    // own mirror image); endpoint: y toward the gap; junction: x = the canonical first arm,
    // y = +-(n x x) so that the canonical arm order proceeds counterclockwise in (x, y).
    std::array<Point3D, 3> axes = {Point3D{}, Point3D{}, n_ref};
    if (site.type == "ConvexCorner" || site.type == "ConcaveCorner")
    {
      signature = CanonicalCornerSignature(site.interfaces, site.boundary_law,
                                           site.angle_degrees, site.corner_radius / R);
      Point3D first = site.arm_directions[0], second = site.arm_directions[1];
      if (Dot(Cross(n_ref, first), second) < 0.0)
      {
        std::swap(first, second);
      }
      axes[0] = Normalize(Sub(first, Scale(Dot(first, n_ref), n_ref)));
      axes[1] = Normalize(Cross(n_ref, axes[0]));
    }
    else if (site.type == "Junction")
    {
      JunctionCanonicalOrder order;
      signature = CanonicalJunctionSignature(site.interfaces, site.boundary_law,
                                             site.arm_angles, site.arm_conductors, &order);
      const Point3D first = site.arm_directions[order.first_arm];
      axes[0] = Normalize(Sub(first, Scale(Dot(first, n_ref), n_ref)));
      axes[1] = Scale(order.reversed ? -1.0 : 1.0, Normalize(Cross(n_ref, axes[0])));
    }
    else if (!site.arm_directions.empty())
    {
      const Point3D first = site.arm_directions.front();
      axes[0] = Normalize(Sub(first, Scale(Dot(first, n_ref), n_ref)));
      axes[1] = Normalize(Cross(n_ref, axes[0]));
      if (Dot(axes[1], site.gap_direction) < 0.0)
      {
        axes[1] = Scale(-1.0, axes[1]);
      }
    }
    const int feature = NewFeature(site.type, signature);
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
    std::sort(run_claims.begin(), run_claims.end(),
              [](const Claim &a, const Claim &b)
              {
                return std::tie(a.priority, a.feature, a.interval) <
                       std::tie(b.priority, b.feature, b.interval);
              });
    // Two claims of one priority by different features overlapping on the run would be
    // decided by the feature id below: the rules (clusters disjoint, windows shortened to
    // abut, pairs / stacks assembled per cross-section) never produce one; counted.
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
    // Pieces below the signature grid are roundoff between the boundaries of two claims (a
    // cluster ball cutting a pair piece, a claim ending next to a run end): each joins the
    // adjacent portion on the run instead of standing as a feature the signature cannot
    // resolve (DS-SCT-001: three CurvedSameConductorStrip features of 1.6e-7 um in total).
    {
      const double sliver = kSignatureLengthQuantumOverRadius * R;
      for (const auto &piece :
           SubtractIntervals({Interval{0.0, runs[r].length}}, taken, Tol()))
      {
        if (piece.second - piece.first <= sliver)
        {
          assigned[r].emplace_back(piece.first, piece.second, -1, 0);  // unclaimed sliver
        }
      }
      std::sort(assigned[r].begin(), assigned[r].end());
      auto &pieces = assigned[r];
      for (std::size_t i = 0; i < pieces.size();)
      {
        auto &[lo, hi, feature, side] = pieces[i];
        (void)side;
        if (hi - lo > sliver && feature >= 0)
        {
          i++;
          continue;
        }
        auto joins = [&](std::size_t j)
        {
          const auto &[jlo, jhi, jfeature, jside] = pieces[j];
          (void)jside;
          return jfeature >= 0 &&
                 (std::abs(jhi - lo) <= Tol() || std::abs(hi - jlo) <= Tol());
        };
        if (i > 0 && joins(i - 1))
        {
          std::get<1>(pieces[i - 1]) = hi;
        }
        else if (i + 1 < pieces.size() && joins(i + 1))
        {
          std::get<0>(pieces[i + 1]) = lo;
        }
        else if (feature >= 0)
        {
          i++;  // a pair piece alone in the remainder keeps its feature
          continue;
        }
        pieces.erase(pieces.begin() + i);
      }
    }
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
    // facing test is the pair tolerance on the separation (a slow taper's wider end is
    // still facing). Recorded as claim resolution: both sides of a pair are then mutual.
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
      const double reach = features[feature].signature["SeparationOverR"].get<double>() *
                           R * (1.0 + kPairSeparationTolerance);
      // The other side's pieces by run: a piece faces only the pieces on runs whose boxes
      // lie within the reach of its own run (the facing intervals of every other piece are
      // empty; the union is sorted, so the candidate order does not matter).
      std::array<std::map<std::size_t, std::vector<const SidePiece *>>, 2> pieces_by_run;
      for (int k = 0; k < 2; k++)
      {
        for (const auto &piece : sides[k])
        {
          pieces_by_run[k][piece.run].push_back(&piece);
        }
      }
      for (int k = 0; k < 2; k++)
      {
        for (const auto &piece : sides[k])
        {
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
        if (it == lengths.end() || it->second <= kSignatureLengthQuantumOverRadius * R)
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
          features[feature].axes = {runs[r].tangent, Cross(n_ref, runs[r].tangent), n_ref};
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
          features[feature].axes = {runs[r].tangent, Cross(n_ref, runs[r].tangent), n_ref};
          it = isolated_features.emplace(key, feature).first;
        }
        for (const auto &piece : straight)
        {
          assigned[r].emplace_back(piece.first, piece.second, it->second, 0);
        }
      }
    }
    std::sort(assigned[r].begin(), assigned[r].end());
    // Curvature annotation of every feature over its assigned portions.
    {
      const Chain &chain = chains[chain_index.at(runs[r].chain)];
      const double offset = chain.run_offset[RunIndexInChain(chain, r)];
      for (const auto &[lo, hi, feature, side] : assigned[r])
      {
        (void)side;
        feature_max_kappa[feature] = std::max(
            feature_max_kappa[feature], MaxCurvature(chain, offset + lo, offset + hi));
      }
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
      for (const auto &[lo, hi, feature, side] : assigned[r])
      {
        const auto portion = SegmentPortion(lo, hi);
        if (!portion)
        {
          continue;
        }
        const auto [s0, s1] = *portion;
        table.portions.push_back({s0, s1, static_cast<double>(feature)});
        features[feature].portions.push_back({rs.segment, s0, s1, side});
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
      continue;
    }
    IdentifiedVertex entry;
    entry.vertex = v;
    entry.point_contact = IsPointContact(v);
    if (!IsFeatureVertex(v))
    {
      entry.type = cross_layer_vertices.find(v) != cross_layer_vertices.end()
                       ? "Excluded"
                       : (vertex.on_truncation_boundary ? "TruncationCut" : "PortCut");
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
    std::size_t corners = 0, bends = 0, absorbed = 0;
    for (const auto &arc : arcs)
    {
      (corners += arc.corner ? 1 : 0), (bends += arc.corner ? 0 : 1);
    }
    for (const int a : vertex_arc)
    {
      absorbed += a >= 0 ? 1 : 0;
    }
    stage.End(std::to_string(corners) + " rounded corners, " + std::to_string(bends) +
              " bends of exact radius, " + std::to_string(absorbed) + " joints absorbed");
  }
  stage.Begin("runs / chains");
  BuildRuns(input, runs, chains, chain_index);
  BuildRunIndex();
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
    stage.Begin("pairs / stacks");
    BuildPairsAndStacks();
  }
  stage.Begin("claim resolution / assignment / tables");
  Assign(result);
  stage.End(std::to_string(result.features.size()) + " features, " +
            std::to_string(result.exclusions.size()) + " exclusion classes");
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
  w.Pod(result.assigned_length);
  w.Pod(result.excluded_length);
  w.String(result.geometry_digest);
  w.Size(result.same_priority_claim_overlaps);
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
    w.Pod(f.matched_model.has_value());
    w.String(f.matched_model.value_or(""));
  }
  w.Size(result.segments.size());
  for (const auto &s : result.segments)
  {
    w.Point(s.key[0]);
    w.Point(s.key[1]);
    w.Pod(s.length);
    w.Pod(s.chain);
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
  return std::move(w.buffer);
}

IdentificationResult DeserializeIdentificationResult(const std::string &buffer)
{
  ByteReader r(buffer);
  IdentificationResult result;
  result.radius = r.Pod<double>();
  result.reference_process_normal = r.Point();
  result.perimeter_length = r.Pod<double>();
  result.assigned_length = r.Pod<double>();
  result.excluded_length = r.Pod<double>();
  result.geometry_digest = r.String();
  result.same_priority_claim_overlaps = r.Size();
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
    const bool has_model = r.Pod<bool>();
    const std::string model = r.String();
    if (has_model)
    {
      f.matched_model = model;
    }
  }
  result.segments.resize(r.Size());
  for (auto &s : result.segments)
  {
    s.key[0] = r.Point();
    s.key[1] = r.Point();
    s.length = r.Pod<double>();
    s.chain = r.Pod<int>();
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
  MFEM_VERIFY(r.Done(), "Identification result buffer was not consumed exactly!");
  return result;
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
  nlohmann::json feature_list = nlohmann::json::array();
  for (const auto &feature : features)
  {
    nlohmann::json portions = nlohmann::json::array();
    nlohmann::json sides = nlohmann::json::array();
    bool multi_sided = false;
    for (const auto &portion : feature.portions)
    {
      portions.push_back({portion.segment, L(portion.s0), L(portion.s1)});
      sides.push_back(portion.side);
      multi_sided = multi_sided || portion.side != 0;
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
        {"Frame",
         {{"Origin", P(feature.origin)},
          {"Axes", {D(feature.axes[0]), D(feature.axes[1]), D(feature.axes[2])}}}},
        {"BendRadiusOverR",
         feature.bend_radius_over_R
             ? nlohmann::json(std::round(*feature.bend_radius_over_R /
                                         kSignatureLengthQuantumOverRadius) *
                              kSignatureLengthQuantumOverRadius)
             : nlohmann::json(nullptr)},
        {"Match", {{"Status", feature.matched_model ? "Matched" : "Missing"}}}};
    if (feature.matched_model)
    {
      entry["Match"]["Model"] = *feature.matched_model;
    }
    if (multi_sided)
    {
      // Side of every portion of a pair / parallel cluster (parallel to Portions): the
      // signature's edge order for chirality +1, reversed for -1.
      entry["Sides"] = sides;
    }
    feature_list.push_back(std::move(entry));
  }
  nlohmann::json segment_list = nlohmann::json::array();
  for (const auto &segment : segments)
  {
    nlohmann::json entry = {{"Key", {P(segment.key[0]), P(segment.key[1])}},
                            {"Length", L(segment.length)},
                            {"Chain", segment.chain}};
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
  return {
      {"Version", 2},
      {"MatchingRadius", scaled_radius},
      {"Conventions",
       {{"CornerTurnToleranceDegrees", kCornerTurnToleranceDegrees},
        {"InteractionDistanceOverR", kInteractionDistanceOverRadius},
        {"ThroughVertexZoneOverR", kThroughVertexZoneOverRadius},
        {"ClusterBallOverR", kClusterBallOverRadius},
        {"VertexJoinsClusterOverR", kVertexJoinsClusterOverRadius},
        {"VertexWindowOverR", kVertexWindowOverRadius},
        {"ParallelCosineTolerance", kParallelCosineTolerance},
        {"RoundedCornerTangentTolerance", kRoundedCornerTangentTolerance},
        {"ArcFitToleranceRelative", kArcFitToleranceRelative},
        {"ArcRule", "joints (>= 3) joined by pieces < 2R, same turn sign, <= 180 deg, one "
                    "circle with the arm tangents within the fit tolerance: radius < R and "
                    "turn > corner threshold = one rounded corner (own chain, arms meet "
                    "through it); radius >= R = a bend of exact radius inside its chain"},
        {"LengthQuantumOverR", kLengthQuantumOverRadius},
        {"DirectionQuantum", kDirectionQuantum},
        {"SignatureLengthQuantumOverR", kSignatureLengthQuantumOverRadius},
        {"SignatureAngleQuantumDegrees", kSignatureAngleQuantumDegrees},
        {"StraightBendRadiusOverR", kStraightBendRadiusOverRadius},
        {"CurvatureWindowOverR", kCurvatureWindowOverRadius},
        {"PairSeparationToleranceRelative", kPairSeparationTolerance},
        {"PairSeparationSamplesPerInterval", kPairSeparationSamplesPerInterval},
        {"PairSeparationEstimate",
         "per sample: chord reading C = min over the two chains of the maximum sampled "
         "closest-point distance within max(R, local chord) of the sample / its foot where "
         "the chain bends, R on straight runs; inscribed reading C / cos(turn / 2) with "
         "the "
         "larger local joint turn; interacting iff both < 2R; feature separation = mean C"},
        {"PairConstancyWindowOverR", 1.0},
        {"PairSampleSpacingOverR", 0.5},
        {"PairCandidateReachOverR",
         kInteractionDistanceOverRadius * (1.0 + kPairSeparationTolerance)},
        {"SamplingMargins", "the bent-pair candidate reach 2R (1 + 0.05) and the mutual-sides "
                            "facing test at the pair separation x 1.05 gather candidates / "
                            "samples only; every interaction decision is the 3D distance "
                            "strictly below 2R on the quantized grid"},
        {"SelfPairNeighbourhoodOverR", kSelfPairNeighbourhoodOverRadius},
        {"StackRule", "links (bent pairs) and translational spans sharing a run over a common "
                      "interval are one cross-section: k >= 3 edges with consecutive "
                      "separations < 2R are one ParallelEdgeCluster / "
                      "CurvedParallelEdgeCluster (offsets from the consecutive links, gap "
                      "pattern, conductors, bend radius), straight and along bends; the "
                      "pairwise candidates inside it are superseded; a chain folding back "
                      "within 2R beyond pi R of arc length pairs with itself"},
        {"CrossLayerReachOverR", kInteractionDistanceOverRadius},
        {"PlaneRule", "features (pairs, clusters, vertex joins, translational classes) "
                      "never span two metal planes; metal of another plane within the "
                      "interaction distance is the CrossLayer exclusion"},
        {"PortRule", "metal perimeter bordering a LumpedPort / WavePort boundary face is "
                     "the Port exclusion; its vertices are PortCut, never features"},
        {"Comparison", "strict less on the quantized grid"}}},
      {"ReferenceProcessNormal", D(reference_process_normal)},
      {"Features", feature_list},
      {"Segments", segment_list},
      {"Vertices", vertex_list},
      {"Exclusions", exclusion_list},
      {"Totals",
       {{"PerimeterLength", L(perimeter_length)},
        {"AssignedLength", L(assigned_length)},
        {"ExcludedLength", L(excluded_length)}}},
      {"Diagnostics",
       {{"SamePriorityClaimOverlaps", same_priority_claim_overlaps},
        {"Rule", "claims of one priority by different features never overlap on a run "
                 "(clusters disjoint, windows abut, pairs / stacks assembled per "
                 "cross-section): claim resolution never decides by feature id"}}},
      {"GeometryDigest", geometry_digest}};
}

}  // namespace palace
