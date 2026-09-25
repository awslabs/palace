// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfaceresponseidentification.hpp"

#include <algorithm>
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
constexpr double kVertexJoinsClusterOverRadius = 3.0;
constexpr double kVertexWindowOverRadius = 1.0;
constexpr double kParallelCosineTolerance = 1.0e-8;
constexpr double kRoundedCornerTangentTolerance = 0.05;

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
// candidate reach 2R (1 + kPairSeparationTolerance) of the other chain, not beyond its ends,
// outside the shared-vertex zones; samples at most R / 2 apart) is constant when the sampled
// distances within R of it along its own chain vary by at most kPairSeparationTolerance x
// their minimum (a polyline of sub-corner turns <= 30 deg at constant width varies by
// 1 / cos(15 deg) - 1 = 3.5 %; the pair response sensitivity d dR/dd is O(1)). The constant
// portions are the pair; the portions that are not (divergence at tees and port ends, fast
// tapers, acute corner arms) keep the event rule, so a slow taper is a pair and a tee is a
// cluster. Whether a constant portion interacts is decided on the separation of the
// underlying curves, so that a discretisation never changes the classification: the chords
// of a polyline inscribed in a curve lie inside it (two concentric inscribed polylines are
// w cos(turn / 2) apart mid-chord and exactly w apart at their vertices), while the exact
// offset polyline of a bent path keeps corresponding chords at the design separation and
// the outer side's samples near the joints project onto the inner vertices at up to
// w / cos(turn / 2). In both constructions the sampled closest-point distance from one
// chain to the other reaches the curve separation w as its maximum on the side whose
// maximum is smaller: the sample's chord reading C = min over the two chains of the maximum
// sampled distance within a window of half-width max(R, local chord) where the chain bends
// and R on straight runs, about the sample on its own chain and about its foot on the other
// chain (a straight taper is read locally). C is exact for an offset polyline; two polylines
// inscribed in the curves at aligned angles are C = w cos(turn / 2) apart everywhere (chords
// and vertex-to-polyline alike) for a curve separation w, and the polyline pair alone cannot
// tell the two constructions apart (they differ at order turn^2): the inscribed reading is
// C / cos(turn / 2) with turn = the larger local joint turn of the two chains. A portion
// interacts iff BOTH readings are below 2R on the quantized grid — the same strict-less
// decision as a straight parallel pair at that separation, taken on the non-interacting
// side of the recorded ambiguity w (1 / cos(turn / 2) - 1) (below 1e-4 w for joints under
// 1.6 deg; a CPW gap of exactly 2R along a bend is isolated edges like a straight one;
// DS-SCT-001's 4 um gaps at R = 2 um read 3.9998 mid-chord and became 3 mm clusters). The
// pair feature's separation is the mean chord reading over its samples; the cross-chord
// interactions of a locally constant portion are never event cores, whether or not it
// interacts.
constexpr double kStraightBendRadiusOverRadius = 10.0;
constexpr double kCurvatureWindowOverRadius = 1.0;
constexpr double kPairSeparationTolerance = 0.05;
constexpr int kPairSeparationSamplesPerInterval = 16;

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
  int next_chain_id =
      segments_by_chain.empty() ? 0 : segments_by_chain.rbegin()->first + 1;
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
      auto [entry, inserted] = id_by_root.try_emplace(
          root, id_by_root.empty() ? chain_id : next_chain_id);
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

nlohmann::json SerializeInFrame(const std::vector<SignaturePortion> &portions,
                                const std::vector<SignatureVertex> &vertices,
                                const Point3D &origin, const Point3D &x, const Point3D &y,
                                double radius)
{
  auto Local = [&](const Point3D &p)
  {
    const Point3D r = Sub(p, origin);
    return std::array<double, 2>{
        RoundTo(Dot(r, x) / radius, kSignatureLengthQuantumOverRadius),
        RoundTo(Dot(r, y) / radius, kSignatureLengthQuantumOverRadius)};
  };
  auto LocalDirection = [&](const Point3D &v)
  {
    return std::array<double, 2>{RoundTo(Dot(v, x), kSignatureLengthQuantumOverRadius),
                                 RoundTo(Dot(v, y), kSignatureLengthQuantumOverRadius)};
  };
  struct Entry
  {
    nlohmann::json geometry;
    int conductor;
  };
  std::vector<Entry> entries;
  for (const auto &portion : portions)
  {
    auto a = Local(portion.p0), b = Local(portion.p1);
    if (b < a)
    {
      std::swap(a, b);
    }
    entries.push_back({nlohmann::json{{"P", {a[0], a[1], b[0], b[1]}},
                                      {"Gap", LocalDirection(portion.gap_direction)},
                                      {"Interfaces", portion.interfaces},
                                      {"Law", portion.boundary_law}},
                       portion.conductor});
  }
  std::sort(entries.begin(), entries.end(), [](const Entry &a, const Entry &b)
            { return a.geometry.dump() < b.geometry.dump(); });
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
                                          std::vector<int> arm_conductors)
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
  std::pair<std::vector<double>, std::vector<int>> best{arm_angles_degrees,
                                                        Relabel(arm_conductors)};
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
      }
    }
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

CanonicalSignature CanonicalClusterSignature(const std::vector<SignaturePortion> &portions,
                                             const std::vector<SignatureVertex> &vertices,
                                             const Point3D &process_normal, double radius)
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
  std::string best_key;
  bool have = false;
  std::set<int> minimal_handedness;
  for (const auto &x : candidates)
  {
    for (const int handedness : {1, -1})
    {
      const Point3D y = Scale(static_cast<double>(handedness), Cross(n, x));
      auto serialized = SerializeInFrame(portions, vertices, origin, x, y, radius);
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
};

struct EventCore
{
  std::size_t run = 0;
  Interval interval;
  Point3D p0{}, p1{};
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
    : input(input_), quantizer(input_.radius), R(input_.radius)
  {
  }

  IdentificationResult Identify();

private:
  const IdentificationInput &input;
  Quantizer quantizer;
  double R;
  std::vector<Run> runs;
  std::vector<Chain> chains;
  std::map<int, std::size_t> chain_index;
  Point3D n_ref{};
  std::vector<std::optional<std::pair<std::string, std::string>>> segment_exclusion;
  std::vector<std::vector<Claim>> claims;  // per run
  // Per run: (other chain, interval) of the pieces of a constant-separation pair along a bend
  // with that chain, interacting (a pair feature) or not; their cross-chord interactions are
  // not event cores.
  std::vector<std::vector<std::pair<int, Interval>>> bent_claims;
  // Decision 73(3) CrossLayer zones: per run, the intervals within 2R of metal off the
  // run's plane (facing layers, walls, staples); vertices within 2R of such metal.
  std::vector<std::vector<Interval>> cross_layer;
  std::set<std::size_t> cross_layer_vertices;
  std::vector<IdentifiedFeature> features;
  std::vector<double> feature_max_kappa;  // per feature, over its assigned portions
  std::vector<VertexFeatureSite> sites;
  std::vector<std::vector<EventCore>> cluster_cores;
  std::vector<std::vector<std::size_t>> cluster_sites;
  std::map<std::size_t, int> vertex_feature;  // mesh vertex -> feature id

  double Tol() const { return quantizer.Quantum(); }

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
        vertex.on_truncation_boundary)
    {
      return false;  // simulation cut
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

  std::vector<std::size_t> RunsAtVertex(std::size_t v) const
  {
    std::vector<std::size_t> result;
    for (std::size_t r = 0; r < runs.size(); r++)
    {
      if (!runs[r].excluded && (runs[r].start_vertex == v || runs[r].end_vertex == v))
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

  void ClassifyPlanes();
  void ClassifyVertices();
  void DetectRoundedCorners();
  void ComputeCurvature();
  void BuildTranslationalFeatures();
  void BuildBentPairs();
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
  ChainPoint ClosestPointOnChain(const Chain &chain, const Point3D &p) const;
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
  if (input.faces.empty())
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
  std::vector<Face> faces;
  faces.reserve(input.faces.size());
  for (std::size_t f = 0; f < input.faces.size(); f++)
  {
    const auto &face = input.faces[f];
    if (face.vertices.size() < 3)
    {
      continue;
    }
    Face entry{f, !DirectionLess(std::abs(Dot(face.normal, n_ref)),
                                 1.0 - kParallelCosineTolerance),
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
    faces.push_back(entry);
  }
  const double reach = 2.0 * R;
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
    std::vector<Interval> zones;
    for (const auto &face : faces)
    {
      if (!OffPlane(face, offset) || !NearBox(face, lower, upper))
      {
        continue;
      }
      const auto &polygon = input.faces[face.index];
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
    for (const auto &face : faces)
    {
      if (!OffPlane(face, *offset) || !NearBox(face, vertex.coordinate, vertex.coordinate))
      {
        continue;
      }
      const auto &polygon = input.faces[face.index];
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
    }
    else if (incident.size() == 2)
    {
      const Point3D da = ArmDirection(incident[0], v);
      const Point3D db = ArmDirection(incident[1], v);
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
      std::vector<std::pair<double, int>> arms;  // angle, conductor
      for (const std::size_t r : incident)
      {
        const Point3D d = ArmDirection(r, v);
        double angle = std::atan2(Dot(d, y), Dot(d, x)) * 180.0 / std::acos(-1.0);
        if (angle < 0.0)
        {
          angle += 360.0;
        }
        arms.emplace_back(angle, runs[r].conductor);
      }
      std::sort(arms.begin(), arms.end());
      for (std::size_t i = 0; i < arms.size(); i++)
      {
        const double next = i + 1 < arms.size() ? arms[i + 1].first : arms[0].first + 360.0;
        site.arm_angles.push_back(next - arms[i].first);
        site.arm_conductors.push_back(arms[i].second);
      }
    }
    sites.push_back(std::move(site));
  }
}

// Fillet arcs: a maximal sequence of consecutive sub-threshold turns inside a chain bounded
// by two straight arms whose tangent distances from the virtual corner are < R, equal
// within 5 %, with a fillet radius in (0, R) (the legacy rounded-corner rule, now on runs
// so that the collinear midpoints inserted by refinement do not break the arc).
void Identifier::DetectRoundedCorners()
{
  for (auto &chain : chains)
  {
    const std::size_t m = chain.runs.size();
    chain.joint_excluded.assign(m, false);
    if (m < 3)
    {
      continue;
    }
    // Boundary k (between run k-1 and run k, cyclic for a closed chain) is a sub-threshold
    // turn; a run shorter than R with turns at both ends is an arc chord, longer runs are
    // arms.
    auto Turn = [&](std::size_t k)
    {
      if (k == 0 || k >= m)
      {
        if (!chain.closed)
        {
          return false;
        }
        k = 0;
      }
      const std::size_t a = chain.runs[(k + m - 1) % m], b = chain.runs[k];
      return !runs[a].excluded && !runs[b].excluded &&
             runs[a].end_vertex == runs[b].start_vertex &&
             DirectionLess(Dot(runs[a].tangent, runs[b].tangent), 1.0 - kDirectionQuantum);
    };
    std::vector<bool> is_arc(m, false);
    for (std::size_t k = 0; k < m; k++)
    {
      is_arc[k] = quantizer.Less(runs[chain.runs[k]].length, R) && Turn(k) && Turn(k + 1);
    }
    const auto first_arm = std::find(is_arc.begin(), is_arc.end(), false);
    if (first_arm == is_arc.end())
    {
      continue;  // a closed curve without straight arms (phase 3: curved class)
    }
    const std::size_t start = static_cast<std::size_t>(first_arm - is_arc.begin());
    std::size_t step = 0;
    while (step < m)
    {
      const std::size_t k = (start + step) % m;
      if (!is_arc[k])
      {
        step++;
        continue;
      }
      std::size_t count = 0;
      while (step + count < m && is_arc[(start + step + count) % m])
      {
        count++;
      }
      const std::size_t before = (k + m - 1) % m, after = (k + count) % m;
      const bool has_arms = (chain.closed || (k > 0 && k + count < m)) && !is_arc[before] &&
                            !is_arc[after] && Turn(k) &&
                            Turn((k + count) % m == 0 ? m : k + count);
      if (has_arms)
      {
        const Run &arm_a = runs[chain.runs[before]];
        const Run &arm_b = runs[chain.runs[after]];
        const Point3D ta = arm_a.tangent, tb = arm_b.tangent;
        const double cos_turn = std::clamp(Dot(ta, tb), -1.0, 1.0);
        const double turn = std::acos(cos_turn);
        if (std::sin(turn) > 1.0e-9)
        {
          // Virtual corner X = arm_a.end + a ta = arm_b.start - b tb, so that
          // w = arm_b.start - arm_a.end = a ta + b tb; solve in the (ta, y) plane.
          const Point3D w = Sub(arm_b.start, arm_a.end);
          const Point3D y = Cross(n_ref, ta);
          const double wx = Dot(w, ta), wy = Dot(w, y);
          const double tbx = Dot(tb, ta), tby = Dot(tb, y);
          if (std::abs(tby) > 1.0e-12)
          {
            const double b = wy / tby;
            const double a = wx - b * tbx;
            if (a > 0.0 && b > 0.0)
            {
              const double radius = 0.5 * (a + b) / std::tan(0.5 * turn);
              if (quantizer.Less(a, R) && quantizer.Less(b, R) &&
                  std::abs(a - b) <= kRoundedCornerTangentTolerance * std::max(a, b) &&
                  radius > 0.0 && quantizer.Less(radius, R))
              {
                VertexFeatureSite site;
                site.point = Add(arm_a.end, Scale(a, ta));
                site.turn_degrees = turn * 180.0 / std::acos(-1.0);
                site.angle_degrees = 180.0 - site.turn_degrees;
                site.corner_radius = radius;
                const Point3D bisector = Add(Scale(-1.0, ta), tb);
                const Point3D gap = Add(arm_a.gap_direction, arm_b.gap_direction);
                site.type = Dot(bisector, gap) > 0.0 ? "ConcaveCorner" : "ConvexCorner";
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
                // Claims: the arc runs completely, R along each arm from the tangent point.
                for (std::size_t c = 0; c < count; c++)
                {
                  const std::size_t r = chain.runs[(k + c) % m];
                  site.window.emplace_back(r, Interval{0.0, runs[r].length});
                  site.runs_at_site.push_back(r);
                }
                site.window.emplace_back(
                    chain.runs[before],
                    Interval{std::max(0.0, arm_a.length - R), arm_a.length});
                site.window.emplace_back(chain.runs[after],
                                         Interval{0.0, std::min(R, arm_b.length)});
                site.runs_at_site.push_back(chain.runs[before]);
                site.runs_at_site.push_back(chain.runs[after]);
                sites.push_back(std::move(site));
                // The rounded corner accounts for the turns of the fillet's joints (both
                // arm joints included): they take no part in the curved-edge chain rule.
                for (std::size_t c = 0; c <= count; c++)
                {
                  chain.joint_excluded[(k + c) % m] = true;
                }
              }
            }
          }
        }
      }
      step += count;
    }
  }
}

// ---------------------------------------------------------------------------------------
// Curved-edge chain rule: windowed curvature along every chain
// ---------------------------------------------------------------------------------------

std::size_t Identifier::RunIndexInChain(const Chain &chain, std::size_t run) const
{
  const auto position = std::find(chain.runs.begin(), chain.runs.end(), run);
  MFEM_VERIFY(position != chain.runs.end(), "Run missing from its chain!");
  return static_cast<std::size_t>(position - chain.runs.begin());
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
  const auto upper = std::upper_bound(nodes.begin(), nodes.end(), std::make_pair(x, 0.0),
                                      [](const auto &a, const auto &b)
                                      { return a.first < b.first; });
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
  for (const auto &[x, kappa] : chain.kappa_nodes)
  {
    if (x > x0 && x < x1)
    {
      best = std::max(best, kappa);
    }
  }
  return best;
}

bool Identifier::IsCurvedAt(const Chain &chain, double x, std::size_t *section) const
{
  for (std::size_t j = 0; j < chain.curved.size(); j++)
  {
    if (x >= chain.curved[j].first - Tol() && x <= chain.curved[j].second + Tol())
    {
      if (section)
      {
        *section = j;
      }
      return true;
    }
  }
  return false;
}

Identifier::ChainPoint Identifier::ClosestPointOnChain(const Chain &chain,
                                                       const Point3D &p) const
{
  ChainPoint best;
  best.distance = std::numeric_limits<double>::infinity();
  for (std::size_t k = 0; k < chain.runs.size(); k++)
  {
    const Run &r = runs[chain.runs[k]];
    if (r.excluded)
    {
      continue;
    }
    const double s = std::clamp(Dot(Sub(p, r.start), r.tangent), 0.0, r.length);
    const double distance = Distance(p, r.At(s));
    if (distance < best.distance)
    {
      best = {chain.runs[k], s, chain.run_offset[k] + s, distance};
    }
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
  for (auto &chain : chains)
  {
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
    for (std::size_t k = 0; k < m; k++)
    {
      if (k == 0 && !chain.closed)
      {
        continue;
      }
      const Run &a = runs[chain.runs[(k + m - 1) % m]];
      const Run &b = runs[chain.runs[k]];
      if (a.excluded || b.excluded || a.end_vertex != b.start_vertex ||
          chain.joint_excluded[k])
      {
        continue;
      }
      chain.joint_turn[k] = std::acos(std::clamp(Dot(a.tangent, b.tangent), -1.0, 1.0));
    }
    // Density breakpoints at the chord midpoints; cumulative turn F(x) = int_0^x density.
    struct Piece
    {
      double x0, x1, density;
    };
    std::vector<Piece> pieces;
    for (std::size_t k = 0; k < m; k++)
    {
      const double half = 0.5 * runs[chain.runs[k]].length;
      const double x0 = chain.run_offset[k];
      // First half of run k carries the joint before it, the second half the joint after.
      const std::size_t next = (k + 1) % m;
      const bool has_next = chain.closed || k + 1 < m;
      auto Density = [&](std::size_t joint)
      {
        const std::size_t prev = (joint + m - 1) % m;
        const double window =
            0.5 * (runs[chain.runs[prev]].length + runs[chain.runs[joint]].length);
        return window > 0.0 ? chain.joint_turn[joint] / window : 0.0;
      };
      pieces.push_back({x0, x0 + half, chain.joint_turn[k] > 0.0 ? Density(k) : 0.0});
      pieces.push_back({x0 + half, x0 + 2.0 * half,
                        has_next && chain.joint_turn[next] > 0.0 ? Density(next) : 0.0});
    }
    std::vector<double> cumulative(pieces.size() + 1, 0.0);
    for (std::size_t i = 0; i < pieces.size(); i++)
    {
      cumulative[i + 1] =
          cumulative[i] + pieces[i].density * (pieces[i].x1 - pieces[i].x0);
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
      double value = 0.0;
      for (std::size_t i = 0; i < pieces.size(); i++)
      {
        if (x <= pieces[i].x0)
        {
          break;
        }
        value = cumulative[i] + pieces[i].density * (std::min(x, pieces[i].x1) - pieces[i].x0);
      }
      return value + shift;
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
    xs.erase(std::unique(xs.begin(), xs.end(),
                         [&](double a, double b) { return std::abs(a - b) <= 1.0e-12 * R; }),
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
        chain.curved_max_kappa.push_back(MaxCurvature(chain, *open_start, x));
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
    // The larger joint turn (radians) at the ends of the sample's run and of the foot's run:
    // the chord reading C and the inscribed-vertex reading C / cos(turn / 2) of the pair
    // separation differ by this discretisation ambiguity.
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
  auto Pieces = [&](const Chain &A, const Chain &B, std::vector<Piece> &pieces)
  {
    std::vector<std::size_t> shared;
    std::set_intersection(A.vertices.begin(), A.vertices.end(), B.vertices.begin(),
                          B.vertices.end(), std::back_inserter(shared));
    for (std::size_t ka = 0; ka < A.runs.size(); ka++)
    {
      const std::size_t a = A.runs[ka];
      const Run &ra = runs[a];
      if (ra.excluded)
      {
        continue;
      }
      std::vector<Interval> within;
      for (const std::size_t b : B.runs)
      {
        const Run &rb = runs[b];
        if (rb.excluded ||
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
        // Beyond either end of B: the half-plane past the end along the outward tangent.
        for (const bool at_end : {false, true})
        {
          const Run &terminal = runs[B.runs[at_end ? B.runs.size() - 1 : 0]];
          const Point3D e = at_end ? terminal.end : terminal.start;
          const Point3D t_out = at_end ? terminal.tangent : Scale(-1.0, terminal.tangent);
          const double g0 = Dot(Sub(ra.start, e), t_out), slope = Dot(ra.tangent, t_out);
          // g(s) = g0 + slope s > 0.
          if (std::abs(slope) <= kDirectionQuantum)
          {
            if (g0 > 0.0)
            {
              excluded.emplace_back(0.0, ra.length);
            }
          }
          else
          {
            const double root = -g0 / slope;
            if (slope > 0.0)
            {
              excluded.emplace_back(std::max(root, 0.0), ra.length);
            }
            else
            {
              excluded.emplace_back(0.0, std::min(root, ra.length));
            }
          }
        }
      }
      for (const std::size_t v : shared)
      {
        const Point3D &p = input.vertices[v].coordinate;
        const auto z = RunIntervalWithin(a, p, p, zone);
        excluded.insert(excluded.end(), z.begin(), z.end());
      }
      const auto paired = SubtractIntervals(within, MergeIntervals(excluded, Tol()), Tol());
      for (const auto &interval : paired)
      {
        // Cuts at A's curved boundaries and at B's curved boundaries mapped onto the run.
        std::vector<double> cuts = {interval.first, interval.second};
        for (const auto &ci : A.curved)
        {
          for (const double x : {ci.first, ci.second})
          {
            const double s = x - A.run_offset[ka];
            if (s > interval.first + Tol() && s < interval.second - Tol())
            {
              cuts.push_back(s);
            }
          }
        }
        for (const auto &ci : B.curved)
        {
          for (const double x : {ci.first, ci.second})
          {
            // Run of B containing chain position x.
            const std::size_t kb = std::min(
                static_cast<std::size_t>(
                    std::upper_bound(B.run_offset.begin(), B.run_offset.end(), x) -
                    B.run_offset.begin()) -
                    1,
                B.runs.size() - 1);
            const std::size_t rb = B.runs[kb];
            const Point3D q =
                runs[rb].At(std::clamp(x - B.run_offset[kb], 0.0, runs[rb].length));
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
          const auto q_mid = ClosestPointOnChain(B, ra.At(0.5 * (piece.first + piece.second)));
          const bool curved = IsCurvedAt(A, mid) || IsCurvedAt(B, q_mid.x);
          const double max_kappa =
              std::max(MaxCurvature(A, x0, x1), WindowedCurvature(B, q_mid.x));
          const int n_samples = std::max(
              kPairSeparationSamplesPerInterval,
              static_cast<int>(std::ceil(2.0 * (piece.second - piece.first) / R)));
          Piece result{a, piece, curved, max_kappa, {}};
          for (int i_s = 0; i_s <= n_samples; i_s++)
          {
            const double s = piece.first + (piece.second - piece.first) * i_s / n_samples;
            const double x = A.run_offset[ka] + s;
            const auto q = ClosestPointOnChain(B, ra.At(s));
            const double half_own = WindowedCurvature(A, x) > 0.0 ? std::max(R, ra.length) : R;
            const double half_other =
                WindowedCurvature(B, q.x) > 0.0 ? std::max(R, runs[q.run].length) : R;
            const double turn =
                std::max(LocalTurn(A, ka), LocalTurn(B, RunIndexInChain(B, q.run)));
            result.samples.push_back({a, s, x, q.distance, q.x, half_own, half_other, turn});
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
    auto WindowMax = [&](const std::vector<const Sample *> &list, double x, double half,
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
        constant.push_back(!quantizer.Less(kPairSeparationTolerance * min_d, max_own - min_d));
        const double w_own = WindowMax(own_samples, sample.x, sample.half_own, nullptr);
        const double w_other = other_samples.empty()
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
        const double s_lo =
            i == 0 ? piece.interval.first
                   : 0.5 * (piece.samples[i - 1].s + piece.samples[i].s);
        const double s_hi = j == n ? piece.interval.second
                                   : 0.5 * (piece.samples[j - 1].s + piece.samples[j].s);
        if (s_hi - s_lo > Tol())
        {
          out.push_back({piece.run, {s_lo, s_hi}, piece.curved, piece.max_kappa, c, inter,
                         weighted, static_cast<int>(j - i)});
        }
        i = j;
      }
      k += n;
    }
  };

  for (std::size_t ca = 0; ca < chains.size(); ca++)
  {
    for (std::size_t cb = ca + 1; cb < chains.size(); cb++)
    {
      const Chain &A = chains[ca];
      const Chain &B = chains[cb];
      if (!usable[ca] || !usable[cb])
      {
        continue;
      }
      if (A.Rigid() && B.Rigid() &&
          !DirectionLess(std::abs(Dot(runs[A.runs.front()].tangent,
                                      runs[B.runs.front()].tangent)),
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
      std::vector<Piece> pieces_a, pieces_b;
      Pieces(A, B, pieces_a);
      Pieces(B, A, pieces_b);
      if (pieces_a.empty() || pieces_b.empty())
      {
        continue;
      }
      std::vector<char> constant_a, constant_b;
      std::vector<double> separation_a, separation_b, upper_a, upper_b;
      Classify(pieces_a, pieces_b, constant_a, separation_a, upper_a);
      Classify(pieces_b, pieces_a, constant_b, separation_b, upper_b);
      std::vector<SubPiece> sub_a, sub_b;
      Split(pieces_a, constant_a, separation_a, upper_a, sub_a);
      Split(pieces_b, constant_b, separation_b, upper_b, sub_b);
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
      // Frame: lateral from A to B at the first interacting piece; gap signs relative to it.
      const SubPiece *lead = nullptr;
      for (const auto &piece : sub_a)
      {
        if (piece.constant && piece.interacting)
        {
          lead = &piece;
          break;
        }
      }
      if (!lead)
      {
        for (const auto &piece : sub_b)
        {
          if (piece.constant && piece.interacting)
          {
            lead = &piece;
            break;
          }
        }
      }
      const bool lead_on_a = lead >= sub_a.data() && lead < sub_a.data() + sub_a.size();
      const Chain &lead_other = lead_on_a ? B : A;
      const Point3D pa = runs[lead->run].At(0.5 * (lead->interval.first + lead->interval.second));
      const auto qb = ClosestPointOnChain(lead_other, pa);
      const Point3D lateral = Normalize(Sub(runs[qb.run].At(qb.s), pa));
      const Run &ra = runs[lead->run];
      const Run &rb = runs[qb.run];
      const int gap_a = Dot(ra.gap_direction, lateral) > 0.0 ? 1 : -1;
      const int gap_b = Dot(rb.gap_direction, lateral) > 0.0 ? 1 : -1;
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
        for (const auto *list : {&sub_a, &sub_b})
        {
          for (const auto &piece : *list)
          {
            if (!piece.constant || !piece.interacting || piece.curved != curved_class)
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
        // One separation per pair and curvature class: the mean chord reading of its
        // samples (exact for a constant pair; a slow taper is described by its mean).
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
        features[feature].origin = pa;
        features[feature].axes = {ra.tangent, lateral, n_ref};
        for (const auto *list : {&sub_a, &sub_b})
        {
          for (const auto &piece : *list)
          {
            if (!piece.constant || !piece.interacting || piece.curved != curved_class)
            {
              continue;
            }
            claims[piece.run].push_back({feature, 2, piece.interval});
          }
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
  std::map<std::array<long long int, 3>, std::vector<std::size_t>> classes;
  for (std::size_t r = 0; r < runs.size(); r++)
  {
    if (runs[r].excluded || !chains[chain_index.at(runs[r].chain)].Rigid())
    {
      continue;  // chains with joints pair through the curved-edge chain rule
    }
    classes[DirectionKey(SignCanonical(runs[r].tangent), 1.0e-9)].push_back(r);
  }
  const double interaction = kInteractionDistanceOverRadius * R;
  for (const auto &[key, members_in_class] : classes)
  {
    (void)key;
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
    // Components per elementary interval, then merge consecutive equal components.
    struct Span
    {
      std::vector<std::size_t> component;  // member indices sorted by w
      double lo, hi;
    };
    std::vector<Span> spans;
    for (std::size_t k = 0; k + 1 < unique_splits.size(); k++)
    {
      const double lo = unique_splits[k], hi = unique_splits[k + 1];
      std::vector<std::size_t> active;
      for (std::size_t m = 0; m < members.size(); m++)
      {
        if (members[m].u0 <= lo + Tol() && members[m].u1 >= hi - Tol())
        {
          active.push_back(m);
        }
      }
      if (active.size() < 2)
      {
        continue;
      }
      UnionFind uf(active.size());
      for (std::size_t i = 0; i < active.size(); i++)
      {
        for (std::size_t j = i + 1; j < active.size(); j++)
        {
          const auto &mi = members[active[i]], &mj = members[active[j]];
          if (runs[mi.run].chain != runs[mj.run].chain &&
              quantizer.Less(std::abs(mi.w - mj.w), interaction))
          {
            uf.Union(i, j);
          }
        }
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
      std::vector<TranslationalEdge> edges;
      for (const std::size_t m : span.component)
      {
        const auto &member = members[m];
        edges.push_back({member.w, member.gap_sign, runs[member.run].conductor,
                         InterfaceNames(runs[member.run].targets),
                         runs[member.run].boundary_law});
      }
      auto translational = CanonicalTranslationalSignature(edges, R);
      nlohmann::json best = std::move(translational.signature);
      const int chirality = translational.chirality;
      std::string type;
      std::optional<std::pair<std::string, std::string>> exclusion;
      if (span.component.size() == 2)
      {
        const auto &lower = members[span.component[0]], &upper = members[span.component[1]];
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
      features[feature].origin = Scale(0.5 * (span.lo + span.hi), axis);
      features[feature].axes = {axis, lateral, n_ref};
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
        if (s1 - s0 <= Tol())
        {
          continue;
        }
        claims[member.run].push_back({feature, 2, {s0, s1}});
      }
    }
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
  const auto position = std::find(chain.runs.begin(), chain.runs.end(), run);
  MFEM_VERIFY(position != chain.runs.end(), "Run missing from its chain!");
  const std::size_t start_index = static_cast<std::size_t>(position - chain.runs.begin());
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
      cores.push_back({a, interval, ra.At(interval.first), ra.At(interval.second)});
    }
  };

  for (std::size_t a = 0; a < runs.size(); a++)
  {
    if (runs[a].excluded)
    {
      continue;
    }
    for (std::size_t b = a + 1; b < runs.size(); b++)
    {
      if (runs[b].excluded || runs[a].chain == runs[b].chain)
      {
        continue;
      }
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
      std::vector<std::size_t> shared;
      std::set_intersection(ca.vertices.begin(), ca.vertices.end(), cb.vertices.begin(),
                            cb.vertices.end(), std::back_inserter(shared));
      std::vector<std::vector<Interval>> zones_a, zones_b;
      for (const std::size_t v : shared)
      {
        const Point3D &p = input.vertices[v].coordinate;
        zones_a.push_back(RunIntervalWithin(a, p, p, zone));
        zones_b.push_back(RunIntervalWithin(b, p, p, zone));
      }
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

  // Two vertex sites within 2R have overlapping radius-R windows (invariant A2): they are
  // an interaction event of their own, with the two points as degenerate cores.
  for (std::size_t i = 0; i < sites.size(); i++)
  {
    for (std::size_t j = i + 1; j < sites.size(); j++)
    {
      if (quantizer.Less(Distance(sites[i].point, sites[j].point), interaction))
      {
        cores.push_back({std::numeric_limits<std::size_t>::max(),
                         {0.0, 0.0},
                         sites[i].point,
                         sites[i].point});
        cores.push_back({std::numeric_limits<std::size_t>::max(),
                         {0.0, 0.0},
                         sites[j].point,
                         sites[j].point});
      }
    }
  }

  // Connected union of the radius-R balls: cores whose distance is below 2R.
  UnionFind uf(cores.size());
  for (std::size_t i = 0; i < cores.size(); i++)
  {
    for (std::size_t j = i + 1; j < cores.size(); j++)
    {
      if (quantizer.Less(
              SegmentSegmentDistance(cores[i].p0, cores[i].p1, cores[j].p0, cores[j].p1),
              interaction))
      {
        uf.Union(i, j);
      }
    }
  }
  // A vertex site whose 2R through-vertex zone reaches into a region joins it (and merges
  // the regions it reaches).
  std::vector<std::optional<std::size_t>> site_root(sites.size());
  for (std::size_t s = 0; s < sites.size(); s++)
  {
    std::optional<std::size_t> root;
    for (std::size_t i = 0; i < cores.size(); i++)
    {
      if (quantizer.Less(PointSegmentDistance(sites[s].point, cores[i].p0, cores[i].p1),
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

  // Claimed portions: run intervals within R of a core, plus the member sites' windows.
  for (std::size_t c = 0; c < cluster_cores.size(); c++)
  {
    std::vector<SignaturePortion> portions;
    std::vector<SignatureVertex> vertices;
    std::vector<std::pair<std::size_t, Interval>> claimed;
    for (std::size_t r = 0; r < runs.size(); r++)
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
    const auto canonical = CanonicalClusterSignature(portions, vertices, n_ref, R);
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
    if (site.type == "ConvexCorner" || site.type == "ConcaveCorner")
    {
      signature = CanonicalCornerSignature(site.interfaces, site.boundary_law,
                                           site.angle_degrees, site.corner_radius / R);
    }
    else if (site.type == "Junction")
    {
      signature = CanonicalJunctionSignature(site.interfaces, site.boundary_law,
                                             site.arm_angles, site.arm_conductors);
    }
    const int feature = NewFeature(site.type, signature);
    features[feature].origin = site.point;
    features[feature].axes = {Point3D{}, Point3D{}, n_ref};
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
  std::vector<std::vector<std::tuple<double, double, int>>> assigned(runs.size());
  for (std::size_t r = 0; r < runs.size(); r++)
  {
    if (runs[r].excluded)
    {
      continue;
    }
    auto &run_claims = claims[r];
    std::sort(run_claims.begin(), run_claims.end(),
              [](const Claim &a, const Claim &b)
              {
                return std::tie(a.priority, a.feature, a.interval) <
                       std::tie(b.priority, b.feature, b.interval);
              });
    // CrossLayer zones are excluded before any feature claims the run.
    std::vector<Interval> taken = cross_layer[r];
    for (const auto &claim : run_claims)
    {
      for (const auto &piece : SubtractIntervals({claim.interval}, taken, Tol()))
      {
        assigned[r].emplace_back(piece.first, piece.second, claim.feature);
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
          assigned[r].emplace_back(piece.first, piece.second, -1);  // unclaimed sliver
        }
      }
      std::sort(assigned[r].begin(), assigned[r].end());
      auto &pieces = assigned[r];
      for (std::size_t i = 0; i < pieces.size();)
      {
        auto &[lo, hi, feature] = pieces[i];
        if (hi - lo > sliver && feature >= 0)
        {
          i++;
          continue;
        }
        auto joins = [&](std::size_t j)
        {
          const auto &[jlo, jhi, jfeature] = pieces[j];
          return jfeature >= 0 && (std::abs(jhi - lo) <= Tol() || std::abs(hi - jlo) <= Tol());
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
      taken = cross_layer[r];
      for (const auto &[lo, hi, feature] : pieces)
      {
        taken.push_back({lo, hi});
      }
      taken = MergeIntervals(taken, Tol());
    }
    const auto remainder = SubtractIntervals({Interval{0.0, runs[r].length}}, taken, Tol());
    if (!remainder.empty())
    {
      // Curved sections of the chain (windowed bend radius below the straight threshold)
      // are CurvedEdge features, one per section; the rest is the chain's isolated edge.
      const Chain &chain = chains[chain_index.at(runs[r].chain)];
      const double offset = chain.run_offset[RunIndexInChain(chain, r)];
      std::vector<Interval> curved_on_run;
      std::vector<std::size_t> curved_section;
      for (std::size_t j = 0; j < chain.curved.size(); j++)
      {
        const Interval local{chain.curved[j].first - offset, chain.curved[j].second - offset};
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
        const auto key = std::make_pair(
            runs[r].chain, signature.dump() + "#" + std::to_string(curved_section[i]));
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
                                 it->second);
      }
      const auto straight = SubtractIntervals(remainder, curved_on_run, Tol());
      if (!straight.empty())
      {
        nlohmann::json signature = {{"Interfaces", InterfaceNames(runs[r].targets)},
                                    {"Law", runs[r].boundary_law}};
        signature["Type"] = "IsolatedEdge";
        const auto key = std::make_pair(runs[r].chain, signature.dump());
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
          assigned[r].emplace_back(piece.first, piece.second, it->second);
        }
      }
    }
    std::sort(assigned[r].begin(), assigned[r].end());
    // Curvature annotation of every feature over its assigned portions.
    {
      const Chain &chain = chains[chain_index.at(runs[r].chain)];
      const double offset = chain.run_offset[RunIndexInChain(chain, r)];
      for (const auto &[lo, hi, feature] : assigned[r])
      {
        feature_max_kappa[feature] =
            std::max(feature_max_kappa[feature], MaxCurvature(chain, offset + lo, offset + hi));
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
  result.exclusions.push_back(
      {"CrossLayer",
       "planar metal edge within 2R of metal off its own plane (facing layer, wall, staple)",
       0, 0.0});
  for (std::size_t r = 0; r < runs.size(); r++)
  {
    if (runs[r].excluded)
    {
      continue;
    }
    for (const auto &rs : runs[r].segments)
    {
      auto &table = result.segments[rs.segment];
      MFEM_VERIFY(!table.exclusion, "An excluded segment received a feature portion!");
      const double segment_length =
          Distance(input.segments[rs.segment].p0, input.segments[rs.segment].p1);
      const double scale = segment_length / (rs.t1 - rs.t0);
      auto SegmentPortion = [&](double lo, double hi) -> std::optional<std::array<double, 2>>
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
      for (const auto &[lo, hi, feature] : assigned[r])
      {
        const auto portion = SegmentPortion(lo, hi);
        if (!portion)
        {
          continue;
        }
        const auto [s0, s1] = *portion;
        table.portions.push_back({s0, s1, static_cast<double>(feature)});
        features[feature].portions.push_back({rs.segment, s0, s1});
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
                       : "TruncationCut";
      result.vertices.push_back(entry);
      continue;
    }
    const auto site =
        std::find_if(sites.begin(), sites.end(), [&](const VertexFeatureSite &s)
                     { return s.vertex && *s.vertex == v; });
    if (site == sites.end())
    {
      entry.type = "Excluded";
      result.vertices.push_back(entry);
      continue;
    }
    entry.type = site->type;
    entry.turn_degrees = site->turn_degrees;
    const auto feature = vertex_feature.find(v);
    entry.feature = feature == vertex_feature.end() ? -1 : feature->second;
    result.vertices.push_back(entry);
  }
  for (const auto &site : sites)
  {
    if (site.vertex)
    {
      continue;
    }
    IdentifiedVertex entry;
    entry.vertex = std::numeric_limits<std::size_t>::max();
    entry.type = "RoundedCorner";
    entry.turn_degrees = site.turn_degrees;
    entry.feature = -1;
    for (const auto &feature : features)
    {
      if (site.cluster < 0 && feature.origin == site.point &&
          (feature.type == "ConvexCorner" || feature.type == "ConcaveCorner"))
      {
        entry.feature = feature.id;
      }
    }
    if (site.cluster >= 0)
    {
      for (const auto &feature : features)
      {
        if (feature.type == "SpatialEdgeCluster" &&
            std::any_of(feature.portions.begin(), feature.portions.end(),
                        [&](const IdentifiedPortion &p)
                        {
                          return std::any_of(site.window.begin(), site.window.end(),
                                             [&](const auto &w)
                                             {
                                               return std::any_of(
                                                   runs[w.first].segments.begin(),
                                                   runs[w.first].segments.end(),
                                                   [&](const RunSegment &rs)
                                                   { return rs.segment == p.segment; });
                                             });
                        }))
        {
          entry.feature = feature.id;
          break;
        }
      }
    }
    result.vertices.push_back(entry);
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
  result.features = features;
  result.radius = R;
  result.reference_process_normal = n_ref;
}

IdentificationResult Identifier::Identify()
{
  IdentificationResult result;
  segment_exclusion.assign(input.segments.size(), std::nullopt);
  BuildRuns(input, runs, chains, chain_index);
  claims.assign(runs.size(), {});
  bent_claims.assign(runs.size(), {});
  if (!runs.empty())
  {
    ClassifyPlanes();
    ClassifyVertices();
    DetectRoundedCorners();
    ComputeCurvature();
    BuildVertexWindows();
    BuildBentPairs();
    BuildClusters();
    BuildTranslationalFeatures();
  }
  Assign(result);
  return result;
}

}  // namespace

IdentificationResult IdentifyMetalPerimeter(const IdentificationInput &input)
{
  Identifier identifier(input);
  return identifier.Identify();
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
    for (const auto &portion : feature.portions)
    {
      portions.push_back({portion.segment, L(portion.s0), L(portion.s1)});
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
        {"BendRadiusOverR", feature.bend_radius_over_R
                                ? nlohmann::json(std::round(*feature.bend_radius_over_R /
                                                            kSignatureLengthQuantumOverRadius) *
                                                 kSignatureLengthQuantumOverRadius)
                                : nlohmann::json(nullptr)},
        {"Match", {{"Status", feature.matched_model ? "Matched" : "Missing"}}}};
    if (feature.matched_model)
    {
      entry["Match"]["Model"] = *feature.matched_model;
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
    if (vertex.type != "TruncationCut" && vertex.type != "Excluded" &&
        vertex.type != "ExclusionCut")
    {
      entry["TurnDegrees"] = std::round(vertex.turn_degrees * 1.0e6) * 1.0e-6;
      entry["Feature"] = vertex.feature;
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
  return {{"Version", 2},
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
             "the chain bends, R on straight runs; inscribed reading C / cos(turn / 2) with the "
             "larger local joint turn; interacting iff both < 2R; feature separation = mean C"},
            {"PairConstancyWindowOverR", 1.0},
            {"PairSampleSpacingOverR", 0.5},
            {"PairCandidateReachOverR",
             kInteractionDistanceOverRadius * (1.0 + kPairSeparationTolerance)},
            {"CrossLayerReachOverR", kInteractionDistanceOverRadius},
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
          {"GeometryDigest", geometry_digest}};
}

}  // namespace palace
