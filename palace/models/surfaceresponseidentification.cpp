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
  auto Less = [](const Point3D &a, const Point3D &b) { return a < b; };
  for (const auto &[chain_id, members] : segments_by_chain)
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
    return std::array<double, 2>{RoundTo(Dot(r, x) / radius, 1.0e-10),
                                 RoundTo(Dot(r, y) / radius, 1.0e-10)};
  };
  auto LocalDirection = [&](const Point3D &v)
  { return std::array<double, 2>{RoundTo(Dot(v, x), 1.0e-8), RoundTo(Dot(v, y), 1.0e-8)}; };
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
    vertex_entries.push_back({{"P", {p[0], p[1]}},
                              {"Type", vertex.type},
                              {"TurnDegrees", RoundTo(vertex.turn_degrees, 1.0e-6)}});
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
      list.push_back(
          {{"OffsetOverR", RoundTo(orientation * (edge.offset - w0) / radius, 1.0e-10)},
           {"GapSide", orientation * edge.gap_sign},
           {"Conductor", it->second},
           {"Interfaces", edge.interfaces},
           {"Law", edge.boundary_law}});
    }
    nlohmann::json candidate = {{"Edges", list}};
    if (edges.size() == 2)
    {
      candidate["SeparationOverR"] =
          RoundTo(std::abs(edges.back().offset - edges.front().offset) / radius, 1.0e-10);
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
  std::vector<IdentifiedFeature> features;
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
    return true;
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
  void BuildTranslationalFeatures();
  void BuildClusters();
  void BuildVertexWindows();
  void Assign(IdentificationResult &result);
  std::vector<Interval> RunIntervalWithin(std::size_t run, const Point3D &a,
                                          const Point3D &b, double distance) const;
  std::vector<Interval>
  ChainWindow(std::size_t run, std::size_t from_vertex, double length,
              std::vector<std::pair<std::size_t, Interval>> &out) const;
};

// Runs whose process normal is not parallel to the reference normal are walls / staples
// (NonPlanar); planar runs off the primary plane (largest perimeter length) are CrossLayer.
void Identifier::ClassifyPlanes()
{
  n_ref = ReferenceProcessNormal(runs);
  std::vector<std::pair<double, std::size_t>> offsets;
  for (std::size_t r = 0; r < runs.size(); r++)
  {
    if (DirectionLess(std::abs(Dot(runs[r].process_normal, n_ref)),
                      1.0 - kParallelCosineTolerance))
    {
      ExcludeRun(r, "NonPlanar",
                 "metal edge whose process normal is not parallel to the reference process "
                 "normal (wall, staple, via)");
      continue;
    }
    offsets.emplace_back(Dot(Scale(0.5, Add(runs[r].start, runs[r].end)), n_ref), r);
  }
  std::sort(offsets.begin(), offsets.end());
  std::vector<std::pair<double, std::vector<std::size_t>>> planes;  // offset, runs
  for (const auto &[offset, r] : offsets)
  {
    if (planes.empty() || !quantizer.Equal(offset, planes.back().first))
    {
      planes.emplace_back(offset, std::vector<std::size_t>{});
    }
    planes.back().second.push_back(r);
  }
  std::size_t primary = 0;
  double primary_length = -1.0;
  for (std::size_t p = 0; p < planes.size(); p++)
  {
    double length = 0.0;
    for (const std::size_t r : planes[p].second)
    {
      length += runs[r].length;
    }
    if (length > primary_length)
    {
      primary_length = length;
      primary = p;
    }
  }
  for (std::size_t p = 0; p < planes.size(); p++)
  {
    if (p == primary)
    {
      continue;
    }
    for (const std::size_t r : planes[p].second)
    {
      ExcludeRun(r, "CrossLayer",
                 "planar metal edge off the primary process plane (facing layer)");
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
      std::vector<double> angles;
      for (const std::size_t r : incident)
      {
        const Point3D d = ArmDirection(r, v);
        double angle = std::atan2(Dot(d, y), Dot(d, x)) * 180.0 / std::acos(-1.0);
        if (angle < 0.0)
        {
          angle += 360.0;
        }
        angles.push_back(angle);
      }
      std::sort(angles.begin(), angles.end());
      for (std::size_t i = 0; i < angles.size(); i++)
      {
        const double next = i + 1 < angles.size() ? angles[i + 1] : angles[0] + 360.0;
        site.arm_angles.push_back(RoundTo(next - angles[i], 1.0e-6));
      }
      // Canonical cyclic order: minimal rotation, both orientations (mirror).
      std::vector<double> best = site.arm_angles;
      for (const bool reverse : {false, true})
      {
        std::vector<double> sequence = site.arm_angles;
        if (reverse)
        {
          std::reverse(sequence.begin(), sequence.end());
        }
        for (std::size_t shift = 0; shift < sequence.size(); shift++)
        {
          std::rotate(sequence.begin(), sequence.begin() + 1, sequence.end());
          if (sequence < best)
          {
            best = sequence;
          }
        }
      }
      site.arm_angles = best;
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
  for (const auto &chain : chains)
  {
    const std::size_t m = chain.runs.size();
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
              }
            }
          }
        }
      }
      step += count;
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
    if (runs[r].excluded)
    {
      continue;
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
      if (!DirectionLess(std::abs(Dot(runs[a].tangent, runs[b].tangent)),
                         1.0 - kParallelCosineTolerance))
      {
        continue;  // parallel: a translational interaction
      }
      if (!quantizer.Less(SegmentSegmentDistance(runs[a].start, runs[a].end, runs[b].start,
                                                 runs[b].end),
                          interaction))
      {
        continue;
      }
      const Chain &ca = chains[chain_index.at(runs[a].chain)];
      const Chain &cb = chains[chain_index.at(runs[b].chain)];
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
      vertices.push_back(
          {sites[s].point, sites[s].type, RoundTo(sites[s].turn_degrees, 1.0e-6)});
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
      signature["AngleDegrees"] = RoundTo(site.angle_degrees, 1.0e-6);
      signature["CornerRadiusOverR"] = RoundTo(site.corner_radius / R, 1.0e-10);
    }
    else if (site.type == "Junction")
    {
      signature["ArmAnglesDegrees"] = site.arm_angles;
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
    std::vector<Interval> taken;
    for (const auto &claim : run_claims)
    {
      for (const auto &piece : SubtractIntervals({claim.interval}, taken, Tol()))
      {
        assigned[r].emplace_back(piece.first, piece.second, claim.feature);
        taken.push_back(piece);
      }
      taken = MergeIntervals(taken, Tol());
    }
    const auto remainder = SubtractIntervals({Interval{0.0, runs[r].length}}, taken, Tol());
    if (!remainder.empty())
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
      for (const auto &piece : remainder)
      {
        assigned[r].emplace_back(piece.first, piece.second, it->second);
      }
    }
    std::sort(assigned[r].begin(), assigned[r].end());
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
      for (const auto &[lo, hi, feature] : assigned[r])
      {
        const double a = std::max(lo, rs.t0), b = std::min(hi, rs.t1);
        if (b - a <= Tol())
        {
          continue;
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
        table.portions.push_back({s0, s1, static_cast<double>(feature)});
        features[feature].portions.push_back({rs.segment, s0, s1});
        features[feature].length += s1 - s0;
        result.assigned_length += s1 - s0;
      }
      std::sort(table.portions.begin(), table.portions.end());
    }
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

  // Vertex table.
  for (std::size_t v = 0; v < input.vertices.size(); v++)
  {
    const auto &vertex = input.vertices[v];
    if (!vertex.physical_type || *vertex.physical_type == MetalEdgeVertexType::REGULAR)
    {
      continue;
    }
    IdentifiedVertex entry;
    entry.vertex = v;
    if (!IsFeatureVertex(v))
    {
      entry.type = "TruncationCut";
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
  if (!runs.empty())
  {
    ClassifyPlanes();
    ClassifyVertices();
    DetectRoundedCorners();
    BuildVertexWindows();
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
    if (vertex.type != "TruncationCut" && vertex.type != "Excluded")
    {
      entry["TurnDegrees"] = std::round(vertex.turn_degrees * 1.0e6) * 1.0e-6;
      entry["Feature"] = vertex.feature;
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
