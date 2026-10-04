// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "models/surfaceresponseidentification.hpp"
#include "models/surfaceresponseoperator.hpp"
#include "utils/metaledge.hpp"

using namespace palace;
using namespace Catch::Matchers;

namespace
{

using Point2 = std::array<double, 2>;

// Plan-view metal loops (counter-clockwise: metal inside) subdivided into mesh segments of
// at most the given length, with the perimeter frames the classifier would supply (gap
// direction = outward normal, process normal = normal_sign x z, one PEC conductor per loop
// or a shared one).
struct LoopSpec
{
  std::vector<Point2> points;
  int conductor = 0;
  double subdivision = 1.0;
  double z = 0.0;  // metal plane of the loop
  // Polygon edges (index i: points[i] -> points[i + 1]) lying on a simulation cut: their
  // mesh segments are truncation segments (excluded, no chain) and their vertices lie on
  // the truncation boundary (ENDPOINT of the physical chain, never a feature), as the
  // perimeter extraction reports a metal polygon clipped by a window wall.
  std::set<std::size_t> truncation_edges;
  // Sign of the loop's process normal (+1: metal facing +z, substrate below; -1: a flipped
  // plane, a flip-chip top chip facing down with its substrate above).
  double normal_sign = 1.0;
};

IdentificationInput MakeInput(const std::vector<LoopSpec> &loops, double radius,
                              bool reverse_segment_order = false)
{
  IdentificationInput input;
  input.radius = radius;
  if (std::getenv("PALACE_IDENTIFICATION_TEST_LOG"))
  {
    input.log = [](const std::string &line) { std::cout << line; };
  }
  std::map<std::array<long long, 3>, std::size_t> vertex_index;
  auto Vertex = [&](const Point2 &p, double z)
  {
    const std::array<long long, 3> key = {
        std::llround(p[0] * 1.0e9), std::llround(p[1] * 1.0e9), std::llround(z * 1.0e9)};
    auto it = vertex_index.find(key);
    if (it == vertex_index.end())
    {
      it = vertex_index.emplace(key, input.vertices.size()).first;
      IdentificationVertex vertex;
      vertex.coordinate = {p[0], p[1], z};
      input.vertices.push_back(vertex);
    }
    return it->second;
  };
  struct Raw
  {
    Point2 a, b, outward;
    int conductor;
    int chain;
    double z;
    bool truncation;
    double normal_sign;
  };
  std::vector<Raw> raw;
  int chain = 0;
  for (const auto &loop : loops)
  {
    const std::size_t n = loop.points.size();
    // Chains break at corners (joints that are not noise under the geometric rule): every
    // polygon edge here is its own chain.
    for (std::size_t i = 0; i < n; i++)
    {
      const Point2 a = loop.points[i], b = loop.points[(i + 1) % n];
      const double length = std::hypot(b[0] - a[0], b[1] - a[1]);
      const Point2 t = {(b[0] - a[0]) / length, (b[1] - a[1]) / length};
      const Point2 outward = {t[1], -t[0]};  // right of a CCW edge
      const int pieces =
          std::max(1, static_cast<int>(std::ceil(length / loop.subdivision - 1.0e-9)));
      for (int k = 0; k < pieces; k++)
      {
        const double s0 = length * k / pieces, s1 = length * (k + 1) / pieces;
        raw.push_back({{a[0] + t[0] * s0, a[1] + t[1] * s0},
                       {a[0] + t[0] * s1, a[1] + t[1] * s1},
                       outward,
                       loop.conductor,
                       chain,
                       loop.z,
                       loop.truncation_edges.count(i) > 0,
                       loop.normal_sign});
      }
      chain++;
    }
  }
  if (reverse_segment_order)
  {
    std::reverse(raw.begin(), raw.end());
  }
  for (const auto &r : raw)
  {
    IdentificationSegment segment;
    segment.p0 = {r.a[0], r.a[1], r.z};
    segment.p1 = {r.b[0], r.b[1], r.z};
    segment.vertices = {Vertex(r.a, r.z), Vertex(r.b, r.z)};
    segment.chain = r.truncation ? -1 : r.chain;
    segment.truncation = r.truncation;
    segment.conductor = r.conductor;
    segment.targets = {{InterfaceDielectric::MS, 1}};
    segment.gap_direction = {r.outward[0], r.outward[1], 0.0};
    segment.process_normal = {0.0, 0.0, r.normal_sign};
    segment.boundary_law = "{\"Type\":\"PEC\"}";
    for (const std::size_t v : segment.vertices)
    {
      input.vertices[v].segments.push_back(input.segments.size());
      input.vertices[v].on_truncation_boundary |= r.truncation;
    }
    input.segments.push_back(segment);
  }
  // The vertex classification and the chains see the physical segments only (metaledge.cpp
  // classifies the physical type on the non-truncation segments).
  auto PhysicalSegments = [&](std::size_t v)
  {
    std::vector<std::size_t> physical;
    for (const std::size_t s : input.vertices[v].segments)
    {
      if (!input.segments[s].truncation)
      {
        physical.push_back(s);
      }
    }
    return physical;
  };
  // Vertex types as metaledge.cpp: two segments -> a collinear continuation is regular,
  // otherwise corner iff the joint is not noise under the geometric rule (the implied
  // sagitta (c / 2) tan(turn / 4) on the shorter adjacent straight piece reaches
  // kJointNoiseSagittaOverRadius x R; the pieces merge the collinear subdivisions).
  auto UnitDirection = [&](std::size_t from, std::size_t to)
  {
    std::array<double, 3> d{};
    double norm = 0.0;
    for (int k = 0; k < 3; k++)
    {
      d[k] = input.vertices[to].coordinate[k] - input.vertices[from].coordinate[k];
      norm += d[k] * d[k];
    }
    for (double &value : d)
    {
      value /= std::sqrt(norm);
    }
    return std::make_pair(d, std::sqrt(norm));
  };
  auto OtherEnd = [&](std::size_t segment, std::size_t vertex)
  {
    const auto &ends = input.segments[segment].vertices;
    return ends[0] == vertex ? ends[1] : ends[0];
  };
  auto PieceLength = [&](std::size_t vertex, std::size_t segment)
  {
    double length = 0.0;
    std::size_t current = vertex, s = segment;
    for (std::size_t steps = 0; steps < input.segments.size(); steps++)
    {
      const std::size_t other = OtherEnd(s, current);
      const auto [in, piece] = UnitDirection(current, other);
      length += piece;
      const auto &next_segments = input.vertices[other].segments;
      if (next_segments.size() != 2 || other == vertex)
      {
        break;
      }
      const std::size_t next = next_segments[0] == s ? next_segments[1] : next_segments[0];
      const auto [out, ignored] = UnitDirection(other, OtherEnd(next, other));
      (void)ignored;
      if (in[0] * out[0] + in[1] * out[1] + in[2] * out[2] < 1.0 - 1.0e-12)
      {
        break;
      }
      current = other;
      s = next;
    }
    return length;
  };
  for (std::size_t v = 0; v < input.vertices.size(); v++)
  {
    auto &vertex = input.vertices[v];
    const auto physical = PhysicalSegments(v);
    if (physical.size() <= 1)
    {
      vertex.physical_type = MetalEdgeVertexType::ENDPOINT;
    }
    else if (physical.size() > 2)
    {
      vertex.physical_type = MetalEdgeVertexType::JUNCTION;
    }
    else
    {
      const auto d0 = UnitDirection(v, OtherEnd(physical[0], v)).first;
      const auto d1 = UnitDirection(v, OtherEnd(physical[1], v)).first;
      const double dot = d0[0] * d1[0] + d0[1] * d1[1] + d0[2] * d1[2];
      if (-dot >= 1.0 - 1.0e-12)
      {
        vertex.physical_type = MetalEdgeVertexType::REGULAR;
        continue;
      }
      const double turn = std::acos(std::clamp(-dot, -1.0, 1.0));
      const double shorter =
          std::min(PieceLength(v, physical[0]), PieceLength(v, physical[1]));
      vertex.physical_type =
          JointIsNoise(turn, shorter, kJointNoiseSagittaOverRadius * radius)
              ? MetalEdgeVertexType::REGULAR
              : MetalEdgeVertexType::CORNER;
    }
  }
  // Chains as metaledge.cpp: maximal paths through regular vertices (a polyline arc with
  // sub-threshold turns stays in one chain).
  std::vector<int> chain_of(input.segments.size(), -1);
  int next_chain = 0;
  for (std::size_t seed = 0; seed < input.segments.size(); seed++)
  {
    if (chain_of[seed] >= 0 || input.segments[seed].truncation)
    {
      continue;
    }
    std::vector<std::size_t> stack = {seed};
    chain_of[seed] = next_chain;
    while (!stack.empty())
    {
      const std::size_t current = stack.back();
      stack.pop_back();
      for (const std::size_t v : input.segments[current].vertices)
      {
        if (input.vertices[v].physical_type != MetalEdgeVertexType::REGULAR)
        {
          continue;
        }
        for (const std::size_t other : input.vertices[v].segments)
        {
          if (chain_of[other] < 0 && !input.segments[other].truncation)
          {
            chain_of[other] = next_chain;
            stack.push_back(other);
          }
        }
      }
    }
    next_chain++;
  }
  for (std::size_t i = 0; i < input.segments.size(); i++)
  {
    input.segments[i].chain = chain_of[i];
  }
  return input;
}

// Rounded rectangle (half extents, fillet radius) as a polyline: chords per quarter circle.
std::vector<Point2> RoundedRectangle(double half_x, double half_y, double radius,
                                     int chords)
{
  std::vector<Point2> points;
  const std::array<Point2, 4> centers = {Point2{half_x - radius, half_y - radius},
                                         Point2{-half_x + radius, half_y - radius},
                                         Point2{-half_x + radius, -half_y + radius},
                                         Point2{half_x - radius, -half_y + radius}};
  for (int corner = 0; corner < 4; corner++)
  {
    for (int k = 0; k <= chords; k++)
    {
      const double angle =
          (corner + static_cast<double>(k) / chords) * 0.5 * std::acos(-1.0);
      points.push_back({centers[corner][0] + radius * std::cos(angle),
                        centers[corner][1] + radius * std::sin(angle)});
    }
  }
  return points;
}

std::vector<Point2> Rectangle(double x0, double y0, double x1, double y1)
{
  return {{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}};
}

// The synthetic arc bar of synthetic_layouts.py: a bar of the given width following a
// circular arc (radius, sweep) discretised with the given turn per vertex, with straight
// leads at both ends, offset exactly along the vertex bisectors (counter-clockwise loop).
// offset shifts the bar across the centreline: it occupies the offsets [offset - width / 2,
// offset + width / 2] of the centreline polyline (two bars of opposite offsets about one
// centreline face each other across a gap whose chords are exactly parallel at the design
// separation, like an offset path).
// Counter-clockwise loop of a bar of the given width around a centreline polyline that
// turns left overall (the smaller-offset side forward, the larger back); `offset` shifts
// the bar sideways. Interior vertices are offset along `normals` (unit, to the left) when
// given — the exact radial direction on an arc, so that both sides are concentric polylines
// of the design circle at the same angles (as a CAD offset of an arc path discretises
// them); otherwise along the bisector of the adjacent chord normals (a mitred offset, whose
// vertices lie 1 / cos(turn / 2) off the concentric circle: a different geometry under the
// exact-parameter arc rule at coarse steps).
std::vector<Point2> BarAroundCentreline(const std::vector<Point2> &centreline, double width,
                                        double offset = 0.0,
                                        const std::vector<Point2> &normals = {})
{
  const double h = 0.5 * width;
  auto Offset = [&](double sign)
  {
    std::vector<Point2> result;
    const std::size_t n = centreline.size();
    auto Normal = [](const Point2 &d)
    {
      const double norm = std::hypot(d[0], d[1]);
      return Point2{-d[1] / norm, d[0] / norm};
    };
    for (std::size_t i = 0; i < n; i++)
    {
      const Point2 &p = centreline[i];
      if (!normals.empty())
      {
        const Point2 &nrm = normals[i];
        result.push_back(
            {p[0] + (offset + sign * h) * nrm[0], p[1] + (offset + sign * h) * nrm[1]});
      }
      else if (i == 0 || i + 1 == n)
      {
        const Point2 d =
            i == 0 ? Point2{centreline[1][0] - p[0], centreline[1][1] - p[1]}
                   : Point2{p[0] - centreline[i - 1][0], p[1] - centreline[i - 1][1]};
        const Point2 nrm = Normal(d);
        result.push_back(
            {p[0] + (offset + sign * h) * nrm[0], p[1] + (offset + sign * h) * nrm[1]});
      }
      else
      {
        const Point2 n0 =
            Normal({p[0] - centreline[i - 1][0], p[1] - centreline[i - 1][1]});
        const Point2 n1 =
            Normal({centreline[i + 1][0] - p[0], centreline[i + 1][1] - p[1]});
        Point2 b = {n0[0] + n1[0], n0[1] + n1[1]};
        const double norm = std::hypot(b[0], b[1]);
        b = {b[0] / norm, b[1] / norm};
        const double scale = (offset + sign * h) / (b[0] * n0[0] + b[1] * n0[1]);
        result.push_back({p[0] + scale * b[0], p[1] + scale * b[1]});
      }
    }
    return result;
  };
  // Counter-clockwise: the smaller offset side forward, the larger offset side back (the
  // centreline turns left, so the right side is the smaller offset).
  std::vector<Point2> points = Offset(-1.0);
  const std::vector<Point2> left = Offset(1.0);
  points.insert(points.end(), left.rbegin(), left.rend());
  return points;
}

// A bar along a left arc of the given radius and sweep (centre (0, radius)), straight leads
// at both ends; both sides are concentric arcs at the centreline's angles.
std::vector<Point2> ArcBar(double width, double radius, double sweep_degrees,
                           double step_degrees, double lead = 6.0, double offset = 0.0)
{
  const int steps =
      std::max(1, static_cast<int>(std::lround(sweep_degrees / step_degrees)));
  const double step = sweep_degrees * std::acos(-1.0) / 180.0 / steps;
  std::vector<Point2> centreline, normals;
  for (int k = 0; k <= steps; k++)
  {
    centreline.push_back(
        {radius * std::sin(k * step), radius - radius * std::cos(k * step)});
    normals.push_back({-std::sin(k * step), std::cos(k * step)});  // toward the centre
  }
  const Point2 d_end = {std::cos(steps * step), std::sin(steps * step)};
  centreline.insert(centreline.begin(), {-lead, 0.0});
  normals.insert(normals.begin(), {0.0, 1.0});
  centreline.push_back(
      {centreline.back()[0] + lead * d_end[0], centreline.back()[1] + lead * d_end[1]});
  normals.push_back({-d_end[1], d_end[0]});
  return BarAroundCentreline(centreline, width, offset, normals);
}

// A bar along an S-bend: a left arc of the given radius and sweep followed by a right arc
// of the same radius and sweep (the tangent returns to +x), straight leads at both ends;
// both sides concentric with the two design circles (the join vertex is offset along their
// common normal).
std::vector<Point2> SBar(double width, double radius, double sweep_degrees,
                         double step_degrees, double lead = 6.0)
{
  const int steps =
      std::max(1, static_cast<int>(std::lround(sweep_degrees / step_degrees)));
  const double step = sweep_degrees * std::acos(-1.0) / 180.0 / steps;
  std::vector<Point2> centreline = {{-lead, 0.0}}, normals = {{0.0, 1.0}};
  for (int k = 0; k <= steps; k++)
  {
    centreline.push_back(
        {radius * std::sin(k * step), radius - radius * std::cos(k * step)});
    normals.push_back({-std::sin(k * step), std::cos(k * step)});  // toward centre 1
  }
  // Second arc: centre at the reflection of the first centre through the join point, the
  // tangent turning back to +x.
  const double sweep = steps * step;
  const Point2 join = centreline.back();
  const Point2 centre2 = {2.0 * join[0], 2.0 * join[1] - radius};
  for (int k = 1; k <= steps; k++)
  {
    const double angle = sweep - k * step;  // tangent angle from +x
    centreline.push_back(
        {centre2[0] - radius * std::sin(angle), centre2[1] + radius * std::cos(angle)});
    normals.push_back({-std::sin(angle), std::cos(angle)});  // away from centre 2 (left)
  }
  centreline.push_back({centreline.back()[0] + lead, centreline.back()[1]});
  normals.push_back({0.0, 1.0});
  return BarAroundCentreline(centreline, width, 0.0, normals);
}

// Every segment is either excluded or covered exactly once; every corner has a feature.
void CheckPartition(const IdentificationInput &input, const IdentificationResult &result)
{
  REQUIRE(result.segments.size() == input.segments.size());
  double assigned = 0.0;
  for (std::size_t i = 0; i < input.segments.size(); i++)
  {
    const auto &segment = result.segments[i];
    const auto &source = input.segments[i];
    const double length =
        std::hypot(source.p1[0] - source.p0[0], source.p1[1] - source.p0[1]);
    if (segment.exclusion)
    {
      continue;
    }
    REQUIRE(!segment.portions.empty());
    double cursor = 0.0, covered = 0.0;
    for (const auto &portion : segment.portions)
    {
      CHECK(portion[0] >= cursor - 1.0e-9);
      CHECK(portion[1] > portion[0]);
      CHECK(portion[2] >= 0);
      CHECK(portion[2] < static_cast<double>(result.features.size()));
      covered += portion[1] - portion[0];
      cursor = portion[1];
    }
    CHECK_THAT(covered, WithinAbs(length, 1.0e-9));
    assigned += covered;
  }
  CHECK_THAT(assigned + result.excluded_length, WithinAbs(result.perimeter_length, 1.0e-9));
  CHECK_THAT(assigned, WithinAbs(result.assigned_length, 1.0e-9));
  double feature_length = 0.0;
  for (const auto &feature : result.features)
  {
    feature_length += feature.length;
  }
  CHECK_THAT(feature_length, WithinAbs(assigned, 1.0e-9));
  for (const auto &vertex : result.vertices)
  {
    if (vertex.type != "TruncationCut" && vertex.type != "PortCut")
    {
      CHECK(vertex.type != "Excluded");
      // A corner vertex absorbed by a bend (BendVertex) belongs to the bend's chain, whose
      // features it takes no part in: it has no feature of its own (no Feature in the
      // manifest); every other vertex maps to exactly one feature.
      if (vertex.type != "BendVertex")
      {
        CHECK(vertex.feature >= 0);
      }
    }
  }
}

// The loop under a rigid motion of the plan view: rotated by the angle about the origin,
// then shifted.
std::vector<Point2> RotatedAndShifted(const std::vector<Point2> &points,
                                      double angle_degrees, const Point2 &shift)
{
  const double pi = std::acos(-1.0);
  const double c = std::cos(angle_degrees * pi / 180.0),
               s = std::sin(angle_degrees * pi / 180.0);
  std::vector<Point2> out;
  for (const auto &p : points)
  {
    out.push_back({c * p[0] - s * p[1] + shift[0], s * p[0] + c * p[1] + shift[1]});
  }
  return out;
}

std::multiset<std::pair<std::string, std::string>>
FeatureHashes(const IdentificationResult &r)
{
  std::multiset<std::pair<std::string, std::string>> hashes;
  for (const auto &feature : r.features)
  {
    hashes.emplace(feature.type, feature.hash);
  }
  return hashes;
}

}  // namespace

TEST_CASE("SurfaceResponseIdentification", "[surfaceresponseidentification][Serial]")
{
  const double R = 2.0;

  // Two 12 x 8 pads with a 2 um gap: the gap ends form two mirror-image clusters with one
  // signature, the gap middle is a same-conductor gap claiming both edges, the outer
  // corners are plain corners, everything else isolated; the partition holds and does not
  // depend on the segment order or on the mesh subdivision.
  {
    const std::vector<LoopSpec> pads = {{Rectangle(-13.0, -6.0, -1.0, 6.0), 0, 1.0},
                                        {Rectangle(1.0, -6.0, 13.0, 6.0), 0, 1.0}};
    const auto input = MakeInput(pads, R);
    const auto result = IdentifyMetalPerimeter(input);
    CheckPartition(input, result);
    std::map<std::string, int> counts;
    for (const auto &feature : result.features)
    {
      counts[feature.type]++;
    }
    CHECK(counts["SpatialEdgeCluster"] == 2);
    CHECK(counts["SameConductorGap"] == 1);
    CHECK(counts["ConvexCorner"] == 4);
    CHECK(counts["IsolatedEdge"] == 6);
    std::set<std::string> cluster_hashes;
    for (const auto &feature : result.features)
    {
      if (feature.type == "SpatialEdgeCluster")
      {
        cluster_hashes.insert(feature.hash);
        CHECK(feature.vertices.size() == 2);
        CHECK(feature.signature["EdgeCount"] == 4);
      }
      if (feature.type == "SameConductorGap")
      {
        CHECK_THAT(feature.signature["SeparationOverR"].get<double>(),
                   WithinAbs(1.0, 1.0e-12));
        // Both partner edges are claimed (design (b) 2).
        std::set<int> chains;
        for (const auto &portion : feature.portions)
        {
          chains.insert(input.segments[portion.segment].chain);
        }
        CHECK(chains.size() == 2);
      }
    }
    CHECK(cluster_hashes.size() == 1);

    // Element ordering (A4) and mesh subdivision (A5) do not change the identification.
    const auto reversed = IdentifyMetalPerimeter(MakeInput(pads, R, true));
    CHECK(reversed.geometry_digest == result.geometry_digest);
    CHECK(FeatureHashes(reversed) == FeatureHashes(result));
    const std::vector<LoopSpec> refined = {{Rectangle(-13.0, -6.0, -1.0, 6.0), 0, 0.5},
                                           {Rectangle(1.0, -6.0, 13.0, 6.0), 0, 0.5}};
    const auto refined_result = IdentifyMetalPerimeter(MakeInput(refined, R));
    CheckPartition(MakeInput(refined, R), refined_result);
    CHECK(refined_result.geometry_digest == result.geometry_digest);
    CHECK(FeatureHashes(refined_result) == FeatureHashes(result));
    CHECK(refined_result.segments.size() == 2 * result.segments.size());
  }

  // Mirror images of an asymmetric cluster (a pad corner near the end of a perpendicular
  // strip, offset to one side) have one signature and opposite chirality; a rotated copy
  // has the same signature and the same chirality.
  {
    auto Scene = [&](double mirror, double rotate)
    {
      auto Transform = [&](std::vector<Point2> points)
      {
        for (auto &p : points)
        {
          p = {mirror * p[0], p[1]};
          const double c = std::cos(rotate), s = std::sin(rotate);
          p = {c * p[0] - s * p[1] + 100.0, s * p[0] + c * p[1] - 50.0};
        }
        if (mirror < 0.0)
        {
          std::reverse(points.begin(), points.end());  // keep counter-clockwise
        }
        return points;
      };
      return std::vector<LoopSpec>{{Transform(Rectangle(0.0, 0.0, 10.0, 10.0)), 0, 1.0},
                                   {Transform(Rectangle(-3.0, 3.0, -1.0, 20.0)), 0, 1.0}};
    };
    const auto base = IdentifyMetalPerimeter(MakeInput(Scene(1.0, 0.0), R));
    const auto mirrored = IdentifyMetalPerimeter(MakeInput(Scene(-1.0, 0.0), R));
    const auto rotated = IdentifyMetalPerimeter(MakeInput(Scene(1.0, 0.7), R));
    // hash -> chiralities of every cluster (the scene has a strip-end cluster at the top
    // and the asymmetric pad-corner cluster at the bottom).
    auto Clusters = [](const IdentificationResult &r)
    {
      std::map<std::string, std::multiset<int>> clusters;
      for (const auto &feature : r.features)
      {
        if (feature.type == "SpatialEdgeCluster")
        {
          clusters[feature.hash].insert(feature.chirality);
        }
      }
      return clusters;
    };
    const auto base_clusters = Clusters(base);
    const auto mirrored_clusters = Clusters(mirrored);
    const auto rotated_clusters = Clusters(rotated);
    REQUIRE(base_clusters.size() >= 2);
    CHECK(rotated_clusters == base_clusters);
    REQUIRE(mirrored_clusters.size() == base_clusters.size());
    bool asymmetric_found = false;
    for (const auto &[hash, chiralities] : base_clusters)
    {
      REQUIRE(mirrored_clusters.count(hash) == 1);
      std::multiset<int> flipped;
      for (const int chirality : mirrored_clusters.at(hash))
      {
        flipped.insert(-chirality);
      }
      CHECK(flipped == chiralities);
      if (chiralities.count(0) < chiralities.size())
      {
        asymmetric_found = true;  // a cluster that is not its own mirror image
      }
    }
    CHECK(asymmetric_found);
    CHECK(base.geometry_digest == mirrored.geometry_digest);
    CHECK(base.geometry_digest == rotated.geometry_digest);
    CheckPartition(MakeInput(Scene(1.0, 0.0), R), base);

    // Chip-scale coordinates: the same scene translated to (4999.123, 7321.789) um and
    // rotated by 37 deg carries coordinate roundoff ~ulp(1e4 um) ~ 1e-12 um into every
    // bisected portion endpoint; the signature grid (1e-6 R) must absorb it so that every
    // feature hash (clusters, corners, strips, isolated edges) is unchanged.
    {
      auto FarScene = [&](double rotate)
      {
        auto Transform = [&](std::vector<Point2> points)
        {
          for (auto &p : points)
          {
            const double c = std::cos(rotate), s = std::sin(rotate);
            p = {c * p[0] - s * p[1] + 4999.123, s * p[0] + c * p[1] + 7321.789};
          }
          return points;
        };
        return std::vector<LoopSpec>{{Transform(Rectangle(0.0, 0.0, 10.0, 10.0)), 0, 1.0},
                                     {Transform(Rectangle(-3.0, 3.0, -1.0, 20.0)), 0, 1.0}};
      };
      const auto far_input = MakeInput(FarScene(37.0 * std::acos(-1.0) / 180.0), R);
      const auto far = IdentifyMetalPerimeter(far_input);
      CheckPartition(far_input, far);
      CHECK(far.geometry_digest == base.geometry_digest);
      CHECK(FeatureHashes(far) == FeatureHashes(base));
      CHECK(Clusters(far) == base_clusters);
    }
  }

  // Two interaction regions merge into one cluster when their radius-R balls overlap (cores
  // closer than 2R) and stay separate otherwise: two thin fingers reaching toward a bar.
  {
    auto Fingers = [&](double spacing)
    {
      return std::vector<LoopSpec>{
          {Rectangle(-30.0, -10.0, 30.0, -3.0), 0, 1.0},
          {Rectangle(-0.5 * spacing - 1.0, 0.0, -0.5 * spacing + 1.0, 30.0), 0, 1.0},
          {Rectangle(0.5 * spacing - 1.0, 0.0, 0.5 * spacing + 1.0, 30.0), 0, 1.0}};
    };
    // Clusters that claim bar perimeter (the bar loop is chains 0-3 of the input).
    auto BarClusters = [](const IdentificationInput &input, const IdentificationResult &r)
    {
      std::vector<const IdentifiedFeature *> clusters;
      for (const auto &feature : r.features)
      {
        if (feature.type == "SpatialEdgeCluster" &&
            std::any_of(feature.portions.begin(), feature.portions.end(),
                        [&](const IdentifiedPortion &p)
                        { return input.segments[p.segment].chain < 4; }))
        {
          clusters.push_back(&feature);
        }
      }
      return clusters;
    };
    // Finger ends 3 um above the bar: each end interacts with the bar. With the fingers 3
    // um apart (inner edges) the two regions' cores overlap on the bar -> one cluster; 20
    // um apart -> two clusters with one signature (translated copies).
    const auto merged_input = MakeInput(Fingers(5.0), R);
    const auto separate_input = MakeInput(Fingers(22.0), R);
    const auto merged = IdentifyMetalPerimeter(merged_input);
    const auto separate = IdentifyMetalPerimeter(separate_input);
    CheckPartition(merged_input, merged);
    CheckPartition(separate_input, separate);
    REQUIRE(BarClusters(merged_input, merged).size() == 1);
    REQUIRE(BarClusters(separate_input, separate).size() == 2);
    CHECK(BarClusters(merged_input, merged).front()->vertices.size() == 4);
    CHECK(BarClusters(separate_input, separate)[0]->hash ==
          BarClusters(separate_input, separate)[1]->hash);
    CHECK(BarClusters(separate_input, separate)[0]->vertices.size() == 2);
  }

  // No knife edge at exactly R: a strip of width 0.975 R, R or 1.025 R reads the same way,
  // one strip pair claiming both long edges and one cluster per end (its two corners are
  // closer than 2R, so their windows overlap: invariant A2), both ends with one signature.
  {
    for (const double width : {1.95, 2.0, 2.05})
    {
      const auto input = MakeInput({{Rectangle(-12.0, 0.0, 12.0, width), 0, 1.0}}, R);
      const auto result = IdentifyMetalPerimeter(input);
      CheckPartition(input, result);
      std::map<std::string, int> counts;
      std::set<std::string> cluster_hashes;
      for (const auto &feature : result.features)
      {
        counts[feature.type]++;
        if (feature.type == "SpatialEdgeCluster")
        {
          cluster_hashes.insert(feature.hash);
          CHECK(feature.vertices.size() == 2);
        }
      }
      CHECK(counts["ConvexCorner"] == 0);
      CHECK(counts["SameConductorStrip"] == 1);
      CHECK(counts["SpatialEdgeCluster"] == 2);
      CHECK(cluster_hashes.size() == 1);
    }
    // Wider than 2R the long edges do not interact and the corners are plain corners.
    const auto wide_input = MakeInput({{Rectangle(-12.0, 0.0, 12.0, 4.5), 0, 1.0}}, R);
    const auto wide = IdentifyMetalPerimeter(wide_input);
    CheckPartition(wide_input, wide);
    CHECK(std::count_if(wide.features.begin(), wide.features.end(),
                        [](const auto &f) { return f.type == "ConvexCorner"; }) == 4);
    CHECK(std::none_of(wide.features.begin(), wide.features.end(),
                       [](const auto &f) { return f.type == "SpatialEdgeCluster"; }));
  }
}

TEST_CASE("SurfaceResponseIdentificationRoundedCorners",
          "[surfaceresponseidentification][Serial]")
{
  // 8 x 6 island with 0.5 um fillets (radius < R): four rounded convex corners separating
  // the four sides (a rounded corner is a corner: one isolated edge per side, as for the
  // sharp rectangle; decision 82(3)), no cluster; the same features when every chord is
  // bisected (refinement).
  const double R = 2.0;
  for (const double subdivision : {1.0, 0.5})
  {
    const auto input = MakeInput({{RoundedRectangle(4.0, 3.0, 0.5, 4), 0, subdivision}}, R);
    const auto result = IdentifyMetalPerimeter(input);
    CheckPartition(input, result);
    std::map<std::string, int> counts;
    for (const auto &feature : result.features)
    {
      counts[feature.type]++;
      if (feature.type == "ConvexCorner")
      {
        CHECK_THAT(feature.signature["CornerRadiusOverR"].get<double>(),
                   WithinAbs(0.25, 1.0e-3));
        CHECK_THAT(feature.signature["AngleDegrees"].get<double>(),
                   WithinAbs(90.0, 1.0e-6));
      }
    }
    CHECK(counts["ConvexCorner"] == 4);
    CHECK(counts["SpatialEdgeCluster"] == 0);
    CHECK(counts["IsolatedEdge"] == 4);
    CHECK(std::count_if(result.vertices.begin(), result.vertices.end(),
                        [](const auto &v) { return v.type == "RoundedCorner"; }) == 4);
  }
}

TEST_CASE("SurfaceResponseIdentificationObtuseCorners",
          "[surfaceresponseidentification][Serial]")
{
  // A trapezoid with 80 and 100 deg corners (taper-10): no arm pair comes within 2R outside
  // the corners' 2R zones (2a sin(40 deg) >= 5.1 R for a, b >= 2R), so all four corners are
  // plain corners and no cluster forms.
  const double R = 2.0;
  const double top = 12.0 - 20.0 * std::tan(10.0 * std::acos(-1.0) / 180.0);
  const std::vector<LoopSpec> taper = {
      {{{-12.0, 0.0}, {12.0, 0.0}, {top, 20.0}, {-top, 20.0}}, 0, 1.0}};
  const auto input = MakeInput(taper, R);
  const auto result = IdentifyMetalPerimeter(input);
  CheckPartition(input, result);
  std::map<std::string, int> counts;
  for (const auto &feature : result.features)
  {
    counts[feature.type]++;
    if (feature.type == "SpatialEdgeCluster")
    {
      for (const auto &portion : feature.portions)
      {
        const auto &s = input.segments[portion.segment];
        UNSCOPED_INFO("cluster portion chain "
                      << s.chain << " (" << s.p0[0] << "," << s.p0[1] << ")-(" << s.p1[0]
                      << "," << s.p1[1] << ") " << portion.s0 << ".." << portion.s1);
      }
    }
  }
  CHECK(counts["SpatialEdgeCluster"] == 0);
  CHECK(counts["ConvexCorner"] == 4);
}

// Whether the concentric offset polylines of ArcBar resolve their design circles under the
// sagitta arc rule (USER decision 117(4)): every chord's sagitta on the wider (outer) side
// below kArcSagittaOverRadius x R (the vertices lie on the circles exactly). A coarser
// polyline is corners: the meshed geometry.
// The concyclicity arc rule (USER decisions 121 / 122): the polyline bar is an arc iff
// every joint turns less than kArcMaxJointTurnDegrees (its vertices are concyclic by
// construction); the chord sagitta is the recorded mesh-coarseness diagnostic.
bool ArcBarIsArc(double sweep_degrees, double step_degrees)
{
  const int steps =
      std::max(1, static_cast<int>(std::lround(sweep_degrees / step_degrees)));
  // One chord has two joints: never an arc (a chamfer).
  return steps >= 2 && sweep_degrees / steps < kArcMaxJointTurnDegrees;
}

bool ArcBarCoarse(double radius, double width, double sweep_degrees, double step_degrees,
                  double R)
{
  const int steps =
      std::max(1, static_cast<int>(std::lround(sweep_degrees / step_degrees)));
  const double step = sweep_degrees * std::acos(-1.0) / 180.0 / steps;
  const double outer = radius + 0.5 * width;
  return outer * (1.0 - std::cos(0.5 * step)) >= kArcSagittaOverRadius * R;
}

TEST_CASE("SurfaceResponseIdentificationCurvedEdges",
          "[surfaceresponseidentification][Serial]")
{
  // Curved-edge chain rule (decision 73(1)) under the concyclicity arc rule (USER decisions
  // 121 / 122): a 3 um bar (1.5 R) along a polyline arc. Where every joint turns less than
  // ArcMaxJointTurnDegrees (the vertices are concyclic by construction) the polyline IS the
  // arc whatever its chord sagitta: the two sides are a pair along the bend (constant
  // separation), never a cluster; the arc is a curved pair when the inner bend radius is
  // below 10 R (decision 75) and a straight-like strip with a curvature annotation
  // otherwise; the leads are a plain strip; the two bar ends are one two-corner cluster
  // each. Chords coarser than the resolution (50 um at 22.5 deg: sagitta 0.5 R, 250 um at
  // 5 deg: 0.12 R) are arcs with the mesh-coarseness diagnostic set (MaxChordSagittaOverR
  // >= SagittaOverR), where the sagitta form of 117(4) read them as strings of corners. A
  // polygon whose joints turn 50 deg or more (a 180 deg sweep in three 60 deg chords) is
  // the meshed geometry: its joints are corners (within 2R of the facing side's corners:
  // clusters), no curved pair and no bend annotation. In both regimes the features do not
  // change under mesh refinement (A5).
  const double R = 2.0;
  struct Case
  {
    double radius, sweep, step;
    bool curved;
  };
  for (const Case &c : {Case{5.0, 90.0, 20.0, true}, Case{5.0, 90.0, 1.0, true},
                        Case{20.0, 90.0, 5.0, true}, Case{50.0, 45.0, 20.0, false},
                        Case{50.0, 45.0, 1.0, false}, Case{250.0, 15.0, 5.0, false},
                        Case{50.0, 180.0, 60.0, false}})
  {
    const auto bar = ArcBar(3.0, c.radius, c.sweep, c.step);
    const auto input = MakeInput({{bar, 0, 1.0}}, R);
    const auto result = IdentifyMetalPerimeter(input);
    const bool resolved = ArcBarIsArc(c.sweep, c.step);
    const bool coarse = ArcBarCoarse(c.radius, 3.0, c.sweep, c.step, R);
    INFO("radius " << c.radius << " step " << c.step << (resolved ? " (arc)" : " (polygon)")
                   << (coarse ? " coarse chords" : ""));
    CheckPartition(input, result);
    std::map<std::string, int> counts;
    int annotated = 0;
    for (const auto &feature : result.features)
    {
      counts[feature.type]++;
      annotated += feature.bend_radius_over_R.has_value();
      if (feature.type == "CurvedSameConductorStrip")
      {
        // The inner side's radius (radius - 1.5): the concentric polyline's exact circle.
        CHECK_THAT(feature.signature["RadiusOverR"].get<double>(),
                   WithinRel((c.radius - 1.5) / R, 1.0e-6));
        CHECK_THAT(feature.signature["SeparationOverR"].get<double>(),
                   WithinRel(1.5, 0.02));
      }
      if (feature.type == "SameConductorStrip" && resolved)
      {
        // (A polygon bar's parallel chords are 1.5 cos(step / 2) R apart: the meshed
        // geometry, not checked.)
        CHECK_THAT(feature.signature["SeparationOverR"].get<double>(),
                   WithinRel(1.5, 0.02));
        if (!c.curved)
        {
          REQUIRE(feature.bend_radius_over_R.has_value());
          CHECK_THAT(*feature.bend_radius_over_R, WithinRel((c.radius - 1.5) / R, 0.03));
        }
      }
    }
    if (resolved)
    {
      CHECK(counts["SpatialEdgeCluster"] == 2);
      CHECK(counts["SameConductorStrip"] >= 1);
      CHECK(counts["CurvedSameConductorStrip"] == (c.curved ? 1 : 0));
      CHECK(counts["IsolatedEdge"] == 0);
      CHECK(counts["CurvedEdge"] == 0);
      CHECK(counts["ConvexCorner"] == 0);
      // The mesh-coarseness diagnostic: the recorded largest chord sagitta of the bar's
      // arcs reaches SagittaOverR exactly when the chords are coarser than the resolution.
      REQUIRE(!result.arcs.empty());
      const double worst =
          std::max_element(result.arcs.begin(), result.arcs.end(),
                           [](const auto &a, const auto &b)
                           { return a.max_sagitta_over_R < b.max_sagitta_over_R; })
              ->max_sagitta_over_R;
      CHECK((worst >= kArcSagittaOverRadius) == coarse);
    }
    else
    {
      // Corners at the polyline joints (joints of 50 deg or more are never arc joints): the
      // polygon is the meshed geometry. The joints of the two sides face each other within
      // 2R (1.5 R apart), so the corner sites form clusters; the vertex table lists every
      // joint as a corner vertex, more than the four bar-end corners, and no bend vertex.
      CHECK(counts["CurvedSameConductorStrip"] == 0);
      CHECK(counts["CurvedEdge"] == 0);
      CHECK(annotated == 0);
      CHECK(counts["ConvexCorner"] + counts["ConcaveCorner"] +
                counts["SpatialEdgeCluster"] >=
            1);
      CHECK(std::count_if(
                result.vertices.begin(), result.vertices.end(), [](const auto &v)
                { return v.type == "ConvexCorner" || v.type == "ConcaveCorner"; }) > 4);
      CHECK(std::count_if(result.vertices.begin(), result.vertices.end(),
                          [](const auto &v) { return v.type == "BendVertex"; }) == 0);
    }
    // Refinement: every chord bisected twice.
    const auto refined_input = MakeInput({{bar, 0, 0.25}}, R);
    const auto refined = IdentifyMetalPerimeter(refined_input);
    CheckPartition(refined_input, refined);
    CHECK(refined.geometry_digest == result.geometry_digest);
    CHECK(FeatureHashes(refined) == FeatureHashes(result));
  }

  // A lone polyline arc edge (no partner within 2R): a 3 um wide bar would pair, so use a
  // wide bar (6 um = 3 R): the inner and outer sides are isolated; the arc part is a
  // CurvedEdge with the side's own bend radius, the leads and the outer straight-like parts
  // are isolated edges.
  {
    const auto bar = ArcBar(6.0, 8.0, 90.0, 5.0);
    const auto input = MakeInput({{bar, 0, 1.0}}, R);
    const auto result = IdentifyMetalPerimeter(input);
    CheckPartition(input, result);
    std::map<std::string, int> counts;
    std::vector<double> radii;
    for (const auto &feature : result.features)
    {
      counts[feature.type]++;
      if (feature.type == "CurvedEdge")
      {
        radii.push_back(feature.signature["RadiusOverR"].get<double>());
      }
    }
    CHECK(counts["SpatialEdgeCluster"] == 0);
    CHECK(counts["ConvexCorner"] == 4);
    CHECK(counts["CurvedEdge"] == 2);
    CHECK(counts["IsolatedEdge"] >= 2);
    std::sort(radii.begin(), radii.end());
    REQUIRE(radii.size() == 2);
    CHECK_THAT(radii[0], WithinRel(5.0 / R, 0.03));   // inner side, radius 8 - 3
    CHECK_THAT(radii[1], WithinRel(11.0 / R, 0.03));  // outer side, radius 8 + 3
  }
}

TEST_CASE("SurfaceResponseIdentificationConvexity",
          "[surfaceresponseidentification][Serial]")
{
  // Convexity (decision 108 / C4): the signed windowed curvature toward the metal. A wide
  // bar (3 R) along a 90 deg bend: the inner side's edge bends around the gap (Concave: a
  // hole-like edge), the outer side's around the metal (Convex: a disk-like edge); every
  // portion carries its signed turn, whose sum over a CurvedEdge is the bend's sweep in the
  // metal-ward sense (within the window's smoothing across the section ends); a mirrored
  // scene keeps the convexities (a mirror does not exchange inside and outside).
  const double R = 2.0;
  const double quarter = 0.5 * std::acos(-1.0);
  auto Convexities = [&](const std::vector<Point2> &loop)
  {
    const auto input = MakeInput({{loop, 0, 1.0}}, R);
    const auto result = IdentifyMetalPerimeter(input);
    CheckPartition(input, result);
    std::map<double, std::pair<std::string, double>>
        by_radius;  // RadiusOverR -> (convexity, turn)
    for (const auto &feature : result.features)
    {
      if (feature.type == "CurvedEdge")
      {
        double turn = 0.0;
        for (const auto &portion : feature.portions)
        {
          turn += portion.turn;
        }
        by_radius[feature.signature["RadiusOverR"].get<double>()] = {
            feature.signature["Convexity"].get<std::string>(), turn};
      }
      else
      {
        CHECK_FALSE(feature.signature.contains("Convexity"));
      }
    }
    return by_radius;
  };
  {
    const auto by_radius = Convexities(ArcBar(6.0, 8.0, 90.0, 5.0));
    REQUIRE(by_radius.size() == 2);
    const auto inner = by_radius.begin(), outer = std::next(inner);
    CHECK(inner->second.first == "Concave");
    CHECK(outer->second.first == "Convex");
    CHECK_THAT(inner->second.second, WithinRel(-quarter, 0.10));
    CHECK_THAT(outer->second.second, WithinRel(quarter, 0.10));
    // Mirrored in y (the bend turns right): the same classes and convexities.
    auto mirrored = ArcBar(6.0, 8.0, 90.0, 5.0);
    for (auto &p : mirrored)
    {
      p[1] = -p[1];
    }
    std::reverse(mirrored.begin(), mirrored.end());  // keep the loop counter-clockwise
    const auto mirrored_by_radius = Convexities(mirrored);
    REQUIRE(mirrored_by_radius.size() == 2);
    CHECK(mirrored_by_radius.begin()->second.first == "Concave");
    CHECK(std::next(mirrored_by_radius.begin())->second.first == "Convex");
    CHECK_THAT(mirrored_by_radius.begin()->second.second, WithinRel(-quarter, 0.10));
  }
  // Straight-like bend (radius 50 = 25 R): no CurvedEdge; the isolated edges carry the
  // bend annotation and their portions the signed turn (the inner side around the gap).
  {
    const auto input = MakeInput({{ArcBar(6.0, 50.0, 45.0, 1.0), 0, 1.0}}, R);
    const auto result = IdentifyMetalPerimeter(input);
    int annotated = 0;
    for (const auto &feature : result.features)
    {
      CHECK(feature.type != "CurvedEdge");
      if (feature.type == "IsolatedEdge" && feature.bend_radius_over_R)
      {
        double turn = 0.0;
        for (const auto &portion : feature.portions)
        {
          turn += portion.turn;
        }
        if (std::abs(turn) > 0.1)
        {
          annotated++;
          const bool inner = *feature.bend_radius_over_R < 25.0;
          CHECK_THAT(turn, WithinRel((inner ? -1.0 : 1.0) * 0.5 * quarter, 0.10));
        }
      }
    }
    CHECK(annotated == 2);
  }
  // A curved strip (1.5 R wide along a 5 um bend) records the convexity of its first side
  // (the edge the coupon's first edge e1 lands on). Its two edges are alike (chirality 0),
  // so the first side is the outermost edge of the bend (decision 214 (i)), whose metal
  // lies inside the bend: a symmetric curved strip is Convex, in every frame.
  {
    const auto input = MakeInput({{ArcBar(3.0, 5.0, 90.0, 5.0), 0, 1.0}}, R);
    const auto result = IdentifyMetalPerimeter(input);
    int curved_strips = 0;
    for (const auto &feature : result.features)
    {
      if (feature.type != "CurvedSameConductorStrip")
      {
        continue;
      }
      curved_strips++;
      const std::string convexity = feature.signature["Convexity"].get<std::string>();
      // The first side (side 0 for chirality >= 0, the last for -1): inner when its
      // portions lie closer to the bend centre (0, 5) than the other side's.
      const int first_side = feature.chirality < 0 ? 1 : 0;
      double first_radius = 0.0, other_radius = 0.0;
      int first_count = 0, other_count = 0;
      for (const auto &portion : feature.portions)
      {
        const auto &p0 = input.segments[portion.segment].p0;
        const double r = std::hypot(p0[0], p0[1] - 5.0);
        (portion.side == first_side ? first_radius : other_radius) += r;
        (portion.side == first_side ? first_count : other_count)++;
      }
      REQUIRE(first_count > 0);
      REQUIRE(other_count > 0);
      CHECK(feature.chirality == 0);
      CHECK(first_radius / first_count > other_radius / other_count + 1.0);
      CHECK(convexity == "Convex");
    }
    CHECK(curved_strips == 1);
  }
  // An S-bend of two 40 deg arcs of radius 6 (3 R): one curved section with both senses in
  // the curved regime -> Convexity "Mixed" (reported, never one convexity).
  {
    const auto by_radius = Convexities(SBar(6.0, 6.0, 40.0, 4.0));
    REQUIRE(!by_radius.empty());
    int mixed = 0;
    for (const auto &[radius, entry] : by_radius)
    {
      (void)radius;
      mixed += entry.first == "Mixed";
    }
    CHECK(mixed >= 1);
  }
}

TEST_CASE("SurfaceResponseIdentificationPairsAtTheThreshold",
          "[surfaceresponseidentification][Serial]")
{
  // A gap of exactly 2R (and 2R +/- 1e-3 R) between two 8 um bars (4 R: the far corners of
  // a bar end are beyond the vertex-join reach of the near corners) that are concentric
  // offsets of one centreline along bends of 50 and 250 um at three discretisations. Where
  // the polylines resolve their circles (sagitta below 0.05 R: 1 deg per vertex, 5 deg on
  // the 50 um bend) the interaction decision uses the separation of the fitted arcs — the
  // radius difference, exact — so exactly 2R and 2R + 1e-3 R are isolated edges at every
  // such discretisation (no cluster, no pair: the mid-chord dips of the inscribed chords
  // below 2R, DS-SCT-001's 4 um gaps at 3.9998 that became 3 mm clusters, create no events)
  // and 2R - 1e-3 R is one DifferentConductorGap along the whole pair (straight-like: the
  // bends are 21 R and 121 R, with the bend annotation). The corner pairs across a
  // 2R - 1e-3 R gap are events (two clusters) in every case. Under the concyclicity rule
  // (USER decisions 121 / 122) every discretisation here (1 / 5 / 15 deg per vertex, all
  // below ArcMaxJointTurnDegrees) is an arc, the coarse ones (15 deg: sagitta 0.23 R on the
  // 50 um bend; 5 deg on the 250 um bend: 0.12 R) with the mesh-coarseness diagnostic set;
  // a polygon of 60 deg joints (a 180 deg sweep in three chords) is corners at its joints:
  // no curved pair, no bend annotation, corner vertices at every joint.
  const double R = 2.0, width = 8.0;
  for (const double radius : {50.0, 250.0})
  {
    for (const double step : {1.0, 5.0, 15.0, 60.0})
    {
      // (The 250 um bend sweeps 30 deg so that the 15 deg discretisation has two chords:
      // one chord's two 7.5 deg end joints on 6 um leads imply 0.049 R and are noise, a
      // straight-like pair, not an arc test.)
      const double sweep = step >= 60.0 ? 180.0 : (radius < 100.0 ? 45.0 : 30.0);
      for (const double gap : {2.0 * R, 2.0 * R - 1.0e-3 * R, 2.0 * R + 1.0e-3 * R})
      {
        // Both bars are concentric offsets of the gap's centreline (radius), so the facing
        // edges are inscribed polylines of two circles gap apart at the same angles.
        const auto inner = ArcBar(width, radius, sweep, step, 6.0, 0.5 * gap + 0.5 * width);
        const auto outer =
            ArcBar(width, radius, sweep, step, 6.0, -0.5 * gap - 0.5 * width);
        const auto input = MakeInput({{inner, 0, 1.0}, {outer, 1, 1.0}}, R);
        const auto result = IdentifyMetalPerimeter(input);
        const bool resolved = ArcBarIsArc(sweep, step);
        const bool coarse = ArcBarCoarse(radius, 2.0 * (0.5 * gap + width), sweep, step, R);
        INFO("radius " << radius << " step " << step << " gap " << gap
                       << (resolved ? " (arc)" : " (polygon)")
                       << (coarse ? " coarse chords" : ""));
        CheckPartition(input, result);
        std::map<std::string, int> counts;
        int annotated = 0;
        for (const auto &feature : result.features)
        {
          counts[feature.type]++;
          annotated += feature.bend_radius_over_R.has_value();
          if (feature.type == "DifferentConductorGap" && resolved)
          {
            // Decision 85(1): the separation along the bends is the radius difference of
            // the two fitted arcs (exact for polylines inscribed in the design circles).
            CHECK_THAT(feature.signature["SeparationOverR"].get<double>(),
                       WithinAbs(gap / R, kSignatureParameterToleranceOverRadius));
          }
        }
        const bool corner_events = gap < 2.0 * R;
        CHECK(counts["CurvedEdge"] == 0);
        CHECK(counts["CurvedDifferentConductorGap"] == 0);
        if (resolved)
        {
          const bool pair = gap < 2.0 * R;
          CHECK(counts["DifferentConductorGap"] == (pair ? 1 : 0));
          CHECK(counts["SpatialEdgeCluster"] == (corner_events ? 2 : 0));
          CHECK(counts["ConvexCorner"] == (corner_events ? 4 : 8));
          CHECK(counts["IsolatedEdge"] >= (pair ? 2 : 4));
          REQUIRE(!result.arcs.empty());
          const double worst =
              std::max_element(result.arcs.begin(), result.arcs.end(),
                               [](const auto &a, const auto &b)
                               { return a.max_sagitta_over_R < b.max_sagitta_over_R; })
                  ->max_sagitta_over_R;
          CHECK((worst >= kArcSagittaOverRadius) == coarse);
        }
        else
        {
          CHECK(annotated == 0);
          CHECK(std::count_if(
                    result.vertices.begin(), result.vertices.end(), [](const auto &v)
                    { return v.type == "ConvexCorner" || v.type == "ConcaveCorner"; }) > 8);
        }
      }
    }
  }
}

TEST_CASE("SurfaceResponseIdentificationOneInteractionDistance",
          "[surfaceresponseidentification][Serial]")
{
  // Decision 82(1): one interaction distance (3D, strictly below 2R) and no feature across
  // two metal planes.
  const double R = 2.0;
  const std::vector<Point2> pad = Rectangle(-10.0, -6.0, 10.0, 6.0);
  const auto single = IdentifyMetalPerimeter(MakeInput({{pad, 0, 1.0}}, R));
  std::map<std::string, int> single_counts;
  for (const auto &feature : single.features)
  {
    single_counts[feature.type]++;
  }
  CHECK(single_counts["ConvexCorner"] == 4);
  CHECK(single_counts["IsolatedEdge"] == 4);

  // Identical pads on planes 2R and 2.4R apart (and a laterally offset one at R): no pair,
  // no cluster, no vertex join between the planes — every feature lies on one plane and
  // the feature multiset is twice the single pad's.
  for (const double z : {2.0 * R, 2.4 * R})
  {
    for (const double shift : {0.0, R})
    {
      std::vector<Point2> upper = pad;
      for (auto &p : upper)
      {
        p[0] += shift;
      }
      const auto input = MakeInput({{pad, 0, 1.0}, {upper, 1, 1.0, z}}, R);
      const auto result = IdentifyMetalPerimeter(input);
      INFO("z " << z << " shift " << shift);
      CheckPartition(input, result);
      std::map<std::string, int> counts;
      for (const auto &feature : result.features)
      {
        counts[feature.type]++;
        std::set<long long> planes;
        for (const auto &portion : feature.portions)
        {
          planes.insert(std::llround(input.segments[portion.segment].p0[2] * 1.0e6));
        }
        CHECK(planes.size() <= 1);
      }
      CHECK(counts["ConvexCorner"] == 8);
      CHECK(counts["IsolatedEdge"] == 8);
      CHECK(counts.size() == 2);
    }
  }

  // Two 3 um slots (corner clusters at both slot ends) on planes 2.4R apart, the upper one
  // shifted by 2 um: its sites are 2.6R from the lower cores (inside the former 3R vertex
  // join) — the clusters of the two planes stay separate, its slot edges above the lower
  // slot's do not pair.
  {
    auto Slot = [](double dy)
    {
      const double hw = 1.5, x_out = 15.5, y_top = 6.0 + dy, y_slot = -6.0 + dy,
                   y_bottom = -12.0 + dy;
      return std::vector<Point2>{{-x_out, y_bottom}, {x_out, y_bottom}, {x_out, y_top},
                                 {hw, y_top},        {hw, y_slot},      {-hw, y_slot},
                                 {-hw, y_top},       {-x_out, y_top}};
    };
    const auto lower_input = MakeInput({{Slot(0.0), 0, 1.0}}, R);
    const auto lower = IdentifyMetalPerimeter(lower_input);
    const auto input = MakeInput({{Slot(0.0), 0, 1.0}, {Slot(2.0), 1, 1.0, 2.4 * R}}, R);
    const auto result = IdentifyMetalPerimeter(input);
    CheckPartition(input, result);
    std::map<std::string, int> counts, lower_counts;
    for (const auto &feature : lower.features)
    {
      lower_counts[feature.type]++;
    }
    for (const auto &feature : result.features)
    {
      counts[feature.type]++;
    }
    CHECK(lower_counts["SpatialEdgeCluster"] == 2);
    CHECK(lower_counts["SameConductorGap"] == 1);
    for (const auto &[type, count] : lower_counts)
    {
      CHECK(counts[type] == 2 * count);
    }
  }
}

TEST_CASE("SurfaceResponseIdentificationPortCutAndBroadcast",
          "[surfaceresponseidentification][Serial]")
{
  // Decision 82(5): the metal perimeter bordering a port face is a Port exclusion and its
  // vertices are PortCut (no corner / endpoint feature); the serialised result round-trips.
  const double R = 2.0;
  auto input = MakeInput({{Rectangle(-20.0, -3.0, -2.0, 3.0), 0, 1.0},
                          {Rectangle(2.0, -3.0, 20.0, 3.0), 1, 1.0}},
                         R);
  // The lead ends at x = +-2 border the port face: excluded, their vertices port cuts.
  for (std::size_t s = 0; s < input.segments.size(); s++)
  {
    auto &segment = input.segments[s];
    if (std::abs(std::abs(segment.p0[0]) - 2.0) < 1.0e-9 &&
        std::abs(std::abs(segment.p1[0]) - 2.0) < 1.0e-9)
    {
      segment.exclusion = std::make_pair("Port", "metal perimeter bordering a port face");
      for (const std::size_t v : segment.vertices)
      {
        input.vertices[v].on_port_boundary = true;
        if (std::abs(std::abs(input.vertices[v].coordinate[1]) - 3.0) < 1.0e-9)
        {
          input.vertices[v].physical_type = MetalEdgeVertexType::ENDPOINT;
        }
      }
    }
  }
  const auto result = IdentifyMetalPerimeter(input);
  CheckPartition(input, result);
  std::map<std::string, int> counts, vertex_types;
  for (const auto &feature : result.features)
  {
    counts[feature.type]++;
  }
  for (const auto &vertex : result.vertices)
  {
    vertex_types[vertex.type]++;
  }
  for (const auto &[type, count] : counts)
  {
    INFO(type << " " << count);
    CHECK((type == "ConvexCorner" || type == "IsolatedEdge"));
  }
  CHECK(counts["ConvexCorner"] == 4);
  CHECK(counts["IsolatedEdge"] == 6);
  CHECK(counts["Endpoint"] == 0);
  CHECK(vertex_types["PortCut"] == 4);
  CHECK(vertex_types["ConvexCorner"] == 4);
  double port_length = 0.0;
  for (const auto &exclusion : result.exclusions)
  {
    if (exclusion.cls == "Port")
    {
      port_length += exclusion.length;
    }
  }
  CHECK_THAT(port_length, WithinAbs(12.0, 1.0e-9));

  // Broadcast form: the deserialised result produces the same manifest object.
  const auto copy = DeserializeIdentificationResult(SerializeIdentificationResult(result));
  CHECK(copy.ToJson(1.0) == result.ToJson(1.0));
  CHECK(copy.geometry_digest == result.geometry_digest);
}

TEST_CASE("SurfaceResponseIdentificationStacks", "[surfaceresponseidentification][Serial]")
{
  // Decision 82(2): a translation-invariant cross-section of k >= 3 edges whose consecutive
  // separations are below 2R is ONE feature (ParallelEdgeCluster, CurvedParallelEdgeCluster
  // along a bend), straight and curved; the pairwise candidates inside it are superseded
  // (no claim of one priority ever overlaps another: Diagnostics); a member taken by a
  // cluster is recomposed out of the cross-section at the stack ends; sides in the
  // canonical order.
  const double R = 2.0;
  auto Offsets = [](const IdentifiedFeature &f)
  {
    std::vector<double> offsets;
    for (const auto &edge : f.signature["Edges"])
    {
      offsets.push_back(edge["OffsetOverR"].get<double>());
    }
    return offsets;
  };
  // Either orientation of the stated offsets may be the canonical one.
  auto Matches = [](std::vector<double> found, std::vector<double> wanted)
  {
    std::sort(found.begin(), found.end());
    std::sort(wanted.begin(), wanted.end());
    if (found.size() != wanted.size())
    {
      return false;
    }
    std::vector<double> mirrored;
    for (const double w : wanted)
    {
      mirrored.push_back(wanted.back() - w);
    }
    std::sort(mirrored.begin(), mirrored.end());
    for (const auto *candidate : {&wanted, &mirrored})
    {
      bool ok = true;
      for (std::size_t i = 0; i < found.size(); i++)
      {
        ok = ok && std::abs(found[i] - (*candidate)[i]) <= 0.02;
      }
      if (ok)
      {
        return true;
      }
    }
    return false;
  };
  SECTION("straight 4-edge stack: ground | 2 | trace 2 | 2 | ground")
  {
    // The bars end inside the domain (end corners): one corner cluster per end. Both
    // grounds are cluster material next to the trace's end edge (events within 2R of it)
    // over a longer reach than the trace edges, exactly R from the cores: the trace strip
    // is recomposed there (the stack-end rule).
    const auto input = MakeInput({{Rectangle(-40.0, -10.0, 40.0, -2.0), 0, 1.0},
                                  {Rectangle(-30.0, 0.0, 30.0, 2.0), 0, 1.0},
                                  {Rectangle(-40.0, 4.0, 40.0, 12.0), 0, 1.0}},
                                 R);
    const auto result = IdentifyMetalPerimeter(input);
    CheckPartition(input, result);
    CHECK(result.same_priority_claim_overlaps == 0);
    std::map<std::string, int> counts;
    const IdentifiedFeature *stack = nullptr;
    for (const auto &feature : result.features)
    {
      counts[feature.type]++;
      if (feature.type == "ParallelEdgeCluster")
      {
        stack = &feature;
      }
    }
    CHECK(counts["ParallelEdgeCluster"] == 1);
    CHECK(counts["SameConductorStrip"] == 1);
    CHECK(counts["SpatialEdgeCluster"] == 2);
    REQUIRE(stack != nullptr);
    CHECK(Matches(Offsets(*stack), {0.0, 1.0, 2.0, 3.0}));
    CHECK(stack->signature["Edges"].size() == 4);
    // Every side present with the same length (a symmetric cross-section), all four sides
    // over the same longitudinal extent: the ground cores within 2R of the trace end edge
    // reach sqrt((2R)^2 - (2 um)^2) = sqrt(12) um past the trace end, their R balls another
    // R, so |x| < 30 - sqrt(12) - R um.
    std::map<int, double> side_length;
    for (const auto &portion : stack->portions)
    {
      side_length[portion.side] += portion.s1 - portion.s0;
    }
    REQUIRE(side_length.size() == 4);
    for (const auto &[side, length] : side_length)
    {
      INFO("side " << side);
      CHECK_THAT(length, WithinAbs(2.0 * (30.0 - std::sqrt(12.0) - R), 0.05));
    }
    CHECK(stack->chirality != -1);
  }
  SECTION("asymmetric straight stack: sides in the canonical order")
  {
    // ground | 1 | trace 1.5 | 3 | ground: offsets 0 / 0.5 / 1.25 / 2.75 R; chirality +1
    // and one side index per physical edge over the whole feature.
    const auto input = MakeInput({{Rectangle(-40.0, -9.0, 40.0, -1.0), 0, 1.0},
                                  {Rectangle(-30.0, 0.0, 30.0, 1.5), 0, 1.0},
                                  {Rectangle(-40.0, 4.5, 40.0, 12.5), 0, 1.0}},
                                 R);
    const auto result = IdentifyMetalPerimeter(input);
    CheckPartition(input, result);
    CHECK(result.same_priority_claim_overlaps == 0);
    int stacks = 0;
    for (const auto &feature : result.features)
    {
      if (feature.type != "ParallelEdgeCluster" || feature.signature["Edges"].size() != 4)
      {
        continue;
      }
      stacks++;
      CHECK(Matches(Offsets(feature), {0.0, 0.5, 1.25, 2.75}));
      CHECK(feature.chirality == 1);
      // The side of a portion is a function of its edge (the segment's y coordinate).
      std::map<int, std::set<int>> sides_of_y;
      for (const auto &portion : feature.portions)
      {
        const double y = input.segments[portion.segment].p0[1];
        sides_of_y[static_cast<int>(std::lround(10.0 * y))].insert(portion.side);
      }
      CHECK(sides_of_y.size() == 4);
      for (const auto &[y, sides] : sides_of_y)
      {
        INFO("y / 10 " << y);
        CHECK(sides.size() == 1);
      }
    }
    CHECK(stacks == 1);
  }
  SECTION("curved 3-edge stack along a 3R bend")
  {
    // A 3 um trace (inner radius 3R) and an 8 um ground band 2 um outside it along a 90 deg
    // bend of the centreline with 12 um leads: a straight ParallelEdgeCluster on the leads
    // and a CurvedParallelEdgeCluster along the bend with RadiusOverR = the trace's inner
    // radius.
    const double centre = 3.0 * R + 1.5;
    const auto trace = ArcBar(3.0, centre, 90.0, 5.0, 12.0, 0.0);
    const auto ground = ArcBar(8.0, centre, 90.0, 5.0, 12.0, -(1.5 + 2.0 + 4.0));
    const auto input = MakeInput({{trace, 0, 1.0}, {ground, 0, 1.0}}, R);
    const auto result = IdentifyMetalPerimeter(input);
    CheckPartition(input, result);
    CHECK(result.same_priority_claim_overlaps == 0);
    std::map<std::string, int> counts;
    for (const auto &feature : result.features)
    {
      counts[feature.type]++;
      if (feature.type == "ParallelEdgeCluster" ||
          feature.type == "CurvedParallelEdgeCluster")
      {
        CHECK(Matches(Offsets(feature), {0.0, 1.5, 2.5}));
      }
      if (feature.type == "CurvedParallelEdgeCluster")
      {
        CHECK_THAT(feature.signature["RadiusOverR"].get<double>(), WithinAbs(3.0, 0.15));
      }
    }
    CHECK(counts["ParallelEdgeCluster"] == 1);
    CHECK(counts["CurvedParallelEdgeCluster"] == 1);
    CHECK(counts["SpatialEdgeCluster"] == 2);
    CHECK(counts["CurvedSameConductorGap"] == 0);
    CHECK(counts["SameConductorGap"] == 0);
  }
  SECTION("curved 3-edge stack: convexity of every side related through the GapSides")
  {
    // The convexity of a curved stack is the signature's first edge's; AssembleStack reads
    // it from the first side, else from the far side (opposite), else from the first
    // interior side carrying curvature (review J/C4 m5): concentric sides bend in one
    // geometric sense, so the convexity of edge k equals the first edge's when their
    // GapSides agree and is the opposite otherwise. Asserted on every side of the 3-edge
    // stack (the ground band outside and inside the bend: the first edge Convex in one
    // scene and Concave in the other), reading each side's own convexity from its portions'
    // signed turns (toward the metal positive). A stack whose outer sides are exactly
    // straight while an interior one bends cannot be composed (a bent side against a
    // straight partner leaves the pair constancy band where its windowed bend radius drops
    // below 10R: the two rules meet at 2R x 5 %), so the interior-side branch itself is
    // unreachable through IdentifyMetalPerimeter and the relation it encodes is what this
    // test pins. Ground band outside the bend (trace centreline radius 3R + 1.5) or inside
    // it (radius 3R + 1.5 + 2 + 8: the band's far edge at 3R).
    std::set<std::string> first_edge_convexities;
    for (const double ground_offset : {-(1.5 + 2.0 + 4.0), 1.5 + 2.0 + 4.0})
    {
      const double centre = 3.0 * R + 1.5 + (ground_offset > 0.0 ? 2.0 + 8.0 : 0.0);
      const auto trace = ArcBar(3.0, centre, 90.0, 5.0, 12.0, 0.0);
      const auto ground = ArcBar(8.0, centre, 90.0, 5.0, 12.0, ground_offset);
      const auto input = MakeInput({{trace, 0, 1.0}, {ground, 0, 1.0}}, R);
      const auto result = IdentifyMetalPerimeter(input);
      CheckPartition(input, result);
      int curved_stacks = 0;
      for (const auto &feature : result.features)
      {
        if (feature.type != "CurvedParallelEdgeCluster")
        {
          continue;
        }
        curved_stacks++;
        const auto &edges = feature.signature["Edges"];
        REQUIRE(edges.size() == 3);
        const std::string convexity = feature.signature["Convexity"].get<std::string>();
        REQUIRE((convexity == "Convex" || convexity == "Concave"));
        first_edge_convexities.insert(convexity);
        std::map<int, double> side_turn;
        for (const auto &portion : feature.portions)
        {
          side_turn[portion.side] += portion.turn;
        }
        REQUIRE(side_turn.size() == 3);
        const int first_gap = edges.at(0)["GapSide"].get<int>();
        for (int k = 0; k < 3; k++)
        {
          // Signature edge k lies on side k for chirality >= 0 and on side 2 - k for -1.
          const int side = feature.chirality < 0 ? 2 - k : k;
          const double turn = side_turn.at(side);
          INFO("ground offset " << ground_offset << " edge " << k << " side " << side
                                << " turn " << turn);
          REQUIRE(std::abs(turn) > 0.1);
          const std::string own = turn > 0.0 ? "Convex" : "Concave";
          const bool same_gap = edges.at(k)["GapSide"].get<int>() * first_gap > 0;
          const std::string related =
              same_gap ? own : (own == "Convex" ? "Concave" : "Convex");
          CHECK(related == convexity);
        }
      }
      CHECK(curved_stacks == 1);
    }
    CHECK(first_edge_convexities == std::set<std::string>{"Concave", "Convex"});
  }
}

TEST_CASE("SurfaceResponseIdentificationExactParametersAndExtension",
          "[surfaceresponseidentification][Serial]")
{
  // Decision 85 (2026-09-26). (1) Exact signature parameters: the curved 3-edge stack of
  // the previous test at two discretisations of the bend (5 and 2.5 deg steps) gives
  // IDENTICAL offsets (0 / 1.5 / 2.5 R exactly: the arc radius differences) and bend
  // radius, hence one signature key; the tolerance API groups near-identical instances and
  // tells the mirror orientation apart from a different topology. (2) Cluster extension:
  // the two 12 x 8 pads of the first test have every single-edge portion within 2R of a
  // cluster's claimed perimeter absorbed (the gap edges between the end clusters are the
  // pair, the remainder isolated only where nothing is within 2R across).
  const double R = 2.0;
  // A band between two concentric circles (vertices ON the circles: the inscribed
  // construction of a CAD polygonisation) with tangent leads, counter-clockwise.
  auto InscribedBand =
      [](double r_in, double r_out, double sweep_degrees, double step_degrees, double lead)
  {
    const int steps = static_cast<int>(std::lround(sweep_degrees / step_degrees));
    const double sweep = sweep_degrees * std::acos(-1.0) / 180.0;
    std::vector<Point2> inner, outer;
    for (int k = 0; k <= steps; k++)
    {
      const double phi = sweep * k / steps;
      inner.push_back({r_in * std::cos(phi), r_in * std::sin(phi)});
      outer.push_back({r_out * std::cos(phi), r_out * std::sin(phi)});
    }
    // Leads along the tangents at phi = 0 (direction -y) and phi = sweep.
    const Point2 t0 = {0.0, -1.0}, t1 = {-std::sin(sweep), std::cos(sweep)};
    std::vector<Point2> loop = {{r_in, -lead}};
    loop.insert(loop.end(), inner.begin(), inner.end());
    loop.push_back({inner.back()[0] + lead * t1[0], inner.back()[1] + lead * t1[1]});
    loop.push_back({outer.back()[0] + lead * t1[0], outer.back()[1] + lead * t1[1]});
    loop.insert(loop.end(), outer.rbegin(), outer.rend());
    loop.push_back({r_out, -lead});
    (void)t0;
    return loop;
  };
  SECTION("exact offsets at two discretisations (inscribed bands)")
  {
    // Trace 3 um from 3R, gap 2 um, ground 8 um: offsets 0 / 1.5 / 2.5 R exactly at 5 and
    // 2.5 deg steps, bend radius exactly 3R, identical keys.
    std::vector<std::string> keys;
    for (const double step : {5.0, 2.5})
    {
      const double r0 = 3.0 * R;
      const auto trace = InscribedBand(r0, r0 + 3.0, 90.0, step, 12.0);
      const auto ground = InscribedBand(r0 + 5.0, r0 + 13.0, 90.0, step, 12.0);
      const auto result =
          IdentifyMetalPerimeter(MakeInput({{trace, 0, 1.0}, {ground, 0, 1.0}}, R));
      int stacks = 0;
      for (const auto &feature : result.features)
      {
        if (feature.type == "ParallelEdgeCluster" ||
            feature.type == "CurvedParallelEdgeCluster")
        {
          INFO("step " << step << " " << feature.type << " " << feature.signature.dump());
          std::vector<double> offsets;
          for (const auto &edge : feature.signature["Edges"])
          {
            offsets.push_back(edge["OffsetOverR"].get<double>());
          }
          std::sort(offsets.begin(), offsets.end());
          REQUIRE(offsets.size() == 3);
          // Either orientation: the consecutive separations are the gap (1R) and the
          // trace width (1.5R), the span 2.5R.
          std::vector<double> gaps = {offsets[1] - offsets[0], offsets[2] - offsets[1]};
          std::sort(gaps.begin(), gaps.end());
          CHECK_THAT(gaps[0], WithinAbs(1.0, 1.0e-9));
          CHECK_THAT(gaps[1], WithinAbs(1.5, 1.0e-9));
          CHECK_THAT(offsets[2], WithinAbs(2.5, 1.0e-9));
          CHECK(feature.exact_parameters);
          if (feature.type == "CurvedParallelEdgeCluster")
          {
            CHECK_THAT(feature.signature["RadiusOverR"].get<double>(),
                       WithinAbs(3.0, 1.0e-9));
          }
          keys.push_back(feature.type + feature.signature_key);
          stacks++;
        }
      }
      CHECK(stacks == 2);
    }
    std::sort(keys.begin(), keys.end());
    REQUIRE(keys.size() == 4);
    CHECK(keys[0] == keys[1]);
    CHECK(keys[2] == keys[3]);
  }
  SECTION("offset polylines agree within the tolerance")
  {
    // The mitre-offset construction of ArcBar (parallel chords, vertices off the design
    // circles): no exact arc reading, the chord reading stays within the signature
    // parameter tolerance of the design offsets at both discretisations (the recorded
    // ambiguity).
    for (const double step : {5.0, 2.5})
    {
      const double centre = 3.0 * R + 1.5;
      const auto trace = ArcBar(3.0, centre, 90.0, step, 12.0, 0.0);
      const auto ground = ArcBar(8.0, centre, 90.0, step, 12.0, -(1.5 + 2.0 + 4.0));
      const auto result =
          IdentifyMetalPerimeter(MakeInput({{trace, 0, 1.0}, {ground, 0, 1.0}}, R));
      for (const auto &feature : result.features)
      {
        if (feature.type == "ParallelEdgeCluster" ||
            feature.type == "CurvedParallelEdgeCluster")
        {
          INFO("step " << step << " " << feature.type << " " << feature.signature.dump());
          std::vector<double> offsets;
          for (const auto &edge : feature.signature["Edges"])
          {
            offsets.push_back(edge["OffsetOverR"].get<double>());
          }
          std::sort(offsets.begin(), offsets.end());
          REQUIRE(offsets.size() == 3);
          std::vector<double> gaps = {offsets[1] - offsets[0], offsets[2] - offsets[1]};
          std::sort(gaps.begin(), gaps.end());
          CHECK_THAT(gaps[0], WithinAbs(1.0, kSignatureParameterToleranceOverRadius));
          CHECK_THAT(gaps[1], WithinAbs(1.5, kSignatureParameterToleranceOverRadius));
          CHECK_THAT(offsets[2], WithinAbs(2.5, kSignatureParameterToleranceOverRadius));
        }
      }
    }
  }
  SECTION("signature tolerance API")
  {
    auto Stack = [&](std::vector<double> offsets, std::vector<int> conductors)
    {
      std::vector<TranslationalEdge> edges;
      const int gaps[4] = {1, -1, 1, -1};
      for (std::size_t i = 0; i < offsets.size(); i++)
      {
        edges.push_back({offsets[i] * R, gaps[i], conductors[i], {"SA"}, "{}"});
      }
      nlohmann::json signature = CanonicalTranslationalSignature(edges, R).signature;
      signature["Type"] = "ParallelEdgeCluster";
      return signature;
    };
    const auto a = Stack({0.0, 1.0004, 2.0006, 3.0008}, {1, 2, 2, 1});
    const auto b = Stack({0.0, 0.9996, 1.9998, 3.0002}, {1, 2, 2, 1});
    const auto far = Stack({0.0, 1.003, 2.004, 3.005}, {1, 2, 2, 1});
    const auto other = Stack({0.0, 1.0, 2.0, 3.0}, {1, 2, 3, 1});
    const auto pa = SplitSignatureParameters(a), pb = SplitSignatureParameters(b);
    CHECK(pa.topology_key == pb.topology_key);
    CHECK(pa.lengths_over_R.size() == 4);
    REQUIRE(SignatureDeviation(a, b).has_value());
    CHECK(*SignatureDeviation(a, b) <= 1.0);
    CHECK(*SignatureDeviation(a, MirrorTranslationalSignature(b)) <= 1.0);
    CHECK(*SignatureDeviation(a, far) > 1.0);
    CHECK(!SignatureDeviation(a, other).has_value());
    const auto representative = RepresentativeSignature({a, b});
    CHECK(*SignatureDeviation(representative, a) <= 1.0);
    CHECK(*SignatureDeviation(representative, b) <= 1.0);
    CHECK(RepresentativeSignature({b, a}) == representative);
    // A corner: angles within 1e-2 deg agree, beyond do not.
    const auto c1 = CanonicalCornerSignature({"SA"}, "{}", 90.0, 0.25);
    const auto c2 = CanonicalCornerSignature({"SA"}, "{}", 90.005, 0.2505);
    const auto c3 = CanonicalCornerSignature({"SA"}, "{}", 90.05, 0.25);
    CHECK(*SignatureDeviation(c1, c2) <= 1.0);
    CHECK(*SignatureDeviation(c1, c3) > 1.0);
  }
  SECTION("cluster extension on the two pads")
  {
    const std::vector<LoopSpec> pads = {{Rectangle(-13.0, -6.0, -1.0, 6.0), 0, 1.0},
                                        {Rectangle(1.0, -6.0, 13.0, 6.0), 0, 1.0}};
    const auto input = MakeInput(pads, R);
    const auto result = IdentifyMetalPerimeter(input);
    CheckPartition(input, result);
    // The gap edges (x = -1 and x = 1, |y| < 6): the end clusters claim R around their
    // cores; the 2 um gap between them is the pair; no isolated portion remains on them
    // (a gap-edge portion facing the cluster across at 1R joins the cluster).
    for (const auto &feature : result.features)
    {
      if (feature.type == "IsolatedEdge")
      {
        for (const auto &portion : feature.portions)
        {
          const auto &s = input.segments[portion.segment];
          INFO("isolated portion (" << s.p0[0] << "," << s.p0[1] << ")-(" << s.p1[0] << ","
                                    << s.p1[1] << ") " << portion.s0 << ".." << portion.s1);
          CHECK(std::abs(std::abs(s.p0[0]) - 1.0) > 1.0e-9);
        }
      }
    }
    CHECK(result.extension.passes >= 1);
    CHECK(result.extension.length >= 0.0);
  }
}

namespace
{

// A hairpin strip (synthetic_layouts.py hairpin): width 2 rho - gap folded through a
// semicircle of centreline radius rho about the origin (the fold below y = 0), the legs up
// to y = half_y and closed across the top; inner fold radius gap / 2, outer 2 rho - gap /
// 2, each fold a polyline of `chords` chords (counter-clockwise loop).
std::vector<Point2> Hairpin(double rho, double gap, double half_y, int chords)
{
  const double r_in = 0.5 * gap, r_out = 2.0 * rho - 0.5 * gap;
  std::vector<Point2> points = {{-r_out, half_y}};
  for (int k = 0; k <= chords; k++)
  {
    const double angle = std::acos(-1.0) * (1.0 + static_cast<double>(k) / chords);
    points.push_back({r_out * std::cos(angle), r_out * std::sin(angle)});
  }
  points.push_back({r_out, half_y});
  points.push_back({r_in, half_y});
  for (int k = chords; k >= 0; k--)
  {
    const double angle = std::acos(-1.0) * (1.0 + static_cast<double>(k) / chords);
    points.push_back({r_in * std::cos(angle), r_in * std::sin(angle)});
  }
  points.push_back({-r_in, half_y});
  return points;
}

}  // namespace

TEST_CASE("SurfaceResponseIdentificationArcClusters",
          "[surfaceresponseidentification][Serial]")
{
  // Option A (decision 91(1)): a cluster holding fitted arcs is described on the design
  // circles, so its signature, extent and edge count do not depend on the chord count.
  const double R = 2.0;
  SECTION("rounded finger ends: two fillets and the end edge")
  {
    // A 3 x 20 um finger (1.5 R wide) with 1 um fillets (0.5 R: rounded corners): the two
    // corner sites of each end are 3 um < 2R apart, an event of their own; each end is one
    // cluster of two arcs, the end edge and R along both sides (5 edges), identical for
    // 2 / 4 / 8 / 16 chords per fillet; the long sides between the clusters are a strip.
    std::optional<std::string> hash;
    std::optional<nlohmann::json> signature;
    for (const int chords : {2, 4, 8, 16})
    {
      const auto input = MakeInput({{RoundedRectangle(10.0, 1.5, 1.0, chords), 0, 0.5}}, R);
      const auto result = IdentifyMetalPerimeter(input);
      CheckPartition(input, result);
      std::map<std::string, int> counts;
      std::set<std::string> hashes;
      for (const auto &feature : result.features)
      {
        counts[feature.type]++;
        if (feature.type == "SpatialEdgeCluster")
        {
          hashes.insert(feature.hash);
          CHECK(feature.signature["EdgeCount"].get<int>() == 5);
          int arcs = 0;
          for (const auto &portion : feature.signature["Portions"])
          {
            arcs += portion.contains("Arc") ? 1 : 0;
            if (portion.contains("Arc"))
            {
              // Metal inside the fillet circle.
              CHECK(portion["GapRadial"].get<int>() == 1);
            }
          }
          CHECK(arcs == 2);
          if (!signature)
          {
            signature = feature.signature;
          }
          else
          {
            INFO("chords " << chords << ": " << feature.signature.dump() << " vs "
                           << signature->dump());
            CHECK(feature.signature == *signature);
          }
        }
      }
      INFO("chords " << chords);
      CHECK(counts["SpatialEdgeCluster"] == 2);
      CHECK(hashes.size() == 1);  // the two ends are congruent
      if (!hash)
      {
        hash = *hashes.begin();
      }
      else
      {
        CHECK(*hashes.begin() == *hash);
      }
      CHECK(counts["ConvexCorner"] == 0);
      CHECK(result.arcs.size() == 4);
      for (const auto &arc : result.arcs)
      {
        CHECK(arc.kind == "RoundedCorner");
        CHECK_THAT(arc.radius, WithinAbs(1.0, 1.0e-9));
        CHECK_THAT(arc.turn_degrees, WithinAbs(90.0, 1.0e-9));
      }
      // Every mesh segment of a fillet (its chords, subdivided at 0.5 um) points at its
      // arc.
      std::size_t on_arcs = 0, arc_segments = 0;
      for (const auto &segment : result.segments)
      {
        on_arcs += segment.arc >= 0 ? 1 : 0;
      }
      for (const auto &arc : result.arcs)
      {
        arc_segments += arc.segments;
      }
      CHECK(on_arcs == arc_segments);
      CHECK(on_arcs >= static_cast<std::size_t>(4 * chords));
    }
  }
  SECTION("hairpin: a rounded corner and a concentric bend")
  {
    // Centreline radius 1.5 R, legs 1.8 R apart: inner fold 0.9 R (a 180 deg rounded
    // corner), outer fold 2.1 R (a bend) 1.2 R away, the strip 1.2 R wide. The outer fold
    // joins the corner across (decision 85(2)): one cluster of two concentric semicircles
    // and the four leg pieces of length R (6 edges), identical for 8 / 16 / 32 chords. The
    // loop closes across the top with sharp corners 1.8 R apart: a second, straight
    // cluster.
    std::optional<nlohmann::json> signature;
    for (const int chords : {8, 16, 32})
    {
      const auto input = MakeInput({{Hairpin(1.5 * R, 1.8 * R, 30.0, chords), 0, 0.5}}, R);
      const auto result = IdentifyMetalPerimeter(input);
      CheckPartition(input, result);
      std::map<std::string, int> counts;
      int arc_clusters = 0;
      for (const auto &feature : result.features)
      {
        counts[feature.type]++;
        if (feature.type != "SpatialEdgeCluster")
        {
          continue;
        }
        int arcs = 0;
        for (const auto &portion : feature.signature["Portions"])
        {
          arcs += portion.contains("Arc") ? 1 : 0;
        }
        if (arcs == 0)
        {
          continue;
        }
        arc_clusters++;
        CHECK(arcs == 2);
        CHECK(feature.signature["EdgeCount"].get<int>() == 6);
        if (!signature)
        {
          signature = feature.signature;
        }
        else
        {
          INFO("chords " << chords << ": " << feature.signature.dump() << " vs "
                         << signature->dump());
          CHECK(feature.signature == *signature);
        }
        CHECK(feature.vertices.empty());  // the rounded corner is a virtual site
        CHECK(feature.signature["Vertices"].size() == 1);
      }
      INFO("chords " << chords);
      CHECK(arc_clusters == 1);
      CHECK(counts["SpatialEdgeCluster"] == 2);
      CHECK(counts["ParallelEdgeCluster"] == 1);  // the 4-edge stack of the legs
      CHECK(counts["ConcaveCorner"] == 0);
    }
  }
  SECTION("closed circles: small round pads inside a cluster")
  {
    // A disc of radius 0.5 R with its centre 1.2 R above a long bar edge lies entirely
    // within 2R of the bar: the cluster claims the WHOLE circle as one closed-circle
    // portion, serialised as centre + radius (its ends at the frame's +x point of the
    // circle, the midpoint at the antipode), identical for 24 / 36 chords and for a mesh
    // starting at another joint (the polygon rotated by a fraction of a chord). Two such
    // discs alone (0.3 R, centres 1.2 R apart) form a cluster of closed circles only, whose
    // frame comes from the centre-to-centre direction.
    auto Disc = [&](double cx, double cy, double r, int chords, double start_degrees)
    {
      std::vector<Point2> disc;
      for (int k = 0; k < chords; k++)
      {
        const double angle = (start_degrees + 360.0 * k / chords) * std::acos(-1.0) / 180.0;
        disc.push_back({cx + r * std::cos(angle), cy + r * std::sin(angle)});
      }
      return disc;
    };
    const std::vector<std::pair<int, double>> meshes = {
        {24, 0.0}, {36, 0.0}, {24, 7.0}, {36, 3.5}};
    auto ClosedCirclePortions = [&](const IdentificationResult &result, int expected_arcs,
                                    std::optional<nlohmann::json> &signature,
                                    std::optional<std::string> &hash)
    {
      int clusters = 0;
      for (const auto &feature : result.features)
      {
        if (feature.type != "SpatialEdgeCluster")
        {
          continue;
        }
        clusters++;
        int closed = 0;
        for (const auto &portion : feature.signature["Portions"])
        {
          if (!portion.contains("Arc"))
          {
            continue;
          }
          const auto P = portion["P"].get<std::array<double, 4>>();
          const auto arc = portion["Arc"].get<std::array<double, 4>>();
          CHECK(P[0] == P[2]);
          CHECK(P[1] == P[3]);
          // Ends at the circle's +x point, midpoint at the antipode.
          const double r = P[0] - arc[0];
          CHECK(r > 0.0);
          CHECK_THAT(P[1], WithinAbs(arc[1], 1.0e-6));
          CHECK_THAT(arc[2], WithinAbs(arc[0] - r, 1.0e-6));
          CHECK_THAT(arc[3], WithinAbs(arc[1], 1.0e-6));
          CHECK(portion["GapRadial"].get<int>() == 1);
          closed++;
        }
        CHECK(closed == expected_arcs);
        CHECK(!feature.hash.empty());
        if (!signature)
        {
          signature = feature.signature;
          hash = feature.hash;
        }
        else
        {
          INFO(feature.signature.dump() << " vs " << signature->dump());
          CHECK(feature.signature == *signature);
          CHECK(feature.hash == *hash);
        }
      }
      return clusters;
    };
    {
      std::optional<nlohmann::json> signature;
      std::optional<std::string> hash;
      for (const auto &[chords, start_degrees] : meshes)
      {
        const auto input =
            MakeInput({{Rectangle(-20.0, -10.0, 20.0, 0.0), 0, 1.0},
                       {Disc(0.0, 1.2 * R, 0.5 * R, chords, start_degrees), 1, 0.25}},
                      R);
        const auto result = IdentifyMetalPerimeter(input);
        CheckPartition(input, result);
        INFO("bar + disc: chords " << chords << " start " << start_degrees);
        CHECK(ClosedCirclePortions(result, 1, signature, hash) == 1);
      }
      REQUIRE(signature);
      CHECK((*signature)["Portions"].size() == 2);  // the bar window and the circle
    }
    {
      std::optional<nlohmann::json> signature;
      std::optional<std::string> hash;
      for (const auto &[chords, start_degrees] : meshes)
      {
        const auto input =
            MakeInput({{Disc(0.0, 0.0, 0.3 * R, chords, start_degrees), 0, 0.2},
                       {Disc(1.2 * R, 0.0, 0.3 * R, chords, start_degrees + 5.0), 1, 0.2}},
                      R);
        const auto result = IdentifyMetalPerimeter(input);
        CheckPartition(input, result);
        INFO("two discs: chords " << chords << " start " << start_degrees);
        CHECK(ClosedCirclePortions(result, 2, signature, hash) == 1);
        CHECK(result.features.size() == 1);
      }
      REQUIRE(signature);
      CHECK((*signature)["Portions"].size() == 2);
    }
  }
}

TEST_CASE("SurfaceResponseIdentificationDecision184Fixes",
          "[surfaceresponseidentification][Serial]")
{
  // USER decision 184 (2026-10-01): the three identification defects of the stage-0 audit
  // (E8-1 / E8-3 / E8-7) and the cluster-composition band of the knife-edge census.
  const double R = 1.9;
  auto Counts = [](const IdentificationResult &result)
  {
    std::map<std::string, int> counts;
    for (const auto &feature : result.features)
    {
      counts[feature.type]++;
    }
    return counts;
  };
  SECTION("near-parallel rigid runs: a mesh-noise tilt of one stack edge (E8-1)")
  {
    // ground | 2 | trace 2 | 2 | ground over 80 um; the lower ground's top edge is tilted
    // by 2e-4 um over its 80 um (2.5e-6 rad: the DS-SCT-002 flux lines' 1.8e-4 um wobble
    // over 186 um), i.e. by more than the 1e-9 DirectionKey grid and less than the parallel
    // cosine tolerance (1.4e-4 rad). Before the one-class rule the translational stage put
    // it in another direction class and the bent-pair stage skipped it as "exactly
    // parallel": a 3-edge stack + an 80 um IsolatedEdge facing it at 1.05 R. Rule: one
    // 4-edge ParallelEdgeCluster, as for the untilted geometry.
    std::map<double, std::map<std::string, double>> lengths;  // per tilt, per type
    for (const double tilt : {0.0, 2.0e-4})
    {
      const auto input = MakeInput(
          {{{{-40.0, -10.0}, {40.0, -10.0}, {40.0, -2.0}, {-40.0, -2.0 - tilt}}, 0, 1.0},
           {Rectangle(-30.0, 0.0, 30.0, 2.0), 0, 1.0},
           {Rectangle(-40.0, 4.0, 40.0, 12.0), 0, 1.0}},
          R);
      const auto result = IdentifyMetalPerimeter(input);
      CheckPartition(input, result);
      INFO("tilt " << tilt);
      int stacks = 0;
      for (const auto &feature : result.features)
      {
        lengths[tilt][feature.type] += feature.length;
        if (feature.type == "ParallelEdgeCluster")
        {
          INFO(feature.signature.dump());
          CHECK(feature.signature["Edges"].size() == 4);
          CHECK(feature.exact_parameters);
          stacks++;
        }
      }
      CHECK(stacks == 1);
    }
    // The tilted geometry reads as the untilted one: the same length per feature type (the
    // outer ground edges are the only isolated edges; the stack has four sides over the
    // whole run between the end clusters).
    REQUIRE(lengths.size() == 2);
    const auto &untilted = lengths.at(0.0), &tilted = lengths.at(2.0e-4);
    CHECK(untilted.size() == tilted.size());
    for (const auto &[type, length] : untilted)
    {
      INFO(type);
      REQUIRE(tilted.count(type) == 1);
      CHECK_THAT(tilted.at(type), WithinAbs(length, 0.01));
    }
    CHECK(untilted.at("IsolatedEdge") >= 2.0 * 80.0);  // the outer ground edges
    CHECK(untilted.at("ParallelEdgeCluster") >
          4.0 * 2.0 * (30.0 - std::sqrt(12.0) - 2.0 * R));
  }
  SECTION("a long straight edge with tiny same-sign end joints is no arc (E8-3)")
  {
    // A 500 um straight island edge whose ends each carry two 1.2 / 2.4 deg joints on 12 um
    // pieces (the first chords of the bends a chip trace enters; noise under the geometric
    // joint rule, so one chain; exactly concyclic by mirror symmetry). Before the
    // every-point clause the four joints alone passed the concyclicity test: a bend of
    // radius ~4 mm bowing 7.8 um (4 R) off the metal toward the ground edge 10 um away ->
    // false 2-edge clusters over most of the run (DS-OSC-003 3.9 mm, DS-SCT-002 2.1 mm).
    // Rule: no arc; the edge and the ground edge are isolated (10 um > 2R apart).
    const double half = 250.0, piece = 12.0, deg = std::acos(-1.0) / 180.0;
    const Point2 m0 = {-half, 0.0}, m1 = {half, 0.0};
    const Point2 q1 = {m1[0] + piece * std::cos(1.2 * deg),
                       m1[1] + piece * std::sin(1.2 * deg)};
    const Point2 q2 = {q1[0] + piece * std::cos(3.6 * deg),
                       q1[1] + piece * std::sin(3.6 * deg)};
    const Point2 p1 = {-q1[0], q1[1]}, p2 = {-q2[0], q2[1]};
    const double top = 20.0;
    const std::vector<Point2> island = {p2, p1, m0, m1, q1, q2, {q2[0], top}, {p2[0], top}};
    for (const bool with_ground : {false, true})
    {
      std::vector<LoopSpec> loops = {{island, 0, 1.0}};
      if (with_ground)
      {
        loops.push_back({Rectangle(-half - 60.0, -30.0, half + 60.0, -10.0), 1, 1.0});
      }
      const auto input = MakeInput(loops, R);
      const auto result = IdentifyMetalPerimeter(input);
      CheckPartition(input, result);
      INFO("with ground " << with_ground);
      for (const auto &arc : result.arcs)
      {
        INFO("arc " << arc.kind << " radius " << arc.radius << " joints " << arc.joints
                    << " sagitta " << arc.max_sagitta_over_R);
        CHECK(arc.kind != "Bend");
      }
      CHECK(result.arcs.empty());
      const auto counts = Counts(result);
      CHECK(counts.count("SpatialEdgeCluster") == 0);
      CHECK(counts.count("CurvedEdge") == 0);
      CHECK(counts.at("IsolatedEdge") >= 4);
    }
  }
  SECTION("a taper strip reads its own width, not its narrow lead's (E8-7)")
  {
    // A 2 um lead (x < 0) widening along a smooth (quadratic) taper to 3.7 um at x = 120
    // (just under 2R = 3.8), chorded every 5 um: every joint is noise, each side is one
    // chain of many runs, the chord readings step by < 5 % between consecutive runs (a slow
    // taper: constant within the pair tolerance). Before the amendment the 5 % steps linked
    // the lead's exact 2.0 um pieces to every taper piece and the link read the exact mean
    // 2.0 um (ExactParameters true) over the whole strip up to 3.7 um (the DS-CTX-003
    // launcher strips keyed 2.0 um at 3.7 um). Rule: a strip feature's separation is the
    // length-weighted mean of the local width over its portions, and it is exact only where
    // the width is constant within the parameter tolerance.
    const double lead = 60.0, taper = 120.0, w0 = 2.0, w1 = 3.7, step = 5.0;
    auto HalfWidth = [&](double x)
    {
      if (x <= 0.0)
      {
        return 0.5 * w0;
      }
      const double u = std::min(x, taper) / taper;
      return 0.5 * (w0 + (w1 - w0) * u * u);
    };
    std::vector<Point2> bar = {{-lead, -HalfWidth(-lead)}};
    const int chords = static_cast<int>(std::lround(taper / step));
    for (int k = 0; k <= chords; k++)
    {
      const double x = taper * k / chords;
      bar.push_back({x, -HalfWidth(x)});
    }
    for (int k = chords; k >= 0; k--)
    {
      const double x = taper * k / chords;
      bar.push_back({x, HalfWidth(x)});
    }
    bar.push_back({-lead, HalfWidth(-lead)});
    const auto input = MakeInput({{bar, 0, 1.0}}, R);
    const auto result = IdentifyMetalPerimeter(input);
    CheckPartition(input, result);
    double strip_length = 0.0;
    for (const auto &feature : result.features)
    {
      if (feature.type != "SameConductorStrip")
      {
        continue;
      }
      const double separation = feature.signature["SeparationOverR"].get<double>() * R;
      double length = 0.0, width_sum = 0.0, width_min = 1.0e300, width_max = 0.0;
      for (const auto &portion : feature.portions)
      {
        const auto &segment = result.segments[portion.segment];
        const double s = 0.5 * (portion.s0 + portion.s1) / segment.length;
        const double x = segment.key[0][0] + s * (segment.key[1][0] - segment.key[0][0]);
        const double width = 2.0 * HalfWidth(x);
        const double l = portion.s1 - portion.s0;
        length += l;
        width_sum += l * width;
        width_min = std::min(width_min, width);
        width_max = std::max(width_max, width);
      }
      strip_length += length;
      INFO("strip " << feature.signature.dump() << " exact " << feature.exact_parameters
                    << " length " << length << " width " << width_min << " .. "
                    << width_max);
      // The feature's separation is its own length-weighted mean width (the lateral
      // chords of a slow taper read the local width within the 5 % pair tolerance).
      CHECK_THAT(separation, WithinRel(width_sum / length, 0.05));
      if (feature.exact_parameters)
      {
        // Exact only on the constant lead (chord pieces within the pair tolerance of the
        // exact 2.0 um may join its link: a width within 5 % of the exact value).
        CHECK(width_max <= 1.05 * separation);
      }
    }
    // The lead and the taper up to the 2R crossing are strip material (the end clusters
    // and the taper end excepted).
    CHECK(strip_length > 2.0 * (lead + taper) * 0.6);
  }
  SECTION("knife-edge census: cluster-composition band")
  {
    // Pads A | 3.0 | B (3 um wide) | 3.82 | C: at R the corners of A's right edge and of B
    // are one cluster (every gap below 2R = 3.8); at R (1 + 1 %) the gap of 3.82 um is
    // below 2R = 3.838 and C's left corners join it (edge count / member vertices change);
    // at R (1 - 1 %) nothing changes. The census reports the cluster's claimed length on
    // the Above side only.
    const auto input = MakeInput({{Rectangle(0.0, 0.0, 8.0, 8.0), 0, 1.0},
                                  {Rectangle(11.0, 0.0, 14.0, 8.0), 1, 1.0},
                                  {Rectangle(17.82, 0.0, 25.82, 8.0), 2, 1.0}},
                                 R);
    const auto result = IdentifyMetalPerimeter(input);
    CheckPartition(input, result);
    const auto counts = Counts(result);
    REQUIRE(counts.at("SpatialEdgeCluster") >= 1);
    const auto census = nlohmann::json::parse(result.knife_edge_census);
    REQUIRE(census.contains("ClusterComposition"));
    const auto &band = census["ClusterComposition"];
    INFO(band.dump());
    CHECK(band["Clusters"].get<int>() == counts.at("SpatialEdgeCluster"));
    CHECK(band["ClustersBelow"].get<int>() == 0);
    CHECK_THAT(band["Below"].get<double>(), WithinAbs(0.0, 1.0e-12));
    CHECK(band["ClustersAbove"].get<int>() >= 1);
    CHECK(band["Above"].get<double>() > 0.0);
    CHECK_THAT(
        band["Total"].get<double>(),
        WithinAbs(band["Below"].get<double>() + band["Above"].get<double>(), 1.0e-9));
    // Stable without a knife edge: the pads 6 um apart report nothing.
    const auto quiet =
        IdentifyMetalPerimeter(MakeInput({{Rectangle(0.0, 0.0, 8.0, 8.0), 0, 1.0},
                                          {Rectangle(11.0, 0.0, 14.0, 8.0), 1, 1.0},
                                          {Rectangle(20.0, 0.0, 28.0, 8.0), 2, 1.0}},
                                         R));
    const auto quiet_band =
        nlohmann::json::parse(quiet.knife_edge_census)["ClusterComposition"];
    INFO(quiet_band.dump());
    CHECK(quiet_band["ClustersBelow"].get<int>() == 0);
    CHECK(quiet_band["ClustersAbove"].get<int>() == 0);
  }
  SECTION("a straight lead next to a slow taper keeps its exact reading (184 (3) locality)")
  {
    // ground | 1 | trace 2 | 1 | ground, straight over 80 um (x < 0) and then one LINEAR
    // taper run of 300 um to 2.5 / 3.5 / 2.5 um (the DS-CTX-003 launcher routes: 0.9-1.3 mm
    // tapers from 1 / 2 / 1 um to ~4 um); every taper joint is noise (0.14-0.43 deg), so
    // every edge is one chain of two runs. Before the locality correction the curvature
    // rule's half-run spreading of the taper joint's turn made kappa > 0 over 40 um of the
    // lead, the lead's samples there had no exact reading and read a window of the whole
    // run length (80 um) deep into the taper: the straight 1 / 2 / 1 um stack read
    // 1 / 2.3 / 1 um and ExactParameters false (the chip's 532 um 6-edge stacks 2.2 / 2.7 /
    // 2.0 um). Rule: the straight part is one exact 4-edge stack at 0 / 1 / 3 / 4 um; the
    // taper reads its own mean separations, not exact.
    const double L = 80.0, T = 300.0, w0 = 1.0, w1 = 1.75, g0 = 1.0, g1 = 2.5, H = 12.0;
    const auto input = MakeInput(
        {{{{-L, -w0}, {0.0, -w0}, {T, -w1}, {T, w1}, {0.0, w0}, {-L, w0}}, 0, 1.0},
         {{{-L, -w0 - g0 - H},
           {T, -w1 - g1 - H},
           {T, -w1 - g1},
           {0.0, -w0 - g0},
           {-L, -w0 - g0}},
          0,
          1.0},
         {{{-L, w0 + g0},
           {0.0, w0 + g0},
           {T, w1 + g1},
           {T, w1 + g1 + H},
           {-L, w0 + g0 + H}},
          0,
          1.0}},
        R);
    const auto result = IdentifyMetalPerimeter(input);
    CheckPartition(input, result);
    const std::vector<double> design = {0.0, g0 / R, (g0 + 2.0 * w0) / R,
                                        (2.0 * g0 + 2.0 * w0) / R};
    double exact_length = 0.0, taper_length = 0.0;
    for (const auto &feature : result.features)
    {
      if (feature.type != "ParallelEdgeCluster")
      {
        continue;
      }
      const auto &edges = feature.signature["Edges"];
      INFO(feature.signature.dump()
           << " exact " << feature.exact_parameters << " length " << feature.length);
      REQUIRE(edges.size() == 4);
      if (feature.exact_parameters)
      {
        // An exact stack here is the straight part: the design offsets.
        for (std::size_t k = 0; k < 4; k++)
        {
          CHECK_THAT(edges[k]["OffsetOverR"].get<double>(), WithinAbs(design[k], 1.0e-5));
        }
        exact_length += feature.length;
      }
      else
      {
        // The taper: wider than the design everywhere beyond its first 5 %.
        CHECK(edges[3]["OffsetOverR"].get<double>() > 1.05 * design[3]);
        taper_length += feature.length;
      }
    }
    // The straight part (less the end cluster at x = -L) is exact; the taper is read.
    CHECK(exact_length >= 4.0 * (L - 4.0 * R));
    CHECK(taper_length >= 4.0 * 0.5 * T);
  }
  SECTION("a straight lead next to a bend keeps its exact reading (184 (3) locality)")
  {
    // Two conductors 2.85 um = 1.5 R apart: straight leads of 60 um, concentric 90 deg
    // bends of radius 19 / 16.15 um in 4 chords (7.4 / 6.3 um, beyond 2R: coarse) or 18
    // chords (1.7 / 1.4 um: fine) whose vertices are NOT aligned (the inner polyline's
    // vertices sit mid-chord of the outer's, as a mesh generator places them), then
    // straight leads again. An outer vertex (on its circle) is farther from the inner chord
    // across it than the design separation (3.16 um coarse, 2.865 um fine), so the bend's
    // sampled distances exceed 1.5 R. Before the correction the lead's samples within half
    // the lead of the bend joint read a window of the whole run (60 um), which held the
    // bend's locally constant samples (fine chords; the coarse chords' readings vary by
    // more than the pair tolerance and were masked): the fine case read the lead's pair
    // non-exact above 1.5 R. Rule: the leads are the exact 1.5 R pair in both cases (their
    // samples' own distances are the line distance).
    const double rho = 19.0, sep = 1.5 * R, lead = 60.0, quarter = 0.5 * std::acos(-1.0);
    for (const int chords : {4, 18})
    {
      std::vector<Point2> outer = {
          {-lead, -10.0}, {30.0, -10.0}, {30.0, 40.0}, {rho, 40.0}};
      for (int k = chords; k >= 0; k--)
      {
        const double theta = quarter * k / chords;
        outer.push_back({rho * std::sin(theta), rho - rho * std::cos(theta)});
      }
      outer.push_back({-lead, 0.0});
      std::vector<Point2> inner = {{-lead, sep}, {0.0, sep}};
      for (int k = 0; k < chords; k++)
      {
        const double theta = quarter * (k + 0.5) / chords;
        inner.push_back(
            {(rho - sep) * std::sin(theta), rho - (rho - sep) * std::cos(theta)});
      }
      inner.push_back({rho - sep, rho});
      inner.push_back({rho - sep, 40.0});
      inner.push_back({-lead, 40.0});
      const auto input = MakeInput({{outer, 0, 1.0}, {inner, 1, 1.0}}, R);
      const auto result = IdentifyMetalPerimeter(input);
      CheckPartition(input, result);
      double exact_length = 0.0;
      std::ostringstream summary;
      for (const auto &feature : result.features)
      {
        summary << feature.type << " " << feature.length << " exact "
                << feature.exact_parameters << " bend "
                << feature.bend_radius_over_R.value_or(0.0) << "; ";
        if (feature.type != "DifferentConductorGap")
        {
          continue;
        }
        const double separation = feature.signature["SeparationOverR"].get<double>();
        INFO("chords " << chords << " " << feature.signature.dump() << " exact "
                       << feature.exact_parameters << " length " << feature.length
                       << " bend " << feature.bend_radius_over_R.value_or(0.0));
        if (feature.exact_parameters)
        {
          // Exact only at the design separation (the straight-class pair of the leads; the
          // bend is the curved type, exact by its concentric fitted arcs).
          CHECK_THAT(separation, WithinAbs(1.5, 1.0e-5));
          exact_length += feature.length;
        }
      }
      // The horizontal leads (less the end cluster at x = -lead) and the vertical ones
      // read exact: more than the two horizontal leads alone.
      INFO("chords " << chords << ": " << summary.str());
      CHECK(exact_length >= 2.0 * (lead - 4.0 * R));
    }
  }
}

TEST_CASE("SurfaceResponseIdentificationDecision203ExactStretches",
          "[surfaceresponseidentification][Serial]")
{
  // USER decision 203 (2026-10-02): the sub-piece exactness is local along the run (exact
  // stretches), completing decision 184 (3). Before, exactness was all-or-nothing per
  // sub-piece: a straight run facing a partner that is parallel over part of its length and
  // then changes separation one-sidedly lost its exact part entirely (one judged sample
  // made the whole run non-exact at its chord mean; the lead's exact group then had no
  // partner piece and the strip read as two IsolatedEdges), or — when every non-exact
  // sample sat within R of a partner joint (finely chorded taper) and was not judged — read
  // exact at the lead's value over the whole taper (the twin mode, E8-7 again).
  const double R = 1.9;
  // A strip island whose LEFT side is one straight run (x = 0) over the whole length and
  // whose RIGHT side is straight (width w0) for y < 0 and then a quadratic one-sided taper
  // to w1 at y = taper, chorded every `step` um (noise joints of < 0.01 deg), optionally
  // followed by a straight part of width w1 over `tail` um. The reviewer's reproducer
  // (straight 150 um, taper 300 um, w0 2.0 um = 0.53 R, w1 3.7 um < 2R, chords 3 um < 2R
  // and 10 um > 2R).
  const double straight = 150.0, taper = 300.0, w0 = 2.0, w1 = 3.7;
  auto Bar = [&](double step, double tail)
  {
    const int chords = static_cast<int>(std::lround(taper / step));
    std::vector<Point2> right = {{w0, -straight}};
    for (int k = 0; k <= chords; k++)
    {
      const double y = taper * k / chords;
      right.push_back({w0 + (w1 - w0) * (y / taper) * (y / taper), y});
    }
    if (tail > 0.0)
    {
      right.push_back({w1, taper + tail});
    }
    right.push_back({0.0, taper + tail});
    right.push_back({0.0, -straight});
    return right;
  };
  // Per feature: the length of its portions on the left edge (x = 0) and on the right side.
  auto SideLengths =
      [&](const IdentificationResult &result, const IdentifiedFeature &feature)
  {
    std::array<double, 2> lengths = {0.0, 0.0};
    for (const auto &portion : feature.portions)
    {
      const auto &segment = result.segments[portion.segment];
      const bool left =
          std::abs(segment.key[0][0]) < 1.0e-9 && std::abs(segment.key[1][0]) < 1.0e-9;
      lengths[left ? 0 : 1] += portion.s1 - portion.s0;
    }
    return lengths;
  };
  // The left edge's y range of a feature's portions.
  auto LeftRange = [&](const IdentificationResult &result, const IdentifiedFeature &feature)
  {
    double lo = 1.0e300, hi = -1.0e300;
    for (const auto &portion : feature.portions)
    {
      const auto &segment = result.segments[portion.segment];
      if (std::abs(segment.key[0][0]) > 1.0e-9 || std::abs(segment.key[1][0]) > 1.0e-9)
      {
        continue;
      }
      const double y0 = segment.key[0][1], y1 = segment.key[1][1];
      for (const double s : {portion.s0, portion.s1})
      {
        const double y = y0 + (y1 - y0) * s / segment.length;
        lo = std::min(lo, y);
        hi = std::max(hi, y);
      }
    }
    return std::make_pair(lo, hi);
  };
  auto CheckOneSided = [&](double step, double tail)
  {
    INFO("step " << step << " tail " << tail);
    const auto input = MakeInput({{Bar(step, tail), 0, 1.0}}, R);
    const auto result = IdentifyMetalPerimeter(input);
    CheckPartition(input, result);
    double exact_left = 0.0, exact_right = 0.0, taper_left = 0.0, taper_right = 0.0,
           tail_left = 0.0, isolated = 0.0;
    std::pair<double, double> exact_range = {0.0, 0.0};
    for (const auto &feature : result.features)
    {
      if (feature.type == "IsolatedEdge")
      {
        isolated += feature.length;
        continue;
      }
      if (feature.type != "SameConductorStrip")
      {
        continue;
      }
      const double separation = feature.signature["SeparationOverR"].get<double>() * R;
      const auto sides = SideLengths(result, feature);
      INFO("step " << step << " tail " << tail << " strip " << feature.signature.dump()
                   << " exact " << feature.exact_parameters << " length " << feature.length
                   << " left " << sides[0] << " right " << sides[1]);
      if (feature.exact_parameters && std::abs(separation - w0) < 0.05 * w0)
      {
        // The exact strip of the straight part: exactly w0.
        CHECK_THAT(separation, WithinAbs(w0, 1.0e-5 * R));
        exact_left += sides[0];
        exact_right += sides[1];
        exact_range = LeftRange(result, feature);
      }
      else if (feature.exact_parameters)
      {
        // Only the tail's straight part can be exact otherwise: exactly w1.
        CHECK_THAT(separation, WithinAbs(w1, 1.0e-5 * R));
        tail_left += sides[0];
      }
      else
      {
        // The taper: its own mean width, wider than w0.
        CHECK(separation > 1.05 * w0);
        CHECK(separation < w1);
        taper_left += sides[0];
        taper_right += sides[1];
      }
    }
    INFO("step " << step << " tail " << tail << ": exact left " << exact_left << " right "
                 << exact_right << " taper left " << taper_left << " right " << taper_right
                 << " tail left " << tail_left << " isolated " << isolated);
    // The exact strip's left portion is the straight part (less the end cluster at
    // y = -straight) plus the taper start whose chord pieces lie within the 5 % pair
    // tolerance of w0 (the established link rule of decision 184 (3): a chord piece joins
    // the exact group within the pair tolerance of its own separation; here the first
    // ~73 um of the taper, 2.0 -> 2.1 um, keyed with the exact 2.0 um strip on BOTH
    // sides — the non-exact stretch is cut at the partner's runs so that both sides are
    // discretised alike) and at most one partner chord + R beyond.
    const double five_percent = taper * std::sqrt(0.05 * w0 / (w1 - w0));
    CHECK(exact_left >= straight - 4.0 * R);
    CHECK(exact_left <= straight + five_percent + step + R);
    CHECK(exact_range.second <= five_percent + step + R);
    CHECK(exact_right >= straight - 4.0 * R);
    CHECK(exact_right <= straight + five_percent + step + R);
    CHECK_THAT(exact_left, WithinAbs(exact_right, 2.0 * R));
    // The taper beyond reads non-exact over the rest of its length on both sides (with a
    // tail, its end within 5 % of w1 keys with the tail's exact w1 strip likewise).
    const double five_percent_end =
        tail > 0.0 ? taper * (1.0 - std::sqrt((0.95 * w1 - w0) / (w1 - w0))) + step : 0.0;
    CHECK(taper_left >= taper - five_percent - step - five_percent_end - 4.0 * R);
    CHECK(taper_right >= taper - five_percent - step - five_percent_end - 4.0 * R);
    CHECK_THAT(taper_left, WithinAbs(taper_right, 2.0 * R));
    if (tail > 0.0)
    {
      CHECK(tail_left >= tail - 4.0 * R - R);
    }
    // No part of the strip (0.53 .. 1.95 R wide) is an isolated edge.
    CHECK(isolated == 0.0);
  };
  SECTION("a straight run facing a partner that is parallel then tapers one-sidedly")
  {
    // Chords 3 um (every taper point within R of a noise joint: the left run's samples
    // facing the taper are near-bend, not judged, except one on the last chord) and 10 um
    // (judged samples mid-chord). Before: the strip's straight part read as two
    // IsolatedEdges (220 + 217 um) and the taper as a 455 um non-exact strip.
    for (const double step : {3.0, 10.0})
    {
      CheckOneSided(step, 0.0);
    }
  }
  SECTION("the twin mode: a finely chorded one-sided taper ending at a joint")
  {
    // Chords 2.5 um: every point of the taper lies within R of a noise joint, the taper
    // ends at the island's top corner, so the left run's samples facing the taper are
    // near-bend, not judged, except the last one at the island's corner (stretch log
    // "exact 0, judged 1, near-bend 2" for the final cells; it lies in the non-exact
    // stretch either way). Before: the whole 450 um left run read exact at 2.0 um (the
    // taper's 300 um keyed at the lead's width, E8-7). Rule: the near-bend samples farther
    // than R from the exact stretch form a non-exact stretch at its chord mean.
    CheckOneSided(2.5, 0.0);
  }
  SECTION("a finely chorded one-sided taper between two straight parts")
  {
    // The taper ends in a straight part of width w1 (1.95 R) over 40 um: two exact
    // stretches (w0 and w1) with the taper's non-exact stretch between them.
    CheckOneSided(2.5, 40.0);
  }
  SECTION("parallel-class stack: the tilt x half-length exactness guard (MINOR-1)")
  {
    // ground | 2 | trace 2 | 2 | ground over 100 um as E8-1, the lower ground's top edge
    // tilted by 1.2e-4 rad (within the parallel cosine tolerance 1.4e-4 rad: one class,
    // one 4-edge stack): its lateral offset is read at the run's midpoint, but the lines'
    // separation at the span's ends deviates by tilt x half-length = 6e-3 um = 3.2e-3 R,
    // beyond the 1e-3 R signature tolerance -> the stack reads the midpoint offsets with
    // ExactParameters false. The E8-1 tilt (2.5e-6 rad, 5e-5 R over the half-length) stays
    // exact.
    for (const double tilt_rad : {2.5e-6, 1.2e-4})
    {
      const double tilt = tilt_rad * 100.0;
      const auto input = MakeInput(
          {{{{-50.0, -10.0}, {50.0, -10.0}, {50.0, -2.0}, {-50.0, -2.0 - tilt}}, 0, 1.0},
           {Rectangle(-40.0, 0.0, 40.0, 2.0), 0, 1.0},
           {Rectangle(-50.0, 4.0, 50.0, 12.0), 0, 1.0}},
          R);
      const auto result = IdentifyMetalPerimeter(input);
      CheckPartition(input, result);
      int stacks = 0;
      for (const auto &feature : result.features)
      {
        if (feature.type != "ParallelEdgeCluster" || feature.length < R)
        {
          // (The tilted class axis projects the members' ends differently: a 0.02 um
          // 3-edge sliver at one stack end, below the signature grid's interest.)
          continue;
        }
        INFO("tilt " << tilt_rad << " " << feature.signature.dump() << " exact "
                     << feature.exact_parameters << " length " << feature.length);
        CHECK(feature.signature["Edges"].size() == 4);
        CHECK(feature.exact_parameters == (tilt_rad < 1.0e-5));
        stacks++;
      }
      CHECK(stacks == 1);
    }
  }
  SECTION("no regression: an oblique straight strip with grid-rounded vertices")
  {
    // A 2.6 um = 1.37 R strip of 120 um at 23 deg whose vertices (every 10 um on both
    // edges) are rounded to a 1 nm grid, so the chords carry noise joints of ~1e-4 rad and
    // facing chords are parallel within the cosine tolerance or not by chance: a mix of
    // exact stretches and judged (non-exact) samples along one straight run. Reads as at
    // e2b148f9bd: ONE exact strip over the whole length (less the end clusters) at the
    // mean of the exact readings, within the parameter tolerance of the design width.
    const double length = 120.0, width = 2.6, angle = 23.0 * std::acos(-1.0) / 180.0;
    const Point2 t = {std::cos(angle), std::sin(angle)}, m = {-t[1], t[0]};
    auto Grid = [](const Point2 &p)
    {
      return Point2{std::round(p[0] * 1000.0) / 1000.0, std::round(p[1] * 1000.0) / 1000.0};
    };
    std::vector<Point2> strip;
    const int chords = 12;
    for (int k = 0; k <= chords; k++)
    {
      const double s = length * k / chords;
      strip.push_back(Grid({t[0] * s, t[1] * s}));
    }
    for (int k = chords; k >= 0; k--)
    {
      const double s = length * k / chords;
      strip.push_back(Grid({t[0] * s + m[0] * width, t[1] * s + m[1] * width}));
    }
    const auto input = MakeInput({{strip, 0, 1.0}}, R);
    const auto result = IdentifyMetalPerimeter(input);
    CheckPartition(input, result);
    int strips = 0;
    for (const auto &feature : result.features)
    {
      INFO(feature.type << " " << feature.signature.dump() << " exact "
                        << feature.exact_parameters << " length " << feature.length);
      CHECK(feature.type != "IsolatedEdge");
      if (feature.type != "SameConductorStrip")
      {
        continue;
      }
      strips++;
      const double separation = feature.signature["SeparationOverR"].get<double>() * R;
      CHECK(feature.exact_parameters);
      CHECK_THAT(separation, WithinAbs(width, 1.0e-3 * R));
      CHECK(feature.length >= 2.0 * (length - 4.0 * R));
    }
    CHECK(strips == 1);
  }
}

TEST_CASE("SurfaceResponseIdentificationCollinearSubdivision",
          "[surfaceresponseidentification][Serial]")
{
  // Collinear-subdivision invariance (decision 212; VALIDATION-PLAN (h)-8): inserting
  // collinear vertices on a chord is a geometric no-op, so the identification must read a
  // joint-only polyline and the same polyline with every chord subdivided (a thin window
  // mesh at LC 4 um, a second-order mid-edge node) identically, on the CURRENT arc / corner
  // rule. The defect: a coarse round pad with a straight lead (DS-CTX-003 C4: r 23 um, 8
  // chords of 43 deg, sagitta 0.86 R; the ground's r 36 hole, 9 chords of 33 deg) read a
  // Bend over its interior joints on the chip mesh (one mesh edge per chord: the
  // decision-122 coarse-chord exception) and sharp corners on the thin window mesh. The
  // lead attaches through two joints above the 50 deg cap, so the bend is the least-squares
  // fit over the interior joints and its end arms are the first and last chords of the
  // circle; their end-joint test read the arm's far vertex off the MESH SEGMENT (on the
  // circle for the design chord, off it for a 4 um sub-chord) instead of the far end of the
  // arm's straight piece (the rigid-run joint).
  const double R = 1.9;
  // A round pad of the given radius with a vertical lead of width 2 w and length L
  // attached at the top (counter-clockwise, metal inside): the arc between the lead's
  // attach points in `chords` equal chords (the chip's design polygon).
  auto RoundPadWithLead = [](double rho, int chords, double w, double L)
  {
    const double theta0 = std::atan2(std::sqrt(rho * rho - w * w), w);
    const double y0 = std::sqrt(rho * rho - w * w);
    const double a0 = std::acos(-1.0) - theta0, a1 = 2.0 * std::acos(-1.0) + theta0;
    std::vector<Point2> points = {{w, y0 + L}, {-w, y0 + L}};
    for (int k = 0; k <= chords; k++)
    {
      const double angle = a0 + (a1 - a0) * k / chords;
      points.push_back({rho * std::cos(angle), rho * std::sin(angle)});
    }
    return points;
  };
  // Collinear vertices inserted on every chord: uniform pieces of at most `spacing`, or an
  // irregular pattern (fractions of the chord, rotated from chord to chord).
  auto Subdivide = [](const std::vector<Point2> &points, double spacing)
  {
    std::vector<Point2> result;
    const std::size_t n = points.size();
    for (std::size_t i = 0; i < n; i++)
    {
      const Point2 a = points[i], b = points[(i + 1) % n];
      const double length = std::hypot(b[0] - a[0], b[1] - a[1]);
      const int pieces =
          std::max(1, static_cast<int>(std::ceil(length / spacing - 1.0e-9)));
      for (int k = 0; k < pieces; k++)
      {
        const double s = static_cast<double>(k) / pieces;
        result.push_back({a[0] + (b[0] - a[0]) * s, a[1] + (b[1] - a[1]) * s});
      }
    }
    return result;
  };
  auto SubdivideIrregularly = [](const std::vector<Point2> &points)
  {
    const std::array<double, 4> pattern = {0.11, 0.37, 0.58, 0.83};
    std::vector<Point2> result;
    const std::size_t n = points.size();
    for (std::size_t i = 0; i < n; i++)
    {
      const Point2 a = points[i], b = points[(i + 1) % n];
      result.push_back(a);
      for (std::size_t k = 0; k < pattern.size(); k++)
      {
        const double s = pattern[(k + i) % pattern.size()];
        result.push_back({a[0] + (b[0] - a[0]) * s, a[1] + (b[1] - a[1]) * s});
      }
      // The pattern is rotated per chord, so sort the inserted vertices along the chord.
      std::sort(result.end() - static_cast<std::ptrdiff_t>(pattern.size()), result.end(),
                [&](const Point2 &p, const Point2 &q)
                {
                  return (p[0] - a[0]) * (b[0] - a[0]) + (p[1] - a[1]) * (b[1] - a[1]) <
                         (q[0] - a[0]) * (b[0] - a[0]) + (q[1] - a[1]) * (b[1] - a[1]);
                });
    }
    return result;
  };
  struct Reading
  {
    std::string digest;
    std::map<std::pair<std::string, std::string>, double> length;    // per (type, hash)
    std::vector<std::tuple<std::string, double, std::size_t>> arcs;  // kind, radius, joints
    std::map<std::string, int> vertex_types;
    int corners = 0;  // ConvexCorner / ConcaveCorner features
  };
  auto Read = [&](const std::vector<LoopSpec> &loops)
  {
    const auto input = MakeInput(loops, R);
    const auto result = IdentifyMetalPerimeter(input);
    CheckPartition(input, result);
    Reading reading;
    reading.digest = result.geometry_digest;
    for (const auto &feature : result.features)
    {
      reading.length[{feature.type, feature.hash}] += feature.length;
      reading.corners += feature.type == "ConvexCorner" || feature.type == "ConcaveCorner";
    }
    for (const auto &arc : result.arcs)
    {
      reading.arcs.emplace_back(arc.kind, std::round(arc.radius / R * 1.0e6), arc.joints);
    }
    std::sort(reading.arcs.begin(), reading.arcs.end());
    for (const auto &vertex : result.vertices)
    {
      reading.vertex_types[vertex.type]++;
    }
    return reading;
  };
  auto CheckSame = [](const Reading &a, const Reading &b)
  {
    CHECK(a.digest == b.digest);
    CHECK(a.arcs == b.arcs);
    CHECK(a.vertex_types == b.vertex_types);
    CHECK(a.length.size() == b.length.size());
    for (const auto &[key, length] : a.length)
    {
      INFO(key.first << " " << key.second);
      CHECK(b.length.count(key) == 1);
      if (b.length.count(key) == 1)
      {
        CHECK_THAT(b.length.at(key), WithinAbs(length, 1.0e-6));
      }
    }
  };
  const double joint_only = 1.0e9;  // no subdivision by MakeInput
  SECTION("coarse polygon arcs read alike joint-only and with subdivided chords")
  {
    // Pads of 10-60 um chords whose interior joints turn 6-44 deg (174-136 deg corners),
    // every chord sagitta above the 0.05 R resolution (0.07-0.86 R: the decision-122
    // coarse-chord regime); the lead's attach joints turn above the 50 deg cap.
    struct Pad
    {
      double radius;
      int chords;
    };
    for (const Pad &pad : {Pad{23.0, 8}, Pad{36.0, 9}, Pad{80.0, 8}, Pad{95.5, 60}})
    {
      INFO("pad radius " << pad.radius << " chords " << pad.chords);
      const auto design = RoundPadWithLead(pad.radius, pad.chords, 2.5, 60.0);
      const Reading plain = Read({{design, 0, joint_only}});
      // Joint-only: the interior joints form least-squares bends of the pad's exact radius
      // (each at most 180 deg of turn under the current rule, the leftover joints corners:
      // the reading this block preserves); the attach joints and the lead's end corners
      // are corners.
      REQUIRE(!plain.arcs.empty());
      for (const auto &arc : plain.arcs)
      {
        CHECK(std::get<0>(arc) == "Bend");
        CHECK(std::get<1>(arc) == std::round(pad.radius / R * 1.0e6));
        CHECK(std::get<2>(arc) >= 4);
      }
      CHECK(plain.vertex_types.at("BendVertex") >= 4);
      CHECK(plain.corners >= 4);  // the lead's two end corners and two attach joints
      const Reading uniform = Read({{Subdivide(design, 4.0), 0, joint_only}});
      const Reading irregular = Read({{SubdivideIrregularly(design), 0, joint_only}});
      const Reading meshed = Read({{design, 0, 4.0}});
      CheckSame(plain, uniform);
      CheckSame(plain, irregular);
      CheckSame(plain, meshed);
    }
  }
  SECTION("a straight edge and a genuine sharp corner with subdivided arms")
  {
    // A rectangle (four 90 deg corners, straight edges), a rhombus (two 60 deg turns, two
    // 120 deg turns) and a polygon whose side carries four sub-cap joints (20 / 25 / 30 /
    // 20 deg on 10-30 um pieces) that are NOT concyclic (least-squares residuals 0.02-0.06
    // um against the 1.9e-3 um fit tolerance): no arc at any discretisation; the corners
    // stay corners and the digest, features and lengths are the same.
    const std::vector<Point2> rhombus = {
        {0.0, 0.0}, {30.0, -17.32}, {60.0, 0.0}, {30.0, 17.32}};
    const std::vector<Point2> sub_cap_corners = {
        {30.0, 0.0},        {39.3969, 3.4202}, {57.0746, 21.0979}, {60.9569, 35.5868},
        {58.3422, 65.4726}, {-20.0, 65.4726},  {-20.0, 0.0}};
    for (const auto &design :
         {Rectangle(-30.0, -10.0, 30.0, 10.0), rhombus, sub_cap_corners})
    {
      INFO("polygon of " << design.size() << " vertices");
      const Reading plain = Read({{design, 0, joint_only}});
      CHECK(plain.arcs.empty());
      CHECK(plain.vertex_types.count("BendVertex") == 0);
      CHECK(plain.vertex_types.count("RoundedCornerVertex") == 0);
      CHECK(plain.corners == static_cast<int>(design.size()));
      const Reading uniform = Read({{Subdivide(design, 4.0), 0, joint_only}});
      const Reading irregular = Read({{SubdivideIrregularly(design), 0, joint_only}});
      CheckSame(plain, uniform);
      CheckSame(plain, irregular);
    }
  }
  SECTION("a closed loop whose longest piece is a chord inside the arc (start rule)")
  {
    // The closed-loop scan starts after the loop's LONGEST piece; when that piece is a
    // chord inside a bend of four or more joints (short leads, long chords) the scan
    // started inside the arc. Before decision 212 the least-squares sub-range anchored at
    // that start was accepted on a joint-only mesh (its first arm, the previous chord, ends
    // on the circle) and the arc was chopped there, while on a subdivided mesh the
    // sub-vertex failed the mesh-segment test and the whole arc was found from its real
    // start: the two discretisations disagreed. Rule: a chord arm whose far joint is itself
    // absorbable is no arm (the range is not maximal at its start), so the scan passes and
    // finds the whole arc from its real start at every discretisation and whatever the
    // loop's start vertex. (i) a bar with TANGENT leads (6 um) along a 250 um bend of six 5
    // deg chords (21.8 um, the longest pieces): one 7-joint tangent bend per side; (ii) a
    // bar whose 6 um leads meet a 120 um bend of seven unequal chords (8-12 deg, 16.7-25.1
    // um) at a 20 deg kink (below the cap, same sign; neither tangent nor a chord: not
    // absorbable): the kink joints stay corners and the six interior joints are one
    // least-squares bend per side; (iii) the same with ten chords whose longest lies
    // mid-arc (both scan directions chopped it).
    auto Rotated = [](const std::vector<Point2> &points, std::size_t start)
    {
      std::vector<Point2> result(points.begin() + static_cast<std::ptrdiff_t>(start),
                                 points.end());
      result.insert(result.end(), points.begin(),
                    points.begin() + static_cast<std::ptrdiff_t>(start));
      return result;
    };
    auto KinkedBar = [](double width, double radius,
                        const std::vector<double> &chord_degrees, double lead,
                        double kink_degrees)
    {
      const double deg = std::acos(-1.0) / 180.0;
      std::vector<Point2> centreline, normals;
      double angle = 0.0;  // tangent angle from +x along the left-turning arc
      std::vector<double> angles = {0.0};
      for (const double chord : chord_degrees)
      {
        angle += chord * deg;
        angles.push_back(angle);
      }
      // Leads rotated by -kink from the tangents at both ends (the junction turns by kink +
      // half the chord angle, in the arc's sense).
      const double a_in = -kink_degrees * deg, a_out = angles.back() + kink_degrees * deg;
      const Point2 p0 = {0.0, 0.0};
      centreline.push_back({p0[0] - lead * std::cos(a_in), p0[1] - lead * std::sin(a_in)});
      normals.push_back({-std::sin(a_in), std::cos(a_in)});
      for (const double a : angles)
      {
        centreline.push_back({radius * std::sin(a), radius - radius * std::cos(a)});
        normals.push_back({-std::sin(a), std::cos(a)});
      }
      const Point2 end = centreline.back();
      centreline.push_back(
          {end[0] + lead * std::cos(a_out), end[1] + lead * std::sin(a_out)});
      normals.push_back({-std::sin(a_out), std::cos(a_out)});
      return BarAroundCentreline(centreline, width, 0.0, normals);
    };
    struct Loop
    {
      std::string name;
      std::vector<Point2> points;
      double radius;
      std::size_t joints_per_arc;
      int bend_vertices;  // absorbed CORNER joints (a 2.5 deg tangent end joint is noise)
      int corners;
    };
    const double width = 8.0;
    for (const Loop &loop :
         {Loop{"tangent leads", ArcBar(width, 250.0, 30.0, 5.0, 6.0), 250.0, 7, 10, 4},
          Loop{"kinked leads",
               KinkedBar(width, 120.0, {8.0, 10.0, 12.0, 10.0, 8.0, 10.0, 11.0}, 6.0, 20.0),
               120.0, 6, 12, 8},
          // The longest chord in the MIDDLE of a 9-joint bend: both scan directions start
          // inside it and both chopped it (5 + 4 joints) before the rule.
          Loop{"kinked leads, longest chord mid-arc",
               KinkedBar(width, 120.0, {6.0, 7.0, 8.0, 7.0, 6.0, 10.0, 6.0, 7.0, 8.0, 7.0},
                         6.0, 20.0),
               120.0, 9, 18, 8}})
    {
      INFO(loop.name);
      const Reading plain = Read({{loop.points, 0, joint_only}});
      // One whole bend per side (radius +/- width / 2), all interior joints absorbed.
      REQUIRE(plain.arcs.size() == 2);
      for (const auto &arc : plain.arcs)
      {
        CHECK(std::get<0>(arc) == "Bend");
        CHECK(std::get<2>(arc) == loop.joints_per_arc);
      }
      CHECK(std::get<1>(plain.arcs[0]) ==
            std::round((loop.radius - 0.5 * width) / R * 1.0e6));
      CHECK(std::get<1>(plain.arcs[1]) ==
            std::round((loop.radius + 0.5 * width) / R * 1.0e6));
      CHECK(plain.vertex_types.at("BendVertex") == loop.bend_vertices);
      CHECK(plain.corners == loop.corners);
      CheckSame(plain, Read({{Subdivide(loop.points, 4.0), 0, joint_only}}));
      CheckSame(plain, Read({{SubdivideIrregularly(loop.points), 0, joint_only}}));
      CheckSame(plain, Read({{loop.points, 0, 1.0}}));
      // Start-vertex invariance: the loop's point list rotated to start at every vertex.
      for (std::size_t start = 1; start < loop.points.size(); start++)
      {
        INFO("start vertex " << start);
        CheckSame(plain, Read({{Rotated(loop.points, start), 0, joint_only}}));
        CheckSame(plain, Read({{Rotated(loop.points, start), 0, 4.0}}));
      }
    }
  }
}

TEST_CASE("SurfaceResponseIdentificationSymmetricCurvedStackOrientation",
          "[surfaceresponseidentification][Serial]")
{
  // Decision 214 (i): the Convexity of a curved stack whose cross-section is its own mirror
  // image (chirality 0: two or three equal traces along one bend) is read on side 0, and
  // side 0 of a symmetric cross-section was the member with the lowest chain id, a
  // coordinate-sorted (frame-dependent) numbering: the subdivision gate's rotation variant
  // flipped the k4 / k6 curved stacks Convex <-> Concave. The bend sense now orients a
  // curved symmetric cross-section (side 0 = the outermost edge of the bend), so the
  // signature is the same under every rigid rotation and translation of the layout and
  // under any numbering of the perimeter (here: the segment order reversed, which gives the
  // outer trace the lower chain ids as a rotation does through the canonical numbering).
  const double R = 2.0;
  struct Stack
  {
    std::string name;
    int traces;
    double rho_over_R;
  };
  struct Reading
  {
    std::string signature;
    int chirality;
    std::string convexity;
    bool first_side_outermost;
  };
  for (const Stack &stack : {Stack{"k4-rho3", 2, 3.0}, Stack{"k4-rho8", 2, 8.0},
                             Stack{"k6-rho3", 3, 3.0}, Stack{"k6-rho8", 3, 8.0}})
  {
    INFO(stack.name);
    // 3 um traces 2 um apart along a 90 deg bend with 12 um leads, the innermost trace's
    // inner edge at rho (as the gate's stack-curved-k4 / k6 layouts); ArcBar offsets run
    // toward the bend centre (0, centre), so the further traces sit outside.
    const double centre = stack.rho_over_R * R + 1.5;
    std::vector<std::vector<Point2>> traces;
    for (int t = 0; t < stack.traces; t++)
    {
      traces.push_back(ArcBar(3.0, centre, 90.0, 5.0, 12.0, -5.0 * t));
    }
    std::vector<double> wanted_offsets;
    for (int t = 0; t < stack.traces; t++)
    {
      wanted_offsets.push_back(2.5 * t);
      wanted_offsets.push_back(2.5 * t + 1.5);
    }
    auto Read = [&](double angle_degrees, const Point2 &shift, bool reverse_segments)
    {
      std::vector<LoopSpec> loops;
      for (const auto &trace : traces)
      {
        loops.push_back({RotatedAndShifted(trace, angle_degrees, shift), 0, 1.0});
      }
      const auto input = MakeInput(loops, R, reverse_segments);
      const auto result = IdentifyMetalPerimeter(input);
      CheckPartition(input, result);
      const Point2 bend_centre =
          RotatedAndShifted({{0.0, centre}}, angle_degrees, shift).front();
      std::optional<Reading> reading;
      int curved_stacks = 0;
      for (const auto &feature : result.features)
      {
        if (feature.type != "CurvedParallelEdgeCluster")
        {
          continue;
        }
        curved_stacks++;
        std::vector<double> offsets;
        for (const auto &edge : feature.signature["Edges"])
        {
          offsets.push_back(edge["OffsetOverR"].get<double>());
        }
        REQUIRE(offsets.size() == wanted_offsets.size());
        for (std::size_t i = 0; i < offsets.size(); i++)
        {
          CHECK_THAT(offsets[i], WithinAbs(wanted_offsets[i], 0.02));
        }
        CHECK_THAT(feature.signature["RadiusOverR"].get<double>(),
                   WithinAbs(stack.rho_over_R, 0.15));
        // Mean distance of every side's portions from the bend centre: side 0 (the first
        // signature edge, chirality 0) must be the outermost edge.
        std::map<int, std::pair<double, int>> side_radius;
        for (const auto &portion : feature.portions)
        {
          const auto &p0 = input.segments[portion.segment].p0;
          side_radius[portion.side].first +=
              std::hypot(p0[0] - bend_centre[0], p0[1] - bend_centre[1]);
          side_radius[portion.side].second++;
        }
        REQUIRE(side_radius.size() == wanted_offsets.size());
        const double first_radius = side_radius.at(0).first / side_radius.at(0).second;
        bool outermost = true;
        for (const auto &[side, sum_count] : side_radius)
        {
          outermost = outermost && (side == 0 || sum_count.first / sum_count.second <
                                                     first_radius - 1.0);
        }
        reading = Reading{feature.signature.dump(), feature.chirality,
                          feature.signature["Convexity"].get<std::string>(), outermost};
      }
      REQUIRE(curved_stacks == 1);
      return *reading;
    };
    const Reading base = Read(0.0, {0.0, 0.0}, false);
    CHECK(base.chirality == 0);
    CHECK(base.convexity == "Convex");  // the outermost trace edge bends around its metal
    CHECK(base.first_side_outermost);
    // The gate's 37 deg and other angles (one per quadrant, a half turn), the rotated
    // layouts shifted off the origin; the perimeter numbering reversed in every frame.
    for (const double angle : {0.0, 37.0, 101.0, 180.0, 253.7})
    {
      const Point2 shift = angle == 0.0 ? Point2{0.0, 0.0} : Point2{-7.25, 3.5};
      for (const bool reverse_segments : {false, true})
      {
        INFO("angle " << angle << " shift (" << shift[0] << ", " << shift[1]
                      << ") reversed segment order " << reverse_segments);
        const Reading other = Read(angle, shift, reverse_segments);
        CHECK(other.signature == base.signature);
        CHECK(other.chirality == base.chirality);
        CHECK(other.convexity == base.convexity);
        CHECK(other.first_side_outermost);
      }
    }
  }
}

TEST_CASE("SurfaceResponseIdentificationSymmetricCurvedPairOrientation",
          "[surfaceresponseidentification][Serial]")
{
  // Decision 214 (i) at k = 2: a curved pair whose two edges are alike (chirality 0) reads
  // the Convexity of its outermost edge in every frame — a symmetric curved gap is Concave
  // (the outer edge's metal lies outside the bend), a symmetric curved strip Convex (the
  // outer edge's metal lies inside it) — under rigid rotations and translations of the
  // layout and under a reversed perimeter numbering. On the (chain, position) key the same
  // strip read Concave or Convex with the chain order (the identity gate's
  // syn-arc-r5-step20 carried both readings of one arc).
  const double R = 2.0;
  struct Pair
  {
    std::string name;
    std::string type;
    std::string convexity;
    double centre;  // the bend centre (0, centre)
    std::vector<std::vector<Point2>> loops;
  };
  // Gap: two 6 um (3 R) bars of one conductor 2 um (1 R) apart along a 90 deg bend, the
  // inner bar's inner edge at 3 R, 12 um leads (ArcBar offsets run toward the bend centre,
  // so the second bar sits outside). Strip: the Convexity test's 3 um bar along a 5 um
  // bend.
  const double gap_centre = 3.0 * R + 3.0;
  const std::vector<Pair> pairs = {
      {"gap",
       "CurvedSameConductorGap",
       "Concave",
       gap_centre,
       {ArcBar(6.0, gap_centre, 90.0, 5.0, 12.0, 0.0),
        ArcBar(6.0, gap_centre, 90.0, 5.0, 12.0, -(3.0 + 2.0 + 3.0))}},
      {"strip", "CurvedSameConductorStrip", "Convex", 5.0, {ArcBar(3.0, 5.0, 90.0, 5.0)}}};
  struct Reading
  {
    std::string signature;
    int chirality;
    std::string convexity;
    bool first_side_outermost;
  };
  for (const Pair &pair : pairs)
  {
    INFO(pair.name);
    auto Read = [&](double angle_degrees, const Point2 &shift, bool reverse_segments)
    {
      std::vector<LoopSpec> loops;
      for (const auto &loop : pair.loops)
      {
        loops.push_back({RotatedAndShifted(loop, angle_degrees, shift), 0, 1.0});
      }
      const auto input = MakeInput(loops, R, reverse_segments);
      const auto result = IdentifyMetalPerimeter(input);
      CheckPartition(input, result);
      const Point2 bend_centre =
          RotatedAndShifted({{0.0, pair.centre}}, angle_degrees, shift).front();
      std::optional<Reading> reading;
      int curved_pairs = 0;
      for (const auto &feature : result.features)
      {
        if (feature.type != pair.type)
        {
          continue;
        }
        curved_pairs++;
        // Mean distance of every side's portions from the bend centre: side 0 (the first
        // signature edge, chirality 0) must be the outermost edge.
        std::map<int, std::pair<double, int>> side_radius;
        for (const auto &portion : feature.portions)
        {
          const auto &p0 = input.segments[portion.segment].p0;
          side_radius[portion.side].first +=
              std::hypot(p0[0] - bend_centre[0], p0[1] - bend_centre[1]);
          side_radius[portion.side].second++;
        }
        REQUIRE(side_radius.size() == 2);
        const double first_radius = side_radius.at(0).first / side_radius.at(0).second;
        const double other_radius = side_radius.at(1).first / side_radius.at(1).second;
        reading = Reading{feature.signature.dump(), feature.chirality,
                          feature.signature["Convexity"].get<std::string>(),
                          first_radius > other_radius + 1.0};
      }
      REQUIRE(curved_pairs == 1);
      return *reading;
    };
    const Reading base = Read(0.0, {0.0, 0.0}, false);
    CHECK(base.chirality == 0);
    CHECK(base.convexity == pair.convexity);
    CHECK(base.first_side_outermost);
    for (const double angle : {0.0, 37.0, 101.0, 180.0, 253.7})
    {
      const Point2 shift = angle == 0.0 ? Point2{0.0, 0.0} : Point2{-7.25, 3.5};
      for (const bool reverse_segments : {false, true})
      {
        INFO("angle " << angle << " shift (" << shift[0] << ", " << shift[1]
                      << ") reversed segment order " << reverse_segments);
        const Reading other = Read(angle, shift, reverse_segments);
        CHECK(other.signature == base.signature);
        CHECK(other.chirality == base.chirality);
        CHECK(other.convexity == base.convexity);
        CHECK(other.first_side_outermost);
      }
    }
  }
}

namespace
{

// Feature portions shorter than the signature parameter tolerance that continue into NO
// portion of the same feature at either end (on the same segment or across the shared
// vertex of the neighbouring segment): the slivers decision 222 forbids. A sub-tolerance
// MESH SEGMENT inside a long edge of one feature is not one.
struct SubTolerancePortion
{
  std::size_t segment;
  double s0, s1;
  int feature;
};

std::vector<SubTolerancePortion> SubTolerancePortions(const IdentificationInput &input,
                                                      const IdentificationResult &result)
{
  const double tolerance = kSignatureParameterToleranceOverRadius * result.radius;
  const double eps = 1.0e-9;
  std::map<std::array<double, 3>, std::vector<std::size_t>> segments_at;
  for (std::size_t i = 0; i < result.segments.size(); i++)
  {
    for (const auto &end : result.segments[i].key)
    {
      segments_at[end].push_back(i);
    }
  }
  // Does the feature own a portion of segment j reaching the given endpoint of j?
  auto ReachesEnd = [&](std::size_t j, const std::array<double, 3> &end, int feature)
  {
    const auto &segment = result.segments[j];
    for (const auto &portion : segment.portions)
    {
      if (static_cast<int>(portion[2]) != feature)
      {
        continue;
      }
      if ((end == segment.key[0] && portion[0] <= eps) ||
          (end == segment.key[1] && portion[1] >= segment.length - eps))
      {
        return true;
      }
    }
    return false;
  };
  std::vector<SubTolerancePortion> slivers;
  for (std::size_t i = 0; i < result.segments.size(); i++)
  {
    const auto &segment = result.segments[i];
    for (const auto &portion : segment.portions)
    {
      const double s0 = portion[0], s1 = portion[1];
      const int feature = static_cast<int>(portion[2]);
      if (s1 - s0 >= tolerance)
      {
        continue;
      }
      bool continues = false;
      for (const auto &other : segment.portions)
      {
        if (static_cast<int>(other[2]) == feature && &other != &portion &&
            (std::abs(other[1] - s0) <= eps || std::abs(other[0] - s1) <= eps))
        {
          continues = true;
        }
      }
      for (int end = 0; end < 2 && !continues; end++)
      {
        if ((end == 0 && s0 > eps) || (end == 1 && s1 < segment.length - eps))
        {
          continue;
        }
        for (const std::size_t j : segments_at.at(segment.key[end]))
        {
          if (j != i && !input.segments[j].truncation &&
              ReachesEnd(j, segment.key[end], feature))
          {
            continues = true;
          }
        }
      }
      if (!continues)
      {
        slivers.push_back({i, s0, s1, feature});
      }
    }
  }
  return slivers;
}

}  // namespace

TEST_CASE("SurfaceResponseIdentificationSubTolerancePortions",
          "[surfaceresponseidentification][Serial]")
{
  // Supervisor decision 222 (the S1p stage-1 window): NO portion shorter than the signature
  // parameter tolerance 1e-3 R exists — a remainder below it joins its adjacent portion on
  // the chain (the longer neighbour; ties: the one before it), replacing the former
  // knife-edge at the 1e-6 R signature grid. The reproducer: a near-parallel 4-edge stack
  // whose members end on a truncation cut with a tilt of ~1e-6 rad between the two bodies,
  // so that the perpendicular foot of one body's cut end lies a few 1e-6 um along the other
  // body's edges and the stack's claim leaves a remainder > 1e-6 R and < 1e-3 R at the cut;
  // on S1p that remainder belonged to the 3-edge stack whose real end lay 64 um away (a
  // sample placed on it pointed its lateral axis along the segment: the placement aborted)
  // and its twin read as a 2e-6 um IsolatedEdge.
  const double R = 2.0;
  const double tolerance = kSignatureParameterToleranceOverRadius * R;
  auto EdgeCount = [](const IdentifiedFeature &f)
  { return f.signature.contains("Edges") ? f.signature["Edges"].size() : 0; };
  SECTION("near-parallel stack members cut by a truncation line")
  {
    // Body A (x in [0, 2]) and body B (x in [4, 5] below y = -20, [4, 9] above) with all
    // four edges exactly parallel, tilted by theta against the normal of the cut y = -60
    // (the S1p configuration: the window set's edges are mutually parallel, the window wall
    // is not perpendicular to them), so that the foot of every member's cut end on its
    // right-hand neighbour lies 2 theta along it (3e-6 um = 1.5e-6 R: above the signature
    // grid, far below the tolerance) and the stack's claim on the inner members starts
    // there. Edges 0 / 2 / 4 / 5 are a 4-edge stack below y = -20, 0 / 2 / 4 a 3-edge stack
    // above it (the x = 9 edge is 5 um = 2.5 R from x = 4). Edge 0 of each loop (its
    // bottom) is the truncation cut y = -60.
    const double theta = 1.5e-6;
    const std::vector<Point2> body_b = {{4.0, -60.0},
                                        {5.0, -60.0},
                                        {5.0 + 40.0 * theta, -20.0},
                                        {9.0 + 40.0 * theta, -20.0},
                                        {9.0 + 60.0 * theta, 0.0},
                                        {4.0 + 60.0 * theta, 0.0}};
    const std::vector<Point2> body_a = {
        {0.0, -60.0}, {2.0, -60.0}, {2.0 + 60.0 * theta, 0.0}, {60.0 * theta, 0.0}};
    const auto input =
        MakeInput({{body_a, 0, 1.0, 0.0, {0}}, {body_b, 0, 1.0, 0.0, {0}}}, R);
    const auto result = IdentifyMetalPerimeter(input);
    CheckPartition(input, result);
    CHECK(result.same_priority_claim_overlaps == 0);
    // The cut itself: both bottom edges excluded as TruncationCut, their vertices cuts.
    int truncation_segments = 0, truncation_vertices = 0;
    for (const auto &segment : result.segments)
    {
      truncation_segments +=
          segment.exclusion && segment.exclusion->first == "TruncationCut" ? 1 : 0;
    }
    for (const auto &vertex : result.vertices)
    {
      truncation_vertices += vertex.type == "TruncationCut" ? 1 : 0;
    }
    CHECK(truncation_segments == 3);  // A's 2 um bottom in 2 pieces, B's 1 um bottom
    CHECK(truncation_vertices == 5);  // their 4 corners and the midpoint of A's bottom
    if (std::getenv("PALACE_IDENTIFICATION_TEST_LOG"))
    {
      for (std::size_t i = 0; i < input.segments.size(); i++)
      {
        const auto &segment = input.segments[i];
        if (std::min(segment.p0[1], segment.p1[1]) > -59.0 || segment.truncation)
        {
          continue;
        }
        std::cout << "segment " << i << " (" << segment.p0[0] << ", " << segment.p0[1]
                  << ") -> (" << segment.p1[0] << ", " << segment.p1[1] << ")\n";
        for (const auto &portion : result.segments[i].portions)
        {
          const auto &feature = result.features[static_cast<int>(portion[2])];
          std::cout << "  [" << std::setprecision(12) << portion[0] << ", " << portion[1]
                    << "] feature " << feature.id << " " << feature.type << " edges "
                    << (feature.signature.contains("Edges")
                            ? feature.signature["Edges"].size()
                            : 0)
                    << " length " << feature.length << "\n";
        }
      }
    }
    const IdentifiedFeature *stack4 = nullptr, *stack3 = nullptr;
    for (const auto &feature : result.features)
    {
      if (feature.type == "ParallelEdgeCluster" && EdgeCount(feature) == 4)
      {
        CHECK(stack4 == nullptr);
        stack4 = &feature;
      }
      else if (feature.type == "ParallelEdgeCluster" && EdgeCount(feature) == 3)
      {
        CHECK(stack3 == nullptr);
        stack3 = &feature;
      }
    }
    REQUIRE(stack4 != nullptr);
    REQUIRE(stack3 != nullptr);
    CHECK(stack4->exact_parameters);
    CHECK(stack3->exact_parameters);
    // The rule: no sub-tolerance sliver anywhere, no feature shorter than the tolerance.
    const auto slivers = SubTolerancePortions(input, result);
    for (const auto &sliver : slivers)
    {
      const auto &p0 = input.segments[sliver.segment].p0;
      INFO("sliver on segment " << sliver.segment << " at (" << p0[0] << ", " << p0[1]
                                << ") [" << sliver.s0 << ", " << sliver.s1 << "] feature "
                                << sliver.feature << " "
                                << result.features[sliver.feature].type);
      CHECK(false);
    }
    CHECK(slivers.empty());
    for (const auto &feature : result.features)
    {
      INFO("feature " << feature.id << " " << feature.type);
      CHECK((feature.length >= tolerance || feature.length == 0.0));
    }
    // The 3-edge stack owns nothing at the cut: every portion of it lies above its stack
    // end at y = -20 less the cluster reach there (2R).
    for (const auto &portion : stack3->portions)
    {
      const auto &segment = input.segments[portion.segment];
      INFO("3-edge stack portion on segment " << portion.segment << " at y "
                                              << std::min(segment.p0[1], segment.p1[1]));
      CHECK(std::min(segment.p0[1], segment.p1[1]) > -20.0 - 2.0 * R);
    }
    // The 4-edge stack reaches the cut on every member: its claim starts at the cut (the
    // sub-tolerance remainders joined it).
    std::map<int, double> side_start;
    for (const auto &portion : stack4->portions)
    {
      const auto &segment = input.segments[portion.segment];
      const double y_low = std::min(segment.p0[1], segment.p1[1]);
      auto it = side_start.find(portion.side);
      side_start[portion.side] =
          it == side_start.end() ? y_low : std::min(it->second, y_low);
    }
    REQUIRE(side_start.size() == 4);
    for (const auto &[side, y_low] : side_start)
    {
      INFO("side " << side);
      CHECK_THAT(y_low, WithinAbs(-60.0, 1.0e-9));
    }
    // The census: the slivers of the former rule (two per body-B edge: the feet of A's two
    // cut ends) joined; nothing left without a neighbour.
    CHECK(result.sub_tolerance_portions.count >= 2);
    CHECK(result.sub_tolerance_portions.max_length < tolerance);
    CHECK(result.sub_tolerance_portions.max_length > kSignatureLengthQuantumOverRadius * R);
    CHECK(result.sub_tolerance_portions.isolated == 0);
  }
  SECTION("a sub-tolerance mesh segment inside a long edge is not a sliver")
  {
    // Two exactly parallel strips with a 1 nm collinear mesh segment (5e-4 R) inside the
    // trace's right edge: the stack's portion on that segment continues into its portions
    // on both neighbouring segments, so the rule has nothing to join.
    const std::vector<Point2> trace = {{0.0, -60.0},          {2.0, -60.0}, {2.0, -30.0},
                                       {2.0, -30.0 + 1.0e-3}, {2.0, 0.0},   {0.0, 0.0}};
    const auto input =
        MakeInput({{trace, 0, 1.0}, {Rectangle(4.0, -60.0, 5.0, 0.0), 0, 1.0}}, R);
    const auto result = IdentifyMetalPerimeter(input);
    CheckPartition(input, result);
    CHECK(SubTolerancePortions(input, result).empty());
    CHECK(result.sub_tolerance_portions.count == 0);
    CHECK(result.sub_tolerance_portions.isolated == 0);
    bool found = false;
    for (std::size_t i = 0; i < input.segments.size(); i++)
    {
      const double length = result.segments[i].length;
      if (length < tolerance)
      {
        found = true;
        REQUIRE(result.segments[i].portions.size() == 1);
        CHECK(result.features[static_cast<int>(result.segments[i].portions[0][2])].type ==
              "ParallelEdgeCluster");
      }
    }
    CHECK(found);
  }
}

TEST_CASE("SurfaceResponseIdentificationStackPieceInsideCluster",
          "[surfaceresponseidentification][Serial]")
{
  // Supervisor decision 224 (the S1p 41-edge loop end): pair / stack stretches that exist
  // only because a cluster's claim boundary cut them are absorbed by that cluster — the one
  // exception to "pairs / stacks are never absorbed" (the spatial coupon's volume and the
  // stack's translational patches would otherwise correct the same surface twice): a
  // stretch of one cross-section's claims, within the cluster ball radius R of the
  // cluster's claims, bounded at both ends by the SAME cluster (two-sided) or adjacent to
  // it at one end and continuing a larger stack at the other (stack-end recomposition).
  // Never between two different clusters, never a genuine pair whose far end is free or a
  // bend. Reproducer: a square loop
  // wire (2 um) attached to the left ground around a hole, with a ground edge 2 um to the
  // right of its right side, so that the loop's right side and the ground edge form a
  // 3-edge stack (0 / 2 / 4 um) between the loop's top and bottom corner clusters; the hole
  // is short enough for the two corner groups to be ONE cluster and for the stack stretch
  // between its claims to be shorter than 2R.
  const double R = 2.0;
  auto Layout = [&](double half_hole)
  {
    // A free square ring wire (outer x in [-2, 8], |y| <= half_hole + 2; hole x in [0, 6],
    // |y| < half_hole, clockwise = metal outside), a straight ground edge at x = 10 and a
    // serrated ground to the left (teeth 1 um wide every 2 um reaching x = -3.5, 1.5 um
    // from the ring): the teeth's corners are events all along the ring's left side, so
    // the ring's top and bottom corner groups are ONE cluster that wraps around the loop,
    // while the ring's right side (edges x = 6 / 8) and the ground edge x = 10.5 form a
    // 3-edge stack (offsets 0 / 1 / 2.25 R) between that cluster's claims.
    const double h = half_hole, H = half_hole + 2.0;
    std::vector<Point2> left = {{-20.0, -20.0}, {-4.0, -20.0}};
    for (double y = -H - 1.0; y < H + 1.0; y += 2.0)
    {
      left.push_back({-4.0, y});
      left.push_back({-3.5, y});
      left.push_back({-3.5, y + 1.0});
      left.push_back({-4.0, y + 1.0});
    }
    left.push_back({-4.0, 20.0});
    left.push_back({-20.0, 20.0});
    std::vector<Point2> ring = {{-2.0, -H}, {8.0, -H}, {8.0, H}, {-2.0, H}};
    std::vector<Point2> hole = {{0.0, -h}, {0.0, h}, {6.0, h}, {6.0, -h}};
    // The ground edge at x = 10.5: 2.5 um from the ring (4.5 um from the hole edge, clear
    // of the exactly-2R knife-edge a ground at x = 10 would sit on).
    std::vector<Point2> ground = Rectangle(10.5, -20.0, 20.0, 20.0);
    return MakeInput({{left, 0, 1.0}, {ring, 0, 1.0}, {hole, 0, 1.0}, {ground, 0, 1.0}}, R);
  };
  // The stack length on the three lead edges between the cluster's claims (x = 6 / 8 /
  // 10.5), the clusters, and whether every lead-edge portion between |y| < half_hole - 2
  // belongs to the single cluster.
  struct Reading
  {
    int clusters = 0;
    double stack_length = 0.0;
    std::size_t stack_features = 0;
    bool leads_owned_by_cluster = true;
  };
  auto Read = [&](const IdentificationInput &input, const IdentificationResult &result,
                  double half_hole)
  {
    Reading reading;
    for (const auto &feature : result.features)
    {
      if (feature.type == "SpatialEdgeCluster")
      {
        reading.clusters++;
      }
      else if (feature.type == "ParallelEdgeCluster")
      {
        reading.stack_features++;
        reading.stack_length += feature.length;
      }
      for (const auto &portion : feature.portions)
      {
        const auto &segment = input.segments[portion.segment];
        const bool lead_edge = std::abs(segment.p0[0] - segment.p1[0]) < 1.0e-9 &&
                               (std::abs(segment.p0[0] - 6.0) < 1.0e-9 ||
                                std::abs(segment.p0[0] - 8.0) < 1.0e-9 ||
                                std::abs(segment.p0[0] - 10.5) < 1.0e-9);
        if (lead_edge && std::abs(segment.p0[1]) < half_hole - 2.0 &&
            std::abs(segment.p1[1]) < half_hole - 2.0 &&
            feature.type != "SpatialEdgeCluster")
        {
          reading.leads_owned_by_cluster = false;
        }
      }
    }
    return reading;
  };
  SECTION("a stack stretch shorter than 2R between the claims of one cluster is absorbed")
  {
    // half_hole 7: the cluster's claims leave a 3.07 um stretch (1.54 R) of the 3-edge
    // stack on every lead; one cluster bounds it on both sides and every point lies within
    // R of one of its claims (two-sided). Before decision 224 the stretch stayed a
    // ParallelEdgeCluster (9.2 um over the three edges) inside the cluster's extent; now
    // the cluster owns the leads entirely.
    const auto input = Layout(7.0);
    const auto result = IdentifyMetalPerimeter(input);
    CheckPartition(input, result);
    CHECK(SubTolerancePortions(input, result).empty());
    const Reading reading = Read(input, result, 7.0);
    CHECK(reading.clusters == 1);
    CHECK(reading.stack_features == 0);
    CHECK(reading.leads_owned_by_cluster);
    // Pass 1 absorbs the one stretch already bounded by the cluster's claims at both ends
    // (the ring edge x = 8, 3.0718 um); the stack recomposed without it leaves the hole
    // edge x = 6 and the ground edge x = 10.5 (4.5 um apart: no pair) as single-edge
    // remainders, which the ordinary extension absorbs in pass 2.
    CHECK(result.extension.translational_pieces == 1);
    CHECK(result.extension.translational_two_sided == 1);
    CHECK_THAT(result.extension.translational_length, WithinAbs(3.0718, 0.001));
    CHECK_THAT(result.extension.translational_two_sided_length, WithinAbs(3.0718, 0.001));
    CHECK_THAT(result.extension.translational_max_length, WithinAbs(3.0718, 0.001));
    CHECK(result.extension.passes >= 2);
    // The stretch is cut into its three runs' intervals, all into the one cluster: its
    // portions on the leads are contiguous stretches from the top claim to the bottom one.
    for (const auto &feature : result.features)
    {
      if (feature.type != "SpatialEdgeCluster")
      {
        continue;
      }
      std::map<int, std::set<int>> stretches_on_lead;
      for (const auto &portion : feature.portions)
      {
        const auto &segment = input.segments[portion.segment];
        if (std::abs(segment.p0[0] - segment.p1[0]) < 1.0e-9 &&
            (std::abs(segment.p0[0] - 6.0) < 1.0e-9 ||
             std::abs(segment.p0[0] - 8.0) < 1.0e-9))
        {
          stretches_on_lead[static_cast<int>(std::lround(segment.p0[0]))].insert(
              portion.stretch);
        }
      }
      REQUIRE(stretches_on_lead.size() == 2);
      for (const auto &[x, stretches] : stretches_on_lead)
      {
        INFO("lead x = " << x);
        CHECK(stretches.size() == 1);
      }
    }
  }
  SECTION("a stretch reaching beyond the ball radius of the cluster's claims stays a stack")
  {
    // half_hole 9: a 7.07 um stretch (3.54 R) on every lead, adjacent to the cluster at
    // both ends but with its middle 3.5 um from either claim: the identification leaves it
    // to the stack (the placement's ownership check judges it against the coupon volume).
    const auto input = Layout(9.0);
    const auto result = IdentifyMetalPerimeter(input);
    CheckPartition(input, result);
    const Reading reading = Read(input, result, 9.0);
    CHECK(reading.clusters == 1);
    CHECK(reading.stack_features >= 1);
    CHECK_THAT(reading.stack_length, WithinAbs(3.0 * 7.0718, 0.01));
    CHECK(!reading.leads_owned_by_cluster);
    CHECK(result.extension.translational_pieces == 0);
  }
  SECTION("a stretch between the claims of two different clusters stays with the stack")
  {
    // The ring attached to a plain left ground (no teeth): the top and bottom corner groups
    // are two clusters, and the 1.07 um stretch of the 3-edge stack between their claims is
    // bounded by claims of DIFFERENT clusters: a genuine stack between two clusters, not a
    // claim-boundary artefact of one of them — never absorbed (the placement's
    // spatial-vs-spatial check covers overlapping coupon boxes); half_hole 6.
    const double h = 6.0, H = 8.0;
    std::vector<Point2> metal = {{-20.0, -20.0}, {-2.0, -20.0}, {-2.0, -H},
                                 {8.0, -H},      {8.0, H},      {-2.0, H},
                                 {-2.0, 20.0},   {-20.0, 20.0}};
    std::vector<Point2> hole = {{0.0, -h}, {0.0, h}, {6.0, h}, {6.0, -h}};
    const auto input = MakeInput(
        {{metal, 0, 1.0}, {hole, 0, 1.0}, {Rectangle(10.5, -20.0, 20.0, 20.0), 0, 1.0}}, R);
    const auto result = IdentifyMetalPerimeter(input);
    CheckPartition(input, result);
    const Reading reading = Read(input, result, h);
    CHECK(reading.clusters == 2);
    CHECK(reading.stack_features >= 1);
    CHECK_THAT(reading.stack_length, WithinAbs(3.0 * 1.0718, 0.01));
    CHECK(!reading.leads_owned_by_cluster);
    CHECK(result.extension.translational_pieces == 0);
  }
  SECTION("a stack-end recomposition piece next to a cluster is absorbed, the stack kept")
  {
    // Decision 230 (ii), class StackEndRecomposition (the gate's stack-k3-1p5-3 /
    // stack-k4-1-1p5-3 layouts): a trace bar between grounds whose edges reach the
    // truncation box (no ground corners). At each trace end the cluster's claims reach
    // unequally far along the members (the ground within 2R of the end corners' cores
    // farther than the trace edges), so between the claim ends the stack is recomposed as a
    // SMALLER cross-section continuing the larger stack toward the cluster: adjacent to the
    // cluster at one end, continuing the k-edge stack at the other, within the 1 R ball —
    // absorbed; the k-edge stack itself (reaching far beyond the ball) is kept.
    auto GroundBar = [&](double y0, double y1)
    {
      LoopSpec ground{Rectangle(-40.0, y0, 40.0, y1), 0, 1.0};
      ground.truncation_edges = {1, 3};  // the x = +-40 sides lie on the window walls
      return ground;
    };
    auto Count = [&](const IdentificationResult &result, const std::string &type)
    {
      int count = 0;
      for (const auto &feature : result.features)
      {
        count += feature.type == type ? 1 : 0;
      }
      return count;
    };
    {
      // k = 3: ground | 1.5 um gap | 3 um trace (offsets 0 / 0.75 / 2.25 R). The trace
      // strip between the ground's claim end and the trace's own (0.75 R from the cores vs
      // the trace's edges) continues the 3-edge stack: one recomposition piece per end.
      const auto input = MakeInput(
          {LoopSpec{Rectangle(-30.0, 0.0, 30.0, 3.0), 0, 1.0}, GroundBar(-9.5, -1.5)}, R);
      const auto result = IdentifyMetalPerimeter(input);
      CheckPartition(input, result);
      CHECK(Count(result, "SpatialEdgeCluster") == 2);
      CHECK(Count(result, "ParallelEdgeCluster") >= 1);
      CHECK(Count(result, "SameConductorStrip") == 0);
      CHECK(result.extension.translational_two_sided == 0);
      CHECK(result.extension.translational_pieces -
                result.extension.translational_two_sided ==
            2);
      CHECK_THAT(result.extension.translational_length, WithinAbs(1.3542, 0.001));
      CHECK_THAT(result.extension.translational_max_length, WithinAbs(0.6771, 0.001));
      double stack_length = 0.0;
      for (const auto &feature : result.features)
      {
        stack_length += feature.type == "ParallelEdgeCluster" ? feature.length : 0.0;
      }
      CHECK(stack_length > 3.0 * 40.0);  // the 3 members over most of the 60 um bar
    }
    {
      // k = 4: ground | 1 | trace 1.5 | 3 | ground (offsets 0 / 0.5 / 1.25 / 2.75 R): the
      // 3-edge stack (trace | 3 | ground) then the DifferentConductorGap (trace top |
      // ground) each continue the larger cross-section toward the cluster: three pieces
      // per end over successive passes (the ball grows with the claims), the 4-edge stack
      // kept.
      const auto input = MakeInput({LoopSpec{Rectangle(-30.0, 0.0, 30.0, 1.5), 0, 1.0},
                                    GroundBar(-9.0, -1.0), GroundBar(4.5, 12.5)},
                                   R);
      const auto result = IdentifyMetalPerimeter(input);
      CheckPartition(input, result);
      CHECK(Count(result, "SpatialEdgeCluster") == 2);
      CHECK(Count(result, "ParallelEdgeCluster") >= 1);
      CHECK(Count(result, "DifferentConductorGap") == 0);
      CHECK(result.extension.translational_two_sided == 0);
      CHECK(result.extension.translational_pieces -
                result.extension.translational_two_sided ==
            6);
      CHECK_THAT(result.extension.translational_length, WithinAbs(2.9904, 0.001));
      CHECK_THAT(result.extension.translational_max_length, WithinAbs(0.9593, 0.001));
      CHECK(result.extension.passes >= 3);
    }
  }
  SECTION("a bent strip's halves next to their end clusters are genuine pairs")
  {
    // A 2.5 um strip (1.25 R) of 7.5 um between its own end corners, bent by 8 deg at
    // mid-length: each half is its own cross-section, adjacent to an end-corner cluster and
    // shorter than R beyond it — but its far end is the bend, not a larger stack, so the
    // halves stay SameConductorStrip features (the translational mortar strip test of
    // test-surfaceresponseoperator.cpp relies on it).
    const double half = 1.25, length = 7.5, tilt = std::tan(8.0 * std::acos(-1.0) / 180.0);
    const std::vector<Point2> strip = {{-half, 0.0},
                                       {half, 0.0},
                                       {half, 0.5 * length},
                                       {half + tilt * 0.5 * length, length},
                                       {-half + tilt * 0.5 * length, length},
                                       {-half, 0.5 * length}};
    const auto input = MakeInput({{strip, 0, 0.5}}, R);
    const auto result = IdentifyMetalPerimeter(input);
    CheckPartition(input, result);
    int strips = 0;
    for (const auto &feature : result.features)
    {
      strips += feature.type == "SameConductorStrip" ? 1 : 0;
    }
    CHECK(strips >= 1);
    CHECK(result.extension.translational_pieces == 0);
  }
}

TEST_CASE("SurfaceResponseIdentificationFlippedPlaneFrames",
          "[surfaceresponseidentification][Serial]")
{
  // Decision 266 (flip-chip orientation): a spatial cluster and a vertex feature are framed
  // with their OWN signed process normal (substrate -> vacuum), not with the device-global
  // sign-canonical reference normal. Two planes carry congruent copies of one asymmetric
  // scene (a pad corner near the end of an offset strip: an asymmetric cluster, plain
  // convex corners elsewhere): the upright plane U at z = 0 facing +z, and the flipped
  // plane F at z = 4.8 facing -z (a flip-chip top chip, its substrate above) = U rotated by
  // 180 deg about the line {x = 30, z = 2.4} parallel to y (a proper rigid motion: (x, y,
  // z) -> (60 - x, y, 4.8 - z), vectors (vx, vy, vz) -> (-vx, vy, -vz)).
  const double R = 2.0;
  const double flipped_z = 4.8, shift_x = 60.0;
  const std::vector<Point2> pad = Rectangle(0.0, 0.0, 10.0, 10.0);
  const std::vector<Point2> strip = Rectangle(-3.0, 3.0, -1.0, 20.0);
  auto Image = [&](std::vector<Point2> points)
  {
    for (auto &p : points)
    {
      p = {shift_x - p[0], p[1]};
    }
    std::reverse(points.begin(), points.end());  // keep counter-clockwise in plan view
    return points;
  };
  auto ImageVector = [](const std::array<double, 3> &v)
  { return std::array<double, 3>{-v[0], v[1], -v[2]}; };
  auto ImagePoint = [&](const std::array<double, 3> &p)
  { return std::array<double, 3>{shift_x - p[0], p[1], flipped_z - p[2]}; };
  auto Close = [](const std::array<double, 3> &a, const std::array<double, 3> &b)
  {
    return std::abs(a[0] - b[0]) < 1.0e-9 && std::abs(a[1] - b[1]) < 1.0e-9 &&
           std::abs(a[2] - b[2]) < 1.0e-9;
  };
  auto Handedness = [](const IdentifiedFeature &f)
  {
    const auto &x = f.axes[0], &y = f.axes[1], &w = f.axes[2];
    const std::array<double, 3> c = {x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2],
                                     x[0] * y[1] - x[1] * y[0]};
    return c[0] * w[0] + c[1] * w[1] + c[2] * w[2];
  };
  const std::vector<LoopSpec> upright = {{pad, 0, 1.0}, {strip, 0, 1.0}};
  std::vector<LoopSpec> flipped = {{Image(pad), 0, 1.0}, {Image(strip), 0, 1.0}};
  for (auto &loop : flipped)
  {
    loop.z = flipped_z;
    loop.normal_sign = -1.0;
  }
  std::vector<LoopSpec> both = upright;
  both.insert(both.end(), flipped.begin(), flipped.end());

  const auto upright_input = MakeInput(upright, R);
  const auto upright_alone = IdentifyMetalPerimeter(upright_input);
  CheckPartition(upright_input, upright_alone);
  const auto both_input = MakeInput(both, R);
  const auto result = IdentifyMetalPerimeter(both_input);
  CheckPartition(both_input, result);
  CHECK(result.reference_process_normal == std::array<double, 3>{0.0, 0.0, 1.0});

  auto PlaneOf = [&](const IdentifiedFeature &f)
  {
    REQUIRE(!f.portions.empty());
    return both_input.segments[f.portions.front().segment].p0[2] > 0.5 * flipped_z ? 1 : 0;
  };
  std::vector<const IdentifiedFeature *> upright_features, flipped_features;
  for (const auto &f : result.features)
  {
    (PlaneOf(f) == 1 ? flipped_features : upright_features).push_back(&f);
  }
  REQUIRE(upright_features.size() == upright_alone.features.size());
  REQUIRE(flipped_features.size() == upright_alone.features.size());

  // The upright plane reads exactly as it does alone: keys, chirality, frames.
  {
    auto Key = [](const IdentifiedFeature &f)
    { return std::make_tuple(f.type, f.hash, f.chirality, f.origin, f.axes); };
    std::multiset<decltype(Key(result.features.front()))> alone, together;
    for (const auto &f : upright_alone.features)
    {
      alone.insert(Key(f));
    }
    for (const auto *f : upright_features)
    {
      together.insert(Key(*f));
    }
    CHECK(alone == together);
  }

  // Every spatial / vertex feature of the flipped plane: w = -z (its own vacuum side), the
  // key AND chirality of its congruent upright counterpart (the hash of a cluster signature
  // is handedness-invariant by design — a plan-view mirror image has the same hash with the
  // opposite chirality — so the identity under test is the pair (hash, chirality)), the
  // frame = the rigid image of the upright frame, and the coupon's substrate half-space
  // (support w in [-1.95, 0)) inside the flipped plane's substrate (z > 4.8).
  int clusters_checked = 0, corners_checked = 0, asymmetric_clusters = 0;
  for (const auto *f : flipped_features)
  {
    if (f->type != "SpatialEdgeCluster" && f->type != "ConvexCorner" &&
        f->type != "ConcaveCorner")
    {
      continue;
    }
    INFO("flipped " << f->type << " at (" << f->origin[0] << ", " << f->origin[1] << ")");
    CHECK(Close(f->axes[2], {0.0, 0.0, -1.0}));
    const double coupon_substrate_z = f->origin[2] - 1.95 * f->axes[2][2];
    CHECK(coupon_substrate_z > flipped_z);
    const IdentifiedFeature *counterpart = nullptr;
    for (const auto &u : upright_alone.features)
    {
      if (u.type == f->type && Close(ImagePoint(u.origin), f->origin))
      {
        counterpart = &u;
      }
    }
    REQUIRE(counterpart != nullptr);
    CHECK(counterpart->hash == f->hash);
    CHECK(counterpart->chirality == f->chirality);
    CHECK(Close(counterpart->axes[2], {0.0, 0.0, 1.0}));
    if (f->type == "SpatialEdgeCluster")
    {
      clusters_checked++;
      if (f->chirality != 0)
      {
        // Unique canonical frame: the flipped frame is the image of the upright one, and
        // its handedness is the chirality (right-handed for +1, the mirror frame for -1).
        asymmetric_clusters++;
        CHECK(Close(ImageVector(counterpart->axes[0]), f->axes[0]));
        CHECK(Close(ImageVector(counterpart->axes[1]), f->axes[1]));
        CHECK_THAT(Handedness(*f), WithinAbs(static_cast<double>(f->chirality), 1.0e-12));
        CHECK_THAT(Handedness(*counterpart),
                   WithinAbs(static_cast<double>(f->chirality), 1.0e-12));
      }
    }
    else
    {
      // A corner frame is right-handed (x = the arm from which the other is
      // counterclockwise about the plane's own normal, y = n x x): the flipped corner's
      // frame is the rigid image of the upright one.
      corners_checked++;
      CHECK_THAT(Handedness(*f), WithinAbs(1.0, 1.0e-12));
      CHECK(Close(ImageVector(counterpart->axes[0]), f->axes[0]));
      CHECK(Close(ImageVector(counterpart->axes[1]), f->axes[1]));
    }
  }
  CHECK(clusters_checked >= 2);
  CHECK(asymmetric_clusters >= 1);
  CHECK(corners_checked >= 2);

  // The old reading (the flipped geometry framed with the sign-canonical +z normal) is the
  // plan-view MIRROR of the congruent upright scene: the same hash with the opposite
  // chirality and w = +z (the coupon upside down). The new frames differ from it exactly
  // there.
  {
    std::vector<LoopSpec> mirror = {{Image(pad), 0, 1.0}, {Image(strip), 0, 1.0}};
    const auto mirror_result = IdentifyMetalPerimeter(MakeInput(mirror, R));
    for (const auto *f : flipped_features)
    {
      if (f->type != "SpatialEdgeCluster" || f->chirality == 0)
      {
        continue;
      }
      const IdentifiedFeature *old_reading = nullptr;
      for (const auto &m : mirror_result.features)
      {
        if (m.type == "SpatialEdgeCluster" &&
            std::abs(m.origin[0] - f->origin[0]) < 1.0e-9 &&
            std::abs(m.origin[1] - f->origin[1]) < 1.0e-9)
        {
          old_reading = &m;
        }
      }
      REQUIRE(old_reading != nullptr);
      CHECK(old_reading->hash == f->hash);
      CHECK(old_reading->chirality == -f->chirality);
      CHECK(Close(old_reading->axes[2], {0.0, 0.0, 1.0}));
      CHECK(Close(old_reading->axes[0], f->axes[0]));
      CHECK(Close(old_reading->axes[1], f->axes[1]));  // the in-plane map is the same
    }
  }

  // A cluster whose runs disagree on the sign of their process normal has no common
  // substrate -> vacuum side: the identification fails closed naming the cluster (the strip
  // of the same plane facing -z while the pad faces +z: the asymmetric cluster spans both).
  {
    std::vector<LoopSpec> mixed = {{pad, 0, 1.0}, {strip, 0, 1.0}};
    mixed[1].normal_sign = -1.0;
    CHECK_THROWS_WITH(IdentifyMetalPerimeter(MakeInput(mixed, R)),
                      ContainsSubstring("disagree on the sign") &&
                          ContainsSubstring("spatial cluster"));
  }
}

namespace
{

const IdentifiedFeature *ClusterContaining(const IdentificationResult &result,
                                           const IdentificationInput &input,
                                           const std::array<double, 2> &point)
{
  // The spatial cluster one of whose claimed portions passes within 1e-6 of the point.
  for (const auto &feature : result.features)
  {
    if (feature.type != "SpatialEdgeCluster")
    {
      continue;
    }
    for (const auto &portion : feature.portions)
    {
      // Portions run from the segment's canonical key[0] (the lexicographically smaller
      // end), not from the input's p0.
      const auto &key = result.segments[portion.segment].key;
      const double length = std::hypot(key[1][0] - key[0][0], key[1][1] - key[0][1]);
      const double dx = (key[1][0] - key[0][0]) / length,
                   dy = (key[1][1] - key[0][1]) / length;
      const double s = std::clamp((point[0] - key[0][0]) * dx + (point[1] - key[0][1]) * dy,
                                  portion.s0, portion.s1);
      if (std::hypot(key[0][0] + s * dx - point[0], key[0][1] + s * dy - point[1]) <=
          1.0e-6)
      {
        return &feature;
      }
    }
  }
  return nullptr;
}

// A signature-frame point (units of R) in device coordinates.
std::array<double, 3> DevicePoint(const IdentifiedFeature &feature, double x, double y,
                                  double radius)
{
  std::array<double, 3> p = feature.origin;
  for (int d = 0; d < 3; d++)
  {
    p[d] += radius * (x * feature.axes[0][d] + y * feature.axes[1][d]);
  }
  return p;
}

// Device-coordinate bounding box of the feature's recorded support box.
std::array<double, 4> DeviceBox(const IdentifiedFeature &feature, const char *key,
                                double radius)
{
  const auto &box = feature.spatial_support.at(key);
  std::array<double, 4> device = {1.0e300, 1.0e300, -1.0e300, -1.0e300};
  for (const double bx : {box[0].get<double>(), box[2].get<double>()})
  {
    for (const double by : {box[1].get<double>(), box[3].get<double>()})
    {
      const auto p = DevicePoint(feature, bx, by, radius);
      device[0] = std::min(device[0], p[0]);
      device[1] = std::min(device[1], p[1]);
      device[2] = std::max(device[2], p[0]);
      device[3] = std::max(device[3], p[1]);
    }
  }
  return device;
}

// Distance of a device point from the metal perimeter of the input (its segments).
double PerimeterDistance(const IdentificationInput &input, const std::array<double, 3> &p)
{
  double best = 1.0e300;
  for (const auto &segment : input.segments)
  {
    const double dx = segment.p1[0] - segment.p0[0], dy = segment.p1[1] - segment.p0[1];
    const double length2 = dx * dx + dy * dy;
    double t = ((p[0] - segment.p0[0]) * dx + (p[1] - segment.p0[1]) * dy) / length2;
    t = std::clamp(t, 0.0, 1.0);
    best = std::min(
        best, std::hypot(p[0] - (segment.p0[0] + t * dx), p[1] - (segment.p0[1] + t * dy)));
  }
  return best;
}

}  // namespace

TEST_CASE("SurfaceResponseIdentificationSpatialSupportContract",
          "[surfaceresponseidentification][Serial]")
{
  // Spatial-support contract v3 (USER decision 281, supervisor decision 282): the coupon
  // metal of a cluster is the device plan clipped to the claims-derived box. Pad A
  // (conductor 0, top edge y = 0, right corner at corner_x) and a 2-um lead B (conductor 1)
  // ending 1 um above it form one cluster whose pad-edge claim is cut at x = +-6.873 and
  // whose box (claims + 3R past the cuts, 2R laterally) spans x in [-12.873, 12.873],
  // y in [-5, 12] in device coordinates.
  const double R = 2.0;
  auto Scene =
      [&](double corner_x, std::vector<LoopSpec> extra = {}, double subdivision = 100.0)
  {
    std::vector<LoopSpec> loops = {{Rectangle(-40.0, -10.0, corner_x, 0.0), 0, subdivision},
                                   {Rectangle(-1.0, 1.0, 1.0, 40.0), 1, subdivision}};
    loops.insert(loops.end(), extra.begin(), extra.end());
    return loops;
  };
  auto Identify =
      [&](const std::vector<LoopSpec> &loops, double mirror = 1.0, double rotate = 0.0)
  {
    std::vector<LoopSpec> transformed = loops;
    for (auto &loop : transformed)
    {
      for (auto &p : loop.points)
      {
        p = {mirror * p[0], p[1]};
      }
      if (mirror < 0.0)
      {
        std::reverse(loop.points.begin(), loop.points.end());
      }
      loop.points = RotatedAndShifted(
          loop.points, rotate, {rotate != 0.0 ? 100.0 : 0.0, rotate != 0.0 ? -50.0 : 0.0});
    }
    const auto input = MakeInput(transformed, R);
    auto result = IdentifyMetalPerimeter(input);
    CheckPartition(input, result);
    return std::make_pair(input, std::move(result));
  };
  const std::array<double, 2> probe = {0.0, 1.0};  // the lead's end edge: in the cluster

  // 1. Empty context (the device IS the straight continuation): the claims-only key byte
  //    for byte, no Box / Context in the signature, the record says Contract 2.
  std::string legacy_key, legacy_hash;
  {
    const auto [input, result] = Identify(Scene(40.0));
    const auto *cluster = ClusterContaining(result, input, probe);
    REQUIRE(cluster != nullptr);
    CHECK(cluster->signature["EdgeCount"] == 4);
    CHECK(!cluster->signature.contains("Box"));
    CHECK(!cluster->signature.contains("Context"));
    CHECK(!cluster->signature.contains("Unboxable"));
    const auto &support = cluster->spatial_support;
    REQUIRE(!support.is_null());
    CHECK(support["Contract"] == 2);
    CHECK(support["LegacyEquivalent"] == true);
    CHECK(support["Growth"]["Grown"] == false);
    CHECK(support["Context"]["ForeignPieces"] == 0);
    CHECK(support["Context"]["ChainPieces"] == 4);  // the four straight continuations
    CHECK_THAT(
        support["LegacyContinuation"]["FictitiousContinuationLengthOverR"].get<double>(),
        WithinAbs(0.0, 1.0e-9));
    CHECK_THAT(
        support["LegacyContinuation"]["StraightContinuationLengthOverR"].get<double>(),
        WithinAbs(12.0, 1.0e-5));
    const auto box = DeviceBox(*cluster, "Box", R);
    CHECK_THAT(box[0], WithinAbs(-12.872984, 1.0e-5));
    CHECK_THAT(box[2], WithinAbs(12.872984, 1.0e-5));
    CHECK_THAT(box[1], WithinAbs(-5.0, 1.0e-5));
    CHECK_THAT(box[3], WithinAbs(12.0, 1.0e-5));
    CHECK(support["Box"] == support["ClaimsBox"]);
    // The key is the claims-only canonicalisation (the unchanged v2 function) of the
    // feature's own portions and vertices: byte-identical to today's.
    std::vector<SignaturePortion> portions;
    for (const auto &portion : cluster->portions)
    {
      const auto &segment = input.segments[portion.segment];
      const auto &key = result.segments[portion.segment].key;
      const double length = std::hypot(key[1][0] - key[0][0], key[1][1] - key[0][1]);
      auto At = [&](double s)
      {
        std::array<double, 3> p;
        for (int d = 0; d < 3; d++)
        {
          p[d] = key[0][d] + s / length * (key[1][d] - key[0][d]);
        }
        return p;
      };
      portions.push_back({At(portion.s0),
                          At(portion.s1),
                          segment.gap_direction,
                          segment.conductor,
                          {"MS"},
                          segment.boundary_law});
    }
    std::vector<SignatureVertex> vertices;
    for (const auto &vertex : result.vertices)
    {
      if (vertex.feature == cluster->id && vertex.type != "BendVertex")
      {
        vertices.push_back(
            {input.vertices[vertex.vertex].coordinate, vertex.type, vertex.turn_degrees});
      }
    }
    REQUIRE(vertices.size() == 2);
    const auto canonical =
        CanonicalClusterSignature(portions, vertices, {0.0, 0.0, 1.0}, R);
    nlohmann::json expected = canonical.signature;
    expected["EdgeCount"] = portions.size();
    CHECK(SignatureKeyAndHash(expected, "SpatialEdgeCluster").first ==
          cluster->signature_key);
    CHECK(canonical.chirality == cluster->chirality);
    legacy_key = cluster->signature_key;
    legacy_hash = cluster->hash;
    // The box rule on the serialised signature is the one the record used.
    const auto from_signature = SupportBoxFromSignature(cluster->signature);
    for (int k = 0; k < 4; k++)
    {
      CHECK_THAT(from_signature[k], WithinAbs(support["Box"][k].get<double>(), 1.0e-12));
    }
    // Broadcast form carries the record.
    const auto copy =
        DeserializeIdentificationResult(SerializeIdentificationResult(result));
    CHECK(copy.ToJson(1.0) == result.ToJson(1.0));
    CHECK(copy.spatial_support.claims_keyed == result.spatial_support.claims_keyed);
    CHECK(result.spatial_support.claims_keyed == 2);  // this cluster + the lead's far end
    CHECK(result.spatial_support.context_keyed == 0);
  }

  // 2. The S1p pattern: the pad corner 1.56 R past the claim cut, inside the box. The
  //    straight continuation would run 1.436 R past the corner over the device's gap
  //    (fictitious metal, the D2 / D3-C mechanism); the v3 geometry follows the device
  //    chain through the corner, every context piece lies on the device perimeter and the
  //    key changes.
  std::string corner_hash;
  int corner_chirality = 0;
  {
    const auto [input, result] = Identify(Scene(10.0));
    const auto *cluster = ClusterContaining(result, input, probe);
    REQUIRE(cluster != nullptr);
    REQUIRE(cluster->signature.contains("Box"));
    REQUIRE(cluster->signature.contains("Context"));
    CHECK(cluster->signature["EdgeCount"] == 4);  // the claims only
    CHECK(cluster->signature_key != legacy_key);
    const auto &support = cluster->spatial_support;
    CHECK(support["Contract"] == 3);
    CHECK(support["LegacyEquivalent"] == false);
    CHECK(support["Growth"]["Grown"] == false);
    CHECK(support["Context"]["ForeignPieces"] == 0);
    CHECK(support["Context"]["ChainPieces"] == 5);
    REQUIRE(support["Context"]["ChainVertices"].size() == 1);
    const auto &chain_vertex = support["Context"]["ChainVertices"][0];
    CHECK(chain_vertex["Type"] == "ConvexCorner");
    // The corner's distance from the nearest face: 2.873 um = 1.436 R from the right face.
    CHECK_THAT(chain_vertex["FaceDistanceOverR"].get<double>(),
               WithinAbs((12.872984 - 10.0) / R, 1.0e-5));
    // The corner is a vertex feature of its own (no cluster): its feature id is recorded
    // (review MINOR-1) and names a ConvexCorner feature at the corner.
    CHECK(chain_vertex["Cluster"].is_null());
    REQUIRE(chain_vertex["Feature"].is_number_integer());
    {
      const auto &corner_feature =
          result.features[chain_vertex["Feature"].get<std::size_t>()];
      CHECK(corner_feature.type == "ConvexCorner");
      CHECK_THAT(corner_feature.origin[0], WithinAbs(10.0, 1.0e-5));
      CHECK_THAT(corner_feature.origin[1], WithinAbs(0.0, 1.0e-5));
    }
    const auto corner =
        DevicePoint(*cluster, chain_vertex["P"][0], chain_vertex["P"][1], R);
    CHECK_THAT(corner[0], WithinAbs(10.0, 1.0e-5));
    CHECK_THAT(corner[1], WithinAbs(0.0, 1.0e-5));
    CHECK_THAT(
        support["LegacyContinuation"]["FictitiousContinuationLengthOverR"].get<double>(),
        WithinAbs((12.872984 - 10.0) / R, 1.0e-5));
    // Every context piece lies on the device perimeter (no fictitious metal boundary) and
    // the pad's right edge below the corner is among them, flagged Chain.
    bool right_edge = false;
    double context_length = 0.0;
    for (const auto &entry : cluster->signature["Context"])
    {
      REQUIRE(entry.contains("Chain"));
      CHECK(entry["Chain"] == true);
      const auto a = DevicePoint(*cluster, entry["P"][0], entry["P"][1], R);
      const auto b = DevicePoint(*cluster, entry["P"][2], entry["P"][3], R);
      for (const auto &p : {a, b})
      {
        CHECK(PerimeterDistance(input, p) <= 1.0e-5);
      }
      context_length += std::hypot(b[0] - a[0], b[1] - a[1]);
      if (std::abs(a[0] - 10.0) < 1.0e-5 && std::abs(b[0] - 10.0) < 1.0e-5)
      {
        right_edge = true;
        CHECK_THAT(std::min(a[1], b[1]), WithinAbs(-5.0, 1.0e-5));
        CHECK_THAT(std::max(a[1], b[1]), WithinAbs(0.0, 1.0e-5));
        CHECK(entry["Conductor"] == cluster->signature["Portions"][3]["Conductor"]);
      }
    }
    CHECK(right_edge);
    CHECK_THAT(context_length / R,
               WithinAbs(support["Context"]["ChainLengthOverR"].get<double>(), 1.0e-5));
    // The legacy straight continuation (12 R) minus the fictitious part plus the right edge
    // (2.5 R): the chain length.
    CHECK_THAT(support["Context"]["ChainLengthOverR"].get<double>(),
               WithinAbs(12.0 - (12.872984 - 10.0) / R + 2.5, 1.0e-4));
    corner_hash = cluster->hash;
    corner_chirality = cluster->chirality;
    CHECK(corner_chirality != 0);  // the context breaks the mirror symmetry of the claims
    CHECK(result.spatial_support.context_keyed == 1);
    const auto copy =
        DeserializeIdentificationResult(SerializeIdentificationResult(result));
    CHECK(copy.ToJson(1.0) == result.ToJson(1.0));
  }

  // 3. Invariance: a mirror image has the same key with the opposite chirality; a rotated
  //    and shifted copy and a finer mesh the same key and chirality (the context comes
  //    from the runs, not from the mesh segments).
  {
    const auto [mirror_input, mirror] = Identify(Scene(10.0), -1.0);
    const auto *m = ClusterContaining(mirror, mirror_input, {0.0, 1.0});
    REQUIRE(m != nullptr);
    CHECK(m->hash == corner_hash);
    CHECK(m->chirality == -corner_chirality);
    const auto [rotated_input, rotated] = Identify(Scene(10.0), 1.0, 37.0);
    const auto rotated_probe = RotatedAndShifted({probe}, 37.0, {100.0, -50.0})[0];
    const auto *r = ClusterContaining(rotated, rotated_input, rotated_probe);
    REQUIRE(r != nullptr);
    CHECK(r->hash == corner_hash);
    CHECK(r->chirality == corner_chirality);
    const auto [fine_input, fine] = Identify(Scene(10.0, {}, 0.7));
    const auto *f = ClusterContaining(fine, fine_input, probe);
    REQUIRE(f != nullptr);
    CHECK(f->hash == corner_hash);
    CHECK(f->chirality == corner_chirality);
    CHECK(f->signature == ClusterContaining(Identify(Scene(10.0)).second,
                                            Identify(Scene(10.0)).first, probe)
                              ->signature);
  }

  // 4. Foreign metal in the box: a third conductor's narrow lead (0.1 R wide; below the
  //    joint noise resolution, so its end has no corner features and it reads as one
  //    chain) entering through the right face carries two crossings 0.1 R apart around
  //    METAL (allowed, recorded) and its three edges inside are foreign context: the key
  //    differs from the empty-context key, the context lists the foreign pieces with the
  //    next conductor label, the mirror image keeps the key.
  {
    const LoopSpec thin = {Rectangle(6.0, 8.0, 30.0, 8.2), 2, 100.0};
    const auto [input, result] = Identify(Scene(40.0, {thin}));
    const auto *cluster = ClusterContaining(result, input, probe);
    REQUIRE(cluster != nullptr);
    CHECK(cluster->hash != legacy_hash);
    CHECK(cluster->hash != corner_hash);
    const auto &support = cluster->spatial_support;
    CHECK(support["Contract"] == 3);
    CHECK(support["Growth"]["Grown"] == false);
    CHECK(support["Context"]["ForeignPieces"] == 3);
    CHECK(support["Context"]["ForeignConductors"] == 1);
    CHECK(support["Context"]["ChainPieces"] == 4);
    CHECK(support["Context"]["ForeignVertices"].empty());
    CHECK(support["FaceRules"]["NarrowCrossSections"] == 1);
    CHECK_THAT(support["FaceRules"]["MinCrossSectionOverR"].get<double>(),
               WithinAbs(0.1, 1.0e-6));
    CHECK_THAT(support["Context"]["ForeignLengthOverR"].get<double>(),
               WithinAbs((2.0 * (12.872984 - 6.0) + 0.2) / R, 1.0e-4));
    std::set<int> labels;
    for (const auto &entry : cluster->signature["Context"])
    {
      if (entry["Chain"] == false)
      {
        labels.insert(entry["Conductor"].get<int>());
      }
    }
    CHECK(labels == std::set<int>{3});
    const auto [mirror_input, mirror] = Identify(Scene(40.0, {thin}), -1.0);
    const auto *m = ClusterContaining(mirror, mirror_input, probe);
    REQUIRE(m != nullptr);
    CHECK(m->hash == cluster->hash);
    CHECK(m->chirality == -cluster->chirality);
  }

  // 5. T2 / T3: the pad corner 0.1 R inside the right face fails the vertex clearance; the
  //    face grows by one 0.25 R step and the context then holds the corner and the right
  //    edge; the box differs from the claims box on that face only.
  {
    const auto [input, result] = Identify(Scene(12.872984 - 0.2));
    const auto *cluster = ClusterContaining(result, input, probe);
    REQUIRE(cluster != nullptr);
    const auto &support = cluster->spatial_support;
    CHECK(support["Contract"] == 3);
    CHECK(support["Growth"]["Grown"] == true);
    int grown_faces = 0, total_steps = 0;
    for (const auto &steps : support["Growth"]["Steps"])
    {
      grown_faces += steps.get<int>() > 0 ? 1 : 0;
      total_steps += steps.get<int>();
    }
    CHECK(grown_faces == 1);
    CHECK(total_steps == 1);
    const auto claims_box = DeviceBox(*cluster, "ClaimsBox", R);
    const auto box = DeviceBox(*cluster, "Box", R);
    CHECK_THAT(box[2] - claims_box[2], WithinAbs(0.25 * R, 1.0e-5));
    CHECK_THAT(box[0], WithinAbs(claims_box[0], 1.0e-9));
    CHECK_THAT(box[1], WithinAbs(claims_box[1], 1.0e-9));
    CHECK_THAT(box[3], WithinAbs(claims_box[3], 1.0e-9));
    CHECK(support["Context"]["ChainVertices"].size() == 1);
    CHECK(support["FaceRules"]["MinClearanceOverR"].get<double>() >= 0.25);
    CHECK(result.spatial_support.grown == 1);
  }

  // 6. T1: the pad corner ON the claims-box face (within the snap) — the pad's right edge
  //    runs along the face (sin theta = 0): the face grows one step and the edge then sits
  //    exactly at the clearance (a threshold-band reading, reported).
  {
    const auto [input, result] = Identify(Scene(12.872984));
    const auto *cluster = ClusterContaining(result, input, probe);
    REQUIRE(cluster != nullptr);
    const auto &support = cluster->spatial_support;
    CHECK(support["Contract"] == 3);
    CHECK(support["Growth"]["Grown"] == true);
    int total_steps = 0;
    for (const auto &steps : support["Growth"]["Steps"])
    {
      total_steps += steps.get<int>();
    }
    CHECK(total_steps == 1);
    CHECK(support["FaceRules"]["ThresholdBandHits"].get<int>() >= 1);
    CHECK(support["FaceRules"]["MinClearanceOverR"].get<double>() >= 0.25 - 1.0e-9);
    CHECK(support["FaceRules"]["MinCrossingSine"].get<double>() >= 0.25);
    CHECK(support["Context"]["ChainVertices"].size() == 1);
  }

  // 7. A grazing face crossing: a foreign rectangle (conductor 2; its own 4-corner
  //    cluster, > 2R from everything else) tilted 10.4 deg whose long edge crosses the left
  //    face at sin(theta) = 0.18 < 0.25; the face grows until the whole rectangle is inside
  //    (4 steps: the grazing crossing, then one corner after another within the clearance
  //    of the moved face), then every crossing is gone and the rectangle is foreign context
  //    with its four corners.
  {
    const std::vector<Point2> slanted = {
        {-12.6, 5.2}, {-13.3, 9.0}, {-14.2838, 8.8188}, {-13.5838, 5.0188}};
    const auto [input, result] = Identify(Scene(40.0, {{slanted, 2, 100.0}}));
    const auto *cluster = ClusterContaining(result, input, probe);
    REQUIRE(cluster != nullptr);
    const auto &support = cluster->spatial_support;
    CHECK(support["Contract"] == 3);
    CHECK(support["Growth"]["Grown"] == true);
    int grown_faces = 0, total_steps = 0;
    for (const auto &steps : support["Growth"]["Steps"])
    {
      grown_faces += steps.get<int>() > 0 ? 1 : 0;
      total_steps += steps.get<int>();
    }
    CHECK(grown_faces == 1);
    CHECK(total_steps == 4);
    const auto claims_box = DeviceBox(*cluster, "ClaimsBox", R);
    const auto box = DeviceBox(*cluster, "Box", R);
    CHECK_THAT(claims_box[0] - box[0], WithinAbs(4 * 0.25 * R, 1.0e-5));
    CHECK(support["Context"]["ForeignPieces"] == 4);
    CHECK(support["Context"]["ForeignVertices"].size() == 4);
    double perimeter = 0.0;
    for (std::size_t k = 0; k < slanted.size(); k++)
    {
      const Point2 &a = slanted[k], &b = slanted[(k + 1) % slanted.size()];
      perimeter += std::hypot(b[0] - a[0], b[1] - a[1]);
    }
    CHECK_THAT(support["Context"]["ForeignLengthOverR"].get<double>(),
               WithinAbs(perimeter / R, 1.0e-5));
    for (const auto &entry : cluster->signature["Context"])
    {
      if (entry["Chain"] == false)
      {
        for (const auto &p : {DevicePoint(*cluster, entry["P"][0], entry["P"][1], R),
                              DevicePoint(*cluster, entry["P"][2], entry["P"][3], R)})
        {
          CHECK(PerimeterDistance(input, p) <= 1.0e-5);
        }
      }
    }
    CHECK(support["FaceRules"]["MinCrossingSine"].get<double>() >= 0.25);
    CHECK(support["FaceRules"]["MinClearanceOverR"].get<double>() >= 0.25 - 1.0e-9);
    CHECK(support["Unboxable"].is_null());
  }

  // 8. Two foreign leads of different conductors entering through the right face around a
  //    0.1 R GAP: no growth resolves it (the channel moves with the face) — after the step
  //    cap the cluster is UNBOXABLE: a Missing placeholder key (the claims-only signature +
  //    Unboxable), the record names the face and the reason.
  {
    const std::vector<LoopSpec> channel = {{Rectangle(6.0, 8.0, 30.0, 8.6), 2, 100.0},
                                           {Rectangle(6.0, 8.8, 30.0, 9.4), 3, 100.0}};
    const auto [input, result] = Identify(Scene(40.0, channel));
    const auto *cluster = ClusterContaining(result, input, probe);
    REQUIRE(cluster != nullptr);
    CHECK(cluster->signature["Unboxable"] == true);
    CHECK(!cluster->signature.contains("Box"));
    CHECK(cluster->signature_key != legacy_key);
    const auto &support = cluster->spatial_support;
    CHECK(support["Contract"] == 0);
    REQUIRE(support["Unboxable"].is_string());
    CHECK_THAT(support["Unboxable"].get<std::string>(), ContainsSubstring("bound a gap"));
    int max_steps = 0;
    for (const auto &steps : support["Growth"]["Steps"])
    {
      max_steps = std::max(max_steps, steps.get<int>());
    }
    CHECK(max_steps == kSupportFaceGrowthMaxSteps);
    // The two leads' ends form clusters of their own whose boxes hold the same channel.
    std::size_t unboxable = 0;
    for (const auto &feature : result.features)
    {
      if (feature.type == "SpatialEdgeCluster" &&
          feature.signature.value("Unboxable", false))
      {
        unboxable++;
        CHECK(feature.spatial_support["Contract"] == 0);
        CHECK(!feature.matched_model);
      }
    }
    CHECK(unboxable >= 1);
    CHECK(result.spatial_support.unboxable == unboxable);
    const auto copy =
        DeserializeIdentificationResult(SerializeIdentificationResult(result));
    CHECK(copy.ToJson(1.0) == result.ToJson(1.0));
  }

  // 9. Two clusters whose boxes overlap each other's claims (the S1p 19-edge pattern;
  //    decision 285 (1) on the R1a review's MAJOR-1): pad A now ends at x = 11.5 (its
  //    corner inside B's box), a second pad D (conductor 3) starts at x = 13.5 and a lead C
  //    (conductor 2) ends above it at x = 16..18. C's cluster claims pad D's edge and,
  //    within 2R of it, pad A's corner (11.5, 0) with the edge next to it: inside B's box
  //    those claims are context of B hashed as foreign (Chain false; C's coupon owns them),
  //    B's chain stops where C's claim begins, the pieces are split there, and A's corner
  //    is listed as a foreign vertex of B naming C's feature, never as B's chain vertex.
  {
    const std::vector<LoopSpec> loops = {{Rectangle(-40.0, -10.0, 11.5, 0.0), 0, 100.0},
                                         {Rectangle(-1.0, 1.0, 1.0, 40.0), 1, 100.0},
                                         {Rectangle(13.5, -10.0, 60.0, 0.0), 3, 100.0},
                                         {Rectangle(16.0, 1.0, 18.0, 40.0), 2, 100.0}};
    const auto [input, result] = Identify(loops);
    const auto *cluster_b = ClusterContaining(result, input, probe);
    const auto *cluster_c = ClusterContaining(result, input, {17.0, 1.0});
    REQUIRE(cluster_b != nullptr);
    REQUIRE(cluster_c != nullptr);
    REQUIRE(cluster_b != cluster_c);
    const auto &support = cluster_b->spatial_support;
    CHECK(support["Contract"] == 3);
    CHECK(support["LegacyEquivalent"] == false);
    CHECK(cluster_b->signature_key != legacy_key);
    CHECK(support["Context"]["ClaimedByOtherFeature"]["Pieces"].get<int>() >= 1);
    CHECK(support["Context"]["ClaimedByOtherFeature"]["LengthOverR"].get<double>() > 0.0);
    for (const auto &entry : support["Context"]["ClaimedByOtherFeature"]["Entries"])
    {
      CHECK(entry["Feature"] == cluster_c->id);
    }
    // Where C's claims begin on pad A's top edge (y = 0): the left-most such x of C.
    double c_start = 1.0e300;
    for (const auto &portion : cluster_c->signature["Portions"])
    {
      for (const auto &p : {DevicePoint(*cluster_c, portion["P"][0], portion["P"][1], R),
                            DevicePoint(*cluster_c, portion["P"][2], portion["P"][3], R)})
      {
        if (std::abs(p[1]) <= 1.0e-5 && p[0] <= 11.5 + 1.0e-5)
        {
          c_start = std::min(c_start, p[0]);
        }
      }
    }
    REQUIRE(c_start < 11.5);
    CHECK(c_start > 1.0);
    std::size_t chain_on_pad_edge = 0;
    for (const auto &entry : cluster_b->signature["Context"])
    {
      const auto a = DevicePoint(*cluster_b, entry["P"][0], entry["P"][1], R);
      const auto b = DevicePoint(*cluster_b, entry["P"][2], entry["P"][3], R);
      const bool on_pad_edge = std::abs(a[1]) <= 1.0e-5 && std::abs(b[1]) <= 1.0e-5;
      if (entry["Chain"] == true)
      {
        // A chain piece never lies on C's claims: on pad A's edge it ends where C's claim
        // begins; otherwise it is a side of lead B.
        if (on_pad_edge)
        {
          chain_on_pad_edge++;
          CHECK(std::max(a[0], b[0]) <= c_start + 1.0e-5);
        }
        else
        {
          CHECK(std::max(std::abs(a[0]), std::abs(b[0])) <= 1.0 + 1.0e-5);
        }
      }
      else if (on_pad_edge && std::min(a[0], b[0]) > 1.0)
      {
        // A pad-edge context piece right of B's chain is C's claim, split at c_start.
        CHECK(std::min(a[0], b[0]) >= c_start - 1.0e-5);
      }
    }
    CHECK(chain_on_pad_edge == 2);  // the left continuation and the right one up to C
    CHECK(support["Context"]["ChainVertices"].empty());
    bool a_corner_listed = false;
    for (const auto &vertex : support["Context"]["ForeignVertices"])
    {
      const auto p = DevicePoint(*cluster_b, vertex["P"][0], vertex["P"][1], R);
      if (std::abs(p[0] - 11.5) <= 1.0e-5 && std::abs(p[1]) <= 1.0e-5)
      {
        a_corner_listed = true;
        CHECK(vertex["Cluster"].is_number_integer());
        CHECK(vertex["Feature"] == cluster_c->id);
      }
    }
    CHECK(a_corner_listed);
    const auto copy =
        DeserializeIdentificationResult(SerializeIdentificationResult(result));
    CHECK(copy.ToJson(1.0) == result.ToJson(1.0));
  }

  // 10. Two-sided T2 (decision 285 (2) on the R1a review's MAJOR-2): the pad corner 0.16 R
  //     OUTSIDE the right face of the claims box. Under the one-sided rule the corner was
  //     invisible (the pad's top edge crosses the face straight: a legacy-equivalent key);
  //     the corner field sits on the face trace all the same. The face fails on the
  //     exterior vertex, grows one step (the corner is then 0.09 R inside: the interior
  //     clearance rule) and one more: two steps, the corner a chain vertex, the key v3.
  {
    const auto [input, result] = Identify(Scene(12.872984 + 0.16 * R));
    const auto *cluster = ClusterContaining(result, input, probe);
    REQUIRE(cluster != nullptr);
    const auto &support = cluster->spatial_support;
    CHECK(support["Contract"] == 3);
    CHECK(support["LegacyEquivalent"] == false);
    CHECK(support["Growth"]["Grown"] == true);
    int grown_faces = 0, total_steps = 0;
    for (const auto &steps : support["Growth"]["Steps"])
    {
      grown_faces += steps.get<int>() > 0 ? 1 : 0;
      total_steps += steps.get<int>();
    }
    CHECK(grown_faces == 1);
    CHECK(total_steps == 2);
    CHECK(support["Growth"]["AttemptedSteps"] == support["Growth"]["Steps"]);
    REQUIRE(support["Growth"]["StepReasons"].size() == 2);
    CHECK_THAT(support["Growth"]["StepReasons"][0].get<std::string>(),
               ContainsSubstring("outside face"));
    const auto claims_box = DeviceBox(*cluster, "ClaimsBox", R);
    const auto box = DeviceBox(*cluster, "Box", R);
    CHECK_THAT(box[2] - claims_box[2], WithinAbs(2 * 0.25 * R, 1.0e-5));
    CHECK(support["Context"]["ChainVertices"].size() == 1);
    // The exterior vertex is read on the FIRST pass (accumulated over every pass, decision
    // 287 (a)); the final pass has no exterior reading left.
    // (the corner is the end of two shell pieces, the pad's top and right edges: 2
    // readings; the canonical frame maps the device's right face to face 3)
    CHECK(support["FaceRules"]["ExteriorVertices"] == 2);
    CHECK(support["FaceRules"]["ExteriorEdges"] == 0);
    REQUIRE(support["Growth"]["Passes"].size() == 3);
    CHECK(support["Growth"]["Passes"][0]["ExteriorVertices"] == 2);
    CHECK(support["Growth"]["Passes"][0]["Failing"][3] == true);
    CHECK(support["Growth"]["Passes"][0]["Failures"].size() >= 1);
    CHECK(support["Growth"]["Passes"][2]["ExteriorVertices"] == 0);
    for (const auto &failing : support["Growth"]["Passes"][2]["Failing"])
    {
      CHECK(failing == false);
    }
    CHECK(support["Growth"]["Passes"][2]["Failures"].empty());
    CHECK(support["FaceRules"]["Passes"] == 3);
    CHECK(support["FaceRules"]["MinClearanceOverR"].get<double>() >= 0.25 - 1.0e-9);
    // The same corner 0.5 R outside the face is beyond the clearance shell: invisible by
    // the rule, the legacy-equivalent key as before.
    const auto [input_far, result_far] = Identify(Scene(12.872984 + 0.5 * R));
    const auto *far = ClusterContaining(result_far, input_far, probe);
    REQUIRE(far != nullptr);
    CHECK(far->spatial_support["Contract"] == 2);
    CHECK(far->signature_key == legacy_key);
  }

  // 11. The knife edge (R1 final review MAJOR-1; decisions 287 (b) / 288): the pad corner
  //     ON the claims-box face is the S1p / C2p pattern — step 1 moves the face by exactly
  //     0.25 R and the pad's right edge then reads EXACTLY the clearance from the moved
  //     face, up to the box quantisation residual (half a 1e-6 R quantum) and noise.
  //     Sub-quantum perturbations (+-4e-7 R, the largest grid-preserving shift) of the
  //     layout (the corner) and of the frame origin (lead B with the pad-edge claims it
  //     cuts) give the SAME growth sequence and the SAME key: the quantised comparison
  //     reads the edge AT the clearance, which passes (the rule's own ">="). Before the
  //     fix `distance + 1e-9 < clearance` gave 2 steps on the minus side, 1 on the plus
  //     side, and two keys. A FULL-quantum move (+-1e-6 R) is a real geometric change: the
  //     hashed context coordinate moves by one quantum (the key changes) and the corner
  //     moved OUTWARD legitimately grows twice (the edge is then one quantum closer than
  //     the clearance to the moved face) — keys are geometry-precise to the quantum,
  //     sub-quantum noise cannot move them. (The canonical frame maps the device's right
  //     face to face 3.)
  {
    const auto [input_far, result_far] = Identify(Scene(40.0));
    const auto *far = ClusterContaining(result_far, input_far, probe);
    REQUIRE(far != nullptr);
    const double face_x = DeviceBox(*far, "ClaimsBox", R)[2];  // the right face
    CHECK_THAT(face_x, WithinAbs(12.872984, 1.0e-5));
    struct Reading
    {
      std::string key;
      std::array<int, 4> steps;
      int band_hits;
      nlohmann::json box;
    };
    auto Read = [&](double corner_x, double lead_shift)
    {
      const std::vector<LoopSpec> loops = {
          {Rectangle(-40.0, -10.0, corner_x, 0.0), 0, 100.0},
          {Rectangle(-1.0 + lead_shift, 1.0, 1.0 + lead_shift, 40.0), 1, 100.0}};
      const auto [input, result] = Identify(loops);
      const auto *cluster = ClusterContaining(result, input, {lead_shift, 1.0});
      REQUIRE(cluster != nullptr);
      REQUIRE(cluster->spatial_support["Contract"] == 3);
      return Reading{cluster->signature_key,
                     cluster->spatial_support["Growth"]["Steps"].get<std::array<int, 4>>(),
                     cluster->spatial_support["FaceRules"]["ThresholdBandHits"].get<int>(),
                     cluster->spatial_support["Box"]};
    };
    const Reading base = Read(face_x, 0.0);
    CHECK(base.steps == std::array<int, 4>{0, 0, 0, 1});
    const double sub_quantum = 4.0e-7 * R;
    for (const double sign : {-1.0, 1.0})
    {
      // The layout: the corner moved relative to the cluster.
      const Reading layout = Read(face_x + sign * sub_quantum, 0.0);
      CHECK(layout.steps == base.steps);
      CHECK(layout.key == base.key);
      CHECK(layout.box == base.box);
      // The frame origin: the lead (its claims cut the pad edge, the canonical frame and
      // the claims box follow it) moved, the corner fixed.
      const Reading frame = Read(face_x, sign * sub_quantum);
      CHECK(frame.steps == base.steps);
      CHECK(frame.key == base.key);
      CHECK(frame.box == base.box);
    }
    // The threshold reading is reported on every pass (decision 287 (a)): the base case
    // reads the edge at the clearance on its second pass.
    CHECK(base.band_hits >= 1);
    // The complement: one full quantum is a real move.
    const double quantum = 1.0e-6 * R;
    const Reading inside = Read(face_x - quantum, 0.0);
    const Reading outside = Read(face_x + quantum, 0.0);
    CHECK(inside.steps == std::array<int, 4>{0, 0, 0, 1});
    CHECK(outside.steps == std::array<int, 4>{0, 0, 0, 2});
    CHECK(inside.key != base.key);
    CHECK(outside.key != base.key);
    CHECK(inside.key != outside.key);
  }

  // 12. Two faces failing on one pass (R1 final review MINOR-8): the pad corner 0.1 R
  // inside
  //     the right face (scene 5) and a foreign strip (conductor 2, 2.5 R above the pad)
  //     whose end sits 0.1 R inside the LEFT face: both faces grow on the first step and
  //     the pass record lists BOTH failures (StepReasons keeps the first only).
  {
    const double left_x = -12.872984;
    const std::vector<LoopSpec> strip = {
        {Rectangle(left_x - 3.0, 5.0, left_x + 0.1 * R, 6.0), 2, 100.0}};
    const auto [input, result] = Identify(Scene(12.872984 - 0.1 * R, strip));
    const auto *cluster = ClusterContaining(result, input, probe);
    REQUIRE(cluster != nullptr);
    const auto &support = cluster->spatial_support;
    CHECK(support["Growth"]["Steps"] == std::array<int, 4>{0, 1, 0, 1});  // frame faces
    REQUIRE(support["Growth"]["StepReasons"].size() == 1);
    REQUIRE(support["Growth"]["Passes"].size() == 2);
    const auto &first = support["Growth"]["Passes"][0];
    CHECK(first["Failing"] == std::array<bool, 4>{false, true, false, true});
    for (const int face : {1, 3})
    {
      bool listed = false;
      for (const auto &failure : first["Failures"])
      {
        listed = listed || failure.get<std::string>().find(
                               "face " + std::to_string(face)) != std::string::npos;
      }
      CHECK(listed);
    }
    CHECK(support["Growth"]["Passes"][1]["Failures"].empty());
  }

  // 13. Another cluster's claim CUT is no device vertex (R1 final review MINOR-1): pad A
  //     continues far past B's box (its corner outside the clearance shell), a lead C ends
  //     above it to the right so that C's claim on the pad edge begins INSIDE B's box
  //     within 0.25 R of B's right face. The cut splits the context (C's piece is
  //     ClaimedByOtherFeature, B's chain ends there) but the end it makes is exempt from
  //     the clearance rule: no growth. Before the fix the cut read as "a device vertex or
  //     face crossing" within the clearance and the face grew.
  {
    auto WithLeadC = [&](double lead_x)
    {
      std::vector<LoopSpec> loops = Scene(60.0);
      loops.push_back({Rectangle(lead_x - 1.0, 1.0, lead_x + 1.0, 40.0), 2, 100.0});
      return loops;
    };
    // C's claim window on the pad edge from a first placement: the left-most claimed x.
    auto ClaimStart = [&](double lead_x)
    {
      const auto [input, result] = Identify(WithLeadC(lead_x));
      const auto *cluster_c = ClusterContaining(result, input, {lead_x, 1.0});
      REQUIRE(cluster_c != nullptr);
      double start = 1.0e300;
      for (const auto &portion : cluster_c->signature["Portions"])
      {
        for (const auto &p : {DevicePoint(*cluster_c, portion["P"][0], portion["P"][1], R),
                              DevicePoint(*cluster_c, portion["P"][2], portion["P"][3], R)})
        {
          if (std::abs(p[1]) <= 1.0e-5)
          {
            start = std::min(start, p[0]);
          }
        }
      }
      REQUIRE(start < 1.0e299);
      return start;
    };
    const double probe_x = 30.0;
    const double window = probe_x - ClaimStart(probe_x);  // C's reach along the pad edge
    const double face_x = 12.872984;
    const double lead_x = face_x - 0.2 * R + window;  // the cut 0.2 R inside B's face
    const auto [input, result] = Identify(WithLeadC(lead_x));
    const auto *cluster_b = ClusterContaining(result, input, probe);
    const auto *cluster_c = ClusterContaining(result, input, {lead_x, 1.0});
    REQUIRE(cluster_b != nullptr);
    REQUIRE(cluster_c != nullptr);
    REQUIRE(cluster_b != cluster_c);
    const auto &support = cluster_b->spatial_support;
    CHECK(support["Contract"] == 3);
    CHECK(support["LegacyEquivalent"] == false);
    CHECK(support["Growth"]["Grown"] == false);
    CHECK(support["Growth"]["Steps"] == std::array<int, 4>{0, 0, 0, 0});
    CHECK(support["Context"]["ClaimedByOtherFeature"]["Pieces"].get<int>() >= 1);
    // The cut lies where C's claim begins: inside the box, within the clearance of face 2.
    double cut_x = 1.0e300;
    for (const auto &entry : support["Context"]["ClaimedByOtherFeature"]["Entries"])
    {
      CHECK(entry["Feature"] == cluster_c->id);
      for (const auto &p : {DevicePoint(*cluster_b, entry["P"][0], entry["P"][1], R),
                            DevicePoint(*cluster_b, entry["P"][2], entry["P"][3], R)})
      {
        cut_x = std::min(cut_x, p[0]);
      }
    }
    CHECK_THAT(cut_x, WithinAbs(face_x - 0.2 * R, 1.0e-3));
    CHECK(support["Box"] == support["ClaimsBox"]);
  }
}

TEST_CASE("SurfaceResponseIdentificationLegacyContractAlias",
          "[surfaceresponseidentification][surfaceresponseoperator][Serial]")
{
  // USER decision 283: a legacy model (built under the decision-236 straight-continuation
  // contract) serves a contract-3 key ONLY through an explicit library alias naming that
  // key and its context digest; the alias resolves to a recorded LegacyContract match and
  // fails closed on a digest or claims mismatch. The scene of
  // SurfaceResponseIdentificationSpatialSupportContract: the legacy model's Signature is
  // the claims-only signature of the corner-far cluster (Contract 2), the aliased key the
  // corner-at-10 cluster's (Contract 3, the same claims, a foreign-free chain context).
  const double R = 2.0;
  auto Scene = [&](double corner_x)
  {
    return std::vector<LoopSpec>{{Rectangle(-40.0, -10.0, corner_x, 0.0), 0, 100.0},
                                 {Rectangle(-1.0, 1.0, 1.0, 40.0), 1, 100.0}};
  };
  const auto legacy_input = MakeInput(Scene(40.0), R);
  const auto legacy = IdentifyMetalPerimeter(legacy_input);
  const auto v3_input = MakeInput(Scene(10.0), R);
  const auto v3 = IdentifyMetalPerimeter(v3_input);
  const auto *legacy_cluster = ClusterContaining(legacy, legacy_input, {0.0, 1.0});
  const auto *v3_cluster = ClusterContaining(v3, v3_input, {0.0, 1.0});
  REQUIRE(legacy_cluster != nullptr);
  REQUIRE(v3_cluster != nullptr);
  REQUIRE(legacy_cluster->spatial_support["Contract"] == 2);
  REQUIRE(v3_cluster->spatial_support["Contract"] == 3);
  const auto *other_cluster = ClusterContaining(v3, v3_input, {0.0, 40.0});  // the lead end
  REQUIRE(other_cluster != nullptr);
  REQUIRE(other_cluster->hash != legacy_cluster->hash);

  // The context digest: of the Box + Context only (the claims-only signature has none), and
  // recorded by the identification for every contract-3 feature.
  const std::string digest = SpatialSupportContextDigest(v3_cluster->signature);
  CHECK(digest.size() == 64);
  CHECK(SpatialSupportContextDigest(legacy_cluster->signature).empty());
  CHECK(v3_cluster->spatial_support["ContextDigest"] == digest);
  CHECK(legacy_cluster->spatial_support["ContextDigest"].is_null());
  // The recorded claims-only key of the v3 cluster IS the legacy key (the v3 Portions are
  // serialised in the Box + Context frame, so the key is recorded, not re-derived).
  CHECK(v3_cluster->spatial_support["ClaimsKey"] == legacy_cluster->hash);
  CHECK(legacy_cluster->spatial_support["ClaimsKey"] == legacy_cluster->hash);
  // Neither feature carries an alias record out of the identification itself.
  CHECK(v3_cluster->legacy_contract.is_null());
  CHECK(legacy_cluster->legacy_contract.is_null());
  // The claims-only frame of the v3 cluster is the legacy cluster's frame (the same claims
  // canonicalised alone): the frame a legacy model is placed in. The v3 frame differs here
  // (the context breaks the mirror symmetry: chirality -+1 vs 0).
  CHECK(v3_cluster->claims_chirality == legacy_cluster->chirality);
  CHECK(v3_cluster->claims_chirality != v3_cluster->chirality);
  for (int d = 0; d < 3; d++)
  {
    CHECK_THAT(v3_cluster->claims_origin[d], WithinAbs(legacy_cluster->origin[d], 1.0e-9));
    for (int k = 0; k < 3; k++)
    {
      CHECK_THAT(v3_cluster->claims_axes[k][d],
                 WithinAbs(legacy_cluster->axes[k][d], 1.0e-12));
    }
  }
  for (int d = 0; d < 3; d++)
  {
    CHECK_THAT(legacy_cluster->claims_origin[d], WithinAbs(legacy_cluster->origin[d], 0.0));
  }

  LegacyContractAlias alias;
  alias.key = v3_cluster->hash;
  alias.context_digest = digest;
  alias.reason = "unit test: USER decision 283";
  const nlohmann::json record = ResolveLegacyContractAlias(
      "legacy-model", legacy_cluster->signature, alias, *v3_cluster);
  CHECK(record["Model"] == "legacy-model");
  CHECK(record["Key"] == v3_cluster->hash);
  CHECK(record["ContextDigest"] == digest);
  CHECK(record["ClaimsKey"] == legacy_cluster->hash);
  CHECK(record["Reason"] == alias.reason);
  CHECK(record["Context"]["Box"] == v3_cluster->signature["Box"]);
  CHECK(record["Context"]["Context"] == v3_cluster->signature["Context"]);

  // Fail closed: a context digest that is not the feature's.
  LegacyContractAlias wrong_digest = alias;
  wrong_digest.context_digest = std::string(64, '0');
  CHECK_THROWS_WITH(ResolveLegacyContractAlias("legacy-model", legacy_cluster->signature,
                                               wrong_digest, *v3_cluster),
                    ContainsSubstring("does not match the alias's ContextDigest"));
  // Fail closed: the alias on a model whose claims are not the feature's.
  CHECK_THROWS_WITH(ResolveLegacyContractAlias("other-model", other_cluster->signature,
                                               alias, *v3_cluster),
                    ContainsSubstring("is not the legacy model's key"));
  // Fail closed: an alias resolved for a feature other than its key (never a fallback).
  CHECK_THROWS_WITH(ResolveLegacyContractAlias("legacy-model", legacy_cluster->signature,
                                               alias, *legacy_cluster),
                    ContainsSubstring("is not the alias key"));
  // A contract-2 feature carries no context: an alias naming it cannot verify a digest.
  LegacyContractAlias claims_alias = alias;
  claims_alias.key = legacy_cluster->hash;
  CHECK_THROWS_WITH(ResolveLegacyContractAlias("legacy-model", legacy_cluster->signature,
                                               claims_alias, *legacy_cluster),
                    ContainsSubstring("(no Box)"));

  // The record survives the broadcast form and reaches the manifest's Match entry.
  IdentificationResult copy_source = v3;
  for (auto &feature : copy_source.features)
  {
    if (feature.id == v3_cluster->id)
    {
      feature.legacy_contract = record;
      feature.matched_model = "legacy-model";
      feature.match_deviation = 0.0;
    }
  }
  const auto copy =
      DeserializeIdentificationResult(SerializeIdentificationResult(copy_source));
  CHECK(copy.ToJson(1.0) == copy_source.ToJson(1.0));
  bool found = false;
  const nlohmann::json manifest = copy.ToJson(1.0);
  for (const auto &entry : manifest["Features"])
  {
    if (entry["Id"].get<int>() == v3_cluster->id)
    {
      found = true;
      CHECK(entry["Match"]["Status"] == "Matched");
      CHECK(entry["Match"]["Model"] == "legacy-model");
      CHECK(entry["Match"]["LegacyContract"]["Key"] == v3_cluster->hash);
      CHECK(entry["ClaimsFrame"]["Chirality"] == legacy_cluster->chirality);
    }
    else
    {
      CHECK(!entry["Match"].contains("LegacyContract"));
    }
  }
  CHECK(found);
}
