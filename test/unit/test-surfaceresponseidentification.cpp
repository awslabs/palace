// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <optional>
#include <set>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "models/surfaceresponseidentification.hpp"
#include "utils/metaledge.hpp"

using namespace palace;
using namespace Catch::Matchers;

namespace
{

using Point2 = std::array<double, 2>;

// Plan-view metal loops (counter-clockwise: metal inside) subdivided into mesh segments of
// at most the given length, with the perimeter frames the classifier would supply (gap
// direction = outward normal, process normal = +z, one PEC conductor per loop or a shared
// one).
struct LoopSpec
{
  std::vector<Point2> points;
  int conductor = 0;
  double subdivision = 1.0;
  double z = 0.0;  // metal plane of the loop
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
  };
  std::vector<Raw> raw;
  int chain = 0;
  for (const auto &loop : loops)
  {
    const std::size_t n = loop.points.size();
    // Chains break at corners (turn > 30 deg): every polygon edge here is its own chain.
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
                       loop.z});
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
    segment.chain = r.chain;
    segment.conductor = r.conductor;
    segment.targets = {{InterfaceDielectric::MS, 1}};
    segment.gap_direction = {r.outward[0], r.outward[1], 0.0};
    segment.process_normal = {0.0, 0.0, 1.0};
    segment.boundary_law = "{\"Type\":\"PEC\"}";
    input.vertices[segment.vertices[0]].segments.push_back(input.segments.size());
    input.vertices[segment.vertices[1]].segments.push_back(input.segments.size());
    input.segments.push_back(segment);
  }
  // Vertex types as metaledge.cpp: two segments -> corner iff turn > 30 deg, else regular.
  for (auto &vertex : input.vertices)
  {
    if (vertex.segments.size() == 1)
    {
      vertex.physical_type = MetalEdgeVertexType::ENDPOINT;
    }
    else if (vertex.segments.size() > 2)
    {
      vertex.physical_type = MetalEdgeVertexType::JUNCTION;
    }
    else
    {
      std::array<std::array<double, 3>, 2> directions;
      for (int i = 0; i < 2; i++)
      {
        const auto &segment = input.segments[vertex.segments[i]];
        const auto &other =
            input.vertices[segment.vertices[0] == (&vertex - &input.vertices[0])
                               ? segment.vertices[1]
                               : segment.vertices[0]];
        double norm = 0.0;
        for (int d = 0; d < 3; d++)
        {
          directions[i][d] = other.coordinate[d] - vertex.coordinate[d];
          norm += directions[i][d] * directions[i][d];
        }
        for (double &value : directions[i])
        {
          value /= std::sqrt(norm);
        }
      }
      double dot = 0.0;
      for (int d = 0; d < 3; d++)
      {
        dot += directions[0][d] * directions[1][d];
      }
      vertex.physical_type = dot <= -std::cos(30.0 * std::acos(-1.0) / 180.0)
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
    if (chain_of[seed] >= 0)
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
          if (chain_of[other] < 0)
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
// the bar sideways.
std::vector<Point2> BarAroundCentreline(const std::vector<Point2> &centreline, double width,
                                        double offset = 0.0)
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
      if (i == 0 || i + 1 == n)
      {
        const Point2 d = i == 0 ? Point2{centreline[1][0] - p[0], centreline[1][1] - p[1]}
                                : Point2{p[0] - centreline[i - 1][0],
                                         p[1] - centreline[i - 1][1]};
        const Point2 nrm = Normal(d);
        result.push_back(
            {p[0] + (offset + sign * h) * nrm[0], p[1] + (offset + sign * h) * nrm[1]});
      }
      else
      {
        const Point2 n0 = Normal({p[0] - centreline[i - 1][0], p[1] - centreline[i - 1][1]});
        const Point2 n1 = Normal({centreline[i + 1][0] - p[0], centreline[i + 1][1] - p[1]});
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

std::vector<Point2> ArcBar(double width, double radius, double sweep_degrees,
                           double step_degrees, double lead = 6.0, double offset = 0.0)
{
  const int steps = std::max(1, static_cast<int>(std::lround(sweep_degrees / step_degrees)));
  const double step = sweep_degrees * std::acos(-1.0) / 180.0 / steps;
  std::vector<Point2> centreline;
  for (int k = 0; k <= steps; k++)
  {
    centreline.push_back({radius * std::sin(k * step), radius - radius * std::cos(k * step)});
  }
  const Point2 d_end = {std::cos(steps * step), std::sin(steps * step)};
  centreline.insert(centreline.begin(), {-lead, 0.0});
  centreline.push_back({centreline.back()[0] + lead * d_end[0],
                        centreline.back()[1] + lead * d_end[1]});
  return BarAroundCentreline(centreline, width, offset);
}

// A bar along an S-bend: a left arc of the given radius and sweep followed by a right arc
// of the same radius and sweep (the tangent returns to +x), straight leads at both ends.
std::vector<Point2> SBar(double width, double radius, double sweep_degrees,
                         double step_degrees, double lead = 6.0)
{
  const int steps = std::max(1, static_cast<int>(std::lround(sweep_degrees / step_degrees)));
  const double step = sweep_degrees * std::acos(-1.0) / 180.0 / steps;
  std::vector<Point2> centreline = {{-lead, 0.0}};
  for (int k = 0; k <= steps; k++)
  {
    centreline.push_back({radius * std::sin(k * step), radius - radius * std::cos(k * step)});
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
  }
  centreline.push_back({centreline.back()[0] + lead, centreline.back()[1]});
  return BarAroundCentreline(centreline, width);
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
      CHECK(vertex.feature >= 0);
    }
  }
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

TEST_CASE("SurfaceResponseIdentificationCurvedEdges",
          "[surfaceresponseidentification][Serial]")
{
  // Curved-edge chain rule (decision 73(1)): a 3 um bar (1.5 R) along a polyline arc. The
  // two sides are a pair along the bend (constant separation), never a cluster; the arc is
  // a curved pair when the inner bend radius is below 10 R (decision 75) and a straight-like strip with a
  // curvature annotation otherwise; the leads are a plain strip; the two bar ends are one
  // two-corner cluster each. The classes do not depend on the discretisation (turn per
  // vertex) and the features do not change under mesh refinement (A5).
  const double R = 2.0;
  struct Case
  {
    double radius, sweep, step;
    bool curved;
  };
  for (const Case &c : {Case{5.0, 90.0, 20.0, true}, Case{5.0, 90.0, 1.0, true},
                        Case{20.0, 90.0, 5.0, true}, Case{50.0, 45.0, 20.0, false},
                        Case{50.0, 45.0, 1.0, false}, Case{250.0, 15.0, 5.0, false}})
  {
    const auto bar = ArcBar(3.0, c.radius, c.sweep, c.step);
    const auto input = MakeInput({{bar, 0, 1.0}}, R);
    const auto result = IdentifyMetalPerimeter(input);
    INFO("radius " << c.radius << " step " << c.step);
    CheckPartition(input, result);
    std::map<std::string, int> counts;
    for (const auto &feature : result.features)
    {
      counts[feature.type]++;
      if (feature.type == "CurvedSameConductorStrip")
      {
        // The inner side's radius (radius - 1.5): the arc rule reads the polyline's own
        // inscribed circle, which for an offset polyline at 20 deg per vertex differs from
        // the design radius within the arc-fit tolerance (5 %).
        CHECK_THAT(feature.signature["RadiusOverR"].get<double>(),
                   WithinRel((c.radius - 1.5) / R, 0.05));
        CHECK_THAT(feature.signature["SeparationOverR"].get<double>(),
                   WithinRel(1.5, 0.02));
      }
      if (feature.type == "SameConductorStrip")
      {
        CHECK_THAT(feature.signature["SeparationOverR"].get<double>(),
                   WithinRel(1.5, 0.02));
        if (!c.curved)
        {
          REQUIRE(feature.bend_radius_over_R.has_value());
          CHECK_THAT(*feature.bend_radius_over_R, WithinRel((c.radius - 1.5) / R, 0.03));
        }
      }
    }
    CHECK(counts["SpatialEdgeCluster"] == 2);
    CHECK(counts["SameConductorStrip"] >= 1);
    CHECK(counts["CurvedSameConductorStrip"] == (c.curved ? 1 : 0));
    CHECK(counts["IsolatedEdge"] == 0);
    CHECK(counts["CurvedEdge"] == 0);
    CHECK(counts["ConvexCorner"] == 0);
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
    std::map<double, std::pair<std::string, double>> by_radius;  // RadiusOverR -> (convexity, turn)
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
  // (the edge the coupon's first edge e1 lands on): the inner side of a strip has its metal
  // outside the bend (Concave), the outer side inside (Convex).
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
      const bool first_inner =
          first_radius / first_count < other_radius / other_count;
      CHECK(convexity == (first_inner ? "Concave" : "Convex"));
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
  // A gap of exactly 2R (and 2R +/- 1e-3 R) between two 8 um bars (4 R: the far corners of a
  // bar end are beyond the 3R vertex-join reach of the near corners) that are offsets of one
  // centreline along bends of 50 and 250 um at three discretisations: the interaction
  // decision uses the separation of the underlying curves (rule at kPairSeparationTolerance):
  // exactly 2R and 2R + 1e-3 R are isolated edges at every discretisation (no cluster, no
  // pair — the mid-chord dips of the chords below 2R, DS-SCT-001's 4 um gaps at 3.9998 that
  // became 3 mm clusters, create no events); 2R - 1e-3 R interacts where both readings of the
  // separation are below 2R, i.e. where the local joint turn satisfies
  // gap / cos(turn / 2) < 2R (the recorded discretisation ambiguity: an inscribed polyline
  // pair at this chord separation could be 2R apart): the whole pair at 1 deg per vertex,
  // the straight leads only at 5 deg (their joint turn is half a step), nothing at 15 deg.
  // The corner pairs across a 2R - 1e-3 R gap are events (two clusters) in every case.
  const double R = 2.0, width = 8.0;
  for (const double radius : {50.0, 250.0})
  {
    const double sweep = radius < 100.0 ? 45.0 : 15.0;
    for (const double step : {1.0, 5.0, 15.0})
    {
      for (const double gap : {2.0 * R, 2.0 * R - 1.0e-3 * R, 2.0 * R + 1.0e-3 * R})
      {
        // Both bars are offsets of the gap's centreline (radius), so the facing edges are
        // exactly parallel chords at the design gap along the bend (an offset path).
        const auto inner = ArcBar(width, radius, sweep, step, 6.0, 0.5 * gap + 0.5 * width);
        const auto outer = ArcBar(width, radius, sweep, step, 6.0, -0.5 * gap - 0.5 * width);
        const auto input = MakeInput({{inner, 0, 1.0}, {outer, 1, 1.0}}, R);
        const auto result = IdentifyMetalPerimeter(input);
        INFO("radius " << radius << " step " << step << " gap " << gap);
        CheckPartition(input, result);
        std::map<std::string, int> counts;
        for (const auto &feature : result.features)
        {
          counts[feature.type]++;
          if (feature.type == "DifferentConductorGap")
          {
            // Decision 85(1): the separation along the bends is the radius difference of the
            // two fitted arcs (exact for polylines inscribed in the design circles); an
            // OFFSET polyline pair (parallel chords at the design gap, this construction)
            // reads gap / cos(turn / 2) at its vertices — the recorded discretisation
            // ambiguity, within the signature parameter tolerance (1e-3 R) at these steps.
            CHECK_THAT(feature.signature["SeparationOverR"].get<double>(),
                       WithinAbs(gap / R, kSignatureParameterToleranceOverRadius));
          }
        }
        const bool corner_events = gap < 2.0 * R;
        const double lead_turn = 0.5 * step * std::acos(-1.0) / 180.0;
        const bool pair = gap / std::cos(0.5 * lead_turn) < 2.0 * R;
        CHECK(counts["DifferentConductorGap"] == (pair ? 1 : 0));
        CHECK(counts["SpatialEdgeCluster"] == (corner_events ? 2 : 0));
        CHECK(counts["ConvexCorner"] == (corner_events ? 4 : 8));
        CHECK(counts["CurvedEdge"] == 0);
        CHECK(counts["CurvedDifferentConductorGap"] == 0);
        CHECK(counts["IsolatedEdge"] >= (pair ? 2 : 4));
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
  // along a bend), straight and curved; the pairwise candidates inside it are superseded (no
  // claim of one priority ever overlaps another: Diagnostics); a member taken by a cluster is
  // recomposed out of the cross-section at the stack ends; sides in the canonical order.
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
    // over a longer reach than the trace edges, exactly R from the cores: the trace strip is
    // recomposed there (the stack-end rule).
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
    // ground | 1 | trace 1.5 | 3 | ground: offsets 0 / 0.5 / 1.25 / 2.75 R; chirality +1 and
    // one side index per physical edge over the whole feature.
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
    // bend of the centreline with 12 um leads: a straight ParallelEdgeCluster on the leads and
    // a CurvedParallelEdgeCluster along the bend with RadiusOverR = the trace's inner radius.
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
      if (feature.type == "ParallelEdgeCluster" || feature.type == "CurvedParallelEdgeCluster")
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
}

TEST_CASE("SurfaceResponseIdentificationExactParametersAndExtension",
          "[surfaceresponseidentification][Serial]")
{
  // Decision 85 (2026-09-26). (1) Exact signature parameters: the curved 3-edge stack of the
  // previous test at two discretisations of the bend (5 and 2.5 deg steps) gives IDENTICAL
  // offsets (0 / 1.5 / 2.5 R exactly: the arc radius differences) and bend radius, hence
  // one signature key; the tolerance API groups near-identical instances and tells the
  // mirror orientation apart from a different topology. (2) Cluster extension: the two
  // 12 x 8 pads of the first test have every single-edge portion within 2R of a cluster's
  // claimed perimeter absorbed (the gap edges between the end clusters are the pair, the
  // remainder isolated only where nothing is within 2R across).
  const double R = 2.0;
  // A band between two concentric circles (vertices ON the circles: the inscribed
  // construction of a CAD polygonisation) with tangent leads, counter-clockwise.
  auto InscribedBand = [](double r_in, double r_out, double sweep_degrees, double step_degrees,
                          double lead)
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
      const auto result = IdentifyMetalPerimeter(MakeInput({{trace, 0, 1.0}, {ground, 0, 1.0}}, R));
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
            CHECK_THAT(feature.signature["RadiusOverR"].get<double>(), WithinAbs(3.0, 1.0e-9));
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
    // circles): no exact arc reading, the chord reading stays within the signature parameter
    // tolerance of the design offsets at both discretisations (the recorded ambiguity).
    for (const double step : {5.0, 2.5})
    {
      const double centre = 3.0 * R + 1.5;
      const auto trace = ArcBar(3.0, centre, 90.0, step, 12.0, 0.0);
      const auto ground = ArcBar(8.0, centre, 90.0, step, 12.0, -(1.5 + 2.0 + 4.0));
      const auto result = IdentifyMetalPerimeter(MakeInput({{trace, 0, 1.0}, {ground, 0, 1.0}}, R));
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
