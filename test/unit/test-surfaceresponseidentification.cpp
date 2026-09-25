// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
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
};

IdentificationInput MakeInput(const std::vector<LoopSpec> &loops, double radius,
                              bool reverse_segment_order = false)
{
  IdentificationInput input;
  input.radius = radius;
  std::map<std::array<long long, 2>, std::size_t> vertex_index;
  auto Vertex = [&](const Point2 &p)
  {
    const std::array<long long, 2> key = {std::llround(p[0] * 1.0e9),
                                          std::llround(p[1] * 1.0e9)};
    auto it = vertex_index.find(key);
    if (it == vertex_index.end())
    {
      it = vertex_index.emplace(key, input.vertices.size()).first;
      IdentificationVertex vertex;
      vertex.coordinate = {p[0], p[1], 0.0};
      input.vertices.push_back(vertex);
    }
    return it->second;
  };
  struct Raw
  {
    Point2 a, b, outward;
    int conductor;
    int chain;
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
                       chain});
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
    segment.p0 = {r.a[0], r.a[1], 0.0};
    segment.p1 = {r.b[0], r.b[1], 0.0};
    segment.vertices = {Vertex(r.a), Vertex(r.b)};
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
std::vector<Point2> ArcBar(double width, double radius, double sweep_degrees,
                           double step_degrees, double lead = 6.0, double offset = 0.0)
{
  const int steps = std::max(1, static_cast<int>(std::lround(sweep_degrees / step_degrees)));
  const double step = sweep_degrees * std::acos(-1.0) / 180.0 / steps;
  const double h = 0.5 * width;
  std::vector<Point2> centreline;
  for (int k = 0; k <= steps; k++)
  {
    centreline.push_back({radius * std::sin(k * step), radius - radius * std::cos(k * step)});
  }
  const Point2 d_end = {std::cos(steps * step), std::sin(steps * step)};
  centreline.insert(centreline.begin(), {-lead, 0.0});
  centreline.push_back({centreline.back()[0] + lead * d_end[0],
                        centreline.back()[1] + lead * d_end[1]});
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
    if (vertex.type != "TruncationCut")
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
  // 8 x 6 island with 0.5 um fillets (radius < R): four rounded convex corners, one
  // isolated edge chain, no cluster; the same features when every chord is bisected
  // (refinement).
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
    CHECK(counts["IsolatedEdge"] == 1);
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
        // The tightest windowed bend radius is the inner side's (radius - 1.5).
        CHECK_THAT(feature.signature["RadiusOverR"].get<double>(),
                   WithinRel((c.radius - 1.5) / R, 0.03));
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

TEST_CASE("SurfaceResponseIdentificationPairsAtTheThreshold",
          "[surfaceresponseidentification][Serial]")
{
  // A gap of exactly 2R (and 2R +/- 1e-3 R) between two concentric 8 um bars (4 R: the far
  // corners of a bar end are beyond the 3R vertex-join reach of the near corners) along
  // bends of 50 and 250 um at three discretisations: the interaction decision uses the separation of
  // the underlying curves (rule at kPairSeparationTolerance), so every discretisation gives
  // the straight-pair answer — exactly 2R and 2R + 1e-3 R: isolated edges, no cluster, no
  // pair; 2R - 1e-3 R: one DifferentConductorGap along the bend (straight-like: inner radius
  // >= 10 R) and the two corner pairs across the gap as clusters. The mid-chord dips of the
  // chords below 2R (DS-SCT-001: 4 um gaps at 3.9999 -> 3 mm clusters) must not create
  // events.
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
            CHECK_THAT(feature.signature["SeparationOverR"].get<double>(),
                       WithinAbs(gap / R, 1.0e-6));
          }
        }
        const bool interacting = gap < 2.0 * R;
        CHECK(counts["DifferentConductorGap"] == (interacting ? 1 : 0));
        CHECK(counts["SpatialEdgeCluster"] == (interacting ? 2 : 0));
        CHECK(counts["ConvexCorner"] == (interacting ? 4 : 8));
        CHECK(counts["CurvedEdge"] == 0);
        CHECK(counts["CurvedDifferentConductorGap"] == 0);
        CHECK(counts["IsolatedEdge"] >= (interacting ? 2 : 4));
      }
    }
  }
}
