// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
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
  return input;
}

std::vector<Point2> Rectangle(double x0, double y0, double x1, double y1)
{
  return {{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}};
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
