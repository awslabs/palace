// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <utility>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include "models/surfaceresponsesegmentindex.hpp"

using namespace palace;

namespace
{

using Point3D = std::array<double, 3>;

// A chorded sub-segment in the placement's form: p0, the unit tangent, the parameter
// range [a, b] = [0, length] (a FramedPortion covering its whole segment).
struct Chord
{
  Point3D p0{}, tangent{};
  double a = 0.0, b = 0.0;
};

// A CPW-like meander in the plane z = 0.5: `turns` periods of a straight run of `run` along
// x (`run_chords` chords) followed by a half-turn of radius `radius` (`turn_chords`
// chords), (run_chords + turn_chords) x turns chords in all (the chain of a long chorded
// resonator edge: QUO-023's class). Coordinates in mm-scale units offset from the origin so
// that rounding is of the record's magnitude.
std::vector<Chord> Meander(int turns, int turn_chords, double run, double radius,
                           int run_chords)
{
  std::vector<Point3D> points;
  const Point3D origin = {12.5, -7.25, 0.5};
  double x = origin[0], y = origin[1];
  int direction = 1;
  points.push_back({x, y, origin[2]});
  for (int turn = 0; turn < turns; turn++)
  {
    for (int k = 1; k <= run_chords; k++)
    {
      points.push_back({x + direction * run * k / run_chords, y, origin[2]});
    }
    x += direction * run;
    const double cx = x, cy = y + radius;
    for (int k = 1; k <= turn_chords; k++)
    {
      const double angle = M_PI * k / turn_chords;
      points.push_back({cx + direction * radius * std::sin(angle),
                        cy - radius * std::cos(angle), origin[2]});
    }
    y += 2.0 * radius;
    direction = -direction;
  }
  std::vector<Chord> chords;
  for (std::size_t i = 1; i < points.size(); i++)
  {
    Chord c;
    c.p0 = points[i - 1];
    Point3D d = {points[i][0] - c.p0[0], points[i][1] - c.p0[1], points[i][2] - c.p0[2]};
    const double length = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    c.tangent = {d[0] / length, d[1] / length, d[2] / length};
    c.a = 0.0;
    c.b = length;
    chords.push_back(c);
  }
  return chords;
}

std::vector<std::pair<Point3D, Point3D>> Pieces(const std::vector<Chord> &chords)
{
  std::vector<std::pair<Point3D, Point3D>> pieces;
  for (const auto &c : chords)
  {
    pieces.emplace_back(Point3D{c.p0[0] + c.a * c.tangent[0], c.p0[1] + c.a * c.tangent[1],
                                c.p0[2] + c.a * c.tangent[2]},
                        Point3D{c.p0[0] + c.b * c.tangent[0], c.p0[1] + c.b * c.tangent[1],
                                c.p0[2] + c.b * c.tangent[2]});
  }
  return pieces;
}

// The placement's foot distance expression (ClosestFoot: the parameter clamped to the
// portion, the point interpolated, the Euclidean distance).
double FootDistance(const Point3D &p, const Chord &c)
{
  const double t =
      std::clamp((p[0] - c.p0[0]) * c.tangent[0] + (p[1] - c.p0[1]) * c.tangent[1] +
                     (p[2] - c.p0[2]) * c.tangent[2],
                 c.a, c.b);
  const Point3D q = {c.p0[0] + t * c.tangent[0], c.p0[1] + t * c.tangent[1],
                     c.p0[2] + t * c.tangent[2]};
  return std::sqrt((p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]) +
                   (p[2] - q[2]) * (p[2] - q[2]));
}

// The linear scan the index replaces: the first chord attaining the minimum distance.
std::size_t ScanNearest(const Point3D &p, const std::vector<Chord> &chords)
{
  std::size_t winner = std::numeric_limits<std::size_t>::max();
  double best = std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < chords.size(); i++)
  {
    const double distance = FootDistance(p, chords[i]);
    if (distance < best)
    {
      best = distance;
      winner = i;
    }
  }
  return winner;
}

}  // namespace

TEST_CASE("NearestSubSegmentIndex returns the linear scan's nearest sub-segment",
          "[surfaceresponse][segmentindex][Serial]")
{
  // 25 periods x (24 run chords + 16 turn chords) = 1,000 chords of a meander at pitch 2 r.
  const auto chords = Meander(25, 16, 0.6, 0.012, 24);
  REQUIRE(chords.size() == 1000);
  const NearestSubSegmentIndex index(Pieces(chords));
  REQUIRE(index.Size() == chords.size());
  auto Query = [&](const Point3D &p)
  { return index.Nearest(p, [&](std::size_t i) { return FootDistance(p, chords[i]); }); };
  SECTION("points on the chain: the patch origins and the end-face feet of every chord")
  {
    std::size_t checked = 0;
    for (const auto &c : chords)
    {
      // Quadrature-like interior points, the chord's own ends (a tie between neighbouring
      // chords at a shared vertex: the scan's first index must win) and a point just off
      // the chord on its +y side.
      for (const double t : {0.1127 * c.b, 0.5 * c.b, 0.8873 * c.b, c.a, c.b})
      {
        const Point3D p = {c.p0[0] + t * c.tangent[0], c.p0[1] + t * c.tangent[1],
                           c.p0[2] + t * c.tangent[2]};
        CHECK(Query(p) == ScanNearest(p, chords));
        const Point3D off = {p[0], p[1] + 1.0e-3, p[2]};
        CHECK(Query(off) == ScanNearest(off, chords));
        checked += 2;
      }
    }
    CHECK(checked == 10 * chords.size());
  }
  SECTION("points across the meander: the other-side feet of a pair at every separation")
  {
    std::mt19937_64 generator(1112);
    std::uniform_real_distribution<double> ux(12.5 - 0.3, 12.5 + 0.9),
        uy(-7.25 - 0.1, -7.25 + 25 * 0.024 + 0.1), uz(0.5 - 0.05, 0.5 + 0.05);
    for (int k = 0; k < 4000; k++)
    {
      const Point3D p = {ux(generator), uy(generator), uz(generator)};
      CHECK(Query(p) == ScanNearest(p, chords));
    }
    // Far points (every chord beyond the first query radius; the doubling reaches them).
    for (const Point3D &p :
         {Point3D{-3.0, 2.0, 0.5}, Point3D{40.0, -7.0, 9.0}, Point3D{12.8, 100.0, 0.5}})
    {
      CHECK(Query(p) == ScanNearest(p, chords));
    }
  }
  SECTION("a chain of one piece and the degenerate point chain")
  {
    const NearestSubSegmentIndex one(Pieces({chords.front()}));
    const Point3D p = {12.4, -7.3, 0.5};
    CHECK(one.Nearest(p, [&](std::size_t i) { return FootDistance(p, chords[i]); }) == 0);
    const std::vector<std::pair<Point3D, Point3D>> points = {
        {{1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}}, {{1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}}};
    const NearestSubSegmentIndex degenerate(points);
    // Equal distances: the first index wins, as in the scan.
    CHECK(degenerate.Nearest(Point3D{2.0, 1.0, 1.0}, [](std::size_t) { return 1.0; }) == 0);
    // No finite distance: the scan's "no winner".
    CHECK(degenerate.Nearest(Point3D{2.0, 1.0, 1.0}, [](std::size_t)
                             { return std::numeric_limits<double>::quiet_NaN(); }) ==
          std::numeric_limits<std::size_t>::max());
  }
}

TEST_CASE("NearestSubSegmentIndex keeps the patch build's foot queries linear in the chain",
          "[surfaceresponse][segmentindex][Serial]")
{
  // The scaling assertion of OI-PLACEMENT-PASS-SCALING (decision 1112): the placement
  // queries ~3 feet per placed patch on the feature's own chain, one patch per quadrature
  // point per chord, so the chain of N chords is queried ~3 Q N times; a scan costs N
  // evaluations per query (3 Q N^2 in all, the 677-s QUO-023 build), the index O(1).
  // Between N and 4 N chords the total evaluations must stay within a factor 6 of linear (4
  // x) - the scan would read 16 x.
  auto Work = [](int turns)
  {
    const auto chords = Meander(turns, 16, 0.6, 0.012, 24);
    const NearestSubSegmentIndex index(Pieces(chords));
    std::size_t evaluated = 0;
    for (const auto &c : chords)
    {
      for (const double t : {0.1127 * c.b, 0.5 * c.b, 0.8873 * c.b})
      {
        const Point3D p = {c.p0[0] + t * c.tangent[0], c.p0[1] + t * c.tangent[1],
                           c.p0[2] + t * c.tangent[2]};
        // The origin and the two end-face points of the cell about p (the probe form).
        for (const double shift : {0.0, -0.4 * c.b / 3.0, 0.4 * c.b / 3.0})
        {
          const Point3D q = {p[0] + shift * c.tangent[0], p[1] + shift * c.tangent[1],
                             p[2] + shift * c.tangent[2]};
          const std::size_t nearest = index.Nearest(
              q, [&](std::size_t i) { return FootDistance(q, chords[i]); }, &evaluated);
          REQUIRE(nearest < chords.size());
        }
      }
    }
    return std::make_pair(chords.size(), evaluated);
  };
  const auto [n_small, work_small] = Work(50);   // 2,000 chords
  const auto [n_large, work_large] = Work(200);  // 8,000 chords
  REQUIRE(n_large == 4 * n_small);
  CHECK(work_large <= 6 * work_small);
  // The per-query work is bounded (a handful of candidates per foot), not O(N).
  CHECK(work_small <= 40 * 9 * n_small);
  CHECK(work_large <= 40 * 9 * n_large);
}
