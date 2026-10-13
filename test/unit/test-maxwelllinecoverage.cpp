// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "models/surfaceresponseoperator.hpp"

namespace palace
{

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace
{

using Interval = std::pair<double, double>;

// A vertical contour line of the given length (mesh units) at the given position.
struct Line
{
  std::array<double, 3> begin;
  std::array<double, 3> end;
  double length;
};

Line VerticalLine(double x, double y, double z0, double z1)
{
  return {{x, y, z0}, {x, y, z1}, std::abs(z1 - z0)};
}

// The two cracks of record (decision 745 (1); the lane contourcover REPORT 1.1): the
// exact element intersections of two tets sharing the face a 158-nm line crosses, whose
// independently rounded clips of that face leave a gap of a few tens of ULP of the
// coordinate magnitude. ctx003-C4 (PBS 62551): line 297, parameter 9.1975876448039284e-01
// to 9.1975876468043716e-01 on a line of length 1.5011e-04 (internal units, Lc = 1055 um;
// coordinate scale 6.75); osc003-O2 (PBS 62556): line 1018, 6.1012003333674247e-01 to
// 6.1012003360164113e-01 on 1.0383e-04 (Lc = 1525 um; coordinate scale 1.81).
struct CrackOfRecord
{
  const char *window;
  int line_index;
  Line line;
  double coordinate_scale;
  double gap_begin;
  double gap_end;
};

const std::array<CrackOfRecord, 2> kCracksOfRecord = {
    CrackOfRecord{"ctx003-C4", 297,
                  VerticalLine(-6.6490033793415808e+00, 3.6334089021739100e+00,
                               -1.3509670079635962e-03, -1.5010744532928826e-03),
                  6.75, 9.1975876448039284e-01, 9.1975876468043716e-01},
    CrackOfRecord{"osc003-O2", 1018,
                  VerticalLine(1.7546010928963933e+00, -1.2027177539972731e+00,
                               1.0382513661202186e-03, 1.1420765027322405e-03),
                  1.81, 6.1012003333674247e-01, 6.1012003360164113e-01}};

}  // namespace

TEST_CASE("SurfaceResponseOperatorMaxwellLineCoverageCrack",
          "[surfaceresponseoperator][maxwell][Serial]")
{
  // T1 of the lane REPORT: the two cracks of record, each a sub-ULP uncovered piece
  // between the intersections of two adjacent tets, are floating-point cracks (in physical
  // length 3.0e-14 / 2.8e-14 mesh units = 20 / 68 ULP of the coordinate scale) and must be
  // closed: the integration intervals are the two element intervals, the crack is counted
  // and its length reported; the run of record aborted here instead.
  for (const auto &crack : kCracksOfRecord)
  {
    INFO(crack.window);
    const double gap = crack.gap_end - crack.gap_begin;
    const double ulps = gap * crack.line.length /
                        (std::numeric_limits<double>::epsilon() * crack.coordinate_scale);
    CHECK(ulps > 10.0);
    CHECK(ulps < 100.0);
    std::vector<Interval> element_intervals = {{0.0, crack.gap_begin},
                                               {crack.gap_end, 1.0}};
    MaxwellLineCoverage coverage;
    REQUIRE_NOTHROW(coverage = ResolveMaxwellLineCoverage(element_intervals,
                                                          crack.line_index,
                                                          crack.line.begin, crack.line.end));
    REQUIRE(coverage.intervals.size() == 2);
    CHECK(coverage.intervals[0] == Interval{0.0, crack.gap_begin});
    CHECK(coverage.intervals[1] == Interval{crack.gap_end, 1.0});
  }
}

TEST_CASE("SurfaceResponseOperatorMaxwellLineCoverageTrueGap",
          "[surfaceresponseoperator][maxwell][Serial]")
{
  // T2: a true gap (the line spends a finite length outside the mesh) still fails closed
  // with the message of record, naming the line and the parameter range.
  const auto line = VerticalLine(-6.6490033793415808e+00, 3.6334089021739100e+00,
                                 -1.3509670079635962e-03, -1.5010744532928826e-03);
  std::vector<Interval> element_intervals = {{0.0, 0.4}, {0.6, 1.0}};
  CHECK_THROWS_WITH(ResolveMaxwellLineCoverage(element_intervals, 297, line.begin, line.end),
                    ContainsSubstring("Maxwell response contour line 297 leaves the "
                                      "finite-element mesh over parameter"));
}

TEST_CASE("SurfaceResponseOperatorMaxwellLineCoverageBitwise",
          "[surfaceresponseoperator][maxwell][Serial]")
{
  // T3: where no uncovered piece exists the resolution is the pre-existing one to the bit:
  // ends within 1e-10 of 0 / 1 snap, breakpoints within 1e-10 of the last kept one are
  // dropped, overlaps and duplicates resolve to the breakpoint pieces.
  const auto line = VerticalLine(0.5, 0.25, 0.0, 2.0);
  std::vector<Interval> element_intervals = {
      {3.0e-11, 0.5 + 2.0e-11}, {0.5 - 4.0e-11, 1.0 - 2.0e-11}, {0.5 - 4.0e-11, 0.75},
      {0.75, 0.75 + 5.0e-11}};
  MaxwellLineCoverage coverage;
  REQUIRE_NOTHROW(coverage =
                      ResolveMaxwellLineCoverage(element_intervals, 7, line.begin, line.end));
  // Breakpoints after snapping / uniquing: 0, 0.5 - 4e-11 (the first of the 0.5 cluster),
  // 0.75, 0.75 + 5e-11 is within 1e-10 of 0.75 (dropped), 1.
  REQUIRE(coverage.intervals.size() == 3);
  CHECK(coverage.intervals[0] == Interval{0.0, 0.5 - 4.0e-11});
  CHECK(coverage.intervals[1] == Interval{0.5 - 4.0e-11, 0.75});
  CHECK(coverage.intervals[2] == Interval{0.75, 1.0});
}

}  // namespace palace
