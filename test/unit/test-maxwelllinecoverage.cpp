// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
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
    REQUIRE_NOTHROW(coverage = ResolveMaxwellLineCoverage(
                        element_intervals, {}, crack.line.length, crack.coordinate_scale,
                        crack.line_index, crack.line.begin, crack.line.end));
    REQUIRE(coverage.intervals.size() == 2);
    CHECK(coverage.intervals[0] == Interval{0.0, crack.gap_begin});
    CHECK(coverage.intervals[1] == Interval{crack.gap_end, 1.0});
    CHECK(coverage.closed_crack_count == 1);
    CHECK_THAT(coverage.max_closed_crack, WithinRel(gap * crack.line.length, 1.0e-12));
    CHECK(coverage.max_closed_crack < kMaxwellLineCrackUlps *
                                          std::numeric_limits<double>::epsilon() *
                                          crack.coordinate_scale);
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
  CHECK_THROWS_WITH(ResolveMaxwellLineCoverage(element_intervals, {}, line.length, 6.75,
                                               297, line.begin, line.end),
                    ContainsSubstring("Maxwell response contour line 297 leaves the "
                                      "finite-element mesh over parameter"));

  // The boundary of the rule in physical length: a gap of 1.5 x the tolerance fails closed
  // (with the gap and the tolerance named), a gap of 0.5 x the tolerance is a crack; and
  // the tolerance is PHYSICAL: the same 1.5 x gap in parameter on a 100x longer line is 150
  // x the tolerance and fails too, while the 0.5 x crack on a 100x shorter line still
  // closes.
  const double coordinate_scale = 6.75;
  const double tolerance =
      kMaxwellLineCrackUlps * std::numeric_limits<double>::epsilon() * coordinate_scale;
  const double split = 0.3;
  auto Gapped = [&](double gap_length, double line_length)
  {
    std::vector<Interval> intervals = {{0.0, split},
                                       {split + gap_length / line_length, 1.0}};
    return intervals;
  };
  CHECK_THROWS_WITH(ResolveMaxwellLineCoverage(Gapped(1.5 * tolerance, line.length), {},
                                               line.length, coordinate_scale, 297,
                                               line.begin, line.end),
                    ContainsSubstring("above the floating-point crack tolerance"));
  MaxwellLineCoverage coverage;
  REQUIRE_NOTHROW(coverage = ResolveMaxwellLineCoverage(
                      Gapped(0.5 * tolerance, line.length), {}, line.length,
                      coordinate_scale, 297, line.begin, line.end));
  CHECK(coverage.closed_crack_count == 1);
  // (1e-6: the fixture's split + gap / L - split rounds at eps x split / (gap / L) ~ 3e-9.)
  CHECK_THAT(coverage.max_closed_crack, WithinRel(0.5 * tolerance, 1.0e-6));
  REQUIRE(coverage.intervals.size() == 2);
  const auto longer = VerticalLine(line.begin[0], line.begin[1], 0.0, 100.0 * line.length);
  CHECK_THROWS_WITH(ResolveMaxwellLineCoverage(Gapped(1.5 * tolerance, line.length), {},
                                               longer.length, coordinate_scale, 297,
                                               longer.begin, longer.end),
                    ContainsSubstring("above the floating-point crack tolerance"));
  const auto shorter = VerticalLine(line.begin[0], line.begin[1], 0.0, 0.01 * line.length);
  REQUIRE_NOTHROW(coverage = ResolveMaxwellLineCoverage(
                      Gapped(0.5 * tolerance, line.length), {}, shorter.length,
                      coordinate_scale, 297, shorter.begin, shorter.end));
  CHECK(coverage.closed_crack_count == 1);
}

TEST_CASE("SurfaceResponseOperatorMaxwellLineCoverageLocalElementDiameter",
          "[surfaceresponseoperator][maxwell][Serial]")
{
  // OI-MX-CONTOUR-COVERAGE-NC-REFINEMENT (idfix DESIGN section 7; decision 1072 Q2): a
  // hanging-face closure crack of a nonconforming refinement is O(eps x h). The tolerance
  // of a piece = max(4096 ULP x scale, 1e-6 x the local element diameter), the diameter =
  // the max over the elements of the intervals bounding the piece.
  const auto line = VerticalLine(-6.6490033793415808e+00, 3.6334089021739100e+00,
                                 -1.3509670079635962e-03, -1.5010744532928826e-03);
  const double coordinate_scale = 6.75;
  const double ulp_tolerance =
      kMaxwellLineCrackUlps * std::numeric_limits<double>::epsilon() * coordinate_scale;
  const double split = 0.3;
  auto Gapped = [&](double gap_length)
  {
    std::vector<Interval> intervals = {{0.0, split},
                                       {split + gap_length / line.length, 1.0}};
    return intervals;
  };
  // A crack of 20x the ULP tolerance (the S1b-rc record: 3.06e-11 vs 5.69e-13 = 54x): the
  // ULP rule alone fails closed; with the bounding elements of diameter h it is closed iff
  // 1e-6 x h >= the gap. The diameter of a NON-bounding interval does not count.
  const double gap = 20.0 * ulp_tolerance;
  CHECK_THROWS_WITH(ResolveMaxwellLineCoverage(Gapped(gap), {}, line.length,
                                               coordinate_scale, 744, line.begin, line.end),
                    ContainsSubstring("above the floating-point crack tolerance"));
  MaxwellLineCoverage coverage;
  const double h_closing = gap / kMaxwellLineCrackLocalDiameterFraction * 1.5;
  REQUIRE_NOTHROW(coverage = ResolveMaxwellLineCoverage(
                      Gapped(gap), {}, line.length, coordinate_scale, 744, line.begin,
                      line.end, {h_closing, 0.1 * h_closing}));
  CHECK(coverage.closed_crack_count == 1);
  CHECK_THAT(coverage.max_closed_crack, WithinRel(gap, 1.0e-6));
  REQUIRE(coverage.intervals.size() == 2);
  // The larger of the two bounding elements sets the scale.
  REQUIRE_NOTHROW(coverage = ResolveMaxwellLineCoverage(
                      Gapped(gap), {}, line.length, coordinate_scale, 744, line.begin,
                      line.end, {0.1 * h_closing, h_closing}));
  CHECK(coverage.closed_crack_count == 1);
  // Elements too small to explain the crack: fails closed, the message naming both terms.
  const double h_failing = gap / kMaxwellLineCrackLocalDiameterFraction * 0.5;
  CHECK_THROWS_WITH(ResolveMaxwellLineCoverage(Gapped(gap), {}, line.length,
                                               coordinate_scale, 744, line.begin, line.end,
                                               {h_failing, h_failing}),
                    ContainsSubstring("x the local element diameter"));
  // A large element NOT bounding the crack (a third interval elsewhere) does not count.
  std::vector<Interval> three = {
      {0.0, 0.1}, {0.1, split}, {split + gap / line.length, 1.0}};
  CHECK_THROWS_WITH(ResolveMaxwellLineCoverage(three, {}, line.length, coordinate_scale,
                                               744, line.begin, line.end,
                                               {h_closing, h_failing, h_failing}),
                    ContainsSubstring("above the floating-point crack tolerance"));
  // Mismatched diameter counts are refused by name.
  CHECK_THROWS_WITH(ResolveMaxwellLineCoverage(Gapped(gap), {}, line.length,
                                               coordinate_scale, 744, line.begin, line.end,
                                               {h_closing}),
                    ContainsSubstring("element diameters must parallel the intervals"));
  // Where no crack exists the diameters change nothing: the covered resolution is bitwise.
  std::vector<Interval> covered = {{0.0, split}, {split, 1.0}};
  const auto without = ResolveMaxwellLineCoverage(
      covered, {}, line.length, coordinate_scale, 744, line.begin, line.end);
  const auto with =
      ResolveMaxwellLineCoverage(covered, {}, line.length, coordinate_scale, 744,
                                 line.begin, line.end, {h_closing, h_closing});
  CHECK(with.intervals == without.intervals);
  CHECK(with.closed_crack_count == 0);
}

TEST_CASE("SurfaceResponseOperatorMaxwellLineCoverageUncoveredReported",
          "[surfaceresponseoperator][maxwell][geometryfallback][h14a][Serial]")
{
  // H14a (hardening-2 DESIGN; idfix DESIGN section 7; decisions 1072 Q2 / 1074): a line
  // that leaves the mesh beyond the crack tolerance is REPORTED (the former abort's text,
  // the gap, the tolerance) when the caller asks, so that the operator can exclude the
  // patch; without the request the abort stands; a covered line reports nothing and
  // resolves bitwise with or without the request.
  const auto line = VerticalLine(-6.6490033793415808e+00, 3.6334089021739100e+00,
                                 -1.3509670079635962e-03, -1.5010744532928826e-03);
  const double coordinate_scale = 6.75;
  const double ulp_tolerance =
      kMaxwellLineCrackUlps * std::numeric_limits<double>::epsilon() * coordinate_scale;
  const double gap = 1.0e-6 * line.length;  // the line pushed 1e-6 of its length outside
  REQUIRE(gap > ulp_tolerance);
  std::vector<Interval> gapped = {{0.0, 0.3}, {0.3 + gap / line.length, 1.0}};
  // Without the request: the abort.
  CHECK_THROWS_WITH(ResolveMaxwellLineCoverage(gapped, {}, line.length, coordinate_scale,
                                               744, line.begin, line.end),
                    ContainsSubstring("leaves the finite-element mesh"));
  // With it: no throw, the report carries the same text and the numbers.
  MaxwellLineUncovered uncovered;
  MaxwellLineCoverage coverage;
  REQUIRE_NOTHROW(coverage = ResolveMaxwellLineCoverage(gapped, {}, line.length,
                                                        coordinate_scale, 744, line.begin,
                                                        line.end, {}, {}, &uncovered));
  CHECK(uncovered.Reported());
  CHECK_THAT(uncovered.message, ContainsSubstring("contour line 744 leaves the "
                                                  "finite-element mesh"));
  CHECK_THAT(uncovered.message, ContainsSubstring("above the floating-point crack "
                                                  "tolerance"));
  CHECK_THAT(uncovered.gap, WithinRel(gap, 1.0e-6));
  CHECK_THAT(uncovered.tolerance, WithinRel(ulp_tolerance, 1.0e-12));
  // The resolution stopped at the uncovered piece: the covered piece before it stands.
  REQUIRE(coverage.intervals.size() == 1);
  CHECK(coverage.intervals.front() == Interval{0.0, 0.3});
  // A covered line: nothing reported, bitwise with and without the request.
  std::vector<Interval> covered = {{0.0, 0.3}, {0.3, 1.0}};
  MaxwellLineUncovered none;
  const auto with = ResolveMaxwellLineCoverage(covered, {}, line.length, coordinate_scale,
                                               744, line.begin, line.end, {}, {}, &none);
  const auto without = ResolveMaxwellLineCoverage(
      covered, {}, line.length, coordinate_scale, 744, line.begin, line.end);
  CHECK(!none.Reported());
  CHECK(with.intervals == without.intervals);
  CHECK(with.closed_crack_count == without.closed_crack_count);
}

TEST_CASE("SurfaceResponseOperatorMaxwellLineCoverageBitwise",
          "[surfaceresponseoperator][maxwell][Serial]")
{
  // T3: where no uncovered piece exists the resolution is the pre-existing one to the bit:
  // ends within 1e-10 of 0 / 1 snap, breakpoints within 1e-10 of the last kept one are
  // dropped, overlaps and duplicates resolve to the breakpoint pieces.
  const auto line = VerticalLine(0.5, 0.25, 0.0, 2.0);
  std::vector<Interval> element_intervals = {{3.0e-11, 0.5 + 2.0e-11},
                                             {0.5 - 4.0e-11, 1.0 - 2.0e-11},
                                             {0.5 - 4.0e-11, 0.75},
                                             {0.75, 0.75 + 5.0e-11}};
  MaxwellLineCoverage coverage;
  REQUIRE_NOTHROW(coverage = ResolveMaxwellLineCoverage(element_intervals, {}, line.length,
                                                        2.0, 7, line.begin, line.end));
  // Breakpoints after snapping / uniquing: 0, 0.5 - 4e-11 (the first of the 0.5 cluster),
  // 0.75, 0.75 + 5e-11 is within 1e-10 of 0.75 (dropped), 1.
  REQUIRE(coverage.intervals.size() == 3);
  CHECK(coverage.intervals[0] == Interval{0.0, 0.5 - 4.0e-11});
  CHECK(coverage.intervals[1] == Interval{0.5 - 4.0e-11, 0.75});
  CHECK(coverage.intervals[2] == Interval{0.75, 1.0});
  CHECK(coverage.closed_crack_count == 0);
  CHECK(coverage.max_closed_crack == 0.0);
}

TEST_CASE("SurfaceResponseOperatorMaxwellLineCoverageInFace",
          "[surfaceresponseoperator][maxwell][Serial]")
{
  // T4 (decisions 829 / 831, lane inplane): the S1 stage-1 line of record (the e713e17dd494
  // anchor line in the trace_331_L2 sheet plane z = 4.80 um; PBS 63913 line 234, locally
  // line 10665): its exact intersections cover [0, 0.2943] and [0.8609, 1] plus five
  // zero-length slivers (the tets with two sheet vertices), the cracked sheet's tets
  // reject the in-crack stretch 0.2943-0.8609 from both sides. The in-face intersections
  // (the same tets with the sheet face taken as coplanar) cover that stretch in element
  // pieces, overlapping from the two sides: the resolution integrates it in the pieces
  // between their breakpoints, counts the stretch and its length, and the exact pieces are
  // untouched.
  const Line line = {{0.54346795625350652, -0.063844533280995167, 0.004550625711035268},
                     {0.54722248990917588, -0.052038301104279205, 0.0045506257110352671},
                     1.2389e-02};
  const double coordinate_scale = 1.0;
  const std::vector<Interval> exact = {{0.0, 0.29428948215455036},
                                       {0.42631914836662643, 0.42631914836662993},
                                       {0.5158780723205124, 0.5158780723205147},
                                       {0.767823476863604, 0.767823476863607},
                                       {0.8609476618148033, 1.0}};
  // The upper (substrate) tets' passages and the lower (vacuum) tets' passages.
  const std::vector<Interval> in_face = {
      {0.2943, 0.4263}, {0.4263, 0.5159}, {0.5159, 0.7678}, {0.7678, 0.8609},
      {0.2741, 0.4263}, {0.4263, 0.5159}, {0.5159, 0.7679}, {0.7678, 0.8610},
      {0.0901, 0.2741}, {0.8609, 1.0}};
  MaxwellLineCoverage coverage;
  REQUIRE_NOTHROW(coverage = ResolveMaxwellLineCoverage(exact, in_face, line.length,
                                                        coordinate_scale, 234, line.begin,
                                                        line.end));
  CHECK(coverage.closed_crack_count == 0);
  CHECK(coverage.in_face_count == 1);
  // The stretch between the last exact breakpoint before the gap and the first after it.
  CHECK_THAT(coverage.in_face_length,
             WithinRel((0.8609476618148033 - 0.29428948215455036) * line.length, 1.0e-12));
  // The exact pieces: [0, 0.2943] and [0.8609, 1] (the slivers are below the parameter
  // tolerance and skipped, their breakpoints kept); the in-face pieces split at the in-face
  // breakpoints inside the gap: 0.2943 | 0.4263 | 0.5159 | 0.7678 | 0.7679 | 0.8609.
  REQUIRE(coverage.intervals.size() >= 6);
  CHECK(coverage.intervals.front() == Interval{0.0, 0.29428948215455036});
  CHECK(coverage.intervals.back() == Interval{0.8609476618148033, 1.0});
  double covered_length = 0.0;
  for (std::size_t i = 0; i < coverage.intervals.size(); i++)
  {
    const auto &[begin, end] = coverage.intervals[i];
    CHECK(end > begin);
    covered_length += end - begin;
    if (i > 0)
    {
      CHECK(begin >= coverage.intervals[i - 1].second);
    }
  }
  // Every piece of the line is integrated exactly once (the slivers excepted).
  CHECK_THAT(covered_length, WithinAbs(1.0, 1.0e-9));
  for (const auto &[begin, end] : coverage.intervals)
  {
    if (begin < 0.2942 || end > 0.8610)
    {
      continue;  // An exact piece.
    }
    // No in-face piece straddles an in-face breakpoint strictly inside it.
    for (const auto &[face_begin, face_end] : in_face)
    {
      for (const double bound : {face_begin, face_end})
      {
        CHECK_FALSE((bound > begin + 1.0e-10 && bound < end - 1.0e-10));
      }
    }
  }

  // T5: a gap the in-face intersections do not cover still fails closed with the message
  // of record (a true gap, here the stretch 0.5159-0.7678 left out).
  std::vector<Interval> partial = in_face;
  partial.erase(std::remove_if(partial.begin(), partial.end(), [](const Interval &interval)
                               { return interval.first > 0.51 && interval.first < 0.52; }),
                partial.end());
  CHECK_THROWS_WITH(ResolveMaxwellLineCoverage(exact, partial, line.length,
                                               coordinate_scale, 234, line.begin, line.end),
                    ContainsSubstring("Maxwell response contour line 234 leaves the "
                                      "finite-element mesh over parameter 5.159"));
  CHECK_THROWS_WITH(ResolveMaxwellLineCoverage(exact, partial, line.length,
                                               coordinate_scale, 234, line.begin, line.end),
                    ContainsSubstring("no element intersection, exact or in a coplanar "
                                      "face, covers it"));
  CHECK_THROWS_WITH(ResolveMaxwellLineCoverage(exact, {}, line.length, coordinate_scale,
                                               234, line.begin, line.end),
                    ContainsSubstring("above the floating-point crack tolerance"));

  // T6: where the exact intervals leave no true gap, in-face intervals change nothing: the
  // T3 line with in-face intersections of its own (an uncracked interior face both tets
  // keep) resolves to the bit as without them.
  const auto bitwise_line = VerticalLine(0.5, 0.25, 0.0, 2.0);
  const std::vector<Interval> bitwise_exact = {{3.0e-11, 0.5 + 2.0e-11},
                                               {0.5 - 4.0e-11, 1.0 - 2.0e-11},
                                               {0.5 - 4.0e-11, 0.75},
                                               {0.75, 0.75 + 5.0e-11}};
  const std::vector<Interval> bitwise_in_face = {{0.1, 0.3}, {0.3, 0.9}, {0.45, 0.6}};
  MaxwellLineCoverage without, with;
  REQUIRE_NOTHROW(
      without = ResolveMaxwellLineCoverage(bitwise_exact, {}, bitwise_line.length, 2.0, 7,
                                           bitwise_line.begin, bitwise_line.end));
  REQUIRE_NOTHROW(with = ResolveMaxwellLineCoverage(bitwise_exact, bitwise_in_face,
                                                    bitwise_line.length, 2.0, 7,
                                                    bitwise_line.begin, bitwise_line.end));
  CHECK(with.intervals == without.intervals);
  CHECK(with.closed_crack_count == without.closed_crack_count);
  CHECK(with.max_closed_crack == without.max_closed_crack);
  CHECK(with.in_face_count == 0);
  CHECK(with.in_face_length == 0.0);
  // A floating-point crack stays a closed crack (not an in-face stretch) with in-face
  // intervals present.
  for (const auto &crack : kCracksOfRecord)
  {
    const std::vector<Interval> element_intervals = {{0.0, crack.gap_begin},
                                                     {crack.gap_end, 1.0}};
    const std::vector<Interval> crack_in_face = {{0.0, 1.0}};
    MaxwellLineCoverage crack_coverage;
    REQUIRE_NOTHROW(crack_coverage = ResolveMaxwellLineCoverage(
                        element_intervals, crack_in_face, crack.line.length,
                        crack.coordinate_scale, crack.line_index, crack.line.begin,
                        crack.line.end));
    CHECK(crack_coverage.closed_crack_count == 1);
    CHECK(crack_coverage.in_face_count == 0);
    REQUIRE(crack_coverage.intervals.size() == 2);
  }
}

}  // namespace palace
