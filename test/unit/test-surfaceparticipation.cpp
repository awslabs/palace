// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include "drivers/basesolver.hpp"
#include "fixtures.hpp"
#include "models/surfaceparticipationtable.hpp"
#include "models/surfaceresponsecoverage.hpp"
#include "models/surfaceresponseresolution.hpp"
#include "utils/communication.hpp"
#include "utils/filesystem.hpp"

// CUT-PLAN step 4 (lane v1table): the per-iteration designer table
// surface-participation.csv. Its verdicts against the 19 test vectors of lane meshrule
// VERDICT-INPUT.md (real values of its adaptive ladders, at full precision), the edge
// cases, the verdict data (the built-in file and a bars file's Verdict member), the
// eigenmode match, the join of levels, the region rows' own verdicts, the file across
// levels and runs, and the log block.

using namespace palace;

namespace
{

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

struct VerdictVector
{
  int id;
  const char *value;
  int order;
  double unknowns, p, m90;
  std::vector<ParticipationHistoryEntry> history;  // (N_j, p_j, m90_j)
  double per_doubling;
  const char *convergence, *resolution;
};

// The 19 test vectors of lane meshrule VERDICT-INPUT.md s5 (final, MESHTEXT=OK;
// reads/verdict-vectors-v2.csv sha256-16 91aa716a5be8cb90; 10 real, 9 synthetic): order,
// N_k, p_k, the cell census m90_k, the value's earlier (N_j, p_j, m90_j), the expected g
// (%/doubling; NaN = not evaluated), convergence and resolution words.
const std::vector<VerdictVector> kVerdictVectors = {
    {1,
     "X3-ES-O3-A it 7 window 1/SA",
     3,
     23394544,
     4.933939224253e-04,
     1.972423,
     {{562169, 3.956801834422e-04, 14.621772},
      {597189, 4.244301615266e-04, 9.885531},
      {730037, 4.426736608846e-04, 8.609938},
      {1149243, 4.607205136411e-04, 6.683439},
      {2524845, 4.753070615050e-04, 4.415704},
      {7145399, 4.865797415173e-04, 2.985383}},
     0.818442,
     "SETTLED",
     "RESOLVED"},
    {2,
     "X3-ES-O4-A it 6 window 1/SA",
     4,
     17075627,
     4.934047321631e-04,
     3.126079,
     {{1314749, 4.417166954086e-04, 14.621772},
      {1403639, 4.570954831982e-04, 10.351422},
      {1725499, 4.677369887120e-04, 9.225714},
      {2729225, 4.787843719915e-04, 6.998420},
      {6012723, 4.872871550515e-04, 4.623810}},
     0.833707,
     "SETTLED",
     "RESOLVED"},
    {3,
     "X2-ES-O2-A it 9 window 2/MS",
     2,
     27460837,
     2.868663113993e-04,
     0.767361,
     {{165078, 9.583618647313e-05, 17.988709},
      {167154, 1.537438122667e-04, 13.645831},
      {176089, 2.019874903626e-04, 7.498942},
      {209548, 2.471624197300e-04, 3.935501},
      {326513, 2.681468976611e-04, 2.162719},
      {734708, 2.782534240316e-04, 1.798871},
      {2168941, 2.830762955051e-04, 1.109175},
      {7504567, 2.854210907633e-04, 0.944061}},
     0.270552,
     "SETTLED",
     "RESOLVED"},
    {4,
     "X1-ES-O2-A it 8 B1P30-L2 -/MS",
     2,
     21201741,
     3.790137614837e-05,
     1.109175,
     {{210096, 2.707056897786e-05, 15.310875},
      {216640, 2.863804576578e-05, 10.115795},
      {241177, 2.927455542412e-05, 8.222426},
      {339873, 3.116845814413e-05, 4.954502},
      {690631, 3.381074577132e-05, 3.427678},
      {1912228, 3.576288769927e-05, 2.317395},
      {6115696, 3.705217006789e-05, 1.640590}},
     1.277836,
     "MOVING",
     "MARGINAL"},
    {5,
     "X1-ES-O3-A it 5 B1P30-L1 -/SA",
     3,
     2878078,
     2.554947047588e-06,
     3.273407,
     {{690521, 2.336842154179e-06, 15.667511},
      {721821, 2.396186419960e-06, 10.115795},
      {849092, 2.449933189045e-06, 7.673615},
      {1307168, 2.518733719200e-06, 4.841724}},
     1.262675,
     "MOVING",
     "MARGINAL"},
    {6,
     "X2-ES-O2-ACF it 6 window 2/MS",
     2,
     14838704,
     2.790409840851e-04,
     0.683912,
     {{165078, 9.583618647313e-05, 17.988709},
      {208531, 1.485109245822e-04, 9.015711},
      {351609, 2.119279679357e-04, 4.518559},
      {859205, 2.581889588515e-04, 2.213095},
      {3124896, 2.712293239836e-04, 1.161449}},
     1.281474,
     "MOVING",
     "MARGINAL"},
    {7,
     "X2-ES-O2-A it 1 window 1/SA",
     2,
     165078,
     7.740901851509e-05,
     17.179084,
     {},
     kNaN,
     "MOVING: FIRST-ITERATION",
     "UNDER-RESOLVED"},
    {8,
     "X1-EIG-O2-A it 2 B1P30-L2 -/SA",
     2,
     960868,
     1.628651844008e-05,
     15.310875,
     {{960850, 1.628643576764e-05, 15.310875}},
     kNaN,
     "MOVING: UNKNOWNS-GREW-<20%",
     "UNDER-RESOLVED"},
    {9,
     "M3-ES-O3-I1 it 1 window 1/SA",
     3,
     1859973,
     4.874319633988e-04,
     1.678804,
     {},
     kNaN,
     "MOVING: FIRST-ITERATION",
     "MARGINAL"},
    {10,
     "X1-EIG-O3-A it 7 B1P30-L2 -/SA",
     3,
     3610755,
     1.889726386086e-05,
     9.015711,
     {{2774733, 1.781989019561e-05, 15.667511},
      {2774805, 1.781995457538e-05, 15.667511},
      {2774949, 1.781996454556e-05, 15.667511},
      {2775987, 1.782019432380e-05, 15.667511},
      {2850381, 1.824882139227e-05, 10.115795},
      {3058158, 1.864224706569e-05, 10.115795}},
     10.415896,
     "MOVING",
     "UNDER-RESOLVED"},
    {11,
     "synthetic: SETTLED with B < m90 <= 2 B: MARGINAL",
     3,
     2000000,
     1.005000000000e-03,
     3.000000,
     {{1000000, 1.000000000000e-03, 6.000000}},
     0.500000,
     "SETTLED",
     "MARGINAL"},
    {12,
     "synthetic: SETTLED with m90 > 2 B: UNDER-RESOLVED",
     3,
     2000000,
     1.005000000000e-03,
     5.000000,
     {{1000000, 1.000000000000e-03, 12.000000}},
     0.500000,
     "SETTLED",
     "UNDER-RESOLVED"},
    {13,
     "synthetic: g = T exactly (1.0 in IEEE double: 25 / 25): MOVING (strict ",
     3,
     134217728,
     1.220703125000e-03,
     1.500000,
     {{4, 9.765625000000e-04, 40.000000}},
     1.000000,
     "MOVING",
     "MARGINAL"},
    {14,
     "synthetic: m90 = B exactly with g < T: RESOLVED (<=)",
     3,
     2000000,
     1.005000000000e-03,
     2.000000,
     {{1000000, 1.000000000000e-03, 4.000000}},
     0.500000,
     "SETTLED",
     "RESOLVED"},
    {15,
     "synthetic: p_j* = 0: NOT EVALUATED, ZERO-ANCHOR-VALUE",
     3,
     2000000,
     1.000000000000e-06,
     1.000000,
     {{1000000, 0.000000000000e+00, 4.000000}},
     kNaN,
     "MOVING: ZERO-ANCHOR-VALUE",
     "MARGINAL"},
    {16,
     "synthetic: unknowns +50 % elsewhere, own cell unchanged: CELL-NOT-REFIN",
     3,
     1500000,
     1.000000000000e-03,
     1.900000,
     {{1000000, 1.000000000000e-03, 1.900000}},
     kNaN,
     "MOVING: CELL-NOT-REFINED",
     "MARGINAL"},
    {17,
     "synthetic: m90_k = 0.95 m90_j* exactly: guard passes (<=), RESOLVED",
     3,
     2000000,
     1.005000000000e-03,
     1.900000,
     {{1000000, 1.000000000000e-03, 2.000000}},
     0.500000,
     "SETTLED",
     "RESOLVED"},
    {18,
     "synthetic: Order 1: NO-BARS, indicator printed (SETTLED)",
     1,
     2000000,
     1.005000000000e-03,
     0.400000,
     {{1000000, 1.000000000000e-03, 0.800000}},
     0.500000,
     "SETTLED",
     "NO-BARS"},
    {19,
     "synthetic: Order 5: NO-BARS, indicator printed (MOVING)",
     5,
     2000000,
     1.030000000000e-03,
     0.400000,
     {{1000000, 1.000000000000e-03, 0.800000}},
     3.000000,
     "MOVING",
     "NO-BARS"},
};

std::string ReadText(const std::filesystem::path &path)
{
  std::ifstream input(path);
  std::stringstream text;
  text << input.rdbuf();
  return text.str();
}

std::size_t CountLines(const std::string &text)
{
  return static_cast<std::size_t>(std::count(text.begin(), text.end(), '\n'));
}

SurfaceParticipationRow Row(int point, int interface, const std::string &type, double p_raw,
                            double p_corrected, bool region = false)
{
  SurfaceParticipationRow row;
  row.point = point;
  row.interface = interface;
  row.type = type;
  row.region = region;
  row.p_raw = p_raw;
  row.p_corrected = p_corrected;
  row.m90 = 0.5;
  row.m50 = 0.4;
  row.coverage = "OK";
  row.regime = "IN";
  row.confidence = "pass";
  return row;
}

ParticipationVerdictBars FinalBars()
{
  auto bars = BuiltInParticipationVerdictBars();
  bars.interim = false;
  return bars;
}

ParticipationVerdictBars InterimBars()
{
  auto bars = BuiltInParticipationVerdictBars();
  bars.interim = true;
  return bars;
}

// The previous point a current mode joined, or -1 (UNMATCHED).
int JoinedPoint(const std::map<int, ParticipationModeJoin> &joins, int point)
{
  const auto &join = joins.at(point);
  return join.joined ? join.nearest : -1;
}

// The five modes of the CTX-002 window at Order 2 (lane v1table, run KA-X1-EIG-O2, level 1:
// Re{f} in GHz, p raw of the six whole interfaces of surface-participation.csv, |p| of the
// six lumped ports of port-EPR.csv). The participation vectors of the same field are nearly
// parallel: modes 1 / 3 / 5 overlap by 0.9976 to 0.9993, modes 2 / 4 by 0.9986 (the others
// 0.82 to 0.86); with the ports, modes 1 / 3 / 5 overlap by 0.02 to 0.33.
const std::vector<ParticipationModeSignature> kCtx002Modes = {
    {1,
     4.374736052950,
     {{1, 5.346838568654e-05},
      {2, 5.012625975262e-07},
      {3, 2.757530528807e-07},
      {4, 3.201838700776e-05},
      {5, 4.294285846984e-06},
      {6, 9.195304920505e-07}},
     {{64, 8.710769076445e-01},
      {65, 9.297285364748e-02},
      {102, 8.806376256590e-03},
      {103, 4.672794159397e-03},
      {106, 6.856868476300e-03},
      {107, 4.110182175216e-03}}},
    {2,
     5.178263320190,
     {{1, 2.916394745738e-04},
      {2, 8.413334801103e-16},
      {3, 3.226557865179e-12},
      {4, 7.690297228123e-11},
      {5, 8.697349620206e-13},
      {6, 7.431932309097e-06}},
     {{64, 7.742537815721e-12},
      {65, 9.463552451136e-13},
      {102, 1.033464678425e-13},
      {103, 6.295862344434e-14},
      {106, 1.068016472258e-13},
      {107, 1.022488263987e-13}}},
    {3,
     8.149317590894,
     {{1, 3.980883957009e-05},
      {2, 9.015407081454e-07},
      {3, 7.042438274160e-08},
      {4, 2.770620737767e-05},
      {5, 3.190179135131e-06},
      {6, 8.906149957107e-07}},
     {{64, 3.524733544933e-03},
      {65, 3.771214526963e-04},
      {102, 5.330668929498e-01},
      {103, 2.832213915417e-01},
      {106, 9.088414177334e-02},
      {107, 5.437473814845e-02}}},
    {4,
     8.401627481141,
     {{1, 2.883555069063e-04},
      {2, 2.850960474239e-12},
      {3, 1.727104660860e-12},
      {4, 2.657675987983e-10},
      {5, 1.314763608845e-11},
      {6, 2.247437334546e-05}},
     {{64, 6.863474646239e-09},
      {65, 7.138848223163e-10},
      {102, 1.193013135173e-09},
      {103, 6.328275548649e-10},
      {106, 1.858444104225e-06},
      {107, 1.111279022132e-06}}},
    {5,
     8.657838609138,
     {{1, 3.917593757516e-05},
      {2, 1.440679618318e-06},
      {3, 8.325456791519e-08},
      {4, 2.496962367196e-05},
      {5, 2.976096390633e-06},
      {6, 9.341052094816e-07}},
     {{64, 1.496896460086e-02},
      {65, 1.596508402998e-03},
      {102, 8.279085554538e-02},
      {103, 4.401219694469e-02},
      {106, 4.998011524404e-01},
      {107, 2.988103910093e-01}}}};

}  // namespace

TEST_CASE("Surface participation verdicts against the VERDICT-INPUT test vectors",
          "[surfaceparticipation][Serial]")
{
  const auto bars = FinalBars();
  const auto interim = InterimBars();
  REQUIRE(interim.built_in);
  for (const auto &v : kVerdictVectors)
  {
    INFO("vector " << v.id << ": " << v.value);
    const auto indicator =
        ComputeParticipationIndicator(v.p, v.m90, v.unknowns, v.history, bars);
    if (std::isnan(v.per_doubling))
    {
      CHECK(std::isnan(indicator.per_doubling));
    }
    else
    {
      CHECK(indicator.per_doubling == Catch::Approx(v.per_doubling).margin(1.0e-6));
    }
    CHECK(ParticipationConvergenceName(indicator, bars, "Electrostatic") == v.convergence);
    CHECK(ParticipationResolutionWord(v.m90, v.order, indicator.settled, false, bars,
                                      "Electrostatic") == v.resolution);
    // The interim data print the same rule with the INTERIM words.
    const auto word = ParticipationResolutionWord(v.m90, v.order, indicator.settled, false,
                                                  interim, "Electrostatic");
    CHECK(word == (std::string(v.resolution) == "RESOLVED" ? "RESOLVED-INTERIM"
                                                           : std::string(v.resolution)));
  }
  // Vector 10: the anchor is the latest earlier level with N_k >= 1.2 N_j (level 5; level 6
  // grew the unknowns by less than 20 %).
  const auto &v10 = kVerdictVectors[9];
  REQUIRE(v10.history.size() == 6);
  CHECK(ComputeParticipationIndicator(v10.p, v10.m90, v10.unknowns, v10.history, bars)
            .anchor == 4);
}

TEST_CASE("Surface participation verdict edge cases", "[surfaceparticipation][Serial]")
{
  using W = ParticipationNotEvaluated;
  using H = std::vector<ParticipationHistoryEntry>;
  const std::string es = "Electrostatic";
  const auto bars = FinalBars();
  CHECK(bars.min_unknowns_ratio == 1.2);
  CHECK(bars.max_census_ratio == 0.95);
  CHECK(bars.max_change_per_doubling == 1.0);
  CHECK(std::isnan(bars.CensusBar(1)));
  CHECK(bars.CensusBar(2) == 1.7);
  CHECK(bars.CensusBar(3) == 2.0);
  CHECK(bars.CensusBar(4) == 4.0);
  CHECK(std::isnan(bars.CensusBar(5)));
  CHECK(bars.under_resolved_factor == 2.0);
  CHECK(bars.tested_problems == std::vector<std::string>{"Electrostatic"});
  // No earlier level; none that much smaller; a non-finite value or p_j = 0; the value's
  // own census not fallen to 0.95 of the anchor's (or not read).
  CHECK(ComputeParticipationIndicator(1.0, 1.0, 1000, {}, bars).not_evaluated ==
        W::FIRST_ITERATION);
  const auto small =
      ComputeParticipationIndicator(1.0, 1.0, 1199, H{{1000, 1.0, 2.0}}, bars);
  CHECK(small.not_evaluated == W::UNKNOWNS_GREW_LITTLE);
  CHECK(ParticipationConvergenceName(small, bars, es) == "MOVING: UNKNOWNS-GREW-<20%");
  CHECK(ComputeParticipationIndicator(1.0, 1.0, 2000, H{{1000, 0.0, 2.0}}, bars)
            .not_evaluated == W::ZERO_ANCHOR_VALUE);
  CHECK(ComputeParticipationIndicator(kNaN, 1.0, 2000, H{{1000, 1.0, 2.0}}, bars)
            .not_evaluated == W::NO_VALUE);
  const auto stayed =
      ComputeParticipationIndicator(1.0, 1.91, 2000, H{{1000, 1.0, 2.0}}, bars);
  CHECK(stayed.not_evaluated == W::CELL_NOT_REFINED);
  CHECK(std::isnan(stayed.per_doubling));
  CHECK(ParticipationConvergenceName(stayed, bars, es) == "MOVING: CELL-NOT-REFINED");
  CHECK(ComputeParticipationIndicator(1.0, kNaN, 2000, H{{1000, 1.0, 2.0}}, bars)
            .not_evaluated == W::CELL_NOT_REFINED);
  CHECK(ComputeParticipationIndicator(1.0, 1.9, 2000, H{{1000, 1.0, 2.0}}, bars).settled);
  // A ratio of exactly G is read; the anchor is the LATEST such level.
  const auto exact =
      ComputeParticipationIndicator(1.0, 1.0, 1200, H{{1000, 1.0, 2.0}}, bars);
  CHECK(exact.not_evaluated == W::NONE);
  CHECK(exact.settled);
  const auto latest = ComputeParticipationIndicator(
      1.02, 1.0, 4000, H{{1000, 2.0, 4.0}, {2000, 1.0, 3.0}, {3900, 9.0, 2.0}}, bars);
  CHECK(latest.anchor == 1);
  CHECK(latest.per_doubling == Catch::Approx(2.0).epsilon(1e-12));
  CHECK(ParticipationConvergenceName(latest, bars, es) == "MOVING");
  // g = T exactly (25 % over 25 doublings) is MOVING (strict).
  const auto at_bar =
      ComputeParticipationIndicator(1.25, 1.0, 33554432, H{{1, 1.0, 2.0}}, bars);
  CHECK(at_bar.per_doubling == 1.0);
  CHECK_FALSE(at_bar.settled);
  CHECK(ComputeParticipationIndicator(1.005, 1.0, 4000, H{{2000, 1.0, 2.0}}, bars).settled);
  // The resolution word: no bar at Orders 1 and >= 5; no reading; the bars' edges (census =
  // bar with g < T is RESOLVED).
  CHECK(ParticipationResolutionWord(0.1, 1, true, false, bars, es) == "NO-BARS");
  CHECK(ParticipationResolutionWord(0.1, 5, true, false, bars, es) == "NO-BARS");
  CHECK(ParticipationResolutionWord(kNaN, 3, true, false, bars, es) == "NOT-EVALUATED");
  CHECK(ParticipationResolutionWord(2.0, 3, true, false, bars, es) == "RESOLVED");
  CHECK(ParticipationResolutionWord(2.0, 3, false, false, bars, es) == "MARGINAL");
  CHECK(ParticipationResolutionWord(4.0, 3, true, false, bars, es) == "MARGINAL");
  CHECK(ParticipationResolutionWord(4.000001, 3, false, false, bars, es) ==
        "UNDER-RESOLVED");
  CHECK(ParticipationResolutionWord(0.1, 2, false, false, bars, es) == "MARGINAL");
  // A single solve never reads RESOLVED.
  const auto first = ComputeParticipationIndicator(1.0, 0.1, 1000, {}, bars);
  CHECK(ParticipationResolutionWord(0.1, 3, first.settled, false, bars, es) == "MARGINAL");
  // The rule is untested on eigenmode and driven values: their words say so.
  CHECK(ParticipationResolutionWord(2.0, 3, true, false, bars, "Eigenmode") ==
        "RESOLVED-UNTESTED");
  // A region (box) value reads its own patches' census, which the bars were not graded
  // with: its graded words read -UNTESTED (once, also on an untested problem); the words
  // without a grading stay.
  CHECK(ParticipationResolutionWord(2.0, 3, true, true, bars, es) == "RESOLVED-UNTESTED");
  CHECK(ParticipationResolutionWord(2.0, 3, false, true, bars, es) == "MARGINAL-UNTESTED");
  CHECK(ParticipationResolutionWord(4.000001, 3, true, true, bars, es) ==
        "UNDER-RESOLVED-UNTESTED");
  CHECK(ParticipationResolutionWord(2.0, 3, true, true, bars, "Eigenmode") ==
        "RESOLVED-UNTESTED");
  CHECK(ParticipationResolutionWord(2.0, 3, true, true, InterimBars(), es) ==
        "RESOLVED-INTERIM-UNTESTED");
  CHECK(ParticipationResolutionWord(0.1, 1, true, true, bars, es) == "NO-BARS");
  CHECK(ParticipationResolutionWord(kNaN, 3, true, true, bars, es) == "NOT-EVALUATED");
  CHECK(ParticipationConvergenceName(exact, bars, "Driven") == "SETTLED-UNTESTED");
  CHECK(ParticipationConvergenceName(exact, bars, es) == "SETTLED");
  CHECK(ParticipationConvergenceName(exact, InterimBars(), es) == "SETTLED-INTERIM");
}

TEST_CASE("Surface participation verdict data from the built-in file and a bars file",
          "[surfaceparticipation][Serial]")
{
  const auto dir = std::filesystem::temp_directory_path() /
                   ("palace-participation-bars-" + std::to_string(Mpi::Rank(Mpi::World())));
  std::filesystem::create_directories(dir);
  const auto built_in = BuiltInParticipationVerdictBars();
  CHECK(built_in.built_in);
  CHECK_FALSE(built_in.interim);  // the supervisor's word once the 19 vectors passed
  CHECK(built_in.WordSuffix("Electrostatic").empty());
  CHECK(built_in.WordSuffix("Eigenmode") == "-UNTESTED");
  // The box values were graded with their cell's census: region words read -UNTESTED.
  CHECK(built_in.tested_scopes == std::vector<std::string>{"Type"});
  CHECK_FALSE(built_in.RegionScopeTested());
  CHECK(built_in.sha256.size() == 64);
  CHECK(built_in.provenance.find("VERDICT-INPUT") != std::string::npos);
  const std::string resolution = R"("RadiusMicrons": 1.9, "Bars": {
    "1": {"b5": null, "b10": 0.3}, "2": {"b5": 0.5, "b10": 1.2}, "3": {"b5": 0.6, "b10": 2.0},
    "4": {"b5": 1.0, "b10": 3.0}})";
  // No Verdict member: the built-in data; the resolution bars read as before.
  const auto plain = dir / "plain.json";
  std::ofstream(plain) << "{" << resolution << "}";
  const auto read_plain = ReadParticipationVerdictBars(plain.string());
  CHECK(read_plain.built_in);
  CHECK(read_plain.sha256 == built_in.sha256);
  // A Verdict member replaces every number; Interim defaults to false.
  const auto full = dir / "full.json";
  std::ofstream(full) << "{" << resolution << R"(, "Verdict": {"MinUnknownsRatio": 1.5,
    "MaxCensusRatio": 0.9, "MaxChangePerDoublingPercent": 0.5, "CensusBarOverR": {"1": 0.5,
    "2": 1.0, "3": 1.5, "4": null}, "UnderResolvedFactor": 3.0, "TestedProblems":
    ["Electrostatic", "Eigenmode"], "TestedScopes": ["Type", "Region"], "Provenance":
    "unit test"}})";
  const auto read = ReadParticipationVerdictBars(full.string());
  CHECK_FALSE(read.built_in);
  CHECK_FALSE(read.interim);
  CHECK(read.min_unknowns_ratio == 1.5);
  CHECK(read.max_change_per_doubling == 0.5);
  CHECK(read.CensusBar(1) == 0.5);
  CHECK(std::isnan(read.CensusBar(4)));
  CHECK(read.under_resolved_factor == 3.0);
  CHECK(read.max_census_ratio == 0.9);
  CHECK(read.WordSuffix("Eigenmode").empty());
  CHECK(read.WordSuffix("Driven") == "-UNTESTED");
  CHECK(read.RegionScopeTested());
  CHECK(read.provenance == "unit test");
  CHECK(read.sha256.size() == 64);
  CHECK(ReadResolutionBars(full.string()).B10(3) == 2.0);
  CHECK(ParticipationResolutionWord(4.4, 3, false, false, read, "Electrostatic") ==
        "MARGINAL");
  CHECK(ParticipationResolutionWord(4.6, 3, false, false, read, "Electrostatic") ==
        "UNDER-RESOLVED");
  // A file naming the Region scope: region words plain.
  CHECK(ParticipationResolutionWord(4.4, 3, false, true, read, "Electrostatic") ==
        "MARGINAL");
  CHECK(ParticipationResolutionWord(1.0, 3, true, true, read, "Eigenmode") == "RESOLVED");
  // Malformed members are refused.
  const std::string good_tail =
      R"("CensusBarOverR": {"1": null, "2": 1.7, "3": 2.0, "4": 4.0}, "UnderResolvedFactor": 2,
      "TestedProblems": ["Electrostatic"], "TestedScopes": ["Type"])";
  const std::string g12 = R"("MinUnknownsRatio": 1.2, )";
  const std::string g12_f_t =
      g12 + R"("MaxCensusRatio": 0.95, )" + R"("MaxChangePerDoublingPercent": 1.0, )";
  // The good member is read.
  {
    const auto good = dir / "good.json";
    std::ofstream(good) << "{" << resolution << R"(, "Verdict": {)" << g12_f_t << good_tail
                        << "}}";
    CHECK_FALSE(ReadParticipationVerdictBars(good.string()).RegionScopeTested());
  }
  for (
      const std::string member :
      {std::string(R"("MinUnknownsRatio": 1.0, "MaxCensusRatio": 0.95, )") +
           R"("MaxChangePerDoublingPercent": 1.0, )" + good_tail,
       g12 + R"("MaxCensusRatio": 0.95, "MaxChangePerDoublingPercent": 0.0, )" + good_tail,
       g12 + R"("MaxCensusRatio": 1.05, "MaxChangePerDoublingPercent": 1.0, )" + good_tail,
       g12 + R"("MaxChangePerDoublingPercent": 1.0, )" + good_tail,
       g12 + R"("MaxCensusRatio": 0.95, "MaxChangePerDoublingPercent": 1.0,
          "CensusBarOverR": {"1": null, "2": 1.7, "3": 2.0}, "UnderResolvedFactor": 2,
          "TestedProblems": ["Electrostatic"])",
       g12 + R"("MaxCensusRatio": 0.95, "MaxChangePerDoublingPercent": 1.0,
          "CensusBarOverR": {"1": null, "2": 1.7, "3": 2.0, "4": 4.0}, "UnderResolvedFactor": 2)",
       g12 +
           R"("MaxCensusRatio": 0.95, "MaxChangePerDoublingPercent": 1.0, "Interim": 1, )" +
           good_tail,
       g12_f_t + R"("CensusBarOverR": {"1": null, "2": 1.7, "3": 2.0, "4": 4.0},
          "UnderResolvedFactor": 2, "TestedProblems": ["Electrostatic"])",
       g12_f_t + R"("CensusBarOverR": {"1": null, "2": 1.7, "3": 2.0, "4": 4.0},
          "UnderResolvedFactor": 2, "TestedProblems": ["Electrostatic"],
          "TestedScopes": ["Type", "Box"])",
       g12_f_t + R"("CensusBarOverR": {"1": null, "2": 1.7, "3": 2.0, "4": 4.0},
          "UnderResolvedFactor": 2, "TestedProblems": ["Electrostatic"],
          "TestedScopes": [])",
       g12_f_t + R"("CensusBarOverR": {"1": null, "2": 1.7, "3": 2.0, "4": 4.0},
          "UnderResolvedFactor": 2, "TestedProblems": ["Electrostatic"],
          "TestedScopes": ["Region"])"})
  {
    const auto bad = dir / "bad.json";
    std::ofstream(bad) << "{" << resolution << R"(, "Verdict": {)" << member << "}}";
    CHECK_THROWS(ReadParticipationVerdictBars(bad.string()));
  }
  std::filesystem::remove_all(dir);
}

TEST_CASE("Surface participation eigenmode match across levels",
          "[surfaceparticipation][Serial]")
{
  using S = ParticipationModeSignature;
  const std::map<int, double> qubit = {{1, 0.9}, {2, 0.1}},
                              resonator = {{1, 0.1}, {2, 0.9}};
  // Two modes that swap their order between levels: each matched to its own partner.
  {
    const std::vector<S> previous = {{1, 5.0, qubit}, {2, 7.0, resonator}};
    const std::vector<S> current = {{1, 6.98, resonator}, {2, 7.4, qubit}};
    const auto match = MatchParticipationModes(previous, current);
    CHECK(JoinedPoint(match, 1) == 2);
    CHECK(JoinedPoint(match, 2) == -1);  // 7.4 GHz is nearest to mode 2, which mode 1 holds
    const std::vector<S> swapped = {{1, 7.01, resonator}, {2, 5.02, qubit}};
    const auto match_swapped = MatchParticipationModes(previous, swapped);
    CHECK(JoinedPoint(match_swapped, 1) == 2);
    CHECK(JoinedPoint(match_swapped, 2) == 1);
  }
  // A near-degenerate pair whose fields swap: the frequency rule pairs them crosswise and
  // the overlap refuses it: UNMATCHED, never the wrong mode.
  {
    const std::vector<S> previous = {{1, 5.000, qubit}, {2, 5.010, resonator}};
    const std::vector<S> current = {{1, 5.002, resonator}, {2, 5.008, qubit}};
    const auto match = MatchParticipationModes(previous, current);
    CHECK(JoinedPoint(match, 1) == -1);
    CHECK(JoinedPoint(match, 2) == -1);
  }
  // A mode that enters the computed window, and one beyond half the gap: UNMATCHED.
  {
    const std::vector<S> previous = {{1, 5.0, qubit}};
    const std::vector<S> current = {{1, 5.01, qubit}, {2, 6.5, resonator}};
    const auto match = MatchParticipationModes(previous, current);
    CHECK(JoinedPoint(match, 1) == 1);
    CHECK(JoinedPoint(match, 2) == -1);
    const std::vector<S> pair = {{1, 5.0, qubit}, {2, 5.2, qubit}};
    CHECK(JoinedPoint(MatchParticipationModes(pair, {{1, 5.35, qubit}}), 1) == -1);
    // One mode of that field left: the remaining one is UNMATCHED for this level (its rank
    // among the previous modes of its field is 2, among the current ones 1).
    const auto left = MatchParticipationModes(pair, {{1, 5.25, qubit}});
    CHECK(JoinedPoint(left, 1) == -1);
    CHECK(left.at(1).nearest == 2);
    CHECK(left.at(1).failed == "same-field rank differs");
    // Not unique: a second previous mode of the same field lies within twice the gap.
    const auto ambiguous = MatchParticipationModes(pair, {{1, 5.12, qubit}});
    CHECK(JoinedPoint(ambiguous, 1) == -1);
    CHECK(ambiguous.at(1).nearest == 2);
    CHECK(ambiguous.at(1).failed == "not unique");
  }
}

TEST_CASE("Surface participation eigenmode join on the real CTX-002 vectors",
          "[surfaceparticipation][Serial]")
{
  auto Shifted = [](std::vector<ParticipationModeSignature> modes, double shift)
  {
    for (auto &mode : modes)
    {
      mode.frequency += shift;
    }
    return modes;
  };
  for (const bool ports : {false, true})
  {
    INFO("fingerprint with the ports: " << ports);
    auto modes = kCtx002Modes;
    for (auto &mode : modes)
    {
      mode.ports = ports ? mode.ports : std::map<int, double>{};
    }
    // No shift and small shifts: every mode joins itself.
    for (const double shift : {0.0, 1.0e-3, -1.0e-3, 0.01, -0.01, 0.05, 0.1})
    {
      INFO("shift " << shift << " GHz");
      const auto joins = MatchParticipationModes(modes, Shifted(modes, shift));
      for (int point = 1; point <= 5; point++)
      {
        CHECK(JoinedPoint(joins, point) == point);
      }
    }
    // +0.40 GHz, beyond half the spacing of modes 3 and 5 (0.51 GHz; mode 4 lies between):
    // the rule without the rank condition joins mode 3 to mode 5 on the participations
    // alone; now mode 3 is UNMATCHED, 4 and 5 too (gap), 1 and 2 join themselves.
    {
      const auto joins = MatchParticipationModes(modes, Shifted(modes, 0.40));
      CHECK(JoinedPoint(joins, 1) == 1);
      CHECK(JoinedPoint(joins, 2) == 2);
      CHECK(joins.at(3).nearest == 5);
      CHECK(JoinedPoint(joins, 3) == -1);
      CHECK(joins.at(3).failed ==
            (ports ? "overlap below 0.9" : "same-field rank differs"));
      CHECK(JoinedPoint(joins, 4) == -1);
      CHECK(JoinedPoint(joins, 5) == -1);
    }
    // The sweep of the lane's tools (fix/tools/joinsim_v2.py, joinsim_ports_v1.py): 801
    // uniform shifts from 0 to 4 GHz, 4005 joins: 0 wrong (the rule without the rank
    // condition: 152 on the participations alone, 50 with the ports), 240 right and 3765
    // UNMATCHED with either fingerprint.
    int right = 0, wrong = 0, unmatched = 0;
    for (int k = 0; k <= 800; k++)
    {
      for (const auto &[point, join] :
           MatchParticipationModes(modes, Shifted(modes, 0.005 * k)))
      {
        (join.joined ? (join.nearest == point ? right : wrong) : unmatched)++;
      }
    }
    CHECK(wrong == 0);
    CHECK(right == 240);
    CHECK(unmatched == 3765);
    // A mode of the field of mode 3 enters above the window: the others keep their joins.
    {
      auto current = modes;
      current.push_back({6, 9.2, modes[2].participation, modes[2].ports});
      const auto joins = MatchParticipationModes(modes, current);
      for (int point = 1; point <= 5; point++)
      {
        CHECK(JoinedPoint(joins, point) == point);
      }
      CHECK(JoinedPoint(joins, 6) == -1);
    }
    // The lowest mode leaves the window (the others renumbered, +1 MHz): on the
    // participations alone modes 3 and 5 (the field of mode 1) are UNMATCHED for this
    // level; with the ports they keep their joins.
    {
      std::vector<ParticipationModeSignature> current(modes.begin() + 1, modes.end());
      for (auto &mode : current)
      {
        mode.point--;
        mode.frequency += 1.0e-3;
      }
      const auto joins = MatchParticipationModes(modes, current);
      CHECK(JoinedPoint(joins, 1) == 2);
      CHECK(JoinedPoint(joins, 2) == (ports ? 3 : -1));
      CHECK(JoinedPoint(joins, 3) == 4);
      CHECK(JoinedPoint(joins, 4) == (ports ? 5 : -1));
    }
    // The reviewer's case: two modes of one field (3 and 5) shifted beyond half their
    // spacing (+0.28 of 0.51 GHz): the old rule joined mode 3 to mode 5.
    {
      const std::vector<ParticipationModeSignature> pair = {modes[2], modes[4]};
      const auto joins = MatchParticipationModes(pair, Shifted(pair, 0.28));
      CHECK(JoinedPoint(joins, 3) == -1);
      CHECK(JoinedPoint(joins, 5) == -1);
      CHECK(joins.at(3).failed == (ports ? "overlap below 0.9" : "not unique"));
    }
  }
}

TEST_CASE("Surface participation join, history, words and flags across levels",
          "[surfaceparticipation][Serial]")
{
  const auto bars = InterimBars();
  auto Level = [](const std::string &problem, long long int unknowns,
                  std::vector<SurfaceParticipationRow> rows)
  {
    for (auto &row : rows)
    {
      row.problem = problem;
      row.unknowns = unknowns;
    }
    return rows;
  };
  // Electrostatic: two terminals; terminal 1 settles, terminal 2 moves; a raw-only row.
  auto first = Level("Electrostatic", 1000,
                     {Row(1, 1, "SA", 1.0e-4, 1.2e-4), Row(2, 1, "SA", 2.0e-4, 2.4e-4),
                      Row(1, 9, "SA", 1.0e-5, kNaN)});
  SetParticipationNotTarget(first[2]);
  JoinParticipationLevel(first, nullptr, bars, 3);
  CHECK(first[0].convergence == "MOVING: FIRST-ITERATION");
  CHECK(first[0].resolution == "MARGINAL");  // a single solve never reads RESOLVED
  CHECK(first[0].flags == "MARGINAL; MOVING: FIRST-ITERATION");
  CHECK(first[2].flags == "NOT-TARGET");
  CHECK(first[2].resolution == "NOT-TARGET");
  CHECK(first[2].reason == "NOT-TARGET: raw only");
  auto second = Level("Electrostatic", 1100,
                      {Row(1, 1, "SA", 1.1e-4, 1.2001e-4), Row(2, 1, "SA", 2.2e-4, 2.6e-4),
                       Row(1, 9, "SA", 1.1e-5, kNaN)});
  SetParticipationNotTarget(second[2]);
  second[1].m90 = 20.0;
  JoinParticipationLevel(second, &first, bars, 3);
  CHECK(second[0].change_corrected == 1.2001e-4 / 1.2e-4 - 1.0);
  CHECK(second[0].change_raw == 1.1e-4 / 1.0e-4 - 1.0);
  CHECK(std::isnan(second[0].per_doubling));
  CHECK(second[0].convergence == "MOVING: UNKNOWNS-GREW-<20%");
  REQUIRE(second[0].history.size() == 1);
  // The third level reads against the latest level with x1.2 fewer unknowns (the second).
  auto third = Level("Electrostatic", 4000,
                     {Row(1, 1, "SA", 1.1e-4, 1.2002e-4), Row(2, 1, "SA", 2.2e-4, 2.8e-4),
                      Row(1, 9, "SA", 1.1e-5, kNaN)});
  SetParticipationNotTarget(third[2]);
  third[0].m90 = 0.3;  // the value's own census fell since the anchor
  third[1].m90 = 9.0;
  third[1].regime = "NO-RANGE";
  third[1].marked = true;
  JoinParticipationLevel(third, &second, bars, 3);
  REQUIRE(third[0].history.size() == 2);
  CHECK(third[0].per_doubling ==
        Catch::Approx(100.0 * (1.2002e-4 / 1.2001e-4 - 1.0) / std::log2(4000.0 / 1100.0)));
  CHECK(third[0].convergence == "SETTLED-INTERIM");
  CHECK(third[0].resolution == "RESOLVED-INTERIM");
  CHECK(third[0].flags == "RESOLVED-INTERIM; SETTLED-INTERIM");
  CHECK(third[1].convergence == "MOVING");
  CHECK(third[1].resolution == "UNDER-RESOLVED");
  CHECK(third[1].flags ==
        "UNDER-RESOLVED; regime NO-RANGE; MOVING; developer: MarkNotApplied");
  CHECK(third[2].convergence == "MOVING: NO-VALUE");
  CHECK(third[2].flags == "NOT-TARGET");
  // Final (not interim) data: a row whose verdicts are all met has no flag.
  auto final_bars = bars;
  final_bars.interim = false;
  JoinParticipationLevel(third, &second, final_bars, 3);
  CHECK(third[0].convergence == "SETTLED");
  CHECK(third[0].resolution == "RESOLVED");
  CHECK(third[0].flags.empty());
  // CUT-PLAN step 3 (lane v1lsdefault): the impedance-metal rule's two flags words follow
  // the verdict words; a NOT-TARGET row keeps its one word.
  third[0].sheet_served = true;
  CHECK(ParticipationFlags(third[0]) == kParticipationSheetServedFlag);
  CHECK(ParticipationFlags(third[0]) == "Ls sheet served as PEC");
  third[0].sheet_not_served = true;
  CHECK(ParticipationFlags(third[0]) ==
        "Ls sheet served as PEC; metal law not served: raw there");
  third[0].sheet_served = false;
  CHECK(ParticipationFlags(third[0]) == kParticipationSheetNotServedFlag);
  third[1].sheet_served = true;
  CHECK(ParticipationFlags(third[1]) ==
        "UNDER-RESOLVED; regime NO-RANGE; MOVING; "
        "developer: MarkNotApplied; Ls sheet served as PEC");
  third[2].sheet_not_served = true;
  CHECK(ParticipationFlags(third[2]) == "NOT-TARGET");
  // The flag S (lane v1verdict1): S >= 5 % on corner and cluster models, more than 10 % on
  // cluster models (known shares of the row's census cell); no value changes.
  auto s_row = third[0];
  s_row.sheet_not_served = false;
  CHECK(ParticipationFlags(s_row).empty());
  ResolutionCell s_cell;
  s_cell.quantiles = {s_row.m50, s_row.m50, s_row.m90};
  s_cell.spatial_share = 0.05;
  s_cell.cluster_share = 0.10;
  SetParticipationResolution(s_row, &s_cell);
  CHECK(s_row.spatial_share == 0.05);
  CHECK(ParticipationFlags(s_row) == kParticipationSpatialShareFlag);
  CHECK(ParticipationFlags(s_row) == "S >= 5% on corner and cluster models");
  s_row.spatial_share = 0.0499;
  CHECK(ParticipationFlags(s_row).empty());
  s_row.spatial_share = 0.2;
  s_row.cluster_share = 0.1001;
  CHECK(ParticipationFlags(s_row) ==
        "S >= 5% on corner and cluster models; over 10% on cluster models");
  CHECK(s_row.p_corrected == third[0].p_corrected);
  SetParticipationNotTarget(s_row);
  CHECK(std::isnan(s_row.spatial_share));
  // A Type total carries its interface rows' words.
  const auto totals = AddParticipationTypeTotals({third[0]});
  REQUIRE(totals.size() == 2);
  CHECK(totals[0].type_total);
  CHECK(totals[0].sheet_not_served);
  third[0].sheet_not_served = false;
  third[1].sheet_served = false;
  third[2].sheet_not_served = false;
  // Driven: the same output-frequency index joins only at the same frequency and
  // excitation.
  auto driven_first = Level("Driven", 1000, {Row(1, 1, "SA", 1.0e-4, 1.2e-4)});
  driven_first[0].frequency = 5.0;
  driven_first[0].excitation = 1;
  JoinParticipationLevel(driven_first, nullptr, bars, 2);
  auto driven_second = Level("Driven", 3000, {Row(1, 1, "SA", 1.0e-4, 1.2e-4)});
  driven_second[0].m90 = 0.4;
  driven_second[0].frequency = 5.0;
  driven_second[0].excitation = 1;
  JoinParticipationLevel(driven_second, &driven_first, bars, 2);
  CHECK(driven_second[0].convergence == "SETTLED-INTERIM");
  driven_second[0].m90 = 0.5;  // the edges were not refined: no change means nothing
  JoinParticipationLevel(driven_second, &driven_first, bars, 2);
  CHECK(driven_second[0].convergence == "MOVING: CELL-NOT-REFINED");
  driven_second[0].m90 = 0.4;
  driven_second[0].frequency = 5.0 * (1.0 + 1.0e-9);
  JoinParticipationLevel(driven_second, &driven_first, bars, 2);
  CHECK(driven_second[0].convergence == "MOVING: UNMATCHED");
  CHECK(std::isnan(driven_second[0].change_corrected));
  CHECK(driven_second[0].history.empty());
  driven_second[0].frequency = 5.0;
  driven_second[0].excitation = 2;
  JoinParticipationLevel(driven_second, &driven_first, bars, 2);
  CHECK(driven_second[0].convergence == "MOVING: UNMATCHED");
  // Eigenmode: the rows of a mode follow its match (two modes that swap order).
  auto Mode = [&](int point, double f, double p1, double p2)
  {
    auto a = Row(point, 1, "SA", p1, p1), b = Row(point, 2, "MS", p2, p2);
    a.frequency = b.frequency = f;
    return std::vector<SurfaceParticipationRow>{a, b};
  };
  std::vector<SurfaceParticipationRow> eig_first, eig_second;
  for (const auto &rows : {Mode(1, 5.0, 0.9, 0.1), Mode(2, 7.0, 0.1, 0.9)})
  {
    eig_first.insert(eig_first.end(), rows.begin(), rows.end());
  }
  for (const auto &rows : {Mode(1, 7.001, 0.1, 0.9), Mode(2, 5.001, 0.99, 0.1)})
  {
    eig_second.insert(eig_second.end(), rows.begin(), rows.end());
  }
  eig_first = Level("Eigenmode", 1000, eig_first);
  eig_second = Level("Eigenmode", 5000, eig_second);
  for (auto &row : eig_second)
  {
    row.m90 = 0.3;
  }
  CHECK(JoinParticipationLevel(eig_first, nullptr, bars, 3).empty());
  const auto joins = JoinParticipationLevel(eig_second, &eig_first, bars, 3);
  CHECK(JoinedPoint(joins, 1) == 2);
  CHECK(JoinedPoint(joins, 2) == 1);
  CHECK(eig_second[0].change_corrected == Catch::Approx(0.0).margin(1e-15));  // 1 <- 2
  CHECK(eig_second[0].convergence == "SETTLED-INTERIM");
  CHECK(eig_second[2].change_corrected == Catch::Approx(0.99 / 0.9 - 1.0));  // 2 <- 1
  CHECK(eig_second[2].convergence == "MOVING");
  // The log block prints every join with the two frequencies and the overlap.
  const auto block = ParticipationLogBlock(eig_second, 2, 5000, 1000, bars, 3, joins);
  CHECK(
      block.find(" mode 1 <- mode 2 of iteration 1: 7 -> 7.001 GHz, overlap 1.000000\n") !=
      std::string::npos);
  CHECK(block.find(" mode 2 <- mode 1 of iteration 1: 5 -> 5.001 GHz, overlap 0.999") !=
        std::string::npos);
  // A third level whose mode 2 moved by more than half the spacing: UNMATCHED, printed so.
  auto eig_third = eig_second;
  for (auto &row : eig_third)
  {
    row.frequency = row.point == 2 ? 3.9 : row.frequency;
  }
  const auto third_joins = JoinParticipationLevel(eig_third, &eig_second, bars, 3);
  CHECK(JoinedPoint(third_joins, 2) == -1);
  CHECK(eig_third[2].convergence == "MOVING: UNMATCHED");
  CHECK(ParticipationLogBlock(eig_third, 3, 5000, 5000, bars, 3, third_joins)
            .find(" mode 2 UNMATCHED (gap not below half the spacing): nearest mode 2 of "
                  "iteration 2: 5.001 -> 3.9 GHz, overlap 1.000000\n") !=
        std::string::npos);
}

TEST_CASE("Surface participation region rows read their own patches",
          "[surfaceparticipation][Serial]")
{
  // The census: the box's patches at h/R = 17 (UNDER-RESOLVED at the built-in bars, Order
  // 3) hold 1 % of the Type's energy; the rest of the Type at h/R = 0.4 (RESOLVED). The
  // Type's cell reads RESOLVED, the box's own cell (the same patches, each with its share)
  // UNDER-RESOLVED: the GATE's 17 R box on a chip resolved elsewhere.
  const auto bars = BuiltInResolutionBars();
  const auto coarse = ResolutionPatchStatus(17.0, 3, bars, true);
  const auto fine = ResolutionPatchStatus(0.4, 3, bars, true);
  REQUIRE(coarse == ResolutionStatus::UNDER_RESOLVED);
  REQUIRE(fine == ResolutionStatus::RESOLVED);
  ResolutionAccumulator sums;
  sums.Reset(1, 1);
  for (int patch = 0; patch < 99; patch++)
  {
    sums.Add(0, 0, fine, 0.4, 1.0);
  }
  // Two box patches: one wholly inside (share 1), one cut by the wall (share 0.5).
  sums.Add(0, 0, coarse, 17.0, 0.6);
  sums.AddRegion(7, coarse, 17.0, 1.0 * 0.6);
  sums.Add(0, 0, coarse, 17.0, 0.8);
  sums.AddRegion(7, coarse, 17.0, 0.5 * 0.8);
  REQUIRE(sums.regions.size() == 1);
  const auto type_cell =
      ComputeResolutionCell(sums.Slot(0, 0), "SA", true, 0.0, true, 1.0, "");
  const auto box_cell =
      ComputeResolutionCell(sums.regions.at(7).data(), "SA", true, 0.0, true, 1.0, "");
  CHECK(type_cell.status == ResolutionStatus::RESOLVED);
  CHECK(box_cell.status == ResolutionStatus::UNDER_RESOLVED);
  CHECK(box_cell.patches == 2);
  CHECK(box_cell.energy == 1.0);
  CHECK(box_cell.shares[2] == 1.0);
  // Reset clears the box slots too.
  auto reset = sums;
  reset.Reset(1, 1);
  CHECK(reset.regions.empty());
  // The coverage: the whole SA interface 1 is served; the box interface 7 (same Type) has
  // half its W uncovered and no Missing key holds it (a method gap). The Type's cell reads
  // OK, the box's own cell NOT-ACHIEVABLE.
  std::map<int, SurfacePostOperator::InterfaceCoveragePartition> partitions;
  partitions[1].within = 1000.0;
  partitions[1].served = {{"m0", 1000.0}};
  partitions[7].within = 10.0;
  partitions[7].served = {{"m0", 5.0}};
  const auto coverage = ComputeCoverageMeasurement(
      partitions, {{1, "SA"}, {7, "SA"}}, {{1, 1000.0}, {7, 10.0}}, {{1, 1100.0}, {7, 6.0}},
      {}, {{0, CoverageModelQualification::QUALIFIED}}, CoverageThresholds{},
      CoverageValidatedRegime{}, CoverageExcitationClass::ELECTROSTATIC_TERMINAL);
  REQUIRE(coverage.interface_cells.size() == 2);
  // The rows: 'whole' from the Type's cells, 'region' from the box's own.
  auto whole = Row(1, 1, "SA", 1.0e-4, 1.1e-4);
  auto box = Row(1, 7, "SA", 1.0e-6, 1.6e-6, true);
  SetParticipationResolution(whole, &type_cell);
  SetParticipationCoverage(whole, coverage.FindCell("SA"));
  SetParticipationResolution(box, &box_cell);
  SetParticipationCoverage(box, &coverage.interface_cells.at(7));
  const auto final_bars = FinalBars();
  CHECK(whole.m90 == ResolutionHistogramCentre(ResolutionHistogramBin(0.4)));
  CHECK(ParticipationResolutionWord(whole.m90, 3, true, false, final_bars,
                                    "Electrostatic") == "RESOLVED");
  CHECK(whole.coverage == "OK");
  CHECK(whole.uncovered_share == Catch::Approx(5.0 / 1010.0));
  CHECK(box.m90 == ResolutionHistogramCentre(ResolutionHistogramBin(17.0)));
  CHECK(ParticipationResolutionWord(box.m90, 3, true, true, final_bars, "Electrostatic") ==
        "UNDER-RESOLVED-UNTESTED");
  CHECK(box.m50 == ResolutionHistogramCentre(ResolutionHistogramBin(17.0)));
  CHECK(box.coverage == "NOT-ACHIEVABLE");
  CHECK(box.uncovered_share == 0.5);
  // A target no model serves: its own cell reads TARGET-NOT-SERVED as its Type's does.
  auto not_served = coverage;
  MarkCoverageTargetsNotServed(not_served, {{7, "SA"}});
  CHECK(not_served.interface_cells.at(7).status == CoverageStatus::TARGET_NOT_SERVED);
  CHECK(not_served.interface_cells.at(1).status == CoverageStatus::OK);
  // No cell: NOT-EVALUATED / NOT-MEASURED.
  SetParticipationResolution(box, nullptr);
  SetParticipationCoverage(box, nullptr);
  CHECK(std::isnan(box.m90));
  CHECK(ParticipationResolutionWord(box.m90, 3, true, true, final_bars, "Electrostatic") ==
        "NOT-EVALUATED");
  CHECK(box.coverage == "NOT-MEASURED");
  CHECK(box.regime == "NO-MODEL");
}

TEST_CASE("Surface participation box rows are guarded by their own census",
          "[surfaceparticipation][Serial]")
{
  // The unknowns double; the Type's cell is refined (m90 3.0 -> 1.5) but the box's own
  // patches are not (17 R both times): the Type total is evaluated, the box value reads
  // CELL-NOT-REFINED although its p did not move (on the Type's census it would SETTLE).
  const auto bars = FinalBars();
  auto Level = [](long long int unknowns, double type_m90, double box_m90)
  {
    auto total = Row(1, 0, "SA", 1.0e-4, 1.2e-4);
    total.type_total = true;
    total.m90 = type_m90;
    auto whole = Row(1, 1, "SA", 1.0e-4, 1.2e-4);
    whole.m90 = type_m90;
    auto box = Row(1, 7, "SA", 1.0e-6, 1.6e-6, true);
    box.m90 = box_m90;
    std::vector<SurfaceParticipationRow> rows = {total, whole, box};
    for (auto &row : rows)
    {
      row.problem = "Electrostatic";
      row.unknowns = unknowns;
    }
    return rows;
  };
  auto first = Level(1000000, 3.0, 17.0);
  JoinParticipationLevel(first, nullptr, bars, 3);
  auto second = Level(2000000, 1.5, 17.0);
  JoinParticipationLevel(second, &first, bars, 3);
  CHECK(second[0].convergence == "SETTLED");
  CHECK(second[0].per_doubling == 0.0);
  CHECK(second[1].convergence == "SETTLED");  // its Type's word
  CHECK(second[2].convergence == "MOVING: CELL-NOT-REFINED");
  CHECK(std::isnan(second[2].per_doubling));
  CHECK(second[2].resolution == "UNDER-RESOLVED-UNTESTED");
  // Once the box's own patches are refined (17 -> 1.8 R), the box is evaluated and
  // RESOLVED on its own census, a scope the bars were not graded with: flagged.
  auto third = Level(4000000, 1.4, 1.8);
  JoinParticipationLevel(third, &second, bars, 3);
  CHECK(third[2].convergence == "SETTLED");
  CHECK(third[2].resolution == "RESOLVED-UNTESTED");
  CHECK(third[2].flags == "RESOLVED-UNTESTED");
  CHECK(third[0].resolution == "RESOLVED");
  CHECK(third[0].flags.empty());
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "Surface participation table file across levels and a restarted run",
                 "[surfaceparticipation][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  const auto path = temp_dir / kSurfaceParticipationFile;
  std::ofstream(path) << "a previous run's table\n";
  auto Run = [&]()
  { return SurfaceParticipationTable(path.string(), "Electrostatic", 3, InterimBars()); };
  auto table = Run();
  CHECK_FALSE(table.Started());
  table.AddRow(Row(1, 1, "SA", 1.0e-4, 1.2e-4));
  table.AddRow(Row(1, 2, "MS", 2.0e-4, 2.4e-4));
  table.FinishLevel(1, 1000, true, comm);
  CHECK(table.Started());
  std::string text = ReadText(path);
  CHECK(CountLines(text) == 1 + 4);  // a Type total per Type, then the interface rows
  CHECK(text.find("previous run") == std::string::npos);
  CHECK(text.find("convergence (INTERIM)") != std::string::npos);
  CHECK(text.find("resolution (INTERIM; RESOLVED needs SETTLED)") != std::string::npos);
  CHECK(text.find("change per doubling (%; an interface row: its Type's)") !=
        std::string::npos);
  CHECK(text.find("change corrected (fraction)") != std::string::npos);
  CHECK(text.find("change raw (fraction)") != std::string::npos);
  CHECK(text.find("h/R m90") != std::string::npos);
  CHECK(text.find("FIRST") != std::string::npos);
  for (auto row : {Row(1, 1, "SA", 1.1e-4, 1.21e-4), Row(1, 2, "MS", 2.2e-4, 2.6e-4)})
  {
    row.m90 = 0.3;
    table.AddRow(row);
  }
  table.FinishLevel(2, 4000, true, comm);
  text = ReadText(path);
  CHECK(CountLines(text) == 1 + 8);
  CHECK(text.find("SETTLED-INTERIM") != std::string::npos);
  CHECK(text.find("MOVING") != std::string::npos);
  REQUIRE(table.LastLevel().size() == 4);
  CHECK(table.LastLevel()[0].type_total);
  CHECK(table.LastLevel()[0].interface == 0);
  CHECK(table.LastLevel()[0].p_corrected == 1.21e-4);
  CHECK(table.LastLevel()[0].convergence == "SETTLED-INTERIM");
  CHECK(table.LastLevel()[2].convergence == "SETTLED-INTERIM");  // its Type's word
  CHECK(table.LastLevel()[3].convergence == "MOVING");
  CHECK(table.LastLevel()[0].iteration == 2);
  CHECK(table.LastLevel()[0].unknowns == 4000);
  CHECK(table.LastLevel()[0].problem == "Electrostatic");
  // A level whose solve did not converge appends nothing and keeps the last level.
  table.AddRow(Row(1, 1, "SA", 9.0, 9.0));
  table.FinishLevel(3, 8000, false, comm);
  CHECK(ReadText(path) == text);
  CHECK(table.LastLevel()[0].iteration == 2);
  const auto meta = table.Metadata();
  CHECK(meta["File"] == kSurfaceParticipationFile);
  CHECK(meta["Iterations"] == 2);
  CHECK(meta["Rows"] == 8);
  CHECK(meta["Verdict"]["BuiltIn"] == true);
  CHECK(meta["Verdict"]["Interim"] == true);
  CHECK(meta["Verdict"]["MinUnknownsRatio"] == 1.2);
  CHECK(meta["Verdict"]["MaxChangePerDoublingPercent"] == 1.0);
  CHECK(meta["Verdict"]["CensusBarOverR"]["3"] == 2.0);
  CHECK(meta["Verdict"]["CensusBarOverR"]["1"].is_null());
  CHECK(meta["Verdict"]["UnderResolvedFactor"] == 2.0);
  // The iteration archive moves every other file into iterationNN/ and leaves the table.
  std::ofstream(temp_dir / "surface-Q-corrected.csv") << "i\n1\n";
  SaveIteration(comm, temp_dir, 1, 1);
  CHECK(fs::is_regular_file(temp_dir / "iteration1" / "surface-Q-corrected.csv"));
  CHECK(fs::is_symlink(temp_dir / "surface-Q-corrected.csv"));
  CHECK_FALSE(fs::is_symlink(path));
  CHECK_FALSE(fs::exists(temp_dir / "iteration1" / kSurfaceParticipationFile));
  CHECK(ReadText(path) == text);
  // A new run in the same directory replaces the file at its first level.
  auto restarted = Run();
  restarted.AddRow(Row(1, 1, "SA", 1.0e-4, 1.2e-4));
  restarted.FinishLevel(1, 1000, true, comm);
  CHECK(CountLines(ReadText(path)) == 1 + 2);
  // A run that never has a row writes no file.
  fs::remove(path);
  auto empty = Run();
  empty.FinishLevel(1, 1000, true, comm);
  CHECK_FALSE(empty.Started());
  CHECK_FALSE(fs::exists(path));
}

// CUT-PLAN step 4b (lane v1step4b): a PreflightOnly run's rows carry the census PREFLIGHT
// word of their Type's all-planes row (a region row and a Type without a census row:
// NOT-EVALUATED), compute no value, and are written with their Type totals without a join.
// ONE verdict (lane v1verdict1): the word is bare, the table's own (no 'PREFLIGHT ' label).
TEST_CASE_METHOD(palace::test::SharedTempDir, "Surface participation PreflightOnly rows",
                 "[surfaceparticipation][preflightonly][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  std::vector<ResolutionPreflightRow> preflight(3);
  preflight[0].type = "SA";
  preflight[0].resolution = "UNDER-RESOLVED";
  preflight[1].type = "SA";
  preflight[1].all_planes = false;  // a plane's word is not the Type's
  preflight[1].resolution = "MARGINAL";
  preflight[2].type = "MS";
  preflight[2].resolution = "NO-BARS";
  const std::map<int, std::string> targets = {{4, "SA"}, {5, "MS"}, {7, "SA"}, {8, "MA"}};
  const auto rows = PreflightParticipationRows(targets, {7}, preflight);
  REQUIRE(rows.size() == 4);
  CHECK(rows[0].resolution == "UNDER-RESOLVED");
  CHECK(rows[1].resolution == "NO-BARS");
  CHECK(rows[2].region);
  CHECK(rows[2].resolution == "NOT-EVALUATED");
  CHECK(rows[3].resolution == "NOT-EVALUATED");
  for (const auto &row : rows)
  {
    INFO("interface " << row.interface);
    CHECK(row.point == 0);
    CHECK(row.excitation == 0);
    CHECK(std::isnan(row.frequency));
    for (const double value :
         {row.p_raw, row.p_corrected, row.q_raw, row.q_corrected, row.m90, row.m50,
          row.uncovered_share, row.change_corrected, row.change_raw, row.per_doubling})
    {
      CHECK(std::isnan(value));
    }
    for (const auto *word : {&row.coverage, &row.regime, &row.confidence, &row.convergence})
    {
      CHECK(*word == kParticipationNotComputed);
    }
    CHECK(row.flags == kParticipationPreflightFlags);
    CHECK(row.reason == kParticipationPreflightReason);
    CHECK(row.reason.find(',') == std::string::npos);  // a text cell holds no comma
  }
  const auto path = temp_dir / (std::string("preflight-") + kSurfaceParticipationFile);
  SurfaceParticipationTable table(path.string(), "Eigenmode", 2, InterimBars());
  table.WritePreflightLevel(rows, 12345, comm);
  CHECK(table.Started());
  const std::string text = ReadText(path);
  CHECK(CountLines(text) == 1 + 3 + 4);  // the SA, MS, MA totals (not the box), the rows
  CHECK(text.find("UNDER-RESOLVED") != std::string::npos);
  CHECK(text.find("PREFLIGHT ") == std::string::npos);
  CHECK(text.find("S corner+cluster share (fraction)") != std::string::npos);
  CHECK(text.find("12345") != std::string::npos);
  CHECK(text.find("MOVING") == std::string::npos);
  CHECK(text.find("FIRST-ITERATION") == std::string::npos);
  std::size_t not_computed = 0;
  for (auto at = text.find(kParticipationNotComputed); at != std::string::npos;
       at = text.find(kParticipationNotComputed, at + 1))
  {
    not_computed++;
  }
  CHECK(not_computed == 4 * 7);  // coverage, regime, confidence, convergence of each row
  const auto meta = table.Metadata();
  CHECK(meta["PreflightOnly"] == true);
  CHECK(meta["Iterations"] == 1);
  CHECK(meta["Rows"] == 7);
  // A solve-form table records no PreflightOnly member.
  SurfaceParticipationTable solved((temp_dir / kSurfaceParticipationFile).string(),
                                   "Eigenmode", 2, InterimBars());
  solved.AddRow(Row(1, 4, "SA", 1.0e-4, 1.2e-4));
  solved.FinishLevel(1, 1000, true, comm);
  CHECK_FALSE(solved.Metadata().contains("PreflightOnly"));
}

TEST_CASE("Surface participation log block", "[surfaceparticipation][Serial]")
{
  const auto bars = InterimBars();
  std::vector<SurfaceParticipationRow> rows = {Row(1, 1, "SA", 1.0e-4, 1.2e-4),
                                               Row(1, 7, "SA", 1.0e-6, 1.6e-6, true)};
  rows[0].convergence = "MOVING";
  rows[0].type_total = true;
  const auto block = ParticipationLogBlock(rows, 2, 3000, 1000, bars, 3);
  CHECK(block.find("iteration 2: 3000 unknowns (x3.000 over iteration 1), Order 3") !=
        std::string::npos);
  CHECK(block.find("Verdicts (INTERIM)") != std::string::npos);
  CHECK(block.find("built-in palace/models/surfaceparticipationverdict.json") !=
        std::string::npos);
  CHECK(block.find("a single solve never reads RESOLVED") != std::string::npos);
  CHECK(block.find("region resolution words -UNTESTED") != std::string::npos);
  CHECK(CountLines(block) == 1 + 2 + 1 + rows.size());
  // Only the targeted values count: 45 untargeted region rows (raw only) leave the block
  // in its per-value form with the same two values.
  auto with_untargeted = rows;
  for (int interface = 10; interface < 55; interface++)
  {
    auto raw_only = Row(1, interface, "SA", 1.0e-7, kNaN, true);
    SetParticipationNotTarget(raw_only);
    with_untargeted.push_back(raw_only);
  }
  CHECK(ParticipationLogBlock(with_untargeted, 2, 3000, 1000, bars, 3) == block);
  // No region value, no region note.
  CHECK(ParticipationLogBlock({rows[0]}, 2, 3000, 1000, bars, 3)
            .find("region resolution words") == std::string::npos);
  // Above kParticipationLogRows values: one line per (Type, whole | region), with the raw
  // and the corrected ranges.
  std::vector<SurfaceParticipationRow> many;
  for (int point = 1; point <= 30; point++)
  {
    many.push_back(Row(point, 0, "SA", 1.0e-4, 1.1e-4));
    many.push_back(Row(point, 0, "MS", 1.0e-4, 1.3e-4));
    many[many.size() - 2].type_total = many.back().type_total = true;
  }
  const auto grouped = ParticipationLogBlock(many, 1, 1000, 0, bars, 1);
  CHECK(grouped.find("NO-BARS") != std::string::npos);
  CHECK(grouped.find("60 values") != std::string::npos);
  CHECK(grouped.find(
            "SA      Type   30 values: p raw 1.00000e-04 .. 1.00000e-04, p corrected "
            "1.10000e-04 .. 1.10000e-04, corr/raw 1.1000 .. 1.1000") != std::string::npos);
  CHECK(CountLines(grouped) == 1 + 2 + 1 + 2);
}
