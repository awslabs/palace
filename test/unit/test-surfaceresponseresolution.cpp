// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include "models/surfaceresponseoperator.hpp"
#include "models/surfaceresponseresolution.hpp"
#include "utils/communication.hpp"

// Lane rescensus-product (decision 1817; ruling 1819): the edge-resolution census. The
// definitions (bars, patch and cell statuses, the quantile histogram, the PREFLIGHT's
// census weighted by length and its word) on hand-made sums, and the knot location on
// meshes whose element sizes are known by construction.

using namespace palace;

namespace
{

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// The longest edge of a linear element from its raw vertices (independent of the census).
double RawLongestEdge(const mfem::Mesh &mesh, int element)
{
  const mfem::Element *el = mesh.GetElement(element);
  double longest = 0.0;
  for (int edge = 0; edge < el->GetNEdges(); edge++)
  {
    const int *ends = el->GetEdgeVertices(edge);
    const double *a = mesh.GetVertex(el->GetVertices()[ends[0]]);
    const double *b = mesh.GetVertex(el->GetVertices()[ends[1]]);
    longest = std::max(longest, std::hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]));
  }
  return longest;
}

}  // namespace

TEST_CASE("Resolution census bars, patch statuses and the share rule",
          "[SurfaceResponseResolution][Serial]")
{
  const auto bars = BuiltInResolutionBars();
  CHECK(std::isnan(bars.B5(1)));
  CHECK(bars.B10(1) == 0.2657095439599277);
  CHECK(bars.B5(2) == 0.540763397421779);
  CHECK(bars.B10(2) == 1.3277577835971868);
  CHECK(bars.B10(3) == 2.597499085375387);
  CHECK(bars.B5(4) == 1.3277577835971868);
  CHECK(bars.B10(4) == 3.2419880492445246);
  CHECK(bars.B5(5) == bars.B5(4));
  CHECK(bars.B10(7) == bars.B10(4));
  // Valid at R = 1.9 um only (ERRATUM-1662 E-3): another radius reads NO-BARS.
  CHECK(ResolutionBarsHold(bars, 1.9e-6));
  CHECK(ResolutionBarsHold(bars, 1.9e-6 * (1.0 + 5.0e-7)));
  CHECK_FALSE(ResolutionBarsHold(bars, 1.9e-6 * (1.0 + 2.0e-6)));
  CHECK_FALSE(ResolutionBarsHold(bars, 2.0e-6));
  CHECK_FALSE(ResolutionBarsHold(bars, 0.0));

  using S = ResolutionStatus;
  CHECK(ResolutionPatchStatus(0.5, 2, bars, true) == S::RESOLVED);
  CHECK(ResolutionPatchStatus(bars.B5(2), 2, bars, true) == S::RESOLVED);
  CHECK(ResolutionPatchStatus(1.0, 2, bars, true) == S::MARGINAL);
  CHECK(ResolutionPatchStatus(bars.B10(2), 2, bars, true) == S::MARGINAL);
  CHECK(ResolutionPatchStatus(2.0, 2, bars, true) == S::UNDER_RESOLVED);
  CHECK(ResolutionPatchStatus(2.0, 3, bars, true) == S::MARGINAL);
  CHECK(ResolutionPatchStatus(0.1, 1, bars, true) == S::MARGINAL);  // no b5 at Order 1
  CHECK(ResolutionPatchStatus(0.3, 1, bars, true) == S::UNDER_RESOLVED);
  CHECK(ResolutionPatchStatus(1.0, 5, bars, true) == S::RESOLVED);  // Order 4's bars
  CHECK(ResolutionPatchStatus(kNaN, 2, bars, true) == S::NOT_EVALUATED);
  CHECK(ResolutionPatchStatus(0.5, 2, bars, false) == S::NO_BARS);
  CHECK(ResolutionPatchStatus(kNaN, 2, bars, false) == S::NOT_EVALUATED);

  CHECK(ResolutionShareStatus({0.9, 0.05, 0.05, 0.0}) == S::RESOLVED);
  CHECK(ResolutionShareStatus({0.89, 0.06, 0.05, 0.0}) == S::MARGINAL);
  CHECK(ResolutionShareStatus({0.2, 0.2, 0.51, 0.09}) == S::UNDER_RESOLVED);
  CHECK(ResolutionShareStatus({0.2, 0.3, 0.5, 0.0}) == S::MARGINAL);
  CHECK(ResolutionShareStatus({0.0, 0.0, 0.49, 0.51}) == S::NOT_EVALUATED);
  CHECK(ResolutionShareStatus({0.0, 0.0, 0.0, 0.5}) == S::MARGINAL);

  CHECK(ResolutionStatusRank(S::UNDER_RESOLVED) > ResolutionStatusRank(S::MARGINAL));
  CHECK(ResolutionStatusRank(S::MARGINAL) > ResolutionStatusRank(S::NO_BARS));
  CHECK(ResolutionStatusRank(S::NO_BARS) > ResolutionStatusRank(S::NOT_EVALUATED));
  CHECK(ResolutionStatusRank(S::NOT_EVALUATED) > ResolutionStatusRank(S::RESOLVED));

  CHECK(ResolutionHistogramBin(1.0) == 300);
  CHECK(ResolutionHistogramBin(1.0e-4) == 0);
  CHECK(ResolutionHistogramBin(1.0e4) == kResolutionHistogramBins - 1);
  CHECK(ResolutionHistogramCentre(300) == std::pow(10.0, -3.0 + 300.5 * 6.0 / 600.0));

  std::vector<double> odd = {3.0, 1.0, 2.0}, even = {4.0, 1.0, 3.0, 2.0}, none;
  CHECK(ResolutionLowerMedian(odd) == 2.0);
  CHECK(ResolutionLowerMedian(even) == 2.0);
  CHECK(std::isnan(ResolutionLowerMedian(none)));
}

TEST_CASE("Resolution census cells from known patch energies",
          "[SurfaceResponseResolution][Serial]")
{
  using S = ResolutionStatus;
  ResolutionAccumulator sums;
  sums.Reset(1, 2);
  sums.Add(0, 1, S::RESOLVED, 0.4, 6.0);
  sums.Add(0, 1, S::MARGINAL, 1.0, 3.0);
  sums.Add(0, 1, S::UNDER_RESOLVED, 2.0, 1.0);
  sums.Add(0, 1, S::NOT_EVALUATED, kNaN, 0.5);
  sums.Add(0, 1, S::RESOLVED, 0.4, 0.0);
  for (const int slot : {0, 1})
  {
    const auto cell =
        ComputeResolutionCell(sums.Slot(0, slot), "SA", slot == 0, 4.8, true, 2.0, "");
    CHECK(cell.patches == 5);
    CHECK(cell.not_evaluated == 1);
    CHECK(cell.energy == 21.0);
    CHECK(cell.shares[0] == 6.0 / 10.5);
    CHECK(cell.shares[1] == 3.0 / 10.5);
    CHECK(cell.shares[2] == 1.0 / 10.5);
    CHECK(cell.shares[3] == 0.5 / 10.5);
    CHECK(cell.status == S::MARGINAL);
    CHECK(cell.quantiles[0] == ResolutionHistogramCentre(ResolutionHistogramBin(0.4)));
    CHECK(cell.quantiles[1] == ResolutionHistogramCentre(ResolutionHistogramBin(0.4)));
    CHECK(cell.quantiles[2] == ResolutionHistogramCentre(ResolutionHistogramBin(1.0)));
  }
  {
    // Every patch UNDER-RESOLVED but one RESOLVED holding 95 % of the energy: RESOLVED.
    ResolutionAccumulator resolved;
    resolved.Reset(1, 1);
    resolved.Add(0, 0, S::RESOLVED, 0.3, 95.0);
    resolved.Add(0, 0, S::UNDER_RESOLVED, 4.0, 5.0);
    CHECK(
        ComputeResolutionCell(resolved.Slot(0, 0), "MS", true, 0.0, true, 1.0, "").status ==
        S::RESOLVED);
    resolved.Add(0, 0, S::NOT_EVALUATED, kNaN, 150.0);
    const auto ne =
        ComputeResolutionCell(resolved.Slot(0, 0), "MS", true, 0.0, true, 1.0, "");
    CHECK(ne.status == S::NOT_EVALUATED);
    CHECK(ne.shares[3] == 150.0 / 250.0);
  }
  {
    // NO-BARS: no R / M / U share; the not-evaluated share and the quantiles remain.
    ResolutionAccumulator no_bars;
    no_bars.Reset(1, 1);
    no_bars.Add(0, 0, S::NO_BARS, 0.4, 2.0);
    no_bars.Add(0, 0, S::NOT_EVALUATED, kNaN, 1.0);
    const auto cell =
        ComputeResolutionCell(no_bars.Slot(0, 0), "MA", true, 0.0, false, 1.0, "");
    CHECK(cell.status == S::NO_BARS);
    CHECK(std::isnan(cell.shares[0]));
    CHECK(std::isnan(cell.shares[2]));
    CHECK(cell.shares[3] == 1.0 / 3.0);
    CHECK(cell.quantiles[1] == ResolutionHistogramCentre(ResolutionHistogramBin(0.4)));
  }
  {
    // Disabled: every patch NOT-EVALUATED (counted), the energy kept.
    ResolutionAccumulator off;
    off.Reset(1, 1);
    off.Add(0, 0, S::NOT_EVALUATED, kNaN, 2.0);
    off.Add(0, 0, S::NOT_EVALUATED, kNaN, 1.0);
    const auto cell =
        ComputeResolutionCell(off.Slot(0, 0), "SA", true, 0.0, true, 1.0, "disabled");
    CHECK(cell.status == S::NOT_EVALUATED);
    CHECK(cell.reason == "disabled");
    CHECK(cell.energy == 3.0);
    CHECK(cell.shares == std::array<double, 4>{0.0, 0.0, 0.0, 1.0});
    CHECK(cell.not_evaluated == 2);
    const auto none = ComputeResolutionCell(off.Slot(0, 0), "SA", true, 0.0, true, 1.0,
                                            "no evaluation-0 energy");
    CHECK(none.status == S::NOT_EVALUATED);
    CHECK(std::isnan(none.energy));
    ResolutionAccumulator empty;
    empty.Reset(1, 1);
    empty.Add(0, 0, S::RESOLVED, 0.3, 0.0);
    const auto zero =
        ComputeResolutionCell(empty.Slot(0, 0), "SA", true, 0.0, true, 1.0, "");
    CHECK(zero.status == S::NOT_EVALUATED);
    CHECK(zero.reason == "no energy");
  }
  std::ostringstream table;
  table << ResolutionTableHeader("i");
  WriteResolutionRows(
      table, "1",
      {ComputeResolutionCell(sums.Slot(0, 0), "SA", true, 0.0, true, 1.0, ""),
       ComputeResolutionCell(sums.Slot(0, 1), "SA", false, 4.8, true, 1.0, "")});
  const std::string text = table.str();
  CHECK(text.rfind("i,Type,plane z (um),patches,patches not evaluated,energy (J),", 0) ==
        0);
  CHECK(text.find("\n1,SA,all,5,1,+1.050000000000e+01,") != std::string::npos);
  CHECK(text.find("\n1,SA,4.8,5,1,") != std::string::npos);
  CHECK(text.find(",MARGINAL,\n") != std::string::npos);
}

// ONE verdict (decision 2403 (3), lane v1verdict1, DESIGN-NOTE 1.5): the PREFLIGHT census
// weighted by length (a spatial patch as edge patches with as many knots) and its word, the
// table's word at a first iteration: known answers.
TEST_CASE("Resolution census PREFLIGHT weighs patches by length and reads the table's word",
          "[SurfaceResponseResolution][Serial]")
{
  const auto bars = BuiltInParticipationVerdictBars();  // B = 2.0 R at Order 3
  const int fine = ResolutionHistogramBin(0.5), coarse = ResolutionHistogramBin(10.0);
  // Edge patches: total length and knots, and the length at m = 0.5; spatial: knots at m.
  auto Row = [&](double edge_length, double edge_knots, double spatial_knots, int bin)
  {
    std::vector<double> slot(kResolutionPreflightStride, 0.0);
    slot[0] = edge_length;
    slot[10] = edge_length > 0.0 ? 9.0 : 0.0;
    slot[12] = edge_knots;
    slot[kResolutionPreflightEdgeHistogram + fine] = edge_length;
    slot[11] = spatial_knots > 0.0 ? 1.0 : 0.0;
    slot[7] = slot[11];
    slot[kResolutionPreflightSpatialHistogram + bin] = spatial_knots;
    return ComputeResolutionPreflightRow(slot.data(), "SA", true, 0.0, true, 1.0e-6);
  };
  auto Word = [&](const ResolutionPreflightRow &row, int order)
  {
    return ParticipationPreflightWord(row.length_quantiles[2], order, bars,
                                      "Electrostatic");
  };
  // Fine edges (9 of length 1, 10 knots each: 0.1 per knot) and ONE coarse spatial patch of
  // 20 knots (2.0 of 11.0, above 10 %): the m90 is the coarse patch's; UNDER-RESOLVED.
  const auto above = Row(9.0, 90.0, 20.0, coarse);
  CHECK(above.length_quantiles[2] == ResolutionHistogramCentre(coarse));
  CHECK(above.length_quantiles[1] == ResolutionHistogramCentre(fine));
  CHECK(Word(above, 3) == "UNDER-RESOLVED");
  CHECK(above.edge_length == Catch::Approx(9.0e-6).epsilon(1e-14));
  CHECK(above.edge_shares[0] == 1.0);  // the census bars' shares stay numbers
  CHECK(above.spatial_shares[2] == 1.0);
  // The same spatial patch at 5 knots (0.5 of 9.5, below 10 %): the edges' m90; MARGINAL,
  // never RESOLVED before a solve although m90 <= B.
  const auto below = Row(9.0, 90.0, 5.0, coarse);
  CHECK(below.length_quantiles[2] == ResolutionHistogramCentre(fine));
  CHECK(Word(below, 3) == "MARGINAL");
  CHECK(Word(below, 3) == ParticipationResolutionWord(below.length_quantiles[2], 3, false,
                                                      false, bars, "Electrostatic"));
  // A spatial-only census by its knots; no reading; the orders without a bar.
  const auto spatial_only = Row(0.0, 0.0, 20.0, ResolutionHistogramBin(3.8));
  CHECK(spatial_only.length_quantiles[2] ==
        ResolutionHistogramCentre(ResolutionHistogramBin(3.8)));
  CHECK(Word(spatial_only, 3) == "MARGINAL");
  CHECK(Word(spatial_only, 2) == "UNDER-RESOLVED");  // 3.76 > 2 x 1.7
  CHECK(Word(Row(0.0, 0.0, 0.0, 0), 3) == "NOT-EVALUATED");
  CHECK(Word(above, 1) == "NO-BARS");
  CHECK(Word(above, 5) == "NO-BARS");
  for (const int order : {1, 2, 3, 4, 5})
  {
    for (const auto &row : {above, below, spatial_only})
    {
      CHECK(Word(row, order).find("RESOLVED") != 0);  // never RESOLVED
    }
  }
  auto worded = above;
  worded.resolution = Word(above, 3);
  const std::string line = ResolutionPreflightLine({worded, below}, true, 3, 2.0, 4.0);
  CHECK(line.find("B = 2 R, UNDER-RESOLVED above 4 R") != std::string::npos);
  CHECK(line.find("SA all planes: UNDER-RESOLVED (length m90 ") != std::string::npos);
  CHECK(line.find("never RESOLVED") != std::string::npos);
  CHECK(ResolutionPreflightLine({above}, true, 1, kNaN, kNaN).find("NO-BARS") !=
        std::string::npos);
  CHECK(ResolutionPreflightLine({above}, false, 3, 2.0, 4.0).find("DISABLED") !=
        std::string::npos);
}

// The PREFLIGHT's quantile is the energy cell's (one function).
TEST_CASE("Resolution census quantiles are one function",
          "[SurfaceResponseResolution][Serial]")
{
  ResolutionAccumulator accumulator;
  accumulator.Reset(1, 1);
  accumulator.Add(0, 0, ResolutionStatus::MARGINAL, 0.7, 2.0);
  accumulator.Add(0, 0, ResolutionStatus::UNDER_RESOLVED, 9.0, 1.0);
  const double *slot = accumulator.Slot(0, 0);
  const auto cell = ComputeResolutionCell(slot, "SA", true, 0.0, true, 1.0, "");
  CHECK(cell.quantiles ==
        ResolutionHistogramQuantiles(slot + ResolutionAccumulator::kStatusBins));
  CHECK(cell.quantiles[2] == ResolutionHistogramCentre(ResolutionHistogramBin(9.0)));
}

TEST_CASE("Resolution bars file: read and checked against the run radius",
          "[SurfaceResponseResolution][Serial]")
{
  const auto dir = std::filesystem::temp_directory_path() /
                   ("palace-resolution-bars-" + std::to_string(Mpi::Rank(Mpi::World())));
  std::filesystem::create_directories(dir);
  const auto good = dir / "bars-r2.json";
  std::ofstream(good) << R"({"RadiusMicrons": 2.0, "Provenance": "unit test", "Bars": {
    "1": {"b5": null, "b10": 0.3}, "2": {"b5": 0.5, "b10": 1.2}, "3": {"b5": 0.6, "b10": 2.0},
    "4": {"b5": 1.0, "b10": 3.0}}})";
  const auto bars = ReadResolutionBars(good.string());
  CHECK(bars.radius == Catch::Approx(2.0e-6).epsilon(1e-15));
  CHECK(std::isnan(bars.B5(1)));
  CHECK(bars.B10(2) == 1.2);
  CHECK_FALSE(bars.built_in);
  CHECK(bars.sha256.size() == 64);
  CHECK(ResolutionBarsHold(bars, 2.0e-6));
  CHECK_FALSE(ResolutionBarsHold(bars, 1.9e-6));
  // The status follows the file's bars (a bars file switches it).
  CHECK(ResolutionPatchStatus(1.0, 2, bars, true) == ResolutionStatus::MARGINAL);
  CHECK(ResolutionPatchStatus(1.25, 2, bars, true) == ResolutionStatus::UNDER_RESOLVED);
  const auto bad = dir / "bars-missing-order.json";
  std::ofstream(bad)
      << R"({"RadiusMicrons": 2.0, "Bars": {"1": {"b5": null, "b10": 0.3}}})";
  CHECK_THROWS(ReadResolutionBars(bad.string()));
  std::filesystem::remove_all(dir);
}

TEST_CASE("Resolution census locates knots on meshes of known element sizes",
          "[SurfaceResponseResolution][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("The census's point locator requires MFEM_USE_GSLIB");
#else
  MPI_Comm comm = Mpi::World();
  const auto bars = BuiltInResolutionBars();
  {
    // Linear tetrahedra of unit boxes (6 per box, every one with the box diagonal as its
    // longest edge): h = sqrt(3) wherever a knot lies.
    mfem::Mesh serial =
        mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::TETRAHEDRON, 4.0, 4.0, 4.0);
    double smallest = mfem::infinity(), largest = 0.0;
    for (int e = 0; e < serial.GetNE(); e++)
    {
      smallest = std::min(smallest, RawLongestEdge(serial, e));
      largest = std::max(largest, RawLongestEdge(serial, e));
    }
    REQUIRE(smallest == Catch::Approx(std::sqrt(3.0)).epsilon(1e-14));
    REQUIRE(largest == Catch::Approx(std::sqrt(3.0)).epsilon(1e-14));
    mfem::ParMesh mesh(comm, serial);
    const std::vector<std::vector<std::array<double, 3>>> knots = {
        {{0.3, 0.4, 0.6},
         {1.7, 2.2, 3.1},
         {3.9, 0.1, 2.5},
         {2.4, 3.3, 0.7},
         {0.5, 1.5, 2.5}},
        {{10.0, 0.5, 0.5}, {-1.0, 0.5, 0.5}, {0.5, 9.0, 0.5}},  // outside: NOT-EVALUATED
        {{1.1, 1.2, 1.3}, {2.6, 2.7, 2.8}, {9.0, 9.0, 9.0}, {-3.0, 0.0, 0.0}},  // 2 of 4
        {{1.1, 1.2, 1.3},
         {2.6, 2.7, 2.8},
         {9.0, 9.0, 9.0},
         {-3.0, 0.0, 0.0},
         {7.0, 7.0, 7.0}},  // 2 of 5: fewer than half located
        {}};
    const auto readings = LocateResolutionPatchKnots(mesh, 3, knots);
    REQUIRE(readings.size() == 5);
    CHECK(readings[0].knots == 5);
    CHECK(readings[0].located == 5);
    CHECK(readings[0].h == Catch::Approx(std::sqrt(3.0)).epsilon(1e-14));
    CHECK(readings[1].located == 0);
    CHECK_FALSE(readings[1].Evaluated());
    CHECK(readings[2].located == 2);
    CHECK(readings[2].Evaluated());
    CHECK(readings[2].h == Catch::Approx(std::sqrt(3.0)).epsilon(1e-14));
    CHECK(readings[3].located == 2);
    CHECK_FALSE(readings[3].Evaluated());
    CHECK(readings[4].knots == 0);
    CHECK_FALSE(readings[4].Evaluated());
    // m = h / R: R = sqrt(3) / 0.5 -> m 0.5: RESOLVED at Orders 2-4, UNDER at Order 1.
    const double radius = std::sqrt(3.0) / 0.5;
    CHECK(ResolutionPatchStatus(readings[0].h / radius, 2, bars, true) ==
          ResolutionStatus::RESOLVED);
    CHECK(ResolutionPatchStatus(readings[0].h / radius, 1, bars, true) ==
          ResolutionStatus::UNDER_RESOLVED);
  }
  {
    // Hexahedra of 1 x 2 x 0.5 boxes, the boxes with x < 2 refined NONCONFORMINGLY (leaves
    // 0.5 x 1 x 0.25), with a Nodes grid function as Palace's meshes: the longest EDGE (not
    // a diagonal) of the LEAF owning the knot: 1 where refined, 2 elsewhere.
    mfem::Mesh serial =
        mfem::Mesh::MakeCartesian3D(4, 2, 4, mfem::Element::HEXAHEDRON, 4.0, 4.0, 2.0);
    mfem::Array<int> refine;
    for (int e = 0; e < serial.GetNE(); e++)
    {
      mfem::Vector center(3);
      serial.GetElementCenter(e, center);
      if (center(0) < 2.0)
      {
        refine.Append(e);
      }
    }
    serial.GeneralRefinement(refine, 1);
    serial.EnsureNodes();
    mfem::ParMesh mesh(comm, serial);
    const std::vector<std::vector<std::array<double, 3>>> knots = {
        {{0.3, 0.7, 0.2}, {1.1, 3.3, 1.6}, {1.9, 2.1, 0.9}},   // refined: h = 1
        {{2.3, 0.7, 0.2}, {3.1, 3.3, 1.6}, {3.9, 2.1, 0.9}},   // coarse: h = 2
        {{0.3, 0.7, 0.2}, {1.1, 3.3, 1.6}, {3.1, 3.3, 1.6}}};  // lower median: 1
    const auto readings = LocateResolutionPatchKnots(mesh, 3, knots);
    REQUIRE(readings.size() == 3);
    CHECK(readings[0].located == 3);
    CHECK(readings[0].h == Catch::Approx(1.0).epsilon(1e-14));
    CHECK(readings[1].located == 3);
    CHECK(readings[1].h == Catch::Approx(2.0).epsilon(1e-14));
    CHECK(readings[2].h == Catch::Approx(1.0).epsilon(1e-14));
    // R = 1: m 1 / 2: Order 2 MARGINAL / UNDER-RESOLVED; Order 4 RESOLVED / MARGINAL.
    CHECK(ResolutionPatchStatus(readings[0].h, 2, bars, true) ==
          ResolutionStatus::MARGINAL);
    CHECK(ResolutionPatchStatus(readings[1].h, 2, bars, true) ==
          ResolutionStatus::UNDER_RESOLVED);
    CHECK(ResolutionPatchStatus(readings[0].h, 4, bars, true) ==
          ResolutionStatus::RESOLVED);
    CHECK(ResolutionPatchStatus(readings[1].h, 4, bars, true) ==
          ResolutionStatus::MARGINAL);
    // Decision 1863 (the owner rule of a knot on a face shared by several elements): knots
    // ON the face x = 2 between the refined leaves (h = 1) and the coarse boxes (h = 2)
    // read the LARGEST container, 2, at every rank count (the two sides may lie on
    // different ranks); knots on the face x = 1 between two refined leaves read 1.
    const std::vector<std::vector<std::array<double, 3>>> face_knots = {
        {{2.0, 0.7, 0.2}, {2.0, 3.3, 1.6}, {2.0, 2.1, 0.9}},
        {{1.0, 0.7, 0.2}, {1.0, 3.3, 1.6}, {1.0, 2.1, 0.9}}};
    const auto face_readings = LocateResolutionPatchKnots(mesh, 3, face_knots);
    REQUIRE(face_readings.size() == 2);
    CHECK(face_readings[0].located == 3);
    CHECK(face_readings[0].h == Catch::Approx(2.0).epsilon(1e-14));
    CHECK(face_readings[1].located == 3);
    CHECK(face_readings[1].h == Catch::Approx(1.0).epsilon(1e-14));
  }
#endif
}

TEST_CASE("Resolution census palace.json block and banner",
          "[SurfaceResponseResolution][Serial]")
{
  using S = ResolutionStatus;
  ResolutionRunRecord record;
  record.order = 2;
  record.radius = 2.0e-6;
  record.bars = BuiltInResolutionBars();
  record.bars_hold = ResolutionBarsHold(record.bars, record.radius);
  ResolutionCell under, no_bars;
  under.type = "SA";
  under.status = S::UNDER_RESOLVED;
  under.quantiles = {2.0, 3.0, 4.0};
  no_bars.type = "MS";
  no_bars.status = S::NO_BARS;
  no_bars.quantiles = {kNaN, kNaN, kNaN};
  record.entries = {{{{"i", 1.0}}, under}, {{{"i", 1.0}}, no_bars}};
  const auto metadata = ResolutionRunMetadata(record, 2, 3);
  CHECK_FALSE(metadata["Bars"]["Hold"].get<bool>());
  CHECK(metadata["Level"] == 2);
  CHECK(metadata["Levels"] == 3);
  CHECK(metadata["Complete"].get<bool>());
  CHECK_FALSE(metadata.contains("Final"));
  CHECK(metadata["WorstByType"]["SA"]["Worst"] == "UNDER-RESOLVED");
  CHECK(metadata["WorstByType"]["SA"]["Flagged"] == 1);
  CHECK(metadata["WorstByType"]["MS"]["Worst"] == "NO-BARS");
  CHECK(metadata["WorstByType"]["MS"]["m50Min"].is_null());
  CHECK(metadata["FlaggedCount"] == 2);
  const std::string banner = ResolutionRunBanner(metadata);
  CHECK(banner.find("SA UNDER-RESOLVED over every cell and plane (m50 3 .. 3;") !=
        std::string::npos);
  CHECK(banner.find("NOT valid for this radius") != std::string::npos);
  CHECK(banner.find("INTERIM in-sample built-in bars, valid for R = 1.9 um ONLY") !=
        std::string::npos);
  CHECK(banner.find("RESOLVED is necessary, not sufficient.") != std::string::npos);
}

TEST_CASE("Resolution census summary reads every plane's cell",
          "[SurfaceResponseResolution][Serial]")
{
  // Review r1817 MJ-2: 85 % of a Type's energy RESOLVED on one plane, 15 % on another plane
  // at 1 / 3 RESOLVED and 2 / 3 UNDER-RESOLVED: the all-planes cell reads RESOLVED (s_R =
  // 0.90), the second plane UNDER-RESOLVED; the summary must read UNDER-RESOLVED and count
  // the plane's cell.
  using S = ResolutionStatus;
  auto Summary = [](S second_status)
  {
    ResolutionAccumulator sums;
    sums.Reset(1, 3);
    sums.Add(0, 1, S::RESOLVED, 0.3, 85.0);
    sums.Add(0, 2, S::RESOLVED, 0.3, 5.0);
    sums.Add(0, 2, second_status, second_status == S::NOT_EVALUATED ? kNaN : 2.0, 10.0);
    ResolutionRunRecord record;
    record.order = 2;
    record.radius = kResolutionBuiltInRadius;
    record.bars = BuiltInResolutionBars();
    record.bars_hold = ResolutionBarsHold(record.bars, record.radius);
    std::vector<ResolutionCell> cells;
    for (int slot = 0; slot < 3; slot++)
    {
      cells.push_back(ComputeResolutionCell(sums.Slot(0, slot), "SA", slot == 0,
                                            slot == 2 ? 4.8 : 0.0, true, 1.0, ""));
      record.entries.push_back({{{"i", 1.0}}, cells.back()});
    }
    return std::make_pair(cells, ResolutionRunMetadata(record, 1, 1));
  };
  {
    const auto [cells, metadata] = Summary(S::UNDER_RESOLVED);
    CHECK(cells[0].status == S::RESOLVED);
    CHECK(cells[0].shares[0] == Catch::Approx(0.90));
    CHECK(cells[1].status == S::RESOLVED);
    CHECK(cells[2].status == S::UNDER_RESOLVED);
    CHECK(cells[2].shares[2] == Catch::Approx(2.0 / 3.0));
    CHECK(metadata["WorstByType"]["SA"]["Worst"] == "UNDER-RESOLVED");
    CHECK(metadata["WorstByType"]["SA"]["Flagged"] == 1);
    CHECK(metadata["WorstByType"]["SA"]["Cells"] == 3);
    CHECK(metadata["WorstByType"]["SA"]["ShareUnderResolvedMax"].get<double>() ==
          Catch::Approx(2.0 / 3.0));
    CHECK(metadata["FlaggedCount"] == 1);
    REQUIRE(metadata["Flagged"].size() == 1);
    CHECK(metadata["Flagged"][0]["Plane"].get<double>() == 4.8);
    CHECK(metadata["Flagged"][0]["Status"] == "UNDER-RESOLVED");
    // The m50 range over every cell: the second plane's m50 lies in the UNDER bin.
    CHECK(metadata["WorstByType"]["SA"]["m50Max"].get<double>() ==
          ResolutionHistogramCentre(ResolutionHistogramBin(2.0)));
    const std::string banner = ResolutionRunBanner(metadata);
    CHECK(banner.find("SA UNDER-RESOLVED over every cell and plane") != std::string::npos);
    CHECK(banner.find("1 of 3 cell(s) not RESOLVED") != std::string::npos);
    CHECK(banner.find("SA RESOLVED") == std::string::npos);
  }
  {
    // The same with the 2 / 3 NOT-EVALUATED: the plane reads NOT-EVALUATED, the Type is
    // never RESOLVED.
    const auto [cells, metadata] = Summary(S::NOT_EVALUATED);
    CHECK(cells[0].status == S::RESOLVED);
    CHECK(cells[2].status == S::NOT_EVALUATED);
    CHECK(metadata["WorstByType"]["SA"]["Worst"] == "NOT-EVALUATED");
    CHECK(metadata["FlaggedCount"] == 1);
    CHECK(ResolutionRunBanner(metadata).find("1 of 3 cell(s) not RESOLVED") !=
          std::string::npos);
  }
}

TEST_CASE("Resolution census summary counts a single plane's repeated cell once",
          "[SurfaceResponseResolution][Serial]")
{
  // Review r1846 mn-A: with ONE plane in the run, each plane cell repeats the all-planes
  // cell of its measurement and Type; the summary counts distinct readings (the table keeps
  // both rows). A plane cell whose reading differs from the all-planes cell (here: a patch
  // without a plane slot enters the all-planes cell only) is a distinct reading, counted.
  using S = ResolutionStatus;
  auto Summary = [](bool patch_without_plane)
  {
    ResolutionRunRecord record;
    record.order = 2;
    record.radius = kResolutionBuiltInRadius;
    record.bars = BuiltInResolutionBars();
    record.bars_hold = ResolutionBarsHold(record.bars, record.radius);
    for (const double mode : {1.0, 2.0})
    {
      ResolutionAccumulator sums;
      sums.Reset(2, 2);
      for (int type = 0; type < 2; type++)
      {
        sums.Add(type, 1, S::UNDER_RESOLVED, 2.0, 7.0);
        sums.Add(type, 1, type == 0 ? S::UNDER_RESOLVED : S::RESOLVED, 0.3, 3.0);
        if (patch_without_plane && type == 0)
        {
          sums.Add(type, 0, S::RESOLVED, 0.3, 1.0);
        }
      }
      for (int type = 0; type < 2; type++)
      {
        for (int slot = 0; slot < 2; slot++)
        {
          record.entries.push_back(
              {{{"m", mode}},
               ComputeResolutionCell(sums.Slot(type, slot), type == 0 ? "SA" : "MS",
                                     slot == 0, 0.0, true, 1.0, "")});
        }
      }
    }
    return std::make_pair(record, ResolutionRunMetadata(record, 1, 1));
  };
  {
    const auto [record, metadata] = Summary(false);
    REQUIRE(record.entries.size() == 8);
    const auto repeated = ResolutionRepeatedPlaneCells(record.entries);
    CHECK(std::count(repeated.begin(), repeated.end(), true) == 4);
    CHECK(metadata["Cells"] == 8);
    CHECK(metadata["RepeatedPlaneCells"] == 4);
    CHECK(metadata["WorstByType"]["SA"]["Cells"] == 2);
    CHECK(metadata["WorstByType"]["SA"]["Flagged"] == 2);
    CHECK(metadata["WorstByType"]["SA"]["Worst"] == "UNDER-RESOLVED");
    CHECK(metadata["WorstByType"]["MS"]["Cells"] == 2);
    CHECK(metadata["WorstByType"]["MS"]["Flagged"] == 2);
    CHECK(metadata["WorstByType"]["MS"]["Worst"] == "UNDER-RESOLVED");
    CHECK(metadata["FlaggedCount"] == 4);
    CHECK(metadata["Flagged"].size() == 4);
    for (const auto &entry : metadata["Flagged"])
    {
      CHECK(entry["Plane"].is_null());
    }
    CHECK(metadata["StatusCounts"]["UNDER-RESOLVED"] == 4);
    const std::string banner = ResolutionRunBanner(metadata);
    CHECK(banner.find("SA UNDER-RESOLVED over every cell and plane") != std::string::npos);
    CHECK(banner.find("2 of 2 cell(s) not RESOLVED") != std::string::npos);
    CHECK(banner.find("4 of 4 cell(s)") == std::string::npos);
  }
  {
    // SA's all-planes cell holds a patch the plane cell does not: SA's plane cells are
    // distinct readings (counted); MS's still repeat.
    const auto [record, metadata] = Summary(true);
    CHECK(metadata["RepeatedPlaneCells"] == 2);
    CHECK(metadata["WorstByType"]["SA"]["Cells"] == 4);
    CHECK(metadata["WorstByType"]["SA"]["Flagged"] == 4);
    CHECK(metadata["WorstByType"]["MS"]["Cells"] == 2);
    CHECK(metadata["FlaggedCount"] == 6);
  }
  {
    // Two planes: nothing repeats (each plane is its own reading).
    ResolutionRunRecord record;
    ResolutionCell all, plane_a, plane_b;
    all.type = plane_a.type = plane_b.type = "SA";
    all.status = plane_a.status = plane_b.status = S::UNDER_RESOLVED;
    plane_a.all_planes = plane_b.all_planes = false;
    plane_b.plane = 4.8;
    record.entries = {
        {{{"m", 1.0}}, all}, {{{"m", 1.0}}, plane_a}, {{{"m", 1.0}}, plane_b}};
    const auto metadata = ResolutionRunMetadata(record, 1, 1);
    CHECK(metadata["RepeatedPlaneCells"] == 0);
    CHECK(metadata["WorstByType"]["SA"]["Cells"] == 3);
  }
}

TEST_CASE("Resolution census Flagged list is capped, FlaggedCount is not",
          "[SurfaceResponseResolution][Serial]")
{
  // Review r1846 mn-B: Flagged keeps the first kResolutionFlaggedListCap flagged cells;
  // FlaggedCount and the per-Type counts count every one.
  ResolutionRunRecord record;
  ResolutionCell cell;
  cell.type = "SA";
  cell.status = ResolutionStatus::UNDER_RESOLVED;
  const int count = kResolutionFlaggedListCap + 50;
  for (int mode = 1; mode <= count; mode++)
  {
    record.entries.push_back({{{"m", static_cast<double>(mode)}}, cell});
  }
  const auto metadata = ResolutionRunMetadata(record, 1, 1);
  CHECK(kResolutionFlaggedListCap == 200);
  CHECK(metadata["Flagged"].size() == static_cast<std::size_t>(kResolutionFlaggedListCap));
  CHECK(metadata["Flagged"].back()["Keys"]["m"].get<double>() ==
        static_cast<double>(kResolutionFlaggedListCap));
  CHECK(metadata["FlaggedCount"] == count);
  CHECK(metadata["WorstByType"]["SA"]["Flagged"] == count);
  CHECK(metadata["WorstByType"]["SA"]["Cells"] == count);
}

TEST_CASE("Resolution census block without a cell names its reason",
          "[SurfaceResponseResolution][Serial]")
{
  // Review r1817 mn-1 (the coverage block's rule): a raw fallback with no applied patch.
  ResolutionRunRecord record;
  record.order = 1;
  record.radius = kResolutionBuiltInRadius;
  record.bars = BuiltInResolutionBars();
  record.bars_hold = true;
  {
    const auto metadata = ResolutionRunMetadata(record, 1, 1);
    CHECK(metadata["Cells"] == 0);
    CHECK(metadata["NoCellsReason"].is_null());
    CHECK_FALSE(metadata["Complete"].get<bool>());
    CHECK(ResolutionRunBanner(metadata).find("no census cell was measured in this run:") !=
          std::string::npos);
  }
  record.no_cells_reason = "RawFallback: NoMatchedFeature: no correction applied";
  const auto metadata = ResolutionRunMetadata(record, 1, 1);
  CHECK(metadata["NoCellsReason"] == record.no_cells_reason);
  CHECK(metadata["Complete"].get<bool>());
  CHECK(ResolutionRunBanner(metadata).find("no census cell was measured in this run "
                                           "(RawFallback: NoMatchedFeature: no "
                                           "correction applied)") != std::string::npos);
}
