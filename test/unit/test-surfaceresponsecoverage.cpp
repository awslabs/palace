// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <fmt/format.h>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "models/surfaceresponsecoverage.hpp"
#include "models/surfaceresponseidentification.hpp"

// Lane L1 (decision 1634; rulings 1641 / 1649 / 1654 (III)): the coverage cells on the
// closed partition. Every energy is a dyadic rational, so every identity holds exactly.

namespace palace
{

namespace
{

using Partition = SurfacePostOperator::InterfaceCoveragePartition;

// Interface 1: W = 16. Kept raw: Missing m3 (Missing key) 1, m5 (Unboxable) 0.5, the joint
// bin m3+m5 0.25, DomainBoundary 0.25, ConductorConsistency 1; served: m0 (qualified) 8,
// m1 (pending) 2, m0+m2 (m2 no status) 1, m2 0.5; UnpatchedPairs 0.25; unclaimed 1.25.
Partition InterfaceOne()
{
  Partition partition;
  partition.within = 16.0;
  partition.classes = {{{"Missing", "m3"}, 1.0},
                       {{"Missing", "m5"}, 0.5},
                       {{"Missing", "m3+m5"}, 0.25},
                       {{"DomainBoundary", "db"}, 0.25},
                       {{"ConductorConsistency", "cc"}, 1.0},
                       {{"UnpatchedPairs", "unpatched"}, 0.25}};
  partition.served = {{"m0", 8.0}, {"m1", 2.0}, {"m0+m2", 1.0}, {"m2", 0.5}};
  partition.unclaimed = 1.25;
  return partition;
}

std::map<int, CoverageKeyRecord> Keys()
{
  return {{3, ClassifyCoverageKey("hashA", false, "", false)},
          {5, ClassifyCoverageKey("hashB", false, "", true)},
          {7, ClassifyCoverageKey("hashA", true, "ConductorRefused", false)}};
}

std::map<int, CoverageModelQualification> Models()
{
  return {{0, CoverageModelQualification::QUALIFIED},
          {1, CoverageModelQualification::PENDING},
          {2, CoverageModelQualification::NO_STATUS}};
}

}  // namespace

TEST_CASE("Coverage cell shares and the sum identities",
          "[SurfaceResponseCoverage][Serial]")
{
  const auto partition = InterfaceOne();
  double bins = partition.unclaimed;
  for (const auto &[key, energy] : partition.classes)
  {
    bins += energy;
  }
  for (const auto &[label, energy] : partition.served)
  {
    bins += energy;
  }
  REQUIRE(bins == partition.within);  // the closed partition (the fixture's own check)

  const auto summary = SummarizeCoverageInterface(partition, 16.0, 13.5, Keys(), Models());
  CHECK(summary.served == 11.5);
  // Missing = the bins with a Missing-class token (m3, m3+m5); m5 alone is method gap.
  CHECK(summary.missing == 1.25);
  CHECK(summary.pending == 2.0);    // m1 whole
  CHECK(summary.no_status == 1.5);  // m0+m2 whole + m2
  CHECK(summary.partition_within == 16.0);
  // Review m-5: without qualification records (an explicit-model configuration) every
  // served label is no-status, never silently qualified.
  const auto unrecorded = SummarizeCoverageInterface(partition, 16.0, 13.5, Keys(), {});
  CHECK(unrecorded.served == 11.5);
  CHECK(unrecorded.no_status == 11.5);
  CHECK(unrecorded.pending == 0.0);

  const auto cell =
      ComputeCoverageCell("SA", {summary}, CoverageThresholds{}, CoverageValidatedRegime{},
                          CoverageExcitationClass::ELECTROSTATIC_TERMINAL);
  CHECK(cell.within == 16.0);
  CHECK(cell.uncovered == 4.5);  // fail-closed: W - served
  CHECK(cell.missing + cell.method_gap == cell.uncovered);
  CHECK(cell.method_gap == 3.25);
  CHECK(cell.method_gap_signed_share == 3.25 / 16.0);  // positive: the clamp is inactive
  CHECK(cell.uncovered_share == 4.5 / 16.0);
  CHECK(cell.pending_share == 2.0 / 16.0);
  CHECK(cell.no_status_share == 1.5 / 16.0);
  // The method-gap share alone (0.203) exceeds MaxUncoveredShare.
  CHECK(cell.status == CoverageStatus::NOT_ACHIEVABLE);
  CHECK(cell.factor == 13.5 / 11.5);
  CHECK(cell.regime == RegimeStatus::NO_RANGE);
  CHECK(std::isnan(cell.validated_min));
  CHECK(std::isnan(cell.validated_max));
  // Decision 1708 (1): served + missing above W (round-off of two integrals, or an overlap
  // term) clamps the method gap at 0 and the signed share keeps the negative value.
  auto over = summary;
  over.served = 15.0;
  over.missing = 1.5;
  const auto clamped =
      ComputeCoverageCell("SA", {over}, CoverageThresholds{}, CoverageValidatedRegime{},
                          CoverageExcitationClass::ELECTROSTATIC_TERMINAL);
  CHECK(clamped.method_gap == 0.0);
  CHECK(clamped.method_gap_share == 0.0);
  CHECK(clamped.method_gap_signed_share == -0.5 / 16.0);
  // Ruling 1714 (B)(3): below -1e-9 W the cell is INCONSISTENT (never OK), U as computed;
  // a round-off excess (-1e-10 W) keeps the band.
  CHECK(clamped.status == CoverageStatus::INCONSISTENT);
  CHECK(clamped.uncovered == 16.0 - 15.0);
  CHECK(CoverageStatusName(clamped.status) == "INCONSISTENT");
  CHECK(CoverageStatusRank(CoverageStatus::INCONSISTENT) >
        CoverageStatusRank(CoverageStatus::NOT_MEASURED));
  auto roundoff = summary;
  roundoff.served = 16.0 - 1.0e-10 * 16.0 - 1.0;
  roundoff.missing = 1.0 + 2.0e-10 * 16.0;
  const auto kept_band =
      ComputeCoverageCell("SA", {roundoff}, CoverageThresholds{}, CoverageValidatedRegime{},
                          CoverageExcitationClass::ELECTROSTATIC_TERMINAL);
  CHECK(kept_band.method_gap_signed_share < 0.0);
  CHECK(kept_band.status != CoverageStatus::INCONSISTENT);
}

TEST_CASE("Coverage claim-end match by the round-off rule",
          "[SurfaceResponseCoverage][Serial]")
{
  // Ruling 1715: 16 ulps per coordinate of max(|v_d|, R); near = within 1e-5 R beyond it.
  using Match = SurfacePostOperator::ClaimEndMatch;
  const double R = 0.04;
  const std::array<double, 3> v = {0.9, 0.1, 0.24};
  CHECK(SurfacePostOperator::ClassifyClaimEnd(v, v, R).first == Match::EXACT);
  auto p = v;
  p[0] = std::nextafter(v[0], 0.0);  // the restart case: 1.1e-16 short of x = 0.9
  CHECK(SurfacePostOperator::ClassifyClaimEnd(p, v, R).first == Match::ROUNDOFF);
  p[0] = v[0] - 1000.0 * std::numeric_limits<double>::epsilon() * v[0];
  const auto [near, distance] = SurfacePostOperator::ClassifyClaimEnd(p, v, R);
  CHECK(near == Match::NEAR);
  CHECK(distance > 0.0);
  p[0] = v[0] - 2.0e-5 * R;
  CHECK(SurfacePostOperator::ClassifyClaimEnd(p, v, R).first == Match::BEYOND);
  // A zero coordinate: the bar is 16 eps R there.
  const std::array<double, 3> zero = {0.0, 0.0, 0.0};
  auto q = zero;
  q[1] = 8.0 * std::numeric_limits<double>::epsilon() * R;
  CHECK(SurfacePostOperator::ClassifyClaimEnd(q, zero, R).first == Match::ROUNDOFF);
  q[1] = 32.0 * std::numeric_limits<double>::epsilon() * R;
  CHECK(SurfacePostOperator::ClassifyClaimEnd(q, zero, R).first == Match::NEAR);
}

TEST_CASE("Coverage status bands on U + pending", "[SurfaceResponseCoverage][Serial]")
{
  // W = 1024, served = 1024 - u, nothing Missing: U = u = method gap; pending p.
  auto Cell = [](double u, double missing, double pending)
  {
    CoverageInterfaceSummary summary;
    summary.within = summary.partition_within = 1024.0;
    summary.served = 1024.0 - u;
    summary.missing = missing;
    summary.pending = pending;
    summary.fabricated = summary.served;
    return ComputeCoverageCell("MS", {summary}, CoverageThresholds{0.02, 0.05},
                               CoverageValidatedRegime{},
                               CoverageExcitationClass::EIGENMODE);
  };
  CHECK(Cell(0.0, 0.0, 0.0).status == CoverageStatus::OK);
  CHECK(Cell(20.0, 20.0, 0.0).status == CoverageStatus::OK);   // 0.0195
  CHECK(Cell(20.0, 20.0, 2.0).status == CoverageStatus::LOW);  // 0.0215
  CHECK(Cell(50.0, 50.0, 1.0).status == CoverageStatus::LOW);  // 0.0498
  CHECK(Cell(50.0, 50.0, 2.0).status == CoverageStatus::PARTIAL);
  CHECK(Cell(100.0, 100.0, 0.0).status == CoverageStatus::PARTIAL);        // all Missing
  CHECK(Cell(100.0, 40.0, 0.0).status == CoverageStatus::NOT_ACHIEVABLE);  // gap 0.0586
  CHECK(Cell(100.0, 60.0, 0.0).status == CoverageStatus::PARTIAL);         // gap 0.039
  // Every eigen / driven cell reads NO-RANGE without a validated range of its class.
  CHECK(Cell(0.0, 0.0, 0.0).regime == RegimeStatus::NO_RANGE);
  CHECK(Cell(0.0, 0.0, 0.0).factor == 1.0);
  // The worse-status order of the matrix entries and the confidence row.
  CHECK(CoverageStatusRank(CoverageStatus::OK) < CoverageStatusRank(CoverageStatus::LOW));
  CHECK(CoverageStatusRank(CoverageStatus::LOW) <
        CoverageStatusRank(CoverageStatus::PARTIAL));
  CHECK(CoverageStatusRank(CoverageStatus::PARTIAL) <
        CoverageStatusRank(CoverageStatus::NOT_ACHIEVABLE));
  CHECK(CoverageStatusRank(CoverageStatus::NOT_ACHIEVABLE) <
        CoverageStatusRank(CoverageStatus::NO_ENERGY));
  CHECK(CoverageStatusRank(CoverageStatus::NO_ENERGY) <
        CoverageStatusRank(CoverageStatus::NOT_MEASURED));
  CHECK(RegimeStatusRank(RegimeStatus::NO_RANGE) < RegimeStatusRank(RegimeStatus::OUT));
}

TEST_CASE("Coverage cell without energy or without a served model",
          "[SurfaceResponseCoverage][Serial]")
{
  // Ruling 1654 (III) Q1: W = 0 is NO-ENERGY, every share nan (never 0), the factor nan
  // with regime NO-MODEL, the validated range nan.
  CoverageInterfaceSummary empty;
  const auto none =
      ComputeCoverageCell("MA", {empty}, CoverageThresholds{}, CoverageValidatedRegime{},
                          CoverageExcitationClass::ELECTROSTATIC_TERMINAL);
  CHECK(none.status == CoverageStatus::NO_ENERGY);
  for (const double share :
       {none.uncovered_share, none.missing_share, none.method_gap_share, none.pending_share,
        none.no_status_share, none.method_gap_signed_share, none.factor, none.validated_min,
        none.validated_max})
  {
    CHECK(std::isnan(share));
  }
  CHECK(none.regime == RegimeStatus::NO_MODEL);
  // Nothing served (W > 0): PARTIAL with U = 1, factor nan, NO-MODEL.
  CoverageInterfaceSummary raw;
  raw.within = raw.partition_within = 2.0;
  const auto unserved =
      ComputeCoverageCell("MA", {raw}, CoverageThresholds{}, CoverageValidatedRegime{},
                          CoverageExcitationClass::ELECTROSTATIC_TERMINAL);
  CHECK(unserved.uncovered_share == 1.0);
  CHECK(unserved.status == CoverageStatus::NOT_ACHIEVABLE);  // all of it method gap
  CHECK(unserved.regime == RegimeStatus::NO_MODEL);
  CHECK(std::isnan(unserved.factor));
  const auto unmeasured = NotMeasuredCoverageCell("SA");
  CHECK(unmeasured.status == CoverageStatus::NOT_MEASURED);
  CHECK(std::isnan(unmeasured.uncovered_share));
  CHECK(CoverageStatusName(unmeasured.status) == "NOT-MEASURED");
}

TEST_CASE("Coverage regime status against a validated range",
          "[SurfaceResponseCoverage][Serial]")
{
  CoverageValidatedRegime regime;
  regime.tolerance = 0.25;
  regime.ranges[{"SA", CoverageExcitationClass::ELECTROSTATIC_TERMINAL}] = {0.875, 1.125,
                                                                            4};
  auto Cell = [&](double fabricated, CoverageExcitationClass excitation)
  {
    CoverageInterfaceSummary summary;
    summary.within = summary.partition_within = summary.served = 8.0;
    summary.fabricated = fabricated;
    return ComputeCoverageCell("SA", {summary}, CoverageThresholds{}, regime, excitation);
  };
  const auto electrostatic = CoverageExcitationClass::ELECTROSTATIC_TERMINAL;
  CHECK(Cell(8.0, electrostatic).regime == RegimeStatus::IN);
  CHECK(Cell(8.0, electrostatic).validated_min == 0.875);
  CHECK(Cell(5.6, electrostatic).regime == RegimeStatus::IN);  // 0.7 = 0.875 / 1.25
  CHECK(Cell(5.5, electrostatic).regime == RegimeStatus::OUT);
  CHECK(Cell(11.25, electrostatic).regime == RegimeStatus::IN);  // 1.40625 = 1.125 x 1.25
  CHECK(Cell(11.5, electrostatic).regime == RegimeStatus::OUT);
  // The range is per excitation class: an eigen cell of the same Type has none.
  CHECK(Cell(8.0, CoverageExcitationClass::EIGENMODE).regime == RegimeStatus::NO_RANGE);
  CHECK(RegimeStatusName(RegimeStatus::OUT) == "OUT");
}

TEST_CASE("Coverage key rows: per key, joint keys, method gap",
          "[SurfaceResponseCoverage][Serial]")
{
  std::map<int, Partition> cell;
  cell[1] = InterfaceOne();
  cell[2].classes = {{{"Missing", "m7"}, 0.5}, {{"Missing", "m3"}, 0.125}};
  const auto rows = CoverageKeyRows(cell, Keys());
  // Keys: hashA (m3 on 1 and 2, m7 on 2: one key), hashA+hashB (m3+m5), hashB (m5).
  REQUIRE(rows.size() == 3);
  CHECK(rows[0].key == "hashA");
  CHECK(rows[0].kind == "Missing+PlacementRejected");
  CHECK(rows[0].reason == "ConductorRefused");
  CHECK_FALSE(rows[0].method_gap);
  CHECK(rows[0].energy.at(1) == 1.0);
  CHECK(rows[0].energy.at(2) == 0.625);
  CHECK(rows[0].features == std::vector<int>{3, 7});
  CHECK(rows[1].key == "hashA+hashB");
  CHECK_FALSE(rows[1].method_gap);  // one token is Missing
  CHECK(rows[2].key == "hashB");
  CHECK(rows[2].kind == "Unboxable");
  CHECK(rows[2].method_gap);
  // The key rows sum to the Missing class per interface.
  double sum = 0.0;
  for (const auto &row : rows)
  {
    sum += row.energy.count(1) ? row.energy.at(1) : 0.0;
  }
  CHECK(sum == 1.75);
}

TEST_CASE("Coverage validated-regime file", "[SurfaceResponseCoverage][Serial]")
{
  namespace fs = std::filesystem;
  const auto dir = fs::temp_directory_path() / "palace-test-coverage-regime";
  fs::create_directories(dir);
  const auto library = dir / "library.json";
  const std::string library_text = "{\"Name\": \"test-library\"}\n";
  std::ofstream(library) << library_text;
  auto Write = [&](const nlohmann::json &data)
  {
    const auto path = dir / "validated-regime.json";
    std::ofstream(path) << data.dump(2);
    return path.string();
  };
  nlohmann::json data = {
      {"Version", 1},
      {"Library", {{"Name", "test-library"}, {"Sha256", Sha256Hex(library_text)}}},
      {"Tolerance", 0.0},
      {"Ranges",
       {{{"Type", "SA"},
         {"ExcitationClass", "ElectrostaticTerminal"},
         {"FactorMin", 0.86},
         {"FactorMax", 1.13},
         {"N", 6},
         {"Provenance", {{"Windows", {"FQ1", "FQ2"}}}}}}}};
  const auto regime = ReadCoverageValidatedRegime(Write(data), library.string());
  REQUIRE(regime.ranges.size() == 1);
  const auto &range =
      regime.ranges.at({"SA", CoverageExcitationClass::ELECTROSTATIC_TERMINAL});
  CHECK(range.factor_min == 0.86);
  CHECK(range.factor_max == 1.13);
  CHECK(range.n == 6);
  // Another library's range is a config error; so is a duplicate or an unknown class.
  auto other = data;
  other["Library"]["Sha256"] = Sha256Hex("another library");
  CHECK_THROWS(ReadCoverageValidatedRegime(Write(other), library.string()));
  auto twice = data;
  twice["Ranges"].push_back(data["Ranges"][0]);
  CHECK_THROWS(ReadCoverageValidatedRegime(Write(twice), library.string()));
  auto unknown = data;
  unknown["Ranges"][0]["ExcitationClass"] = "JunctionPort";
  CHECK_THROWS(ReadCoverageValidatedRegime(Write(unknown), library.string()));
  // Erratum 1784 (review m-6): a Driven range is refused by name (the class lumps every
  // driven excitation), also beside a valid range.
  auto driven = data;
  driven["Ranges"].push_back(data["Ranges"][0]);
  driven["Ranges"][1]["ExcitationClass"] = "Driven";
  CHECK_THROWS_WITH(ReadCoverageValidatedRegime(Write(driven), library.string()),
                    Catch::Matchers::ContainsSubstring(
                        "driven ranges are not supported in this version"));
  // Review m-1: N >= 1 and a Provenance per range, a bounded Tolerance.
  CHECK(range.provenance == data["Ranges"][0]["Provenance"].dump());
  for (const auto &n : {nlohmann::json(0), nlohmann::json(1.5), nlohmann::json("6")})
  {
    auto bad_n = data;
    bad_n["Ranges"][0]["N"] = n;
    CHECK_THROWS(ReadCoverageValidatedRegime(Write(bad_n), library.string()));
  }
  auto no_n = data;
  no_n["Ranges"][0].erase("N");
  CHECK_THROWS(ReadCoverageValidatedRegime(Write(no_n), library.string()));
  for (const auto &provenance :
       {nlohmann::json(nullptr), nlohmann::json(""), nlohmann::json::object()})
  {
    auto bad_provenance = data;
    bad_provenance["Ranges"][0]["Provenance"] = provenance;
    CHECK_THROWS(ReadCoverageValidatedRegime(Write(bad_provenance), library.string()));
  }
  auto no_provenance = data;
  no_provenance["Ranges"][0].erase("Provenance");
  CHECK_THROWS(ReadCoverageValidatedRegime(Write(no_provenance), library.string()));
  auto text_provenance = data;
  text_provenance["Ranges"][0]["Provenance"] = "windows FQ1 / FQ2";
  CHECK_NOTHROW(ReadCoverageValidatedRegime(Write(text_provenance), library.string()));
  auto edge_tolerance = data;
  edge_tolerance["Tolerance"] = kCoverageRegimeToleranceMax;
  CHECK_NOTHROW(ReadCoverageValidatedRegime(Write(edge_tolerance), library.string()));
  auto wide_tolerance = data;
  wide_tolerance["Tolerance"] = 0.5;
  CHECK_THROWS(ReadCoverageValidatedRegime(Write(wide_tolerance), library.string()));
  fs::remove_all(dir);
}

TEST_CASE("Coverage model qualification from the library records",
          "[SurfaceResponseCoverage][Serial]")
{
  // Note 1642: the library of record carries QualificationStatus / LibraryQualified only
  // on SpatialEdgeCluster records; LibraryQualified false with no status is PENDING.
  namespace fs = std::filesystem;
  const auto dir = fs::temp_directory_path() / "palace-test-coverage-qualification";
  fs::create_directories(dir);
  const auto library = dir / "library.json";
  const nlohmann::json data = {{"Models",
                                {{{"Name", "qualified"},
                                  {"QualificationStatus", "Qualified"},
                                  {"LibraryQualified", true}},
                                 {{"Name", "pending"},
                                  {"QualificationStatus", "PendingQualification"},
                                  {"LibraryQualified", false}},
                                 {{"Name", "unqualified"}, {"LibraryQualified", false}},
                                 {{"Name", "plain"}, {"Topology", "ConvexCorner"}}}}};
  std::ofstream(library) << data.dump();
  const auto qualifications = ReadCoverageModelQualifications(library.string());
  REQUIRE(qualifications.size() == 4);
  CHECK(qualifications.at("qualified") == CoverageModelQualification::QUALIFIED);
  CHECK(qualifications.at("pending") == CoverageModelQualification::PENDING);
  CHECK(qualifications.at("unqualified") == CoverageModelQualification::PENDING);
  CHECK(qualifications.at("plain") == CoverageModelQualification::NO_STATUS);
  fs::remove_all(dir);
}

}  // namespace palace

namespace palace
{

TEST_CASE("Coverage measurement per Type and the run record",
          "[SurfaceResponseCoverage][Serial]")
{
  // Interfaces 1 (SA) and 3 (SA) share a cell, 2 (MA) is its own; 4 is not a target.
  std::map<int, SurfacePostOperator::InterfaceCoveragePartition> partitions;
  partitions[1] = InterfaceOne();
  partitions[3] = InterfaceOne();
  partitions[2].within = 4.0;
  partitions[2].served = {{"m0", 4.0}};
  const std::map<int, std::string> types = {{1, "SA"}, {2, "MA"}, {3, "SA"}};
  const std::map<int, double> within = {{1, 16.0}, {2, 4.0}, {3, 16.0}};
  const std::map<int, double> fabricated = {{1, 13.5}, {2, 8.0}, {3, 13.5}};
  const auto measurement = ComputeCoverageMeasurement(
      partitions, types, within, fabricated, Keys(), Models(), CoverageThresholds{},
      CoverageValidatedRegime{}, CoverageExcitationClass::EIGENMODE);
  CHECK(measurement.measured);
  REQUIRE(measurement.cells.size() == 2);
  // First appearance over the ascending targets: SA (1), MA (2).
  CHECK(measurement.cells[0].type == "SA");
  CHECK(measurement.cells[1].type == "MA");
  const auto *sa = measurement.FindCell("SA");
  REQUIRE(sa);
  CHECK(sa->within == 32.0);
  CHECK(sa->served == 23.0);
  CHECK(sa->uncovered == 9.0);
  CHECK(sa->uncovered_share == 9.0 / 32.0);
  CHECK(sa->missing == 2.5);
  CHECK(sa->method_gap == 6.5);
  CHECK(sa->factor == 27.0 / 23.0);
  CHECK(sa->status == CoverageStatus::NOT_ACHIEVABLE);  // method gap 6.5 / 32 > 0.05
  CHECK(sa->regime == RegimeStatus::NO_RANGE);
  const auto *ma = measurement.FindCell("MA");
  REQUIRE(ma);
  CHECK(ma->uncovered_share == 0.0);
  CHECK(ma->status == CoverageStatus::OK);
  CHECK(ma->factor == 2.0);
  CHECK(measurement.WorstCoverageCell() == sa);
  CHECK(measurement.served_pending_models == std::set<int>{1});
  CHECK(measurement.served_no_status_models == std::set<int>{2});
  // The key rows: hashA (m3), hashB (m5, Unboxable), hashA+hashB (m3+m5), twice over.
  REQUIRE(measurement.keys.size() == 3);
  CHECK(measurement.keys[0].key == "hashA");
  CHECK(measurement.keys[0].energy.at(1) == 1.0);
  CHECK(measurement.keys[0].energy.at(3) == 1.0);

  const auto unmeasured =
      NotMeasuredCoverageMeasurement(types, CoverageExcitationClass::EIGENMODE);
  CHECK(!unmeasured.measured);
  REQUIRE(unmeasured.cells.size() == 2);
  CHECK(unmeasured.cells[0].status == CoverageStatus::NOT_MEASURED);
  CHECK(std::isnan(unmeasured.cells[0].uncovered_share));

  // The run record: a flagged entry per non-OK or OUT cell; nan as null; the banner names
  // the interim thresholds, every flagged cell and the status-less served models once.
  CoverageRunRecord record;
  record.Add({{"m", 1.0}}, measurement);
  record.Add({{"m", 2.0}}, unmeasured);
  const auto metadata = CoverageRunMetadata(record, CoverageThresholds{},
                                            CoverageValidatedRegime{}, 2.0, 2, 3);
  CHECK(metadata.at("Policy") == "Warn");
  // Erratum 1784 E-3: the block names the level that wrote it (2 of 3 here) and never
  // carries "Final" (written by BaseSolver after the AMR loop only).
  CHECK(metadata.at("Level") == 2);
  CHECK(metadata.at("Levels") == 3);
  CHECK_FALSE(metadata.contains("Final"));
  CHECK(metadata.at("Thresholds").at("Interim") == true);
  CHECK(metadata.at("ValidatedRegime").is_null());
  CHECK(metadata.at("Cells") == 4);
  CHECK(metadata.at("StatusCounts").at("OK") == 1);
  CHECK(metadata.at("StatusCounts").at("NOT-ACHIEVABLE") == 1);
  CHECK(metadata.at("StatusCounts").at("NOT-MEASURED") == 2);
  CHECK(metadata.at("RegimeCounts").at("NO-RANGE") == 2);
  CHECK(metadata.at("RegimeCounts").at("NO-MODEL") == 2);
  // Erratum 1784 E-2: the OK + NO-RANGE cell (m = 1, MA) is flagged too.
  const auto &flagged = metadata.at("Flagged");
  REQUIRE(flagged.size() == 4);
  CHECK(flagged[0].at("Cell").at("m") == 1.0);
  CHECK(flagged[0].at("Type") == "SA");
  CHECK(flagged[0].at("W within-R raw (J)") == 64.0);
  CHECK(flagged[0].at("validated factor min").is_null());
  CHECK(flagged[1].at("Type") == "MA");
  CHECK(flagged[1].at("coverage status") == "OK");
  CHECK(flagged[1].at("regime status") == "NO-RANGE");
  CHECK(flagged[2].at("uncovered raw share of W").is_null());
  CHECK(metadata.at("ServedPendingModels") == 1);
  CHECK(metadata.at("ServedNoStatusModels") == 1);
  const std::string banner = CoverageRunBanner(metadata);
  CHECK(banner.find("interim thresholds") != std::string::npos);
  CHECK(banner.find("*** COVERAGE NOT-ACHIEVABLE (interim thresholds): m = 1 SA: U = "
                    "28.125 % of W") != std::string::npos);
  CHECK(banner.find("*** COVERAGE NOT-MEASURED (interim thresholds): m = 2 SA") !=
        std::string::npos);
  CHECK(banner.find("2 of 4 cell(s) NO-RANGE (no validated range of the covered correction "
                    "factor for this library and excitation class: computed, NOT VALIDATED "
                    "in this regime)") != std::string::npos);
  CHECK(banner.find("1 served model(s) carry no QualificationStatus") != std::string::npos);
  // Erratum 1784 E-2 / review m-5: an OK + NO-RANGE cell has its own line, with the words
  // of ERRATUM-1581 MJ-2.
  CHECK(banner.find("*** COVERAGE OK (interim thresholds): m = 1 MA: U = 0.000 % of W") !=
        std::string::npos);
  CHECK(banner.find("covered factor 2.0000 (no validated range: NO-RANGE: computed, NOT "
                    "VALIDATED in this regime)") != std::string::npos);
  // An OK + IN cell is not flagged.
  CoverageValidatedRegime ma_range;
  ma_range.ranges[{"MA", CoverageExcitationClass::EIGENMODE}] = {1.5, 2.5, 1, "{}"};
  const auto validated = ComputeCoverageMeasurement(
      partitions, types, within, fabricated, Keys(), Models(), CoverageThresholds{},
      ma_range, CoverageExcitationClass::EIGENMODE);
  REQUIRE(validated.FindCell("MA"));
  CHECK(validated.FindCell("MA")->regime == RegimeStatus::IN);
  CoverageRunRecord validated_record;
  validated_record.Add({{"m", 1.0}}, validated);
  const auto validated_metadata =
      CoverageRunMetadata(validated_record, CoverageThresholds{}, ma_range, 2.0, 1, 1);
  REQUIRE(validated_metadata.at("Flagged").size() == 1);
  CHECK(validated_metadata.at("Flagged")[0].at("Type") == "SA");
  CHECK(CoverageRunBanner(validated_metadata).find("COVERAGE OK") == std::string::npos);
  // Review m-8: the banner prints at most 40 cell lines, then one line counting the rest
  // (palace.json's Flagged lists every flagged cell).
  CoverageRunRecord many;
  for (int m = 1; m <= 21; m++)
  {
    many.Add({{"m", double(m)}}, measurement);  // 2 flagged cells each
  }
  const auto many_metadata =
      CoverageRunMetadata(many, CoverageThresholds{}, CoverageValidatedRegime{}, 2.0, 1, 1);
  CHECK(many_metadata.at("Flagged").size() == 42);
  const std::string many_banner = CoverageRunBanner(many_metadata);
  int cell_lines = 0;
  for (std::size_t at = many_banner.find("(interim thresholds): m = ");
       at != std::string::npos; at = many_banner.find("(interim thresholds): m = ", at + 1))
  {
    cell_lines++;
  }
  CHECK(cell_lines == 40);
  CHECK(many_banner.find("*** COVERAGE: 2 more flagged cell(s): see "
                         "SurfaceResponse.Coverage.Flagged of palace.json ***") !=
        std::string::npos);
}

TEST_CASE("Coverage table rows", "[SurfaceResponseCoverage][Serial]")
{
  std::map<int, SurfacePostOperator::InterfaceCoveragePartition> partitions;
  partitions[1] = InterfaceOne();
  const auto measurement = ComputeCoverageMeasurement(
      partitions, {{1, "SA"}}, {{1, 16.0}}, {{1, 13.5}}, Keys(), Models(),
      CoverageThresholds{}, CoverageValidatedRegime{}, CoverageExcitationClass::DRIVEN);
  const Units units(1.0, 1.0);
  std::ostringstream rows, keys;
  WriteCoverageRows(rows, "1,0", measurement, units);
  const std::string header = CoverageTableHeader("m,exc");
  CHECK(header.rfind("m,exc,Type,W within-R raw (J),", 0) == 0);
  CHECK(header.find(",coverage status,regime status,partition within-R raw (J),untargeted "
                    "raw energy of the Type (J),untargeted raw share of the Type interface "
                    "energy,loop remainder unserved (J),loop remainder of excluded patches "
                    "kept raw twice (J),method-gap raw share of W signed (before the clamp "
                    "at 0)\n") != std::string::npos);
  const std::string row = rows.str();
  CHECK(std::count(header.begin(), header.end(), ',') ==
        std::count(row.begin(), row.end(), ','));
  CHECK(row.find(fmt::format(
            ",nan,nan,NOT-ACHIEVABLE,NO-RANGE,{:+.12e},+0.000000000000e+00,"
            "+0.000000000000e+00,+0.000000000000e+00,+0.000000000000e+00,{:+.12e}\n",
            units.Dimensionalize<Units::ValueType::ENERGY>(16.0),
            measurement.cells.at(0).method_gap_signed_share)) != std::string::npos);
  CHECK(measurement.cells.at(0).method_gap_signed_share ==
        measurement.cells.at(0).method_gap_share);
  WriteCoverageByKeyRows(keys, "1,0", measurement, {1},
                         {{3, {2, 0.5, false}}, {5, {1, 0.25, true}}}, units);
  const std::string by_key_header = CoverageByKeyTableHeader("m,exc", {1});
  CHECK(by_key_header ==
        "m,exc,evaluation,key hash,class,reason,method gap,features,portions,length "
        "(m),raw within-R energy[1] (J),raw share of W[1]\n");
  std::istringstream lines(keys.str());
  std::string line;
  std::vector<std::string> parsed;
  while (std::getline(lines, line))
  {
    parsed.push_back(line);
  }
  REQUIRE(parsed.size() == 3);
  CHECK(parsed[0].rfind("1,0,0,hashA,Missing,-,0,3,2,", 0) == 0);
  CHECK(parsed[1].rfind("1,0,0,hashA+hashB,Missing+Unboxable,MirrorFormed,0,3+5,3,", 0) ==
        0);
  CHECK(parsed[2].rfind("1,0,0,hashB,Unboxable,MirrorFormed,1,5,1,", 0) == 0);
  CHECK(parsed[2].find(fmt::format(",{:+.12e},+3.125000000000e-02",
                                   units.Dimensionalize<Units::ValueType::ENERGY>(0.5))) !=
        std::string::npos);
}

// Ruling 1694 (1) + review M-4 / m-6: the untargeted interfaces of a Type (their raw
// energy and its share of the Type's raw interface energy; a NOT-TARGET cell for a Type
// with untargeted interfaces only: W = 0, shares / factor nan, NO-MODEL, never the worst
// cell, never flagged, named once), the loop-remainder terms summed per Type, the run
// record's cell set with "Complete", a null factor in the banner, an empty record's line.
TEST_CASE("Coverage untargeted interfaces, the cell set and the banner edge cases",
          "[SurfaceResponseCoverage][Serial]")
{
  std::map<int, Partition> partitions;
  partitions[1] = InterfaceOne();
  partitions[1].loop_remainder_unserved = 0.5;
  auto measurement =
      ComputeCoverageMeasurement(partitions, {{1, "SA"}}, {{1, 16.0}}, {{1, 13.5}}, Keys(),
                                 Models(), CoverageThresholds{}, CoverageValidatedRegime{},
                                 CoverageExcitationClass::EIGENMODE, {{1, 0.25}});
  AddCoverageUntargetedInterfaces(measurement, {{2, "SA"}, {3, "MS"}},
                                  {{1, 24.0}, {2, 8.0}, {3, 4.0}});
  REQUIRE(measurement.cells.size() == 2);
  const auto &sa = measurement.cells[0];
  CHECK(sa.type == "SA");
  CHECK(sa.untargeted_raw == 8.0);
  CHECK(sa.untargeted_share == 0.25);  // 8 / (24 + 8)
  CHECK(sa.loop_remainder_unserved == 0.5);
  CHECK(sa.loop_remainder_excluded == 0.25);
  CHECK(sa.partition_within == 16.0);
  CHECK(sa.status == CoverageStatus::NOT_ACHIEVABLE);
  const auto &ms = measurement.cells[1];
  CHECK(ms.type == "MS");
  CHECK(ms.status == CoverageStatus::NOT_TARGET);
  CHECK(ms.regime == RegimeStatus::NO_MODEL);
  CHECK(ms.within == 0.0);
  for (const double value :
       {ms.uncovered_share, ms.missing_share, ms.method_gap_share, ms.pending_share,
        ms.no_status_share, ms.factor, ms.validated_min, ms.validated_max})
  {
    CHECK(std::isnan(value));
  }
  CHECK(ms.untargeted_raw == 4.0);
  CHECK(ms.untargeted_share == 1.0);
  CHECK(measurement.untargeted.size() == 2);
  CHECK(measurement.WorstCoverageCell() == &sa);
  CHECK(measurement.WorstRegimeCell() == &sa);

  CoverageRunRecord record;
  record.Add({{"m", 1.0}}, measurement);
  record.Add({{"m", 2.0}}, measurement);
  const auto metadata = CoverageRunMetadata(record, CoverageThresholds{},
                                            CoverageValidatedRegime{}, 1.0, 1, 1);
  CHECK(metadata.at("Complete") == true);
  const auto &cell_set = metadata.at("CellSet");
  CHECK(cell_set.at("Kind") == "m");
  CHECK(cell_set.at("Measurements") == 2);
  CHECK(cell_set.at("Cells") == 4);
  CHECK(cell_set.at("Types") == nlohmann::json::array({"SA", "MS"}));
  CHECK(metadata.at("StatusCounts").at("NOT-TARGET") == 2);
  REQUIRE(metadata.at("UntargetedInterfaces").size() == 2);
  CHECK(metadata.at("UntargetedInterfaces")[0].at("Interface") == 2);
  CHECK(metadata.at("UntargetedInterfaces")[1].at("Type") == "MS");
  REQUIRE(metadata.at("Flagged").size() == 2);
  for (const auto &entry : metadata.at("Flagged"))
  {
    CHECK(entry.at("coverage status") == "NOT-ACHIEVABLE");
    CHECK(entry.at("untargeted raw share of the Type interface energy") == 0.25);
  }
  const std::string banner = CoverageRunBanner(metadata);
  CHECK(banner.find("2 untargeted postprocessing interface(s) 2 (SA), 3 (MS) carry their "
                    "raw energy; not corrected") != std::string::npos);
  // Review m-6: a null factor / range prints nan (never a JSON type error).
  auto with_null = metadata;
  with_null["Flagged"][0]["covered correction factor"] = nullptr;
  with_null["Flagged"][0]["regime status"] = "OUT";
  std::string null_banner;
  CHECK_NOTHROW(null_banner = CoverageRunBanner(with_null));
  CHECK(null_banner.find("covered factor nan (validated [nan, nan] by 0 window(s): OUT)") !=
        std::string::npos);
  // Review m-6: an empty record is never silent.
  const auto empty = CoverageRunMetadata(CoverageRunRecord{}, CoverageThresholds{},
                                         CoverageValidatedRegime{}, 1.0, 1, 1);
  CHECK(empty.at("Cells") == 0);
  // Ruling 1738: zero cells without a reason is NOT complete; with the reason it is, and
  // the reason is printed.
  CHECK(empty.at("Complete") == false);
  CHECK(empty.at("NoCellsReason").is_null());
  CHECK(empty.at("ConfiguredTargetsNotServed").empty());
  CHECK(CoverageRunBanner(empty).find("no coverage cell was measured") !=
        std::string::npos);
  CoverageRunRecord fallback;
  fallback.no_cells_reason = "RawFallback: NoMatchedFeature: no correction applied";
  const auto reasoned = CoverageRunMetadata(fallback, CoverageThresholds{},
                                            CoverageValidatedRegime{}, 1.0, 1, 1);
  CHECK(reasoned.at("Complete") == true);
  CHECK(reasoned.at("NoCellsReason") == fallback.no_cells_reason);
  CHECK(CoverageRunBanner(reasoned).find("(RawFallback: NoMatchedFeature: no correction "
                                         "applied)") != std::string::npos);
}

// Ruling 1738 (OI-CONFIGURED-TARGET-NOT-SERVED): a configured target no loaded model serves
// is a coverage interface of its Type (its within-R raw energy in W, uncovered: U counts
// it), the Type's cell reads TARGET-NOT-SERVED whatever its band (ranked above every band,
// below NO-ENERGY), a Type without one keeps its band; palace.json lists the interface with
// its reason, the banner names it, the startup warning names it.
TEST_CASE("Coverage configured target not served", "[SurfaceResponseCoverage][Serial]")
{
  CHECK(CoverageStatusName(CoverageStatus::TARGET_NOT_SERVED) == "TARGET-NOT-SERVED");
  CHECK(CoverageStatusRank(CoverageStatus::TARGET_NOT_SERVED) >
        CoverageStatusRank(CoverageStatus::NOT_ACHIEVABLE));
  CHECK(CoverageStatusRank(CoverageStatus::TARGET_NOT_SERVED) <
        CoverageStatusRank(CoverageStatus::NO_ENERGY));
  std::map<int, Partition> partitions;
  partitions[1] = InterfaceOne();  // SA, served
  Partition unserved;              // SA, configured, no model surface: nothing served
  unserved.within = 4.0;
  unserved.classes = {{{"Missing", "m3"}, 3.0}};
  unserved.unclaimed = 1.0;
  partitions[2] = unserved;
  Partition full;  // MS, served whole
  full.within = 10.0;
  full.served = {{"m0", 10.0}};
  partitions[3] = full;
  auto measurement = ComputeCoverageMeasurement(
      partitions, {{1, "SA"}, {2, "SA"}, {3, "MS"}}, {{1, 16.0}, {2, 4.0}, {3, 10.0}},
      {{1, 13.5}, {3, 12.0}}, Keys(), Models(), CoverageThresholds{},
      CoverageValidatedRegime{}, CoverageExcitationClass::EIGENMODE, {});
  MarkCoverageTargetsNotServed(measurement, {{2, "SA"}});
  CHECK(measurement.not_served == std::map<int, std::string>{{2, "SA"}});
  const auto *sa = measurement.FindCell("SA");
  const auto *ms = measurement.FindCell("MS");
  REQUIRE(sa);
  REQUIRE(ms);
  CHECK(sa->within == 20.0);
  CHECK(std::abs(sa->uncovered_share - 8.5 / 20.0) <= 1.0e-12);  // 4.5 + all of 4
  CHECK(sa->status == CoverageStatus::TARGET_NOT_SERVED);
  CHECK(ms->status == CoverageStatus::OK);
  CHECK(measurement.WorstCoverageCell() == sa);
  CoverageRunRecord record;
  record.Add({{"m", 1.0}}, measurement);
  const auto metadata = CoverageRunMetadata(record, CoverageThresholds{},
                                            CoverageValidatedRegime{}, 1.0, 1, 1);
  const auto &listed = metadata.at("ConfiguredTargetsNotServed");
  REQUIRE(listed.size() == 1);
  CHECK(listed[0].at("Interface") == 2);
  CHECK(listed[0].at("Type") == "SA");
  CHECK(listed[0].at("Reason") == kCoverageTargetNotServedReason);
  CHECK(metadata.at("StatusCounts").at("TARGET-NOT-SERVED") == 1);
  CHECK(metadata.at("Complete") == true);
  CHECK(CoverageRunBanner(metadata).find(
            "*** COVERAGE TARGET-NOT-SERVED: configured target interface 2 (SA): no loaded "
            "model has a fabricated surface on it; it is NOT corrected") !=
        std::string::npos);
  CHECK(CoverageTargetsNotServedWarning({{2, "SA"}}) ==
        "Configured response-correction target interface 2 (SA) is NOT SERVED: no loaded "
        "model has a fabricated surface on it; it is not corrected (its corrected values "
        "are its raw values) and its within-R raw energy counts as uncovered in its Type's "
        "coverage (status TARGET-NOT-SERVED)\n");
  CHECK(CoverageTargetsNotServedWarning({}).empty());
}

// Ruling 1682: the compact claim list returns every appended value BIT FOR BIT in the
// appended order; end points are shared only by exact bit equality (-0.0 and +0.0 are two
// points; a value 1 ulp away is another point); each label is held once; appending after
// FinishAppending continues the same tables; a NaN end is refused.
TEST_CASE("Coverage compact claim list is exact", "[SurfaceResponseCoverage][Serial]")
{
  using Portion = SurfacePostOperator::UncoveredPerimeterPortion;
  const double next = std::nextafter(0.3, 1.0);
  const std::vector<Portion> plain = {
      {{0.1, 0.2, 0.3}, {0.4, 0.5, 0.6}, "m1", 7},
      {{0.4, 0.5, 0.6}, {0.7, 0.8, 0.9}, "m1", 7},          // shares its first end
      {{0.1, 0.2, next}, {0.0, 0.0, 0.0}, "m1+m2", 8},      // 1 ulp from the first point
      {{-0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, "unpatched", 9},  // -0.0 kept apart from +0.0
      {{0.7, 0.8, 0.9}, {0.1, 0.2, 0.3}, "m1", -1}};
  SurfacePostOperator::CompactPerimeterPortions compact;
  for (std::size_t i = 0; i < 3; i++)
  {
    compact.Append(plain[i].p0, plain[i].p1, plain[i].type, plain[i].feature);
  }
  compact.FinishAppending();
  for (std::size_t i = 3; i < plain.size(); i++)
  {
    compact.Append(plain[i].p0, plain[i].p1, plain[i].type, plain[i].feature);
  }
  compact.FinishAppending();
  REQUIRE(compact.Size() == plain.size());
  CHECK(compact.LabelCount() == 3);
  // Distinct points: (0.1,0.2,0.3), (0.4,0.5,0.6), (0.7,0.8,0.9), (0.1,0.2,next), (0,0,0),
  // (-0,0,0).
  CHECK(compact.PointCount() == 6);
  auto Bits = [](const std::array<double, 3> &p)
  {
    std::array<std::uint64_t, 3> bits;
    std::memcpy(bits.data(), p.data(), sizeof(bits));
    return bits;
  };
  const auto expanded = compact.Expand();
  const SurfacePostOperator::PerimeterPortionView view(compact), plain_view(plain);
  REQUIRE(expanded.size() == plain.size());
  for (std::size_t i = 0; i < plain.size(); i++)
  {
    INFO("portion " << i);
    CHECK(Bits(compact.P0(i)) == Bits(plain[i].p0));
    CHECK(Bits(compact.P1(i)) == Bits(plain[i].p1));
    CHECK(compact.Label(i) == plain[i].type);
    CHECK(compact.Feature(i) == plain[i].feature);
    CHECK(Bits(expanded[i].p0) == Bits(plain[i].p0));
    CHECK(Bits(expanded[i].p1) == Bits(plain[i].p1));
    CHECK(expanded[i].type == plain[i].type);
    CHECK(expanded[i].feature == plain[i].feature);
    CHECK(Bits(view.P0(i)) == Bits(plain_view.P0(i)));
    CHECK(Bits(view.P1(i)) == Bits(plain_view.P1(i)));
    CHECK(view.Label(i) == plain_view.Label(i));
    CHECK(view.Feature(i) == plain_view.Feature(i));
  }
  CHECK(std::signbit(compact.P0(3)[0]));
  CHECK_FALSE(std::signbit(compact.P1(3)[0]));
  // The trimmed buffers: 6 points x 24 B + 5 portions x 16 B + 3 labels (in the string
  // objects' small buffers).
  CHECK(compact.HeapBytes() == 6 * 24 + 5 * 16 + 3 * sizeof(std::string));
  CHECK_THROWS(compact.Append({0.0, std::nan(""), 0.0}, {0.0, 0.0, 0.0}, "m1", 1));
}

}  // namespace palace
