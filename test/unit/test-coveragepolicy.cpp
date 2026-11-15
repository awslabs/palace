// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include "models/surfaceresponsecoverage.hpp"

// Lane L4 (decision 1829; DESIGN-L4-v2 s2; PREREG-L4-IMPL s2): the mark-not-applied
// policy's marked set, the measurement marks, the run metadata and banner, the
// terminal-matrix marks (ERRATUM-1829 E-1) and the CTX-002 replay cells (prereg P9) through
// the C++ status rule.

namespace palace
{

namespace
{

constexpr CoverageStatus kAllStatuses[] = {CoverageStatus::OK,
                                           CoverageStatus::LOW,
                                           CoverageStatus::PARTIAL,
                                           CoverageStatus::NOT_ACHIEVABLE,
                                           CoverageStatus::NO_ENERGY,
                                           CoverageStatus::NOT_MEASURED,
                                           CoverageStatus::NOT_TARGET,
                                           CoverageStatus::INCONSISTENT,
                                           CoverageStatus::TARGET_NOT_SERVED};

CoverageCell Cell(const std::string &type, CoverageStatus status, double within,
                  double uncovered)
{
  CoverageCell cell;
  cell.type = type;
  cell.status = status;
  cell.within = within;
  cell.uncovered = uncovered;
  cell.uncovered_share = within > 0.0 ? uncovered / within : std::nan("");
  cell.regime = RegimeStatus::NO_RANGE;
  return cell;
}

CoverageMeasurement Measurement(std::vector<CoverageCell> cells)
{
  CoverageMeasurement measurement;
  measurement.measured = true;
  measurement.cells = std::move(cells);
  return measurement;
}

}  // namespace

TEST_CASE("Coverage policy marked set", "[CoveragePolicy][Serial]")
{
  // DESIGN-L4-v2 s2.1 (Q2 = the recommended "yes"): five statuses mark, four never do.
  for (const auto status : kAllStatuses)
  {
    const bool expected = status == CoverageStatus::PARTIAL ||
                          status == CoverageStatus::NOT_ACHIEVABLE ||
                          status == CoverageStatus::INCONSISTENT ||
                          status == CoverageStatus::TARGET_NOT_SERVED ||
                          status == CoverageStatus::NOT_MEASURED;
    INFO(std::string(CoverageStatusName(status)));
    CHECK(IsCoverageStatusMarked(status) == expected);
  }
  // A regime status never marks: the same coverage status marks alike under every regime.
  for (const auto regime : {RegimeStatus::IN, RegimeStatus::OUT, RegimeStatus::NO_RANGE,
                            RegimeStatus::NO_MODEL})
  {
    auto ok = Cell("SA", CoverageStatus::OK, 1.0, 0.0);
    auto partial = Cell("MS", CoverageStatus::PARTIAL, 1.0, 0.5);
    ok.regime = partial.regime = regime;
    const auto measurement = Measurement({ok, partial});
    CHECK_FALSE(measurement.IsTypeMarked("SA"));
    CHECK(measurement.IsTypeMarked("MS"));
  }
  CHECK(std::string(kCoveragePartialColumnSuffix) == " partial (uncovered part raw)");
  CHECK_FALSE(kCoverageMarkWholeMeasurement);
}

TEST_CASE("Coverage policy measurement marks", "[CoveragePolicy][Serial]")
{
  // An OK and a LOW Type: nothing marked, nothing uncorrected.
  const auto covered = Measurement({Cell("SA", CoverageStatus::OK, 4.0, 0.0625),
                                    Cell("MS", CoverageStatus::LOW, 8.0, 0.25)});
  CHECK_FALSE(covered.AnyMarked());
  CHECK_FALSE(covered.MarksWholeMeasurement());
  CHECK_FALSE(covered.IsTypeMarked("SA"));
  CHECK(covered.MarkedUncoveredEnergy() == 0.0);
  // A PARTIAL Type beside an OK one (Q3 keep + flag): only its own Type's interfaces are
  // marked; the whole values are; the exposure = its U W.
  const auto mixed = Measurement({Cell("SA", CoverageStatus::OK, 4.0, 0.0625),
                                  Cell("MS", CoverageStatus::PARTIAL, 8.0, 2.0),
                                  Cell("MA", CoverageStatus::NO_ENERGY, 0.0, 0.0),
                                  Cell("D", CoverageStatus::NOT_TARGET, 0.0, 0.0)});
  CHECK(mixed.AnyMarked());
  CHECK(mixed.MarksWholeMeasurement());
  CHECK_FALSE(mixed.IsTypeMarked("SA"));
  CHECK(mixed.IsTypeMarked("MS"));
  CHECK_FALSE(mixed.IsTypeMarked("MA"));
  CHECK_FALSE(mixed.IsTypeMarked("D"));
  CHECK_FALSE(mixed.IsTypeMarked("absent"));
  CHECK(mixed.MarkedUncoveredEnergy() == 2.0);
  // Two marked Types sum; a NOT-MEASURED cell (no U) makes the exposure unknown (nan).
  const auto two = Measurement({Cell("SA", CoverageStatus::NOT_ACHIEVABLE, 4.0, 1.0),
                                Cell("MS", CoverageStatus::TARGET_NOT_SERVED, 8.0, 0.5)});
  CHECK(two.MarkedUncoveredEnergy() == 1.5);
  const auto unmeasured = NotMeasuredCoverageMeasurement(
      {{1, "SA"}, {2, "MS"}}, CoverageExcitationClass::EIGENMODE);
  CHECK(unmeasured.AnyMarked());
  CHECK(unmeasured.IsTypeMarked("SA"));
  CHECK(std::isnan(unmeasured.MarkedUncoveredEnergy()));
  // No cell at all: the whole values are marked (fail-closed), no Type is.
  const CoverageMeasurement empty;
  CHECK_FALSE(empty.AnyMarked());
  CHECK(empty.MarksWholeMeasurement());
}

TEST_CASE("Coverage policy run metadata and banner", "[CoveragePolicy][Serial]")
{
  CoverageRunRecord record;
  record.Add({{"m", 1.0}}, Measurement({Cell("SA", CoverageStatus::OK, 4.0, 0.0),
                                        Cell("MS", CoverageStatus::PARTIAL, 8.0, 2.0)}));
  record.Add({{"m", 2.0}}, Measurement({Cell("SA", CoverageStatus::OK, 4.0, 0.0),
                                        Cell("MS", CoverageStatus::LOW, 8.0, 0.25)}));
  record.Add({{"m", 3.0}}, NotMeasuredCoverageMeasurement(
                               {{1, "SA"}}, CoverageExcitationClass::EIGENMODE));
  // Warn: the block and the banner of the default argument, without any new key or line.
  const auto warn = CoverageRunMetadata(record, CoverageThresholds{},
                                        CoverageValidatedRegime{}, 1.0, 1, 1);
  const auto warn_explicit = CoverageRunMetadata(
      record, CoverageThresholds{}, CoverageValidatedRegime{}, 1.0, 1, 1, false);
  CHECK(warn == warn_explicit);
  CHECK(warn.at("Policy") == "Warn");
  for (const char *key : {"MarkedStatuses", "MarkedCells", "MarkedMeasurements"})
  {
    CHECK_FALSE(warn.contains(key));
  }
  for (const auto &entry : warn.at("Flagged"))
  {
    CHECK_FALSE(entry.contains("Marked"));
  }
  const std::string warn_banner = CoverageRunBanner(warn);
  CHECK(warn_banner.find("MARK-NOT-APPLIED") == std::string::npos);
  CHECK(warn_banner.find("[NOT APPLIED]") == std::string::npos);
  CHECK(warn_banner.find("Coverage labels") == std::string::npos);

  // MarkNotApplied: the same block plus the policy keys; the marked cells are the PARTIAL
  // one (m = 1) and the NOT-MEASURED one (m = 3).
  const auto marked = CoverageRunMetadata(record, CoverageThresholds{},
                                          CoverageValidatedRegime{}, 1.0, 1, 1, true);
  CHECK(marked.at("Policy") == "MarkNotApplied");
  CHECK(marked.at("MarkedStatuses") ==
        nlohmann::json({"PARTIAL", "NOT-ACHIEVABLE", "TARGET-NOT-SERVED", "NOT-MEASURED",
                        "INCONSISTENT"}));
  CHECK(marked.at("MarkedCells") == 2);
  CHECK(marked.at("MarkedMeasurements") == 2);
  auto stripped = marked;
  for (const char *key : {"MarkedStatuses", "MarkedCells", "MarkedMeasurements"})
  {
    stripped.erase(key);
  }
  stripped["Policy"] = "Warn";
  int marked_flags = 0;
  for (auto &entry : stripped.at("Flagged"))
  {
    REQUIRE(entry.contains("Marked"));
    const bool expected = IsCoverageStatusMarked(
        entry.at("coverage status") == "PARTIAL"
            ? CoverageStatus::PARTIAL
            : (entry.at("coverage status") == "NOT-MEASURED" ? CoverageStatus::NOT_MEASURED
                                                             : CoverageStatus::OK));
    CHECK(entry.at("Marked") == expected);
    marked_flags += entry.at("Marked").get<bool>();
    entry.erase("Marked");
  }
  CHECK(marked_flags == 2);
  CHECK(stripped == warn);  // nothing else in the block moves
  const std::string banner = CoverageRunBanner(marked);
  CHECK(banner.find("*** MARK-NOT-APPLIED: 2 cell(s) print nan in their corrected columns "
                    "(statuses PARTIAL, NOT-ACHIEVABLE, TARGET-NOT-SERVED, NOT-MEASURED, "
                    "INCONSISTENT); the values are in the '... partial (uncovered part "
                    "raw)' columns; 2 measurement(s)' corrected totals") !=
        std::string::npos);
  CHECK(banner.find("Coverage labels: COVERAGE = ") != std::string::npos);
  CHECK(banner.find("REGIME (validated-range status)") != std::string::npos);
  CHECK(banner.find("*** COVERAGE PARTIAL (interim thresholds): m = 1 MS") !=
        std::string::npos);
  CHECK(banner.find(" *** [NOT APPLIED]\n") != std::string::npos);
  // The LOW line (m = 2) is not marked; every line of the Warn banner is in this one.
  const auto low = banner.find("*** COVERAGE LOW (interim thresholds): m = 2 MS");
  REQUIRE(low != std::string::npos);
  CHECK(banner.substr(low, banner.find('\n', low) - low).find("[NOT APPLIED]") ==
        std::string::npos);
  std::size_t from = 0;
  int count = 0;
  while ((from = banner.find("[NOT APPLIED]", from)) != std::string::npos)
  {
    count++;
    from++;
  }
  CHECK(count == 2);
  std::string unmarked = banner;
  for (std::size_t at; (at = unmarked.find(" [NOT APPLIED]")) != std::string::npos;)
  {
    unmarked.erase(at, std::string(" [NOT APPLIED]").size());
  }
  // The banner starts with an empty line; its first line is the summary line.
  const auto first = unmarked.find('\n', unmarked.find("Surface-response coverage")) + 1;
  const auto extra_end = unmarked.find("Coverage: the share denominators");
  REQUIRE(extra_end != std::string::npos);
  unmarked.erase(first, unmarked.find('\n', extra_end) + 1 - first);
  CHECK(unmarked == warn_banner);
}

TEST_CASE("Coverage policy terminal matrix marks", "[CoveragePolicy][Serial]")
{
  // ERRATUM-1829 E-1: source 1 OK, source 2 PARTIAL + NO-ENERGY, source 3 PARTIAL +
  // TARGET-NOT-SERVED. The coverage table's worse-of-two status of source 2 is NO-ENERGY
  // (rank 5 > PARTIAL 2): unmarked by the status table. The any-cell rule marks it.
  const std::vector<CoverageMeasurement> sources = {
      Measurement({Cell("SA", CoverageStatus::OK, 1.0, 0.0),
                   Cell("MS", CoverageStatus::OK, 1.0, 0.0)}),
      Measurement({Cell("SA", CoverageStatus::PARTIAL, 1.0, 0.5),
                   Cell("MS", CoverageStatus::NO_ENERGY, 0.0, 0.0)}),
      Measurement({Cell("SA", CoverageStatus::PARTIAL, 1.0, 0.5),
                   Cell("MS", CoverageStatus::TARGET_NOT_SERVED, 1.0, 1.0)})};
  REQUIRE(sources[1].WorstCoverageCell());
  CHECK(sources[1].WorstCoverageCell()->status == CoverageStatus::NO_ENERGY);
  CHECK_FALSE(IsCoverageStatusMarked(sources[1].WorstCoverageCell()->status));
  CHECK(sources[2].WorstCoverageCell()->status == CoverageStatus::TARGET_NOT_SERVED);
  std::vector<bool> source_marked;
  for (const auto &source : sources)
  {
    source_marked.push_back(source.MarksWholeMeasurement());
  }
  CHECK(source_marked == std::vector<bool>{false, true, true});
  const auto marks = ComputeCoverageTerminalMatrixMarks(source_marked);
  const std::vector<std::vector<bool>> C = {
      {false, true, true}, {true, true, true}, {true, true, true}};
  CHECK(marks.capacitance == C);
  // C_m: off-diagonal as C; the diagonal of row 1 sums C(1, 2) and C(1, 3): marked.
  CHECK(marks.mutual == std::vector<std::vector<bool>>{
                            {true, true, true}, {true, true, true}, {true, true, true}});
  CHECK(marks.inverse == std::vector<std::vector<bool>>(3, std::vector<bool>(3, true)));
  // Nothing marked: no entry of any matrix is.
  const auto none = ComputeCoverageTerminalMatrixMarks({false, false});
  CHECK(none.capacitance == std::vector<std::vector<bool>>(2, std::vector<bool>(2, false)));
  CHECK(none.mutual == none.capacitance);
  CHECK(none.inverse == none.capacitance);
  // One source of two: C(1, 1) unmarked, the rest marked; the inverse wholly marked.
  const auto one = ComputeCoverageTerminalMatrixMarks({false, true});
  CHECK(one.capacitance == std::vector<std::vector<bool>>{{false, true}, {true, true}});
  CHECK(one.mutual == std::vector<std::vector<bool>>{{true, true}, {true, true}});
  CHECK(one.inverse == std::vector<std::vector<bool>>(2, std::vector<bool>(2, true)));
}

TEST_CASE("Coverage policy CTX-002 replay cells", "[CoveragePolicy][Serial]")
{
  // PREREG-L4-IMPL P9 (ERRATUM-1829 m-2): the 33 cells of PBS 111274 (W, U, Missing,
  // pending per mode x Type from the frozen-reader reports on the tables of record; U on
  // the reader basis) through the branch's ComputeCoverageCell and IsCoverageStatusMarked.
  // The expected status and mark are the offline replay's; they are not inputs of the
  // computation.
  struct Row
  {
    int mode;
    const char *type;
    double within, uncovered_share, missing_share, pending_share;
    const char *status;
    bool marked;
  };
  const std::vector<Row> rows = {
      {1, "SA", 2.829418e-15, 0.000182, 0.000000, 0.000e+00, "OK", false},
      {1, "MS", 2.830970e-15, 0.003880, 0.000000, 0.000e+00, "OK", false},
      {1, "MA", 8.448309e-17, 0.000307, 0.000000, 0.000e+00, "OK", false},
      {2, "SA", 5.646148e-16, 0.594594, 0.594594, 0.000e+00, "PARTIAL", true},
      {2, "MS", 3.756847e-16, 0.230197, 0.229764, 0.000e+00, "PARTIAL", true},
      {2, "MA", 1.568160e-17, 0.430407, 0.430406, 0.000e+00, "PARTIAL", true},
      {3, "SA", 5.593909e-16, 0.586515, 0.586514, 0.000e+00, "PARTIAL", true},
      {3, "MS", 3.842151e-16, 0.215832, 0.215381, 0.000e+00, "PARTIAL", true},
      {3, "MA", 1.529191e-17, 0.457135, 0.457135, 0.000e+00, "PARTIAL", true},
      {4, "SA", 5.596392e-16, 0.597282, 0.597282, 0.000e+00, "PARTIAL", true},
      {4, "MS", 3.636122e-16, 0.216991, 0.216534, 0.000e+00, "PARTIAL", true},
      {4, "MA", 1.549891e-17, 0.430325, 0.430324, 0.000e+00, "PARTIAL", true},
      {5, "SA", 5.531573e-16, 0.587550, 0.587550, 0.000e+00, "PARTIAL", true},
      {5, "MS", 3.757004e-16, 0.224215, 0.223785, 0.000e+00, "PARTIAL", true},
      {5, "MA", 1.544733e-17, 0.449544, 0.449543, 0.000e+00, "PARTIAL", true},
      {6, "SA", 5.533853e-16, 0.589579, 0.589578, 0.000e+00, "PARTIAL", true},
      {6, "MS", 3.784115e-16, 0.221261, 0.220803, 0.000e+00, "PARTIAL", true},
      {6, "MA", 1.535452e-17, 0.452100, 0.452099, 0.000e+00, "PARTIAL", true},
      {7, "SA", 5.513645e-16, 0.579312, 0.578408, 0.000e+00, "PARTIAL", true},
      {7, "MS", 3.932189e-16, 0.237672, 0.237274, 0.000e+00, "PARTIAL", true},
      {7, "MA", 1.547279e-17, 0.440696, 0.440695, 0.000e+00, "PARTIAL", true},
      {8, "SA", 5.679854e-16, 0.589433, 0.589433, 0.000e+00, "PARTIAL", true},
      {8, "MS", 4.080793e-16, 0.232483, 0.232072, 0.000e+00, "PARTIAL", true},
      {8, "MA", 1.576534e-17, 0.428096, 0.428096, 0.000e+00, "PARTIAL", true},
      {9, "SA", 2.281645e-16, 0.006136, 0.005601, 0.000e+00, "OK", false},
      {9, "MS", 2.139152e-16, 0.005760, 0.004827, 0.000e+00, "OK", false},
      {9, "MA", 1.382392e-17, 0.002004, 0.001067, 0.000e+00, "OK", false},
      {10, "SA", 2.283482e-16, 0.006173, 0.005624, 0.000e+00, "OK", false},
      {10, "MS", 2.157045e-16, 0.004363, 0.003344, 0.000e+00, "OK", false},
      {10, "MA", 1.387104e-17, 0.001997, 0.001115, 0.000e+00, "OK", false},
      {11, "SA", 2.284281e-16, 0.006157, 0.005612, 0.000e+00, "OK", false},
      {11, "MS", 2.159420e-16, 0.004691, 0.003711, 0.000e+00, "OK", false},
      {11, "MA", 1.395154e-17, 0.001921, 0.001026, 0.000e+00, "OK", false},
  };
  REQUIRE(rows.size() == 33);
  int marked = 0;
  for (const auto &row : rows)
  {
    CoverageInterfaceSummary interface;
    interface.within = row.within;
    interface.partition_within = row.within;
    interface.served = row.within * (1.0 - row.uncovered_share);
    interface.missing = row.within * row.missing_share;
    interface.pending = row.within * row.pending_share;
    const auto cell =
        ComputeCoverageCell(row.type, {interface}, CoverageThresholds{},
                            CoverageValidatedRegime{}, CoverageExcitationClass::EIGENMODE);
    INFO("mode " << row.mode << " " << row.type);
    CHECK(std::string(CoverageStatusName(cell.status)) == row.status);
    const bool cell_marked = IsCoverageStatusMarked(cell.status);
    CHECK(cell_marked == row.marked);
    CHECK(cell_marked == (row.mode >= 2 && row.mode <= 8));
    marked += cell_marked;
  }
  CHECK(marked == 21);
}

}  // namespace palace
