// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_RESPONSE_COVERAGE_HPP
#define PALACE_MODELS_SURFACE_RESPONSE_COVERAGE_HPP

#include <map>
#include <ostream>
#include <set>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include <nlohmann/json.hpp>
#include "models/surfacepostoperator.hpp"
#include "utils/units.hpp"

namespace palace
{

// The coverage outputs of the surface-response correction (lane L1, decision 1634; the
// coverage design DESIGN-v2 s1 / s4 / s6.2; rulings 1641 / 1649 / 1654 (III)). Pure
// bookkeeping on the closed coverage partition of one field
// (SurfacePostOperator::GetInterfaceCoveragePartition): nothing here changes a corrected
// value, and no corrected value reads it.
//
// W = the RAW within-R surface energy of an interface (E_surf raw - E_out at the matching
// radius R). U = the part of W where NO correction is applied in the cell (fail-closed,
// ruling 1641: U = W - served). U = missing + method gap: missing = the energy on the
// uncovered features whose key the identification reports Missing (a library could hold
// a model for it); method gap = everything else left raw (the structural classes, the
// Unboxable and placement-rejected keys, the unclaimed remainder): a LOWER BOUND of what
// no library can cover today. The pending share = the served energy on labels holding a
// model admitted with a non-Qualified QualificationStatus (or LibraryQualified false and
// no status), counted whole; the no-status share = the served energy on labels holding
// a model whose record carries no QualificationStatus (information only).

// The coverage status of a cell on U + pending, with the config's interim bars (L1b,
// pre-registered: PREREG-L1B). NO-ENERGY: W = 0 (nothing to cover: the shares print
// nan). NOT-MEASURED: no partition exists for the cell. NOT-TARGET: an interface in no
// Type (per-interface rows only). Every one but OK is a non-OK status of the banner.
enum class CoverageStatus : char
{
  OK,
  LOW,
  PARTIAL,
  NOT_ACHIEVABLE,
  NO_ENERGY,
  NOT_MEASURED,
  NOT_TARGET,
  // Ruling 1714 (B)(3): the signed method gap is below -kCoverageInconsistentShare x W: the
  // partition serves MORE than the corrected total removes (a defect to report; never OK).
  INCONSISTENT,
  // Ruling 1738 (OI-CONFIGURED-TARGET-NOT-SERVED): a Type holding an interface the
  // configuration targets (TargetInterfaces, or every typed interface with EdgeDistances
  // when the list is absent) that no loaded model has a fabricated surface on: that
  // interface is not corrected at all (its corrected value is its raw value); worse than
  // any band, never OK.
  TARGET_NOT_SERVED
};
// Ruling 1714 (B)(3): the round-off margin of the signed method-gap share.
inline constexpr double kCoverageInconsistentShare = 1.0e-9;

// The regime status of a cell (DESIGN-v2 s6.2; decision 1636 (1)): whether the covered
// part's correction factor lies in the range validated for this library and excitation
// class. A status, never a correction: no value depends on it.
enum class RegimeStatus : char
{
  IN,
  OUT,
  NO_RANGE,
  NO_MODEL
};

enum class CoverageExcitationClass : char
{
  ELECTROSTATIC_TERMINAL,
  EIGENMODE,
  DRIVEN,
  NONE  // a solve type with no validated range class (boundary mode): NO-RANGE
};

// How a served model's library record states its qualification (item (5); note 1642).
enum class CoverageModelQualification : char
{
  QUALIFIED,  // QualificationStatus "Qualified"
  PENDING,    // a non-Qualified status admitted by AllowPendingModels, or LibraryQualified
              // false with no status
  NO_STATUS   // no QualificationStatus field
};

std::string_view CoverageStatusName(CoverageStatus status);
std::string_view RegimeStatusName(RegimeStatus status);
std::string_view CoverageExcitationClassName(CoverageExcitationClass excitation);

// The order of "the worse status" (a matrix entry, the confidence row's worst over
// Types): OK < LOW < PARTIAL < NOT-ACHIEVABLE < NO-ENERGY < NOT-MEASURED (a cell that
// says nothing about its coverage counts as the worst); NOT-TARGET never enters.
// Regime: IN < NO-MODEL < NO-RANGE < OUT.
int CoverageStatusRank(CoverageStatus status);
int RegimeStatusRank(RegimeStatus status);

// Lane L4 (decision 1829; USER 1636 (2); DESIGN-L4-v2 s2.1): THE MARKED SET of the
// mark-not-applied policy, one table (IsCoverageStatusMarked): a corrected value built from
// a cell whose coverage status is marked prints nan under Policy "MarkNotApplied" and its
// value moves to the appended partial column. PARTIAL and NOT-ACHIEVABLE (USER 1636 (2)),
// INCONSISTENT, TARGET-NOT-SERVED and NOT-MEASURED (question Q2 to the USER, the
// recommended answer "yes": this table is the switch) are marked; OK, LOW, NO-ENERGY (W =
// 0: no within-R raw energy near identified edges to correct) and NOT-TARGET are not. A
// regime status never marks.
bool IsCoverageStatusMarked(CoverageStatus status);
// Lane L4 (question Q4, the recommended name): the header suffix of the appended column
// that holds a marked corrected value (the covered part corrected PLUS the uncovered
// within-R part raw: today's number). One constant.
inline constexpr const char *kCoveragePartialColumnSuffix = " partial (uncovered part raw)";
// Lane L4 (question Q3, the recommended answer "keep + flag" = false): whether a marked
// cell marks the per-interface corrected values of EVERY Type of its measurement (true) or
// only of its own Type (false: a covered Type keeps its p / Q, flagged by the nan corrected
// normalization of its row and the exposure column). One constant.
inline constexpr bool kCoverageMarkWholeMeasurement = false;
// Lane L4 round 2 (decision 1914): the EXPLICIT per-cell marks. Under Policy
// "MarkNotApplied" every corrected table appends, LAST, integer 0 / 1 columns = the
// writer's own mark of the cells it covers (1 = the cells print nan, their values in the
// partial columns), so that a reader keys on a stated mark, never on a token; under "Warn"
// none exists. The readers and the guard match these headers EXACTLY (docs/src/reference.md
// lists them). Maxwell surface-Q-corrected.csv: the E_norm corrected / corrected fixed-flux
// / self-consistent cells of the row's measurement, and the row interface's E / p / Q
// corrected, corrected fixed-flux and self-consistent cells; electrostatic
// surface-Q-corrected.csv: the source's E_elec postprocessed fixed-trace / fixed-flux /
// corrected cells, and per interface k its E / p / Q postprocessed fixed-trace / fixed-flux
// / corrected cells; surface-response-confidence.csv: the row's measurement cells; the
// corrected terminal matrices: entry (i, j), fixed-trace and corrected.
inline constexpr const char *kCoverageMarkNormalizationColumn = "coverage mark E_norm";
inline constexpr const char *kCoverageMarkInterfaceColumn = "coverage mark E_surf";
inline constexpr const char *kCoverageMarkSourceColumn = "coverage mark E_elec";
inline constexpr const char *kCoverageMarkIndexedInterfaceColumn = "coverage mark[{}]";
inline constexpr const char *kCoverageMarkMeasurementColumn = "coverage mark measurement";
inline constexpr const char *kCoverageMarkMatrixEntryColumn = "coverage mark[i][{}]";

// The interim bars (SurfaceResponseCorrection.Coverage; the gate in force's
// BuildTarget 0.02 and MaxUncoveredShare 0.05 as placeholders).
struct CoverageThresholds
{
  double build_target = 0.02;
  double max_uncovered_share = 0.05;
};

// The validated range of the covered correction factor per Type and excitation class (a
// sidecar file named in the config, SurfaceResponseCorrection.Coverage.ValidatedRegime;
// read by ReadCoverageValidatedRegime). Empty = every cell NO-RANGE.
struct CoverageValidatedRange
{
  double factor_min = 0.0, factor_max = 0.0;
  // The number of validation windows behind the range (>= 1) and their provenance (the
  // file's Provenance entry as JSON text); review m-1: both required.
  int n = 0;
  std::string provenance;
};
// The largest regime tolerance a ValidatedRegime file may carry (review m-1: a sanity bound
// on a data file, not a validated value; tau itself is set by the validation track).
inline constexpr double kCoverageRegimeToleranceMax = 0.1;
struct CoverageValidatedRegime
{
  std::string path, library_name, library_sha256;
  double tolerance = 0.0;
  std::map<std::pair<std::string, CoverageExcitationClass>, CoverageValidatedRange> ranges;
};

// Parses the validated-regime file and checks it against the loaded library: its
// Library.Sha256 must equal the SHA-256 of the library file's bytes (a range of another
// library would print a wrong status: a config error).
CoverageValidatedRegime ReadCoverageValidatedRegime(const std::string &path,
                                                    const std::string &library_path);

// The qualification of every model of a library file by Name, read from the model records
// top-level fields (item (5); note 1642): QualificationStatus "Qualified" = QUALIFIED,
// another QualificationStatus = PENDING (admitted by AllowPendingModels at load), none with
// LibraryQualified false = PENDING, none otherwise = NO_STATUS. A read only: the library's
// load (ReadProcessLibrary) is unchanged.
std::map<std::string, CoverageModelQualification>
ReadCoverageModelQualifications(const std::string &library_path);

// The classification of an uncovered feature's key record (ruling 1649 Q2).
struct CoverageKeyRecord
{
  std::string hash;
  std::string kind;    // "Missing" | "Unboxable" | "PlacementRejected"
  std::string reason;  // the placement rejection's class, when any
  bool method_gap = false;
};
CoverageKeyRecord ClassifyCoverageKey(const std::string &hash, bool matched,
                                      const std::string &placement_rejection,
                                      bool unboxable);

// Ruling 2277 (lane dropcover): the NEAR-EDGE ENERGY BUDGET of an interface / cell: every
// EXCLUSIVE bin of the coverage partition (J), each within-R point once by the partition's
// first-list rule (the kept classes' own tables integrate each class alone and can count a
// point in two windows). `served` = by a claim (after the unpatched and loop-remainder
// moves); `context` / `context_overlap` = held by one / two or more matched boxes as
// context (ruling 2276); `missing` = the Missing class on keys a library could hold (the
// coverage's missing part), `missing_method_gap` on Unboxable / placement-rejected keys;
// `loop_remainder` = the part the partition moved out of the served claims;
// `inside_excluded_box` = a dropped window's points that only a box excluded at this
// excitation contains (ruling 2295 (A): such a box holds nothing; dropped);
// `kept_inside_box` (the mirror) is a PART of `spatial_cluster_priority`, not a bin, and
// `context_pending` / `context_no_status` the PARTS of the two context bins that the
// coverage counts pending / no-status (by the holders' models, labels counted whole).
struct CoverageNearEdgeBudget
{
  double served = 0.0, context = 0.0, context_overlap = 0.0;
  double missing = 0.0, missing_method_gap = 0.0, domain_boundary = 0.0,
         spatial_cluster_priority = 0.0, conductor_consistency = 0.0, loop_remainder = 0.0;
  double unpatched_pairs = 0.0, excluded = 0.0, beyond_box = 0.0, unclaimed_claim_end = 0.0,
         unclaimed_no_claims = 0.0, unclaimed_elsewhere = 0.0, inside_excluded_box = 0.0;
  double kept_inside_box = 0.0, context_pending = 0.0, context_no_status = 0.0;
  // Per holder label ("m4+m5"), the energy held in a box overlap.
  std::map<std::string, double> overlaps;
  // The sum of the exclusive bins (every member but the parts and overlaps).
  double Sum() const;
  // The energy the corrected total neither replaces by a model nor keeps raw.
  double Dropped() const;
  void Add(const CoverageNearEdgeBudget &other);
};
// One interface of one cell: W and the partition's sums (J).
struct CoverageInterfaceSummary
{
  double within = 0.0;            // W (raw - outside at R): the corrected bookkeeping's
  double partition_within = 0.0;  // the partition pass's own within-R total
  double served = 0.0, missing = 0.0, pending = 0.0, no_status = 0.0;
  double fabricated = 0.0;  // the correction's fabricated energy (evaluation 0)
  // Ruling 1689 M-2: the loop remainder no served label held (the partition's
  // loop_remainder_unserved). Ruling 1694 (i): the loop remainder of the conductor-
  // consistency-EXCLUDED patches by its own integral (OI-CORRECTED-TOTAL-CC-LOOPREMAINDER-
  // OVERLAP: the corrected total keeps it raw twice, in the ConductorConsistency term and
  // in the loop-remainder term). Maxwell only (0 on the electrostatic path).
  double loop_remainder_unserved = 0.0, loop_remainder_excluded = 0.0;
  // Rulings 1714 / 1715: the partition's end census (served claims; kept classes).
  SurfacePostOperator::InterfaceCoveragePartition::EndCensus served_ends, kept_ends;
  // Lane ccpx (ruling 2229 D4): the ConductorConsistency kept class's energy (inside U and
  // the method gap, which keep their definitions; printed only under PerExcitation).
  double conductor_consistency = 0.0;
  // Ruling 2277: the near-edge energy budget (`served` above = budget.served +
  // budget.context + budget.context_overlap: ruling 2276 counts the context as served).
  CoverageNearEdgeBudget budget;
};
CoverageInterfaceSummary
SummarizeCoverageInterface(const SurfacePostOperator::InterfaceCoveragePartition &partition,
                           double within, double fabricated,
                           const std::map<int, CoverageKeyRecord> &keys_by_feature,
                           const std::map<int, CoverageModelQualification> &models);

// One cell x Type (the rows of surface-response-coverage.csv). Energies in J; shares are
// fractions of W (nan when the status is NO-ENERGY or NOT-MEASURED); the factor is nan
// exactly when the regime is NO-MODEL, the validated min / max exactly when it is
// NO-RANGE or NO-MODEL (rulings 1641 (4) / 1654 (III) Q1).
struct CoverageCell
{
  std::string type;
  double within = 0.0, uncovered = 0.0, missing = 0.0, method_gap = 0.0, pending = 0.0,
         no_status = 0.0, served = 0.0, fabricated = 0.0;
  // The partition pass's own within-R total (its difference from W is the round-off of the
  // two integrals; review m-3: printed, never used for a status).
  double partition_within = 0.0;
  double uncovered_share = 0.0, missing_share = 0.0, method_gap_share = 0.0,
         pending_share = 0.0, no_status_share = 0.0;
  // (U - missing) / W BEFORE the clamp at 0 (decision 1708 (1)): printed, never a status.
  double method_gap_signed_share = 0.0;
  double factor = 0.0, validated_min = 0.0, validated_max = 0.0;
  // The validated range's window count (0 unless the regime is IN or OUT).
  int validated_n = 0;
  CoverageStatus status = CoverageStatus::NOT_MEASURED;
  RegimeStatus regime = RegimeStatus::NO_MODEL;
  // Ruling 1694 (1): the raw interface energy (whole interface, J) of this Type's
  // UNTARGETED postprocessing interfaces (they carry their raw energy, not corrected) and
  // its share of the Type's total raw interface energy (targeted + untargeted; nan when
  // that total is 0). The shares above are over the TARGETED interfaces only. A Type with
  // untargeted interfaces only is a NOT-TARGET cell (W = 0, shares nan, NO-MODEL).
  double untargeted_raw = 0.0, untargeted_share = 0.0;
  // The Type's sums of CoverageInterfaceSummary's two loop-remainder terms (J).
  double loop_remainder_unserved = 0.0, loop_remainder_excluded = 0.0;
  // The Type's sums of the end census (counts and energies summed, distances the largest).
  SurfacePostOperator::InterfaceCoveragePartition::EndCensus served_ends, kept_ends;
  // Lane ccpx (ruling 2229 D4): the Type's ConductorConsistency kept-class energy (J) and
  // its share of W (nan where the other shares are): the part of U (inside the method gap)
  // held by the patches excluded at this excitation.
  double conductor_consistency = 0.0, conductor_consistency_share = 0.0;
  // Ruling 2277: the Type's near-edge energy budget (the interfaces' sums, J).
  CoverageNearEdgeBudget budget;
  // The edge-resolution census (decision 1817) is not a column here: it writes its own
  // table, block, PREFLIGHT line and banner (models/surfaceresponseresolution.hpp).
};
CoverageCell ComputeCoverageCell(const std::string &type,
                                 const std::vector<CoverageInterfaceSummary> &interfaces,
                                 const CoverageThresholds &thresholds,
                                 const CoverageValidatedRegime &regime,
                                 CoverageExcitationClass excitation);
// A cell with no partition (Q2 (b) of ruling 1654 (III)): NOT-MEASURED, no shares.
CoverageCell NotMeasuredCoverageCell(const std::string &type);

// The rows of surface-response-uncovered-energy-by-key.csv: per key (the hash of a
// Missing-class bin's features, '+'-joined in ascending order when one bin joins several
// keys: such a row cannot be attributed to one key), its kind / reason / method-gap flag
// (a joint row is method gap only when every token is) and the raw within-R energy per
// interface.
struct CoverageKeyRow
{
  std::string key, kind, reason;
  bool method_gap = false;
  std::vector<int> features;
  std::map<int, double> energy;
};
std::vector<CoverageKeyRow>
CoverageKeyRows(const std::map<int, SurfacePostOperator::InterfaceCoveragePartition> &cell,
                const std::map<int, CoverageKeyRecord> &keys_by_feature);

// The coverage of ONE measurement (a Maxwell (m | f, exc), an electrostatic source): per
// target interface its summary and its Type, the cells per Type (in order of first
// appearance over the ascending target interfaces) and the key rows. Computed in the
// measurement step BEFORE any corrected column or matrix is written; every writer reads
// the cell status from it (the hook of a coverage policy: lane L4's MarkNotApplied reads
// `cells` there, the accumulators never read a policy). `measured` false = the cells are
// NOT-MEASURED (no partition: ruling 1654 (III) Q2 (b)).
struct CoverageMeasurement
{
  bool measured = false;
  CoverageExcitationClass excitation = CoverageExcitationClass::NONE;
  std::map<int, CoverageInterfaceSummary> interfaces;
  std::map<int, std::string> types;
  std::vector<CoverageCell> cells;
  // CUT-PLAN step 4 (lane v1table): per target interface a cell of that interface ALONE
  // (its own summary; the surface-participation.csv verdicts of a Region target's row). No
  // other writer and no policy reads it.
  std::map<int, CoverageCell> interface_cells;
  std::vector<CoverageKeyRow> keys;
  // The served models (runtime index) with served energy whose library record carries no
  // QualificationStatus (item (5): information) and the pending ones.
  std::set<int> served_no_status_models, served_pending_models;
  // Ruling 1694 (1): the untargeted postprocessing interfaces and their Types
  // (AddCoverageUntargetedInterfaces).
  std::map<int, std::string> untargeted;
  // Ruling 1738: the configured targets no loaded model serves and their Types (in `types`
  // too: their within-R raw energy is in their Type's W, uncovered;
  // MarkCoverageTargetsNotServed).
  std::map<int, std::string> not_served;
  const CoverageCell *FindCell(const std::string &type) const;
  // The worst cell over Types (CoverageStatusRank / RegimeStatusRank); NOT-TARGET cells
  // never enter; nullptr when no other cell exists.
  const CoverageCell *WorstCoverageCell() const;
  const CoverageCell *WorstRegimeCell() const;
  // Lane L4 (DESIGN-L4-v2 s2.3): whether any cell of the measurement is marked
  // (IsCoverageStatusMarked; the whole-measurement values: the corrected normalizations and
  // totals, the frequency shift, the self-consistent frequency, the terminal-matrix entries
  // of a source); whether the per-interface corrected values of a Type are marked (its own
  // cell, or AnyMarked under kCoverageMarkWholeMeasurement); and the raw within-R energy
  // the marked cells leave uncorrected (the sum of their U W; nan when a marked cell has no
  // measured U, e.g. NOT-MEASURED; 0 when no cell is marked).
  bool AnyMarked() const;
  bool IsTypeMarked(const std::string &type) const;
  double MarkedUncoveredEnergy() const;
  // The mark of the measurement's whole values: AnyMarked, or no cell at all (fail-closed:
  // a corrected measurement without cells says nothing about its coverage).
  bool MarksWholeMeasurement() const { return AnyMarked() || cells.empty(); }
  // Lane L4 (ERRATUM-1847 m-2): the Types of the coverage interfaces (`types`) that hold no
  // cell; the writers mark such an interface's corrected values (fail-closed, like a
  // measurement without cells).
  std::vector<std::string> TypesWithoutCell() const;
};
// Lane L4 (ERRATUM-1829 E-1; DESIGN-L4-v2 s2.3 E2): the marks of the corrected terminal
// matrices from the per-source marks (MarksWholeMeasurement of each source's FINAL
// measurement): C(i, j) iff source i or source j is marked (never the coverage table's
// worse-of-two status: NO-ENERGY / TARGET-NOT-SERVED / NOT-MEASURED outrank PARTIAL there);
// C_m off-diagonal as C(i, j), diagonal iff any entry of row i of C is; C⁻¹ every entry iff
// any C entry is.
struct CoverageTerminalMatrixMarks
{
  std::vector<std::vector<bool>> capacitance, mutual, inverse;
};
CoverageTerminalMatrixMarks
ComputeCoverageTerminalMatrixMarks(const std::vector<bool> &source_marked);
// `within` (W) and `fabricated` per target interface (the corrected bookkeeping's raw -
// outside at R and the correction's fabricated energy, evaluation 0); `types` the target
// interfaces' Types; `loop_remainder_excluded` per target interface (J; Maxwell; absent =
// 0).
CoverageMeasurement ComputeCoverageMeasurement(
    const std::map<int, SurfacePostOperator::InterfaceCoveragePartition> &partitions,
    const std::map<int, std::string> &types, const std::map<int, double> &within,
    const std::map<int, double> &fabricated,
    const std::map<int, CoverageKeyRecord> &keys_by_feature,
    const std::map<int, CoverageModelQualification> &models,
    const CoverageThresholds &thresholds, const CoverageValidatedRegime &regime,
    CoverageExcitationClass excitation,
    const std::map<int, double> &loop_remainder_excluded = {});
CoverageMeasurement NotMeasuredCoverageMeasurement(const std::map<int, std::string> &types,
                                                   CoverageExcitationClass excitation);
// Ruling 1694 (1): records the untargeted postprocessing interfaces (`untargeted`: index ->
// Type) in the measurement and, per cell, its Type's untargeted raw energy and share from
// `raw_energy` (the whole-interface raw energy of every postprocessing interface, J); a
// Type with untargeted interfaces only gets a NOT-TARGET cell (W = 0, shares / factor /
// range nan, NO-MODEL), after the target Types' cells.
void AddCoverageUntargetedInterfaces(CoverageMeasurement &measurement,
                                     const std::map<int, std::string> &untargeted,
                                     const std::map<int, double> &raw_energy);
// Ruling 1738: records the configured targets no loaded model serves (`not_served`: index
// -> Type; each must be in the measurement's `types`) and gives every cell of their Types
// the worse of its status and TARGET-NOT-SERVED.
void MarkCoverageTargetsNotServed(CoverageMeasurement &measurement,
                                  const std::map<int, std::string> &not_served);
// Ruling 1738: the reason a configured target is not served (the startup warning, the
// banner and palace.json's ConfiguredTargetsNotServed).
inline constexpr const char *kCoverageTargetNotServedReason =
    "no loaded model has a fabricated surface on it";
// Ruling 1738: the startup warning naming each configured target no loaded model serves
// (empty when there is none).
std::string CoverageTargetsNotServedWarning(const std::map<int, std::string> &not_served);

struct CoverageFeaturePortions
{
  int portions = 0;
  double length = 0.0;
  bool mirror_formed = false;
};
using CoverageFeaturePortionsMap = std::map<int, CoverageFeaturePortions>;

// The plain-text rows of surface-response-coverage.csv and
// surface-response-uncovered-energy-by-key.csv for one measurement (both solve families):
// `keys` = the row's leading cells ("m,exc" | "f (GHz),exc" | "i" values, comma-joined),
// `key_header` their header; energies in J and lengths in m through `units`; nan prints
// "nan". `portions_by_feature`: CoveragePortionsByFeature.
std::string CoverageTableHeader(const std::string &key_header);
std::string CoverageByKeyTableHeader(const std::string &key_header,
                                     const std::set<int> &targets);
void WriteCoverageRows(std::ostream &output, const std::string &keys,
                       const CoverageMeasurement &measurement, const Units &units);
// Ruling 2277: the near-edge energy budget table (one row per coverage cell, every
// exclusive bin of the partition as a raw share of W, the closure W - the bins; beside
// surface-response-coverage.csv, on every corrected run).
inline constexpr const char *kNearEdgeBudgetTableName =
    "surface-response-near-edge-energy-budget.csv";
std::string NearEdgeBudgetTableHeader(const std::string &key_header);
void WriteNearEdgeBudgetRows(std::ostream &output, const std::string &keys,
                             const CoverageMeasurement &measurement, const Units &units);
void WriteCoverageByKeyRows(std::ostream &output, const std::string &keys,
                            const CoverageMeasurement &measurement,
                            const std::set<int> &targets,
                            const CoverageFeaturePortionsMap &portions_by_feature,
                            const Units &units);
// Per uncovered feature: its portions' count and summed length (the by-key rows' columns)
// and whether the identification formed it by the mirror rule (its portions' type carries
// ":MirrorFormed"; ruling 1649 Q2: such a key stays class Missing, the reason printed).
CoverageFeaturePortionsMap CoveragePortionsByFeature(
    const std::vector<SurfacePostOperator::UncoveredPerimeterPortion> &uncovered);

// The run's coverage record (root only): the cells of every measurement, for the end-of-run
// banner and palace.json SurfaceResponse.Coverage (BaseSolver::SaveMetadata). `keys` are
// the cell's table keys ("m" | "f (GHz)", "exc" | "i").
struct CoverageRunRecord
{
  struct Entry
  {
    std::vector<std::pair<std::string, double>> keys;
    CoverageCell cell;
  };
  std::vector<Entry> entries;
  std::set<int> served_no_status_models, served_pending_models;
  std::map<int, std::string> untargeted;  // ruling 1694 (1): named once in the banner
  std::map<int, std::string>
      not_served;  // ruling 1738: named in the banner and palace.json
  // Ruling 1738: why no cell was measured (a raw fallback: its geometry fallback class);
  // with no cell, palace.json writes "Complete": true only when this is set.
  std::string no_cells_reason;
  // Lane L4 (ERRATUM-1847 m-2): the measurements marked without a marked-status cell (no
  // cell at all, or a coverage interface Type without a cell) and why; counted in
  // MarkedMeasurements and listed as MarkedWithoutCell under MarkNotApplied.
  struct MarkWithoutCell
  {
    std::vector<std::pair<std::string, double>> keys;
    std::string reason;
  };
  std::vector<MarkWithoutCell> marked_without_cell;
  // Lane ccpx (ruling 2229 D4): Exclusion PerExcitation writes each flagged cell's
  // conductor-consistency raw share of W (palace.json and the banner); false: as before.
  bool conductor_consistency_share = false;
  void Add(const std::vector<std::pair<std::string, double>> &keys,
           const CoverageMeasurement &measurement);
};
// palace.json SurfaceResponse.Coverage: the policy, the interim bars, the regime file, the
// status counts and one entry per cell whose coverage status is not OK or whose regime is
// not IN (erratum 1784 E-2: an OK + NO-RANGE cell is flagged as an OK + OUT cell is;
// energies in J by `energy_scale`; nan as null). Erratum 1784 E-3: `level` of
// `levels` = the solve of the AMR sequence that wrote the block (1 = the initial mesh) of
// Refinement.MaxIts + 1 (1 without AMR), written as Level / Levels; the block never holds
// "Final" (BaseSolver adds it once the AMR loop has ended).
// Lane L4 (DESIGN-L4-v2 s2.4): `mark_not_applied` (Policy "MarkNotApplied") writes "Policy"
// "MarkNotApplied", "MarkedStatuses", "MarkedCells", "MarkedMeasurements" (a measurement
// with a marked cell or marked without one: ERRATUM-1847 m-2), "MarkedWithoutCell" (when
// any) and each Flagged entry's "Marked"; false (Policy "Warn") writes the block exactly as
// before.
nlohmann::json CoverageRunMetadata(const CoverageRunRecord &record,
                                   const CoverageThresholds &thresholds,
                                   const CoverageValidatedRegime &regime,
                                   double energy_scale, int level, int levels,
                                   bool mark_not_applied = false);
// The end-of-run banner of palace.json's block (a line per flagged cell: coverage not OK or
// regime not IN, a NO-RANGE cell's line saying "computed, NOT VALIDATED in this regime";
// at most 40 such lines, then one line counting the rest, which palace.json's Flagged
// lists; one NO-RANGE summary line; one line naming the status-less served models). Lane
// L4: a block with "Policy" "MarkNotApplied" adds the MARK-NOT-APPLIED header and label
// lines after the first line, ends each flagged line of a marked cell with " [NOT
// APPLIED]" and adds the share-denominator line; a "Warn" block prints the banner exactly
// as before.
std::string CoverageRunBanner(const nlohmann::json &metadata);
// Ruling 2276 R7 (b): the run's near-edge warnings (each line once; none: empty): energy
// held in a box overlap (its holders; the fabricated energy there may be counted by each
// holder, not measured) and energy the corrected total neither models nor keeps raw.
std::vector<std::string> CoverageNearEdgeWarnings(const CoverageRunRecord &record);

// Ruling 1688 (the stage-0 Missing LENGTH line; a log line only): the uncovered portions
// split by their feature's key class: first = Missing (and a feature without a key record),
// second = method gap (PlacementRejected, Unboxable).
std::pair<std::vector<SurfacePostOperator::UncoveredPerimeterPortion>,
          std::vector<SurfacePostOperator::UncoveredPerimeterPortion>>
SplitCoverageUncoveredByKeyClass(
    const std::vector<SurfacePostOperator::UncoveredPerimeterPortion> &uncovered,
    const std::map<int, CoverageKeyRecord> &keys_by_feature);
// The line: per interface Type (of `types`, every interface) the Missing length over the
// edge-tree length of the Type's TARGET interfaces (`lengths` per target: lists = {Missing,
// method gap}), the method-gap length beside it; "n/a" for a Type without a target
// interface or with a zero length; then the identification's all-feature figure
// unmatched / (matched + unmatched) ("n/a" when no feature pass ran in this process).
std::string CoverageStageZeroLine(
    const std::map<int, std::string> &types, const std::set<int> &targets,
    const std::map<int, SurfacePostOperator::InterfacePerimeterLengths> &lengths,
    std::pair<double, double> identification_lengths, const Units &units);

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_RESPONSE_COVERAGE_HPP
