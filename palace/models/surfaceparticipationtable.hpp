// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_PARTICIPATION_TABLE_HPP
#define PALACE_MODELS_SURFACE_PARTICIPATION_TABLE_HPP

#include <array>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <nlohmann/json_fwd.hpp>
#include "utils/communication.hpp"
#include "utils/labels.hpp"
#include "utils/tablecsv.hpp"

namespace palace
{

struct CoverageCell;
struct ResolutionCell;
struct ResolutionPreflightRow;

// CUT-PLAN step 4 (decisions 2378 / 2379; lane v1table DESIGN-NOTE): the designer's table
// surface-participation.csv. One row per (adaptive iteration, point, excitation, interface)
// with the raw and the fixed-trace corrected participation, the row's verdicts and the
// change since the previous iteration; one file per run at the top level, appended at the
// end of every converged level. Output only: no value, solve or marking reads it.
inline constexpr const char *kSurfaceParticipationFile = "surface-participation.csv";

// The words of the impedance-metal rule (CUT-PLAN step 3; SurfaceParticipationRow
// sheet_served / sheet_not_served): a value with Ls sheets served by the perfect-conductor
// models (an information word, lane v1wordscode), and a value raw at a metal the rule does
// not serve (an Ls sheet above the bound, a sheet with Rs or Cs, another law; a flags
// word).
inline constexpr const char *kParticipationSheetServedFlag = "Ls sheet served as PEC";
inline constexpr const char *kParticipationSheetNotServedFlag =
    "metal law not served: raw there";

// The flags word of a value that a mismatch between the configuration and its library
// affects and the run tolerates (surfaceresponselibrarycheck.hpp; SurfaceParticipationRow
// library_mismatch): the mismatches, their sizes and estimated effects are in palace.json
// SurfaceResponse.LibraryMismatch.
inline constexpr const char *kParticipationLibraryMismatchFlag =
    "library mismatch tolerated";
// Whether a tolerated mismatch acts on a row: an interface row of a listed target
// interface, a region row of a listed Type (a Type total: through its members).
inline bool LibraryMismatchActsOn(const std::set<int> &interfaces,
                                  const std::set<std::string> &types, int interface,
                                  const std::string &type, bool region)
{
  return interfaces.count(interface) > 0 || (region && types.count(type) > 0);
}

// The flag S (lane v1verdict1; ROBUSTNESS-REGISTER G19; bench-heldout round 4, DESIGNER
// RULE): S = the share of a value's within-R corrected energy (its census cell's
// fixed-trace fabricated energy) served by spatial-basis models, corner and spatial-cluster
// models, whose trace map onto 16-knot rings is not within 1 %. At S >= 5 % SA can be off
// by up to 5 S points and MS by up to 6 S: the first word; more than 10 % on
// spatial-cluster models: the second. Information words (the column `information`, lane
// v1wordscode): no value changes and no summary word reads them.
inline constexpr double kParticipationSpatialShareFlagAt = 0.05;
inline constexpr double kParticipationClusterShareFlagAbove = 0.10;
inline constexpr const char *kParticipationSpatialShareFlag =
    "S >= 5% on corner and cluster models";
inline constexpr const char *kParticipationClusterShareFlag =
    "S over 10% on cluster models";

// The Maxwell gates that concern the fixed-trace (v1) value (SUPERVISOR-RULINGS R2, lane
// v1wordscode): the operator's closure-independent gate reads these limits
// (surfaceresponseoperator.cpp) and the table names a failed one in `flags` with its number
// (the loop residual is the response-weighted one); R/rho beyond its limit is an
// information word (not validated as a verdict). An MA value is never certified in v1.
inline constexpr double kMaxwellGateMaximumKR = 0.1;
inline constexpr double kMaxwellGateMaximumWeightedLoopResidual = 0.05;
inline constexpr double kMaxwellGateMaximumLoopResponseFailureFraction = 0.01;
inline constexpr double kMaxwellGateMaximumLibraryDistance = 0.8;
inline constexpr double kMaxwellGateMaximumCurvatureRatio = 0.25;
inline constexpr const char *kParticipationBoundaryLawFlag = "boundary law not verified";
inline constexpr const char *kParticipationMaNotCertified = "MA not certified";

// The summary word of a row (column `summary`; SUPERVISOR-RULINGS R5, lane v1wordscode).
inline constexpr const char *kParticipationSummaryOk = "OK";
inline constexpr const char *kParticipationSummaryCheck = "CHECK";
inline constexpr const char *kParticipationSummaryDoNotUse = "DO-NOT-USE";

// The verdict data (lane meshrule VERDICT-INPUT.md and its check; DATA, not code: the
// built-in values are the repository file palace/models/surfaceparticipationverdict.json,
// embedded at build time, and the optional "Verdict" member of a Resolution.BarsFile file
// replaces them). G = the growth of the unknowns from an earlier level for the change to be
// read; F = the largest ratio m90_k / m90_j of the value's own cell census for the change
// to be read (its edges were refined since j); T = the largest change in % per doubling of
// the unknowns that reads SETTLED; B = per FE order 1-4 the largest cell census m90 (h / R)
// of a RESOLVED value (NaN = none; orders >= 5 have none); a census above
// `under_resolved_factor` x B reads UNDER-RESOLVED. `tested_problems`: the problem types
// the rule was tested on; the words of any other read -UNTESTED. `tested_scopes`: the
// values whose census the bars were graded with ("Type": a Type's cell census; "Region": a
// box's own patches); the resolution words of a region row read -UNTESTED unless "Region"
// is named (the bars were graded with the cell's census, which reads coarser than a
// junction box's own).
struct ParticipationVerdictBars
{
  double min_unknowns_ratio = 0.0;
  double max_census_ratio = 0.0;
  double max_change_per_doubling = 0.0;  // % per doubling
  std::array<double, 4> census_bar{};
  double under_resolved_factor = 0.0;
  std::vector<std::string> tested_problems, tested_scopes;
  bool built_in = true;
  bool interim = true;
  std::string source, sha256, provenance;
  // The bar at `order` (NaN at orders >= 5: no bar was measured there).
  double CensusBar(int order) const;
  // The suffix of the SETTLED / RESOLVED words: -INTERIM while the data are interim,
  // -UNTESTED for a problem type the rule was not tested on, else none.
  std::string WordSuffix(const std::string &problem) const;
  // Whether the bars were graded on the census of region (box) values.
  bool RegionScopeTested() const;
};
// The "Verdict" member of a JSON text: {"MinUnknownsRatio": G > 1, "MaxCensusRatio": F in
// (0, 1], "MaxChangePerDoublingPercent": T > 0, "CensusBarOverR": {"1": null | B > 0, ...,
// "4": ...}, "UnderResolvedFactor": >= 1, "TestedProblems": [names], "TestedScopes":
// ["Type"] or ["Type", "Region"], "Interim": bool
// (default false), "Provenance": text}.
ParticipationVerdictBars ParseParticipationVerdictBars(const std::string &text,
                                                       const std::string &source,
                                                       bool built_in);
ParticipationVerdictBars BuiltInParticipationVerdictBars();
// The "Verdict" member of a bars file; the built-in data when the file has no such member.
ParticipationVerdictBars ReadParticipationVerdictBars(const std::string &path);
// The verdict data of a run (its Resolution.BarsFile `bars_file`; built-in when empty): the
// one selection of the table and of the census operator's PREFLIGHT word.
ParticipationVerdictBars ParticipationVerdictBarsFor(const std::string &bars_file);
// The problem name of the table's rows and words (Electrostatic, Eigenmode, Driven; empty
// for a problem type without the table).
std::string ParticipationProblemName(ProblemType type);

// Why a value's indicator is not evaluated (printed after MOVING), the first that applies
// in the order of VERDICT-INPUT.md s2: no earlier level (FIRST-ITERATION: the first level,
// every fixed mesh), no earlier level with N_k >= G N_j (UNKNOWNS-GREW-<20% for G = 1.2),
// p_j* = 0 (ZERO-ANCHOR-VALUE), then this table's NO-VALUE (p_k or p_j* not finite), then
// the value's own cell census not fallen to F m90_j* or not read (CELL-NOT-REFINED: on a
// chip the unknowns can grow elsewhere while its edges stay); and UNMATCHED (an eigenmode,
// a driven frequency or an interface with no row at the previous level).
enum class ParticipationNotEvaluated : char
{
  NONE,
  FIRST_ITERATION,
  UNKNOWNS_GREW_LITTLE,
  UNMATCHED,
  ZERO_ANCHOR_VALUE,
  NO_VALUE,
  CELL_NOT_REFINED
};

// One earlier level of a value: its unknowns, its value and its cell census m90.
struct ParticipationHistoryEntry
{
  double unknowns = 0.0, p = 0.0, m90 = 0.0;
};

// The indicator of a value p with cell census m90 at a level of N unknowns, from its own
// earlier levels `history`, oldest first: j* = the LATEST with N >= G N_j (the guard
// compares with j* only); read only when p_j* != 0, p and p_j* are finite and m90 <= F
// m90_j* (exact comparisons); g = 100 |p / p_j* - 1| / log2(N / N_j*) (% per doubling);
// SETTLED iff g < T (strict). A region value's m90 is its own patches' census, in the guard
// as in the bar.
struct ParticipationIndicator
{
  double per_doubling = std::numeric_limits<double>::quiet_NaN();
  int anchor = -1;  // the index of j* in the history
  ParticipationNotEvaluated not_evaluated = ParticipationNotEvaluated::FIRST_ITERATION;
  bool settled = false;
};
ParticipationIndicator
ComputeParticipationIndicator(double p, double m90, double unknowns,
                              const std::vector<ParticipationHistoryEntry> &history,
                              const ParticipationVerdictBars &bars);
// SETTLED (+ WordSuffix), MOVING, or MOVING: <why not evaluated>.
std::string ParticipationConvergenceName(const ParticipationIndicator &indicator,
                                         const ParticipationVerdictBars &bars,
                                         const std::string &problem);
// The resolution word of a value from its cell census m90 (NaN: no reading) at `order`:
// RESOLVED (+ WordSuffix) iff m90 <= B AND the value is SETTLED, UNDER-RESOLVED iff
// m90 > factor x B, MARGINAL otherwise; NO-BARS without a B at the order (orders 1 and >=
// 5), NOT-EVALUATED without a reading. A single solve is never SETTLED, so never RESOLVED.
// A `region` value's RESOLVED, MARGINAL and UNDER-RESOLVED read -UNTESTED unless the bars
// name the "Region" scope.
std::string ParticipationResolutionWord(double m90, int order, bool settled, bool region,
                                        const ParticipationVerdictBars &bars,
                                        const std::string &problem);
// ONE verdict (decision 2403 (3), lane v1verdict1): the resolution word of a Type BEFORE a
// solve, from its census weighted by length (ResolutionPreflightRow::length_quantiles): the
// word above with no convergence evaluated (settled false, a Type's census), as the table
// prints it at a first iteration. Never RESOLVED.
std::string ParticipationPreflightWord(double length_m90, int order,
                                       const ParticipationVerdictBars &bars,
                                       const std::string &problem);

// The Maxwell gates of a row's fixed-trace estimate (R2 above; the numbers of
// SurfaceResponseOperator::MaxwellResponse). `evaluated` false on an electrostatic row (its
// path computes none of them) and on a row without the estimate: no word. A gate fails
// iff NOT (value <= limit), as the operator's own gate (a NaN fails).
struct ParticipationMaxwellGates
{
  bool evaluated = false;
  double kR = std::numeric_limits<double>::quiet_NaN();
  double loop_residual = std::numeric_limits<double>::quiet_NaN();  // response-weighted
  double loop_fraction = std::numeric_limits<double>::quiet_NaN();
  double library_distance = std::numeric_limits<double>::quiet_NaN();
  double curvature_ratio = std::numeric_limits<double>::quiet_NaN();  // max R/rho
  bool boundary_law_verified = true;
};
// The flags words of the failed gates, in the order kR, loop residual, loop fraction,
// library distance, boundary law ('kR 0.153 above 0.1'; empty when none failed or not
// evaluated), and the information word of R/rho beyond its limit (empty otherwise).
std::vector<std::string> ParticipationGateFlags(const ParticipationMaxwellGates &gates);
std::string ParticipationCurvatureInformation(const ParticipationMaxwellGates &gates);

// One row of surface-participation.csv (the DESIGN-NOTE's 27 columns; lane v1wordscode:
// regime and confidence are not written, summary and information are).
struct SurfaceParticipationRow
{
  int iteration = 0;
  long long int unknowns = 0;
  std::string problem;
  int point = 0;
  double frequency = std::numeric_limits<double>::quiet_NaN();  // GHz
  int excitation = 0;
  int interface = 0;
  std::string type;
  bool region = false;
  double p_raw = std::numeric_limits<double>::quiet_NaN();
  double p_corrected = std::numeric_limits<double>::quiet_NaN();
  double q_raw = std::numeric_limits<double>::quiet_NaN();
  double q_corrected = std::numeric_limits<double>::quiet_NaN();
  std::string resolution = "NOT-EVALUATED";
  double m90 = std::numeric_limits<double>::quiet_NaN();
  double m50 = std::numeric_limits<double>::quiet_NaN();
  std::string coverage = "NOT-MEASURED";
  double uncovered_share = std::numeric_limits<double>::quiet_NaN();
  // The coverage cell's regime word and the estimate's confidence word: not written to this
  // table and read by none of its words (R1, R3; lane v1wordscode); the developer tables
  // (surface-response-coverage.csv, surface-Q-corrected.csv) print their own.
  std::string regime = "NO-MODEL";
  std::string confidence = "none";
  double change_corrected = std::numeric_limits<double>::quiet_NaN();
  double change_raw = std::numeric_limits<double>::quiet_NaN();
  double per_doubling = std::numeric_limits<double>::quiet_NaN();
  std::string convergence;
  std::string flags;
  // The summary word (ParticipationSummary) and the information words
  // (ParticipationInformation) of the row.
  std::string summary;
  std::string information;
  std::string reason;
  // The flag S (a column) and the share on spatial-cluster models (a flags word only) of
  // the row's census cell.
  double spatial_share = std::numeric_limits<double>::quiet_NaN();
  double cluster_share = std::numeric_limits<double>::quiet_NaN();
  // A Type total: the sum over the Type's whole target interfaces of the point (interface
  // 0, scope Type): THE value of the verdicts, which its interface rows carry.
  bool type_total = false;
  // Not columns: the developer policy MarkNotApplied marked the row's cell (a flag); the
  // value's earlier levels through its join chain, oldest first; an eigenmode row's |p| of
  // the mode's inductive lumped ports (lumped port index -> |EPR|; the join fingerprint).
  bool marked = false;
  // Not columns (CUT-PLAN step 3): the run's impedance-metal rule "Apply" corrected a
  // feature on an impedance sheet by a perfect-conductor model, and a feature of the run
  // lies on a metal the rule does not serve (raw there). A Type's value sums every feature
  // of the run, so both are set on every row of the run; flags print their words.
  bool sheet_served = false, sheet_not_served = false;
  // Not a column: a tolerated library mismatch acts on this value (an interface row of a
  // listed target interface; a region row or a Type total of an affected Type).
  bool library_mismatch = false;
  // Not columns: the Maxwell gates of the row's fixed-trace estimate (flags and summary).
  ParticipationMaxwellGates gates;
  std::vector<ParticipationHistoryEntry> history;
  std::map<int, double> mode_ports;
};

// CUT-PLAN step 4b (lane v1step4b; the supervisor's ruling of 2026-11-25): the level-1 rows
// of a PreflightOnly run, which solves nothing. Per target interface (`targets`: interface
// -> Type; `regions`: the targets with a Region) one row of point 0 and excitation 0 whose
// values (p, Q, h/R m90 and m50, the share, the changes) are NaN, whose resolution reads
// the Type's all-planes PREFLIGHT word (ONE verdict, decision 2403 (3), lane v1verdict1:
// the table's word with no convergence evaluated, on the census weighted by length; a
// region row NOT-EVALUATED: a box has no census of its own before a solve), whose
// coverage, regime, confidence, convergence and summary read kParticipationNotComputed,
// whose information is ParticipationInformation's, whose flags
// read kParticipationPreflightFlags and whose reason is kParticipationPreflightReason.
// Where the Type's all-planes PREFLIGHT row has cells beyond twice the bar
// (ResolutionPreflightRow::UnderResolvedCells, the larger set: a box's patches are among
// them), every row of that Type, region rows too, adds the flags word
// kParticipationPreflightUnderResolvedCellsFlag; without it the table cannot read
// UNDER-RESOLVED on that mesh after the solve.
inline constexpr const char *kParticipationNotComputed = "NOT-COMPUTED";
inline constexpr const char *kParticipationPreflightFlags = "PREFLIGHT-ONLY";
inline constexpr const char *kParticipationPreflightUnderResolvedCellsFlag =
    "cells beyond twice the bar: may read UNDER-RESOLVED after the solve";
inline constexpr const char *kParticipationPreflightReason =
    "PREFLIGHT-ONLY: no solve; values not computed; resolution = the table's rule on the "
    "census weighted by length (no field before a solve); no convergence evaluated: "
    "never RESOLVED";
std::vector<SurfaceParticipationRow>
PreflightParticipationRows(const std::map<int, std::string> &targets,
                           const std::set<int> &regions,
                           const std::vector<ResolutionPreflightRow> &preflight,
                           const std::set<int> &mismatch_interfaces = {},
                           const std::set<std::string> &mismatch_types = {});

// The verdict inputs of a target row from its cells (nullptr = no cell: no census reading /
// NOT-MEASURED): m90 and m50 of the census cell (the resolution word is set by the join,
// which knows the convergence), coverage and regime of the coverage cell.
void SetParticipationResolution(SurfaceParticipationRow &row, const ResolutionCell *cell);
void SetParticipationCoverage(SurfaceParticipationRow &row, const CoverageCell *cell);
// An untargeted interface: raw only (p corrected NaN, verdicts NOT-TARGET).
void SetParticipationNotTarget(SurfaceParticipationRow &row);
// The census cell of `type` over every plane (nullptr when the census has none).
const ResolutionCell *FindParticipationCensusCell(const std::vector<ResolutionCell> &cells,
                                                  const std::string &type);

// The eigenmode join (DESIGN-NOTE s3). A mode's fingerprint is the raw participations of
// its whole interfaces with the |EPR| of its inductive lumped ports (the EPR's sign follows
// the eigenvector's phase); the overlap of two modes is the cosine of their fingerprints.
// Mode m of the current level joins mode j of the previous level iff (1) j is the previous
// mode nearest to f_m, (2) |f_m - f_j| is below half the distance from f_j to its nearest
// other previous mode, (3) m is the current mode nearest to f_j, (4) their overlap is at
// least kParticipationModeMinOverlap, (5) the join is unique: no other previous mode whose
// overlap with m passes lies within kParticipationModeAmbiguityRatio |f_m - f_j| of f_m,
// and (6) the same-field ranks agree: j's rank in frequency among the previous modes whose
// overlap with m passes equals m's rank among the current modes whose overlap with j
// passes (same-field fingerprints are nearly parallel, so a shift of the spectrum would
// otherwise join a mode to its neighbour).
inline constexpr double kParticipationModeMinOverlap = 0.9;
inline constexpr double kParticipationModeAmbiguityRatio = 2.0;
struct ParticipationModeSignature
{
  int point = 0;
  double frequency = 0.0;
  std::map<int, double> participation;  // interface -> p raw of the mode's whole rows
  std::map<int, double> ports;          // inductive lumped port -> |EPR|
};
// The join of a current mode: the previous mode nearest in frequency (-1: none), the two
// frequencies and the overlap, and whether conditions (1)-(6) hold (else the first that
// failed, in words).
struct ParticipationModeJoin
{
  int point = 0, nearest = -1;
  double frequency = std::numeric_limits<double>::quiet_NaN();
  double nearest_frequency = std::numeric_limits<double>::quiet_NaN();
  double overlap = std::numeric_limits<double>::quiet_NaN();
  bool joined = false;
  std::string failed;
};
// For every current mode (by point), its join.
std::map<int, ParticipationModeJoin>
MatchParticipationModes(const std::vector<ParticipationModeSignature> &previous,
                        const std::vector<ParticipationModeSignature> &current);

// The rows of a level with, per (point, excitation), one Type total per Type with whole
// target interfaces first (interface 0: p raw and p corrected summed, Q = 1 / sum 1 / Q_i,
// the Type's cells), then the given rows in their order.
std::vector<SurfaceParticipationRow>
AddParticipationTypeTotals(const std::vector<SurfaceParticipationRow> &rows);
// Joins the rows of a level to the previous level's (nullptr: the first level) and sets
// their history, change, indicator, convergence, resolution word and flags; a whole
// interface row then carries its Type total's indicator and words. ES and driven
// rows join by (point, excitation, interface; a driven frequency must agree to 1e-12
// relative); eigenmode rows through MatchParticipationModes, whose joins it returns (empty
// for the other problems and the first level).
std::map<int, ParticipationModeJoin>
JoinParticipationLevel(std::vector<SurfaceParticipationRow> &rows,
                       const std::vector<SurfaceParticipationRow> *previous,
                       const ParticipationVerdictBars &bars, int order);
// The flags column: the bad words only (lane v1wordscode), in the order resolution,
// coverage, the failed Maxwell gates (R2: kR, loop residual, loop fraction, library
// distance, boundary law), convergence, then the developer mark, the metal law not served
// and a tolerated library mismatch; empty only when every verdict is met with final (not
// interim) data. No regime or confidence word (R1, R3).
std::string ParticipationFlags(const SurfaceParticipationRow &row);
// The information column (R2, R4, R5): 'Ls sheet served as PEC', the two S words, R/rho
// beyond its limit, 'MA not certified' on an MA row; empty on a NOT-TARGET row.
std::string ParticipationInformation(const SurfaceParticipationRow &row);
// The summary column (R5 with the supervisor's two additions), in this order: NOT-TARGET on
// a NOT-TARGET row; DO-NOT-USE iff the metal law or the target is not served, the
// resolution is UNDER-RESOLVED(-UNTESTED) or NOT-EVALUATED, p corrected is not a finite
// number, or the coverage policy MarkNotApplied marked the row; OK iff the Type is SA or
// MS, the resolution RESOLVED or RESOLVED-UNTESTED, the coverage OK, no Maxwell gate of R2
// failed, no library mismatch is tolerated on the row, and the row is not a box (region)
// row while the bars' region scope is untested (review MAJOR 1, ruled (a): no box value has
// read OK on the benchmark, and a box census reads its own patches); CHECK otherwise
// (NO-BARS, every MA row). A PreflightOnly row reads NOT-COMPUTED
// (PreflightParticipationRows).
std::string ParticipationSummary(const SurfaceParticipationRow &row,
                                 const ParticipationVerdictBars &bars);
// The log's one information line per run on the trace-closure gates (R1): the largest
// fixed-trace / fixed-flux spread over interfaces at the first solve, not a verdict.
std::string ParticipationClosureGapLine(double maximum_spread);
// The log block of a finished level (DESIGN-NOTE s5): a header line, the eigenmode joins
// (one line per mode: joined or UNMATCHED, the two frequencies and the overlap), then the
// values (the Type totals and the targeted region rows): one line per value up to
// kParticipationLogRows values, else one line per (Type, whole | region).
inline constexpr std::size_t kParticipationLogRows = 40;
std::string ParticipationLogBlock(const std::vector<SurfaceParticipationRow> &rows,
                                  int iteration, long long int unknowns,
                                  long long int previous_unknowns,
                                  const ParticipationVerdictBars &bars, int order,
                                  const std::map<int, ParticipationModeJoin> &joins = {});

class SurfaceParticipationTable
{
  std::string problem;
  int order = 0;
  ParticipationVerdictBars bars;
  TableWithCSVFile file;
  std::vector<SurfaceParticipationRow> current, previous;
  long long int previous_unknowns = 0;
  int iterations = 0;
  std::size_t rows_written = 0;
  bool started = false;
  // CUT-PLAN step 4b (lane v1step4b): the file holds a PreflightOnly run's rows.
  bool preflight_only = false;
  // Lane v1wordscode (R1): the run's one closure-gap information line was printed.
  bool closure_gap_line_printed = false;
  // The rows into the file (root): the 28 columns, one line per row.
  void AppendRows(const std::vector<SurfaceParticipationRow> &rows);

public:
  // `path` the table file; `problem` Electrostatic / Eigenmode / Driven; `order` the FE
  // order of the run (the bar's order).
  SurfaceParticipationTable(std::string path, std::string problem, int order,
                            ParticipationVerdictBars bars);
  const std::string &Problem() const { return problem; }
  const ParticipationVerdictBars &Bars() const { return bars; }
  // A row of the level being solved (root; iteration, unknowns and the join are set by
  // FinishLevel).
  void AddRow(SurfaceParticipationRow row) { current.push_back(std::move(row)); }
  // The end of a level (collective call; root work): with `converged`, joins the level's
  // rows to the previous level's, appends them to the file and prints the log block;
  // otherwise drops them and prints one line.
  void FinishLevel(int iteration, long long int unknowns, bool converged, MPI_Comm comm);
  // CUT-PLAN step 4b (lane v1step4b): a PreflightOnly run's only level (collective call;
  // root work): the PreflightParticipationRows of iteration 1 with their Type totals
  // (AddParticipationTypeTotals), written without a join, an indicator or a log block of
  // values (one log line).
  void WritePreflightLevel(std::vector<SurfaceParticipationRow> rows,
                           long long int unknowns, MPI_Comm comm);
  // The rows of the last finished level (root).
  const std::vector<SurfaceParticipationRow> &LastLevel() const { return previous; }
  // Whether the file exists: a run that never has a row writes no file (root).
  bool Started() const { return started; }
  // True the first time it is asked in the run (every rank): the run prints its one
  // ParticipationClosureGapLine then (R1).
  bool TakeClosureGapLine() { return !std::exchange(closure_gap_line_printed, true); }
  // The palace.json block SurfaceResponse.ParticipationTable.
  nlohmann::json Metadata() const;
};

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_PARTICIPATION_TABLE_HPP
