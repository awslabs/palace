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
  NOT_TARGET
};

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
  int n = 0;
};
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

// One interface of one cell: W and the partition's sums (J).
struct CoverageInterfaceSummary
{
  double within = 0.0;            // W (raw - outside at R): the corrected bookkeeping's
  double partition_within = 0.0;  // the partition pass's own within-R total
  double served = 0.0, missing = 0.0, pending = 0.0, no_status = 0.0;
  double fabricated = 0.0;  // the correction's fabricated energy (evaluation 0)
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
  double uncovered_share = 0.0, missing_share = 0.0, method_gap_share = 0.0,
         pending_share = 0.0, no_status_share = 0.0;
  double factor = 0.0, validated_min = 0.0, validated_max = 0.0;
  CoverageStatus status = CoverageStatus::NOT_MEASURED;
  RegimeStatus regime = RegimeStatus::NO_MODEL;
  // Lane rescensus (decision 1637 (e)): its per-cell edge-resolution status is appended
  // here and printed after `regime status`; its stage-0 line beside the Missing length
  // share. Not part of this lane.
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
  std::vector<CoverageKeyRow> keys;
  // The served models (runtime index) with served energy whose library record carries no
  // QualificationStatus (item (5): information) and the pending ones.
  std::set<int> served_no_status_models, served_pending_models;
  const CoverageCell *FindCell(const std::string &type) const;
  // The worst cell over Types (CoverageStatusRank / RegimeStatusRank); nullptr when no
  // cell exists.
  const CoverageCell *WorstCoverageCell() const;
  const CoverageCell *WorstRegimeCell() const;
};
// `within` (W) and `fabricated` per target interface (the corrected bookkeeping's raw -
// outside at R and the correction's fabricated energy, evaluation 0); `types` the target
// interfaces' Types.
CoverageMeasurement ComputeCoverageMeasurement(
    const std::map<int, SurfacePostOperator::InterfaceCoveragePartition> &partitions,
    const std::map<int, std::string> &types, const std::map<int, double> &within,
    const std::map<int, double> &fabricated,
    const std::map<int, CoverageKeyRecord> &keys_by_feature,
    const std::map<int, CoverageModelQualification> &models,
    const CoverageThresholds &thresholds, const CoverageValidatedRegime &regime,
    CoverageExcitationClass excitation);
CoverageMeasurement NotMeasuredCoverageMeasurement(const std::map<int, std::string> &types,
                                                   CoverageExcitationClass excitation);

// The plain-text rows of surface-response-coverage.csv and
// surface-response-uncovered-energy-by-key.csv for one measurement (both solve families):
// `keys` = the row's leading cells ("m,exc" | "f (GHz),exc" | "i" values, comma-joined),
// `key_header` their header; energies in J and lengths in m through `units`; nan prints
// "nan". `portions_by_feature` = per uncovered feature its portions' count and length (mesh
// units). Lane rescensus (decision 1637 (e)) appends its per-cell edge-resolution status
// column after `regime status`.
std::string CoverageTableHeader(const std::string &key_header);
std::string CoverageByKeyTableHeader(const std::string &key_header,
                                     const std::set<int> &targets);
void WriteCoverageRows(std::ostream &output, const std::string &keys,
                       const CoverageMeasurement &measurement, const Units &units);
void WriteCoverageByKeyRows(
    std::ostream &output, const std::string &keys, const CoverageMeasurement &measurement,
    const std::set<int> &targets,
    const std::map<int, std::pair<int, double>> &portions_by_feature, const Units &units);
// Per uncovered feature: its portions' count and summed length (the by-key rows' columns).
std::map<int, std::pair<int, double>> CoveragePortionsByFeature(
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
  void Add(const std::vector<std::pair<std::string, double>> &keys,
           const CoverageMeasurement &measurement);
};
// palace.json SurfaceResponse.Coverage: the policy, the interim bars, the regime file, the
// status counts and one entry per cell whose coverage status is not OK or whose regime is
// OUT (energies in J by `energy_scale`; nan as null).
nlohmann::json CoverageRunMetadata(const CoverageRunRecord &record,
                                   const CoverageThresholds &thresholds,
                                   const CoverageValidatedRegime &regime,
                                   double energy_scale);
// The end-of-run banner of palace.json's block (a line per non-OK / OUT cell, capped; one
// NO-RANGE summary line; one line naming the status-less served models).
std::string CoverageRunBanner(const nlohmann::json &metadata);

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_RESPONSE_COVERAGE_HPP
