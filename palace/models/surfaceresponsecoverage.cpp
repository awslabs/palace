// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfaceresponsecoverage.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iterator>
#include <limits>
#include <set>
#include <sstream>
#include <fmt/format.h>
#include <mfem.hpp>
#include <nlohmann/json.hpp>
#include "models/surfaceresponseidentification.hpp"

namespace palace
{

namespace
{

// The '+'-joined tokens of a partition label.
std::vector<std::string> LabelTokens(const std::string &label)
{
  std::vector<std::string> tokens;
  std::size_t begin = 0;
  while (begin <= label.size())
  {
    const std::size_t end = std::min(label.find('+', begin), label.size());
    tokens.push_back(label.substr(begin, end - begin));
    begin = end + 1;
  }
  return tokens;
}

// The index of an "m<k>" token, or -1.
int TokenIndex(const std::string &token)
{
  if (token.size() < 2 || token[0] != 'm' ||
      !std::all_of(std::next(token.begin()), token.end(),
                   [](char c) { return c >= '0' && c <= '9'; }))
  {
    return -1;
  }
  return std::stoi(token.substr(1));
}

std::string JoinSorted(const std::set<std::string> &values)
{
  std::string joined;
  for (const auto &value : values)
  {
    joined += (joined.empty() ? "" : "+") + value;
  }
  return joined;
}

}  // namespace

std::string_view CoverageStatusName(CoverageStatus status)
{
  switch (status)
  {
    case CoverageStatus::OK:
      return "OK";
    case CoverageStatus::LOW:
      return "LOW";
    case CoverageStatus::PARTIAL:
      return "PARTIAL";
    case CoverageStatus::NOT_ACHIEVABLE:
      return "NOT-ACHIEVABLE";
    case CoverageStatus::NO_ENERGY:
      return "NO-ENERGY";
    case CoverageStatus::NOT_MEASURED:
      return "NOT-MEASURED";
    case CoverageStatus::NOT_TARGET:
      return "NOT-TARGET";
    case CoverageStatus::INCONSISTENT:
      return "INCONSISTENT";
    case CoverageStatus::TARGET_NOT_SERVED:
      return "TARGET-NOT-SERVED";
  }
  return "NOT-MEASURED";
}

std::string_view RegimeStatusName(RegimeStatus status)
{
  switch (status)
  {
    case RegimeStatus::IN:
      return "IN";
    case RegimeStatus::OUT:
      return "OUT";
    case RegimeStatus::NO_RANGE:
      return "NO-RANGE";
    case RegimeStatus::NO_MODEL:
      return "NO-MODEL";
  }
  return "NO-RANGE";
}

std::string_view CoverageExcitationClassName(CoverageExcitationClass excitation)
{
  switch (excitation)
  {
    case CoverageExcitationClass::ELECTROSTATIC_TERMINAL:
      return "ElectrostaticTerminal";
    case CoverageExcitationClass::EIGENMODE:
      return "Eigenmode";
    case CoverageExcitationClass::DRIVEN:
      return "Driven";
    case CoverageExcitationClass::NONE:
      return "None";
  }
  return "None";
}

int CoverageStatusRank(CoverageStatus status)
{
  switch (status)
  {
    case CoverageStatus::OK:
      return 0;
    case CoverageStatus::LOW:
      return 1;
    case CoverageStatus::PARTIAL:
      return 2;
    case CoverageStatus::NOT_ACHIEVABLE:
      return 3;
    case CoverageStatus::TARGET_NOT_SERVED:
      return 4;  // ruling 1738: worse than any band (the existing order unchanged)
    case CoverageStatus::NO_ENERGY:
      return 5;
    case CoverageStatus::NOT_MEASURED:
      return 6;
    case CoverageStatus::NOT_TARGET:
      return -1;
    case CoverageStatus::INCONSISTENT:
      return 7;  // a defect of the partition: worse than any coverage reading
  }
  return 6;
}

int RegimeStatusRank(RegimeStatus status)
{
  switch (status)
  {
    case RegimeStatus::IN:
      return 0;
    case RegimeStatus::NO_MODEL:
      return 1;
    case RegimeStatus::NO_RANGE:
      return 2;
    case RegimeStatus::OUT:
      return 3;
  }
  return 3;
}

CoverageValidatedRegime ReadCoverageValidatedRegime(const std::string &path,
                                                    const std::string &library_path)
{
  auto ReadBytes = [](const std::string &file, const char *what)
  {
    std::ifstream input(file, std::ios::binary);
    MFEM_VERIFY(input, "Unable to open the " << what << " \"" << file << "\"!");
    std::ostringstream bytes;
    bytes << input.rdbuf();
    return bytes.str();
  };
  const std::string text = ReadBytes(path, "coverage ValidatedRegime file");
  const auto data = nlohmann::json::parse(text, nullptr, false);
  MFEM_VERIFY(data.is_object(), "Coverage ValidatedRegime file \"" << path
                                                                   << "\" is not a JSON "
                                                                      "object!");
  MFEM_VERIFY(data.value("Version", 0) == 1,
              "Coverage ValidatedRegime file \"" << path << "\" must have Version 1!");
  CoverageValidatedRegime regime;
  regime.path = path;
  const auto library = data.find("Library");
  MFEM_VERIFY(library != data.end() && library->is_object() &&
                  library->contains("Sha256") && library->at("Sha256").is_string(),
              "Coverage ValidatedRegime file \"" << path << "\" requires Library.Sha256!");
  regime.library_name = library->value("Name", std::string{});
  regime.library_sha256 = library->at("Sha256").get<std::string>();
  const std::string library_sha256 =
      Sha256Hex(ReadBytes(library_path, "surface-response library"));
  MFEM_VERIFY(regime.library_sha256 == library_sha256,
              "Coverage ValidatedRegime file \""
                  << path << "\" was written for the library with SHA-256 "
                  << regime.library_sha256 << ", not the loaded library \"" << library_path
                  << "\" (SHA-256 " << library_sha256
                  << "): a range of another library would print a wrong regime status!");
  regime.tolerance = data.value("Tolerance", 0.0);
  MFEM_VERIFY(std::isfinite(regime.tolerance) && regime.tolerance >= 0.0 &&
                  regime.tolerance <= kCoverageRegimeToleranceMax,
              "Coverage ValidatedRegime \"Tolerance\" must be a number in [0, "
                  << kCoverageRegimeToleranceMax << "]!");
  for (const auto &entry : data.value("Ranges", nlohmann::json::array()))
  {
    MFEM_VERIFY(entry.is_object() && entry.contains("Type") &&
                    entry.contains("ExcitationClass") && entry.contains("FactorMin") &&
                    entry.contains("FactorMax"),
                "Coverage ValidatedRegime Ranges entries require Type, ExcitationClass, "
                "FactorMin and FactorMax!");
    const std::string type = entry.at("Type").get<std::string>();
    const std::string name = entry.at("ExcitationClass").get<std::string>();
    CoverageExcitationClass excitation = CoverageExcitationClass::NONE;
    for (const auto candidate :
         {CoverageExcitationClass::ELECTROSTATIC_TERMINAL,
          CoverageExcitationClass::EIGENMODE, CoverageExcitationClass::DRIVEN})
    {
      if (name == CoverageExcitationClassName(candidate))
      {
        excitation = candidate;
      }
    }
    MFEM_VERIFY(
        excitation != CoverageExcitationClass::NONE,
        "Coverage ValidatedRegime ExcitationClass must be \"ElectrostaticTerminal\", "
        "\"Eigenmode\" or \"Driven\", not \""
            << name << "\"!");
    CoverageValidatedRange range;
    range.factor_min = entry.at("FactorMin").get<double>();
    range.factor_max = entry.at("FactorMax").get<double>();
    // Review m-1: a range names how many windows validated it and where they come from.
    MFEM_VERIFY(entry.contains("N") && entry.at("N").is_number_integer() &&
                    entry.at("N").get<int>() >= 1,
                "Coverage ValidatedRegime range for "
                    << type << " / " << name
                    << " requires N (the number of validation windows) >= 1!");
    range.n = entry.at("N").get<int>();
    const auto provenance = entry.find("Provenance");
    MFEM_VERIFY(provenance != entry.end() &&
                    ((provenance->is_object() && !provenance->empty()) ||
                     (provenance->is_string() && !provenance->get<std::string>().empty())),
                "Coverage ValidatedRegime range for "
                    << type << " / " << name
                    << " requires a non-empty Provenance (object or text)!");
    range.provenance = provenance->dump();
    MFEM_VERIFY(std::isfinite(range.factor_min) && std::isfinite(range.factor_max) &&
                    range.factor_min <= range.factor_max,
                "Coverage ValidatedRegime range for " << type << " / " << name
                                                      << " needs FactorMin <= FactorMax!");
    MFEM_VERIFY(regime.ranges.emplace(std::make_pair(type, excitation), range).second,
                "Coverage ValidatedRegime lists " << type << " / " << name << " twice!");
  }
  return regime;
}

std::map<std::string, CoverageModelQualification>
ReadCoverageModelQualifications(const std::string &library_path)
{
  std::ifstream input(library_path);
  MFEM_VERIFY(input,
              "Unable to open the surface-response library \"" << library_path << "\"!");
  nlohmann::json data;
  input >> data;
  std::map<std::string, CoverageModelQualification> qualifications;
  for (const auto &entry : data.value("Models", nlohmann::json::array()))
  {
    if (!entry.is_object() || !entry.contains("Name") || !entry.at("Name").is_string())
    {
      continue;  // the load refuses such a record by name
    }
    auto qualification = CoverageModelQualification::NO_STATUS;
    const auto status = entry.find("QualificationStatus");
    const auto library_qualified = entry.find("LibraryQualified");
    if (status != entry.end() && status->is_string())
    {
      qualification = status->get<std::string>() == "Qualified"
                          ? CoverageModelQualification::QUALIFIED
                          : CoverageModelQualification::PENDING;
    }
    else if (library_qualified != entry.end() && library_qualified->is_boolean() &&
             !library_qualified->get<bool>())
    {
      qualification = CoverageModelQualification::PENDING;
    }
    qualifications[entry.at("Name").get<std::string>()] = qualification;
  }
  return qualifications;
}

CoverageKeyRecord ClassifyCoverageKey(const std::string &hash, bool matched,
                                      const std::string &placement_rejection,
                                      bool unboxable)
{
  // Ruling 1649 Q2: Unboxable and placement-rejected keys are method gap; every other
  // uncovered feature's key is Missing (a library could hold a model for it).
  CoverageKeyRecord record;
  record.hash = hash;
  if (unboxable)
  {
    record.kind = "Unboxable";
    record.method_gap = true;
  }
  else if (matched)
  {
    record.kind = "PlacementRejected";
    record.reason = placement_rejection;
    record.method_gap = true;
  }
  else
  {
    record.kind = "Missing";
  }
  return record;
}

CoverageInterfaceSummary
SummarizeCoverageInterface(const SurfacePostOperator::InterfaceCoveragePartition &partition,
                           double within, double fabricated,
                           const std::map<int, CoverageKeyRecord> &keys_by_feature,
                           const std::map<int, CoverageModelQualification> &models)
{
  CoverageInterfaceSummary summary;
  summary.within = within;
  summary.fabricated = fabricated;
  summary.partition_within = partition.within;
  summary.loop_remainder_unserved = partition.loop_remainder_unserved;
  summary.served_ends = partition.served_ends;
  summary.kept_ends = partition.kept_ends;
  for (const auto &[label, energy] : partition.served)
  {
    summary.served += energy;
    bool pending = false, no_status = false;
    for (const auto &token : LabelTokens(label))
    {
      // Review m-5: a served model without a qualification record (an explicit-model
      // configuration has none at all) is NO_STATUS, never silently qualified.
      const auto model = models.find(TokenIndex(token));
      const auto qualification =
          model != models.end() ? model->second : CoverageModelQualification::NO_STATUS;
      pending = pending || qualification == CoverageModelQualification::PENDING;
      no_status = no_status || qualification == CoverageModelQualification::NO_STATUS;
    }
    // Shared labels count WHOLE (DESIGN-v2 s1: an upper bound); a label holding a pending
    // model is pending, not no-status.
    summary.pending += pending ? energy : 0.0;
    summary.no_status += (!pending && no_status) ? energy : 0.0;
  }
  for (const auto &[key, energy] : partition.classes)
  {
    if (key.first != "Missing")
    {
      continue;
    }
    // A bin joining several features is Missing unless EVERY token is method gap; a
    // feature with no key record counts as Missing (method gap stays a lower bound).
    bool method_gap = true;
    for (const auto &token : LabelTokens(key.second))
    {
      const auto record = keys_by_feature.find(TokenIndex(token));
      method_gap =
          method_gap && record != keys_by_feature.end() && record->second.method_gap;
    }
    summary.missing += method_gap ? 0.0 : energy;
  }
  return summary;
}

CoverageCell ComputeCoverageCell(const std::string &type,
                                 const std::vector<CoverageInterfaceSummary> &interfaces,
                                 const CoverageThresholds &thresholds,
                                 const CoverageValidatedRegime &regime,
                                 CoverageExcitationClass excitation)
{
  const double nan = std::numeric_limits<double>::quiet_NaN();
  CoverageCell cell;
  cell.type = type;
  for (const auto &interface : interfaces)
  {
    cell.within += interface.within;
    cell.served += interface.served;
    cell.missing += interface.missing;
    cell.pending += interface.pending;
    cell.no_status += interface.no_status;
    cell.fabricated += interface.fabricated;
    cell.partition_within += interface.partition_within;
    cell.loop_remainder_unserved += interface.loop_remainder_unserved;
    cell.loop_remainder_excluded += interface.loop_remainder_excluded;
    for (const auto &[sum, add] :
         {std::make_pair(&cell.served_ends, &interface.served_ends),
          std::make_pair(&cell.kept_ends, &interface.kept_ends)})
    {
      sum->roundoff_count += add->roundoff_count;
      sum->near_count += add->near_count;
      sum->roundoff_distance = std::max(sum->roundoff_distance, add->roundoff_distance);
      sum->near_distance = std::max(sum->near_distance, add->near_distance);
      sum->roundoff_energy += add->roundoff_energy;
      sum->near_energy += add->near_energy;
      sum->roundoff_beyond_count += add->roundoff_beyond_count;
      sum->roundoff_beyond_distance =
          std::max(sum->roundoff_beyond_distance, add->roundoff_beyond_distance);
    }
  }
  // Ruling 1641: U is fail-closed, whatever is not served is uncovered. The method gap
  // (U - missing) is clamped at 0: a negative value is round-off of the two integrals
  // (review m-3); the partition's own within-R total is printed beside W (the residual).
  // Decision 1708 (1): the clamp must not hide a negative gap, so the signed share before
  // the clamp is printed beside it (its own column; never used for a status).
  cell.uncovered = cell.within - cell.served;
  cell.method_gap = std::max(0.0, cell.uncovered - cell.missing);
  if (cell.within > 0.0)
  {
    cell.uncovered_share = cell.uncovered / cell.within;
    cell.missing_share = cell.missing / cell.within;
    cell.method_gap_share = cell.method_gap / cell.within;
    cell.method_gap_signed_share = (cell.uncovered - cell.missing) / cell.within;
    cell.pending_share = cell.pending / cell.within;
    cell.no_status_share = cell.no_status / cell.within;
    if (cell.method_gap_signed_share < -kCoverageInconsistentShare)
    {
      // Ruling 1714 (B)(3): served + Missing above W beyond round-off: the partition and
      // the corrected total disagree (U stays as computed).
      cell.status = CoverageStatus::INCONSISTENT;
    }
    else if (cell.method_gap_share > thresholds.max_uncovered_share)
    {
      cell.status = CoverageStatus::NOT_ACHIEVABLE;
    }
    else
    {
      const double share = cell.uncovered_share + cell.pending_share;
      cell.status = share <= thresholds.build_target          ? CoverageStatus::OK
                    : share <= thresholds.max_uncovered_share ? CoverageStatus::LOW
                                                              : CoverageStatus::PARTIAL;
    }
  }
  else
  {
    // Ruling 1654 (III) Q1: nothing to cover; the shares print nan, never 0.
    cell.uncovered_share = cell.missing_share = cell.method_gap_share = cell.pending_share =
        cell.no_status_share = cell.method_gap_signed_share = nan;
    cell.status = CoverageStatus::NO_ENERGY;
  }
  const auto range = regime.ranges.find({type, excitation});
  if (!(cell.served > 0.0) || cell.status == CoverageStatus::NO_ENERGY)
  {
    cell.factor = cell.validated_min = cell.validated_max = nan;
    cell.regime = RegimeStatus::NO_MODEL;
  }
  else if (range == regime.ranges.end())
  {
    cell.factor = cell.fabricated / cell.served;
    cell.validated_min = cell.validated_max = nan;
    cell.regime = RegimeStatus::NO_RANGE;
  }
  else
  {
    // DESIGN-v2 s6.2: IN when F_lo / (1 + tau) <= F_cov <= F_hi (1 + tau).
    cell.factor = cell.fabricated / cell.served;
    cell.validated_min = range->second.factor_min;
    cell.validated_max = range->second.factor_max;
    cell.validated_n = range->second.n;
    const double tau = regime.tolerance;
    cell.regime = (cell.factor >= cell.validated_min / (1.0 + tau) &&
                   cell.factor <= cell.validated_max * (1.0 + tau))
                      ? RegimeStatus::IN
                      : RegimeStatus::OUT;
  }
  return cell;
}

CoverageCell NotMeasuredCoverageCell(const std::string &type)
{
  const double nan = std::numeric_limits<double>::quiet_NaN();
  CoverageCell cell;
  cell.type = type;
  cell.within = cell.uncovered = cell.missing = cell.method_gap = cell.pending =
      cell.no_status = cell.served = cell.fabricated = cell.partition_within =
          cell.loop_remainder_unserved = cell.loop_remainder_excluded = nan;
  cell.uncovered_share = cell.missing_share = cell.method_gap_share = cell.pending_share =
      cell.no_status_share = cell.method_gap_signed_share = nan;
  cell.factor = cell.validated_min = cell.validated_max = nan;
  cell.status = CoverageStatus::NOT_MEASURED;
  cell.regime = RegimeStatus::NO_MODEL;
  return cell;
}

std::vector<CoverageKeyRow>
CoverageKeyRows(const std::map<int, SurfacePostOperator::InterfaceCoveragePartition> &cell,
                const std::map<int, CoverageKeyRecord> &keys_by_feature)
{
  std::map<std::string, CoverageKeyRow> rows;
  std::map<std::string, std::pair<std::set<std::string>, std::set<std::string>>> labels;
  for (const auto &[interface, partition] : cell)
  {
    for (const auto &[key, energy] : partition.classes)
    {
      if (key.first != "Missing")
      {
        continue;
      }
      std::set<std::string> hashes, kinds, reasons;
      std::vector<int> features;
      bool method_gap = true;
      for (const auto &token : LabelTokens(key.second))
      {
        const int feature = TokenIndex(token);
        features.push_back(feature);
        const auto record = keys_by_feature.find(feature);
        if (record == keys_by_feature.end())
        {
          hashes.insert("unrecorded-" + token);
          kinds.insert("Missing");
          method_gap = false;
          continue;
        }
        hashes.insert(record->second.hash);
        kinds.insert(record->second.kind);
        if (!record->second.reason.empty())
        {
          reasons.insert(record->second.reason);
        }
        method_gap = method_gap && record->second.method_gap;
      }
      // Bins of one key set merge: the kinds and reasons by union, the method-gap flag
      // only when every bin's is.
      const std::string joined = JoinSorted(hashes);
      auto &row = rows[joined];
      auto &[row_kinds, row_reasons] = labels[joined];
      row.method_gap = (row.key.empty() || row.method_gap) && method_gap;
      row.key = joined;
      row_kinds.insert(kinds.begin(), kinds.end());
      row_reasons.insert(reasons.begin(), reasons.end());
      row.kind = JoinSorted(row_kinds);
      row.reason = JoinSorted(row_reasons);
      for (const int feature : features)
      {
        if (std::find(row.features.begin(), row.features.end(), feature) ==
            row.features.end())
        {
          row.features.push_back(feature);
        }
      }
      row.energy[interface] += energy;
    }
  }
  std::vector<CoverageKeyRow> result;
  for (auto &[key, row] : rows)
  {
    (void)key;
    std::sort(row.features.begin(), row.features.end());
    result.push_back(std::move(row));
  }
  return result;
}

namespace
{

// The Types of the target interfaces in order of first appearance (ascending targets).
std::vector<std::string> CoverageTypeOrder(const std::map<int, std::string> &types)
{
  std::vector<std::string> order;
  for (const auto &[interface, type] : types)
  {
    (void)interface;
    if (std::find(order.begin(), order.end(), type) == order.end())
    {
      order.push_back(type);
    }
  }
  return order;
}

}  // namespace

const CoverageCell *CoverageMeasurement::FindCell(const std::string &type) const
{
  const auto cell = std::find_if(cells.begin(), cells.end(),
                                 [&](const auto &entry) { return entry.type == type; });
  return cell != cells.end() ? &*cell : nullptr;
}

const CoverageCell *CoverageMeasurement::WorstCoverageCell() const
{
  const CoverageCell *worst = nullptr;
  for (const auto &cell : cells)
  {
    if (cell.status == CoverageStatus::NOT_TARGET)
    {
      continue;
    }
    if (!worst || CoverageStatusRank(cell.status) > CoverageStatusRank(worst->status))
    {
      worst = &cell;
    }
  }
  return worst;
}

const CoverageCell *CoverageMeasurement::WorstRegimeCell() const
{
  const CoverageCell *worst = nullptr;
  for (const auto &cell : cells)
  {
    if (cell.status == CoverageStatus::NOT_TARGET)
    {
      continue;
    }
    if (!worst || RegimeStatusRank(cell.regime) > RegimeStatusRank(worst->regime))
    {
      worst = &cell;
    }
  }
  return worst;
}

CoverageMeasurement ComputeCoverageMeasurement(
    const std::map<int, SurfacePostOperator::InterfaceCoveragePartition> &partitions,
    const std::map<int, std::string> &types, const std::map<int, double> &within,
    const std::map<int, double> &fabricated,
    const std::map<int, CoverageKeyRecord> &keys_by_feature,
    const std::map<int, CoverageModelQualification> &models,
    const CoverageThresholds &thresholds, const CoverageValidatedRegime &regime,
    CoverageExcitationClass excitation,
    const std::map<int, double> &loop_remainder_excluded)
{
  CoverageMeasurement measurement;
  measurement.measured = true;
  measurement.excitation = excitation;
  measurement.types = types;
  const SurfacePostOperator::InterfaceCoveragePartition empty;
  for (const auto &[interface, type] : types)
  {
    (void)type;
    const auto partition = partitions.find(interface);
    const auto w = within.find(interface);
    const auto fab = fabricated.find(interface);
    // An interface with no partition entry (no within-R point on any rank) has nothing
    // served: it is uncovered by its W (fail-closed).
    measurement.interfaces[interface] = SummarizeCoverageInterface(
        partition != partitions.end() ? partition->second : empty,
        w != within.end() ? w->second : 0.0, fab != fabricated.end() ? fab->second : 0.0,
        keys_by_feature, models);
    const auto excluded = loop_remainder_excluded.find(interface);
    measurement.interfaces[interface].loop_remainder_excluded =
        excluded != loop_remainder_excluded.end() ? excluded->second : 0.0;
  }
  for (const auto &type : CoverageTypeOrder(types))
  {
    std::vector<CoverageInterfaceSummary> summaries;
    for (const auto &[interface, interface_type] : types)
    {
      if (interface_type == type)
      {
        summaries.push_back(measurement.interfaces.at(interface));
      }
    }
    measurement.cells.push_back(
        ComputeCoverageCell(type, summaries, thresholds, regime, excitation));
  }
  measurement.keys = CoverageKeyRows(partitions, keys_by_feature);
  for (const auto &[interface, partition] : partitions)
  {
    if (types.count(interface) == 0)
    {
      continue;
    }
    for (const auto &[label, energy] : partition.served)
    {
      if (!(energy > 0.0))
      {
        continue;
      }
      for (const auto &token : LabelTokens(label))
      {
        const auto model = models.find(TokenIndex(token));
        if (model != models.end() && model->second == CoverageModelQualification::PENDING)
        {
          measurement.served_pending_models.insert(model->first);
        }
        else if (model != models.end() &&
                 model->second == CoverageModelQualification::NO_STATUS)
        {
          measurement.served_no_status_models.insert(model->first);
        }
      }
    }
  }
  return measurement;
}

CoverageMeasurement NotMeasuredCoverageMeasurement(const std::map<int, std::string> &types,
                                                   CoverageExcitationClass excitation)
{
  CoverageMeasurement measurement;
  measurement.excitation = excitation;
  measurement.types = types;
  for (const auto &type : CoverageTypeOrder(types))
  {
    measurement.cells.push_back(NotMeasuredCoverageCell(type));
  }
  return measurement;
}

void AddCoverageUntargetedInterfaces(CoverageMeasurement &measurement,
                                     const std::map<int, std::string> &untargeted,
                                     const std::map<int, double> &raw_energy)
{
  const double nan = std::numeric_limits<double>::quiet_NaN();
  measurement.untargeted = untargeted;
  auto Raw = [&](int interface)
  {
    const auto raw = raw_energy.find(interface);
    return raw != raw_energy.end() ? raw->second : 0.0;
  };
  auto Fill = [&](CoverageCell &cell)
  {
    double targeted = 0.0, untargeted_raw = 0.0;
    for (const auto &[interface, type] : measurement.types)
    {
      targeted += type == cell.type ? Raw(interface) : 0.0;
    }
    for (const auto &[interface, type] : untargeted)
    {
      untargeted_raw += type == cell.type ? Raw(interface) : 0.0;
    }
    cell.untargeted_raw = untargeted_raw;
    cell.untargeted_share = targeted + untargeted_raw > 0.0
                                ? untargeted_raw / (targeted + untargeted_raw)
                                : nan;
  };
  for (auto &cell : measurement.cells)
  {
    Fill(cell);
  }
  for (const auto &type : CoverageTypeOrder(untargeted))
  {
    if (measurement.FindCell(type))
    {
      continue;
    }
    // Ruling 1694 (A): no target = no correction radius = an empty within-R region.
    CoverageCell cell;
    cell.type = type;
    cell.uncovered_share = cell.missing_share = cell.method_gap_share = cell.pending_share =
        cell.no_status_share = cell.method_gap_signed_share = nan;
    cell.partition_within = 0.0;
    cell.factor = cell.validated_min = cell.validated_max = nan;
    cell.status = CoverageStatus::NOT_TARGET;
    cell.regime = RegimeStatus::NO_MODEL;
    Fill(cell);
    measurement.cells.push_back(cell);
  }
}

void MarkCoverageTargetsNotServed(CoverageMeasurement &measurement,
                                  const std::map<int, std::string> &not_served)
{
  measurement.not_served = not_served;
  for (auto &cell : measurement.cells)
  {
    const bool holds =
        std::any_of(not_served.begin(), not_served.end(),
                    [&](const auto &entry) { return entry.second == cell.type; });
    if (holds && CoverageStatusRank(CoverageStatus::TARGET_NOT_SERVED) >
                     CoverageStatusRank(cell.status))
    {
      cell.status = CoverageStatus::TARGET_NOT_SERVED;
    }
  }
}

std::string CoverageTargetsNotServedWarning(const std::map<int, std::string> &not_served)
{
  std::string text;
  for (const auto &[interface, type] : not_served)
  {
    text += fmt::format("Configured response-correction target interface {:d} ({}) is NOT "
                        "SERVED: {}; it is not corrected (its corrected values are its raw "
                        "values) and its within-R raw energy counts as uncovered in its "
                        "Type's coverage (status TARGET-NOT-SERVED)\n",
                        interface, type, kCoverageTargetNotServedReason);
  }
  return text;
}

namespace
{

std::string CoverageNumber(double value)
{
  return std::isnan(value) ? std::string("nan") : fmt::format("{:+.12e}", value);
}

}  // namespace

std::string CoverageTableHeader(const std::string &key_header)
{
  return key_header +
         ",Type,W within-R raw (J),uncovered raw share of W,missing raw share "
         "of W,method-gap raw share of W,pending-model share of W,no-status-"
         "model share of W,covered correction factor,validated factor "
         "min,validated factor max,coverage status,regime status,partition "
         "within-R raw (J),untargeted raw energy of the Type (J),untargeted "
         "raw share of the Type interface energy,loop remainder unserved "
         "(J),loop remainder of excluded patches kept raw twice (J),method-gap "
         "raw share of W signed (before the clamp at 0)\n";
}

std::string CoverageByKeyTableHeader(const std::string &key_header,
                                     const std::set<int> &targets)
{
  std::string header =
      key_header +
      ",evaluation,key hash,class,reason,method gap,features,portions,length (m)";
  for (const int interface : targets)
  {
    header += fmt::format(",raw within-R energy[{}] (J),raw share of W[{}]", interface,
                          interface);
  }
  return header + "\n";
}

void WriteCoverageRows(std::ostream &output, const std::string &keys,
                       const CoverageMeasurement &measurement, const Units &units)
{
  using VT = Units::ValueType;
  for (const auto &cell : measurement.cells)
  {
    output << keys << ',' << cell.type << ','
           << CoverageNumber(units.Dimensionalize<VT::ENERGY>(cell.within)) << ','
           << CoverageNumber(cell.uncovered_share) << ','
           << CoverageNumber(cell.missing_share) << ','
           << CoverageNumber(cell.method_gap_share) << ','
           << CoverageNumber(cell.pending_share) << ','
           << CoverageNumber(cell.no_status_share) << ',' << CoverageNumber(cell.factor)
           << ',' << CoverageNumber(cell.validated_min) << ','
           << CoverageNumber(cell.validated_max) << ',' << CoverageStatusName(cell.status)
           << ',' << RegimeStatusName(cell.regime) << ','
           << CoverageNumber(units.Dimensionalize<VT::ENERGY>(cell.partition_within)) << ','
           << CoverageNumber(units.Dimensionalize<VT::ENERGY>(cell.untargeted_raw)) << ','
           << CoverageNumber(cell.untargeted_share) << ','
           << CoverageNumber(units.Dimensionalize<VT::ENERGY>(cell.loop_remainder_unserved))
           << ','
           << CoverageNumber(units.Dimensionalize<VT::ENERGY>(cell.loop_remainder_excluded))
           << ',' << CoverageNumber(cell.method_gap_signed_share) << '\n';
  }
}

void WriteCoverageByKeyRows(std::ostream &output, const std::string &keys,
                            const CoverageMeasurement &measurement,
                            const std::set<int> &targets,
                            const CoverageFeaturePortionsMap &portions_by_feature,
                            const Units &units)
{
  using VT = Units::ValueType;
  for (const auto &row : measurement.keys)
  {
    std::string features, reason = row.reason;
    int portions = 0;
    double length = 0.0;
    bool mirror_formed = false;
    for (const int feature : row.features)
    {
      features += (features.empty() ? "" : "+") + std::to_string(feature);
      const auto entry = portions_by_feature.find(feature);
      if (entry != portions_by_feature.end())
      {
        portions += entry->second.portions;
        length += entry->second.length;
        mirror_formed = mirror_formed || entry->second.mirror_formed;
      }
    }
    if (mirror_formed)
    {
      reason += (reason.empty() ? "" : "+") + std::string("MirrorFormed");
    }
    std::replace(reason.begin(), reason.end(), ',', ';');
    // Evaluation 0: the raw field (the partition's field; the fixed-trace evaluation's).
    output << keys << ",0," << row.key << ',' << row.kind << ','
           << (reason.empty() ? "-" : reason) << ',' << (row.method_gap ? 1 : 0) << ','
           << features << ',' << portions << ','
           << CoverageNumber(units.Dimensionalize<VT::LENGTH>(length));
    for (const int interface : targets)
    {
      const auto energy = row.energy.find(interface);
      const double value = energy != row.energy.end() ? energy->second : 0.0;
      const auto summary = measurement.interfaces.find(interface);
      const double within =
          summary != measurement.interfaces.end() ? summary->second.within : 0.0;
      output << ',' << CoverageNumber(units.Dimensionalize<VT::ENERGY>(value)) << ','
             << CoverageNumber(within > 0.0 ? value / within
                                            : std::numeric_limits<double>::quiet_NaN());
    }
    output << '\n';
  }
}

CoverageFeaturePortionsMap CoveragePortionsByFeature(
    const std::vector<SurfacePostOperator::UncoveredPerimeterPortion> &uncovered)
{
  constexpr std::string_view mirror = ":MirrorFormed";
  CoverageFeaturePortionsMap portions_by_feature;
  for (const auto &portion : uncovered)
  {
    double length2 = 0.0;
    for (int d = 0; d < 3; d++)
    {
      length2 += (portion.p1[d] - portion.p0[d]) * (portion.p1[d] - portion.p0[d]);
    }
    auto &entry = portions_by_feature[portion.feature];
    entry.portions++;
    entry.length += std::sqrt(length2);
    entry.mirror_formed =
        entry.mirror_formed || (portion.type.size() >= mirror.size() &&
                                portion.type.compare(portion.type.size() - mirror.size(),
                                                     mirror.size(), mirror) == 0);
  }
  return portions_by_feature;
}

void CoverageRunRecord::Add(const std::vector<std::pair<std::string, double>> &keys,
                            const CoverageMeasurement &measurement)
{
  for (const auto &cell : measurement.cells)
  {
    entries.push_back({keys, cell});
  }
  served_no_status_models.insert(measurement.served_no_status_models.begin(),
                                 measurement.served_no_status_models.end());
  served_pending_models.insert(measurement.served_pending_models.begin(),
                               measurement.served_pending_models.end());
  untargeted.insert(measurement.untargeted.begin(), measurement.untargeted.end());
  not_served.insert(measurement.not_served.begin(), measurement.not_served.end());
}

nlohmann::json CoverageRunMetadata(const CoverageRunRecord &record,
                                   const CoverageThresholds &thresholds,
                                   const CoverageValidatedRegime &regime,
                                   double energy_scale, int level, int levels)
{
  auto Number = [](double value)
  { return std::isfinite(value) ? nlohmann::json(value) : nlohmann::json(nullptr); };
  nlohmann::json metadata = nlohmann::json::object();
  metadata["Policy"] = "Warn";
  metadata["Thresholds"] = {{"BuildTarget", thresholds.build_target},
                            {"MaxUncoveredShare", thresholds.max_uncovered_share},
                            {"Interim", true}};
  metadata["ValidatedRegime"] =
      regime.path.empty()
          ? nlohmann::json(nullptr)
          : nlohmann::json{{"Path", regime.path},
                           {"LibraryName", regime.library_name},
                           {"LibrarySha256", regime.library_sha256},
                           {"Tolerance", regime.tolerance},
                           {"Ranges", static_cast<int>(regime.ranges.size())}};
  metadata["Cells"] = static_cast<int>(record.entries.size());
  // Ruling 1689 M-4: the measured cell set, written at the END of each solve with
  // "Complete": true (gate v2 fails without it or when the coverage table's row count
  // differs; erratum 1784 E-3: a run killed during a later AMR level keeps this block of
  // an earlier level, without "Final"). DESIGN-v2 s4.5's expected set "Of" / "Certificate"
  // and the stage-0 PortAdjacent count are DEFERRED to the library-tool lane.
  {
    std::set<std::vector<double>> measurements;
    std::vector<std::string> types;
    std::string kind;
    for (const auto &entry : record.entries)
    {
      std::vector<double> key;
      for (const auto &[name, value] : entry.keys)
      {
        kind = kind.empty() ? name : kind;
        key.push_back(value);
      }
      measurements.insert(key);
      if (std::find(types.begin(), types.end(), entry.cell.type) == types.end())
      {
        types.push_back(entry.cell.type);
      }
    }
    metadata["CellSet"] = {{"Kind", kind},
                           {"Measurements", static_cast<int>(measurements.size())},
                           {"Types", types},
                           {"Cells", static_cast<int>(record.entries.size())}};
  }
  nlohmann::json untargeted = nlohmann::json::array();
  for (const auto &[interface, type] : record.untargeted)
  {
    untargeted.push_back({{"Interface", interface}, {"Type", type}});
  }
  metadata["UntargetedInterfaces"] = std::move(untargeted);
  // Ruling 1738: the configured targets no loaded model serves (in their Type's W,
  // uncovered; the Type's cells TARGET-NOT-SERVED).
  nlohmann::json not_served = nlohmann::json::array();
  for (const auto &[interface, type] : record.not_served)
  {
    not_served.push_back({{"Interface", interface},
                          {"Type", type},
                          {"Reason", kCoverageTargetNotServedReason}});
  }
  metadata["ConfiguredTargetsNotServed"] = std::move(not_served);
  nlohmann::json status_counts = nlohmann::json::object(),
                 regime_counts = nlohmann::json::object();
  nlohmann::json flagged = nlohmann::json::array();
  // Rulings 1714 / 1715: the end census per population over the cells (the largest count,
  // distance and energy share of W of any cell; the cells with a nonzero count listed).
  using EndCensus = SurfacePostOperator::InterfaceCoveragePartition::EndCensus;
  nlohmann::json end_cells = nlohmann::json::array();
  std::map<std::string, std::array<double, 4>> end_totals;  // count, distance, J, share
  auto EndPopulation =
      [&](const std::string &name, int count, double distance, double energy, double within)
  {
    auto &total = end_totals[name];
    total[0] = std::max(total[0], double(count));
    total[1] = std::max(total[1], distance);
    total[2] = std::max(total[2], energy * energy_scale);
    total[3] = std::max(total[3], within > 0.0 ? energy / within : 0.0);
    return nlohmann::json{
        {"Count", count}, {"MaxDistance", distance}, {"Energy (J)", energy * energy_scale}};
  };
  for (const auto &entry : record.entries)
  {
    const auto &cell = entry.cell;
    {
      auto Populations = [&](const std::string &prefix, const EndCensus &census)
      {
        return nlohmann::json{
            {"RoundoffShort",
             EndPopulation(prefix + "RoundoffShort", census.roundoff_count,
                           census.roundoff_distance, census.roundoff_energy, cell.within)},
            {"NearNotRoundoff",
             EndPopulation(prefix + "NearNotRoundoff", census.near_count,
                           census.near_distance, census.near_energy, cell.within)},
            {"RoundoffBeyond",
             EndPopulation(prefix + "RoundoffBeyond", census.roundoff_beyond_count,
                           census.roundoff_beyond_distance, 0.0, cell.within)}};
      };
      auto served = Populations("Served", cell.served_ends);
      auto kept = Populations("Kept", cell.kept_ends);
      if (cell.served_ends.roundoff_count + cell.served_ends.near_count +
              cell.served_ends.roundoff_beyond_count + cell.kept_ends.roundoff_count +
              cell.kept_ends.near_count + cell.kept_ends.roundoff_beyond_count >
          0)
      {
        nlohmann::json keys = nlohmann::json::object();
        for (const auto &[name, value] : entry.keys)
        {
          keys[name] = value;
        }
        end_cells.push_back({{"Cell", std::move(keys)},
                             {"Type", cell.type},
                             {"Served", std::move(served)},
                             {"Kept", std::move(kept)}});
      }
    }
    const std::string status(CoverageStatusName(cell.status));
    const std::string regime_status(RegimeStatusName(cell.regime));
    status_counts[status] = status_counts.value(status, 0) + 1;
    regime_counts[regime_status] = regime_counts.value(regime_status, 0) + 1;
    // Erratum 1784 E-2 (ERRATUM-1581 MJ-2: NO-RANGE is never OK, a contact like OUT): a
    // cell is flagged unless its coverage is OK AND its regime IN; an OK + NO-RANGE cell is
    // flagged per cell as an OK + OUT cell is.
    if ((cell.status == CoverageStatus::OK && cell.regime == RegimeStatus::IN) ||
        cell.status == CoverageStatus::NOT_TARGET)
    {
      continue;  // a NOT-TARGET Type: named once in UntargetedInterfaces
    }
    nlohmann::json keys = nlohmann::json::object();
    for (const auto &[name, value] : entry.keys)
    {
      keys[name] = value;
    }
    flagged.push_back({{"Cell", std::move(keys)},
                       {"Type", cell.type},
                       {"W within-R raw (J)", Number(cell.within * energy_scale)},
                       {"uncovered raw share of W", Number(cell.uncovered_share)},
                       {"missing raw share of W", Number(cell.missing_share)},
                       {"method-gap raw share of W", Number(cell.method_gap_share)},
                       {"pending-model share of W", Number(cell.pending_share)},
                       {"no-status-model share of W", Number(cell.no_status_share)},
                       {"covered correction factor", Number(cell.factor)},
                       {"validated factor min", Number(cell.validated_min)},
                       {"validated factor max", Number(cell.validated_max)},
                       {"validated N", cell.validated_n},
                       {"coverage status", status},
                       {"regime status", regime_status},
                       {"untargeted raw share of the Type interface energy",
                        Number(cell.untargeted_share)},
                       {"loop remainder of excluded patches kept raw twice (J)",
                        Number(cell.loop_remainder_excluded * energy_scale)}});
  }
  nlohmann::json end_census = {
      {"Rule",
       fmt::format("by direction (ruling 1731): a served claim end within {:g} ulps per "
                   "coordinate of a tree-segment end point whose clamped foot falls short "
                   "of the segment end (round-off short) is snapped to it; its energy = "
                   "the within-R raw energy at that end the snap gives back; kept-class "
                   "ends are counted, never changed (OI-FOOT-TEST-VERTEX-ROUNDOFF), their "
                   "energy = the within-R raw energy at that end no kept class holds; an "
                   "end within the same ulps whose clamped foot reaches the segment end "
                   "exactly (round-off beyond) is counted only (nothing at stake); near = "
                   "beyond the ulps, within {:g} R (never snapped; its energy an upper "
                   "bound)",
                   kClaimEndRoundoffUlps, kSupportEndCoincidenceOverRadius)},
      {"Cells", std::move(end_cells)}};
  for (const auto &[name, total] : end_totals)
  {
    end_census[name] = {{"MaxCount", static_cast<int>(total[0])},
                        {"MaxDistance", total[1]},
                        {"MaxEnergy (J)", total[2]},
                        {"MaxShareOfW", total[3]}};
  }
  metadata["ClaimEndCensus"] = std::move(end_census);
  metadata["StatusCounts"] = std::move(status_counts);
  metadata["RegimeCounts"] = std::move(regime_counts);
  metadata["Flagged"] = std::move(flagged);
  metadata["ServedPendingModels"] = static_cast<int>(record.served_pending_models.size());
  metadata["ServedNoStatusModels"] =
      static_cast<int>(record.served_no_status_models.size());
  // Ruling 1738: a run with no cell says why; without a reason it is not complete.
  if (record.entries.empty())
  {
    metadata["NoCellsReason"] = record.no_cells_reason.empty()
                                    ? nlohmann::json(nullptr)
                                    : nlohmann::json(record.no_cells_reason);
  }
  metadata["Complete"] = !record.entries.empty() || !record.no_cells_reason.empty();
  // Erratum 1784 E-3: "Complete" says the block holds every cell of its level; under AMR
  // the root holds the last COMPLETED level's block until the next level writes, so the
  // block names its level and only BaseSolver's write after the AMR loop adds "Final".
  metadata["Level"] = level;
  metadata["Levels"] = levels;
  return metadata;
}

std::string CoverageRunBanner(const nlohmann::json &metadata)
{
  constexpr int maximum_lines = 40;
  auto Percent = [](const nlohmann::json &value)
  {
    return value.is_number() ? fmt::format("{:.3f} %", 100.0 * value.get<double>())
                             : std::string("nan");
  };
  // A factor / range entry is null when not finite (review m-6: never .get<double>() on
  // it).
  auto Factor = [](const nlohmann::json &value)
  { return value.is_number() ? fmt::format("{:.4f}", value.get<double>()) : "nan"; };
  std::string banner;
  const auto &thresholds = metadata.at("Thresholds");
  banner += fmt::format(
      "\nSurface-response coverage (interim thresholds: BuildTarget {:g}, "
      "MaxUncoveredShare {:g}; a coverage status is about coverage only): {:d} cell(s);",
      thresholds.at("BuildTarget").get<double>(),
      thresholds.at("MaxUncoveredShare").get<double>(), metadata.at("Cells").get<int>());
  for (const auto &[status, count] : metadata.at("StatusCounts").items())
  {
    banner += fmt::format(" {} {:d}", status, count.get<int>());
  }
  banner += "; regime";
  for (const auto &[status, count] : metadata.at("RegimeCounts").items())
  {
    banner += fmt::format(" {} {:d}", status, count.get<int>());
  }
  banner += "\n";
  if (metadata.at("Cells").get<int>() == 0)
  {
    const auto reason = metadata.find("NoCellsReason");
    banner += fmt::format(
        "*** COVERAGE: no coverage cell was measured in this run (no corrected "
        "measurement{}): nothing of its corrected output is covered by this record ***\n",
        reason != metadata.end() && reason->is_string()
            ? " (" + reason->get<std::string>() + ")"
            : std::string());
  }
  int lines = 0;
  const auto &flagged = metadata.at("Flagged");
  for (const auto &entry : flagged)
  {
    if (lines == maximum_lines)
    {
      banner += fmt::format("*** COVERAGE: {:d} more flagged cell(s): see "
                            "SurfaceResponse.Coverage.Flagged of palace.json ***\n",
                            static_cast<int>(flagged.size()) - maximum_lines);
      break;
    }
    std::string cell;
    for (const auto &[name, value] : entry.at("Cell").items())
    {
      cell +=
          fmt::format("{}{} = {:g}", cell.empty() ? "" : ", ", name, value.get<double>());
    }
    const std::string status = entry.at("coverage status").get<std::string>();
    const std::string regime_status = entry.at("regime status").get<std::string>();
    std::string factor;
    if (regime_status == "NO-MODEL")
    {
      factor = "no covered part (regime NO-MODEL)";
    }
    else if (regime_status == "NO-RANGE")
    {
      factor = fmt::format("covered factor {} (no validated range: NO-RANGE: computed, NOT "
                           "VALIDATED in this regime)",
                           Factor(entry.at("covered correction factor")));
    }
    else
    {
      factor = fmt::format("covered factor {} (validated [{}, {}] by {} window(s): {})",
                           Factor(entry.at("covered correction factor")),
                           Factor(entry.at("validated factor min")),
                           Factor(entry.at("validated factor max")),
                           entry.value("validated N", 0), regime_status);
    }
    if (status == "NO-ENERGY" || status == "NOT-MEASURED")
    {
      banner += fmt::format("*** COVERAGE {} (interim thresholds): {} {}: {} ***\n", status,
                            cell, entry.at("Type").get<std::string>(),
                            status == "NO-ENERGY"
                                ? "W = 0, nothing within R to cover (shares nan)"
                                : "no coverage partition for this cell (shares nan)");
    }
    else
    {
      banner += fmt::format(
          "*** COVERAGE {} (interim thresholds): {} {}: U = {} of W (Missing {}, method "
          "gap {}, pending {}); {} ***\n",
          status, cell, entry.at("Type").get<std::string>(),
          Percent(entry.at("uncovered raw share of W")),
          Percent(entry.at("missing raw share of W")),
          Percent(entry.at("method-gap raw share of W")),
          Percent(entry.at("pending-model share of W")), factor);
    }
    lines++;
  }
  const auto &untargeted = metadata.at("UntargetedInterfaces");
  if (!untargeted.empty())
  {
    std::string names;
    for (const auto &entry : untargeted)
    {
      names += fmt::format("{}{:d} ({})", names.empty() ? "" : ", ",
                           entry.at("Interface").get<int>(),
                           entry.at("Type").get<std::string>());
    }
    banner += fmt::format("Coverage: {:d} untargeted postprocessing interface(s) {} carry "
                          "their raw energy; not corrected (the Type shares are over the "
                          "targeted interfaces only; see `untargeted raw share of the Type "
                          "interface energy`)\n",
                          static_cast<int>(untargeted.size()), names);
  }
  for (const auto &entry : metadata.at("ConfiguredTargetsNotServed"))
  {
    banner +=
        fmt::format("*** COVERAGE TARGET-NOT-SERVED: configured target interface {:d} "
                    "({}): {}; it is NOT corrected (its raw within-R energy is "
                    "uncovered in its Type's W) ***\n",
                    entry.at("Interface").get<int>(), entry.at("Type").get<std::string>(),
                    entry.at("Reason").get<std::string>());
  }
  const int no_range = metadata.at("RegimeCounts").value("NO-RANGE", 0);
  if (no_range > 0)
  {
    banner +=
        fmt::format("Coverage: {:d} of {:d} cell(s) NO-RANGE (no validated range of the "
                    "covered correction factor for this library and excitation class: "
                    "computed, NOT VALIDATED in this regime)\n",
                    no_range, metadata.at("Cells").get<int>());
  }
  if (const auto census = metadata.find("ClaimEndCensus"); census != metadata.end())
  {
    // Rulings 1714 / 1715: one line per run (zeros included).
    auto Population = [&](const char *name)
    {
      const auto entry = census->find(name);
      return entry == census->end() ? std::string("0")
                                    : fmt::format("{:d} (max share of W {:.3e})",
                                                  entry->at("MaxCount").get<int>(),
                                                  entry->at("MaxShareOfW").get<double>());
    };
    banner += fmt::format(
        "Coverage claim-end census: served ends snapped (round-off short) {}, near {}; "
        "kept-class ends round-off short {} (OI-FOOT-TEST-VERTEX-ROUNDOFF: counted, not "
        "changed), near {}; round-off beyond (the foot reaches the segment end: nothing "
        "at stake) served {}, kept {}\n",
        Population("ServedRoundoffShort"), Population("ServedNearNotRoundoff"),
        Population("KeptRoundoffShort"), Population("KeptNearNotRoundoff"),
        Population("ServedRoundoffBeyond"), Population("KeptRoundoffBeyond"));
  }
  const int no_status = metadata.at("ServedNoStatusModels").get<int>();
  if (no_status > 0)
  {
    banner += fmt::format("Coverage: {:d} served model(s) carry no QualificationStatus "
                          "record (information: the `no-status-model share of W` column of "
                          "surface-response-coverage.csv)\n",
                          no_status);
  }
  return banner;
}

std::pair<std::vector<SurfacePostOperator::UncoveredPerimeterPortion>,
          std::vector<SurfacePostOperator::UncoveredPerimeterPortion>>
SplitCoverageUncoveredByKeyClass(
    const std::vector<SurfacePostOperator::UncoveredPerimeterPortion> &uncovered,
    const std::map<int, CoverageKeyRecord> &keys_by_feature)
{
  std::pair<std::vector<SurfacePostOperator::UncoveredPerimeterPortion>,
            std::vector<SurfacePostOperator::UncoveredPerimeterPortion>>
      split;
  for (const auto &portion : uncovered)
  {
    const auto record = keys_by_feature.find(portion.feature);
    const bool method_gap = record != keys_by_feature.end() && record->second.method_gap;
    (method_gap ? split.second : split.first).push_back(portion);
  }
  return split;
}

std::string CoverageStageZeroLine(
    const std::map<int, std::string> &types, const std::set<int> &targets,
    const std::map<int, SurfacePostOperator::InterfacePerimeterLengths> &lengths,
    std::pair<double, double> identification_lengths, const Units &units)
{
  struct TypeLengths
  {
    int targets = 0;
    double total = 0.0, missing = 0.0, method_gap = 0.0;
  };
  std::vector<std::pair<std::string, TypeLengths>> by_type;  // order of first appearance
  for (const auto &[interface, type] : types)
  {
    auto it = std::find_if(by_type.begin(), by_type.end(), [&type = type](const auto &entry)
                           { return entry.first == type; });
    if (it == by_type.end())
    {
      by_type.emplace_back(type, TypeLengths{});
      it = std::prev(by_type.end());
    }
    const auto length = lengths.find(interface);
    if (targets.count(interface) == 0 || length == lengths.end())
    {
      continue;
    }
    it->second.targets++;
    it->second.total += length->second.total;
    it->second.missing += length->second.lists.size() > 0 ? length->second.lists[0] : 0.0;
    it->second.method_gap +=
        length->second.lists.size() > 1 ? length->second.lists[1] : 0.0;
  }
  std::string line =
      "Coverage stage 0: Missing LENGTH share per Type (a LENGTH figure: it is "
      "not the energy share and a small value does not mean good coverage):";
  for (std::size_t i = 0; i < by_type.size(); i++)
  {
    const auto &[type, entry] = by_type[i];
    line += fmt::format("{} {} ", i == 0 ? "" : ",", type);
    if (entry.targets == 0)
    {
      line += "n/a (no target interface)";
    }
    else if (!(entry.total > 0.0))
    {
      line += "n/a (zero perimeter length)";
    }
    else
    {
      line += fmt::format("{:.3f} % (Missing {:.6e} m, method gap {:.6e} m, of {:.6e} m)",
                          100.0 * entry.missing / entry.total,
                          units.Dimensionalize<Units::ValueType::LENGTH>(entry.missing),
                          units.Dimensionalize<Units::ValueType::LENGTH>(entry.method_gap),
                          units.Dimensionalize<Units::ValueType::LENGTH>(entry.total));
    }
  }
  const auto [matched, unmatched] = identification_lengths;
  line += "; identification (all features): unmatched / (matched + unmatched) = ";
  if (matched >= 0.0 && unmatched >= 0.0 && matched + unmatched > 0.0)
  {
    line += fmt::format("{:.3f} %", 100.0 * unmatched / (matched + unmatched));
  }
  else
  {
    line += "n/a (no feature pass in this process: a response-geometry cache load or "
            "explicit models)";
  }
  // ERRATUM-1662 HOOK (lane rescensus; no code here): the length-weighted resolution
  // census of the target edges adds its per-Type stage-0 line beside this one.
  return line + "\n";
}

}  // namespace palace
