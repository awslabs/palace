// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfaceresponsecoverage.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iterator>
#include <limits>
#include <set>
#include <sstream>
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
    case CoverageStatus::NO_ENERGY:
      return 4;
    case CoverageStatus::NOT_MEASURED:
      return 5;
    case CoverageStatus::NOT_TARGET:
      return -1;
  }
  return 5;
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
  MFEM_VERIFY(std::isfinite(regime.tolerance) && regime.tolerance >= 0.0,
              "Coverage ValidatedRegime \"Tolerance\" must be a nonnegative number!");
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
    range.n = entry.value("N", 0);
    MFEM_VERIFY(std::isfinite(range.factor_min) && std::isfinite(range.factor_max) &&
                    range.factor_min <= range.factor_max,
                "Coverage ValidatedRegime range for " << type << " / " << name
                                                      << " needs FactorMin <= FactorMax!");
    MFEM_VERIFY(regime.ranges.emplace(std::make_pair(type, excitation), range).second,
                "Coverage ValidatedRegime lists " << type << " / " << name << " twice!");
  }
  return regime;
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
  for (const auto &[label, energy] : partition.served)
  {
    summary.served += energy;
    bool pending = false, no_status = false;
    for (const auto &token : LabelTokens(label))
    {
      const auto model = models.find(TokenIndex(token));
      if (model != models.end())
      {
        pending = pending || model->second == CoverageModelQualification::PENDING;
        no_status = no_status || model->second == CoverageModelQualification::NO_STATUS;
      }
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
  }
  // Ruling 1641: U is fail-closed, whatever is not served is uncovered.
  cell.uncovered = cell.within - cell.served;
  cell.method_gap = cell.uncovered - cell.missing;
  if (cell.within > 0.0)
  {
    cell.uncovered_share = cell.uncovered / cell.within;
    cell.missing_share = cell.missing / cell.within;
    cell.method_gap_share = cell.method_gap / cell.within;
    cell.pending_share = cell.pending / cell.within;
    cell.no_status_share = cell.no_status / cell.within;
    if (cell.method_gap_share > thresholds.max_uncovered_share)
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
        cell.no_status_share = nan;
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
      cell.no_status = cell.served = cell.fabricated = nan;
  cell.uncovered_share = cell.missing_share = cell.method_gap_share = cell.pending_share =
      cell.no_status_share = nan;
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

}  // namespace palace
