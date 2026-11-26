// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfaceresponselibrarycheck.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <set>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include "models/surfaceparticipationtable.hpp"
#include "utils/enum_string.hpp"
#include "utils/iodata.hpp"

namespace palace
{

namespace
{

std::string Percent(double fraction)
{
  return fmt::format("{:+.3g} %", 100.0 * fraction);
}

std::string Number(double value)
{
  return fmt::format("{:.10g}", value);
}

bool Equal(double config, double library)
{
  return std::abs(config - library) <=
         kLibraryMismatchEqualTolerance * std::max(std::abs(config), std::abs(library));
}

}  // namespace

// THE MEASURED FACTORS (lane v1loadcheck, reads/eps-cl-moves.tsv; cluster runs on the
// container of main fe67940f30): the merge tier's electrostatic member (the CTX-002 flat
// window FC1v2cg3, Order 3, 3 terminals x 6 target interfaces) with ONE value of the
// configuration moved by +0.175 % and +1 %, the library's models unchanged; effect = |the
// move of the corrected value - the move of the raw value| (the raw value carries the
// configuration's value everywhere, as a library made for it would), its largest over the
// rows of the Type, over the mismatch, rounded up to two decimals; linear in the mismatch.
double LibraryMismatchFactor(LibraryMismatchQuantity quantity, InterfaceDielectric layer,
                             InterfaceDielectric affected)
{
  const auto Of = [affected](double ma, double ms, double sa)
  {
    return affected == InterfaceDielectric::MA   ? ma
           : affected == InterfaceDielectric::MS ? ms
           : affected == InterfaceDielectric::SA ? sa
                                                 : 0.0;
  };
  switch (quantity)
  {
    case LibraryMismatchQuantity::LAYER_PERMITTIVITY:
      return layer == affected ? Of(0.90, 0.94, 0.87) : 0.0;  // 0.896 / 0.939 / 0.865
    case LibraryMismatchQuantity::LAYER_THICKNESS:
      return layer == affected ? Of(0.91, 0.95, 0.94) : 0.0;  // 0.905 / 0.940 / 0.934
    case LibraryMismatchQuantity::SUBSTRATE_PERMITTIVITY:
      return Of(0.16, 1.86, 0.11);  // 0.157 / 1.858 / 0.102
  }
  return 0.0;
}

namespace
{

std::string TargetEntry(const LibraryCheckTarget &target)
{
  return fmt::format("Boundaries.Postprocessing.Dielectric Index {} ({}; attributes {})",
                     target.index, ToString(target.type),
                     fmt::join(target.attributes, ", "));
}

// The tolerance of one quantity alone as the rule reads it (the limit over the largest
// factor of that quantity), printed beside each row.
std::string AloneTolerance(LibraryMismatchQuantity quantity, InterfaceDielectric layer)
{
  double factor = 0.0;
  for (const auto affected :
       {InterfaceDielectric::MA, InterfaceDielectric::MS, InterfaceDielectric::SA})
  {
    factor = std::max(factor, LibraryMismatchFactor(quantity, layer, affected));
  }
  return fmt::format("effect = {:.2f} x |mismatch| on {}; alone within {:.3g} %", factor,
                     quantity == LibraryMismatchQuantity::SUBSTRATE_PERMITTIVITY
                         ? std::string("the most affected Type")
                         : std::string(ToString(layer)),
                     100.0 * kLibraryMismatchEffectLimit / factor);
}

}  // namespace

std::vector<LibraryCheckTarget>
LibraryCheckTargets(const IoData &iodata,
                    const config::ElectrostaticSolverData::ResponseCorrectionData &request)
{
  const double scale = iodata.units.GetMeshLengthRelativeScale();
  const std::set<int> filter(request.target_interfaces.begin(),
                             request.target_interfaces.end());
  std::vector<LibraryCheckTarget> targets;
  for (const auto &[index, dielectric] : iodata.boundaries.postpro.dielectric)
  {
    if ((!filter.empty() && filter.count(index) == 0) ||
        dielectric.type == InterfaceDielectric::DEFAULT ||
        dielectric.edge_distances.empty())
    {
      continue;
    }
    targets.push_back({index, dielectric.type, dielectric.attributes, dielectric.t * scale,
                       dielectric.epsilon_r, dielectric.edge_distances.back() * scale});
  }
  return targets;
}

namespace
{

// The metres of a Fabrication.LengthUnit word (0: not a unit this check knows).
double LengthUnitMetres(const std::string &unit)
{
  const std::map<std::string, double> units = {
      {"m", 1.0}, {"cm", 1.0e-2}, {"mm", 1.0e-3}, {"um", 1.0e-6}, {"nm", 1.0e-9}};
  const auto it = units.find(unit);
  return it != units.end() ? it->second : 0.0;
}

LibraryCheckRow StopRow(std::string quantity, std::string entry, std::string config_value,
                        std::string library_value, std::string tolerance = "exact")
{
  LibraryCheckRow row;
  row.kind = LibraryCheckRow::Kind::STOP;
  row.quantity = std::move(quantity);
  row.entry = std::move(entry);
  row.config_value = std::move(config_value);
  row.library_value = std::move(library_value);
  row.tolerance = std::move(tolerance);
  return row;
}

// A row judged by the rule: the relative mismatch and its effect on the listed targets.
LibraryCheckRow EffectRow(LibraryMismatchQuantity quantity, InterfaceDielectric layer,
                          std::string name, std::string entry, double config,
                          double library, const std::string &library_key,
                          const std::vector<LibraryCheckTarget> &affected)
{
  LibraryCheckRow row;
  row.kind = LibraryCheckRow::Kind::EFFECT;
  row.quantity = std::move(name);
  row.entry = std::move(entry);
  row.config_number = config;
  row.library_number = library;
  row.config_value = Number(config);
  row.library_value = Number(library) + " (" + library_key + ")";
  row.mismatch = config / library - 1.0;
  row.tolerance = AloneTolerance(quantity, layer);
  for (const auto &target : affected)
  {
    row.effect[target.index] =
        LibraryMismatchFactor(quantity, layer, target.type) * std::abs(row.mismatch);
  }
  return row;
}

}  // namespace

std::vector<LibraryCheckRow>
FindLibraryLoadRows(const IoData &iodata,
                    const config::ElectrostaticSolverData::ResponseCorrectionData &request,
                    const std::vector<LibraryCheckTarget> &targets,
                    const LibraryCheckLibrary &library)
{
  std::vector<LibraryCheckRow> rows;
  const std::string library_entry = "Library \"" + request.library + "\"";
  if (!library.refusal.empty())
  {
    rows.push_back(StopRow("library file", library_entry, "-", library.refusal, "-"));
    return rows;
  }
  // The target set (the matching's own refusals, unchanged in action).
  if (targets.empty())
  {
    rows.push_back(StopRow("target set", "Boundaries.Postprocessing.Dielectric",
                           "no entry of an explicit Type with EdgeDistances",
                           "at least one target interface", "-"));
  }
  std::set<int> missing(request.target_interfaces.begin(), request.target_interfaces.end());
  // (AutomaticEdges / EdgeAttributes depend on the mesh's dimension: checked with it.)
  for (const auto &target : targets)
  {
    missing.erase(target.index);
  }
  for (const int index : missing)
  {
    rows.push_back(StopRow("target set", fmt::format("TargetInterfaces entry {}", index),
                           "missing, untyped, or without EdgeDistances",
                           "a target interface", "-"));
  }
  // R: the library's models replace the energy inside it, so any other R is another model.
  for (const auto &target : targets)
  {
    if (!Equal(target.radius, library.matching_radius))
    {
      rows.push_back(StopRow("R (the largest EdgeDistances)", TargetEntry(target),
                             Number(target.radius) + " mesh length units",
                             Number(library.matching_radius) + " (MatchingRadius)"));
    }
  }
  if (library.length_unit)
  {
    const double metres = LengthUnitMetres(*library.length_unit);
    if (!Equal(iodata.model.L0, metres))
    {
      rows.push_back(StopRow("length unit", "Model.L0", Number(iodata.model.L0) + " m",
                             "\"" + *library.length_unit + "\" (Fabrication.LengthUnit)"));
    }
  }
  if (!request.developer)
  {
    for (const auto &model : library.parked_models)
    {
      rows.push_back(StopRow(
          "metal law", "model \"" + model + "\"", "the device carries its law",
          "a conductivity or rational-impedance coupon (parked)", "\"Developer\": true"));
    }
  }
  // The layers (decision 1009's InterfaceLayers): a library without them cannot be checked.
  if (library.interface_layers.empty())
  {
    LibraryCheckRow row = StopRow("InterfaceLayers", library_entry, "-",
                                  "no InterfaceLayers metadata: the target thicknesses "
                                  "and permittivities cannot be verified",
                                  "-");
    row.kind = LibraryCheckRow::Kind::WARNING;
    rows.push_back(std::move(row));
    return rows;
  }
  for (const auto &target : targets)
  {
    const std::string type = ToString(target.type);
    const auto layer = library.interface_layers.find(target.type);
    if (layer == library.interface_layers.end())
    {
      rows.push_back(StopRow(type + " InterfaceLayers entry", TargetEntry(target), "-",
                             "none for " + type, "-"));
      continue;
    }
    if (!Equal(target.thickness, layer->second.thickness))
    {
      rows.push_back(EffectRow(
          LibraryMismatchQuantity::LAYER_THICKNESS, target.type, type + " layer thickness",
          TargetEntry(target), target.thickness, layer->second.thickness,
          "Fabrication.InterfaceLayers." + type + ".Thickness", {target}));
    }
    if (!Equal(target.permittivity, layer->second.permittivity))
    {
      rows.push_back(EffectRow(LibraryMismatchQuantity::LAYER_PERMITTIVITY, target.type,
                               type + " layer permittivity", TargetEntry(target),
                               target.permittivity, layer->second.permittivity,
                               "Fabrication.InterfaceLayers." + type + ".Permittivity",
                               {target}));
    }
  }
  return rows;
}

namespace
{

// A material's permittivity against the library's substrate: vacuum (every component 1),
// else the component farthest from the substrate, as config / library - 1.
std::optional<double> SubstrateMismatch(const std::array<double, 3> &s, double substrate)
{
  if (std::all_of(s.begin(), s.end(), [](double e) { return Equal(e, 1.0); }))
  {
    return std::nullopt;
  }
  double mismatch = 0.0;
  for (const double e : s)
  {
    if (!Equal(e, substrate) && std::abs(e / substrate - 1.0) > std::abs(mismatch))
    {
      mismatch = e / substrate - 1.0;
    }
  }
  return mismatch;
}

std::string MaterialValue(const std::array<double, 3> &s)
{
  return s[0] == s[1] && s[1] == s[2]
             ? Number(s[0])
             : fmt::format("{} / {} / {}", Number(s[0]), Number(s[1]), Number(s[2]));
}

double LargestSubstrateFactor()
{
  double factor = 0.0;
  for (const auto affected :
       {InterfaceDielectric::MA, InterfaceDielectric::MS, InterfaceDielectric::SA})
  {
    factor = std::max(factor,
                      LibraryMismatchFactor(LibraryMismatchQuantity::SUBSTRATE_PERMITTIVITY,
                                            InterfaceDielectric::DEFAULT, affected));
  }
  return factor;
}

}  // namespace

std::optional<LibraryCheckRow>
FindSubstrateLoadRow(const IoData &iodata, const LibraryCheckLibrary &library,
                     const std::vector<LibraryCheckTarget> &targets,
                     std::string *candidates)
{
  candidates->clear();
  if (!library.substrate_permittivity || targets.empty())
  {
    return std::nullopt;
  }
  const double substrate = *library.substrate_permittivity;
  // A material could be the substrate if its own effect alone is within the rule.
  const double tolerance = kLibraryMismatchEffectLimit / LargestSubstrateFactor();
  const auto &materials = iodata.domains.materials;
  std::vector<std::string> names;
  std::optional<std::size_t> closest;
  double closest_mismatch = 0.0;
  for (std::size_t m = 0; m < materials.size(); m++)
  {
    const auto mismatch = SubstrateMismatch(materials[m].epsilon_r.s, substrate);
    if (!mismatch)
    {
      continue;
    }
    if (std::abs(*mismatch) <= tolerance)
    {
      names.push_back(fmt::format("Domains.Materials[{}] (attributes {}) {} ({})", m,
                                  fmt::join(materials[m].attributes, ", "),
                                  MaterialValue(materials[m].epsilon_r.s),
                                  Percent(*mismatch)));
    }
    if (!closest || std::abs(*mismatch) < std::abs(closest_mismatch))
    {
      closest = m;
      closest_mismatch = *mismatch;
    }
  }
  *candidates = fmt::format(
      "substrate (checked after the mesh, on the domains next to the corrected metal): "
      "library {} (Fabrication.SubstratePermittivity); could be the substrate: {}",
      Number(substrate),
      names.empty() ? std::string("none") : fmt::format("{}", fmt::join(names, "; ")));
  if (!closest || !names.empty())
  {
    return std::nullopt;
  }
  // No material but vacuum can pass the check after the mesh: judged now, on the closest.
  const auto &s = materials[*closest].epsilon_r.s;
  auto row = EffectRow(
      LibraryMismatchQuantity::SUBSTRATE_PERMITTIVITY, InterfaceDielectric::DEFAULT,
      "substrate permittivity (no material of the configuration within "
      "the rule: the check after the mesh cannot pass)",
      fmt::format("Domains.Materials[{}] (attributes {}; the closest)", *closest,
                  fmt::join(materials[*closest].attributes, ", ")),
      (1.0 + closest_mismatch) * substrate, substrate, "Fabrication.SubstratePermittivity",
      targets);
  row.config_value = MaterialValue(s);
  return row;
}

std::vector<LibraryCheckRow>
SubstrateMeshRows(const IoData &iodata, const std::map<int, std::vector<int>> &mismatches,
                  double substrate_permittivity,
                  const std::vector<LibraryCheckTarget> &targets)
{
  std::vector<LibraryCheckRow> rows;
  for (const auto &[material, attributes] : mismatches)
  {
    const auto &s = iodata.domains.materials[material].epsilon_r.s;
    const auto mismatch = SubstrateMismatch(s, substrate_permittivity);
    if (!mismatch)
    {
      continue;
    }
    auto row = EffectRow(
        LibraryMismatchQuantity::SUBSTRATE_PERMITTIVITY, InterfaceDielectric::DEFAULT,
        "substrate permittivity",
        fmt::format("Domains.Materials[{}] (attributes next to the corrected metal: {})",
                    material, fmt::join(attributes, ", ")),
        (1.0 + *mismatch) * substrate_permittivity, substrate_permittivity,
        "Fabrication.SubstratePermittivity", targets);
    row.config_value = MaterialValue(s);
    rows.push_back(std::move(row));
  }
  return rows;
}

LibraryCheckReport JudgeLibraryRows(std::vector<LibraryCheckRow> rows,
                                    const std::vector<LibraryCheckTarget> &targets,
                                    bool developer)
{
  LibraryCheckReport report;
  report.developer = developer;
  for (const auto &target : targets)
  {
    report.target_type[target.index] = target.type;
  }
  for (const auto &row : rows)
  {
    report.stop = report.stop || row.kind == LibraryCheckRow::Kind::STOP;
    for (const auto &[index, effect] : row.effect)
    {
      report.target_effect[index] += effect;
    }
  }
  for (const auto &[index, effect] : report.target_effect)
  {
    report.beyond_rule = report.beyond_rule || effect > kLibraryMismatchEffectLimit;
  }
  report.stop = report.stop || (report.beyond_rule && !developer);
  report.rows = std::move(rows);
  return report;
}

namespace
{

std::string RowAction(const LibraryCheckRow &row, const LibraryCheckReport &report)
{
  if (row.kind != LibraryCheckRow::Kind::EFFECT)
  {
    return row.kind == LibraryCheckRow::Kind::STOP ? "STOP" : "WARNING";
  }
  const bool beyond = std::any_of(
      row.effect.begin(), row.effect.end(), [&](const auto &entry)
      { return report.target_effect.at(entry.first) > kLibraryMismatchEffectLimit; });
  return !beyond            ? "TOLERATED"
         : report.developer ? "TOLERATED by \"Developer\": true (beyond the rule)"
                            : "STOP (beyond the rule)";
}

std::string RowLine(const LibraryCheckRow &row, const LibraryCheckReport &report)
{
  std::string line =
      fmt::format("    {}: {}; {}; config {}; library {}", RowAction(row, report),
                  row.quantity, row.entry, row.config_value, row.library_value);
  if (row.kind == LibraryCheckRow::Kind::EFFECT)
  {
    double largest = 0.0;
    for (const auto &[index, effect] : row.effect)
    {
      largest = std::max(largest, effect);
    }
    line += fmt::format("; mismatch {}; estimated effect {:.3g} %", Percent(row.mismatch),
                        100.0 * largest);
  }
  return line + "; tolerance: " + row.tolerance + "\n";
}

}  // namespace

std::string LibraryCheckMessage(const LibraryCheckReport &report,
                                const LibraryCheckLibrary &library, const std::string &when,
                                const std::string &substrate_line)
{
  std::string message = fmt::format(
      "Response-correction library check {}: library \"{}\" ({})\n  the rule: a mismatch "
      "is tolerated while its estimated effect on a corrected value stays within {:.3g} % "
      "(effect = the measured factor x |config / library - 1|, summed over the mismatches "
      "acting on one target interface); R, a missing InterfaceLayers entry, the library "
      "file, the length unit, a parked metal law and the target set stop at any "
      "difference\n",
      when, library.name, library.path, 100.0 * kLibraryMismatchEffectLimit);
  if (!library.profile.empty())
  {
    std::vector<std::string> profile;
    for (const auto &[key, value] : library.profile)
    {
      profile.push_back(key + " " + (value ? Number(*value) : std::string("absent")));
    }
    message +=
        fmt::format("  the library's profile (Fabrication; length unit {}, "
                    "degrees): {}\n",
                    library.length_unit.value_or("not stated"), fmt::join(profile, ", "));
  }
  if (!substrate_line.empty())
  {
    message += "  " + substrate_line + "\n";
  }
  if (report.Empty())
  {
    return message + "  no mismatch\n";
  }
  message += fmt::format("  {} row(s):\n", report.rows.size());
  for (const auto &row : report.rows)
  {
    message += RowLine(row, report);
  }
  // Per Type: the largest summed effect of one of its target interfaces.
  std::vector<std::string> per_type;
  for (const auto type :
       {InterfaceDielectric::MA, InterfaceDielectric::MS, InterfaceDielectric::SA})
  {
    std::optional<std::pair<int, double>> largest;
    for (const auto &[index, effect] : report.target_effect)
    {
      if (report.target_type.at(index) == type && (!largest || effect > largest->second))
      {
        largest = std::make_pair(index, effect);
      }
    }
    if (largest && largest->second > 0.0)
    {
      per_type.push_back(fmt::format("{} {:.3g} % (Index {})", ToString(type),
                                     100.0 * largest->second, largest->first));
    }
  }
  if (!per_type.empty())
  {
    message += fmt::format("  estimated effect per Type (the largest sum over its target "
                           "interfaces; the limit {:.3g} %): {}\n",
                           100.0 * kLibraryMismatchEffectLimit, fmt::join(per_type, ", "));
  }
  const bool effect_rows =
      std::any_of(report.rows.begin(), report.rows.end(), [](const LibraryCheckRow &row)
                  { return row.kind == LibraryCheckRow::Kind::EFFECT; });
  if (report.stop)
  {
    message += "  -> STOP: set each entry above to the library's value, or use a library "
               "made for the configuration's values";
    message += report.beyond_rule ? " (\"Developer\": true tolerates an effect beyond the "
                                    "rule, never a STOP row)\n"
                                  : "\n";
  }
  else if (effect_rows)
  {
    message +=
        "  -> TOLERATED: the run goes on; the outputs name it (the flags word \"" +
        std::string(kParticipationLibraryMismatchFlag) +
        "\" of surface-participation.csv, palace.json SurfaceResponse.LibraryMismatch)\n";
  }
  else
  {
    message += "  -> WARNING: the run goes on unverified\n";
  }
  return message;
}

void RecordLibraryMismatch(
    const LibraryCheckReport &report,
    const config::ElectrostaticSolverData::ResponseCorrectionData &request)
{
  auto &record = request.library_mismatch;
  record = {};
  if (report.stop)
  {
    return;
  }
  for (const auto &row : report.rows)
  {
    if (row.kind != LibraryCheckRow::Kind::EFFECT)
    {
      continue;
    }
    config::ElectrostaticSolverData::ResponseCorrectionData::LibraryMismatchRecord::Entry
        entry;
    entry.quantity = row.quantity;
    entry.entry = row.entry;
    entry.config_value = row.config_number;
    entry.library_value = row.library_number;
    entry.mismatch = row.mismatch;
    for (const auto &[index, effect] : row.effect)
    {
      auto &value = entry.effect[ToString(report.target_type.at(index))];
      value = std::max(value, effect);
    }
    record.entries.push_back(std::move(entry));
  }
  for (const auto &[index, effect] : report.target_effect)
  {
    if (effect > 0.0)
    {
      auto &value = record.effect_per_type[ToString(report.target_type.at(index))];
      value = std::max(value, effect);
      record.interfaces.push_back(index);
    }
  }
  record.beyond_rule = report.beyond_rule;
}

nlohmann::json LibraryMismatchJson(
    const config::ElectrostaticSolverData::ResponseCorrectionData::LibraryMismatchRecord
        &record)
{
  if (record.Empty())
  {
    return nlohmann::json::object();
  }
  nlohmann::json json;
  json["Rule"] = fmt::format(
      "tolerated while the estimated effect on a corrected value stays within {:.3g} % "
      "(effect = the measured factor x |config / library - 1|, summed per target "
      "interface); BeyondRule: tolerated beyond it by \"Developer\": true",
      100.0 * kLibraryMismatchEffectLimit);
  json["BeyondRule"] = record.beyond_rule;
  json["Mismatches"] = nlohmann::json::array();
  for (const auto &entry : record.entries)
  {
    json["Mismatches"].push_back({{"Quantity", entry.quantity},
                                  {"Entry", entry.entry},
                                  {"Config", entry.config_value},
                                  {"Library", entry.library_value},
                                  {"Mismatch", entry.mismatch},
                                  {"EstimatedEffect", entry.effect}});
  }
  json["EstimatedEffectPerType"] = record.effect_per_type;
  json["Interfaces"] = record.interfaces;
  return json;
}

config::ElectrostaticSolverData::ResponseCorrectionData::LibraryMismatchRecord
LibraryMismatchFromJson(const nlohmann::json &json)
{
  config::ElectrostaticSolverData::ResponseCorrectionData::LibraryMismatchRecord record;
  if (!json.contains("Mismatches"))
  {
    return record;
  }
  for (const auto &item : json.at("Mismatches"))
  {
    record.entries.push_back(
        {item.at("Quantity").get<std::string>(), item.at("Entry").get<std::string>(),
         item.at("Config").get<double>(), item.at("Library").get<double>(),
         item.at("Mismatch").get<double>(),
         item.at("EstimatedEffect").get<std::map<std::string, double>>()});
  }
  record.effect_per_type =
      json.at("EstimatedEffectPerType").get<std::map<std::string, double>>();
  record.interfaces = json.at("Interfaces").get<std::vector<int>>();
  record.beyond_rule = json.at("BeyondRule").get<bool>();
  return record;
}

std::set<int> LibraryMismatchInterfaces(
    const config::ElectrostaticSolverData::ResponseCorrectionData::LibraryMismatchRecord
        &record)
{
  return {record.interfaces.begin(), record.interfaces.end()};
}

std::set<std::string> LibraryMismatchTypes(
    const config::ElectrostaticSolverData::ResponseCorrectionData::LibraryMismatchRecord
        &record)
{
  std::set<std::string> types;
  for (const auto &[type, effect] : record.effect_per_type)
  {
    types.insert(type);
  }
  return types;
}

}  // namespace palace
