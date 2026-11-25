// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfaceresponsefilmrule.hpp"

#include <algorithm>
#include <set>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <mfem.hpp>
#include <nlohmann/json.hpp>
#include "utils/communication.hpp"
#include "utils/units.hpp"

namespace palace
{

const char *FilmServiceName(FilmService service)
{
  switch (service)
  {
    case FilmService::FILM_PEC:
      return "PEC";
    case FilmService::FILM_INDUCTIVE:
      return "FILM";
    case FilmService::NOT_DECLARED:
      return "NOT-DECLARED";
    case FilmService::SERVED_BY_OVERRIDE:
      return "SERVED-BY-OVERRIDE";
    case FilmService::NOT_SERVED_BY_OVERRIDE:
      return "NOT-SERVED-BY-OVERRIDE";
    case FilmService::INDUCTIVE_ABOVE_RANGE:
      return "ABOVE-RANGE";
  }
  return "NOT-DECLARED";
}

namespace
{

// The sheet form of a metal attribute in the mesh (an inductive sheet is sometimes cracked
// and sometimes not): the per-face law differs between the forms for one configured Ls
// (surfaceimpedanceoperator.cpp: each boundary element of the attribute carries 1 / (s Ls),
// s = 2 on a legacy split): uncracked, one element carries 1 / Ls once (one sheet of total
// Ls); cracked, the two crack faces carry 1 / Ls each; the legacy split, 1 / (2 Ls) each.
std::string SheetFormName(const config::BoundaryData &boundaries, int attribute)
{
  if (boundaries.split_impedance_attributes.contains(attribute))
  {
    return "cracked, legacy split (two faces of 2 Ls each)";
  }
  if (boundaries.cracked_attributes.contains(attribute))
  {
    return "cracked (two faces of Ls each)";
  }
  return "uncracked (one sheet of total Ls)";
}

}  // namespace

const char *ImpedanceLawStatusName(ImpedanceLawStatus status)
{
  switch (status)
  {
    case ImpedanceLawStatus::PEC:
      return "PEC";
    case ImpedanceLawStatus::EXACT:
      return "EXACT";
    case ImpedanceLawStatus::SUBSTITUTED:
      return "SUBSTITUTED";
    case ImpedanceLawStatus::MIXED_LAW:
      return "MIXED-LAW";
    case ImpedanceLawStatus::LAW_JOINT:
      return "LAW-JOINT";
    case ImpedanceLawStatus::NOT_DECLARED:
      return "NOT-DECLARED";
    case ImpedanceLawStatus::NOT_SERVED_BY_OVERRIDE:
      return "NOT-SERVED-BY-OVERRIDE";
    case ImpedanceLawStatus::SERVED_BY_OVERRIDE:
      return "SERVED-BY-OVERRIDE";
  }
  return "NOT-DECLARED";
}

FilmService
FilmRuleResolution::SegmentService(const std::vector<int> &metal_attributes) const
{
  bool all_served = !metal_attributes.empty();
  bool any_override_refusal = false;
  bool any_override_service = false;
  bool any_inductive = false;
  bool all_pec = !metal_attributes.empty();
  for (const int attribute : metal_attributes)
  {
    const auto it = attributes.find(attribute);
    const FilmService service =
        it == attributes.end() ? FilmService::NOT_DECLARED : it->second.service;
    all_served = all_served && FilmServiceServes(service);
    all_pec = all_pec && service == FilmService::FILM_PEC;
    any_override_refusal =
        any_override_refusal || service == FilmService::NOT_SERVED_BY_OVERRIDE;
    any_override_service =
        any_override_service || service == FilmService::SERVED_BY_OVERRIDE;
    any_inductive = any_inductive || service == FilmService::FILM_INDUCTIVE;
  }
  if (any_override_refusal)
  {
    return FilmService::NOT_SERVED_BY_OVERRIDE;
  }
  if (!all_served)
  {
    return FilmService::NOT_DECLARED;
  }
  if (any_override_service)
  {
    return FilmService::SERVED_BY_OVERRIDE;
  }
  return all_pec ? FilmService::FILM_PEC : FilmService::FILM_INDUCTIVE;
}

FilmRuleResolution ResolveImpedanceMetalFilm(
    const config::BoundaryData &boundaries,
    const config::ElectrostaticSolverData::ResponseCorrectionData &request,
    const Units &units, bool nondimensionalized)
{
  using ImpedanceMetalRule =
      config::ElectrostaticSolverData::ResponseCorrectionData::ImpedanceMetalRule;
  FilmRuleResolution resolution;
  const bool apply = request.impedance_metal_rule == ImpedanceMetalRule::APPLY;
  resolution.apply = apply;
  resolution.max_inductance = request.impedance_metal_max_inductance;
  resolution.max_inductance_is_default = !request.impedance_metal_max_inductance_given;
  const double impedance_scale =
      nondimensionalized ? units.GetScaleFactor<Units::ValueType::IMPEDANCE>() : 1.0;
  const double capacitance_scale =
      nondimensionalized ? units.GetScaleFactor<Units::ValueType::CAPACITANCE>() : 1.0;
  auto Add = [&](int attribute, const std::string &law, FilmService service)
  {
    auto &entry = resolution.attributes[attribute];
    entry.attribute = attribute;
    entry.law = law;
    entry.service = service;
  };
  for (const int attribute : boundaries.pec.attributes)
  {
    Add(attribute, "PEC", FilmService::FILM_PEC);
  }
  for (const int attribute : boundaries.auxpec.attributes)
  {
    Add(attribute, "PEC", FilmService::FILM_PEC);
  }
  for (const auto &data : boundaries.conductivity)
  {
    for (const int attribute : data.attributes)
    {
      Add(attribute, "Conductivity", FilmService::NOT_DECLARED);
    }
  }
  for (const auto &data : boundaries.impedance)
  {
    // The configured Ls (as the designer wrote it, H/sq): the scaled Ls times its scale can
    // land one unit in the last place above it (bench-driven-ls NOTES 21).
    const double inductance = nondimensionalized ? data.configured_Ls : data.Ls;
    // Apply (decision 2206 (a); 2208 (3); 2215): a purely inductive sheet is served when
    // its CONFIGURED Ls (the number of the Impedance boundary, whatever the sheet form; the
    // per-face meaning is printed with the form) lies in (0, ImpedanceMetalMaxInductance];
    // above it it is not served; Ls = 0 with Rs = Cs = 0 is no sheet (the solve refuses
    // such a boundary): not served, named. Record serves any Ls (the measuring mode of
    // stage A).
    const bool no_sheet = apply && data.Rs == 0.0 && data.Cs == 0.0 && data.Ls == 0.0;
    const bool inductive = data.Rs == 0.0 && data.Cs == 0.0 && (!apply || inductance > 0.0);
    const bool in_range = !apply || inductance <= resolution.max_inductance;
    const std::string law =
        fmt::format("Impedance Rs {:.6g} Ls {:.6g} Cs {:.6g}", data.Rs * impedance_scale,
                    inductance, data.Cs * capacitance_scale);
    for (const int attribute : data.attributes)
    {
      Add(attribute, law,
          inductive ? (in_range ? FilmService::FILM_INDUCTIVE
                                : FilmService::INDUCTIVE_ABOVE_RANGE)
                    : FilmService::NOT_DECLARED);
      auto &entry = resolution.attributes[attribute];
      entry.inductance = inductive ? inductance : 0.0;
      entry.form = SheetFormName(boundaries, attribute);
      entry.no_sheet = no_sheet;
    }
  }
  for (const auto &data : boundaries.rational_impedance)
  {
    for (const int attribute : data.attributes)
    {
      Add(attribute, "RationalImpedance", FilmService::NOT_DECLARED);
    }
  }

  std::vector<std::string> override_lines;
  for (const auto &entry : request.impedance_metal_film_overrides)
  {
    for (const int attribute : entry.attributes)
    {
      auto it = resolution.attributes.find(attribute);
      MFEM_VERIFY(it != resolution.attributes.end(),
                  "ImpedanceMetalFilmOverrides: attribute "
                      << attribute
                      << " is not a metal boundary of this run (PEC, AuxPEC, Conductivity, "
                         "Impedance or RationalImpedance)!");
      auto &record = it->second;
      const bool film_by_default = record.service == FilmService::FILM_PEC ||
                                   record.service == FilmService::FILM_INDUCTIVE;
      record.overridden = true;
      if (entry.serve)
      {
        if (film_by_default)
        {
          override_lines.push_back(fmt::format(
              "  override: attribute {:d} ({}) served: it is the library's film by default "
              "(no change)",
              attribute, record.law));
        }
        else
        {
          if (apply)
          {
            // Decision 2215: the refusal follows from the configuration alone, so every
            // rank takes this branch and stops; the root prints the reason, the other ranks
            // one short line.
            if (Mpi::Root(Mpi::World()))
            {
              MFEM_ABORT(
                  "ImpedanceMetalFilmOverrides: attribute "
                  << attribute << " (" << record.law
                  << ") cannot be served by an override under \"ImpedanceMetalRule\": "
                     "\"Apply\": "
                  << (record.service == FilmService::INDUCTIVE_ABOVE_RANGE
                          ? fmt::format("its Ls is above the Ls range served "
                                        "(ImpedanceMetalMaxInductance = {:.6g} H/sq); a "
                                        "larger ImpedanceMetalMaxInductance serves it",
                                        resolution.max_inductance)
                          : std::string(record.no_sheet
                                            ? "it is no sheet (Rs, Ls and Cs are all 0)"
                                            : "Apply serves perfect conductors and purely "
                                              "inductive sheets (Rs = 0, Cs = 0) only"))
                  << "!");
            }
            MFEM_ABORT("ImpedanceMetalFilmOverrides: an override is refused under "
                       "\"ImpedanceMetalRule\": \"Apply\" (the reason is printed by the "
                       "root process)!");
          }
          record.service = FilmService::SERVED_BY_OVERRIDE;
          override_lines.push_back(fmt::format(
              "  override: attribute {:d} ({}) served by override: not a purely inductive "
              "sheet; only the measured check stands behind it",
              attribute, record.law));
        }
      }
      else
      {
        override_lines.push_back(fmt::format(
            "  override: attribute {:d} ({}) NOT served: marked by the designer as not the "
            "library's film",
            attribute, record.law));
        record.service = FilmService::NOT_SERVED_BY_OVERRIDE;
      }
    }
  }

  // The printed lines: the film line (the laws treated as the library's film and their
  // attributes; under Apply with each sheet's form, the Ls range and the Ls values met),
  // the overrides, then one line per law not served.
  std::map<std::string, std::vector<int>> film_laws, not_served_laws;
  std::set<double> served_inductances;
  for (const auto &[attribute, record] : resolution.attributes)
  {
    const bool sheet = record.law.rfind("Impedance", 0) == 0;
    const std::string law = apply && sheet && !record.form.empty()
                                ? record.law + ", " + record.form
                                : record.law;
    if (FilmServiceServes(record.service))
    {
      film_laws[law + (record.service == FilmService::SERVED_BY_OVERRIDE
                           ? " (served by override)"
                           : "")]
          .push_back(attribute);
      if (record.service == FilmService::FILM_INDUCTIVE)
      {
        served_inductances.insert(record.inductance);
      }
    }
    else
    {
      not_served_laws[law + (record.service == FilmService::NOT_SERVED_BY_OVERRIDE
                                 ? " (not served by override)"
                                 : (record.service == FilmService::INDUCTIVE_ABOVE_RANGE
                                        ? fmt::format(": ABOVE the Ls range <= {:.6g} H/sq",
                                                      resolution.max_inductance)
                                        : (record.no_sheet
                                               ? ": no sheet (Rs, Ls and Cs are all 0)"
                                               : "")))]
          .push_back(attribute);
    }
    if (record.service == FilmService::FILM_INDUCTIVE ||
        record.service == FilmService::SERVED_BY_OVERRIDE ||
        record.service == FilmService::NOT_SERVED_BY_OVERRIDE)
    {
      resolution.candidate_differs = true;
    }
  }
  resolution.served_inductances.assign(served_inductances.begin(),
                                       served_inductances.end());
  std::string film;
  for (const auto &[law, attributes] : film_laws)
  {
    film += fmt::format("{}{} [{}]", film.empty() ? "" : "; ", law,
                        fmt::join(attributes, ", "));
  }
  if (apply)
  {
    std::vector<std::string> met;
    for (const double inductance : resolution.served_inductances)
    {
      met.push_back(fmt::format("{:.6g}", inductance));
    }
    resolution.lines.push_back(fmt::format(
        "Impedance-metal rule (APPLY): treated as the library's perfect-conductor film (an "
        "ASSUMPTION): {}; Ls range served: 0 < configured Ls <= {:.6g} H/sq in every "
        "sheet "
        "form (ImpedanceMetalMaxInductance, {}); served Ls met: {}",
        film.empty() ? "none" : film, resolution.max_inductance,
        resolution.max_inductance_is_default ? "the default" : "configured",
        met.empty() ? "none" : fmt::format("{} H/sq", fmt::join(met, ", "))));
    // CUT-PLAN step 3 (decision 2208 (1) as ruled for the default): the regime label that
    // replaced the line "experimental and not validated": the bound and the validated
    // range.
    resolution.lines.push_back(fmt::format(
        "Impedance-metal rule (APPLY): regime: purely inductive sheets at configured Ls <= "
        "{:.6g} H/sq are served by the perfect-conductor models; {}",
        resolution.max_inductance, kImpedanceApplyValidatedRange));
  }
  else
  {
    resolution.lines.push_back(fmt::format(
        "Impedance-metal rule (RECORD): film = the library's film: an ASSUMPTION: "
        "{}",
        film.empty() ? "none" : film));
  }
  resolution.lines.insert(resolution.lines.end(), override_lines.begin(),
                          override_lines.end());
  for (const auto &[law, attributes] : not_served_laws)
  {
    // Not served by the perfect-conductor models; raw unless the library holds a model of
    // the device law itself (the feature's law status EXACT: corrected by that model). The
    // line precedes the matching, so it names both outcomes.
    resolution.lines.push_back(fmt::format(
        "  not served ({}): attributes [{}]{}", law, fmt::join(attributes, ", "),
        apply ? fmt::format(": not served by the perfect-conductor models: raw unless the "
                            "library holds a model of this law (law status EXACT, "
                            "corrected by that model); the perfect-conductor models serve "
                            "Ls sheets at or below {:.6g} H/sq ({})",
                            resolution.max_inductance, kImpedanceApplyValidatedRange)
              : std::string()));
  }
  return resolution;
}

bool ImpedanceApplyActs(
    const config::BoundaryData &boundaries,
    const config::ElectrostaticSolverData::ResponseCorrectionData &request)
{
  using ImpedanceMetalRule =
      config::ElectrostaticSolverData::ResponseCorrectionData::ImpedanceMetalRule;
  return request.impedance_metal_rule == ImpedanceMetalRule::APPLY &&
         (!boundaries.impedance.empty() || !request.impedance_metal_film_overrides.empty());
}

bool ImpedanceLawStatusServesSheet(const std::string &status)
{
  return status == ImpedanceLawStatusName(ImpedanceLawStatus::SUBSTITUTED) ||
         status == ImpedanceLawStatusName(ImpedanceLawStatus::MIXED_LAW) ||
         status == ImpedanceLawStatusName(ImpedanceLawStatus::SERVED_BY_OVERRIDE);
}

ImpedanceApplyCounts CountImpedanceApplyFeatures(
    const std::vector<
        config::ElectrostaticSolverData::ResponseCorrectionData::ImpedanceLawFeatureData>
        &features)
{
  ImpedanceApplyCounts counts;
  for (const auto &feature : features)
  {
    counts.sheet += feature.impedance_sheet ? 1 : 0;
    if (feature.status == ImpedanceLawStatusName(ImpedanceLawStatus::NOT_DECLARED) ||
        feature.status ==
            ImpedanceLawStatusName(ImpedanceLawStatus::NOT_SERVED_BY_OVERRIDE) ||
        feature.status == ImpedanceLawStatusName(ImpedanceLawStatus::LAW_JOINT))
    {
      counts.not_served++;
    }
    if (!ImpedanceLawStatusServesSheet(feature.status))
    {
      continue;
    }
    counts.substituted++;
    if (feature.model.empty())
    {
      counts.missing++;
    }
    else
    {
      counts.corrected++;
      counts.corrected_length += feature.length;
    }
  }
  return counts;
}

std::string ImpedanceApplyCountText(const ImpedanceApplyCounts &counts)
{
  return fmt::format("{:d} of {:d} feature(s) on impedance sheets are corrected by "
                     "perfect-conductor models ({:d} Missing, as for a perfect conductor)",
                     counts.corrected, counts.sheet, counts.missing);
}

nlohmann::json ImpedanceApplyRecord(
    const FilmRuleResolution &film,
    const std::vector<
        config::ElectrostaticSolverData::ResponseCorrectionData::ImpedanceLawFeatureData>
        &features,
    const Units &units, bool with_features)
{
  std::map<std::string, int> counts;
  nlohmann::json rows = nlohmann::json::array();
  for (const auto &feature : features)
  {
    counts[feature.status]++;
    if (with_features)
    {
      rows.push_back({{"Id", feature.feature},
                      {"Type", feature.type},
                      {"DeviceLaw", feature.device_law},
                      {"LawStatus", feature.status},
                      {"Override", feature.override_status},
                      {"DeviceKeyHash", feature.device_hash},
                      {"ServedKeyHash", feature.served_hash},
                      {"Model", feature.model}});
    }
  }
  const auto sheet_counts = CountImpedanceApplyFeatures(features);
  nlohmann::json record = {
      {"Mode", "Apply"},
      {"Lines", film.lines},
      {"MaxInductance", film.max_inductance},
      {"MaxInductanceIsDefault", film.max_inductance_is_default},
      {"ServedInductances", film.served_inductances},
      {"LawStatusCounts", counts},
      {"ImpedanceSheetFeatures", sheet_counts.sheet},
      {"SubstitutedImpedanceFeatures", sheet_counts.substituted},
      {"ServedImpedanceFeatures", sheet_counts.corrected},
      {"MissingImpedanceFeatures", sheet_counts.missing},
      {"ServedImpedanceLength",
       units.Dimensionalize<Units::ValueType::LENGTH>(sheet_counts.corrected_length)}};
  if (with_features)
  {
    record["Features"] = std::move(rows);
  }
  return record;
}

}  // namespace palace
