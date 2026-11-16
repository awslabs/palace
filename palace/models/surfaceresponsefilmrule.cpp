// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfaceresponsefilmrule.hpp"

#include <algorithm>
#include <set>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <mfem.hpp>
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
  }
  return "NOT-DECLARED";
}

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
  FilmRuleResolution resolution;
  const double impedance_scale =
      nondimensionalized ? units.GetScaleFactor<Units::ValueType::IMPEDANCE>() : 1.0;
  const double inductance_scale =
      nondimensionalized ? units.GetScaleFactor<Units::ValueType::INDUCTANCE>() : 1.0;
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
    const bool inductive = data.Rs == 0.0 && data.Cs == 0.0;
    const std::string law =
        fmt::format("Impedance Rs {:.6g} Ls {:.6g} Cs {:.6g}", data.Rs * impedance_scale,
                    data.Ls * inductance_scale, data.Cs * capacitance_scale);
    for (const int attribute : data.attributes)
    {
      Add(attribute, law,
          inductive ? FilmService::FILM_INDUCTIVE : FilmService::NOT_DECLARED);
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
  // attributes), the overrides, then one line per law not served.
  std::map<std::string, std::vector<int>> film_laws, not_served_laws;
  for (const auto &[attribute, record] : resolution.attributes)
  {
    if (FilmServiceServes(record.service))
    {
      film_laws[record.law + (record.service == FilmService::SERVED_BY_OVERRIDE
                                  ? " (served by override)"
                                  : "")]
          .push_back(attribute);
    }
    else
    {
      not_served_laws[record.law + (record.service == FilmService::NOT_SERVED_BY_OVERRIDE
                                        ? " (not served by override)"
                                        : "")]
          .push_back(attribute);
    }
    if (record.service == FilmService::FILM_INDUCTIVE ||
        record.service == FilmService::SERVED_BY_OVERRIDE ||
        record.service == FilmService::NOT_SERVED_BY_OVERRIDE)
    {
      resolution.candidate_differs = true;
    }
  }
  std::string film;
  for (const auto &[law, attributes] : film_laws)
  {
    film += fmt::format("{}{} [{}]", film.empty() ? "" : "; ", law,
                        fmt::join(attributes, ", "));
  }
  resolution.lines.push_back(fmt::format(
      "Impedance-metal rule (RECORD): film = the library's film: an ASSUMPTION: "
      "{}",
      film.empty() ? "none" : film));
  resolution.lines.insert(resolution.lines.end(), override_lines.begin(),
                          override_lines.end());
  for (const auto &[law, attributes] : not_served_laws)
  {
    resolution.lines.push_back(fmt::format("  not served ({}): attributes [{}]", law,
                                           fmt::join(attributes, ", ")));
  }
  return resolution;
}

}  // namespace palace
