// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfaceparticipationtable.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iterator>
#include <set>
#include <tuple>
#include <fmt/format.h>
#include <mfem.hpp>
#include <nlohmann/json.hpp>
#include "embedded_verdict_bars.hpp"
#include "models/surfaceresponsecoverage.hpp"
#include "models/surfaceresponseidentification.hpp"
#include "models/surfaceresponseresolution.hpp"

namespace palace
{

namespace
{

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// A driven row joins the previous level's row of the same output-frequency index only when
// the two frequencies agree to this relative tolerance.
constexpr double kParticipationFrequencyTolerance = 1.0e-12;

}  // namespace

double ParticipationVerdictBars::CensusBar(int order) const
{
  MFEM_VERIFY(order >= 1, "A census bar needs an order >= 1!");
  return order <= 4 ? census_bar[order - 1] : kNaN;
}

std::string ParticipationVerdictBars::WordSuffix(const std::string &problem) const
{
  if (interim)
  {
    return "-INTERIM";
  }
  return std::find(tested_problems.begin(), tested_problems.end(), problem) !=
                 tested_problems.end()
             ? std::string()
             : std::string("-UNTESTED");
}

bool ParticipationVerdictBars::RegionScopeTested() const
{
  return std::find(tested_scopes.begin(), tested_scopes.end(), "Region") !=
         tested_scopes.end();
}

ParticipationVerdictBars ParseParticipationVerdictBars(const std::string &text,
                                                       const std::string &source,
                                                       bool built_in)
{
  nlohmann::json data;
  try
  {
    data = nlohmann::json::parse(text);
  }
  catch (const std::exception &error)
  {
    MFEM_ABORT("The verdict data \"" << source << "\" is not JSON: " << error.what());
  }
  MFEM_VERIFY(data.is_object() && data.contains("Verdict") && data["Verdict"].is_object(),
              "The verdict data \"" << source << "\" has no \"Verdict\" object!");
  const auto &verdict = data["Verdict"];
  auto Number = [&](const char *key, double low, bool open)
  {
    const auto entry = verdict.find(key);
    MFEM_VERIFY(entry != verdict.end() && entry->is_number() &&
                    (open ? entry->get<double>() > low : entry->get<double>() >= low),
                "The verdict data \"" << source << "\" requires Verdict." << key
                                      << " to be a "
                                      << "number " << (open ? "> " : ">= ") << low << "!");
    return entry->get<double>();
  };
  ParticipationVerdictBars bars;
  bars.min_unknowns_ratio = Number("MinUnknownsRatio", 1.0, true);
  bars.max_census_ratio = Number("MaxCensusRatio", 0.0, true);
  MFEM_VERIFY(bars.max_census_ratio <= 1.0,
              "The verdict data \"" << source
                                    << "\" requires Verdict.MaxCensusRatio <= 1!");
  bars.max_change_per_doubling = Number("MaxChangePerDoublingPercent", 0.0, true);
  bars.under_resolved_factor = Number("UnderResolvedFactor", 1.0, false);
  const auto census = verdict.find("CensusBarOverR");
  MFEM_VERIFY(census != verdict.end() && census->is_object(),
              "The verdict data \"" << source << "\" requires Verdict.CensusBarOverR!");
  for (int order = 1; order <= 4; order++)
  {
    const auto entry = census->find(std::to_string(order));
    MFEM_VERIFY(
        entry != census->end() &&
            (entry->is_null() || (entry->is_number() && entry->get<double>() > 0.0)),
        "The verdict data \"" << source << "\" requires Verdict.CensusBarOverR." << order
                              << " = null or a positive number!");
    bars.census_bar[order - 1] = entry->is_null() ? kNaN : entry->get<double>();
  }
  const auto tested = verdict.find("TestedProblems");
  MFEM_VERIFY(tested != verdict.end() && tested->is_array() &&
                  std::all_of(tested->begin(), tested->end(),
                              [](const nlohmann::json &name) { return name.is_string(); }),
              "The verdict data \""
                  << source << "\" requires Verdict.TestedProblems (problem names)!");
  bars.tested_problems = tested->get<std::vector<std::string>>();
  const auto scopes = verdict.find("TestedScopes");
  // The verdict words of Type rows carry no scope suffix: the data must name "Type" among
  // its tested scopes (and may name "Region").
  MFEM_VERIFY(scopes != verdict.end() && scopes->is_array() &&
                  std::all_of(scopes->begin(), scopes->end(),
                              [](const nlohmann::json &name)
                              {
                                return name.is_string() &&
                                       (name.get<std::string>() == "Type" ||
                                        name.get<std::string>() == "Region");
                              }) &&
                  std::any_of(scopes->begin(), scopes->end(), [](const nlohmann::json &name)
                              { return name.get<std::string>() == "Type"; }),
              "The verdict data \""
                  << source
                  << "\" requires Verdict.TestedScopes naming \"Type\" (and optionally "
                     "\"Region\")!");
  bars.tested_scopes = scopes->get<std::vector<std::string>>();
  const auto interim = verdict.find("Interim");
  MFEM_VERIFY(interim == verdict.end() || interim->is_boolean(),
              "The verdict data \"" << source
                                    << "\" requires Verdict.Interim to be a boolean!");
  bars.interim = interim != verdict.end() && interim->get<bool>();
  const auto provenance = verdict.find("Provenance");
  bars.provenance = provenance != verdict.end() && provenance->is_string()
                        ? provenance->get<std::string>()
                        : std::string();
  bars.built_in = built_in;
  bars.source = source;
  bars.sha256 = Sha256Hex(text);
  return bars;
}

ParticipationVerdictBars BuiltInParticipationVerdictBars()
{
  return ParseParticipationVerdictBars(
      verdict::surface_participation_verdict_json,
      "built-in palace/models/surfaceparticipationverdict.json", true);
}

ParticipationVerdictBars ReadParticipationVerdictBars(const std::string &path)
{
  std::ifstream input(path, std::ios::binary);
  MFEM_VERIFY(input, "Unable to open the resolution bars file \"" << path << "\"!");
  const std::string text((std::istreambuf_iterator<char>(input)),
                         std::istreambuf_iterator<char>());
  nlohmann::json data;
  try
  {
    data = nlohmann::json::parse(text);
  }
  catch (const std::exception &error)
  {
    MFEM_ABORT("The resolution bars file \"" << path << "\" is not JSON: " << error.what());
  }
  if (!data.is_object() || !data.contains("Verdict"))
  {
    return BuiltInParticipationVerdictBars();
  }
  return ParseParticipationVerdictBars(text, path, false);
}

ParticipationVerdictBars ParticipationVerdictBarsFor(const std::string &bars_file)
{
  return bars_file.empty() ? BuiltInParticipationVerdictBars()
                           : ReadParticipationVerdictBars(bars_file);
}

std::string ParticipationProblemName(ProblemType type)
{
  switch (type)
  {
    case ProblemType::ELECTROSTATIC:
      return "Electrostatic";
    case ProblemType::EIGENMODE:
      return "Eigenmode";
    case ProblemType::DRIVEN:
      return "Driven";
    default:
      return {};
  }
}

ParticipationIndicator
ComputeParticipationIndicator(double p, double m90, double unknowns,
                              const std::vector<ParticipationHistoryEntry> &history,
                              const ParticipationVerdictBars &bars)
{
  using W = ParticipationNotEvaluated;
  ParticipationIndicator indicator;
  if (history.empty())
  {
    indicator.not_evaluated = W::FIRST_ITERATION;
    return indicator;
  }
  for (int j = static_cast<int>(history.size()) - 1; j >= 0; j--)
  {
    if (unknowns >= bars.min_unknowns_ratio * history[j].unknowns)
    {
      indicator.anchor = j;
      break;
    }
  }
  if (indicator.anchor < 0)
  {
    indicator.not_evaluated = W::UNKNOWNS_GREW_LITTLE;
    return indicator;
  }
  const auto &anchor = history[indicator.anchor];
  if (anchor.p == 0.0)
  {
    indicator.not_evaluated = W::ZERO_ANCHOR_VALUE;
    return indicator;
  }
  if (!std::isfinite(p) || !std::isfinite(anchor.p))
  {
    indicator.not_evaluated = W::NO_VALUE;
    return indicator;
  }
  if (!(m90 <= bars.max_census_ratio * anchor.m90))
  {
    indicator.not_evaluated = W::CELL_NOT_REFINED;  // also without a census reading
    return indicator;
  }
  indicator.not_evaluated = W::NONE;
  indicator.per_doubling =
      100.0 * std::abs(p / anchor.p - 1.0) / std::log2(unknowns / anchor.unknowns);
  indicator.settled = indicator.per_doubling < bars.max_change_per_doubling;
  return indicator;
}

std::string ParticipationConvergenceName(const ParticipationIndicator &indicator,
                                         const ParticipationVerdictBars &bars,
                                         const std::string &problem)
{
  using W = ParticipationNotEvaluated;
  if (indicator.settled)
  {
    return "SETTLED" + bars.WordSuffix(problem);
  }
  switch (indicator.not_evaluated)
  {
    case W::NONE:
      return "MOVING";
    case W::FIRST_ITERATION:
      return "MOVING: FIRST-ITERATION";
    case W::UNKNOWNS_GREW_LITTLE:
      return fmt::format("MOVING: UNKNOWNS-GREW-<{:g}%",
                         100.0 * (bars.min_unknowns_ratio - 1.0));
    case W::UNMATCHED:
      return "MOVING: UNMATCHED";
    case W::ZERO_ANCHOR_VALUE:
      return "MOVING: ZERO-ANCHOR-VALUE";
    case W::NO_VALUE:
      return "MOVING: NO-VALUE";
    case W::CELL_NOT_REFINED:
      return "MOVING: CELL-NOT-REFINED";
  }
  return "MOVING";
}

std::string ParticipationResolutionWord(double m90, int order, bool settled, bool region,
                                        const ParticipationVerdictBars &bars,
                                        const std::string &problem)
{
  const double bar = bars.CensusBar(order);
  if (!std::isfinite(m90))
  {
    return std::string(ResolutionStatusName(ResolutionStatus::NOT_EVALUATED));
  }
  if (!std::isfinite(bar))
  {
    return std::string(ResolutionStatusName(ResolutionStatus::NO_BARS));
  }
  // A box value's census is its own patches', which the bars were not graded with.
  auto Scoped = [&](std::string word)
  {
    const bool untested =
        region && !bars.RegionScopeTested() && word.find("-UNTESTED") == std::string::npos;
    return untested ? word + "-UNTESTED" : word;
  };
  if (m90 > bars.under_resolved_factor * bar)
  {
    return Scoped(std::string(ResolutionStatusName(ResolutionStatus::UNDER_RESOLVED)));
  }
  if (m90 <= bar && settled)
  {
    return Scoped("RESOLVED" + bars.WordSuffix(problem));
  }
  return Scoped(std::string(ResolutionStatusName(ResolutionStatus::MARGINAL)));
}

std::string ParticipationPreflightWord(double length_m90, int order,
                                       const ParticipationVerdictBars &bars,
                                       const std::string &problem)
{
  return ParticipationResolutionWord(length_m90, order, false, false, bars, problem);
}

void SetParticipationResolution(SurfaceParticipationRow &row, const ResolutionCell *cell)
{
  row.m90 = cell ? cell->quantiles[2] : kNaN;
  row.m50 = cell ? cell->quantiles[1] : kNaN;
  row.spatial_share = cell ? cell->spatial_share : kNaN;
  row.cluster_share = cell ? cell->cluster_share : kNaN;
}

void SetParticipationCoverage(SurfaceParticipationRow &row, const CoverageCell *cell)
{
  if (!cell)
  {
    row.coverage = std::string(CoverageStatusName(CoverageStatus::NOT_MEASURED));
    row.uncovered_share = kNaN;
    row.regime = std::string(RegimeStatusName(RegimeStatus::NO_MODEL));
    return;
  }
  row.coverage = std::string(CoverageStatusName(cell->status));
  row.uncovered_share = cell->uncovered_share + cell->pending_share;
  row.regime = std::string(RegimeStatusName(cell->regime));
}

const ResolutionCell *FindParticipationCensusCell(const std::vector<ResolutionCell> &cells,
                                                  const std::string &type)
{
  const auto cell = std::find_if(cells.begin(), cells.end(), [&](const ResolutionCell &c)
                                 { return c.all_planes && c.type == type; });
  return cell != cells.end() ? &*cell : nullptr;
}

void SetParticipationNotTarget(SurfaceParticipationRow &row)
{
  const std::string not_target(CoverageStatusName(CoverageStatus::NOT_TARGET));
  row.p_corrected = row.q_corrected = kNaN;
  row.resolution = row.coverage = not_target;
  row.m90 = row.m50 = row.uncovered_share = row.spatial_share = row.cluster_share = kNaN;
  row.regime = std::string(RegimeStatusName(RegimeStatus::NO_MODEL));
  row.reason = not_target + ": raw only";
}

std::vector<SurfaceParticipationRow>
PreflightParticipationRows(const std::map<int, std::string> &targets,
                           const std::set<int> &regions,
                           const std::vector<ResolutionPreflightRow> &preflight)
{
  std::vector<SurfaceParticipationRow> rows;
  for (const auto &[interface, type] : targets)
  {
    SurfaceParticipationRow row;
    row.point = 0;
    row.excitation = 0;
    row.interface = interface;
    row.type = type;
    row.region = regions.count(interface) > 0;
    std::string word(ResolutionStatusName(ResolutionStatus::NOT_EVALUATED));
    const auto census = std::find_if(preflight.begin(), preflight.end(),
                                     [&](const ResolutionPreflightRow &r)
                                     { return r.all_planes && r.type == type; });
    if (!row.region && census != preflight.end())
    {
      word = census->resolution;
    }
    row.resolution = word;
    row.coverage = row.regime = row.confidence = row.convergence =
        kParticipationNotComputed;
    row.flags = kParticipationPreflightFlags;
    if (census != preflight.end() && census->UnderResolvedCells())
    {
      row.flags += std::string("; ") + kParticipationPreflightUnderResolvedCellsFlag;
    }
    row.reason = kParticipationPreflightReason;
    rows.push_back(std::move(row));
  }
  return rows;
}

std::map<int, ParticipationModeJoin>
MatchParticipationModes(const std::vector<ParticipationModeSignature> &previous,
                        const std::vector<ParticipationModeSignature> &current)
{
  auto Nearest =
      [](const std::vector<ParticipationModeSignature> &modes, double f, std::size_t skip)
  {
    std::size_t best = modes.size();
    for (std::size_t i = 0; i < modes.size(); i++)
    {
      if (i != skip && (best == modes.size() || std::abs(modes[i].frequency - f) <
                                                    std::abs(modes[best].frequency - f)))
      {
        best = i;
      }
    }
    return best;
  };
  // The cosine of two fingerprints: the interfaces' raw participations and the ports'
  // |EPR|.
  auto Overlap =
      [](const ParticipationModeSignature &a, const ParticipationModeSignature &b)
  {
    double ab = 0.0, aa = 0.0, bb = 0.0;
    auto Add = [&](const std::map<int, double> &x, const std::map<int, double> &y)
    {
      for (const auto &[key, value] : x)
      {
        aa += value * value;
        const auto other = y.find(key);
        ab += other != y.end() ? value * other->second : 0.0;
      }
      for (const auto &[key, value] : y)
      {
        bb += value * value;
      }
    };
    Add(a.participation, b.participation);
    Add(a.ports, b.ports);
    return aa > 0.0 && bb > 0.0 ? ab / std::sqrt(aa * bb) : 0.0;
  };
  auto SameField =
      [&](const ParticipationModeSignature &a, const ParticipationModeSignature &b)
  { return Overlap(a, b) >= kParticipationModeMinOverlap; };
  std::map<int, ParticipationModeJoin> joins;
  for (std::size_t m = 0; m < current.size(); m++)
  {
    auto &join = joins[current[m].point];
    join.point = current[m].point;
    join.frequency = current[m].frequency;
    const std::size_t j = Nearest(previous, current[m].frequency, previous.size());
    if (j == previous.size())
    {
      join.failed = "no previous mode";
      continue;
    }
    join.nearest = previous[j].point;
    join.nearest_frequency = previous[j].frequency;
    join.overlap = Overlap(current[m], previous[j]);
    const double gap = std::abs(current[m].frequency - previous[j].frequency);
    const std::size_t neighbour = Nearest(previous, previous[j].frequency, j);
    if (neighbour != previous.size() &&
        !(gap < 0.5 * std::abs(previous[neighbour].frequency - previous[j].frequency)))
    {
      join.failed = "gap not below half the spacing";
      continue;
    }
    if (Nearest(current, previous[j].frequency, current.size()) != m)
    {
      join.failed = "not mutual";
      continue;
    }
    if (!(join.overlap >= kParticipationModeMinOverlap))
    {
      join.failed = fmt::format("overlap below {:g}", kParticipationModeMinOverlap);
      continue;
    }
    bool unique = true;
    int previous_rank = 0, current_rank = 0;
    for (std::size_t k = 0; k < previous.size(); k++)
    {
      if (k != j && SameField(current[m], previous[k]))
      {
        unique = unique && std::abs(current[m].frequency - previous[k].frequency) >
                               kParticipationModeAmbiguityRatio * gap;
        previous_rank += previous[k].frequency < previous[j].frequency;
      }
    }
    for (std::size_t n = 0; n < current.size(); n++)
    {
      if (n != m && SameField(current[n], previous[j]))
      {
        current_rank += current[n].frequency < current[m].frequency;
      }
    }
    if (!unique)
    {
      join.failed = "not unique";
      continue;
    }
    if (previous_rank != current_rank)
    {
      join.failed = "same-field rank differs";
      continue;
    }
    join.joined = true;
  }
  return joins;
}

std::vector<SurfaceParticipationRow>
AddParticipationTypeTotals(const std::vector<SurfaceParticipationRow> &rows)
{
  const std::string not_target(CoverageStatusName(CoverageStatus::NOT_TARGET));
  // Per (point, excitation) in order of first appearance: its Type totals, then its rows.
  std::vector<std::pair<int, int>> groups;
  std::map<std::pair<int, int>, std::vector<SurfaceParticipationRow>> totals, members;
  for (const auto &row : rows)
  {
    const std::pair<int, int> group = {row.point, row.excitation};
    if (members.find(group) == members.end())
    {
      groups.push_back(group);
    }
    members[group].push_back(row);
    if (row.region || row.type_total || row.coverage == not_target)
    {
      continue;
    }
    auto &list = totals[group];
    auto total = std::find_if(list.begin(), list.end(),
                              [&](const auto &t) { return t.type == row.type; });
    if (total == list.end())
    {
      SurfaceParticipationRow sum = row;
      sum.interface = 0;
      sum.type_total = true;
      sum.p_raw = sum.p_corrected = 0.0;
      sum.q_raw = sum.q_corrected = 0.0;  // the sums of 1 / Q until the end
      sum.reason.clear();
      list.push_back(std::move(sum));
      total = list.end() - 1;
    }
    total->p_raw += row.p_raw;
    total->p_corrected += row.p_corrected;
    total->q_raw += 1.0 / row.q_raw;
    total->q_corrected += 1.0 / row.q_corrected;
    total->reason = total->reason.empty() ? row.reason : total->reason;
  }
  std::vector<SurfaceParticipationRow> out;
  out.reserve(rows.size() + 3 * groups.size());
  for (const auto &group : groups)
  {
    for (auto &total : totals[group])
    {
      total.q_raw = total.q_raw == 0.0 ? mfem::infinity() : 1.0 / total.q_raw;
      total.q_corrected =
          total.q_corrected == 0.0 ? mfem::infinity() : 1.0 / total.q_corrected;
      out.push_back(std::move(total));
    }
    for (auto &row : members[group])
    {
      out.push_back(std::move(row));
    }
  }
  return out;
}

std::map<int, ParticipationModeJoin>
JoinParticipationLevel(std::vector<SurfaceParticipationRow> &rows,
                       const std::vector<SurfaceParticipationRow> *previous,
                       const ParticipationVerdictBars &bars, int order)
{
  using Key = std::tuple<int, int, int, std::string>;  // point, excitation, interface, Type
  std::map<Key, const SurfaceParticipationRow *> previous_rows;
  if (previous)
  {
    for (const auto &row : *previous)
    {
      previous_rows[{row.point, row.excitation, row.interface, row.type}] = &row;
    }
  }
  const bool eigenmode = !rows.empty() && rows.front().problem == "Eigenmode";
  std::map<int, ParticipationModeJoin> mode_joins;
  if (previous && eigenmode)
  {
    auto Signatures = [](const std::vector<SurfaceParticipationRow> &level)
    {
      std::map<int, ParticipationModeSignature> by_point;
      for (const auto &row : level)
      {
        auto &signature = by_point[row.point];
        signature.point = row.point;
        signature.frequency = row.frequency;
        signature.ports = row.mode_ports;
        if (!row.region && !row.type_total && std::isfinite(row.p_raw))
        {
          signature.participation[row.interface] = row.p_raw;
        }
      }
      std::vector<ParticipationModeSignature> signatures;
      for (auto &[point, signature] : by_point)
      {
        (void)point;
        signatures.push_back(std::move(signature));
      }
      return signatures;
    };
    mode_joins = MatchParticipationModes(Signatures(*previous), Signatures(rows));
  }
  const std::string not_target(CoverageStatusName(CoverageStatus::NOT_TARGET));
  std::map<std::tuple<int, int, std::string>, const SurfaceParticipationRow *> type_rows;
  for (auto &row : rows)
  {
    row.change_corrected = row.change_raw = row.per_doubling = kNaN;
    row.history.clear();
    bool matched = false;
    if (previous)
    {
      const SurfaceParticipationRow *match = nullptr;
      const auto mode = mode_joins.find(row.point);
      const bool point_matched =
          !eigenmode || (mode != mode_joins.end() && mode->second.joined);
      if (point_matched)
      {
        const int point = eigenmode ? mode->second.nearest : row.point;
        const auto found =
            previous_rows.find({point, row.excitation, row.interface, row.type});
        match = found != previous_rows.end() ? found->second : nullptr;
      }
      if (match && row.problem == "Driven" &&
          !(std::abs(row.frequency - match->frequency) <=
            kParticipationFrequencyTolerance *
                std::max(std::abs(row.frequency), std::abs(match->frequency))))
      {
        match = nullptr;
      }
      if (match)
      {
        matched = true;
        row.history = match->history;
        row.history.push_back(
            {static_cast<double>(match->unknowns), match->p_corrected, match->m90});
        if (std::isfinite(row.p_corrected) && std::isfinite(match->p_corrected) &&
            match->p_corrected != 0.0)
        {
          row.change_corrected = row.p_corrected / match->p_corrected - 1.0;
        }
        if (std::isfinite(row.p_raw) && std::isfinite(match->p_raw) && match->p_raw != 0.0)
        {
          row.change_raw = row.p_raw / match->p_raw - 1.0;
        }
      }
    }
    auto indicator = ComputeParticipationIndicator(
        row.p_corrected, row.m90, static_cast<double>(row.unknowns), row.history, bars);
    if (previous && !matched)
    {
      indicator.not_evaluated = ParticipationNotEvaluated::UNMATCHED;
    }
    row.per_doubling = indicator.per_doubling;
    row.convergence = ParticipationConvergenceName(indicator, bars, row.problem);
    if (row.coverage != not_target)
    {
      row.resolution = ParticipationResolutionWord(row.m90, order, indicator.settled,
                                                   row.region, bars, row.problem);
    }
    if (row.type_total)
    {
      type_rows[{row.point, row.excitation, row.type}] = &row;
    }
  }
  // A whole interface row carries its Type total's verdicts (the verdict is of that value;
  // its own change stays in its change columns).
  for (auto &row : rows)
  {
    const auto total = type_rows.find({row.point, row.excitation, row.type});
    if (!row.type_total && !row.region && row.coverage != not_target &&
        total != type_rows.end())
    {
      row.per_doubling = total->second->per_doubling;
      row.convergence = total->second->convergence;
      row.resolution = total->second->resolution;
    }
    row.flags = ParticipationFlags(row);
  }
  return mode_joins;
}

std::string ParticipationFlags(const SurfaceParticipationRow &row)
{
  const std::string not_target(CoverageStatusName(CoverageStatus::NOT_TARGET));
  if (row.coverage == not_target)
  {
    return not_target;
  }
  std::vector<std::string> parts;
  if (row.resolution != "RESOLVED")
  {
    parts.push_back(row.resolution);
  }
  if (row.coverage != CoverageStatusName(CoverageStatus::OK))
  {
    parts.push_back(row.coverage);
  }
  if (row.regime != RegimeStatusName(RegimeStatus::IN))
  {
    parts.push_back("regime " + row.regime);
  }
  if (row.confidence != "pass")
  {
    parts.push_back("confidence " + row.confidence);
  }
  if (row.convergence != "SETTLED")
  {
    parts.push_back(row.convergence);
  }
  if (row.marked)
  {
    parts.push_back("developer: MarkNotApplied");
  }
  if (row.sheet_served)
  {
    parts.push_back(kParticipationSheetServedFlag);
  }
  if (row.sheet_not_served)
  {
    parts.push_back(kParticipationSheetNotServedFlag);
  }
  if (row.spatial_share >= kParticipationSpatialShareFlagAt)
  {
    parts.push_back(kParticipationSpatialShareFlag);
  }
  if (row.cluster_share > kParticipationClusterShareFlagAbove)
  {
    parts.push_back(kParticipationClusterShareFlag);
  }
  std::string flags;
  for (const auto &part : parts)
  {
    flags += (flags.empty() ? "" : "; ") + part;
  }
  return flags;
}

std::string ParticipationLogBlock(const std::vector<SurfaceParticipationRow> &level,
                                  int iteration, long long int unknowns,
                                  long long int previous_unknowns,
                                  const ParticipationVerdictBars &bars, int order,
                                  const std::map<int, ParticipationModeJoin> &joins)
{
  // The values of the verdicts: the Type totals and the targeted region values.
  const std::string not_target(CoverageStatusName(CoverageStatus::NOT_TARGET));
  std::vector<SurfaceParticipationRow> rows;
  std::copy_if(level.begin(), level.end(), std::back_inserter(rows),
               [&](const SurfaceParticipationRow &row)
               { return row.type_total || (row.region && row.coverage != not_target); });
  const std::string growth =
      previous_unknowns > 0
          ? fmt::format(" (x{:.3f} over iteration {})",
                        static_cast<double>(unknowns) / previous_unknowns, iteration - 1)
          : std::string();
  const double bar = bars.CensusBar(order);
  const std::string bar_text =
      std::isfinite(bar) ? fmt::format("RESOLVED iff h/R m90 <= {:g} AND SETTLED, "
                                       "UNDER-RESOLVED iff m90 > {:g}",
                                       bar, bars.under_resolved_factor * bar)
                         : std::string("no census bar at this order (NO-BARS)");
  const std::string problem = level.empty() ? std::string() : level.front().problem;
  const std::string tested =
      bars.interim || bars.WordSuffix(problem).empty()
          ? std::string()
          : fmt::format("; NOT YET TESTED on {} values (words -UNTESTED)", problem);
  const bool regions =
      std::any_of(rows.begin(), rows.end(),
                  [](const SurfaceParticipationRow &row) { return row.region; });
  const std::string region_text =
      regions && !bars.RegionScopeTested()
          ? std::string(
                "; region resolution words -UNTESTED (the bars were graded with the "
                "cell's census, not a box's own)")
          : std::string();
  std::string out = fmt::format(
      "\nSurface participation (surface-participation.csv; p corrected = the fixed-trace "
      "value), iteration {:d}: {:d} unknowns{}, Order {:d}\n Verdicts{}: SETTLED iff the "
      "change per doubling of the unknowns (against the latest level with x{:g} fewer, "
      "the value's own census m90 fallen to x{:g} of it) < {:g} %; {}; a single solve "
      "never reads RESOLVED{}{} (data: {}{})\n",
      iteration, unknowns, growth, order, bars.interim ? " (INTERIM)" : "",
      bars.min_unknowns_ratio, bars.max_census_ratio, bars.max_change_per_doubling,
      bar_text, tested, region_text, bars.source,
      bars.sha256.empty() ? "" : ", sha256 " + bars.sha256.substr(0, 16));
  for (const auto &[point, join] : joins)
  {
    out += join.joined
               ? fmt::format(" mode {:d} <- mode {:d} of iteration {:d}: {:.9g} -> {:.9g} "
                             "GHz, overlap {:.6f}\n",
                             point, join.nearest, iteration - 1, join.nearest_frequency,
                             join.frequency, join.overlap)
           : join.nearest < 0
               ? fmt::format(" mode {:d} UNMATCHED ({}): {:.9g} GHz\n", point, join.failed,
                             join.frequency)
               : fmt::format(" mode {:d} UNMATCHED ({}): nearest mode {:d} of iteration "
                             "{:d}: {:.9g} -> {:.9g} GHz, overlap {:.6f}\n",
                             point, join.failed, join.nearest, iteration - 1,
                             join.nearest_frequency, join.frequency, join.overlap);
  }
  auto Percent = [](double value, bool sign)
  {
    return std::isfinite(value) ? (sign ? fmt::format("{:+.3f}", 100.0 * value)
                                        : fmt::format("{:.3f}", value))
                                : std::string("-");
  };
  if (rows.size() <= kParticipationLogRows)
  {
    out += fmt::format(" {:>5} {:>3} {:<7} {:<6} {:>9} {:>12} {:>12} {:>9} {:>9} {:>9}  "
                       "{:<16} {:>8} {:<17} {:<24} {}\n",
                       "point", "exc", "Type", "scope", "interface", "p raw", "p corrected",
                       "corr/raw", "change %", "%/dbl", "resolution", "h/R m90", "coverage",
                       "convergence", "flags");
    for (const auto &row : rows)
    {
      out += fmt::format(
          " {:>5d} {:>3d} {:<7} {:<6} {:>9d} {:>12.5e} {:>12.5e} {:>9.4f} "
          "{:>9} {:>9}  {:<16} {:>8.3f} {:<17} {:<24} {}\n",
          row.point, row.excitation, row.type, row.region ? "region" : "Type",
          row.interface, row.p_raw, row.p_corrected, row.p_corrected / row.p_raw,
          Percent(row.change_corrected, true), Percent(row.per_doubling, false),
          row.resolution, row.m90, row.coverage, row.convergence, row.flags);
    }
    return out;
  }
  // Above kParticipationLogRows values: one line per (Type, whole | region).
  struct Group
  {
    std::size_t rows = 0;
    double raw_min = kNaN, raw_max = kNaN, p_min = kNaN, p_max = kNaN, ratio_min = kNaN,
           ratio_max = kNaN, g_max = kNaN;
    std::map<std::string, int> resolution, coverage, convergence;
  };
  std::map<std::pair<std::string, bool>, Group> groups;
  auto Widen = [](double &low, double &high, double value)
  {
    if (std::isfinite(value))
    {
      low = std::isfinite(low) ? std::min(low, value) : value;
      high = std::isfinite(high) ? std::max(high, value) : value;
    }
  };
  for (const auto &row : rows)
  {
    auto &group = groups[{row.type, row.region}];
    group.rows++;
    Widen(group.raw_min, group.raw_max, row.p_raw);
    Widen(group.p_min, group.p_max, row.p_corrected);
    Widen(group.ratio_min, group.ratio_max, row.p_corrected / row.p_raw);
    double unused = kNaN;
    Widen(unused, group.g_max, row.per_doubling);
    group.resolution[row.resolution]++;
    group.coverage[row.coverage]++;
    group.convergence[row.convergence]++;
  }
  auto Counts = [](const std::map<std::string, int> &counts)
  {
    std::string text;
    for (const auto &[name, count] : counts)
    {
      text += fmt::format("{}{} x{:d}", text.empty() ? "" : ", ", name, count);
    }
    return text;
  };
  out += fmt::format(
      " {:d} values (the full rows are in the file), per Type and Type | region:\n",
      rows.size());
  for (const auto &[key, group] : groups)
  {
    out += fmt::format(
        " {:<7} {:<6} {:d} values: p raw {:.5e} .. {:.5e}, p corrected {:.5e} .. {:.5e}, "
        "corr/raw {:.4f} .. {:.4f}, largest %/dbl {}; resolution {}; coverage {}; "
        "convergence {}\n",
        key.first, key.second ? "region" : "Type", group.rows, group.raw_min, group.raw_max,
        group.p_min, group.p_max, group.ratio_min, group.ratio_max,
        Percent(group.g_max, false), Counts(group.resolution), Counts(group.coverage),
        Counts(group.convergence));
  }
  return out;
}

SurfaceParticipationTable::SurfaceParticipationTable(std::string path, std::string problem_,
                                                     int order_,
                                                     ParticipationVerdictBars bars_)
  : problem(std::move(problem_)), order(order_), bars(std::move(bars_)),
    file(std::move(path))
{
  auto &table = file.table;
  table.insert(Column::Integer("iteration", "iteration"));
  table.insert(Column::Integer("unknowns", "unknowns"));
  table.insert(Column::Text("problem", "problem"));
  table.insert(Column::Integer("point", "point"));
  table.insert("frequency", "f (GHz)");
  table.insert(Column::Integer("excitation", "exc"));
  table.insert(Column::Integer("interface", "interface"));
  table.insert(Column::Text("type", "Type"));
  table.insert(Column::Text("scope", "scope"));
  table.insert("p_raw", "p raw");
  table.insert("p_corrected", "p corrected (fixed-trace)");
  table.insert("ratio", "corrected / raw");
  table.insert("q_raw", "Q raw");
  table.insert("q_corrected", "Q corrected");
  // The label of the words' state: INTERIM data, or a problem type the rule was not tested
  // on (the eigenmode and driven rows), in both headers.
  const std::string label = bars.interim                        ? "INTERIM; "
                            : !bars.WordSuffix(problem).empty() ? "NOT YET TESTED here; "
                                                                : "";
  table.insert(
      Column::Text("resolution", "resolution (" + label + "RESOLVED needs SETTLED)"));
  table.insert("m90", "h/R m90");
  table.insert("m50", "h/R m50");
  table.insert(Column::Text("coverage", "coverage"));
  table.insert("uncovered_share", "uncovered+pending share of W");
  table.insert(Column::Text("regime", "regime"));
  table.insert(Column::Text("confidence", "confidence"));
  table.insert("change_corrected", "change corrected (fraction)");
  table.insert("change_raw", "change raw (fraction)");
  table.insert("per_doubling", "change per doubling (%; an interface row: its Type's)");
  table.insert(Column::Text(
      "convergence", label.empty()
                         ? "convergence"
                         : "convergence (" + label.substr(0, label.size() - 2) + ")"));
  table.insert("spatial_share", "S corner+cluster share (fraction)");
  table.insert(Column::Text("flags", "flags"));
  table.insert(Column::Text("reason", "reason"));
}

void SurfaceParticipationTable::FinishLevel(int iteration, long long int unknowns,
                                            bool converged, MPI_Comm comm)
{
  if (!Mpi::Root(comm))
  {
    current.clear();
    return;
  }
  if (!converged)
  {
    Mpi::Print(comm,
               "\nSurface participation (surface-participation.csv), iteration {:d}: solve "
               "not converged; no rows appended\n",
               iteration);
    current.clear();
    return;
  }
  if (current.empty() && !started)
  {
    return;  // no row yet (e.g. a library-build run): no file
  }
  for (auto &row : current)
  {
    row.iteration = iteration;
    row.unknowns = unknowns;
    row.problem = problem;
  }
  current = AddParticipationTypeTotals(current);
  const auto joins =
      JoinParticipationLevel(current, iterations > 0 ? &previous : nullptr, bars, order);
  AppendRows(current);
  started = true;
  Mpi::Print(comm, "{}",
             ParticipationLogBlock(current, iteration, unknowns,
                                   iterations > 0 ? previous_unknowns : 0, bars, order,
                                   joins));
  rows_written += current.size();
  iterations++;
  previous = std::move(current);
  current.clear();
  previous_unknowns = unknowns;
}

void SurfaceParticipationTable::WritePreflightLevel(
    std::vector<SurfaceParticipationRow> rows, long long int unknowns, MPI_Comm comm)
{
  if (!Mpi::Root(comm))
  {
    return;
  }
  for (auto &row : rows)
  {
    row.iteration = 1;
    row.unknowns = unknowns;
    row.problem = problem;
  }
  rows = AddParticipationTypeTotals(rows);
  AppendRows(rows);
  started = true;
  preflight_only = true;
  rows_written += rows.size();
  iterations = 1;
  Mpi::Print(
      comm,
      "\nSurface participation (surface-participation.csv), PreflightOnly: {:d} "
      "row(s) of iteration 1 ({:d} unknowns), the resolution cells the table's word "
      "of the census PREFLIGHT (never RESOLVED); no solve, no value computed (NaN)\n",
      rows.size(), unknowns);
}

void SurfaceParticipationTable::AppendRows(const std::vector<SurfaceParticipationRow> &rows)
{
  auto &table = file.table;
  for (const auto &row : rows)
  {
    table["iteration"] << static_cast<double>(row.iteration);
    table["unknowns"] << static_cast<double>(row.unknowns);
    table["problem"] << row.problem;
    table["point"] << static_cast<double>(row.point);
    table["frequency"] << row.frequency;
    table["excitation"] << static_cast<double>(row.excitation);
    table["interface"] << static_cast<double>(row.interface);
    table["type"] << row.type;
    table["scope"] << std::string(row.type_total ? "Type"
                                  : row.region   ? "region"
                                                 : "interface");
    table["p_raw"] << row.p_raw;
    table["p_corrected"] << row.p_corrected;
    table["ratio"] << row.p_corrected / row.p_raw;
    table["q_raw"] << row.q_raw;
    table["q_corrected"] << row.q_corrected;
    table["resolution"] << row.resolution;
    table["m90"] << row.m90;
    table["m50"] << row.m50;
    table["coverage"] << row.coverage;
    table["uncovered_share"] << row.uncovered_share;
    table["regime"] << row.regime;
    table["confidence"] << row.confidence;
    table["change_corrected"] << row.change_corrected;
    table["change_raw"] << row.change_raw;
    table["per_doubling"] << row.per_doubling;
    table["convergence"] << row.convergence;
    table["spatial_share"] << row.spatial_share;
    table["flags"] << row.flags;
    table["reason"] << row.reason;
  }
  file.AppendRowsAndClear();
}

nlohmann::json SurfaceParticipationTable::Metadata() const
{
  nlohmann::json block;
  block["File"] = kSurfaceParticipationFile;
  block["Iterations"] = iterations;
  block["Rows"] = rows_written;
  if (preflight_only)
  {
    block["PreflightOnly"] = true;
  }
  auto &verdict = block["Verdict"];
  verdict["MinUnknownsRatio"] = bars.min_unknowns_ratio;
  verdict["MaxCensusRatio"] = bars.max_census_ratio;
  verdict["MaxChangePerDoublingPercent"] = bars.max_change_per_doubling;
  verdict["TestedProblems"] = bars.tested_problems;
  verdict["TestedScopes"] = bars.tested_scopes;
  for (int order_ = 1; order_ <= 4; order_++)
  {
    const double bar = bars.census_bar[order_ - 1];
    verdict["CensusBarOverR"][std::to_string(order_)] =
        std::isfinite(bar) ? nlohmann::json(bar) : nlohmann::json(nullptr);
  }
  verdict["UnderResolvedFactor"] = bars.under_resolved_factor;
  verdict["BuiltIn"] = bars.built_in;
  verdict["Interim"] = bars.interim;
  verdict["Source"] = bars.source;
  verdict["Sha256"] = bars.sha256;
  verdict["Provenance"] = bars.provenance;
  return block;
}

}  // namespace palace
