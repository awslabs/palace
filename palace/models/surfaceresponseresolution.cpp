// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfaceresponseresolution.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <fmt/format.h>
#include "models/surfaceresponseidentification.hpp"

namespace palace
{

namespace
{

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

std::string ResolutionNumber(double value)
{
  return std::isnan(value) ? std::string("nan") : fmt::format("{:+.12e}", value);
}

nlohmann::json JsonNumber(double value)
{
  return std::isfinite(value) ? nlohmann::json(value) : nlohmann::json(nullptr);
}

// The status of a group's four shares, or "EMPTY" when the group carries no weight.
std::string GroupStatusName(const double *bins, double total)
{
  if (!(total > 0.0))
  {
    return "EMPTY";
  }
  std::array<double, 4> shares;
  for (int k = 0; k < 4; k++)
  {
    shares[k] = bins[k] / total;
  }
  return std::string(ResolutionStatusName(ResolutionShareStatus(shares)));
}

ResolutionStatus StatusOfName(const std::string &name)
{
  for (const auto status : {ResolutionStatus::RESOLVED, ResolutionStatus::MARGINAL,
                            ResolutionStatus::UNDER_RESOLVED, ResolutionStatus::NO_BARS,
                            ResolutionStatus::NOT_EVALUATED})
  {
    if (ResolutionStatusName(status) == name)
    {
      return status;
    }
  }
  return ResolutionStatus::NOT_EVALUATED;
}

}  // namespace

std::string_view ResolutionStatusName(ResolutionStatus status)
{
  switch (status)
  {
    case ResolutionStatus::RESOLVED:
      return "RESOLVED";
    case ResolutionStatus::MARGINAL:
      return "MARGINAL";
    case ResolutionStatus::UNDER_RESOLVED:
      return "UNDER-RESOLVED";
    case ResolutionStatus::NO_BARS:
      return "NO-BARS";
    case ResolutionStatus::NOT_EVALUATED:
      return "NOT-EVALUATED";
  }
  return "NOT-EVALUATED";
}

int ResolutionStatusRank(ResolutionStatus status)
{
  switch (status)
  {
    case ResolutionStatus::RESOLVED:
      return 0;
    case ResolutionStatus::NOT_EVALUATED:
      return 1;
    case ResolutionStatus::NO_BARS:
      return 2;
    case ResolutionStatus::MARGINAL:
      return 3;
    case ResolutionStatus::UNDER_RESOLVED:
      return 4;
  }
  return 1;
}

double ResolutionBars::B5(int order) const
{
  return b5[std::clamp(order, 1, 4) - 1];
}

double ResolutionBars::B10(int order) const
{
  return b10[std::clamp(order, 1, 4) - 1];
}

ResolutionBars BuiltInResolutionBars()
{
  // rescensus/tools/rescensus_bars_v1.json (sha256 4d6eae3e6c539192...), K0_bars_h_over_R.
  ResolutionBars bars;
  bars.b5 = {kNaN, 0.540763397421779, 0.540763397421779, 1.3277577835971868};
  bars.b10 = {0.2657095439599277, 1.3277577835971868, 2.597499085375387,
              3.2419880492445246};
  bars.radius = kResolutionBuiltInRadius;
  bars.source = "built-in INTERIM (lane rescensus rescensus_bars_v1.json 4d6eae3e; the 2D "
                "map at R = 1.9 um)";
  bars.built_in = true;
  return bars;
}

ResolutionBars ReadResolutionBars(const std::string &path)
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
  MFEM_VERIFY(data.is_object() && data.contains("RadiusMicrons") &&
                  data["RadiusMicrons"].is_number() && data["RadiusMicrons"] > 0.0 &&
                  data.contains("Bars") && data["Bars"].is_object(),
              "The resolution bars file \""
                  << path << "\" requires a positive RadiusMicrons and a Bars object!");
  ResolutionBars bars;
  bars.radius = data["RadiusMicrons"].get<double>() * 1.0e-6;
  for (int order = 1; order <= 4; order++)
  {
    const auto entry = data["Bars"].find(std::to_string(order));
    MFEM_VERIFY(entry != data["Bars"].end() && entry->is_object() &&
                    entry->contains("b10") && (*entry)["b10"].is_number() &&
                    (*entry)["b10"] > 0.0 && entry->contains("b5") &&
                    ((*entry)["b5"].is_null() ||
                     ((*entry)["b5"].is_number() && (*entry)["b5"] > 0.0 &&
                      (*entry)["b5"] <= (*entry)["b10"])),
                "The resolution bars file \""
                    << path << "\" requires Bars." << order
                    << " = {\"b5\": null or 0 < b5 <= b10, \"b10\": a positive number}!");
    bars.b10[order - 1] = (*entry)["b10"].get<double>();
    bars.b5[order - 1] = (*entry)["b5"].is_null() ? kNaN : (*entry)["b5"].get<double>();
  }
  bars.source = path;
  bars.sha256 = Sha256Hex(text);
  bars.built_in = false;
  return bars;
}

bool ResolutionBarsHold(const ResolutionBars &bars, double radius)
{
  return radius > 0.0 && bars.radius > 0.0 &&
         std::abs(radius / bars.radius - 1.0) <= kResolutionRadiusTolerance;
}

double ResolutionLowerMedian(std::vector<double> &values)
{
  if (values.empty())
  {
    return kNaN;
  }
  const auto middle = values.begin() + static_cast<std::ptrdiff_t>((values.size() - 1) / 2);
  std::nth_element(values.begin(), middle, values.end());
  return *middle;
}

double ResolutionElementLongestEdge(mfem::Mesh &mesh, int element)
{
  const mfem::Element *el = mesh.GetElement(element);
  const auto geometry = el->GetGeometryType();
  const int vertex_count = mfem::Geometry::NumVerts[geometry];
  const int dimension = mesh.SpaceDimension();
  std::array<std::array<double, 3>, 8> x{};
  if (mesh.GetNodes())
  {
    // The element transformation at the reference vertices (the offline reader's form).
    const mfem::IntegrationRule *vertices = mfem::Geometries.GetVertices(geometry);
    mfem::ElementTransformation *transformation = mesh.GetElementTransformation(element);
    mfem::Vector point(dimension);
    for (int k = 0; k < vertex_count; k++)
    {
      transformation->Transform(vertices->IntPoint(k), point);
      for (int d = 0; d < dimension; d++)
      {
        x[k][d] = point(d);
      }
    }
  }
  else
  {
    const int *v = el->GetVertices();
    for (int k = 0; k < vertex_count; k++)
    {
      const double *coordinates = mesh.GetVertex(v[k]);
      for (int d = 0; d < dimension; d++)
      {
        x[k][d] = coordinates[d];
      }
    }
  }
  double longest = 0.0;
  for (int edge = 0; edge < el->GetNEdges(); edge++)
  {
    const int *ends = el->GetEdgeVertices(edge);
    double length2 = 0.0;
    for (int d = 0; d < dimension; d++)
    {
      const double delta = x[ends[0]][d] - x[ends[1]][d];
      length2 += delta * delta;
    }
    longest = std::max(longest, std::sqrt(length2));
  }
  return longest;
}

ResolutionStatus ResolutionPatchStatus(double m, int order, const ResolutionBars &bars,
                                       bool bars_hold)
{
  if (std::isnan(m))
  {
    return ResolutionStatus::NOT_EVALUATED;
  }
  if (!bars_hold)
  {
    return ResolutionStatus::NO_BARS;
  }
  const double b5 = bars.B5(order);
  if (!std::isnan(b5) && m <= b5)
  {
    return ResolutionStatus::RESOLVED;
  }
  return m <= bars.B10(order) ? ResolutionStatus::MARGINAL
                              : ResolutionStatus::UNDER_RESOLVED;
}

int ResolutionHistogramBin(double m)
{
  // The offline census's int((log10(max(m, 1e-300)) - LMIN) / (LMAX - LMIN) * NB), clamped.
  const int bin = static_cast<int>(
      (std::log10(std::max(m, 1.0e-300)) - kResolutionHistogramMin) /
      (kResolutionHistogramMax - kResolutionHistogramMin) * kResolutionHistogramBins);
  return std::min(std::max(bin, 0), kResolutionHistogramBins - 1);
}

double ResolutionHistogramCentre(int bin)
{
  return std::pow(10.0, kResolutionHistogramMin +
                            (bin + 0.5) *
                                (kResolutionHistogramMax - kResolutionHistogramMin) /
                                kResolutionHistogramBins);
}

void ResolutionAccumulator::Reset(int types_, int slots_)
{
  types = types_;
  slots = slots_;
  values.assign(static_cast<std::size_t>(types) * slots * kStride, 0.0);
}

void ResolutionAccumulator::Add(int type, int plane_slot, ResolutionStatus status, double m,
                                double energy)
{
  // The status bins in the order R, M, U, NO-BARS, NE; the histogram after them; then the
  // patch count and the not-evaluated count.
  const int bin = status == ResolutionStatus::RESOLVED         ? 0
                  : status == ResolutionStatus::MARGINAL       ? 1
                  : status == ResolutionStatus::UNDER_RESOLVED ? 2
                  : status == ResolutionStatus::NO_BARS        ? 3
                                                               : 4;
  // Slot 0 (every plane), then the patch's plane slot (none when <= 0: no located census).
  for (const int slot : {0, plane_slot > 0 ? plane_slot : -1})
  {
    if (slot < 0)
    {
      continue;
    }
    double *values = Slot(type, slot);
    values[bin] += energy;
    if (status != ResolutionStatus::NOT_EVALUATED && energy > 0.0)
    {
      values[kStatusBins + ResolutionHistogramBin(m)] += energy;
    }
    values[kStride - 2] += 1.0;
    values[kStride - 1] += status == ResolutionStatus::NOT_EVALUATED ? 1.0 : 0.0;
  }
}

ResolutionStatus ResolutionShareStatus(const std::array<double, 4> &shares)
{
  if (shares[3] > kResolutionShareNotEvaluated)
  {
    return ResolutionStatus::NOT_EVALUATED;
  }
  if (shares[0] >= 1.0 - kResolutionShareResolvedMargin)
  {
    return ResolutionStatus::RESOLVED;
  }
  if (shares[2] > kResolutionShareUnder)
  {
    return ResolutionStatus::UNDER_RESOLVED;
  }
  return ResolutionStatus::MARGINAL;
}

ResolutionCell ComputeResolutionCell(const double *slot, const std::string &type,
                                     bool all_planes, double plane, bool bars_hold,
                                     double energy_scale, const std::string &reason)
{
  constexpr int bins = ResolutionAccumulator::kStatusBins;
  constexpr int stride = ResolutionAccumulator::kStride;
  ResolutionCell cell;
  cell.type = type;
  cell.all_planes = all_planes;
  cell.plane = plane;
  cell.patches = std::llround(slot[stride - 2]);
  cell.not_evaluated = std::llround(slot[stride - 1]);
  cell.shares.fill(kNaN);
  cell.quantiles.fill(kNaN);
  cell.energy = kNaN;
  cell.status = ResolutionStatus::NOT_EVALUATED;
  if (!reason.empty())
  {
    cell.reason = reason;
    if (reason == "disabled")
    {
      // No locate: every patch NOT-EVALUATED (the census off counts, never drops).
      double total = 0.0;
      for (int k = 0; k < bins; k++)
      {
        total += slot[k];
      }
      cell.energy = total * energy_scale;
      cell.shares = {0.0, 0.0, 0.0, 1.0};
      cell.not_evaluated = cell.patches;
    }
    return cell;
  }
  // The offline census sums the bins R, M, U, NE in that order (NO-BARS is empty there).
  const double total = slot[0] + slot[1] + slot[2] + slot[4] + slot[3];
  cell.energy = total * energy_scale;
  if (!(total > 0.0))
  {
    cell.reason = "no energy";
    return cell;
  }
  cell.shares[3] = slot[4] / total;
  if (bars_hold)
  {
    cell.shares[0] = slot[0] / total;
    cell.shares[1] = slot[1] / total;
    cell.shares[2] = slot[2] / total;
  }
  // The energy quantiles of m (the offline census's bin-centre rule).
  const double *histogram = slot + bins;
  double histogram_total = 0.0;
  for (int k = 0; k < kResolutionHistogramBins; k++)
  {
    histogram_total += histogram[k];
  }
  if (histogram_total > 0.0)
  {
    const std::array<double, 3> levels = {0.1 * histogram_total, 0.5 * histogram_total,
                                          0.9 * histogram_total};
    double cumulative = 0.0;
    std::size_t next = 0;
    for (int k = 0; k < kResolutionHistogramBins && next < levels.size(); k++)
    {
      cumulative += histogram[k];
      while (next < levels.size() && cumulative >= levels[next])
      {
        cell.quantiles[next++] = ResolutionHistogramCentre(k);
      }
    }
    for (; next < levels.size(); next++)
    {
      cell.quantiles[next] = std::pow(10.0, kResolutionHistogramMax);
    }
  }
  if (cell.shares[3] > kResolutionShareNotEvaluated)
  {
    cell.status = ResolutionStatus::NOT_EVALUATED;
    cell.reason = "more than half of the energy on patches without a reading";
  }
  else if (!bars_hold)
  {
    cell.status = ResolutionStatus::NO_BARS;
    cell.reason = "no bars for this matching radius";
  }
  else
  {
    cell.status = ResolutionShareStatus(cell.shares);
  }
  return cell;
}

std::string ResolutionTableHeader(const std::string &key_header)
{
  return key_header +
         ",Type,plane z (um),patches,patches not evaluated,energy (J),share RESOLVED,share "
         "MARGINAL,share UNDER-RESOLVED,share NOT-EVALUATED,m10 h/R,m50 h/R,m90 h/R,"
         "status,reason\n";
}

void WriteResolutionRows(std::ostream &output, const std::string &keys,
                         const std::vector<ResolutionCell> &cells)
{
  for (const auto &cell : cells)
  {
    output << keys << ',' << cell.type << ','
           << (cell.all_planes ? std::string("all") : fmt::format("{:.1f}", cell.plane))
           << ',' << cell.patches << ',' << cell.not_evaluated << ','
           << ResolutionNumber(cell.energy);
    for (const double share : cell.shares)
    {
      output << ',' << ResolutionNumber(share);
    }
    for (const double quantile : cell.quantiles)
    {
      output << ',' << ResolutionNumber(quantile);
    }
    output << ',' << ResolutionStatusName(cell.status) << ',' << cell.reason << '\n';
  }
}

ResolutionPreflightRow ComputeResolutionPreflightRow(const double *slot,
                                                     const std::string &type,
                                                     bool all_planes, double plane,
                                                     bool bars_hold, double length_scale)
{
  ResolutionPreflightRow row;
  row.type = type;
  row.all_planes = all_planes;
  row.plane = plane;
  row.edge_patches = std::llround(slot[10]);
  row.spatial_patches = std::llround(slot[11]);
  // The group bins in the order R, M, U, NO-BARS, NE; the shares of R, M, U, NE (NO-BARS in
  // the total), the status by the share rule or NO-BARS.
  auto Group = [&](const double *bins, std::array<double, 4> &shares, std::string &status)
  {
    const double total = bins[0] + bins[1] + bins[2] + bins[4] + bins[3];
    shares.fill(kNaN);
    if (!(total > 0.0))
    {
      status = "EMPTY";
      return total;
    }
    const std::array<double, 4> four = {bins[0], bins[1], bins[2], bins[4]};
    shares[3] = bins[4] / total;
    if (bars_hold)
    {
      for (int k = 0; k < 3; k++)
      {
        shares[k] = four[k] / total;
      }
      status = GroupStatusName(four.data(), total);
    }
    else
    {
      status = shares[3] > kResolutionShareNotEvaluated
                   ? std::string(ResolutionStatusName(ResolutionStatus::NOT_EVALUATED))
                   : std::string(ResolutionStatusName(ResolutionStatus::NO_BARS));
    }
    return total;
  };
  row.edge_length = Group(slot, row.edge_shares, row.edge_status) * length_scale;
  Group(slot + 5, row.spatial_shares, row.spatial_status);
  if (row.edge_status == "EMPTY" && row.spatial_status == "EMPTY")
  {
    row.status = ResolutionStatus::NOT_EVALUATED;
  }
  else if (row.edge_status == "EMPTY" || row.spatial_status == "EMPTY")
  {
    row.status =
        StatusOfName(row.edge_status == "EMPTY" ? row.spatial_status : row.edge_status);
  }
  else
  {
    const auto edge = StatusOfName(row.edge_status),
               spatial = StatusOfName(row.spatial_status);
    row.status =
        ResolutionStatusRank(edge) >= ResolutionStatusRank(spatial) ? edge : spatial;
  }
  return row;
}

std::string ResolutionPreflightLine(const std::vector<ResolutionPreflightRow> &rows,
                                    bool enabled)
{
  if (!enabled)
  {
    return "Resolution census PREFLIGHT: the census is DISABLED "
           "(SurfaceResponseCorrection.Resolution.Enabled false): no knot located, every "
           "cell NOT-EVALUATED\n";
  }
  std::string line = "Resolution census PREFLIGHT (geometry only, an indicator; edge "
                     "patches by cell length, spatial patches by count; the worse of the "
                     "two; NOT the offline length-weighted form of record; the energy "
                     "census follows the solve):";
  for (const auto &row : rows)
  {
    line += fmt::format(
        " {} {}: {} (edge {:d} patches, {:.3e} m, beyond b10 {}; spatial {:d} patches, "
        "beyond b10 {});",
        row.type,
        row.all_planes ? std::string("all planes") : fmt::format("z {:.1f} um", row.plane),
        ResolutionStatusName(row.status), row.edge_patches, row.edge_length,
        std::isnan(row.edge_shares[2]) ? std::string("n/a")
                                       : fmt::format("{:.3f}", row.edge_shares[2]),
        row.spatial_patches,
        std::isnan(row.spatial_shares[2]) ? std::string("n/a")
                                          : fmt::format("{:.3f}", row.spatial_shares[2]));
  }
  return line + " RESOLVED is necessary, not sufficient.\n";
}

std::vector<bool>
ResolutionRepeatedPlaneCells(const std::vector<ResolutionRunRecord::Entry> &entries)
{
  std::vector<bool> repeated(entries.size(), false);
  std::set<double> planes;
  for (const auto &entry : entries)
  {
    if (!entry.cell.all_planes)
    {
      planes.insert(entry.cell.plane);
    }
  }
  if (planes.size() != 1)
  {
    return repeated;
  }
  // The same reading: every count, the energy, the shares, the quantiles (NaN equal to
  // NaN), the status and the reason; anything else is a distinct reading (counted).
  auto Same = [](double a, double b) { return (std::isnan(a) && std::isnan(b)) || a == b; };
  auto SameReading = [&](const ResolutionCell &a, const ResolutionCell &b)
  {
    bool same = a.status == b.status && a.reason == b.reason && a.patches == b.patches &&
                a.not_evaluated == b.not_evaluated && Same(a.energy, b.energy);
    for (std::size_t k = 0; k < a.shares.size(); k++)
    {
      same = same && Same(a.shares[k], b.shares[k]);
    }
    for (std::size_t k = 0; k < a.quantiles.size(); k++)
    {
      same = same && Same(a.quantiles[k], b.quantiles[k]);
    }
    return same;
  };
  using CellKey = std::pair<std::vector<std::pair<std::string, double>>, std::string>;
  std::map<CellKey, std::size_t> all_planes_cell;
  for (std::size_t i = 0; i < entries.size(); i++)
  {
    if (entries[i].cell.all_planes)
    {
      all_planes_cell.emplace(CellKey{entries[i].keys, entries[i].cell.type}, i);
    }
  }
  for (std::size_t i = 0; i < entries.size(); i++)
  {
    if (entries[i].cell.all_planes)
    {
      continue;
    }
    const auto it = all_planes_cell.find(CellKey{entries[i].keys, entries[i].cell.type});
    repeated[i] = it != all_planes_cell.end() &&
                  SameReading(entries[i].cell, entries[it->second].cell);
  }
  return repeated;
}

nlohmann::json ResolutionRunMetadata(const ResolutionRunRecord &record, int level,
                                     int levels)
{
  nlohmann::json metadata;
  metadata["Enabled"] = record.enabled;
  metadata["Measure"] =
      "K3: per patch the lower median over its model's basis-point knots of the longest "
      "edge of the leaf element owning each knot (a knot on a face shared by several "
      "elements: the largest of their longest edges), over the matching radius; per cell "
      "the "
      "energy-weighted shares (the fixed-trace fabricated surface energy of the "
      "postprocessed response) of the patch statuses. RESOLVED is necessary, not "
      "sufficient.";
  metadata["Order"] = record.order;
  metadata["Radius"] = JsonNumber(record.radius);
  metadata["Bars"] = {{"Source", record.bars.source},
                      {"Sha256", record.bars.sha256},
                      {"BuiltIn", record.bars.built_in},
                      {"Radius", JsonNumber(record.bars.radius)},
                      {"Hold", record.bars_hold},
                      {"b5", JsonNumber(record.bars.B5(record.order))},
                      {"b10", JsonNumber(record.bars.B10(record.order))}};
  metadata["ShareRule"] = {{"NotEvaluatedAbove", kResolutionShareNotEvaluated},
                           {"ResolvedAtLeast", 1.0 - kResolutionShareResolvedMargin},
                           {"UnderResolvedAbove", kResolutionShareUnder},
                           {"Interim", true}};
  metadata["Knots"] = {{"Patches", record.patches},
                       {"Total", record.knots},
                       {"Located", record.located},
                       {"PatchesNotEvaluated", record.patches_not_evaluated},
                       {"MeshElements", record.mesh_elements}};
  nlohmann::json preflight = nlohmann::json::array();
  for (const auto &row : record.preflight)
  {
    auto shares = [](const std::array<double, 4> &values)
    {
      nlohmann::json list = nlohmann::json::array();
      for (const double value : values)
      {
        list.push_back(JsonNumber(value));
      }
      return list;
    };
    preflight.push_back(
        {{"Type", row.type},
         {"Plane", row.all_planes ? nlohmann::json(nullptr) : nlohmann::json(row.plane)},
         {"EdgePatches", row.edge_patches},
         {"EdgeLength", row.edge_length},
         {"EdgeShares", shares(row.edge_shares)},
         {"EdgeStatus", row.edge_status},
         {"SpatialPatches", row.spatial_patches},
         {"SpatialShares", shares(row.spatial_shares)},
         {"SpatialStatus", row.spatial_status},
         {"Status", std::string(ResolutionStatusName(row.status))}});
  }
  metadata["Preflight"] = std::move(preflight);
  // The summary over EVERY cell (each plane and the all-planes cell; review r1817 MJ-2): a
  // plane carrying a small part of a Type's energy can be UNDER-RESOLVED or NOT-EVALUATED
  // under an all-planes cell that reads RESOLVED.
  struct TypeSummary
  {
    ResolutionStatus worst = ResolutionStatus::RESOLVED;
    std::array<double, 2> m50{kNaN, kNaN}, share_under{kNaN, kNaN};
    int cells = 0, flagged = 0;
  };
  auto Widen = [](std::array<double, 2> &range, double value)
  {
    if (!std::isnan(value))
    {
      range[0] = std::isnan(range[0]) ? value : std::min(range[0], value);
      range[1] = std::isnan(range[1]) ? value : std::max(range[1], value);
    }
  };
  std::map<std::string, int> counts;
  std::map<std::string, TypeSummary> by_type_summary;
  nlohmann::json flagged = nlohmann::json::array();
  int flagged_count = 0;
  // Review r1846 mn-A: with ONE plane in the run, each plane cell repeats the all-planes
  // cell of its measurement and Type (the same patches): the summary counts DISTINCT
  // readings, so such a plane cell (when its reading equals the all-planes cell's) is not
  // counted again. The table keeps every row; RepeatedPlaneCells says how many were not
  // counted.
  const auto repeated = ResolutionRepeatedPlaneCells(record.entries);
  int repeated_count = 0;
  for (std::size_t i = 0; i < record.entries.size(); i++)
  {
    const auto &entry = record.entries[i];
    if (repeated[i])
    {
      repeated_count++;
      continue;
    }
    const auto &cell = entry.cell;
    counts[std::string(ResolutionStatusName(cell.status))]++;
    auto [it, inserted] = by_type_summary.try_emplace(cell.type);
    auto &summary = it->second;
    if (inserted || ResolutionStatusRank(cell.status) > ResolutionStatusRank(summary.worst))
    {
      summary.worst = cell.status;
    }
    summary.cells++;
    Widen(summary.m50, cell.quantiles[1]);
    Widen(summary.share_under, cell.shares[2]);
    if (cell.status != ResolutionStatus::RESOLVED)
    {
      summary.flagged++;
      if (++flagged_count <= kResolutionFlaggedListCap)
      {
        nlohmann::json keys = nlohmann::json::object();
        for (const auto &[name, value] : entry.keys)
        {
          keys[name] = value;
        }
        flagged.push_back({{"Keys", keys},
                           {"Type", cell.type},
                           {"Plane", cell.all_planes ? nlohmann::json(nullptr)
                                                     : nlohmann::json(cell.plane)},
                           {"Status", std::string(ResolutionStatusName(cell.status))},
                           {"ShareUnderResolved", JsonNumber(cell.shares[2])},
                           {"ShareNotEvaluated", JsonNumber(cell.shares[3])},
                           {"m50", JsonNumber(cell.quantiles[1])},
                           {"Reason", cell.reason}});
      }
    }
  }
  metadata["StatusCounts"] = counts;
  nlohmann::json by_type = nlohmann::json::object();
  for (const auto &[type, summary] : by_type_summary)
  {
    by_type[type] = {{"Worst", std::string(ResolutionStatusName(summary.worst))},
                     {"m50Min", JsonNumber(summary.m50[0])},
                     {"m50Max", JsonNumber(summary.m50[1])},
                     {"ShareUnderResolvedMin", JsonNumber(summary.share_under[0])},
                     {"ShareUnderResolvedMax", JsonNumber(summary.share_under[1])},
                     {"Cells", summary.cells},
                     {"Flagged", summary.flagged}};
  }
  metadata["WorstByType"] = std::move(by_type);
  metadata["Flagged"] = std::move(flagged);
  metadata["FlaggedCount"] = flagged_count;
  metadata["RepeatedPlaneCells"] = repeated_count;
  metadata["Cells"] = static_cast<int>(record.entries.size());
  // Review r1817 mn-1 (the coverage block's rule, ruling 1738): a run with no cell says
  // why; without a reason it is not complete.
  if (record.entries.empty())
  {
    metadata["NoCellsReason"] = record.no_cells_reason.empty()
                                    ? nlohmann::json(nullptr)
                                    : nlohmann::json(record.no_cells_reason);
  }
  metadata["Complete"] = !record.entries.empty() || !record.no_cells_reason.empty();
  metadata["Level"] = level;
  metadata["Levels"] = levels;
  return metadata;
}

std::string ResolutionRunBanner(const nlohmann::json &metadata)
{
  if (!metadata.is_object() || !metadata.contains("WorstByType"))
  {
    return "";
  }
  std::string banner;
  const bool built_in = metadata["Bars"].value("BuiltIn", true);
  const bool hold = metadata["Bars"].value("Hold", false);
  for (const auto &[type, value] : metadata["WorstByType"].items())
  {
    const auto number = [&](const char *key)
    {
      return value[key].is_number() ? fmt::format("{:.3g}", value[key].get<double>())
                                    : std::string("n/a");
    };
    banner += fmt::format(
        "\n*** Resolution census (knot-owner size / R, Order {:d}, {}{}): {} {} over every "
        "cell and plane (m50 {} .. {}; share beyond b10 {} .. {}; {:d} of {:d} cell(s) not "
        "RESOLVED: surface-response-resolution.csv, SurfaceResponse.Resolution of "
        "palace.json). RESOLVED is necessary, not sufficient. ***\n",
        metadata.value("Order", 0),
        built_in ? "INTERIM in-sample built-in bars, valid for R = 1.9 um ONLY"
                 : "bars file",
        hold ? "" : ": NOT valid for this radius", type,
        value.value("Worst", std::string("NOT-EVALUATED")), number("m50Min"),
        number("m50Max"), number("ShareUnderResolvedMin"), number("ShareUnderResolvedMax"),
        value.value("Flagged", 0), value.value("Cells", 0));
  }
  if (metadata.value("Cells", 0) == 0)
  {
    const auto reason = metadata.find("NoCellsReason");
    banner += fmt::format(
        "\n*** Resolution census: no census cell was measured in this run{}: no resolution "
        "reading exists for its output ***\n",
        reason != metadata.end() && reason->is_string()
            ? " (" + reason->get<std::string>() + ")"
            : std::string());
  }
  if (!metadata.value("Enabled", true))
  {
    banner += "\n*** Resolution census DISABLED (SurfaceResponseCorrection.Resolution."
              "Enabled false): every cell NOT-EVALUATED ***\n";
  }
  return banner;
}

}  // namespace palace
