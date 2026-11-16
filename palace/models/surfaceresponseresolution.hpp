// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_RESPONSE_RESOLUTION_HPP
#define PALACE_MODELS_SURFACE_RESPONSE_RESOLUTION_HPP

#include <array>
#include <ostream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include <mfem.hpp>
#include <nlohmann/json.hpp>

namespace palace
{

// The edge-resolution census of the surface-response correction (lane rescensus-product,
// decision 1817; the offline measure K3 of lane rescensus, ERRATUM-1662; ruling 1819). Per
// applied patch: the knots = every basis point of its runtime model, placed as the
// correction places them; the owner of a knot = the leaf element containing it (the
// correction's distributed point locator, no fallback search; a knot on a face shared by
// several elements: the one of the LARGEST longest edge, decision 1863, rank-count
// independent); h = the owner's LONGEST EDGE; the patch's h = the LOWER median of h over
// its located knots (NOT-EVALUATED when it has no knot or fewer than half of them located);
// m = h / R, R the library matching radius. Per cell (measurement x Type x plane) the
// energy-weighted shares of the patch statuses, the energy quantiles of m and a status on
// INTERIM bars (the 2D map at R = 1.9 um only). RESOLVED IS NECESSARY, NOT SUFFICIENT: it
// says the mesh at the trace contour is as fine as the interim bars ask, never that a
// corrected value is accurate. Outputs only: no value, no coverage status and no policy
// reads a census status (Warn only).

enum class ResolutionStatus : char
{
  RESOLVED,
  MARGINAL,
  UNDER_RESOLVED,
  NO_BARS,  // no bar applies: the run's radius has no validated bars (fail-closed)
  NOT_EVALUATED
};
std::string_view ResolutionStatusName(ResolutionStatus status);
// The order of "the worse status" (the offline census: UNDER > MARGINAL > NOT-EVALUATED >
// RESOLVED, with NO-BARS between MARGINAL and NOT-EVALUATED).
int ResolutionStatusRank(ResolutionStatus status);

// The INTERIM share rule of a cell (named placeholders of the offline census, fixed):
// NOT-EVALUATED iff s_N > kResolutionShareNotEvaluated; else RESOLVED iff s_R >= 1 -
// kResolutionShareResolvedMargin; else UNDER-RESOLVED iff s_U > kResolutionShareUnder; else
// MARGINAL.
inline constexpr double kResolutionShareNotEvaluated = 0.50;
inline constexpr double kResolutionShareResolvedMargin = 0.10;
inline constexpr double kResolutionShareUnder = 0.50;
// The quantile histogram of the offline census: log10(h / R) in [-3, 3], 600 bins (0.01
// decade: the bin centre lies within a factor 10^0.005 = 1.0116 of any value of its bin).
inline constexpr int kResolutionHistogramBins = 600;
inline constexpr double kResolutionHistogramMin = -3.0, kResolutionHistogramMax = 3.0;
// The radius of the built-in bars (m) and the relative tolerance of a radius match.
inline constexpr double kResolutionBuiltInRadius = 1.9e-6;
inline constexpr double kResolutionRadiusTolerance = 1.0e-6;

// The bars on m = h / R per FE order 1-4 (an order >= 5 reads Order 4's: a placeholder); a
// b5 of NaN = no RESOLVED bar at that order.
struct ResolutionBars
{
  std::array<double, 4> b5{}, b10{};
  double radius = 0.0;  // m: the radius the bars hold for
  std::string source, sha256;
  bool built_in = true;
  double B5(int order) const;
  double B10(int order) const;
};
// The built-in INTERIM bars (lane rescensus rescensus_bars_v1.json 4d6eae3e, from the 2D
// map rvsh map-v1.csv at R = 1.9 um; the doubles of that file).
ResolutionBars BuiltInResolutionBars();
// A bars file: {"RadiusMicrons": R, "Bars": {"1": {"b5": null | x, "b10": x}, ..., "4":
// ...}, "Provenance": text}; every order 1-4 required, b10 > 0, b5 null or in (0, b10].
ResolutionBars ReadResolutionBars(const std::string &path);
// Whether the bars hold for a run of matching radius `radius` (m).
bool ResolutionBarsHold(const ResolutionBars &bars, double radius);

// One patch's reading (knot counts and the lower-median owner edge, mesh units).
struct ResolutionPatchReading
{
  int knots = 0, located = 0;
  double h = 0.0;
  bool Evaluated() const { return knots > 0 && 2 * located >= knots; }
};
// The lower median of `values` (reordered; NaN when empty).
double ResolutionLowerMedian(std::vector<double> &values);
// The longest edge of a mesh element (its geometry's edges between its own vertices; with a
// Nodes grid function, the element transformation at the reference vertices).
double ResolutionElementLongestEdge(mfem::Mesh &mesh, int element);
// The status of a patch of ratio m (NaN = not evaluated) at `order`; `bars_hold` false =
// NO-BARS for every evaluated patch.
ResolutionStatus ResolutionPatchStatus(double m, int order, const ResolutionBars &bars,
                                       bool bars_hold);
int ResolutionHistogramBin(double m);
double ResolutionHistogramCentre(int bin);

// The per-cell sums of one measurement: per (Type, plane slot; slot 0 = every plane) the
// energy of the patches of each status (R, M, U, NO-BARS, NE), the m histogram (patches
// with a reading and a positive energy), the patch count and the not-evaluated count.
// Formed per rank inside the response's per-patch loop, summed once (collective).
struct ResolutionAccumulator
{
  static constexpr int kStatusBins = 5;
  static constexpr int kStride = kStatusBins + kResolutionHistogramBins + 2;
  int types = 0, slots = 0;
  std::vector<double> values;
  void Reset(int types_, int slots_);
  double *Slot(int type, int slot)
  {
    return values.data() + (type * slots + slot) * kStride;
  }
  const double *Slot(int type, int slot) const
  {
    return values.data() + (type * slots + slot) * kStride;
  }
  // One patch's energy `energy` of Type `type` on plane slot `plane_slot` (>= 1).
  void Add(int type, int plane_slot, ResolutionStatus status, double m, double energy);
};

// One row of surface-response-resolution.csv.
struct ResolutionCell
{
  std::string type;
  bool all_planes = true;
  double plane = 0.0;  // um
  long long int patches = 0, not_evaluated = 0;
  double energy = 0.0;  // J
  std::array<double, 4>
      shares{};  // R, M, U, NE of the cell's energy (nan: no bars / energy)
  std::array<double, 3> quantiles{};  // m10, m50, m90 (nan: no reading)
  ResolutionStatus status = ResolutionStatus::NOT_EVALUATED;
  std::string reason;
};
// The cell of one accumulator slot (`energy_scale` J per energy unit); `reason` non-empty =
// every cell of the measurement NOT-EVALUATED for that reason (the census disabled, no
// evaluation-0 energy).
ResolutionCell ComputeResolutionCell(const double *slot, const std::string &type,
                                     bool all_planes, double plane, bool bars_hold,
                                     double energy_scale, const std::string &reason);
// The status of four shares by the interim share rule.
ResolutionStatus ResolutionShareStatus(const std::array<double, 4> &shares);

std::string ResolutionTableHeader(const std::string &key_header);
void WriteResolutionRows(std::ostream &output, const std::string &keys,
                         const std::vector<ResolutionCell> &cells);

// The PREFLIGHT (stage-0, geometry only; ruling 1819 F1): per Type and plane, the EDGE
// patches (a longitudinal cell of positive length) by their cell length and the SPATIAL
// patches (cell {0, 0}) by count; each group's shares and status by the same share rule;
// the printed status = the WORSE of the two (an empty group does not enter; both empty =
// NOT-EVALUATED). Not the offline length-weighted form of record; an indicator only.
struct ResolutionPreflightRow
{
  std::string type;
  bool all_planes = true;
  double plane = 0.0;
  long long int edge_patches = 0, spatial_patches = 0;
  double edge_length = 0.0;  // m
  std::array<double, 4> edge_shares{}, spatial_shares{};
  std::string edge_status, spatial_status;  // a status name, or "EMPTY"
  ResolutionStatus status = ResolutionStatus::NOT_EVALUATED;
};
// Per (Type, slot): [5 edge-length bins][5 spatial-count bins][edge count][spatial count].
inline constexpr int kResolutionPreflightStride = 12;
ResolutionPreflightRow ComputeResolutionPreflightRow(const double *slot,
                                                     const std::string &type,
                                                     bool all_planes, double plane,
                                                     bool bars_hold, double length_scale);
std::string ResolutionPreflightLine(const std::vector<ResolutionPreflightRow> &rows,
                                    bool enabled);

// The run's census record (root only): the measurement cells and the preflight of the
// level, for palace.json SurfaceResponse.Resolution and the end-of-run banner.
struct ResolutionRunRecord
{
  struct Entry
  {
    std::vector<std::pair<std::string, double>> keys;
    ResolutionCell cell;
  };
  bool enabled = true;
  int order = 0;
  double radius = 0.0;  // m
  bool bars_hold = false;
  ResolutionBars bars;
  long long int patches = 0, knots = 0, located = 0, patches_not_evaluated = 0;
  long long int mesh_elements = 0;
  std::vector<ResolutionPreflightRow> preflight;
  std::vector<Entry> entries;
  // Why the run measured no census cell (the coverage record's reason: a raw fallback
  // applies no patch); empty when it has cells. With no cell and no reason the block does
  // not read Complete (the coverage block's rule, ruling 1738).
  std::string no_cells_reason;
};
// palace.json SurfaceResponse.Resolution (Level / Levels as the coverage block; "Final" is
// added by BaseSolver once the AMR loop has ended). StatusCounts, WorstByType, Flagged,
// FlaggedCount and the m50 / share ranges are taken over EVERY distinct cell reading (each
// plane and the all-planes cell): a Type whose energy on one plane is not RESOLVED never
// reads RESOLVED. With one plane in the run a plane cell that repeats its all-planes cell
// is counted once (RepeatedPlaneCells). Flagged lists at most kResolutionFlaggedListCap
// cells; FlaggedCount counts every distinct flagged reading (the table keeps every row, so
// with one plane its flagged rows can exceed FlaggedCount).
inline constexpr int kResolutionFlaggedListCap = 200;
nlohmann::json ResolutionRunMetadata(const ResolutionRunRecord &record, int level,
                                     int levels);
// Per entry: whether it is a plane cell that repeats the all-planes cell of its measurement
// and Type (the run has ONE plane and the two readings are equal in every field).
std::vector<bool>
ResolutionRepeatedPlaneCells(const std::vector<ResolutionRunRecord::Entry> &entries);
// The end-of-run banner: one line per Type (its worst status over every cell, the count of
// cells not RESOLVED, the m50 range and the range of the share beyond b10, the bars'
// source and scope) and "RESOLVED is necessary, not sufficient"; a run without a census
// cell says so, with its reason.
std::string ResolutionRunBanner(const nlohmann::json &metadata);

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_RESPONSE_RESOLUTION_HPP
