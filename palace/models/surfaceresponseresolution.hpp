// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_RESPONSE_RESOLUTION_HPP
#define PALACE_MODELS_SURFACE_RESPONSE_RESOLUTION_HPP

#include <array>
#include <limits>
#include <map>
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
// The m10 / m50 / m90 of a histogram of kResolutionHistogramBins weights (the offline
// census's bin-centre rule: the centre of the first bin where the cumulative weight reaches
// 0.1 / 0.5 / 0.9 of the total; 10^kResolutionHistogramMax if rounding never reaches it;
// NaN without weight). The energy cell's and the PREFLIGHT's one quantile.
std::array<double, 3> ResolutionHistogramQuantiles(const double *histogram);

// The per-cell sums of one measurement: per (Type, plane slot; slot 0 = every plane) the
// energy of the patches of each status (R, M, U, NO-BARS, NE), the m histogram (patches
// with a reading and a positive energy), the patch count and the not-evaluated count, and
// (the flag S, lane v1verdict1) the energy of the patches of spatial-basis models (corner
// and spatial-cluster models) and of spatial-cluster models alone.
// Formed per rank inside the response's per-patch loop, summed once (collective).
struct ResolutionAccumulator
{
  static constexpr int kStatusBins = 5;
  static constexpr int kPatchCount = kStatusBins + kResolutionHistogramBins;
  static constexpr int kNotEvaluatedCount = kPatchCount + 1;
  static constexpr int kSpatialEnergy = kPatchCount + 2;
  static constexpr int kClusterEnergy = kPatchCount + 3;
  static constexpr int kStride = kPatchCount + 4;
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
  // One patch's energy `energy` of Type `type` on plane slot `plane_slot` (>= 1); `spatial`
  // its model has a spatial basis, `cluster` it is a spatial-cluster model.
  void Add(int type, int plane_slot, ResolutionStatus status, double m, double energy,
           bool spatial = false, bool cluster = false);
  // CUT-PLAN step 4 (lane v1table): per target with a Region (its interface index) ONE
  // slot of the same sums over the patches that add to the box, each patch's fixed-trace
  // fabricated energy on that interface counted with its region share s_p (ruling 2251),
  // so that a box's census is its own patches'. Summed apart from `values`
  // (SurfaceResponseOperator::MeasureRegionResolution); empty without Region targets.
  std::map<int, std::vector<double>> regions;
  void AddRegion(int interface, ResolutionStatus status, double m, double energy,
                 bool spatial = false, bool cluster = false);
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
  // The flag S (lane v1verdict1; G19): the shares of the cell's energy on spatial-basis
  // (corner and spatial-cluster) models and on spatial-cluster models (nan: no energy). Not
  // a column of surface-response-resolution.csv.
  double spatial_share = std::numeric_limits<double>::quiet_NaN();
  double cluster_share = std::numeric_limits<double>::quiet_NaN();
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

// The PREFLIGHT (stage-0, geometry only; ruling 1819 F1; ONE verdict, decision 2403 (3),
// lane v1verdict1): per Type and plane the census WEIGHTED BY LENGTH, since no field exists
// before a solve: an EDGE patch (a longitudinal cell of positive length) weighs its cell
// length, a SPATIAL patch (cell {0, 0}: corner and cluster models) as much as edge patches
// with as many knots (its knot count times the slot's edge length per edge knot; by its
// knots alone in a slot without edge patches). length_quantiles = m10 / m50 / m90 of h / R
// over the patches with a reading by that weight (ResolutionHistogramQuantiles; NaN without
// a reading). The word (`resolution`) is the per-iteration table's word of that m90 with no
// convergence evaluated (ParticipationResolutionWord, set by the census operator): never
// RESOLVED. The edge and spatial shares by the census bars stay as numbers.
struct ResolutionPreflightRow
{
  std::string type;
  bool all_planes = true;
  double plane = 0.0;
  long long int edge_patches = 0, spatial_patches = 0;
  double edge_length = 0.0;  // m
  std::array<double, 4> edge_shares{}, spatial_shares{};
  std::array<double, 3> length_quantiles{};
  std::string resolution = "NOT-EVALUATED";
};
// Per (Type, slot): [5 edge-length bins][5 spatial-count bins][edge count][spatial count]
// [edge knots][the edge-length histogram][the spatial-knot histogram] (histograms of
// kResolutionHistogramBins over ResolutionHistogramBin(m), patches with a reading).
inline constexpr int kResolutionPreflightEdgeHistogram = 13;
inline constexpr int kResolutionPreflightSpatialHistogram =
    kResolutionPreflightEdgeHistogram + kResolutionHistogramBins;
inline constexpr int kResolutionPreflightStride =
    kResolutionPreflightSpatialHistogram + kResolutionHistogramBins;
ResolutionPreflightRow ComputeResolutionPreflightRow(const double *slot,
                                                     const std::string &type,
                                                     bool all_planes, double plane,
                                                     bool bars_hold, double length_scale);
// The PREFLIGHT line: per row its word and its length m90, with the rule's bar `bar` (NaN:
// none at `order`) and its UNDER-RESOLVED limit `under_resolved_bar` in the header.
std::string ResolutionPreflightLine(const std::vector<ResolutionPreflightRow> &rows,
                                    bool enabled, int order, double bar,
                                    double under_resolved_bar);

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
