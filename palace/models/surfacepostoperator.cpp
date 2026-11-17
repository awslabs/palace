// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfacepostoperator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstring>
#include <fstream>
#include <limits>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "linalg/vector.hpp"
#include "models/materialoperator.hpp"
#include "models/strattonchu.hpp"
#include "models/surfaceresponseidentification.hpp"
#include "utils/communication.hpp"
#include "utils/edgedistance.hpp"
#include "utils/geodata.hpp"
#include "utils/interfaceownership.hpp"
#include "utils/interfaceregion.hpp"
#include "utils/iodata.hpp"
#include "utils/metaledge.hpp"
#include "utils/ownershipquadrature.hpp"
#include "utils/prettyprint.hpp"
#include "utils/timer.hpp"

namespace palace
{

using namespace std::complex_literals;

namespace
{

bool IsSupportedSurfaceFunctionalDimension(const mfem::ParMesh &mesh)
{
  return (mesh.Dimension() == 2 && mesh.SpaceDimension() == 2) ||
         (mesh.Dimension() == 3 && mesh.SpaceDimension() == 3);
}

bool IsSupportedSurfaceFluxDimension(const mfem::ParMesh &mesh)
{
  // The libCEED surface-flux kernels currently implement 3D surface integrals. 2D line
  // flux semantics remain on the legacy path explicitly rather than by silent invalid
  // fallback.
  return mesh.Dimension() == 3 && mesh.SpaceDimension() == 3;
}

void RequireCeedSurfaceFunctional(const SurfaceFunctional *func, const std::string &what)
{
  MFEM_VERIFY(func && func->IsValid(), "libCEED postprocessing was expected for " + what +
                                           ", but SurfaceFunctional could not assemble!");
}

template <typename T>
mfem::Array<int> SetUpBoundaryProperties(const T &data,
                                         const mfem::Array<int> &bdr_attr_marker)
{
  mfem::Array<int> attr_list;
  attr_list.Reserve(static_cast<int>(data.attributes.size()));
  std::set<int> bdr_warn_list;
  for (auto attr : data.attributes)
  {
    // MFEM_VERIFY(attr > 0 && attr <= bdr_attr_max,
    //             "Boundary postprocessing attribute tags must be non-negative and "
    //             "correspond to attributes in the mesh!");
    // MFEM_VERIFY(bdr_attr_marker[attr - 1],
    //             "Unknown boundary postprocessing attribute " << attr << "!");
    if (attr <= 0 || attr > bdr_attr_marker.Size() || !bdr_attr_marker[attr - 1])
    {
      bdr_warn_list.insert(attr);
    }
    else
    {
      attr_list.Append(attr);
    }
  }
  if (!bdr_warn_list.empty())
  {
    Mpi::Print("\n");
    Mpi::Warning(
        "Unknown boundary postprocessing attributes!\nSolver will just ignore them!");
    utils::PrettyPrint(bdr_warn_list, "Boundary attribute list:");
    Mpi::Print("\n");
  }
  return attr_list;
}

double EdgeDistanceOutsideWeight(double distance, double radius, double smoothing)
{
  if (!std::isfinite(radius))
  {
    return 0.0;
  }
  if (smoothing == 0.0)
  {
    return distance >= radius ? 1.0 : 0.0;
  }

  const double lower = radius * (1.0 - smoothing);
  const double upper = radius * (1.0 + smoothing);
  if (distance <= lower)
  {
    return 0.0;
  }
  if (distance >= upper)
  {
    return 1.0;
  }
  const double x = (distance - lower) / (upper - lower);
  return x * x * (3.0 - 2.0 * x);
}

double EdgeDistanceWindowWeight(double distance, double distance_min, double distance_max,
                                double smoothing)
{
  return std::max(0.0, EdgeDistanceOutsideWeight(distance, distance_min, smoothing) -
                           EdgeDistanceOutsideWeight(distance, distance_max, smoothing));
}

// Physical quadrature weight of a point, including the revolution measure 2 pi x of an
// axisymmetric (r, z) mesh (config::ModelData::axisymmetric).
double PhysicalQuadratureWeight(mfem::ElementTransformation &T,
                                const mfem::IntegrationPoint &ip, bool axisymmetric)
{
  double weight = ip.weight * T.Weight();
  if (axisymmetric)
  {
    mfem::Vector point(T.GetSpaceDim());
    T.Transform(ip, point);
    weight *= 2.0 * M_PI * point[0];
  }
  return weight;
}

// Wraps a coefficient with the revolution measure 2 pi x for MFEM integrators on an
// axisymmetric (r, z) mesh.
class AxisymmetricMeasureCoefficient : public mfem::Coefficient
{
private:
  mfem::Coefficient &coefficient;

public:
  AxisymmetricMeasureCoefficient(mfem::Coefficient &coefficient) : coefficient(coefficient)
  {
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    mfem::Vector point(T.GetSpaceDim());
    T.Transform(ip, point);
    return 2.0 * M_PI * point[0] * coefficient.Eval(T, ip);
  }
};

class OwnershipCoefficient : public mfem::Coefficient
{
private:
  std::unique_ptr<mfem::Coefficient> coefficient;
  std::shared_ptr<const InterfaceOwnershipPartition> ownership;
  int slot;

public:
  OwnershipCoefficient(std::unique_ptr<mfem::Coefficient> coefficient_,
                       std::shared_ptr<const InterfaceOwnershipPartition> ownership_,
                       int slot_)
    : coefficient(std::move(coefficient_)), ownership(std::move(ownership_)), slot(slot_)
  {
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    mfem::Vector point(T.GetSpaceDim());
    T.Transform(ip, point);
    return ownership->SelectSlot(point) == slot ? coefficient->Eval(T, ip) : 0.0;
  }
};

// The quadrature-level spatial filter of an interface entry (decision 352 follow-up (2)):
// the coefficient inside the region, zero outside.
class RegionCoefficient : public mfem::Coefficient
{
private:
  std::unique_ptr<mfem::Coefficient> coefficient;
  std::shared_ptr<const InterfaceRegion> region;

public:
  RegionCoefficient(std::unique_ptr<mfem::Coefficient> coefficient_,
                    std::shared_ptr<const InterfaceRegion> region_)
    : coefficient(std::move(coefficient_)), region(std::move(region_))
  {
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    mfem::Vector point(T.GetSpaceDim());
    T.Transform(ip, point);
    return region->Contains(point) ? coefficient->Eval(T, ip) : 0.0;
  }
};

class EdgeDistanceCoefficient : public mfem::Coefficient
{
private:
  mfem::Coefficient &coefficient;
  const EdgeDistanceTree &edge_distance_tree;
  const double distance_min;
  const double distance_max;
  const double smoothing;

public:
  EdgeDistanceCoefficient(mfem::Coefficient &coefficient,
                          const EdgeDistanceTree &edge_distance_tree, double distance_min,
                          double distance_max, double smoothing)
    : coefficient(coefficient), edge_distance_tree(edge_distance_tree),
      distance_min(distance_min), distance_max(distance_max), smoothing(smoothing)
  {
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    mfem::Vector point(T.GetSpaceDim());
    T.Transform(ip, point);
    const double distance = std::sqrt(edge_distance_tree.DistanceSquared(point));
    return EdgeDistanceWindowWeight(distance, distance_min, distance_max, smoothing) *
           coefficient.Eval(T, ip);
  }
};

// The uncovered intervals of one interface's edge-distance tree (decision 394 F2): per tree
// segment, the closed foot-parameter intervals [a, b] (from the segment's first end) of the
// uncovered portions lying on it, each with its type index into `types`.
// Rulings 1714 / 1715: the end census of one list (and, with `snap`, the served claims'
// round-off-short ends moved to exactly 0 / L). `ends` maps (segment, end 0 | 1) to a mask:
// 1 round-off short, 2 near, 4 (ruling 1731; set only where a round-off-short end of a
// list sharing this check lies) an UNSNAPPED portion of this list has its clamped foot at
// exactly that segment end (the list holds the vertex without the snap).
struct ClaimEndCheck
{
  double radius = 0.0;
  bool snap = false;
  std::map<std::pair<std::size_t, int>, int> ends;
  SurfacePostOperator::InterfaceCoveragePartition::EndCensus census;
};

struct UncoveredEdgeIntervals
{
  struct Interval
  {
    double a = 0.0, b = 0.0;
    int type = -1;
  };
  std::map<std::size_t, std::vector<Interval>> by_segment;
  std::vector<std::string> types;

  UncoveredEdgeIntervals(const EdgeDistanceTree &tree,
                         SurfacePostOperator::PerimeterPortionView portions,
                         double tolerance, ClaimEndCheck *check = nullptr)
  {
    // A portion on this interface's perimeter: its midpoint's tree segment, the segment's
    // length and the clamped foot parameters (from the segment's first end) of its two
    // ends with their distances off the segment.
    struct Located
    {
      bool on_tree = false;
      std::size_t segment = 0;
      double length = 0.0, a0 = 0.0, a1 = 0.0, off0 = 0.0, off1 = 0.0;
    };
    auto Locate = [&](std::size_t i)
    {
      Located located;
      const auto &portion_p0 = portions.P0(i), &portion_p1 = portions.P1(i);
      mfem::Vector midpoint(3), p0(3), p1(3);
      for (int d = 0; d < 3; d++)
      {
        p0(d) = portion_p0[d];
        p1(d) = portion_p1[d];
        midpoint(d) = 0.5 * (portion_p0[d] + portion_p1[d]);
      }
      const auto nearest = tree.Nearest(midpoint);
      if (std::sqrt(nearest.distance_squared) > tolerance)
      {
        return located;  // not on this interface's perimeter (another plane)
      }
      const auto &segment = tree.GetSegment(nearest.segment);
      mfem::Vector origin(3), chord(3);
      for (int d = 0; d < 3; d++)
      {
        origin(d) = segment.p0[d];
        chord(d) = segment.p1[d] - segment.p0[d];
      }
      const double length2 = chord * chord;
      MFEM_VERIFY(length2 > 0.0, "An edge-distance tree segment has zero length!");
      const double length = std::sqrt(length2);
      auto Foot = [&](const mfem::Vector &q)
      {
        mfem::Vector r(q);
        r -= origin;
        const double t = (r * chord) / length2;
        mfem::Vector residual(r);
        residual.Add(-t, chord);
        return std::make_pair(std::clamp(t, 0.0, 1.0) * length, residual.Norml2());
      };
      located.on_tree = true;
      located.segment = nearest.segment;
      located.length = length;
      std::tie(located.a0, located.off0) = Foot(p0);
      std::tie(located.a1, located.off1) = Foot(p1);
      return located;
    };
    bool short_inside = false;
    for (std::size_t i = 0; i < portions.Size(); i++)
    {
      const auto &portion_p0 = portions.P0(i), &portion_p1 = portions.P1(i);
      const std::string &portion_type = portions.Label(i);
      const Located located = Locate(i);
      if (!located.on_tree)
      {
        continue;
      }
      const auto &segment = tree.GetSegment(located.segment);
      const double length = located.length, a0 = located.a0, a1 = located.a1;
      MFEM_VERIFY(located.off0 <= tolerance && located.off1 <= tolerance,
                  "An uncovered perimeter portion of feature "
                      << portions.Feature(i) << " (" << portion_type
                      << ") lies within the tolerance of a perimeter segment at its "
                         "midpoint but not at its ends (ends "
                      << located.off0 << " / " << located.off1
                      << " off the segment): the identification's perimeter differs from "
                         "the postprocessing perimeter!");
      auto type = std::find(types.begin(), types.end(), portion_type);
      if (type == types.end())
      {
        types.push_back(portion_type);
        type = std::prev(types.end());
      }
      double e0 = a0, e1 = a1;
      if (check)
      {
        for (int k = 0; k < 2; k++)
        {
          for (int j = 0; j < 2; j++)
          {
            const auto [match, distance] = SurfacePostOperator::ClassifyClaimEnd(
                k == 0 ? portion_p0 : portion_p1, j == 0 ? segment.p0 : segment.p1,
                check->radius);
            auto &census = check->census;
            if (match == SurfacePostOperator::ClaimEndMatch::ROUNDOFF)
            {
              // Ruling 1731: the direction, on the clamped foot before any snap.
              if ((k == 0 ? a0 : a1) == (j == 0 ? 0.0 : length))
              {
                census.roundoff_beyond_count++;
                census.roundoff_beyond_distance =
                    std::max(census.roundoff_beyond_distance, distance);
                continue;
              }
              census.roundoff_count++;
              census.roundoff_distance = std::max(census.roundoff_distance, distance);
              check->ends[{located.segment, j}] |= 1;
              short_inside = true;
              if (check->snap)
              {
                (k == 0 ? e0 : e1) = j == 0 ? 0.0 : length;
              }
            }
            else if (match == SurfacePostOperator::ClaimEndMatch::NEAR)
            {
              census.near_count++;
              census.near_distance = std::max(census.near_distance, distance);
              check->ends[{located.segment, j}] |= 2;
            }
          }
        }
      }
      by_segment[located.segment].push_back(
          {std::min(e0, e1), std::max(e0, e1), static_cast<int>(type - types.begin())});
    }
    // Ruling 1731: bit 4 at the round-off-short ends this list also reaches unsnapped (a
    // second pass, so that the order of the portions does not matter; nothing is stored
    // for the other ends).
    if (check && (short_inside ||
                  std::any_of(check->ends.begin(), check->ends.end(),
                              [](const auto &end) { return (end.second & 1) != 0; })))
    {
      for (std::size_t i = 0; i < portions.Size(); i++)
      {
        const Located located = Locate(i);
        if (!located.on_tree)
        {
          continue;
        }
        for (int j = 0; j < 2; j++)
        {
          const double segment_end = j == 0 ? 0.0 : located.length;
          const auto found = check->ends.find({located.segment, j});
          if (found != check->ends.end() && (found->second & 1) &&
              (located.a0 == segment_end || located.a1 == segment_end))
          {
            found->second |= 4;
          }
        }
      }
    }
  }
};

// The within-R energy density of the quadrature points whose nearest perimeter point lies
// on an uncovered portion of the given type (every type when `type` < 0): the complement
// of EdgeDistanceCoefficient's outside window times the coefficient, restricted by the
// foot test above.
class UncoveredEdgeCoefficient : public mfem::Coefficient
{
private:
  mfem::Coefficient &coefficient;
  const EdgeDistanceTree &edge_distance_tree;
  const UncoveredEdgeIntervals &intervals;
  const double radius;
  const double smoothing;
  const int type;

public:
  UncoveredEdgeCoefficient(mfem::Coefficient &coefficient,
                           const EdgeDistanceTree &edge_distance_tree,
                           const UncoveredEdgeIntervals &intervals, double radius,
                           double smoothing, int type)
    : coefficient(coefficient), edge_distance_tree(edge_distance_tree),
      intervals(intervals), radius(radius), smoothing(smoothing), type(type)
  {
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    mfem::Vector point(T.GetSpaceDim());
    T.Transform(ip, point);
    const auto nearest = edge_distance_tree.Nearest(point);
    const double within = 1.0 - EdgeDistanceOutsideWeight(
                                    std::sqrt(nearest.distance_squared), radius, smoothing);
    if (within <= 0.0)
    {
      return 0.0;
    }
    const auto it = intervals.by_segment.find(nearest.segment);
    if (it == intervals.by_segment.end())
    {
      return 0.0;
    }
    const auto &segment = edge_distance_tree.GetSegment(nearest.segment);
    double length2 = 0.0, projection = 0.0;
    for (int d = 0; d < point.Size(); d++)
    {
      const double chord = segment.p1[d] - segment.p0[d];
      length2 += chord * chord;
      projection += (point(d) - segment.p0[d]) * chord;
    }
    const double foot = std::clamp(projection / length2, 0.0, 1.0) * std::sqrt(length2);
    for (const auto &interval : it->second)
    {
      if (foot >= interval.a && foot <= interval.b && (type < 0 || interval.type == type))
      {
        return within * coefficient.Eval(T, ip);
      }
    }
    return 0.0;
  }
};

}  // namespace

namespace
{

// PatchEnergy (decision 922 B-1; the 851 (C) export): the exact partition of one
// interface's perimeter into elementary pieces, each labelled by the UNION of the labels of
// the portions covering it (labels are '+'-joined "m<idx>" tokens; the union is the
// ascending token set). Pieces are the sub-intervals between consecutive distinct interval
// ends on a tree segment; a piece is labelled by the portions containing its midpoint
// (closed intervals), so overlapping portions of different labels yield a shared piece
// instead of a double count.
struct PartitionedEdgePieces
{
  struct Piece
  {
    double a = 0.0, b = 0.0;
    int group = -1;
  };
  std::map<std::size_t, std::vector<Piece>> by_segment;  // ascending, non-overlapping
  std::vector<std::string> groups;

  PartitionedEdgePieces(const UncoveredEdgeIntervals &intervals)
  {
    auto Tokens = [](const std::string &label)
    {
      std::set<std::string> tokens;
      std::size_t begin = 0;
      while (begin <= label.size())
      {
        const auto end = label.find('+', begin);
        tokens.insert(label.substr(begin, end == std::string::npos ? std::string::npos
                                                                   : end - begin));
        if (end == std::string::npos)
        {
          break;
        }
        begin = end + 1;
      }
      return tokens;
    };
    auto Model = [](const std::string &token)
    { return token.size() > 1 && token[0] == 'm' ? std::stoi(token.substr(1)) : -1; };
    for (const auto &[segment, list] : intervals.by_segment)
    {
      std::vector<double> ends;
      for (const auto &interval : list)
      {
        ends.push_back(interval.a);
        ends.push_back(interval.b);
      }
      std::sort(ends.begin(), ends.end());
      ends.erase(std::unique(ends.begin(), ends.end()), ends.end());
      auto &pieces = by_segment[segment];
      for (std::size_t i = 0; i + 1 < ends.size(); i++)
      {
        const double midpoint = 0.5 * (ends[i] + ends[i + 1]);
        std::set<std::string> tokens;
        for (const auto &interval : list)
        {
          if (midpoint >= interval.a && midpoint <= interval.b)
          {
            const auto portion_tokens = Tokens(intervals.types[interval.type]);
            tokens.insert(portion_tokens.begin(), portion_tokens.end());
          }
        }
        if (tokens.empty())
        {
          continue;
        }
        std::vector<std::string> ordered(tokens.begin(), tokens.end());
        std::sort(ordered.begin(), ordered.end(),
                  [&](const auto &x, const auto &y) { return Model(x) < Model(y); });
        std::string label;
        for (const auto &token : ordered)
        {
          label += (label.empty() ? "" : "+") + token;
        }
        auto group = std::find(groups.begin(), groups.end(), label);
        if (group == groups.end())
        {
          groups.push_back(label);
          group = std::prev(groups.end());
        }
        pieces.push_back({ends[i], ends[i + 1], static_cast<int>(group - groups.begin())});
      }
    }
  }
};

// The within-R energy density of the quadrature points whose nearest perimeter foot lies
// on a piece of the requested group (`group` = the groups' size: on NO piece, the unclaimed
// within-R remainder). Pieces are half-open [a, b); the last piece of a segment is closed.
class PartitionedEdgeCoefficient : public mfem::Coefficient
{
private:
  mfem::Coefficient &coefficient;
  const EdgeDistanceTree &edge_distance_tree;
  const PartitionedEdgePieces &pieces;
  const double radius;
  const double smoothing;
  const int group;

public:
  PartitionedEdgeCoefficient(mfem::Coefficient &coefficient,
                             const EdgeDistanceTree &edge_distance_tree,
                             const PartitionedEdgePieces &pieces, double radius,
                             double smoothing, int group)
    : coefficient(coefficient), edge_distance_tree(edge_distance_tree), pieces(pieces),
      radius(radius), smoothing(smoothing), group(group)
  {
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    mfem::Vector point(T.GetSpaceDim());
    T.Transform(ip, point);
    const auto nearest = edge_distance_tree.Nearest(point);
    const double within = 1.0 - EdgeDistanceOutsideWeight(
                                    std::sqrt(nearest.distance_squared), radius, smoothing);
    if (within <= 0.0)
    {
      return 0.0;
    }
    int found = static_cast<int>(pieces.groups.size());
    const auto it = pieces.by_segment.find(nearest.segment);
    if (it != pieces.by_segment.end() && !it->second.empty())
    {
      const auto &segment = edge_distance_tree.GetSegment(nearest.segment);
      double length2 = 0.0, projection = 0.0;
      for (int d = 0; d < point.Size(); d++)
      {
        const double chord = segment.p1[d] - segment.p0[d];
        length2 += chord * chord;
        projection += (point(d) - segment.p0[d]) * chord;
      }
      const double foot = std::clamp(projection / length2, 0.0, 1.0) * std::sqrt(length2);
      const auto &list = it->second;
      auto piece = std::upper_bound(list.begin(), list.end(), foot,
                                    [](double x, const auto &p) { return x < p.a; });
      if (piece != list.begin())
      {
        piece = std::prev(piece);
        const bool last = std::next(piece) == list.end();
        if (foot >= piece->a && (foot < piece->b || (last && foot == piece->b)))
        {
          found = piece->group;
        }
      }
    }
    return found == group ? within * coefficient.Eval(T, ip) : 0.0;
  }
};

// OI-MX-LOOPGATE (decision 955): the exact piece partition of one interface's perimeter by
// WEIGHTED portions. The intervals are built with one type per portion index so that each
// interval knows its portion; a piece takes the largest weight among the portions covering
// its midpoint and that portion's label.
struct WeightedEdgePieces
{
  struct Piece
  {
    double a = 0.0, b = 0.0;
    int group = -1;
    double weight = 0.0;
  };
  std::map<std::size_t, std::vector<Piece>> by_segment;  // ascending, non-overlapping
  std::vector<std::string> groups;

  WeightedEdgePieces(
      const EdgeDistanceTree &tree,
      const std::vector<SurfacePostOperator::UncoveredPerimeterPortion> &portions,
      const std::vector<double> &weights, double tolerance)
  {
    MFEM_VERIFY(portions.size() == weights.size(),
                "Weighted perimeter portions require one weight per portion!");
    std::vector<SurfacePostOperator::UncoveredPerimeterPortion> indexed(portions);
    for (std::size_t i = 0; i < indexed.size(); i++)
    {
      indexed[i].type = std::to_string(i);
    }
    const UncoveredEdgeIntervals intervals(tree, indexed, tolerance);
    for (const auto &[segment, list] : intervals.by_segment)
    {
      std::vector<double> ends;
      for (const auto &interval : list)
      {
        ends.push_back(interval.a);
        ends.push_back(interval.b);
      }
      std::sort(ends.begin(), ends.end());
      ends.erase(std::unique(ends.begin(), ends.end()), ends.end());
      auto &pieces = by_segment[segment];
      for (std::size_t i = 0; i + 1 < ends.size(); i++)
      {
        const double midpoint = 0.5 * (ends[i] + ends[i + 1]);
        int best = -1;
        for (const auto &interval : list)
        {
          if (midpoint >= interval.a && midpoint <= interval.b)
          {
            const int portion = std::stoi(intervals.types[interval.type]);
            if (best < 0 || weights[portion] > weights[best])
            {
              best = portion;
            }
          }
        }
        if (best < 0)
        {
          continue;
        }
        auto group = std::find(groups.begin(), groups.end(), portions[best].type);
        if (group == groups.end())
        {
          groups.push_back(portions[best].type);
          group = std::prev(groups.end());
        }
        pieces.push_back({ends[i], ends[i + 1], static_cast<int>(group - groups.begin()),
                          weights[best]});
      }
    }
  }
};

// rho x the within-R energy density of the quadrature points whose nearest perimeter foot
// lies on a weighted piece of the requested group. Pieces are half-open [a, b); the last
// piece of a segment is closed.
class WeightedEdgeCoefficient : public mfem::Coefficient
{
private:
  mfem::Coefficient &coefficient;
  const EdgeDistanceTree &edge_distance_tree;
  const WeightedEdgePieces &pieces;
  const double radius;
  const double smoothing;
  const int group;

public:
  WeightedEdgeCoefficient(mfem::Coefficient &coefficient,
                          const EdgeDistanceTree &edge_distance_tree,
                          const WeightedEdgePieces &pieces, double radius, double smoothing,
                          int group)
    : coefficient(coefficient), edge_distance_tree(edge_distance_tree), pieces(pieces),
      radius(radius), smoothing(smoothing), group(group)
  {
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    mfem::Vector point(T.GetSpaceDim());
    T.Transform(ip, point);
    const auto nearest = edge_distance_tree.Nearest(point);
    const double within = 1.0 - EdgeDistanceOutsideWeight(
                                    std::sqrt(nearest.distance_squared), radius, smoothing);
    if (within <= 0.0)
    {
      return 0.0;
    }
    const auto it = pieces.by_segment.find(nearest.segment);
    if (it == pieces.by_segment.end() || it->second.empty())
    {
      return 0.0;
    }
    const auto &segment = edge_distance_tree.GetSegment(nearest.segment);
    double length2 = 0.0, projection = 0.0;
    for (int d = 0; d < point.Size(); d++)
    {
      const double chord = segment.p1[d] - segment.p0[d];
      length2 += chord * chord;
      projection += (point(d) - segment.p0[d]) * chord;
    }
    const double foot = std::clamp(projection / length2, 0.0, 1.0) * std::sqrt(length2);
    const auto &list = it->second;
    auto piece = std::upper_bound(list.begin(), list.end(), foot,
                                  [](double x, const auto &p) { return x < p.a; });
    if (piece == list.begin())
    {
      return 0.0;
    }
    piece = std::prev(piece);
    const bool last = std::next(piece) == list.end();
    if (piece->group != group ||
        !(foot >= piece->a && (foot < piece->b || (last && foot == piece->b))))
    {
      return 0.0;
    }
    return piece->weight * within * coefficient.Eval(T, ip);
  }
};

}  // namespace

SurfacePostOperator::SurfaceFluxData::SurfaceFluxData(
    const config::SurfaceFluxData &data, const mfem::ParMesh &mesh,
    const mfem::Array<int> &bdr_attr_marker)
{
  // Store boundary attributes for this postprocessing boundary.
  attr_list = SetUpBoundaryProperties(data, bdr_attr_marker);

  // Store the type of flux.
  switch (data.type)
  {
    case SurfaceFlux::ELECTRIC:
      type = SurfaceFlux::ELECTRIC;
      break;
    case SurfaceFlux::MAGNETIC:
      type = SurfaceFlux::MAGNETIC;
      break;
    case SurfaceFlux::POWER:
      type = SurfaceFlux::POWER;
      break;
  }

  // Store information about the global direction for orientation. Note the true boundary
  // normal is used in calculating the flux, this is just used to determine the sign.
  two_sided = data.two_sided;
  if (!two_sided)
  {
    center.SetSize(mesh.SpaceDimension());
    if (data.no_center)
    {
      // Compute the center as the bounding box centroid for all boundary elements making up
      // this postprocessing boundary.
      mfem::Vector bbmin, bbmax;
      mesh::GetAxisAlignedBoundingBox(
          mesh, mesh::AttrToMarker(bdr_attr_marker.Size(), attr_list), true, bbmin, bbmax);
      for (int d = 0; d < mesh.SpaceDimension(); d++)
      {
        center(d) = 0.5 * (bbmin(d) + bbmax(d));
      }
    }
    else
    {
      std::copy_n(data.center.begin(), mesh.SpaceDimension(), center.begin());
    }
  }
}

std::unique_ptr<mfem::Coefficient> SurfacePostOperator::SurfaceFluxData::GetCoefficient(
    const mfem::ParGridFunction *E, const mfem::ParGridFunction *B,
    const MaterialOperator &mat_op, bool imag_permittivity) const
{
  switch (type)
  {
    case SurfaceFlux::ELECTRIC:
      return std::make_unique<
          RestrictedCoefficient<BdrSurfaceFluxCoefficient<SurfaceFlux::ELECTRIC>>>(
          attr_list, E, nullptr, mat_op, two_sided, center, 1.0, imag_permittivity);
    case SurfaceFlux::MAGNETIC:
      return std::make_unique<
          RestrictedCoefficient<BdrSurfaceFluxCoefficient<SurfaceFlux::MAGNETIC>>>(
          attr_list, nullptr, B, mat_op, two_sided, center);
    case SurfaceFlux::POWER:
      return std::make_unique<
          RestrictedCoefficient<BdrSurfaceFluxCoefficient<SurfaceFlux::POWER>>>(
          attr_list, E, B, mat_op, two_sided, center);
  }
  return {};
}

SurfacePostOperator::InterfaceDielectricData::InterfaceDielectricData(
    const config::InterfaceDielectricData &data, const mfem::Array<int> &bdr_attr_marker)
{
  // Store boundary attributes for this postprocessing boundary.
  attr_list = SetUpBoundaryProperties(data, bdr_attr_marker);
  edge_distances = data.edge_distances;
  edge_distance_smoothing = data.edge_distance_smoothing;
  localize_edge_energy = data.localize_edge_energy;
  save_local_edge_energy = data.save_local_edge_energy;
  edge_frame_normal = data.edge_frame_normal.value_or(std::array<double, 3>{});
  ownership_slot = data.ownership_slot;
  if (!data.ownership_data_file.empty())
  {
    ownership = std::make_shared<InterfaceOwnershipPartition>(
        data.ownership_data_file, data.ownership_group, data.ownership_coordinate_scale);
    MFEM_VERIFY(ownership->Slots().count(ownership_slot),
                "Interface ownership group does not contain the selected slot!");
    ownership_rule = std::make_shared<OwnershipQuadrature>(data.ownership_quadrature_order);
  }
  ownership_quadrature_order = data.ownership_quadrature_order;
  if (data.region)
  {
    region = std::make_shared<InterfaceRegion>(data.region->box_min, data.region->box_max,
                                               data.region->segments, data.region->distance,
                                               data.region->normal);
  }
  flux_recovery = data.flux_recovery;

  // Calculate surface dielectric loss according to the formulas from J. Wenner et al.,
  // Surface loss simulations of superconducting coplanar waveguide resonators, Appl. Phys.
  // Lett. (2011). If only a general layer permittivity is specified and not any special
  // metal-air (MA), metal-substrate (MS), or substrate-air (SA) permittivity, compute the
  // numerator of the participation ratio according to the regular formula
  //                       p * E_elec = 1/2 t Re{∫ (ε E)ᴴ E_m dS} .
  switch (data.type)
  {
    case InterfaceDielectric::DEFAULT:
      type = InterfaceDielectric::DEFAULT;
      break;
    case InterfaceDielectric::MA:
      type = InterfaceDielectric::MA;
      break;
    case InterfaceDielectric::MS:
      type = InterfaceDielectric::MS;
      break;
    case InterfaceDielectric::SA:
      type = InterfaceDielectric::SA;
      break;
  }
  t = data.t;
  epsilon = data.epsilon_r;
  tandelta = data.tandelta;
}

std::unique_ptr<mfem::Coefficient>
SurfacePostOperator::InterfaceDielectricData::GetCoefficient(
    const GridFunction &E, const GridFunction *D, const MaterialOperator &mat_op,
    InterfaceDielectricComponent component) const
{
  const GridFunction *recovered_flux = flux_recovery ? D : nullptr;
  MFEM_VERIFY(!flux_recovery || recovered_flux,
              "Interface dielectric flux recovery was requested but no recovered electric "
              "flux was supplied!");
  std::unique_ptr<mfem::Coefficient> coefficient;
  switch (type)
  {
    case InterfaceDielectric::DEFAULT:
      coefficient = std::make_unique<RestrictedCoefficient<
          InterfaceDielectricCoefficient<InterfaceDielectric::DEFAULT>>>(
          attr_list, E, mat_op, t, epsilon, &sum_sides_marker, recovered_flux, component);
      break;
    case InterfaceDielectric::MA:
      coefficient = std::make_unique<
          RestrictedCoefficient<InterfaceDielectricCoefficient<InterfaceDielectric::MA>>>(
          attr_list, E, mat_op, t, epsilon, &sum_sides_marker, recovered_flux, component);
      break;
    case InterfaceDielectric::MS:
      coefficient = std::make_unique<
          RestrictedCoefficient<InterfaceDielectricCoefficient<InterfaceDielectric::MS>>>(
          attr_list, E, mat_op, t, epsilon, &sum_sides_marker, recovered_flux, component);
      break;
    case InterfaceDielectric::SA:
      coefficient = std::make_unique<
          RestrictedCoefficient<InterfaceDielectricCoefficient<InterfaceDielectric::SA>>>(
          attr_list, E, mat_op, t, epsilon, nullptr, recovered_flux, component);
      break;
  }
  if (ownership)
  {
    coefficient = std::make_unique<OwnershipCoefficient>(std::move(coefficient), ownership,
                                                         ownership_slot);
  }
  if (region)
  {
    coefficient = std::make_unique<RegionCoefficient>(std::move(coefficient), region);
  }
  return coefficient;
}

bool SurfacePostOperator::InterfaceDielectricData::Includes(const mfem::Vector &point) const
{
  return (!ownership || ownership->SelectSlot(point) == ownership_slot) &&
         (!region || region->Contains(point));
}

SurfacePostOperator::FarFieldData::FarFieldData(const config::FarFieldPostData &data,
                                                const mfem::ParMesh &mesh,
                                                const mfem::Array<int> &bdr_attr_marker)
  : thetaphis(data.thetaphis)
{
  // Store boundary attributes for this postprocessing boundary.
  attr_list = SetUpBoundaryProperties(data, bdr_attr_marker);
  if (attr_list.Size() == 0)
  {
    return;
  }

  // An interior integration surface (for example around a radiator or scatterer inside a
  // PML) is traced from the side of each face facing the center of its bounding box (see
  // SurfaceFunctional), which only the libCEED path supports.
  const auto marker = mesh::AttrToMarker(bdr_attr_marker.Size(), attr_list);
  for (int i = 0; i < mesh.GetNBE() && !interior; i++)
  {
    if (!marker[mesh.GetBdrAttribute(i) - 1])
    {
      continue;
    }
    int f, o;
    mesh.GetBdrElementFace(i, &f, &o);
    interior = !mesh.GetFaceInformation(f).IsBoundary();
  }
  Mpi::GlobalOr(1, &interior, mesh.GetComm());
}

bool ExclusionEnergyRequested(const config::SolverData &solver)
{
  return (solver.electrostatic.response_correction &&
          solver.electrostatic.response_correction->exclusion_energy) ||
         (solver.surface_response_correction &&
          solver.surface_response_correction->exclusion_energy);
}

SurfacePostOperator::SurfacePostOperator(
    const config::BoundaryPostData &postpro, ProblemType problem_type,
    const MaterialOperator &mat_op, mfem::ParFiniteElementSpace &h1_fespace,
    mfem::ParFiniteElementSpace &nd_fespace, const std::vector<int> &sheet_attributes,
    const config::BoundaryData *boundaries,
    std::shared_ptr<const SurfacePostGeometry> *automatic_geometry,
    bool classified_edge_trees)
  : mat_op(mat_op), h1_fespace(h1_fespace), nd_fespace(nd_fespace),
    axisymmetric(mat_op.GetMesh().IsAxisymmetric()),
    classified_edge_trees(classified_edge_trees)
{
  BlockTimer setup_timer(Timer::CONSTRUCT_SURFACE_POST);

  // Check that boundary attributes have been specified correctly.
  const auto &mesh = *h1_fespace.GetParMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> bdr_attr_marker;
  if (!postpro.flux.empty() || !postpro.dielectric.empty() || !postpro.farfield.empty())
  {
    bdr_attr_marker.SetSize(bdr_attr_max);
    bdr_attr_marker = 0;
    for (auto attr : mesh.bdr_attributes)
    {
      bdr_attr_marker[attr - 1] = 1;
    }
  }

  // Surface flux postprocessing.
  for (const auto &[idx, data] : postpro.flux)
  {
    MFEM_VERIFY(problem_type != ProblemType::ELECTROSTATIC ||
                    data.type == SurfaceFlux::ELECTRIC,
                "Magnetic field or power surface flux postprocessing are not available "
                "for electrostatic problems!");
    MFEM_VERIFY(problem_type != ProblemType::MAGNETOSTATIC ||
                    data.type == SurfaceFlux::MAGNETIC,
                "Electric field or power surface flux postprocessing are not available "
                "for magnetostatic problems!");
    flux_surfs.try_emplace(idx, data, *h1_fespace.GetParMesh(), bdr_attr_marker);
  }

  // Interface dielectric postprocessing.
  MFEM_VERIFY(postpro.dielectric.empty() || problem_type != ProblemType::MAGNETOSTATIC,
              "Interface dielectric loss postprocessing is not available for "
              "magnetostatic problems!");
  using EdgeDistanceTreeKey =
      std::tuple<std::vector<int>, std::vector<int>, std::optional<std::array<double, 3>>,
                 std::vector<std::array<double, 6>>, std::vector<double>>;
  struct EdgeDistanceTreeEntry
  {
    std::shared_ptr<const EdgeDistanceTree> tree;
    std::optional<EdgeExclusionRecord> edge_exclusion;
  };
  std::map<EdgeDistanceTreeKey, EdgeDistanceTreeEntry> edge_distance_trees;
  // The leaf / root flag first (decision 904 A3.2): on a refined nonconforming mesh the
  // retained trees are built from the root-level perimeter and the localized ones from the
  // leaf perimeter, whose segment indices are different spaces.
  using AutomaticEdgeDistanceTreeKey =
      std::tuple<bool, std::vector<std::size_t>, std::optional<std::array<double, 3>>,
                 std::vector<std::array<double, 6>>, std::vector<double>>;
  std::map<AutomaticEdgeDistanceTreeKey, std::shared_ptr<const EdgeDistanceTree>>
      automatic_edge_distance_trees;
  std::map<AutomaticEdgeDistanceTreeKey, std::shared_ptr<const EdgeDistanceTree>>
      classified_edge_distance_trees;
  const SurfacePostGeometry *cached_automatic_geometry =
      automatic_geometry && *automatic_geometry ? automatic_geometry->get() : nullptr;
  std::shared_ptr<SurfacePostGeometry> new_automatic_geometry;
  if (automatic_geometry && !*automatic_geometry)
  {
    new_automatic_geometry = std::make_shared<SurfacePostGeometry>();
  }
  // The perimeter of the RETAINED trees (root-level, decision 904 A3.2) and, where it
  // differs, the LEAF perimeter of the LocalizeEdgeEnergy trees.
  MetalEdgeGeometry metal_edges;
  MetalEdgeGeometry leaf_metal_edges;
  bool leaf_extraction = false;
  const bool use_automatic_edges = std::any_of(
      postpro.dielectric.begin(), postpro.dielectric.end(),
      [cached_automatic_geometry](const auto &entry)
      {
        const auto &[idx, data] = entry;
        if (!data.automatic_edges)
        {
          return false;
        }
        return data.localize_edge_energy || !cached_automatic_geometry ||
               cached_automatic_geometry->automatic_edge_distance_trees.find(idx) ==
                   cached_automatic_geometry->automatic_edge_distance_trees.end();
      });
  if (use_automatic_edges)
  {
    MFEM_VERIFY(boundaries,
                "Automatic metal edge extraction requires complete boundary data!");
    // The joint noise rule at the largest edge distance of the automatic-edge interfaces.
    // The retained (non-localized) trees are, under AMR, the producing run's INITIAL-mesh
    // trees (SurfacePostGeometry); a run on a SAVED nonconforming mesh builds the same
    // trees from the ROOT-LEVEL perimeter (decision 904 A3.2; the identification's rule of
    // decision 899: the refinement-generated vertices merged out, fail closed), so that the
    // uncovered / DomainBoundary / ConductorConsistency portions - root segments - are
    // credited over their whole segments exactly as the producer credits them (on the leaf
    // perimeter a root portion spanning refined sub-segments was credited on its midpoint's
    // sub-segment only). Without refinement-generated vertices the extraction is the leaf
    // perimeter's byte for byte.
    MetalSurfaceExtraction retained_surface = JointNoiseExtractionFor(*boundaries);
    retained_surface.root_level_perimeter = true;
    metal_edges = ExtractMetalEdgeGeometry(mesh, *boundaries, retained_surface);
    // The LocalizeEdgeEnergy trees are rebuilt on every mesh from the LEAF perimeter (the
    // producing run's own form on its last cycle; OI-LOCALEDGE, decision 854): a second
    // extraction only where the two perimeters differ.
    const bool localized = std::any_of(postpro.dielectric.begin(), postpro.dielectric.end(),
                                       [](const auto &entry)
                                       {
                                         const auto &[idx, data] = entry;
                                         (void)idx;
                                         return data.automatic_edges &&
                                                data.localize_edge_energy &&
                                                !data.edge_distances.empty();
                                       });
    if (localized && !metal_edges.refinement_points.empty())
    {
      leaf_metal_edges =
          ExtractMetalEdgeGeometry(mesh, *boundaries, JointNoiseExtractionFor(*boundaries));
      leaf_extraction = true;
    }
  }
  using OwnershipKey = std::tuple<std::string, int, InterfaceDielectric, std::vector<int>,
                                  double, double, bool>;
  struct OwnershipCoverage
  {
    std::set<int> expected, configured;
    int order;
  };
  std::map<OwnershipKey, OwnershipCoverage> ownership_coverage;
  for (const auto &[idx, data] : postpro.dielectric)
  {
    MFEM_VERIFY(data.ownership_data_file.empty() ||
                    problem_type == ProblemType::ELECTROSTATIC,
                "Quadrature ownership is currently implemented for electrostatic "
                "postprocessing only!");
    // FluxRecovery is refused on the cracked internal boundaries (the mesh crack set, as at
    // OURS) and on the INTERIOR faces of the sheet attributes (uncracked by default): the
    // volume L2 projection does not provide a controlled normal trace on a zero-thickness
    // PEC surface. An exterior sheet boundary keeps it, as at OURS.
    if (data.flux_recovery)
    {
      constexpr const char *refusal =
          "Interface dielectric \"FluxRecovery\" is not supported on internal boundary "
          "sheets, cracked or not: the volume L2 projection does not provide a controlled "
          "normal trace on a zero-thickness PEC surface!";
      std::vector<int> entry_sheet_attributes;
      for (const int attr : data.attributes)
      {
        MFEM_VERIFY(!(boundaries && boundaries->cracked_attributes.contains(attr)),
                    refusal);
        if (std::ranges::find(sheet_attributes, attr) != sheet_attributes.end())
        {
          entry_sheet_attributes.push_back(attr);
        }
      }
      int interior_sheet = !entry_sheet_attributes.empty() &&
                           mat_op.GetMesh().GetCrackSides(entry_sheet_attributes).Any();
      Mpi::GlobalMax(1, &interior_sheet, mesh.GetComm());
      MFEM_VERIFY(!interior_sheet, refusal);
    }
    auto it = eps_surfs.try_emplace(idx, data, bdr_attr_marker).first;
    it->second.sum_sides_marker = mesh::AttrToMarker(bdr_attr_max, sheet_attributes, true);
    if (it->second.ownership)
    {
      MFEM_VERIFY(!data.attributes.empty(),
                  "Ownership partition requires boundary attributes!");
      MFEM_VERIFY(
          std::all_of(
              data.attributes.begin(), data.attributes.end(), [&](int a)
              { return a > 0 && a <= bdr_attr_marker.Size() && bdr_attr_marker[a - 1]; }),
          "Ownership partition references an absent boundary attribute!");
      const OwnershipKey key{data.ownership_data_file,
                             data.ownership_group,
                             data.type,
                             data.attributes,
                             data.t,
                             data.epsilon_r,
                             data.flux_recovery};
      auto [coverage, inserted] = ownership_coverage.try_emplace(
          key, OwnershipCoverage{
                   it->second.ownership->Slots(), {}, data.ownership_quadrature_order});
      MFEM_VERIFY(
          coverage->second.order == data.ownership_quadrature_order,
          "All slots of one ownership partition must use the same quadrature order!");
      MFEM_VERIFY(coverage->second.configured.insert(data.ownership_slot).second,
                  "Duplicate slot in an interface ownership partition!");
    }
    if (data.edge_distances.empty())
    {
      continue;
    }

    if (data.automatic_edges)
    {
      if (!data.localize_edge_energy && cached_automatic_geometry)
      {
        const auto cached_tree =
            cached_automatic_geometry->automatic_edge_distance_trees.find(idx);
        if (cached_tree != cached_automatic_geometry->automatic_edge_distance_trees.end())
        {
          it->second.edge_distance_tree = cached_tree->second;
          const auto cached_classified =
              cached_automatic_geometry->classified_edge_distance_trees.find(idx);
          if (cached_classified !=
              cached_automatic_geometry->classified_edge_distance_trees.end())
          {
            it->second.classified_edge_distance_tree = cached_classified->second;
          }
          continue;
        }
      }

      const bool leaf = data.localize_edge_energy && leaf_extraction;
      const MetalEdgeGeometry &edges = leaf ? leaf_metal_edges : metal_edges;
      auto segment_indices = GetInterfaceMetalEdgeSegmentIndices(edges, idx, data.type);
      ExcludeMetalEdgeSegmentIndices(mesh, edges, data.edge_exclude_attributes,
                                     segment_indices);
      EdgeExclusionRecord exclusion_record;
      ExcludeCoincidentMetalEdgeSegmentIndices(
          edges, data.edge_exclude_segments, data.edge_exclude_segment_sagittas,
          data.edge_exclude_segment_tolerance, mesh.SpaceDimension(), segment_indices,
          &exclusion_record);
      if (!data.edge_exclude_segments.empty())
      {
        it->second.edge_exclusion = exclusion_record;
      }
      const AutomaticEdgeDistanceTreeKey tree_key{
          leaf, segment_indices, data.edge_frame_normal, data.edge_exclude_segments,
          data.edge_exclude_segment_sagittas};
      auto tree_it = automatic_edge_distance_trees.find(tree_key);
      if (tree_it == automatic_edge_distance_trees.end())
      {
        std::vector<std::array<double, 3>> process_normals;
        if (data.localize_edge_energy)
        {
          process_normals = BuildMetalEdgeProcessNormals(
              mesh, edges, segment_indices, [this](int attribute)
              { return this->mat_op.GetLightSpeedMax(attribute); }, data.edge_frame_normal);
        }
        tree_it =
            automatic_edge_distance_trees
                .try_emplace(tree_key, BuildEdgeDistanceTree(edges, segment_indices,
                                                             std::move(process_normals)))
                .first;
      }
      it->second.edge_distance_tree = tree_it->second;
      if (!data.localize_edge_energy && new_automatic_geometry)
      {
        new_automatic_geometry->automatic_edge_distance_trees.emplace(idx, tree_it->second);
      }
      // The exclusion-energy export's CLASSIFIED tree (decision 1076): every segment of the
      // interface, physical or not, under the same exclusion rules; cached alike. Built
      // only when the export is requested (d1073 MINOR-5): every other table reads the
      // physical tree above, so an unrequested run holds no second tree.
      if (!classified_edge_trees)
      {
        continue;
      }
      auto classified_indices =
          GetInterfaceMetalEdgeSegmentIndices(edges, idx, data.type, false);
      ExcludeMetalEdgeSegmentIndices(mesh, edges, data.edge_exclude_attributes,
                                     classified_indices);
      ExcludeCoincidentMetalEdgeSegmentIndices(
          edges, data.edge_exclude_segments, data.edge_exclude_segment_sagittas,
          data.edge_exclude_segment_tolerance, mesh.SpaceDimension(), classified_indices);
      const AutomaticEdgeDistanceTreeKey classified_key{
          leaf, classified_indices, data.edge_frame_normal, data.edge_exclude_segments,
          data.edge_exclude_segment_sagittas};
      auto classified_it = classified_edge_distance_trees.find(classified_key);
      if (classified_it == classified_edge_distance_trees.end())
      {
        // No process normals (a truncation / fold segment has none to infer): the
        // classified tree reads the plain distance, under LocalizeEdgeEnergy too.
        classified_it =
            classified_edge_distance_trees
                .try_emplace(classified_key,
                             BuildEdgeDistanceTree(edges, classified_indices, {}, false))
                .first;
      }
      it->second.classified_edge_distance_tree = classified_it->second;
      if (!data.localize_edge_energy && new_automatic_geometry)
      {
        new_automatic_geometry->classified_edge_distance_trees.emplace(
            idx, classified_it->second);
      }
    }
    else
    {
      const EdgeDistanceTreeKey tree_key{data.edge_attributes, data.edge_exclude_attributes,
                                         data.edge_frame_normal, data.edge_exclude_segments,
                                         data.edge_exclude_segment_sagittas};
      auto tree_it = edge_distance_trees.find(tree_key);
      if (tree_it == edge_distance_trees.end())
      {
        EdgeDistanceTreeEntry entry;
        EdgeExclusionRecord exclusion_record;
        entry.tree =
            BuildEdgeDistanceTree(mesh, data.edge_attributes, data.edge_exclude_attributes,
                                  data.edge_frame_normal, data.edge_exclude_segments,
                                  data.edge_exclude_segment_tolerance,
                                  data.edge_exclude_segment_sagittas, &exclusion_record);
        if (!data.edge_exclude_segments.empty())
        {
          entry.edge_exclusion = exclusion_record;
        }
        tree_it = edge_distance_trees.try_emplace(tree_key, std::move(entry)).first;
      }
      it->second.edge_distance_tree = tree_it->second.tree;
      it->second.edge_exclusion = tree_it->second.edge_exclusion;
    }
  }
  for (const auto &[key, coverage] : ownership_coverage)
  {
    MFEM_VERIFY(
        coverage.configured == coverage.expected,
        "Incomplete interface ownership partition: every file slot must be configured!");
  }
  if (automatic_geometry && !*automatic_geometry && new_automatic_geometry &&
      !new_automatic_geometry->automatic_edge_distance_trees.empty())
  {
    *automatic_geometry = std::move(new_automatic_geometry);
  }

  // FarField postprocessing.
  MFEM_VERIFY(postpro.farfield.empty() || problem_type == ProblemType::DRIVEN ||
                  problem_type == ProblemType::EIGENMODE,
              "Far-field extraction is only available for driven and eigenmode problems!");

  // Check the materials on the far-field surface. The Stratton-Chu integrand assumes an
  // isotropic medium, and on an interior surface it uses the medium of the traced side for
  // both sides, which must then be the same (same speed of light). Both sides of each
  // boundary element are checked, including face neighbor elements on parallel interfaces,
  // and the result is combined across processes so that all of them fail consistently.
  if (!postpro.farfield.empty())
  {
    const auto &mesh = *nd_fespace.GetParMesh();
    const int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
    const int attr_max = mesh.attributes.Size() ? mesh.attributes.Max() : 0;
    mfem::Array<int> farfield_marker =
        mesh::AttrToMarker(bdr_attr_max, postpro.farfield.attributes);
    // A process only has the materials of its elements and face neighbor elements: check
    // them there, and reduce the result.
    mfem::Array<int> anisotropic(attr_max);
    anisotropic = 0;
    int mismatch[2] = {attr_max + 1, attr_max + 1};
    mfem::FaceElementTransformations FET;
    mfem::IsoparametricTransformation T1, T2;
    for (int i = 0; i < mesh.GetNBE(); i++)
    {
      if (!farfield_marker[mesh.GetBdrAttribute(i) - 1])
      {
        continue;
      }
      BdrGridFunctionCoefficient::GetBdrElementNeighborTransformations(i, mesh, FET, T1,
                                                                       T2);
      const int attr1 = FET.Elem1->Attribute;
      anisotropic[attr1 - 1] = !mat_op.IsIsotropic(attr1);
      if (FET.Elem2)
      {
        const int attr2 = FET.Elem2->Attribute;
        anisotropic[attr2 - 1] = !mat_op.IsIsotropic(attr2);
        const double c1 = mat_op.GetLightSpeedMax(attr1),
                     c2 = mat_op.GetLightSpeedMax(attr2);
        if (std::abs(c1 - c2) > 1.0e-12 * std::max(c1, c2))
        {
          // Keep the pair with the smallest attributes, for a deterministic message.
          const auto [a, b] = std::minmax(attr1, attr2);
          if (a < mismatch[0] || (a == mismatch[0] && b < mismatch[1]))
          {
            mismatch[0] = a;
            mismatch[1] = b;
          }
        }
      }
    }
    Mpi::GlobalMax(attr_max, anisotropic.HostReadWrite(), mesh.GetComm());
    {
      // Global pair with the smallest attributes: the smallest first attribute, then the
      // smallest second attribute among the processes which have it.
      const int a_loc = mismatch[0];
      Mpi::GlobalMin(1, &mismatch[0], mesh.GetComm());
      if (a_loc != mismatch[0])
      {
        mismatch[1] = attr_max + 1;
      }
      Mpi::GlobalMin(1, &mismatch[1], mesh.GetComm());
    }
    for (int attr = 1; attr <= attr_max; attr++)
    {
      MFEM_VERIFY(!anisotropic[attr - 1],
                  "FarField requires isotropic materials, but attribute " +
                      std::to_string(attr) + " is not.");
    }
    MFEM_VERIFY(mismatch[0] > attr_max,
                "FarField requires the same material on both sides of an interior surface, "
                "but attributes " +
                    std::to_string(mismatch[0]) + " and " + std::to_string(mismatch[1]) +
                    " have a different speed of light.");
  }

  farfield = FarFieldData(postpro.farfield, *nd_fespace.GetParMesh(), bdr_attr_marker);
}

SurfacePostOperator::SurfacePostOperator(
    const IoData &iodata, const MaterialOperator &mat_op,
    mfem::ParFiniteElementSpace &h1_fespace, mfem::ParFiniteElementSpace &nd_fespace,
    std::shared_ptr<const SurfacePostGeometry> *automatic_geometry)
  : SurfacePostOperator(iodata.boundaries.postpro, iodata.problem.type, mat_op, h1_fespace,
                        nd_fespace, iodata.boundaries.GetSheetAttributes(),
                        &iodata.boundaries, automatic_geometry,
                        ExclusionEnergyRequested(iodata.solver))
{
}

std::complex<double> SurfacePostOperator::GetSurfaceFlux(int idx, const GridFunction *E,
                                                         const GridFunction *B) const
{
  // For complex-valued fields, output the separate real and imaginary parts for the time-
  // harmonic quantity. For power flux (Poynting vector), output only the stationary real
  // part and not the part which has double the frequency. For the electric flux of a
  // complex field in a lossy dielectric, D = ε E with the complex permittivity
  // ε = Re{ε} + i Im{ε} (Im{ε} = -Re{ε} tan δ): the Im{ε}-weighted flux F'' of the same
  // field adds as Φ = F'(E) + i F''(E), i.e. Φ_r = F'(E_r) - F''(E_i), Φ_i = F'(E_i) +
  // F''(E_r).
  auto it = flux_surfs.find(idx);
  MFEM_VERIFY(it != flux_surfs.end(),
              "Unknown surface flux postprocessing index requested!");
  const bool has_imag = (E) ? E->HasImag() : B->HasImag();
  const bool complex_permittivity =
      (it->second.type == SurfaceFlux::ELECTRIC && has_imag && mat_op.HasLossTangent());
  const auto &mesh = *h1_fespace.GetParMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> attr_marker = mesh::AttrToMarker(bdr_attr_max, it->second.attr_list);

  // Use the libCEED surface functional path when supported, avoiding per-call
  // boundary LinearForm assembly in the legacy coefficient path.
  auto &func = flux_funcs[idx];
  if (!func)
  {
    func = std::make_unique<SurfaceFunctional>(
        mat_op.GetMesh(), attr_marker, E ? E->ParFESpace() : nullptr,
        B ? B->ParFESpace() : nullptr, mat_op, it->second.type, it->second.two_sided,
        it->second.center);
  }
  if (func && func->IsValid())
  {
    std::complex<double> flux = func->EvalFlux(E, B);
    if (complex_permittivity)
    {
      auto &func_imag = flux_funcs_imag[idx];
      if (!func_imag)
      {
        func_imag = std::make_unique<SurfaceFunctional>(
            mat_op.GetMesh(), attr_marker, E->ParFESpace(), nullptr, mat_op,
            it->second.type, it->second.two_sided, it->second.center,
            /*imag_permittivity*/ true);
      }
      MFEM_VERIFY(func_imag && func_imag->IsValid(),
                  "libCEED postprocessing was expected for the lossy electric surface "
                  "flux, but SurfaceFunctional could not assemble!");
      flux += 1i * func_imag->EvalFlux(E, B);
    }
    return flux;
  }
  if (IsSupportedSurfaceFluxDimension(mesh))
  {
    RequireCeedSurfaceFunctional(func.get(), "3D surface flux");
  }

  auto Integrate = [&](const mfem::ParGridFunction *Ei, const mfem::ParGridFunction *Bi,
                       bool imag_permittivity)
  {
    auto f = it->second.GetCoefficient(Ei, Bi, mat_op, imag_permittivity);
    return GetLocalSurfaceIntegral(*f, attr_marker);
  };
  std::complex<double> dot(
      Integrate(E ? &E->Real() : nullptr, B ? &B->Real() : nullptr, false), 0.0);
  if (has_imag)
  {
    double doti = Integrate(E ? &E->Imag() : nullptr, B ? &B->Imag() : nullptr, false);
    if (it->second.type == SurfaceFlux::POWER)
    {
      dot += doti;
    }
    else
    {
      dot.imag(doti);
    }
    if (complex_permittivity)
    {
      dot += 1i * std::complex<double>(Integrate(&E->Real(), nullptr, true),
                                       Integrate(&E->Imag(), nullptr, true));
    }
  }
  Mpi::GlobalSum(1, &dot, (E) ? E->GetComm() : B->GetComm());
  return dot;
}

double SurfacePostOperator::GetInterfaceLossTangent(int idx) const
{
  auto it = eps_surfs.find(idx);
  MFEM_VERIFY(it != eps_surfs.end(),
              "Unknown interface dielectric postprocessing index requested!");
  return it->second.tandelta;
}

double SurfacePostOperator::GetInterfaceElectricFieldEnergy(int idx, const GridFunction &E,
                                                            const GridFunction *D) const
{
  auto it = eps_surfs.find(idx);
  MFEM_VERIFY(it != eps_surfs.end(),
              "Unknown interface dielectric postprocessing index requested!");
  const auto &mesh = *h1_fespace.GetParMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> attr_marker = mesh::AttrToMarker(bdr_attr_max, it->second.attr_list);

  // Use the libCEED surface functional path when supported, avoiding per-call
  // boundary LinearForm assembly in the legacy coefficient path.
  // It evaluates the plain interface energy density with the same two-sided sheet rule
  // (sum_sides_marker) as the coefficient; an entry with recovered flux, an ownership
  // slot, a region, its own ownership quadrature order, or on an axisymmetric mesh (the
  // revolution measure) keeps the coefficient path, which the edge / partitioned /
  // response energies of the entry share. So does an entry with edge distances (the
  // entries the response correction and its coverage accounting act on): its within-R
  // energy W (this energy minus the edge energy outside R) and the coverage partition's
  // within-R total must come from the same quadrature, else their closure residual is
  // quadrature error instead of round-off (P4 FM3 merge interaction; ruling 2036 (2)).
  const bool coefficient_path = it->second.flux_recovery || it->second.ownership ||
                                it->second.region ||
                                it->second.ownership_quadrature_order > 0 ||
                                !it->second.edge_distances.empty() || axisymmetric;
  if (!coefficient_path)
  {
    auto &func = eps_funcs[idx];
    if (!func)
    {
      func = std::make_unique<SurfaceFunctional>(
          mat_op.GetMesh(), attr_marker, *E.ParFESpace(), mat_op, it->second.type,
          it->second.t, it->second.epsilon, &it->second.sum_sides_marker);
    }
    if (func && func->IsValid())
    {
      return func->Eval(E);
    }
    if (IsSupportedSurfaceFunctionalDimension(mesh))
    {
      RequireCeedSurfaceFunctional(func.get(), "interface dielectric postprocessing");
    }
  }

  auto f = it->second.GetCoefficient(E, D, mat_op);
  double dot =
      GetLocalSurfaceIntegral(*f, attr_marker, it->second.ownership_quadrature_order);
  Mpi::GlobalSum(1, &dot, E.GetComm());
  return dot;
}

std::vector<SurfacePostOperator::InterfaceEdgeEnergy>
SurfacePostOperator::GetInterfaceEdgeElectricFieldEnergies(int idx, const GridFunction &E,
                                                           const GridFunction *D) const
{
  auto it = eps_surfs.find(idx);
  MFEM_VERIFY(it != eps_surfs.end(),
              "Unknown interface dielectric postprocessing index requested!");
  const auto &data = it->second;
  if (data.edge_distances.empty())
  {
    return {};
  }

  const auto &mesh = *h1_fespace.GetParMesh();
  const int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  const auto attr_marker = mesh::AttrToMarker(bdr_attr_max, data.attr_list);
  auto coefficient = data.GetCoefficient(E, D, mat_op);

  std::vector<InterfaceEdgeEnergy> energies;
  energies.reserve(data.edge_distances.size());
  std::vector<double> local_energy(2 * data.edge_distances.size());
  for (std::size_t i = 0; i < data.edge_distances.size(); i++)
  {
    const double distance = data.edge_distances[i];
    EdgeDistanceCoefficient outside(*coefficient, *data.edge_distance_tree, distance,
                                    mfem::infinity(), data.edge_distance_smoothing);
    EdgeDistanceCoefficient annulus(*coefficient, *data.edge_distance_tree, distance,
                                    2.0 * distance, data.edge_distance_smoothing);
    local_energy[2 * i] =
        GetLocalSurfaceIntegral(outside, attr_marker, data.ownership_quadrature_order);
    local_energy[2 * i + 1] =
        GetLocalSurfaceIntegral(annulus, attr_marker, data.ownership_quadrature_order);
  }
  Mpi::GlobalSum(static_cast<int>(local_energy.size()), local_energy.data(), E.GetComm());
  for (std::size_t i = 0; i < data.edge_distances.size(); i++)
  {
    energies.push_back(
        {data.edge_distances[i], local_energy[2 * i], local_energy[2 * i + 1]});
  }
  return energies;
}

std::map<int, SurfacePostOperator::UncoveredEdgeEnergy>
SurfacePostOperator::GetInterfaceUncoveredEdgeEnergies(
    const std::set<int> &indices, const GridFunction &E, const GridFunction *D,
    const std::vector<UncoveredPerimeterPortion> &portions) const
{
  const auto &mesh = *h1_fespace.GetParMesh();
  const int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  std::map<int, UncoveredEdgeEnergy> energies;
  std::vector<double> local_energies;
  std::vector<std::pair<int, std::string>> entries;  // (interface, type) per local energy
  for (const int idx : indices)
  {
    auto it = eps_surfs.find(idx);
    MFEM_VERIFY(it != eps_surfs.end(),
                "Unknown interface dielectric postprocessing index requested!");
    const auto &data = it->second;
    MFEM_VERIFY(!data.edge_distances.empty() && data.edge_distance_tree,
                "Uncovered-requirement energies require a matching edge distance!");
    energies.emplace(idx, UncoveredEdgeEnergy{});
    if (portions.empty())
    {
      continue;
    }
    const double radius = data.edge_distances.back();
    const UncoveredEdgeIntervals intervals(*data.edge_distance_tree, portions,
                                           kSignatureParameterToleranceOverRadius * radius);
    if (intervals.types.empty())
    {
      continue;
    }
    const auto attr_marker = mesh::AttrToMarker(bdr_attr_max, data.attr_list);
    auto coefficient = data.GetCoefficient(E, D, mat_op);
    for (std::size_t type = 0; type < intervals.types.size(); type++)
    {
      UncoveredEdgeCoefficient uncovered(*coefficient, *data.edge_distance_tree, intervals,
                                         radius, data.edge_distance_smoothing,
                                         static_cast<int>(type));
      local_energies.push_back(
          GetLocalSurfaceIntegral(uncovered, attr_marker, data.ownership_quadrature_order));
      entries.emplace_back(idx, intervals.types[type]);
    }
  }
  if (!local_energies.empty())
  {
    Mpi::GlobalSum(static_cast<int>(local_energies.size()), local_energies.data(),
                   E.GetComm());
  }
  for (std::size_t i = 0; i < entries.size(); i++)
  {
    auto &energy = energies.at(entries[i].first);
    energy.by_type[entries[i].second] += local_energies[i];
    energy.energy += local_energies[i];
  }
  return energies;
}

std::map<int, std::map<std::string, double>>
SurfacePostOperator::GetInterfacePartitionedEdgeEnergies(
    const std::set<int> &indices, const GridFunction &E, const GridFunction *D,
    PerimeterPortionView portions) const
{
  // PatchEnergy (decision 922 B-1): the physical tree, the remainder "unclaimed".
  return PartitionInterfaceEdgeEnergies(indices, E, D, portions, false, "unclaimed");
}

std::map<int, std::map<std::string, double>>
SurfacePostOperator::GetInterfaceClassifiedEdgeEnergies(
    const std::set<int> &indices, const GridFunction &E, const GridFunction *D,
    const std::vector<UncoveredPerimeterPortion> &portions) const
{
  // The exclusion-energy export (decision 1076): the classified tree, the remainder
  // "Claimed".
  return PartitionInterfaceEdgeEnergies(indices, E, D, portions, true, "Claimed");
}

std::map<int, std::map<std::string, double>>
SurfacePostOperator::PartitionInterfaceEdgeEnergies(
    const std::set<int> &indices, const GridFunction &E, const GridFunction *D,
    PerimeterPortionView portions, bool classified,
    const std::string &remainder_label) const
{
  // One integral per (interface, group) and one for the remainder; the groups come from
  // the exact piece partition of the portions on the selected tree.
  const auto &mesh = *h1_fespace.GetParMesh();
  const int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  std::map<int, std::map<std::string, double>> energies;
  std::vector<double> local_energies;
  std::vector<std::pair<int, std::string>> entries;
  for (const int idx : indices)
  {
    auto it = eps_surfs.find(idx);
    MFEM_VERIFY(it != eps_surfs.end(),
                "Unknown interface dielectric postprocessing index requested!");
    const auto &data = it->second;
    MFEM_VERIFY(!data.edge_distances.empty() && data.edge_distance_tree,
                "Partitioned edge energies require a matching edge distance!");
    MFEM_VERIFY(!classified || data.classified_edge_distance_tree,
                "The exclusion-energy export requires the classified edge tree of an "
                "automatic-edge interface (interface "
                    << idx << " has none)!");
    const EdgeDistanceTree &tree =
        classified ? *data.classified_edge_distance_tree : *data.edge_distance_tree;
    energies.emplace(idx, std::map<std::string, double>{});
    const double radius = data.edge_distances.back();
    const UncoveredEdgeIntervals intervals(tree, portions,
                                           kSignatureParameterToleranceOverRadius * radius);
    const PartitionedEdgePieces pieces(intervals);
    const auto attr_marker = mesh::AttrToMarker(bdr_attr_max, data.attr_list);
    auto coefficient = data.GetCoefficient(E, D, mat_op);
    for (std::size_t group = 0; group <= pieces.groups.size(); group++)
    {
      PartitionedEdgeCoefficient partitioned(*coefficient, tree, pieces, radius,
                                             data.edge_distance_smoothing,
                                             static_cast<int>(group));
      local_energies.push_back(GetLocalSurfaceIntegral(partitioned, attr_marker,
                                                       data.ownership_quadrature_order));
      entries.emplace_back(idx, group < pieces.groups.size() ? pieces.groups[group]
                                                             : remainder_label);
    }
  }
  if (!local_energies.empty())
  {
    Mpi::GlobalSum(static_cast<int>(local_energies.size()), local_energies.data(),
                   E.GetComm());
  }
  for (std::size_t i = 0; i < entries.size(); i++)
  {
    energies.at(entries[i].first)[entries[i].second] += local_energies[i];
  }
  return energies;
}

std::map<int, std::map<std::pair<std::string, std::string>, double>>
SurfacePostOperator::GetInterfaceUncoveredFeatureEdgeEnergies(
    const std::set<int> &indices, const GridFunction &E, const GridFunction *D,
    const std::vector<UncoveredPerimeterPortion> &portions) const
{
  // Per type (in the order of first appearance) the type's portions relabelled by feature:
  // the token "m<feature>" makes PartitionedEdgePieces order a union label by feature id.
  std::vector<std::string> types;
  std::vector<std::vector<UncoveredPerimeterPortion>> portions_by_type;
  for (const auto &portion : portions)
  {
    auto type = std::find(types.begin(), types.end(), portion.type);
    if (type == types.end())
    {
      types.push_back(portion.type);
      portions_by_type.emplace_back();
      type = std::prev(types.end());
    }
    auto relabelled = portion;
    relabelled.type = "m" + std::to_string(portion.feature);
    portions_by_type[type - types.begin()].push_back(std::move(relabelled));
  }
  auto FeatureLabel = [](const std::string &label)
  {
    // "m12+m40" -> "12+40".
    std::string features;
    for (std::size_t i = 0; i < label.size(); i++)
    {
      if (label[i] == 'm' && (i == 0 || label[i - 1] == '+'))
      {
        continue;
      }
      features += label[i];
    }
    return features;
  };
  auto &mesh = *h1_fespace.GetParMesh();
  const int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  std::map<int, std::map<std::pair<std::string, std::string>, double>> energies;
  std::vector<double> local_energies;
  std::vector<std::tuple<int, std::size_t, std::string>>
      entries;  // (interface, type, label)
  for (const int idx : indices)
  {
    auto it = eps_surfs.find(idx);
    MFEM_VERIFY(it != eps_surfs.end(),
                "Unknown interface dielectric postprocessing index requested!");
    const auto &data = it->second;
    MFEM_VERIFY(!data.edge_distances.empty() && data.edge_distance_tree,
                "Uncovered feature energies require a matching edge distance!");
    energies.emplace(idx, std::map<std::pair<std::string, std::string>, double>{});
    const EdgeDistanceTree &tree = *data.edge_distance_tree;
    const double radius = data.edge_distances.back();
    std::vector<PartitionedEdgePieces> pieces;
    std::vector<std::size_t> offsets;  // the type's first entry in local_energies
    for (std::size_t type = 0; type < types.size(); type++)
    {
      const UncoveredEdgeIntervals intervals(
          tree, portions_by_type[type], kSignatureParameterToleranceOverRadius * radius);
      pieces.emplace_back(intervals);
      offsets.push_back(local_energies.size());
      for (const auto &label : pieces.back().groups)
      {
        local_energies.push_back(0.0);
        entries.emplace_back(idx, type, FeatureLabel(label));
      }
    }
    if (entries.empty() || std::get<0>(entries.back()) != idx)
    {
      continue;  // no portion on this interface's perimeter
    }
    const auto attr_marker = mesh::AttrToMarker(bdr_attr_max, data.attr_list);
    auto coefficient = data.GetCoefficient(E, D, mat_op);
    const int quadrature_order = data.ownership_quadrature_order;
    const std::optional<OwnershipQuadrature> ownership_rules =
        quadrature_order > 0 ? std::optional<OwnershipQuadrature>(quadrature_order)
                             : std::nullopt;
    mfem::Vector point(mesh.SpaceDimension());
    for (int be = 0; be < mesh.GetNBE(); be++)
    {
      const int attr = mesh.GetBdrAttribute(be);
      if (attr <= 0 || attr > attr_marker.Size() || !attr_marker[attr - 1])
      {
        continue;
      }
      auto *T = mesh.GetBdrElementTransformation(be);
      const auto geometry = h1_fespace.GetBE(be)->GetGeomType();
      const mfem::IntegrationRule &rule =
          ownership_rules
              ? ownership_rules->Get(geometry)
              : mfem::IntRules.Get(geometry, fem::DefaultIntegrationOrder::Get(*T));
      for (int q = 0; q < rule.GetNPoints(); q++)
      {
        const auto &ip = rule.IntPoint(q);
        T->SetIntPoint(&ip);
        T->Transform(ip, point);
        const auto nearest = tree.Nearest(point);
        const double within =
            1.0 - EdgeDistanceOutsideWeight(std::sqrt(nearest.distance_squared), radius,
                                            data.edge_distance_smoothing);
        if (within <= 0.0)
        {
          continue;
        }
        // PartitionedEdgeCoefficient's foot test, per type.
        const auto &segment = tree.GetSegment(nearest.segment);
        double length2 = 0.0, projection = 0.0;
        for (int d = 0; d < point.Size(); d++)
        {
          const double chord = segment.p1[d] - segment.p0[d];
          length2 += chord * chord;
          projection += (point(d) - segment.p0[d]) * chord;
        }
        const double foot = std::clamp(projection / length2, 0.0, 1.0) * std::sqrt(length2);
        std::optional<double> value;
        for (std::size_t type = 0; type < types.size(); type++)
        {
          const auto list = pieces[type].by_segment.find(nearest.segment);
          if (list == pieces[type].by_segment.end() || list->second.empty())
          {
            continue;
          }
          auto piece = std::upper_bound(list->second.begin(), list->second.end(), foot,
                                        [](double x, const auto &p) { return x < p.a; });
          if (piece == list->second.begin())
          {
            continue;
          }
          piece = std::prev(piece);
          const bool last = std::next(piece) == list->second.end();
          if (!(foot >= piece->a && (foot < piece->b || (last && foot == piece->b))))
          {
            continue;
          }
          if (!value)
          {
            value = PhysicalQuadratureWeight(*T, ip, axisymmetric) *
                    (within * coefficient->Eval(*T, ip));
          }
          local_energies[offsets[type] + piece->group] += *value;
        }
      }
    }
  }
  if (!local_energies.empty())
  {
    Mpi::GlobalSum(static_cast<int>(local_energies.size()), local_energies.data(),
                   E.GetComm());
  }
  for (std::size_t i = 0; i < entries.size(); i++)
  {
    const auto &[idx, type, label] = entries[i];
    energies.at(idx)[{types[type], label}] += local_energies[i];
  }
  return energies;
}

std::pair<SurfacePostOperator::ClaimEndMatch, double>
SurfacePostOperator::ClassifyClaimEnd(const std::array<double, 3> &p,
                                      const std::array<double, 3> &v, double radius)
{
  bool exact = true, roundoff = true;
  double distance2 = 0.0;
  for (int d = 0; d < 3; d++)
  {
    const double difference = p[d] - v[d];
    exact = exact && difference == 0.0;
    roundoff =
        roundoff && std::abs(difference) <= kClaimEndRoundoffUlps *
                                                std::numeric_limits<double>::epsilon() *
                                                std::max(std::abs(v[d]), radius);
    distance2 += difference * difference;
  }
  const double distance = std::sqrt(distance2);
  const ClaimEndMatch match = exact      ? ClaimEndMatch::EXACT
                              : roundoff ? ClaimEndMatch::ROUNDOFF
                              : distance <= kSupportEndCoincidenceOverRadius * radius
                                  ? ClaimEndMatch::NEAR
                                  : ClaimEndMatch::BEYOND;
  return {match, distance};
}

std::map<int, SurfacePostOperator::InterfacePerimeterLengths>
SurfacePostOperator::GetInterfacePerimeterLengths(
    const std::set<int> &indices, const std::vector<PerimeterPortionView> &lists) const
{
  std::map<int, InterfacePerimeterLengths> lengths;
  for (const int idx : indices)
  {
    auto it = eps_surfs.find(idx);
    MFEM_VERIFY(it != eps_surfs.end(),
                "Unknown interface dielectric postprocessing index requested!");
    const auto &data = it->second;
    auto &entry = lengths[idx];
    entry.lists.assign(lists.size(), 0.0);
    if (data.edge_distances.empty() || !data.edge_distance_tree)
    {
      continue;
    }
    const EdgeDistanceTree &tree = *data.edge_distance_tree;
    for (std::size_t i = 0; i < tree.Size(); i++)
    {
      const auto &segment = tree.GetSegment(i);
      double length2 = 0.0;
      for (int d = 0; d < 3; d++)
      {
        length2 += (segment.p1[d] - segment.p0[d]) * (segment.p1[d] - segment.p0[d]);
      }
      entry.total += std::sqrt(length2);
    }
    const double radius = data.edge_distances.back();
    for (std::size_t k = 0; k < lists.size(); k++)
    {
      if (lists[k].Empty())
      {
        continue;
      }
      const UncoveredEdgeIntervals intervals(
          tree, lists[k], kSignatureParameterToleranceOverRadius * radius);
      for (const auto &[segment, pieces] : intervals.by_segment)
      {
        std::vector<std::pair<double, double>> sorted;
        for (const auto &piece : pieces)
        {
          sorted.emplace_back(piece.a, piece.b);
        }
        std::sort(sorted.begin(), sorted.end());
        double a = sorted.front().first, b = sorted.front().second;
        for (const auto &[a1, b1] : sorted)
        {
          if (a1 > b)
          {
            entry.lists[k] += b - a;
            a = a1;
            b = b1;
          }
          else
          {
            b = std::max(b, b1);
          }
        }
        entry.lists[k] += b - a;
      }
    }
  }
  return lengths;
}

std::map<int, SurfacePostOperator::InterfaceCoveragePartition>
SurfacePostOperator::GetInterfaceCoveragePartition(
    const std::set<int> &indices, const GridFunction &E, const GridFunction *D,
    const std::vector<CoveragePartitionClass> &kept_classes, PerimeterPortionView served,
    const std::vector<CoveragePartitionClass> &remainder_classes) const
{
  // Per interface the local bins in a fixed order (identical on every rank: the portion
  // lists are): per kept class its piece groups, the served groups, the remainder and the
  // within-R total; ONE reduction for all of them.
  auto &mesh = *h1_fespace.GetParMesh();
  const int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  struct InterfaceBins
  {
    int idx = 0;
    std::vector<PartitionedEdgePieces> pieces;  // the kept classes, then the served claims
    std::vector<std::size_t> offsets;           // each one's first bin
    std::size_t unclaimed = 0, within = 0;
    // Rulings 1714 / 1715: the served claims' and the kept classes' end census; `stake`
    // = the first of the 4 bins of the energy at those segment ends (served round-off /
    // near, kept round-off / near).
    ClaimEndCheck served_check, kept_check;
    std::size_t stake = 0;
  };
  std::vector<InterfaceBins> bins;
  std::vector<double> local_energies;
  for (const int idx : indices)
  {
    auto it = eps_surfs.find(idx);
    MFEM_VERIFY(it != eps_surfs.end(),
                "Unknown interface dielectric postprocessing index requested!");
    const auto &data = it->second;
    MFEM_VERIFY(!data.edge_distances.empty() && data.edge_distance_tree,
                "The coverage partition requires a matching edge distance!");
    const EdgeDistanceTree &tree = *data.edge_distance_tree;
    const double radius = data.edge_distances.back();
    const double tolerance = kSignatureParameterToleranceOverRadius * radius;
    auto &interface = bins.emplace_back();
    interface.idx = idx;
    interface.served_check.radius = interface.kept_check.radius = radius;
    interface.served_check.snap = true;  // ruling 1714 (B): the served claims only
    auto AddPieces = [&](PerimeterPortionView portions, ClaimEndCheck *check)
    {
      interface.pieces.emplace_back(
          UncoveredEdgeIntervals(tree, portions, tolerance, check));
      interface.offsets.push_back(local_energies.size());
      local_energies.resize(local_energies.size() + interface.pieces.back().groups.size(),
                            0.0);
    };
    for (const auto &kept : kept_classes)
    {
      MFEM_VERIFY(kept.portions, "A coverage-partition class requires its portions!");
      AddPieces(*kept.portions, &interface.kept_check);
    }
    AddPieces(served, &interface.served_check);
    for (const auto &remainder : remainder_classes)
    {
      MFEM_VERIFY(remainder.portions, "A coverage-partition class requires its portions!");
      AddPieces(*remainder.portions, nullptr);
    }
    interface.unclaimed = local_energies.size();
    local_energies.push_back(0.0);
    interface.within = local_energies.size();
    local_energies.push_back(0.0);
    interface.stake = local_energies.size();
    local_energies.resize(local_energies.size() + 4, 0.0);

    const auto attr_marker = mesh::AttrToMarker(bdr_attr_max, data.attr_list);
    auto coefficient = data.GetCoefficient(E, D, mat_op);
    const int quadrature_order = data.ownership_quadrature_order;
    const std::optional<OwnershipQuadrature> ownership_rules =
        quadrature_order > 0 ? std::optional<OwnershipQuadrature>(quadrature_order)
                             : std::nullopt;
    mfem::Vector point(mesh.SpaceDimension());
    for (int be = 0; be < mesh.GetNBE(); be++)
    {
      const int attr = mesh.GetBdrAttribute(be);
      if (attr <= 0 || attr > attr_marker.Size() || !attr_marker[attr - 1])
      {
        continue;
      }
      auto *T = mesh.GetBdrElementTransformation(be);
      const auto geometry = h1_fespace.GetBE(be)->GetGeomType();
      const mfem::IntegrationRule &rule =
          ownership_rules
              ? ownership_rules->Get(geometry)
              : mfem::IntRules.Get(geometry, fem::DefaultIntegrationOrder::Get(*T));
      for (int q = 0; q < rule.GetNPoints(); q++)
      {
        const auto &ip = rule.IntPoint(q);
        T->SetIntPoint(&ip);
        T->Transform(ip, point);
        const auto nearest = tree.Nearest(point);
        const double within =
            1.0 - EdgeDistanceOutsideWeight(std::sqrt(nearest.distance_squared), radius,
                                            data.edge_distance_smoothing);
        if (within <= 0.0)
        {
          continue;
        }
        // PartitionedEdgeCoefficient's foot test; the first piece list holding the foot
        // takes the point.
        const auto &segment = tree.GetSegment(nearest.segment);
        double length2 = 0.0, projection = 0.0;
        for (int d = 0; d < point.Size(); d++)
        {
          const double chord = segment.p1[d] - segment.p0[d];
          length2 += chord * chord;
          projection += (point(d) - segment.p0[d]) * chord;
        }
        const double t = projection / length2;
        const double foot = std::clamp(t, 0.0, 1.0) * std::sqrt(length2);
        std::size_t bin = interface.unclaimed;
        std::size_t held = interface.pieces.size();  // the holding list (none: unclaimed)
        for (std::size_t list = 0; list < interface.pieces.size(); list++)
        {
          const auto found = interface.pieces[list].by_segment.find(nearest.segment);
          if (found == interface.pieces[list].by_segment.end() || found->second.empty())
          {
            continue;
          }
          auto piece = std::upper_bound(found->second.begin(), found->second.end(), foot,
                                        [](double x, const auto &p) { return x < p.a; });
          if (piece == found->second.begin())
          {
            continue;
          }
          piece = std::prev(piece);
          const bool last = std::next(piece) == found->second.end();
          if (foot >= piece->a && (foot < piece->b || (last && foot == piece->b)))
          {
            bin = interface.offsets[list] + piece->group;
            held = list;
            break;
          }
        }
        const double value = PhysicalQuadratureWeight(*T, ip, axisymmetric) *
                             (within * coefficient->Eval(*T, ip));
        local_energies[bin] += value;
        local_energies[interface.within] += value;
        if (t <= 0.0 || t >= 1.0)
        {
          // The clamped foot at a segment end: the vertex quadrant (census only). Ruling
          // 1731: round-off short = EXACT energies: served, what the snap gives back (the
          // served list holds the point and no unsnapped served claim reaches the end);
          // kept, what no kept class holds. Near: every point there (an upper bound).
          const std::pair<std::size_t, int> end{nearest.segment, t <= 0.0 ? 0 : 1};
          const std::size_t served_list = kept_classes.size();
          for (const auto &[check, first, short_inside_held] :
               {std::make_tuple(&interface.served_check, interface.stake,
                                held == served_list),
                std::make_tuple(&interface.kept_check, interface.stake + 2,
                                held >= served_list)})
          {
            const auto found = check->ends.find(end);
            if (found != check->ends.end())
            {
              local_energies[first] +=
                  ((found->second & 1) && !(found->second & 4) && short_inside_held) ? value
                                                                                     : 0.0;
              local_energies[first + 1] += (found->second & 2) ? value : 0.0;
            }
          }
        }
      }
    }
  }
  if (!local_energies.empty())
  {
    Mpi::GlobalSum(static_cast<int>(local_energies.size()), local_energies.data(),
                   E.GetComm());
  }
  std::map<int, InterfaceCoveragePartition> partitions;
  for (const auto &interface : bins)
  {
    auto &partition = partitions[interface.idx];
    for (std::size_t list = 0; list < interface.pieces.size(); list++)
    {
      const auto &groups = interface.pieces[list].groups;
      for (std::size_t group = 0; group < groups.size(); group++)
      {
        const double energy = local_energies[interface.offsets[list] + group];
        if (list < kept_classes.size())
        {
          partition.classes[{kept_classes[list].name, groups[group]}] += energy;
        }
        else if (list == kept_classes.size())
        {
          partition.served[groups[group]] += energy;
        }
        else
        {
          partition.classes[{remainder_classes[list - kept_classes.size() - 1].name,
                             groups[group]}] += energy;
        }
      }
    }
    partition.unclaimed = local_energies[interface.unclaimed];
    partition.within = local_energies[interface.within];
    partition.served_ends = interface.served_check.census;
    partition.served_ends.roundoff_energy = local_energies[interface.stake];
    partition.served_ends.near_energy = local_energies[interface.stake + 1];
    partition.kept_ends = interface.kept_check.census;
    partition.kept_ends.roundoff_energy = local_energies[interface.stake + 2];
    partition.kept_ends.near_energy = local_energies[interface.stake + 3];
  }
  return partitions;
}

std::vector<SurfacePostOperator::UncoveredPerimeterPortion>
SurfacePostOperator::CoverageMissingByFeature(
    const std::vector<UncoveredPerimeterPortion> &uncovered)
{
  std::vector<UncoveredPerimeterPortion> by_feature;
  by_feature.reserve(uncovered.size());
  for (const auto &portion : uncovered)
  {
    by_feature.push_back(
        {portion.p0, portion.p1, "m" + std::to_string(portion.feature), portion.feature});
  }
  return by_feature;
}

void SurfacePostOperator::MoveCoverageServedLabel(
    std::map<int, InterfaceCoveragePartition> &partitions, const std::string &label,
    const std::string &name)
{
  for (auto &[interface, partition] : partitions)
  {
    (void)interface;
    const auto served = partition.served.find(label);
    if (served != partition.served.end())
    {
      partition.classes[{name, label}] += served->second;
      partition.served.erase(served);
    }
  }
}

void SurfacePostOperator::MoveCoverageLoopRemainder(
    std::map<int, InterfaceCoveragePartition> &partitions,
    const std::map<int, std::map<std::string, double>> &remainder)
{
  // A served label holds a model token when one of its '+'-joined tokens equals it.
  auto HoldsToken = [](const std::string &label, const std::string &token)
  {
    std::size_t begin = 0;
    while (begin <= label.size())
    {
      const std::size_t end = std::min(label.find('+', begin), label.size());
      if (label.compare(begin, end - begin, token) == 0 && end - begin == token.size())
      {
        return true;
      }
      begin = end + 1;
    }
    return false;
  };
  for (const auto &[interface, rows] : remainder)
  {
    auto found = partitions.find(interface);
    if (found == partitions.end())
    {
      continue;
    }
    auto &partition = found->second;
    for (const auto &[row, energy] : rows)
    {
      partition.classes[{"LoopRemainder", row}] += energy;
      double holding = 0.0;
      for (const auto &[label, served] : partition.served)
      {
        if (served > 0.0 && HoldsToken(label, row))
        {
          holding += served;
        }
      }
      const double taken = std::min(energy, holding);
      if (holding > 0.0)
      {
        for (auto &[label, served] : partition.served)
        {
          if (served > 0.0 && HoldsToken(label, row))
          {
            served -= taken * (served / holding);
          }
        }
      }
      partition.loop_remainder_unserved += energy - taken;
    }
  }
}

namespace
{

// The three coordinates' bit patterns: the exact-point lookup's key (no tolerance; -0.0 and
// +0.0 differ).
std::array<std::uint64_t, 3> PointBits(const std::array<double, 3> &point)
{
  std::array<std::uint64_t, 3> bits;
  static_assert(sizeof(bits) == sizeof(point));
  std::memcpy(bits.data(), point.data(), sizeof(bits));
  return bits;
}

std::uint64_t PointHash(const std::array<std::uint64_t, 3> &bits)
{
  // splitmix64 finalizer over the combined patterns.
  std::uint64_t h = 0x9E3779B97F4A7C15ull;
  for (const std::uint64_t b : bits)
  {
    std::uint64_t z = h ^ (b + 0x9E3779B97F4A7C15ull + (h << 6) + (h >> 2));
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    h = z ^ (z >> 31);
  }
  return h;
}

}  // namespace

std::uint32_t SurfacePostOperator::CompactPerimeterPortions::FindOrAddPoint(
    const std::array<double, 3> &point)
{
  MFEM_VERIFY(!std::isnan(point[0]) && !std::isnan(point[1]) && !std::isnan(point[2]),
              "A perimeter portion end point has a NaN coordinate!");
  const auto bits = PointBits(point);
  // Keep the table at most half full (power-of-two size); rebuilt after FinishAppending.
  if (2 * (points.size() + 1) > point_slots.size())
  {
    std::size_t size = 16;
    while (size < 4 * (points.size() + 1))
    {
      size *= 2;
    }
    point_slots.assign(size, kEmptySlot);
    for (std::size_t i = 0; i < points.size(); i++)
    {
      std::size_t slot = PointHash(PointBits(points[i])) & (size - 1);
      while (point_slots[slot] != kEmptySlot)
      {
        slot = (slot + 1) & (size - 1);
      }
      point_slots[slot] = static_cast<std::uint32_t>(i);
    }
  }
  const std::size_t mask = point_slots.size() - 1;
  std::size_t slot = PointHash(bits) & mask;
  while (point_slots[slot] != kEmptySlot)
  {
    if (PointBits(points[point_slots[slot]]) == bits)
    {
      return point_slots[slot];
    }
    slot = (slot + 1) & mask;
  }
  MFEM_VERIFY(points.size() < kEmptySlot, "Too many perimeter portion end points!");
  point_slots[slot] = static_cast<std::uint32_t>(points.size());
  points.push_back(point);
  return point_slots[slot];
}

void SurfacePostOperator::CompactPerimeterPortions::Append(const std::array<double, 3> &p0,
                                                           const std::array<double, 3> &p1,
                                                           const std::string &label,
                                                           int feature)
{
  if (label_lookup.empty() && !labels.empty())
  {
    for (std::size_t i = 0; i < labels.size(); i++)
    {
      label_lookup.emplace(labels[i], static_cast<std::int32_t>(i));
    }
  }
  const auto [entry, added] =
      label_lookup.emplace(label, static_cast<std::int32_t>(labels.size()));
  if (added)
  {
    labels.push_back(label);
  }
  Portion portion;
  portion.p0 = FindOrAddPoint(p0);
  portion.p1 = FindOrAddPoint(p1);
  portion.label = entry->second;
  portion.feature = feature;
  portions.push_back(portion);
}

void SurfacePostOperator::CompactPerimeterPortions::FinishAppending()
{
  std::vector<std::uint32_t>().swap(point_slots);
  std::unordered_map<std::string, std::int32_t>().swap(label_lookup);
  points.shrink_to_fit();
  labels.shrink_to_fit();
  portions.shrink_to_fit();
}

std::size_t SurfacePostOperator::CompactPerimeterPortions::HeapBytes() const
{
  std::size_t bytes = points.capacity() * sizeof(points[0]) +
                      labels.capacity() * sizeof(std::string) +
                      portions.capacity() * sizeof(Portion) +
                      point_slots.capacity() * sizeof(std::uint32_t);
  for (const auto &label : labels)
  {
    // A label beyond the small-string buffer holds its characters on the heap (the memory
    // census's rule).
    const char *data = label.data(), *self = reinterpret_cast<const char *>(&label);
    if (!(data >= self && data < self + sizeof(label)))
    {
      bytes += label.capacity() + 1;
    }
  }
  // The label lookup while appending: one node (key, value, link, cached hash) per label
  // and the buckets.
  if (!label_lookup.empty())
  {
    bytes += label_lookup.size() * (sizeof(std::string) + 4 * sizeof(void *)) +
             label_lookup.bucket_count() * sizeof(void *);
  }
  return bytes;
}

std::vector<SurfacePostOperator::UncoveredPerimeterPortion>
SurfacePostOperator::CompactPerimeterPortions::Expand() const
{
  std::vector<UncoveredPerimeterPortion> expanded;
  expanded.reserve(portions.size());
  for (std::size_t i = 0; i < portions.size(); i++)
  {
    expanded.push_back({P0(i), P1(i), Label(i), Feature(i)});
  }
  return expanded;
}

std::map<int, std::map<std::string, double>>
SurfacePostOperator::GetInterfaceWeightedPortionEdgeEnergies(
    const std::set<int> &indices, const GridFunction &E, const GridFunction *D,
    const std::vector<UncoveredPerimeterPortion> &portions,
    const std::vector<double> &weights) const
{
  // OI-MX-LOOPGATE (decision 955). One integral per (interface, label).
  const auto &mesh = *h1_fespace.GetParMesh();
  const int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  std::map<int, std::map<std::string, double>> energies;
  std::vector<double> local_energies;
  std::vector<std::pair<int, std::string>> entries;
  for (const int idx : indices)
  {
    auto it = eps_surfs.find(idx);
    MFEM_VERIFY(it != eps_surfs.end(),
                "Unknown interface dielectric postprocessing index requested!");
    const auto &data = it->second;
    MFEM_VERIFY(!data.edge_distances.empty() && data.edge_distance_tree,
                "Weighted portion edge energies require a matching edge distance!");
    energies.emplace(idx, std::map<std::string, double>{});
    if (portions.empty())
    {
      continue;
    }
    const double radius = data.edge_distances.back();
    const WeightedEdgePieces pieces(*data.edge_distance_tree, portions, weights,
                                    kSignatureParameterToleranceOverRadius * radius);
    if (pieces.groups.empty())
    {
      continue;
    }
    const auto attr_marker = mesh::AttrToMarker(bdr_attr_max, data.attr_list);
    auto coefficient = data.GetCoefficient(E, D, mat_op);
    for (std::size_t group = 0; group < pieces.groups.size(); group++)
    {
      WeightedEdgeCoefficient weighted(*coefficient, *data.edge_distance_tree, pieces,
                                       radius, data.edge_distance_smoothing,
                                       static_cast<int>(group));
      local_energies.push_back(
          GetLocalSurfaceIntegral(weighted, attr_marker, data.ownership_quadrature_order));
      entries.emplace_back(idx, pieces.groups[group]);
    }
  }
  if (!local_energies.empty())
  {
    Mpi::GlobalSum(static_cast<int>(local_energies.size()), local_energies.data(),
                   E.GetComm());
  }
  for (std::size_t i = 0; i < entries.size(); i++)
  {
    energies.at(entries[i].first)[entries[i].second] += local_energies[i];
  }
  return energies;
}

SurfacePostOperator::InterfaceEdgeEnergy
SurfacePostOperator::GetInterfaceOuterElectricFieldEnergy(int idx, const GridFunction &E,
                                                          const GridFunction *D) const
{
  auto energies = GetInterfaceOuterElectricFieldEnergies({idx}, E, D);
  MFEM_ASSERT(energies.size() == 1 && energies.find(idx) != energies.end(),
              "Missing response-corrected interface energy!");
  return energies.at(idx);
}

std::map<int, SurfacePostOperator::InterfaceEdgeEnergy>
SurfacePostOperator::GetInterfaceOuterElectricFieldEnergies(const std::set<int> &indices,
                                                            const GridFunction &E,
                                                            const GridFunction *D) const
{
  const auto &mesh = *h1_fespace.GetParMesh();
  const int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  std::map<int, InterfaceEdgeEnergy> energies;
  std::vector<double> local_energies;
  local_energies.reserve(indices.size());
  for (const int idx : indices)
  {
    auto it = eps_surfs.find(idx);
    MFEM_VERIFY(it != eps_surfs.end(),
                "Unknown interface dielectric postprocessing index requested!");
    const auto &data = it->second;
    MFEM_VERIFY(!data.edge_distances.empty(),
                "Response-corrected interface requires a matching edge distance!");

    const auto attr_marker = mesh::AttrToMarker(bdr_attr_max, data.attr_list);
    auto coefficient = data.GetCoefficient(E, D, mat_op);
    const double distance = data.edge_distances.back();
    EdgeDistanceCoefficient outside(*coefficient, *data.edge_distance_tree, distance,
                                    mfem::infinity(), data.edge_distance_smoothing);
    local_energies.push_back(
        GetLocalSurfaceIntegral(outside, attr_marker, data.ownership_quadrature_order));
    energies.emplace(idx, InterfaceEdgeEnergy{distance, 0.0, 0.0});
  }
  Mpi::GlobalSum(static_cast<int>(local_energies.size()), local_energies.data(),
                 E.GetComm());
  std::size_t i = 0;
  for (auto &[idx, energy] : energies)
  {
    (void)idx;
    energy.energy_outside = local_energies[i++];
  }
  return energies;
}

std::size_t SurfacePostOperator::GetNInterfaceEdgeEntries() const
{
  std::size_t size = 0;
  for (const auto &[idx, data] : eps_surfs)
  {
    size += data.edge_distances.size();
  }
  return size;
}

std::vector<SurfacePostOperator::InterfaceLocalEdgeEnergy>
SurfacePostOperator::GetInterfaceLocalEdgeElectricFieldEnergies(int idx,
                                                                const GridFunction &E,
                                                                const GridFunction *D,
                                                                bool include_volume) const
{
  auto it = eps_surfs.find(idx);
  MFEM_VERIFY(it != eps_surfs.end(),
              "Unknown interface dielectric postprocessing index requested!");
  const auto &data = it->second;
  if (!data.localize_edge_energy)
  {
    return {};
  }
  MFEM_VERIFY(
      !data.ownership || !include_volume,
      "Quadrature ownership partitions surfaces, not localized volume diagnostics!");
  MFEM_VERIFY(!data.region || !include_volume,
              "An interface region filters surfaces, not localized volume diagnostics!");

  const auto &mesh = *h1_fespace.GetParMesh();
  const int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  const auto attr_marker = mesh::AttrToMarker(bdr_attr_max, data.attr_list);
  const bool polarized = data.type != InterfaceDielectric::DEFAULT;
  auto coefficient = polarized ? nullptr : data.GetCoefficient(E, D, mat_op);
  auto normal_coefficient =
      polarized ? data.GetCoefficient(E, D, mat_op, InterfaceDielectricComponent::NORMAL)
                : nullptr;
  auto tangential_coefficient =
      polarized
          ? data.GetCoefficient(E, D, mat_op, InterfaceDielectricComponent::TANGENTIAL)
          : nullptr;
  const std::size_t edge_count = data.edge_distance_tree->Size();
  const std::size_t radius_count = data.edge_distances.size();
  std::vector<double> local_total_energy(edge_count, 0.0);
  std::vector<double> local_total_polarized_energy(2 * edge_count, 0.0);
  std::vector<double> local_energy(2 * radius_count * edge_count, 0.0);
  std::vector<double> local_vertex_energy(radius_count * edge_count, 0.0);
  std::vector<double> local_polarized_energy(4 * radius_count * edge_count, 0.0);
  const double max_distance =
      2.0 * data.edge_distances.back() * (1.0 + data.edge_distance_smoothing);

  mfem::Vector point(mesh.SpaceDimension());
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    const int attr = mesh.GetBdrAttribute(be);
    if (attr <= 0 || attr > attr_marker.Size() || !attr_marker[attr - 1])
    {
      continue;
    }
    auto *T = const_cast<mfem::ParMesh &>(mesh).GetBdrElementTransformation(be);
    const auto *fe = h1_fespace.GetBE(be);
    const auto &ir =
        data.ownership_rule
            ? data.ownership_rule->Get(fe->GetGeomType())
            : mfem::IntRules.Get(fe->GetGeomType(), fem::DefaultIntegrationOrder::Get(*T));
    for (int q = 0; q < ir.GetNPoints(); q++)
    {
      const auto &ip = ir.IntPoint(q);
      T->SetIntPoint(&ip);
      T->Transform(ip, point);
      if (!data.Includes(point))
      {
        continue;
      }
      const auto nearest = data.edge_distance_tree->Nearest(point);
      const double distance = std::sqrt(nearest.distance_squared);
      const double vertex_distance =
          data.edge_distance_tree->DistanceAlongEdgeToNonregularVertex(point,
                                                                       nearest.segment);
      const double integration_weight = PhysicalQuadratureWeight(*T, ip, axisymmetric);
      const double normal_energy =
          polarized ? integration_weight * normal_coefficient->Eval(*T, ip) : 0.0;
      const double tangential_energy =
          polarized ? integration_weight * tangential_coefficient->Eval(*T, ip) : 0.0;
      const double energy = polarized ? normal_energy + tangential_energy
                                      : integration_weight * coefficient->Eval(*T, ip);
      local_total_energy[nearest.segment] += energy;
      local_total_polarized_energy[2 * nearest.segment] += normal_energy;
      local_total_polarized_energy[2 * nearest.segment + 1] += tangential_energy;
      if (distance >= max_distance)
      {
        continue;
      }

      for (std::size_t radius = 0; radius < radius_count; radius++)
      {
        const std::size_t entry = radius * edge_count + nearest.segment;
        const double matching_distance = data.edge_distances[radius];
        const double outside_weight = EdgeDistanceOutsideWeight(
            distance, matching_distance, data.edge_distance_smoothing);
        const double annulus_weight =
            EdgeDistanceWindowWeight(distance, matching_distance, 2.0 * matching_distance,
                                     data.edge_distance_smoothing);
        const double vertex_weight =
            1.0 - EdgeDistanceOutsideWeight(vertex_distance, matching_distance,
                                            data.edge_distance_smoothing);
        local_energy[2 * entry] += (1.0 - outside_weight) * energy;
        local_energy[2 * entry + 1] += annulus_weight * energy;
        local_vertex_energy[entry] += (1.0 - outside_weight) * vertex_weight * energy;
        local_polarized_energy[4 * entry] += (1.0 - outside_weight) * normal_energy;
        local_polarized_energy[4 * entry + 1] += annulus_weight * normal_energy;
        local_polarized_energy[4 * entry + 2] += (1.0 - outside_weight) * tangential_energy;
        local_polarized_energy[4 * entry + 3] += annulus_weight * tangential_energy;
      }
    }
  }
  Mpi::GlobalSum(static_cast<int>(local_total_energy.size()), local_total_energy.data(),
                 E.GetComm());
  Mpi::GlobalSum(static_cast<int>(local_total_polarized_energy.size()),
                 local_total_polarized_energy.data(), E.GetComm());
  Mpi::GlobalSum(static_cast<int>(local_energy.size()), local_energy.data(), E.GetComm());
  Mpi::GlobalSum(static_cast<int>(local_vertex_energy.size()), local_vertex_energy.data(),
                 E.GetComm());
  Mpi::GlobalSum(static_cast<int>(local_polarized_energy.size()),
                 local_polarized_energy.data(), E.GetComm());
  const LocalVolumeEdgeEnergyCache *local_volume_energy = nullptr;
  if (include_volume)
  {
    local_volume_energy = &GetLocalVolumeEdgeElectricFieldEnergies(data, E);
  }

  std::vector<InterfaceLocalEdgeEnergy> energies;
  energies.reserve(local_energy.size());
  for (std::size_t radius = 0; radius < radius_count; radius++)
  {
    for (std::size_t edge = 0; edge < edge_count; edge++)
    {
      const auto &segment = data.edge_distance_tree->GetSegment(edge);
      double length_squared = 0.0;
      for (int d = 0; d < mesh.SpaceDimension(); d++)
      {
        const double delta = segment.p1[d] - segment.p0[d];
        length_squared += delta * delta;
      }
      const std::size_t entry = radius * edge_count + edge;
      const auto frame =
          BuildEdgeFrame(segment,
                         data.edge_distance_tree->HasProcessNormals()
                             ? data.edge_distance_tree->GetProcessNormal(edge)
                             : data.edge_frame_normal,
                         mesh.SpaceDimension());
      int component = -1;
      int chain = -1;
      std::array<int, 2> vertex_types{-1, -1};
      if (data.edge_distance_tree->HasMetadata())
      {
        const auto &metadata = data.edge_distance_tree->GetMetadata(edge);
        component = metadata.component;
        chain = metadata.chain;
        vertex_types = metadata.vertex_types;
      }
      InterfaceLocalEdgeEnergy energy{
          static_cast<int>(edge + 1),
          data.edge_distance_tree->HasMetadata(),
          component,
          chain,
          vertex_types,
          frame.normal,
          segment.p0,
          segment.p1,
          std::sqrt(length_squared),
          data.edge_distances[radius],
          local_total_energy[edge],
          local_energy[2 * entry],
          local_energy[2 * entry + 1],
          local_vertex_energy[entry],
          {local_total_polarized_energy[2 * edge],
           local_total_polarized_energy[2 * edge + 1]},
          {local_polarized_energy[4 * entry], local_polarized_energy[4 * entry + 2]},
          {local_polarized_energy[4 * entry + 1], local_polarized_energy[4 * entry + 3]},
          local_volume_energy ? local_volume_energy->energies[entry] : 0.0,
          local_volume_energy ? local_volume_energy->vertex_energies[entry] : 0.0,
          {}};
      if (local_volume_energy)
      {
        for (std::size_t component = 0; component < 6; component++)
        {
          energy.energy_volume_annulus_polarized[component] =
              local_volume_energy->polarized_energies[component][entry];
        }
      }
      energies.push_back(std::move(energy));
    }
  }
  return energies;
}

template <InterfaceDielectric Type>
std::vector<SurfacePostOperator::InterfaceResponseMatrix>
SurfacePostOperator::GetInterfaceElectricFieldEnergyMatricesImpl(
    const InterfaceDielectricData &data, const std::vector<const GridFunction *> &E,
    const std::vector<const GridFunction *> &D, int interface_index, int quadrature_extra,
    std::vector<InterfaceQuadratureRule> *quadrature_rules) const
{
  MFEM_VERIFY(quadrature_extra >= 0 && quadrature_extra <= 12,
              "Surface quadrature extra order must be in [0,12]!");
  MFEM_VERIFY(quadrature_extra == 0 || !data.ownership_rule,
              "Surface quadrature override is unsupported with an ownership rule!");
  MFEM_VERIFY(!E.empty() && (D.empty() || D.size() == E.size()),
              "Invalid batched interface response fields!");
  MFEM_VERIFY(!data.flux_recovery || !D.empty(),
              "Batched interface response requires recovered electric flux!");
  const auto &mesh = *h1_fespace.GetParMesh();
  const int sdim = mesh.SpaceDimension();
  const int basis_size = static_cast<int>(E.size());
  const int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  const auto attr_marker = mesh::AttrToMarker(bdr_attr_max, data.attr_list);

  std::vector<std::unique_ptr<InterfaceDielectricCoefficient<Type>>> evaluators;
  evaluators.reserve(E.size());
  for (std::size_t i = 0; i < E.size(); i++)
  {
    evaluators.push_back(std::make_unique<InterfaceDielectricCoefficient<Type>>(
        *E[i], mat_op, data.t, data.epsilon, &data.sum_sides_marker,
        data.flux_recovery ? D[i] : nullptr));
  }

  std::vector<double> weights, distances, amplitudes;
  mfem::Vector point(sdim), field;
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    const int attr = mesh.GetBdrAttribute(be);
    if (attr <= 0 || attr > attr_marker.Size() || !attr_marker[attr - 1])
    {
      continue;
    }
    auto *T = const_cast<mfem::ParMesh &>(mesh).GetBdrElementTransformation(be);
    const auto *fe = h1_fespace.GetBE(be);
    const int requested =
        data.ownership_rule ? -1 : fem::DefaultIntegrationOrder::Get(*T) + quadrature_extra;
    const auto &ir = data.ownership_rule ? data.ownership_rule->Get(fe->GetGeomType())
                                         : mfem::IntRules.Get(fe->GetGeomType(), requested);
    if (quadrature_extra > 0 || quadrature_rules)
    {
      double minimum_weight = std::numeric_limits<double>::infinity();
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
        minimum_weight = std::min(minimum_weight, ir.IntPoint(q).weight);
        MFEM_VERIFY(std::isfinite(ir.IntPoint(q).weight) && ir.IntPoint(q).weight > 0.0,
                    "Nonpositive surface matrix quadrature weight for geometry "
                        << fe->GetGeomType() << ", requested order " << requested
                        << "; no rule substitution is permitted!");
      }
      if (quadrature_rules)
      {
        auto found = std::find_if(quadrature_rules->begin(), quadrature_rules->end(),
                                  [&](const auto &r)
                                  {
                                    return r.interface_index == interface_index &&
                                           r.geometry == fe->GetGeomType() &&
                                           r.requested_order == requested;
                                  });
        if (found == quadrature_rules->end())
        {
          quadrature_rules->push_back({interface_index, fe->GetGeomType(), requested,
                                       ir.GetOrder(), ir.GetNPoints(), 1, minimum_weight,
                                       bool(data.ownership_rule)});
        }
        else
        {
          found->local_faces++;
        }
      }
    }
    for (int q = 0; q < ir.GetNPoints(); q++)
    {
      const auto &ip = ir.IntPoint(q);
      T->SetIntPoint(&ip);
      T->Transform(ip, point);
      if (!data.Includes(point))
      {
        continue;
      }
      const auto nearest = data.edge_distance_tree->Nearest(point);
      distances.push_back(std::sqrt(nearest.distance_squared));
      MFEM_VERIFY(!data.ownership || ip.weight >= 0.0,
                  "Ownership matrix quadrature requires nonnegative weights!");
      const double weight = PhysicalQuadratureWeight(*T, ip, axisymmetric);
      MFEM_VERIFY(!(quadrature_extra > 0 || quadrature_rules) ||
                      (std::isfinite(weight) && weight > 0.0),
                  "Invalid physical surface quadrature weight!");
      weights.push_back(weight);
      for (const auto &evaluator : evaluators)
      {
        evaluator->EvalEnergyField(*T, ip, field);
        MFEM_ASSERT(field.Size() == sdim + 1,
                    "Invalid batched interface energy field size!");
        amplitudes.insert(amplitudes.end(), field.begin(), field.end());
      }
    }
  }

  const int sample_count = static_cast<int>(weights.size());
  auto Assemble = [&](const std::vector<double> *inside_weights, mfem::DenseMatrix &normal,
                      mfem::DenseMatrix &tangential)
  {
    normal.SetSize(basis_size);
    normal = 0.0;
    tangential.SetSize(basis_size);
    tangential = 0.0;
    if (sample_count == 0)
    {
      return;
    }
    mfem::DenseMatrix normal_samples(basis_size, sample_count);
    mfem::DenseMatrix tangential_samples(basis_size, sdim * sample_count);
    for (int sample = 0; sample < sample_count; sample++)
    {
      const double window = inside_weights ? (*inside_weights)[sample] : 1.0;
      const double scale = std::sqrt(0.5 * weights[sample] * window);
      for (int basis = 0; basis < basis_size; basis++)
      {
        const std::size_t offset =
            (static_cast<std::size_t>(sample) * basis_size + basis) * (sdim + 1);
        normal_samples(basis, sample) = scale * amplitudes[offset];
        for (int d = 0; d < sdim; d++)
        {
          tangential_samples(basis, sdim * sample + d) = scale * amplitudes[offset + d + 1];
        }
      }
    }
    mfem::AddMult_a_AAt(1.0, normal_samples, normal);
    mfem::AddMult_a_AAt(1.0, tangential_samples, tangential);
  };
  auto Reduce = [&](mfem::DenseMatrix &matrix)
  {
    Mpi::GlobalSum(matrix.Height() * matrix.Width(), matrix.GetData(),
                   E.front()->GetComm());
  };

  mfem::DenseMatrix total_normal, total_tangential;
  Assemble(nullptr, total_normal, total_tangential);
  Reduce(total_normal);
  Reduce(total_tangential);

  std::vector<InterfaceResponseMatrix> result;
  result.reserve(data.edge_distances.size());
  std::vector<double> inside_weights(sample_count);
  for (const double radius : data.edge_distances)
  {
    for (int sample = 0; sample < sample_count; sample++)
    {
      inside_weights[sample] =
          1.0 - EdgeDistanceOutsideWeight(distances[sample], radius,
                                          data.edge_distance_smoothing);
    }
    InterfaceResponseMatrix entry;
    entry.distance = radius;
    entry.energy_total_normal = total_normal;
    entry.energy_total_tangential = total_tangential;
    Assemble(&inside_weights, entry.energy_inside_normal, entry.energy_inside_tangential);
    Reduce(entry.energy_inside_normal);
    Reduce(entry.energy_inside_tangential);
    entry.energy_total = entry.energy_total_normal;
    entry.energy_total += entry.energy_total_tangential;
    entry.energy_inside = entry.energy_inside_normal;
    entry.energy_inside += entry.energy_inside_tangential;
    result.push_back(std::move(entry));
  }
  return result;
}

std::map<int, std::vector<SurfacePostOperator::InterfaceResponseMatrix>>
SurfacePostOperator::GetInterfaceElectricFieldEnergyMatrices(
    const std::vector<const GridFunction *> &E, const std::vector<const GridFunction *> &D,
    int quadrature_extra, std::vector<InterfaceQuadratureRule> *quadrature_rules) const
{
  std::map<int, std::vector<InterfaceResponseMatrix>> result;
  for (const auto &[idx, data] : eps_surfs)
  {
    if (!data.localize_edge_energy)
    {
      continue;
    }
    switch (data.type)
    {
      case InterfaceDielectric::DEFAULT:
        result.emplace(
            idx, GetInterfaceElectricFieldEnergyMatricesImpl<InterfaceDielectric::DEFAULT>(
                     data, E, D, idx, quadrature_extra, quadrature_rules));
        break;
      case InterfaceDielectric::MA:
        result.emplace(idx,
                       GetInterfaceElectricFieldEnergyMatricesImpl<InterfaceDielectric::MA>(
                           data, E, D, idx, quadrature_extra, quadrature_rules));
        break;
      case InterfaceDielectric::MS:
        result.emplace(idx,
                       GetInterfaceElectricFieldEnergyMatricesImpl<InterfaceDielectric::MS>(
                           data, E, D, idx, quadrature_extra, quadrature_rules));
        break;
      case InterfaceDielectric::SA:
        result.emplace(idx,
                       GetInterfaceElectricFieldEnergyMatricesImpl<InterfaceDielectric::SA>(
                           data, E, D, idx, quadrature_extra, quadrature_rules));
        break;
    }
  }
  return result;
}

namespace
{

// Cache-blocked symmetric rank-k update of the response Gram: for basis rows i >= j,
// normal(i, j) += sum_s w_s a_i(s, 0) a_j(s, 0) and tangential(i, j) += sum_s w_s
// sum_{c >= 1} a_i(s, c) a_j(s, c), a_i = rows + i x row_size (sample-major, `components`
// amplitudes per sample), w_s = sample_weights[s] or 1. Rows are visited in blocks of
// basis fields over chunks of samples so every product streams contiguous memory.
void AccumulateInterfaceResponseGram(const double *rows, int basis_size,
                                     std::size_t sample_count, int components,
                                     const double *sample_weights,
                                     mfem::DenseMatrix &normal,
                                     mfem::DenseMatrix &tangential)
{
  constexpr int basis_block = 8;
  constexpr std::size_t sample_chunk = 2048;
  const std::size_t row_size = sample_count * static_cast<std::size_t>(components);
  for (int ib = 0; ib < basis_size; ib += basis_block)
  {
    const int ie = std::min(ib + basis_block, basis_size);
    for (int jb = 0; jb <= ib; jb += basis_block)
    {
      const int je = std::min(jb + basis_block, basis_size);
      for (std::size_t sb = 0; sb < sample_count; sb += sample_chunk)
      {
        const std::size_t se = std::min(sb + sample_chunk, sample_count);
        for (int i = ib; i < ie; i++)
        {
          const double *a = rows + static_cast<std::size_t>(i) * row_size;
          for (int j = jb; j < std::min(je, i + 1); j++)
          {
            const double *b = rows + static_cast<std::size_t>(j) * row_size;
            double n = 0.0, t = 0.0;
            for (std::size_t sample = sb; sample < se; sample++)
            {
              const double *as = a + sample * components;
              const double *bs = b + sample * components;
              double tangential_product = 0.0;
              for (int c = 1; c < components; c++)
              {
                tangential_product += as[c] * bs[c];
              }
              const double w = sample_weights ? sample_weights[sample] : 1.0;
              n += w * (as[0] * bs[0]);
              t += w * tangential_product;
            }
            normal(i, j) += n;
            tangential(i, j) += t;
          }
        }
      }
    }
  }
  for (int i = 0; i < basis_size; i++)
  {
    for (int j = 0; j < i; j++)
    {
      normal(j, i) = normal(i, j);
      tangential(j, i) = tangential(i, j);
    }
  }
}

}  // namespace

std::vector<SurfacePostOperator::InterfaceResponseSamples>
SurfacePostOperator::CacheInterfaceResponseSamples() const
{
  const auto &mesh = *h1_fespace.GetParMesh();
  const int sdim = mesh.SpaceDimension();
  const int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  std::vector<InterfaceResponseSamples> result;
  mfem::Vector point(sdim);
  for (const auto &[idx, data] : eps_surfs)
  {
    if (!data.localize_edge_energy)
    {
      continue;
    }
    InterfaceResponseSamples samples;
    samples.interface_index = idx;
    samples.data = &data;
    samples.components = sdim + 1;
    const auto attr_marker = mesh::AttrToMarker(bdr_attr_max, data.attr_list);
    for (int be = 0; be < mesh.GetNBE(); be++)
    {
      const int attr = mesh.GetBdrAttribute(be);
      if (attr <= 0 || attr > attr_marker.Size() || !attr_marker[attr - 1])
      {
        continue;
      }
      auto *T = const_cast<mfem::ParMesh &>(mesh).GetBdrElementTransformation(be);
      const auto *fe = h1_fespace.GetBE(be);
      const auto &ir = data.ownership_rule
                           ? data.ownership_rule->Get(fe->GetGeomType())
                           : mfem::IntRules.Get(fe->GetGeomType(),
                                                fem::DefaultIntegrationOrder::Get(*T));
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
        const auto &ip = ir.IntPoint(q);
        T->SetIntPoint(&ip);
        T->Transform(ip, point);
        if (!data.Includes(point))
        {
          continue;
        }
        const auto nearest = data.edge_distance_tree->Nearest(point);
        MFEM_VERIFY(!data.ownership || ip.weight >= 0.0,
                    "Ownership matrix quadrature requires nonnegative weights!");
        samples.elements.push_back(be);
        samples.rules.push_back(&ir);
        samples.points.push_back(q);
        samples.weights.push_back(PhysicalQuadratureWeight(*T, ip, axisymmetric));
        samples.distances.push_back(std::sqrt(nearest.distance_squared));
      }
    }
    result.push_back(std::move(samples));
  }
  return result;
}

template <InterfaceDielectric Type>
void SurfacePostOperator::EvaluateInterfaceResponseRowImpl(
    const InterfaceResponseSamples &samples, const GridFunction &E, const GridFunction *D,
    double *row) const
{
  const auto &data = *samples.data;
  MFEM_VERIFY(!data.flux_recovery || D,
              "Streamed interface response requires recovered electric flux!");
  const auto &mesh = *h1_fespace.GetParMesh();
  InterfaceDielectricCoefficient<Type> evaluator(E, mat_op, data.t, data.epsilon,
                                                 &data.sum_sides_marker,
                                                 data.flux_recovery ? D : nullptr);
  mfem::Vector field;
  mfem::ElementTransformation *T = nullptr;
  int element = -1;
  for (std::size_t sample = 0; sample < samples.Count(); sample++)
  {
    if (samples.elements[sample] != element)
    {
      element = samples.elements[sample];
      T = const_cast<mfem::ParMesh &>(mesh).GetBdrElementTransformation(element);
    }
    const auto &ip = samples.rules[sample]->IntPoint(samples.points[sample]);
    T->SetIntPoint(&ip);
    evaluator.EvalEnergyField(*T, ip, field);
    MFEM_ASSERT(field.Size() == samples.components,
                "Invalid streamed interface energy field size!");
    const double scale = std::sqrt(0.5 * samples.weights[sample]);
    double *out = row + sample * static_cast<std::size_t>(samples.components);
    for (int c = 0; c < samples.components; c++)
    {
      out[c] = scale * field[c];
    }
  }
}

void SurfacePostOperator::EvaluateInterfaceResponseRow(
    const InterfaceResponseSamples &samples, const GridFunction &E, const GridFunction *D,
    double *row) const
{
  switch (samples.data->type)
  {
    case InterfaceDielectric::DEFAULT:
      EvaluateInterfaceResponseRowImpl<InterfaceDielectric::DEFAULT>(samples, E, D, row);
      break;
    case InterfaceDielectric::MA:
      EvaluateInterfaceResponseRowImpl<InterfaceDielectric::MA>(samples, E, D, row);
      break;
    case InterfaceDielectric::MS:
      EvaluateInterfaceResponseRowImpl<InterfaceDielectric::MS>(samples, E, D, row);
      break;
    case InterfaceDielectric::SA:
      EvaluateInterfaceResponseRowImpl<InterfaceDielectric::SA>(samples, E, D, row);
      break;
  }
}

template <InterfaceDielectric Type>
int SurfacePostOperator::CountSummedInterfaceResponseSamples(
    const InterfaceResponseSamples &samples, const GridFunction &E) const
{
  const auto &data = *samples.data;
  const auto &mesh = *h1_fespace.GetParMesh();
  InterfaceDielectricCoefficient<Type> probe(E, mat_op, data.t, data.epsilon,
                                             &data.sum_sides_marker);
  mfem::ElementTransformation *T = nullptr;
  int element = -1, count = 0;
  for (std::size_t sample = 0; sample < samples.Count(); sample++)
  {
    if (samples.elements[sample] != element)
    {
      element = samples.elements[sample];
      T = const_cast<mfem::ParMesh &>(mesh).GetBdrElementTransformation(element);
    }
    const auto &ip = samples.rules[sample]->IntPoint(samples.points[sample]);
    T->SetIntPoint(&ip);
    count += probe.IsSummedSample(*T, ip);
  }
  return count;
}

void SurfacePostOperator::VerifyBatchedInterfaceResponse(const GridFunction &E) const
{
  // The neighbouring elements of shared sheet faces (the other side on another process).
  auto &mesh = const_cast<mfem::ParMesh &>(*h1_fespace.GetParMesh());
  mesh.ExchangeFaceNbrData();
  for (const auto &samples : CacheInterfaceResponseSamples())
  {
    int summed = 0;
    switch (samples.data->type)
    {
      case InterfaceDielectric::DEFAULT:
        summed =
            CountSummedInterfaceResponseSamples<InterfaceDielectric::DEFAULT>(samples, E);
        break;
      case InterfaceDielectric::MA:
        summed = CountSummedInterfaceResponseSamples<InterfaceDielectric::MA>(samples, E);
        break;
      case InterfaceDielectric::MS:
        summed = CountSummedInterfaceResponseSamples<InterfaceDielectric::MS>(samples, E);
        break;
      case InterfaceDielectric::SA:
        break;
    }
    Mpi::GlobalSum(1, &summed, mesh.GetComm());
    MFEM_VERIFY(summed == 0,
                "Batched interface response fields (the aggregated response matrix, the "
                "response archive) cannot represent the summed energies of both sides of "
                "an interior sheet: interface dielectric "
                    << samples.interface_index << " has " << summed
                    << " quadrature samples on interior sheet faces whose two sides both "
                       "qualify (refused before the solve; SD-401 (d))!");
  }
}

std::vector<SurfacePostOperator::InterfaceResponseMatrix>
SurfacePostOperator::AssembleInterfaceResponseMatrices(
    const InterfaceResponseSamples &samples, const double *rows, int basis_size,
    MPI_Comm comm) const
{
  MFEM_VERIFY(basis_size > 0, "Streamed interface response requires basis fields!");
  const auto &data = *samples.data;
  auto Assemble = [&](const std::vector<double> *inside_weights, mfem::DenseMatrix &normal,
                      mfem::DenseMatrix &tangential)
  {
    normal.SetSize(basis_size);
    normal = 0.0;
    tangential.SetSize(basis_size);
    tangential = 0.0;
    if (samples.Count() > 0)
    {
      AccumulateInterfaceResponseGram(rows, basis_size, samples.Count(), samples.components,
                                      inside_weights ? inside_weights->data() : nullptr,
                                      normal, tangential);
    }
    Mpi::GlobalSum(normal.Height() * normal.Width(), normal.GetData(), comm);
    Mpi::GlobalSum(tangential.Height() * tangential.Width(), tangential.GetData(), comm);
  };

  mfem::DenseMatrix total_normal, total_tangential;
  Assemble(nullptr, total_normal, total_tangential);

  std::vector<InterfaceResponseMatrix> result;
  result.reserve(data.edge_distances.size());
  std::vector<double> inside_weights(samples.Count());
  for (const double radius : data.edge_distances)
  {
    for (std::size_t sample = 0; sample < samples.Count(); sample++)
    {
      inside_weights[sample] =
          1.0 - EdgeDistanceOutsideWeight(samples.distances[sample], radius,
                                          data.edge_distance_smoothing);
    }
    InterfaceResponseMatrix entry;
    entry.distance = radius;
    entry.energy_total_normal = total_normal;
    entry.energy_total_tangential = total_tangential;
    Assemble(&inside_weights, entry.energy_inside_normal, entry.energy_inside_tangential);
    entry.energy_total = entry.energy_total_normal;
    entry.energy_total += entry.energy_total_tangential;
    entry.energy_inside = entry.energy_inside_normal;
    entry.energy_inside += entry.energy_inside_tangential;
    result.push_back(std::move(entry));
  }
  return result;
}

const SurfacePostOperator::LocalVolumeEdgeEnergyCache &
SurfacePostOperator::GetLocalVolumeEdgeElectricFieldEnergies(
    const InterfaceDielectricData &data, const GridFunction &E) const
{
  for (const auto &entry : local_volume_edge_energy_cache)
  {
    if (entry.edge_distance_tree == data.edge_distance_tree.get() &&
        entry.edge_distances == data.edge_distances &&
        entry.edge_distance_smoothing == data.edge_distance_smoothing &&
        entry.edge_frame_normal == data.edge_frame_normal)
    {
      return entry;
    }
  }

  const auto &mesh = *E.ParFESpace()->GetParMesh();
  const std::size_t edge_count = data.edge_distance_tree->Size();
  const std::size_t radius_count = data.edge_distances.size();
  auto &entry = local_volume_edge_energy_cache.emplace_back(
      LocalVolumeEdgeEnergyCache{data.edge_distance_tree.get(),
                                 data.edge_distances,
                                 data.edge_distance_smoothing,
                                 data.edge_frame_normal,
                                 std::vector<double>(radius_count * edge_count, 0.0),
                                 std::vector<double>(radius_count * edge_count, 0.0),
                                 {}});
  for (auto &component : entry.polarized_energies)
  {
    component.assign(radius_count * edge_count, 0.0);
  }
  const double max_distance =
      2.0 * data.edge_distances.back() * (1.0 + data.edge_distance_smoothing);

  std::vector<EdgeFrame> frames;
  frames.reserve(edge_count);
  for (std::size_t edge = 0; edge < edge_count; edge++)
  {
    const auto &normal = data.edge_distance_tree->HasProcessNormals()
                             ? data.edge_distance_tree->GetProcessNormal(edge)
                             : data.edge_frame_normal;
    frames.push_back(BuildEdgeFrame(data.edge_distance_tree->GetSegment(edge), normal,
                                    mesh.SpaceDimension()));
  }

  mfem::Vector point(mesh.SpaceDimension());
  mfem::Vector field(mesh.SpaceDimension()), displacement(mesh.SpaceDimension());
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    mfem::IsoparametricTransformation T;
    mesh.GetElementTransformation(e, &T);
    const auto &ir =
        mfem::IntRules.Get(T.GetGeometryType(), fem::DefaultIntegrationOrder::Get(T));
    for (int q = 0; q < ir.GetNPoints(); q++)
    {
      const auto &ip = ir.IntPoint(q);
      T.SetIntPoint(&ip);
      T.Transform(ip, point);
      const auto nearest = data.edge_distance_tree->Nearest(point);
      const double distance = std::sqrt(nearest.distance_squared);
      if (distance >= max_distance)
      {
        continue;
      }

      const auto &frame = frames[nearest.segment];
      const auto &segment = data.edge_distance_tree->GetSegment(nearest.segment);
      const double vertex_distance =
          data.edge_distance_tree->DistanceAlongEdgeToNonregularVertex(point,
                                                                       nearest.segment);
      std::array<double, 6> component_energy{};
      auto AddFieldEnergy = [&](const mfem::ParGridFunction &grid_function)
      {
        grid_function.GetVectorValue(T, ip, field);
        mat_op.GetPermittivityReal(T.Attribute).Mult(field, displacement);
        const auto polarized =
            GetPolarizedEdgeEnergyDensity(point, segment, frame, field, displacement);
        for (std::size_t component = 0; component < 6; component++)
        {
          component_energy[component] += polarized[component];
        }
      };
      AddFieldEnergy(E.Real());
      if (E.HasImag())
      {
        AddFieldEnergy(E.Imag());
      }
      const double integration_weight = PhysicalQuadratureWeight(T, ip, axisymmetric);
      for (std::size_t radius = 0; radius < radius_count; radius++)
      {
        const double matching_distance = data.edge_distances[radius];
        const double annulus_weight =
            EdgeDistanceWindowWeight(distance, matching_distance, 2.0 * matching_distance,
                                     data.edge_distance_smoothing);
        const double vertex_weight =
            1.0 - EdgeDistanceOutsideWeight(vertex_distance, matching_distance,
                                            data.edge_distance_smoothing);
        const std::size_t result = radius * edge_count + nearest.segment;
        for (std::size_t component = 0; component < 6; component++)
        {
          const double energy =
              annulus_weight * integration_weight * component_energy[component];
          entry.polarized_energies[component][result] += energy;
          entry.energies[result] += energy;
          entry.vertex_energies[result] += vertex_weight * energy;
        }
      }
    }
  }
  Mpi::GlobalSum(static_cast<int>(entry.energies.size()), entry.energies.data(),
                 E.GetComm());
  Mpi::GlobalSum(static_cast<int>(entry.vertex_energies.size()),
                 entry.vertex_energies.data(), E.GetComm());
  for (auto &component : entry.polarized_energies)
  {
    Mpi::GlobalSum(static_cast<int>(component.size()), component.data(), E.GetComm());
  }
  return entry;
}

std::size_t SurfacePostOperator::GetNInterfaceLocalEdgeEntries() const
{
  std::size_t size = 0;
  for (const auto &[idx, data] : eps_surfs)
  {
    if (data.localize_edge_energy && data.save_local_edge_energy)
    {
      size += data.edge_distances.size() * data.edge_distance_tree->Size();
    }
  }
  return size;
}

bool SurfacePostOperator::NeedsFluxRecovery() const
{
  return std::any_of(eps_surfs.begin(), eps_surfs.end(),
                     [](const auto &entry) { return entry.second.flux_recovery; });
}

void SurfacePostOperator::ResetInterfaceLocalEdgeEnergyCache() const
{
  local_volume_edge_energy_cache.clear();
}

double SurfacePostOperator::GetLocalSurfaceIntegral(mfem::Coefficient &f,
                                                    const mfem::Array<int> &attr_marker,
                                                    int quadrature_order) const
{
  // Integrate the coefficient over the boundary attributes making up this surface index.
  if (quadrature_order > 0)
  {
    const OwnershipQuadrature rules(quadrature_order);
    auto &mesh = *h1_fespace.GetParMesh();
    double integral = 0.0;
    for (int be = 0; be < mesh.GetNBE(); be++)
    {
      const int attr = mesh.GetBdrAttribute(be);
      if (attr <= 0 || attr > attr_marker.Size() || !attr_marker[attr - 1])
      {
        continue;
      }
      auto *T = mesh.GetBdrElementTransformation(be);
      const auto &rule = rules.Get(h1_fespace.GetBE(be)->GetGeomType());
      for (int q = 0; q < rule.GetNPoints(); q++)
      {
        const auto &point = rule.IntPoint(q);
        MFEM_VERIFY(point.weight >= 0.0,
                    "Ownership quadrature requires nonnegative weights!");
        T->SetIntPoint(&point);
        integral += PhysicalQuadratureWeight(*T, point, axisymmetric) * f.Eval(*T, point);
      }
    }
    return integral;
  }
  AxisymmetricMeasureCoefficient axisymmetric_f(f);
  mfem::LinearForm s(&h1_fespace);
  s.AddBoundaryIntegrator(
      new BoundaryLFIntegrator(
          axisymmetric ? static_cast<mfem::Coefficient &>(axisymmetric_f) : f),
      const_cast<mfem::Array<int> &>(attr_marker));
  s.UseFastAssembly(false);
  s.UseDevice(false);
  s.Assemble();
  s.UseDevice(true);
  return linalg::LocalSum(s);
}

std::vector<std::array<std::complex<double>, 3>> SurfacePostOperator::GetFarFieldrE(
    const std::vector<std::pair<double, double>> &theta_phi_pairs, const GridFunction &E,
    const GridFunction &B, std::complex<double> omega) const
{
  if (theta_phi_pairs.empty())
    return {};
  MFEM_VERIFY(nd_fespace.GetParMesh()->SpaceDimension() == 3,
              "Far-field computation is only available for 3D simulations!");
  BlockTimer bt0(Timer::POSTPRO_FARFIELD);

  // Compute target unit vectors from the given theta and phis.
  std::vector<std::array<double, 3>> r_naughts;
  r_naughts.reserve(theta_phi_pairs.size());
  for (const auto &[theta, phi] : theta_phi_pairs)
  {
    r_naughts.emplace_back(std::array<double, 3>{
        std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)});
  }

  const auto &mesh = *nd_fespace.GetParMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> attr_marker = mesh::AttrToMarker(bdr_attr_max, farfield.attr_list);

  // Use the libCEED surface functional path when supported, avoiding per-point host
  // coefficient evaluation in the MFEM path. Observation directions are part of the
  // assembled operator and must therefore invalidate the cached functional when changed.
  const bool use_ceed = E.HasImag() && B.HasImag();
  if (use_ceed && (!farfield_func || farfield_func_dirs != r_naughts))
  {
    farfield_func = std::make_unique<SurfaceFunctional>(
        mat_op.GetMesh(), attr_marker, *E.ParFESpace(), *B.ParFESpace(), mat_op, r_naughts);
    farfield_func_dirs = r_naughts;
  }
  if (use_ceed && farfield_func && farfield_func->IsValid())
  {
    return farfield_func->EvalFarField(E, B, omega);
  }
  if (use_ceed)
  {
    RequireCeedSurfaceFunctional(farfield_func.get(), "3D far-field postprocessing");
  }

  MFEM_VERIFY(!farfield.interior, "Far-field extraction on interior surfaces requires "
                                  "complex-valued fields (libCEED path)!");

  // Integrate. Each MPI process computes its contribution and we will reduce
  // everything at the end. We make them std::vector<std::array<double, 3>>
  // because we want a very simple memory layout so that we can reduce
  // everything with two MPI calls.
  std::vector<std::array<double, 3>> integrals_r(theta_phi_pairs.size());
  std::vector<std::array<double, 3>> integrals_i(theta_phi_pairs.size());

  for (int i = 0; i < mesh.GetNBE(); i++)
  {
    if (!attr_marker[mesh.GetBdrAttribute(i) - 1])
      continue;

    auto *T = const_cast<mfem::ParMesh &>(mesh).GetBdrElementTransformation(i);
    const auto *fe = nd_fespace.GetBE(i);
    const auto *ir =
        &mfem::IntRules.Get(fe->GetGeomType(), fem::DefaultIntegrationOrder::Get(*T));

    AddStrattonChuIntegrandAtElement(E, B, mat_op, omega.real(), omega.imag(), r_naughts,
                                     *T, *ir, integrals_r, integrals_i);
  }

  double *data_r_ptr = integrals_r.data()->data();
  double *data_i_ptr = integrals_i.data()->data();
  std::size_t total_elements = integrals_r.size() * 3;
  Mpi::GlobalSum(total_elements, data_i_ptr, E.GetComm());
  Mpi::GlobalSum(total_elements, data_r_ptr, E.GetComm());

  // Finally, we apply cross product to reduced integrals and package the result
  // in a neatly accessible vector of arrays of complex numbers.
  std::vector<std::array<std::complex<double>, 3>> result(theta_phi_pairs.size());
  StaticVector<3> tmp_r, tmp_i;
  for (std::size_t k = 0; k < theta_phi_pairs.size(); k++)
  {
    linalg::Cross3(r_naughts[k], integrals_r[k], tmp_r);
    linalg::Cross3(r_naughts[k], integrals_i[k], tmp_i);
    for (std::size_t d = 0; d < 3; d++)
    {
      result[k][d] = std::complex<double>{tmp_r[d], tmp_i[d]};
    }
  }
  return result;
}

}  // namespace palace
