// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_RESPONSE_SEGMENT_INDEX_HPP
#define PALACE_MODELS_SURFACE_RESPONSE_SEGMENT_INDEX_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>
#include <mfem.hpp>

namespace palace
{

//
// Deterministic broad-phase index for nearby segment queries (an AABB tree over straight
// segments, split at the median along the widest axis, leaves of 16). Exact distance and
// topology predicates remain at the call sites; this index only removes pairs whose
// expanded AABBs cannot interact. Shared by the operator's segment proximity passes
// (surfaceresponseoperator.cpp) and the nearest-sub-segment query below.
//
class SegmentBoxIndex
{
private:
  using Point3D = std::array<double, 3>;
  static constexpr std::size_t leaf_size = 16;
  struct Box
  {
    Point3D min{};
    Point3D max{};
  };
  struct Node
  {
    Box box;
    std::size_t begin = 0;
    std::size_t end = 0;
    int left = -1;
    int right = -1;
    bool IsLeaf() const { return left < 0; }
  };

  std::vector<Box> boxes;
  std::vector<std::size_t> indices;
  std::vector<Node> nodes;

  static bool Intersects(const Box &a, const Box &b)
  {
    for (int d = 0; d < 3; d++)
    {
      if (a.max[d] < b.min[d] || b.max[d] < a.min[d])
      {
        return false;
      }
    }
    return true;
  }

  int Build(std::size_t begin, std::size_t end)
  {
    Node node;
    node.begin = begin;
    node.end = end;
    node.box.min.fill(mfem::infinity());
    node.box.max.fill(-mfem::infinity());
    for (std::size_t i = begin; i < end; i++)
    {
      const auto &box = boxes[indices[i]];
      for (int d = 0; d < 3; d++)
      {
        node.box.min[d] = std::min(node.box.min[d], box.min[d]);
        node.box.max[d] = std::max(node.box.max[d], box.max[d]);
      }
    }
    const int node_index = static_cast<int>(nodes.size());
    nodes.push_back(node);
    if (end - begin <= leaf_size)
    {
      return node_index;
    }
    int axis = 0;
    for (int d = 1; d < 3; d++)
    {
      if (node.box.max[d] - node.box.min[d] > node.box.max[axis] - node.box.min[axis])
      {
        axis = d;
      }
    }
    const std::size_t middle = begin + (end - begin) / 2;
    std::nth_element(indices.begin() + begin, indices.begin() + middle,
                     indices.begin() + end,
                     [&](std::size_t first, std::size_t second)
                     {
                       return boxes[first].min[axis] + boxes[first].max[axis] <
                                  boxes[second].min[axis] + boxes[second].max[axis] ||
                              (boxes[first].min[axis] + boxes[first].max[axis] ==
                                   boxes[second].min[axis] + boxes[second].max[axis] &&
                               first < second);
                     });
    nodes[node_index].left = Build(begin, middle);
    nodes[node_index].right = Build(middle, end);
    return node_index;
  }

  void Query(int node_index, const Box &query, std::vector<std::size_t> &matches) const
  {
    const auto &node = nodes[node_index];
    if (!Intersects(node.box, query))
    {
      return;
    }
    if (node.IsLeaf())
    {
      for (std::size_t i = node.begin; i < node.end; i++)
      {
        const std::size_t index = indices[i];
        if (Intersects(boxes[index], query))
        {
          matches.push_back(index);
        }
      }
      return;
    }
    Query(node.left, query, matches);
    Query(node.right, query, matches);
  }

public:
  explicit SegmentBoxIndex(const std::vector<std::pair<Point3D, Point3D>> &segments)
    : boxes(segments.size()), indices(segments.size())
  {
    MFEM_VERIFY(!segments.empty(), "Cannot build an empty segment proximity index!");
    for (std::size_t i = 0; i < segments.size(); i++)
    {
      for (int d = 0; d < 3; d++)
      {
        boxes[i].min[d] = std::min(segments[i].first[d], segments[i].second[d]);
        boxes[i].max[d] = std::max(segments[i].first[d], segments[i].second[d]);
      }
    }
    std::iota(indices.begin(), indices.end(), 0);
    nodes.reserve(2 * (segments.size() / leaf_size + 1));
    Build(0, segments.size());
  }

  std::size_t Size() const { return boxes.size(); }

  // The root box (every segment's box) and the largest coordinate magnitude in it.
  const Point3D &Min() const { return nodes.front().box.min; }
  const Point3D &Max() const { return nodes.front().box.max; }

  // The segments whose boxes meet the box of the segment (p0, p1) widened by `distance`
  // along every axis, in increasing index order.
  std::vector<std::size_t> Query(const Point3D &p0, const Point3D &p1,
                                 double distance) const
  {
    Box query;
    for (int d = 0; d < 3; d++)
    {
      query.min[d] = std::min(p0[d], p1[d]) - distance;
      query.max[d] = std::max(p0[d], p1[d]) + distance;
    }
    std::vector<std::size_t> matches;
    Query(0, query, matches);
    std::sort(matches.begin(), matches.end());
    return matches;
  }

  std::vector<std::pair<std::size_t, std::size_t>> CandidatePairs(double distance) const
  {
    std::vector<std::pair<std::size_t, std::size_t>> result;
    for (std::size_t first = 0; first < boxes.size(); first++)
    {
      Point3D p0 = boxes[first].min;
      Point3D p1 = boxes[first].max;
      for (const std::size_t second : Query(p0, p1, distance))
      {
        if (second > first)
        {
          result.emplace_back(first, second);
        }
      }
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
  }
};

//
// The nearest sub-segment of a chain to a point, with the broad phase on a SegmentBoxIndex
// and the SAME verdict as a linear scan of the chain (OI-PLACEMENT-PASS-SCALING, decision
// 1112): the placement's ClosestFoot (surfaceresponseoperator.cpp BuildFeaturePatches)
// reads, for every placed translational patch, the foot of the patch origin and of the two
// end-face points on the feature's side chains (the R0-a conductor probes of the
// identification fix round, DESIGN A.1 (c), and the pair / stack placers' feet); scanning a
// chain of S sub-segments once per patch costs O(S) per patch and O(S^2) per feature, which
// is what made the patch build of a chip of long chorded resonators (QUO-023: 381,527
// patches) take 677 s.
//
// THE CONTRACT. `sub_segments` are the chain's straight pieces by their end points (the
// index i of a piece = its position in the chain). `Nearest(p, distance_of)` returns the
// index the linear scan `for i in 0..S-1: if (distance_of(i) < best) best = distance_of(i),
// winner = i` would return, i.e. the FIRST index attaining the minimum of the computed
// distances, provided `distance_of(i)` is the Euclidean distance from p to piece i up to
// floating-point rounding (the caller's own expression; the index never computes a distance
// itself). It never changes a result: the broad phase only excludes pieces whose box lies
// farther than the best computed distance plus a margin that exceeds every rounding term
// (kNearestMarginRelative x the coordinate scale), and the Chebyshev distance from p to a
// piece's box never exceeds the Euclidean distance to the piece, so every piece whose
// computed distance equals the minimum is among the candidates; the candidates are
// evaluated in increasing index order with the scan's strict `<`. The query: a box query at
// a doubling radius (from the mean piece extent) until a candidate is found; its
// candidates' best computed distance d1 bounds the minimum; when d1 + margin exceeds the
// query radius the candidates of radius d1 + margin are evaluated instead (a superset of
// the first set). `evaluated` (optional) counts the distance evaluations (the unit test's
// scaling reading).
//
class NearestSubSegmentIndex
{
public:
  using Point3D = std::array<double, 3>;
  static constexpr double kNearestMarginRelative = 1.0e-9;

  explicit NearestSubSegmentIndex(
      const std::vector<std::pair<Point3D, Point3D>> &sub_segments)
    : tree(sub_segments)
  {
    double extent_sum = 0.0;
    for (const auto &[a, b] : sub_segments)
    {
      double extent = 0.0;
      for (int d = 0; d < 3; d++)
      {
        extent = std::max(extent, std::abs(b[d] - a[d]));
        scale = std::max({scale, std::abs(a[d]), std::abs(b[d])});
      }
      extent_sum += extent;
    }
    mean_extent = extent_sum / static_cast<double>(sub_segments.size());
    double diagonal = 0.0;
    for (int d = 0; d < 3; d++)
    {
      diagonal = std::max(diagonal, tree.Max()[d] - tree.Min()[d]);
    }
    reach = diagonal;
  }

  std::size_t Size() const { return tree.Size(); }

  // The result is std::numeric_limits<std::size_t>::max() when no piece has a computed
  // distance below +infinity (the scan's "no winner").
  template <typename DistanceOf>
  std::size_t Nearest(const Point3D &p, DistanceOf &&distance_of,
                      std::size_t *evaluated = nullptr) const
  {
    constexpr std::size_t none = std::numeric_limits<std::size_t>::max();
    double point_scale = scale;
    for (int d = 0; d < 3; d++)
    {
      point_scale = std::max(point_scale, std::abs(p[d]));
    }
    const double margin = kNearestMarginRelative * point_scale;
    auto Evaluate = [&](const std::vector<std::size_t> &candidates, double &best_distance)
    {
      std::size_t winner = none;
      best_distance = mfem::infinity();
      for (const std::size_t index : candidates)
      {
        const double distance = distance_of(index);
        if (evaluated)
        {
          (*evaluated)++;
        }
        if (distance < best_distance)
        {
          best_distance = distance;
          winner = index;
        }
      }
      return winner;
    };
    // Every piece lies within `reach` + the point's distance to the root box; a radius
    // beyond that returns the whole chain (the scan itself).
    double root_distance = 0.0;
    for (int d = 0; d < 3; d++)
    {
      root_distance = std::max({root_distance, tree.Min()[d] - p[d], p[d] - tree.Max()[d]});
    }
    const double everything = reach + root_distance + margin;
    double radius = std::max(mean_extent, margin);
    if (!(radius > 0.0))
    {
      radius = everything > 0.0 ? everything : 1.0;
    }
    for (;;)
    {
      const auto candidates = tree.Query(p, p, radius);
      if (!candidates.empty())
      {
        double best_distance;
        const std::size_t winner = Evaluate(candidates, best_distance);
        if (winner == none)
        {
          // No finite distance among the candidates: widen to the whole chain once.
          if (radius >= everything)
          {
            return none;
          }
          radius = everything;
          continue;
        }
        const double bound = best_distance + margin;
        if (bound <= radius || radius >= everything)
        {
          return winner;
        }
        return Evaluate(tree.Query(p, p, std::min(bound, everything)), best_distance);
      }
      if (radius >= everything)
      {
        return none;
      }
      radius = std::min(2.0 * radius, everything);
    }
  }

private:
  SegmentBoxIndex tree;
  double scale = 0.0;
  double mean_extent = 0.0;
  double reach = 0.0;
};

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_RESPONSE_SEGMENT_INDEX_HPP
