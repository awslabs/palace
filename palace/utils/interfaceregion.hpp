// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_INTERFACEREGION_HPP
#define PALACE_UTILS_INTERFACEREGION_HPP

#include <array>
#include <cmath>
#include <optional>
#include <vector>
#include <mfem.hpp>

namespace palace
{

// Opt-in quadrature-level spatial filter of one interface dielectric postprocessing entry
// (decision 352 follow-up (2)): a point belongs to the region when it lies in the optional
// axis-aligned box AND, when segments are given, in the translational cell of at least one
// segment. The translational cell of a segment AB is the references lane's: with the
// normal n projected out (in-plane geometry), the along-coordinate s = (P - A) . t lies in
// [0, |AB|] (flat cuts perpendicular to the segment at its ends) and the in-plane
// transverse distance to the line AB is at most `distance`. All intervals are closed; a
// quadrature point exactly on a cut is a measure-zero event. A region partition of one
// interface (adjacent cells, complementary boxes) therefore reproduces the unfiltered
// energy to roundoff, independent of the mesh faces.
class InterfaceRegion
{
private:
  struct Segment
  {
    std::array<double, 3> first;
    std::array<double, 3> direction;  // unit, in-plane
    double length;
  };
  std::optional<std::array<double, 3>> box_min, box_max;
  std::vector<Segment> segments;
  double distance = 0.0;
  std::array<double, 3> normal = {0.0, 0.0, 1.0};

  std::array<double, 3> ProjectInPlane(const std::array<double, 3> &v) const
  {
    const double along_normal = v[0] * normal[0] + v[1] * normal[1] + v[2] * normal[2];
    return {v[0] - along_normal * normal[0], v[1] - along_normal * normal[1],
            v[2] - along_normal * normal[2]};
  }

public:
  // Box bounds, segments [x0, y0, z0, x1, y1, z1], the transverse distance and the plane
  // normal in one consistent (nondimensional) length unit; a box without segments, segments
  // without a box, or both. Segments of zero in-plane length are refused.
  InterfaceRegion(const std::optional<std::array<double, 3>> &box_min_,
                  const std::optional<std::array<double, 3>> &box_max_,
                  const std::vector<std::array<double, 6>> &segments_, double distance_,
                  const std::array<double, 3> &normal_)
    : box_min(box_min_), box_max(box_max_), distance(distance_)
  {
    MFEM_VERIFY(box_min.has_value() == box_max.has_value(),
                "An interface region box requires both of its bounds!");
    MFEM_VERIFY(box_min || !segments_.empty(),
                "An interface region requires a box or segments!");
    MFEM_VERIFY(segments_.empty() || (std::isfinite(distance) && distance > 0.0),
                "Interface region segments require a positive transverse distance!");
    const double norm = std::sqrt(normal_[0] * normal_[0] + normal_[1] * normal_[1] +
                                  normal_[2] * normal_[2]);
    MFEM_VERIFY(std::isfinite(norm) && norm > 0.0,
                "An interface region normal must be nonzero!");
    for (int d = 0; d < 3; d++)
    {
      normal[d] = normal_[d] / norm;
      MFEM_VERIFY(!box_min ||
                      (std::isfinite((*box_min)[d]) && std::isfinite((*box_max)[d]) &&
                       (*box_min)[d] <= (*box_max)[d]),
                  "Interface region box bounds must be finite and ordered!");
    }
    segments.reserve(segments_.size());
    for (const auto &segment : segments_)
    {
      const std::array<double, 3> first = {segment[0], segment[1], segment[2]};
      const auto chord = ProjectInPlane(
          {segment[3] - segment[0], segment[4] - segment[1], segment[5] - segment[2]});
      const double length =
          std::sqrt(chord[0] * chord[0] + chord[1] * chord[1] + chord[2] * chord[2]);
      MFEM_VERIFY(std::isfinite(length) && length > 0.0,
                  "Interface region segments must have a finite, nonzero in-plane length!");
      segments.push_back(
          {first, {chord[0] / length, chord[1] / length, chord[2] / length}, length});
    }
  }

  bool HasBox() const { return box_min.has_value(); }
  std::size_t SegmentCount() const { return segments.size(); }

  bool Contains(const double *point) const
  {
    if (box_min)
    {
      for (int d = 0; d < 3; d++)
      {
        if (point[d] < (*box_min)[d] || point[d] > (*box_max)[d])
        {
          return false;
        }
      }
    }
    if (segments.empty())
    {
      return true;
    }
    for (const auto &segment : segments)
    {
      const auto relative =
          ProjectInPlane({point[0] - segment.first[0], point[1] - segment.first[1],
                          point[2] - segment.first[2]});
      const double along = relative[0] * segment.direction[0] +
                           relative[1] * segment.direction[1] +
                           relative[2] * segment.direction[2];
      if (along < 0.0 || along > segment.length)
      {
        continue;
      }
      double transverse_squared = 0.0;
      for (int d = 0; d < 3; d++)
      {
        const double transverse = relative[d] - along * segment.direction[d];
        transverse_squared += transverse * transverse;
      }
      if (transverse_squared <= distance * distance)
      {
        return true;
      }
    }
    return false;
  }

  // A two-dimensional point lies in the plane z = 0.
  bool Contains(const mfem::Vector &point) const
  {
    MFEM_ASSERT(point.Size() == 2 || point.Size() == 3,
                "Interface regions are two- or three-dimensional!");
    double coordinates[3] = {0.0, 0.0, 0.0};
    for (int d = 0; d < point.Size(); d++)
    {
      coordinates[d] = point(d);
    }
    return Contains(coordinates);
  }
};

}  // namespace palace

#endif  // PALACE_UTILS_INTERFACEREGION_HPP
