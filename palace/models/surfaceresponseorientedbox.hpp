// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_RESPONSE_ORIENTED_BOX_HPP
#define PALACE_MODELS_SURFACE_RESPONSE_ORIENTED_BOX_HPP

#include <algorithm>
#include <array>
#include <cmath>

namespace palace
{

//
// The oriented (frame-aligned) box of a surface-response coupon domain and the
// separating-axis overlap test on two of them — ONE definition shared by the placement
// (surfaceresponseoperator: the same-class overlap test of hardening-2 B2 / H4a) and the
// identification (surfaceresponseidentification: the vertex-vertex interaction rule of the
// identification fix round, decision 1058 MAJOR-3 / 1062 A.5 (a)): two vertex sites are
// clustered iff THIS test says their +-R coupon cubes overlap, so the identification's
// clustering decision and the placement's overlap verdict are literally one criterion.
//
// The box: `axes` = the frame (rows = unit axes), `min` / `max` = the box's extent along
// each axis in that frame's coordinates (the frame-projected bounding box of a support's
// basis points; for a vertex coupon the +-R cube about the site: min = p . axis - R,
// max = p . axis + R). For mesh-aligned frames the test is the mesh-aligned
// InteriorOverlaps. The separating-axis test: the boxes are disjoint iff some axis among
// the face normals of both (and, in 3D, the 9 cross products of their frame axes)
// separates the projected centres by at least the sum of the projected half extents,
// within the tolerance 1e-12 x the coordinate scale (the ElementBox tolerance of record).
// Two R-cubes with parallel frames overlap up to 2 R apart along an axis and 2 sqrt(2) R
// along a face diagonal; frames 45 deg apart overlap up to (1 + sqrt 2) R = 2.414 R along
// a frame axis and 2.613 R along the 22.5-deg bisector (decision 1062 C.2); no pair of
// R-cubes overlaps beyond 2 sqrt(3) R (the probe / query domain).
//
struct OrientedSupportBox
{
  std::array<std::array<double, 3>, 3> axes = {
      {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
  std::array<double, 3> min{}, max{};
};

// The +-R cube of a vertex coupon about `origin` in the frame `axes` (the object the
// placement builds for a corner coupon's SpatialSupportBounds: the library models' basis
// points span exactly [-R, R] on every axis of the vertex frame).
inline OrientedSupportBox VertexCubeBox(const std::array<double, 3> &origin,
                                        const std::array<std::array<double, 3>, 3> &axes,
                                        double R)
{
  OrientedSupportBox box;
  box.axes = axes;
  for (int k = 0; k < 3; k++)
  {
    const double c =
        origin[0] * axes[k][0] + origin[1] * axes[k][1] + origin[2] * axes[k][2];
    box.min[k] = c - R;
    box.max[k] = c + R;
  }
  return box;
}

inline bool OrientedBoxesInteriorOverlap(const OrientedSupportBox &a,
                                         const OrientedSupportBox &b, int dimension)
{
  using Vec = std::array<double, 3>;
  auto Dot = [](const Vec &u, const Vec &v)
  { return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]; };
  auto Cross = [](const Vec &u, const Vec &v)
  {
    return Vec{u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
               u[0] * v[1] - u[1] * v[0]};
  };
  double scale = 1.0;
  for (int k = 0; k < dimension; k++)
  {
    scale = std::max({scale, std::abs(a.min[k]), std::abs(a.max[k]), std::abs(b.min[k]),
                      std::abs(b.max[k])});
  }
  const double tolerance = 1.0e-12 * scale;
  Vec centre_a{}, centre_b{};
  std::array<double, 3> half_a{}, half_b{};
  for (int k = 0; k < dimension; k++)
  {
    half_a[k] = 0.5 * (a.max[k] - a.min[k]);
    half_b[k] = 0.5 * (b.max[k] - b.min[k]);
    for (int d = 0; d < 3; d++)
    {
      centre_a[d] += 0.5 * (a.min[k] + a.max[k]) * a.axes[k][d];
      centre_b[d] += 0.5 * (b.min[k] + b.max[k]) * b.axes[k][d];
    }
  }
  const Vec d = {centre_a[0] - centre_b[0], centre_a[1] - centre_b[1],
                 centre_a[2] - centre_b[2]};
  auto Separates = [&](Vec axis)
  {
    const double length = std::sqrt(Dot(axis, axis));
    if (length <= 1.0e-12)
    {
      return false;  // parallel frame axes: their cross product is no direction
    }
    for (double &c : axis)
    {
      c /= length;
    }
    double reach = 0.0;
    for (int k = 0; k < dimension; k++)
    {
      reach += half_a[k] * std::abs(Dot(a.axes[k], axis)) +
               half_b[k] * std::abs(Dot(b.axes[k], axis));
    }
    return std::abs(Dot(d, axis)) >= reach - tolerance;
  };
  for (int k = 0; k < dimension; k++)
  {
    if (Separates(a.axes[k]) || Separates(b.axes[k]))
    {
      return false;
    }
  }
  if (dimension == 3)
  {
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        if (Separates(Cross(a.axes[i], b.axes[j])))
        {
          return false;
        }
      }
    }
  }
  return true;
}

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_RESPONSE_ORIENTED_BOX_HPP
