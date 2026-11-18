// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "lumpedelement.hpp"

#include <array>
#include <cmath>
#include <numbers>
#include "fem/coefficient.hpp"
#include "fem/integrator.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"

namespace palace
{

UniformElementData::UniformElementData(const std::array<double, 3> &input_dir,
                                       const mfem::Array<int> &attr_list,
                                       const mfem::ParMesh &mesh)
  : LumpedElementData(attr_list)
{
  const int sdim = mesh.SpaceDimension();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> attr_marker = mesh::AttrToMarker(bdr_attr_max, attr_list);
  direction.SetSize(sdim);
  for (int i = 0; i < sdim; i++)
  {
    direction[i] = input_dir[i];
  }
  direction /= direction.Norml2();

  // The element is a strip along the direction, between conductors at its two ends. Its
  // width w is the measure of its cross-section normal to the direction half way along it,
  // and its length l = area / w is the mean distance between its ends. The strip may be
  // folded or bent about the direction, have faces on both sides (which doubles w), and
  // have ends which are not straight. In 2D, the boundary is a 1D edge with a unit
  // out-of-plane depth: the cross-section counts the edges crossing the line, and the
  // "area" is the edge length.
  const double area = mesh::GetSurfaceArea(mesh, attr_marker);
  MFEM_VERIFY(area > 0.0, "Uniform lumped element has zero area!");
  constexpr double rel_tol = 1.0e-6;
  mfem::Vector bbmin, bbmax;
  mesh::GetAxisAlignedBoundingBox(mesh, attr_marker, true, bbmin, bbmax);
  bbmax -= bbmin;
  const double size = bbmax.Normlinf();
  const auto [x_min, x_max] =
      mesh::GetProjectedInterval(mesh, attr_marker, true, direction);
  const double extent = x_max - x_min;
  MFEM_VERIFY(extent > rel_tol * size,
              fmt::format("Uniform lumped element has no extent along its direction "
                          "(extent {:.3e}, element size {:.3e})!",
                          extent, size));
  std::array<double, 3> section;
  std::array<mfem::Vector, 3> centroid;
  for (int i = 0; i < 3; i++)
  {
    section[i] = mesh::GetCrossSection(mesh, attr_marker, direction,
                                       x_min + 0.25 * (i + 1) * extent, centroid[i]);
    MFEM_VERIFY(section[i] > 0.0,
                fmt::format("Uniform lumped element has no cross-section normal to its "
                            "direction at {:d}% of its extent: the element is not a strip "
                            "along its direction!",
                            25 * (i + 1)));
  }
  w = section[1];
  l = area / w;

  // Check that the element is a strip along the direction. Its axis, the line through the
  // centroids of the cross-sections at 1/4 and 3/4 of the extent, must run along the
  // direction, which catches a direction at an angle to the strip. The mean length must
  // match the extent of the element along its axis, which catches area beyond the ends of
  // the strip, as for an L-shaped element or a square with the direction along its
  // diagonal. Measured along the axis, not the direction, the extent of a strip at a small
  // angle to the direction matches its mean length to second order in the angle, even when
  // the strip is wider than it is long.
  constexpr double angle_warning_deg = 0.1, angle_error_deg = 1.0;
  constexpr double length_warning = 1.0e-2, length_error = 5.0e-2;
  mfem::Vector axis(centroid[2]);
  axis -= centroid[0];
  const double along = axis * direction;
  mfem::Vector across(axis);
  across.Add(-along, direction);
  const double angle_deg = std::atan2(across.Norml2(), along) * 180.0 / std::numbers::pi;
  MFEM_VERIFY(
      angle_deg <= angle_error_deg,
      fmt::format("Uniform lumped element direction is {:.3e} degrees from the axis "
                  "of the element (tolerance {:g} degrees): the element is not a "
                  "strip along its direction!",
                  angle_deg, angle_error_deg));
  const auto [a_min, a_max] = mesh::GetProjectedInterval(mesh, attr_marker, true, axis);
  const double axis_extent = a_max - a_min;
  const double length_dev = std::abs(l / axis_extent - 1.0);
  MFEM_VERIFY(length_dev <= length_error,
              fmt::format("Uniform lumped element mean length {:.3e} differs from its "
                          "extent {:.3e} along its axis by {:.2f}% (tolerance {:g}%): the "
                          "element is not a strip along its direction!",
                          l, axis_extent, 100.0 * length_dev, 100.0 * length_error));
  if (angle_deg > angle_warning_deg)
  {
    Mpi::Warning(mesh.GetComm(),
                 "Uniform lumped element at attributes {} has its direction {:.3e} degrees "
                 "from the axis of the element!\n",
                 fmt::join(attr_list, ", "), angle_deg);
  }
  if (length_dev > length_warning)
  {
    Mpi::Warning(mesh.GetComm(),
                 "Uniform lumped element at attributes {} has a mean length {:.3e} which "
                 "differs from its extent {:.3e} along its axis by {:.2f}%: its ends are "
                 "not straight and normal to its axis!\n",
                 fmt::join(attr_list, ", "), l, axis_extent, 100.0 * length_dev);
  }
}

std::unique_ptr<mfem::VectorCoefficient>
UniformElementData::GetModeCoefficient(double coeff) const
{
  mfem::Vector source = direction;
  source *= coeff;
  return std::make_unique<RestrictedVectorCoefficient<mfem::VectorConstantCoefficient>>(
      attr_list, source);
}

CoaxialElementData::CoaxialElementData(const std::array<double, 3> &input_dir,
                                       const mfem::Array<int> &attr_list,
                                       const mfem::ParMesh &mesh)
  : LumpedElementData(attr_list)
{
  MFEM_VERIFY(mesh.SpaceDimension() == 3,
              "Coaxial lumped elements are only supported for 3D simulations!");
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> attr_marker = mesh::AttrToMarker(bdr_attr_max, attr_list);
  auto bounding_ball = mesh::GetBoundingBall(mesh, attr_marker, true);
  MFEM_VERIFY(bounding_ball.planar,
              "Boundary elements must be coplanar to define a coaxial lumped element!");

  // Direction of the excitation as +/-r̂.
  direction = (input_dir[0] > 0 ? +1 : -1);
  origin = bounding_ball.center;

  // Get outer and inner radius of the annulus, assuming full 2π circumference.
  r_outer = 0.5 * bounding_ball.Lengths()(0);
  r_inner = mesh::GetDistanceFromPoint(mesh, attr_marker, true, bounding_ball.center);
  MFEM_VERIFY(r_inner > 0.0,
              "Coaxial element annulus has should have positive inner radius!");
  MFEM_VERIFY(r_outer > r_inner, "Coaxial element annulus has unexpected outer radius "
                                     << r_outer << " <= inner radius " << r_inner << "!");
}

std::unique_ptr<mfem::VectorCoefficient>
CoaxialElementData::GetModeCoefficient(double coeff) const
{
  coeff *= direction;
  mfem::Vector x0(origin);
  auto Source = [coeff, x0](const mfem::Vector &x, mfem::Vector &f) -> void
  {
    f = x;
    f -= x0;
    double oor = 1.0 / f.Norml2();
    f *= coeff * oor * oor;
  };
  return std::make_unique<RestrictedVectorCoefficient<mfem::VectorFunctionCoefficient>>(
      attr_list, x0.Size(), Source);
}

}  // namespace palace
