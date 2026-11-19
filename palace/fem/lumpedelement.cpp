// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "lumpedelement.hpp"

#include <array>
#include <cmath>
#include <numbers>
#include <string>
#include <Eigen/Dense>
#include "fem/coefficient.hpp"
#include "fem/integrator.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"

namespace palace
{

namespace
{

// Find the axis of a uniform lumped element from the covariance of its surface. For a flat
// element (planar in 3D, straight in 2D) with distinct principal axes in its plane, this
// is the principal axis closest to the direction: cuts normal to it then cross the element
// from side to side, even for a strip much wider than it is long at an angle to the
// direction. Otherwise, as for a folded or cylindrical element, whose axis of extrusion is
// the direction, or one close to square, whose principal axes are sensitive to small
// asymmetries, it is the direction itself. Returns whether the axis is a principal axis.
bool GetStripAxis(const mfem::DenseMatrix &covariance, const mfem::Vector &direction,
                  mfem::Vector &axis)
{
  constexpr double flat_tol = 1.0e-10, distinct_tol = 0.5;
  const int sdim = direction.Size();
  Eigen::MatrixXd C(sdim, sdim);
  for (int i = 0; i < sdim; i++)
  {
    for (int j = 0; j < sdim; j++)
    {
      C(i, j) = covariance(i, j);
    }
  }
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(C);
  const auto &lambda = eig.eigenvalues();  // Ascending
  axis = direction;
  if (!(lambda(sdim - 1) > 0.0) || lambda(0) > flat_tol * lambda(sdim - 1))
  {
    return false;
  }
  // The principal axes in the plane of the element are all but the first.
  int k = sdim - 1;
  if (sdim == 3)
  {
    if (lambda(2) - lambda(1) < distinct_tol * lambda(2))
    {
      return false;
    }
    double dot1 = 0.0, dot2 = 0.0;
    for (int i = 0; i < sdim; i++)
    {
      dot1 += eig.eigenvectors()(i, 1) * direction(i);
      dot2 += eig.eigenvectors()(i, 2) * direction(i);
    }
    k = (std::abs(dot1) > std::abs(dot2)) ? 1 : 2;
  }
  double dot = 0.0;
  for (int i = 0; i < sdim; i++)
  {
    dot += eig.eigenvectors()(i, k) * direction(i);
  }
  for (int i = 0; i < sdim; i++)
  {
    axis(i) = (dot < 0.0 ? -1.0 : 1.0) * eig.eigenvectors()(i, k);
  }
  return true;
}

// Angle in degrees between two vectors.
double AngleDeg(const mfem::Vector &u, const mfem::Vector &v)
{
  const double v_norm = v.Norml2(), u_along = (u * v) / v_norm;
  mfem::Vector u_across(u);
  u_across.Add(-u_along / v_norm, v);
  return std::atan2(u_across.Norml2(), u_along) * 180.0 / std::numbers::pi;
}

}  // namespace

UniformElementData::UniformElementData(const std::array<double, 3> &input_dir,
                                       const mfem::Array<int> &attr_list,
                                       const mfem::ParMesh &mesh)
  : LumpedElementData(attr_list)
{
  const int sdim = mesh.SpaceDimension();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> attr_marker = mesh::AttrToMarker(bdr_attr_max, attr_list);
  const std::string name =
      fmt::format("Uniform lumped element at attributes {}", fmt::join(attr_list, ", "));
  direction.SetSize(sdim);
  for (int i = 0; i < sdim; i++)
  {
    direction[i] = input_dir[i];
  }
  const double direction_norm = direction.Norml2();
  MFEM_VERIFY(
      direction_norm > 0.0,
      fmt::format("{} has a direction with no component in the {:d}D mesh!", name, sdim));
  direction /= direction_norm;

  // The element is a strip along its axis, between conductors at its two ends. Its width w
  // is the measure of its cross-section normal to its axis half way along it, and its
  // length l = area / w is the mean distance between its ends. The axis is the direction,
  // or for a flat element its principal axis closest to the direction, so that l and w
  // differ from those normal to the direction only at second order in the angle between
  // them. The strip may be folded or bent about its axis, have faces on both sides (which
  // doubles w), and have ends which are not straight. In 2D, the boundary is a 1D edge with
  // a unit out-of-plane depth: the cross-section counts the edges crossing the line, and
  // the "area" is the edge length.
  mfem::Vector center;
  mfem::DenseMatrix covariance;
  const double area = mesh::GetSurfaceMoments(mesh, attr_marker, center, covariance);
  MFEM_VERIFY(area > 0.0, name + " has zero area!");
  constexpr double rel_tol = 1.0e-6;
  const double size = (sdim == 3) ? std::sqrt(area) : area;
  const auto [d_min, d_max] = mesh::GetProjectedInterval(mesh, attr_marker, direction);
  MFEM_VERIFY(d_max - d_min > rel_tol * size,
              fmt::format("{} has no extent along its direction (extent {:.3e}, element "
                          "size {:.3e})!",
                          name, d_max - d_min, size));
  mfem::Vector axis;
  const bool principal = GetStripAxis(covariance, direction, axis);
  const auto [x_min, x_max] = principal
                                  ? mesh::GetProjectedInterval(mesh, attr_marker, axis)
                                  : std::pair<double, double>{d_min, d_max};
  const double extent = x_max - x_min;
  std::array<double, 3> section;
  std::array<mfem::Vector, 3> centroid;
  for (int i = 0; i < 3; i++)
  {
    section[i] = mesh::GetCrossSection(mesh, attr_marker, axis,
                                       x_min + 0.25 * (i + 1) * extent, centroid[i]);
    MFEM_VERIFY(section[i] > 0.0,
                fmt::format("{} has no cross-section normal to its axis at {:d}% of its "
                            "extent: the element is not a strip along its direction!",
                            name, 25 * (i + 1)));
  }
  w = section[1];
  l = area / w;

  // Check that the element is a strip along the direction:
  //  - The axis must be close to the direction.
  //  - The cross-sections at 1/4 and 3/4 of the extent must be centered on a line along the
  //    axis, which catches a direction at an angle to a strip which is not flat or close to
  //    square, and elements which are not strips, such as L-shaped ones.
  //  - The mean length must match the extent of the element along its axis (the principal
  //    axis, or the line through the centroids of the cross-sections), which catches area
  //    beyond the ends of the strip, as for a square with the direction along its diagonal.
  constexpr double angle_warning_deg = 0.1, angle_error_deg = 1.0;
  constexpr double length_warning = 1.0e-2, length_error = 5.0e-2;
  const char *axis_name = principal ? "its principal axis" : "its direction";
  mfem::Vector line(centroid[2]);
  line -= centroid[0];
  const double misalign_deg = AngleDeg(axis, direction);
  const double straight_deg = AngleDeg(line, axis);
  const auto [a_min, a_max] = principal
                                  ? std::pair<double, double>{x_min, x_max}
                                  : mesh::GetProjectedInterval(mesh, attr_marker, line);
  const double length_dev = std::abs(l / (a_max - a_min) - 1.0);
  MFEM_VERIFY(misalign_deg <= angle_error_deg,
              fmt::format("{} has its direction {:.3e} degrees from its principal axis "
                          "(tolerance {:g} degrees): the element is not a strip along its "
                          "direction!",
                          name, misalign_deg, angle_error_deg));
  MFEM_VERIFY(straight_deg <= angle_error_deg,
              fmt::format("{} has cross-sections centered on a line {:.3e} degrees from {} "
                          "(tolerance {:g} degrees): the element is not a strip along its "
                          "direction!",
                          name, straight_deg, axis_name, angle_error_deg));
  MFEM_VERIFY(
      length_dev <= length_error,
      fmt::format("{} has a mean length {:.3e} which differs from its extent {:.3e} "
                  "along its axis by {:.2f}% (tolerance {:g}%): the element is not "
                  "a strip along its direction!",
                  name, l, a_max - a_min, 100.0 * length_dev, 100.0 * length_error));
  if (misalign_deg > angle_warning_deg)
  {
    Mpi::Warning(mesh.GetComm(),
                 "{} has its direction {:.3e} degrees from its principal axis!\n", name,
                 misalign_deg);
  }
  if (straight_deg > angle_warning_deg)
  {
    Mpi::Warning(mesh.GetComm(),
                 "{} has cross-sections centered on a line {:.3e} degrees from {}!\n", name,
                 straight_deg, axis_name);
  }
  if (length_dev > length_warning)
  {
    Mpi::Warning(mesh.GetComm(),
                 "{} has a mean length {:.3e} which differs from its extent {:.3e} along "
                 "its axis by {:.2f}%: its ends are not straight and normal to its axis!\n",
                 name, l, a_max - a_min, 100.0 * length_dev);
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
