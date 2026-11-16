// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "lumpedelement.hpp"

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

  // The length is the extent of the element along the direction, which is zero for a
  // direction normal to the element (or, in 2D, to the boundary edge).
  constexpr double rel_tol = 1.0e-6;
  // Relative tolerance of the rectangle checks and of whether an element is planar.
  constexpr double shape_tol = 1.0e-3;
  mfem::Vector bbmin, bbmax;
  mesh::GetAxisAlignedBoundingBox(mesh, attr_marker, true, bbmin, bbmax);
  bbmax -= bbmin;
  const double size = bbmax.Normlinf();
  l = mesh::GetProjectedLength(mesh, attr_marker, true, direction);
  MFEM_VERIFY(l > rel_tol * size,
              "Uniform lumped element has no extent along its direction (length "
                  << l << ", element size " << size << ")!");

  // Compute the width as area / length in 3D. This allows the lumped element to be
  // non-planar, for example folded or bent about the direction. In 2D, the boundary is a 1D
  // edge with a unit out-of-plane depth, and the "area" is its length.
  const double area = mesh::GetSurfaceArea(mesh, attr_marker);
  MFEM_VERIFY(area > 0.0, "Uniform lumped element has zero area!");
  w = (sdim == 3) ? area / l : 1.0;

  // An element swept along the direction has an area of the length times its unfolded
  // width, which is at least its extent across the direction. A planar element must match
  // it: it is a rectangle with sides along and across the direction (in 2D, a straight edge
  // along the direction). The direction of an element which is not planar is assumed to
  // lie in its surface.
  mfem::Vector normal = mesh::GetSurfaceNormal(mesh, attr_marker);
  double w_across = 1.0;
  if (sdim == 3)
  {
    mfem::Vector across(sdim);
    linalg::Cross3(normal, direction, across);
    w_across = mesh::GetProjectedLength(mesh, attr_marker, true, across);
  }
  MFEM_VERIFY(area >= (1.0 - shape_tol) * l * w_across,
              "Uniform lumped element area ("
                  << area << ") is less than its length times width (" << l * w_across
                  << "): port geometry is not rectangular!");
  const bool planar =
      mesh::GetProjectedLength(mesh, attr_marker, true, normal) <= shape_tol * size;
  MFEM_VERIFY(!planar || area <= (1.0 + shape_tol) * l * w_across,
              "Planar uniform lumped element area ("
                  << area << ") exceeds its length times width (" << l * w_across
                  << "): port geometry is not rectangular!");
  if (!planar)
  {
    Mpi::Warning(mesh.GetComm(),
                 "Uniform lumped element at attributes {} is not planar; its direction is "
                 "not checked!\n",
                 fmt::join(attr_list, ", "));
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
