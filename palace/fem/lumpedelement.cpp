// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "lumpedelement.hpp"

#include <string>
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
  auto bounding_box = mesh::GetBoundingBox(mesh, attr_marker, true);

  // Convert the input direction to mfem::Vector for use with BoundingBox methods.
  const int dim = bounding_box.Dim();
  mfem::Vector dir_vec(dim);
  for (int i = 0; i < dim; i++)
  {
    dir_vec(i) = input_dir[i];
  }

  // Check the user specified direction aligns with an axis direction.
  constexpr double angle_warning_deg = 0.1;
  constexpr double angle_error_deg = 1.0;
  auto lengths = bounding_box.Lengths();
  auto deviations_deg = bounding_box.Deviations(dir_vec);
  if (std::none_of(deviations_deg.begin(), deviations_deg.end(),
                   [](double x) { return x < angle_warning_deg; }))
  {
    auto normals = bounding_box.Normals();
    if (dim == 3)
    {
      mfem::Vector n0(normals.GetColumn(0), 3);
      mfem::Vector n1(normals.GetColumn(1), 3);
      mfem::Vector n2(normals.GetColumn(2), 3);
      Mpi::Warning("User specified direction {} does not align with either bounding "
                   "box axis up to {:.3e} degrees!\n"
                   "Axis 1: ({}) ({:.3e} degrees)\nAxis 2: ({}) ({:.3e} degrees)\n"
                   "Axis 3: ({}) ({:.3e} degrees)!\n",
                   input_dir, angle_warning_deg, fmt::join(n0, ", "), deviations_deg(0),
                   fmt::join(n1, ", "), deviations_deg(1), fmt::join(n2, ", "),
                   deviations_deg(2));
    }
    else
    {
      mfem::Vector n0(normals.GetColumn(0), 2);
      mfem::Vector n1(normals.GetColumn(1), 2);
      Mpi::Warning("User specified direction does not align with either bounding box "
                   "axis up to {:.3e} degrees!\n"
                   "Axis 1: ({}) ({:.3e} degrees)\n"
                   "Axis 2: ({}) ({:.3e} degrees)!\n",
                   angle_warning_deg, fmt::join(n0, ", "), deviations_deg(0),
                   fmt::join(n1, ", "), deviations_deg(1));
    }
  }
  if (std::none_of(deviations_deg.begin(), deviations_deg.end(),
                   [](double x) { return x < angle_error_deg; }))
  {
    Mpi::Barrier(mesh.GetComm());
    std::string dev_str;
    for (int i = 0; i < dim; i++)
    {
      if (i > 0)
      {
        dev_str += ", ";
      }
      dev_str += std::to_string(deviations_deg(i));
    }
    MFEM_ABORT("Specified direction does not align sufficiently with bounding box axes ("
               << dev_str << " vs. tolerance " << angle_error_deg << ")!");
  }
  direction.SetSize(sdim);
  for (int i = 0; i < sdim; i++)
  {
    direction[i] = input_dir[i];
  }
  direction /= direction.Norml2();

  // Compute the length from the most aligned normal direction.
  constexpr double rel_tol = 1.0e-6;
  int l_component = 0;
  for (int i = 1; i < dim; i++)
  {
    if (deviations_deg(i) < deviations_deg(l_component))
    {
      l_component = i;
    }
  }
  l = lengths(l_component);
  // Cross-check the bounding box length against the projected length. In 2D, the
  // direction may be perpendicular to a collinear boundary (e.g., +Y for a horizontal
  // edge), in which case the projected length is near zero and the check does not apply.
  {
    double proj_l = mesh::GetProjectedLength(mesh, attr_marker, true, dir_vec);
    MFEM_VERIFY(proj_l < rel_tol * l || std::abs(l - proj_l) < rel_tol * l,
                "Bounding box discovered length ("
                    << l << ") should match projected length (" << proj_l << "!");
  }

  if (sdim == 3)
  {
    // Compute the width as area / length. This allows the lumped element to be non-planar,
    // and generalizes nicely to the case for an infinitely thin rectangular lumped element
    // with elements on both sides (for which the width computed from the bounding box
    // would be incorrect by a factor of 2).
    double area = mesh::GetSurfaceArea(mesh, attr_marker);
    MFEM_VERIFY(area > 0.0, "Uniform lumped element has zero area!");
    w = area / l;
  }
  else
  {
    // In 2D, the boundary is a 1D edge. Use a unit out-of-plane depth convention so that
    // impedance and admittance calculations reduce to per-unit-length quantities.
    w = 1.0;
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
