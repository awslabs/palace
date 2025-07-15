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
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> attr_marker = mesh::AttrToMarker(bdr_attr_max, attr_list);
  auto bounding_box = mesh::GetBoundingBox(mesh, attr_marker, true);

  // Check the user specified direction aligns with an axis direction.
  constexpr double angle_warning_deg = 0.1;
  constexpr double angle_error_deg = 1.0;
  auto lengths = bounding_box.Lengths();
  auto deviations_deg = bounding_box.Deviations(input_dir);
  if (std::none_of(deviations_deg.begin(), deviations_deg.end(),
                   [](double x) { return x < angle_warning_deg; }))
  {
    auto normals = bounding_box.Normals();
    Mpi::Warning("User specified direction {} does not align with either bounding box "
                 "axis up to {:.3e} degrees!\n"
                 "Axis 1: {} ({:.3e} degrees)\nAxis 2: {} ({:.3e} degrees)\nAxis 3: "
                 "{} ({:.3e} degrees)!\n",
                 input_dir, angle_warning_deg, normals[0], deviations_deg[0], normals[1],
                 deviations_deg[1], normals[2], deviations_deg[2]);
  }
  if (std::none_of(deviations_deg.begin(), deviations_deg.end(),
                   [](double x) { return x < angle_error_deg; }))
  {
    Mpi::Barrier(mesh.GetComm());
    MFEM_ABORT("Specified direction does not align sufficiently with bounding box axes ("
               << deviations_deg[0] << ", " << deviations_deg[1] << ", "
               << deviations_deg[2] << " vs. tolerance " << angle_error_deg << ")!");
  }
  direction.SetSize(input_dir.size());
  std::copy(input_dir.begin(), input_dir.end(), direction.begin());
  direction /= direction.Norml2();

  // Compute the length from the most aligned normal direction.
  constexpr double rel_tol = 1.0e-6;
  auto l_component =
      std::distance(deviations_deg.begin(),
                    std::min_element(deviations_deg.begin(), deviations_deg.end()));
  l = lengths[l_component];
  MFEM_VERIFY(std::abs(l - mesh::GetProjectedLength(mesh, attr_marker, true, input_dir)) <
                  rel_tol * l,
              "Bounding box discovered length ("
                  << l << ") should match projected length ("
                  << mesh::GetProjectedLength(mesh, attr_marker, true, input_dir) << "!");

  // Compute the width as area / length. This allows the lumped element to be non-planar,
  // and generalizes nicely to the case for an infinitely thin rectangular lumped element
  // with elements on both sides (for which the width computed from the bounding box would
  // be incorrect by a factor of 2).
  double area = mesh::GetSurfaceArea(mesh, attr_marker);
  MFEM_VERIFY(area > 0.0, "Uniform lumped element has zero area!");
  w = area / l;
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
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> attr_marker = mesh::AttrToMarker(bdr_attr_max, attr_list);
  auto bounding_ball = mesh::GetBoundingBall(mesh, attr_marker, true);
  MFEM_VERIFY(bounding_ball.planar,
              "Boundary elements must be coplanar to define a coaxial lumped element!");

  // Direction of the excitation as +/-r̂.
  direction = (input_dir[0] > 0 ? +1 : -1);
  origin.SetSize(bounding_ball.center.size());
  std::copy(bounding_ball.center.begin(), bounding_ball.center.end(), origin.begin());

  // Get outer and inner radius of the annulus, assuming full 2π circumference.
  r_outer = 0.5 * bounding_ball.Lengths()[0];
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
