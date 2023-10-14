// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "lumpedelement.hpp"

#include "fem/integrator.hpp"
#include "utils/communication.hpp"

namespace palace
{

double LumpedElementData::GetArea(mfem::ParFiniteElementSpace &fespace)
{
  mfem::ConstantCoefficient one_func(1.0);
  mfem::LinearForm s(&fespace);
  s.AddBoundaryIntegrator(new BoundaryLFIntegrator(one_func), attr_marker);
  s.UseFastAssembly(false);
  s.Assemble();

  mfem::GridFunction ones(&fespace);
  ones = 1.0;
  double dot = s * ones;
  Mpi::GlobalSum(1, &dot, fespace.GetComm());
  return dot;
}

UniformElementData::UniformElementData(const std::array<double, 3> &input_dir,
                                       const mfem::Array<int> &marker,
                                       mfem::ParFiniteElementSpace &fespace)
  : LumpedElementData(fespace.GetParMesh()->SpaceDimension(), marker),
    bounding_box(mesh::GetBoundingBox(*fespace.GetParMesh(), marker, true)), direction(3)
{
  // Check that the bounding box discovered matches the area. This validates that the
  // boundary elements form a right angled quadrilateral port.
  constexpr double rel_tol = 1.0e-6;
  double A = GetArea(fespace);
  MFEM_VERIFY((!bounding_box.planar || (std::abs(A - bounding_box.Area()) / A < rel_tol)),
              "Discovered bounding box area "
                  << bounding_box.Area() << " and integrated area " << A
                  << " do not match: Planar port geometry is not a quadrilateral!");

  // Check the user specified direction aligns with an axis direction.
  constexpr double angle_warning_deg = 0.1;
  constexpr double angle_error_deg = 1.0;
  auto lengths = bounding_box.Lengths();
  auto deviation_deg = bounding_box.Deviation(input_dir);
  if (std::none_of(deviation_deg.begin(), deviation_deg.end(),
                   [](double x) { return x < angle_warning_deg; }))
  {
    auto normal_0 = bounding_box.normals[0];
    for (auto &x : normal_0)
    {
      x /= lengths[0];
    }
    auto normal_1 = bounding_box.normals[1];
    for (auto &x : normal_1)
    {
      x /= lengths[1];
    }
    auto normal_2 = bounding_box.normals[2];
    for (auto &x : normal_2)
    {
      x /= lengths[2];
    }
    Mpi::Warning(
        "User specified direction {:.3e} does not align with either bounding box "
        "axis up to {:.3e} degrees!\n"
        "Axis 1: {:.3e} ({:.3e} degrees)\nAxis 2: {:.3e} ({:.3e} degrees)\nAxis 3: "
        "{:.3e} ({:.3e} degrees)!\n",
        input_dir, angle_warning_deg, normal_0, deviation_deg[0], normal_1,
        deviation_deg[1], normal_2, deviation_deg[2]);
  }
  MFEM_VERIFY(std::any_of(deviation_deg.begin(), deviation_deg.end(),
                          [](double x) { return x < angle_error_deg; }),
              "Specified direction does not align sufficiently with bounding box axes: "
                  << deviation_deg[0] << ' ' << deviation_deg[1] << ' ' << deviation_deg[2]
                  << " tolerance " << angle_error_deg << '!');
  std::copy(input_dir.begin(), input_dir.end(), direction.begin());
  direction /= direction.Norml2();

  // Compute the length from the most aligned normal direction.
  l = lengths[std::distance(deviation_deg.begin(),
                            std::min_element(deviation_deg.begin(), deviation_deg.end()))];
  MFEM_ASSERT(
      (l - mesh::GetProjectedLength(*fespace.GetParMesh(), marker, true, input_dir)) / l <
          rel_tol,
      "Bounding box discovered length should match projected length!");
  w = A / l;
}

std::unique_ptr<mfem::VectorCoefficient>
UniformElementData::GetModeCoefficient(double coef) const
{
  mfem::Vector source = direction;
  source *= coef;
  return std::make_unique<mfem::VectorConstantCoefficient>(source);
}

CoaxialElementData::CoaxialElementData(const std::array<double, 3> &direction,
                                       const mfem::Array<int> &marker,
                                       mfem::ParFiniteElementSpace &fespace)
  : LumpedElementData(fespace.GetParMesh()->SpaceDimension(), marker),
    bounding_ball(mesh::GetBoundingBall(*fespace.GetParMesh(), marker, true)),
    sign(direction[0] > 0)
{
  MFEM_VERIFY(bounding_ball.planar,
              "Boundary elements must be coplanar to define a coaxial lumped element!");

  // Get inner radius of annulus assuming full 2π circumference.
  double A = GetArea(fespace);
  MFEM_VERIFY(bounding_ball.radius > 0.0 &&
                  std::pow(bounding_ball.radius, 2) - A / M_PI > 0.0,
              "Coaxial element boundary is not defined correctly: Radius "
                  << bounding_ball.radius << ", area " << A << "!");
  ra = std::sqrt(std::pow(bounding_ball.radius, 2) - A / M_PI);
}

std::unique_ptr<mfem::VectorCoefficient>
CoaxialElementData::GetModeCoefficient(double coef) const
{
  double scoef = (sign ? 1.0 : -1.0) * coef;
  mfem::Vector x0(3);
  std::copy(bounding_ball.center.begin(), bounding_ball.center.end(), x0.begin());
  auto Source = [scoef, x0](const mfem::Vector &x, mfem::Vector &f) -> void
  {
    f = x;
    f -= x0;
    double oor = 1.0 / f.Norml2();
    f *= scoef * oor * oor;
  };
  return std::make_unique<mfem::VectorFunctionCoefficient>(dim, Source);
}

}  // namespace palace
