// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_LUMPED_ELEMENT_HPP
#define PALACE_FEM_LUMPED_ELEMENT_HPP

#include <memory>
#include <string>
#include <mfem.hpp>
#include "fem/integrator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"

namespace palace
{

//
// Base class handling geometry of lumped elements for uniform and coaxial lumped port and
// surface current source boundaries.
//
class LumpedElementData
{
protected:
  // Spatial dimension.
  const int dim;

  // Marker for all boundary attributes making up this lumped element boundary.
  mfem::Array<int> attr_marker;

  double GetArea(mfem::ParFiniteElementSpace &fespace)
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

public:
  LumpedElementData(int d, const mfem::Array<int> &marker) : dim(d), attr_marker(marker) {}
  virtual ~LumpedElementData() = default;

  mfem::Array<int> &GetMarker() { return attr_marker; }
  const mfem::Array<int> &GetMarker() const { return attr_marker; }

  virtual double GetGeometryLength() const = 0;
  virtual double GetGeometryWidth() const = 0;

  virtual std::unique_ptr<mfem::VectorCoefficient>
  GetModeCoefficient(double coef = 1.0) const = 0;
};

class UniformElementData : public LumpedElementData
{
protected:
  mesh::BoundingBox bounding_box;  // Bounding box defining the rectangular lumped port
  mfem::Vector direction;  // Cartesian vector specifying signed direction of incident field
  double l, w;             // Lumped element length and width
public:
  UniformElementData(const std::array<double, 3> &input_dir, const mfem::Array<int> &marker,
                     mfem::ParFiniteElementSpace &fespace)
    : LumpedElementData(fespace.GetParMesh()->SpaceDimension(), marker),
      bounding_box(mesh::GetBoundingBox(*fespace.GetParMesh(), marker, true)), direction(3)
  {
    MFEM_VERIFY(bounding_box.planar,
                "Boundary elements must be coplanar to define a lumped element!");

    // Check that the bounding box discovered matches the area. This validates that the
    // boundary elements form a right angled quadrilateral port.
    constexpr double rel_tol = 1.0e-6;
    double A = GetArea(fespace);
    MFEM_VERIFY(std::abs(A - bounding_box.Area()) / A < rel_tol,
                "Assumed bounding box area "
                    << bounding_box.Area() << " and integrated area " << A
                    << " do not match: Port geometry is not rectangular!");

    // Check the user specified direction aligns with an axis direction.
    constexpr double angle_warning_deg = 0.1;
    constexpr double angle_error_deg = 1;
    auto lengths = bounding_box.Lengths();
    auto deviation_deg = bounding_box.Deviation(input_dir);
    if ((deviation_deg[0] > angle_warning_deg && deviation_deg[1] > angle_warning_deg) ||
        std::isnan(deviation_deg[0]) || std::isnan(deviation_deg[1]))
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
      Mpi::Warning("User specified direction {} does not align with either bounding box "
                   "axis up to {} degrees.\n"
                   "Axis 1: {} ({} degrees)\nAxis 2: {} ({} degrees)!\n",
                   input_dir, angle_warning_deg, normal_0, deviation_deg[0], normal_1,
                   deviation_deg[1]);
    }
    MFEM_VERIFY(deviation_deg[0] <= angle_error_deg || deviation_deg[1] <= angle_error_deg,
                "Specified direction does not align sufficiently with bounding box axes: "
                    << deviation_deg[0] << ' ' << deviation_deg[1] << ' ' << angle_error_deg
                    << '!');
    std::copy(input_dir.begin(), input_dir.end(), direction.begin());
    direction /= direction.Norml2();

    // Compute the length from the most aligned normal direction.
    l = lengths[deviation_deg[0] < deviation_deg[1] ? 0 : 1];
    w = A / l;
  }

  double GetGeometryLength() const override { return l; }
  double GetGeometryWidth() const override { return w; }

  std::unique_ptr<mfem::VectorCoefficient>
  GetModeCoefficient(double coef = 1.0) const override
  {
    mfem::Vector source = direction;
    source *= coef;
    return std::make_unique<mfem::VectorConstantCoefficient>(source);
  }
};

class CoaxialElementData : public LumpedElementData
{
protected:
  mesh::BoundingBall bounding_ball;  // Bounding ball defined by boundary element
  bool sign;                         // Sign of incident field, +r̂ if true
  double ra;                         // Inner radius of coaxial annulus

public:
  CoaxialElementData(const std::array<double, 3> &direction, const mfem::Array<int> &marker,
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

  double GetGeometryLength() const override { return std::log(bounding_ball.radius / ra); }
  double GetGeometryWidth() const override { return 2.0 * M_PI; }

  std::unique_ptr<mfem::VectorCoefficient>
  GetModeCoefficient(double coef = 1.0) const override
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
};

}  // namespace palace

#endif  // PALACE_FEM_LUMPED_ELEMENT_HPP
