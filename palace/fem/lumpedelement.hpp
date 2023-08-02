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
  mfem::Vector direction;  // Cartesian vector specifying signed direction of incident field
  mesh::BoundingBox bounding_box;  // Bounding box defining the rectangular lumped port
  double l, w;                     // Lumped element length and width
public:
  UniformElementData(const std::array<double, 3> &input_dir, const mfem::Array<int> &marker,
                     mfem::ParFiniteElementSpace &fespace)
    : LumpedElementData(fespace.GetParMesh()->SpaceDimension(), marker), direction(3),
      bounding_box(mesh::GetBoundingBox(*fespace.GetParMesh(), marker, true))
  {
    std::copy(input_dir.begin(), input_dir.end(), direction.begin());

    MFEM_VERIFY(bounding_box.planar, "The set of boundary elements must be coplanar");

    double A = GetArea(fespace);

    // Check that the bounding box discovered matches the area. This validates
    // that the boundary elements form a right angled quadrilateral port.
    constexpr double rel_tol = 1e-6;
    MFEM_VERIFY(std::abs(A - bounding_box.Area()) / A < rel_tol,
                "Assumed bounding box area "
                    << bounding_box.Area() << " and integrated area " << A
                    << " do not match: Port geometry is not rectangular");

    // Pick normal most aligned with direction -> should be forgiving to errors
    // in specification of direction.
    auto Dot = [](const auto &x_1, const auto &x_2)
    { return x_1[0] * x_2[0] + x_1[1] * x_2[1] + x_1[2] * x_2[2]; };

    const auto &normal = std::abs(Dot(input_dir, bounding_box.normals[0])) >
                                 std::abs(Dot(input_dir, bounding_box.normals[1]))
                             ? bounding_box.normals[0]
                             : bounding_box.normals[1];
    std::copy(normal.begin(), normal.end(), direction.begin());

    if ((direction[0] * input_dir[0] + direction[1] * input_dir[1] +
         direction[2] * input_dir[2]) < 0)
    {
      // Ensure the correct orientation of the direction was chosen.
      direction *= -1.0;
    }

    // Get the lumped element length and width.
    l = std::sqrt(Dot(direction, direction)) * 2;
    w = A / l;

    // Normalize the direction vector
    direction /= direction.Norml2();
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
  bool sign;                         // Sign of incident field, +r̂ if true
  double ra;                         // Inner radius of coaxial annulus
  mesh::BoundingBall bounding_ball;  // Bounding ball defined by boundary element

public:
  CoaxialElementData(const std::array<double, 3> &direction, const mfem::Array<int> &marker,
                     mfem::ParFiniteElementSpace &fespace)
    : LumpedElementData(fespace.GetParMesh()->SpaceDimension(), marker),
      sign(direction[0] > 0),
      bounding_ball(mesh::GetBoundingBall(*fespace.GetParMesh(), marker, true))
  {
    MFEM_VERIFY(std::abs(bounding_ball.planar_normal[0]) == 0 &&
                    std::abs(bounding_ball.planar_normal[1]) == 0 &&
                    std::abs(bounding_ball.planar_normal[2]) == 0,
                "Boundary elements must be coplanar to define a coaxial port.");

    double A = GetArea(fespace);
    // Get inner radius of annulus assuming full 2π circumference.
    MFEM_VERIFY(bounding_ball.radius > 0.0 &&
                    std::pow(bounding_ball.radius, 2) - A / M_PI > 0.0,
                "Coaxial element boundary is not defined correctly!");
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
    auto Source = [scoef, x0](const mfem::Vector &x, mfem::Vector &f)
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
