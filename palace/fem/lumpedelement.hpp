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
  double l, w;             // Lumped element length and width

public:
  UniformElementData(const std::array<double, 3> &input_dir, const mfem::Array<int> &marker,
                     mfem::ParFiniteElementSpace &fespace)
    : LumpedElementData(fespace.GetParMesh()->SpaceDimension(), marker), direction(3)
  {
    std::copy(input_dir.begin(), input_dir.end(), direction.begin());

    // Get the lumped element length and width.
    mfem::Vector bbmin, bbmax;
    mesh::GetBoundingBox(*fespace.GetParMesh(), marker, true, bbmin, bbmax);
    double A = GetArea(fespace);
    bbmax -= bbmin;
    l = std::abs(bbmax * direction);
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
  bool sign;       // Sign of incident field, +r̂ if true
  double ra, rb;   // Inner and outer radius of coaxial annulus
  mfem::Vector c;  // Center coordinates of coaxial annulus

public:
  CoaxialElementData(const std::array<double, 3> &direction, const mfem::Array<int> &marker,
                     mfem::ParFiniteElementSpace &fespace)
    : LumpedElementData(fespace.GetParMesh()->SpaceDimension(), marker), sign(direction[0] > 0)
  {
    // Get the outer annulus radius.
    mfem::Vector bbmin, bbmax;
    mesh::GetBoundingBox(*fespace.GetParMesh(), marker, true, bbmin, bbmax);
    double A = GetArea(fespace);
    rb = 0.0;
    for (int d = 0; d < dim; d++)
    {
      double diff = 0.5 * (bbmax(d) - bbmin(d));
      if (diff > rb)
      {
        rb = diff;
      }
    }

    // Get inner radius of annulus assuming full 2π circumference.
    MFEM_VERIFY(rb > 0.0 && rb * rb - A / M_PI > 0.0,
                "Coaxial element boundary is not defined correctly!");
    ra = std::sqrt(rb * rb - A / M_PI);
    c.SetSize(dim);
    for (int d = 0; d < dim; d++)
    {
      c(d) = 0.5 * (bbmax(d) + bbmin(d));
    }
  }

  double GetGeometryLength() const override { return std::log(rb / ra); }
  double GetGeometryWidth() const override { return 2.0 * M_PI; }

  std::unique_ptr<mfem::VectorCoefficient>
  GetModeCoefficient(double coef = 1.0) const override
  {
    double scoef = (sign ? 1.0 : -1.0) * coef;
    mfem::Vector x0(c);
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
