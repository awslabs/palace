// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LUMPED_ELEMENT_HPP
#define PALACE_LUMPED_ELEMENT_HPP

#include <memory>
#include <string>
#include <mfem.hpp>
#include "utils/geodata.hpp"
#include "utils/mfemintegrators.hpp"

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
    mfem::ParGridFunction ones(&fespace);
    ones.mfem::Vector::operator=(1.0);
    mfem::ParLinearForm s(&fespace);
    mfem::ConstantCoefficient one_func(1.0);
    // NOLINTNEXTLINE(cppcoreguidelines-owning-memory)
    s.AddBoundaryIntegrator(new BoundaryLFIntegrator(one_func), attr_marker);
    s.UseFastAssembly(true);
    s.Assemble();
    return s(ones);
  }

public:
  LumpedElementData(int d, mfem::Array<int> marker) : dim(d), attr_marker(std::move(marker))
  {
  }
  virtual ~LumpedElementData() = default;

  LumpedElementData() = delete;
  LumpedElementData(const LumpedElementData &) = default;
  LumpedElementData(LumpedElementData &&) = default;
  LumpedElementData &operator=(const LumpedElementData &) = delete;
  LumpedElementData &operator=(LumpedElementData &&) = delete;

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
  bool sign;      // Sign of incident field, +x̂ / ŷ / ẑ if true
  int component;  // Lumped element direction (0: x, 1: y, 2: z)
  double l, w;    // Lumped element length and width

public:
  UniformElementData(const std::string &direction, const mfem::Array<int> &marker,
                     mfem::ParFiniteElementSpace &fespace)
    : LumpedElementData(fespace.GetParMesh()->SpaceDimension(), marker),
      sign(direction[0] == '+')
  {
    switch (direction[1])
    {
      case 'x':
        component = 0;
        break;
      case 'y':
        component = 1;
        break;
      case 'z':
        component = 2;
        break;
      default:
        MFEM_ABORT("Lumped element direction is not correctly formatted!");
        component = 0;  // For compiler warning
        break;
    }

    // Get the lumped element length and width assuming axis-aligned rectangle.
    mfem::Vector bbmin, bbmax;
    mesh::GetBoundingBox(*fespace.GetParMesh(), marker, true, bbmin, bbmax);
    double A = GetArea(fespace);
    l = bbmax(component) - bbmin(component);
    w = A / l;
  }

  double GetGeometryLength() const override { return l; }
  double GetGeometryWidth() const override { return w; }

  std::unique_ptr<mfem::VectorCoefficient>
  GetModeCoefficient(double coef = 1.0) const override
  {
    mfem::Vector source(dim);
    source = 0.0;
    source(component) = (sign ? 1.0 : -1.0) * coef;
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
  CoaxialElementData(const std::string &direction, const mfem::Array<int> &marker,
                     mfem::ParFiniteElementSpace &fespace)
    : LumpedElementData(fespace.GetParMesh()->SpaceDimension(), marker),
      sign(direction[0] == '+')
  {
    // Get the outer annulus radius.
    MFEM_VERIFY(direction[1] == 'r',
                "Lumped element direction is not correctly formatted!");
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

#endif  // PALACE_LUMPED_ELEMENT_HPP
