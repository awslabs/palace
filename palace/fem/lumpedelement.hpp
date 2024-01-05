// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_LUMPED_ELEMENT_HPP
#define PALACE_FEM_LUMPED_ELEMENT_HPP

#include <memory>
#include <mfem.hpp>
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
  // List of all boundary attributes making up this lumped element boundary.
  mfem::Array<int> attr_list;

public:
  LumpedElementData(const mfem::Array<int> &attr_list) : attr_list(attr_list) {}
  virtual ~LumpedElementData() = default;

  const auto &GetAttrList() const { return attr_list; }

  virtual double GetGeometryLength() const = 0;
  virtual double GetGeometryWidth() const = 0;

  virtual std::unique_ptr<mfem::VectorCoefficient>
  GetModeCoefficient(double coef = 1.0) const = 0;
};

class UniformElementData : public LumpedElementData
{
private:
  // Bounding box defining the rectangular lumped port.
  mesh::BoundingBox bounding_box;

  // Cartesian vector specifying signed direction of incident field.
  mfem::Vector direction;

  // Lumped element length and width.
  double l, w;

public:
  UniformElementData(const std::array<double, 3> &input_dir,
                     const mfem::Array<int> &attr_list,
                     mfem::ParFiniteElementSpace &fespace);

  double GetGeometryLength() const override { return l; }
  double GetGeometryWidth() const override { return w; }

  std::unique_ptr<mfem::VectorCoefficient>
  GetModeCoefficient(double coef = 1.0) const override;
};

class CoaxialElementData : public LumpedElementData
{
private:
  // Bounding ball defined by boundary element.
  mesh::BoundingBall bounding_ball;

  // Sign of incident field, +r̂ if true.
  bool sign;

  // Inner radius of coaxial annulus.
  double ra;

public:
  CoaxialElementData(const std::array<double, 3> &direction,
                     const mfem::Array<int> &attr_list,
                     mfem::ParFiniteElementSpace &fespace);

  double GetGeometryLength() const override { return std::log(bounding_ball.radius / ra); }
  double GetGeometryWidth() const override { return 2.0 * M_PI; }

  std::unique_ptr<mfem::VectorCoefficient>
  GetModeCoefficient(double coef = 1.0) const override;
};

}  // namespace palace

#endif  // PALACE_FEM_LUMPED_ELEMENT_HPP
