// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_LUMPED_ELEMENT_HPP
#define PALACE_FEM_LUMPED_ELEMENT_HPP

#include <memory>
#include <mfem.hpp>

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
  GetModeCoefficient(double coeff = 1.0) const = 0;
};

class UniformElementData : public LumpedElementData
{
private:
  // Cartesian vector specifying signed direction of incident field.
  mfem::Vector direction;

  // Lumped element length and width.
  double l, w;

public:
  UniformElementData(const std::array<double, 3> &input_dir,
                     const mfem::Array<int> &attr_list, const mfem::ParMesh &mesh);

  double GetGeometryLength() const override { return l; }
  double GetGeometryWidth() const override { return w; }

  std::unique_ptr<mfem::VectorCoefficient>
  GetModeCoefficient(double coeff = 1.0) const override;
};

class CoaxialElementData : public LumpedElementData
{
private:
  // Sign of incident field, +1 for +r̂, -1 for -r̂.
  double direction;

  // Origin of the coaxial annulus.
  mfem::Vector origin;

  // Outer and inner radii of coaxial annulus.
  double r_outer, r_inner;

public:
  CoaxialElementData(const std::array<double, 3> &input_dir,
                     const mfem::Array<int> &attr_list, const mfem::ParMesh &mesh);

  double GetGeometryLength() const override { return std::log(r_outer / r_inner); }
  double GetGeometryWidth() const override { return 2.0 * M_PI; }

  std::unique_ptr<mfem::VectorCoefficient>
  GetModeCoefficient(double coeff = 1.0) const override;
};

}  // namespace palace

#endif  // PALACE_FEM_LUMPED_ELEMENT_HPP
