// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FARFIELD_BOUNDARY_OPERATOR_HPP
#define PALACE_FARFIELD_BOUNDARY_OPERATOR_HPP

#include <mfem.hpp>

namespace palace
{

class IoData;
class MaterialOperator;
class SumMatrixCoefficient;

//
// A class handling farfield, or absorbing, boundaries.
//
class FarfieldBoundaryOperator
{
private:
  // Reference to input data (not owned).
  const MaterialOperator &mat_op;

  // First- or second-order absorbing boundary condition.
  int order;

  // Marker for all absorbing boundary condition attributes.
  mfem::Array<int> farfield_marker;
  void SetUpBoundaryProperties(const IoData &iodata, const mfem::ParMesh &mesh);

public:
  FarfieldBoundaryOperator(const IoData &iodata, const MaterialOperator &mat,
                           const mfem::ParMesh &mesh);

  // Returns order of absorbing BC approximation.
  int GetOrder() const { return order; }

  // Returns array marking farfield BC attributes.
  const mfem::Array<int> &GetMarker() const { return farfield_marker; }

  // Add contributions to system matrices from first- or second-order absorbing boundary
  // condition.
  void AddDampingBdrCoefficients(double coef, SumMatrixCoefficient &fb);
  void AddExtraSystemBdrCoefficients(double omega, SumMatrixCoefficient &dfbr,
                                     SumMatrixCoefficient &dfbi);
};

}  // namespace palace

#endif  // PALACE_FARFIELD_BOUNDARY_OPERATOR_HPP
