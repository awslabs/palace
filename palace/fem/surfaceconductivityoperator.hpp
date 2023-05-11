// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_SURF_CONDUCTIVITY_OPERATOR_HPP
#define PALACE_FEM_SURF_CONDUCTIVITY_OPERATOR_HPP

#include <mfem.hpp>

namespace palace
{

class IoData;
class SumMatrixCoefficient;

//
// A class handling finite conductivity boundaries.
//
class SurfaceConductivityOperator
{
private:
  // Surface properties for finite conductivity boundary attributes: conductor conductivity
  // and permeability, and (optionally) thickness.
  mfem::Vector bdr_sigma, bdr_mu, bdr_h;
  mfem::Array<int> conductivity_marker;
  void SetUpBoundaryProperties(const IoData &iodata, const mfem::ParMesh &mesh);
  void PrintBoundaryInfo(const IoData &iodata, mfem::ParMesh &mesh);

public:
  SurfaceConductivityOperator(const IoData &iodata, mfem::ParMesh &mesh);

  // Returns array marking finite conductivity boundary attributes.
  const mfem::Array<int> &GetMarker() const { return conductivity_marker; }

  // Add contributions to system matrix for a finite conductivity boundary condition.
  void AddExtraSystemBdrCoefficients(double omega, SumMatrixCoefficient &fbr,
                                     SumMatrixCoefficient &fbi);
};

}  // namespace palace

#endif  // PALACE_FEM_SURF_CONDUCTIVITY_OPERATOR_HPP
