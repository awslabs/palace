// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_CONDUCTIVITY_OPERATOR_HPP
#define PALACE_MODELS_SURFACE_CONDUCTIVITY_OPERATOR_HPP

#include <vector>
#include <mfem.hpp>

namespace palace
{

class IoData;
class MaterialOperator;
class MaterialPropertyCoefficient;

//
// A class handling finite conductivity boundaries.
//
class SurfaceConductivityOperator
{
private:
  // Reference to material property data (not owned).
  const MaterialOperator &mat_op;

  // Surface properties for finite conductivity boundary attributes: conductor conductivity
  // and permeability, and (optionally) thickness.
  struct ConductivityData
  {
    double sigma, mu, h;
    mfem::Array<int> attr_list;
  };
  std::vector<ConductivityData> boundaries;

  void SetUpBoundaryProperties(const IoData &iodata, const mfem::ParMesh &mesh);
  void PrintBoundaryInfo(const IoData &iodata, const mfem::ParMesh &mesh);

public:
  SurfaceConductivityOperator(const IoData &iodata, const MaterialOperator &mat_op,
                              const mfem::ParMesh &mesh);

  // Returns array of finite conductivity boundary attributes.
  mfem::Array<int> GetAttrList() const;

  // Add contributions to system matrix for a finite conductivity boundary condition.
  void AddExtraSystemBdrCoefficients(double omega, MaterialPropertyCoefficient &fbr,
                                     MaterialPropertyCoefficient &fbi);
};

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_CONDUCTIVITY_OPERATOR_HPP
