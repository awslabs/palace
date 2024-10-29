// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_PERIODIC_BOUNDARY_OPERATOR_HPP
#define PALACE_MODELS_PERIODIC_BOUNDARY_OPERATOR_HPP

#include <mfem.hpp>

namespace palace
{

class IoData;
class MaterialOperator;
class MaterialPropertyCoefficient;

//
// A class handling periodic boundaries.
//
class PeriodicBoundaryOperator
{
private:
  // Reference to material property data (not owned).
  const MaterialOperator &mat_op;

  // List of all periodic boundary condition attributes.
  mfem::Array<int> periodic_attr;

  // Bloch wave vector for Floquet boundary conditions.
  mfem::Vector wave_vector;

  // Matrix representation of cross product with the wave_vector;
  mfem::DenseMatrix wave_vector_cross;

  mfem::Array<int> SetUpBoundaryProperties(const IoData &iodata, const mfem::ParMesh &mesh);

public:
  PeriodicBoundaryOperator(const IoData &iodata, const MaterialOperator &mat_op,
                           const mfem::ParMesh &mesh);

  // Returns array of periodic BC attributes.
  const auto &GetAttrList() const { return periodic_attr; }

  // Add contributions to system matrices
  void AddRealMassCoefficients(double coeff, MaterialPropertyCoefficient &f);
  void AddWeakCurlCoefficients(double coeff, MaterialPropertyCoefficient &f);
  void AddCurlCoefficients(double coeff, MaterialPropertyCoefficient &f);
};

}  // namespace palace

#endif  // PALACE_MODELS_PERIODIC_BOUNDARY_OPERATOR_HPP
