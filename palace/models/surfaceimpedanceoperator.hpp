// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_IMPEDANCE_OPERATOR_HPP
#define PALACE_MODELS_SURFACE_IMPEDANCE_OPERATOR_HPP

#include <vector>
#include <mfem.hpp>

namespace palace
{

class IoData;
class MaterialOperator;
class MaterialPropertyCoefficient;

//
// A class handling impedance boundaries.
//
class SurfaceImpedanceOperator
{
private:
  // Reference to material property data (not owned).
  const MaterialOperator &mat_op;

  // Surface properties for impedance boundary attributes: surface resistance, capacitance,
  // and inductance.
  struct ImpedanceData
  {
    double Rs, Ls, Cs;
    mfem::Array<int> attr_list;
    double scaling = 1.0;
  };
  std::vector<ImpedanceData> boundaries;

  void SetUpBoundaryProperties(const IoData &iodata, const mfem::ParMesh &mesh);
  void PrintBoundaryInfo(const IoData &iodata, const mfem::ParMesh &mesh);

public:
  SurfaceImpedanceOperator(const IoData &iodata, const MaterialOperator &mat_op,
                           const mfem::ParMesh &mesh);

  // Returns array of surface impedance attributes.
  mfem::Array<int> GetAttrList() const;
  mfem::Array<int> GetRsAttrList() const;
  mfem::Array<int> GetLsAttrList() const;
  mfem::Array<int> GetCsAttrList() const;

  // Add contributions to system matrices from impedance boundaries with nonzero inductance,
  // resistance, and/or capacitance. For boundaries with more than R/L/C, impedances add in
  // parallel.
  void AddStiffnessBdrCoefficients(double coeff, MaterialPropertyCoefficient &fb);
  void AddDampingBdrCoefficients(double coeff, MaterialPropertyCoefficient &fb);
  void AddMassBdrCoefficients(double coeff, MaterialPropertyCoefficient &fb);
};

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_IMPEDANCE_OPERATOR_HPP
