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

  // List of all impedance boundary attributes.
  mfem::Array<int> impedance_attr, impedance_Rs_attr, impedance_Ls_attr, impedance_Cs_attr;

  // Surface properties for impedance boundary attributes: surface resistance, capacitance,
  // and inductance.
  struct ImpedanceData
  {
    double Rsinv, Lsinv, Cs;
    mfem::Array<int> attr;
  };
  std::vector<ImpedanceData> impedance_data;

  mfem::Array<int> SetUpBoundaryProperties(const IoData &iodata, const mfem::ParMesh &mesh);
  void PrintBoundaryInfo(const IoData &iodata, mfem::ParMesh &mesh);

public:
  SurfaceImpedanceOperator(const IoData &iodata, const MaterialOperator &mat_op,
                           mfem::ParMesh &mesh);

  // Returns array of surface impedance attributes.
  const auto &GetAttrList() const { return impedance_attr; }
  const auto &GetRsAttrList() const { return impedance_Rs_attr; }
  const auto &GetLsAttrList() const { return impedance_Ls_attr; }
  const auto &GetCsAttrList() const { return impedance_Cs_attr; }

  // Add contributions to system matrices from impedance boundaries with nonzero inductance,
  // capacitance, and/or resistance. For boundaries with more than R/L/C, impedances add in
  // parallel.
  void AddStiffnessBdrCoefficients(double coef, MaterialPropertyCoefficient &fb);
  void AddDampingBdrCoefficients(double coef, MaterialPropertyCoefficient &fb);
  void AddMassBdrCoefficients(double coef, MaterialPropertyCoefficient &fb);
};

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_IMPEDANCE_OPERATOR_HPP
