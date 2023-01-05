// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_SURF_IMPEDANCE_OPERATOR_HPP
#define PALACE_SURF_IMPEDANCE_OPERATOR_HPP

#include <mfem.hpp>

namespace palace
{

class IoData;
class SumMatrixCoefficient;

//
// A class handling impedance boundaries.
//
class SurfaceImpedanceOperator
{
private:
  // Surface properties for impedance boundary attributes: surface resistance, capacitance,
  // and inductance.
  mfem::Vector Z_Rsinv, Z_Lsinv, Z_Cs;
  mfem::Array<int> impedance_marker, impedance_Rs_marker, impedance_Ls_marker,
      impedance_Cs_marker;
  void SetUpBoundaryProperties(const IoData &iodata, const mfem::ParMesh &mesh);
  void PrintBoundaryInfo(const IoData &iodata, mfem::ParMesh &mesh);

public:
  SurfaceImpedanceOperator(const IoData &iodata, mfem::ParMesh &mesh);

  // Returns array marking surface impedance attributes.
  const mfem::Array<int> &GetMarker() const { return impedance_marker; }
  const mfem::Array<int> &GetRsMarker() const { return impedance_Rs_marker; }
  const mfem::Array<int> &GetLsMarker() const { return impedance_Ls_marker; }
  const mfem::Array<int> &GetCsMarker() const { return impedance_Cs_marker; }

  // Add contributions to system matrices from impedance boundaries with nonzero inductance,
  // capacitance, and/or resistance. For boundaries with more than R/L/C, impedances add in
  // parallel.
  void AddStiffnessBdrCoefficients(double coef, SumMatrixCoefficient &fb);
  void AddMassBdrCoefficients(double coef, SumMatrixCoefficient &fb);
  void AddDampingBdrCoefficients(double coef, SumMatrixCoefficient &fb);
};

}  // namespace palace

#endif  // PALACE_SURF_IMPEDANCE_OPERATOR_HPP
