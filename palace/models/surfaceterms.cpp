// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfaceterms.hpp"

#include "utils/configfile.hpp"

namespace palace
{

SurfaceTerms::SurfaceTerms(const config::BoundaryData &boundaries, ProblemType problem_type,
                           const Units &units, const MaterialOperator &mat_op,
                           const mfem::ParMesh &mesh)
  : farfield_op(boundaries.farfield, problem_type, mat_op, mesh),
    surf_sigma_op(boundaries.conductivity, problem_type, units, mat_op, mesh),
    surf_z_op(boundaries.impedance, boundaries.split_impedance_attributes, units, mat_op,
              mesh),
    surf_rz_op(boundaries.rational_impedance, boundaries.split_impedance_attributes,
               problem_type, units, mat_op, mesh)
{
}

mfem::Array<int> SurfaceTerms::GetAttrList() const
{
  mfem::Array<int> attr_list;
  attr_list.Append(farfield_op.GetAttrList());
  attr_list.Append(surf_sigma_op.GetAttrList());
  attr_list.Append(surf_z_op.GetAttrList());
  attr_list.Append(surf_rz_op.GetAttrList());
  return attr_list;
}

void SurfaceTerms::AddStiffnessBdrCoefficients(double coeff,
                                               MaterialPropertyCoefficient &fb)
{
  surf_z_op.AddStiffnessBdrCoefficients(coeff, fb);
}

void SurfaceTerms::AddDampingBdrCoefficients(double coeff, MaterialPropertyCoefficient &fb)
{
  farfield_op.AddDampingBdrCoefficients(coeff, fb);
  surf_z_op.AddDampingBdrCoefficients(coeff, fb);
}

void SurfaceTerms::AddMassBdrCoefficients(double coeff, MaterialPropertyCoefficient &fb)
{
  surf_z_op.AddMassBdrCoefficients(coeff, fb);
}

template <typename OmegaType>
void SurfaceTerms::AddExtraSystemBdrMassImpl(OmegaType omega,
                                             MaterialPropertyCoefficient &fbr,
                                             MaterialPropertyCoefficient &fbi)
{
  surf_sigma_op.AddExtraSystemBdrCoefficients(omega, fbr, fbi);
  surf_rz_op.AddExtraSystemBdrCoefficients(omega, fbr, fbi);
}

void SurfaceTerms::AddExtraSystemBdrCoefficients(
    double omega, MaterialPropertyCoefficient &dfbr, MaterialPropertyCoefficient &dfbi,
    MaterialPropertyCoefficient &fbr, MaterialPropertyCoefficient &fbi,
    MaterialPropertyCoefficient &frbr, MaterialPropertyCoefficient &frbi,
    MaterialPropertyCoefficient *fpbr, MaterialPropertyCoefficient *fpbi)
{
  farfield_op.AddExtraSystemBdrCoefficients(omega, dfbr, dfbi, frbr, frbi, fpbr, fpbi);
  AddExtraSystemBdrMassImpl(omega, fbr, fbi);
}

void SurfaceTerms::AddExtraSystemBdrCoefficients(
    std::complex<double> omega, MaterialPropertyCoefficient &dfbr,
    MaterialPropertyCoefficient &dfbi, MaterialPropertyCoefficient &fbr,
    MaterialPropertyCoefficient &fbi, MaterialPropertyCoefficient &frbr,
    MaterialPropertyCoefficient &frbi, MaterialPropertyCoefficient *fpbr,
    MaterialPropertyCoefficient *fpbi)
{
  farfield_op.AddExtraSystemBdrCoefficients(omega, dfbr, dfbi, frbr, frbi, fpbr, fpbi);
  AddExtraSystemBdrMassImpl(omega, fbr, fbi);
}

}  // namespace palace
