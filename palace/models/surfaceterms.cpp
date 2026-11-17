// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfaceterms.hpp"

#include <fmt/format.h>
#include "utils/configfile.hpp"
#include "utils/iodata.hpp"

namespace palace
{

namespace
{

// Surface conductivity groups with nonzero conductivity.
std::vector<std::size_t> GetActiveGroups(const SurfaceConductivityOperator &surf_sigma_op)
{
  std::vector<std::size_t> groups;
  for (std::size_t g = 0; g < surf_sigma_op.Size(); g++)
  {
    if (surf_sigma_op.IsActive(g))
    {
      groups.push_back(g);
    }
  }
  return groups;
}

}  // namespace

SurfaceTerms::SurfaceTerms(const config::BoundaryData &boundaries, ProblemType problem_type,
                           const Units &units, const MaterialOperator &mat_op,
                           const mfem::ParMesh &mesh)
  : farfield_op(boundaries.farfield, problem_type, mat_op, mesh),
    surf_sigma_op(boundaries.conductivity, problem_type, units, mat_op, mesh),
    surf_z_op(boundaries.impedance, boundaries.split_impedance_attributes, units, mat_op,
              mesh),
    surf_rz_op(boundaries.rational_impedance, boundaries.split_impedance_attributes,
               problem_type, units, mat_op, mesh),
    sigma_terms(GetActiveGroups(surf_sigma_op))
{
}

SurfaceTerms::SurfaceTerms(const IoData &iodata, const MaterialOperator &mat_op,
                           const mfem::ParMesh &mesh)
  : farfield_op(iodata, mat_op, mesh), surf_sigma_op(iodata, mat_op, mesh),
    surf_z_op(iodata, mat_op, mesh), surf_rz_op(iodata, mat_op, mesh),
    sigma_terms(GetActiveGroups(surf_sigma_op))
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

bool SurfaceTerms::HasDamping() const
{
  return surf_z_op.GetRsAttrList().Size() > 0 || farfield_op.GetAttrList().Size() > 0;
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

void SurfaceTerms::AddScalarDampingBdrCoefficients(double coeff,
                                                   MaterialPropertyCoefficient &fb)
{
  farfield_op.AddScalarDampingBdrCoefficients(coeff, fb);
  surf_z_op.AddDampingBdrCoefficients(coeff, fb);
}

template <typename OmegaType>
void SurfaceTerms::AddExtraSystemBdrMassImpl(OmegaType omega,
                                             MaterialPropertyCoefficient &fbr,
                                             MaterialPropertyCoefficient &fbi)
{
  surf_sigma_op.AddExtraSystemBdrCoefficients(omega, fbr, fbi);
  surf_rz_op.AddExtraSystemBdrCoefficients(omega, fbr, fbi);
}

void SurfaceTerms::AddExtraSystemBdrMassCoefficients(double omega,
                                                     MaterialPropertyCoefficient &fbr,
                                                     MaterialPropertyCoefficient &fbi)
{
  AddExtraSystemBdrMassImpl(omega, fbr, fbi);
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

int SurfaceTerms::GetNumExtraSystemBdrMassTerms() const
{
  return static_cast<int>(sigma_terms.size()) + surf_rz_op.GetNumBoundaries();
}

void SurfaceTerms::AddExtraSystemBdrMassTermCoefficients(
    int i, MaterialPropertyCoefficient &fb) const
{
  const int num_sigma = static_cast<int>(sigma_terms.size());
  if (i < num_sigma)
  {
    surf_sigma_op.AddBoundaryMassBdrCoefficients(sigma_terms[i], fb);
  }
  else
  {
    surf_rz_op.AddUnitBdrCoefficients(i - num_sigma, fb);
  }
}

void SurfaceTerms::AddScalarExtraSystemBdrMassTermCoefficients(
    int i, MaterialPropertyCoefficient &fb) const
{
  // The coefficients of both types are scalar, the same for both fields.
  AddExtraSystemBdrMassTermCoefficients(i, fb);
}

std::complex<double>
SurfaceTerms::GetExtraSystemBdrMassTermScale(int i, std::complex<double> omega) const
{
  const int num_sigma = static_cast<int>(sigma_terms.size());
  if (i < num_sigma)
  {
    return surf_sigma_op.EvaluateScalar(sigma_terms[i], omega);
  }
  return surf_rz_op.EvalRobinCoefficient(i - num_sigma,
                                         std::complex<double>(0.0, 1.0) * omega);
}

std::string SurfaceTerms::GetExtraSystemBdrMassTermLabel(int i) const
{
  const int num_sigma = static_cast<int>(sigma_terms.size());
  if (i < num_sigma)
  {
    return fmt::format("surfsigma_{:d}", sigma_terms[i]);
  }
  return fmt::format("rationalz_{:d}", i - num_sigma);
}

std::optional<SurfaceTerms::ExactSplit>
SurfaceTerms::GetExtraSystemBdrMassTermSplit(int i) const
{
  const int num_sigma = static_cast<int>(sigma_terms.size());
  if (i < num_sigma)
  {
    return std::nullopt;
  }
  const int b = i - num_sigma;
  return ExactSplit{
      surf_rz_op.GetRobinQuotient(b), surf_rz_op.GetRobinQuotientDegree(b),
      [this, b](std::complex<double> s) { return surf_rz_op.EvalRobinCoefficient(b, s); },
      [this, b](std::complex<double> s) { return surf_rz_op.EvalRobinRemainder(b, s); }};
}

}  // namespace palace
