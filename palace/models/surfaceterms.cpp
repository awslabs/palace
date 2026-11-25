// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfaceterms.hpp"

#include <fmt/format.h>
#include <fmt/ranges.h>
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

void SurfaceTerms::AddExtraSystemBdrMassImpl(std::complex<double> omega, bool scalar,
                                             bool preconditioner,
                                             MaterialPropertyCoefficient &fbr,
                                             MaterialPropertyCoefficient &fbi)
{
  // Each term is its unit coefficient times its scale at ω: the real part in fbr, then the
  // imaginary part in fbi, also when one is zero (except that the vector coefficient of an
  // absorbing boundary skips a zero part, FarfieldBoundaryOperator::
  // AddDispersiveDampingBdrCoefficients). At a real frequency, a rational impedance
  // boundary is checked for passivity. In a preconditioner, the term of the absorbing
  // boundaries next to a frequency-dependent material is the magnitude of its scale in fbi.
  const int num_sigma = static_cast<int>(sigma_terms.size());
  const int num_rz = surf_rz_op.GetNumBoundaries();
  auto add = [this, scalar](int i, MaterialPropertyCoefficient &fb, double coeff)
  {
    if (scalar)
    {
      AddScalarExtraSystemBdrMassTermCoefficients(i, fb, coeff);
    }
    else
    {
      AddExtraSystemBdrMassTermCoefficients(i, fb, coeff);
    }
  };
  for (int i = 0; i < GetNumExtraSystemBdrMassTerms(); i++)
  {
    const std::complex<double> scale = GetExtraSystemBdrMassTermScale(i, omega);
    if (i >= num_sigma && i < num_sigma + num_rz && omega.imag() == 0.0)
    {
      surf_rz_op.CheckPassivity(i - num_sigma, omega.real());
    }
    if (preconditioner && i >= num_sigma + num_rz)
    {
      add(i, fbi, std::abs(scale));
    }
    else
    {
      add(i, fbr, scale.real());
      add(i, fbi, scale.imag());
    }
  }
}

void SurfaceTerms::AddExtraSystemBdrMassCoefficients(double omega,
                                                     MaterialPropertyCoefficient &fbr,
                                                     MaterialPropertyCoefficient &fbi)
{
  AddExtraSystemBdrMassImpl(omega, false, false, fbr, fbi);
}

void SurfaceTerms::AddScalarExtraSystemBdrMassCoefficients(double omega,
                                                           MaterialPropertyCoefficient &fbr,
                                                           MaterialPropertyCoefficient &fbi)
{
  AddExtraSystemBdrMassImpl(omega, true, false, fbr, fbi);
}

void SurfaceTerms::AddPreconditionerExtraSystemBdrMassCoefficients(
    double omega, MaterialPropertyCoefficient &fbr, MaterialPropertyCoefficient &fbi)
{
  AddExtraSystemBdrMassImpl(omega, false, true, fbr, fbi);
}

void SurfaceTerms::AddScalarPreconditionerExtraSystemBdrMassCoefficients(
    double omega, MaterialPropertyCoefficient &fbr, MaterialPropertyCoefficient &fbi)
{
  AddExtraSystemBdrMassImpl(omega, true, true, fbr, fbi);
}

void SurfaceTerms::AddExtraSystemBdrCoefficients(
    double omega, MaterialPropertyCoefficient &dfbr, MaterialPropertyCoefficient &dfbi,
    MaterialPropertyCoefficient &fbr, MaterialPropertyCoefficient &fbi,
    MaterialPropertyCoefficient &frbr, MaterialPropertyCoefficient &frbi,
    MaterialPropertyCoefficient *fpbr, MaterialPropertyCoefficient *fpbi)
{
  farfield_op.AddExtraSystemBdrCoefficients(omega, dfbr, dfbi, frbr, frbi, fpbr, fpbi);
  AddExtraSystemBdrMassImpl(omega, false, false, fbr, fbi);
}

void SurfaceTerms::AddExtraSystemBdrCoefficients(
    std::complex<double> omega, MaterialPropertyCoefficient &dfbr,
    MaterialPropertyCoefficient &dfbi, MaterialPropertyCoefficient &fbr,
    MaterialPropertyCoefficient &fbi, MaterialPropertyCoefficient &frbr,
    MaterialPropertyCoefficient &frbi, MaterialPropertyCoefficient *fpbr,
    MaterialPropertyCoefficient *fpbi)
{
  farfield_op.AddExtraSystemBdrCoefficients(omega, dfbr, dfbi, frbr, frbi, fpbr, fpbi);
  AddExtraSystemBdrMassImpl(omega, false, false, fbr, fbi);
}

void SurfaceTerms::AddStaticExtraSystemBdrMassCoefficients(MaterialPropertyCoefficient &fb,
                                                           double nu2) const
{
  const int num_sigma = static_cast<int>(sigma_terms.size());
  for (int b = 0; b < surf_rz_op.GetNumBoundaries(); b++)
  {
    // The quotient q of the Robin coefficient is stored with its highest coefficient first.
    const auto &q = surf_rz_op.GetRobinQuotient(b);
    const double q2 = (q.size() >= 3) ? q[q.size() - 3] : 0.0;
    const double g = surf_rz_op.GetStaticRobinCoefficient(b) + nu2 * q2;
    if (g != 0.0)
    {
      AddExtraSystemBdrMassTermCoefficients(num_sigma + b, fb, g);
    }
  }
}

void SurfaceTerms::AddRealPreconditionerExtraSystemBdrCoefficients(
    double omega, MaterialPropertyCoefficient &dfb, MaterialPropertyCoefficient &fb,
    MaterialPropertyCoefficient &frb)
{
  farfield_op.AddPreconditionerBdrCoefficients(omega, dfb, frb);
  AddExtraSystemBdrMassImpl(omega, false, true, fb, fb);
}

void SurfaceTerms::AddRealPreconditionerExtraSystemBdrCoefficients(
    std::complex<double> omega, MaterialPropertyCoefficient &dfb,
    MaterialPropertyCoefficient &fb, MaterialPropertyCoefficient &frb)
{
  farfield_op.AddPreconditionerBdrCoefficients(omega, dfb, frb);
  AddExtraSystemBdrMassImpl(omega, false, true, fb, fb);
}

void SurfaceTerms::AddFlippedPreconditionerBdrCoefficients(
    std::complex<double> omega, MaterialPropertyCoefficient &dfbr,
    MaterialPropertyCoefficient &dfbi) const
{
  farfield_op.AddFlippedPreconditionerBdrCoefficients(omega, dfbr, dfbi);
}

int SurfaceTerms::GetNumExtraSystemBdrMassTerms() const
{
  return static_cast<int>(sigma_terms.size()) + surf_rz_op.GetNumBoundaries() +
         static_cast<int>(farfield_op.NumDispersiveMaterials());
}

void SurfaceTerms::AddExtraSystemBdrMassTermCoefficients(int i,
                                                         MaterialPropertyCoefficient &fb,
                                                         double coeff) const
{
  const int num_sigma = static_cast<int>(sigma_terms.size());
  const int num_rz = surf_rz_op.GetNumBoundaries();
  if (i < num_sigma)
  {
    surf_sigma_op.AddBoundaryMassBdrCoefficients(sigma_terms[i], fb, coeff);
  }
  else if (i < num_sigma + num_rz)
  {
    surf_rz_op.AddUnitBdrCoefficients(i - num_sigma, fb, coeff);
  }
  else
  {
    farfield_op.AddDispersiveDampingBdrCoefficients(
        static_cast<std::size_t>(i - num_sigma - num_rz), coeff, fb);
  }
}

void SurfaceTerms::AddScalarExtraSystemBdrMassTermCoefficients(
    int i, MaterialPropertyCoefficient &fb, double coeff) const
{
  // The coefficients of finite conductivity and rational impedance boundaries are scalar,
  // the same for both fields. An absorbing boundary has the scalar inverse impedance of
  // its material.
  const int first_ff = static_cast<int>(sigma_terms.size()) + surf_rz_op.GetNumBoundaries();
  if (i < first_ff)
  {
    AddExtraSystemBdrMassTermCoefficients(i, fb, coeff);
  }
  else
  {
    farfield_op.AddDispersiveScalarDampingBdrCoefficients(
        static_cast<std::size_t>(i - first_ff), coeff, fb);
  }
}

std::complex<double>
SurfaceTerms::GetExtraSystemBdrMassTermScale(int i, std::complex<double> omega) const
{
  const int num_sigma = static_cast<int>(sigma_terms.size());
  const int num_rz = surf_rz_op.GetNumBoundaries();
  if (i < num_sigma)
  {
    return surf_sigma_op.EvaluateScalar(sigma_terms[i], omega);
  }
  if (i < num_sigma + num_rz)
  {
    return surf_rz_op.EvalRobinCoefficient(i - num_sigma,
                                           std::complex<double>(0.0, 1.0) * omega);
  }
  return std::complex<double>(0.0, 1.0) * omega *
         farfield_op.GetDispersiveScale(static_cast<std::size_t>(i - num_sigma - num_rz),
                                        omega);
}

std::string SurfaceTerms::GetExtraSystemBdrMassTermLabel(int i) const
{
  const int num_sigma = static_cast<int>(sigma_terms.size());
  const int num_rz = surf_rz_op.GetNumBoundaries();
  if (i < num_sigma)
  {
    return fmt::format("surfsigma_{:d}", sigma_terms[i]);
  }
  if (i < num_sigma + num_rz)
  {
    return fmt::format("rationalz_{:d}", i - num_sigma);
  }
  return fmt::format("farfield_permittivity_{:d}", i - num_sigma - num_rz);
}

std::string SurfaceTerms::GetExtraSystemBdrMassTermName(int i) const
{
  const int num_sigma = static_cast<int>(sigma_terms.size());
  const int num_rz = surf_rz_op.GetNumBoundaries();
  if (i < num_sigma)
  {
    return fmt::format("finite conductivity boundary (attribute {})",
                       fmt::join(surf_sigma_op.GetAttrList(sigma_terms[i]), ", "));
  }
  if (i < num_sigma + num_rz)
  {
    return fmt::format("rational impedance boundary (attribute {})",
                       fmt::join(surf_rz_op.GetAttrList(i - num_sigma), ", "));
  }
  return fmt::format(
      "absorbing boundary next to material {:d}",
      farfield_op.GetDispersiveMaterial(static_cast<std::size_t>(i - num_sigma - num_rz)) +
          1);
}

std::optional<SurfaceTerms::ExactSplit>
SurfaceTerms::GetExtraSystemBdrMassTermSplit(int i) const
{
  const int num_sigma = static_cast<int>(sigma_terms.size());
  if (i < num_sigma || i >= num_sigma + surf_rz_op.GetNumBoundaries())
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
