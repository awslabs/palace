// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_CONDUCTIVITY_OPERATOR_HPP
#define PALACE_MODELS_SURFACE_CONDUCTIVITY_OPERATOR_HPP

#include <complex>
#include <vector>
#include <mfem.hpp>

namespace palace
{

class IoData;
class MaterialOperator;
class MaterialPropertyCoefficient;
class Units;

namespace config
{

struct ConductivityData;

}  // namespace config

enum class ProblemType : char;

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

  void SetUpBoundaryProperties(const std::vector<config::ConductivityData> &conductivity,
                               ProblemType problem_type, const mfem::ParMesh &mesh);
  void PrintBoundaryInfo(const Units &units, const mfem::ParMesh &mesh);

  // Per-group complex surface-admittance coefficient i·ω/Z(ω) for the boundary BC,
  // evaluated at (possibly complex) ω. Z is the skin-depth surface impedance with the
  // optional finite-thickness HFSS correction. The skin-depth √ and the cosh/cos/sinh/sin
  // thickness terms accept complex ω directly (analytic continuation ω = -i·λ); for real
  // ω this returns the same scalar the real-ω stamping always used. The formula lives
  // here so both AddExtraSystemBdrCoefficients overloads share it.
  std::complex<double> EvaluateScalar(std::size_t group_idx,
                                      std::complex<double> omega) const;

public:
  SurfaceConductivityOperator(const std::vector<config::ConductivityData> &conductivity,
                              ProblemType problem_type, const Units &units,
                              const MaterialOperator &mat_op, const mfem::ParMesh &mesh);
  SurfaceConductivityOperator(const IoData &iodata, const MaterialOperator &mat_op,
                              const mfem::ParMesh &mesh);

  // Returns array of finite conductivity boundary attributes.
  mfem::Array<int> GetAttrList() const;

  // Add contributions to system matrix for a finite conductivity boundary condition.
  void AddExtraSystemBdrCoefficients(double omega, MaterialPropertyCoefficient &fbr,
                                     MaterialPropertyCoefficient &fbi);

  // Complex-ω overload, used by the eigenmode nonlinear solve where the surface-admittance
  // coefficient i·ω/Z(ω) is evaluated at the genuinely complex eigenvalue (ω = -i·λ). For
  // real ω (imag = 0) this produces the same (fbr, fbi) split as the double overload.
  void AddExtraSystemBdrCoefficients(std::complex<double> omega,
                                     MaterialPropertyCoefficient &fbr,
                                     MaterialPropertyCoefficient &fbi);
};

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_CONDUCTIVITY_OPERATOR_HPP
