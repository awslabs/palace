// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_FARFIELD_BOUNDARY_OPERATOR_HPP
#define PALACE_MODELS_FARFIELD_BOUNDARY_OPERATOR_HPP

#include <complex>
#include <utility>
#include <vector>
#include <mfem.hpp>

namespace palace
{

class IoData;
class MaterialOperator;
class MaterialPropertyCoefficient;

namespace config
{

struct FarfieldBoundaryData;

}  // namespace config

enum class ProblemType : char;

//
// A class handling farfield, or absorbing, boundaries.
//
class FarfieldBoundaryOperator
{
private:
  // Reference to material property data (not owned).
  const MaterialOperator &mat_op;

  // List of all absorbing boundary condition attributes.
  mfem::Array<int> farfield_attr;

  // First- or second-order absorbing boundary condition.
  int order;

  // libCEED boundary attributes of the absorbing boundaries next to no frequency-dependent
  // material, and those next to each frequency-dependent material (index and attributes).
  // The boundary condition uses the neighboring material's permittivity, so on the latter
  // it depends on frequency through r(ω) = √(ε_w(iω)/ε∞), with the wave permittivity
  // ε_w = ε + σ/s of the material
  // (MaterialOperator::EvaluateFrequencyDependentWavePermittivity).
  mfem::Array<int> ceed_attr_other;
  std::vector<std::pair<std::size_t, mfem::Array<int>>> ceed_attr_dispersive;

  mfem::Array<int> SetUpBoundaryProperties(const config::FarfieldBoundaryData &farfield,
                                           ProblemType problem_type,
                                           const mfem::ParMesh &mesh);

  void AddDispersiveExtraSystemBdrCoefficients(std::complex<double> omega,
                                               MaterialPropertyCoefficient &dfbr,
                                               MaterialPropertyCoefficient &dfbi,
                                               MaterialPropertyCoefficient &fbr,
                                               MaterialPropertyCoefficient &fbi) const;

public:
  FarfieldBoundaryOperator(const config::FarfieldBoundaryData &farfield,
                           ProblemType problem_type, const MaterialOperator &mat_op,
                           const mfem::ParMesh &mesh);
  FarfieldBoundaryOperator(const IoData &iodata, const MaterialOperator &mat_op,
                           const mfem::ParMesh &mesh);

  // Returns array of farfield BC attributes.
  const auto &GetAttrList() const { return farfield_attr; }

  // Returns order of absorbing BC approximation.
  int GetOrder() const { return order; }

  // Add contributions to system matrices from first- or second-order absorbing boundary
  // condition. The damping and curl-curl coefficients cover the boundaries next to no
  // frequency-dependent material; the extra system coefficients include the rest.
  void AddDampingBdrCoefficients(double coeff, MaterialPropertyCoefficient &fb);
  void AddExtraSystemBdrCoefficients(double omega, MaterialPropertyCoefficient &dfbr,
                                     MaterialPropertyCoefficient &dfbi,
                                     MaterialPropertyCoefficient &fbr,
                                     MaterialPropertyCoefficient &fbi);

  // Complex-ω overload, used by the eigenmode nonlinear solve where the second-order
  // absorbing-BC coefficient i·(0.5/ω) is evaluated at the genuinely complex eigenvalue
  // (ω = -i·λ). For real ω (imag = 0) this produces the same (dfbr, dfbi) split as the
  // double overload (the whole term on dfbi).
  void AddExtraSystemBdrCoefficients(std::complex<double> omega,
                                     MaterialPropertyCoefficient &dfbr,
                                     MaterialPropertyCoefficient &dfbi,
                                     MaterialPropertyCoefficient &fbr,
                                     MaterialPropertyCoefficient &fbi);

  // Absorbing boundaries next to frequency-dependent material k (of
  // NumDispersiveMaterials): the first-order term is iω·r_k(ω) times the damping
  // coefficient, and the second-order term i·(0.5/ω)/r_k(ω) times the curl-curl
  // coefficient, both with the high-frequency permittivity.
  std::size_t NumDispersiveMaterials() const { return ceed_attr_dispersive.size(); }
  std::complex<double> GetDispersiveScale(std::size_t k, std::complex<double> omega) const;
  void AddDispersiveDampingBdrCoefficients(std::size_t k, double coeff,
                                           MaterialPropertyCoefficient &fb) const;
  void AddDispersiveCurlCurlBdrCoefficients(std::size_t k, double coeff,
                                            MaterialPropertyCoefficient &df) const;
  // The scalar version of the damping coefficient for the H1 block of the mode problem, on
  // each libCEED boundary attribute of material k.
  void AddDispersiveScalarDampingBdrCoefficients(std::size_t k, double coeff,
                                                 MaterialPropertyCoefficient &fb) const;
  const mfem::Array<int> &GetCeedAttrListOther() const { return ceed_attr_other; }

  // Add the ω-independent boundary curl-curl coefficient of the second-order absorbing
  // BC (the μ⁻¹·c₀ normal-projected term that the real-ω stamping scales by 0.5/ω),
  // scaled by `coeff`. Pulled out so both the real-ω and complex-ω overloads share the
  // material-coefficient construction. No-op when order < 2 or there are no farfield
  // attributes.
  void AddExtraSystemBoundaryCurlCurlBdrCoefficients(double coeff,
                                                     MaterialPropertyCoefficient &df) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_FARFIELD_BOUNDARY_OPERATOR_HPP
