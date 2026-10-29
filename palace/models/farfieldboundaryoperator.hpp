// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_FARFIELD_BOUNDARY_OPERATOR_HPP
#define PALACE_MODELS_FARFIELD_BOUNDARY_OPERATOR_HPP

#include <complex>
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

  mfem::Array<int> SetUpBoundaryProperties(const config::FarfieldBoundaryData &farfield,
                                           ProblemType problem_type,
                                           const mfem::ParMesh &mesh);

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
  // condition.
  void AddDampingBdrCoefficients(double coeff, MaterialPropertyCoefficient &fb);

  // Add the second-order absorbing BC term at the frequency ω (see
  // AddCurlCurlBdrCoefficients with a = β = 1/(2iω)) to the real and imaginary parts of the
  // system matrix. Without a Floquet wave vector, only the boundary curl-curl coefficients
  // dfbr, dfbi are modified.
  void AddExtraSystemBdrCoefficients(double omega, MaterialPropertyCoefficient &dfbr,
                                     MaterialPropertyCoefficient &dfbi,
                                     MaterialPropertyCoefficient &frbr,
                                     MaterialPropertyCoefficient &frbi,
                                     MaterialPropertyCoefficient *fpbr,
                                     MaterialPropertyCoefficient *fpbi);

  // Complex-ω overload, used by the eigenmode nonlinear solve where the second-order
  // absorbing-BC coefficient β = 1/(2iω) is evaluated at the genuinely complex eigenvalue
  // (ω = -i·λ). For real ω (imag = 0) this produces the same split as the double overload.
  void AddExtraSystemBdrCoefficients(std::complex<double> omega,
                                     MaterialPropertyCoefficient &dfbr,
                                     MaterialPropertyCoefficient &dfbi,
                                     MaterialPropertyCoefficient &frbr,
                                     MaterialPropertyCoefficient &frbi,
                                     MaterialPropertyCoefficient *fpbr,
                                     MaterialPropertyCoefficient *fpbi);

  // With a Floquet wave vector k = ω k₀ which depends on frequency
  // (FloquetReferenceFrequency), the Floquet terms of the second-order absorbing BC (see
  // AddCurlCurlBdrCoefficients) scale as iω (rotated boundary mass) and 1 (mixed boundary
  // curl), and are added with coefficient coeff to the damping and stiffness matrices
  // instead of the extra system matrix. No-op without a frequency dependent Floquet wave
  // vector.
  void AddDampingFloquetBdrCoefficients(double coeff, MaterialPropertyCoefficient &frb);
  void AddStiffnessFloquetBdrCoefficients(double coeff, MaterialPropertyCoefficient &fpb);

  // Add a·(μ⁻¹c₀ n·∇ₖ×u, n·∇₋ₖ×v) on the absorbing boundaries for the complex scalar a,
  // with μ⁻¹c₀ normal projected. For a = β = 1/(2iω), this is the second-order absorbing BC
  // term (see docs/src/reference.md). With the Floquet wave vector k, the surface curl of
  // the field u exp(-i k·x) is exp(-i k·x) n·∇ₖ×u, where n·∇ₖ×u = n·∇×u - i n·(k x u),
  // which splits the Hermitian form into
  //   (μ⁻¹c₀ n·∇ₖ×u, n·∇₋ₖ×v) = (μ⁻¹c₀ n·∇×u, n·∇×v) + (μ⁻¹c₀ k kᵀ n x u, n x v)
  //       + i [(μ⁻¹c₀ [k x]ᵀ ∇×u, v) - (μ⁻¹c₀ [k x] u, ∇×v)],
  // with the boundary curl ∇×u = (n·∇×u) n of the tangential trace. The boundary curl-curl
  // term is added to dfr, dfi, the rotated boundary mass term
  // (RotatedVectorFEMassIntegrator) to frr, fri, and the mixed boundary curl terms
  // (transposed MixedVectorCurlIntegrator and MixedVectorWeakCurlIntegrator) to fpr, fpi,
  // unless null. The Floquet terms are only included with a Floquet wave vector which does
  // not depend on frequency (see AddDampingFloquetBdrCoefficients). No-op when order < 2 or
  // there are no farfield attributes.
  void AddCurlCurlBdrCoefficients(std::complex<double> a, MaterialPropertyCoefficient &dfr,
                                  MaterialPropertyCoefficient &dfi,
                                  MaterialPropertyCoefficient &frr,
                                  MaterialPropertyCoefficient &fri,
                                  MaterialPropertyCoefficient *fpr,
                                  MaterialPropertyCoefficient *fpi) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_FARFIELD_BOUNDARY_OPERATOR_HPP
