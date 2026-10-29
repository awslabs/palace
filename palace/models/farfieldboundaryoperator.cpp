// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "farfieldboundaryoperator.hpp"

#include <array>
#include "linalg/densematrix.hpp"
#include "models/boundaryattributes.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"

namespace palace
{

FarfieldBoundaryOperator::FarfieldBoundaryOperator(const IoData &iodata,
                                                   const MaterialOperator &mat_op,
                                                   const mfem::ParMesh &mesh)
  : FarfieldBoundaryOperator(iodata.boundaries.farfield, iodata.problem.type, mat_op, mesh)
{
}

FarfieldBoundaryOperator::FarfieldBoundaryOperator(
    const config::FarfieldBoundaryData &farfield, ProblemType problem_type,
    const MaterialOperator &mat_op, const mfem::ParMesh &mesh)
  : mat_op(mat_op), farfield_attr(SetUpBoundaryProperties(farfield, problem_type, mesh))
{
  // Print out BC info for all farfield attributes.
  if (farfield_attr.Size())
  {
    Mpi::Print("\nConfiguring Robin absorbing BC (order {:d}) at attributes:\n", order);
    std::ranges::sort(farfield_attr);
    utils::PrettyPrint(farfield_attr);
  }
}

mfem::Array<int> FarfieldBoundaryOperator::SetUpBoundaryProperties(
    const config::FarfieldBoundaryData &farfield, ProblemType problem_type,
    const mfem::ParMesh &mesh)
{
  // Check that impedance boundary attributes have been specified correctly.
  MeshBoundaryAttributes mesh_attrs(mesh);
  if (!farfield.empty())
  {
    CheckBoundaryAttributes(mesh_attrs, std::array{farfield}, "absorbing", false);
  }

  // Set the order of the farfield boundary condition.
  order = farfield.order;

  // Mark selected boundary attributes from the mesh as farfield.
  mfem::Array<int> farfield_bcs;
  farfield_bcs.Reserve(static_cast<int>(farfield.attributes.size()));
  for (auto attr : farfield.attributes)
  {
    if (!mesh_attrs.Contains(attr))
    {
      continue;  // Can just ignore if wrong
    }
    farfield_bcs.Append(attr);
  }
  MFEM_VERIFY(farfield_bcs.Size() == 0 || order < 2 ||
                  problem_type == ProblemType::DRIVEN ||
                  problem_type == ProblemType::EIGENMODE ||
                  problem_type == ProblemType::BOUNDARYMODE,
              "Second-order farfield boundaries are only available for frequency "
              "domain simulations!");
  return farfield_bcs;
}

void FarfieldBoundaryOperator::AddDampingBdrCoefficients(double coeff,
                                                         MaterialPropertyCoefficient &fb)
{
  // First-order absorbing boundary condition.
  if (farfield_attr.Size())
  {
    MaterialPropertyCoefficient invz0_func(mat_op.GetBdrAttributeToMaterial(),
                                           mat_op.GetInvImpedance());
    invz0_func.RestrictCoefficient(mat_op.GetCeedBdrAttributes(farfield_attr));
    fb.AddCoefficient(invz0_func.GetAttributeToMaterial(),
                      invz0_func.GetMaterialProperties(), coeff);
  }
}

namespace
{

// Coefficient μ⁻¹c₀ of the second-order absorbing BC on the given boundary attributes.
MaterialPropertyCoefficient GetCurlCurlCoefficient(const MaterialOperator &mat_op,
                                                   const mfem::Array<int> &attr_list)
{
  mfem::DenseTensor muinvc0 =
      linalg::Mult(mat_op.GetInvPermeability(), mat_op.GetLightSpeed());
  MaterialPropertyCoefficient muinvc0_func(mat_op.GetBdrAttributeToMaterial(), muinvc0);
  muinvc0_func.RestrictCoefficient(mat_op.GetCeedBdrAttributes(attr_list));

  // Instead getting the correct normal of farfield boundary elements, just pick the
  // the first element normal. This is fine as long as the farfield material properties
  // are not anisotropic.
  mfem::Vector normal(mat_op.SpaceDimension());
  normal = 0.0;
  normal(0) = 1.0;
  muinvc0_func.NormalProjectedCoefficient(normal);
  return muinvc0_func;
}

// Material properties c M, for the scalar material properties c and the matrix M.
mfem::DenseTensor ScaledMaterialProperties(const mfem::DenseTensor &c,
                                           const mfem::DenseMatrix &M)
{
  MFEM_ASSERT(c.SizeI() == 1 && c.SizeJ() == 1,
              "Expected scalar material properties for the absorbing boundary!");
  mfem::DenseTensor mat_coeff(M.Height(), M.Width(), c.SizeK());
  for (int k = 0; k < c.SizeK(); k++)
  {
    mat_coeff(k).Set(c(0, 0, k), M);
  }
  return mat_coeff;
}

// Matrices k kᵀ and [k x] of the Floquet wave vector k, for the rotated boundary mass and
// mixed boundary curl terms of the second-order absorbing BC.
void GetFloquetMatrices(const MaterialOperator &mat_op, mfem::DenseMatrix &kkT,
                        mfem::DenseMatrix &kx)
{
  const auto &k = mat_op.GetWaveVector();
  MFEM_VERIFY(k.Size() == 3, "Floquet periodicity is only available in 3D!");
  kkT.SetSize(3);
  mfem::MultVVt(k, kkT);
  kx.SetSize(3);
  kx = 0.0;
  kx(0, 1) = -k(2);
  kx(0, 2) = k(1);
  kx(1, 0) = k(2);
  kx(1, 2) = -k(0);
  kx(2, 0) = -k(1);
  kx(2, 1) = k(0);
}

}  // namespace

void FarfieldBoundaryOperator::AddCurlCurlBdrCoefficients(
    std::complex<double> a, MaterialPropertyCoefficient &dfr,
    MaterialPropertyCoefficient &dfi, MaterialPropertyCoefficient &frr,
    MaterialPropertyCoefficient &fri, MaterialPropertyCoefficient *fpr,
    MaterialPropertyCoefficient *fpi) const
{
  // Contribution for second-order absorbing BC. See Jin Section 9.3 for reference. The β
  // coefficient for the second-order ABC is 1/(2ik+2/r). Taking the radius of curvature
  // as infinity (plane wave scattering), the r-dependence vanishes and the contribution
  // is purely imaginary: β = 1/(2ik) = -i/(2k). Multiplying through by μ⁻¹ we get the
  // material coefficient to ω as 1 / (μ √(με)), and the weak form term is
  // +β (μ⁻¹c₀ ∇ₛ×E, ∇ₛ×F), see docs/src/reference.md. Also, this implementation ignores the
  // divergence term ∇⋅Eₜ, as COMSOL does as well. The ω-dependence enters only through the
  // scalar a, so the real-ω, complex-ω, and ω-independent (M_ff) matrices share the
  // material coefficient construction: the sign of β is in a = β = 0.5 / (iω) for the
  // system matrix (AddExtraSystemBdrCoefficients), and in a = -1 for the matrix M_ff, whose
  // term is i·(0.5/ω)·M_ff (SpaceOperator::GetFarfieldBoundaryCurlCurlMatrix).
  if (!farfield_attr.Size() || order < 2)
  {
    return;
  }
  const auto muinvc0_func = GetCurlCurlCoefficient(mat_op, farfield_attr);

  // a (C + T + i (X - Xᵀ)) = a_r (C + T) - a_i (X - Xᵀ) + i (a_i (C + T) + a_r (X - Xᵀ)),
  // for the curl-curl term C, the rotated boundary mass term T, and the mixed boundary curl
  // term X.
  const auto &attr_mat = muinvc0_func.GetAttributeToMaterial();
  const auto &mat_c = muinvc0_func.GetMaterialProperties();
  if (a.real() != 0.0)
  {
    dfr.AddCoefficient(attr_mat, mat_c, a.real());
  }
  if (a.imag() != 0.0)
  {
    dfi.AddCoefficient(attr_mat, mat_c, a.imag());
  }

  // Floquet periodicity, E = E_p exp(-i k·x): the surface curl of E is
  // exp(-i k·x) (n·∇×E_p - i n·(k x E_p)), where n·(k x E_p) = -k·(n x E_p) only depends on
  // the tangential trace of E_p. With a frequency dependent wave vector, these terms do not
  // scale with 1/ω and are added to the damping and stiffness matrices instead.
  if (!mat_op.HasWaveVector() || mat_op.HasFloquetFrequencyScaling())
  {
    return;
  }
  mfem::DenseMatrix kkT, kx;
  GetFloquetMatrices(mat_op, kkT, kx);
  const auto mat_ckkT = ScaledMaterialProperties(mat_c, kkT);
  const auto mat_ckx = ScaledMaterialProperties(mat_c, kx);
  if (a.real() != 0.0)
  {
    frr.AddCoefficient(attr_mat, mat_ckkT, a.real());
    if (fpi)
    {
      fpi->AddCoefficient(attr_mat, mat_ckx, a.real());
    }
  }
  if (a.imag() != 0.0)
  {
    fri.AddCoefficient(attr_mat, mat_ckkT, a.imag());
    if (fpr)
    {
      fpr->AddCoefficient(attr_mat, mat_ckx, -a.imag());
    }
  }
}

void FarfieldBoundaryOperator::AddExtraSystemBdrCoefficients(
    double omega, MaterialPropertyCoefficient &dfbr, MaterialPropertyCoefficient &dfbi,
    MaterialPropertyCoefficient &frbr, MaterialPropertyCoefficient &frbi,
    MaterialPropertyCoefficient *fpbr, MaterialPropertyCoefficient *fpbi)
{
  // Real-ω second-order ABC: β = 1/(2iω) = -i·(0.5/ω) is imaginary, so the curl-curl and
  // rotated boundary mass terms land on the imaginary part and the mixed boundary curl
  // terms on the real part.
  AddCurlCurlBdrCoefficients(0.5 / (std::complex<double>(0.0, 1.0) * omega), dfbr, dfbi,
                             frbr, frbi, fpbr, fpbi);
}

void FarfieldBoundaryOperator::AddExtraSystemBdrCoefficients(
    std::complex<double> omega, MaterialPropertyCoefficient &dfbr,
    MaterialPropertyCoefficient &dfbi, MaterialPropertyCoefficient &frbr,
    MaterialPropertyCoefficient &frbi, MaterialPropertyCoefficient *fpbr,
    MaterialPropertyCoefficient *fpbi)
{
  // Complex-ω second-order ABC for the eigenmode nonlinear solve: with ω complex,
  // β = 1/(2iω) has a real and an imaginary part, so every term contributes to both parts.
  // For real ω (imag = 0) this reduces to the double overload.
  AddCurlCurlBdrCoefficients(0.5 / (std::complex<double>(0.0, 1.0) * omega), dfbr, dfbi,
                             frbr, frbi, fpbr, fpbi);
}

void FarfieldBoundaryOperator::AddDampingFloquetBdrCoefficients(
    double coeff, MaterialPropertyCoefficient &frb)
{
  // With k = ω k₀, the rotated boundary mass term of the second-order ABC is
  // β ω² (μ⁻¹c₀ k₀ k₀ᵀ n x E, n x F) = iω (-0.5) (μ⁻¹c₀ k₀ k₀ᵀ n x E, n x F).
  if (!farfield_attr.Size() || order < 2 || !mat_op.HasWaveVector() ||
      !mat_op.HasFloquetFrequencyScaling())
  {
    return;
  }
  mfem::DenseMatrix kkT, kx;
  GetFloquetMatrices(mat_op, kkT, kx);
  const auto muinvc0_func = GetCurlCurlCoefficient(mat_op, farfield_attr);
  frb.AddCoefficient(muinvc0_func.GetAttributeToMaterial(),
                     ScaledMaterialProperties(muinvc0_func.GetMaterialProperties(), kkT),
                     -0.5 * coeff);
}

void FarfieldBoundaryOperator::AddStiffnessFloquetBdrCoefficients(
    double coeff, MaterialPropertyCoefficient &fpb)
{
  // With k = ω k₀, the mixed boundary curl terms of the second-order ABC are
  // iβ ω [(μ⁻¹c₀ [k₀ x]ᵀ ∇×E, F) - (μ⁻¹c₀ [k₀ x] E, ∇×F)], with iβ ω = 0.5.
  if (!farfield_attr.Size() || order < 2 || !mat_op.HasWaveVector() ||
      !mat_op.HasFloquetFrequencyScaling())
  {
    return;
  }
  mfem::DenseMatrix kkT, kx;
  GetFloquetMatrices(mat_op, kkT, kx);
  const auto muinvc0_func = GetCurlCurlCoefficient(mat_op, farfield_attr);
  fpb.AddCoefficient(muinvc0_func.GetAttributeToMaterial(),
                     ScaledMaterialProperties(muinvc0_func.GetMaterialProperties(), kx),
                     0.5 * coeff);
}

}  // namespace palace
