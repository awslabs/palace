// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "farfieldboundaryoperator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <string>
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
  CheckFloquetPassivity(iodata.solver, iodata.problem.type, iodata.units, mesh);
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
  SetUpCurlCurlCoefficients();
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

}  // namespace

void FarfieldBoundaryOperator::SetUpCurlCurlCoefficients()
{
  if (!farfield_attr.Size() || order < 2)
  {
    return;
  }
  const auto muinvc0_func = GetCurlCurlCoefficient(mat_op, farfield_attr);
  curlcurl_attr_mat = muinvc0_func.GetAttributeToMaterial();
  curlcurl_mat = muinvc0_func.GetMaterialProperties();

  // Matrices k kᵀ and [k x] of the Floquet wave vector k (k₀ with frequency scaling), for
  // the rotated boundary mass and mixed boundary curl terms. The cross product matrix is
  // the one of the domain Floquet terms.
  if (!mat_op.HasWaveVector())
  {
    return;
  }
  const auto &k = mat_op.GetWaveVector();
  MFEM_VERIFY(k.Size() == 3, "Floquet periodicity is only available in 3D!");
  mfem::DenseMatrix kkT(3);
  mfem::MultVVt(k, kkT);
  curlcurl_kkT_mat = ScaledMaterialProperties(curlcurl_mat, kkT);
  curlcurl_kx_mat = ScaledMaterialProperties(curlcurl_mat, mat_op.GetWaveVectorCross());
}

void FarfieldBoundaryOperator::CheckFloquetPassivity(const config::SolverData &solver,
                                                     ProblemType problem_type,
                                                     const Units &units,
                                                     const mfem::ParMesh &mesh) const
{
  if (!farfield_attr.Size() || order < 2 || !mat_op.HasWaveVector() ||
      mesh.SpaceDimension() != 3)
  {
    return;
  }

  // Lowest angular frequency of the simulation, where the wavenumber k is smallest.
  double omega_min = 0.0;
  if (problem_type == ProblemType::DRIVEN && !solver.driven.sample_f.empty())
  {
    omega_min = *std::ranges::min_element(solver.driven.sample_f);
  }
  else if (problem_type == ProblemType::EIGENMODE)
  {
    omega_min = solver.eigenmode.target;
  }
  const bool scaling = mat_op.HasFloquetFrequencyScaling();
  if (!scaling && omega_min <= 0.0)
  {
    return;
  }

  // Largest tangential part of the wave vector on the absorbing boundary elements, and
  // largest light speed of the materials next to them (with frequency scaling, the wave
  // vector is k₀ = k / ω).
  const auto &k = mat_op.GetWaveVector();
  const auto marker = mesh::AttrToMarker(mesh.bdr_attributes.Max(), farfield_attr);
  mfem::IsoparametricTransformation T;
  mfem::Vector normal(3);
  double kt_c_max[2] = {0.0, 0.0};  // |k_t|, c
  for (int i = 0; i < mesh.GetNBE(); i++)
  {
    if (!marker[mesh.GetBdrAttribute(i) - 1])
    {
      continue;
    }
    mesh.GetBdrElementTransformation(i, &T);
    T.SetIntPoint(&mfem::Geometries.GetCenter(T.GetGeometryType()));
    mfem::CalcOrtho(T.Jacobian(), normal);
    const double k_n = (k * normal) / normal.Norml2();
    kt_c_max[0] = std::max(kt_c_max[0], std::sqrt(std::max(k * k - k_n * k_n, 0.0)));
  }
  const auto bdr_attr_mat = mat_op.GetBdrAttributeToMaterial();
  for (auto attr : mat_op.GetCeedBdrAttributes(farfield_attr))
  {
    if (const int mat = bdr_attr_mat[attr - 1]; mat >= 0)
    {
      kt_c_max[1] = std::max(kt_c_max[1], mat_op.GetLightSpeedMax()[mat]);
    }
  }
  Mpi::GlobalMax(2, kt_c_max, mesh.GetComm());

  // k_t / k, with k = ω / c.
  const double ratio = kt_c_max[0] * kt_c_max[1] / (scaling ? 1.0 : omega_min);
  if (ratio > std::sqrt(2.0))
  {
    const std::string where =
        scaling ? std::string("at every frequency")
                : fmt::format("at {:.3e} GHz, the lowest frequency of the simulation",
                              units.Dimensionalize<Units::ValueType::FREQUENCY>(omega_min) /
                                  (2.0 * std::numbers::pi));
    Mpi::Warning(
        "The second-order absorbing boundary condition is not passive with this Floquet "
        "wave vector: its part tangential to the absorbing boundary is {:.3f} times the "
        "wavenumber k in the material at the boundary {}, above √2. For TE fields, the "
        "condition imposes the effective wavenumber k - k_t²/(2k), which is then "
        "negative, and the boundary injects energy (for example, into growing eigenmodes). "
        "The first-order absorbing boundary condition is passive.\n",
        ratio, where);
  }
}

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

  // a (C + T + i (X - Xᵀ)) = a_r (C + T) - a_i (X - Xᵀ) + i (a_i (C + T) + a_r (X - Xᵀ)),
  // for the curl-curl term C, the rotated boundary mass term T, and the mixed boundary curl
  // term X.
  const auto &attr_mat = curlcurl_attr_mat;
  if (a.real() != 0.0)
  {
    dfr.AddCoefficient(attr_mat, curlcurl_mat, a.real());
  }
  if (a.imag() != 0.0)
  {
    dfi.AddCoefficient(attr_mat, curlcurl_mat, a.imag());
  }

  // Floquet periodicity, E = E_p exp(-i k·x): the surface curl of E is
  // exp(-i k·x) (n·∇×E_p - i n·(k x E_p)), where n·(k x E_p) = -k·(n x E_p) only depends on
  // the tangential trace of E_p. With a frequency dependent wave vector, these terms do not
  // scale with 1/ω and are added to the damping and stiffness matrices instead.
  if (!mat_op.HasWaveVector() || mat_op.HasFloquetFrequencyScaling())
  {
    return;
  }
  if (a.real() != 0.0)
  {
    frr.AddCoefficient(attr_mat, curlcurl_kkT_mat, a.real());
    if (fpi)
    {
      fpi->AddCoefficient(attr_mat, curlcurl_kx_mat, a.real());
    }
  }
  if (a.imag() != 0.0)
  {
    fri.AddCoefficient(attr_mat, curlcurl_kkT_mat, a.imag());
    if (fpr)
    {
      fpr->AddCoefficient(attr_mat, curlcurl_kx_mat, -a.imag());
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
  frb.AddCoefficient(curlcurl_attr_mat, curlcurl_kkT_mat, -0.5 * coeff);
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
  fpb.AddCoefficient(curlcurl_attr_mat, curlcurl_kx_mat, 0.5 * coeff);
}

}  // namespace palace
