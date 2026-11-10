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
  {
    auto split =
        mat_op.SplitCeedBdrAttributesByFrequencyDependentPermittivity(farfield_attr);
    ceed_attr_other = std::move(split.other);
    ceed_attr_dispersive = std::move(split.materials);
  }
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
    invz0_func.RestrictCoefficient(ceed_attr_other);
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
  AddCurlCurlBdrCoefficients(ceed_attr_other, a,
                             mat_op.HasFloquetFrequencyScaling() ? 0.0 : 1.0, dfr, dfi, frr,
                             fri, fpr, fpi);
}

void FarfieldBoundaryOperator::AddCurlCurlBdrCoefficients(
    const mfem::Array<int> &attrs, std::complex<double> a,
    std::complex<double> wave_vector_scale, MaterialPropertyCoefficient &dfr,
    MaterialPropertyCoefficient &dfi, MaterialPropertyCoefficient &frr,
    MaterialPropertyCoefficient &fri, MaterialPropertyCoefficient *fpr,
    MaterialPropertyCoefficient *fpi) const
{
  // β = 1/(2ik) = -i/(2k). The form C + T + i(X - Xᵀ) uses the surface curl of
  // E exp(-i k·x), as in docs/src/reference.md. Here k = wave_vector_scale k₀ for
  // frequency scaling, or k₀ otherwise. The sign of β is in a for system stamping,
  // and in a = -1 for M_ff, whose scalar is i·0.5/ω.
  if (attrs.Size() == 0 || order < 2)
  {
    return;
  }
  auto add =
      [&](MaterialPropertyCoefficient &f, const mfem::DenseTensor &tensor, double coeff)
  {
    if (coeff != 0.0)
    {
      MaterialPropertyCoefficient c(curlcurl_attr_mat, tensor);
      c.RestrictCoefficient(attrs);
      f.AddCoefficient(c.GetAttributeToMaterial(), c.GetMaterialProperties(), coeff);
    }
  };
  add(dfr, curlcurl_mat, a.real());
  add(dfi, curlcurl_mat, a.imag());
  if (!mat_op.HasWaveVector() || wave_vector_scale == 0.0)
  {
    return;
  }
  const auto rotated = a * wave_vector_scale * wave_vector_scale;
  const auto mixed = std::complex<double>(0.0, 1.0) * a * wave_vector_scale;
  add(frr, curlcurl_kkT_mat, rotated.real());
  add(fri, curlcurl_kkT_mat, rotated.imag());
  if (fpr)
  {
    add(*fpr, curlcurl_kx_mat, mixed.real());
  }
  if (fpi)
  {
    add(*fpi, curlcurl_kx_mat, mixed.imag());
  }
}

void FarfieldBoundaryOperator::AddExtraSystemBdrCoefficients(
    double omega, MaterialPropertyCoefficient &dfbr, MaterialPropertyCoefficient &dfbi,
    MaterialPropertyCoefficient &fbr, MaterialPropertyCoefficient &fbi,
    MaterialPropertyCoefficient &frbr, MaterialPropertyCoefficient &frbi,
    MaterialPropertyCoefficient *fpbr, MaterialPropertyCoefficient *fpbi)
{
  AddExtraSystemBdrCoefficients(std::complex<double>(omega, 0.0), dfbr, dfbi, fbr, fbi,
                                frbr, frbi, fpbr, fpbi);
}

void FarfieldBoundaryOperator::AddExtraSystemBdrCoefficients(
    std::complex<double> omega, MaterialPropertyCoefficient &dfbr,
    MaterialPropertyCoefficient &dfbi, MaterialPropertyCoefficient &fbr,
    MaterialPropertyCoefficient &fbi, MaterialPropertyCoefficient &frbr,
    MaterialPropertyCoefficient &frbi, MaterialPropertyCoefficient *fpbr,
    MaterialPropertyCoefficient *fpbi)
{
  constexpr std::complex<double> i(0.0, 1.0);
  AddCurlCurlBdrCoefficients(0.5 / (i * omega), dfbr, dfbi, frbr, frbi, fpbr, fpbi);
  for (std::size_t k = 0; k < ceed_attr_dispersive.size(); k++)
  {
    const auto r = GetDispersiveScale(k, omega);
    const auto z1 = i * omega * r;
    AddDispersiveDampingBdrCoefficients(k, z1.real(), fbr);
    AddDispersiveDampingBdrCoefficients(k, z1.imag(), fbi);
    // All three second-order terms carry 1/r. With k = ω k₀, none belongs to the
    // frequency-independent damping or stiffness matrix on these attributes.
    AddCurlCurlBdrCoefficients(ceed_attr_dispersive[k].second, 0.5 / (i * omega * r),
                               mat_op.HasFloquetFrequencyScaling() ? omega : 1.0, dfbr,
                               dfbi, frbr, frbi, fpbr, fpbi);
  }
}

void FarfieldBoundaryOperator::AddPreconditionerBdrCoefficients(
    std::complex<double> omega, MaterialPropertyCoefficient &df,
    MaterialPropertyCoefficient &f, MaterialPropertyCoefficient &fr) const
{
  // As in the upstream real-valued PC, make only the boundary curl-curl term
  // positive semidefinite and drop the antisymmetric mixed terms. The rotated
  // Floquet trace keeps its physical sign. Dispersive scalars use magnitudes.
  const auto a = 0.5 / (std::complex<double>(0.0, 1.0) * omega);
  MaterialPropertyCoefficient c(curlcurl_attr_mat, curlcurl_mat);
  if (order > 1 && ceed_attr_other.Size())
  {
    c.RestrictCoefficient(ceed_attr_other);
    df.AddCoefficient(c.GetAttributeToMaterial(), c.GetMaterialProperties(),
                      -(a.real() + a.imag()));
  }
  if (order > 1 && ceed_attr_other.Size() && mat_op.HasWaveVector() &&
      !mat_op.HasFloquetFrequencyScaling())
  {
    MaterialPropertyCoefficient t(curlcurl_attr_mat, curlcurl_kkT_mat);
    t.RestrictCoefficient(ceed_attr_other);
    fr.AddCoefficient(t.GetAttributeToMaterial(), t.GetMaterialProperties(),
                      a.real() + a.imag());
  }
  constexpr std::complex<double> i(0.0, 1.0);
  for (std::size_t k = 0; k < ceed_attr_dispersive.size(); k++)
  {
    const auto r = GetDispersiveScale(k, omega);
    AddDispersiveDampingBdrCoefficients(k, std::abs(i * omega * r), f);
    AddDispersiveCurlCurlBdrCoefficients(k, std::abs(0.5 / (omega * r)), df);
    if (order > 1 && mat_op.HasWaveVector())
    {
      MaterialPropertyCoefficient t(curlcurl_attr_mat, curlcurl_kkT_mat);
      t.RestrictCoefficient(ceed_attr_dispersive[k].second);
      const auto scale = mat_op.HasFloquetFrequencyScaling() ? omega : 1.0;
      fr.AddCoefficient(t.GetAttributeToMaterial(), t.GetMaterialProperties(),
                        -std::abs(0.5 * scale * scale / (omega * r)));
    }
  }
}

std::complex<double>
FarfieldBoundaryOperator::GetDispersiveScale(std::size_t k,
                                             std::complex<double> omega) const
{
  // Take the passive lossy limit on the negative-real branch cut: the outgoing root
  // decays below the plasma frequency. User conductivity is included once in ε_wave.
  const auto m = ceed_attr_dispersive.at(k).first;
  auto eps = mat_op.EvaluateFrequencyDependentWavePermittivity(
                 m, std::complex<double>(0.0, 1.0) * omega) /
             mat_op.GetFrequencyDependentPermittivityHighFrequency(m);
  if (eps.imag() == 0.0)
  {
    eps.imag(-0.0);
  }
  return std::sqrt(eps);
}

void FarfieldBoundaryOperator::AddDispersiveDampingBdrCoefficients(
    std::size_t k, double coeff, MaterialPropertyCoefficient &fb) const
{
  const auto &attrs = ceed_attr_dispersive.at(k).second;
  if (coeff != 0.0 && attrs.Size())
  {
    MaterialPropertyCoefficient invz0_func(mat_op.GetBdrAttributeToMaterial(),
                                           mat_op.GetInvImpedance());
    invz0_func.RestrictCoefficient(attrs);
    fb.AddCoefficient(invz0_func.GetAttributeToMaterial(),
                      invz0_func.GetMaterialProperties(), coeff);
  }
}

void FarfieldBoundaryOperator::AddDispersiveCurlCurlBdrCoefficients(
    std::size_t k, double coeff, MaterialPropertyCoefficient &df) const
{
  const auto &attrs = ceed_attr_dispersive.at(k).second;
  if (coeff != 0.0 && order > 1 && attrs.Size())
  {
    MaterialPropertyCoefficient c(curlcurl_attr_mat, curlcurl_mat);
    c.RestrictCoefficient(attrs);
    df.AddCoefficient(c.GetAttributeToMaterial(), c.GetMaterialProperties(), coeff);
  }
}

void FarfieldBoundaryOperator::AddDispersiveCurlCurlBdrCoefficients(
    std::size_t k, std::complex<double> a, MaterialPropertyCoefficient &dfr,
    MaterialPropertyCoefficient &dfi, MaterialPropertyCoefficient &frr,
    MaterialPropertyCoefficient &fri, MaterialPropertyCoefficient *fpr,
    MaterialPropertyCoefficient *fpi) const
{
  AddCurlCurlBdrCoefficients(ceed_attr_dispersive.at(k).second, a,
                             mat_op.HasFloquetFrequencyScaling() ? 0.0 : 1.0, dfr, dfi, frr,
                             fri, fpr, fpi);
}

void FarfieldBoundaryOperator::AddDispersiveScalarDampingBdrCoefficients(
    std::size_t k, double coeff, MaterialPropertyCoefficient &fb) const
{
  const auto &inv_z = mat_op.GetInvImpedance();
  const auto bdr_attr_to_mat = mat_op.GetBdrAttributeToMaterial();
  for (auto ceed_attr : ceed_attr_dispersive.at(k).second)
  {
    const int mat_idx = bdr_attr_to_mat[ceed_attr - 1];
    MFEM_ASSERT(mat_idx >= 0, "Missing material for farfield boundary attribute!");
    fb.AddMaterialProperty(ceed_attr, inv_z(0, 0, mat_idx), coeff);
  }
}

void FarfieldBoundaryOperator::AddDispersiveDampingFloquetBdrCoefficients(
    std::size_t k, double coeff, MaterialPropertyCoefficient &fr) const
{
  if (order < 2 || !mat_op.HasWaveVector() || !mat_op.HasFloquetFrequencyScaling())
  {
    return;
  }
  MaterialPropertyCoefficient t(curlcurl_attr_mat, curlcurl_kkT_mat);
  t.RestrictCoefficient(ceed_attr_dispersive.at(k).second);
  fr.AddCoefficient(t.GetAttributeToMaterial(), t.GetMaterialProperties(), -0.5 * coeff);
}

void FarfieldBoundaryOperator::AddDispersiveStiffnessFloquetBdrCoefficients(
    std::size_t k, double coeff, MaterialPropertyCoefficient &fp) const
{
  if (order < 2 || !mat_op.HasWaveVector() || !mat_op.HasFloquetFrequencyScaling())
  {
    return;
  }
  MaterialPropertyCoefficient x(curlcurl_attr_mat, curlcurl_kx_mat);
  x.RestrictCoefficient(ceed_attr_dispersive.at(k).second);
  fp.AddCoefficient(x.GetAttributeToMaterial(), x.GetMaterialProperties(), 0.5 * coeff);
}

void FarfieldBoundaryOperator::AddDampingFloquetBdrCoefficients(
    double coeff, MaterialPropertyCoefficient &frb)
{
  // The rotated trace scales as iω when k = ω k₀, on nondispersive attributes only.
  if (!farfield_attr.Size() || order < 2 || !mat_op.HasWaveVector() ||
      !mat_op.HasFloquetFrequencyScaling())
  {
    return;
  }
  MaterialPropertyCoefficient t(curlcurl_attr_mat, curlcurl_kkT_mat);
  t.RestrictCoefficient(ceed_attr_other);
  frb.AddCoefficient(t.GetAttributeToMaterial(), t.GetMaterialProperties(), -0.5 * coeff);
}

void FarfieldBoundaryOperator::AddStiffnessFloquetBdrCoefficients(
    double coeff, MaterialPropertyCoefficient &fpb)
{
  // The mixed boundary curl scales as 1 when k = ω k₀, on nondispersive attributes only.
  if (!farfield_attr.Size() || order < 2 || !mat_op.HasWaveVector() ||
      !mat_op.HasFloquetFrequencyScaling())
  {
    return;
  }
  MaterialPropertyCoefficient x(curlcurl_attr_mat, curlcurl_kx_mat);
  x.RestrictCoefficient(ceed_attr_other);
  fpb.AddCoefficient(x.GetAttributeToMaterial(), x.GetMaterialProperties(), 0.5 * coeff);
}

}  // namespace palace
