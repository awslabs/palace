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

void FarfieldBoundaryOperator::AddExtraSystemBoundaryCurlCurlBdrCoefficients(
    double coeff, MaterialPropertyCoefficient &df) const
{
  // Contribution for second-order absorbing BC. See Jin Section 9.3 for reference. The β
  // coefficient for the second-order ABC is 1/(2ik+2/r). Taking the radius of curvature
  // as infinity (plane wave scattering), the r-dependence vanishes and the contribution
  // is purely imaginary. Multiplying through by μ⁻¹ we get the material coefficient to ω
  // as 1 / (μ √(με)). Also, this implementation ignores the divergence term ∇⋅Eₜ, as
  // COMSOL does as well. The ω-dependence enters only through the scalar `coeff` (= 0.5/ω
  // for the imaginary slot in the real-ω stamping), so the material coefficient is pulled
  // out here and reused by both the real-ω and complex-ω overloads.
  if (farfield_attr.Size() && order > 1)
  {
    mfem::DenseTensor muinvc0 =
        linalg::Mult(mat_op.GetInvPermeability(), mat_op.GetLightSpeed());
    MaterialPropertyCoefficient muinvc0_func(mat_op.GetBdrAttributeToMaterial(), muinvc0);
    muinvc0_func.RestrictCoefficient(ceed_attr_other);

    // Instead getting the correct normal of farfield boundary elements, just pick the
    // the first element normal. This is fine as long as the farfield material properties
    // are not anisotropic.
    mfem::Vector normal(mat_op.SpaceDimension());
    normal = 0.0;
    normal(0) = 1.0;
    muinvc0_func.NormalProjectedCoefficient(normal);

    df.AddCoefficient(muinvc0_func.GetAttributeToMaterial(),
                      muinvc0_func.GetMaterialProperties(), coeff);
  }
}

void FarfieldBoundaryOperator::AddExtraSystemBdrCoefficients(
    double omega, MaterialPropertyCoefficient &dfbr, MaterialPropertyCoefficient &dfbi,
    MaterialPropertyCoefficient &fbr, MaterialPropertyCoefficient &fbi)
{
  // Real-ω second-order ABC: the term is i·(0.5/ω)·(curl-curl). At real ω the scalar
  // 0.5/ω is real, so the whole contribution is imaginary and lands on dfbi.
  AddExtraSystemBoundaryCurlCurlBdrCoefficients(0.5 / omega, dfbi);
  AddDispersiveExtraSystemBdrCoefficients(omega, dfbr, dfbi, fbr, fbi);
}

void FarfieldBoundaryOperator::AddExtraSystemBdrCoefficients(
    std::complex<double> omega, MaterialPropertyCoefficient &dfbr,
    MaterialPropertyCoefficient &dfbi, MaterialPropertyCoefficient &fbr,
    MaterialPropertyCoefficient &fbi)
{
  // Complex-ω second-order ABC for the eigenmode nonlinear solve. The contribution is
  // i·(0.5/ω)·(curl-curl); with ω complex, s = 0.5/ω is complex and i·s splits as
  // i·(s.real() + i·s.imag()) = -s.imag() + i·s.real(), so the real part -s.imag() goes
  // on dfbr and the imaginary part s.real() on dfbi. For real ω (s.imag() == 0) this
  // reduces to the double overload (whole term on dfbi).
  const std::complex<double> s = 0.5 / omega;
  AddExtraSystemBoundaryCurlCurlBdrCoefficients(s.real(), dfbi);
  if (s.imag() != 0.0)
  {
    AddExtraSystemBoundaryCurlCurlBdrCoefficients(-s.imag(), dfbr);
  }
  AddDispersiveExtraSystemBdrCoefficients(omega, dfbr, dfbi, fbr, fbi);
}

void FarfieldBoundaryOperator::AddDispersiveExtraSystemBdrCoefficients(
    std::complex<double> omega, MaterialPropertyCoefficient &dfbr,
    MaterialPropertyCoefficient &dfbi, MaterialPropertyCoefficient &fbr,
    MaterialPropertyCoefficient &fbi) const
{
  // The first-order term iω·r(ω) and second-order term i·(0.5/ω)/r(ω), split into their
  // real and imaginary parts.
  constexpr std::complex<double> i(0.0, 1.0);
  for (std::size_t k = 0; k < ceed_attr_dispersive.size(); k++)
  {
    const auto r = GetDispersiveScale(k, omega);
    const auto z1 = i * omega * r;
    AddDispersiveDampingBdrCoefficients(k, z1.real(), fbr);
    AddDispersiveDampingBdrCoefficients(k, z1.imag(), fbi);
    if (order > 1)
    {
      const auto z2 = i * 0.5 / (omega * r);
      AddDispersiveCurlCurlBdrCoefficients(k, z2.real(), dfbr);
      AddDispersiveCurlCurlBdrCoefficients(k, z2.imag(), dfbi);
    }
  }
}

std::complex<double>
FarfieldBoundaryOperator::GetDispersiveScale(std::size_t k,
                                             std::complex<double> omega) const
{
  // The principal root, with a nonnegative real part, is the outgoing wave for a passive
  // material (Im ε ≤ 0). A lossless ε < 0 (Drude below the plasma frequency) lies on the
  // branch cut: take the lossy limit Im ε → 0⁻, the decaying wave.
  const auto m = ceed_attr_dispersive.at(k).first;
  // The wave permittivity, with the user conductivity of the material.
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
    mfem::DenseTensor muinvc0 =
        linalg::Mult(mat_op.GetInvPermeability(), mat_op.GetLightSpeed());
    MaterialPropertyCoefficient muinvc0_func(mat_op.GetBdrAttributeToMaterial(), muinvc0);
    muinvc0_func.RestrictCoefficient(attrs);
    mfem::Vector normal(mat_op.SpaceDimension());
    normal = 0.0;
    normal(0) = 1.0;
    muinvc0_func.NormalProjectedCoefficient(normal);
    df.AddCoefficient(muinvc0_func.GetAttributeToMaterial(),
                      muinvc0_func.GetMaterialProperties(), coeff);
  }
}

void FarfieldBoundaryOperator::AddDispersiveScalarDampingBdrCoefficients(
    std::size_t k, double coeff, MaterialPropertyCoefficient &fb) const
{
  // Scalar inverse impedance of the material adjacent to each boundary element, for the H1
  // block of the mode problem (assumes an isotropic material at the boundary).
  const auto &inv_z = mat_op.GetInvImpedance();
  const auto bdr_attr_to_mat = mat_op.GetBdrAttributeToMaterial();
  for (auto ceed_attr : ceed_attr_dispersive.at(k).second)
  {
    const int mat_idx = bdr_attr_to_mat[ceed_attr - 1];
    MFEM_ASSERT(mat_idx >= 0, "Missing material for farfield boundary attribute!");
    fb.AddMaterialProperty(ceed_attr, inv_z(0, 0, mat_idx), coeff);
  }
}

}  // namespace palace
