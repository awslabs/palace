// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "domainpostoperator.hpp"

#include "fem/bilinearform.hpp"
#include "fem/coefficient.hpp"
#include "fem/integrator.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

DomainPostOperator::DomainPostOperator(const IoData &iodata, const MaterialOperator &mat_op,
                                       const mfem::ParFiniteElementSpace *nd_fespace,
                                       const mfem::ParFiniteElementSpace *rt_fespace)
{
  // Mass operators are always partially assembled.
  constexpr auto MatTypeEps = MaterialPropertyType::PERMITTIVITY_REAL;
  constexpr auto MatTypeMuInv = MaterialPropertyType::INV_PERMEABILITY;
  if (nd_fespace)
  {
    // Construct ND mass matrix to compute the electric field energy integral as:
    //              E_elec = 1/2 Re{∫_Ω Dᴴ E dV} as (M_eps * e)ᴴ e.
    // Only the real part of the permeability contributes to the energy (imaginary part
    // cancels out in the inner product due to symmetry).
    MaterialPropertyCoefficient<MatTypeEps> epsilon_func(mat_op);
    BilinearForm m_nd(*nd_fespace);
    m_nd.AddDomainIntegrator<VectorFEMassIntegrator>(epsilon_func);
    M_ND = m_nd.Assemble();
    D.SetSize(M_ND->Height());
    D.UseDevice(true);
  }

  if (rt_fespace)
  {
    // Construct RT mass matrix to compute the magnetic field energy integral as:
    //              E_mag = 1/2 Re{∫_Ω Bᴴ H dV} as (M_muinv * b)ᴴ b.
    MaterialPropertyCoefficient<MatTypeMuInv> muinv_func(mat_op);
    BilinearForm m_rt(*rt_fespace);
    m_rt.AddDomainIntegrator<VectorFEMassIntegrator>(muinv_func);
    M_RT = m_rt.Assemble();
    H.SetSize(M_RT->Height());
    H.UseDevice(true);
  }

  // Use the provided domain postprocessing indices for postprocessing the electric and
  // magnetic field energy in specific regions of the domain.
  const auto &mesh = nd_fespace ? *nd_fespace->GetParMesh() : *rt_fespace->GetParMesh();
  int attr_max = mesh.attributes.Max();
  for (const auto &[idx, data] : iodata.domains.postpro.energy)
  {
    mfem::Array<int> attr_marker(attr_max);
    attr_marker = 0;
    for (auto attr : data.attributes)
    {
      attr_marker[attr - 1] = 1;
    }
    std::unique_ptr<Operator> M_ND_i, M_RT_i;
    if (nd_fespace)
    {
      SumMatrixCoefficient epsilon_func_i(nd_fespace->GetParMesh()->SpaceDimension());
      epsilon_func_i.AddCoefficient(
          std::make_unique<MaterialPropertyCoefficient<MatTypeEps>>(mat_op), attr_marker);
      BilinearForm m_nd_i(*nd_fespace);
      m_nd_i.AddDomainIntegrator<VectorFEMassIntegrator>(epsilon_func_i);
      M_ND_i = m_nd_i.Assemble();
    }
    if (rt_fespace)
    {
      SumMatrixCoefficient muinv_func_i(rt_fespace->GetParMesh()->SpaceDimension());
      muinv_func_i.AddCoefficient(
          std::make_unique<MaterialPropertyCoefficient<MatTypeMuInv>>(mat_op), attr_marker);
      BilinearForm m_rt_i(*rt_fespace);
      m_rt_i.AddDomainIntegrator<VectorFEMassIntegrator>(muinv_func_i);
      M_RT_i = m_rt_i.Assemble();
    }
    M_i.emplace(idx, std::make_pair(std::move(M_ND_i), std::move(M_RT_i)));
  }
}

double
DomainPostOperator::GetElectricFieldEnergy(const mfem::ParComplexGridFunction &E) const
{
  if (M_ND)
  {
    M_ND->Mult(E.real(), D);
    double res = mfem::InnerProduct(E.real(), D);
    M_ND->Mult(E.imag(), D);
    res += mfem::InnerProduct(E.imag(), D);
    Mpi::GlobalSum(1, &res, E.ParFESpace()->GetComm());
    return 0.5 * res;
  }
  MFEM_ABORT(
      "Domain postprocessing is not configured for electric field energy calculation!");
  return 0.0;
}

double DomainPostOperator::GetElectricFieldEnergy(const mfem::ParGridFunction &E) const
{
  if (M_ND)
  {
    M_ND->Mult(E, D);
    double res = mfem::InnerProduct(E, D);
    Mpi::GlobalSum(1, &res, E.ParFESpace()->GetComm());
    return 0.5 * res;
  }
  MFEM_ABORT(
      "Domain postprocessing is not configured for electric field energy calculation!");
  return 0.0;
}

double
DomainPostOperator::GetMagneticFieldEnergy(const mfem::ParComplexGridFunction &B) const
{
  if (M_RT)
  {
    M_RT->Mult(B.real(), H);
    double res = mfem::InnerProduct(B.real(), H);
    M_RT->Mult(B.imag(), H);
    res += mfem::InnerProduct(B.imag(), H);
    Mpi::GlobalSum(1, &res, B.ParFESpace()->GetComm());
    return 0.5 * res;
  }
  MFEM_ABORT(
      "Domain postprocessing is not configured for magnetic field energy calculation!");
  return 0.0;
}

double DomainPostOperator::GetMagneticFieldEnergy(const mfem::ParGridFunction &B) const
{
  if (M_RT)
  {
    M_RT->Mult(B, H);
    double res = mfem::InnerProduct(B, H);
    Mpi::GlobalSum(1, &res, B.ParFESpace()->GetComm());
    return 0.5 * res;
  }
  MFEM_ABORT(
      "Domain postprocessing is not configured for magnetic field energy calculation!");
  return 0.0;
}

double DomainPostOperator::GetDomainElectricFieldEnergy(
    int idx, const mfem::ParComplexGridFunction &E) const
{
  // Compute the electric field energy integral for only a portion of the domain.
  auto it = M_i.find(idx);
  MFEM_VERIFY(it != M_i.end() && it->second.first,
              "Invalid domain index when postprocessing domain electric field energy!");
  if (!it->second.first)
  {
    return 0.0;
  }
  it->second.first->Mult(E.real(), D);
  double res = mfem::InnerProduct(E.real(), D);
  it->second.first->Mult(E.imag(), D);
  res += mfem::InnerProduct(E.imag(), D);
  Mpi::GlobalSum(1, &res, E.ParFESpace()->GetComm());
  return 0.5 * res;
}

double
DomainPostOperator::GetDomainElectricFieldEnergy(int idx,
                                                 const mfem::ParGridFunction &E) const
{
  auto it = M_i.find(idx);
  MFEM_VERIFY(it != M_i.end() && it->second.first,
              "Invalid domain index when postprocessing domain electric field energy!");
  if (!it->second.first)
  {
    return 0.0;
  }
  it->second.first->Mult(E, D);
  double res = mfem::InnerProduct(E, D);
  Mpi::GlobalSum(1, &res, E.ParFESpace()->GetComm());
  return 0.5 * res;
}

double DomainPostOperator::GetDomainMagneticFieldEnergy(
    int idx, const mfem::ParComplexGridFunction &B) const
{
  // Compute the magnetic field energy integral for only a portion of the domain.
  auto it = M_i.find(idx);
  MFEM_VERIFY(it != M_i.end(),
              "Invalid domain index when postprocessing domain magnetic field energy!");
  if (!it->second.second)
  {
    return 0.0;
  }
  it->second.second->Mult(B.real(), H);
  double res = mfem::InnerProduct(B.real(), H);
  it->second.second->Mult(B.imag(), H);
  res += mfem::InnerProduct(B.imag(), H);
  Mpi::GlobalSum(1, &res, B.ParFESpace()->GetComm());
  return 0.5 * res;
}

double
DomainPostOperator::GetDomainMagneticFieldEnergy(int idx,
                                                 const mfem::ParGridFunction &B) const
{
  auto it = M_i.find(idx);
  MFEM_VERIFY(it != M_i.end() && it->second.second,
              "Invalid domain index when postprocessing domain magnetic field energy!");
  if (!it->second.second)
  {
    return 0.0;
  }
  it->second.second->Mult(B, H);
  double res = mfem::InnerProduct(B, H);
  Mpi::GlobalSum(1, &res, B.ParFESpace()->GetComm());
  return 0.5 * res;
}

}  // namespace palace
