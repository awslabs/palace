// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "domainpostoperator.hpp"

#include <mfem.hpp>
#include "fem/bilinearform.hpp"
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

DomainPostOperator::DomainPostOperator(const IoData &iodata, const MaterialOperator &mat_op,
                                       const FiniteElementSpace *nd_fespace,
                                       const FiniteElementSpace *rt_fespace)
{
  // Mass operators are always partially assembled.
  if (nd_fespace)
  {
    // Construct ND mass matrix to compute the electric field energy integral as:
    //              E_elec = 1/2 Re{∫_Ω Dᴴ E dV} as (M_eps * e)ᴴ e.
    // Only the real part of the permeability contributes to the energy (imaginary part
    // cancels out in the inner product due to symmetry).
    MaterialPropertyCoefficient epsilon_func(mat_op.GetAttributeToMaterial(),
                                             mat_op.GetPermittivityReal());
    BilinearForm m_nd(*nd_fespace);
    m_nd.AddDomainIntegrator<VectorFEMassIntegrator>(epsilon_func);
    M_ND = m_nd.PartialAssemble();
    D.SetSize(M_ND->Height());
    D.UseDevice(true);
  }

  if (rt_fespace)
  {
    // Construct RT mass matrix to compute the magnetic field energy integral as:
    //              E_mag = 1/2 Re{∫_Ω Bᴴ H dV} as (M_muinv * b)ᴴ b.
    MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                           mat_op.GetInvPermeability());
    BilinearForm m_rt(*rt_fespace);
    m_rt.AddDomainIntegrator<VectorFEMassIntegrator>(muinv_func);
    M_RT = m_rt.PartialAssemble();
    H.SetSize(M_RT->Height());
    H.UseDevice(true);
  }

  // Use the provided domain postprocessing indices for postprocessing the electric and
  // magnetic field energy in specific regions of the domain.
  for (const auto &[idx, data] : iodata.domains.postpro.energy)
  {
    std::unique_ptr<Operator> M_ND_i, M_RT_i;
    if (nd_fespace)
    {
      MaterialPropertyCoefficient epsilon_func(mat_op.GetAttributeToMaterial(),
                                               mat_op.GetPermittivityReal());
      epsilon_func.RestrictCoefficient(mat_op.GetCeedAttributes(data.attributes));
      BilinearForm m_nd_i(*nd_fespace);
      m_nd_i.AddDomainIntegrator<VectorFEMassIntegrator>(epsilon_func);
      M_ND_i = m_nd_i.PartialAssemble();
    }
    if (rt_fespace)
    {
      MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                             mat_op.GetInvPermeability());
      muinv_func.RestrictCoefficient(mat_op.GetCeedAttributes(data.attributes));
      BilinearForm m_rt_i(*rt_fespace);
      m_rt_i.AddDomainIntegrator<VectorFEMassIntegrator>(muinv_func);
      M_RT_i = m_rt_i.PartialAssemble();
    }
    M_i.emplace(idx, std::make_pair(std::move(M_ND_i), std::move(M_RT_i)));
  }
}

double DomainPostOperator::GetElectricFieldEnergy(const GridFunction &E) const
{
  if (M_ND)
  {
    M_ND->Mult(E.Real(), D);
    double dot = linalg::LocalDot(E.Real(), D);
    if (E.HasImag())
    {
      M_ND->Mult(E.Imag(), D);
      dot += linalg::LocalDot(E.Imag(), D);
    }
    Mpi::GlobalSum(1, &dot, E.GetComm());
    return 0.5 * dot;
  }
  MFEM_ABORT(
      "Domain postprocessing is not configured for electric field energy calculation!");
  return 0.0;
}

double DomainPostOperator::GetMagneticFieldEnergy(const GridFunction &B) const
{
  if (M_RT)
  {
    M_RT->Mult(B.Real(), H);
    double dot = linalg::LocalDot(B.Real(), H);
    if (B.HasImag())
    {
      M_RT->Mult(B.Imag(), H);
      dot += linalg::LocalDot(B.Imag(), H);
    }
    Mpi::GlobalSum(1, &dot, B.GetComm());
    return 0.5 * dot;
  }
  MFEM_ABORT(
      "Domain postprocessing is not configured for magnetic field energy calculation!");
  return 0.0;
}

double DomainPostOperator::GetDomainElectricFieldEnergy(int idx,
                                                        const GridFunction &E) const
{
  // Compute the electric field energy integral for only a portion of the domain.
  auto it = M_i.find(idx);
  MFEM_VERIFY(it != M_i.end() && it->second.first,
              "Invalid domain index when postprocessing domain electric field energy!");
  if (!it->second.first)
  {
    return 0.0;
  }
  it->second.first->Mult(E.Real(), D);
  double dot = linalg::LocalDot(E.Real(), D);
  if (E.HasImag())
  {
    it->second.first->Mult(E.Imag(), D);
    dot += linalg::LocalDot(E.Imag(), D);
  }
  Mpi::GlobalSum(1, &dot, E.GetComm());
  return 0.5 * dot;
}

double DomainPostOperator::GetDomainMagneticFieldEnergy(int idx,
                                                        const GridFunction &B) const
{
  // Compute the magnetic field energy integral for only a portion of the domain.
  auto it = M_i.find(idx);
  MFEM_VERIFY(it != M_i.end(),
              "Invalid domain index when postprocessing domain magnetic field energy!");
  if (!it->second.second)
  {
    return 0.0;
  }
  it->second.second->Mult(B.Real(), H);
  double dot = linalg::LocalDot(B.Real(), H);
  if (B.HasImag())
  {
    it->second.second->Mult(B.Imag(), H);
    dot += linalg::LocalDot(B.Imag(), H);
  }
  Mpi::GlobalSum(1, &dot, B.GetComm());
  return 0.5 * dot;
}

}  // namespace palace
