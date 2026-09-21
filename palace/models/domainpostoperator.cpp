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
#include "utils/constants.hpp"
#include "utils/iodata.hpp"

namespace palace
{

DomainPostOperator::DomainPostOperator(const config::DomainPostData &postpro,
                                       const MaterialOperator &mat_op,
                                       const FiniteElementSpace &nd_fespace,
                                       const FiniteElementSpace &rt_fespace)
{
  // Mass operators are always partially assembled.
  const int dim = nd_fespace.Dimension();
  const auto curl_map = rt_fespace.GetFEColl().GetMapType(dim);
  MFEM_VERIFY(nd_fespace.GetFEColl().GetMapType(dim) == mfem::FiniteElement::H_CURL &&
                  (curl_map == mfem::FiniteElement::H_DIV ||
                   curl_map == mfem::FiniteElement::INTEGRAL),
              "Unexpected finite element space types for domain energy postprocessing!");
  {
    // Construct ND mass matrix to compute the electric field energy integral as:
    //              E_elec = 1/2 Re{∫_Ω Dᴴ E dV} as (M_eps * e)ᴴ e.
    // Only the real part of the permeability contributes to the energy (imaginary part
    // cancels out in the inner product due to symmetry).
    MaterialPropertyCoefficient epsilon_func(mat_op.GetAttributeToMaterial(),
                                             mat_op.GetPermittivityReal());
    BilinearForm m(nd_fespace);
    m.AddDomainIntegrator<VectorFEMassIntegrator>(epsilon_func);
    M_elec = m.PartialAssemble();
    D.SetSize(M_elec->Height());
    D.UseDevice(true);
  }
  {
    // Construct mass matrix for B-field to compute the magnetic field energy integral as:
    //              E_mag = 1/2 Re{∫_Ω Hᴴ B dV} as (M_muinv * b)ᴴ b.
    // In 2D, B is scalar (L2 space) so use scalar mass integrator with scalar μ⁻¹.
    // In 3D, B is a vector (RT space) so use vector FE mass integrator.
    if (curl_map == mfem::FiniteElement::INTEGRAL)
    {
      // Scalar curl (2D): need scalar μ⁻¹ (z-z component).
      // GetCurlCurlInvPermeability() returns 1x1 for 2D MaterialOperator.
      MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                             mat_op.GetCurlCurlInvPermeability());
      BilinearForm m(rt_fespace);
      m.AddDomainIntegrator<MassIntegrator>(muinv_func);
      M_mag = m.PartialAssemble();
    }
    else
    {
      MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                             mat_op.GetInvPermeability());
      BilinearForm m(rt_fespace);
      m.AddDomainIntegrator<VectorFEMassIntegrator>(muinv_func);
      M_mag = m.PartialAssemble();
    }
    H.SetSize(M_mag->Height());
    H.UseDevice(true);
  }

  // Use the provided domain postprocessing indices for postprocessing the electric and
  // magnetic field energy in specific regions of the domain.
  for (const auto &[idx, data] : postpro.energy)
  {
    std::unique_ptr<Operator> M_elec_i, M_mag_i;
    {
      MaterialPropertyCoefficient epsilon_func(mat_op.GetAttributeToMaterial(),
                                               mat_op.GetPermittivityReal());
      epsilon_func.RestrictCoefficient(mat_op.GetCeedAttributes(data.attributes));
      BilinearForm m(nd_fespace);
      m.AddDomainIntegrator<VectorFEMassIntegrator>(epsilon_func);
      M_elec_i = m.PartialAssemble();
    }
    {
      if (curl_map == mfem::FiniteElement::INTEGRAL)
      {
        MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                               mat_op.GetCurlCurlInvPermeability());
        muinv_func.RestrictCoefficient(mat_op.GetCeedAttributes(data.attributes));
        BilinearForm m(rt_fespace);
        m.AddDomainIntegrator<MassIntegrator>(muinv_func);
        M_mag_i = m.PartialAssemble();
      }
      else
      {
        MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                               mat_op.GetInvPermeability());
        muinv_func.RestrictCoefficient(mat_op.GetCeedAttributes(data.attributes));
        BilinearForm m(rt_fespace);
        m.AddDomainIntegrator<VectorFEMassIntegrator>(muinv_func);
        M_mag_i = m.PartialAssemble();
      }
    }
    M_i.emplace(idx, std::make_pair(std::move(M_elec_i), std::move(M_mag_i)));
  }
}

DomainPostOperator::DomainPostOperator(const IoData &iodata, const MaterialOperator &mat_op,
                                       const FiniteElementSpace &nd_fespace,
                                       const FiniteElementSpace &rt_fespace)
  : DomainPostOperator(iodata.domains.postpro, mat_op, nd_fespace, rt_fespace)
{
}

DomainPostOperator::DomainPostOperator(const config::DomainPostData &postpro,
                                       const MaterialOperator &mat_op,
                                       const FiniteElementSpace &fespace)
{
  const auto map_type = fespace.GetFEColl().GetMapType(fespace.Dimension());
  if (map_type == mfem::FiniteElement::VALUE)
  {
    // H1 space for voltage and electric field energy.
    {
      MaterialPropertyCoefficient epsilon_func(mat_op.GetAttributeToMaterial(),
                                               mat_op.GetPermittivityReal());
      BilinearForm m(fespace);
      m.AddDomainIntegrator<DiffusionIntegrator>(epsilon_func);
      M_elec = m.PartialAssemble();
      D.SetSize(M_elec->Height());
      D.UseDevice(true);
    }

    for (const auto &[idx, data] : postpro.energy)
    {
      std::unique_ptr<Operator> M_elec_i;
      {
        MaterialPropertyCoefficient epsilon_func(mat_op.GetAttributeToMaterial(),
                                                 mat_op.GetPermittivityReal());
        epsilon_func.RestrictCoefficient(mat_op.GetCeedAttributes(data.attributes));
        BilinearForm m(fespace);
        m.AddDomainIntegrator<DiffusionIntegrator>(epsilon_func);
        M_elec_i = m.PartialAssemble();
      }
      M_i.emplace(idx, std::make_pair(std::move(M_elec_i), nullptr));
    }
  }
  else if (map_type == mfem::FiniteElement::H_CURL)
  {
    // H(curl) space for magnetic vector potential and magnetic field energy.
    // (This is the magnetostatic case — creates curl-curl mass for B-field energy.)
    {
      MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                             mat_op.GetCurlCurlInvPermeability());
      BilinearForm m(fespace);
      m.AddDomainIntegrator<CurlCurlIntegrator>(muinv_func);
      M_mag = m.PartialAssemble();
      H.SetSize(M_mag->Height());
      H.UseDevice(true);
    }

    for (const auto &[idx, data] : postpro.energy)
    {
      std::unique_ptr<Operator> M_mag_i;
      {
        MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                               mat_op.GetCurlCurlInvPermeability());
        muinv_func.RestrictCoefficient(mat_op.GetCeedAttributes(data.attributes));
        BilinearForm m(fespace);
        m.AddDomainIntegrator<CurlCurlIntegrator>(muinv_func);
        M_mag_i = m.PartialAssemble();
      }
      M_i.emplace(idx, std::make_pair(nullptr, std::move(M_mag_i)));
    }
  }
  else
  {
    MFEM_ABORT("Unexpected finite element space type for domain energy postprocessing!");
  }
}

DomainPostOperator::DomainPostOperator(const IoData &iodata, const MaterialOperator &mat_op,
                                       const FiniteElementSpace &fespace)
  : DomainPostOperator(iodata.domains.postpro, mat_op, fespace)
{
}

DomainPostOperator::DomainPostOperator(const IoData &iodata, const MaterialOperator &mat_op,
                                       const FiniteElementSpace &nd_fespace,
                                       bool electric_energy_only)
{
  // Mode analysis: ND space for electric field energy only (VectorFE mass with ε).
  MFEM_VERIFY(nd_fespace.GetFEColl().GetMapType(nd_fespace.Dimension()) ==
                  mfem::FiniteElement::H_CURL,
              "Electric energy only constructor requires H(curl) space!");
  {
    MaterialPropertyCoefficient epsilon_func(mat_op.GetAttributeToMaterial(),
                                             mat_op.GetPermittivityReal());
    BilinearForm m(nd_fespace);
    m.AddDomainIntegrator<VectorFEMassIntegrator>(epsilon_func);
    M_elec = m.PartialAssemble();
    D.SetSize(M_elec->Height());
    D.UseDevice(true);
  }

  for (const auto &[idx, data] : iodata.domains.postpro.energy)
  {
    std::unique_ptr<Operator> M_elec_i;
    {
      MaterialPropertyCoefficient epsilon_func(mat_op.GetAttributeToMaterial(),
                                               mat_op.GetPermittivityReal());
      epsilon_func.RestrictCoefficient(mat_op.GetCeedAttributes(data.attributes));
      BilinearForm m(nd_fespace);
      m.AddDomainIntegrator<VectorFEMassIntegrator>(epsilon_func);
      M_elec_i = m.PartialAssemble();
    }
    M_i.emplace(idx, std::make_pair(std::move(M_elec_i), nullptr));
  }
}

double DomainPostOperator::GetFieldEnergy(const Operator &M, const GridFunction &u,
                                          Vector &Mu)
{
  // Instantaneous energy 1/2 uᵀ M u for a real field, or time-averaged energy
  // 1/4 (u_rᵀ M u_r + u_iᵀ M u_i) for a complex peak phasor.
  M.Mult(u.Real(), Mu);
  double dot = linalg::LocalDot(u.Real(), Mu);
  if (u.HasImag())
  {
    M.Mult(u.Imag(), Mu);
    dot += linalg::LocalDot(u.Imag(), Mu);
  }
  Mpi::GlobalSum(1, &dot, u.GetComm());
  return 0.5 * electromagnetics::TimeAverageWeight(u.HasImag()) * dot;
}

double DomainPostOperator::GetElectricFieldEnergy(const GridFunction &E) const
{
  if (M_elec)
  {
    return GetFieldEnergy(*M_elec, E, D);
  }
  MFEM_ABORT(
      "Domain postprocessing is not configured for electric field energy calculation!");
  return 0.0;
}

double DomainPostOperator::GetMagneticFieldEnergy(const GridFunction &B) const
{
  if (M_mag)
  {
    return GetFieldEnergy(*M_mag, B, H);
  }
  return 0.0;
}

double DomainPostOperator::GetDomainElectricFieldEnergy(int idx,
                                                        const GridFunction &E) const
{
  // Compute the electric field energy integral for only a portion of the domain.
  auto it = M_i.find(idx);
  MFEM_VERIFY(it != M_i.end(),
              "Invalid domain index when postprocessing domain electric field energy!");
  if (!it->second.first)
  {
    return 0.0;
  }
  return GetFieldEnergy(*it->second.first, E, D);
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
  return GetFieldEnergy(*it->second.second, B, H);
}

}  // namespace palace
