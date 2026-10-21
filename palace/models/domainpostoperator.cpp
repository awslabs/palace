// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "domainpostoperator.hpp"

#include <algorithm>
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
      // The coefficient is exactly zero outside of the postprocessing domain and the
      // operator is only ever applied, never assembled as a matrix, so a postprocessing
      // domain inside the shared lossy-domain set is assembled over that subset only.
      m.SkipZeroCoefficientElements();
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
        m.SkipZeroCoefficientElements();
        m.AddDomainIntegrator<MassIntegrator>(muinv_func);
        M_mag_i = m.PartialAssemble();
      }
      else
      {
        MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                               mat_op.GetInvPermeability());
        muinv_func.RestrictCoefficient(mat_op.GetCeedAttributes(data.attributes));
        BilinearForm m(rt_fespace);
        m.SkipZeroCoefficientElements();
        m.AddDomainIntegrator<VectorFEMassIntegrator>(muinv_func);
        M_mag_i = m.PartialAssemble();
      }
    }
    M_i.emplace(idx, std::make_pair(std::move(M_elec_i), std::move(M_mag_i)));
  }
  SetUpFrequencyDependentForms(postpro, mat_op, nd_fespace);
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
        m.SkipZeroCoefficientElements();
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
      // Volume magnetic energy only; the London kinetic sheet energy is added explicitly in
      // the driver, so GetMagneticFieldEnergy stays the physical field energy.
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
        m.SkipZeroCoefficientElements();
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
      m.SkipZeroCoefficientElements();
      m.AddDomainIntegrator<VectorFEMassIntegrator>(epsilon_func);
      M_elec_i = m.PartialAssemble();
    }
    M_i.emplace(idx, std::make_pair(std::move(M_elec_i), nullptr));
  }
  SetUpFrequencyDependentForms(iodata.domains.postpro, mat_op, nd_fespace);
}

void DomainPostOperator::SetUpFrequencyDependentForms(const config::DomainPostData &postpro,
                                                      const MaterialOperator &mat_op,
                                                      const FiniteElementSpace &nd_fespace)
{
  // Every rank builds the same forms for the globally supported materials, empty where a
  // rank has none of their elements.
  this->mat_op = &mat_op;
  auto MakeForms = [&](const std::vector<int> *domain_attributes)
  {
    FrequencyDependentForms forms;
    for (std::size_t m = 0; m < mat_op.NumFrequencyDependentPermittivityMaterials(); m++)
    {
      if (!mat_op.HasFrequencyDependentPermittivitySupport(m))
      {
        continue;
      }
      std::vector<int> attributes = mat_op.GetFrequencyDependentPermittivityAttributes(m);
      if (domain_attributes)
      {
        std::erase_if(attributes,
                      [domain_attributes](int attr)
                      {
                        return std::ranges::find(*domain_attributes, attr) ==
                               domain_attributes->end();
                      });
      }
      MaterialPropertyCoefficient f(mat_op.MaxCeedAttribute());
      if (!attributes.empty())
      {
        f.AddMaterialProperty(mat_op.GetCeedAttributes(attributes), 1.0);
      }
      BilinearForm a(nd_fespace);
      if (!f.empty())
      {
        a.AddDomainIntegrator<VectorFEMassIntegrator>(f);
      }
      forms.emplace_back(m, a.PartialAssemble());
    }
    return forms;
  };
  M_elec_fd = MakeForms(nullptr);
  for (const auto &[idx, data] : postpro.energy)
  {
    M_elec_fd_i.emplace(idx, MakeForms(&data.attributes));
  }
}

double DomainPostOperator::GetFrequencyDependentEnergyWeight(std::size_t material_idx,
                                                             double omega,
                                                             ElectricEnergy type) const
{
  MFEM_VERIFY(omega > 0.0, "The electric energy of a frequency-dependent material requires "
                           "a positive frequency!");
  return (type == ElectricEnergy::STORED)
             ? mat_op->GetFrequencyDependentStoredEnergySusceptibility(material_idx, omega)
             : mat_op->GetFrequencyDependentSusceptibilityReal(material_idx, omega);
}

double DomainPostOperator::ElectricFieldEnergy(const Operator &M,
                                               const FrequencyDependentForms &M_fd,
                                               const GridFunction &E, double omega,
                                               ElectricEnergy type) const
{
  auto dot = [&](const Operator &op)
  {
    op.Mult(E.Real(), D);
    double d = linalg::LocalDot(E.Real(), D);
    if (E.HasImag())
    {
      op.Mult(E.Imag(), D);
      d += linalg::LocalDot(E.Imag(), D);
    }
    return d;
  };
  double energy = dot(M);
  for (const auto &[m, op] : M_fd)
  {
    energy += GetFrequencyDependentEnergyWeight(m, omega, type) * dot(*op);
  }
  Mpi::GlobalSum(1, &energy, E.GetComm());
  return 0.5 * energy;
}

double DomainPostOperator::GetElectricFieldEnergy(const GridFunction &E, double omega,
                                                  ElectricEnergy type) const
{
  MFEM_VERIFY(M_elec, "Domain postprocessing is not configured for electric field energy "
                      "calculation!");
  return ElectricFieldEnergy(*M_elec, M_elec_fd, E, omega, type);
}

double DomainPostOperator::GetMagneticFieldEnergy(const GridFunction &B) const
{
  if (M_mag)
  {
    M_mag->Mult(B.Real(), H);
    double dot = linalg::LocalDot(B.Real(), H);
    if (B.HasImag())
    {
      M_mag->Mult(B.Imag(), H);
      dot += linalg::LocalDot(B.Imag(), H);
    }
    Mpi::GlobalSum(1, &dot, B.GetComm());
    return 0.5 * dot;
  }
  return 0.0;
}

double DomainPostOperator::GetDomainElectricFieldEnergy(int idx, const GridFunction &E,
                                                        double omega,
                                                        ElectricEnergy type) const
{
  // Compute the electric field energy integral for only a portion of the domain.
  auto it = M_i.find(idx);
  MFEM_VERIFY(it != M_i.end(),
              "Invalid domain index when postprocessing domain electric field energy!");
  if (!it->second.first)
  {
    return 0.0;
  }
  static const FrequencyDependentForms none;
  auto fd = M_elec_fd_i.find(idx);
  return ElectricFieldEnergy(*it->second.first,
                             (fd != M_elec_fd_i.end()) ? fd->second : none, E, omega, type);
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
