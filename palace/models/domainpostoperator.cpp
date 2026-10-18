// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "domainpostoperator.hpp"

#include <mfem.hpp>
#include "fem/bilinearform.hpp"
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "models/dispersiveconductivityoperator.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

namespace
{

// Kinetic energy and loss as weighted sums of the forms 1/2 eᴴ M e, with M = get_M(form).
template <typename GetOperator>
std::pair<double, double>
KineticEnergyAndLoss(const std::vector<DomainPostOperator::ConductionForm> &forms,
                     GetOperator &&get_M, const GridFunction &E, double omega, Vector &D)
{
  if (forms.empty())
  {
    return {0.0, 0.0};
  }
  std::vector<double> dots(forms.size(), 0.0);
  for (std::size_t k = 0; k < forms.size(); k++)
  {
    const Operator &M = get_M(forms[k]);
    M.Mult(E.Real(), D);
    dots[k] = linalg::LocalDot(E.Real(), D);
    if (E.HasImag())
    {
      M.Mult(E.Imag(), D);
      dots[k] += linalg::LocalDot(E.Imag(), D);
    }
  }
  Mpi::GlobalSum(static_cast<int>(dots.size()), dots.data(), E.GetComm());
  double kinetic = 0.0, loss = 0.0;
  for (std::size_t k = 0; k < forms.size(); k++)
  {
    const auto [w_kinetic, w_loss] = forms[k].weights(omega);
    kinetic += 0.5 * w_kinetic * dots[k];
    loss += 0.5 * w_loss * dots[k];
  }
  return {kinetic, loss};
}

}  // namespace

DomainPostOperator::DomainPostOperator(const config::DomainPostData &postpro,
                                       const MaterialOperator &mat_op,
                                       const FiniteElementSpace &nd_fespace,
                                       const FiniteElementSpace &rt_fespace,
                                       const DispersiveConductivityOperator *disp_op)
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

  // Kinetic energy and conduction loss forms.
  if (!disp_op)
  {
    return;
  }
  auto add_form = [&](const MaterialPropertyCoefficient &coeff,
                      std::function<std::pair<double, double>(double)> weights)
  {
    auto assemble = [&](const MaterialPropertyCoefficient &q)
    {
      // Processes without elements of the attributes have an empty coefficient.
      BilinearForm m(nd_fespace);
      if (!q.empty())
      {
        m.AddDomainIntegrator<VectorFEMassIntegrator>(q);
      }
      return m.PartialAssemble();
    };
    auto &form = cond_forms.emplace_back();
    form.M = assemble(coeff);
    for (const auto &[idx, data] : postpro.energy)
    {
      MaterialPropertyCoefficient coeff_i(coeff);
      coeff_i.RestrictCoefficient(mat_op.GetCeedAttributes(data.attributes));
      form.M_i.emplace(idx, assemble(coeff_i));
    }
    form.weights = std::move(weights);
  };
  if (mat_op.HasLondonDepth())
  {
    has_kinetic = true;
    add_form({mat_op.GetAttributeToMaterial(), mat_op.GetInvLondonDepth()},
             [](double omega) { return std::make_pair(1.0 / (omega * omega), 0.0); });
  }
  if (mat_op.HasConductivity())
  {
    has_loss = true;
    add_form({mat_op.GetAttributeToMaterial(), mat_op.GetConductivity()},
             [](double) { return std::make_pair(0.0, 1.0); });
  }
  for (std::size_t g = 0; g < disp_op->Size(); g++)
  {
    if (disp_op->IsActive(g))
    {
      has_kinetic = has_loss = true;
      MaterialPropertyCoefficient unit(mat_op.MaxCeedAttribute());
      disp_op->AddUnitCoefficients(g, unit);
      add_form(unit,
               [disp_op, g](double omega)
               {
                 return std::make_pair(disp_op->EvaluateKineticScalar(g, omega),
                                       disp_op->EvaluateLossScalar(g, omega));
               });
    }
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

double DomainPostOperator::GetElectricFieldEnergy(const GridFunction &E) const
{
  if (M_elec)
  {
    M_elec->Mult(E.Real(), D);
    double dot = linalg::LocalDot(E.Real(), D);
    if (E.HasImag())
    {
      M_elec->Mult(E.Imag(), D);
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

std::pair<double, double> DomainPostOperator::GetKineticEnergyAndLoss(const GridFunction &E,
                                                                      double omega) const
{
  return KineticEnergyAndLoss(
      cond_forms, [](const ConductionForm &form) -> const Operator & { return *form.M; }, E,
      omega, D);
}

std::pair<double, double>
DomainPostOperator::GetDomainKineticEnergyAndLoss(int idx, const GridFunction &E,
                                                  double omega) const
{
  MFEM_VERIFY(M_i.contains(idx),
              "Invalid domain index when postprocessing domain kinetic energy and loss!");
  return KineticEnergyAndLoss(
      cond_forms, [idx](const ConductionForm &form) -> const Operator &
      { return *form.M_i.at(idx); }, E, omega, D);
}

}  // namespace palace
