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

// Local part of uᴴ M u for the real and imaginary parts of u.
double LocalQuadratic(const Operator &M, const GridFunction &u, Vector &w)
{
  M.Mult(u.Real(), w);
  double dot = linalg::LocalDot(u.Real(), w);
  if (u.HasImag())
  {
    M.Mult(u.Imag(), w);
    dot += linalg::LocalDot(u.Imag(), w);
  }
  return dot;
}

// Kinetic energy and loss as weighted sums of the forms 1/2 eᴴ M e, with (M, M_n) =
// get_M(form) the forms of E and, for mode analysis, En.
template <typename GetOperators>
std::pair<double, double>
KineticEnergyAndLoss(const std::vector<DomainPostOperator::ConductionForm> &forms,
                     GetOperators &&get_M, const GridFunction &E, const GridFunction *En,
                     double omega, Vector &D, Vector &D_n)
{
  if (forms.empty())
  {
    return {0.0, 0.0};
  }
  std::vector<double> dots(forms.size(), 0.0);
  for (std::size_t k = 0; k < forms.size(); k++)
  {
    const auto [M, M_n] = get_M(forms[k]);
    dots[k] = LocalQuadratic(*M, E, D) + (En ? LocalQuadratic(*M_n, *En, D_n) : 0.0);
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

// Partially assembled mass form of a coefficient: vector on ND and scalar on H1 spaces.
// Processes without elements of the attributes have an empty coefficient.
std::unique_ptr<Operator> AssembleMass(const FiniteElementSpace &fespace,
                                       const MaterialPropertyCoefficient &q)
{
  BilinearForm m(fespace);
  if (!q.empty())
  {
    if (fespace.GetFEColl().GetMapType(fespace.Dimension()) == mfem::FiniteElement::H_CURL)
    {
      m.AddDomainIntegrator<VectorFEMassIntegrator>(q);
    }
    else
    {
      m.AddDomainIntegrator<MassIntegrator>(q);
    }
  }
  return m.PartialAssemble();
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

  if (disp_op)
  {
    AddConductionForms(postpro, mat_op, nd_fespace, nullptr, *disp_op);
  }
}

DomainPostOperator::DomainPostOperator(const config::DomainPostData &postpro,
                                       const MaterialOperator &mat_op,
                                       const FiniteElementSpace &nd_fespace,
                                       const FiniteElementSpace &h1_fespace,
                                       const FiniteElementSpace &l2_fespace,
                                       const DispersiveConductivityOperator &disp_op)
  : DomainPostOperator(postpro, mat_op, nd_fespace, l2_fespace)
{
  // Out-of-plane permittivity on the H1 space, and in-plane inverse permeability on the ND
  // space.
  const MaterialPropertyCoefficient eps_n(mat_op.GetAttributeToMaterial(),
                                          mat_op.GetPermittivityScalar()),
      muinv_t(mat_op.GetAttributeToMaterial(), mat_op.GetInvPermeability());
  M_elec_n = AssembleMass(h1_fespace, eps_n);
  M_mag_t = AssembleMass(nd_fespace, muinv_t);
  D_n.SetSize(M_elec_n->Height());
  D_n.UseDevice(true);
  for (const auto &[idx, data] : postpro.energy)
  {
    const auto attr_list = mat_op.GetCeedAttributes(data.attributes);
    MaterialPropertyCoefficient eps_n_i(eps_n), muinv_t_i(muinv_t);
    eps_n_i.RestrictCoefficient(attr_list);
    muinv_t_i.RestrictCoefficient(attr_list);
    M_n_i.emplace(idx, std::make_pair(AssembleMass(h1_fespace, eps_n_i),
                                      AssembleMass(nd_fespace, muinv_t_i)));
  }
  AddConductionForms(postpro, mat_op, nd_fespace, &h1_fespace, disp_op);
}

void DomainPostOperator::AddConductionForms(const config::DomainPostData &postpro,
                                            const MaterialOperator &mat_op,
                                            const FiniteElementSpace &nd_fespace,
                                            const FiniteElementSpace *h1_fespace,
                                            const DispersiveConductivityOperator &disp_op)
{
  // Forms of a coefficient on the ND space and, for mode analysis, of its out-of-plane part
  // make_n() on the H1 space, for the entire domain and the postprocessing domains.
  auto add_form = [&](const MaterialPropertyCoefficient &coeff, auto &&make_n,
                      std::function<std::pair<double, double>(double)> weights)
  {
    auto add = [&](const FiniteElementSpace &fespace, const MaterialPropertyCoefficient &q,
                   std::unique_ptr<Operator> &M,
                   std::map<int, std::unique_ptr<Operator>> &M_i)
    {
      M = AssembleMass(fespace, q);
      for (const auto &[idx, data] : postpro.energy)
      {
        MaterialPropertyCoefficient q_i(q);
        q_i.RestrictCoefficient(mat_op.GetCeedAttributes(data.attributes));
        M_i.emplace(idx, AssembleMass(fespace, q_i));
      }
    };
    auto &form = cond_forms.emplace_back();
    add(nd_fespace, coeff, form.M, form.M_i);
    if (h1_fespace)
    {
      add(*h1_fespace, make_n(), form.M_n, form.M_n_i);
    }
    form.weights = std::move(weights);
  };
  const auto &attr_mat = mat_op.GetAttributeToMaterial();
  if (mat_op.HasLondonDepth())
  {
    has_kinetic = true;
    add_form(
        {attr_mat, mat_op.GetInvLondonDepth()}, [&]() -> MaterialPropertyCoefficient
        { return {attr_mat, mat_op.GetInvLondonDepthScalar()}; },
        [](double omega) { return std::make_pair(1.0 / (omega * omega), 0.0); });
  }
  if (mat_op.HasConductivity())
  {
    has_loss = true;
    add_form(
        {attr_mat, mat_op.GetConductivity()}, [&]() -> MaterialPropertyCoefficient
        { return {attr_mat, mat_op.GetConductivityScalar()}; },
        [](double) { return std::make_pair(0.0, 1.0); });
  }
  for (std::size_t g = 0; g < disp_op.Size(); g++)
  {
    if (disp_op.IsActive(g))
    {
      has_kinetic = has_loss = true;
      MaterialPropertyCoefficient unit(mat_op.MaxCeedAttribute());
      disp_op.AddUnitCoefficients(g, unit);
      add_form(
          unit, [&]() { return unit; },
          [&disp_op, g](double omega)
          {
            return std::make_pair(disp_op.EvaluateKineticScalar(g, omega),
                                  disp_op.EvaluateLossScalar(g, omega));
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

double DomainPostOperator::GetElectricFieldEnergy(const GridFunction &E,
                                                  const GridFunction *En) const
{
  MFEM_VERIFY(
      M_elec,
      "Domain postprocessing is not configured for electric field energy calculation!");
  double dot =
      LocalQuadratic(*M_elec, E, D) + (En ? LocalQuadratic(*M_elec_n, *En, D_n) : 0.0);
  Mpi::GlobalSum(1, &dot, E.GetComm());
  return 0.5 * dot;
}

double DomainPostOperator::GetMagneticFieldEnergy(const GridFunction &B,
                                                  const GridFunction *Bt) const
{
  if (!M_mag)
  {
    return 0.0;
  }
  double dot = LocalQuadratic(*M_mag, B, H) + (Bt ? LocalQuadratic(*M_mag_t, *Bt, D) : 0.0);
  Mpi::GlobalSum(1, &dot, B.GetComm());
  return 0.5 * dot;
}

double DomainPostOperator::GetDomainElectricFieldEnergy(int idx, const GridFunction &E,
                                                        const GridFunction *En) const
{
  // Compute the electric field energy integral for only a portion of the domain.
  auto it = M_i.find(idx);
  MFEM_VERIFY(it != M_i.end(),
              "Invalid domain index when postprocessing domain electric field energy!");
  if (!it->second.first)
  {
    return 0.0;
  }
  double dot = LocalQuadratic(*it->second.first, E, D) +
               (En ? LocalQuadratic(*M_n_i.at(idx).first, *En, D_n) : 0.0);
  Mpi::GlobalSum(1, &dot, E.GetComm());
  return 0.5 * dot;
}

double DomainPostOperator::GetDomainMagneticFieldEnergy(int idx, const GridFunction &B,
                                                        const GridFunction *Bt) const
{
  // Compute the magnetic field energy integral for only a portion of the domain.
  auto it = M_i.find(idx);
  MFEM_VERIFY(it != M_i.end(),
              "Invalid domain index when postprocessing domain magnetic field energy!");
  if (!it->second.second)
  {
    return 0.0;
  }
  double dot = LocalQuadratic(*it->second.second, B, H) +
               (Bt ? LocalQuadratic(*M_n_i.at(idx).second, *Bt, D) : 0.0);
  Mpi::GlobalSum(1, &dot, B.GetComm());
  return 0.5 * dot;
}

std::pair<double, double>
DomainPostOperator::GetKineticEnergyAndLoss(const GridFunction &E, double omega,
                                            const GridFunction *En) const
{
  return KineticEnergyAndLoss(
      cond_forms, [](const ConductionForm &form)
      { return std::make_pair(form.M.get(), form.M_n.get()); }, E, En, omega, D, D_n);
}

std::pair<double, double> DomainPostOperator::GetDomainKineticEnergyAndLoss(
    int idx, const GridFunction &E, double omega, const GridFunction *En) const
{
  MFEM_VERIFY(M_i.contains(idx),
              "Invalid domain index when postprocessing domain kinetic energy and loss!");
  return KineticEnergyAndLoss(
      cond_forms,
      [idx, En](const ConductionForm &form)
      {
        return std::make_pair(form.M_i.at(idx).get(),
                              En ? form.M_n_i.at(idx).get() : nullptr);
      },
      E, En, omega, D, D_n);
}

}  // namespace palace
