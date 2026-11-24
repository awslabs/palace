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
#include "utils/constants.hpp"
#include "utils/iodata.hpp"

namespace palace
{

namespace
{

// Local part of uᴴ M u for the real and imaginary parts of u.
double LocalDot(const Operator &M, const GridFunction &u, Vector &Mu)
{
  M.Mult(u.Real(), Mu);
  double dot = linalg::LocalDot(u.Real(), Mu);
  if (u.HasImag())
  {
    M.Mult(u.Imag(), Mu);
    dot += linalg::LocalDot(u.Imag(), Mu);
  }
  return dot;
}

}  // namespace

DomainPostOperator::DomainPostOperator(const config::DomainPostData &postpro,
                                       const MaterialOperator &mat_op,
                                       const FiniteElementSpace &nd_fespace,
                                       const FiniteElementSpace &rt_fespace,
                                       const std::vector<const Operator *> &unit_mass)
{
  // Mass operators are always partially assembled. The total energies exclude the PML
  // regions (see MaterialOperator::GetBulkAttributeToMaterial), and the energies of the
  // domains and of the PML regions use the material properties of all attributes, with the
  // background material in the PML regions (see MaterialOperator::GetAttributeToMaterial).
  const int dim = nd_fespace.Dimension();
  const auto curl_map = rt_fespace.GetFEColl().GetMapType(dim);
  MFEM_VERIFY(nd_fespace.GetFEColl().GetMapType(dim) == mfem::FiniteElement::H_CURL &&
                  (curl_map == mfem::FiniteElement::H_DIV ||
                   curl_map == mfem::FiniteElement::INTEGRAL),
              "Unexpected finite element space types for domain energy postprocessing!");
  {
    // Construct the unweighted electric quadratic form eᴴ M_eps e = Re{∫_Ω Dᴴ E dV}.
    // ElectricFieldEnergy applies EnergyWeight: 1/2 for instantaneous real fields and 1/4
    // for complex peak phasors, to its sum with the polarization forms of the
    // frequency-dependent materials. Only the real, symmetric permittivity contributes to
    // the real part of this quadratic form.
    MaterialPropertyCoefficient epsilon_func(mat_op.GetBulkAttributeToMaterial(),
                                             mat_op.GetPermittivityReal());
    BilinearForm m(nd_fespace);
    m.AddDomainIntegrator<VectorFEMassIntegrator>(epsilon_func);
    M_elec = m.PartialAssemble();
    D.SetSize(M_elec->Height());
    D.UseDevice(true);
  }
  {
    // Construct the unweighted magnetic quadratic form bᴴ M_muinv b = Re{∫_Ω Hᴴ B dV}.
    // GetFieldEnergy applies the same instantaneous or time-average weight as above.
    // In 2D, B is scalar (L2 space) so use scalar mass integrator with scalar μ⁻¹.
    // In 3D, B is a vector (RT space) so use vector FE mass integrator.
    if (curl_map == mfem::FiniteElement::INTEGRAL)
    {
      // Scalar curl (2D): need scalar μ⁻¹ (z-z component).
      // GetCurlCurlInvPermeability() returns 1x1 for 2D MaterialOperator.
      MaterialPropertyCoefficient muinv_func(mat_op.GetBulkAttributeToMaterial(),
                                             mat_op.GetCurlCurlInvPermeability());
      BilinearForm m(rt_fespace);
      m.AddDomainIntegrator<MassIntegrator>(muinv_func);
      M_mag = m.PartialAssemble();
    }
    else
    {
      MaterialPropertyCoefficient muinv_func(mat_op.GetBulkAttributeToMaterial(),
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

  // Energies in the PML regions.
  if (mat_op.HasPML())
  {
    const auto pml_attributes = mat_op.GetCeedAttributes(mat_op.GetPMLAttributes());
    {
      MaterialPropertyCoefficient epsilon_func(mat_op.GetAttributeToMaterial(),
                                               mat_op.GetPermittivityReal());
      epsilon_func.RestrictCoefficient(pml_attributes);
      BilinearForm m(nd_fespace);
      m.AddDomainIntegrator<VectorFEMassIntegrator>(epsilon_func);
      M_elec_pml = m.PartialAssemble();
    }
    {
      MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                             mat_op.GetInvPermeability());
      muinv_func.RestrictCoefficient(pml_attributes);
      BilinearForm m(rt_fespace);
      m.AddDomainIntegrator<VectorFEMassIntegrator>(muinv_func);
      M_mag_pml = m.PartialAssemble();
    }
  }
  SetUpFrequencyDependentForms(postpro, mat_op, nd_fespace, unit_mass, M_elec_fd,
                               M_elec_fd_i);
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
      MaterialPropertyCoefficient epsilon_func(mat_op.GetBulkAttributeToMaterial(),
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
      MaterialPropertyCoefficient muinv_func(mat_op.GetBulkAttributeToMaterial(),
                                             mat_op.GetCurlCurlInvPermeability());
      BilinearForm m(fespace);
      m.AddDomainIntegrator<CurlCurlIntegrator>(muinv_func);
      // Volume magnetic energy only; the London kinetic sheet energy is added explicitly in
      // the driver, so GetMagneticFieldEnergy stays the physical field energy.
      M_mag = m.PartialAssemble();
      H.SetSize(M_mag->Height());
      H.UseDevice(true);
    }
    if (mat_op.HasLondonDepth())
    {
      MaterialPropertyCoefficient f(mat_op.GetAttributeToMaterial(),
                                    mat_op.GetInvLondonDepth());
      BilinearForm m(fespace);
      m.SkipZeroCoefficientElements();
      m.AddDomainIntegrator<VectorFEMassIntegrator>(f);
      M_london = m.PartialAssemble();
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

DomainPostOperator::DomainPostOperator(const config::DomainPostData &postpro,
                                       const MaterialOperator &mat_op,
                                       const FiniteElementSpace &nd_fespace,
                                       const FiniteElementSpace &h1_fespace,
                                       const FiniteElementSpace &l2_fespace)
  : DomainPostOperator(postpro, mat_op, nd_fespace, l2_fespace)
{
  // Out-of-plane permittivity on the H1 space and in-plane inverse permeability on the ND
  // space, as in the mode eigenproblem.
  auto Assemble = [&](const FiniteElementSpace &fespace, const mfem::DenseTensor &mat,
                      const std::vector<int> *attributes)
  {
    MaterialPropertyCoefficient q(attributes ? mat_op.GetAttributeToMaterial()
                                             : mat_op.GetBulkAttributeToMaterial(),
                                  mat);
    BilinearForm m(fespace);
    if (attributes)
    {
      q.RestrictCoefficient(mat_op.GetCeedAttributes(*attributes));
      m.SkipZeroCoefficientElements();
    }
    if (fespace.GetFEColl().GetMapType(fespace.Dimension()) == mfem::FiniteElement::H_CURL)
    {
      m.AddDomainIntegrator<VectorFEMassIntegrator>(q);
    }
    else
    {
      m.AddDomainIntegrator<MassIntegrator>(q);
    }
    return m.PartialAssemble();
  };
  M_elec_n = Assemble(h1_fespace, mat_op.GetPermittivityScalar(), nullptr);
  M_mag_t = Assemble(nd_fespace, mat_op.GetInvPermeability(), nullptr);
  D_n.SetSize(M_elec_n->Height());
  D_n.UseDevice(true);
  for (const auto &[idx, data] : postpro.energy)
  {
    M_n_i.emplace(
        idx, std::make_pair(
                 Assemble(h1_fespace, mat_op.GetPermittivityScalar(), &data.attributes),
                 Assemble(nd_fespace, mat_op.GetInvPermeability(), &data.attributes)));
  }
  SetUpFrequencyDependentForms(postpro, mat_op, h1_fespace, {}, M_elec_n_fd, M_elec_n_fd_i);
}

DomainPostOperator::DomainPostOperator(const IoData &iodata, const MaterialOperator &mat_op,
                                       const FiniteElementSpace &fespace)
  : DomainPostOperator(iodata.domains.postpro, mat_op, fespace)
{
}

void DomainPostOperator::SetUpFrequencyDependentForms(
    const config::DomainPostData &postpro, const MaterialOperator &mat_op,
    const FiniteElementSpace &fespace, const std::vector<const Operator *> &unit_mass,
    FrequencyDependentForms &M_fd, std::map<int, FrequencyDependentForms> &M_fd_i)
{
  // Every rank builds the same forms for the globally supported materials, empty where a
  // rank has none of their elements.
  this->mat_op = &mat_op;
  auto MakeForms = [&](const std::vector<int> *domain_attributes)
  {
    FrequencyDependentForms forms;
    for (std::size_t m = 0; m < mat_op.NumFrequencyDependentPermittivityMaterials(); m++)
    {
      if (!mat_op.HasFrequencyDependentPermittivityModel(m))
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
        // A material outside the postprocessing domain (on every process: the attributes
        // are global) does not enter its energy, nor whether that energy is available.
        if (attributes.empty())
        {
          continue;
        }
      }
      // The form on every attribute of the material is its unit mass, borrowed if given.
      if (attributes.size() ==
              mat_op.GetFrequencyDependentPermittivityAttributes(m).size() &&
          m < unit_mass.size() && unit_mass[m])
      {
        MFEM_VERIFY(unit_mass[m]->Width() == fespace.GetVSize(),
                    "Invalid unit mass of a frequency-dependent material!");
        forms.emplace_back(m, unit_mass[m]);
        continue;
      }
      MaterialPropertyCoefficient f(mat_op.MaxCeedAttribute());
      if (!attributes.empty())
      {
        f.AddMaterialProperty(mat_op.GetCeedAttributes(attributes), 1.0);
      }
      // Only ever applied, as M_elec_i: assembled on the lossy element set, which contains
      // the material (SpaceOperator::SetUpCeedElementSets).
      BilinearForm a(fespace);
      a.SkipZeroCoefficientElements();
      if (!f.empty())
      {
        if (fespace.GetFEColl().GetMapType(fespace.Dimension()) ==
            mfem::FiniteElement::H_CURL)
        {
          a.AddDomainIntegrator<VectorFEMassIntegrator>(f);
        }
        else
        {
          a.AddDomainIntegrator<MassIntegrator>(f);
        }
      }
      forms.emplace_back(m, M_elec_fd_owned.emplace_back(a.PartialAssemble()).get());
    }
    return forms;
  };
  M_fd = MakeForms(nullptr);
  for (const auto &[idx, data] : postpro.energy)
  {
    M_fd_i.emplace(idx, MakeForms(&data.attributes));
  }
}

bool DomainPostOperator::HasStoredEnergy(const FrequencyDependentForms &M_fd) const
{
  return std::ranges::all_of(M_fd, [this](const auto &form)
                             { return mat_op->HasStoredElectricEnergy(form.first); });
}

bool DomainPostOperator::HasStoredElectricFieldEnergy() const
{
  return HasStoredEnergy(M_elec_fd);
}

bool DomainPostOperator::HasStoredDomainElectricFieldEnergy(int idx) const
{
  MFEM_VERIFY(M_i.contains(idx),
              "Invalid domain index when postprocessing domain electric field energy!");
  auto fd = M_elec_fd_i.find(idx);
  return fd == M_elec_fd_i.end() || HasStoredEnergy(fd->second);
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
                                               const GridFunction &E, const Operator *M_n,
                                               const FrequencyDependentForms *M_n_fd,
                                               const GridFunction *En, double omega,
                                               ElectricEnergy type) const
{
  auto Part = [&](const Operator &M, const FrequencyDependentForms &M_fd,
                  const GridFunction &u, Vector &Mu)
  {
    double dot = LocalDot(M, u, Mu);
    for (const auto &[m, op] : M_fd)
    {
      dot += GetFrequencyDependentEnergyWeight(m, omega, type) * LocalDot(*op, u, Mu);
    }
    return dot;
  };
  MFEM_VERIFY(!En || M_n, "Domain postprocessing is not configured for the out-of-plane "
                          "electric field energy!");
  double energy = Part(M, M_fd, E, D) + (En ? Part(*M_n, *M_n_fd, *En, D_n) : 0.0);
  Mpi::GlobalSum(1, &energy, E.GetComm());
  return EnergyWeight(E.HasImag()) * energy;
}

double DomainPostOperator::GetFieldEnergy(const Operator &M, const GridFunction &u,
                                          Vector &Mu)
{
  // Instantaneous energy 1/2 uᵀ M u for a real field, or time-averaged energy
  // 1/4 (u_rᵀ M u_r + u_iᵀ M u_i) for a complex peak phasor.
  double dot = LocalDot(M, u, Mu);
  Mpi::GlobalSum(1, &dot, u.GetComm());
  return EnergyWeight(u.HasImag()) * dot;
}

double DomainPostOperator::GetPMLElectricFieldEnergy(const GridFunction &E) const
{
  return M_elec_pml ? GetFieldEnergy(*M_elec_pml, E, D) : 0.0;
}

double DomainPostOperator::GetPMLMagneticFieldEnergy(const GridFunction &B) const
{
  return M_mag_pml ? GetFieldEnergy(*M_mag_pml, B, H) : 0.0;
}

double DomainPostOperator::GetElectricFieldEnergy(const GridFunction &E, double omega,
                                                  ElectricEnergy type,
                                                  const GridFunction *En) const
{
  MFEM_VERIFY(M_elec, "Domain postprocessing is not configured for electric field energy "
                      "calculation!");
  return ElectricFieldEnergy(*M_elec, M_elec_fd, E, M_elec_n.get(), &M_elec_n_fd, En, omega,
                             type);
}

DomainPostOperator::ElectricFieldEnergies
DomainPostOperator::GetElectricFieldEnergies(const GridFunction &E, double omega,
                                             const GridFunction *En) const
{
  MFEM_VERIFY(M_elec, "Domain postprocessing is not configured for electric field energy "
                      "calculation!");
  MFEM_VERIFY(!En || M_elec_n, "Domain postprocessing is not configured for the "
                               "out-of-plane electric field energy!");
  const bool stored = HasStoredElectricFieldEnergy();
  double energy[2] = {0.0, 0.0};  // STORED, MATERIAL
  auto Part = [&](const Operator &M, const FrequencyDependentForms &M_fd,
                  const GridFunction &u, Vector &Mu)
  {
    const double dot = LocalDot(M, u, Mu);
    energy[0] += dot;
    energy[1] += dot;
    for (const auto &[m, op] : M_fd)
    {
      const double d = LocalDot(*op, u, Mu);
      if (stored)
      {
        energy[0] +=
            GetFrequencyDependentEnergyWeight(m, omega, ElectricEnergy::STORED) * d;
      }
      energy[1] +=
          GetFrequencyDependentEnergyWeight(m, omega, ElectricEnergy::MATERIAL) * d;
    }
  };
  Part(*M_elec, M_elec_fd, E, D);
  if (En)
  {
    Part(*M_elec_n, M_elec_n_fd, *En, D_n);
  }
  Mpi::GlobalSum(2, energy, E.GetComm());
  const double weight = EnergyWeight(E.HasImag());
  return {stored ? std::optional<double>(weight * energy[0]) : std::nullopt,
          weight * energy[1]};
}

double DomainPostOperator::GetHighFrequencyElectricFieldEnergy(const GridFunction &E,
                                                               const GridFunction *En) const
{
  MFEM_VERIFY(M_elec, "Domain postprocessing is not configured for electric field energy "
                      "calculation!");
  static const FrequencyDependentForms none;
  return ElectricFieldEnergy(*M_elec, none, E, M_elec_n.get(), &none, En, 0.0,
                             ElectricEnergy::STORED);
}

double DomainPostOperator::GetMagneticFieldEnergy(const GridFunction &B,
                                                  const GridFunction *Bt) const
{
  if (!M_mag)
  {
    return 0.0;
  }
  MFEM_VERIFY(!Bt || M_mag_t, "Domain postprocessing is not configured for the in-plane "
                              "magnetic field energy!");
  double dot = LocalDot(*M_mag, B, H) + (Bt ? LocalDot(*M_mag_t, *Bt, D) : 0.0);
  Mpi::GlobalSum(1, &dot, B.GetComm());
  return EnergyWeight(B.HasImag()) * dot;
}

double DomainPostOperator::GetDomainElectricFieldEnergy(int idx, const GridFunction &E,
                                                        double omega, ElectricEnergy type,
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
  static const FrequencyDependentForms none;
  auto Forms = [idx](const std::map<int, FrequencyDependentForms> &M_fd_i)
  {
    auto fd = M_fd_i.find(idx);
    return (fd != M_fd_i.end()) ? &fd->second : &none;
  };
  const auto n = M_n_i.find(idx);
  return ElectricFieldEnergy(*it->second.first, *Forms(M_elec_fd_i), E,
                             (n != M_n_i.end()) ? n->second.first.get() : nullptr,
                             Forms(M_elec_n_fd_i), En, omega, type);
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
  double dot = LocalDot(*it->second.second, B, H);
  if (Bt)
  {
    const auto n = M_n_i.find(idx);
    MFEM_VERIFY(n != M_n_i.end(), "Domain postprocessing is not configured for the "
                                  "in-plane magnetic field energy!");
    dot += LocalDot(*n->second.second, *Bt, D);
  }
  Mpi::GlobalSum(1, &dot, B.GetComm());
  return EnergyWeight(B.HasImag()) * dot;
}

double DomainPostOperator::GetMagneticFieldEnergy(const Vector &a_l, MPI_Comm comm) const
{
  if (!M_mag)
  {
    return 0.0;
  }
  M_mag->Mult(a_l, H);
  double dot = linalg::LocalDot(a_l, H);
  Mpi::GlobalSum(1, &dot, comm);
  return EnergyWeight(false) * dot;
}

double DomainPostOperator::GetLondonKineticEnergy(const Vector &a_l, MPI_Comm comm) const
{
  if (!M_london)
  {
    return 0.0;
  }
  M_london->Mult(a_l, H);
  double dot = linalg::LocalDot(a_l, H);
  Mpi::GlobalSum(1, &dot, comm);
  return EnergyWeight(false) * dot;
}

double DomainPostOperator::GetDomainMagneticFieldEnergy(int idx, const Vector &a_l,
                                                        MPI_Comm comm) const
{
  auto it = M_i.find(idx);
  MFEM_VERIFY(it != M_i.end(),
              "Invalid domain index when postprocessing domain magnetic field energy!");
  if (!it->second.second)
  {
    return 0.0;
  }
  it->second.second->Mult(a_l, H);
  double dot = linalg::LocalDot(a_l, H);
  Mpi::GlobalSum(1, &dot, comm);
  return EnergyWeight(false) * dot;
}

}  // namespace palace
