// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "domainpostoperator.hpp"

#include "fem/coefficient.hpp"
#include "fem/integrator.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

DomainPostOperatorMF::DomainPostOperatorMF(const IoData &iodata,
                                           const MaterialOperator &mat,
                                           mfem::ParFiniteElementSpace &h1_fespace)
  : mat_op(mat), ones(&h1_fespace)
{
  // Define a constant 1 function on the scalar finite element space for computing volume
  // integrals.
  ones.mfem::Vector::operator=(1.0);

  // Use the provided domain postprocessing indices to group for postprocessing bulk
  // dielectric loss.
  for (const auto &[idx, data] : iodata.domains.postpro.dielectric)
  {
    mfem::Array<int> &attr_marker =
        attr_markers.emplace(idx, h1_fespace.GetParMesh()->attributes.Max()).first->second;
    attr_marker = 0;
    for (auto attr : data.attributes)
    {
      attr_marker[attr - 1] = 1;
    }
  }
}

double
DomainPostOperatorMF::GetElectricFieldEnergy(const mfem::ParComplexGridFunction &E) const
{
  // Compute the electric field energy integral as: E_elec = 1/2 Re{∫_Ω Dᴴ E dV}.
  std::map<int, int> dummy_l2s;
  EnergyDensityCoefficient<EnergyDensityType::ELECTRIC, EnergyDensityValueType::COMPLEX>
      ue_func(E, mat_op, dummy_l2s);
  return GetVolumeIntegral(ue_func);
}

double DomainPostOperatorMF::GetElectricFieldEnergy(const mfem::ParGridFunction &E) const
{
  std::map<int, int> dummy_l2s;
  EnergyDensityCoefficient<EnergyDensityType::ELECTRIC, EnergyDensityValueType::REAL>
      ue_func(E, mat_op, dummy_l2s);
  return GetVolumeIntegral(ue_func);
}

double
DomainPostOperatorMF::GetMagneticFieldEnergy(const mfem::ParComplexGridFunction &B) const
{
  // Compute the magnetic field energy integral as: E_mag = 1/2 Re{∫_Ω Bᴴ H dV}.
  std::map<int, int> dummy_l2s;
  EnergyDensityCoefficient<EnergyDensityType::MAGNETIC, EnergyDensityValueType::COMPLEX>
      um_func(B, mat_op, dummy_l2s);
  return GetVolumeIntegral(um_func);
}

double DomainPostOperatorMF::GetMagneticFieldEnergy(const mfem::ParGridFunction &B) const
{
  std::map<int, int> dummy_l2s;
  EnergyDensityCoefficient<EnergyDensityType::MAGNETIC, EnergyDensityValueType::REAL>
      um_func(B, mat_op, dummy_l2s);
  return GetVolumeIntegral(um_func);
}

double DomainPostOperatorMF::GetDomainElectricFieldEnergy(
    int idx, const mfem::ParComplexGridFunction &E) const
{
  // Compute the electric field energy integral for only a portion of the domain.
  auto it = attr_markers.find(idx);
  MFEM_VERIFY(it != attr_markers.end(),
              "Invalid domain index when postprocessing bulk dielectric loss!");
  std::map<int, int> dummy_l2s;
  EnergyDensityCoefficient<EnergyDensityType::ELECTRIC, EnergyDensityValueType::COMPLEX>
      ue_func(E, mat_op, dummy_l2s);
  return GetVolumeIntegral(ue_func, it->second);
}

double
DomainPostOperatorMF::GetDomainElectricFieldEnergy(int idx,
                                                   const mfem::ParGridFunction &E) const
{
  auto it = attr_markers.find(idx);
  MFEM_VERIFY(it != attr_markers.end(),
              "Invalid domain index when postprocessing bulk dielectric loss!");
  std::map<int, int> dummy_l2s;
  EnergyDensityCoefficient<EnergyDensityType::ELECTRIC, EnergyDensityValueType::REAL>
      ue_func(E, mat_op, dummy_l2s);
  return GetVolumeIntegral(ue_func, it->second);
}

double DomainPostOperatorMF::GetDomainElectricFieldEnergyLoss(
    int idx, const mfem::ParComplexGridFunction &E) const
{
  // Compute the electric field energy integral for only a portion of the domain.
  auto it = attr_markers.find(idx);
  MFEM_VERIFY(it != attr_markers.end(),
              "Invalid domain index when postprocessing bulk dielectric loss!");
  std::map<int, int> dummy_l2s;
  EnergyDensityCoefficient<EnergyDensityType::ELECTRIC_LOSS,
                           EnergyDensityValueType::COMPLEX>
      uei_func(E, mat_op, dummy_l2s);
  return GetVolumeIntegral(uei_func, it->second);
}

double
DomainPostOperatorMF::GetDomainElectricFieldEnergyLoss(int idx,
                                                       const mfem::ParGridFunction &E) const
{
  auto it = attr_markers.find(idx);
  MFEM_VERIFY(it != attr_markers.end(),
              "Invalid domain index when postprocessing bulk dielectric loss!");
  std::map<int, int> dummy_l2s;
  EnergyDensityCoefficient<EnergyDensityType::ELECTRIC_LOSS, EnergyDensityValueType::REAL>
      uei_func(E, mat_op, dummy_l2s);
  return GetVolumeIntegral(uei_func, it->second);
}

double DomainPostOperatorMF::GetVolumeIntegral(mfem::Coefficient &f) const
{
  // Integrate the coefficient over the entire domain.
  mfem::ParLinearForm s(ones.ParFESpace());
  s.AddDomainIntegrator(new DomainLFIntegrator(f));
  s.UseFastAssembly(true);
  s.Assemble();
  return s(ones);
}

double DomainPostOperatorMF::GetVolumeIntegral(mfem::Coefficient &f,
                                               mfem::Array<int> &attr_marker) const
{
  // Integrate the coefficient over the domain attributes making up this domain index.
  mfem::ParLinearForm s(ones.ParFESpace());
  s.AddDomainIntegrator(new DomainLFIntegrator(f), attr_marker);
  s.UseFastAssembly(true);
  s.Assemble();
  return s(ones);
}

DomainPostOperator::DomainPostOperator(const IoData &iodata, const MaterialOperator &mat_op,
                                       mfem::ParFiniteElementSpace *nd_fespace,
                                       mfem::ParFiniteElementSpace *rt_fespace)
  : m0ND(nd_fespace ? std::optional<mfem::ParBilinearForm>(nd_fespace) : std::nullopt),
    m0RT(rt_fespace ? std::optional<mfem::ParBilinearForm>(rt_fespace) : std::nullopt)
{
  if (m0ND.has_value())
  {
    // Construct ND mass matrix to compute the electric field energy integral as:
    //              E_elec = 1/2 Re{∫_Ω Dᴴ E dV} as (M_eps * e)ᴴ e.
    // Only the real part of the permeability contributes to the energy (imaginary part
    // cancels out in the inner product due to symmetry).
    MaterialPropertyCoefficient<MaterialPropertyType::PERMITTIVITY_REAL> epsilon_func(
        mat_op);
    m0ND->AddDomainIntegrator(new mfem::MixedVectorMassIntegrator(epsilon_func));
    // m0ND->SetAssemblyLevel(mfem::AssemblyLevel::FULL);
    m0ND->Assemble();
    m0ND->Finalize();

    // Use the provided domain postprocessing indices to group for postprocessing bulk
    // dielectric loss.
    int attr_max = nd_fespace->GetParMesh()->attributes.Max();
    for (const auto &[idx, data] : iodata.domains.postpro.dielectric)
    {
      mfem::Array<int> attr_marker(attr_max);
      attr_marker = 0;
      for (auto attr : data.attributes)
      {
        attr_marker[attr - 1] = 1;
      }
      SumMatrixCoefficient epsilon_func_r(nd_fespace->GetParMesh()->SpaceDimension());
      SumMatrixCoefficient epsilon_func_i(nd_fespace->GetParMesh()->SpaceDimension());
      epsilon_func_r.AddCoefficient(
          std::make_unique<
              MaterialPropertyCoefficient<MaterialPropertyType::PERMITTIVITY_REAL>>(mat_op),
          attr_marker);
      epsilon_func_i.AddCoefficient(
          std::make_unique<
              MaterialPropertyCoefficient<MaterialPropertyType::PERMITTIVITY_IMAG>>(mat_op,
                                                                                    -1.0),
          attr_marker);
      auto &m0 = m0NDi.emplace(idx, std::make_pair(nd_fespace, nd_fespace)).first->second;
      mfem::ParBilinearForm &m0r = m0.first;
      mfem::ParBilinearForm &m0i = m0.second;
      m0r.AddDomainIntegrator(new mfem::MixedVectorMassIntegrator(epsilon_func_r));
      m0i.AddDomainIntegrator(new mfem::MixedVectorMassIntegrator(epsilon_func_i));
      // m0r.SetAssemblyLevel(mfem::AssemblyLevel::FULL);
      // m0i.SetAssemblyLevel(mfem::AssemblyLevel::FULL);
      m0r.Assemble();
      m0i.Assemble();
      m0r.Finalize();
      m0i.Finalize();
    }
  }

  if (m0RT.has_value())
  {
    // Construct RT mass matrix to compute the magnetic field energy integral as:
    //              E_mag = 1/2 Re{∫_Ω Bᴴ H dV} as (M_muinv * b)ᴴ b.
    MaterialPropertyCoefficient<MaterialPropertyType::INV_PERMEABILITY> muinv_func(mat_op);
    m0RT->AddDomainIntegrator(new mfem::MixedVectorMassIntegrator(muinv_func));
    // m0RT->SetAssemblyLevel(mfem::AssemblyLevel::FULL);
    m0RT->Assemble();
    m0RT->Finalize();
  }
}

double
DomainPostOperator::GetElectricFieldEnergy(const mfem::ParComplexGridFunction &E) const
{
  if (m0ND.has_value())
  {
    double res = m0ND->InnerProduct(E.real(), E.real());
    res += m0ND->InnerProduct(E.imag(), E.imag());
    Mpi::GlobalSum(1, &res, E.ParFESpace()->GetComm());
    return 0.5 * res;
  }
  MFEM_ABORT(
      "Domain postprocessing is not configured for electric field energy calculation!");
  return 0.0;
}

double DomainPostOperator::GetElectricFieldEnergy(const mfem::ParGridFunction &E) const
{
  if (m0ND.has_value())
  {
    double res = m0ND->InnerProduct(E, E);
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
  if (m0RT.has_value())
  {
    double res = m0RT->InnerProduct(B.real(), B.real());
    res += m0RT->InnerProduct(B.imag(), B.imag());
    Mpi::GlobalSum(1, &res, B.ParFESpace()->GetComm());
    return 0.5 * res;
  }
  MFEM_ABORT(
      "Domain postprocessing is not configured for magnetic field energy calculation!");
  return 0.0;
}

double DomainPostOperator::GetMagneticFieldEnergy(const mfem::ParGridFunction &B) const
{
  if (m0RT.has_value())
  {
    double res = m0RT->InnerProduct(B, B);
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
  auto it = m0NDi.find(idx);
  MFEM_VERIFY(it != m0NDi.end(),
              "Invalid domain index when postprocessing bulk dielectric loss!");
  double res = it->second.first.InnerProduct(E.real(), E.real());
  res += it->second.first.InnerProduct(E.imag(), E.imag());
  Mpi::GlobalSum(1, &res, E.ParFESpace()->GetComm());
  return 0.5 * res;
}

double
DomainPostOperator::GetDomainElectricFieldEnergy(int idx,
                                                 const mfem::ParGridFunction &E) const
{
  auto it = m0NDi.find(idx);
  MFEM_VERIFY(it != m0NDi.end(),
              "Invalid domain index when postprocessing bulk dielectric loss!");
  double res = it->second.first.InnerProduct(E, E);
  Mpi::GlobalSum(1, &res, E.ParFESpace()->GetComm());
  return 0.5 * res;
}

double DomainPostOperator::GetDomainElectricFieldEnergyLoss(
    int idx, const mfem::ParComplexGridFunction &E) const
{
  // Compute the electric field energy integral for only a portion of the domain.
  auto it = m0NDi.find(idx);
  MFEM_VERIFY(it != m0NDi.end(),
              "Invalid domain index when postprocessing bulk dielectric loss!");
  double res = it->second.second.InnerProduct(E.real(), E.real());
  res += it->second.second.InnerProduct(E.imag(), E.imag());
  Mpi::GlobalSum(1, &res, E.ParFESpace()->GetComm());
  return 0.5 * res;
}

double
DomainPostOperator::GetDomainElectricFieldEnergyLoss(int idx,
                                                     const mfem::ParGridFunction &E) const
{
  auto it = m0NDi.find(idx);
  MFEM_VERIFY(it != m0NDi.end(),
              "Invalid domain index when postprocessing bulk dielectric loss!");
  double res = it->second.second.InnerProduct(E, E);
  Mpi::GlobalSum(1, &res, E.ParFESpace()->GetComm());
  return 0.5 * res;
}

}  // namespace palace
