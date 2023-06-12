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

DomainPostOperator::DomainPostOperator(const IoData &iodata, const MaterialOperator &mat_op,
                                       mfem::ParFiniteElementSpace *nd_fespace,
                                       mfem::ParFiniteElementSpace *rt_fespace)
  : M_ND(nd_fespace ? std::optional<mfem::BilinearForm>(nd_fespace) : std::nullopt),
    M_RT(rt_fespace ? std::optional<mfem::BilinearForm>(rt_fespace) : std::nullopt)
{
  if (M_ND.has_value())
  {
    // Construct ND mass matrix to compute the electric field energy integral as:
    //              E_elec = 1/2 Re{∫_Ω Dᴴ E dV} as (M_eps * e)ᴴ e.
    // Only the real part of the permeability contributes to the energy (imaginary part
    // cancels out in the inner product due to symmetry).
    constexpr MaterialPropertyType MatTypeEpsReal = MaterialPropertyType::PERMITTIVITY_REAL;
    constexpr MaterialPropertyType MatTypeEpsImag = MaterialPropertyType::PERMITTIVITY_IMAG;
    MaterialPropertyCoefficient<MatTypeEpsReal> epsilon_func(mat_op);
    M_ND->AddDomainIntegrator(new mfem::MixedVectorMassIntegrator(epsilon_func));
    // XX TODO: Partial assembly option?
    M_ND->SetAssemblyLevel(mfem::AssemblyLevel::LEGACY);
    M_ND->Assemble(0);
    M_ND->Finalize(0);
    D.SetSize(M_ND->Height());

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
          std::make_unique<MaterialPropertyCoefficient<MatTypeEpsReal>>(mat_op),
          attr_marker);
      epsilon_func_i.AddCoefficient(
          std::make_unique<MaterialPropertyCoefficient<MatTypeEpsImag>>(mat_op, -1.0),
          attr_marker);
      auto &M = M_NDi.emplace(idx, std::make_pair(nd_fespace, nd_fespace)).first->second;
      mfem::BilinearForm &Mr = M.first;
      mfem::BilinearForm &Mi = M.second;
      Mr.AddDomainIntegrator(new mfem::MixedVectorMassIntegrator(epsilon_func_r));
      Mi.AddDomainIntegrator(new mfem::MixedVectorMassIntegrator(epsilon_func_i));
      // XX TODO: Partial assembly option?
      Mr.SetAssemblyLevel(mfem::AssemblyLevel::LEGACY);
      Mi.SetAssemblyLevel(mfem::AssemblyLevel::LEGACY);
      Mr.Assemble();
      Mi.Assemble();
      Mr.Finalize();
      Mi.Finalize();
    }
  }

  if (M_RT.has_value())
  {
    // Construct RT mass matrix to compute the magnetic field energy integral as:
    //              E_mag = 1/2 Re{∫_Ω Bᴴ H dV} as (M_muinv * b)ᴴ b.
    constexpr MaterialPropertyType MatTypeMuInv = MaterialPropertyType::INV_PERMEABILITY;
    MaterialPropertyCoefficient<MatTypeMuInv> muinv_func(mat_op);
    M_RT->AddDomainIntegrator(new mfem::MixedVectorMassIntegrator(muinv_func));
    // XX TODO: Partial assembly option?
    M_RT->SetAssemblyLevel(mfem::AssemblyLevel::LEGACY);
    M_RT->Assemble(0);
    M_RT->Finalize(0);
    H.SetSize(M_RT->Height());
  }
}

double
DomainPostOperator::GetElectricFieldEnergy(const mfem::ParComplexGridFunction &E) const
{
  if (M_ND.has_value())
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
  if (M_ND.has_value())
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
  if (M_RT.has_value())
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
  if (M_RT.has_value())
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
  auto it = M_NDi.find(idx);
  MFEM_VERIFY(it != M_NDi.end(),
              "Invalid domain index when postprocessing bulk dielectric loss!");
  it->second.first.Mult(E.real(), D);
  double res = mfem::InnerProduct(E.real(), D);
  it->second.first.Mult(E.imag(), D);
  res += mfem::InnerProduct(E.imag(), D);
  Mpi::GlobalSum(1, &res, E.ParFESpace()->GetComm());
  return 0.5 * res;
}

double
DomainPostOperator::GetDomainElectricFieldEnergy(int idx,
                                                 const mfem::ParGridFunction &E) const
{
  auto it = M_NDi.find(idx);
  MFEM_VERIFY(it != M_NDi.end(),
              "Invalid domain index when postprocessing bulk dielectric loss!");
  it->second.first.Mult(E, D);
  double res = mfem::InnerProduct(E, D);
  Mpi::GlobalSum(1, &res, E.ParFESpace()->GetComm());
  return 0.5 * res;
}

double DomainPostOperator::GetDomainElectricFieldEnergyLoss(
    int idx, const mfem::ParComplexGridFunction &E) const
{
  // Compute the electric field energy integral for only a portion of the domain.
  auto it = M_NDi.find(idx);
  MFEM_VERIFY(it != M_NDi.end(),
              "Invalid domain index when postprocessing bulk dielectric loss!");
  it->second.second.Mult(E.real(), D);
  double res = mfem::InnerProduct(E.real(), D);
  it->second.second.Mult(E.imag(), D);
  res += mfem::InnerProduct(E.imag(), D);
  Mpi::GlobalSum(1, &res, E.ParFESpace()->GetComm());
  return 0.5 * res;
}

double
DomainPostOperator::GetDomainElectricFieldEnergyLoss(int idx,
                                                     const mfem::ParGridFunction &E) const
{
  auto it = M_NDi.find(idx);
  MFEM_VERIFY(it != M_NDi.end(),
              "Invalid domain index when postprocessing bulk dielectric loss!");
  it->second.second.Mult(E, D);
  double res = mfem::InnerProduct(E, D);
  Mpi::GlobalSum(1, &res, E.ParFESpace()->GetComm());
  return 0.5 * res;
}

}  // namespace palace
