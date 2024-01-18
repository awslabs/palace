// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfacepostoperator.hpp"

#include <complex>
#include "fem/integrator.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

namespace palace
{

SurfacePostOperator::InterfaceDielectricData::InterfaceDielectricData(
    const config::InterfaceDielectricData &data, const mfem::ParMesh &mesh)
  : ts(data.ts), tandelta(data.tandelta)
{
  // Calculate surface dielectric loss according to the formulas from J. Wenner et al.,
  // Surface loss simulations of superconducting coplanar waveguide resonators, Appl. Phys.
  // Lett. (2011). If only a general layer permittivity is specified and not any special
  // metal-air (MA), metal-substrate (MS), or substrate-air (SA) permittivity, compute the
  // numerator of the participation ratio according to the regular formula
  //                       p * E_elec = 1/2 t Re{∫ (ε E)ᴴ E_m dS} .
  bool has_eps = (std::abs(data.epsilon_r) > 0.0);
  bool has_eps_ma = (std::abs(data.epsilon_r_ma) > 0.0);
  bool has_eps_ms = (std::abs(data.epsilon_r_ms) > 0.0);
  bool has_eps_sa = (std::abs(data.epsilon_r_sa) > 0.0);
  MFEM_VERIFY(has_eps + has_eps_ma + has_eps_ms + has_eps_sa == 1,
              "Surface dielectric loss postprocessing should only be specialized as one of "
              "metal-air, metal-substrate, or substrate-air, or not specialized at all!");
  if (has_eps)
  {
    type = DielectricInterfaceType::DEFAULT;
    epsilon = data.epsilon_r;
  }
  else if (has_eps_ma)
  {
    type = DielectricInterfaceType::MA;
    epsilon = data.epsilon_r_ma;
  }
  else if (has_eps_ms)
  {
    type = DielectricInterfaceType::MS;
    epsilon = data.epsilon_r_ms;
  }
  else if (has_eps_sa)
  {
    type = DielectricInterfaceType::SA;
    epsilon = data.epsilon_r_sa;
  }
  MFEM_VERIFY(data.ts > 0.0,
              "Surface dielectric loss postprocessing requires positive thickness!");

  // Construct the postprocessing data allowing for multiple groups of attribute with
  // different side values.
  for (const auto &elem : data.elements)
  {
    // Store information about the surface side to consider.
    mfem::Vector &side = sides.emplace_back();
    if (elem.direction[0] == 0 && elem.direction[1] == 0 && elem.direction[2] == 0)
    {
      // This is OK if surface is single sided, just push back an empty Vector.
    }
    else
    {
      side.SetSize(mesh.SpaceDimension());
      std::copy(elem.direction.begin(), elem.direction.end(), side.begin());
      side /= side.Norml2();
    }

    // Store boundary attributes for this element of the postprocessing boundary.
    auto &attr_list = attr_lists.emplace_back();
    attr_list.Append(elem.attributes.data(), elem.attributes.size());
  }
}

std::unique_ptr<mfem::Coefficient>
SurfacePostOperator::InterfaceDielectricData::GetCoefficient(
    std::size_t i, const mfem::ParGridFunction &U, const MaterialOperator &mat_op) const
{
  switch (type)
  {
    case DielectricInterfaceType::MA:
      return std::make_unique<RestrictedCoefficient<
          DielectricInterfaceCoefficient<DielectricInterfaceType::MA>>>(
          attr_lists[i], U, mat_op, ts, epsilon, sides[i]);
    case DielectricInterfaceType::MS:
      return std::make_unique<RestrictedCoefficient<
          DielectricInterfaceCoefficient<DielectricInterfaceType::MS>>>(
          attr_lists[i], U, mat_op, ts, epsilon, sides[i]);
    case DielectricInterfaceType::SA:
      return std::make_unique<RestrictedCoefficient<
          DielectricInterfaceCoefficient<DielectricInterfaceType::SA>>>(
          attr_lists[i], U, mat_op, ts, epsilon, sides[i]);
    case DielectricInterfaceType::DEFAULT:
      return std::make_unique<RestrictedCoefficient<
          DielectricInterfaceCoefficient<DielectricInterfaceType::DEFAULT>>>(
          attr_lists[i], U, mat_op, ts, epsilon, sides[i]);
  }
  return {};  // For compiler warning
}

SurfacePostOperator::SurfaceChargeData::SurfaceChargeData(
    const config::CapacitanceData &data, const mfem::ParMesh &mesh)
{
  // Store boundary attributes for this element of the postprocessing boundary.
  auto &attr_list = attr_lists.emplace_back();
  attr_list.Append(data.attributes.data(), data.attributes.size());
}

std::unique_ptr<mfem::Coefficient> SurfacePostOperator::SurfaceChargeData::GetCoefficient(
    std::size_t i, const mfem::ParGridFunction &U, const MaterialOperator &mat_op) const
{
  return std::make_unique<RestrictedCoefficient<BdrChargeCoefficient>>(attr_lists[i], U,
                                                                       mat_op);
}

SurfacePostOperator::SurfaceFluxData::SurfaceFluxData(const config::InductanceData &data,
                                                      const mfem::ParMesh &mesh)
{
  // Store information about the global direction for orientation. Note the true boundary
  // normal is used in calculating the flux, this is just used to determine the sign.
  direction.SetSize(mesh.SpaceDimension());
  std::copy(data.direction.begin(), data.direction.end(), direction.begin());
  direction /= direction.Norml2();

  // Store boundary attributes for this element of the postprocessing boundary.
  auto &attr_list = attr_lists.emplace_back();
  attr_list.Append(data.attributes.data(), data.attributes.size());
}

std::unique_ptr<mfem::Coefficient> SurfacePostOperator::SurfaceFluxData::GetCoefficient(
    std::size_t i, const mfem::ParGridFunction &U, const MaterialOperator &mat_op) const
{
  return std::make_unique<RestrictedCoefficient<BdrFluxCoefficient>>(attr_lists[i], U,
                                                                     mat_op, direction);
}

SurfacePostOperator::SurfacePostOperator(const IoData &iodata,
                                         const MaterialOperator &mat_op,
                                         mfem::ParFiniteElementSpace &h1_fespace)
  : mat_op(mat_op), ones(&h1_fespace)
{
  // Define a constant 1 function on the scalar finite element space for computing surface
  // integrals.
  ones = 1.0;

  // Surface dielectric loss postprocessing.
  for (const auto &[idx, data] : iodata.boundaries.postpro.dielectric)
  {
    eps_surfs.try_emplace(idx, data, *h1_fespace.GetParMesh());
  }

  // Surface capacitance postprocessing.
  for (const auto &[idx, data] : iodata.boundaries.postpro.capacitance)
  {
    charge_surfs.try_emplace(idx, data, *h1_fespace.GetParMesh());
  }

  // Surface inductance postprocessing.
  for (const auto &[idx, data] : iodata.boundaries.postpro.inductance)
  {
    flux_surfs.try_emplace(idx, data, *h1_fespace.GetParMesh());
  }
}

double SurfacePostOperator::GetInterfaceLossTangent(int idx) const
{
  auto it = eps_surfs.find(idx);
  MFEM_VERIFY(it != eps_surfs.end(),
              "Unknown dielectric loss postprocessing surface index requested!");
  return it->second.tandelta;
}

double SurfacePostOperator::GetInterfaceElectricFieldEnergy(
    int idx, const mfem::ParComplexGridFunction &E) const
{
  auto it = eps_surfs.find(idx);
  MFEM_VERIFY(it != eps_surfs.end(),
              "Unknown dielectric loss postprocessing surface index requested!");
  double dot = GetLocalSurfaceIntegral(it->second, E.real()) +
               GetLocalSurfaceIntegral(it->second, E.imag());
  Mpi::GlobalSum(1, &dot, E.ParFESpace()->GetComm());
  return dot;
}

double
SurfacePostOperator::GetInterfaceElectricFieldEnergy(int idx,
                                                     const mfem::ParGridFunction &E) const
{
  auto it = eps_surfs.find(idx);
  MFEM_VERIFY(it != eps_surfs.end(),
              "Unknown dielectric loss postprocessing surface index requested!");
  double dot = GetLocalSurfaceIntegral(it->second, E);
  Mpi::GlobalSum(1, &dot, E.ParFESpace()->GetComm());
  return dot;
}

double
SurfacePostOperator::GetSurfaceElectricCharge(int idx,
                                              const mfem::ParComplexGridFunction &E) const
{
  auto it = charge_surfs.find(idx);
  MFEM_VERIFY(it != charge_surfs.end(),
              "Unknown capacitance postprocessing surface index requested!");
  std::complex<double> dot(GetLocalSurfaceIntegral(it->second, E.real()),
                           GetLocalSurfaceIntegral(it->second, E.imag()));
  Mpi::GlobalSum(1, &dot, E.ParFESpace()->GetComm());
  return std::copysign(std::abs(dot), dot.real());
}

double SurfacePostOperator::GetSurfaceElectricCharge(int idx,
                                                     const mfem::ParGridFunction &E) const
{
  auto it = charge_surfs.find(idx);
  MFEM_VERIFY(it != charge_surfs.end(),
              "Unknown capacitance postprocessing surface index requested!");
  double dot = GetLocalSurfaceIntegral(it->second, E);
  Mpi::GlobalSum(1, &dot, E.ParFESpace()->GetComm());
  return dot;
}

double
SurfacePostOperator::GetSurfaceMagneticFlux(int idx,
                                            const mfem::ParComplexGridFunction &B) const
{
  auto it = flux_surfs.find(idx);
  MFEM_VERIFY(it != flux_surfs.end(),
              "Unknown inductance postprocessing surface index requested!");
  std::complex<double> dot(GetLocalSurfaceIntegral(it->second, B.real()),
                           GetLocalSurfaceIntegral(it->second, B.imag()));
  Mpi::GlobalSum(1, &dot, B.ParFESpace()->GetComm());
  return std::copysign(std::abs(dot), dot.real());
}

double SurfacePostOperator::GetSurfaceMagneticFlux(int idx,
                                                   const mfem::ParGridFunction &B) const
{
  auto it = flux_surfs.find(idx);
  MFEM_VERIFY(it != flux_surfs.end(),
              "Unknown inductance postprocessing surface index requested!");
  double dot = GetLocalSurfaceIntegral(it->second, B);
  Mpi::GlobalSum(1, &dot, B.ParFESpace()->GetComm());
  return dot;
}

double SurfacePostOperator::GetLocalSurfaceIntegral(const SurfaceData &data,
                                                    const mfem::ParGridFunction &U) const
{
  // Integrate the coefficient over the boundary attributes making up this surface index.
  const auto &mesh = *U.ParFESpace()->GetParMesh();
  SumCoefficient fb;
  mfem::Array<int> attr_list;
  for (std::size_t i = 0; i < data.attr_lists.size(); i++)
  {
    fb.AddCoefficient(data.GetCoefficient(i, U, mat_op));
    attr_list.Append(data.attr_lists[i]);
  }
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> attr_marker = mesh::AttrToMarker(bdr_attr_max, attr_list);
  mfem::LinearForm s(ones.FESpace());
  s.AddBoundaryIntegrator(new BoundaryLFIntegrator(fb), attr_marker);
  s.UseFastAssembly(false);
  s.Assemble();
  return s * ones;
}

}  // namespace palace
