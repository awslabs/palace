// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfacepostoperator.hpp"

#include <string>
#include "fem/integrator.hpp"
#include "models/materialoperator.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

namespace palace
{

SurfacePostOperator::InterfaceDielectricData::InterfaceDielectricData(
    const config::InterfaceDielectricData &data, mfem::ParMesh &mesh)
  : type(DielectricInterfaceType::INVALID), epsilon(0.0), ts(data.ts),
    tandelta(data.tandelta)
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
  for (const auto &node : data.nodes)
  {
    // Check inputs.
    MFEM_VERIFY(node.side.length() == 0 ||
                    (node.side.length() == 2 &&
                     (node.side[0] == '-' || node.side[0] == '+') &&
                     (node.side[1] == 'x' || node.side[1] == 'y' || node.side[1] == 'z')),
                "Postprocessing surface side is not correctly formatted!");

    // Store information about the surface side to consider.
    int component;
    mfem::Vector &side = sides.emplace_back();
    if (node.side.length() == 0)
    {
      // This is OK if surface is single sided, just push back an empty Vector.
    }
    else
    {
      side.SetSize(mesh.SpaceDimension());
      side = 0.0;
      switch (node.side[1])
      {
        case 'x':
          component = 0;
          break;
        case 'y':
          component = 1;
          break;
        case 'z':
          component = 2;
          break;
        default:
          MFEM_ABORT("Invalid side for surface boundary!");
          component = 0;  // For compiler warning
          break;
      }
      side(component) = (node.side[0] == '-') ? -1.0 : 1.0;
    }

    // Store markers for this element of the postprocessing boundary.
    mesh::AttrToMarker(mesh.bdr_attributes.Max(), node.attributes,
                       attr_markers.emplace_back());
  }
}

std::unique_ptr<mfem::Coefficient>
SurfacePostOperator::InterfaceDielectricData::GetCoefficient(
    int i, const mfem::ParGridFunction &U, const MaterialOperator &mat_op,
    const std::map<int, int> &local_to_shared) const
{
  switch (type)
  {
    case DielectricInterfaceType::MA:
      return std::make_unique<DielectricInterfaceCoefficient<DielectricInterfaceType::MA>>(
          U, mat_op, ts, epsilon, sides[i], local_to_shared);
    case DielectricInterfaceType::MS:
      return std::make_unique<DielectricInterfaceCoefficient<DielectricInterfaceType::MS>>(
          U, mat_op, ts, epsilon, sides[i], local_to_shared);
    case DielectricInterfaceType::SA:
      return std::make_unique<DielectricInterfaceCoefficient<DielectricInterfaceType::SA>>(
          U, mat_op, ts, epsilon, sides[i], local_to_shared);
    case DielectricInterfaceType::DEFAULT:
    default:
      return std::make_unique<
          DielectricInterfaceCoefficient<DielectricInterfaceType::DEFAULT>>(
          U, mat_op, ts, epsilon, sides[i], local_to_shared);
  }
}

SurfacePostOperator::SurfaceChargeData::SurfaceChargeData(
    const config::CapacitanceData &data, mfem::ParMesh &mesh)
{
  mesh::AttrToMarker(mesh.bdr_attributes.Max(), data.attributes,
                     attr_markers.emplace_back());
}

std::unique_ptr<mfem::Coefficient> SurfacePostOperator::SurfaceChargeData::GetCoefficient(
    int i, const mfem::ParGridFunction &U, const MaterialOperator &mat_op,
    const std::map<int, int> &local_to_shared) const
{
  return std::make_unique<BdrChargeCoefficient>(U, mat_op, local_to_shared);
}

SurfacePostOperator::SurfaceFluxData::SurfaceFluxData(const config::InductanceData &data,
                                                      mfem::ParMesh &mesh)
{
  // Check inputs.
  MFEM_VERIFY(data.direction.length() == 2 &&
                  (data.direction[0] == '-' || data.direction[0] == '+') &&
                  (data.direction[1] == 'x' || data.direction[1] == 'y' ||
                   data.direction[1] == 'z'),
              "Inductance postprocessing surface direction is not correctly formatted!");

  // Store information about the global direction for orientation. Note the true boundary
  // normal is used in calculating the flux, this is just used to determine the sign.
  int component;
  direction.SetSize(mesh.SpaceDimension());
  direction = 0.0;
  switch (data.direction[1])
  {
    case 'x':
      component = 0;
      break;
    case 'y':
      component = 1;
      break;
    case 'z':
      component = 2;
      break;
    default:
      MFEM_ABORT("Invalid direction for surface boundary!");
      component = 0;  // For compiler warning
      break;
  }
  direction(component) = (data.direction[0] == '-') ? -1.0 : 1.0;

  // Construct the coefficient for this postprocessing boundary (copies the direction
  // vector).
  mesh::AttrToMarker(mesh.bdr_attributes.Max(), data.attributes,
                     attr_markers.emplace_back());
}

std::unique_ptr<mfem::Coefficient> SurfacePostOperator::SurfaceFluxData::GetCoefficient(
    int i, const mfem::ParGridFunction &U, const MaterialOperator &mat_op,
    const std::map<int, int> &local_to_shared) const
{
  return std::make_unique<BdrFluxCoefficient>(U, direction, local_to_shared);
}

SurfacePostOperator::SurfacePostOperator(const IoData &iodata, const MaterialOperator &mat,
                                         const std::map<int, int> &l2s,
                                         mfem::ParFiniteElementSpace &h1_fespace)
  : mat_op(mat), local_to_shared(l2s), ones(&h1_fespace)
{
  // Define a constant 1 function on the scalar finite element space for computing surface
  // integrals.
  ones.mfem::Vector::operator=(1.0);

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

double
SurfacePostOperator::GetInterfaceElectricFieldEnergy(int idx,
                                                     const mfem::ParGridFunction &E) const
{
  auto it = eps_surfs.find(idx);
  MFEM_VERIFY(it != eps_surfs.end(),
              "Unknown dielectric loss postprocessing surface index requested!");
  return GetSurfaceIntegral(it->second, E);
}

double SurfacePostOperator::GetInterfaceLossTangent(int idx) const
{
  auto it = eps_surfs.find(idx);
  MFEM_VERIFY(it != eps_surfs.end(),
              "Unknown dielectric loss postprocessing surface index requested!");
  return it->second.tandelta;
}

double SurfacePostOperator::GetSurfaceElectricCharge(int idx,
                                                     const mfem::ParGridFunction &E) const
{
  auto it = charge_surfs.find(idx);
  MFEM_VERIFY(it != charge_surfs.end(),
              "Unknown capacitance postprocessing surface index requested!");
  return GetSurfaceIntegral(it->second, E);
}

double SurfacePostOperator::GetSurfaceMagneticFlux(int idx,
                                                   const mfem::ParGridFunction &B) const
{
  auto it = flux_surfs.find(idx);
  MFEM_VERIFY(it != flux_surfs.end(),
              "Unknown inductance postprocessing surface index requested!");
  return GetSurfaceIntegral(it->second, B);
}

double SurfacePostOperator::GetSurfaceIntegral(const SurfaceData &data,
                                               const mfem::ParGridFunction &U) const
{
  // Integrate the coefficient over the boundary attributes making up this surface index.
  std::vector<std::unique_ptr<mfem::Coefficient>> fb;
  mfem::ParLinearForm s(ones.ParFESpace());
  for (int i = 0; i < static_cast<int>(data.attr_markers.size()); i++)
  {
    fb.emplace_back(data.GetCoefficient(i, U, mat_op, local_to_shared));
    s.AddBoundaryIntegrator(new BoundaryLFIntegrator(*fb.back()), data.attr_markers[i]);
  }
  s.UseFastAssembly(false);
  s.Assemble();
  return s(ones);
}

}  // namespace palace
