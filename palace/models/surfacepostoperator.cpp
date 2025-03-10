// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfacepostoperator.hpp"

#include <complex>
#include <set>
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "linalg/vector.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"

namespace palace
{

namespace
{

template <typename T>
mfem::Array<int> SetUpBoundaryProperties(const T &data,
                                         const mfem::Array<int> &bdr_attr_marker)
{
  mfem::Array<int> attr_list;
  attr_list.Reserve(static_cast<int>(data.attributes.size()));
  std::set<int> bdr_warn_list;
  for (auto attr : data.attributes)
  {
    // MFEM_VERIFY(attr > 0 && attr <= bdr_attr_max,
    //             "Boundary postprocessing attribute tags must be non-negative and "
    //             "correspond to attributes in the mesh!");
    // MFEM_VERIFY(bdr_attr_marker[attr - 1],
    //             "Unknown boundary postprocessing attribute " << attr << "!");
    if (attr <= 0 || attr > bdr_attr_marker.Size() || !bdr_attr_marker[attr - 1])
    {
      bdr_warn_list.insert(attr);
    }
    else
    {
      attr_list.Append(attr);
    }
  }
  if (!bdr_warn_list.empty())
  {
    Mpi::Print("\n");
    Mpi::Warning(
        "Unknown boundary postprocessing attributes!\nSolver will just ignore them!");
    utils::PrettyPrint(bdr_warn_list, "Boundary attribute list:");
    Mpi::Print("\n");
  }
  return attr_list;
}

}  // namespace

SurfacePostOperator::SurfaceFluxData::SurfaceFluxData(
    const config::SurfaceFluxData &data, const mfem::ParMesh &mesh,
    const mfem::Array<int> &bdr_attr_marker)
{
  // Store boundary attributes for this postprocessing boundary.
  attr_list = SetUpBoundaryProperties(data, bdr_attr_marker);

  // Store the type of flux.
  switch (data.type)
  {
    case config::SurfaceFluxData::Type::ELECTRIC:
      type = SurfaceFluxType::ELECTRIC;
      break;
    case config::SurfaceFluxData::Type::MAGNETIC:
      type = SurfaceFluxType::MAGNETIC;
      break;
    case config::SurfaceFluxData::Type::POWER:
      type = SurfaceFluxType::POWER;
      break;
  }

  // Store information about the global direction for orientation. Note the true boundary
  // normal is used in calculating the flux, this is just used to determine the sign.
  two_sided = data.two_sided;
  if (!two_sided)
  {
    center.SetSize(mesh.SpaceDimension());
    if (data.no_center)
    {
      // Compute the center as the bounding box centroid for all boundary elements making up
      // this postprocessing boundary.
      mfem::Vector bbmin, bbmax;
      mesh::GetAxisAlignedBoundingBox(
          mesh, mesh::AttrToMarker(bdr_attr_marker.Size(), attr_list), true, bbmin, bbmax);
      for (int d = 0; d < mesh.SpaceDimension(); d++)
      {
        center(d) = 0.5 * (bbmin(d) + bbmax(d));
      }
    }
    else
    {
      std::copy(data.center.begin(), data.center.end(), center.begin());
    }
  }
}

std::unique_ptr<mfem::Coefficient>
SurfacePostOperator::SurfaceFluxData::GetCoefficient(const mfem::ParGridFunction *E,
                                                     const mfem::ParGridFunction *B,
                                                     const MaterialOperator &mat_op) const
{
  switch (type)
  {
    case SurfaceFluxType::ELECTRIC:
      return std::make_unique<
          RestrictedCoefficient<BdrSurfaceFluxCoefficient<SurfaceFluxType::ELECTRIC>>>(
          attr_list, E, nullptr, mat_op, two_sided, center);
    case SurfaceFluxType::MAGNETIC:
      return std::make_unique<
          RestrictedCoefficient<BdrSurfaceFluxCoefficient<SurfaceFluxType::MAGNETIC>>>(
          attr_list, nullptr, B, mat_op, two_sided, center);
    case SurfaceFluxType::POWER:
      return std::make_unique<
          RestrictedCoefficient<BdrSurfaceFluxCoefficient<SurfaceFluxType::POWER>>>(
          attr_list, E, B, mat_op, two_sided, center);
  }
  return {};
}

SurfacePostOperator::InterfaceDielectricData::InterfaceDielectricData(
    const config::InterfaceDielectricData &data, const mfem::ParMesh &mesh,
    const mfem::Array<int> &bdr_attr_marker)
{
  // Store boundary attributes for this postprocessing boundary.
  attr_list = SetUpBoundaryProperties(data, bdr_attr_marker);

  // Calculate surface dielectric loss according to the formulas from J. Wenner et al.,
  // Surface loss simulations of superconducting coplanar waveguide resonators, Appl. Phys.
  // Lett. (2011). If only a general layer permittivity is specified and not any special
  // metal-air (MA), metal-substrate (MS), or substrate-air (SA) permittivity, compute the
  // numerator of the participation ratio according to the regular formula
  //                       p * E_elec = 1/2 t Re{∫ (ε E)ᴴ E_m dS} .
  switch (data.type)
  {
    case config::InterfaceDielectricData::Type::DEFAULT:
      type = InterfaceDielectricType::DEFAULT;
      break;
    case config::InterfaceDielectricData::Type::MA:
      type = InterfaceDielectricType::MA;
      break;
    case config::InterfaceDielectricData::Type::MS:
      type = InterfaceDielectricType::MS;
      break;
    case config::InterfaceDielectricData::Type::SA:
      type = InterfaceDielectricType::SA;
      break;
  }
  t = data.t;
  epsilon = data.epsilon_r;
  tandelta = data.tandelta;
}

std::unique_ptr<mfem::Coefficient>
SurfacePostOperator::InterfaceDielectricData::GetCoefficient(
    const GridFunction &E, const MaterialOperator &mat_op) const
{
  switch (type)
  {
    case InterfaceDielectricType::DEFAULT:
      return std::make_unique<RestrictedCoefficient<
          InterfaceDielectricCoefficient<InterfaceDielectricType::DEFAULT>>>(
          attr_list, E, mat_op, t, epsilon);
    case InterfaceDielectricType::MA:
      return std::make_unique<RestrictedCoefficient<
          InterfaceDielectricCoefficient<InterfaceDielectricType::MA>>>(attr_list, E,
                                                                        mat_op, t, epsilon);
    case InterfaceDielectricType::MS:
      return std::make_unique<RestrictedCoefficient<
          InterfaceDielectricCoefficient<InterfaceDielectricType::MS>>>(attr_list, E,
                                                                        mat_op, t, epsilon);
    case InterfaceDielectricType::SA:
      return std::make_unique<RestrictedCoefficient<
          InterfaceDielectricCoefficient<InterfaceDielectricType::SA>>>(attr_list, E,
                                                                        mat_op, t, epsilon);
  }
  return {};  // For compiler warning
}

SurfacePostOperator::SurfacePostOperator(const IoData &iodata,
                                         const MaterialOperator &mat_op,
                                         mfem::ParFiniteElementSpace &h1_fespace)
  : mat_op(mat_op), h1_fespace(h1_fespace)
{
  // Check that boundary attributes have been specified correctly.
  const auto &mesh = *h1_fespace.GetParMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> bdr_attr_marker;
  if (!iodata.boundaries.postpro.flux.empty() ||
      !iodata.boundaries.postpro.dielectric.empty())
  {
    bdr_attr_marker.SetSize(bdr_attr_max);
    bdr_attr_marker = 0;
    for (auto attr : mesh.bdr_attributes)
    {
      bdr_attr_marker[attr - 1] = 1;
    }
  }

  // Surface flux postprocessing.
  for (const auto &[idx, data] : iodata.boundaries.postpro.flux)
  {
    MFEM_VERIFY(iodata.problem.type != config::ProblemData::Type::ELECTROSTATIC ||
                    data.type == config::SurfaceFluxData::Type::ELECTRIC,
                "Magnetic field or power surface flux postprocessing are not available "
                "for electrostatic problems!");
    MFEM_VERIFY(iodata.problem.type != config::ProblemData::Type::MAGNETOSTATIC ||
                    data.type == config::SurfaceFluxData::Type::MAGNETIC,
                "Electric field or power surface flux postprocessing are not available "
                "for magnetostatic problems!");
    flux_surfs.try_emplace(idx, data, *h1_fespace.GetParMesh(), bdr_attr_marker);
  }

  // Interface dielectric postprocessing.
  MFEM_VERIFY(iodata.boundaries.postpro.dielectric.empty() ||
                  iodata.problem.type != config::ProblemData::Type::MAGNETOSTATIC,
              "Interface dielectric loss postprocessing is not available for "
              "magnetostatic problems!");
  for (const auto &[idx, data] : iodata.boundaries.postpro.dielectric)
  {
    eps_surfs.try_emplace(idx, data, *h1_fespace.GetParMesh(), bdr_attr_marker);
  }
}

std::complex<double> SurfacePostOperator::GetSurfaceFlux(int idx, const GridFunction *E,
                                                         const GridFunction *B) const
{
  // For complex-valued fields, output the separate real and imaginary parts for the time-
  // harmonic quantity. For power flux (Poynting vector), output only the stationary real
  // part and not the part which has double the frequency.
  auto it = flux_surfs.find(idx);
  MFEM_VERIFY(it != flux_surfs.end(),
              "Unknown surface flux postprocessing index requested!");
  const bool has_imag = (E) ? E->HasImag() : B->HasImag();
  const auto &mesh = *h1_fespace.GetParMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> attr_marker = mesh::AttrToMarker(bdr_attr_max, it->second.attr_list);
  auto f =
      it->second.GetCoefficient(E ? &E->Real() : nullptr, B ? &B->Real() : nullptr, mat_op);
  std::complex<double> dot(GetLocalSurfaceIntegral(*f, attr_marker), 0.0);
  if (has_imag)
  {
    f = it->second.GetCoefficient(E ? &E->Imag() : nullptr, B ? &B->Imag() : nullptr,
                                  mat_op);
    double doti = GetLocalSurfaceIntegral(*f, attr_marker);
    if (it->second.type == SurfaceFluxType::POWER)
    {
      dot += doti;
    }
    else
    {
      dot.imag(doti);
    }
  }
  Mpi::GlobalSum(1, &dot, (E) ? E->GetComm() : B->GetComm());
  return dot;
}

double SurfacePostOperator::GetInterfaceLossTangent(int idx) const
{
  auto it = eps_surfs.find(idx);
  MFEM_VERIFY(it != eps_surfs.end(),
              "Unknown interface dielectric postprocessing index requested!");
  return it->second.tandelta;
}

double SurfacePostOperator::GetInterfaceElectricFieldEnergy(int idx,
                                                            const GridFunction &E) const
{
  auto it = eps_surfs.find(idx);
  MFEM_VERIFY(it != eps_surfs.end(),
              "Unknown interface dielectric postprocessing index requested!");
  const auto &mesh = *h1_fespace.GetParMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> attr_marker = mesh::AttrToMarker(bdr_attr_max, it->second.attr_list);
  auto f = it->second.GetCoefficient(E, mat_op);
  double dot = GetLocalSurfaceIntegral(*f, attr_marker);
  Mpi::GlobalSum(1, &dot, E.GetComm());
  return dot;
}

double
SurfacePostOperator::GetLocalSurfaceIntegral(mfem::Coefficient &f,
                                             const mfem::Array<int> &attr_marker) const
{
  // Integrate the coefficient over the boundary attributes making up this surface index.
  mfem::LinearForm s(&h1_fespace);
  s.AddBoundaryIntegrator(new BoundaryLFIntegrator(f),
                          const_cast<mfem::Array<int> &>(attr_marker));
  s.UseFastAssembly(false);
  s.UseDevice(false);
  s.Assemble();
  s.UseDevice(true);
  return linalg::LocalSum(s);
}

}  // namespace palace
