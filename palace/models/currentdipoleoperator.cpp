// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "currentdipoleoperator.hpp"

#include "fem/coefficient.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

namespace palace
{

CurrentDipoleData::CurrentDipoleData(const config::CurrentDipoleData &data,
                                     const mfem::ParMesh &mesh)
{
  // Set up dipole moment vector
  moment.SetSize(mesh.SpaceDimension());
  MFEM_VERIFY(data.moment.size() == static_cast<std::size_t>(mesh.SpaceDimension()),
              "Current dipole moment dimension must match mesh space dimension!");
  for (int d = 0; d < mesh.SpaceDimension(); d++)
  {
    moment[d] = data.moment[d];
  }

  // Set up dipole center position
  center.SetSize(mesh.SpaceDimension());
  MFEM_VERIFY(data.center.size() == static_cast<std::size_t>(mesh.SpaceDimension()),
              "Current dipole center dimension must match mesh space dimension!");
  for (int d = 0; d < mesh.SpaceDimension(); d++)
  {
    center[d] = data.center[d];
  }

  // Create the VectorDeltaCoefficient
  dipole_coeff = std::make_unique<mfem::VectorDeltaCoefficient>(moment);
  dipole_coeff->SetDeltaCenter(center);
}

mfem::VectorDeltaCoefficient *CurrentDipoleData::GetDeltaCoefficient(double scale) const
{
  // Scale the dipole coefficient by the given factor
  dipole_coeff->SetScale(scale);
  return dipole_coeff.get();
}

CurrentDipoleOperator::CurrentDipoleOperator(const IoData &iodata,
                                             const mfem::ParMesh &mesh)
{
  // Set up current dipole source properties.
  SetUpDipoleProperties(iodata, mesh);
  PrintDipoleInfo(iodata, mesh);
}

void CurrentDipoleOperator::SetUpDipoleProperties(const IoData &iodata,
                                                  const mfem::ParMesh &mesh)
{
  // Set up current dipole data structures.
  for (const auto &[idx, data] : iodata.domains.current_dipole)
  {
    dipoles.emplace(idx, CurrentDipoleData(data, mesh));
  }
}

void CurrentDipoleOperator::PrintDipoleInfo(const IoData &iodata, const mfem::ParMesh &mesh)
{
  if (dipoles.empty())
  {
    return;
  }

  Mpi::Print("\nConfiguring electrical current dipole source terms:\n");
  for (const auto &[idx, data] : dipoles)
  {
    // Convert center coordinates back to physical units for display
    mfem::Vector physical_center = data.center;
    iodata.units.DimensionalizeInPlace<Units::ValueType::LENGTH>(physical_center);
    Mpi::Print(" Dipole {:d}: \n"
               " \tCenter = ({:.3e}, {:.3e}, {:.3e})\n"
               " \tMoment = ({:.3e}, {:.3e}, {:.3e})\n",
               idx, physical_center[0], physical_center[1],
               (mesh.SpaceDimension() > 2) ? physical_center[2] : 0.0, data.moment[0],
               data.moment[1], (mesh.SpaceDimension() > 2) ? data.moment[2] : 0.0);
  }
}

const CurrentDipoleData &CurrentDipoleOperator::GetDipole(int idx) const
{
  return dipoles.at(idx);
}

std::vector<mfem::Vector> CurrentDipoleOperator::GetMoments() const
{
  std::vector<mfem::Vector> moments;
  moments.reserve(dipoles.size());
  for (const auto &[idx, dipole] : dipoles)
  {
    moments.push_back(dipole.moment);
  }
  return moments;
}

std::vector<mfem::Vector> CurrentDipoleOperator::GetCenters() const
{
  std::vector<mfem::Vector> centers;
  centers.reserve(dipoles.size());
  for (const auto &[idx, dipole] : dipoles)
  {
    centers.push_back(dipole.center);
  }
  return centers;
}

void CurrentDipoleOperator::AddExcitationDomainIntegrators(int idx, mfem::LinearForm &rhs)
{
  // Add only the specific dipole with the given index
  auto it = dipoles.find(idx);
  if (it != dipoles.end())
  {
    const auto &dipole = it->second;

    // Create a VectorDeltaCoefficient for this specific dipole (RHS = J_source)
    auto dipole_coeff = std::make_unique<mfem::VectorDeltaCoefficient>(dipole.moment);
    dipole_coeff->SetDeltaCenter(dipole.center);
    dipole_coeff->SetScale(1.0);

    // Add as domain integrator
    rhs.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(*dipole_coeff));

    // Store the coefficient to prevent deletion.
    dipole_integrator_coeffs.push_back(std::move(dipole_coeff));
  }
}

void CurrentDipoleOperator::AddExcitationDomainIntegrators(mfem::LinearForm &rhs)
{
  // Add each dipole as a separate integrator
  for (const auto &[idx, dipole] : dipoles)
  {
    AddExcitationDomainIntegrators(idx, rhs);
  }
}

}  // namespace palace
