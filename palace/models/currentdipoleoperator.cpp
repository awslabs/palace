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
    dipoles.emplace_back(data, mesh);
  }
}

void CurrentDipoleOperator::PrintDipoleInfo(const IoData &iodata, const mfem::ParMesh &mesh)
{
  if (dipoles.empty())
  {
    return;
  }

  Mpi::Print("\nConfiguring current dipole excitation source terms:\n");
  for (std::size_t i = 0; i < dipoles.size(); ++i)
  {
    const auto &data = dipoles[i];
    Mpi::Print(" Dipole {:d}: Center = ({:.3e}, {:.3e}, {:.3e}), "
               "Moment = ({:.3e}, {:.3e}, {:.3e}) A⋅m\n",
               static_cast<int>(i), data.center[0], data.center[1],
               (mesh.SpaceDimension() > 2) ? data.center[2] : 0.0, data.moment[0],
               data.moment[1], (mesh.SpaceDimension() > 2) ? data.moment[2] : 0.0);
  }
}

const CurrentDipoleData &CurrentDipoleOperator::GetDipole(std::size_t idx) const
{
  MFEM_VERIFY(idx < dipoles.size(), "Current dipole index out of range!");
  return dipoles[idx];
}

std::vector<mfem::Vector> CurrentDipoleOperator::GetMoments() const
{
  std::vector<mfem::Vector> moments;
  moments.reserve(dipoles.size());
  for (const auto &dipole : dipoles)
  {
    moments.push_back(dipole.moment);
  }
  return moments;
}

std::vector<mfem::Vector> CurrentDipoleOperator::GetCenters() const
{
  std::vector<mfem::Vector> centers;
  centers.reserve(dipoles.size());
  for (const auto &dipole : dipoles)
  {
    centers.push_back(dipole.center);
  }
  return centers;
}

void CurrentDipoleOperator::AddExcitationDomainIntegrators(mfem::LinearForm &rhs)
{
  // Add each dipole as a separate integrator
  for (const auto &data : dipoles)
  {
    // Create a VectorDeltaCoefficient for this dipole with negative sign (RHS = -J_source)
    auto dipole_coeff = std::make_unique<mfem::VectorDeltaCoefficient>(data.moment);
    dipole_coeff->SetDeltaCenter(data.center);
    dipole_coeff->SetScale(-1.0);

    // Add as domain integrator
    rhs.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(*dipole_coeff));

    // Store the coefficient to prevent deletion (LinearForm doesn't take ownership of coefficient)
    dipole_integrator_coeffs.push_back(std::move(dipole_coeff));
  }
}

void CurrentDipoleOperator::AddExcitationDomainIntegrators(int idx, mfem::LinearForm &rhs)
{
  // For current dipoles, the index doesn't matter for now. Just call the main function that
  // adds all dipoles
  AddExcitationDomainIntegrators(rhs);
}

}  // namespace palace
