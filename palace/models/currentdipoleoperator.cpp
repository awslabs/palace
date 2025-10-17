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

double CurrentDipoleData::GetDipoleMoment() const
{
  return moment.Norml2();
}

double CurrentDipoleData::GetExcitationPower() const
{
  // The current dipole excitation is normalized such that the power is 1
  return 1.0;
}

mfem::VectorCoefficient* CurrentDipoleData::GetModeCoefficient(double scale) const
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
  for (const auto &data : iodata.domains.current_dipole.dipoles)
  {
    dipoles.emplace_back(data, mesh);
  }
}

void CurrentDipoleOperator::PrintDipoleInfo(const IoData &iodata,
                                            const mfem::ParMesh &mesh)
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
               (mesh.SpaceDimension() > 2) ? data.center[2] : 0.0,
               data.moment[0], data.moment[1],
               (mesh.SpaceDimension() > 2) ? data.moment[2] : 0.0);
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

void CurrentDipoleOperator::AddExcitationDomainCoefficients(SumVectorCoefficient &fd)
{
  // Construct the RHS source term for current dipole sources, which looks like
  // -J_dipole for a current dipole source. The chosen current dipole J_dipole corresponds
  // to a unit current excitation.
  for (const auto &data : dipoles)
  {
    // Create a copy of the VectorDeltaCoefficient for this dipole
    auto dipole_coeff = std::make_unique<mfem::VectorDeltaCoefficient>(data.moment);
    dipole_coeff->SetDeltaCenter(data.center);

    // Add to the sum with negative sign (RHS = -J_source)
    fd.AddCoefficient(std::move(dipole_coeff), -1.0);
  }
}

}  // namespace palace
