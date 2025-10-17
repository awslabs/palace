// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "currentdipoleoperator.hpp"

#include "fem/coefficient.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

namespace palace
{

CurrentDipoleData::CurrentDipoleData(const config::CurrentDipoleData &data,
                                     const mfem::ParMesh &mesh, const Units &units)
{
  // Set up normalized direction
  direction.SetSize(mesh.SpaceDimension());
  MFEM_VERIFY(data.direction.size() == static_cast<std::size_t>(mesh.SpaceDimension()),
              "Current dipole direction dimension must match mesh space dimension!");
  for (int d = 0; d < mesh.SpaceDimension(); d++)
  {
    direction[d] = data.direction[d];
  }

  double dir_norm = direction.Norml2();
  MFEM_VERIFY(dir_norm > 0.0, "Current dipole direction magnitude must be positive!");
  direction /= dir_norm;  // Normalize to unit vector

  // Nondimensionalize moment [A·m]
  double current_scale = units.Nondimensionalize<Units::ValueType::CURRENT>(1.0);
  double length_scale = units.Nondimensionalize<Units::ValueType::LENGTH>(1.0);
  moment = data.moment * current_scale * length_scale;

  // Set up dipole center position
  center.SetSize(mesh.SpaceDimension());
  MFEM_VERIFY(data.center.size() == static_cast<std::size_t>(mesh.SpaceDimension()),
              "Current dipole center dimension must match mesh space dimension!");
  for (int d = 0; d < mesh.SpaceDimension(); d++)
  {
    center[d] = data.center[d];
  }
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
    dipoles.emplace(idx, CurrentDipoleData(data, mesh, iodata.units));
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

    // Convert moment back to physical units for display
    double physical_moment = iodata.units.Dimensionalize<Units::ValueType::CURRENT>(
        iodata.units.Dimensionalize<Units::ValueType::LENGTH>(data.moment));

    Mpi::Print(" Dipole {:d}: \n"
               " \tMoment = {:.3e} A·m\n"
               " \tCenter = ({:.3e}, {:.3e}, {:.3e})\n"
               " \tDirection = ({:.3e}, {:.3e}, {:.3e})\n",
               idx, physical_moment, physical_center[0], physical_center[1],
               (mesh.SpaceDimension() > 2) ? physical_center[2] : 0.0, data.direction[0],
               data.direction[1], (mesh.SpaceDimension() > 2) ? data.direction[2] : 0.0);
  }
}

void CurrentDipoleOperator::AddExcitationDomainIntegrators(int idx, mfem::LinearForm &rhs)
{
  auto it = dipoles.find(idx);
  if (it != dipoles.end())
  {
    const auto &dipole = it->second;

    // Create a VectorDeltaCoefficient for this specific dipole (RHS = -iω J_e)
    // RHS will be scaled by iω later in SpaceOperator.
    auto dipole_coeff =
        std::make_unique<mfem::VectorDeltaCoefficient>(dipole.direction.Size());
    dipole_coeff->SetDirection(dipole.direction);
    dipole_coeff->SetDeltaCenter(dipole.center);
    dipole_coeff->SetScale(-dipole.moment);

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
