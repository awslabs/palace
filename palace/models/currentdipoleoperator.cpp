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
  : coef(mesh.SpaceDimension())
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

  // Assemble the coefficient.
  coef.SetDeltaCenter(center);
  coef.SetDirection(direction);
  coef.SetScale(-moment);
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
    auto physical_center =
        iodata.units.Dimensionalize<Units::ValueType::LENGTH>(data.center);

    // Convert moment back to physical units for display
    double physical_moment = iodata.units.Dimensionalize<Units::ValueType::CURRENT>(
        iodata.units.Dimensionalize<Units::ValueType::LENGTH>(data.moment));

    Mpi::Print(" Dipole {:d}: \n"
               " \tMoment = {:.3e} A·m\n"
               " \tCenter = ({:.3e}) m\n"
               " \tDirection = ({:.3e}) m\n",
               idx, physical_moment, fmt::join(physical_center, ", "),
               fmt::join(data.direction, ", "));
  }
}

void CurrentDipoleOperator::AddExcitationDomainIntegrator(int idx, mfem::LinearForm &rhs)
{
  auto it = dipoles.find(idx);
  MFEM_VERIFY(it != dipoles.end(), "Invalid dipole index!");
  rhs.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(it->second.coef));
}

void CurrentDipoleOperator::AddExcitationDomainIntegrators(mfem::LinearForm &rhs)
{
  // Add each dipole as a separate integrator
  for (auto &[idx, dipole] : dipoles)
  {
    // Create a VectorDeltaCoefficient for this specific dipole (RHS = -iω J_e),
    // where J_e = moment × δ(x-x₀). RHS will be scaled by iω later in SpaceOperator.
    rhs.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(dipole.coef));
  }
}

}  // namespace palace
