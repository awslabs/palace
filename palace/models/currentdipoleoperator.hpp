// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_CURRENT_DIPOLE_OPERATOR_HPP
#define PALACE_MODELS_CURRENT_DIPOLE_OPERATOR_HPP

#include <map>
#include <memory>
#include <vector>
#include <mfem.hpp>

namespace palace
{

class IoData;
class SumVectorCoefficient;

namespace config
{

struct CurrentDipoleData;

}  // namespace config

//
// Helper class for current dipole sources in a model.
//
class CurrentDipoleData
{
public:
  // Current dipole properties.
  mfem::Vector moment;    // Dipole moment vector [A⋅m]
  mfem::Vector center;    // Dipole center position [m]

  // Internal MFEM coefficient for the dipole
  std::unique_ptr<mfem::VectorDeltaCoefficient> dipole_coeff;

public:
  CurrentDipoleData(const config::CurrentDipoleData &data, const mfem::ParMesh &mesh);

  // Get the dipole moment magnitude
  double GetDipoleMoment() const;

  // Get excitation power (normalized to unit power)
  double GetExcitationPower() const;

  // Get the VectorDeltaCoefficient for this dipole
  mfem::VectorCoefficient* GetModeCoefficient(double scale = 1.0) const;
};

//
// A class handling current dipole sources and their excitation.
//
class CurrentDipoleOperator
{
private:
  // Vector of current dipole data structures.
  std::vector<CurrentDipoleData> dipoles;

  void SetUpDipoleProperties(const IoData &iodata, const mfem::ParMesh &mesh);
  void PrintDipoleInfo(const IoData &iodata, const mfem::ParMesh &mesh);

public:
  CurrentDipoleOperator(const IoData &iodata, const mfem::ParMesh &mesh);

  // Access data structures for current dipoles.
  const CurrentDipoleData &GetDipole(std::size_t idx) const;
  auto begin() const { return dipoles.begin(); }
  auto end() const { return dipoles.end(); }
  auto rbegin() const { return dipoles.rbegin(); }
  auto rend() const { return dipoles.rend(); }
  auto Size() const { return dipoles.size(); }
  bool Empty() const { return dipoles.empty(); }

  // Get moment and center arrays separately (like donor_attributes and receiver_attributes)
  std::vector<mfem::Vector> GetMoments() const;
  std::vector<mfem::Vector> GetCenters() const;

  // Add contributions to the right-hand side source term vector for current dipole
  // excitations in the domain.
  void AddExcitationDomainCoefficients(SumVectorCoefficient &fd);
};

}  // namespace palace

#endif  // PALACE_MODELS_CURRENT_DIPOLE_OPERATOR_HPP
