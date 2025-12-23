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

// Helper class for current dipole sources in a model.
class CurrentDipoleData
{
public:
  // Current dipole properties.
  mfem::Vector direction;  // Normalized direction vector (unitless)
  double moment;           // Dipole moment magnitude [A·m]
  mfem::Vector center;     // Dipole center position [m, m, m]

  // Internal MFEM coefficient for the dipole
  std::unique_ptr<mfem::VectorDeltaCoefficient> dipole_coeff;

public:
  CurrentDipoleData(const config::CurrentDipoleData &data, const mfem::ParMesh &mesh);
};

//
// A class handling current dipole sources and their excitation.
//
class CurrentDipoleOperator
{
private:
  // Mapping from dipole index to data structure containing dipole information.
  std::map<int, CurrentDipoleData> dipoles;

  // Storage for integrator coefficients to manage their lifetime
  mutable std::vector<std::unique_ptr<mfem::VectorDeltaCoefficient>>
      dipole_integrator_coeffs;

  void SetUpDipoleProperties(const IoData &iodata, const mfem::ParMesh &mesh);
  void PrintDipoleInfo(const IoData &iodata, const mfem::ParMesh &mesh);

public:
  CurrentDipoleOperator(const IoData &iodata, const mfem::ParMesh &mesh);

  // Access data structures for current dipoles.
  auto begin() const { return dipoles.begin(); }
  auto end() const { return dipoles.end(); }
  auto rbegin() const { return dipoles.rbegin(); }
  auto rend() const { return dipoles.rend(); }
  auto Size() const { return dipoles.size(); }
  bool Empty() const { return dipoles.empty(); }

  // Add integrators directly to a LinearForm
  void AddExcitationDomainIntegrators(mfem::LinearForm &rhs);
  void AddExcitationDomainIntegrators(int idx, mfem::LinearForm &rhs);
};

}  // namespace palace

#endif  // PALACE_MODELS_CURRENT_DIPOLE_OPERATOR_HPP
