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
class Mesh;
class SumVectorCoefficient;

namespace config
{

struct CurrentDipoleData;

}  // namespace config

// Helper class for current dipole sources in a model.
class CurrentDipoleData
{
public:
  // Current dipole properties (stored in nondimensional form).
  mfem::Vector direction;  // Normalized direction vector
  double moment;           // Dipole moment magnitude
  mfem::Vector center;     // Dipole center position

  // Coefficient stored with data given unique status of delta coefficient.
  mfem::VectorDeltaCoefficient coef;  // Vector coefficient for dipole source
public:
  CurrentDipoleData(const config::CurrentDipoleData &data, const Mesh &mesh,
                    const class Units &units);
};

//
// A class handling current dipole sources and their excitation.
//
class CurrentDipoleOperator
{
private:
  // Mapping from dipole index to data structure containing dipole information.
  std::map<int, CurrentDipoleData> dipoles;

  void SetUpDipoleProperties(const IoData &iodata, const Mesh &mesh);
  void PrintDipoleInfo(const IoData &iodata, const Mesh &mesh);

public:
  CurrentDipoleOperator(const IoData &iodata, const Mesh &mesh);

  // Access data structures for current dipoles.
  auto begin() const { return dipoles.begin(); }
  auto end() const { return dipoles.end(); }
  auto rbegin() const { return dipoles.rbegin(); }
  auto rend() const { return dipoles.rend(); }
  auto Size() const { return dipoles.size(); }
  bool Empty() const { return dipoles.empty(); }

  // Add integrators directly to a LinearForm
  void AddExcitationDomainIntegrator(int idx, mfem::LinearForm &rhs);
  void AddExcitationDomainIntegrators(mfem::LinearForm &rhs);
};

}  // namespace palace

#endif  // PALACE_MODELS_CURRENT_DIPOLE_OPERATOR_HPP
