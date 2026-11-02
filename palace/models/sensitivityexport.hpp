// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SENSITIVITY_EXPORT_HPP
#define PALACE_MODELS_SENSITIVITY_EXPORT_HPP

#include <array>
#include <filesystem>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include <nlohmann/json.hpp>
#include "linalg/vector.hpp"

namespace palace
{

class IoData;
class Mesh;

//
// Export of the derivatives of objectives: their gradients with respect to the positions
// of the mesh nodes in a region, and other derivatives as JSON entries. The format is
// specified in scripts/julia/PalaceSensitivity/README.md.
//
class SensitivityExport
{
private:
  const IoData &iodata;
  MPI_Comm comm;

  // Owned true DOFs of the mesh nodes in the export region: the true DOF of each component,
  // and the position in mesh length units.
  std::vector<std::array<int, 3>> tdofs;
  std::vector<double> points;

  // The boundary attributes of each exported node, if requested.
  std::vector<std::vector<int>> boundary;

  // Objectives, with their gradients in their unit per mesh length unit.
  std::vector<nlohmann::json> objectives;
  std::vector<std::vector<double>> gradients;

  bool shape;
  int space_dim;

public:
  SensitivityExport(const IoData &iodata, const Mesh &mesh);

  // Whether gradients with respect to the node positions are exported.
  bool HasShape() const { return shape; }

  // Add an objective, given its JSON entry ("Index", "Type", its configuration, "Value" and
  // "Unit") and, if HasShape(), its gradient g on the true DOFs of the mesh nodal space in
  // nondimensional units. The gradient is multiplied by scale, the dimensional unit of the
  // objective.
  void AddObjective(nlohmann::json entry, const Vector *g, double scale);

  // Write the export to dir. Collective.
  void Write(const std::filesystem::path &dir) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_SENSITIVITY_EXPORT_HPP
