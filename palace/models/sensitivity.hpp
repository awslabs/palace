// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SENSITIVITY_HPP
#define PALACE_MODELS_SENSITIVITY_HPP

#include <filesystem>
#include <functional>
#include <memory>
#include <string>
#include <vector>
#include <mfem.hpp>
#include <nlohmann/json_fwd.hpp>
#include "linalg/vector.hpp"

namespace palace
{

class IoData;
class Mesh;

//
// The velocities of the design parameters at the mesh nodes: the input mesh's interpolant
// of a table of their values at its nodes, which moves the mesh nodes as morphing the input
// mesh does.
//
class DesignVelocities
{
private:
  const mfem::ParMesh &mesh;
  std::vector<std::string> names;

  // The velocity of each parameter on the true DOFs of the mesh nodal space, in
  // nondimensional length units per unit of the parameter.
  std::vector<Vector> velocities;

  // A copy of the mesh, whose nodes move for the derivatives of its geometry.
  mutable std::unique_ptr<mfem::ParMesh> moved_mesh;

public:
  DesignVelocities(const IoData &iodata, const Mesh &mesh);
  ~DesignVelocities();

  const auto &GetNames() const { return names; }
  auto Size() const { return names.size(); }
  const auto &GetVelocity(std::size_t p) const { return velocities[p]; }

  // The derivatives with respect to each parameter of a quantity with gradient g with
  // respect to the node positions, on the true DOFs of the mesh nodal space.
  std::vector<double> Contract(const Vector &g) const;

  // The derivatives [p][k] with respect to parameter p of geometric quantities f_k of the
  // mesh, by central differences on the mesh moved along the velocity. f is collective.
  std::vector<std::vector<double>>
  Differentiate(const std::function<std::vector<double>(const mfem::ParMesh &)> &f) const;
};

// The products u_jᵀ op u_i of L-vectors, for a symmetric operator.
mfem::DenseMatrix GramMatrix(const mfem::Operator &op, const std::vector<Vector> &u,
                             MPI_Comm comm);

// A dense matrix as a JSON array of its rows, times a scale.
nlohmann::json ToJson(const mfem::DenseMatrix &A, double scale);

// The names of the derivatives with respect to a relative change of a property of material
// m (0-based) or of an element (R, L or C) of lumped port idx.
std::string MaterialDerivativeName(int m, const std::string &property);
std::string LumpedPortDerivativeName(int idx, const std::string &element);

// Write the sensitivities to sensitivity.json in the given directory, from the root
// process.
void WriteSensitivities(const nlohmann::json &sensitivities,
                        const std::filesystem::path &post_dir, MPI_Comm comm);

}  // namespace palace

#endif  // PALACE_MODELS_SENSITIVITY_HPP
