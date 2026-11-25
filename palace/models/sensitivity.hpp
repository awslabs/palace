// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SENSITIVITY_HPP
#define PALACE_MODELS_SENSITIVITY_HPP

#include <array>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <mfem.hpp>
#include <nlohmann/json_fwd.hpp>
#include "linalg/vector.hpp"

namespace palace
{

class IoData;
class Mesh;

//
// The velocities of the design parameters on the input mesh: the interpolant of a table of
// their values at the nodes of the mesh file. The table is read and matched to the nodes
// on the root process, which keeps the fields and a point search for the evaluations.
//
class InputVelocities
{
private:
  std::vector<std::string> names;

  // On the root: the input mesh, the velocity of each parameter on its nodal space, its
  // nodes by their cells in a grid of spacing tol, and a point search, set up when first
  // needed.
  std::unique_ptr<mfem::Mesh> input_mesh;
  std::vector<mfem::GridFunction> fields;
  std::vector<std::array<double, 3>> nodes;
  std::unordered_map<std::uint64_t, std::vector<int>> cells;
  double tol = 0.0;
#if defined(MFEM_USE_GSLIB)
  mutable std::unique_ptr<mfem::FindPointsGSLIB> finder;
#endif

  void IndexNodes();
  std::vector<int> FindNodes(const std::array<double, 3> &x) const;

public:
  InputVelocities(const IoData &iodata, MPI_Comm comm);
  ~InputVelocities();

  const auto &GetNames() const { return names; }

  // The velocities [p][3 k + c] of the parameters at the points x_k, in mesh length units,
  // on the root process.
  std::vector<std::vector<double>> Evaluate(const std::vector<double> &x) const;
};

//
// The velocities of the design parameters at the mesh nodes: the input velocities, which
// move the mesh nodes as morphing the input mesh does.
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
  DesignVelocities(const InputVelocities *input, const IoData &iodata, const Mesh &mesh);
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

// Write the sensitivities to sensitivity.json in the given directory, from the root.
// Collective.
void WriteSensitivities(const nlohmann::json &sensitivities,
                        const std::filesystem::path &post_dir, MPI_Comm comm);

}  // namespace palace

#endif  // PALACE_MODELS_SENSITIVITY_HPP
