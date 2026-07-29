// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_ELECTROSTATIC_SOLVER_HPP
#define PALACE_DRIVERS_ELECTROSTATIC_SOLVER_HPP

#include <map>
#include <memory>
#include <vector>
#include "drivers/basesolver.hpp"
#include "drivers/singularsolver.hpp"
#include "fem/singularfeatures.hpp"
#include "linalg/vector.hpp"
#include "utils/configfile.hpp"

namespace mfem
{

class DenseMatrix;

template <typename T>
class Array;

}  // namespace mfem

namespace palace
{

class ErrorIndicator;
class LaplaceOperator;
class Mesh;
template <ProblemType>
class PostOperator;

//
// Driver class for electrostatic simulations.
//
class ElectrostaticSolver : public BaseSolver
{
private:
  // Singular feature identity is extracted from the immutable source mesh and then
  // reconstructed on the partition using exact source-serial entity IDs.
  mutable fem::singular::FeatureTopology serial_singular_features;
  mutable fem::singular::FeatureTopology local_singular_features;
  mutable fem::singular::TriangleFeatureTopology serial_triangle_singular_features;
  mutable fem::singular::TriangleFeatureTopology local_triangle_singular_features;
  mutable std::vector<fem::singular::GlobalVertexId> source_vertex_ids;
  mutable std::vector<fem::singular::GlobalVertexId> source_element_ids;
  mutable NonconformingVertexIdentity vertex_identity;

  void PostprocessTerminals(PostOperator<ProblemType::ELECTROSTATIC> &post_op,
                            const std::map<int, mfem::Array<int>> &terminal_sources,
                            const std::vector<Vector> &V) const;
  void PostprocessSingularTerminals(const LaplaceOperator &laplace_op,
                                    const std::map<int, mfem::Array<int>> &terminal_sources,
                                    const std::vector<Vector> &V) const;
  void WriteTerminalMatrices(const std::map<int, mfem::Array<int>> &terminal_sources,
                             const mfem::DenseMatrix &C) const;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

public:
  using BaseSolver::BaseSolver;

  void Preprocess(IoData &iodata, std::unique_ptr<mfem::Mesh> &smesh,
                  MPI_Comm comm) const override;
  bool RequiresSourceSerialMeshMetadata() const override;
  void ProcessPartitionedMesh(const mfem::ParMesh &parallel_mesh,
                              const mesh::PartitionMetadata &metadata) const override;
  mesh::PartitionMetadata GetSourceEntityMetadata() const override;
  mfem::Array<int>
  GetRefinementProtection(const mfem::ParMesh &mesh, bool *conforming = nullptr,
                          mfem::Array<int> *repair = nullptr) const override;
  void ProcessRefinedMesh(const mfem::ParMesh &mesh) const override;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_ELECTROSTATIC_SOLVER_HPP
