// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_EIGEN_SOLVER_HPP
#define PALACE_DRIVERS_EIGEN_SOLVER_HPP

#include <memory>
#include <vector>
#include "drivers/basesolver.hpp"
#include "drivers/singularsolver.hpp"

namespace palace
{

class ErrorIndicator;
class Mesh;
class ComplexOperator;
class SpaceOperator;

//
// Driver class for eigenmode simulations.
//
class EigenSolver : public BaseSolver
{
private:
  mutable FullWaveSingularFeatures singular_features;

  std::pair<ErrorIndicator, long long int>
  SolveSingular(SpaceOperator &space_op, std::unique_ptr<ComplexOperator> K,
                std::unique_ptr<ComplexOperator> M) const;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

public:
  using BaseSolver::BaseSolver;

  void Preprocess(IoData &iodata, std::unique_ptr<mfem::Mesh> &smesh,
                  MPI_Comm comm) const override;
  bool RequiresSourceSerialMeshMetadata() const override;
  void ProcessPartitionedMesh(const mfem::ParMesh &parallel_mesh,
                              const mesh::PartitionMetadata &metadata) const override;
  mfem::Array<int> GetRefinementProtection(const mfem::ParMesh &mesh) const override;
  void ProcessRefinedMesh(const mfem::ParMesh &mesh) const override;
  bool RebalanceRefinedMesh() const override;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_EIGEN_SOLVER_HPP
