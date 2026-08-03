// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_DRIVEN_SOLVER_HPP
#define PALACE_DRIVERS_DRIVEN_SOLVER_HPP

#include <memory>
#include <vector>
#include "drivers/basesolver.hpp"
#include "drivers/singularsolver.hpp"
#include "utils/configfile.hpp"

namespace palace
{

class ErrorIndicator;
class Mesh;
class SpaceOperator;

//
// Driver class for driven terminal simulations.
//
class DrivenSolver : public BaseSolver
{
private:
  mutable FullWaveSingularFeatures singular_features;

  ErrorIndicator SweepUniformSingular(SpaceOperator &space_op) const;

  ErrorIndicator SweepUniform(SpaceOperator &space_op) const;

  ErrorIndicator SweepAdaptive(SpaceOperator &space_op) const;

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
  void ObserveRefinementAncestry(const mfem::ParMesh &mesh) const override;
  mfem::Array<int> GetEnrichedElements(const mfem::ParMesh &mesh) const override;
  void ReportTraceComponents(const mfem::ParMesh &mesh,
                             const mfem::Array<int> &primary_marks) const override;
  void ProcessRefinedMesh(const mfem::ParMesh &mesh) const override;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_DRIVEN_SOLVER_HPP
