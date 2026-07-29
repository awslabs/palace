// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_BASE_SOLVER_HPP
#define PALACE_DRIVERS_BASE_SOLVER_HPP

#include <memory>
#include <string_view>
#include <vector>
#include <fmt/os.h>
#include <nlohmann/json_fwd.hpp>
#include "fem/errorindicator.hpp"
#include "utils/filesystem.hpp"
#include "utils/memoryreporting.hpp"

namespace mfem
{
class Mesh;
class ParMesh;
}  // namespace mfem

namespace palace
{

class FiniteElementSpaceHierarchy;
class IoData;
class Mesh;
class Timer;
class PortExcitations;

namespace mesh
{
struct PartitionMetadata;
}  // namespace mesh

//
// Base driver class for all simulation types.
//
class BaseSolver
{
protected:
  // Reference to configuration file data (not owned).
  // TODO(C++20): Update to reference wrapper of incomplete type.
  const IoData &iodata;

  // Parameters for writing postprocessing outputs.
  fs::path post_dir;
  bool root;

  // Performs a solve using the mesh sequence, then reports error indicators and the number
  // of global true dofs.
  virtual std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const = 0;

public:
  BaseSolver(const IoData &iodata, bool root, int size = 0, int num_thread = 0,
             const char *git_tag = nullptr);
  virtual ~BaseSolver() = default;

  // Applying driver specific transformations to config data or mesh, including
  // nondimensionalization.
  virtual void Preprocess(IoData &iodata, std::unique_ptr<mfem::Mesh> &smesh,
                          MPI_Comm comm) const;

  // Driver-specific access to exact source-serial entity IDs transported through mesh
  // partitioning. Most solvers do not require this metadata.
  virtual bool RequiresSourceSerialMeshMetadata() const { return false; }
  virtual void ProcessPartitionedMesh(const mfem::ParMesh &mesh,
                                      const mesh::PartitionMetadata &metadata) const;
  virtual mesh::PartitionMetadata GetSourceEntityMetadata() const;

  // Optional AMR hooks for discretizations with mesh-dependent auxiliary topology.
  // A nonempty protection marker excludes local elements from refinement. The post-refine
  // hook runs before optional rebalancing and before the next operator reconstruction.
  virtual mfem::Array<int> GetRefinementProtection(const mfem::ParMesh &mesh) const;
  virtual void ProcessRefinedMesh(const mfem::ParMesh &mesh) const;
  virtual bool RebalanceRefinedMesh() const { return true; }

  // Performs adaptive mesh refinement using the solve-estimate-mark-refine paradigm.
  // Dispatches to the Solve method for the driver specific calculations.
  void SolveEstimateMarkRefine(std::vector<std::unique_ptr<Mesh>> &mesh) const;

  // These methods write different simulation metadata to a JSON file in post_dir.
  void SaveMetadata(const FiniteElementSpaceHierarchy &fespaces) const;
  template <typename SolverType>
  void SaveMetadata(const SolverType &ksp) const;
  void SaveMetadata(const Timer &timer) const;
  void SaveMetadata(const memory_reporting::MemoryStats &peak_memory) const;
  void SaveMetadata(const PortExcitations &excitation_helper) const;
  void SaveMetadata(std::string_view section, const nlohmann::json &data) const;
};

// Archive the current postprocessing output for an AMR iteration. Creates a subfolder
// "iterationXX" inside output_dir and moves all files and directories into it, leaving
// relative symlinks behind so that the output directory always has accessible results.
// The palace.json metadata file is copied (not moved) since it is read and updated by
// subsequent iterations.
void SaveIteration(MPI_Comm comm, const fs::path &output_dir, int step, int width);

}  // namespace palace

#endif  // PALACE_DRIVERS_BASE_SOLVER_HPP
