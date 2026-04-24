// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_BASE_SOLVER_HPP
#define PALACE_DRIVERS_BASE_SOLVER_HPP

#include <memory>
#include <vector>
#include <fmt/os.h>
#include "fem/errorindicator.hpp"
#include "utils/filesystem.hpp"
#include "utils/geodata.hpp"
#include "utils/memoryreporting.hpp"

namespace palace
{

class FiniteElementSpaceHierarchy;
class IoData;
class Mesh;
class Timer;
class PortExcitations;

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

  // Problem-type-specific serial-stage mesh preprocessing hook. Called in the read-mesh
  // pipeline between mesh::Load and mesh::Partition. The default is a no-op; BoundaryMode
  // overrides this to extract a 2D submesh from a 3D parent before partitioning, so the
  // rest of the pipeline (and AMR) operates on the mesh the problem actually solves on.
  virtual void PreprocessMesh(mesh::SerialMesh &smesh, MPI_Comm comm) const {}

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
};

// Archive the current postprocessing output for an AMR iteration. Creates a subfolder
// "iterationXX" inside output_dir and moves all files and directories into it, leaving
// relative symlinks behind so that the output directory always has accessible results.
// The palace.json metadata file is copied (not moved) since it is read and updated by
// subsequent iterations.
void SaveIteration(MPI_Comm comm, const fs::path &output_dir, int step, int width);

}  // namespace palace

#endif  // PALACE_DRIVERS_BASE_SOLVER_HPP
