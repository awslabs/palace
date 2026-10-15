// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_BASE_SOLVER_HPP
#define PALACE_DRIVERS_BASE_SOLVER_HPP

#include <memory>
#include <string>
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
class SurfacePostOperator;
class SurfaceResponseOperator;

namespace config
{
struct RefinementData;
}  // namespace config

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

  // Performs adaptive mesh refinement using the solve-estimate-mark-refine paradigm.
  // Dispatches to the Solve method for the driver specific calculations.
  void SolveEstimateMarkRefine(std::vector<std::unique_ptr<Mesh>> &mesh) const;

  // These methods write different simulation metadata to a JSON file in post_dir.
  void SaveMetadata(const FiniteElementSpaceHierarchy &fespaces) const;
  template <typename SolverType>
  void SaveMetadata(const SolverType &ksp) const;
  void SaveLinearSolverMetadata(MPI_Comm comm, long long int solves,
                                long long int iterations) const;
  void SaveMetadata(const Timer &timer) const;
  void SaveMetadata(const memory_reporting::MemoryStats &peak_memory) const;
  void SaveMetadata(const PortExcitations &excitation_helper) const;
  void SaveMetadata(const SurfaceResponseOperator &response) const;
  // Root only: palace.json SurfaceResponse.Diagnostics.EdgeExclusion.<interface index> =
  // the EdgeExcludeSegments coverage record of every interface dielectric carrying the
  // list (decision 657 (b): ExcludeSegmentsLength, ExcludedPerimeterLength and their ratio
  // Coverage, mesh length units; the qualification's EXCLUSION-COVERAGE gate reads it).
  // Nothing is written when no interface lists EdgeExcludeSegments.
  void SaveMetadata(const SurfacePostOperator &surf_post_op) const;
  void SaveSurfaceResponseSolverMetadata(MPI_Comm comm, const std::string &name,
                                         long long int solves,
                                         long long int iterations) const;
  // Root only: palace.json SurfaceResponse.<key> = value (the record of a response-derived
  // postprocessing result next to the operator statistics).
  void SaveSurfaceResponseMetadata(const std::string &key,
                                   const nlohmann::json &value) const;
  // Root only: palace.json AdaptiveMeshRefinement.MaxSizeStop = record (the AMR step at
  // which the predicted size of the refined mesh exceeded Refinement.MaxSize and the
  // solve was skipped).
  void SaveMaxSizeStopMetadata(const nlohmann::json &record) const;
};

// Predicted number of global true dofs of a refined mesh: the last solved size scaled by
// the element-count growth of the refinement step and rounded up (for a fixed polynomial
// order the dof count scales with the element count; the exact count depends on the space
// of the driver and is only known once that space is built on the refined mesh).
long long int PredictRefinedSize(long long int ntdof, long long int initial_elem_count,
                                 long long int final_elem_count);

// Refinement.MaxSize rule applied after a refinement step and before its solve: true when
// MaxSize is set and the predicted size of the refined mesh exceeds it. The legacy rule
// (MaxSizePredicted = false) never stops here and leaves MaxSize to the check against the
// last solved size at the head of the AMR loop.
bool ExceedsPredictedMaxSize(const config::RefinementData &refinement, long long int ntdof,
                             long long int initial_elem_count,
                             long long int final_elem_count);

// Refinement.MinElementSize applied to the local elements of a mesh: eligible[i] is true
// when the volume-based size of element i (mfem::Mesh::GetElementSize(i, 0), the h whose
// minimum and maximum the mesh statistics report) is at least 2 x min_element_size, so one
// refinement of a marked element stays at or above the floor. Both lengths are
// nondimensional, as the mesh coordinates.
std::vector<bool> ElementsEligibleForMarking(mfem::Mesh &mesh, double min_element_size);

// Global number of elements whose volume-based size is below min_element_size (after a
// refinement step: the elements the conformity closure created below the floor).
long long int CountElementsBelowSize(mfem::ParMesh &mesh, double min_element_size);

// Dörfler marking over the eligible elements only: the threshold is computed on the
// eligible indicators and an element is marked when it is eligible and its indicator
// reaches the threshold. The ineligible elements are counted with their share of the total
// squared indicator (global sums). When no element is eligible (globally) nothing is marked
// and the threshold and marked_error are 0.
struct FloorMarking
{
  mfem::Array<int> marked_elements;       // Local indices.
  long long int ineligible_elements = 0;  // Global count.
  double ineligible_error = 0.0;          // Share of the global sum of e^2.
  double threshold = 0.0;                 // Dörfler threshold over the eligible set.
  double marked_error = 0.0;              // Marked share of the eligible sum of e^2.
};
FloorMarking MarkEligibleElements(MPI_Comm comm, const Vector &e,
                                  const std::vector<bool> &eligible, double fraction);

// Archive the current postprocessing output for an AMR iteration. Creates a subfolder
// "iterationXX" inside output_dir and moves all files and directories into it, leaving
// relative symlinks behind so that the output directory always has accessible results.
// The palace.json metadata file is copied (not moved) since it is read and updated by
// subsequent iterations.
void SaveIteration(MPI_Comm comm, const fs::path &output_dir, int step, int width);

}  // namespace palace

#endif  // PALACE_DRIVERS_BASE_SOLVER_HPP
