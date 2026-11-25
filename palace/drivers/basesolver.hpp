// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_BASE_SOLVER_HPP
#define PALACE_DRIVERS_BASE_SOLVER_HPP

#include <memory>
#include <optional>
#include <string>
#include <vector>
#include <fmt/os.h>
#include <nlohmann/json_fwd.hpp>
#include "fem/errorindicator.hpp"
#include "models/surfaceparticipationtable.hpp"
#include "utils/filesystem.hpp"
#include "utils/geodata.hpp"
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
class SurfaceResponseGeometry;
class SurfaceResponseOperator;

namespace config
{
struct RefinementData;
}  // namespace config

//
// Base driver class for all simulation types.
//
// CUT-PLAN step 4b (lane v1step4b): what PreflightResolutionCensus did: whether it
// evaluated the census (false: the problem type has none, or the correction falls back to
// the raw run) and the features of its operator's recorded geometry fallbacks (raw recorded
// for them).
struct PreflightCensusOutcome
{
  bool evaluated = false;
  std::size_t fallback_features = 0;
};

class BaseSolver
{
protected:
  // Reference to configuration file data (not owned).
  const IoData &iodata;

  // Parameters for writing postprocessing outputs.
  fs::path post_dir;
  bool root;

  // Erratum 1784 E-3: the solve of the AMR sequence in progress (1 = the initial mesh) and
  // the largest number of solves (Refinement.MaxIts + 1 with AMR, 1 without), set by
  // SolveEstimateMarkRefine before each Solve and written as Level / Levels in the
  // coverage block of palace.json; whether a coverage block was written; whether the AMR
  // loop has ended with one written (set after the loop, never before): the end-of-run
  // timer write SaveMetadata(const Timer &) then adds "Final": true to the root's block in
  // its own load of palace.json (review mn-3: no extra full load). A run killed before that
  // write keeps the last written block without Final.
  mutable int solve_level = 0, solve_levels = 0;
  mutable bool coverage_block_written = false, coverage_final = false;

  // Performs a solve using the mesh sequence, then reports error indicators and the number
  // of global true dofs.
  virtual std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const = 0;

  // Set false by Solve when a linear solve fails to converge, or when the London kinetic
  // penalty energy spikes (a stiff λ→0 penalty the solver cannot resolve on a refined
  // mesh). Lets SolveEstimateMarkRefine halt adaptation gracefully, keeping the last
  // converged iteration, instead of aborting or reporting a collapsed inductance.
  mutable bool solve_converged_ = true;

  // CUT-PLAN step 4 (lane v1table): the per-iteration table surface-participation.csv of a
  // corrected electrostatic, eigenmode or driven run (null otherwise); made before the
  // first solve, alive across the adaptive levels; the drivers add each level's rows to it
  // (root) and SolveEstimateMarkRefine finishes each level after its solve.
  mutable std::unique_ptr<SurfaceParticipationTable> surface_participation;
  void FinishSurfaceParticipationLevel(long long int ntdof, MPI_Comm comm) const;

  // The first-level mesh of SolveEstimateMarkRefine (the a priori sequence flattened under
  // AMR or edge-tube refinement, then the configured correction tubes refined); returns
  // whether the run is adaptive. CUT-PLAN step 4b (lane v1step4b): PreflightOnly's census
  // reads the same mesh.
  bool PrepareFirstLevelMesh(std::vector<std::unique_ptr<Mesh>> &mesh) const;

  // CUT-PLAN step 4b (lane v1step4b): the tail of PreflightResolutionCensus after the
  // response operator's construction on the first mesh (`unknowns` the first solve's true
  // dofs): the census Types and the PREFLIGHT line, surface-response-resolution-patches.csv
  // (PatchEnergy), the Resolution block of palace.json and the level-1 PREFLIGHT rows of
  // surface-participation.csv. Collective.
  PreflightCensusOutcome WritePreflightResolutionCensus(SurfaceResponseOperator &response,
                                                        long long int unknowns,
                                                        MPI_Comm comm) const;

  // CUT-PLAN step 4b (ruling R1, G3): `fail_closed` = the configuration sets
  // GeometryFailurePolicy FailClosed and the caller built the census operator under
  // RecordRaw (the same branches; only the deferred abort differs). True, with one line,
  // when that operator recorded a geometry fallback: a solve-form run stops at its deferred
  // FailClosed abort there, so no census is written. Collective.
  bool PreflightCensusStopsFailClosed(const SurfaceResponseOperator &response,
                                      bool fail_closed, MPI_Comm comm) const;

  // PreflightResolutionCensus of the eigenmode and driven solvers: the space operator (no
  // system matrix) and the Maxwell response operator as their PostOperator builds it.
  PreflightCensusOutcome PreflightMaxwellResolutionCensus(
      const std::vector<std::unique_ptr<Mesh>> &mesh,
      std::shared_ptr<const SurfaceResponseGeometry> *response_geometry,
      bool fail_closed) const;

  // Previous adaptation iteration's total London kinetic energy S(A_t − a_h), for detecting
  // an anomalous spike; negative until the first London solve sets it.
  mutable double prev_london_kinetic_ = -1.0;

  // Called after each modification (refinement, rebalancing) of the finest mesh during
  // adaptation, for example to transfer solution data to the modified mesh.
  virtual void UpdateAfterMeshChange(Mesh &mesh) const {}

  // Whether the adaptive mesh refinement goes on after the solve of adaptation iteration it
  // (0 for the initial solve), with ntdof unknowns and the error estimate err.
  bool ContinueAdaptation(int it, long long int ntdof, double err) const;

public:
  BaseSolver(const IoData &iodata, bool root, int size = 0, int num_thread = 0,
             const char *git_tag = nullptr);
  virtual ~BaseSolver() = default;

  // Applying driver specific transformations to config data or mesh, including
  // nondimensionalization.
  virtual void Preprocess(IoData &iodata, std::unique_ptr<mfem::Mesh> &smesh,
                          MPI_Comm comm) const;

  // CUT-PLAN step 4b (lane v1step4b): PreflightOnly's resolution census. On the mesh the
  // first solve would read, builds the FE space and the response operator as the solve form
  // does before its first solve (no system matrix is assembled and nothing is solved; a
  // wave port's boundary-mode matrices are built with the space operator) and writes what
  // that solve form computes there (WritePreflightResolutionCensus; `fail_closed` as in
  // PreflightCensusStopsFailClosed).
  // Electrostatic, eigenmode and driven; any other problem type prints one line.
  // Collective.
  virtual PreflightCensusOutcome
  PreflightResolutionCensus(std::vector<std::unique_ptr<Mesh>> &mesh,
                            bool fail_closed) const;

  // Performs adaptive mesh refinement using the solve-estimate-mark-refine paradigm.
  // Dispatches to the Solve method for the driver specific calculations.
  void SolveEstimateMarkRefine(std::vector<std::unique_ptr<Mesh>> &mesh) const;

  // These methods write different simulation metadata to a JSON file in post_dir.
  void SaveMetadata(const FiniteElementSpaceHierarchy &fespaces) const;
  // Writes the given (1-based) AMR iteration index to post_dir's palace.json as
  // "Problem"/"Iteration", and first to archive_dir's copy if given. The index matches the
  // "iterationXX" archive subdirectory (the initial solve is iteration 1). It is written
  // once the iteration's output is in its final location, archived or not.
  void SaveAdaptationIteration(int iteration, const fs::path &archive_dir = {}) const;
  template <typename SolverType>
  void SaveMetadata(const SolverType &ksp) const;
  void SaveLinearSolverMetadata(MPI_Comm comm, long long int solves,
                                long long int iterations) const;
  void SaveMetadata(const Timer &timer) const;
  void SaveMetadata(const memory_reporting::MemoryStats &peak_memory) const;
  void SaveMetadata(const PortExcitations &excitation_helper) const;
  // The response operator's FINAL metadata write (round 2a R1, decision 1475 (1)): its
  // diagnostics record is moved into palace.json (SurfaceResponseOperator::
  // ReleaseStatistics), so the operator is not read after it.
  void SaveMetadata(SurfaceResponseOperator &response) const;
  // The impedance-metal rule's candidate block (ImpedanceMetalRule "Record"; root):
  // palace.json SurfaceResponse.ImpedanceRuleCandidate = `block` (marked "CANDIDATE: NOT OF
  // RECORD"); called only in RECORD mode, after the operator's metadata.
  void SaveImpedanceRuleCandidateMetadata(const nlohmann::json &block) const;
  // Root only: palace.json SurfaceResponse.Diagnostics.EdgeExclusion.<interface index> =
  // the EdgeExcludeSegments coverage record of every interface dielectric carrying the
  // list (decision 657 (b): ExcludeSegmentsLength, ExcludedPerimeterLength and their ratio
  // Coverage, mesh length units; the qualification's EXCLUSION-COVERAGE gate reads it).
  // Nothing is written when no interface lists EdgeExcludeSegments.
  void SaveMetadata(const SurfacePostOperator &surf_post_op) const;
  // Root only, after a rank-consistency check (collective): palace.json
  // SurfaceResponse.CorrectedSolver.<name>.LinearSolves / LinearIterations, and
  // NonConvergedLinearSolves when `unconverged` is given (the corrected solves whose Krylov
  // solver did not reach its tolerance; the self-consistent eigenproblem, decision 2159).
  void SaveSurfaceResponseSolverMetadata(
      MPI_Comm comm, const std::string &name, long long int solves,
      long long int iterations,
      std::optional<long long int> unconverged = std::nullopt) const;
  // Root only: palace.json SurfaceResponse.<key> = value (the record of a response-derived
  // postprocessing result next to the operator statistics).
  void SaveSurfaceResponseMetadata(const std::string &key,
                                   const nlohmann::json &value) const;
  // Root only: palace.json AdaptiveMeshRefinement.MaxSizeStop = record (the AMR step at
  // which the predicted size of the refined mesh exceeded Refinement.MaxSize and the
  // solve was skipped).
  void SaveMaxSizeStopMetadata(const nlohmann::json &record) const;
  void SaveMetadata(const mesh::MeshEntityCounts &counts) const;
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
// squared indicator (global sums). When no element is eligible (globally), or every
// eligible indicator is zero, markable is false, nothing is marked and the threshold and
// marked_error are 0.
struct FloorMarking
{
  mfem::Array<int> marked_elements;       // Local indices.
  long long int ineligible_elements = 0;  // Global count.
  double ineligible_error = 0.0;          // Share of the global sum of e^2.
  bool markable = false;                  // Some eligible element has a nonzero indicator.
  double threshold = 0.0;                 // Dörfler threshold over the eligible set.
  double marked_error = 0.0;              // Marked share of the eligible sum of e^2.
};
FloorMarking MarkEligibleElements(MPI_Comm comm, const Vector &e,
                                  const std::vector<bool> &eligible, double fraction);

// Archive the current postprocessing output for an AMR iteration. Creates a subfolder
// "iterationXX" inside output_dir and moves all files and directories into it, leaving
// relative symlinks behind so that the output directory always has accessible results.
// The palace.json metadata file is copied (not moved) since it is read and updated by
// subsequent iterations. Returns the path of the "iterationXX" subfolder.
fs::path SaveIteration(MPI_Comm comm, const fs::path &output_dir, int step, int width);

}  // namespace palace

#endif  // PALACE_DRIVERS_BASE_SOLVER_HPP
