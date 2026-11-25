// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_ELECTROSTATIC_SOLVER_HPP
#define PALACE_DRIVERS_ELECTROSTATIC_SOLVER_HPP

#include <filesystem>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>
#include "drivers/basesolver.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "utils/configfile.hpp"

namespace mfem
{

template <typename T>
class Array;

}  // namespace mfem

namespace palace
{

// Where a resolved setting came from: its default, its key, its alias variable, or both
// (equal).
enum class ResponseArchiveSource
{
  DEFAULT,
  KEY,
  VARIABLE,
  KEY_AND_VARIABLE
};

// The response-archive run mode of coupon generation: Solver.Electrostatic.ResponseArchive
// with its alias variables (Mode "Archive" / "Reduce" / "Estimate" =
// PALACE_RESPONSE_ARCHIVE_ONLY / _REDUCE_ONLY / _ESTIMATE_ONLY = "1"; Directory =
// PALACE_RESPONSE_ARCHIVE_DIR; BlockSize = PALACE_RESPONSE_BLOCK_SIZE; SourceTiming =
// PALACE_RESPONSE_SOURCE_TIMING; RecycleInitialGuess =
// PALACE_RESPONSE_RECYCLE_INITIAL_GUESS; EstimateRequest =
// PALACE_RESPONSE_ESTIMATE_REQUEST). Per setting: the key when only the key is given, the
// variable when only the variable is set, and an error when both are given with different
// values. The variables are read collectively (the root's value on every rank:
// utils/collectiveenvironment.hpp).
struct ResponseArchiveSettings
{
  bool archive_only = false;
  bool reduce_only = false;
  bool estimate_only = false;
  // Absolute and lexically normal; relative paths are taken against the working directory.
  std::optional<std::filesystem::path> directory;
  int block_size = 3;
  bool source_timing = false;
  bool recycle_initial_guess = false;
  std::optional<std::string> estimate_request;
  // The given settings an archive estimate refuses (an alias variable counts when it is
  // set at all, as before the keys existed).
  std::vector<std::string> estimate_conflicts;
  ResponseArchiveSource mode_source = ResponseArchiveSource::DEFAULT;
  ResponseArchiveSource directory_source = ResponseArchiveSource::DEFAULT;
  ResponseArchiveSource block_size_source = ResponseArchiveSource::DEFAULT;
  ResponseArchiveSource source_timing_source = ResponseArchiveSource::DEFAULT;
  ResponseArchiveSource recycle_initial_guess_source = ResponseArchiveSource::DEFAULT;
  ResponseArchiveSource estimate_request_source = ResponseArchiveSource::DEFAULT;

  // Whether a mode or an archive Directory is active (an ordinary solve otherwise).
  bool Active() const
  {
    return archive_only || reduce_only || estimate_only || directory.has_value();
  }
};
ResponseArchiveSettings
ResolveResponseArchive(const config::ElectrostaticSolverData &electrostatic, MPI_Comm comm);

// One log line naming every resolved setting with its source, e.g. "Response-archive
// settings: Mode=Archive (variable), Directory=/a (key), ...".
std::string ResponseArchiveSummary(const ResponseArchiveSettings &settings);

// Called before creating any output; the experimental reader must never reuse output.
// Also refuses, on the key path only, a Solver.Electrostatic.ResponseArchive key that has
// no effect in the resolved mode and a key-given Archive / Reduce mode without a
// Directory.
void ValidateArchiveEstimateOptions(const IoData &iodata, MPI_Comm comm,
                                    bool check_output = true);

class ErrorIndicator;
class LaplaceOperator;
class Mesh;
class SurfacePostGeometry;
class SurfaceResponseGeometry;
class SurfaceResponseOperator;
template <ProblemType>
class PostOperator;
template <OperatorType OperType>
class BaseKspSolver;
using KspSolver = BaseKspSolver<Operator>;

// The self-consistent response-corrected solve: PCG on K + Pᵀ (Q_fab - Q_thin) P with the
// thin-metal solver's preconditioner, at the correction's own tolerance, from the raw field
// as the initial guess. The corrected field is ACCEPTED only when PCG converged AND met no
// search direction of non-positive curvature (Ap, p) <= 0: a direction of negative
// curvature proves the corrected operator indefinite, and CG then "converges" to a
// meaningless field (sc-closure diagnostics 2026-09-29: the 120-degree corner model at p5
// gave a corrected domain correction 300x the fixed-trace one with five such directions).
// The record carries the CG diagnostics; ritz_min / ritz_max are the extreme Ritz values of
// the preconditioned corrected operator from the CG coefficients (a negative ritz_min is
// the same proof); they are the eigenvalues of the CG Lanczos tridiagonal of the m
// iterations (typically 3-30) and are NaN (not evaluated) when m exceeds
// corrected_solve_ritz_max_iterations. The curvature record exists for CgSolver only: with
// another Krylov solver (GMRES, FGMRES) the record stays empty (iterations 0, count 0, Ritz
// values NaN), the curvature check is a no-op and acceptance reduces to convergence, as
// before this check.
constexpr int corrected_solve_ritz_max_iterations = 1000;
struct CorrectedSolveRecord
{
  bool converged = false;
  int iterations = 0;
  double relative_residual = 0.0;
  int negative_curvature_count = 0;
  int first_negative_curvature_iteration = -1;
  double ritz_min = 0.0;
  double ritz_max = 0.0;
  bool accepted = false;
};
CorrectedSolveRecord SolveCorrectedField(KspSolver &ksp, const Operator &K,
                                         const Operator &corrected_K, double solve_tol,
                                         bool initial_guess, const Vector &rhs, Vector &x);

//
// Driver class for electrostatic simulations.
//
class ElectrostaticSolver : public BaseSolver
{
private:
  mutable std::shared_ptr<const SurfacePostGeometry> surface_post_geometry;
  mutable std::shared_ptr<const SurfaceResponseGeometry> response_geometry;

  void PostprocessTerminals(PostOperator<ProblemType::ELECTROSTATIC> &post_op,
                            const std::map<int, mfem::Array<int>> &terminal_sources,
                            const std::vector<Vector> &V) const;
  // The response-corrected capacitance matrices (terminal-C-corrected.csv with the C_m and
  // C⁻¹ variants; the raw terminal files are untouched): the fixed-trace matrix
  // C_ft(i, j) = V_jᵀ M V_i + V_jᵀ Pᵀ W (Q_fab,dom - Q_thin,dom) P V_i on the raw fields
  // (NaN unless fixed_trace: the configuration evaluates the fixed-trace response) and
  // the self-consistent matrix, the same fixed-trace bilinear form on the corrected fields
  // (the capacitance of the corrected operator when the corrected solve couples the
  // domain in its fixed-trace form; NaN unless self_consistent: every source's corrected
  // solve was accepted; fail closed). W are the final patch weights of the run. Lane L4
  // (decision 1829): `source_marked` (Coverage.Policy "MarkNotApplied"; root) marks the
  // entries built from a source whose final coverage holds a marked cell (nan + the
  // appended partial columns; palace.json nulls + the Partial matrices); nullptr = "Warn".
  void PostprocessCorrectedTerminals(
      PostOperator<ProblemType::ELECTROSTATIC> &post_op,
      const std::map<int, mfem::Array<int>> &terminal_sources, const std::vector<Vector> &V,
      const std::vector<Vector> &V_corrected, const SurfaceResponseOperator &response,
      bool fixed_trace, bool self_consistent,
      const std::vector<bool> *source_marked = nullptr) const;
  void PostprocessResponseMatrix(PostOperator<ProblemType::ELECTROSTATIC> &post_op,
                                 const LaplaceOperator &laplace_op, const Operator &Grad,
                                 const std::vector<Vector> &V,
                                 const std::vector<Vector> &D) const;
  void PostprocessArchivedResponseMatrix(PostOperator<ProblemType::ELECTROSTATIC> &post_op,
                                         const LaplaceOperator &laplace_op,
                                         const Operator &Grad,
                                         const std::filesystem::path &archive,
                                         int block_size) const;

  ErrorIndicator EstimateArchivedFields(LaplaceOperator &laplace_op, const Operator &K,
                                        const std::filesystem::path &archive,
                                        const std::string &estimate_request) const;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

public:
  using BaseSolver::BaseSolver;

  // CUT-PLAN step 4b (lane v1step4b): the Laplace operator's FE spaces (no stiffness
  // matrix) and the response operator as Solve constructs them, then the census.
  PreflightCensusOutcome PreflightResolutionCensus(std::vector<std::unique_ptr<Mesh>> &mesh,
                                                   bool fail_closed) const override;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_ELECTROSTATIC_SOLVER_HPP
