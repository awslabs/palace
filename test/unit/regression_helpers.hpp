// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_TEST_UNIT_REGRESSION_HELPERS_HPP
#define PALACE_TEST_UNIT_REGRESSION_HELPERS_HPP

#include <filesystem>
#include <functional>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace palace
{
class Table;
}  // namespace palace

namespace palace::test
{

// Callback signature for per-file custom checks. Tables are passed by
// non-const reference because palace::Table's indexing operators are
// non-const; the data itself is never mutated through them.
using CustomCheck = std::function<void(Table &actual, Table &reference)>;

// How a regression case consumes the global CLI solver override.
enum class SolverOverridePolicy
{
  // Use the global value (`--palace-linear-solver` / `--palace-eigensolver`).
  // Empty global means "no override; honour the case JSON".
  UseGlobalOverride,
  // Ignore the global and inject "Default" for this case.
  ForceDefault
};

// Options controlling a single regression case.
struct RegressionOptions
{
  double rtol = 1e-6;
  double atol = 1e-18;
  // Partial-match substrings on header_text; columns whose header
  // contains any of these are dropped from both tables before diffing.
  std::vector<std::string> excluded_columns;
  // Partial-match substrings on header_text; columns whose header contains
  // any of these are compared by magnitude only (|a| vs |b|). Used for
  // gauge-dependent signed quantities whose sign is not reproducible
  // (e.g. near-zero κ_ext for high-Q eigenmodes).
  std::vector<std::string> abs_columns;
  // Allow row-count mismatch (eigen / adaptive cases).
  bool skip_rowcount = false;
  // Compare at most this many leading rows. If unset for an Eigenmode case,
  // RunRegressionCase derives it from Solver.Eigenmode.N so extra converged
  // modes beyond the requested count are ignored.
  std::optional<std::size_t> max_rows;
  // Expected volumetric-output directories.
  bool paraview_fields = true;
  bool gridfunction_fields = false;
  SolverOverridePolicy linear_solver_policy = SolverOverridePolicy::UseGlobalOverride;
  SolverOverridePolicy eigen_solver_policy = SolverOverridePolicy::UseGlobalOverride;
  // Custom per-file checks keyed by relative path under postpro/ (e.g.
  // "farfield-rE.csv", "iteration1/port-S.csv").
  std::unordered_map<std::string, CustomCheck> custom_checks;
};

// End-to-end case runner. Calls palace::Run on every rank; rank 0
// diffs the live postpro/<subdir>/ tree against the reference tree.
void RunRegressionCase(std::string_view case_dir, std::string_view config_json,
                       std::string_view postpro_subdir, const RegressionOptions &opts = {});

// Directory overrides consulted in this order:
//   1. Catch2 CLI flag (--examples-dir / --regression-ref-dir /
//      --regression-run-dir).
//   2. PALACE_EXAMPLES_DIR / PALACE_REGRESSION_REF_DIR environment variables.
//   3. Local checkout paths next to the current build directory
//      (../examples and ../test/examples/ref), when present.
//   4. Compile-time defaults for examples/ref; for the run dir the
//      fallback is std::filesystem::temp_directory_path() /
//      "palace-regression" so the source tree stays read-only.
void SetExamplesDirOverride(std::string value);
void SetRegressionRefDirOverride(std::string value);
void SetRegressionRunDirOverride(std::string value);
std::filesystem::path GetExamplesDir();
std::filesystem::path GetRegressionRefDir();
std::filesystem::path GetRegressionRunDir();

// Solver overrides. Empty = no override; the config JSON value is used.
void SetSolverOverride(std::string value);
void SetEigenSolverOverride(std::string value);

}  // namespace palace::test

#endif  // PALACE_TEST_UNIT_REGRESSION_HELPERS_HPP
