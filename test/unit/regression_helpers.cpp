// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "regression_helpers.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <set>
#include <sstream>
#include <string>
#include <system_error>
#include <utility>

#include <mpi.h>
#include <fmt/format.h>
#include <nlohmann/json.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "driver.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/omp.hpp"
#include "utils/tablecsv.hpp"

namespace palace::test
{

namespace
{

// In-TU override storage. Set once at test startup from the Catch2 CLI
// flag callbacks (see test/unit/main.cpp), read per-case.
std::string g_examples_dir_override;
std::string g_regression_ref_dir_override;
std::string g_regression_run_dir_override;
// These two default to the old Julia ArgConfig defaults. An omitted
// Catch2 flag should still inject "Default" for cases whose Julia port
// used `linear_solver=solver` / `eigen_solver=eigensolver`.
std::string g_solver_override = "Default";
std::string g_eigensolver_override = "Default";
std::string g_device_override;

// Recursively collect CSV and non-CSV metadata filenames relative to
// `root`. Folders named in `skip_dirs` are recorded at their top level
// but not descended into (they hold volumetric field output we don't
// diff). Returns two sorted sets so per-case file lists are stable.
struct FileListing
{
  std::set<std::string> dirs;
  std::set<std::string> csv_files;
  std::set<std::string> meta_files;
};

FileListing ListRegressionFiles(const std::filesystem::path &root,
                                const std::set<std::string> &skip_dirs)
{
  FileListing out;
  if (!std::filesystem::is_directory(root))
  {
    return out;
  }
  for (auto it = std::filesystem::recursive_directory_iterator(root);
       it != std::filesystem::recursive_directory_iterator(); ++it)
  {
    if (it->is_directory())
    {
      out.dirs.insert(std::filesystem::relative(it->path(), root).generic_string());
      if (skip_dirs.count(it->path().filename().string()))
      {
        it.disable_recursion_pending();
      }
      continue;
    }
    if (!it->is_regular_file())
    {
      continue;
    }
    const auto rel = std::filesystem::relative(it->path(), root).generic_string();
    if (it->path().extension() == ".csv")
    {
      out.csv_files.insert(rel);
    }
    else
    {
      out.meta_files.insert(rel);
    }
  }
  return out;
}

// Does `header` partial-match any excluded-column substring? Mirrors
// Julia's `Not(Cols(contains(col)))` behaviour. The match is on the
// column's human-readable `header_text`; Table::has() keys off the
// opaque "col_N" name assigned at parse time, which isn't what we want.
bool HeaderExcluded(const std::string &header, const std::vector<std::string> &excl)
{
  for (const auto &e : excl)
  {
    if (!e.empty() && header.find(e) != std::string::npos)
    {
      return true;
    }
  }
  return false;
}

std::string FormatSet(const std::set<std::string> &values)
{
  if (values.empty())
  {
    return "  <none>";
  }
  std::stringstream ss;
  for (const auto &value : values)
  {
    ss << "  " << value << '\n';
  }
  return ss.str();
}

void CheckSetEqual(const char *label, const std::set<std::string> &actual,
                   const std::set<std::string> &expected)
{
  INFO(label << "\nactual:\n" << FormatSet(actual) << "expected:\n" << FormatSet(expected));
  CHECK(actual == expected);
}

std::set<std::string> ExpectedDirectories(int max_its, const RegressionOptions &opts)
{
  std::set<std::string> expected;
  if (opts.gridfunction_fields)
  {
    expected.insert("gridfunction");
  }
  for (int i = 1; i <= max_its; ++i)
  {
    const std::string iteration = "iteration" + std::to_string(i);
    expected.insert(iteration);
    if (opts.gridfunction_fields)
    {
      expected.insert(iteration + "/gridfunction");
    }
    if (opts.paraview_fields)
    {
      expected.insert(iteration + "/paraview");
    }
  }
  if (opts.paraview_fields)
  {
    expected.insert("paraview");
  }
  return expected;
}

std::set<std::string> ExpectedMetadataFiles(int max_its)
{
  std::set<std::string> expected{"palace.json"};
  for (int i = 1; i <= max_its; ++i)
  {
    expected.insert("iteration" + std::to_string(i) + "/palace.json");
  }
  return expected;
}

nlohmann::json LoadConfigJson(const std::filesystem::path &config_path)
{
  std::stringstream buffer = PreprocessFile(config_path.string().c_str());
  return nlohmann::json::parse(buffer);
}

int GetMaxRefinementIterations(const nlohmann::json &config)
{
  const auto model_it = config.find("Model");
  if (model_it == config.end() || !model_it->is_object())
  {
    return 0;
  }
  const auto refinement_it = model_it->find("Refinement");
  if (refinement_it == model_it->end() || !refinement_it->is_object())
  {
    return 0;
  }
  const auto max_its_it = refinement_it->find("MaxIts");
  if (max_its_it == refinement_it->end())
  {
    return 0;
  }
  return max_its_it->get<int>();
}

// Load `path` via palace::TableWithCSVFile (load_existing_file=true),
// which reads the whole file and parses it through Table's string
// constructor.
Table LoadTable(const std::filesystem::path &path)
{
  TableWithCSVFile wrapped(path.string(), /*load_existing_file=*/true);
  return std::move(wrapped.table);
}

// Shared structural CSV checks. Julia first asserted raw column count,
// optionally asserted row count, dropped excluded columns, then asserted
// names before either custom_tests or the default isapprox broadcast.
// We keep the raw table intact for custom callbacks, but run the same
// shape/header validation for every CSV first so custom checks cannot
// silently bypass malformed output.
bool ValidateCSVTables(Table &actual, Table &reference, const RegressionOptions &opts)
{
  CHECK(actual.n_cols() == reference.n_cols());
  if (!opts.skip_rowcount)
  {
    CHECK(actual.n_rows() == reference.n_rows());
  }

  const bool comparable_columns = actual.n_cols() == reference.n_cols();
  const std::size_t n_cols = std::min(actual.n_cols(), reference.n_cols());
  for (std::size_t c = 0; c < n_cols; ++c)
  {
    const Column &a_col = actual[c];
    const Column &r_col = reference[c];
    if (HeaderExcluded(r_col.header_text, opts.excluded_columns))
    {
      continue;
    }

    INFO("column " << c << ": actual='" << a_col.header_text << "' reference='"
                   << r_col.header_text << "'");
    CHECK(a_col.header_text == r_col.header_text);
  }
  return comparable_columns;
}

bool SkipNumericComparison(const RegressionOptions &opts)
{
  return std::isinf(opts.rtol) && std::isinf(opts.atol);
}

// Column-wise diff of two CSV files using Catch2's WithinRel + WithinAbs
// matchers. Julia's isapprox (|a-b| <= max(atol, rtol*max(|a|,|b|))) is
// the "either condition" form, which is exactly
//   WithinRel(ref, rtol) || WithinAbs(ref, atol).
void CompareCSVFiles(Table &actual, Table &reference, const RegressionOptions &opts)
{
  if (!ValidateCSVTables(actual, reference, opts) || SkipNumericComparison(opts))
  {
    return;
  }

  const std::size_t n_cols = std::min(actual.n_cols(), reference.n_cols());
  const std::size_t n_rows = std::min(actual.n_rows(), reference.n_rows());
  for (std::size_t c = 0; c < n_cols; ++c)
  {
    const Column &a_col = actual[c];
    const Column &r_col = reference[c];
    if (HeaderExcluded(r_col.header_text, opts.excluded_columns))
    {
      continue;
    }

    for (std::size_t r = 0; r < n_rows; ++r)
    {
      const double av = a_col.data[r];
      const double rv = r_col.data[r];
      INFO("row " << r + 1 << " column '" << r_col.header_text << "'");
      CHECK_THAT(av, Catch::Matchers::WithinRel(rv, opts.rtol) ||
                         Catch::Matchers::WithinAbs(rv, opts.atol));
    }
  }
}

// RAII: per-case staging directory. Creates a fresh run directory
// under `stage_root / example_dir.basename() / case_label` (where
// `case_label` is the postpro subdirectory name, or just
// `example_dir.basename()` again if empty), symlinks each entry of
// the example directory (other than `postpro` / `log`) into it, then
// chdirs into the stage so palace::Run writes postpro/ + log/ there
// and not into the source tree. cwd is restored on destruction.
//
// The label suffix is what makes the stage path unique per TEST_CASE
// rather than per example directory. Without it, sibling cases that
// share an example directory (e.g. cylinder/cavity_pec and
// cylinder/cavity_impedance) would race for the same
// `<stage_root>/cylinder` path under `ctest -j` and the second
// starter's `remove_all` would pull the first's CWD out from under
// it (PETSc reports this as "Error in getcwd()").
//
// Rank 0 does the filesystem setup; a barrier before the chdir makes
// sure every rank sees a consistent stage. Rank 0 is also responsible
// for wiping a prior attempt, so each case runs against a clean
// layout regardless of what a previous invocation left behind.
class ScopedExampleStage
{
public:
  ScopedExampleStage(const std::filesystem::path &example_dir,
                     const std::filesystem::path &stage_root, std::string_view label,
                     MPI_Comm comm)
    : prev_(std::filesystem::current_path()),
      stage_(label.empty() ? stage_root / example_dir.filename()
                           : stage_root / example_dir.filename() / std::string{label})
  {
    if (Mpi::Root(comm))
    {
      std::error_code ec;
      std::filesystem::remove_all(stage_, ec);
      std::filesystem::create_directories(stage_);
      for (const auto &entry : std::filesystem::directory_iterator(example_dir))
      {
        const auto name = entry.path().filename().string();
        // postpro/ and log/ are Palace outputs; don't stage stale copies.
        if (name == "postpro" || name == "log")
        {
          continue;
        }
        std::filesystem::create_symlink(entry.path(), stage_ / name);
      }
    }
    Mpi::Barrier(comm);
    std::filesystem::current_path(stage_);
  }
  ~ScopedExampleStage()
  {
    std::error_code ec;
    std::filesystem::current_path(prev_, ec);  // never throw from dtor
    // The stage is left on disk for post-mortem inspection; the next
    // invocation of this case wipes and recreates it.
  }
  const std::filesystem::path &path() const { return stage_; }
  ScopedExampleStage(const ScopedExampleStage &) = delete;
  ScopedExampleStage &operator=(const ScopedExampleStage &) = delete;

private:
  std::filesystem::path prev_;
  std::filesystem::path stage_;
};

std::string ResolveSolverOverride(SolverOverridePolicy policy,
                                  const std::string &global_override)
{
  return (policy == SolverOverridePolicy::ForceDefault) ? "Default" : global_override;
}

// Inject any override knobs and return a fully constructed IoData.
// Matches the pre-processing Julia's testcase.jl did on the temp config
// before launching Palace.
IoData LoadCaseIoData(nlohmann::json config, const RegressionOptions &opts)
{
  auto &solver = config["Solver"];  // creates if missing
  if (!g_device_override.empty())
  {
    solver["Device"] = g_device_override;
  }

  const std::string linear_solver =
      ResolveSolverOverride(opts.linear_solver_policy, g_solver_override);
  if (!linear_solver.empty())
  {
    if (!solver.contains("Linear"))
    {
      solver["Linear"] = nlohmann::json::object();
    }
    solver["Linear"]["Type"] = linear_solver;
  }

  const std::string eigen_solver =
      ResolveSolverOverride(opts.eigen_solver_policy, g_eigensolver_override);
  if (!eigen_solver.empty() && solver.contains("Eigenmode"))
  {
    solver["Eigenmode"]["Type"] = eigen_solver;
  }
  return IoData(std::move(config), /*print=*/false);
}

}  // namespace

void SetExamplesDirOverride(std::string value)
{
  g_examples_dir_override = std::move(value);
}
void SetRegressionRefDirOverride(std::string value)
{
  g_regression_ref_dir_override = std::move(value);
}
void SetRegressionRunDirOverride(std::string value)
{
  g_regression_run_dir_override = std::move(value);
}
void SetSolverOverride(std::string value)
{
  g_solver_override = value.empty() ? "Default" : std::move(value);
}
void SetEigenSolverOverride(std::string value)
{
  g_eigensolver_override = value.empty() ? "Default" : std::move(value);
}
void SetDeviceOverride(std::string value)
{
  g_device_override = std::move(value);
}

std::filesystem::path GetExamplesDir()
{
  if (!g_examples_dir_override.empty())
  {
    return g_examples_dir_override;
  }
#ifdef PALACE_EXAMPLES_DIR_DEFAULT
  return PALACE_EXAMPLES_DIR_DEFAULT;
#else
  return {};
#endif
}

std::filesystem::path GetRegressionRefDir()
{
  if (!g_regression_ref_dir_override.empty())
  {
    return g_regression_ref_dir_override;
  }
#ifdef PALACE_REGRESSION_REF_DIR_DEFAULT
  return PALACE_REGRESSION_REF_DIR_DEFAULT;
#else
  return {};
#endif
}

std::filesystem::path GetRegressionRunDir()
{
  if (!g_regression_run_dir_override.empty())
  {
    return g_regression_run_dir_override;
  }
  return std::filesystem::temp_directory_path() / "palace-regression";
}

void RunRegressionCase(std::string_view case_dir, std::string_view config_json,
                       std::string_view postpro_subdir, const RegressionOptions &opts)
{
  MPI_Comm comm = Mpi::World();
  const auto examples_dir = GetExamplesDir();
  const auto ref_dir = GetRegressionRefDir();
  const auto run_root = GetRegressionRunDir();

  const std::filesystem::path example_path = examples_dir / case_dir;
  const std::filesystem::path config_path = example_path / config_json;
  const std::filesystem::path ref_postpro =
      ref_dir / std::filesystem::path(case_dir) / postpro_subdir;

  INFO("Regression case: " << case_dir << " / " << postpro_subdir << "\nexample:  "
                           << example_path << "\nconfig:   " << config_path
                           << "\nrun root: " << run_root << "\nref:      " << ref_postpro);

  REQUIRE(std::filesystem::exists(config_path));
  REQUIRE(std::filesystem::is_directory(ref_postpro));

  // Stage the example inputs under run_root / <case_dir> and chdir into
  // it. Palace writes postpro/ and log/ there instead of into the source
  // tree. The stage is torn down and recreated per invocation by
  // ScopedExampleStage.
  ScopedExampleStage stage(example_path, run_root, postpro_subdir, comm);
  const std::filesystem::path postpro_path = stage.path() / "postpro" / postpro_subdir;

  nlohmann::json config = LoadConfigJson(config_path);
  const int max_refinement_iterations = GetMaxRefinementIterations(config);

  IoData iodata = LoadCaseIoData(std::move(config), opts);
  const int omp_threads = palace::utils::ConfigureOmp();
  palace::Run(iodata, comm, omp_threads, /*git_tag=*/nullptr);

  // Post-Run comparisons run on rank 0 only. Invariants that keep this
  // safe:
  //   * palace::Run is collective and its last action (BlockTimer
  //     reductions + metadata write) synchronises every rank, so
  //     rank 0 can start reading without a barrier.
  //   * Every assertion below is CHECK (non-throwing). A mismatch
  //     records against the TEST_CASE and keeps going; nothing
  //     unwinds out of the rank-0 block to leave other ranks
  //     stranded.
  //   * Non-root ranks exit RunRegressionCase and move to the next
  //     TEST_CASE (or to session shutdown) at their own pace.
  if (Mpi::Root(comm))
  {
    const std::set<std::string> skip_dirs{"paraview", "gridfunction"};
    const FileListing got = ListRegressionFiles(postpro_path, skip_dirs);
    const FileListing want = ListRegressionFiles(ref_postpro, skip_dirs);

    CheckSetEqual("regression output directories", got.dirs,
                  ExpectedDirectories(max_refinement_iterations, opts));
    CheckSetEqual("regression CSV files", got.csv_files, want.csv_files);
    CheckSetEqual("regression metadata files", got.meta_files,
                  ExpectedMetadataFiles(max_refinement_iterations));

    for (const auto &rel : want.csv_files)
    {
      const auto actual = postpro_path / rel;
      const auto reference = ref_postpro / rel;
      const std::string basename = std::filesystem::path(rel).filename().string();
      INFO("CSV: " << rel << "\nactual:    " << actual << "\nreference: " << reference);
      if (!std::filesystem::is_regular_file(actual) ||
          !std::filesystem::is_regular_file(reference))
      {
        CHECK(std::filesystem::is_regular_file(actual));
        CHECK(std::filesystem::is_regular_file(reference));
        continue;
      }

      Table a = LoadTable(actual);
      Table r = LoadTable(reference);

      auto custom_it = opts.custom_checks.find(basename);
      if (custom_it != opts.custom_checks.end())
      {
        if (ValidateCSVTables(a, r, opts) && !SkipNumericComparison(opts))
        {
          custom_it->second(a, r);
        }
        continue;
      }

      CompareCSVFiles(a, r, opts);
    }
  }
}

}  // namespace palace::test
