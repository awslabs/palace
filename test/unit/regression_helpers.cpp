// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "regression_helpers.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <initializer_list>
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
#include "utils/outputdir.hpp"
#include "utils/tablecsv.hpp"
#include "utils/timer.hpp"

namespace palace::test
{

namespace
{

// In-TU override storage. Set once at test startup from the Catch2 CLI
// flag callbacks (see test/unit/main.cpp), read per-case.
std::string g_examples_dir_override;
std::string g_regression_ref_dir_override;
std::string g_regression_run_dir_override;
std::string g_solver_override;
std::string g_eigensolver_override;

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

// Substring match on the human-readable `header_text` (not Table's
// opaque "col_N" key). Empty patterns are ignored.
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

// Shape and header validation: column count, optional row count, then
// header names. Run for every CSV before either the default comparator
// or a custom callback so malformed output can't slip past.
bool ValidateCSVTables(Table &actual, Table &reference, const RegressionOptions &opts)
{
  CHECK(actual.n_cols() == reference.n_cols());
  if (!opts.skip_rowcount)
  {
    CHECK(actual.n_rows() == reference.n_rows());
  }
  else
  {
    // Eigen/adaptive cases may have a different row count, but empty output
    // against a non-empty reference is still a failure.
    CHECK((actual.n_rows() > 0) == (reference.n_rows() > 0));
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

// Column-wise diff using Catch2's WithinRel + WithinAbs matchers.
//   |a - b| <= max(atol, rtol * max(|a|, |b|))
// is the "either condition" form: WithinRel(ref, rtol) || WithinAbs(ref, atol).
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

    const bool by_magnitude = HeaderExcluded(r_col.header_text, opts.abs_columns);
    for (std::size_t r = 0; r < n_rows; ++r)
    {
      double av = a_col.data[r];
      double rv = r_col.data[r];
      if (by_magnitude)
      {
        // Gauge-dependent signed column: compare magnitude only.
        av = std::abs(av);
        rv = std::abs(rv);
      }
      INFO("row " << r + 1 << " column '" << r_col.header_text << "'");
      CHECK_THAT(av, Catch::Matchers::WithinRel(rv, opts.rtol) ||
                         Catch::Matchers::WithinAbs(rv, opts.atol));
    }
  }
}

// RAII: per-case staging directory. Symlinks the example inputs into
// `stage_root / <example>/ <label>` and chdirs in so palace::Run writes
// postpro/ and log/ there. The label suffix keeps sibling cases that
// share an example directory (e.g. cylinder/cavity_pec and
// cylinder/cavity_impedance) on disjoint paths under ctest -j.
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
    int ok = 1;
    if (Mpi::Root(comm))
    {
      std::error_code ec;
      std::filesystem::remove_all(stage_, ec);
      std::filesystem::create_directories(stage_, ec);
      ok = ec ? 0 : 1;
      for (const auto &entry : std::filesystem::directory_iterator(example_dir, ec))
      {
        const auto name = entry.path().filename().string();
        // postpro/ and log/ are Palace outputs; don't stage stale copies.
        if (name == "postpro" || name == "log")
        {
          continue;
        }
        std::filesystem::create_symlink(entry.path(), stage_ / name, ec);
        if (ec)
        {
          ok = 0;
          break;
        }
      }
    }
    // Broadcast the staging result so all ranks fail together instead of
    // blocking at the barrier if root errored.
    Mpi::Broadcast(1, &ok, 0, comm);
    REQUIRE(ok);
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
IoData LoadCaseIoData(nlohmann::json config, const RegressionOptions &opts)
{
  auto &solver = config["Solver"];  // creates if missing

  if (const std::string linear_solver =
          ResolveSolverOverride(opts.linear_solver_policy, g_solver_override);
      !linear_solver.empty())
  {
    if (!solver.contains("Linear"))
    {
      solver["Linear"] = nlohmann::json::object();
    }
    solver["Linear"]["Type"] = linear_solver;
  }

  if (const std::string eigen_solver =
          ResolveSolverOverride(opts.eigen_solver_policy, g_eigensolver_override);
      !eigen_solver.empty() && solver.contains("Eigenmode"))
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
  g_solver_override = std::move(value);
}
void SetEigenSolverOverride(std::string value)
{
  g_eigensolver_override = std::move(value);
}

namespace
{

std::filesystem::path LocalCheckoutPath(std::initializer_list<std::string_view> parts)
{
  auto path = std::filesystem::current_path();
  for (auto part : parts)
  {
    path /= std::string(part);
  }

  std::error_code ec;
  if (std::filesystem::is_directory(path, ec))
  {
    return path.lexically_normal();
  }
  return {};
}

}  // namespace

std::filesystem::path GetExamplesDir()
{
  if (!g_examples_dir_override.empty())
  {
    return g_examples_dir_override;
  }
  if (const char *env = std::getenv("PALACE_EXAMPLES_DIR"))
  {
    return env;
  }
  if (auto path = LocalCheckoutPath({"..", "examples"}); !path.empty())
  {
    return path;
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
  if (const char *env = std::getenv("PALACE_REGRESSION_REF_DIR"))
  {
    return env;
  }
  if (auto path = LocalCheckoutPath({"..", "test", "examples", "ref"}); !path.empty())
  {
    return path;
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

namespace
{

// Fail the case (on every rank) if a required input is missing, telling the
// user how to point the harness at a checkout that vends the data.
void RequireArtifact(const std::filesystem::path &path, std::string_view what)
{
  INFO("Missing regression "
       << what << ": " << path
       << "\nThe regression suite needs the example config, its mesh, and the "
          "reference data. Vend them from a Palace checkout via the "
          "--examples-dir / --regression-ref-dir flags or the PALACE_EXAMPLES_DIR "
          "/ PALACE_REGRESSION_REF_DIR environment variables.");
  REQUIRE(std::filesystem::exists(path));
}

}  // namespace

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

  // Fail clearly if the user has not vended the config, mesh, or reference data.
  RequireArtifact(config_path, "config file");
  RequireArtifact(ref_postpro, "reference data directory");
  nlohmann::json config = LoadConfigJson(config_path);
  if (auto model = config.find("Model"); model != config.end())
  {
    if (auto mesh = model->find("Mesh"); mesh != model->end() && mesh->is_string())
    {
      std::filesystem::path mesh_path = mesh->get<std::string>();
      if (mesh_path.is_relative())
      {
        mesh_path = example_path / mesh_path;
      }
      RequireArtifact(mesh_path, "mesh file");
    }
  }

  ScopedExampleStage stage(example_path, run_root, postpro_subdir, comm);
  const std::filesystem::path postpro_path = stage.path() / "postpro" / postpro_subdir;

  const int max_refinement_iterations = GetMaxRefinementIterations(config);

  IoData iodata = LoadCaseIoData(std::move(config), opts);
  MakeOutputFolder(iodata, comm);
  const int omp_threads = palace::utils::ConfigureOmp();
  // Wipe BlockTimer state so timings/peak-memory don't accumulate across cases.
  BlockTimer::Reset();
  palace::Run(iodata, comm, omp_threads, /*git_tag=*/nullptr);

  // palace::Run ends with a collective metadata write, so rank 0 can read
  // without a barrier. Use CHECK (non-throwing) so a mismatch on rank 0
  // doesn't strand the other ranks.
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

      auto custom_it = opts.custom_checks.find(rel);
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
