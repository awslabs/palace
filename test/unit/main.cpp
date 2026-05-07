// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <iostream>
#include <string>
#include <mfem.hpp>
#include <catch2/catch_session.hpp>
#include "fem/libceed/ceed.hpp"
#include "linalg/hypre.hpp"
#include "linalg/slepc.hpp"
#include "regression_helpers.hpp"
#include "utils/communication.hpp"
#include "utils/device.hpp"
#include "utils/omp.hpp"

using namespace palace;

// Global test options configurable from command line.
int benchmark_ref_levels = 0;
int benchmark_order = 4;
bool benchmark_assemble_q_data = false;
bool benchmark_no_fa = false;
bool benchmark_no_mfem_pa = false;

int main(int argc, char *argv[])
{
  // Initialize MPI.
  Mpi::Init(argc, argv);

  // See https://github.com/catchorg/Catch2/blob/devel/docs/own-main.md.
  Catch::Session session;

  // Extra command line arguments, mostly used for the test-libceed test suite.
  std::string device_str("cpu");          // MFEM device
  std::string ceed_backend("/cpu/self");  // libCEED backend

  // Regression-suite overrides. Each mirrors a Catch2 CLI flag;
  // examples-dir / regression-ref-dir fall back to their compile-time
  // defaults (PALACE_EXAMPLES_DIR_DEFAULT /
  // PALACE_REGRESSION_REF_DIR_DEFAULT wired in via
  // test/unit/CMakeLists.txt) when not overridden. Empty = "no
  // override".
  std::string examples_dir;        // --examples-dir
  std::string regression_ref_dir;  // --regression-ref-dir
  std::string regression_run_dir;  // --regression-run-dir
  std::string palace_solver;       // --palace-solver
  std::string palace_eigensolver;  // --palace-eigensolver
  std::string palace_device;       // --palace-device

  // Build a new parser on top of Catch2's.
  using namespace Catch::Clara;
  auto cli = session.cli() |
             Opt(device_str, "device")["--device"]("MFEM device (default: \"cpu\")") |
             Opt(ceed_backend,
                 "backend")["--backend"]("libCEED backend (default: \"/cpu/self\")") |
             Opt(examples_dir, "path")["--examples-dir"](
                 "Override for the source-tree examples/ directory used by "
                 "[Regression] cases (default: compile-time "
                 "PALACE_EXAMPLES_DIR_DEFAULT)") |
             Opt(regression_ref_dir, "path")["--regression-ref-dir"](
                 "Override for the test/examples/ref/ directory used by [Regression] "
                 "cases (default: compile-time PALACE_REGRESSION_REF_DIR_DEFAULT)") |
             Opt(regression_run_dir, "path")["--regression-run-dir"](
                 "Staging root under which each [Regression] case gets its own "
                 "subdirectory (inputs symlinked from examples/, outputs written "
                 "there). Default: std::filesystem::temp_directory_path() / "
                 "\"palace-regression\"") |
             Opt(palace_solver, "type")["--palace-solver"](
                 "Override Solver.Linear.Type for [Regression] cases (e.g. "
                 "SuperLU, STRUMPACK)") |
             Opt(palace_eigensolver, "type")["--palace-eigensolver"](
                 "Override Solver.Eigenmode.Type for [Regression] cases (e.g. "
                 "SLEPc, ARPACK)") |
             Opt(palace_device, "name")["--palace-device"](
                 "Override Solver.Device for [Regression] cases (CPU or GPU)") |
             Opt(benchmark_ref_levels, "levels")["--benchmark-ref-levels"](
                 "Levels of uniform mesh refinement for benchmarks (default: 0)") |
             Opt(benchmark_order, "order")["--benchmark-order"](
                 "Element order for benchmarks (default: 4)") |
             Opt(benchmark_assemble_q_data)["--benchmark-assemble-q-data"](
                 "Assemble quadrature data for benchmark operators") |
             Opt(benchmark_no_fa)["--benchmark-skip-full-assembly"](
                 "Skip full assembly tests in benchmarks") |
             Opt(benchmark_no_mfem_pa)["--benchmark-skip-mfem-partial-assembly"](
                 "Skip MFEM partial assembly tests in benchmarks");

  // Now pass the new composite back to Catch2 so it uses that.
  session.cli(cli);

  // Let Catch2 (using Clara) parse the command line.
  int result = session.applyCommandLine(argc, argv);
  if (result != 0)
  {
    return result;
  }

  int num_threads = palace::utils::ConfigureOmp();

  if (num_threads > 0)
    device_str += ",omp";

  // Initialize the device and assign GPUs to MPI ranks.
  mfem::Device device(
      device_str.c_str(),
      palace::utils::GetDeviceId(Mpi::World(), palace::utils::GetDeviceCount()));

  // Initialize HYPRE with correct memory location based on device.
  // TODO: Create a palace::Device class that takes care of all of this.
  hypre::Initialize();

  // Initialize SLEPc/PETSc (needed for eigenvalue solver tests).
#if defined(PALACE_WITH_SLEPC)
  slepc::Initialize();
  if (PETSC_COMM_WORLD != Mpi::World())
  {
    Mpi::Print(Mpi::World(), "Error: Problem during MPI initialization!\n\n");
    return 1;
  }
#endif

  // The Palace test suite defines three key tags:

  // - [Serial], for tests that are meaningful when run on a single process
  // - [Parallel], for tests that are meaningful when run on multiple processes
  // - [GPU], for tests that are meaningful when run on GPUs
  //
  // The tags are additive, meaning that a test can be tagged with all of them
  // (this also means that these tags cannot be used to filter out tests, only
  // to filter in).
  //
  // Here, we automatically add the relevant tags depending on the device/number
  // of MPI processes we detect.

  auto cfg = session.configData();

  // Whether the user passed any explicit test selector on the command
  // line (name, wildcard, or tag spec). Bare invocation leaves this
  // empty and triggers our convenience defaults below; anything
  // explicit runs exactly what was asked for. This matters for
  // regression cases in particular: they're tagged `[Regression]`
  // only, so auto-adding `[Serial]` or `~[Regression]` would either
  // fold in unrelated unit tests (separate positive filters are
  // OR'd in Catch2) or silently exclude the very test the user
  // selected by name.
  const bool user_selected_tests = !cfg.testsOrTags.empty();

  // Check if device is GPU capable, if yes, add the [GPU] tag.
  if (device.Allows(mfem::Backend::CUDA_MASK | mfem::Backend::HIP_MASK))
  {
    cfg.testsOrTags.emplace_back("[GPU]");
    if (ceed_backend == "/cpu/self")
    {
      // TODO: We pick magma because this is what Palace main does. We might
      // want to double check if this is the best default backend to choose.
      // Note that magma is a non-deterministic backend.
      ceed_backend =
          device.Allows(mfem::Backend::CUDA_MASK) ? "/gpu/cuda/magma" : "/gpu/hip/magma";
    }
  }
  // Bare invocation: auto-select by MPI world size, and defensively
  // exclude [Regression] (which is slow and has its own explicit
  // entry point; long regression cases are tagged [Regression][Long]
  // so they're caught here too). With any explicit user selector,
  // run exactly what was asked for.
  if (!user_selected_tests)
  {
    if (Mpi::Size(Mpi::World()) > 1)
    {
      cfg.testsOrTags.emplace_back("[Parallel]");
    }
    else
    {
      cfg.testsOrTags.emplace_back("[Serial]");
    }
    cfg.testsOrTags.emplace_back("~[Regression]");
  }
  session.useConfigData(cfg);

  // Forward regression-harness overrides into the helpers. Empty strings
  // leave the chain intact (env var, then compile-time default).
  palace::test::SetExamplesDirOverride(examples_dir);
  palace::test::SetRegressionRefDirOverride(regression_ref_dir);
  palace::test::SetRegressionRunDirOverride(regression_run_dir);
  palace::test::SetSolverOverride(palace_solver);
  palace::test::SetEigenSolverOverride(palace_eigensolver);
  palace::test::SetDeviceOverride(palace_device);

  // Only print from the root process.
  // TODO: Print errors from other processes as well.
  if (Mpi::Rank(Mpi::World()) != 0)
  {
    std::cout.rdbuf(NULL);
  }

  // Run the tests.
  ceed::Initialize(ceed_backend.c_str(), PALACE_LIBCEED_JIT_SOURCE_DIR);

  // Only print device info if not listing tests (for JSON output
  // compatibility), needed for test discovery.
  if (!cfg.listTests)
  {
    std::ostringstream resource(std::stringstream::out);
#ifdef PALACE_WITH_COVERAGE
    resource << "Built with code coverage for " << PALACE_WITH_COVERAGE << "\n";
#endif
    device.Print(resource);
    resource << "libCEED backend: " << ceed::Print();
    Mpi::Print("{}\n", resource.str());
  }

  result = session.run();
  ceed::Finalize();

  // Finalize SLEPc/PETSc.
#if defined(PALACE_WITH_SLEPC)
  slepc::Finalize();
#endif

  return result;
}
