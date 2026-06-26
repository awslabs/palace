```@raw html
<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
```

# Testing

*Palace* comes with two types of tests:

  - Unit tests in `test/unit/` test individual components in isolation
  - Regression tests, registered in `test/unit/regression/cases.cpp`, run
    full Palace solves on fixtures under `test/data/regression/input/`
    and compare generated CSVs against `test/data/regression/ref/`

Both types of tests are run automatically as part of the project's continuous
integration (CI) workflows.

## Building and running unit tests

The simplest way to build unit tests is using `CMake`.

First, follow the [quick start](../install.md#Quick-start) to build *Palace*
with `CMake`. From the *Palace* root directory:

```bash
mkdir build && cd build
cmake -DPALACE_MFEM_USE_EXCEPTIONS=yes ..
make -j $(nproc) palace-tests
```

The `PALACE_MFEM_USE_EXCEPTIONS` option is necessary to properly capture and
test assertions. If you want to also measure test coverage, turn
`PALACE_BUILD_WITH_COVERAGE` on. See [Unit test coverage](#Unit-test-coverage)
for more details on this. For Spack, these are the `mfem+exceptions` and the
`palace+coverage` variants.

Once the build completes, the `palace-unit-tests` executable will be installed
in the same `bin/` directory as the main `palace` executable, and you can run
tests in two ways:

 1. Using CTest (recommended for running all tests): Automatically runs all test
    categories in parallel with proper environment setup
 2. Using the test executable directly: Provides fine-grained control over which
    tests run using Catch2's filtering syntax

### Running tests with CTest

CTest provides automated test execution with parallel support and proper categorization:

```bash
cd palace-build
# If Palace was built with Spack: spack cd -b palace
ctest
```

This discovers and runs all registered tests. CTest automatically:

  - Runs serial tests with a single process
  - Runs MPI tests with 2 processes
  - Runs GPU tests with proper device configuration
  - Prevents concurrent execution of tests that require exclusive resources

To see all the tests available, call:

```bash
ctest -N
```

This will print something like:

```
Test project palace/build
  Test  #1: serial-Config Boundary Ports
  Test  #2: serial-Config Driven Solver
  Test  #3: serial-FarField
  Test  #4: serial-EM Constant Check
  ....
```

If you want to run all tests in parallel (e.g., with 8 processes):

```bash
ctest -j8 --output-on-failure
```

CTest will schedule execution of tests trying to use all the 8 processes. In
this, CTest handles MPI processes correctly, and ensures that only one GPU test
is being run at a time. If you compiled *Palace*, the total number of
processes that could end up used in this example is `8 * OMP_NUM_THREADS`.

To run specific tests categories, you can use regex matching, for example

```bash
ctest -R mpi-            # Run only MPI tests
ctest -R postoperator    # Run tests with postoperator in the name
```

### Running tests directly with the executable

CTest is powerful, but sometimes we want to have more direct control of how the
tests are being run. When this happens, we can use directly the `unit-tests`
executable.

To run all the `Serial` tests (including the benchmarks):

```bash
bin/palace-unit-tests
```

This runs all `Serial` tests. For `Parallel` tests:

```bash
mpirun -np 2 bin/palace-unit-tests
```

For `GPU` tests:

```bash
bin/palace-unit-tests --device cuda
```

You can use Catch2 filters to run specific subsets:

```bash
bin/palace-unit-tests [vector]              # Tests tagged [vector]
bin/palace-unit-tests "Vector Sum - Real"   # Specific test by name
bin/palace-unit-tests [^postoperator]       # Tests not tagged with [postoperator]
```

#### Benchmarks

The unit test application also includes a small number of benchmarks to compare
performance between MFEM's legacy assembly backend, MFEM's partial assembly
backend, and the specified libCEED backend (specified with the `--backend`
option, use `-h`/`--help` to list all command line options for the `palace-unit-tests`
executable). These can be run using, for example:

```bash
bin/palace-unit-tests "[Benchmark]" --benchmark-samples 10
```

These benchmarks can be accelerated using MPI and/or OpenMP parallelism (when
configured with `PALACE_WITH_OPENMP=ON`), but in all cases they are only testing
the local operator assembly on each process.

CTest is set up to skip all the benchmarks. More specifically, CTest uses Catch2's
`--skip-benchmarks` flag which skips the `BENCHMARK` sections within tests but still
runs the test logic itself (unlike the `~[Benchmark]` tag filter, which would exclude
those tests entirely).

### Adding unit tests

Palace uses [Catch2](https://github.com/catchorg/Catch2) for unit testing. Refer
to the `Catch2` documentation for detailed usage information. If you have never
added a test to *Palace*, start by reading our [tutorial on adding a new unit
test](tutorial_add_new_unit_test.md).

The most important `Catch2` feature  in the *Palace* test suite is
_tags_. When creating a test, you provide a name and a series of tags

```cpp
TEST_CASE("My name", "[MyFirstTag][MyOtherTag][Serial]"){ ... }
```

`Catch2` tags are typically used for filtering tests. *Palace* defines two
orthogonal axes of special tags:

**Execution-style tags** — control when a test runs based on the runtime
environment.

  - `[Serial]` tests run only with a single MPI process. Use this for tests that
    verify single-process behavior.
  - `[Parallel]` tests run only when multiple MPI processes are available. Use
    this for tests that verify MPI operations, distributed computations, or
    inter-process communication.
  - `[GPU]` tests run only when GPU devices are available. Use this for tests that
    are meaningful and interesting on GPU hardware.

**Category tags** — say what kind of test this is. Used by CTest registrations
to group cases into separate sweeps.

  - `[Regression]` tests are end-to-end Palace solves diffed against the
    reference data under `test/data/regression/ref/`. Slow enough to deserve
    their own ctest registration (`regression-*`, label `regression`)
    and skipped from the default unit-test sweep. Opt in with
    `ctest -L "^regression$"` or
    `palace-unit-tests "[Regression]~[Long]"`. Each case lives in
    `test/unit/regression/cases.cpp`.
  - `[Long]` is a *modifier* on a regression case for solves that take
    long enough (~10 minutes) that they shouldn't run on every PR.
    These cases are tagged `[Regression][Long]` and registered as
    `long-*` ctest entries (label `long`); the long-tests CI workflow
    runs them via `ctest -L "^long$"` when the `trigger-long-tests`
    PR label is applied.

Category tags are orthogonal to execution-style tags: a regression case is
still `[Serial][Parallel][GPU][Regression]` because the same case is valid
at any rank count and on either CPU or GPU.

These tags are inclusive, meaning that a test can be marked with multiple
special tags, if the test is meaningful in different contexts (e.g., if a test
supports CPU and GPU implementation at the same time). For example,
`[Serial][Parallel][GPU]` indicates a test that should run in every case.

!!! warning

    Tests without any of our special tags will never execute. Make sure to add at
    least one.

!!! note "Tests should use tags for which the test is meaningful"

    Consider testing this code:

    ```cpp
    int a = 1, b = 2;
    Assert(a + b == 3);
    ```

    This test will pass if run with MPI, but this is uninteresting since
    every MPI process executes the same code that the `Serial` test would
    run (same with GPUs). Therefore, this test should only be marked as
    `[Serial]`.

For the other tags, we recommend grouping related tests using descriptive tag
names like `[vector]` or `[materialoperator]` (typically named after files or
classes). This enables effective filtering, as described in [Building and
running unit tests](#Building-and-running-unit-tests).

#### Accessing files

Files required for tests (e.g., meshes or configurations) need to be saved
inside the `test/data` folder. This ensures that the files are accessible
when *Palace* is installed in a folder that is not the source folder (e.g., with
Spack). The path to the content of `test/data` maps to
`PALACE_TEST_DATA_DIR`. For example, if you want to access the `banana.txt` file
in `test/data`, refer to it as

```cpp
auto path_to_banana = fs::path(PALACE_TEST_DATA_DIR) / "banana.txt"
```

Accessing files without using `PALACE_TEST_DATA_DIR` will likely result in
failing tests on Spack builds (unless you know what you are doing).

#### Test configuration

The `CMakeLists.txt` in
[`test/unit/`](https://github.com/awslabs/palace/blob/main/test/unit/) contains
important settings:

  - The list of files compiled as part of the tests
  - The list of files that need to be compiled with a GPU compiler (`TARGET_SOURCES_DEVICE`)
  - Path variables used in tests, such as `MFEM_DATA_PATH`
    `MFEM_DATA_PATH` contains 2D and 3D sample meshes for testing (these meshes come
    from the [MFEM repository](https://github.com/mfem/mfem/tree/master/data)).

### Unit test coverage

When *Palace* is built with `PALACE_BUILD_WITH_COVERAGE`, running the unit tests
generates coverage information that tracks which parts of the codebase are
exercised by the tests. The `scripts/measure-test-coverage` script automates
the entire coverage measurement process, handling compiler detection, test
execution, and report generation.

!!! warning "Before continuing"

    Make sure you have `lcov` installed. To check that it works, run

    ```sh
    lcov --version
    ```

    `lcov` 2.0 or newer is recommended (for `--ignore-errors inconsistent`
    support); older versions will still work but may report non-fatal
    inconsistency warnings. For LLVM, you will also need `llvm-profdata`
    and `llvm-cov`. For GCC, you will need `gcov-tool` (which comes with GCC).

    Also make sure that *Palace* was build with `PALACE_BUILD_WITH_COVERAGE`.
    If not, compile it with

    ```sh
    cmake -DPALACE_MFEM_USE_EXCEPTIONS=yes -DPALACE_BUILD_WITH_COVERAGE=yes ..
    ```

Assuming you correctly built *Palace* in the `build` folder, the simplest
approach runs tests and generates an HTML report in one command:

```sh
cd scripts
./measure-test-coverage report
# With Spack
./measure-test-coverage report $(spack location -b palace)
```

The script takes the build directory as argument (defaulting to `../build`).
This automatically:

 1. Detects your compiler type (GCC or LLVM)
 2. Runs all tests via CTest (serial, MPI, GPU if available)
 3. Collects and merges coverage data
 4. Filters to *Palace* source code only
 5. Generates an HTML report at `build/coverage_html/index.html`

You can control test parallelism with the `-j` flag:

```sh
./measure-test-coverage report -j8
```

Note that MPI tests run with 2 processes each, so `-j8` allows up to 4 MPI tests
simultaneously. As noted in [Running tests with CTest](#Running-tests-with-CTest),
if *Palace* was built with OpenMP, account for the additional threads when choosing
the parallelism level.

The `measure-test-coverage` allows for more control, for example:

```sh
# Generate coverage data only
./measure-test-coverage generate /path/to/build
# With spack
./measure-test-coverage generate $(spack location -b palace)

# Generate HTML report from existing coverage data
./measure-test-coverage report coverage_filtered.info

# Merge coverage from multiple builds
./measure-test-coverage generate /path/to/build1
./measure-test-coverage generate /path/to/build2
./measure-test-coverage merge build1/coverage_filtered.info build2/coverage_filtered.info
./measure-test-coverage report coverage_filtered.info
```

#### Understanding the coverage system

*Palace* supports two coverage implementations:

  - LLVM source-based coverage (Clang, AppleClang, IntelLLVM): Provides accurate
    instrumentation through compiler integration. Coverage data is written to
    `.profraw` files during test execution, then merged and converted to LCOV
    format.
  - GCC gcov coverage: Uses compile-time annotations (`.gcno`) and runtime data
    (`.gcda`). The `measure-test-coverage` script handles the complexity of
    merging parallel test runs and unmangling filenames.

Both produce LCOV-compatible output for standardized processing and HTML
visualization.

Note that `measure-test-coverage` runs all the tests and produces only one
coverage report, merging the results from the serial, parallel, and gpu tests.

!!! note "Why two different coverage systems?"

    LLVM's source-based coverage provides more accurate instrumentation than the
    gcov-compatible mode, which attempts to emulate gcov behavior but has known
    reliability issues. The source-based approach offers better precision
    for complex codebases.

## AddressSanitizer and UndefinedBehaviorSanitizer

*Palace* can be built with [AddressSanitizer
(ASan)](https://clang.llvm.org/docs/AddressSanitizer.html) and
[UndefinedBehaviorSanitizer
(UBSan)](https://clang.llvm.org/docs/UndefinedBehaviorSanitizer.html) to detect
memory errors and undefined behavior at runtime. Both sanitizers are enabled
together with a single CMake option:

```bash
mkdir build && cd build
cmake -DPALACE_BUILD_WITH_SANITIZERS=ON ..
make -j palace-tests
```

This adds `-fsanitize=address,undefined -fno-omit-frame-pointer` to the compile
and link flags for both the `palace` library and the unit test executable.

Run the unit tests as usual:

```bash
bin/palace-unit-tests --skip-benchmarks
mpirun -np 2 bin/palace-unit-tests --skip-benchmarks
```

The following environment variables are useful when running under sanitizers:

  - `ASAN_OPTIONS=detect_leaks=0`: Disables leak detection, which can produce
    false positives from MPI runtimes and third-party libraries.
  - `LSAN_OPTIONS=suppressions=test/unit/lsan_suppressions.txt`: Suppresses
    known third-party leaks from MPI, hwloc, libevent, and libCEED instead of
    disabling leak detection entirely.
  - `UBSAN_OPTIONS=print_stacktrace=1`: Prints a full stack trace when undefined
    behavior is detected.

!!! note "Performance overhead"

    AddressSanitizer introduces singificant slowdowns. Use it only for debugging purposes.

!!! warning "GCC on macOS (Apple Silicon)"

    GCC does not support AddressSanitizer on macOS with Apple Silicon.
    Use Apple Clang or Homebrew LLVM/Clang instead when building
    with sanitizers on macOS.

## Regression tests

In addition to unit tests, *Palace* comes with a series of regression tests.
Regression tests use self-contained fixtures under
[`test/data/regression/input/`](https://github.com/awslabs/palace/blob/main/test/data/regression/input)
and verify that the code reproduces reference CSVs under
[`test/data/regression/ref/`](https://github.com/awslabs/palace/blob/main/test/data/regression/ref).
When a fixture is also a user-facing example, the source-tree regression input
can be a symlink to `examples/` to avoid duplication. Installation dereferences
those symlinks, so the test suite still runs from self-contained installed test
data rather than reading from the source-tree `examples/` directory.

## Tests in CI

*Palace* runs three tiers of CI checks on every pull request:

 1. Static analysis (`style.yml`): Code formatting (clang-format, Julia
    formatter), JSON schema validation, and test tag checks.
 2. Unit + regression tests (`build-and-test-linux.yml`,
    `build-and-test-macos.yml`, `spack.yml`, `containers.yml`): Builds *Palace*
    using CMake, Spack, and Singularity with a matrix of compilers (GCC, Clang,
    Intel), MPI implementations (Open MPI, MPICH, Intel MPI), math libraries
    (OpenBLAS, ARMPL, MKL), and link/integer options across x86, ARM, macOS, and
    GPU runners. Each configuration runs the full unit test suite (serial +
    parallel) and all default regression tests.
 3. Long tests (`long-tests.yml`): Expensive test cases (e.g.,
    `transmon/transmon_amr`) that are too slow to run on every push.

Long tests are a required check for merging non-trivial pull requests and are
not run by default. To run the Long tests, add the `trigger-long-tests` GitHub
label. The typical workflow is:

 1. Develop and iterate using the default (short) tests that run automatically.
 2. When the PR is ready to merge, add the **`trigger-long-tests`** label.
 3. The long tests run and the `long-tests` status is updated.
 4. If you push new commits while the label is present, the long tests re-run
    automatically.

Two special cases bypass the long test requirement:

  - Trivial changes: If a PR only touches documentation, README files, or
    other non-code files, the long test status is automatically set to success.
  - `no-long-tests` label: Adding this label bypasses the long test requirement
    entirely.

### Building and running regression tests

Regression cases live in `test/unit/regression/cases.cpp` and are
exercised by the same `palace-unit-tests` binary as the unit tests. They
carry the `[Regression]` category tag (orthogonal to `[Serial]`,
`[Parallel]`, `[GPU]`), and they are excluded from the default
`palace-unit-tests` invocation and the default `ctest` sweep.

#### Direct invocation

For a one-shot run, invoke the binary explicitly:

```bash
mpirun -n $NUM_PROC_TEST ./palace-unit-tests "[Regression]~[Long]"
```

Reference CSVs are rank-independent; CI rank counts are chosen only to
keep runner occupancy reasonable.

#### CTest invocation

Each regression case is also registered as an individual `regression-*`
CTest entry with label `regression`. The whole suite (serial, MPI, and
regression cases) shares two build-time cache variables:
`PALACE_TESTS_NUMPROC` (default 2) sets the MPI rank count and
`PALACE_TESTS_OMP_THREADS` (default 1) sets the OpenMP threads per process.
Each case reserves `PALACE_TESTS_NUMPROC * PALACE_TESTS_OMP_THREADS` CTest
slots (`PROCESSORS`) and the harness pins `OMP_NUM_THREADS` to the configured
value, so a runtime `OMP_NUM_THREADS` cannot desync the slot accounting. CI
passes both values explicitly; `ctest -j N` then only bounds how many cases
run concurrently.

```bash
# CI-style: let CTest schedule the regression cases
ctest -L "^regression$" -j "$(nproc)" --output-on-failure

# Just the cpw eigen cases
ctest -L "^regression$" -R cpw_.*_eigen -j "$(nproc)" --output-on-failure

# Long regression cases (transmon eigenmodes)
ctest -L "^long$" -j "$(nproc)" --output-on-failure
```

Catch2 tag combinations work as you'd expect:

```bash
mpirun -n 4 ./palace-unit-tests "[Regression]"         # every regression
mpirun -n 2 ./palace-unit-tests "[Regression]~*cpw*"   # skip cpw cases
./palace-unit-tests rings                              # single case, 1 rank
```

When using the CTest registrations, prefer `ctest -j N` over a manual loop;
CTest uses each case's `PROCESSORS` property to avoid oversubscribing
rank/thread slots.

Install-only runners (e.g. Spack jobs that install from a buildcache) have no
build tree, so the build-tree CTestTestfiles are unavailable. A relocatable
CTestTestfile is installed at `<prefix>/share/palace/test/CTestTestfile.cmake`
for exactly this case; with `palace-unit-tests` on `PATH` (e.g. `spack load palace`), run `ctest --test-dir <prefix>/share/palace/test -L "^regression$"`.
The linear solver and device are read from `PALACE_TEST_LINEAR_SOLVER` and
`PALACE_TEST_DEVICE` so one installed build serves the whole solver matrix.

#### Overrides

Regression input fixtures and references are normal test data under
`PALACE_TEST_DATA_DIR/regression` and are read-only. The only directory override
is `--regression-run-dir`, which changes where live outputs are staged; by
default this is `std::filesystem::temp_directory_path() / "palace-regression"`.

Per-case solver knobs:

| CLI flag                 | Effect                                     |
|:------------------------ |:------------------------------------------ |
| `--palace-linear-solver` | Overrides `Solver.Linear.Type` per case    |
| `--palace-eigensolver`   | Overrides `Solver.Eigenmode.Type` per case |

Run on GPU with `--device cuda`, which also adds the `[GPU]` tag.

#### Re-baselining

When Palace behaviour changes legitimately (algorithm improvements,
schema changes), regenerate the reference CSVs by running the case and
refreshing only the CSVs in the reference tree:

```bash
case=cpw config=cpw_lumped_uniform.json subdir=lumped_uniform
input="test/data/regression/input/$case"
mpirun -n "$NUM_PROC_TEST" palace "$input/$config"
dst="test/data/regression/ref/$case/$subdir"
rm -rf "$dst"
rsync -am --include='*/' --include='*.csv' --exclude='*' \
  "$input/postpro/$subdir/" "$dst/"
```

Repeat (or loop) over `case`/`config`/`subdir` to re-baseline several cases
at once, as the old `baseline` script did.

#### Adding a new regression case

 1. Drop the config and mesh/input files under
    `test/data/regression/input/<name>/` and the reference postpro tree
    under `test/data/regression/ref/<name>/<subdir>`. If an input file is also
    a user-facing example, make the regression input a symlink to the file under
    `examples/`; installed test data dereferences symlinks and remains
    self-contained.
 2. Add a `TEST_CASE("<name>", "[Serial][Parallel][GPU][Regression]")`
    to `test/unit/regression/cases.cpp`. Tack on `[Long]` if the case
    is too slow for the always-on regression job. Set `rtol`, `atol`,
    `excluded_columns`, `abs_columns`, `skip_rowcount`, expected
    output-directory flags, solver policies, and any `custom_checks`
    callbacks.
 3. Build `unit-tests`, then run
    `mpirun -n $NUM_PROC_TEST ./palace-unit-tests "<name>"` to validate.
