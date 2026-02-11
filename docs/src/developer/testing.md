```@raw html
<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
```

# Testing

*Palace* comes with two types of tests:

  - Unit tests in `test/unit/` test individual components in isolation
  - Regression tests in `test/examples/` compare code output against saved references

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
for more details on this.

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
is being run at the time. If you compiled *Palace*, the total number of
processes that could end up used in this example is `8 * OMP_NUM_THREADS`.

If you run specific tests categories, you can use regex matching, for example

```bash
ctest -R mpi-            # Run only MPI tests
ctest -R postoperator    # Run tests with postoperator in the name
```

### Running tests directly with the executable

CTest is powerful, but sometimes we want to have more direct control of how the
tests are being run. When this happens, we can use directly the `unit-tests`
executable.

To run all the `Serial` tests (including the benchmarks):
>>>>>>> 7675acc5 (Use ctest, parallelize tests, report coverage)

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

CTest skips all the benchmarks.

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

`Catch2` tags are typically used for filtering tests. *Palace* defines four
special tags that control when tests execute based on the runtime environment:

  - `[Serial]` tests run only with a single MPI process. Use this for tests that
    verify single-process behavior.
  - `[Parallel]` tests run only when multiple MPI processes are available. Use
    this for tests that verify MPI operations, distributed computations, or
    inter-process communication.
  - `[GPU]` tests run only when GPU devices are available. Use this for tests that
    are meaningful and interesting on GPU hardware.
  - `[NoConcurrent]` tests run serially without any other tests executing
    simultaneously. Use this for tests that modify shared resources or have
    timing-sensitive behavior that could be affected by concurrent execution.

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

!!! note "Try improving your test instead of relying on `[NoConcurrent]`"
    
    `[NoConcurrent]` serializes tests to ensure that shared resources are not
    being used/modified at the same time. Sometimes, this is needed (e.g.,
    you have tests that require large amounts of memory, or that talk through
    a socket to a GLVis server). Often times, howoever, it is possible to modify
    the test to remove potential race conditions and ensure isolation. A common
    case where you might be using shared resources is file-system operations. For
    that, you can create and destory temporary directories. `TempDirFixture` in
    `test-fixtures.hpp` helps you with that.

For the other tags, we recommend grouping related tests using descriptive tag
names like `[vector]` or `[materialoperator]` (typically named after files or
classes). This enables effective filtering, as described in [Building and
running unit tests](#Building-and-running-unit-tests).

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
    
    You need a version newer than 1.15. For LLVM, you will also need `llvm-profdata`
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
```

This automatically:

 1. Detects your compiler type (GCC or LLVM)
 2. Runs all tests via CTest (serial, MPI, GPU if available)
 3. Collects and merges coverage data
 4. Filters to *Palace* source code only
 5. Generates an HTML report at `../build/coverage_html/index.html`

You can control test parallelism with the `-j` flag:

```sh
./measure-test-coverage report -j8
```

Note that MPI tests run with 2 processes each, so `-j8` allows up to 4 MPI tests
simultaneously.

The `measure-test-coverage` allows for more control, for example:

```sh
# Generate coverage data only
./measure-test-coverage generate /path/to/unit-test

# Generate HTML report from existing coverage data
./measure-test-coverage report coverage_filtered.info

# Merge coverage from multiple builds
./measure-test-coverage generate /path/to/build1/unit-tests
./measure-test-coverage generate /path/to/build2/unit-tests
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
coverage report, merging the results from the serial, parallel, gpu, and
noconcurrent tests.

!!! note "Why two different coverage systems?"

    LLVM's source-based coverage provides more accurate instrumentation than the
    gcov-compatible mode, which attempts to emulate gcov behavior but has known
    reliability issues. The source-based approach offers better precision
    for complex codebases.

## Regression tests

In addition to unit tests, *Palace* comes with a series of regression tests.
Regression tests based on the provided example applications in the
[`examples/`](https://github.com/awslabs/palace/blob/main/examples/) directory
and verify that the code reproduces results in reference files stored in
[`test/examples/`](https://github.com/awslabs/palace/blob/main/test/examples/ref).

### Building and running example tests

#### Prerequisites

  - Julia
  - Palace executable in PATH or specified via environment variable/command-line argument

#### Setup

First, instantiate the Julia environment:

```bash
julia --project -e "using Pkg; Pkg.instantiate()"
```

You need to do this step only the very first time.

#### Command Line Arguments

The test runner supports command line arguments for configuration. Each argument can also be set via environment variables as fallbacks.

**Key Options:**

  - `--palace-test`: Path to *Palace* executable and optional arguments (default: "`palace`")
  - `--num-proc-test`: Number of MPI processes (default: number of physical cores)
  - `--test-cases`: Space-separated list of test cases to run (default: all examples)

Run `julia --project runtests.jl --help` to see all available options with descriptions and defaults.

#### Execution

Run all tests:

```bash
julia --project runtests.jl
```

Run specific test cases:

```bash
julia --project runtests.jl --test-cases "spheres rings"
```

Run with custom *Palace* executable:

```bash
julia --project runtests.jl --palace-test "../../build/bin/palace"
```

Run with custom number of processes:

```bash
julia --project runtests.jl --num-proc-test 4
```

You can also use environment variables as fallbacks:

```bash
TEST_CASES="spheres rings" julia --project runtests.jl
PALACE_TEST="../../build/bin/palace" julia --project runtests.jl
```

Each test case runs Palace simulations and compares generated CSV files against
reference data using configurable tolerances. When Palace behavior changes
legitimately (e.g., algorithm improvements), reference data can be updated using
the baseline script:

```bash
./baseline                    # Update all reference data
./baseline -e spheres         # Update specific example
./baseline --dry-run          # Test without updating files
./baseline -np 4              # Use 4 MPI processes
```
