```@raw html
<!--- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. --->
<!--- SPDX-License-Identifier: Apache-2.0 --->
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
make -j palace-tests
```

The `PALACE_MFEM_USE_EXCEPTIONS` option is necessary to properly capture and
test assertions. If you want to also measure test coverage, turn
`PALACE_BUILD_WITH_COVERAGE` on. See [Unit test coverage](#Unit-test-coverage)
for more details on this.

Once the build completes, the `palace-unit-tests` executable will be installed in the
same `bin/` directory as the main `palace` executable. Run the tests with:

```bash
bin/palace-unit-tests
```

This runs all the `Serial` tests. To run `Parallel` tests, call

```bash
mpirun -np 2 bin/palace-unit-tests
```

To run the `GPU` tests (this can be combined with `mpirun`):

```bash
bin/palace-unit-tests --device cuda
```

If you are interested in running only a subset of tests, you can use filters.
For example, to run the tests tagged with `[vector]`, call

```bash
bin/palace-unit-tests [vector]
```

To run a specific test by name:

```bash
bin/palace-unit-tests "Vector Sum - Real"
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

`Catch2` tags are typically used for filtering tests. *Palace* defines three
special tags that control when tests execute based on the runtime environment:

  - `[Serial]` tests run only with a single MPI process. Use this for tests that
    verify single-process behavior.
  - `[Parallel]` tests run only when multiple MPI processes are available. Use
    this for tests that verify MPI operations, distributed computations, or
    inter-process communication.
  - `[GPU]` tests run only when GPU devices are available. Use this for tests that
    are meaningful and interesting on GPU hardware.

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
exercised by the tests. This feature uses "[source-based code
coverage](https://clang.llvm.org/docs/SourceBasedCodeCoverage.html)" for
LLVM-based compilers (Clang, AppleClang, IntelLLVM) and
[`gcov`](https://gcc.gnu.org/onlinedocs/gcc/Gcov.html) for GNU compilers.

Both methods can produce coverage data compatible with
[LCOV](https://github.com/linux-test-project/lcov), which provides standardized
processing and HTML report generation. In the sections below, we discuss how to
obtain a `coverage.info` that can be further processed by `lcov`. Before
proceeding, ensure `lcov` is installed on your system. If not, you can obtain it
from your package manager.

!!! warning "Before continuing"
    
    Make sure you have `lcov` installed. To check that it works, run
    
    ```sh
    lcov --version
    ```
    
    Also make sure that *Palace* was build with `PALACE_BUILD_WITH_COVERAGE`.
    If not, compile it with
    
    ```sh
    cmake -DPALACE_MFEM_USE_EXCEPTIONS=yes -DPALACE_BUILD_WITH_COVERAGE=yes ..
    ```

To obtain comprehensive coverage data, tests must be executed in three
configurations: `Serial`, `Parallel`, and `GPU` (if available). Each
configuration produces separate coverage files that are later merged into a
unified report.

*Palace* is a complex codebase that involves inlined and just-in-time compiled
functions. This can lead to minor inconsistencies in the coverage files. These
are not fatal, and ignoring them still produces good reports. For this reason,
most `lcov` commands you will see in the following sections have the
`--ignore-errors inconsistent` flag.

#### Coverage with GNU compilers

GNU compilers automatically generate raw coverage files during program
execution. The `lcov` tool can directly process these files without additional
conversion steps.

Run the test suite in all configurations and capture coverage data:

```sh
# Serial tests
bin/palace-unit-tests
lcov --capture --directory palace-build --ignore-errors inconsistent --output-file coverage1.info

# Parallel tests
mpirun -np 2 bin/palace-unit-tests
lcov --capture --directory palace-build --ignore-errors inconsistent --output-file coverage2.info

# GPU tests (if available)
# (command -v nvidia-smi is for your convenience, so that you can copy this block of code
#  and use it as-is even when a GPU is not available)
command -v nvidia-smi && bin/palace-unit-tests --device cuda
command -v nvidia-smi && lcov --capture --directory palace-build --ignore-errors inconsistent --output-file coverage3.info

# Merge all coverage files (works whether coverage3.info exists or not)
# The $() expands to --add-tracefile coverage1.info --add-tracefile coverage2.info --add-tracefile coverage3.info 
lcov $(for f in coverage[1-3].info; do echo "--add-tracefile $f"; done) --ignore-errors inconsistent --output-file coverage.info
```

#### Coverage with LLVM-based compilers

LLVM compilers use a different coverage format that requires additional
processing steps. The `LLVM_PROFILE_FILE` environment variable controls the
naming pattern for generated coverage files.

First, clean any existing coverage files and configure the profile naming:

```sh
rm -f coverage*.profraw
export LLVM_PROFILE_FILE="coverage-%p.profraw"
```

Execute the test suite in all configurations:

```sh
# Serial tests
bin/palace-unit-tests

# Parallel tests
mpirun -np 2 bin/palace-unit-tests

# GPU tests (if available)
# (command -v nvidia-smi is for your convenience, so that you can copy this block of code
#  and use it as-is even when a GPU is not available)
command -v nvidia-smi && bin/palace-unit-tests --device cuda
```

This generates multiple `coverage-<pid>.profraw` files, where `<pid>` represents
each process ID. These raw files must be merged and converted to LCOV format:

```sh
# Merge raw coverage files
llvm-profdata merge -output=coverage.profdata coverage*.profraw

# Export to JSON format
llvm-cov export -format text -instr-profile coverage.profdata bin/palace-unit-tests > coverage_llvm.json

# Convert to LCOV format
llvm2lcov --ignore-errors inconsistent --output coverage.info coverage_llvm.json
```

!!! note "Why two different coverage systems?"
    
    LLVM's source-based coverage provides more accurate instrumentation than the
    gcov-compatible mode, which attempts to emulate gcov behavior but has known
    reliability issues. The source-based approach offers better precision
    for complex codebases.

#### Generating coverage reports

Once you have produced `coverage.info` using either method, filter the data to
include only *Palace* source code:

```sh
lcov --extract coverage.info \
  '*/palace/drivers/*' \
  '*/palace/fem/*' \
  '*/palace/linalg/*' \
  '*/palace/models/*' \
  '*/palace/utils/*' \
  '*/test/unit/*' \
  --output-file coverage_filtered.info
```

Generate an HTML report for visualization:

```sh
genhtml coverage_filtered.info --output-directory coverage_html
```

Open `coverage_html/index.html` in your web browser to explore the interactive
coverage report, which shows line-by-line coverage statistics and identifies
untested code paths.

## Regression tests

In addition to unit tests, *Palace* comes with a series of regression tests.
Regression tests based on the provided example applications in the
[`examples/`](https://github.com/awslabs/palace/blob/main/examples/) directory
and verify that the code reproduces results in reference files stored in
[`test/examples/`](https://github.com/awslabs/palace/blob/main/test/examples/ref).

### Building and running example tests

#### Prerequisites

  - Julia
  - Palace executable in PATH or specified via environment variable

#### Setup

First, instantiate the Julia environment:

```bash
julia --project -e "using Pkg; Pkg.instantiate()"
```

You need to do this step only the very first time.

#### Environment Variables

  - `PALACE_TEST`: Path to *Palace* executable and optional arguments (default: "`palace`")
  - `NUM_PROC_TEST`: Number of MPI processes (default: number of physical cores)
  - `OMP_NUM_THREADS`: Number of OpenMP threads (default: 1)
  - `TEST_CASES`: Space-separated list of test cases to run (default: all examples)

#### Execution

Run all tests:

```bash
julia --project runtests.jl
```

Run specific test cases:

```bash
TEST_CASES="spheres rings" julia --project runtests.jl
```

Run with custom *Palace* executable:

```bash
PALACE_TEST="../../build/bin/palace" julia --project runtests.jl
```

Run with custom Palace executable and custom MPI launcher:

```bash
PALACE_TEST="../../build/bin/palace --launcher mpiexec" julia --project runtests.jl
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
