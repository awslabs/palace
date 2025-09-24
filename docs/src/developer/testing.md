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

First, follow the [quick start](../install.md#Quick-start) to build *Palace* with `CMake`. From the *Palace* root directory:

```bash
mkdir build && cd build
cmake -DPALACE_MFEM_USE_EXCEPTIONS=yes .. 
make -j palace-tests
```

The `PALACE_MFEM_USE_EXCEPTIONS` option is necessary to properly capture and test assertions.

Once the build completes, run the tests with:

```bash
palace-build/test/unit/unit-tests
```

This runs all the `Serial` tests. To run `Parallel` tests, call

```bash
mpirun -np 2 palace-build/test/unit/unit-tests
```

To run the `GPU` tests (this can be combined with `mpirun`):

```bash
palace-build/test/unit/unit-tests --device cuda
```

If you are interested in running only a subset of tests, you can use filters.
For example, to run the tests tagged with `[vector]`, call

```bash
palace-build/test/unit/unit-tests [vector]
```

To run a specific test by name:

```bash
palace-build/test/unit/unit-tests "Vector Sum - Real"
```

#### Benchmarks

The unit test application also includes a small number of benchmarks to compare
performance between MFEM's legacy assembly backend, MFEM's partial assembly
backend, and the specified libCEED backend (specified with the `--backend`
option, use `-h`/`--help` to list all command line options for the `unit-tests`
executable). These can be run using, for example:

```bash
palace-build/test/unit/unit-tests "[Benchmark]" --benchmark-samples 10
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
