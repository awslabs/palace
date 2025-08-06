```@raw html
<!--- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. --->
<!--- SPDX-License-Identifier: Apache-2.0 --->
```

# Developer Notes

## Style guide

Automated source code formatting is performed using
[`clang-format`](https://clang.llvm.org/docs/ClangFormat.html). Run:

```bash
./scripts/format_source
```

in the repository root directory to automatically use `clang-format` to format `C++` source
as well as [`JuliaFormatter.jl`](https://github.com/domluna/JuliaFormatter.jl) for Julia and
Markdown files. The script can be viewed
[in the repository](https://github.com/awslabs/palace/blob/main/scripts/format_source).

The following conventions also apply:

  - `PascalCase` for classes and function names.
  - Follow 'include what you use' (IWYU), with the include order dictated by the
    [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html#Names_and_Order_of_Includes).
    This order should be automatically enforced by the `clang-format`
    [style file](https://github.com/awslabs/palace/blob/main/.clang-format).
  - Code comments should be full sentences, with punctuation. At this time, no Doxygen API
    reference is generated and so comments generally do not need to conform to Doxygen
    syntax.

## Static analysis

During the `cmake` configuration step, defining the variables `ANALYZE_SOURCES_CLANG_TIDY`
and `ANALYZE_SOURCES_CPPCHECK` to `ON` will turn on static analysis using
[`clang-tidy`](https://clang.llvm.org/extra/clang-tidy/) and
[`cppcheck`](https://cppcheck.sourceforge.io/), respectively, during the build step. This
requires the executables to be installed and findable by CMake on your system.

## JSON Schema for configuration files

A JSON format [configuration file](config/config.md), for example named `config.json`, can
be validated against the provided Schema using:

```bash
./scripts/validate_config config.json
```

[This script](https://github.com/awslabs/palace/blob/main/scripts/validate_config) uses
Julia's [`JSONSchema.jl`](https://github.com/fredo-dedup/JSONSchema.jl) and the Schema
provided in [`scripts/schema/`](https://github.com/awslabs/palace/blob/main/scripts/schema)
to parse the configuration file and check that the fields are correctly specified. This
script and the associated Schema are also installed and can be accessed in
`<INSTALL_DIR>/bin`.

## Timing

Timing facilities are provided by the `Timer` and `BlockTimer` classes.

Creating a block as `BlockTimer b(idx)` where `idx` is a category like `CONSTRUCT`, `SOLVE`,
etc. will record time so long as `b` is in scope; however, timing may be interrupted by
creation of another `BlockTimer` object. It will resume whenever the new block is destroyed.
Only one category is timed at once. This enables functions to declare how calls within them
are timed without needing to know how timing may be done by the calling function.

The `BlockTimer` implementation relies upon a static member object of the `Timer` class,
which behaves as a stopwatch with some memory functions. It is the responsibility of this
`BlockTimer::timer` object to record time spent in each recorded category. Other `Timer`
objects may be created for local timing purposes, but these will not count toward time
reported at the end of a log file or in the metadata JSON.

For each *Palace* simulation, a table of timing statistics is printed out before the program
terminates. The minimum, maximum, and average times across all MPI processes are included in
the table. Timer categories are exclusive, they do not overlap. The categories for breaking
down the total simulation time are:

```
Initialization                  // < Time spent parsing the configuration file,
                                //   preprocessing and partitioning the mesh
Operator Construction           // < Time spent constructing finite element spaces and
                                //   setting up linear and bilinear forms on those spaces
  Wave Ports                    // < Time spent configuring and computing the modes
                                //   associated with wave port boundaries
Linear Solve                    // < Linear solver time
  Setup                         // < Setup time for linear solver and preconditioner
  Preconditioner                // < Preconditioner application time for linear solve
  Coarse Solve                  // < Coarse solve time for geometric multigrid
                                //   preconditioners
Time Stepping                   // < Time spent in the time integrator for transient
                                //   simulation types
Eigenvalue Solve                // < Time spent in the eigenvalue solver (not including
                                //   linear solves for the spectral transformation, or
                                //   divergence-free projection)
Div.-Free Projection            // < Time spent performing divergence-free projection (used
                                //   for eigenmode simulations)
PROM Construction               // < Offline PROM construction time for adaptive fast
                                //   frequency sweep
PROM Solve                      // < Online PROM solve and solution prolongation time for
                                //   adaptive fast frequency sweep
Estimation                      // < Time spent computing element-wise error estimates
  Construction                  // < Construction time for error estimation global linear
                                //   solve
  Solve                         // < Solve time for error estimation global linear solve
Adaptation                      // < Time spent performing adaptive mesh refinement (AMR)
  Rebalancing                   // < Rebalancing time for AMR simulations with rebalancing
Postprocessing                  // < Time spent in postprocessing once the field solution
                                //   has been computed
Disk IO                         // < Disk read/write time for loading the mesh file and
                                //   writing fields for visualization
-----------------------
Total                           // < Total simulation time
```

## Testing

We use [Catch2](https://github.com/catchorg/Catch2) to perform unit testing of the
[libCEED](https://libceed.org/en/latest/) integration in *Palace* against the legacy MFEM
assembly routines. The unit tests source code is located in the
[`test/unit/`](https://github.com/awslabs/palace/blob/main/test/unit/) directory, and can be
built from within the *Palace* build directory using `make unit-tests`, or from the
superbuild as `make palace-tests`. The unit tests can be accelerated using MPI and/or
OpenMP parallelism (when configured with `PALACE_WITH_OPENMP=ON`), but in all cases they are
only testing the local operator assembly on each process. The 2D and 3D sample meshes in
[`test/unit/mesh/`](https://github.com/awslabs/palace/blob/main/test/unit/mesh/) come from
the [MFEM repository](https://github.com/mfem/mfem/tree/master/data).

The unit test application also includes a small number of benchmarks to compare performance
between MFEM's legacy assembly backend, MFEM's partial assembly backend, and the specified
libCEED backend (specified with the `--backend` option, use `-h`/`--help` to list all
command line options for the `unit-tests` executable). These can be run using, for
example:

```bash
./unit-tests "[Benchmark]" --benchmark-samples 10
```

The unit tests are run automatically as part of the project's continuous integration (CI)
workflows. Also run as part of the CI are regression tests based on the provided example
applications in the [`examples/`](https://github.com/awslabs/palace/blob/main/examples/)
directory. These are executed based on the code in
[`test/examples/`](https://github.com/awslabs/palace/blob/main/test/examples/).
See [`README`](https://github.com/awslabs/palace/tree/main/test/examples/) in the test
folder for more information about how to run the reproducibility tests.

## Changelog

Code contributions should generally be accompanied by an entry in the
[changelog](https://github.com/awslabs/palace/blob/main/CHANGELOG.md).
