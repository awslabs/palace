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

During the `cmake` configuration step, definining the variables `ANALYZE_SOURCES_CLANG_TIDY`
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
Only one category is timed at once. This enables functions can declare how calls within them
are timed without needing to know how timing may be done by the calling function.

The `BlockTimer` implementation relies upon a static member object of the `Timer` class,
which behaves as a stopwatch with some memory functions. It is the responsibility of this
`BlockTimer::timer` object to record time spent in each recorded category. Other `Timer`
objects may be created for local timing purposes, but these will not count toward time
reported at the end of a log file or in the metadata JSON.

## Testing

We use [Catch2](https://github.com/catchorg/Catch2) to perform unit testing of the
[libCEED](https://libceed.org/en/latest/) integration in Palace against the legacy MFEM
assembly routines. The unit tests source code is located in the
[`test/unit/`](https://github.com/awslabs/palace/blob/main/test/unit/) directory, and can be
built from within the *Palace* build directory using `make unit-tests`, or from the
superbuild as `make palace-tests`. The unit tests can be accelerated using MPI and/or
OpenMP parallelism(when configured with `PALACE_WITH_OPENMP=ON`), but in all cases they are
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

## Changelog

Code contributions should generally be accompanied by an entry in the
[changelog](https://github.com/awslabs/palace/blob/main/CHANGELOG.md).
