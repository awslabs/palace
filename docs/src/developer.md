```@raw html
<!--- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. --->
<!--- SPDX-License-Identifier: Apache-2.0 --->
```

# Developer Notes

## Style guide

Automated source code formatting is performed using [`clang-format`]
(https://clang.llvm.org/docs/ClangFormat.html). Run:

```bash
./scripts/format_source
```

in the repository root directory to automatically use `clang-format` to format `C++` source
as well as [`JuliaFormatter.jl`](https://github.com/domluna/JuliaFormatter.jl) for Julia and
Markdown files. The script can be viewed [in the repository]
(https://github.com/awslabs/palace/blob/main/scripts/format_source).

The following conventions also apply:

  - `PascalCase` for classes and function names.
  - Follow 'include what you use' (IWYU), with the include order dictated by the
    [Google C++ Style Guide]
    (https://google.github.io/styleguide/cppguide.html#Names_and_Order_of_Includes). This
    order should be automatically enforced by the `clang-format` [style file]
    (https://github.com/awslabs/palace/blob/main/.clang-format).
  - Code comments should be full sentences, with punctuation. At this time, no Doxygen API
    reference is generated and so comments generally do not need to conform to Doxygen
    syntax.

## Static analysis

During the `cmake` configuration step, definining the variables `ANALYZE_SOURCES_CLANG_TIDY`
and `ANALYZE_SOURCES_CPPCHECK` to `ON` will turn on static analysis using [`clang-tidy`]
(https://clang.llvm.org/extra/clang-tidy/) and [`cppcheck`]
(https://cppcheck.sourceforge.io/), respectively, during the build step. This requires the
executables to be installed and findable by CMake on your system.

## JSON Schema for configuration files

A JSON format [configuration file](config/config.md), for example named `config.json`, can
be validated against the provided Schema using:

```bash
./scripts/validate_config config.json
```

[This script](https://github.com/awslabs/palace/blob/main/scripts/validate_config) uses
Julia's [`JSONSchema.jl`](https://github.com/fredo-dedup/JSONSchema.jl) and the Schema
provided in [`scripts/schema/`]
(https://github.com/awslabs/palace/blob/main/scripts/schema) to parse the configuration
file and check that the fields are correctly specified. This script and the associated
Schema are also installed and can be accessed in `<INSTALL_DIR>/bin`.

## Timing

Timing facilities are provided by `BlockTimer`.

Creating a block as `BlockTimer b(idx)` where `idx` is a category like `CONSTRUCT`, `SOLVE`,
etc. will record time so long as `b` is in scope; however, timing may be interrupted by
creation of another `BlockTimer` object. It will resume whenever the new block is destroyed.
Only one category is timed at once. This enables functions can declare how calls within them
are timed without needing to know how timing may be done by the calling function.

The `BlockTimer` implementation relies upon a static member object of the `Timer` class,
which behaves as a stopwatch with some memory functions. It is the responsibility of this
`Timer BlockTimer::timer` object to record time spent in each recorded category. Other
`Timer` objects may be created for local timing purposes, but these will not count toward
time reported at the end of a log file or in the metadata JSON.

Internally, the number of times `BlockTimer b(idx)` was called for each possible value of
`idx` (timing category) can be queried with `BlockTimer::Timer().GetCounts()`.

## Changelog

Code contributions should generally be accompanied by an entry in the [changelog]
(https://github.com/awslabs/palace/blob/main/CHANGELOG.md).
