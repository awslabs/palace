```@raw html
<!--- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. --->
<!--- SPDX-License-Identifier: Apache-2.0 --->
```

# Developer Notes

## Units

For interactions with users (inputs and outputs), *Palace* mostly uses SI units,
and where it doesn't it uses metric units (SI with prefix modifiers e.g. km, um,
nJ). However, internally, there are three unit systems that are used within any
given *Palace* simulation.

First, the mesh comes with its own unit of length. Users specify the conversion
between one unit of length on the mesh to meters through the
[`config["Model"]["L0"]`](../config/problem.md#config%5B%Model%22%5D) parameter.

Second, *Palace* defines a system of non-dimensional units where the actual
solver operates.

Third, postprocessing occurs in a custom variant of SI units where the unit of
time is nanoseconds instead of seconds. This is more suitable for the problems
*Palace* is designed to solve.

The second and third systems are described by the `Units` class defined in
`units.hpp`.

### Non-dimensional unit system

The non-dimensional unit system is constructed as follows:

Length: Lengths are measured in units of the characteristic length `Lc`, which
is the mesh size in mesh units, unless `Lc` is manually specified in
[`config["Model"]["Lc"]`](../config/problem.md#config%5B%Model%22%5D). The extent
of the mesh is 1 in these units, with extent defined as the largest single
dimension (e.g., for a cuboid of size 10m x 20m x 40m, this would be 40m).

Time: One unit of non-dimensional time is the time it takes light to travel one
`Lc` (equivalently, the speed of light is 1 in these units).

Electromagnetic quantities: *Palace* uses
[Heaviside-Lorentz](https://en.wikipedia.org/wiki/Heaviside%E2%80%93Lorentz_units)
units, essentially meaning that `ε₀ = μ₀ = 1` (so that in addition to setting `c = 1`, the impedance of free space `Z₀` is also 1).

The final unit in the system is chosen so that the power flowing through a
surface of area `Lc²` perpendicular to the direction of propagation in vacuum
with a uniform magnetic field of strength 1 is exactly 1 (`P = H^2 * Z * L^2`,
with `Z` impedance). This choice is related to the normalization of power across
ports.

An equivalent mental model about non-dimensional quantities in *Palace* is to
think in terms of characteristic values instead of units: the non-dimensional
quantities are physical quantities normalized by reference values. These
reference values are effectively chosen to satisfy the constraints on space,
time and power given above. For example, the magnetic field strength `H` solved
for is actually `H/Hc`, where `Hc` satisfies `Hc² × Z₀ × Lc² = 1 W`, with `Z₀`
impedance of free space. Similarly, the non-dimensional electric field `E` is
`E/Ec` with `Ec = Z₀ × Hc`.

!!! note "Example: converting times and distances"
    
    Consider a spherical mesh with `L0` of 1 micron and radius of 3 in mesh units.
    Then:
    
      - Mesh size in mesh units: 6
      - Mesh size in non-dimensional units: 1
      - Mesh size in dimensional units: 6e-6
      - Mesh size in SI: 6e-6 meters
      - Speed of light in non-dimensional units: 1
      - Speed of light in dimensional units: 0.3
      - Speed of light in SI: 3e8 meters/seconds
      - Time for light to cross the sphere in non-dimensional units: 1
      - Time for light to cross the sphere in dimensional units: 2e-5
      - Time for light to cross the sphere in SI: 2e-14 seconds

!!! warning "Non-SI time units"
    
    The dimensional unit system uses nanoseconds for time, not seconds. This affects
    unit combinations: `units::POWER * units::TIME` is not equivalent to `units::ENERGY`.
    This choice of nanoseconds reflects a prioritization of the microwave engineering
    domain within *Palace*, as the typical frequencies in GHz (which correspond to time
    scales of nanoseconds).

!!! warning "Omega instead of Frequency"
    
    Internally, *Palace* uses the angular frequency `omega` instead of the regular
    frequency, with both names `omega` and `freq` being used.

!!! warning "Finite element considerations"
    
    Changing mesh units requires additional rescaling for fields defined on
    H(curl) and H(div) spaces to convert them to proper physical units. This is a
    function of the length scales inherent to the reference element transformations
    for Nédélec and Raviart Thomas spaces, used for H(curl) and H(div) respectively.
    See the `ScaleGridFunctions` function in `PostOperator` for details.

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

## Extending GPU code in *Palace*

*Palace* supports [GPU parallelization](../guide/parallelism.md), but not all the
files in *Palace* are compiled with a GPU-compatible compiler (e.g., `nvcc`).
The reason for this is that such compilers might not support all the features
*Palace* needs (e.g., support for `constexpr std::sqrt`), and also because they
tend to be slower. The list of files that contains code that has to run on the
device is defined by the `TARGET_SOURCES_DEVICE` CMake variable.

## JSON Schema for configuration files

A JSON format [configuration file](../config/config.md), for example named
`config.json`, can be validated against the provided Schema using:

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
  Far Fields                    // < Time spent computing surface integrals to extrapolate
                                //   near fields to far fields
  Paraview                      // < Processing and writing Paraview fields for visualization
  Grid function                 // < Processing and writing grid function fields for visualization
Disk IO                         // < Disk read/write time for loading the mesh file and
                                //   writing CSV fields
-----------------------
Total                           // < Total simulation time
```

## Profiling *Palace* on CPUs

A typical *Palace* simulation spends most of its time in libCEED kernels, which, in turn, executed `libsxmm` code on CPUs. Libsxmm generates code just-in-time to ensure it is the most performant on the given architecture and for the given problem. This code generation confuses most profilers. Luckily, [libsxmm](https://libxsmm.readthedocs.io/en/latest/libxsmm_prof/) can integrate with the VTune APIs to enable profiling of jitted functions as well.

!!! note "Using a different libCEED backend?"
    
    The following notes assume that you are using the `libxsmm` backend for libCEED. If this is not the case and you have a preference for other profilers, such as [perf](https://perfwiki.github.io/main/) or [HPCToolkit](https://hpctoolkit.org/), you can consider using them as well.

### Requirements

 1. Verify your system meets the [VTune system requirements](https://www.intel.com/content/www/us/en/developer/articles/system-requirements/vtune-profiler/2025-1.html)
 2. [Install VTune Profiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/vtune-profiler-download.html)
 3. Set VTune's root directory during Palace build (see below)
 4. Use `RelWithDebInfo` build type for optimal profiling results

### Building *Palace* with VTune Support

From the *Palace* root directory:

```bash
$ source /path/to/vtune-vars.sh # Replace with your VTune installation path
$ mkdir build && cd build
$ cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo
$ make -j $(nproc)
```

### Running a Basic Profile

A flamegraph can be created using either:

  - The VTune graphical user interface
  - The command line interface

Follow [Intel's official instructions](https://www.intel.com/content/www/us/en/docs/vtune-profiler/user-guide/2024-2/launch.html) for detailed usage information.

Below is a sample flamegraph resulted from profiling the rings example located at [`/examples/rings/rings.json`](https://github.com/awslabs/palace/blob/main/examples/rings/rings.json) with polynomial order P=6 (`Solver["Order"]: 6`).

```@raw html
<br/><p align="center">
  <img src="../../assets/profiling/vtune-flamegraph.png" width="100%" />
</p><br/>
```

**Note:** Set `Verbose: 0` in the JSON file to reduce terminal output. For profiling, use the binary executable (`build/bin/palace-*.bin`) rather than the `build/bin/palace` wrapper script.

For [MPI profiling](https://www.intel.com/content/www/us/en/docs/vtune-profiler/cookbook/2023-0/profiling-mpi-applications.html) with multiple processes:

```bash
# Run with 10 MPI processes
$ cd /path/to/palace/examples/rings
$ mpirun -n 10 -gtool "vtune -collect hpc-performance --app-working-dir=$(pwd) -result-dir $(pwd)/vtune-mpi:0-9" /path/to/palace/build/bin/palace-x86_64.bin rings.json

# Open the results from VTune GUI
$ vtune-gui vtune-mpi.*
```

## Changelog

Code contributions should generally be accompanied by an entry in the
[changelog](https://github.com/awslabs/palace/blob/main/CHANGELOG.md).
