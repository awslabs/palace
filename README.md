<!--- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. --->
<!--- SPDX-License-Identifier: Apache-2.0 --->
# Palace: 3D Finite Element Solver for Computational Electromagnetics

[![CI (Linux)](https://github.com/awslabs/palace/actions/workflows/build-and-test-linux.yml/badge.svg)](https://github.com/awslabs/palace/actions/workflows/build-and-test-linux.yml)
[![CI (macOS)](https://github.com/awslabs/palace/actions/workflows/build-and-test-macos.yml/badge.svg)](https://github.com/awslabs/palace/actions/workflows/build-and-test-macos.yml)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://awslabs.github.io/palace/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://awslabs.github.io/palace/dev)

*Palace*, for **PA**rallel **LA**rge-scale **C**omputational **E**lectromagnetics, is an
open-source, parallel finite element code for full-wave 3D electromagnetic simulations in
the frequency or time domain, using the
[MFEM finite element discretization library](http://mfem.org) and
[libCEED library](https://github.com/CEED/libCEED) for efficient exascale discretizations.

## Key features

  - Eigenmode calculations with optional material or radiative loss including lumped
    impedance boundaries. Automatic postprocessing of energy-participation ratios (EPRs) for
    [circuit quantization](https://www.nature.com/articles/s41534-021-00461-8) and
    interface or bulk participation ratios for predicting dielectric loss.
  - Frequency domain driven simulations with surface current excitation and lumped or
    numeric wave port boundaries. Wideband frequency response calculation using uniform
    frequency space sampling or an adaptive fast frequency sweep algorithm.
  - Explicit or fully-implicit time domain solver for transient electromagnetic analysis.
  - Lumped capacitance and inductance matrix extraction via electrostatic and magnetostatic
    problem formulations.
  - Support for a wide range of mesh file formats for structured and unstructured meshes,
    with built-in uniform or region-based parallel mesh refinement.
  - Solution-based Adaptive Mesh Refinement (AMR) for all simulation types aside from
    transient. Nonconformal refinement is supported for all mesh types, and conformal
    refinement for simplex meshes.
  - Arbitrary high-order finite element spaces and curvilinear mesh support thanks to the
    [MFEM library](https://mfem.org/features/).
  - Scalable algorithms for the solution of linear systems of equations, including
    matrix-free $p$-multigrid utilizing
    [high-order operator partial assembly](https://mfem.org/performance/), parallel sparse
    direct solvers, and algebraic multigrid (AMG) preconditioners, for fast performance on
    platforms ranging from laptops to HPC systems.
  - Support for hardware acceleration using NVIDIA or AMD GPUs, including multi-GPU
    parallelism, using pure CUDA and HIP code as well as [MAGMA](https://icl.utk.edu/magma/)
    and other libraries.

## Getting started

*Palace* can be installed using the [Spack HPC package manager](https://spack.io/), with the
command `spack install palace`. Run `spack info palace` to get more information about the
available configuration options and dependencies.

Those wishing to work in a containerized environment may use the Singularity/Apptainer
recipe for *Palace* in [`singularity/`](./singularity) to build a container containing
*Palace* and all its dependencies.

Finally, instructions for obtaining *Palace* and building from source can be found in the
[documentation](https://awslabs.github.io/palace/dev/install/). As part of the CMake build
process, most dependencies are downloaded and installed automatically and thus an internet
connection is required.

System requirements:

  - CMake version 3.21 or later
  - C++17 compatible C++ compiler
  - C and Fortran (optional) compilers for dependency builds
  - MPI distribution
  - BLAS, LAPACK libraries
  - CUDA Toolkit or ROCm installation (optional, for GPU support only)

## Documentation

[https://awslabs.github.io/palace/](https://awslabs.github.io/palace/)

The documentation for *Palace* provides full instructions for building the solver and
running electromagnetic simulations.

### Building a local copy of the documentation

[Julia](https://julialang.org) with
[Documenter](https://documenter.juliadocs.org/) is required to build a local
version of the documentation. Obtain Julia following the [official
instructions](https://julialang.org/install/) and install Documenter by
instantiating the `docs` environment

```sh
julia --project=docs -e "using Pkg; Pkg.instantiate()"
```

Then, generate the documentation with `julia --project make.jl` from within the
[`docs/`](./docs) directory. An HTTP server is needed to visualize the
rendered documentation. The simple way is to start a server is with python:

```sh
cd docs/build && python -m http.server 8000
```

Then, navigate to `localhost:8000` with your browser.

## Examples

Some example applications including configuration files and meshes can be found in the
[`examples/`](./examples) directory. Complete tutorials for each example are available in
the [documentation](https://awslabs.github.io/palace/dev/examples/examples/).

## Changelog

Check out the [changelog](./CHANGELOG.md).

## Contributing

We welcome contributions to *Palace* including bug fixes, feature requests, etc. To get
started, check out our [contributing guidelines](CONTRIBUTING.md).

## Contact

*Palace* is developed by the Design and Simulation group in the AWS Center for Quantum
Computing (CQC). Please contact the development team at
[palace-maint@amazon.com](mailto:palace-maint@amazon.com) with any questions or comments, or
[open an issue](https://github.com/awslabs/palace/issues).

## License

This project is licensed under the [Apache-2.0 License](./LICENSE).

See [THIRD-PARTY-LICENSES](./THIRD-PARTY-LICENSES) and
[THIRD-PARTY-NOTICES](./THIRD-PARTY-NOTICES) for licenses and notices of
third-party software in this repository.
