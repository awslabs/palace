<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
# Palace: 3D Finite Element Solver for Computational Electromagnetics

[![CI (Linux)](https://github.com/awslabs/palace/actions/workflows/build-and-test-linux.yml/badge.svg)](https://github.com/awslabs/palace/actions/workflows/build-and-test-linux.yml)
[![CI (macOS)](https://github.com/awslabs/palace/actions/workflows/build-and-test-macos.yml/badge.svg)](https://github.com/awslabs/palace/actions/workflows/build-and-test-macos.yml)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://awslabs.github.io/palace/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://awslabs.github.io/palace/dev)

*Palace*, for **PA**rallel **LA**rge-scale **C**omputational **E**lectromagnetics, is an
open-source parallel finite element code for full-wave 3D electromagnetic simulations in
the frequency or time domain, using the [MFEM](http://mfem.org) finite element
discretization library.

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
  - Arbitrary high-order finite element spaces and curvilinear mesh support thanks to the
    [MFEM library](https://mfem.org/features/).
  - Scalable algorithms for the solution of linear systems of equations, including geometric
    multigrid (GMG), parallel sparse direct solvers, and algebraic multigrid
    (AMG) preconditioners, for fast performance on platforms ranging from laptops to HPC
    systems.

## Getting started

*Palace* can be installed using the [Spack](https://spack.io/) HPC package manager, with the
command `spack install palace`. Run `spack info palace` to get more information about the
available configuration options and dependencies.

Instructions for obtaining *Palace* and building from source can be found in the
[documentation](https://awslabs.github.io/palace/stable/install/). Git submodules are used
to manage the builds of most dependencies, so the repository should be cloned using the
`--recurse-submodules` flag, or alternatively `git submodule update --init` can be run in
order to download the dependency source code for an existing clone.

System requirements:

  - CMake version 3.13 or later
  - C++17 compatible C++ compiler
  - C and (optionally) Fortran compilers for dependency builds
  - MPI distribution
  - BLAS, LAPACK libraries

## Documentation

[https://awslabs.github.io/palace/](https://awslabs.github.io/palace/)

The documentation for *Palace* provides full instructions for building the solver and
running electromagnetic simulations.

To build a local version of the documentation, run `julia make.jl` from within the
[`docs/`](./docs) directory.

## Examples

Some example applications including configuration files and meshes can be found in the
[`examples/`](./examples) directory. Complete tutorials for each example are available in
the [documentation](https://awslabs.github.io/palace/stable/examples/examples/).

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
