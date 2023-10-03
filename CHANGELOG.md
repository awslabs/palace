<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
# Changelog

> Note: *Palace* is under active initial development, pre-v1.0. Functionality and interfaces
> may change rapidly as development progresses.

The format of this changelog is based on
[Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to
[Semantic Versioning](https://semver.org/).

## In progress

  - Changed implementation of complex-valued linear algebra to use new `ComplexVector` and
    `ComplexOperator` types, which are based on the underlying `mfem::Vector` and
    `mfem::Operator` classes, instead of PETSc. PETSc is now fully optional and only
    required when SLEPc eigenvalue solver support is requested. Krylov solvers for real- and
    complex-valued linear systems are implemented via the built-in `IterativeSolver`
    classes.
  - Changed implementation of PROMs for adaptive fast frequency sweep to use the Eigen
    library for sequential dense linear algebra.
  - Changed implementation of numeric wave ports to use MFEM's `SubMesh` functionality. As
    of [#3379](https://github.com/mfem/mfem/pull/3379) in MFEM, this has full ND and RT
    basis support. For now, support for nonconforming mesh boundaries is limited.
  - Added support for operator partial assembly for high-order finite element spaces based
    on libCEED for non-tensor product element meshes. This option is disabled by default,
    but can be activated using `config["Solver"]["PartialAssemblyOrder"]` set to some number
    less than `"Order"` and `config["Solver"]["Device"]: "ceed-cpu"`.
  - Added support for non axis aligned lumped ports and current sources. Key words `"X"`,
    `"Y"`, `"Z"` and `"R"`, with optional prefix `"+"` or `"-"` still work, but now
    directions can be specified as vectors with 3 components. Users will be warned, and
    ultimately errored, if the specified directions do not agree with axis directions
    discovered from the geometry.
  - Fixed bugs for simulations using tetrahedral meshes associated with unexpected mesh
    toplogy changes during parallel mesh construction.
  - Added improved `Timer` and `BlockTimer` classes with more timing categories for
    reporting simulation runtime.
  - Added build dependencies on [libCEED](https://github.com/CEED/libCEED) and
    [LIBXSMM](https://github.com/libxsmm/libxsmm) to support operator partial assembly (CPU-
    based for now).
  - Added improved OpenMP support in `palace` wrapper script and CI tests.
  - Added Apptainer/Singularity container build definition for Palace.
  - Added flux-based error estimation, reported in `error-estimate.csv`. This computes the
    difference between the numerical gradient (electrostatics) or curl (otherwise) of the
    solution, and a smoother approximation obtained through a global mass matrix inversion.
    The results are reported in `error-estimates.csv` within the `"Output"` folder.

## [0.11.2] - 2023-07-14

  - Fixed a regression bug affecting meshes which have domain elements which are not
    assigned material properties in the configuration file.
  - Changed layout and names of `palace/` source directory for better organization.
  - Added many updates to build system: Removed use of Git submodules to download
    dependencies relying instead directly on CMake's ExternalProject, patch GSLIB dependency
    for shared library builds, add CI tests with ARPACK-NG instead of SLEPc, update all
    dependency versions including MFEM to incorporate bug fixes and improvements. This
    affects the Spack package recipe, so a new recipe is distributed as part of Palace in
    in `spack/` which will keep compatibility with `main` while changes are upstreamed to
    the built-in Spack repository.
  - Added a basic Makefile with some useful targets for development.

## [0.11.1] - 2023-05-03

  - Fixed a bug for interface dielectric loss postprocessing including when using
    `--dry-run`.
  - Fixed a regression bug affecting second-order absorbing boundary conditions.
  - Fixed a bug when the number of processors was large enough that some subdomains own no
    true dofs.
  - Fixed memory bug in destruction of SuperLU solver.
  - Fixed some typos in the documentation.
  - Fixed a bug in the cavity convergence study example for the coarsest hexahedral mesh
    case, as well as some other minor Julia fixes in the `examples/` and `scripts/`
    directories.
  - Added updates to superbuild including better ARPACK support, though still experimental
    compared to SLEPc.
  - Added updated submodules for superbuild.

## [0.11.0] - 2023-01-26

  - Initial public release on GitHub.
  - Added support for anisotropic materials through the use of `"MaterialAxes"` under
    `"Materials"`. These material axes allow for specifying symmetric material property
    tensors through a sum of weighted outer products of normal vectors. This can be used in
    conjunction with scalar material properties or multiple different anisotropic
    properties, though all properties are required to use the same set of basis vectors.
    Demonstration given for the coplanar waveguide example which now utilizes a sapphire
    substrate.
  - Added to default postprocessing outputs: Bulk and interface dielectric loss writes
    energy-participation ratios in addition to quality factors associated with lossy
    regions, IO coupling $\kappa$ in addition to quality factor for eigenmode
    simulations, lumped element energy for transient simulations similar to eigenmode and
    frequency domain driven simulations.
  - Changed configuration file syntax to simplify support for multielement lumped ports,
    where each index for a multielement port should now use the `"Elements"` object to
    describe the attributes for each port element.
  - Changed configuration file keywords to better reflect meaning:
    `config["Boundaries"]["Current"]` => `config["Boundaries"]["SurfaceCurrent"]`, and
    `"UseGMG"`, `"UsePCShifted"`, `"MGCycleIts"`, and `"MGSmoothIts"` under
    `config["Solver"]["Linear"]`.
  - Changed geometric multigrid implementation to generalize across problem types (added
    for electrostatics) and use Jacobi-preconditioned Chebyshev smoothing with optional
    auxiliary space smoothing at each multigrid level as well. The auxiliary space matrices
    are constructed directly in H1 now to allow for matrix-free/partial-assembly support
    eventually. Geometric multigrid is turned on by default.
  - Added structured simulation metadata output in the form of a JSON file `palace.json`
    written to the directory specified by `config["Problem"]["Output"]`.
  - Added JSON Schema files for the configuration file format as well as a helper script
    to check a provided configuration file against the schema.
  - Added optional interface to GSLIB library which enables `"Probe"` functionality to
    sample the computed electric and magnetic fields at points in space.
  - Added preliminary support for builds using Spack. This includes an option in the build
    system for user-supplied dependencies as an alternative to the superbuild
    configuration.
  - Added updated submodules for superbuild, fixing a bug in SuperLU_DIST solver causing
    communication hangs for certain numbers of processors.

## [0.10.0] - 2022-10-04

  - Added interfaces to ARPACK and FEAST eigenvalue solvers.
  - Added divergence-free projection for eliminating spurious DC modes for eigenmode
    solves.
  - Added option for visualizing fields on mesh boundaries, output alongside the full 3D
    solution fields for ParaVeiew visualization.
  - Added real and imaginary fields output for complex-valued phasors, and electric and
    magnetic energy density fields.
  - Added convergence study for the cavity example application, and example Julia code for
    automated mesh generation, running of the solver, and postprocessing.
  - Fixed bugs in mesh preprocessing related to precision of nodal coordinates when
    distributing a serial mesh to each process.
  - Added option for `"Driven"` and `"Transient"` simulations to accelerate postprocessing
    by only considering port boundaries.
  - Added `-dry-run`/`--dry-run` command line option to check configuration file for errors
    and exit.
  - Changed dependency build process to rely less on PETSc build system, no longer give
    option to build BLAS/LAPACK libraries from source, added detection for OpenBLAS and AMD
    Optimizing CPU Libraries (AOCL) blis/libflame.
  - Fixed issues when compiling with Intel or Clang compilers.
  - Added CI test builds to test various build options and added easier default build
    options in CMake configuration.
  - Added `clang-format` and `JuliaFormatter` configurations for automated C++, Julia, and
    Markdown code formatting.

## [0.9.0] - 2022-07-11

  - Added region-based mesh refinement with box or sphere refinement regions.
  - Added new sparse direct solver interfaces to SuperLU_DIST and STRUMPACK libraries, to
    complement the existing MUMPS interface.
  - Changed configuration file keywords for linear solver and preconditioner parameters to
    make more options available to the user and clarify them.
  - Added check that all external boundary surface attributes must have an associated
    boundary condition specified in the configuration file (adds new `"PMC"` and
    `"ZeroCharge"` boundaries for natural boundary conditions).
  - Added advanced configuration options for GMRES/FGMRES/SLEPc eigenvalue solver
    orthogonalization, making CGS2 (classical GS with iterative refinement) the default.

## [0.8.0] - 2022-06-06

  - Added new comprehensive example applications and detailed tutorials in the
    documentation.
  - Added CI pipeline for build and regression testing.
  - Fixed bugs for shared library builds, which are now the default for dependency builds.
  - Added complex lumped port voltage and current output for computing impedance and
    admittance parameters from frequency domain driven simulations.
  - Fixed a bug in mesh parsers for COMSOL and Nastran formats for high-order hexahedra,
    prism, and pyramid element types.
  - Changed adaptive frequency sweep implementation to introduce optimizations which, in
    particular, speed up driven simulations with wave ports.
  - Fixed bug related to domain dielectric postprocessing and surface dielectric loss
    postprocessing now automatically detects correct side for substrate-air interfaces
    (always evaluates on vacuum side).

## [0.7.0] - 2022-03-25

  - Added adaptive fast frequency sweep implementation for frequency domain driven
    simulations, activated with `config["Solver"]["Driven"]["AdaptiveTol"] > 0.0`.
  - Changed ASCII file postprocessing to write out .csv files instead of generic whitespace
    delimited .txt for easier parsing.
  - Changed field output postprocessing for visualization to save the mesh in the original
    mesh length units (rather than nondimensionalized ones).
  - Added electric and magnetic field energy output to files as function of mode, frequency,
    time, or solution step, in `domain-E.csv`.
  - Added option to specify circuit R/L/C parameters for lumped port boundaries as an
    alternative to the surface Rs/Ls/Cs ones.
  - Added improved and expanded documentation.
  - Added new MFEM features, including support for ND and RT spaces on wedge or prism
    elements and bug fixes.

## [0.6.0] - 2022-02-17

  - Changed implementation of boundary conditions and excitations to do a better job
    separating spatial discretization and algebraic operators. This is also in preparation
    for future additions of boundary conditions which are not quadratic in the frequency
    such as wave ports, as well as linear algebra improvements.
  - Added a `config["Boundaries"]["Current"]` configuration file parameter to specify a
    surface current excitation, and removed the lumped port `"Termination"` option as it is
    no longer relevant. Magnetostatics simulations use this `"Current"` excitation rather
    than `"Port"` now.
  - Changed lumped port and surface current excitations to now use a unit circuit voltage
    and circuit current excitation, rather than unit electric field or surface current
    density.
  - Added numeric wave port boundaries for frequency domain driven simulations. To
    differentiate, lumped port boundaries are now specified using the keyword
    `"LumpedPort"` instead of just `"Port"` in the configuration file.
  - Fixed farfield boundaries to account for non-vacuum materials at the boundary and added
    second-order option for absorbing boundary condition (ABC).
  - Added finite conductivity surface impedance boundary condition, with optional thickness
    correction.
  - Added dielectric loss calculation specializations for metal-air, metal-substrate, and
    substrate-air interfaces in accordance with
    [this paper](https://aip.scitation.org/doi/10.1063/1.3637047).
  - Changed implementation of frequency domain linear algebra to eliminate the use of PETSc
    matrix types, relying instead on wrappers of the real and imaginary Hypre matrices
    making up a complex operator.
  - Added geometric multigrid AMS solvers (in two variants) using h- or p-coarsening, which
    employs either a specialized smoother for Nedelec space coarsening or geometric
    coarsening for the auxiliary space solves.
  - Changed build process to include PETSc as a submodule and configure as part the standard
    CMake build process.

## [0.5.0] - 2021-11-08

  - Changed postprocessing implementation which leads to some configuration file changes,
    namely `config["Domains"]["Postprocessing"]` and
    `config["Boundaries"]["Postprocessing"]`.
  - Changed some filenames for postprocessed quantities to be more intuitive.
  - Added dielectric (bulk and interface) loss calculations for all simulation types.
    Correctly handles two-sided internal surfaces using `"Side"` specification.
  - Added mutual capacitance matrix extraction to electrostatic solver.
  - Added built-in mesh support for COMSOL (`.mphbin` and `.mphtxt`) and Nastran (`.nas` and
    `.bdf`) mesh formats.
  - Removed Gmsh dependency for mesh format conversion.
  - Added some improvements to the default build flags and dependency build options.

## [0.4.0] - 2021-09-27

  - Added `"Capacitance"` and `"Inductance"` options to extract surface charges and fluxes
    in frequency and time domain simulation types.
  - Changed `"Transient"` solver implementation to use real-valued scalars (user must use a
    domain electrical conductivity as opposed to a loss tangent for loss modeling), for a
    roughly 2x performance improvement.
  - Added new MFEM features, including upgrades to high-order Nedelec element spaces on
    tetrahedra which no longer require global dof reordering.
  - Changed system matrix assembly to use Hypre P^T A P within MFEM and keep the full dense
    element matrices when the coefficient is nonzero (`skip_zeros = 0` always). The
    stiffness and mass matrices thus always share the same sparsity pattern now.

## [0.3.0] - 2021-09-02

  - Added electrostatic solver for capacitance matrix extraction.
  - Added `"Termination"` option for ports to allow for open or short circuiting but still
    enabling port postprocessing.
  - Changed timer to separate out time spent for disk I/O during initialization and
    postprocessing.
  - Added options to job wrappers for new instance types and subnet ID option.
  - Added changes from Hypre and MFEM dependencies to include improvements and bug fixes
    from those libraries.

## [0.2.0] - 2021-06-25

  - Added a proper [changelog](./CHANGELOG.md) for the project.
  - Added time domain solver with various excitation options to accompany existing frequency
    domain driven and eigenmode simulation types.
  - Added wrapper scripts in [`bin/`](./bin) for launching parallel simulations either
    interactively or as batch jobs using PBS.
  - Changed JSON configuration file keywords for frequency and time intervals.
  - Added domain conductivity models for normal and superconducting metals.
  - Fixed a bug for coaxial lumped port postprocessing.
  - Added visualization of surface currents on domain boundaries to accompany existing
    $\bm {E}$- and $\bm{B}$-field visualization.
  - Changed default length scale for nondimensionalization.
  - Changed default filenames for port-related postprocessed quantities.
  - Fixed many smaller issues and added various other improvements.

## [0.1.0] - 2020-11-03

  - Initial development release.
