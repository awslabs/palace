<!--- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. --->
<!--- SPDX-License-Identifier: Apache-2.0 --->
# Changelog

> Note: *Palace* is under active initial development, pre-v1.0. Functionality and interfaces
> may change rapidly as development progresses.

The format of this changelog is based on
[Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to
[Semantic Versioning](https://semver.org/).

## In progress

#### New Features

  - Added support for nonlinear eigenvalue problems arising from frequency-dependent
    boundary conditions. Two nonlinear eigensolvers are now available and can be specified
    by setting the `config["Solver"]["Eigenmode"]["NonlinearType"]` option to `"Hybrid"`
    (default) or `"SLP"`. The nonlinear eigensolver will automatically be used if
    frequency-dependent boundary conditions are used. [PR
    467](https://github.com/awslabs/palace/pull/467).
  - Added support for extraction of electric fields in the radiative zone.
    Consult the
    [documentation](https://awslabs.github.io/palace/dev/features/farfield) for
    additional information. [PR
    449](https://github.com/awslabs/palace/pull/449).
  - Added support for more granular tests and for computing the testing coverage.
    Consult the
    [documentation](https://awslabs.github.io/palace/dev/developer/testing/) for
    additional information. [PR
    398](https://github.com/awslabs/palace/pull/398), [PR
    480](https://github.com/awslabs/palace/pull/480).
  - Added option to export field solutions as MFEM grid functions to visualize with
    [GLVis](https://glvis.org). Output formats can be specified in `config["Problem"]["OutputFormats"]`.
    [PR 518](https://github.com/awslabs/palace/pull/518).
  - Added an option to drop small entries (below machine epsilon) from the matrix used in the sparse
    direct solver. This can be specified with `config["Solver"]["Linear"]["DropSmallEntries"]`.
    [PR 476](https://github.com/awslabs/palace/pull/476).

#### Interface Changes

  - `config["Boundaries"]["Periodic"]` is now a dictionary where all periodic boundary
    pairs, built into the mesh file or not, should be specified in
    `config["Boundaries"]["Periodic"]["BoundaryPairs"]` and a single global Floquet wave
    vector can be specified in `config["Boundaries"]["Periodic"]["FloquetWaveVector"]`. [PR
    471](https://github.com/awslabs/palace/pull/471).
  - Changed Paraview and GLVis output fields from nondimensional to SI units
    [PR 532](https://github.com/awslabs/palace/pull/532).
  - The name of the unit test executable is now `palace-unit-tests` (instead of `unit-tests`)
    [PR 549](https://github.com/awslabs/palace/pull/549).

#### Bug Fixes

  - Changed wave port eigenproblem shift and sorting to fix an issue with the mode ordering.
    The first mode now has the largest propagation constant, closest to the TEM limit, and
    subsequent modes are ordered by decreasing propagation constant. [PR
    448](https://github.com/awslabs/palace/pull/448).
  - Fixed an issue where Gmsh meshes with built-in periodicity (specified in the mesh file)
    were failing. [PR
    471](https://github.com/awslabs/palace/pull/471).
  - Fixed bug where a mesh from a previous nonconformal adaptation could not be loaded to
    use in a non-amr simulation. [PR
    497](https://github.com/awslabs/palace/pull/497).
  - Fixed bug where `"CrackInternalBoundaryElements"` would result in incorrect results for
    some lumped port boundary conditions. [PR
    505](https://github.com/awslabs/palace/pull/505).
  - Changed the sign of the Floquet phase factor from $\exp(+ik\cdot x)$ to $\exp(-ik\cdot x)$.
    [PR 510](https://github.com/awslabs/palace/pull/510) and [PR 526](https://github.com/awslabs/palace/pull/526).
  - Update EM constants to CODATA Recommended Values of the Fundamental Physical Constants 2022
    [PR 525](https://github.com/awslabs/palace/pull/525).
  - Scaled Rs/Ls/Cs of impedance boundary conditions affected by mesh cracking, fixing bug where
    `"CrackInternalBoundaryElements"` would lead to incorrect results. [PR 544](https://github.com/awslabs/palace/pull/544).

## [0.14.0] - 2025-08-20

  - Added `--version` command line flag for displaying Palace version information.
  - Fixed a small regression bug for boundary postprocessing when specifying
    `"Side": "LargerRefractiveIndex"`, introduced as part of v0.13.0.
  - Added an improvement to numeric wave ports to avoid targeting evanescent modes at
    higher operating frequencies. Also finite conductivity boundaries
    (`config["Boundaries"]["Conductivity"]`) are automatically marked as PEC for the wave
    port mode solve (previously these were marked as PMC unless specified under
    `"WavePortPEC"`).
  - Fixed a bug in divergence-free projection for problems without essential or mixed
    boundary conditions.
  - Added `"MakeSimplex"` and `"MakeHexahedral"` mesh options to convert an input mesh to
    all tetrahedra or all hexahedra. Also adds `"SerialUniformLevels"` option to
    `config["Model"]["Refinement"]` for testing or debugging.
  - Added `config["Model"]["CrackInternalBoundaryElements"]` which will separate or "crack" the mesh
    along all internal boundaries. This improves the performance of error estimation and AMR
    as the recovered smooth fields do not enforce additional erroneous continuity at
    internal boundaries. This will change the default behaviour in the case of internal
    impedance boundary conditions, and can be disabled by setting this option to false.
  - Added support for exact periodic boundary conditions, these can be specified as part of
    the mesh file (where supported by the format) or by specification of `"DonorAttributes"`
    and `"ReceiverAttributes"` in `config["Boundaries"]["Periodic"]` which will attempt to
    match the mesh on the boundaries specified by the donor and receiver attributes. This is
    only possible the mesh on the donor and receiver match exactly, non-matching meshes are
    not supported.
  - Exposed configuration linear solver and eigen solver options for the wave port
    subproblem. These can now be specified as part of the `config["Boundaries"]["WavePort"]`
    configuration. The defaults align with the previously hardcoded values.
  - Nonconformal adaptation is now supported for WavePort boundary conditions. This was
    achieved through a patch applied to MFEM to support `mfem::ParSubMesh` on external
    nonconformal surface subdomains.
  - Added adaptive time-stepping capability for transient simulations. The new ODE integrators
    rely on the SUNDIALS library and can be specified by setting the
    `config["Solver"]["Transient"]["Type"]` option to `"CVODE"` or `"ARKODE"`.
  - Added an option to use the complex-valued system matrix for the coarse level solve (sparse
    direct solve) instead of the real-valued approximation. This can be specified with
    `config["Solver"]["Linear"]["ComplexCoarseSolve"]`.
  - Fix bug in London equation implementation where a curl-curl term was added to the
    stiffness operator instead of a mass term.
  - Added support for Floquet periodic boundary conditions with phase-delay constraints.
    The Floquet wave vector can be specified along with periodic boundaries in the
    `config["Boundaries"]["Periodic"]` configuration.
  - Fixed bug where the default `config["Model"]["CrackDisplacementFactor"]` would cause
    lumped port discovery to fail by reducing default from `1e-3` to `1e-12`.
  - Fixed bug in Nastran mesh reader where carriage returns (`\r`) in the mesh file could
    cause a failure to read the mesh.
  - Changed post-processing, so that appropriate measurements are always written to disk as CSV
    files. This is a breaking change, since users are now required to specify a valid output
    folder in `config["Problem"]["Output"]`. Previously an empty string `""` would suppress file
    printing. ParaView printing is still controlled and suppressed by the `"Save"` or `"SaveStep"`
    options in the `config["Solver"][...]`.
  - Changed measurement and printing during post-processing substantially. This change is not
    user-facing, but will enable future multi-excitation support. All measurements are now
    performed as part of the `PostOperator` class, which is templated on solver type. The public
    interface of the `PostOperator` has been simplified: measurements can no longer be called
    individually, but are grouped in a `MeasurePrintAll` function. The actual printing of
    measurements has moved out of individual solvers and is orchestrated by the `PostOperatorCSV`
    class. Other helper classes include a small `Table` class (for data storage and formatting)
    and a light wrapper `TableWithCSVFile` (with file interaction).
  - Changed unit conversion interface: this was moved out of `IOData` into a separate `Units` class.
  - Changed csv output to print more digits.
  - Changed boundary indexing to enforce that the value of `Index` is unique between `LumpedPort`,
    `WavePort`, `SurfaceCurrent` and `Terminal`. This is a breaking change — e.g. currently a lumped
    port and a wave port could share an index.
  - Added multi-excitation support. `Excitation` in `LumpedPort` and `WavePort` boundary conditions
    can now be specified either as bools or excitation integers. Allow driven simulations to have
    multiple excitations indices in same configuration file to be simulated consecutively during the
    same Palace run. The measurement output of multiple excitations is printed to the same csv files
    with distinguish columns being post-indexed by the excitation index.
  - Changed `RomOperator` to split out rational interpolation component into the separate
    `MinimalRationalInterpolation` class to support multiple excitations. The MRI is unique to each
    excitation, but the PROM is shared between excitations.
  - Added `PortExcitations` to manage excitation pattern and print it to json metadata.
  - Update spack configuration to address `0.12` API breakage between libxsmm and libCEED.
  - Support spack vended libCEED and GSLIB builds.
  - Fix bug where `MakeSimplex` would fail for higher order meshes.
  - Fix bug in STRUMPACK build where dependency SCALAPACK library install directory could
    change based on build environment.
  - Introduce `"Samples"` mechanism for driven simulation, which allows for specifying a
    range of sample frequencies by combining `"Linear"`, `"Log"` and `"Point"`
    specifications. These are joined together effectively allowing for variable resolution
    sampling (for examples a coarse sampling supported by a fine sampling in a subset
    region). This is in addition to the prior existing `"MinFreq"`, `"MaxFreq"`,
    `"FreqStep"` interface which is now implicitly converted to a `"Linear"` sample scheme
    internally.
  - Introduce `"Save"` keyword for `"Driven"` simulation type, which will save specific
    frequency choices as opposed to the existing (and maintained `"SaveStep"` based on a
    regular sampling). Only frequencies contained within the set of `"Samples"` are
    supported, as no interpolation is performed.
  - Fix bugs in post processing where `E_cap`, `E_ind`, and `mode_port_kappa` and other
    dependent quantities were not dimensionalized correctly.
  - Refactor `PostOperator` usage of `Measurement` to be entirely non-dimensional, until
    `PostOperatorCSV` which dimensionalizes all measurements before writing to file. Reduces
    the risk of mixed unit bugs throughout `PostOperator`.
  - Fix bug in using `"MakeSimplex"` which would cause undefined behaviour for higher order
    meshes.
  - Fix bug when combining OpenMP and GPU builds in reduction operations over `Vector`.
  - Fix race condition that would affect OpenMP parallelism with periodic boundaries (exact
    and Floquet).
  - Fix race condition in `mfem::DenseTensor::operator()` with OpenMP, due to class member
    variable access.
  - Fix race condition in `DofToQuad` methods within mfem.
  - Normalize eigenmodes so their mean phase is a positive real number.
  - Added `scnlib` as a dependency to Palace.
  - Added parsing ability of existing csv files into the `Table` class. Fix `"Restart"`
    behaviour with the multi-excitation feature of driven solver.

## [0.13.0] - 2024-05-20

  - Changed default value of `config["Solver"]["PartialAssemblyOrder"]` in order to activate
    operator partial assembly by default for all operators in all simulation types.
  - Changed the normalization of computed eigenmodes for consistency across different domain
    decompositions. Eigenvectors are now normalized with respect to the mass matrix for unit
    domain electric field energy.
  - Added documentation for various timer categories and improved timing breakdown of
    various sections of a simulation.
  - Changed mesh files for the cavity and CPW examples, including prism, hexahedral, and
    tetrahedral meshes for the cylindrical cavity and correcting the wave port dimensions
    for the coplanar wave guide.
  - Fixed a few bugs and issues in the implementation of numeric wave ports for driven
    simulations.
  - Added GPU support for *Palace* via its dependencies, and added the
    `config["Solver"]["Device"]` and `config["Solver"]["Backend"]` options for runtime
    configuration of the MFEM device (`"CPU"` or `"GPU"`) and libCEED backend, with suitable
    defaults for users.
  - Added a new section to the documentation on
    [Parallelism and GPU support](https://awslabs.github.io/palace/dev/guide/parallelism/).
  - Removed use of `mfem::SparseMatrix` and replaced with HYPRE's `hypre_CSRMatrix` when
    needed for full assembly, wrapped as `palace::hypre::HypreCSRMatrix`.
  - Added `"Active"` configuration file parameter for lumped and wave port boundaries to
    disable the associated boundary condition and only use the surface for postprocessing.
  - Changed the smooth flux space for the electrostatic error estimator to fix performance
    on problems with material interfaces.
  - Fixed error estimation bug affecting time-dependent simulation types (driven, transient,
    eigenmode) where the recovery of the electric flux density also needs to be taken into
    account in addition to the magnetic field.
  - Fixed a bug related to mesh cleaning for unspecified domains and mesh partitioning.
  - Change computation of domain energy postprocessing for electrostatic and magnetostatic
    simulation types in order to improve performance.
  - Fixed a bug when computing the energy associated with lumped elements with more than
    one nonzero R, L, or C. This also affects the inductive EPR for lumped inductors with
    and associated parallel capacitance.
  - Fixed a bug for coaxial lumped ports which led to incorrect extraction of the geometric
    parameters, especially when coarsely-meshed or non-axis-aligned.
  - Added boundary postprocessing functionality for surface flux including electric,
    magnetic, and power given by the Poynting vector. This results in some breaking changes
    to the configuration file interface, see
    `config["Boundaries"]["Postprocessing"]["SurfaceFlux"]` and
    `config["Boundaries"]["Postprocessing"]["Dielectric"]`. In addition, related
    configuration file keyword changes to for consistency were made to
    `config["Domains"]["Postprocessing"]["Probe"]` and
    `config["Model"]["Refinement"]["Boxes"]`.
  - Fixed a bug in MFEM for nonconformal AMR meshes with internal boundaries affecting
    non-homogeneous Dirichlet boundary conditions for electrostatic simulations (see
    [#236](https://github.com/awslabs/palace/pull/236)).

## [0.12.0] - 2023-12-21

  - Added support for operator partial assembly for high-order finite element spaces based
    on libCEED for mixed and non-tensor product element meshes. This option is disabled by
    default, but can be activated using `config["Solver"]["PartialAssemblyOrder"]` set to
    some number less than or equal to `"Order"`.
  - Added flux-based error estimation, reported in `error-estimate.csv`. This computes the
    difference between the numerical gradient (electrostatics) or curl (otherwise) of the
    solution, and a smoother approximation obtained through a global mass matrix inversion.
    The results are reported in `error-estimates.csv` within the `"Output"` folder.
  - Added Adaptive Mesh Refinement (AMR), specified in the `config["Model"]["Refinement"]`,
    for all problem types aside from transient. To enable AMR, a user must specify
    `"MaxIts"`, while all other options have reasonable defaults. Nonconformal (all mesh
    types) and conformal (simplex meshes) refinement are supported.
  - Added support for non-axis-aligned lumped ports and current sources. Key words `"X"`,
    `"Y"`, `"Z"` and `"R"`, with optional prefix `"+"` or `"-"` still work, but now
    directions can be specified as vectors with 3 components. Users will be warned, and
    ultimately errored, if the specified directions do not agree with axis directions
    discovered from the geometry.
  - Added output of lumped port voltage and current for eigenmode simulations.
  - Added dimensionalized output for energies, voltages, currents, and field values based on
    a choice of the characteristic magnetic field strength used for nondimensionalization.
  - Added output of electric and magnetic field energies and participation ratios in regions
    of the domain, specified with `config["Domains"]["Postprocessing"]["Energy"]` and
    written to `domain-E.csv`. This replaces
    `config["Domains"]["Postprocessing"]["Dielectric"]` and `domain-Q.csv`.
  - Added improved `Timer` and `BlockTimer` classes with more timing categories for
    reporting simulation runtime.
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
  - Added build dependencies on [libCEED](https://github.com/CEED/libCEED), including
    [LIBXSMM](https://github.com/libxsmm/libxsmm) and [MAGMA](https://icl.utk.edu/magma/)
    to support CPU- and GPU-based operator partial assembly.
  - Added unit test framework for all integrators based on
    [Catch2](https://github.com/catchorg/Catch2), which also includes some automated
    benchmarking capabilities for operator assembly and application.
  - Added improved OpenMP support in `palace` wrapper script and CI tests.
  - Added Apptainer/Singularity container build definition for *Palace*.
  - Fixed bugs related to thread-safety for OpenMP builds and parallel tetrahedral meshes in
    the upstream MFEM library.

## [0.11.2] - 2023-07-14

  - Fixed a regression bug affecting meshes which have domain elements which are not
    assigned material properties in the configuration file.
  - Changed layout and names of `palace/` source directory for better organization.
  - Added many updates to build system: Removed use of Git submodules to download
    dependencies relying instead directly on CMake's ExternalProject, patch GSLIB dependency
    for shared library builds, add CI tests with ARPACK-NG instead of SLEPc, update all
    dependency versions including MFEM to incorporate bug fixes and improvements. This
    affects the Spack package recipe, so a new recipe is distributed as part of *Palace* in
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
