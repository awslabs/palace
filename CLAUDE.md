# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is Palace

Palace (**PA**rallel **LA**rge-scale **C**omputational **E**lectromagnetics) is a C++17 parallel finite element solver for 3D electromagnetic simulations. It uses MFEM for finite element discretization and libCEED for efficient operator evaluation. All code lives in the `palace` namespace.

## Build Commands

Palace uses a CMake superbuild that downloads and builds most dependencies automatically.

```bash
# Configure and build (from repo root)
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j

# The binary is installed to build/bin/palace-<arch>.bin
# A wrapper script is at build/bin/palace
```

Key CMake options: `-DPALACE_WITH_CUDA=ON`, `-DPALACE_WITH_HIP=ON`, `-DPALACE_WITH_SUPERLU=ON/OFF`, `-DPALACE_WITH_SLEPC=ON/OFF`, `-DPALACE_WITH_STRUMPACK=ON/OFF`, `-DPALACE_WITH_MUMPS=ON/OFF`.

## Running Tests

### Unit tests (C++, Catch2)
```bash
# Build unit tests (from build directory)
make palace-tests    # builds the unit-tests target via superbuild

# Run all unit tests
cd build && ctest --test-dir palace-build/test/unit

# Run a single test by name
cd build && ctest --test-dir palace-build/test/unit -R "test-name-pattern"

# Run only serial or parallel tests
cd build && ctest --test-dir palace-build/test/unit -R "serial-"
cd build && ctest --test-dir palace-build/test/unit -R "mpi-"
```

Every `TEST_CASE` must have one of these tags: `[Serial]`, `[Parallel]`, `[GPU]`. CI checks this. Tests are discovered automatically by Catch2; MPI tests run with 2 processes.

### Integration tests (Julia)
```bash
julia --project=test/examples -e 'using Pkg; Pkg.instantiate()'
julia --project=test/examples --color=yes test/examples/runtests.jl
```

## Code Formatting and Style

```bash
# Format all C++ and Julia source
./scripts/format-source

# Format only C++ or Julia
./scripts/format-source -cpp
./scripts/format-source -jl
```

- C++ uses `clang-format` (v19 in CI). Config is in `.clang-format`: Allman braces, 92-column limit, 2-space indent.
- Header guards use `#ifndef PALACE_<MODULE>_<FILE>_HPP` pattern.
- Every source file requires this license header:
  ```
  // Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
  // SPDX-License-Identifier: Apache-2.0
  ```
- CMake/Julia files use `#` comment version of the same header.

## CI Checks

The `style.yml` workflow runs on every PR:
1. **clang-format** check on all C++ source
2. **JSON Schema validation** of all example config files (`scripts/validate-config`)
3. **Test tag check** — ensures every `TEST_CASE` has `[Serial]`, `[Parallel]`, or `[GPU]`
4. **License header check** on all source files

## Architecture

### Source Layout (`palace/`)

- **`drivers/`** — Problem-type solvers. `BaseSolver` is the abstract base; concrete implementations: `DrivenSolver` (frequency domain), `EigenSolver` (eigenmode), `ElectrostaticSolver`, `MagnetostaticSolver`, `TransientSolver`. The main solve loop uses `SolveEstimateMarkRefine` for adaptive mesh refinement.

- **`models/`** — Physics operators that assemble the FEM problem. `SpaceOperator` builds the core Maxwell operators (curl-curl, mass). `MaterialOperator` handles material properties. Other operators handle boundary conditions (lumped ports, wave ports, surface impedance, far-field) and postprocessing (domain/surface quantities, Stratton-Chu).

- **`fem/`** — Finite element infrastructure wrapping MFEM. `Mesh` wraps `mfem::ParMesh`. `BilinearForm`/`Integrator` provide libCEED-accelerated operator assembly. `fem/libceed/` contains the libCEED backend (basis, restriction, coefficient, operator). `fem/integ/` has specific integrator implementations. `fem/qfunctions/` has libCEED Q-functions (JIT-compiled at runtime).

- **`linalg/`** — Linear algebra solvers and preconditioners. `KspSolver` is the main iterative solver interface. Includes wrappers for HYPRE (AMG/AMS), PETSc/SLEPc (eigensolvers), SuperLU/STRUMPACK/MUMPS (direct solvers). `Vector` extends `mfem::Vector` with complex number support. `RAP` handles operator assembly (R^T A P pattern).

- **`utils/`** — Configuration parsing (`IoData`, `ConfigFile`), mesh I/O (`GeoData`, `MeshIO`), MPI communication helpers, unit system (`Units`), and diagnostics.

- **`main.cpp`** — Entry point. Initializes MPI, parses config, selects driver based on problem type, runs solve.

### Superbuild Structure

The top-level `CMakeLists.txt` is a superbuild that builds dependencies via `ExternalProject_Add` (in `cmake/External*.cmake` files), then builds Palace itself as a subproject. The inner project CMake is `palace/CMakeLists.txt`. Dependencies are fetched from Git (tags defined in `cmake/ExternalGitTags.cmake`).

### Configuration

Palace is driven by JSON configuration files. JSON schemas are in `scripts/schema/` and embedded into the binary at build time. The `IoData` class parses and validates configuration, then nondimensionalizes all input values.

### Test Data

Unit test data files live in `test/unit/data/`. Example simulation configs are in `examples/` (one directory per example: `coaxial`, `cpw`, `rings`, `spheres`, `cylinder`, `antenna`, `adapter`, `transmon`).
