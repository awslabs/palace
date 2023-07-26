```@raw html
<!--- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. --->
<!--- SPDX-License-Identifier: Apache-2.0 --->
```

# Installation

*Palace* can be built and installed using the [Spack HPC package manager]
(https://spack.io/), following the instructions in the [Build using Spack]
(#Build-using-Spack) section. Alternatively, compiling from source using [CMake]
(https://cmake.org/download) is described in [Build from source](#Build-from-source).

## Build using Spack

*Palace* is a registered package in the built-in Spack package repository. To install the
solver, follow the [instructions for setting up Spack on your system]
(https://spack.readthedocs.io/en/latest/getting_started.html) and run:

```bash
spack install palace
```

More information about about the available configuration options and dependencies can be
found using `spack info palace`.

## Build using Singularity / Apptainer

*Palace* can be built in a [Singularity / Apptainer]
(https://apptainer.org/docs/user/main/introduction.html) container often used in HPC environments.

To build the container, first [set up Singularity or Apptainer]
(https://github.com/apptainer/apptainer/blob/main/INSTALL.md) and subsequently run:

```bash
singularity build palace.sif singularity.def
```

For more information about Singularity / Apptainer,
see the [Quick Start guide](https://apptainer.org/docs/user/main/quick_start.html).

## Build from source

A build from source requires the following prerequisites installed on your system:

  - [CMake](https://cmake.org/download) version 3.18.1 or later
  - C++ compiler supporting C++17
  - C and (optionally) Fortran compilers for dependency builds
  - MPI distribution
  - BLAS, LAPACK libraries (described below in [Math libraries](#Math-libraries))

In addition, builds from source require the following system packages which are typically
already installed and are available from most package managers (`apt`, `dnf`, `brew`, etc.):

  - Python 3
  - [`pkg-config`](https://www.freedesktop.org/wiki/Software/pkg-config/)
  - [`libunwind`](https://www.nongnu.org/libunwind/) (optional)
  - [`zlib`](https://zlib.net/) (optional)

### Quick start

To start, clone the code using

```bash
git clone https://github.com/awslabs/palace.git
```

Then, a build using the default options can be performed by running the following from
within the directory where the repository was cloned:

```bash
mkdir build && cd build
cmake ..
make -j
```

This installs the binary executable in `build/bin/`.

### Configuration options

To configure a *Palace* build in `<BUILD_DIR>` using the source code in `<SOURCE_DIR>`,
run:

```bash
mkdir <BUILD_DIR> && cd <BUILD_DIR>
cmake [OPTIONS] <SOURCE_DIR>
```

Here, `[OPTIONS]` is a list of options passed to `cmake` of the form `-D<VARIABLE>=<VALUE>`.
The *Palace* build respects standard CMake variables, including:

  - `CMAKE_CXX_COMPILER`, `CMAKE_C_COMPILER`, and `CMAKE_Fortran_COMPILER` which define the
    desired compilers.
  - `CMAKE_CXX_FLAGS`, `CMAKE_C_FLAGS`, and `CMAKE_Fortran_FLAGS` which define the
    corresponding compiler flags.
  - `CMAKE_INSTALL_PREFIX` which specifies the path for installation (if none is provided,
    defaults to `<BUILD_DIR>`).
  - `CMAKE_BUILD_TYPE` which defines the build type such as `Release`, `Debug`,
    `RelWithDebInfo`, and `MinSizeRel` (`Release` if not otherwise specified).
  - `BUILD_SHARED_LIBS` which is a flag to create shared libraries for dependency library
    builds instead of static libraries (`OFF` by default).
  - `CMAKE_PREFIX_PATH` which lists directories specifying installation prefixes to be
    searched for dependencies.
  - `CMAKE_INSTALL_RPATH` and `CMAKE_INSTALL_RPATH_USE_LINK_PATH` which configure the rpath
    for installed library and executable targets.

Additional build options are (with default values in brackets):

  - `PALACE_WITH_64BIT_INT [OFF]` :  Build with 64-bit integer support
  - `PALACE_WITH_OPENMP [OFF]` :  Use OpenMP
  - `PALACE_WITH_GSLIB [ON]` :  Build with GSLIB library for high-order field interpolation
  - `PALACE_WITH_SUPERLU [ON]` :  Build with SuperLU_DIST sparse direct solver
  - `PALACE_WITH_STRUMPACK [OFF]` :  Build with STRUMPACK sparse direct solver
  - `PALACE_WITH_MUMPS [OFF]` :  Build with MUMPS sparse direct solver
  - `PALACE_WITH_SLEPC [ON]` :  Build with SLEPc eigenvalue solver
  - `PALACE_WITH_ARPACK [OFF]` :  Build with ARPACK eigenvalue solver

The build step is invoked by running (for example with 4 `make` threads)

```bash
make -j 4
```

or

```bash
cmake --build . -- -j 4
```

which installs the binary executable in `${CMAKE_INSTALL_PREFIX}/bin/`.

### Math libraries

During the configure step, the build system will try to detect system installations of BLAS
and LAPACK libraries depending on the system architecture according to the following
procedure:

  - For `x86_64` systems:
    
      + If the `MKLROOT` environment variable is set, looks for an [Intel MKL]
        (https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html)
        installation.
      + If the `AOCL_DIR` or `AOCLROOT` environment variables are set, looks for an
        [AMD Optimizing CPU Libraries (AOCL)](https://developer.amd.com/amd-aocl)
        installation of BLIS and libFLAME.
      + Otherwise, tries to locate an installation of [OpenBLAS](https://www.openblas.net/)
        which is permissively licensed and available from most package managers.

  - For `aarch64`/`arm64` systems:
    
      + If the `ARMPL_DIR` environment variable is set, looks for an
        [Arm Performance Libraries (PL)]
        (https://www.arm.com/products/development-tools/server-and-hpc/allinea-studio/performance-libraries)
        installation.
      + Otherwise, tries to locate an installation of [OpenBLAS](https://www.openblas.net/).

If the installation path of OpenBLAS is non-standard or is not found by default, it can be
set using the `OPENBLAS_DIR` or `OPENBLASROOT` environment variables, or added to
`CMAKE_PREFIX_PATH` when calling CMake.

It is recommended in most cases to use a serial BLAS and LAPACK builds (not multithreaded),
as the standard parallelization in approach in *Palace* is to use pure MPI parallelism.

## Dependencies

*Palace* leverages the [MFEM finite element discretization library](http://mfem.org). It
always configures and builds its own installation of MFEM internally in order to support
the most up to date features and patches. Likewise, Palace will always build its own
installation of [GSLIB](https://github.com/Nek5000/gslib), when `PALACE_WITH_GSLIB=ON`.

As part of the [Build from source](#Build-from-source), the CMake build will automatically
build and install a small number of third-party dependencies before building *Palace*. The
source code for these dependencies is downloaded using using [Git submodules]
(https://git-scm.com/book/en/v2/Git-Tools-Submodules). These libraries include:

  - [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) and [ParMETIS]
    (http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview)
  - [Hypre](https://github.com/hypre-space/hypre)
  - [SuperLU_DIST](https://github.com/xiaoyeli/superlu_dist) (optional, when
    `PALACE_WITH_SUPERLU=ON`)
  - [STRUMPACK](https://portal.nersc.gov/project/sparse/strumpack) (optional, when
    `PALACE_WITH_STRUMPACK=ON`), including [ButterflyPACK]
    (https://github.com/liuyangzhuan/ButterflyPACK) and [zfp](https://github.com/LLNL/zfp)
    support
  - [MUMPS](http://mumps.enseeiht.fr/) (optional, when `PALACE_WITH_MUMPS=ON`)
  - [SLEPc](https://slepc.upv.es/) (optional, when `PALACE_WITH_SLEPC=ON`), including
    [PETSc](https://petsc.org/release/)
  - [ARPACK-NG](https://github.com/opencollab/arpack-ng) (optional, when
    `PALACE_WITH_ARPACK=ON`)
  - [nlohmann/json](https://github.com/nlohmann/json)
  - [fmt](https://fmt.dev/latest)
  - [Eigen](https://eigen.tuxfamily.org)

For solving eigenvalue problems, at least one of SLEPc or ARPACK-NG must be specified.
Typically only one of the SuperLU_DIST, STRUMPACK, and MUMPS dependencies is required but
all can be built so the user can decide at runtime which solver to use.
