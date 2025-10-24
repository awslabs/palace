```@raw html
<!--- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. --->
<!--- SPDX-License-Identifier: Apache-2.0 --->
```

# Installation

*Palace* can be built and installed using the
[Spack HPC package manager](https://spack.io/), following the instructions in
the [Build using Spack](#Build-using-Spack) section. Containerized builds are possible with
Singularity/Apptainer, described in
[Build using Singularity/Apptainer](#Build-using-Singularity/Apptainer). Alternatively,
compiling from source using [CMake](https://cmake.org/download) is described in
[Build from source](#Build-from-source).

If you are a user, we recommend you install [*Palace* with
Spack](#Build-using-Spack). If you intend to develop *Palace*, [build from
source](#Build-from-source) instead.

## Build using Spack

*Palace* is a registered package in the built-in Spack package repository. To install the
solver, follow the
[instructions for setting up Spack on your system](https://spack.readthedocs.io/en/latest/getting_started.html)
and run:

```bash
spack install palace
```

This will install the default version of *Palace*. Spack supports installing
_variants_ of *Palace*. For instance, if you want to install *Palace* with CUDA,
MUMPS and SLEPc, call

```bash
spack install palace +mumps +slepc +cuda cuda_arch=90
```

where `cuda_arch` is determined by the [generation of your
GPU](https://developer.nvidia.com/cuda-gpus). More information about the
available configuration options and dependencies can be found using `spack info palace`. See the [official
tutorial](https://spack-tutorial.readthedocs.io/en/latest/tutorial_basics.html)
for an introduction.

## Build using Singularity/Apptainer

*Palace* can be built in a
[Singularity/Apptainer](https://apptainer.org/docs/user/main/introduction.html) container
for HPC environments
supporting the Singularity/Apptainer container system. To build the container using the
provided definition file in the
[singularity/](https://github.com/awslabs/palace/blob/main/singularity) directory, first
[set up Singularity/Apptainer on your system](https://github.com/apptainer/apptainer/blob/main/INSTALL.md)
and subsequently run:

```bash
singularity build palace.sif <SOURCE_DIR>/singularity/singularity.def
```

where the repository source code has been cloned to `<SOURCE_DIR>`. For more information
about Singularity/Apptainer, see the
[Quick Start](https://apptainer.org/docs/user/main/quick_start.html) guide in the
Singularity/Apptainer documentation.

## Build from source

A build from source requires the following prerequisites installed on your system:

  - [CMake](https://cmake.org/download) version 3.21 or later
  - C++17 compatible C++ compiler
  - C and Fortran (optional) compilers for dependency builds
  - MPI distribution
  - BLAS, LAPACK libraries (described below in [Math libraries](#Math-libraries))
  - [CUDA Toolkit](https://developer.nvidia.com/cuda-toolkit) or
    [ROCm](https://rocm.docs.amd.com/en/latest/) installation (optional, for GPU support
    only)

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
  - `CMAKE_CUDA_COMPILER`, `CMAKE_CUDA_FLAGS`, `CMAKE_CUDA_ARCHITECTURES`, and the
    corresponding `CMAKE_HIP_COMPILER`, `CMAKE_HIP_FLAGS`, and `CMAKE_HIP_ARCHITECTURES` for
    GPU-accelerated builds with CUDA or HIP.
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
  - `PALACE_WITH_OPENMP [OFF]` :  Use OpenMP for shared-memory parallelism
  - `PALACE_WITH_CUDA [OFF]` :  Use CUDA for NVIDIA GPU support
  - `PALACE_WITH_HIP [OFF]` :  Use HIP for AMD or NVIDIA GPU support
  - `PALACE_WITH_GPU_AWARE_MPI [OFF]` :  Option to set if MPI distribution is GPU aware
  - `PALACE_WITH_SUPERLU [ON]` :  Build with SuperLU_DIST sparse direct solver
  - `PALACE_WITH_STRUMPACK [OFF]` :  Build with STRUMPACK sparse direct solver
  - `PALACE_WITH_MUMPS [OFF]` :  Build with MUMPS sparse direct solver
  - `PALACE_WITH_SLEPC [ON]` :  Build with SLEPc eigenvalue solver
  - `PALACE_WITH_ARPACK [OFF]` :  Build with ARPACK eigenvalue solver
  - `PALACE_WITH_LIBXSMM [ON]` :  Build with LIBXSMM backend for libCEED
  - `PALACE_WITH_MAGMA [ON]` :  Build with MAGMA backend for libCEED
  - `PALACE_WITH_GSLIB [ON]` :  Build with GSLIB library for high-order field interpolation
  - `PALACE_WITH_SUNDIALS [ON]` : Build with SUNDIALS ODE solver library

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
    
      + If the `MKLROOT` environment variable is set, looks for an
        [Intel MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html)
        installation.
      + If the `AOCL_DIR` or `AOCLROOT` environment variables are set, looks for an
        [AMD Optimizing CPU Libraries (AOCL)](https://developer.amd.com/amd-aocl)
        installation of BLIS and libFLAME.
      + Otherwise, tries to locate an installation of [OpenBLAS](https://www.openblas.net/)
        which is permissively licensed and available from most package managers.

  - For `aarch64`/`arm64` systems:
    
      + If the `ARMPL_DIR` environment variable is set, looks for an
        [Arm Performance Libraries (PL)](https://www.arm.com/products/development-tools/server-and-hpc/allinea-studio/performance-libraries)
        installation.
      + Otherwise, tries to locate an installation of [OpenBLAS](https://www.openblas.net/).

If the installation path of OpenBLAS is non-standard or is not found by default, it can be
set using the `OPENBLAS_DIR` or `OPENBLASROOT` environment variables, or added to
`CMAKE_PREFIX_PATH` when calling CMake.

It is recommended in most cases to use serial BLAS and LAPACK builds (not multithreaded),
as the standard parallelization approach in *Palace* is to use pure MPI parallelism.

## Dependencies

*Palace* leverages the [MFEM finite element discretization library](http://mfem.org). It
always configures and builds its own installation of MFEM internally in order to support
the most up to date features and patches. Likewise, *Palace* will always build its own
installation of [libCEED](https://github.com/CEED/libCEED), and
[GSLIB](https://github.com/Nek5000/gslib), when `PALACE_WITH_GSLIB=ON`.

As part of the [Build from source](#Build-from-source), the CMake build will automatically
build and install a small number of third-party dependencies before building *Palace*. The
source code for these dependencies is downloaded during the build process:

  - [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) and
    [ParMETIS](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview)
  - [Hypre](https://github.com/hypre-space/hypre)
  - [SuperLU_DIST](https://github.com/xiaoyeli/superlu_dist) (optional, when
    `PALACE_WITH_SUPERLU=ON`)
  - [STRUMPACK](https://portal.nersc.gov/project/sparse/strumpack) (optional, when
    `PALACE_WITH_STRUMPACK=ON`), including
    [ButterflyPACK](https://github.com/liuyangzhuan/ButterflyPACK) and
    [zfp](https://github.com/LLNL/zfp) support
  - [MUMPS](http://mumps.enseeiht.fr/) (optional, when `PALACE_WITH_MUMPS=ON`)
  - [SLEPc](https://slepc.upv.es/) (optional, when `PALACE_WITH_SLEPC=ON`), including
    [PETSc](https://petsc.org/release/)
  - [ARPACK-NG](https://github.com/opencollab/arpack-ng) (optional, when
    `PALACE_WITH_ARPACK=ON`)
  - [LIBXSMM](https://github.com/libxsmm/libxsmm) (optional, when `PALACE_WITH_LIBXSMM=ON`)
  - [MAGMA](https://icl.utk.edu/magma/) (optional, when `PALACE_WITH_MAGMA=ON`)
  - [SUNDIALS](https://github.com/LLNL/sundials) (optional, when `PALACE_WITH_SUNDIALS=ON`)
  - [nlohmann/json](https://github.com/nlohmann/json)
  - [fmt](https://fmt.dev/latest)
  - [Eigen](https://eigen.tuxfamily.org)

For solving eigenvalue problems, at least one of SLEPc or ARPACK-NG must be specified.
Typically only one of the SuperLU_DIST, STRUMPACK, and MUMPS dependencies is required but
all can be built so the user can decide at runtime which solver to use.

For unit testing, *Palace* relies on the [Catch2
library](https://github.com/catchorg/Catch2), which is automatically downloaded
and built when building the `unit-tests` target. See the [Developer
Notes](developer/testing.md) for more information.
