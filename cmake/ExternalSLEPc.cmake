# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Build PETSc and SLEPc
#

# Force build order
set(PETSC_DEPENDENCIES)
set(SLEPC_DEPENDENCIES petsc)

# First build PETSc
set(PETSC_OPTIONS
  "COPTFLAGS=${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${BUILD_TYPE_UPPER}}"
  "CXXOPTFLAGS=${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${BUILD_TYPE_UPPER}}"
  "--prefix=${CMAKE_INSTALL_PREFIX}"
  "--with-cc=${CMAKE_C_COMPILER}"
  "--with-cxx=${CMAKE_CXX_COMPILER}"
  "--with-fc=0"
  "--with-scalar-type=complex"
  "--with-precision=double"
  "--with-clanguage=c"
  "--with-x=0"
  # "--with-petsc4py=1"
)
if(CMAKE_BUILD_TYPE MATCHES "Debug|debug|DEBUG")
  list(APPEND PETSC_OPTIONS "--with-debugging=1")
else()
  list(APPEND PETSC_OPTIONS "--with-debugging=0")
endif()
if(BUILD_SHARED_LIBS)
  list(APPEND PETSC_OPTIONS "--with-shared-libraries=1")
else()
  list(APPEND PETSC_OPTIONS "--with-shared-libraries=0")
endif()
if(PALACE_WITH_64BIT_INT)
  list(APPEND PETSC_OPTIONS "--with-64-bit-indices")
endif()
if(PALACE_WITH_64BIT_BLAS_INT)
  list(APPEND PETSC_OPTIONS "--known-64-bit-blas-indices=1")
  list(APPEND PETSC_OPTIONS "--with-64-bit-blas-indices")
else()
  list(APPEND PETSC_OPTIONS "--known-64-bit-blas-indices=0")
endif()
if(PALACE_WITH_OPENMP)
  list(APPEND PETSC_OPTIONS "--with-openmp")
endif()

# User might specify the MPI compiler wrappers directly, otherwise we need to supply MPI
# as found from the CMake module
if(NOT MPI_FOUND)
  message(FATAL_ERROR "MPI is not found when trying to build PETSc")
endif()
if(NOT CMAKE_CXX_COMPILER STREQUAL MPI_CXX_COMPILER)
  # For OpenMPI at least, when given a C++ compiler, PETSc needs the C++ MPI libraries for
  # its CxxMPICheck
  string(REPLACE ";" "," PETSC_MPI_LIBRARIES "${MPI_CXX_LIBRARIES}")
  string(REPLACE ";" "," PETSC_MPI_INCLUDE_DIRS "${MPI_CXX_INCLUDE_DIRS}")
  list(APPEND PETSC_OPTIONS
    "--with-mpi-lib=[${PETSC_MPI_LIBRARIES}]"
    "--with-mpi-include=[${PETSC_MPI_INCLUDE_DIRS}]"
  )
endif()

# Configure BLAS/LAPACK
if(NOT "${BLAS_LAPACK_LIBRARIES}" STREQUAL "")
  string(REPLACE "$<SEMICOLON>" "," PETSC_BLAS_LAPACK_LIBRARIES "${BLAS_LAPACK_LIBRARIES}")
  string(REPLACE "$<SEMICOLON>" "," PETSC_BLAS_LAPACK_INCLUDE_DIRS "${BLAS_LAPACK_INCLUDE_DIRS}")
  list(APPEND PETSC_OPTIONS
    "--with-blaslapack-lib=[${PETSC_BLAS_LAPACK_LIBRARIES}]"
    "--with-blaslapack-include=[${BLAS_LAPACK_INCLUDE_DIRS}]"
  )
endif()

# Configure GPU support
if(PALACE_WITH_CUDA)
  list(APPEND PETSC_OPTIONS "--with-cuda")
  if(NOT "${CMAKE_CUDA_ARCHITECTURES}" STREQUAL "")
    list(GET CMAKE_CUDA_ARCHITECTURES 0 PETSC_CUDA_ARCH)
    list(APPEND PETSC_OPTIONS
      "--with-cuda-arch=${PETSC_CUDA_ARCH}"
    )
  endif()
endif()
if(PALACE_WITH_HIP)
  list(APPEND PETSC_OPTIONS "--with-hip")
  if(NOT "${CMAKE_HIP_ARCHITECTURES}" STREQUAL "")
    list(GET CMAKE_HIP_ARCHITECTURES 0 PETSC_HIP_ARCH)
    list(APPEND PETSC_OPTIONS
      "--with-hip-arch=${PETSC_HIP_ARCH}"
    )
  endif()
endif()

string(REPLACE ";" "; " PETSC_OPTIONS_PRINT "${PETSC_OPTIONS}")
message(STATUS "PETSC_OPTIONS: ${PETSC_OPTIONS_PRINT}")

include(ExternalProject)
ExternalProject_Add(petsc
  DEPENDS             ${PETSC_DEPENDENCIES}
  GIT_REPOSITORY      ${EXTERN_PETSC_URL}
  GIT_TAG             ${EXTERN_PETSC_GIT_TAG}
  SOURCE_DIR          ${CMAKE_BINARY_DIR}/extern/petsc
  INSTALL_DIR         ${CMAKE_INSTALL_PREFIX}
  PREFIX              ${CMAKE_BINARY_DIR}/extern/petsc-cmake
  BUILD_IN_SOURCE     TRUE
  UPDATE_COMMAND      ""
  CONFIGURE_COMMAND   ./configure ${PETSC_OPTIONS}
  # TEST_COMMAND        ${CMAKE_MAKE_PROGRAM} check  # Use auto-detected PETSC_DIR/PETSC_ARCH
  TEST_COMMAND        ""
  TEST_BEFORE_INSTALL TRUE
)

# Configure SLEPc eigenvalue solver (most options come from PETSc)
set(SLEPC_OPTIONS
  "--prefix=${CMAKE_INSTALL_PREFIX}"
  "--with-feast=0"
  "--with-arpack=0"
  # "--with-slepc4py=1
)

string(REPLACE ";" "; " SLEPC_OPTIONS_PRINT "${SLEPC_OPTIONS}")
message(STATUS "SLEPC_OPTIONS: ${SLEPC_OPTIONS_PRINT}")

include(ExternalProject)
ExternalProject_Add(slepc
  DEPENDS             ${SLEPC_DEPENDENCIES}
  GIT_REPOSITORY      ${EXTERN_SLEPC_URL}
  GIT_TAG             ${EXTERN_SLEPC_GIT_TAG}
  SOURCE_DIR          ${CMAKE_BINARY_DIR}/extern/slepc
  INSTALL_DIR         ${CMAKE_INSTALL_PREFIX}
  PREFIX              ${CMAKE_BINARY_DIR}/extern/slepc-cmake
  BUILD_IN_SOURCE     TRUE
  UPDATE_COMMAND      ""
  CONFIGURE_COMMAND   SLEPC_DIR=<SOURCE_DIR> PETSC_DIR=${CMAKE_INSTALL_PREFIX} PETSC_ARCH= ./configure ${SLEPC_OPTIONS}
  BUILD_COMMAND       SLEPC_DIR=<SOURCE_DIR> PETSC_DIR=${CMAKE_INSTALL_PREFIX} PETSC_ARCH= ${CMAKE_MAKE_PROGRAM}
  INSTALL_COMMAND     SLEPC_DIR=<SOURCE_DIR> PETSC_DIR=${CMAKE_INSTALL_PREFIX} PETSC_ARCH= ${CMAKE_MAKE_PROGRAM} install
  # TEST_COMMAND        SLEPC_DIR=<SOURCE_DIR> PETSC_DIR=${CMAKE_INSTALL_PREFIX} PETSC_ARCH= ${CMAKE_MAKE_PROGRAM} check
  TEST_COMMAND        ""
  TEST_BEFORE_INSTALL TRUE
)
