# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Build ARPACK/PARPACK (from ARPACK-NG)
#

# Force build order
set(ARPACK_DEPENDENCIES)

# We always build the 32-bit integer ARPACK interface and link with LP64 BLAS/LAPACK
# For PARPACK, this strategy is only not feasible when matrix sizes PER MPI PROCESS exceed
# 2B nonzeros, which is very unlikely
if(PALACE_WITH_64BIT_BLAS_INT)
  message(FATAL_ERROR "ARPACK has not been tested with INTERFACE64 and ILP64 BLAS/LAPACK")
endif()

set(ARPACK_OPTIONS ${PALACE_SUPERBUILD_DEFAULT_ARGS})
list(APPEND ARPACK_OPTIONS
  "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
  "-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
  "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
  "-DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}"
  "-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}"
  "-DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}"
  "-DMPI=ON"
  "-DICB=ON"
  "-DINTERFACE64=OFF"
  "-DTESTS=OFF"
)

# Configure BLAS/LAPACK
if(NOT "${BLAS_LAPACK_LIBRARIES}" STREQUAL "")
  list(APPEND ARPACK_OPTIONS
    "-DLAPACK_LIBRARIES=${BLAS_LAPACK_LIBRARIES}"
    "-DBLAS_LIBRARIES=${BLAS_LAPACK_LIBRARIES}"
  )
endif()

string(REPLACE ";" "; " ARPACK_OPTIONS_PRINT "${ARPACK_OPTIONS}")
message(STATUS "ARPACK_OPTIONS: ${ARPACK_OPTIONS_PRINT}")

# ARPACK-NG patches zdotc to a custom zzdotc, which unfortunately conflicts with a similar
# patch from the reference ScaLAPACK, so we patch the patch
set(ARPACK_PATCH_FILES
  "${CMAKE_SOURCE_DIR}/extern/patch/arpack-ng/patch_build.diff"
  "${CMAKE_SOURCE_DIR}/extern/patch/arpack-ng/patch_zdotc.diff"
  "${CMAKE_SOURCE_DIR}/extern/patch/arpack-ng/patch_pzneupd.diff"
)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  list(APPEND ARPACK_PATCH_FILES
    "${CMAKE_SOURCE_DIR}/extern/patch/arpack-ng/patch_second.diff"
  )
endif()

include(ExternalProject)
ExternalProject_Add(arpack-ng
  DEPENDS           ${ARPACK_DEPENDENCIES}
  GIT_REPOSITORY    ${EXTERN_ARPACK_URL}
  GIT_TAG           ${EXTERN_ARPACK_GIT_TAG}
  SOURCE_DIR        ${CMAKE_BINARY_DIR}/extern/arpack-ng
  BINARY_DIR        ${CMAKE_BINARY_DIR}/extern/arpack-ng-build
  INSTALL_DIR       ${CMAKE_INSTALL_PREFIX}
  PREFIX            ${CMAKE_BINARY_DIR}/extern/arpack-ng-cmake
  UPDATE_COMMAND    ""
  PATCH_COMMAND     git apply "${ARPACK_PATCH_FILES}"
  CONFIGURE_COMMAND ${CMAKE_COMMAND} <SOURCE_DIR> "${ARPACK_OPTIONS}"
  TEST_COMMAND      ""
)
