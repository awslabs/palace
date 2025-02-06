# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Build MUMPS (from scivision, with CMake)
#

# Force build order
set(MUMPS_DEPENDENCIES scalapack parmetis metis)

set(MUMPS_OPTIONS ${PALACE_SUPERBUILD_DEFAULT_ARGS})
list(APPEND MUMPS_OPTIONS
  "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
  "-DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}"
  "-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}"
  "-DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}"
  "-DMUMPS_parallel=ON"
  "-Dopenmp=${PALACE_WITH_OPENMP}"
  "-Dintsize64=OFF"
  "-DBUILD_SINGLE=OFF"
  "-DBUILD_DOUBLE=ON"
  "-DBUILD_COMPLEX=OFF"
  "-DBUILD_COMPLEX16=OFF"
  "-DMUMPS_BUILD_TESTING=OFF"
  "-Dmetis=ON"
  "-Dparmetis=ON"
  "-Dscotch=OFF"
  "-Dscalapack=ON"
  "-DPARMETIS_LIBRARY=${PARMETIS_LIBRARIES}"
  "-DMETIS_LIBRARY=${METIS_LIBRARIES}"
  "-DMETIS_INCLUDE_DIR=${CMAKE_INSTALL_PREFIX}/include"
  "-DSCALAPACK_LIBRARIES=${SCALAPACK_LIBRARIES}"
  "-DSCALAPACK_INCLUDE_DIRS=${CMAKE_INSTALL_PREFIX}/include"
)

# Configure LAPACK dependency
if(NOT "${BLAS_LAPACK_LIBRARIES}" STREQUAL "")
  list(APPEND MUMPS_OPTIONS
    "-DLAPACK_LIBRARIES=${BLAS_LAPACK_LIBRARIES}"
    "-DLAPACK_INCLUDE_DIRS=${BLAS_LAPACK_INCLUDE_DIRS}"
  )
endif()

string(REPLACE ";" "; " MUMPS_OPTIONS_PRINT "${MUMPS_OPTIONS}")
message(STATUS "MUMPS_OPTIONS: ${MUMPS_OPTIONS_PRINT}")

# Fix FindLAPACK and FindScaLAPACK in configuration
set(MUMPS_PATCH_FILES
  "${CMAKE_SOURCE_DIR}/extern/patch/mumps/patch_build.diff"
)

include(ExternalProject)
ExternalProject_Add(mumps
  DEPENDS           ${MUMPS_DEPENDENCIES}
  GIT_REPOSITORY    ${EXTERN_MUMPS_URL}
  GIT_TAG           ${EXTERN_MUMPS_GIT_TAG}
  SOURCE_DIR        ${CMAKE_BINARY_DIR}/extern/mumps
  BINARY_DIR        ${CMAKE_BINARY_DIR}/extern/mumps-build
  INSTALL_DIR       ${CMAKE_INSTALL_PREFIX}
  PREFIX            ${CMAKE_BINARY_DIR}/extern/mumps-cmake
  UPDATE_COMMAND    ""
  PATCH_COMMAND     git apply "${MUMPS_PATCH_FILES}"
  CONFIGURE_COMMAND ${CMAKE_COMMAND} <SOURCE_DIR> "${MUMPS_OPTIONS}"
  TEST_COMMAND      ""
)
