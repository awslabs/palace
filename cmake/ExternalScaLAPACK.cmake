# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Build ScaLAPACK (from scivision, with CMake)
#

# Force build order
set(SCALAPACK_DEPENDENCIES)

set(SCALAPACK_OPTIONS ${PALACE_SUPERBUILD_DEFAULT_ARGS})
list(APPEND SCALAPACK_OPTIONS
  "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
  "-DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}"
  "-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}"
  "-DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}"
  "-DBUILD_SINGLE=ON"
  "-DBUILD_DOUBLE=ON"
  "-DBUILD_COMPLEX=ON"
  "-DBUILD_COMPLEX16=ON"
  "-DBUILD_TESTING=OFF"
)

# Configure LAPACK dependency
if(NOT "${BLAS_LAPACK_LIBRARIES}" STREQUAL "")
  list(APPEND SCALAPACK_OPTIONS
    "-DLAPACK_LIBRARIES=${BLAS_LAPACK_LIBRARIES}"
    "-DBLAS_LIBRARIES=${BLAS_LAPACK_LIBRARIES}"
  )
endif()

string(REPLACE ";" "; " SCALAPACK_OPTIONS_PRINT "${SCALAPACK_OPTIONS}")
message(STATUS "SCALAPACK_OPTIONS: ${SCALAPACK_OPTIONS_PRINT}")

# Fix build
set(SCALAPACK_PATCH_FILES
  "${CMAKE_CURRENT_SOURCE_DIR}/patch/scalapack/patch_build.diff"
  "${CMAKE_CURRENT_SOURCE_DIR}/patch/scalapack/patch_version.diff"
)

include(ExternalProject)
ExternalProject_Add(scalapack
  DEPENDS           ${SCALAPACK_DEPENDENCIES}
  GIT_REPOSITORY    ${EXTERN_SCALAPACK_URL}
  GIT_TAG           ${EXTERN_SCALAPACK_GIT_TAG}
  SOURCE_DIR        ${CMAKE_CURRENT_BINARY_DIR}/scalapack
  BINARY_DIR        ${CMAKE_CURRENT_BINARY_DIR}/scalapack-build
  INSTALL_DIR       ${CMAKE_INSTALL_PREFIX}
  PREFIX            ${CMAKE_CURRENT_BINARY_DIR}/scalapack-cmake
  UPDATE_COMMAND    ""
  PATCH_COMMAND     git apply "${SCALAPACK_PATCH_FILES}"
  CONFIGURE_COMMAND cmake <SOURCE_DIR> "${SCALAPACK_OPTIONS}"
  TEST_COMMAND      ""
)

include(GNUInstallDirs)
if(BUILD_SHARED_LIBS)
  set(_SCALAPACK_LIB_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
else()
  set(_SCALAPACK_LIB_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()
set(_SCALAPACK_LIBRARIES ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/libblacs${_SCALAPACK_LIB_SUFFIX})
set(_SCALAPACK_LIBRARIES ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/libscalapack${_SCALAPACK_LIB_SUFFIX}$<SEMICOLON>${_SCALAPACK_LIBRARIES})
set(SCALAPACK_LIBRARIES ${_SCALAPACK_LIBRARIES} CACHE STRING "List of library files for ScaLAPACK")
set(SCALAPACK_INCLUDE_DIRS ${CMAKE_INSTALL_PREFIX}/include CACHE STRING "Path to ScaLAPACK include directories")