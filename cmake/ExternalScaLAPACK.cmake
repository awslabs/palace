# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Build ScaLAPACK
#

# Force build order
set(SCALAPACK_DEPENDENCIES)

# Silence compiler error
include(CheckCCompilerFlag)
set(SCALAPACK_CFLAGS "${CMAKE_C_FLAGS}")
check_c_compiler_flag(-Wno-implicit-function-declaration SUPPORTS_NOIMPLICITFUNC_WARNING)
if(SUPPORTS_NOIMPLICITFUNC_WARNING)
  set(SCALAPACK_CFLAGS "${SCALAPACK_CFLAGS} -Wno-implicit-function-declaration")
endif()

set(SCALAPACK_OPTIONS ${PALACE_SUPERBUILD_DEFAULT_ARGS})
list(APPEND SCALAPACK_OPTIONS
  "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
  "-DCMAKE_C_FLAGS=${SCALAPACK_CFLAGS}"
  "-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}"
  "-DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}"
  "-DSCALAPACK_BUILD_TESTS=OFF"
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

# With MPICH and the MPI implementations derived from it (MPICH_VERSION is defined), receive
# the point-to-point messages of BLACS with derived datatypes packed and unpack them
# locally. BLACS sends these messages packed and receives them with the datatype of the
# message, and MPICH (ch4) fails to unpack a packed message received in several shared
# memory cells into a noncontiguous derived datatype when a cell boundary splits one of
# its basic elements, as it does for double precision complex elements on macOS arm64
# (https://github.com/pmodels/mpich/issues/8169, fixed on the main branch of MPICH by
# https://github.com/pmodels/mpich/pull/8172, which no release up to 5.0.2 and 5.1.0b1
# has). This corrupts or aborts complex-valued ScaLAPACK factorizations (STRUMPACK and
# MUMPS).
# The patch changes nothing with other MPI implementations. Only the point-to-point
# receives need it: broadcasts and combine operations never receive a packed message with
# a derived datatype.
set(SCALAPACK_PATCH_FILES
  "${CMAKE_SOURCE_DIR}/extern/patch/scalapack/patch_blacs_recv.diff"
)

include(ExternalProject)
ExternalProject_Add(scalapack
  DEPENDS           ${SCALAPACK_DEPENDENCIES}
  GIT_REPOSITORY    ${EXTERN_SCALAPACK_URL}
  GIT_TAG           ${EXTERN_SCALAPACK_GIT_TAG}
  SOURCE_DIR        ${CMAKE_BINARY_DIR}/extern/scalapack
  BINARY_DIR        ${CMAKE_BINARY_DIR}/extern/scalapack-build
  INSTALL_DIR       ${CMAKE_INSTALL_PREFIX}
  PREFIX            ${CMAKE_BINARY_DIR}/extern/scalapack-cmake
  UPDATE_COMMAND    ""
  PATCH_COMMAND     git apply "${SCALAPACK_PATCH_FILES}"
  CONFIGURE_COMMAND ${CMAKE_COMMAND} <SOURCE_DIR> "${SCALAPACK_OPTIONS}"
  TEST_COMMAND      ""
)

include(GNUInstallDirs)
# Save variables to cache
if(BUILD_SHARED_LIBS)
  set(_SCALAPACK_LIB_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
else()
  set(_SCALAPACK_LIB_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()
set(SCALAPACK_LIBRARIES ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/libscalapack${_SCALAPACK_LIB_SUFFIX}
  CACHE STRING "List of library files for ScaLAPACK"
)
