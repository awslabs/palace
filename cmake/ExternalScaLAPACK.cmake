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
