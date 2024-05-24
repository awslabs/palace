# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Build LIBXSMM (for libCEED)
#

# Force build order
set(LIBXSMM_DEPENDENCIES)

set(LIBXSMM_OPTIONS
  # "PREFIX=${CMAKE_INSTALL_PREFIX}"  # Don't use install step, see comment below
  "OUTDIR=${CMAKE_INSTALL_PREFIX}/lib"
  "DIRSTATE=."
  "CC=${CMAKE_C_COMPILER}"
  "CXX=${CMAKE_CXX_COMPILER}"
  "FC="
  "FORTRAN=0"
  "BLAS=0"  # For now, no BLAS linkage (like PyFR)
  "SYM=1"   # Always build with symbols
  "VERBOSE=1"
)

# Always build LIBXSMM as a shared library
list(APPEND LIBXSMM_OPTIONS
  "STATIC=0"
)

# Configure debugging
if(CMAKE_BUILD_TYPE MATCHES "Debug|debug|DEBUG")
  list(APPEND LIBXSMM_OPTIONS
    "DBG=1"
    "TRACE=1"
  )
endif()

# Fix libxsmmext library linkage on macOS
if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
  list(APPEND LIBXSMM_OPTIONS
    "LDFLAGS=-undefined dynamic_lookup"
  )
endif()

string(REPLACE ";" "; " LIBXSMM_OPTIONS_PRINT "${LIBXSMM_OPTIONS}")
message(STATUS "LIBXSMM_OPTIONS: ${LIBXSMM_OPTIONS_PRINT}")

# Don't use LIBXSMM install step, since it just copies shared libraries and doesn't modify
# the dependency locations directly (doesn't use RPATH). Just build directly into the
# installation directory instead. See https://github.com/libxsmm/libxsmm/issues/883.
set(LIBXSMM_INSTALL_HEADERS
  libxsmm.h
  libxsmm_config.h
  libxsmm_version.h
  libxsmm_cpuid.h
  libxsmm_fsspmdm.h
  libxsmm_generator.h
  libxsmm_intrinsics_x86.h
  libxsmm_macros.h
  libxsmm_math.h
  libxsmm_malloc.h
  libxsmm_memory.h
  libxsmm_sync.h
  libxsmm_typedefs.h
)
list(TRANSFORM LIBXSMM_INSTALL_HEADERS PREPEND <SOURCE_DIR>/include/)
set(LIBXSMM_INSTALL_PKGCONFIG
  libxsmm.pc
  libxsmmext.pc
  libxsmmnoblas.pc
  libxsmm-shared.pc
  libxsmmext-shared.pc
  libxsmmnoblas-shared.pc
  libxsmm.env
)
list(TRANSFORM LIBXSMM_INSTALL_PKGCONFIG PREPEND ${CMAKE_INSTALL_PREFIX}/lib/)

include(ExternalProject)
ExternalProject_Add(libxsmm
  DEPENDS           ${LIBXSMM_DEPENDENCIES}
  GIT_REPOSITORY    ${EXTERN_LIBXSMM_URL}
  GIT_TAG           ${EXTERN_LIBXSMM_GIT_TAG}
  SOURCE_DIR        ${CMAKE_BINARY_DIR}/extern/libxsmm
  INSTALL_DIR       ${CMAKE_INSTALL_PREFIX}
  PREFIX            ${CMAKE_BINARY_DIR}/extern/libxsmm-cmake
  BUILD_IN_SOURCE   TRUE
  UPDATE_COMMAND    ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND     ${CMAKE_MAKE_PROGRAM} ${LIBXSMM_OPTIONS}
  INSTALL_COMMAND
    ${CMAKE_COMMAND} -E echo "LIBXSMM installing interface..." &&
    ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/include &&
    ${CMAKE_COMMAND} -E copy ${LIBXSMM_INSTALL_HEADERS} ${CMAKE_INSTALL_PREFIX}/include &&
    ${CMAKE_COMMAND} -E echo "LIBXSMM installing pkg-config and module files..." &&
    ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/lib/pkgconfig &&
    ${CMAKE_COMMAND} -E copy ${LIBXSMM_INSTALL_PKGCONFIG} ${CMAKE_INSTALL_PREFIX}/lib/pkgconfig ||
    ${CMAKE_COMMAND} -E true &&  # No error if files don't exist
    ${CMAKE_COMMAND} -E rm -f ${LIBXSMM_INSTALL_PKGCONFIG} &&
    ${CMAKE_COMMAND} -E rm -f ${CMAKE_INSTALL_PREFIX}/lib/.make
  TEST_COMMAND      ""
)
