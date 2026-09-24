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
  # LIBXSMM 2.x writes its pkg-config files directly to OUTDIR/pkgconfig with includedir
  # defaulting to include/libxsmm; keep the flat include/ layout installed below
  "PINCDIR=include"
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

# Build directly into the installation directory instead of using the LIBXSMM install step,
# which before 2.x left build-tree paths in the installed shared libraries
# (https://github.com/libxsmm/libxsmm/issues/883). The 2.x install step fixes that but also
# moves the headers to include/libxsmm.
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
    ${CMAKE_COMMAND} -E rm -f ${CMAKE_INSTALL_PREFIX}/lib/.make
      ${CMAKE_INSTALL_PREFIX}/lib/pkgconfig/.make
  TEST_COMMAND      ""
)
