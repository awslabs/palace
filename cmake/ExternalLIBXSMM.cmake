# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Build LIBXSMM (for libCEED)
#

# Force build order
set(LIBXSMM_DEPENDENCIES)

set(LIBXSMM_OPTIONS
  "PREFIX=${CMAKE_INSTALL_PREFIX}"
  "OUTDIR=${CMAKE_INSTALL_PREFIX}/lib"
  "DIRSTATE=."
  "PPKGDIR=lib/pkgconfig"
  "PMODDIR=lib/pkgconfig"
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
  BUILD_COMMAND     ""
  INSTALL_COMMAND
    ${CMAKE_MAKE_PROGRAM} ${LIBXSMM_OPTIONS} install-minimal &&
    ${CMAKE_COMMAND} -E remove -f "${CMAKE_INSTALL_PREFIX}/lib/.make" &&
    ${CMAKE_COMMAND} -E remove -f "${CMAKE_INSTALL_PREFIX}/lib/libxsmm.pc" &&
    ${CMAKE_COMMAND} -E remove -f "${CMAKE_INSTALL_PREFIX}/lib/libxsmmext.pc" &&
    ${CMAKE_COMMAND} -E remove -f "${CMAKE_INSTALL_PREFIX}/lib/libxsmmnoblas.pc" &&
    ${CMAKE_COMMAND} -E remove -f "${CMAKE_INSTALL_PREFIX}/lib/libxsmm-shared.pc" &&
    ${CMAKE_COMMAND} -E remove -f "${CMAKE_INSTALL_PREFIX}/lib/libxsmmext-shared.pc" &&
    ${CMAKE_COMMAND} -E remove -f "${CMAKE_INSTALL_PREFIX}/lib/libxsmmnoblas-shared.pc" &&
    ${CMAKE_COMMAND} -E remove -f "${CMAKE_INSTALL_PREFIX}/lib/libxsmm.env"
  TEST_COMMAND      ""
)
