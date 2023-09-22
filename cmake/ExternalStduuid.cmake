# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Configure UUID library from mariusbancila/stduuid (header-only)
#

set(STDUUID_OPTIONS ${PALACE_SUPERBUILD_DEFAULT_ARGS})
list(APPEND STDUUID_OPTIONS
  "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
  "-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
  "-DUUID_ENABLE_INSTALL=ON"
  "-DUUID_BUILD_TESTS=OFF"
  "-DUUID_SYSTEM_GENERATOR=ON"
)

string(REPLACE ";" "; " STDUUID_OPTIONS_PRINT "${STDUUID_OPTIONS}")
message(STATUS "STDUUID_OPTIONS: ${STDUUID_OPTIONS_PRINT}")

include(ExternalProject)
ExternalProject_Add(stduuid
  URL               ${EXTERN_STDUUID_URL}
  SOURCE_DIR        ${CMAKE_BINARY_DIR}/extern/stduuid
  BINARY_DIR        ${CMAKE_BINARY_DIR}/extern/stduuid-build
  INSTALL_DIR       ${CMAKE_INSTALL_PREFIX}
  PREFIX            ${CMAKE_BINARY_DIR}/extern/stduuid-cmake
  UPDATE_COMMAND    ""
  CONFIGURE_COMMAND ${CMAKE_COMMAND} <SOURCE_DIR> "${STDUUID_OPTIONS}"
  TEST_COMMAND      ""
)
