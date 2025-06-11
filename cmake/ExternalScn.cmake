# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Configure scn library
#

set(SCN_OPTIONS ${PALACE_SUPERBUILD_DEFAULT_ARGS})
list(APPEND SCN_OPTIONS
  "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
  "-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
  "-DSCN_INSTALL=ON"
  "-DSCN_REGEX_BACKEND=std"
  "-DSCN_DISABLE_TOP_PROJECT=ON" # disable tests, docs, benchmarks
)

string(REPLACE ";" "; " SCN_OPTIONS_PRINT "${SCN_OPTIONS}")
message(STATUS "SCN_OPTIONS: ${SCN_OPTIONS_PRINT}")

include(ExternalProject)
ExternalProject_Add(scn
  URL               ${EXTERN_SCN_URL}
  SOURCE_DIR        ${CMAKE_BINARY_DIR}/extern/scn
  BINARY_DIR        ${CMAKE_BINARY_DIR}/extern/scn-build
  INSTALL_DIR       ${CMAKE_INSTALL_PREFIX}
  PREFIX            ${CMAKE_BINARY_DIR}/extern/scn-cmake
  UPDATE_COMMAND    ""
  CONFIGURE_COMMAND ${CMAKE_COMMAND} <SOURCE_DIR> "${SCN_OPTIONS}"
  TEST_COMMAND      ""
)
