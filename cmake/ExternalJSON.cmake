# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Configure JSON library from nlohamnn/json (header-only)
#

set(JSON_OPTIONS ${PALACE_SUPERBUILD_DEFAULT_ARGS})
list(APPEND JSON_OPTIONS
  "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
  "-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
  "-DJSON_Install=ON"
  "-DJSON_BuildTests=OFF"
)

string(REPLACE ";" "; " JSON_OPTIONS_PRINT "${JSON_OPTIONS}")
message(STATUS "JSON_OPTIONS: ${JSON_OPTIONS_PRINT}")

include(ExternalProject)
ExternalProject_Add(json
  URL               ${EXTERN_JSON_URL}
  SOURCE_DIR        ${CMAKE_BINARY_DIR}/extern/json
  BINARY_DIR        ${CMAKE_BINARY_DIR}/extern/json-build
  INSTALL_DIR       ${CMAKE_INSTALL_PREFIX}
  PREFIX            ${CMAKE_BINARY_DIR}/extern/json-cmake
  UPDATE_COMMAND    ""
  CONFIGURE_COMMAND ${CMAKE_COMMAND} <SOURCE_DIR> "${JSON_OPTIONS}"
  TEST_COMMAND      ""
)
