# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Configure fmt library
#

set(FMT_OPTIONS ${PALACE_SUPERBUILD_DEFAULT_ARGS})
list(APPEND FMT_OPTIONS
  "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
  "-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
  "-DFMT_INSTALL=ON"
  "-DFMT_DOC=OFF"
  "-DFMT_TEST=OFF"
)

string(REPLACE ";" "; " FMT_OPTIONS_PRINT "${FMT_OPTIONS}")
message(STATUS "FMT_OPTIONS: ${FMT_OPTIONS_PRINT}")

include(ExternalProject)
ExternalProject_Add(fmt
  URL               ${EXTERN_FMT_URL}
  SOURCE_DIR        ${CMAKE_CURRENT_BINARY_DIR}/fmt
  BINARY_DIR        ${CMAKE_CURRENT_BINARY_DIR}/fmt-build
  INSTALL_DIR       ${CMAKE_INSTALL_PREFIX}
  PREFIX            ${CMAKE_CURRENT_BINARY_DIR}/fmt-cmake
  UPDATE_COMMAND    ""
  CONFIGURE_COMMAND cmake <SOURCE_DIR> "${FMT_OPTIONS}"
  TEST_COMMAND      ""
)
