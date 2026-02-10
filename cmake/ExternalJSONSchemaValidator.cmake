# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Configure JSON Schema Validator library (depends on nlohmann/json)
#

set(JSON_SCHEMA_VALIDATOR_OPTIONS ${PALACE_SUPERBUILD_DEFAULT_ARGS})
list(APPEND JSON_SCHEMA_VALIDATOR_OPTIONS
  "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
  "-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
  "-DJSON_VALIDATOR_INSTALL=ON"
  "-DJSON_VALIDATOR_BUILD_TESTS=OFF"
  "-DJSON_VALIDATOR_BUILD_EXAMPLES=OFF"
  "-DJSON_VALIDATOR_SHARED_LIBS=OFF"
  "-Dnlohmann_json_DIR=${CMAKE_INSTALL_PREFIX}/share/cmake/nlohmann_json"
)

string(REPLACE ";" "; " JSON_SCHEMA_VALIDATOR_OPTIONS_PRINT "${JSON_SCHEMA_VALIDATOR_OPTIONS}")
message(STATUS "JSON_SCHEMA_VALIDATOR_OPTIONS: ${JSON_SCHEMA_VALIDATOR_OPTIONS_PRINT}")

include(ExternalProject)
ExternalProject_Add(json_schema_validator
  DEPENDS           json
  URL               ${EXTERN_JSON_SCHEMA_VALIDATOR_URL}
  SOURCE_DIR        ${CMAKE_BINARY_DIR}/extern/json-schema-validator
  BINARY_DIR        ${CMAKE_BINARY_DIR}/extern/json-schema-validator-build
  INSTALL_DIR       ${CMAKE_INSTALL_PREFIX}
  PREFIX            ${CMAKE_BINARY_DIR}/extern/json-schema-validator-cmake
  UPDATE_COMMAND    ""
  CONFIGURE_COMMAND ${CMAKE_COMMAND} <SOURCE_DIR> "${JSON_SCHEMA_VALIDATOR_OPTIONS}"
  TEST_COMMAND      ""
)
