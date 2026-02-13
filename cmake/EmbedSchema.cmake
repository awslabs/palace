# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Embed JSON schema files as C++ string literals at build time.
# Uses add_custom_command so the header regenerates when schema files change,
# without triggering a full CMake reconfigure.

set(SCHEMA_DIR "${CMAKE_SOURCE_DIR}/../scripts/schema")
set(SCHEMA_OUTPUT_DIR "${CMAKE_BINARY_DIR}/generated")
set(SCHEMA_HEADER "${SCHEMA_OUTPUT_DIR}/embedded_schema.hpp")
set(SCHEMA_GENERATOR "${CMAKE_SOURCE_DIR}/../cmake/GenerateSchemaHeader.cmake")

if(NOT EXISTS "${SCHEMA_DIR}")
  message(FATAL_ERROR "Schema directory not found: ${SCHEMA_DIR}")
endif()

file(MAKE_DIRECTORY ${SCHEMA_OUTPUT_DIR})

# Collect all schema files for dependency tracking
file(GLOB_RECURSE SCHEMA_FILES "${SCHEMA_DIR}/*.json")

# Generate the header at build time via cmake -P, re-running only when schema
# files change. The DEPENDS list ensures make/ninja tracks modifications.
add_custom_command(
  OUTPUT "${SCHEMA_HEADER}"
  COMMAND ${CMAKE_COMMAND}
    -DSCHEMA_DIR=${SCHEMA_DIR}
    -DSCHEMA_OUTPUT_DIR=${SCHEMA_OUTPUT_DIR}
    -P ${SCHEMA_GENERATOR}
  DEPENDS ${SCHEMA_FILES} ${SCHEMA_GENERATOR}
  COMMENT "Generating embedded schema header from JSON schema files"
  VERBATIM
)

# Custom target so other targets can depend on the generated header
add_custom_target(generate_schema_header DEPENDS "${SCHEMA_HEADER}")
