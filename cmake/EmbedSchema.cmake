# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Embed JSON schema files as C++ string literals
# This allows schema validation without runtime file dependencies

set(SCHEMA_DIR "${CMAKE_SOURCE_DIR}/../scripts/schema")
set(SCHEMA_OUTPUT_DIR "${CMAKE_BINARY_DIR}/generated")
set(SCHEMA_HEADER "${SCHEMA_OUTPUT_DIR}/embedded_schema.hpp")

if(NOT EXISTS "${SCHEMA_DIR}")
  message(FATAL_ERROR "Schema directory not found: ${SCHEMA_DIR}")
endif()

file(MAKE_DIRECTORY ${SCHEMA_OUTPUT_DIR})

# Collect all schema files
file(GLOB_RECURSE SCHEMA_FILES "${SCHEMA_DIR}/*.json")

# Regenerate the header during build (not configure) when schema files change
add_custom_command(
  OUTPUT ${SCHEMA_HEADER}
  COMMAND ${CMAKE_COMMAND}
    -DSCHEMA_DIR=${SCHEMA_DIR}
    -DSCHEMA_HEADER=${SCHEMA_HEADER}
    -P ${CMAKE_CURRENT_LIST_DIR}/embed_schema.cmake
  DEPENDS ${SCHEMA_FILES}
  COMMENT "Regenerating embedded JSON schemas"
)
add_custom_target(embedded_schemas DEPENDS ${SCHEMA_HEADER})

# Make the library target depend on the generated header
add_dependencies(${LIB_TARGET_NAME} embedded_schemas)
