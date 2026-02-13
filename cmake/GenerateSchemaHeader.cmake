# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Build-time script invoked by cmake -P to generate embedded_schema.hpp.
# Expects -DSCHEMA_DIR=... -DSCHEMA_OUTPUT_DIR=... on the command line.

if(NOT DEFINED SCHEMA_DIR OR NOT DEFINED SCHEMA_OUTPUT_DIR)
  message(FATAL_ERROR "SCHEMA_DIR and SCHEMA_OUTPUT_DIR must be defined")
endif()

set(SCHEMA_HEADER "${SCHEMA_OUTPUT_DIR}/embedded_schema.hpp")

file(GLOB_RECURSE SCHEMA_FILES "${SCHEMA_DIR}/*.json")

set(SCHEMA_HEADER_CONTENT
"// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Auto-generated file - do not edit
// Generated from scripts/schema/*.json
// Provides embedded JSON schema strings for runtime config validation without file dependencies.

#ifndef PALACE_EMBEDDED_SCHEMA_HPP
#define PALACE_EMBEDDED_SCHEMA_HPP

#include <string>
#include <unordered_map>

namespace palace::schema
{

")

foreach(SCHEMA_FILE ${SCHEMA_FILES})
  file(RELATIVE_PATH REL_PATH "${SCHEMA_DIR}" "${SCHEMA_FILE}")

  string(REPLACE "/" "_" VAR_NAME "${REL_PATH}")
  string(REPLACE "." "_" VAR_NAME "${VAR_NAME}")
  string(REPLACE "-" "_" VAR_NAME "${VAR_NAME}")

  file(READ "${SCHEMA_FILE}" FILE_CONTENT)

  # Escape for C++ raw string (avoid )json" sequence)
  string(REPLACE ")json\"" ")json_\"" FILE_CONTENT "${FILE_CONTENT}")

  string(APPEND SCHEMA_HEADER_CONTENT
"inline const char* ${VAR_NAME} = R\"json(
${FILE_CONTENT})json\";

")
endforeach()

string(APPEND SCHEMA_HEADER_CONTENT
"// Look up embedded schema strings by their relative path (e.g. \"config-schema.json\").
inline const std::unordered_map<std::string, const char*>& GetSchemaMap()
{
  static const std::unordered_map<std::string, const char*> schema_map = {
")

foreach(SCHEMA_FILE ${SCHEMA_FILES})
  file(RELATIVE_PATH REL_PATH "${SCHEMA_DIR}" "${SCHEMA_FILE}")
  string(REPLACE "/" "_" VAR_NAME "${REL_PATH}")
  string(REPLACE "." "_" VAR_NAME "${VAR_NAME}")
  string(REPLACE "-" "_" VAR_NAME "${VAR_NAME}")

  string(APPEND SCHEMA_HEADER_CONTENT
"    {\"${REL_PATH}\", ${VAR_NAME}},
")
endforeach()

string(APPEND SCHEMA_HEADER_CONTENT
"  };
  return schema_map;
}

}  // namespace palace::schema

#endif  // PALACE_EMBEDDED_SCHEMA_HPP
")

file(WRITE "${SCHEMA_HEADER}" "${SCHEMA_HEADER_CONTENT}")

message(STATUS "Generated embedded schema header: ${SCHEMA_HEADER}")
