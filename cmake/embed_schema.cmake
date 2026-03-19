# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Helper script for generating embedded_schema.hpp
# Called via: cmake -DSCHEMA_DIR=<dir> -DSCHEMA_HEADER=<output> -P embed_schema.cmake

if(NOT DEFINED SCHEMA_DIR OR NOT DEFINED SCHEMA_HEADER)
  message(FATAL_ERROR "SCHEMA_DIR and SCHEMA_HEADER must be defined")
endif()

if(NOT EXISTS "${SCHEMA_DIR}")
  message(FATAL_ERROR "Schema directory not found: ${SCHEMA_DIR}")
endif()

# Ensure output directory exists
get_filename_component(SCHEMA_OUTPUT_DIR "${SCHEMA_HEADER}" DIRECTORY)
file(MAKE_DIRECTORY "${SCHEMA_OUTPUT_DIR}")

# Collect all schema files
file(GLOB_RECURSE SCHEMA_FILES "${SCHEMA_DIR}/*.json")

# Generate the header file
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
  # Get relative path from schema dir
  file(RELATIVE_PATH REL_PATH "${SCHEMA_DIR}" "${SCHEMA_FILE}")

  # Create valid C++ identifier from path
  string(REPLACE "/" "_" VAR_NAME "${REL_PATH}")
  string(REPLACE "." "_" VAR_NAME "${VAR_NAME}")
  string(REPLACE "-" "_" VAR_NAME "${VAR_NAME}")

  # Read file content
  file(READ "${SCHEMA_FILE}" FILE_CONTENT)

  # Escape for C++ raw string (just need to avoid )json" sequence)
  string(REPLACE ")json\"" ")json_\"" FILE_CONTENT "${FILE_CONTENT}")

  string(APPEND SCHEMA_HEADER_CONTENT
"inline const char* ${VAR_NAME} = R\"json(
${FILE_CONTENT})json\";

")
endforeach()

# Add a map for looking up schemas by path
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
