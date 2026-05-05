# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Generate and embed the Palace JSON schema.
#
# Pipeline:
#   1. The `palace_emit_schema` helper (palace/schema/CMakeLists.txt) is built
#      from the annotated palace::schema::* types. It prints a single draft
#      2020-12 schema with all sections under $defs to stdout.
#   2. A custom command runs that helper, redirecting stdout to
#      build/generated/schema/config-schema.json.
#   3. The existing cmake/embed_schema.cmake helper script reads every *.json
#      under that directory and emits build/generated/embedded_schema.hpp
#      exposing `palace::schema::GetSchemaMap()`.
#   4. libpalace depends on the `embedded_schemas` custom target so the header
#      is regenerated whenever the emit binary's output changes.
#
# Unlike the previous layer-over-the-source-tree pipeline, the schema is a
# build artifact — the source tree no longer carries scripts/schema/config/*.json.

set(SCHEMA_OUTPUT_DIR    "${CMAKE_BINARY_DIR}/generated")
set(SCHEMA_JSON_DIR      "${SCHEMA_OUTPUT_DIR}/schema")
set(SCHEMA_JSON          "${SCHEMA_JSON_DIR}/config-schema.json")
set(SCHEMA_HEADER        "${SCHEMA_OUTPUT_DIR}/embedded_schema.hpp")
set(SCHEMA_COMMITTED_DIR "${CMAKE_SOURCE_DIR}/../scripts/schema")
set(SCHEMA_COMMITTED     "${SCHEMA_COMMITTED_DIR}/config-schema.json")

file(MAKE_DIRECTORY ${SCHEMA_JSON_DIR})

# Emit the schema JSON by running palace_emit_schema. The command re-runs
# whenever the emit binary is rebuilt (any annotated type header change, any
# reflect-cpp change) — CMake tracks the executable as a dependency.
add_custom_command(
  OUTPUT ${SCHEMA_JSON}
  COMMAND palace_emit_schema > ${SCHEMA_JSON}
  DEPENDS palace_emit_schema
  COMMENT "Emitting Palace config schema via palace_emit_schema"
  VERBATIM
)

# Sync the emitted schema into the source tree at scripts/schema so it can be
# committed alongside the annotated types. `copy_if_different` leaves the
# file's mtime alone when the contents match, which keeps no-op rebuilds
# from dirtying git status. Downstream tooling (docs generator, Julia
# validator) consumes this stable path rather than the build directory.
add_custom_command(
  OUTPUT ${SCHEMA_COMMITTED}
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SCHEMA_JSON} ${SCHEMA_COMMITTED}
  DEPENDS ${SCHEMA_JSON}
  COMMENT "Syncing config-schema.json to scripts/schema/"
  VERBATIM
)

add_custom_target(palace_schema_file ALL DEPENDS ${SCHEMA_COMMITTED})

# Run the existing embed helper against the single generated JSON directory.
# Embedding reads from the build-tree copy — the source-tree copy exists for
# humans and external tooling, the build-tree copy is what gets inlined.
add_custom_command(
  OUTPUT ${SCHEMA_HEADER}
  COMMAND ${CMAKE_COMMAND}
    -DSCHEMA_DIR=${SCHEMA_JSON_DIR}
    -DSCHEMA_HEADER=${SCHEMA_HEADER}
    -P ${CMAKE_CURRENT_LIST_DIR}/embed_schema.cmake
  DEPENDS ${SCHEMA_JSON}
  COMMENT "Regenerating embedded JSON schema header"
)

add_custom_target(embedded_schemas DEPENDS ${SCHEMA_HEADER})

# Make libpalace depend on the embedded header so it exists before compile.
add_dependencies(${LIB_TARGET_NAME} embedded_schemas)
