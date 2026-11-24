# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Embed the verdict data of the per-iteration table surface-participation.csv
# (palace/models/surfaceparticipationverdict.json) as a C++ string literal: the built-in
# thresholds are data, changed in that file only, and parsed by the same reader as the
# "Verdict" member of a Resolution.BarsFile file.

set(VERDICT_BARS_FILE "${CMAKE_SOURCE_DIR}/models/surfaceparticipationverdict.json")
set(VERDICT_BARS_HEADER "${CMAKE_BINARY_DIR}/generated/embedded_verdict_bars.hpp")

if(NOT EXISTS "${VERDICT_BARS_FILE}")
  message(FATAL_ERROR "Verdict data file not found: ${VERDICT_BARS_FILE}")
endif()

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/generated)

# Regenerate the header during build (not configure) when the data file changes
add_custom_command(
  OUTPUT ${VERDICT_BARS_HEADER}
  COMMAND ${CMAKE_COMMAND}
    -DVERDICT_BARS_FILE=${VERDICT_BARS_FILE}
    -DVERDICT_BARS_HEADER=${VERDICT_BARS_HEADER}
    -P ${CMAKE_CURRENT_LIST_DIR}/embed_verdict_bars.cmake
  DEPENDS ${VERDICT_BARS_FILE} ${CMAKE_CURRENT_LIST_DIR}/embed_verdict_bars.cmake
  COMMENT "Regenerating the embedded verdict data"
)
add_custom_target(embedded_verdict_bars DEPENDS ${VERDICT_BARS_HEADER})

# Make the library target depend on the generated header
add_dependencies(${LIB_TARGET_NAME} embedded_verdict_bars)
