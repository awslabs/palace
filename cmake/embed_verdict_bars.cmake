# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Helper script for generating embedded_verdict_bars.hpp
# Called via: cmake -DVERDICT_BARS_FILE=<json> -DVERDICT_BARS_HEADER=<output> -P embed_verdict_bars.cmake

if(NOT DEFINED VERDICT_BARS_FILE OR NOT DEFINED VERDICT_BARS_HEADER)
  message(FATAL_ERROR "VERDICT_BARS_FILE and VERDICT_BARS_HEADER must be defined")
endif()

file(READ "${VERDICT_BARS_FILE}" FILE_CONTENT)

# Escape for C++ raw string (just need to avoid )json" sequence)
string(REPLACE ")json\"" ")json_\"" FILE_CONTENT "${FILE_CONTENT}")

file(WRITE "${VERDICT_BARS_HEADER}"
"// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Auto-generated file - do not edit
// Generated from palace/models/surfaceparticipationverdict.json: the built-in verdict data
// of the per-iteration table surface-participation.csv.

#ifndef PALACE_EMBEDDED_VERDICT_BARS_HPP
#define PALACE_EMBEDDED_VERDICT_BARS_HPP

namespace palace::verdict
{

inline const char *surface_participation_verdict_json = R\"json(
${FILE_CONTENT})json\";

}  // namespace palace::verdict

#endif  // PALACE_EMBEDDED_VERDICT_BARS_HPP
")
