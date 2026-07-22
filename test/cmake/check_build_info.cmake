# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.21)

foreach(_required IN ITEMS PALACE_WRAPPER PALACE_EXECUTABLE PALACE_PROJECT_VERSION
                           TEST_BINARY_DIR)
  if(NOT DEFINED ${_required})
    message(FATAL_ERROR "Missing required test argument ${_required}")
  endif()
endforeach()
if(NOT DEFINED PALACE_BUILD_DEPENDENCIES)
  set(PALACE_BUILD_DEPENDENCIES "")
endif()

# Prefer the installed public wrapper. Falling back is only valid when installation has
# not made that wrapper available (for example, a direct target-only developer build).
if(EXISTS "${PALACE_WRAPPER}" AND NOT IS_DIRECTORY "${PALACE_WRAPPER}")
  set(_palace_command "${PALACE_WRAPPER}")
  get_filename_component(_test_binary_dir "${TEST_BINARY_DIR}" ABSOLUTE)
  set(_nonexistent_launcher "${_test_binary_dir}/build-info-nonexistent-launcher")
  if(EXISTS "${_nonexistent_launcher}")
    message(FATAL_ERROR "The nonexistent launcher test path unexpectedly exists")
  endif()
  set(_palace_command_arguments --launcher "${_nonexistent_launcher}")
else()
  if(NOT EXISTS "${PALACE_EXECUTABLE}")
    message(FATAL_ERROR
      "Installed Palace wrapper '${PALACE_WRAPPER}' and direct executable "
      "'${PALACE_EXECUTABLE}' are both unavailable"
    )
  endif()
  set(_palace_command "${PALACE_EXECUTABLE}")
  set(_palace_command_arguments)
endif()

foreach(_flag IN ITEMS -V --version)
  execute_process(
    COMMAND "${_palace_command}" ${_palace_command_arguments} "${_flag}"
    RESULT_VARIABLE _result
    OUTPUT_VARIABLE _output
    ERROR_VARIABLE _error
  )
  if(NOT _result EQUAL 0)
    message(FATAL_ERROR
      "'${_palace_command} ${_flag}' failed with ${_result}:\n${_output}${_error}"
    )
  endif()

  string(REPLACE "\r\n" "\n" _payload "${_output}")
  string(REPLACE "\r" "\n" _payload "${_payload}")
  # The wrapper logs the normalized command before exec; it is not version payload.
  string(REGEX REPLACE "^>> [^\n]*\n" "" _payload "${_payload}")
  if(_flag STREQUAL "-V")
    set(_short_payload "${_payload}")
  else()
    set(_long_payload "${_payload}")
  endif()
endforeach()

if(NOT _short_payload STREQUAL _long_payload)
  message(FATAL_ERROR
    "-V and --version returned different normalized payloads:\n"
    "--- -V ---\n${_short_payload}--- --version ---\n${_long_payload}"
  )
endif()

string(REGEX MATCH "Palace version: [^\n]*" _reported_version "${_short_payload}")
if(NOT _reported_version STREQUAL "Palace version: ${PALACE_PROJECT_VERSION}")
  message(FATAL_ERROR
    "Expected exact project version '${PALACE_PROJECT_VERSION}', got '${_reported_version}'"
  )
endif()

# Derive the expected display independently of BuildDependencyMetadata.cmake. This is deliberately
# a second implementation of the small public wire contract used to test the CLI boundary.
set(_expected_rows)
if(NOT PALACE_BUILD_DEPENDENCIES STREQUAL "")
  string(REPLACE "|" ";" _records "${PALACE_BUILD_DEPENDENCIES}")
  foreach(_record IN LISTS _records)
    string(FIND "${_record}" "=" _equals_position)
    if(_equals_position LESS 1)
      message(FATAL_ERROR "Test received malformed metadata record '${_record}'")
    endif()
    math(EXPR _value_position "${_equals_position} + 1")
    string(SUBSTRING "${_record}" 0 ${_equals_position} _label)
    string(SUBSTRING "${_record}" ${_value_position} -1 _value)
    list(APPEND _expected_rows "  ${_label}: ${_value}")
  endforeach()
endif()
list(JOIN _expected_rows "\n" _expected_dependencies)

string(FIND "${_short_payload}" "\nBuild dependencies:\n" _block_position)
if(PALACE_BUILD_DEPENDENCIES STREQUAL "")
  if(NOT _block_position EQUAL -1)
    message(FATAL_ERROR "Empty metadata unexpectedly printed a dependency block")
  endif()
else()
  if(_block_position EQUAL -1)
    message(FATAL_ERROR "Version output is missing its dependency block")
  endif()
  string(SUBSTRING "${_short_payload}" ${_block_position} -1 _actual_block)
  set(_expected_block "\nBuild dependencies:\n${_expected_dependencies}\n")
  if(NOT _actual_block STREQUAL _expected_block)
    message(FATAL_ERROR
      "Dependency block mismatch:\n--- actual ---${_actual_block}"
      "--- expected ---${_expected_block}"
    )
  endif()
endif()
