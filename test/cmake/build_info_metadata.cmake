# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.21)

foreach(_required IN ITEMS MODULE_FILE HEADER_TEMPLATE TEST_BINARY_DIR)
  if(NOT DEFINED ${_required})
    message(FATAL_ERROR "Missing required test argument ${_required}")
  endif()
endforeach()

include("${MODULE_FILE}")

# Additional equals signs belong to the value. Quotes, backslashes, and the old raw-string
# terminator must all survive formatting and produce escaped ordinary C++ string contents.
set(_metadata [=[alpha=first=second|raw=)palace"|slash=C:\path\"quoted]=])
palace_format_build_dependencies("${_metadata}" _formatted _cpp_string)
set(_expected_formatted [=[  alpha: first=second
  raw: )palace"
  slash: C:\path\"quoted]=])
set(_expected_cpp_string [=[  alpha: first=second\n  raw: )palace\"\n  slash: C:\\path\\\"quoted]=])
if(NOT _formatted STREQUAL _expected_formatted)
  message(FATAL_ERROR
    "Unexpected dependency formatting:\n--- actual ---\n${_formatted}\n"
    "--- expected ---\n${_expected_formatted}"
  )
endif()
if(NOT _cpp_string STREQUAL _expected_cpp_string)
  message(FATAL_ERROR
    "Unexpected escaped C++ string:\n--- actual ---\n${_cpp_string}\n"
    "--- expected ---\n${_expected_cpp_string}"
  )
endif()
string(LENGTH "${_formatted}" _formatted_length)
math(EXPR _last_character_position "${_formatted_length} - 1")
string(SUBSTRING "${_formatted}" ${_last_character_position} 1 _last_character)
if(_last_character STREQUAL "\n")
  message(FATAL_ERROR "Formatted dependency metadata has a trailing newline")
endif()

palace_format_build_dependencies("" _empty_formatted _empty_cpp_string)
if(NOT _empty_formatted STREQUAL "" OR NOT _empty_cpp_string STREQUAL "")
  message(FATAL_ERROR "The empty wire value did not produce empty output")
endif()

# configure_file must replace an existing header when either version or metadata changes.
file(MAKE_DIRECTORY "${TEST_BINARY_DIR}")
set(_generated_header "${TEST_BINARY_DIR}/build_info.hpp")
set(PROJECT_VERSION "1.2.3")
set(PALACE_BUILD_DEPENDENCIES_CPP_STRING "${_cpp_string}")
configure_file("${HEADER_TEMPLATE}" "${_generated_header}" @ONLY)
file(READ "${_generated_header}" _first_header)
if(NOT _first_header MATCHES "palace_version = \"1\\.2\\.3\"")
  message(FATAL_ERROR "The first generated header has the wrong project version")
endif()
if(NOT _first_header MATCHES "first=second")
  message(FATAL_ERROR "The first generated header is missing dependency metadata")
endif()

set(PROJECT_VERSION "9.8.7")
palace_format_build_dependencies("updated=2.0" _formatted
  PALACE_BUILD_DEPENDENCIES_CPP_STRING)
configure_file("${HEADER_TEMPLATE}" "${_generated_header}" @ONLY)
file(READ "${_generated_header}" _second_header)
if(NOT _second_header MATCHES "palace_version = \"9\\.8\\.7\"")
  message(FATAL_ERROR "The regenerated header has the wrong project version")
endif()
if(NOT _second_header MATCHES "formatted_dependency_versions =\n    \"  updated: 2\\.0\";")
  message(FATAL_ERROR "The regenerated header has the wrong dependency string")
endif()
if(_second_header MATCHES "1\\.2\\.3|first=second|palace\\\"")
  message(FATAL_ERROR "The regenerated header retained old version or metadata content")
endif()

# Run invalid inputs in child CMake processes so this test can require every case to fail.
foreach(_case IN ITEMS no_equals empty_record leading_empty_record trailing_empty_record
                       empty_label empty_value invalid_label duplicate semicolon carriage_return
                       line_feed producer_pipe producer_semicolon producer_carriage_return
                       producer_line_feed producer_invalid_label producer_empty_label
                       producer_empty_value producer_duplicate)
  execute_process(
    COMMAND "${CMAKE_COMMAND}" "-DMODULE_FILE=${MODULE_FILE}" "-DINVALID_CASE=${_case}"
      -P "${CMAKE_CURRENT_LIST_DIR}/invalid_build_metadata.cmake"
    RESULT_VARIABLE _invalid_result
    OUTPUT_VARIABLE _invalid_output
    ERROR_VARIABLE _invalid_error
  )
  if(_invalid_result EQUAL 0)
    message(FATAL_ERROR "Invalid metadata case '${_case}' unexpectedly succeeded")
  endif()
endforeach()
