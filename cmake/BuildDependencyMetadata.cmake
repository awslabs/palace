# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

include_guard(GLOBAL)

# Append one producer-owned field to a list of dependency metadata records. The wire
# format itself uses pipe separators, so producer fields must not contain reserved
# delimiters. Values may contain '='; consumers split each record on its first '='.
function(palace_append_build_dependency output_variable label value)
  if(NOT ARGC EQUAL 3)
    message(FATAL_ERROR
      "Build dependency metadata fields must not contain semicolons; expected a label and value"
    )
  endif()

  if(NOT label MATCHES "^[A-Za-z][A-Za-z0-9_]*$")
    message(FATAL_ERROR
      "Invalid build dependency label '${label}'; expected [A-Za-z][A-Za-z0-9_]*"
    )
  endif()
  if(value STREQUAL "")
    message(FATAL_ERROR "Build dependency '${label}' has an empty value")
  endif()
  foreach(_palace_field_name IN ITEMS label value)
    foreach(_palace_reserved IN ITEMS ";" "|" "\r" "\n")
      string(FIND "${${_palace_field_name}}" "${_palace_reserved}"
        _palace_reserved_position)
      if(NOT _palace_reserved_position EQUAL -1)
        message(FATAL_ERROR
          "Build dependency ${_palace_field_name} '${${_palace_field_name}}' "
          "contains a reserved delimiter"
        )
      endif()
    endforeach()
  endforeach()

  set(_palace_entries "${${output_variable}}")
  foreach(_palace_entry IN LISTS _palace_entries)
    string(FIND "${_palace_entry}" "=" _palace_equals_position)
    if(_palace_equals_position GREATER 0)
      string(SUBSTRING "${_palace_entry}" 0 ${_palace_equals_position}
        _palace_existing_label)
      if(_palace_existing_label STREQUAL label)
        message(FATAL_ERROR "Duplicate build dependency label '${label}'")
      endif()
    endif()
  endforeach()
  list(APPEND _palace_entries "${label}=${value}")
  set(${output_variable} "${_palace_entries}" PARENT_SCOPE)
endfunction()

# Append metadata only when the corresponding superbuild ExternalProject target exists.
function(palace_append_external_build_dependency output_variable target label value)
  if(NOT ARGC EQUAL 4)
    message(FATAL_ERROR
      "External build dependency metadata fields must not contain semicolons"
    )
  endif()
  if(NOT TARGET "${target}")
    return()
  endif()

  set(_palace_entries "${${output_variable}}")
  palace_append_build_dependency(_palace_entries "${label}" "${value}")
  set(${output_variable} "${_palace_entries}" PARENT_SCOPE)
endfunction()

# Validate the public label=value|... wire value and produce its display form plus
# escaped C++ string contents. Every display row owns its two-space indent; the final row
# deliberately has no trailing newline.
function(palace_format_build_dependencies wire_value formatted_output
         cpp_string_output)
  if(NOT ARGC EQUAL 3)
    message(FATAL_ERROR
      "PALACE_BUILD_DEPENDENCIES must not contain semicolons; expected one wire value"
    )
  endif()

  foreach(_palace_reserved IN ITEMS ";" "\r" "\n")
    string(FIND "${wire_value}" "${_palace_reserved}" _palace_reserved_position)
    if(NOT _palace_reserved_position EQUAL -1)
      message(FATAL_ERROR
        "PALACE_BUILD_DEPENDENCIES contains a reserved semicolon or line break"
      )
    endif()
  endforeach()

  if(wire_value STREQUAL "")
    set(${formatted_output} "" PARENT_SCOPE)
    set(${cpp_string_output} "" PARENT_SCOPE)
    return()
  endif()
  if(wire_value MATCHES "^\\|" OR wire_value MATCHES "\\|\\|" OR
     wire_value MATCHES "\\|$")
    message(FATAL_ERROR "PALACE_BUILD_DEPENDENCIES contains an empty record")
  endif()

  string(REPLACE "|" ";" _palace_records "${wire_value}")
  set(_palace_labels)
  set(_palace_rows)
  foreach(_palace_record IN LISTS _palace_records)
    if(_palace_record STREQUAL "")
      message(FATAL_ERROR "PALACE_BUILD_DEPENDENCIES contains an empty record")
    endif()

    # Values can contain '='. Only the first one separates the label and value.
    string(FIND "${_palace_record}" "=" _palace_equals_position)
    if(_palace_equals_position LESS 1)
      message(FATAL_ERROR
        "Malformed build dependency record '${_palace_record}'; expected label=value"
      )
    endif()
    string(LENGTH "${_palace_record}" _palace_record_length)
    math(EXPR _palace_value_position "${_palace_equals_position} + 1")
    if(_palace_value_position GREATER_EQUAL _palace_record_length)
      message(FATAL_ERROR "Build dependency record '${_palace_record}' has an empty value")
    endif()

    string(SUBSTRING "${_palace_record}" 0 ${_palace_equals_position} _palace_label)
    string(SUBSTRING "${_palace_record}" ${_palace_value_position} -1 _palace_value)
    if(NOT _palace_label MATCHES "^[A-Za-z][A-Za-z0-9_]*$")
      message(FATAL_ERROR
        "Invalid build dependency label '${_palace_label}'; expected [A-Za-z][A-Za-z0-9_]*"
      )
    endif()
    list(FIND _palace_labels "${_palace_label}" _palace_duplicate_position)
    if(NOT _palace_duplicate_position EQUAL -1)
      message(FATAL_ERROR "Duplicate build dependency label '${_palace_label}'")
    endif()
    list(APPEND _palace_labels "${_palace_label}")
    list(APPEND _palace_rows "  ${_palace_label}: ${_palace_value}")
  endforeach()

  list(JOIN _palace_rows "\n" _palace_formatted)

  set(_palace_cpp_string "${_palace_formatted}")
  string(REPLACE "\\" "\\\\" _palace_cpp_string "${_palace_cpp_string}")
  string(REPLACE "\"" "\\\"" _palace_cpp_string "${_palace_cpp_string}")
  string(REPLACE "\n" "\\n" _palace_cpp_string "${_palace_cpp_string}")

  set(${formatted_output} "${_palace_formatted}" PARENT_SCOPE)
  set(${cpp_string_output} "${_palace_cpp_string}" PARENT_SCOPE)
endfunction()
