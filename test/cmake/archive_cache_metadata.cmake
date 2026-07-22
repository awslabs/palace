# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.21)

if(NOT DEFINED EXTERNAL_GIT_TAGS)
  message(FATAL_ERROR "Missing required test argument EXTERNAL_GIT_TAGS")
endif()

if(DEFINED TEST_MODE)
  set(PALACE_BUILD_EXTERNAL_DEPS ON)
  if(TEST_MODE STREQUAL "coherent")
    set(EXTERN_FMT_VERSION "99.0" CACHE STRING "" FORCE)
    set(EXTERN_FMT_URL
      "https://github.com/fmtlib/fmt/releases/download/99.0/fmt-99.0.zip"
      CACHE STRING "" FORCE)
  elseif(TEST_MODE STREQUAL "opaque_override")
    set(EXTERN_FMT_VERSION "99.0" CACHE STRING "" FORCE)
    set(EXTERN_FMT_URL "file:///cache/fmt.zip" CACHE STRING "" FORCE)
  elseif(TEST_MODE STREQUAL "mismatch")
    set(EXTERN_FMT_VERSION "99.0" CACHE STRING "" FORCE)
    set(EXTERN_FMT_URL
      "https://github.com/fmtlib/fmt/releases/download/12.1.0/fmt-12.1.0.zip"
      CACHE STRING "" FORCE)
  elseif(TEST_MODE STREQUAL "overlap")
    set(EXTERN_FMT_VERSION "2.1.0" CACHE STRING "" FORCE)
    set(EXTERN_FMT_URL
      "https://github.com/fmtlib/fmt/releases/download/12.1.0/fmt-12.1.0.zip"
      CACHE STRING "" FORCE)
  elseif(TEST_MODE STREQUAL "external_off")
    set(PALACE_BUILD_EXTERNAL_DEPS OFF)
    set(EXTERN_FMT_VERSION "99.0" CACHE STRING "" FORCE)
    set(EXTERN_FMT_URL
      "https://github.com/fmtlib/fmt/releases/download/12.1.0/fmt-12.1.0.zip"
      CACHE STRING "" FORCE)
  else()
    message(FATAL_ERROR "Unknown archive cache test mode '${TEST_MODE}'")
  endif()
  include("${EXTERNAL_GIT_TAGS}")
  return()
endif()

foreach(_mode IN ITEMS coherent opaque_override external_off mismatch overlap)
  execute_process(
    COMMAND "${CMAKE_COMMAND}" "-DEXTERNAL_GIT_TAGS=${EXTERNAL_GIT_TAGS}"
      "-DTEST_MODE=${_mode}" -P "${CMAKE_CURRENT_LIST_FILE}"
    RESULT_VARIABLE _result
    OUTPUT_VARIABLE _output
    ERROR_VARIABLE _error
  )
  if(_mode STREQUAL "mismatch" OR _mode STREQUAL "overlap")
    if(_result EQUAL 0)
      message(FATAL_ERROR
        "Stale archive version/URL cache mode '${_mode}' unexpectedly succeeded"
      )
    endif()
    if(NOT "${_output}${_error}" MATCHES "Update both EXTERN_FMT_VERSION" OR
       NOT "${_output}${_error}" MATCHES "EXTERN_FMT_URL")
      message(FATAL_ERROR
        "Stale archive cache failure did not explain how to fix both values:\n${_output}${_error}"
      )
    endif()
  elseif(NOT _result EQUAL 0)
    message(FATAL_ERROR "Archive cache test mode '${_mode}' failed:\n${_output}${_error}")
  endif()
endforeach()
