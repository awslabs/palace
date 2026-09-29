# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

set(PALACE_BUILD_INFO_MODULE_DIR "${CMAKE_CURRENT_LIST_DIR}")

# Record one locally discovered dependency. Producers can instead supply the complete
# PALACE_DEPENDENCY_MANIFEST when they know the concrete dependency set and versions.
function(palace_add_build_info_dependency name version)
  if("${name}" STREQUAL "" OR "${name}" MATCHES "[=|;]" OR "${version}" MATCHES "[=|;]")
    message(FATAL_ERROR "Invalid build-information dependency: ${name}=${version}")
  endif()
  if("${version}" STREQUAL "")
    set(version "unknown")
  endif()
  set_property(GLOBAL APPEND PROPERTY PALACE_BUILD_INFO_LOCAL_MANIFEST "${name}=${version}")
endfunction()

function(palace_escape_build_info_string value output_variable)
  string(REPLACE "\\" "\\\\" escaped "${value}")
  string(REPLACE "\"" "\\\"" escaped "${escaped}")
  string(REPLACE "\n" "\\n" escaped "${escaped}")
  string(REPLACE "\r" "\\r" escaped "${escaped}")
  string(REPLACE "\t" "\\t" escaped "${escaped}")
  set("${output_variable}" "${escaped}" PARENT_SCOPE)
endfunction()

function(palace_configure_build_info target)
  if("${PALACE_DEPENDENCY_MANIFEST}" STREQUAL "")
    get_property(manifest_entries GLOBAL PROPERTY PALACE_BUILD_INFO_LOCAL_MANIFEST)
    list(REMOVE_DUPLICATES manifest_entries)
    string(JOIN "|" PALACE_DEPENDENCY_MANIFEST ${manifest_entries})
  else()
    string(REPLACE "|" ";" manifest_entries "${PALACE_DEPENDENCY_MANIFEST}")
  endif()

  set(PALACE_BUILD_INFO_DEPENDENCIES "")
  foreach(entry IN LISTS manifest_entries)
    string(FIND "${entry}" "=" separator)
    if(separator LESS 1)
      message(FATAL_ERROR "Invalid PALACE_DEPENDENCY_MANIFEST entry: ${entry}")
    endif()
    string(SUBSTRING "${entry}" 0 ${separator} name)
    math(EXPR value_start "${separator} + 1")
    string(SUBSTRING "${entry}" ${value_start} -1 version)
    if("${name}" MATCHES "[=|;]" OR "${version}" MATCHES "[=|;]")
      message(FATAL_ERROR "Unsupported PALACE_DEPENDENCY_MANIFEST delimiter in: ${entry}")
    endif()
    palace_escape_build_info_string("${name}" escaped_name)
    palace_escape_build_info_string("${version}" escaped_version)
    string(APPEND PALACE_BUILD_INFO_DEPENDENCIES "  {\"${escaped_name}\", \"${escaped_version}\"},\n")
  endforeach()

  if("${PALACE_BUILD_SYSTEM}" STREQUAL "")
    set(PALACE_BUILD_SYSTEM "CMake")
  endif()
  set(generator "${CMAKE_BINARY_DIR}/generated/GenerateBuildInfo.cmake")
  configure_file("${PALACE_BUILD_INFO_MODULE_DIR}/GenerateBuildInfo.cmake.in" "${generator}" @ONLY)
  add_custom_target(palace_build_info
    COMMAND ${CMAKE_COMMAND} -P "${generator}"
    BYPRODUCTS "${CMAKE_BINARY_DIR}/generated/BuildInfo.hpp"
    VERBATIM)
  add_dependencies(${target} palace_build_info)
endfunction()
