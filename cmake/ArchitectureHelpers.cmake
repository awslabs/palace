# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

function(palace_escape_external_project_list output values)
  string(REPLACE ";" "$<SEMICOLON>" escaped "${values}")
  set(${output} "${escaped}" PARENT_SCOPE)
endfunction()

# Append -DCMAKE_CUDA_ARCHITECTURES=<list> to an ExternalProject options list.
# ExternalProject_Add splits its command on ";", so the semicolon-separated
# architecture list is escaped or only the first entry survives. No-op when no
# architectures are set, so callers can invoke this unconditionally from within
# their CUDA branch.
function(palace_append_cuda_architectures options_var)
  if("${CMAKE_CUDA_ARCHITECTURES}" STREQUAL "")
    return()
  endif()
  palace_escape_external_project_list(escaped "${CMAKE_CUDA_ARCHITECTURES}")
  set(options "${${options_var}}")
  list(APPEND options "-DCMAKE_CUDA_ARCHITECTURES=${escaped}")
  set(${options_var} "${options}" PARENT_SCOPE)
endfunction()

function(palace_parse_cuda_architecture architecture output_number output_kind)
  # Accept bare compute capabilities (80), the arch-conditional feature suffixes
  # nvcc uses for Hopper/Blackwell (90a, 100f), and CMake's -real/-virtual
  # qualifiers. Special values such as native/all/all-major are intentionally
  # rejected here: they cannot be expanded into explicit gencode flags without a
  # GPU at configure time, which the Makefile-based dependencies require.
  if(NOT architecture MATCHES "^([0-9]+[af]?)(-(real|virtual))?$")
    message(FATAL_ERROR
      "Unsupported CUDA architecture: ${architecture}. Expected forms like "
      "80, 90a, 90-real, or 100f-virtual."
    )
  endif()
  set(${output_number} "${CMAKE_MATCH_1}" PARENT_SCOPE)
  set(${output_kind} "${CMAKE_MATCH_3}" PARENT_SCOPE)
endfunction()

function(palace_cuda_gencode_flags output architectures)
  # Unqualified targets generate cubins only; request -virtual explicitly for PTX.
  set(flags)
  foreach(arch IN LISTS architectures)
    palace_parse_cuda_architecture("${arch}" arch_number arch_kind)
    if(arch_kind STREQUAL "virtual")
      string(APPEND flags
        " --generate-code arch=compute_${arch_number},code=compute_${arch_number}"
      )
    else()
      string(APPEND flags
        " --generate-code arch=compute_${arch_number},code=sm_${arch_number}"
      )
    endif()
  endforeach()
  string(STRIP "${flags}" flags)
  set(${output} "${flags}" PARENT_SCOPE)
endfunction()

# Convert CMake architectures to Make-style targets for dependencies that
# control their own PTX generation.
function(palace_cuda_sm_targets output architectures)
  set(targets)
  foreach(arch IN LISTS architectures)
    palace_parse_cuda_architecture("${arch}" arch_number arch_kind)
    list(APPEND targets "sm_${arch_number}")
  endforeach()
  list(REMOVE_DUPLICATES targets)
  list(JOIN targets " " targets)
  set(${output} "${targets}" PARENT_SCOPE)
endfunction()
