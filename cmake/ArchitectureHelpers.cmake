# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

function(palace_escape_external_project_list output values)
  string(REPLACE ";" "$<SEMICOLON>" escaped "${values}")
  set(${output} "${escaped}" PARENT_SCOPE)
endfunction()

function(palace_parse_cuda_architecture architecture output_number output_kind)
  if(NOT architecture MATCHES "^([0-9]+)(-(real|virtual))?$")
    message(FATAL_ERROR "Unsupported CUDA architecture: ${architecture}")
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
