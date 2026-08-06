# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/..")
include(ArchitectureHelpers)

set(architectures "80;86;89;90")

palace_escape_external_project_list(escaped "${architectures}")
set(expected_escaped "80$<SEMICOLON>86$<SEMICOLON>89$<SEMICOLON>90")
if(NOT escaped STREQUAL expected_escaped)
  message(FATAL_ERROR "Unexpected escaped architecture list: ${escaped}")
endif()

palace_cuda_gencode_flags(cuda_flags "${architectures}")
foreach(arch IN LISTS architectures)
  if(NOT cuda_flags MATCHES "arch=compute_${arch},code=sm_${arch}")
    message(FATAL_ERROR "Missing sm_${arch} in CUDA flags: ${cuda_flags}")
  endif()
  if(cuda_flags MATCHES "arch=compute_${arch},code=compute_${arch}")
    message(FATAL_ERROR "Unexpected compute_${arch} in CUDA flags: ${cuda_flags}")
  endif()
endforeach()

palace_cuda_gencode_flags(cuda_qualified_flags "80-real;90-virtual")
if(NOT cuda_qualified_flags MATCHES "code=sm_80")
  message(FATAL_ERROR "Missing real sm_80 code: ${cuda_qualified_flags}")
endif()
if(cuda_qualified_flags MATCHES "code=compute_80")
  message(FATAL_ERROR "Unexpected virtual compute_80 code: ${cuda_qualified_flags}")
endif()
if(cuda_qualified_flags MATCHES "code=sm_90")
  message(FATAL_ERROR "Unexpected real sm_90 code: ${cuda_qualified_flags}")
endif()
if(NOT cuda_qualified_flags MATCHES "code=compute_90")
  message(FATAL_ERROR "Missing virtual compute_90 code: ${cuda_qualified_flags}")
endif()

palace_parse_cuda_architecture("90-virtual" cuda_arch cuda_kind)
if(NOT cuda_arch STREQUAL "90" OR NOT cuda_kind STREQUAL "virtual")
  message(FATAL_ERROR "Unexpected parsed CUDA architecture: ${cuda_arch} ${cuda_kind}")
endif()
