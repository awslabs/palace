# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Build libCEED
#

# Force build order
if(PALACE_BUILD_EXTERNAL_DEPS AND PALACE_WITH_LIBXSMM)
  set(LIBCEED_DEPENDENCIES libxsmm)
else()
  set(LIBCEED_DEPENDENCIES)
endif()

# Note on recommended flags for libCEED (from Makefile, Spack):
#   OPT: -O3 -g -march=native -ffp-contract=fast [-fopenmp-simd/-qopenmp-simd]
include(CheckCCompilerFlag)
set(LIBCEED_C_FLAGS "${CMAKE_C_FLAGS}")
if(CMAKE_C_COMPILER_ID MATCHES "Intel|IntelLLVM")
  set(OMP_SIMD_FLAG -qopenmp-simd)
else()
  set(OMP_SIMD_FLAG -fopenmp-simd)
endif()
check_c_compiler_flag(${OMP_SIMD_FLAG} SUPPORTS_OMP_SIMD)
if(SUPPORTS_OMP_SIMD)
  set(LIBCEED_C_FLAGS "${LIBCEED_C_FLAGS} ${OMP_SIMD_FLAG}")
endif()

# Build libCEED (always as a shared library)
set(LIBCEED_OPTIONS
  "prefix=${CMAKE_INSTALL_PREFIX}"
  "CC=${CMAKE_C_COMPILER}"
  "OPT=${LIBCEED_C_FLAGS}"
  "STATIC="
)

# Configure OpenMP
if(PALACE_WITH_OPENMP)
  list(APPEND LIBCEED_OPTIONS
    "OPENMP=1"
  )
endif()

# Configure libCEED backends (disable CUDA for now, AVX handled automatically)
list(APPEND LIBCEED_OPTIONS
  "CUDA_DIR=/disable-cuda"
)
if(PALACE_WITH_LIBXSMM)
  if(PALACE_BUILD_EXTERNAL_DEPS)
    list(APPEND LIBCEED_OPTIONS
      "XSMM_DIR=${CMAKE_INSTALL_PREFIX}"
    )
  else()
    list(APPEND LIBCEED_OPTIONS
      "XSMM_DIR=${LIBXSMM_DIR}"
    )
  endif()
  # LIBXSMM can require linkage with BLAS for fallback
  if(NOT "${BLAS_LAPACK_LIBRARIES}" STREQUAL "")
    string(REPLACE "$<SEMICOLON>" " " LIBCEED_BLAS_LAPACK_LIBRARIES "${BLAS_LAPACK_LIBRARIES}")
    list(APPEND LIBCEED_OPTIONS
      "BLAS_LIB=${LIBCEED_BLAS_LAPACK_LIBRARIES}"
    )
  endif()
endif()

string(REPLACE ";" "; " LIBCEED_OPTIONS_PRINT "${LIBCEED_OPTIONS}")
message(STATUS "LIBCEED_OPTIONS: ${LIBCEED_OPTIONS_PRINT}")

# Add OpenMP support to libCEED
set(LIBCEED_PATCH_FILES
  "${CMAKE_SOURCE_DIR}/extern/patch/libCEED/patch_gpu_restriction_dev.diff"
  "${CMAKE_SOURCE_DIR}/extern/patch/libCEED/patch_hcurl_hdiv_basis.diff"
)

include(ExternalProject)
ExternalProject_Add(libCEED
  DEPENDS           ${LIBCEED_DEPENDENCIES}
  GIT_REPOSITORY    ${EXTERN_LIBCEED_URL}
  GIT_TAG           ${EXTERN_LIBCEED_GIT_TAG}
  SOURCE_DIR        ${CMAKE_BINARY_DIR}/extern/libCEED
  INSTALL_DIR       ${CMAKE_INSTALL_PREFIX}
  PREFIX            ${CMAKE_BINARY_DIR}/extern/libCEED-cmake
  BUILD_IN_SOURCE   TRUE
  UPDATE_COMMAND    ""
  PATCH_COMMAND     git apply "${LIBCEED_PATCH_FILES}"
  CONFIGURE_COMMAND ""
  BUILD_COMMAND     ""
  INSTALL_COMMAND   ${CMAKE_MAKE_PROGRAM} ${LIBCEED_OPTIONS} install
  TEST_COMMAND      ""
)
