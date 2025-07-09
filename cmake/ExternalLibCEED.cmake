# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Build libCEED
#

# Force build order
set(LIBCEED_DEPENDENCIES)
if(PALACE_WITH_LIBXSMM)
  list(APPEND LIBCEED_DEPENDENCIES libxsmm)
endif()
if(PALACE_WITH_MAGMA)
  list(APPEND LIBCEED_DEPENDENCIES magma)
endif()

# Note on recommended flags for libCEED (from Makefile, Spack):
#   OPT: -O3 -g -march=native -ffp-contract=fast [-fopenmp-simd/-qopenmp-simd]
include(CheckCCompilerFlag)
set(LIBCEED_OPT_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${BUILD_TYPE_UPPER}}")
if(CMAKE_C_COMPILER_ID MATCHES "Intel|IntelLLVM")
  set(OMP_SIMD_FLAG -qopenmp-simd)
else()
  set(OMP_SIMD_FLAG -fopenmp-simd)
endif()
check_c_compiler_flag(${OMP_SIMD_FLAG} SUPPORTS_OMP_SIMD)
if(SUPPORTS_OMP_SIMD)
  set(LIBCEED_OPT_FLAGS "${LIBCEED_OPT_FLAGS} ${OMP_SIMD_FLAG}")
endif()

# Silence some CUDA/HIP include file warnings
if(PALACE_WITH_CUDA)
  set(LIBCEED_OPT_FLAGS "${LIBCEED_OPT_FLAGS} -isystem ${CUDAToolkit_INCLUDE_DIRS}")
endif()
if(PALACE_WITH_HIP)
  set(LIBCEED_OPT_FLAGS "${LIBCEED_OPT_FLAGS} -isystem ${ROCM_DIR}/include")
endif()

# Configure -pedantic flag if specified (don't want to enable for GPU code)
if(LIBCEED_OPT_FLAGS MATCHES "-pedantic")
  string(REGEX REPLACE "-pedantic" "" LIBCEED_OPT_FLAGS ${LIBCEED_OPT_FLAGS})
  set(LIBCEED_PEDANTIC "1")
else()
  set(LIBCEED_PEDANTIC "")
endif()

# Build libCEED (always as a shared library)
set(LIBCEED_OPTIONS
  "prefix=${CMAKE_INSTALL_PREFIX}"
  "LDFLAGS=${CMAKE_EXE_LINKER_FLAGS}"
  "CC=${CMAKE_C_COMPILER}"
  "CXX=${CMAKE_CXX_COMPILER}"
  "FC="
  "OPT=${LIBCEED_OPT_FLAGS}"
  "STATIC="
  "PEDANTIC=${LIBCEED_PEDANTIC}"
)

# Configure OpenMP
if(PALACE_WITH_OPENMP)
  list(APPEND LIBCEED_OPTIONS
    "OPENMP=1"
  )
endif()

# Configure libCEED backends (nvcc, hipcc flags are configured by libCEED)
if(PALACE_WITH_LIBXSMM)
  list(APPEND LIBCEED_OPTIONS
    "XSMM_DIR=${CMAKE_INSTALL_PREFIX}"
  )
  # LIBXSMM can require linkage with BLAS for fallback
  if(NOT "${BLAS_LAPACK_LIBRARIES}" STREQUAL "")
    string(REPLACE "$<SEMICOLON>" " " LIBCEED_BLAS_LAPACK_LIBRARIES "${BLAS_LAPACK_LIBRARIES}")
    list(APPEND LIBCEED_OPTIONS
      "BLAS_LIB=${LIBCEED_BLAS_LAPACK_LIBRARIES}"
    )
  endif()
endif()
if(PALACE_WITH_CUDA)
  list(APPEND LIBCEED_OPTIONS
    "CUDA_DIR=${CUDAToolkit_LIBRARY_ROOT}"
  )
  if(NOT "${CMAKE_CUDA_ARCHITECTURES}" STREQUAL "")
    list(GET CMAKE_CUDA_ARCHITECTURES 0 LIBCEED_CUDA_ARCH)
    list(APPEND LIBCEED_OPTIONS
      "CUDA_ARCH=sm_${LIBCEED_CUDA_ARCH}"
    )
  endif()
endif()
if(PALACE_WITH_HIP)
  list(APPEND LIBCEED_OPTIONS
    "ROCM_DIR=${ROCM_DIR}"
  )
  if(NOT "${CMAKE_HIP_ARCHITECTURES}" STREQUAL "")
    list(GET CMAKE_HIP_ARCHITECTURES 0 LIBCEED_HIP_ARCH)
    list(APPEND LIBCEED_OPTIONS
      "HIP_ARCH=${LIBCEED_HIP_ARCH}"
    )
  endif()
endif()
if(PALACE_WITH_MAGMA)
  list(APPEND LIBCEED_OPTIONS
    "MAGMA_DIR=${CMAKE_INSTALL_PREFIX}"
  )
endif()

string(REPLACE ";" "; " LIBCEED_OPTIONS_PRINT "${LIBCEED_OPTIONS}")
message(STATUS "LIBCEED_OPTIONS: ${LIBCEED_OPTIONS_PRINT}")

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
  CONFIGURE_COMMAND ""
  BUILD_COMMAND     ""
  INSTALL_COMMAND   ${CMAKE_MAKE_PROGRAM} ${LIBCEED_OPTIONS} install
  TEST_COMMAND      ""
)
