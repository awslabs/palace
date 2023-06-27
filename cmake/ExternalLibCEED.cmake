# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Build libCEED
#

# Force build order
set(LIBCEED_DEPENDENCIES)

# Build LIBXSMM dependency for CPU-based backends
set(PALACE_LIBCEED_WITH_LIBXSMM ON)
if(PALACE_LIBCEED_WITH_LIBXSMM)
  set(LIBXSMM_DEPENDENCIES)

  set(LIBXSMM_OPTIONS
    "PREFIX=${CMAKE_INSTALL_PREFIX}"
    "CC=${CMAKE_C_COMPILER}"
    "CXX=${CMAKE_CXX_COMPILER}"
    "FC=0"
    "FORTRAN=0"
    "BLAS=0"  # For now, no BLAS linkage (like PyFR)
    "SYM=1"   # Always build with symbols
    "VERBOSE=1"
    "PPKGDIR=lib/pkgconfig"
    "PMODDIR=lib/pkgconfig"
  )

  # Always build LIBXSMM as a shared library
  list(APPEND LIBXSMM_OPTIONS
    "STATIC=0"
  )

  # Configure debugging
  if(CMAKE_BUILD_TYPE MATCHES "Debug|debug|DEBUG")
    list(APPEND LIBXSMM_OPTIONS
      "DBG=1"
      "TRACE=1"
    )
  endif()

  # Fix libxsmmext library linkage on macOS
  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    list(APPEND LIBXSMM_OPTIONS
      "LDFLAGS=-undefined dynamic_lookup"
    )
  endif()

  string(REPLACE ";" "; " LIBXSMM_OPTIONS_PRINT "${LIBXSMM_OPTIONS}")
  message(STATUS "LIBXSMM_OPTIONS: ${LIBXSMM_OPTIONS_PRINT}")

  include(ExternalProject)
  ExternalProject_Add(libxsmm
    DEPENDS           ${LIBXSMM_DEPENDENCIES}
    GIT_REPOSITORY    ${EXTERN_LIBXSMM_URL}
    GIT_TAG           ${EXTERN_LIBXSMM_GIT_TAG}
    SOURCE_DIR        ${CMAKE_BINARY_DIR}/extern/libxsmm
    INSTALL_DIR       ${CMAKE_INSTALL_PREFIX}
    PREFIX            ${CMAKE_BINARY_DIR}/extern/libxsmm-cmake
    BUILD_IN_SOURCE   TRUE
    UPDATE_COMMAND    ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND     ""
    INSTALL_COMMAND   ${CMAKE_MAKE_PROGRAM} ${LIBXSMM_OPTIONS} install-minimal
    TEST_COMMAND      ""
  )
  list(APPEND LIBCEED_DEPENDENCIES libxsmm)
endif()

# Note on recommended flags for libCEED (from Makefile, Spack):
#   gcc/clang/icx: -O3 -g -march=native -ffp-contract=fast -fopenmp-simd
#   icc:           -O3 -g -qopenmp-simd
include(CheckCCompilerFlag)
set(LIBCEED_C_FLAGS "${CMAKE_C_FLAGS}")
check_c_compiler_flag(-fopenmp-simd SUPPORTS_OMP_SIMD)
if(SUPPORTS_OMP_SIMD)
  set(LIBCEED_C_FLAGS "${LIBCEED_C_FLAGS} -fopenmp-simd")
endif()

# Build libCEED
set(LIBCEED_OPTIONS
  "prefix=${CMAKE_INSTALL_PREFIX}"
  "CC=${CMAKE_C_COMPILER}"
  "OPT=${LIBCEED_C_FLAGS}"
  "VERBOSE=1"
)

# Always build libCEED as a shared library
list(APPEND LIBCEED_OPTIONS
  "STATIC="
)

# Configure libCEED backends (disable CUDA for now, AVX handled automatically)
list(APPEND LIBCEED_OPTIONS
  "CUDA_DIR=/disable-cuda"
)
if(PALACE_LIBCEED_WITH_LIBXSMM)
  # LIBXSMM requires linkage with BLAS for fallback
  list(APPEND LIBCEED_OPTIONS
    "XSMM_DIR=${CMAKE_INSTALL_PREFIX}"
    "BLAS_LIB="
  )
  # if(NOT "${BLAS_LAPACK_LIBRARIES}" STREQUAL "")
  #   string(REPLACE "$<SEMICOLON>" " " LIBCEED_BLAS_LAPACK_LIBRARIES "${BLAS_LAPACK_LIBRARIES}")
  #   list(APPEND LIBCEED_OPTIONS
  #     "BLAS_LIB=${LIBCEED_BLAS_LAPACK_LIBRARIES}"
  #   )
  # endif()
endif()

string(REPLACE ";" "; " LIBCEED_OPTIONS_PRINT "${LIBCEED_OPTIONS}")
message(STATUS "LIBCEED_OPTIONS: ${LIBCEED_OPTIONS_PRINT}")

# Add H(curl) and H(div) element support to libCEED
set(LIBCEED_PATCH_FILES
  "${CMAKE_SOURCE_DIR}/extern/patch/libCEED/patch_hcurl_hdiv.diff"
  "${CMAKE_SOURCE_DIR}/extern/patch/libCEED/patch_install.diff"
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

if(PALACE_LIBCEED_WITH_LIBXSMM)
  include(GNUInstallDirs)
  set(_LIBCEED_EXTRA_LIBRARIES ${CMAKE_INSTALL_PREFIX}/lib/libxsmm${CMAKE_SHARED_LIBRARY_SUFFIX})
  set(LIBCEED_EXTRA_LIBRARIES ${_LIBCEED_EXTRA_LIBRARIES} CACHE STRING "List of extra library files for libCEED")
endif()
