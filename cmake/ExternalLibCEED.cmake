# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Build libCEED
#

# Force build order
set(LIBCEED_DEPENDENCIES)

set(LIBCEED_CFLAGS ${CMAKE_C_FLAGS})

# Build LIBXSMM dependency for CPU-based backends (header-only)
set(PALACE_LIBCEED_WITH_LIBXSMM ON)
if(PALACE_LIBCEED_WITH_LIBXSMM)
  set(LIBXSMM_DEPENDENCIES)

  # Configure debugging in LIBXSMM
  if(CMAKE_BUILD_TYPE MATCHES "Debug|debug|DEBUG")
    set(LIBCEED_CFLAGS "${LIBCEED_CFLAGS} -D_DEBUG -D__TRACE=1")
  else()
    set(LIBCEED_CFLAGS "${LIBCEED_CFLAGS} -DNDEBUG")
  endif()

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
    INSTALL_COMMAND
      ${CMAKE_COMMAND} -E rm -rf <INSTALL_DIR>/include/libxsmm &&
      ${CMAKE_COMMAND} -E copy_directory <SOURCE_DIR>/src <INSTALL_DIR>/include/libxsmm &&
      ${CMAKE_COMMAND} -E echo "file(GLOB LIBXSMM_INCLUDE_FILES \"<SOURCE_DIR>/include/libxsmm*.h\")" > <SOURCE_DIR>/install-includes.cmake &&
      ${CMAKE_COMMAND} -E echo "file(INSTALL \${LIBXSMM_INCLUDE_FILES} DESTINATION <INSTALL_DIR>/include)" >> <SOURCE_DIR>/install-includes.cmake &&
      ${CMAKE_COMMAND} -P <SOURCE_DIR>/install-includes.cmake &&
      <SOURCE_DIR>/scripts/libxsmm_source.sh libxsmm > <INSTALL_DIR>/include/libxsmm_source.h
    TEST_COMMAND      ""
  )
  list(APPEND LIBCEED_DEPENDENCIES libxsmm)
endif()

# Build libCEED (libCEED's Makefile will add some flags, but we pass CMake's flags just to
# make sure they are included)
set(LIBCEED_OPTIONS
  "prefix=${CMAKE_INSTALL_PREFIX}"
  "CC=${CMAKE_C_COMPILER}"
  "PEDANTIC=1"
  "PEDANTICFLAGS=${LIBCEED_CFLAGS}"
  "VERBOSE=1"
)


#XX TODO EXTRACT USEFUL ONES FROM FROM LIBCEED_CFLAGS AND RELY ON LIBCEED MAKEFILE?
#XX TODO OPT FLAGS LIKE SPACK? OR JUST RELY ON MAKEFILE (like -qopenmp-simd)


# Always build libCEED as a shared library
list(APPEND LIBCEED_OPTIONS
  "STATIC="
)
# if(BUILD_SHARED_LIBS)
#   list(APPEND LIBCEED_OPTIONS
#     "STATIC="
#   )
# else()
#   list(APPEND LIBCEED_OPTIONS
#     "STATIC=1"
#   )
# endif()

# Configure libCEED backends (disable CUDA for now, AVX handled automatically)
list(APPEND LIBCEED_OPTIONS
  "CUDA_DIR=/disable-cuda"
)
if(PALACE_LIBCEED_WITH_LIBXSMM)
  # LIBXSMM requires linkage with BLAS for fallback
  list(APPEND LIBCEED_OPTIONS
    "XSMM_DIR=${CMAKE_INSTALL_PREFIX}"
  )
  if(NOT "${BLAS_LAPACK_LIBRARIES}" STREQUAL "")
    string(REPLACE "$<SEMICOLON>" " " LIBCEED_BLAS_LAPACK_LIBRARIES "${BLAS_LAPACK_LIBRARIES}")
    list(APPEND LIBCEED_OPTIONS
      "BLAS_LIB=${LIBCEED_BLAS_LAPACK_LIBRARIES}"
    )
  endif()
endif()

string(REPLACE ";" "; " LIBCEED_OPTIONS_PRINT "${LIBCEED_OPTIONS}")
message(STATUS "LIBCEED_OPTIONS: ${LIBCEED_OPTIONS_PRINT}")

# Add H(curl) and H(div) element support to libCEED
set(LIBCEED_PATCH_FILES
  "${CMAKE_SOURCE_DIR}/extern/patch/libCEED/patch_hcurl_hdiv.diff"
  "${CMAKE_SOURCE_DIR}/extern/patch/libCEED/patch_install.diff"
  "${CMAKE_SOURCE_DIR}/extern/patch/libCEED/patch_xsmm.diff"
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
