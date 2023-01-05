# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Build GSLIB
#

# Force build order
set(GSLIB_DEPENDENCIES petsc)

set(GSLIB_OPTIONS
  "INSTALL_ROOT=${CMAKE_INSTALL_PREFIX}"
  "CC=${CMAKE_C_COMPILER}"
  "MPI=1"
)
if(CMAKE_BUILD_TYPE MATCHES "Debug|debug|DEBUG")
  list(APPEND GSLIB_OPTIONS
    "DEBUG=1"
  )
endif()

set(GSLIB_LDFLAGS)
set(GSLIB_CFLAGS ${CMAKE_C_FLAGS})

# GSLIB is always built as a static library, even when BUILD_SHARED_LIBS=ON, but we need
# -fPIC
if(BUILD_SHARED_LIBS)
  set(GSLIB_CFLAGS "${GSLIB_CFLAGS} -fPIC")
endif()

# User might specify the MPI compiler wrappers directly, otherwise we need to supply MPI
# as found from the CMake module
if(NOT MPI_FOUND)
  message(FATAL_ERROR "MPI is not found when trying to build PETSc")
endif()
if(NOT "${CMAKE_C_COMPILER}" STREQUAL "${MPI_C_COMPILER}")
  string(REPLACE ";" " " GSLIB_MPI_LIBRARIES "${MPI_C_LIBRARIES}")
  set(GSLIB_LDFLAGS "${GSLIB_LDFLAGS} ${GSLIB_MPI_LIBRARIES}")
  foreach(INCLUDE_DIR IN LISTS MPI_C_INCLUDE_DIRS)
    set(GSLIB_CFLAGS "${GSLIB_CFLAGS} -I${INCLUDE_DIR}")
  endforeach()
endif()

# Configure BLAS dependency
if(NOT "${BLAS_LAPACK_LIBRARIES}" STREQUAL "")
  string(REPLACE "$<SEMICOLON>" " " GSLIB_BLAS_LAPACK_LIBRARIES "${BLAS_LAPACK_LIBRARIES}")
  set(GSLIB_LDFLAGS "${GSLIB_LDFLAGS} ${GSLIB_BLAS_LAPACK_LIBRARIES}")
  foreach(INCLUDE_DIR IN LISTS BLAS_LAPACK_INCLUDE_DIRS)
    set(GSLIB_CFLAGS "${GSLIB_CFLAGS} -I${INCLUDE_DIR}")
  endforeach()
  if("${BLA_VENDOR}" MATCHES "Intel")
    list(APPEND GSLIB_OPTIONS
      "BLAS=1"
      "MKL=1"
    )
  else()
    list(APPEND GSLIB_OPTIONS
      "BLAS=1"
    )
  endif()
else()
  list(APPEND GSLIB_OPTIONS
    "BLAS=0"
  )
endif()

list(APPEND GSLIB_OPTIONS
  "CFLAGS=${GSLIB_CFLAGS}"
  "LDFLAGS=${GSLIB_LDFLAGS}"
)

string(REPLACE ";" "; " GSLIB_OPTIONS_PRINT "${GSLIB_OPTIONS}")
message(STATUS "GSLIB_OPTIONS: ${GSLIB_OPTIONS_PRINT}")

# Build fix for MKL
set(GSLIB_PATCH_FILES
  "${CMAKE_CURRENT_SOURCE_DIR}/patch/gslib/patch_mkl.diff"
)

include(ExternalProject)
ExternalProject_Add(gslib
  DEPENDS           ${GSLIB_DEPENDENCIES}
  GIT_REPOSITORY    ${CMAKE_CURRENT_SOURCE_DIR}/gslib
  GIT_TAG           ${EXTERN_GSLIB_GIT_TAG}
  SOURCE_DIR        ${CMAKE_CURRENT_BINARY_DIR}/gslib
  INSTALL_DIR       ${CMAKE_INSTALL_PREFIX}
  PREFIX            ${CMAKE_CURRENT_BINARY_DIR}/gslib-cmake
  BUILD_IN_SOURCE   TRUE
  UPDATE_COMMAND    ""
  PATCH_COMMAND     git apply "${GSLIB_PATCH_FILES}"
  CONFIGURE_COMMAND ""
  BUILD_COMMAND     ""
  INSTALL_COMMAND   ${CMAKE_MAKE_PROGRAM} ${GSLIB_OPTIONS} install
  TEST_COMMAND      ""
)
