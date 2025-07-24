# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Configure Eigen library (header-only)
#

set(EIGEN_OPTIONS ${PALACE_SUPERBUILD_DEFAULT_ARGS})
list(APPEND EIGEN_OPTIONS
  "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
  "-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
  "-DEIGEN_BUILD_DOC=OFF"
  "-DBUILD_TESTING=OFF"
  "-DCMAKE_Fortran_COMPILER=" # Unneeded
)

string(REPLACE ";" "; " EIGEN_OPTIONS_PRINT "${EIGEN_OPTIONS}")
message(STATUS "EIGEN_OPTIONS: ${EIGEN_OPTIONS_PRINT}")

include(ExternalProject)
ExternalProject_Add(eigen
  URL               ${EXTERN_EIGEN_URL}
  SOURCE_DIR        ${CMAKE_BINARY_DIR}/extern/eigen
  BINARY_DIR        ${CMAKE_BINARY_DIR}/extern/eigen-build
  INSTALL_DIR       ${CMAKE_INSTALL_PREFIX}
  PREFIX            ${CMAKE_BINARY_DIR}/extern/eigen-cmake
  UPDATE_COMMAND    ""
  CONFIGURE_COMMAND ${CMAKE_COMMAND} <SOURCE_DIR> "${EIGEN_OPTIONS}"
  TEST_COMMAND      ""
)
