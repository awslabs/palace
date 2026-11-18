# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Build Scotch and PT-Scotch
#

# Force build order
set(SCOTCH_DEPENDENCIES)

set(SCOTCH_OPTIONS ${PALACE_SUPERBUILD_DEFAULT_ARGS})
list(APPEND SCOTCH_OPTIONS
  "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
  "-DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}"
  "-DBUILD_PTSCOTCH=ON"
  "-DBUILD_FORTRAN=OFF"
  "-DBUILD_LIBESMUMPS=OFF"
  "-DBUILD_LIBSCOTCHMETIS=OFF"  # Its METIS API would clash with METIS
  "-DINSTALL_METIS_HEADERS=OFF"
  "-DENABLE_TESTS=OFF"
  "-DUSE_ZLIB=OFF"
  "-DUSE_BZ2=OFF"
  "-DUSE_LZMA=OFF"
  "-DTHREADS=OFF"
  "-DMPI_THREAD_MULTIPLE=OFF"
  "-DSCOTCH_DETERMINISTIC=FULL"
  "-DIDXSIZE=64"
)

# Configure 64-bit indices
if(PALACE_WITH_64BIT_INT)
  list(APPEND SCOTCH_OPTIONS "-DINTSIZE=64")
else()
  list(APPEND SCOTCH_OPTIONS "-DINTSIZE=32")
endif()

string(REPLACE ";" "; " SCOTCH_OPTIONS_PRINT "${SCOTCH_OPTIONS}")
message(STATUS "SCOTCH_OPTIONS: ${SCOTCH_OPTIONS_PRINT}")

include(ExternalProject)
ExternalProject_Add(scotch
  DEPENDS           ${SCOTCH_DEPENDENCIES}
  GIT_REPOSITORY    ${EXTERN_SCOTCH_URL}
  GIT_TAG           ${EXTERN_SCOTCH_GIT_TAG}
  SOURCE_DIR        ${CMAKE_BINARY_DIR}/extern/scotch
  BINARY_DIR        ${CMAKE_BINARY_DIR}/extern/scotch-build
  INSTALL_DIR       ${CMAKE_INSTALL_PREFIX}
  PREFIX            ${CMAKE_BINARY_DIR}/extern/scotch-cmake
  UPDATE_COMMAND    ""
  CONFIGURE_COMMAND ${CMAKE_COMMAND} <SOURCE_DIR> "${SCOTCH_OPTIONS}"
  # Install the libraries, without the Scotch programs and manual pages
  INSTALL_COMMAND
    ${CMAKE_COMMAND} --install <BINARY_DIR> --component libscotch
    COMMAND ${CMAKE_COMMAND} --install <BINARY_DIR> --component libptscotch
  TEST_COMMAND      ""
)
