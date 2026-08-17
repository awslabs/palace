# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Configure Enzyme compiler plugin
#

if("${LLVM_DIR}" STREQUAL "")
  message(FATAL_ERROR "PALACE_WITH_ENZYME requires LLVM_DIR")
endif()
if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  message(FATAL_ERROR
    "PALACE_WITH_ENZYME requires upstream Clang built against the same LLVM as Enzyme "
    "(got ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION})"
  )
endif()

find_package(LLVM REQUIRED CONFIG HINTS "${LLVM_DIR}")
if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL LLVM_PACKAGE_VERSION)
  message(FATAL_ERROR
    "Enzyme requires matching Clang and LLVM versions: compiler is "
    "${CMAKE_CXX_COMPILER_VERSION}, LLVM is ${LLVM_PACKAGE_VERSION}"
  )
endif()

set(ENZYME_OPTIONS ${PALACE_SUPERBUILD_DEFAULT_ARGS})
list(APPEND ENZYME_OPTIONS
  "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
  "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
  "-DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}"
  "-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
  "-DLLVM_DIR=${LLVM_DIR}"
  "-DENZYME_ENABLE_PLUGINS=ON"
  "-DENZYME_CLANG=ON"
  "-DENZYME_EXTERNAL_SHARED_LIB=OFF"
  "-DENZYME_STATIC_LIB=OFF"
)

string(REPLACE ";" "; " ENZYME_OPTIONS_PRINT "${ENZYME_OPTIONS}")
message(STATUS "ENZYME_OPTIONS: ${ENZYME_OPTIONS_PRINT}")

include(ExternalProject)
ExternalProject_Add(enzyme
  GIT_REPOSITORY    ${EXTERN_ENZYME_URL}
  GIT_TAG           ${EXTERN_ENZYME_GIT_TAG}
  SOURCE_DIR        ${CMAKE_BINARY_DIR}/extern/enzyme-src
  BINARY_DIR        ${CMAKE_BINARY_DIR}/extern/enzyme-build
  INSTALL_DIR       ${CMAKE_INSTALL_PREFIX}
  PREFIX            ${CMAKE_BINARY_DIR}/extern/enzyme-superbuild
  UPDATE_COMMAND    ""
  CONFIGURE_COMMAND ${CMAKE_COMMAND} <SOURCE_DIR>/enzyme "${ENZYME_OPTIONS}"
  TEST_COMMAND      ""
)
