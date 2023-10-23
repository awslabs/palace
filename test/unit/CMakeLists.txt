# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# CMake configuration for Palace unit tests
#

# Add Catch2
Include(FetchContent)
FetchContent_Declare(Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.4.0
)
FetchContent_MakeAvailable(Catch2)

# Add executable target
add_executable(unit-tests
  ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-libceed.cpp
)
target_link_libraries(unit-tests PRIVATE ${LIB_TARGET_NAME} Catch2::Catch2)

# Add JIT source file path definition for libCEED
set_property(
  SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
  APPEND PROPERTY COMPILE_DEFINITIONS "PALACE_LIBCEED_JIT_SOURCE_DIR=\"${CMAKE_INSTALL_PREFIX}/include/palace\""
)

# Add unit test mesh file path definition
set_property(
  SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/test-libceed.cpp
  APPEND PROPERTY COMPILE_DEFINITIONS "PALACE_TEST_MESH_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/mesh\""
)

# Status messages for test target
message(STATUS "Configured unit-tests target for unit tests using Catch2")
