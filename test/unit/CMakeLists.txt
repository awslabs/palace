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
  ${CMAKE_CURRENT_SOURCE_DIR}/test-tablecsv.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-config.cpp
)
target_link_libraries(unit-tests PRIVATE ${LIB_TARGET_NAME} Catch2::Catch2)

# Handle device source code
set(TARGET_SOURCES_DEVICE
  ${CMAKE_SOURCE_DIR}/fem/libceed/operator.cpp
)
if(PALACE_WITH_CUDA OR PALACE_WITH_HIP)
  if(PALACE_WITH_CUDA)
    set(LANGUAGE_PROPERTY CUDA)
  elseif(PALACE_WITH_HIP)
    set(LANGUAGE_PROPERTY HIP)
  endif()
  set(COMPILE_OPTIONS_PROPERTY "-Wno-pedantic")
  if(PALACE_WITH_OPENMP)
    set(COMPILE_OPTIONS_PROPERTY "${COMPILE_OPTIONS_PROPERTY} ${OpenMP_CXX_FLAGS}")
  endif()
  set_property(
    SOURCE ${TARGET_SOURCES_DEVICE}
    PROPERTY LANGUAGE ${LANGUAGE_PROPERTY}
  )
  set_property(
    SOURCE ${TARGET_SOURCES_DEVICE}
    APPEND PROPERTY COMPILE_OPTIONS "${COMPILE_OPTIONS_PROPERTY}"
  )
endif()

# Resolve linker error with CUDA and static linkage
if(PALACE_WITH_CUDA)
  set_target_properties(unit-tests PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
endif()

# Add JIT source file path definition for libCEED
set_property(
  SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
  APPEND PROPERTY COMPILE_DEFINITIONS "PALACE_LIBCEED_JIT_SOURCE_DIR=\"${CMAKE_INSTALL_PREFIX}/include/palace/\""
)

# Add unit test mesh file path definition
set_property(
  SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/test-libceed.cpp
  APPEND PROPERTY COMPILE_DEFINITIONS "PALACE_TEST_MESH_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/mesh\""
)

target_compile_definitions(
  unit-tests
  PUBLIC "PALACE_TEST_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}\""
)

# Status messages for test target
message(STATUS "Configured unit-tests target for unit tests using Catch2")
