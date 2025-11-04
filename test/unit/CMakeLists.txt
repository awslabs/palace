# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# CMake configuration for Palace unit tests
#

# Add Catch2
if(PALACE_BUILD_EXTERNAL_DEPS)
  Include(FetchContent)
  FetchContent_Declare(Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.8.1
  )
  FetchContent_MakeAvailable(Catch2)
else()
  find_package(Catch2 3 QUIET)
  if(NOT Catch2_FOUND)
    message(WARNING "Catch2 not found, unit-tests target will not be available")
    return()
  endif()
endif()

# Add executable target
add_executable(unit-tests
  ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-config.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-constants.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-geodata.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-libceed.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-materialoperator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-postoperator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-postoperatorcsv.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-rap.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-strattonchu.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-tablecsv.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-vector.cpp
)

# Set output name for installed binary
set_target_properties(unit-tests PROPERTIES OUTPUT_NAME palace-unit-tests)
# Below we link using WHOLE_ARCHIVE with coverage, otherwise LLVM won't report coverage of
# all palace files. The actual linker commands to do this by hand is delicate. (see
# https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_LINK_LIBRARY_USING_FEATURE.html#loading-a-whole-static-library)
target_link_libraries(unit-tests 
  PRIVATE Catch2::Catch2
          $<$<BOOL:${PALACE_BUILD_WITH_COVERAGE}>:coverage_flags>
          $<IF:$<BOOL:${PALACE_BUILD_WITH_COVERAGE}>,$<LINK_LIBRARY:WHOLE_ARCHIVE,${LIB_TARGET_NAME}>,${LIB_TARGET_NAME}>
          )

# Handle device source code
set(TARGET_SOURCES_DEVICE
  ${CMAKE_SOURCE_DIR}/fem/libceed/operator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-vector.cpp
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
         ${CMAKE_CURRENT_SOURCE_DIR}/test-materialoperator.cpp
         ${CMAKE_CURRENT_SOURCE_DIR}/test-strattonchu.cpp
  APPEND PROPERTY COMPILE_DEFINITIONS "PALACE_TEST_MESH_DIR=\"${CMAKE_INSTALL_PREFIX}/share/palace/test/mesh\""
)
set_property(
  SOURCE
  ${CMAKE_CURRENT_SOURCE_DIR}/test-geodata.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-rap.cpp
  APPEND PROPERTY COMPILE_DEFINITIONS "MFEM_DATA_PATH=\"${CMAKE_INSTALL_PREFIX}/share/palace/test/mfem-data\""
)

target_compile_definitions(
  unit-tests
  PUBLIC "PALACE_TEST_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}\""
  PRIVATE $<TARGET_PROPERTY:${LIB_TARGET_NAME},COMPILE_DEFINITIONS>
)

# Status messages for test target
message(STATUS "Configured unit-tests target for unit tests using Catch2")

# Install Catch2 libraries
if(PALACE_BUILD_EXTERNAL_DEPS)
  install(
    TARGETS Catch2
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
  )
endif()

# Install unit-tests executable to bin directory
install(
  TARGETS unit-tests
  RUNTIME DESTINATION bin
)

# Install mesh files to share directory
install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/mesh
  DESTINATION share/palace/test
)

# Install MFEM data files if available
if(MFEM_DATA_PATH)
  install(
    DIRECTORY ${MFEM_DATA_PATH}/
    DESTINATION share/palace/test/mfem-data
  )
endif()

install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/mesh/
  DESTINATION share/palace/test/mesh
)
