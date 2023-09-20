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
add_executable(unit-tests "")

# Add source files
target_sources(unit-tests
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-libceed.cpp
  ${CMAKE_SOURCE_DIR}/fem/bilinearform.cpp
  ${CMAKE_SOURCE_DIR}/fem/integrator.cpp
  ${CMAKE_SOURCE_DIR}/fem/integ/curlcurl.cpp
  ${CMAKE_SOURCE_DIR}/fem/integ/curlcurlmass.cpp
  ${CMAKE_SOURCE_DIR}/fem/integ/diffusion.cpp
  ${CMAKE_SOURCE_DIR}/fem/integ/diffusionmass.cpp
  ${CMAKE_SOURCE_DIR}/fem/integ/divdiv.cpp
  ${CMAKE_SOURCE_DIR}/fem/integ/divdivmass.cpp
  ${CMAKE_SOURCE_DIR}/fem/integ/mass.cpp
  ${CMAKE_SOURCE_DIR}/fem/integ/mixedveccurl.cpp
  ${CMAKE_SOURCE_DIR}/fem/integ/mixedvecgrad.cpp
  ${CMAKE_SOURCE_DIR}/fem/integ/vecfemass.cpp
  ${CMAKE_SOURCE_DIR}/fem/libceed/basis.cpp
  ${CMAKE_SOURCE_DIR}/fem/libceed/operator.cpp
  ${CMAKE_SOURCE_DIR}/fem/libceed/restriction.cpp
  ${CMAKE_SOURCE_DIR}/fem/libceed/utils.cpp
)
target_include_directories(unit-tests
  PRIVATE ${CMAKE_SOURCE_DIR}
)

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

# Link with third-party dependencies
if(PALACE_WITH_OPENMP)
  target_link_libraries(unit-tests
    PRIVATE OpenMP::OpenMP_CXX
  )
endif()
target_link_libraries(unit-tests
  PRIVATE mfem ${LIBCEED_TARGET} Catch2::Catch2
)

# Status messages for test target
message(STATUS "Configured unit-tests target for unit tests using Catch2")
