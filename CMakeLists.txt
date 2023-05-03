# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# CMake superbuild for Palace and its dependencies
#

# CMake 3.13 was released in Nov. 2018
cmake_minimum_required(VERSION 3.13)

# Prohibit in-source builds
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are prohibited")
endif()

# Initialize the project
project(palace-superbuild LANGUAGES CXX C VERSION 0.11.0)

# Define build settings and defaults
set(PALACE_WITH_64BIT_INT OFF CACHE BOOL "Use 64 bit integers")
set(PALACE_WITH_64BIT_BLAS_INT OFF CACHE BOOL "Use ILP64 BLAS/LAPACK interface instead of LP64 (not recommended)")
set(PALACE_WITH_OPENMP OFF CACHE BOOL "Use OpenMP")
set(PALACE_WITH_GSLIB ON CACHE BOOL "Build with GSLIB library for high-order field interpolation")
set(PALACE_WITH_SUPERLU ON CACHE BOOL "Build with SuperLU_DIST sparse direct solver")
set(PALACE_WITH_STRUMPACK OFF CACHE BOOL "Build with STRUMPACK sparse direct solver")
set(PALACE_WITH_MUMPS OFF CACHE BOOL "Build with MUMPS sparse direct solver")
set(PALACE_WITH_SLEPC ON CACHE BOOL "Build with SLEPc eigenvalue solver")
set(PALACE_WITH_ARPACK OFF CACHE BOOL "Build with ARPACK eigenvalue solver")
set(ANALYZE_SOURCES_CLANG_TIDY OFF CACHE BOOL "Run static analysis checks using clang-tidy")
set(ANALYZE_SOURCES_CPPCHECK OFF CACHE BOOL "Run static analysis checks using cppcheck")

# Enable Fortran if required
if(PALACE_WITH_STRUMPACK OR PALACE_WITH_MUMPS OR PALACE_WITH_ARPACK)
  enable_language(Fortran)
endif()

# Set a default installation location if none was provided
if(NOT CMAKE_INSTALL_PREFIX OR CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  message(STATUS "Setting CMAKE_INSTALL_PREFIX to '${CMAKE_BINARY_DIR}' as none was specified")
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}" CACHE STRING "Install directory used by install()" FORCE)
endif()

# Set a default build type if none was provided
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting CMAKE_BUILD_TYPE to 'Release' as none was specified")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Specifies the build type ('Debug' or 'Release', for example)" FORCE)
endif()

# Set a default for dependency library builds if none was provided
if(NOT DEFINED BUILD_SHARED_LIBS)
  message(STATUS "Setting BUILD_SHARED_LIBS to 'OFF' as it was not specified")
  set(BUILD_SHARED_LIBS OFF CACHE BOOL "Global flag to cause add_library() to create shared libraries if ON")
endif()

# Add extra CMake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/extern")

# Intel classic compilers no longer supported (prior to CMake 3.20, the new compilers are
# detected as Clang instead of IntelLLVM)
if(CMAKE_C_COMPILER_ID STREQUAL "Intel" OR CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  message(WARNING "The Intel classic compilers (icc/icpc) are no longer supported and \
have been replaced by the newer Clang-based icx/icpx from Intel oneAPI")
endif()

# Default arguments for all external CMake builds
set(PALACE_SUPERBUILD_DEFAULT_ARGS
  "-DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}"
  "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
  "-DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}"
)

# MPI is required for most dependency builds
message(STATUS "====================== Configuring MPI dependency ======================")
find_package(MPI REQUIRED)

# Add BLAS/LAPACK libraries
message(STATUS "================= Configuring BLAS/LAPACK dependencies =================")
include(ExternalBLASLAPACK)

# Add other third-party dependency builds
add_subdirectory("extern")

# Add the main Palace project
message(STATUS "========================== Configuring Palace ==========================")
include(ExternalPalace)

# Finished with superbuild configuration
message(STATUS "======================= Configure stage complete =======================")
