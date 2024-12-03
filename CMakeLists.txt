# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# CMake superbuild for Palace and its dependencies
#

# CMake 3.21 was released in Jul. 2021 (required for HIP support)
cmake_minimum_required(VERSION 3.21)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Prohibit in-source builds
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR "In-source builds are prohibited")
endif()

# Initialize the project
project(palace-superbuild LANGUAGES CXX C VERSION 0.13.0)

# Define build settings and defaults
set(PALACE_WITH_64BIT_INT OFF CACHE BOOL "Use 64 bit integers")
set(PALACE_WITH_64BIT_BLAS_INT OFF CACHE BOOL "Use ILP64 BLAS/LAPACK interface instead of LP64 (experimental, not recommended)")
set(PALACE_WITH_OPENMP OFF CACHE BOOL "Use OpenMP for shared-memory parallelism")
set(PALACE_WITH_CUDA OFF CACHE BOOL "Use CUDA for NVIDIA GPU support")
set(PALACE_WITH_HIP OFF CACHE BOOL "Use HIP for AMD or NVIDIA GPU support")
set(PALACE_WITH_GPU_AWARE_MPI OFF CACHE BOOL "Option to set if MPI distribution is GPU aware")

set(PALACE_BUILD_EXTERNAL_DEPS ON CACHE BOOL "Build external third-party dependency libraries")
set(PALACE_WITH_SUPERLU ON CACHE BOOL "Build with SuperLU_DIST sparse direct solver")
set(PALACE_WITH_STRUMPACK OFF CACHE BOOL "Build with STRUMPACK sparse direct solver")
set(PALACE_WITH_MUMPS OFF CACHE BOOL "Build with MUMPS sparse direct solver")
set(PALACE_WITH_SLEPC ON CACHE BOOL "Build with SLEPc eigenvalue solver")
set(PALACE_WITH_ARPACK OFF CACHE BOOL "Build with ARPACK eigenvalue solver")
set(PALACE_WITH_LIBXSMM ON CACHE BOOL "Build with LIBXSMM backend for libCEED")
set(PALACE_WITH_MAGMA ON CACHE BOOL "Build with MAGMA backend for libCEED")
set(PALACE_WITH_GSLIB ON CACHE BOOL "Build with GSLIB library for high-order field interpolation")

set(PALACE_WITH_STRUMPACK_BUTTERFLYPACK OFF CACHE BOOL "Build with ButterflyPACK support for STRUMPACK solver")
set(PALACE_WITH_STRUMPACK_ZFP OFF CACHE BOOL "Build with ZFP support for STRUMPACK solver")

set(PALACE_WITH_SUNDIALS ON CACHE BOOL "Build with SUNDIALS differential/algebraic equations solver")

set(ANALYZE_SOURCES_CLANG_TIDY OFF CACHE BOOL "Run static analysis checks using clang-tidy")
set(ANALYZE_SOURCES_CPPCHECK OFF CACHE BOOL "Run static analysis checks using cppcheck")

# Enable Fortran if required
if(PALACE_WITH_STRUMPACK OR PALACE_WITH_MUMPS OR PALACE_WITH_ARPACK)
  enable_language(Fortran)
endif()

# Enable CUDA/HIP if required
if(PALACE_WITH_CUDA AND PALACE_WITH_HIP)
  message(FATAL_ERROR "PALACE_WITH_CUDA is not compatible with PALACE_WITH_HIP")
endif()
if(PALACE_WITH_CUDA)
  # Note: The new behavior of CMake policy CMP0104 will initialize CMAKE_CUDA_ARCHITECTURES
  # to an (old) compatible value even when not set by the user.
  enable_language(CUDA)
  get_filename_component(NVCC_DIR ${CMAKE_CUDA_COMPILER} DIRECTORY)
  get_filename_component(BIN_NVCC_DIR ${NVCC_DIR} DIRECTORY)
  set(CUDA_DIR "${BIN_NVCC_DIR}" CACHE STRING
    "CUDA installation directory, typically /usr/local/cuda or /opt/cuda"
  )
  message(STATUS "Found CUDA: ${CUDA_DIR}")
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    message(STATUS "Setting CMAKE_CUDA_ARCHITECTURES to '70' as none were specified")
    set(CMAKE_CUDA_ARCHITECTURES "70" CACHE STRING
      "Specifies the list of NVIDIA GPU architectures to generate device code for"
    )
  endif()
elseif(PALACE_WITH_HIP)
  enable_language(HIP)
  get_filename_component(HIPCC_DIR ${CMAKE_HIP_COMPILER} DIRECTORY)
  get_filename_component(BIN_HIPCC_DIR ${HIPCC_DIR} DIRECTORY)
  set(ROCM_DIR "${BIN_HIPCC_DIR}" CACHE STRING
    "CUDA installation directory, typically /opt/rocm"
  )
  message(STATUS "Found HIP: ${ROCM_DIR}")
  if(NOT DEFINED CMAKE_HIP_ARCHITECTURES)
    message(STATUS "Setting CMAKE_HIP_ARCHITECTURES to 'gfx900' as none were specified")
    set(CMAKE_HIP_ARCHITECTURES "gfx900" CACHE STRING
      "Specifies the list of AMD GPU architectures to generate device code for"
    )
  endif()
endif()

# Set a default build type if none was provided
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting CMAKE_BUILD_TYPE to 'Release' as none was specified")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Specifies the build type ('Debug' or 'Release', for example)" FORCE
  )
endif()

# Set a default installation location if none was provided
if(NOT CMAKE_INSTALL_PREFIX OR CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  message(STATUS "Setting CMAKE_INSTALL_PREFIX to '${CMAKE_BINARY_DIR}' as none was specified")
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}" CACHE STRING
    "Install directory used by install()" FORCE
  )
endif()

# Set a default for dependency library builds if none was provided
if(NOT DEFINED BUILD_SHARED_LIBS)
  message(STATUS "Setting BUILD_SHARED_LIBS to 'OFF' as it was not specified")
  set(BUILD_SHARED_LIBS OFF CACHE BOOL
    "Global flag to cause add_library() to create shared libraries if ON"
  )
endif()

# Configure default RPATH for installed targets if not provided
if(NOT DEFINED CMAKE_INSTALL_RPATH)
  message(STATUS "Setting CMAKE_INSTALL_RPATH to '\${CMAKE_INSTALL_PREFIX}/lib;\${CMAKE_INSTALL_PREFIX}/lib64' as it was not specified")
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib$<SEMICOLON>${CMAKE_INSTALL_PREFIX}/lib64" CACHE STRING
    "Global RPATH to use for installed targets"
  )
endif()
if(NOT DEFINED CMAKE_INSTALL_RPATH_USE_LINK_PATH)
  message(STATUS "Setting CMAKE_INSTALL_RPATH_USE_LINK_PATH to 'ON' as it was not specified")
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON CACHE BOOL
    "Global flag to append to the RPATH of installed binaries any directories which linked library files if ON"
  )
endif()

# Add extra CMake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Intel classic compilers no longer supported (prior to CMake 3.20, the new compilers are
# detected as Clang instead of IntelLLVM)
if(CMAKE_C_COMPILER_ID STREQUAL "Intel" OR CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  message(WARNING "The Intel classic compilers (icc/icpc) are no longer supported and \
have been replaced by the newer Clang-based icx/icpx from Intel oneAPI")
endif()

# MAGMA is only for GPU builds
if(PALACE_WITH_MAGMA AND NOT (PALACE_WITH_CUDA OR PALACE_WITH_HIP))
  message(STATUS "Disabling MAGMA due to lack of CUDA or HIP support")
  set(PALACE_WITH_MAGMA OFF CACHE BOOL
    "Build with MAGMA backend when libCEED is enabled" FORCE
  )
endif()

# MPI is required for most dependency builds
message(STATUS "====================== Configuring MPI dependency ======================")
find_package(MPI REQUIRED)

# Add BLAS/LAPACK libraries
message(STATUS "================= Configuring BLAS/LAPACK dependencies =================")
include(ExternalBLASLAPACK)

# Default arguments for all external CMake builds (needs to happen after BLAS/LAPACK
# detection in case CMAKE_PREFIX_PATH is modified). Compilers and flags are added on a
# package by package basis.
# See Spack docs: http://tinyurl.com/43t2recp
set(PALACE_SUPERBUILD_DEFAULT_ARGS
  "-DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}"
  "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
  "-DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}"
  "-DCMAKE_EXE_LINKER_FLAGS=${CMAKE_EXE_LINKER_FLAGS}"
  "-DCMAKE_MODULE_LINKER_FLAGS=${CMAKE_MODULE_LINKER_FLAGS}"
  "-DCMAKE_SHARED_LINKER_FLAGS=${CMAKE_SHARED_LINKER_FLAGS}"
  "-DCMAKE_STATIC_LINKER_FLAGS=${CMAKE_STATIC_LINKER_FLAGS}"
)
if(NOT "${CMAKE_PREFIX_PATH}" STREQUAL "")
  string(REPLACE ";" "$<SEMICOLON>" PALACE_CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}")
  list(APPEND PALACE_SUPERBUILD_DEFAULT_ARGS
    "-DCMAKE_PREFIX_PATH=${PALACE_CMAKE_PREFIX_PATH}$<SEMICOLON>${CMAKE_INSTALL_PREFIX}"
  )
else()
  list(APPEND PALACE_SUPERBUILD_DEFAULT_ARGS
    "-DCMAKE_PREFIX_PATH=${CMAKE_INSTALL_PREFIX}"
  )
endif()
if(NOT "${CMAKE_INSTALL_RPATH}" STREQUAL "")
  string(REPLACE ";" "$<SEMICOLON>" PALACE_CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH}")
  list(APPEND PALACE_SUPERBUILD_DEFAULT_ARGS
    "-DCMAKE_INSTALL_RPATH=${PALACE_CMAKE_INSTALL_RPATH}"
  )
endif()
if(NOT "${CMAKE_INSTALL_RPATH_USE_LINK_PATH}" STREQUAL "")
  list(APPEND PALACE_SUPERBUILD_DEFAULT_ARGS
    "-DCMAKE_INSTALL_RPATH_USE_LINK_PATH=${CMAKE_INSTALL_RPATH_USE_LINK_PATH}"
  )
endif()

include(ExternalGitTags)
# Avoid DOWNLOAD_EXTRACT_TIMESTAMP warning
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)
endif()

# Add other third-party dependency builds
include(ExternalGitTags)
if(PALACE_BUILD_EXTERNAL_DEPS)
  add_subdirectory("extern")
endif()

# Add GSLIB (always built as part of Palace)
if(PALACE_WITH_GSLIB)
  message(STATUS "===================== Configuring GSLIB dependency =====================")
  include(ExternalGSLIB)
endif()

# Add libCEED (always built as part of Palace)
message(STATUS "==================== Configuring libCEED dependency ====================")
include(ExternalLibCEED)

# Add MFEM (always built as part of Palace)
message(STATUS "====================== Configuring MFEM dependency =====================")
include(ExternalMFEM)

# Add the main Palace project
message(STATUS "========================== Configuring Palace ==========================")
include(ExternalPalace)

# Finished with superbuild configuration
message(STATUS "======================= Configure stage complete =======================")
