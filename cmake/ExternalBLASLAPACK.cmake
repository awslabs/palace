# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Configure BLAS/LAPACK libraries
#

# Configuring 64-bit BLAS/LAPACK integer interface is not supported by older CMake.
if(NOT PALACE_WITH_64BIT_INT AND PALACE_WITH_64BIT_BLAS_INT)
  message(FATAL_ERROR "ILP64 BLAS/LAPACK interface requires PALACE_WITH_64BIT_INT")
endif()
if(NOT ${CMAKE_VERSION} VERSION_LESS "3.22.0")
  if(PALACE_WITH_64BIT_BLAS_INT)
    set(BLA_SIZEOF_INTEGER 8)
  else()
    set(BLA_SIZEOF_INTEGER 4)
  endif()
endif()

# Defines cache variables BLAS_LAPACK_LIBRARIES/BLAS_LAPACK_INCLUDE_DIRS for building
# dependencies on top of BLAS and LAPACK
if(DEFINED ENV{ARMPL_DIR} OR DEFINED ENV{ARMPLROOT} OR DEFINED ENV{ARMPL_ROOT})
  # Arm Performance Libraries for arm64 builds when available
  if(NOT CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm")
    message(WARNING "Arm PL math libraries are only intended for arm64 architecture builds")
  endif()
  if(DEFINED ENV{ARMPL_DIR})
    set(ARMPL_DIR $ENV{ARMPL_DIR})
  elseif(DEFINED ENV{ARMPLROOT})
    set(ARMPL_DIR $ENV{ARMPLROOT})
  elseif(DEFINED ENV{ARMPL_ROOT})
    set(ARMPL_DIR $ENV{ARMPL_ROOT})
  else()
    set(ARMPL_DIR)
  endif()
  if(PALACE_WITH_64BIT_BLAS_INT)
    if(PALACE_WITH_OPENMP)
      set(ARMPL_LIB_SUFFIX "_ilp64_mp")
    else()
      set(ARMPL_LIB_SUFFIX "_ilp64")
    endif()
  else()
    if(PALACE_WITH_OPENMP)
      set(ARMPL_LIB_SUFFIX "_mp")
    else()
      set(ARMPL_LIB_SUFFIX "")
    endif()
  endif()
  list(APPEND CMAKE_PREFIX_PATH ${ARMPL_DIR})
  set(BLA_VENDOR "Arm${ARMPL_LIB_SUFFIX}")
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)

  # Locate include directory
  find_path(_BLAS_LAPACK_INCLUDE_DIRS
    NAMES armpl.h
    PATHS ${ARMPL_DIR}
    PATH_SUFFIXES include
    NO_DEFAULT_PATH
    REQUIRED
  )
  message(STATUS "Using BLAS/LAPACK from Arm Performance Libraries (Arm PL)")
elseif(DEFINED ENV{AOCL_DIR} OR DEFINED ENV{AOCLROOT} OR DEFINED ENV{AOCL_ROOT})
  # AOCL for x86_64 builds when available (part of CMake's FindBLAS/FindLAPACK as of v3.27
  # but unnecessarily adds -fopenmp flag)
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm")
    message(WARNING "AOCL math libraries are not intended for arm64 architecture builds")
  endif()
  if(DEFINED ENV{AOCL_DIR})
    set(AOCL_DIR $ENV{AOCL_DIR})
  elseif(DEFINED ENV{AOCLROOT})
    set(AOCL_DIR $ENV{AOCLROOT})
  elseif(DEFINED ENV{AOCL_ROOT})
    set(AOCL_DIR $ENV{AOCL_ROOT})
  else()
    set(AOCL_DIR)
  endif()
  if(PALACE_WITH_64BIT_BLAS_INT)
    set(AOCL_DIR_SUFFIX "_ILP64")
  else()
    set(AOCL_DIR_SUFFIX "_LP64")
  endif()
  if(PALACE_WITH_OPENMP)
    set(AOCL_LIB_SUFFIX "-mt")
  else()
    set(AOCL_LIB_SUFFIX "")
  endif()
  find_library(BLAS_LIBRARIES
    NAMES blis${AOCL_LIB_SUFFIX} blis
    PATHS ${AOCL_DIR}
    PATH_SUFFIXES lib${AOCL_DIR_SUFFIX} lib lib64
    NO_DEFAULT_PATH
    REQUIRED
  )
  find_library(LAPACK_LIBRARIES
    NAMES flame FLAME
    PATHS ${AOCL_DIR}
    PATH_SUFFIXES lib${AOCL_DIR_SUFFIX} lib lib64
    NO_DEFAULT_PATH
    REQUIRED
  )

  # Locate include directory
  find_path(_BLAS_LAPACK_INCLUDE_DIRS
    NAMES cblas.h
    PATHS ${AOCL_DIR}
    PATH_SUFFIXES include${AOCL_DIR_SUFFIX} include/blis include
    NO_DEFAULT_PATH
    REQUIRED
  )
  message(STATUS "Using BLAS/LAPACK from AMD BLIS/libFLAME")
elseif(DEFINED ENV{MKL_DIR} OR DEFINED ENV{MKLROOT} OR DEFINED ENV{MKL_ROOT})
  # MKL for x86_64 builds when available
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm")
    message(WARNING "MKL math libraries are not intended for arm64 architecture builds")
  endif()
  if(DEFINED ENV{MKL_DIR})
    set(MKL_DIR $ENV{MKL_DIR})
  elseif(DEFINED ENV{MKLROOT})
    set(MKL_DIR $ENV{MKLROOT})
  elseif(DEFINED ENV{MKL_ROOT})
    set(MKL_DIR $ENV{MKL_ROOT})
  else()
    set(MKL_DIR)
  endif()
  if(PALACE_WITH_64BIT_BLAS_INT)
    if(PALACE_WITH_OPENMP)
      set(MKL_LIB_SUFFIX "_64ilp")
    else()
      set(MKL_LIB_SUFFIX "_64ilp_seq")
    endif()
  else()
    if(PALACE_WITH_OPENMP)
      set(MKL_LIB_SUFFIX "_64lp")
    else()
      set(MKL_LIB_SUFFIX "_64lp_seq")
    endif()
  endif()
  list(APPEND CMAKE_PREFIX_PATH ${MKL_DIR})
  set(BLA_VENDOR "Intel10${MKL_LIB_SUFFIX}")
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)

  # Locate include directories
  find_path(_BLAS_LAPACK_INCLUDE_DIRS
    NAMES mkl_cblas.h
    PATHS ${MKL_DIR}
    PATH_SUFFIXES include
    NO_DEFAULT_PATH
    REQUIRED
  )
  message(STATUS "Using BLAS/LAPACK from Intel MKL")
else()
  # Try to find OpenBLAS installation on the system
  # Warning: This does NOT automatically configure for OpenMP support
  if(DEFINED ENV{OPENBLAS_DIR})
    set(OPENBLAS_DIR $ENV{OPENBLAS_DIR})
  elseif(DEFINED ENV{OPENBLASROOT})
    set(OPENBLAS_DIR $ENV{OPENBLASROOT})
  elseif(DEFINED ENV{OPENBLAS_ROOT})
    set(OPENBLAS_DIR $ENV{OPENBLAS_ROOT})
  else()
    set(OPENBLAS_DIR)
    message(STATUS "Using BLAS/LAPACK located by CMake")
  endif()

  if(NOT OPENBLAS_DIR STREQUAL "")
    # If OpenBLAS was found set the vendor to avoid conflict with Accelerate on Darwin
    set(BLA_VENDOR "OpenBLAS")
    message(STATUS "Using BLAS/LAPACK from OpenBLAS")
  endif()

  list(APPEND CMAKE_PREFIX_PATH ${OPENBLAS_DIR})
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)

  # Locate include directory
  set(_BLAS_LAPACK_DIRS)
  foreach(LIB IN LISTS LAPACK_LIBRARIES BLAS_LIBRARIES)
    cmake_path(GET LIB PARENT_PATH LIB_DIR)
    cmake_path(GET LIB_DIR PARENT_PATH LIB_DIR)
    list(APPEND _BLAS_LAPACK_DIRS ${LIB_DIR})
  endforeach()
  list(REMOVE_DUPLICATES _BLAS_LAPACK_DIRS)
  find_path(_BLAS_LAPACK_INCLUDE_DIRS
    NAMES cblas.h
    HINTS ${_BLAS_LAPACK_DIRS}
    PATH_SUFFIXES include include/openblas include/blis
    REQUIRED
  )
endif()
set(LAPACK_LIBRARIES "${LAPACK_LIBRARIES};-lm")

# Save variables to cache
set(_BLAS_LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
list(REMOVE_DUPLICATES _BLAS_LAPACK_LIBRARIES)
string(REPLACE ";" "$<SEMICOLON>" _BLAS_LAPACK_LIBRARIES "${_BLAS_LAPACK_LIBRARIES}")
set(BLAS_LAPACK_LIBRARIES ${_BLAS_LAPACK_LIBRARIES} CACHE STRING
  "List of library files for BLAS/LAPACK"
)
set(BLAS_LAPACK_INCLUDE_DIRS ${_BLAS_LAPACK_INCLUDE_DIRS} CACHE STRING
  "Path to BLAS/LAPACK include directories"
)