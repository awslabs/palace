# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Configure BLAS/LAPACK libraries
#

if(NOT PALACE_WITH_64BIT_INT AND PALACE_WITH_64BIT_BLAS_INT)
  message(FATAL_ERROR "ILP64 BLAS/LAPACK interface requires PALACE_WITH_64BIT_INT")
endif()
if(PALACE_WITH_64BIT_BLAS_INT)
  set(BLA_SIZEOF_INTEGER 8)
else()
  set(BLA_SIZEOF_INTEGER 4)
endif()

# Installation roots of the supported vendors, from the environment
function(palace_blas_lapack_env_dir _out)
  foreach(_var IN LISTS ARGN)
    if(DEFINED ENV{${_var}})
      set(${_out} $ENV{${_var}} PARENT_SCOPE)
      return()
    endif()
  endforeach()
  set(${_out} "" PARENT_SCOPE)
endfunction()
palace_blas_lapack_env_dir(ARMPL_DIR ARMPL_DIR ARMPLROOT ARMPL_ROOT)
palace_blas_lapack_env_dir(AOCL_DIR AOCL_DIR AOCLROOT AOCL_ROOT)
palace_blas_lapack_env_dir(MKL_DIR MKL_DIR MKLROOT MKL_ROOT)
palace_blas_lapack_env_dir(OPENBLAS_DIR OPENBLAS_DIR OPENBLASROOT OPENBLAS_ROOT)

# A user-provided BLA_VENDOR takes precedence, otherwise pick the vendor from the
# environment
if(DEFINED BLA_VENDOR)
  message(STATUS "Using BLAS/LAPACK from BLA_VENDOR=${BLA_VENDOR}")
elseif(NOT ARMPL_DIR STREQUAL "")
  if(PALACE_WITH_64BIT_BLAS_INT)
    set(ARMPL_LIB_SUFFIX "_ilp64")
  else()
    set(ARMPL_LIB_SUFFIX "")
  endif()
  if(PALACE_WITH_OPENMP)
    set(ARMPL_LIB_SUFFIX "${ARMPL_LIB_SUFFIX}_mp")
  endif()
  set(BLA_VENDOR "Arm${ARMPL_LIB_SUFFIX}")
  message(STATUS "Using BLAS/LAPACK from Arm Performance Libraries (Arm PL)")
elseif(NOT AOCL_DIR STREQUAL "")
  if(PALACE_WITH_OPENMP)
    set(BLA_VENDOR "AOCL_mt")
  else()
    set(BLA_VENDOR "AOCL")
  endif()
  message(STATUS "Using BLAS/LAPACK from AMD BLIS/libFLAME")
elseif(NOT MKL_DIR STREQUAL "")
  if(PALACE_WITH_64BIT_BLAS_INT)
    set(MKL_LIB_SUFFIX "_64ilp")
  else()
    set(MKL_LIB_SUFFIX "_64lp")
  endif()
  if(NOT PALACE_WITH_OPENMP)
    set(MKL_LIB_SUFFIX "${MKL_LIB_SUFFIX}_seq")
  endif()
  set(BLA_VENDOR "Intel10${MKL_LIB_SUFFIX}")
  message(STATUS "Using BLAS/LAPACK from Intel MKL")
elseif(NOT OPENBLAS_DIR STREQUAL "")
  # Warning: This does NOT automatically configure for OpenMP support
  # Setting the vendor avoids a conflict with Accelerate on Darwin
  set(BLA_VENDOR "OpenBLAS")
  message(STATUS "Using BLAS/LAPACK from OpenBLAS")
else()
  message(STATUS "Using BLAS/LAPACK located by CMake")
endif()

# Vendor specific search paths and include file
set(_BLAS_LAPACK_ROOT)
set(_BLAS_LAPACK_HEADER cblas.h)
set(_BLAS_LAPACK_INCLUDE_SUFFIXES include include/openblas include/blis)
if(BLA_VENDOR MATCHES "^Arm")
  if(NOT CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm")
    message(WARNING "Arm PL math libraries are only intended for arm64 architecture builds")
  endif()
  set(_BLAS_LAPACK_ROOT ${ARMPL_DIR})
  set(_BLAS_LAPACK_HEADER armpl.h)
elseif(BLA_VENDOR MATCHES "^AOCL")
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm")
    message(WARNING "AOCL math libraries are not intended for arm64 architecture builds")
  endif()
  if(PALACE_WITH_64BIT_BLAS_INT)
    set(AOCL_DIR_SUFFIX "_ILP64")
  else()
    set(AOCL_DIR_SUFFIX "_LP64")
  endif()
  set(_BLAS_LAPACK_ROOT ${AOCL_DIR})
  if(NOT AOCL_DIR STREQUAL "")
    list(APPEND CMAKE_LIBRARY_PATH ${AOCL_DIR}/lib${AOCL_DIR_SUFFIX})
  endif()
  list(PREPEND _BLAS_LAPACK_INCLUDE_SUFFIXES include${AOCL_DIR_SUFFIX})
elseif(BLA_VENDOR MATCHES "^Intel")
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm")
    message(WARNING "MKL math libraries are not intended for arm64 architecture builds")
  endif()
  set(_BLAS_LAPACK_ROOT ${MKL_DIR})
  set(_BLAS_LAPACK_HEADER mkl_cblas.h)
elseif(BLA_VENDOR STREQUAL "OpenBLAS")
  set(_BLAS_LAPACK_ROOT ${OPENBLAS_DIR})
endif()

list(APPEND CMAKE_PREFIX_PATH ${_BLAS_LAPACK_ROOT})
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# FindLAPACK links AOCL libFLAME with -fopenmp, which is not needed
list(REMOVE_ITEM LAPACK_LIBRARIES "-fopenmp")

# Locate include directory
set(_BLAS_LAPACK_DIRS ${_BLAS_LAPACK_ROOT})
foreach(LIB IN LISTS LAPACK_LIBRARIES BLAS_LIBRARIES)
  if(IS_ABSOLUTE "${LIB}")
    cmake_path(GET LIB PARENT_PATH LIB_DIR)
    cmake_path(GET LIB_DIR PARENT_PATH LIB_DIR)
    list(APPEND _BLAS_LAPACK_DIRS ${LIB_DIR})
  endif()
endforeach()
list(REMOVE_DUPLICATES _BLAS_LAPACK_DIRS)
# On macOS frameworks are searched first by default, which finds the Accelerate headers
# for any vendor
set(_CMAKE_FIND_FRAMEWORK ${CMAKE_FIND_FRAMEWORK})
if("${LAPACK_LIBRARIES};${BLAS_LIBRARIES}" MATCHES "\\.framework")
  set(CMAKE_FIND_FRAMEWORK FIRST)
else()
  set(CMAKE_FIND_FRAMEWORK LAST)
endif()
find_path(_BLAS_LAPACK_INCLUDE_DIRS
  NAMES ${_BLAS_LAPACK_HEADER}
  HINTS ${_BLAS_LAPACK_DIRS}
  PATH_SUFFIXES ${_BLAS_LAPACK_INCLUDE_SUFFIXES}
  REQUIRED
)
set(CMAKE_FIND_FRAMEWORK ${_CMAKE_FIND_FRAMEWORK})
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
