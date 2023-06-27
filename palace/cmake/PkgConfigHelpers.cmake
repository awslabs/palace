# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Helper functions to import and test PETSc and SLEPc targets using pkg-config
#

if(__pkg_config_helpers)
  return()
endif()
set(__pkg_config_helpers YES)

function(set_if_empty _variable _arg)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${_arg} PARENT_SCOPE)
  endif()
endfunction()

set(_PETSC_DIR ${PETSC_DIR})
set(_SLEPC_DIR ${SLEPC_DIR})
set_if_empty(_PETSC_DIR "$ENV{PETSC_DIR}")
set_if_empty(_SLEPC_DIR "$ENV{SLEPC_DIR}")
set_if_empty(_SLEPC_DIR "${_PETSC_DIR}")
if(NOT "${_PETSC_DIR}" STREQUAL "")
  set(ENV{PKG_CONFIG_PATH} "${_PETSC_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
  set(ENV{PKG_CONFIG_PATH} "${_PETSC_DIR}/lib64/pkgconfig:$ENV{PKG_CONFIG_PATH}")
endif()
if(NOT "${_SLEPC_DIR}" STREQUAL "" AND NOT "${SLEPC_DIR}" STREQUAL "${_PETSC_DIR}")
  set(ENV{PKG_CONFIG_PATH} "${_SLEPC_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
  set(ENV{PKG_CONFIG_PATH} "${_SLEPC_DIR}/lib64/pkgconfig:$ENV{PKG_CONFIG_PATH}")
endif()
find_package(PkgConfig REQUIRED)

function(check_petsc_build _petsc_target _petsc_test_success)
  set(PETSC_LIB_TEST_DIR ${CMAKE_BINARY_DIR}/CMakeFiles/try_run)
  set(PETSC_LIB_TEST_CPP ${PETSC_LIB_TEST_DIR}/petsc_test_lib.cpp)
  file(WRITE ${PETSC_LIB_TEST_CPP}
"#include <petsc.h>
int main()
{
  TS ts;
  int argc = 0;
  char **argv = NULL;
  PetscCall(PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL));
  PetscCall(TSCreate(PETSC_COMM_WORLD, &ts));
  PetscCall(TSSetFromOptions(ts));
  PetscCall(TSDestroy(&ts));
  PetscCall(PetscFinalize());
  return 0;
}
")
  try_run(
    PETSC_TEST_EXITCODE
    PETSC_TEST_COMPILED
    ${PETSC_LIB_TEST_DIR}
    ${PETSC_LIB_TEST_CPP}
    LINK_LIBRARIES ${_petsc_target}
    COMPILE_OUTPUT_VARIABLE PETSC_TEST_COMPILE_OUTPUT
    RUN_OUTPUT_VARIABLE PETSC_TEST_OUTPUT
  )
  # message(STATUS "PETSC_TEST_COMPILE_OUTPUT: ${PETSC_TEST_COMPILE_OUTPUT}")
  # message(STATUS "PETSC_TEST_OUTPUT: ${PETSC_TEST_OUTPUT}")
  if(PETSC_TEST_COMPILED AND PETSC_TEST_EXITCODE EQUAL 0)
    # message(STATUS "PETSc test program - Successful")
    set(${_petsc_test_success} TRUE PARENT_SCOPE)
  else()
    # message(STATUS "PETSc test program - Failed")
    set(${_petsc_test_success} FALSE PARENT_SCOPE)
  endif()
endfunction()

function(check_slepc_build _slepc_target _slepc_test_success)
  set(SLEPC_LIB_TEST_DIR ${CMAKE_BINARY_DIR}/CMakeFiles/try_run)
  set(SLEPC_LIB_TEST_CPP ${SLEPC_LIB_TEST_DIR}/slepc_test_lib.cpp)
  file(WRITE ${SLEPC_LIB_TEST_CPP}
"#include <petsc.h>
#include <slepc.h>
int main()
{
  EPS eps;
  int argc = 0;
  char **argv = NULL;
  PetscCall(SlepcInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL));
  PetscCall(EPSCreate(PETSC_COMM_SELF, &eps));
  PetscCall(EPSDestroy(&eps));
  PetscCall(SlepcFinalize());
  return 0;
}
")
  try_run(
    SLEPC_TEST_EXITCODE
    SLEPC_TEST_COMPILED
    ${SLEPC_LIB_TEST_DIR}
    ${SLEPC_LIB_TEST_CPP}
    LINK_LIBRARIES ${_slepc_target}
    COMPILE_OUTPUT_VARIABLE SLEPC_TEST_COMPILE_OUTPUT
    RUN_OUTPUT_VARIABLE SLEPC_TEST_OUTPUT
  )
  # message(STATUS "SLEPC_TEST_COMPILE_OUTPUT: ${SLEPC_TEST_COMPILE_OUTPUT}")
  # message(STATUS "SLEPC_TEST_OUTPUT: ${SLEPC_TEST_OUTPUT}")
  if(SLEPC_TEST_COMPILED AND SLEPC_TEST_EXITCODE EQUAL 0)
    # message(STATUS "SLEPc test program - Successful")
    set(${_slepc_test_success} TRUE PARENT_SCOPE)
  else()
    # message(STATUS "SLEPc test program - Failed")
    set(${_slepc_test_success} FALSE PARENT_SCOPE)
  endif()
endfunction()

function(find_petsc_pkgconfig _petsc_deps _petsc_target)
  pkg_check_modules(PETSc IMPORTED_TARGET GLOBAL PETSc)
  if(NOT PETSc_FOUND)
    pkg_check_modules(PETSc IMPORTED_TARGET GLOBAL petsc)
  endif()
  if(NOT PETSc_FOUND)
    set(${_petsc_target} "" PARENT_SCOPE)
    return()
  endif()
  message(STATUS "Found PETSc: ${PETSc_VERSION}")
  check_petsc_build("PkgConfig::PETSc;${_petsc_deps}" PETSC_TEST_SUCCESS)
  if(PETSC_TEST_SUCCESS)
    message(STATUS "PETSc test program - Successful")
    set(${_petsc_target} PkgConfig::PETSc PARENT_SCOPE)
    return()
  endif()

  # Try with --static libraries
  message(STATUS "PETSc test program - Failed")
  set(PKG_CONFIG_EXECUTABLE_BACKUP ${PKG_CONFIG_EXECUTABLE})
  list(APPEND PKG_CONFIG_EXECUTABLE "--static")
  pkg_check_modules(PETSc_STATIC QUIET IMPORTED_TARGET GLOBAL PETSc)
  if(NOT PETSc_STATIC_FOUND)
    pkg_check_modules(PETSc_STATIC QUIET IMPORTED_TARGET GLOBAL petsc)
  endif()
  set(PKG_CONFIG_EXECUTABLE ${PKG_CONFIG_EXECUTABLE_BACKUP})
  if(NOT PETSc_STATIC_FOUND)
    set(${_petsc_target} "" PARENT_SCOPE)
    return()
  endif()
  check_petsc_build("PkgConfig::PETSc_STATIC;${_petsc_deps}" PETSC_TEST_SUCCESS)
  if(PETSC_TEST_SUCCESS)
    message(STATUS "PETSc test program with static linkage - Success")
    set(${_petsc_target} PkgConfig::PETSc_STATIC PARENT_SCOPE)
    return()
  endif()

  # Not able to build a PETSc test program
  message(STATUS "PETSc test program with static linkage - Failed")
  set(${_petsc_target} "" PARENT_SCOPE)
endfunction()

function(find_slepc_pkgconfig _slepc_deps _slepc_target)
  pkg_check_modules(SLEPc IMPORTED_TARGET GLOBAL SLEPc)
  if(NOT SLEPc_FOUND)
    pkg_check_modules(SLEPc IMPORTED_TARGET GLOBAL slepc)
  endif()
  if(NOT SLEPc_FOUND)
    set(${_slepc_target} "" PARENT_SCOPE)
    return()
  endif()
  message(STATUS "Found SLEPc: ${SLEPc_VERSION}")
  check_slepc_build("PkgConfig::SLEPc;${_slepc_deps}" SLEPC_TEST_SUCCESS)
  if(SLEPC_TEST_SUCCESS)
    message(STATUS "SLEPc test program - Success")
    set(${_slepc_target} PkgConfig::SLEPc PARENT_SCOPE)
    return()
  endif()

  # Try with --static libraries
  message(STATUS "SLEPc test program - Failed")
  set(PKG_CONFIG_EXECUTABLE_BACKUP ${PKG_CONFIG_EXECUTABLE})
  list(APPEND PKG_CONFIG_EXECUTABLE "--static")
  pkg_check_modules(SLEPc_STATIC QUIET IMPORTED_TARGET GLOBAL SLEPc)
  if(NOT SLEPc_STATIC_FOUND)
    pkg_check_modules(SLEPc_STATIC QUIET IMPORTED_TARGET GLOBAL slepc)
  endif()
  set(PKG_CONFIG_EXECUTABLE ${PKG_CONFIG_EXECUTABLE_BACKUP})
  if(NOT SLEPc_STATIC_FOUND)
    set(${_slepc_target} "" PARENT_SCOPE)
    return()
  endif()
  check_slepc_build("PkgConfig::SLEPc_STATIC;${_slepc_deps}" SLEPC_TEST_SUCCESS)
  if(SLEPC_TEST_SUCCESS)
    message(STATUS "SLEPc test program with static linkage - Successful")
    set(${_slepc_target} PkgConfig::SLEPc_STATIC PARENT_SCOPE)
    return()
  endif()

  # Not able to build a SLEPc test program
  message(STATUS "SLEPc test program with static linkage - Failed")
  set(${_slepc_target} "" PARENT_SCOPE)
endfunction()
