# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Check compiler support for various features by compiling and running a test program
#

if(__check_compiler_feature_support)
  return()
endif()
set(__check_compiler_feature_support YES)

function(check_constexpr_sqrt_support _has_constexpr_sqrt)
  set(CONSTEXPR_SQRT_TEST_DIR ${CMAKE_BINARY_DIR}/CMakeFiles/try_run)
  set(CONSTEXPR_SQRT_TEST_CPP ${CONSTEXPR_SQRT_TEST_DIR}/constexpr_sqrt_test.cpp)
  file(WRITE ${CONSTEXPR_SQRT_TEST_CPP}
"#include <cmath>
int main()
{
  constexpr double two = 2.0;
  constexpr double four = two*two;
  constexpr double sqrtfour = std::sqrt(four);
  return 0;
}
")
  try_run(
    CONSTEXPR_SQRT_TEST_EXITCODE
    CONSTEXPR_SQRT_TEST_COMPILED
    ${CONSTEXPR_SQRT_TEST_DIR}
    ${CONSTEXPR_SQRT_TEST_CPP}
    COMPILE_OUTPUT_VARIABLE CONSTEXPR_SQRT_TEST_COMPILE_OUTPUT
    RUN_OUTPUT_VARIABLE CONSTEXPR_SQRT_TEST_OUTPUT
  )
  if(CONSTEXPR_SQRT_TEST_COMPILED AND CONSTEXPR_SQRT_TEST_EXITCODE EQUAL 0)
    message(STATUS "CXX compiler supports constexpr std::sqrt")
    set(${_has_constexpr_sqrt} TRUE PARENT_SCOPE)
  else()
    message(STATUS "CXX compiler does not support constexpr std::sqrt")
    set(${_has_constexpr_sqrt} FALSE PARENT_SCOPE)
  endif()
endfunction()

function(check_std_fs_support _has_std_fs_support _extra_fs_libraries)
  set(STD_FS_TEST_DIR ${CMAKE_BINARY_DIR}/CMakeFiles/try_run)
  set(STD_FS_TEST_CPP ${STD_FS_TEST_DIR}/std_fs_test.cpp)
  file(WRITE ${STD_FS_TEST_CPP}
"#include <iostream>
#if defined(__cpp_lib_filesystem) || \
    defined(__has_include) && __has_include(<filesystem>)
#include <filesystem>
namespace fs = std::filesystem;
#else
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif
int main()
{
  std::cout << \"Current path is \" << fs::current_path() << '\\n';
  return 0;
}
")
  try_run(
    STD_FS_TEST_EXITCODE
    STD_FS_TEST_COMPILED
    ${STD_FS_TEST_DIR}
    ${STD_FS_TEST_CPP}
    COMPILE_OUTPUT_VARIABLE STD_FS_TEST_COMPILE_OUTPUT
    RUN_OUTPUT_VARIABLE STD_FS_TEST_OUTPUT
  )
  if(STD_FS_TEST_COMPILED AND STD_FS_TEST_EXITCODE EQUAL 0)
    message(STATUS "CXX compiler supports std::filesystem")
    set(${_has_std_fs_support} TRUE PARENT_SCOPE)
    set(${_extra_fs_libraries} "" PARENT_SCOPE)
    return()
  endif()

  # Try with -lstdc++fs
  try_run(
    STD_FS_TEST_EXITCODE
    STD_FS_TEST_COMPILED
    ${STD_FS_TEST_DIR}
    ${STD_FS_TEST_CPP}
    CMAKE_FLAGS
      "-DLINK_LIBRARIES=stdc++fs"
    COMPILE_OUTPUT_VARIABLE STD_FS_TEST_COMPILE_OUTPUT
    RUN_OUTPUT_VARIABLE STD_FS_TEST_OUTPUT
  )
  if(STD_FS_TEST_COMPILED AND STD_FS_TEST_EXITCODE EQUAL 0)
    message(STATUS "CXX compiler supports std::filesystem with -lstdc++fs")
    set(${_has_std_fs_support} TRUE PARENT_SCOPE)
    set(${_extra_fs_libraries} stdc++fs PARENT_SCOPE)
  else()
    set(${_has_std_fs_support} FALSE PARENT_SCOPE)
    set(${_extra_fs_libraries} "" PARENT_SCOPE)
  endif()
endfunction()
