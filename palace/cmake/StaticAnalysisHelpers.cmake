# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Helper functions configure static source code analysis with clang-tidy or cppcheck
#

if(__static_analysis_helpers)
  return()
endif()
set(__static_analysis_helpers YES)

function(configure_clang_tidy)
  find_program(CLANG_TIDY_EXE
    NAMES clang-tidy
  )
  if(CLANG_TIDY_EXE)
    # If not explicitly specified, clang-tidy will recurse parent folders and find closest
    # .clang-tidy file (for explicit path specification, add "--config-file=
    # ${CMAKE_SOURCE_DIR}/../.clang-tidy")
    message(STATUS "Found clang-tidy for static analysis: ${CLANG_TIDY_EXE}")
    set(CLANG_TIDY_COMMAND "${CLANG_TIDY_EXE}")

    # Try to extract MPI compiler wrapper include paths from compile command line if not
    # found already (clang-tidy will error about not finding mpi.h otherwise)
    if(MPI_FOUND)
      if(NOT MPI_CXX_INCLUDE_DIRS)
        execute_process(
          COMMAND          ${MPI_CXX_COMPILER} -show
          OUTPUT_VARIABLE  MPI_COMPILE_CMDLINE OUTPUT_STRIP_TRAILING_WHITESPACE
          ERROR_VARIABLE   MPI_COMPILE_CMDLINE ERROR_STRIP_TRAILING_WHITESPACE
          RESULT_VARIABLE  MPI_COMPILER_RETURN
        )
        if(MPI_COMPILER_RETURN EQUAL 0)
          string(REGEX
            MATCHALL "(^| )-I([^\" ]+|\"[^\"]+\")"
            MPI_ALL_INCLUDE_PATHS
            "${MPI_COMPILE_CMDLINE}"
          )
          foreach(IPATH IN LISTS MPI_ALL_INCLUDE_PATHS)
            string(REGEX REPLACE "^ ?-I" "" IPATH ${IPATH})
            string(REGEX REPLACE "//" "/" IPATH ${IPATH})
            list(APPEND MPI_CXX_INCLUDE_DIRS ${IPATH})
          endforeach()
        endif()
      endif()
      if(MPI_CXX_INCLUDE_DIRS)
        set(CLANG_TIDY_EXTRA_ARG)
        foreach(INCLUDE_DIR IN LISTS MPI_CXX_INCLUDE_DIRS)
          set(CLANG_TIDY_EXTRA_ARG "${CLANG_TIDY_EXTRA_ARG} -I${INCLUDE_DIR}")
        endforeach()
        string(STRIP "${CLANG_TIDY_EXTRA_ARG}" CLANG_TIDY_EXTRA_ARG)
        list(APPEND CLANG_TIDY_COMMAND
          "-extra-arg=${CLANG_TIDY_EXTRA_ARG}"
        )
      endif()
    endif()
    set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_COMMAND}" CACHE STRING "" FORCE)
  else()
    message(WARNING "Static analysis with clang-tidy requested, but skipped because the \
executable clang-tidy was not found")
  endif()
endfunction()

function(configure_cppcheck)
  find_program(CPPCHECK_EXE
    NAMES cppcheck
  )
  if(CPPCHECK_EXE)
    message(STATUS "Found cppcheck for static analysis: ${CPPCHECK_EXE}")
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/CMakeFiles/cppcheck)
    execute_process(
      COMMAND     ${CMAKE_COMMAND} -E echo "*:${CMAKE_BINARY_DIR}/_deps/*"
      OUTPUT_FILE ${CMAKE_BINARY_DIR}/CMakeFiles/cppcheck/suppressions.txt
    )
    set(CPPCHECK_COMMAND
      "${CPPCHECK_EXE}"
      "--quiet"
      "--force"
      "--template=gcc"
      "--std=c++17"
      "--enable=warning,style,performance,portability"
      "--suppressions-list=${CMAKE_BINARY_DIR}/CMakeFiles/cppcheck/suppressions.txt"
    )
    set(CMAKE_CXX_CPPCHECK "${CPPCHECK_COMMAND}" CACHE STRING "" FORCE)
  else()
    message(WARNING "Static analysis with cppcheck requested, but skipped because the \
executable cppcheck was not found")
  endif()
endfunction()
