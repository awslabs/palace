# Copyright 2009-2013, Iowa State University
# Copyright 2013-2020, Ryan Pavlik
# Copyright 2013-2020, Contributors
# SPDX-License-Identifier: BSL-1.0
# Distributed under the Boost Software License, Version 1.0
# See copy at http://www.boost.org/LICENSE_1_0.txt
# SPDX-License-Identifier: BSL-1.0

#
# Returns the refspec and sha hash of the current head revision of the results of git
# describe on the source tree. These functions force a re-configure on each git commit so
# that you can trust the values of the variables in your build system.
#
# Original author: 2009-2020 Ryan Pavlik <ryan.pavlik@gmail.com> <abiryan@ryand.net>
#

if(__get_git_description)
  return()
endif()
set(__get_git_description YES)

set(CURRENT_LIST_DIR ${CMAKE_CURRENT_LIST_DIR})

function(git_find_closest_git_dir _start_dir _git_dir_var)
  set(cur_dir "${_start_dir}")
  set(git_dir "${_start_dir}/.git")
  while(NOT EXISTS "${git_dir}")
    # .git dir not found, search parent directories
    set(git_previous_parent "${cur_dir}")
    get_filename_component(cur_dir ${cur_dir} DIRECTORY)
    if(cur_dir STREQUAL git_previous_parent)
      # We have reached the root directory, we are not in git
      set(${_git_dir_var} "" PARENT_SCOPE)
      return()
    endif()
    set(git_dir "${cur_dir}/.git")
  endwhile()
  set(${_git_dir_var} "${git_dir}" PARENT_SCOPE)
endfunction()

function(get_git_head_revision _refspecvar _hashvar)
  git_find_closest_git_dir("${CMAKE_CURRENT_SOURCE_DIR}" GIT_DIR)
  if("${GIT_DIR}" STREQUAL "" OR NOT IS_DIRECTORY ${GIT_DIR})
    set(${_refspecvar} "GITDIR-NOTFOUND" PARENT_SCOPE)
    set(${_hashvar} "GITDIR-NOTFOUND" PARENT_SCOPE)
    return()
  endif()

  set(GIT_DATA "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/git-data")
  if(NOT EXISTS "${GIT_DATA}")
    file(MAKE_DIRECTORY "${GIT_DATA}")
  endif()
  set(HEAD_FILE "${GIT_DATA}/HEAD")
  set(HEAD_SOURCE_FILE "${GIT_DIR}/HEAD")
  if(NOT EXISTS "${HEAD_SOURCE_FILE}")
    set(${_refspecvar} "GITHEAD-NOTFOUND" PARENT_SCOPE)
    set(${_hashvar} "GITHEAD-NOTFOUND" PARENT_SCOPE)
    return()
  endif()

  configure_file("${HEAD_SOURCE_FILE}" "${HEAD_FILE}" COPYONLY)
  configure_file("${CURRENT_LIST_DIR}/GetGitDescription.cmake.in"
                 "${GIT_DATA}/GetGitRefHash.cmake" @ONLY)
  include("${GIT_DATA}/GetGitRefHash.cmake")
  set(${_refspecvar} "${HEAD_REF}" PARENT_SCOPE)
  set(${_hashvar} "${HEAD_HASH}" PARENT_SCOPE)
endfunction()

function(git_describe _var)
  if(NOT GIT_FOUND)
    find_package(Git QUIET)
  endif()
  if(NOT GIT_FOUND)
    set(${_var} "GIT-NOTFOUND" PARENT_SCOPE)
    return()
  endif()

  # Force rerun only if Git status has changed
  get_git_head_revision(refspec hash)
  if(NOT hash)
    set(${_var} "HEAD-HASH-NOTFOUND" PARENT_SCOPE)
    return()
  endif()

  # message(STATUS "Git head revision: ${refspec} ${hash}")
  # message(STATUS "Arguments to execute_process: ${ARGN}")

  execute_process(
    COMMAND "${GIT_EXECUTABLE}" describe --tags --always --dirty ${ARGN}
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    RESULT_VARIABLE res
    OUTPUT_VARIABLE out
    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if(NOT res EQUAL 0)
    set(${_var} "${res}-NOTFOUND" PARENT_SCOPE)
  else()
    set(${_var} "${out}" PARENT_SCOPE)
  endif()
endfunction()
