# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Build dependencies with source managed by Git submodules
#

# # Update Git submodules
# find_package(Git QUIET)
# if(GIT_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/.git")
#   option(GIT_SUBMODULE "Check submodules during build" ON)
#   if(GIT_SUBMODULE)
#     message(STATUS "Updating Git submodules")
#     execute_process(
#       COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
#       WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
#       RESULT_VARIABLE GIT_SUBMOD_RESULT
#     )
#     if(NOT GIT_SUBMOD_RESULT EQUAL "0")
#       message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
#     endif()
#   endif()
# endif()

# Git tags for third-party dependency builds
include(ExternalGitTags)

# Add METIS/ParMETIS libraries
message(STATUS "================ Configuring METIS/ParMETIS dependencies ===============")
if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/METIS OR
   NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ParMETIS)
  message(FATAL_ERROR
    "Could not locate source directories for METIS or ParMETIS, maybe try `git submodule update --init`"
  )
endif()
include(ExternalMETIS)

# Add SuperLU_DIST
if(PALACE_WITH_SUPERLU)
  message(STATUS "================= Configuring SuperLU_DIST dependency ==================")
  if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/superlu_dist)
    message(FATAL_ERROR
      "Could not locate source directory for SuperLU_DIST, maybe try `git submodule update --init`"
    )
  endif()
  include(ExternalSuperLU_DIST)
endif()

# Add ScaLAPACK
if(PALACE_WITH_STRUMPACK OR PALACE_WITH_MUMPS)
  message(STATUS "=================== Configuring ScaLAPACK dependency ===================")
  if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/scalapack)
    message(FATAL_ERROR
      "Could not locate source directory for ScaLAPACK, maybe try `git submodule update --init`"
    )
  endif()
  include(ExternalScaLAPACK)
endif()

# Add STRUMPACK
if(PALACE_WITH_STRUMPACK)
  message(STATUS "=================== Configuring STRUMPACK dependency ===================")
  if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/STRUMPACK)
    message(FATAL_ERROR
      "Could not locate source directory for STRUMPACK, maybe try `git submodule update --init`"
    )
  endif()
  include(ExternalSTRUMPACK)
endif()

# Add MUMPS
if(PALACE_WITH_MUMPS)
  message(STATUS "===================== Configuring MUMPS dependency =====================")
  if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/mumps)
    message(FATAL_ERROR
      "Could not locate source directory for MUMPS, maybe try `git submodule update --init`"
    )
  endif()
  include(ExternalMUMPS)
endif()

# Add ARPACK
if(PALACE_WITH_ARPACK)
  message(STATUS "==================== Configuring ARPACK dependency =====================")
  if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/arpack-ng)
    message(FATAL_ERROR
      "Could not locate source directory for ARPACK, maybe try `git submodule update --init`"
    )
  endif()
  include(ExternalARPACK)
endif()

# Add SLEPc (requires PETSc)
if(PALACE_WITH_SLEPC)
  message(STATUS "=============== Configuring PETSc and SLEPc dependencies ===============")
  if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/petsc OR
     NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/slepc)
    message(FATAL_ERROR
      "Could not locate source directory for PETSc or SLEPc, maybe try `git submodule update --init`"
    )
  endif()
  include(ExternalSLEPc)
endif()

# Add GSLIB
if(PALACE_WITH_GSLIB)
  message(STATUS "===================== Configuring GSLIB dependency =====================")
  if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/gslib)
    message(FATAL_ERROR
      "Could not locate source directories for GSLIB, maybe try `git submodule update --init`"
    )
  endif()
  include(ExternalGSLIB)
endif()

# Add HYPRE
message(STATUS "===================== Configuring HYPRE dependency =====================")
if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/hypre)
  message(FATAL_ERROR
    "Could not locate source directory for HYPRE, maybe try `git submodule update --init`"
  )
endif()
include(ExternalHYPRE)
