# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Build MFEM
#

# Force build order
if(PALACE_BUILD_EXTERNAL_DEPS)
  set(MFEM_DEPENDENCIES hypre scotch)
  if(PALACE_WITH_GSLIB)
    list(APPEND MFEM_DEPENDENCIES gslib)
  endif()
  if(PALACE_WITH_MUMPS)
    list(APPEND MFEM_DEPENDENCIES mumps)
  endif()
  if(PALACE_WITH_STRUMPACK)
    list(APPEND MFEM_DEPENDENCIES strumpack)
  endif()
  if(PALACE_WITH_SUPERLU)
    list(APPEND MFEM_DEPENDENCIES superlu_dist)
  endif()
else()
  set(MFEM_DEPENDENCIES)
endif()

# Silence #pragma omp warnings when not building with OpenMP
set(MFEM_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
if(PALACE_WITH_STRUMPACK AND NOT PALACE_WITH_OPENMP)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag(-Wno-unknown-pragmas SUPPORTS_NOPRAGMA_WARNING)
  if(SUPPORTS_NOPRAGMA_WARNING)
    set(MFEM_CXX_FLAGS "${MFEM_CXX_FLAGS} -Wno-unknown-pragmas")
  endif()
endif()

# Find optional MFEM dependencies with CMake because once passed to MFEM, they are
# required
set(PALACE_MFEM_WITH_ZLIB NO)
set(PALACE_MFEM_WITH_LIBUNWIND NO)
find_package(ZLIB)
if(ZLIB_FOUND)
  message(STATUS "Building MFEM with zlib support for binary output compression")
  set(PALACE_MFEM_WITH_ZLIB YES)
endif()
if(CMAKE_BUILD_TYPE MATCHES "Debug|debug|DEBUG")
  find_path(LIBUNWIND_INCLUDE_DIR
    NAMES libunwind.h
    HINTS /usr
  )
  if(LIBUNWIND_INCLUDE_DIR)
    message(STATUS "Building MFEM with libunwind support")
    set(PALACE_MFEM_WITH_LIBUNWIND YES)
  endif()
endif()

set(MFEM_OPTIONS ${PALACE_SUPERBUILD_DEFAULT_ARGS})
list(APPEND MFEM_OPTIONS
  "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
  "-DCMAKE_CXX_FLAGS=${MFEM_CXX_FLAGS}"
  "-DMFEM_USE_MPI=YES"
  "-DMFEM_USE_OPENMP=${PALACE_WITH_OPENMP}"
  "-DMFEM_USE_GSLIB=${PALACE_WITH_GSLIB}"
  "-DMFEM_USE_SUPERLU=${PALACE_WITH_SUPERLU}"
  "-DMFEM_USE_STRUMPACK=${PALACE_WITH_STRUMPACK}"
  "-DMFEM_USE_MUMPS=${PALACE_WITH_MUMPS}"
  "-DMFEM_USE_ZLIB=${PALACE_MFEM_WITH_ZLIB}"
  "-DMFEM_USE_LIBUNWIND=${PALACE_MFEM_WITH_LIBUNWIND}"
  "-DMFEM_USE_METIS_5=YES"
)
if(PALACE_WITH_STRUMPACK OR PALACE_WITH_MUMPS OR PALACE_WITH_ARPACK)
  list(APPEND MFEM_OPTIONS
    "-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}"
    "-DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}"
  )
endif()

# Configure MFEM dependencies
if(PALACE_BUILD_EXTERNAL_DEPS)
  list(APPEND MFEM_OPTIONS
    "-DMETIS_LIBRARIES=${METIS_LIBRARIES}"
    "-DMETIS_INCLUDE_DIRS=${METIS_INCLUDE_DIRS}"
    "-DParMETIS_LIBRARIES=${PARMETIS_LIBRARIES}"
    "-DParMETIS_INCLUDE_DIRS=${PARMETIS_INCLUDE_DIRS}"
    "-DHYPRE_DIR=${CMAKE_INSTALL_PREFIX}"
    "-DHYPRE_REQUIRED_LIBRARIES=${BLAS_LAPACK_LIBRARIES}"
  )

  # Configure GSLIB
  if(PALACE_WITH_GSLIB)
    list(APPEND MFEM_OPTIONS
      "-DGSLIB_DIR=${CMAKE_INSTALL_PREFIX}"
    )
  endif()

  # Need to pass gfortran (or similar) dependency to C++ linker for MFEM link line
  if(PALACE_WITH_STRUMPACK OR PALACE_WITH_MUMPS)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(STRUMPACK_MUMPS_GFORTRAN_LIBRARY gfortran)
      else()
        find_library(STRUMPACK_MUMPS_GFORTRAN_LIBRARY
          NAMES gfortran
          PATHS ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}
          NO_DEFAULT_PATH
          REQUIRED
        )
      endif()
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel|IntelLLVM")
      if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Intel|IntelLLVM")
        message(FATAL_ERROR "Intel Fortran compiler detected but not compatible without \
Intel C++ compiler for MUMPS and STRUMPACK dependencies")
      endif()
      set(STRUMPACK_MUMPS_GFORTRAN_LIBRARY ifport$<SEMICOLON>ifcore)
    endif()
  endif()

  # Configure SuperLU_DIST
  if(PALACE_WITH_SUPERLU)
    set(SUPERLU_REQUIRED_PACKAGES "ParMETIS" "METIS" "MPI")
    if(PALACE_WITH_OPENMP)
      list(APPEND SUPERLU_REQUIRED_PACKAGES "OpenMP")
    endif()
    string(REPLACE ";" "$<SEMICOLON>" SUPERLU_REQUIRED_PACKAGES "${SUPERLU_REQUIRED_PACKAGES}")
    list(APPEND MFEM_OPTIONS
      "-DSuperLUDist_DIR=${CMAKE_INSTALL_PREFIX}"
      "-DSuperLUDist_REQUIRED_PACKAGES=${SUPERLU_REQUIRED_PACKAGES}"
      "-DSuperLUDist_REQUIRED_LIBRARIES=${BLAS_LAPACK_LIBRARIES}"
    )
  endif()

  # Configure STRUMPACK
  if(PALACE_WITH_STRUMPACK)
    set(STRUMPACK_REQUIRED_PACKAGES "ParMETIS" "METIS" "MPI" "MPI_Fortran")
    if(PALACE_WITH_OPENMP)
      list(APPEND STRUMPACK_REQUIRED_PACKAGES "OpenMP")
    endif()
    string(REPLACE ";" "$<SEMICOLON>" STRUMPACK_REQUIRED_PACKAGES "${STRUMPACK_REQUIRED_PACKAGES}")
    set(STRUMPACK_REQUIRED_LIBRARIES)
    if(NOT "${STRUMPACK_EXTRA_LIBRARIES}" STREQUAL "")
      list(APPEND STRUMPACK_REQUIRED_LIBRARIES ${STRUMPACK_EXTRA_LIBRARIES})
    endif()
    list(APPEND STRUMPACK_REQUIRED_LIBRARIES ${SCALAPACK_LIBRARIES} ${BLAS_LAPACK_LIBRARIES} ${STRUMPACK_MUMPS_GFORTRAN_LIBRARY})
    string(REPLACE ";" "$<SEMICOLON>" STRUMPACK_REQUIRED_LIBRARIES "${STRUMPACK_REQUIRED_LIBRARIES}")
    list(APPEND MFEM_OPTIONS
      "-DSTRUMPACK_DIR=${CMAKE_INSTALL_PREFIX}"
      "-DSTRUMPACK_REQUIRED_PACKAGES=${STRUMPACK_REQUIRED_PACKAGES}"
      "-DSTRUMPACK_REQUIRED_LIBRARIES=${STRUMPACK_REQUIRED_LIBRARIES}"
    )
  endif()

  # Configure MUMPS
  if(PALACE_WITH_MUMPS)
    set(MUMPS_REQUIRED_PACKAGES "METIS" "MPI" "MPI_Fortran" "Threads")
    if(PALACE_WITH_OPENMP)
      list(APPEND MUMPS_REQUIRED_PACKAGES "OpenMP")
    endif()
    string(REPLACE ";" "$<SEMICOLON>" MUMPS_REQUIRED_PACKAGES "${MUMPS_REQUIRED_PACKAGES}")
    list(APPEND MFEM_OPTIONS
      "-DMUMPS_DIR=${CMAKE_INSTALL_PREFIX}"
      "-DMUMPS_REQUIRED_PACKAGES=${MUMPS_REQUIRED_PACKAGES}"
      "-DMUMPS_REQUIRED_LIBRARIES=${SCALAPACK_LIBRARIES}$<SEMICOLON>${BLAS_LAPACK_LIBRARIES}$<SEMICOLON>${STRUMPACK_MUMPS_GFORTRAN_LIBRARY}"
    )
  endif()
else()
  # Help find dependencies for the internal MFEM build
  # If we trust MFEM's Find<PACKAGE>.cmake module, we can just set <PACKAGE>_DIR and, if
  # needed, <PACKAGE>_REQUIRED_PACKAGES. The extra <PACKAGE>_REQUIRED_LIBRARIES can be used
  # to add any additional dependency libraries.
  set(PALACE_MFEM_DEPS
    "METIS"
    "ParMETIS"
    "HYPRE"
    "GSLIB"
    "SuperLUDist"
    "STRUMPACK"
    "MUMPS"
  )
  foreach(DEP IN LISTS PALACE_MFEM_DEPS)
    set(${DEP}_DIR "" CACHE STRING "Path to ${DEP} build or installation directory")
    set(${DEP}_REQUIRED_PACKAGES "" CACHE STRING "List of additional required packages for ${DEP}")
    set(${DEP}_REQUIRED_LIBRARIES "" CACHE STRING "List of additional required libraries for ${DEP}")
    # set(${DEP}_LIBRARIES "" CACHE STRING "List of library files for ${DEP}")
    # set(${DEP}_INCLUDE_DIRS "" CACHE STRING "Path to ${DEP} include directories")
    if(NOT "${${DEP}_DIR}" STREQUAL "")
      string(REPLACE ";" "$<SEMICOLON>" DEP_DIR "${${DEP}_DIR}")
      list(APPEND MFEM_OPTIONS
        "-D${DEP}_DIR=${DEP_DIR}"
      )
    endif()
    if(NOT "${${DEP}_REQUIRED_PACKAGES}" STREQUAL "")
      string(REPLACE ";" "$<SEMICOLON>" DEP_REQUIRED_PACKAGES "${${DEP}_REQUIRED_PACKAGES}")
      list(APPEND MFEM_OPTIONS
        "-D${DEP}_REQUIRED_PACKAGES=${DEP_REQUIRED_PACKAGES}"
      )
    endif()
    if(NOT "${${DEP}_REQUIRED_LIBRARIES}" STREQUAL "")
      string(REPLACE ";" "$<SEMICOLON>" DEP_REQUIRED_LIBRARIES "${${DEP}_REQUIRED_LIBRARIES}")
      list(APPEND MFEM_OPTIONS
        "-D${DEP}_REQUIRED_LIBRARIES=${DEP_REQUIRED_LIBRARIES}"
      )
    endif()
    # if(NOT "${${DEP}_LIBRARIES}" STREQUAL "")
    # string(REPLACE ";" "$<SEMICOLON>" DEP_LIBRARIES "${${DEP}_LIBRARIES}")
    #   list(APPEND MFEM_OPTIONS
    #     "-D${DEP}_LIBRARIES=${DEP_LIBRARIES}"
    #   )
    # endif()
    # if(NOT "${${DEP}_INCLUDE_DIRS}" STREQUAL "")
    # string(REPLACE ";" "$<SEMICOLON>" DEP_INCLUDE_DIRS "${${DEP}_INCLUDE_DIRS}")
    #   list(APPEND MFEM_OPTIONS
    #     "-D${DEP}_INCLUDE_DIRS=${DEP_INCLUDE_DIRS}"
    #   )
    # endif()
  endforeach()
endif()

string(REPLACE ";" "; " MFEM_OPTIONS_PRINT "${MFEM_OPTIONS}")
message(STATUS "MFEM_OPTIONS: ${MFEM_OPTIONS_PRINT}")

# A number of patches to MFEM for our use cases
set(MFEM_PATCH_FILES
  "${CMAKE_SOURCE_DIR}/extern/patch/mfem/patch_mesh_part.diff"
  "${CMAKE_SOURCE_DIR}/extern/patch/mfem/patch_mesh_vis.diff"
  "${CMAKE_SOURCE_DIR}/extern/patch/mfem/patch_submesh.diff"
  "${CMAKE_SOURCE_DIR}/extern/patch/mfem/patch_direct_solvers.diff"
)

include(ExternalProject)
ExternalProject_Add(mfem
  DEPENDS           ${MFEM_DEPENDENCIES}
  GIT_REPOSITORY    ${EXTERN_MFEM_URL}
  GIT_TAG           ${EXTERN_MFEM_GIT_TAG}
  SOURCE_DIR        ${CMAKE_BINARY_DIR}/extern/mfem
  BINARY_DIR        ${CMAKE_BINARY_DIR}/extern/mfem-build
  INSTALL_DIR       ${CMAKE_INSTALL_PREFIX}
  PREFIX            ${CMAKE_BINARY_DIR}/extern/mfem-cmake
  UPDATE_COMMAND    ""
  PATCH_COMMAND
    git reset --hard && git clean -fd && git apply "${MFEM_PATCH_FILES}"
  CONFIGURE_COMMAND cmake <SOURCE_DIR> "${MFEM_OPTIONS}"
  TEST_COMMAND      ${CMAKE_MAKE_PROGRAM} ex1 ex1p
)
