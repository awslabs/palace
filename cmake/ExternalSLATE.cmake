# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Build SLATE
#

# Force build order
set(SLATE_DEPENDENCIES)

# Silence #pragma omp warnings when not building with OpenMP
set(SLATE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
if(PALACE_WITH_STRUMPACK AND NOT PALACE_WITH_OPENMP)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag(-Wno-unknown-pragmas SUPPORTS_NOPRAGMA_WARNING)
  if(SUPPORTS_NOPRAGMA_WARNING)
    set(SLATE_CXX_FLAGS "${SLATE_CXX_FLAGS} -Wno-unknown-pragmas")
  endif()
endif()

set(SLATE_OPTIONS ${PALACE_SUPERBUILD_DEFAULT_ARGS})
list(APPEND SLATE_OPTIONS
  "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
  "-DCMAKE_CXX_FLAGS=${SLATE_CXX_FLAGS}"
  "-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}"
  "-DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}"
  "-Duse_mpi=TRUE"
  "-Duse_openmp=${PALACE_WITH_OPENMP}"
  "-Dc_api=FALSE"
  "-Dslate_is_project=FALSE"
  "-Dbuild_tests=FALSE"
)

# Configure BLAS/LAPACK
if(NOT "${BLAS_LAPACK_LIBRARIES}" STREQUAL "")
  list(APPEND SLATE_OPTIONS
    "-DLAPACK_LIBRARIES=${BLAS_LAPACK_LIBRARIES}"
    "-DBLAS_LIBRARIES=${BLAS_LAPACK_LIBRARIES}"
    "-Duse_cmake_find_lapack=FALSE"
    "-Duse_cmake_find_blas=FALSE"
  )
endif()

# Configure GPU support
if(PALACE_WITH_CUDA)
  list(APPEND SLATE_OPTIONS
    "-Dgpu_backend=cuda"
    "-DCMAKE_CUDA_COMPILER=${CMAKE_CUDA_COMPILER}"
    "-DCMAKE_CUDA_FLAGS=${CMAKE_CUDA_FLAGS}"
  )
  if(NOT "${CMAKE_CUDA_ARCHITECTURES}" STREQUAL "")
    list(APPEND SLATE_OPTIONS
      "-DCMAKE_CUDA_ARCHITECTURES=${CMAKE_CUDA_ARCHITECTURES}"
    )
  endif()
endif()
if(PALACE_WITH_HIP)
  list(APPEND SLATE_OPTIONS
    "-Dgpu_backend=hip"
    "-DCMAKE_HIP_COMPILER=${CMAKE_HIP_COMPILER}"
    "-DCMAKE_HIP_FLAGS=${CMAKE_HIP_FLAGS}"
  )
  if(NOT "${CMAKE_HIP_ARCHITECTURES}" STREQUAL "")
    list(APPEND SLATE_OPTIONS
      "-DCMAKE_HIP_ARCHITECTURES=${CMAKE_HIP_ARCHITECTURES}"
    )
  endif()
endif()

string(REPLACE ";" "; " SLATE_OPTIONS_PRINT "${SLATE_OPTIONS}")
message(STATUS "SLATE_OPTIONS: ${SLATE_OPTIONS_PRINT}")

include(ExternalProject)
ExternalProject_Add(slate
  DEPENDS                ${SLATE_DEPENDENCIES}
  GIT_REPOSITORY         ${EXTERN_SLATE_URL}
  GIT_TAG                ${EXTERN_SLATE_GIT_TAG}
  SOURCE_DIR             ${CMAKE_BINARY_DIR}/extern/slate
  BINARY_DIR             ${CMAKE_BINARY_DIR}/extern/slate-build
  INSTALL_DIR            ${CMAKE_INSTALL_PREFIX}
  PREFIX                 ${CMAKE_BINARY_DIR}/extern/slate-cmake
  GIT_SUBMODULES_RECURSE TRUE
  UPDATE_COMMAND         ""
  CONFIGURE_COMMAND      ${CMAKE_COMMAND} <SOURCE_DIR> "${SLATE_OPTIONS}"
  TEST_COMMAND           ""
)
