# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Build SuperLU_DIST
#

# Force build order
set(SUPERLU_DEPENDENCIES parmetis)

set(SUPERLU_OPTIONS ${PALACE_SUPERBUILD_DEFAULT_ARGS})
list(APPEND SUPERLU_OPTIONS
  "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
  "-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
  "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
  "-DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}"
  "-DXSDK_ENABLE_Fortran=OFF"
  "-Denable_tests=OFF"
  "-Denable_examples=OFF"
  "-Denable_double=ON"
  "-Denable_single=ON"
  "-Denable_complex16=ON"
  "-Denable_python=OFF"
  "-DTPL_ENABLE_PARMETISLIB=ON"
  "-DTPL_PARMETIS_LIBRARIES=${PARMETIS_LIBRARIES}$<SEMICOLON>${METIS_LIBRARIES}"
  "-DTPL_PARMETIS_INCLUDE_DIRS=${CMAKE_INSTALL_PREFIX}/include"
  "-DTPL_ENABLE_COMBBLASLIB=OFF"
)

# Always disable OpenMP (seems slower in all cases, just link to threaded BLAS/LAPACK)
# if(PALACE_WITH_OPENMP)
#   list(APPEND SUPERLU_OPTIONS
#     "-Denable_openmp=ON"
#   )
# else()
  list(APPEND SUPERLU_OPTIONS
    "-Denable_openmp=OFF"
  )
# endif()

# SuperLU_DIST has a BUILD_STATIC_LIBS option which defaults to ON
if(BUILD_SHARED_LIBS)
  list(APPEND SUPERLU_OPTIONS
    "-DBUILD_STATIC_LIBS=OFF"
  )
endif()

# Configure 64-bit indices
if(PALACE_WITH_64BIT_INT)
  list(APPEND SUPERLU_OPTIONS
    "-DXSDK_INDEX_SIZE=64"
  )
endif()

# Configure LAPACK dependency
if(NOT "${BLAS_LAPACK_LIBRARIES}" STREQUAL "")
  list(APPEND SUPERLU_OPTIONS
    "-DTPL_ENABLE_LAPACKLIB=ON"
    "-DTPL_ENABLE_INTERNAL_BLASLIB=OFF"
    "-DLAPACK_LIBRARIES=${BLAS_LAPACK_LIBRARIES}"
    "-DBLAS_LIBRARIES=${BLAS_LAPACK_LIBRARIES}"
  )
endif()

# Configure GPU support (for now, disable since faster on CPU)
list(APPEND SUPERLU_OPTIONS
  "-DTPL_ENABLE_CUDALIB=OFF"
  "-DTPL_ENABLE_HIPLIB=OFF"
)
# if(PALACE_WITH_CUDA)
#   list(APPEND SUPERLU_OPTIONS
#     "-DTPL_ENABLE_CUDALIB=ON"
#     "-DCMAKE_CUDA_COMPILER=${CMAKE_CUDA_COMPILER}"
#     "-DCMAKE_CUDA_FLAGS=${CMAKE_CUDA_FLAGS}"
#   )
#   if(NOT "${CMAKE_CUDA_ARCHITECTURES}" STREQUAL "")
#     list(APPEND SUPERLU_OPTIONS
#       "-DCMAKE_CUDA_ARCHITECTURES=${CMAKE_CUDA_ARCHITECTURES}"
#     )
#   endif()
# else()
#   list(APPEND SUPERLU_OPTIONS
#     "-DTPL_ENABLE_CUDALIB=OFF"
#   )
# endif()
# if(PALACE_WITH_HIP)
#   list(APPEND SUPERLU_OPTIONS
#     "-DTPL_ENABLE_HIPLIB=ON"
#     "-DHIP_ROOT_DIR=${ROCM_DIR}"
#     "-DCMAKE_HIP_COMPILER=${CMAKE_HIP_COMPILER}"
#     "-DCMAKE_HIP_FLAGS=${CMAKE_HIP_FLAGS}"
#   )
#   if(NOT "${CMAKE_HIP_ARCHITECTURES}" STREQUAL "")
#     list(APPEND SUPERLU_OPTIONS
#       "-DCMAKE_HIP_ARCHITECTURES=${CMAKE_HIP_ARCHITECTURES}"
#     )
#   endif()
# else()
#   list(APPEND SUPERLU_OPTIONS
#     "-DTPL_ENABLE_HIPLIB=OFF"
#   )
# endif()

string(REPLACE ";" "; " SUPERLU_OPTIONS_PRINT "${SUPERLU_OPTIONS}")
message(STATUS "SUPERLU_OPTIONS: ${SUPERLU_OPTIONS_PRINT}")

# Fix column permutations
set(SUPERLU_PATCH_FILES
  "${CMAKE_SOURCE_DIR}/extern/patch/superlu_dist/patch_parmetis.diff"
)

include(ExternalProject)
ExternalProject_Add(superlu_dist
  DEPENDS           ${SUPERLU_DEPENDENCIES}
  GIT_REPOSITORY    ${EXTERN_SUPERLU_URL}
  GIT_TAG           ${EXTERN_SUPERLU_GIT_TAG}
  SOURCE_DIR        ${CMAKE_BINARY_DIR}/extern/superlu_dist
  BINARY_DIR        ${CMAKE_BINARY_DIR}/extern/superlu_dist-build
  INSTALL_DIR       ${CMAKE_INSTALL_PREFIX}
  PREFIX            ${CMAKE_BINARY_DIR}/extern/superlu_dist-cmake
  UPDATE_COMMAND    ""
  PATCH_COMMAND     git apply "${SUPERLU_PATCH_FILES}"
  CONFIGURE_COMMAND ${CMAKE_COMMAND} <SOURCE_DIR> "${SUPERLU_OPTIONS}"
  TEST_COMMAND      ""
)
