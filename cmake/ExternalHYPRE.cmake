# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Build HYPRE
#

# Force build order
set(HYPRE_DEPENDENCIES)
if(PALACE_WITH_UMPIRE)
  list(APPEND HYPRE_DEPENDENCIES umpire)
endif()

# Silence CUDA/HIP include warnings with -pedantic
set(HYPRE_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${BUILD_TYPE_UPPER}}")
set(HYPRE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${BUILD_TYPE_UPPER}}")
if(HYPRE_C_FLAGS MATCHES "-pedantic" AND (PALACE_WITH_CUDA OR PALACE_WITH_HIP))
  string(REGEX REPLACE "-pedantic" "" HYPRE_C_FLAGS "${HYPRE_C_FLAGS}")
  string(REGEX REPLACE "-pedantic" "" HYPRE_CXX_FLAGS "${HYPRE_CXX_FLAGS}")
endif()

set(HYPRE_OPTIONS
  -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
  -DCMAKE_C_FLAGS=${HYPRE_C_FLAGS}
  -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
  -DHYPRE_ENABLE_FORTRAN=OFF
  -DHYPRE_ENABLE_MPI=ON
)
if(BUILD_SHARED_LIBS)
  list(APPEND HYPRE_OPTIONS -DBUILD_SHARED_LIBS=ON)
else()
  list(APPEND HYPRE_OPTIONS -DBUILD_SHARED_LIBS=OFF)
endif()
if(PALACE_WITH_64BIT_INT)
  list(APPEND HYPRE_OPTIONS
    -DHYPRE_ENABLE_MIXEDINT=ON
    -DHYPRE_ENABLE_BIGINT=OFF
  )
else()
  list(APPEND HYPRE_OPTIONS
    -DHYPRE_ENABLE_MIXEDINT=OFF
    -DHYPRE_ENABLE_BIGINT=OFF
  )
endif()
if(PALACE_WITH_OPENMP)
  list(APPEND HYPRE_OPTIONS -DHYPRE_ENABLE_OPENMP=ON)
endif()

# Configure Umpire support
if(PALACE_WITH_UMPIRE)
  list(APPEND HYPRE_OPTIONS
    -DHYPRE_ENABLE_UMPIRE=ON
    -Dumpire_DIR=${CMAKE_INSTALL_PREFIX}/lib/cmake/umpire
  )
endif()

# Configure BLAS/LAPACK
if(NOT "${BLAS_LAPACK_LIBRARIES}" STREQUAL "")
  list(APPEND HYPRE_OPTIONS
    -DHYPRE_ENABLE_HYPRE_BLAS=OFF
    -DHYPRE_ENABLE_HYPRE_LAPACK=OFF
    -DTPL_BLAS_LIBRARIES=${BLAS_LAPACK_LIBRARIES}
    -DTPL_LAPACK_LIBRARIES=${BLAS_LAPACK_LIBRARIES}
  )
endif()

# Configure GPU support
if(PALACE_WITH_CUDA OR PALACE_WITH_HIP)
  # Palace enables GPU-aware MPI at runtime when MPI supports it, this only sets hypre's
  # default (and lets MFEM compile its GPU-aware hypre communication)
  list(APPEND HYPRE_OPTIONS
    -DHYPRE_ENABLE_UNIFIED_MEMORY=OFF
    -DHYPRE_ENABLE_GPU_AWARE_MPI=ON
  )
endif()
if(PALACE_WITH_CUDA)
  list(APPEND HYPRE_OPTIONS
    -DHYPRE_ENABLE_CUDA=ON
    # https://github.com/hypre-space/hypre/issues/1404
    -DCUDAToolkit_ROOT=${CUDAToolkit_LIBRARY_ROOT}
  )
  palace_append_cuda_architectures(HYPRE_OPTIONS)
endif()
if(PALACE_WITH_HIP)
  list(APPEND HYPRE_OPTIONS
    -DHYPRE_ENABLE_HIP=ON
  )
  if(NOT "${CMAKE_HIP_ARCHITECTURES}" STREQUAL "")
    list(APPEND HYPRE_OPTIONS -DCMAKE_HIP_ARCHITECTURES=${CMAKE_HIP_ARCHITECTURES})
  endif()
endif()

string(REPLACE ";" "; " HYPRE_OPTIONS_PRINT "${HYPRE_OPTIONS}")
message(STATUS "HYPRE_OPTIONS: ${HYPRE_OPTIONS_PRINT}")

# Sync compute stream before GPU-aware MPI in device matvec, using the patch from the Spack
# build (https://github.com/hypre-space/hypre/pull/1629,
# https://github.com/spack/spack-packages/pull/7198)
include(PatchHelpers)
set(HYPRE_PATCH_FILES)
download_patch(HYPRE_PATCH_FILES "${CMAKE_BINARY_DIR}/extern/hypre-patches"
  hypre_pr1629.patch
  "https://github.com/hypre-space/hypre/commit/5d0282f44490f406ea532ea6778fe914df677346.patch?full_index=1"
  8b11d35c29215d462b8e727775e8771d808806ced107e1675082ed80b9ba84fb
)

include(ExternalProject)
ExternalProject_Add(hypre
  DEPENDS           ${HYPRE_DEPENDENCIES}
  GIT_REPOSITORY    ${EXTERN_HYPRE_URL}
  GIT_TAG           ${EXTERN_HYPRE_GIT_TAG}
  SOURCE_DIR        ${CMAKE_BINARY_DIR}/extern/hypre
  BINARY_DIR        ${CMAKE_BINARY_DIR}/extern/hypre-build
  INSTALL_DIR       ${CMAKE_INSTALL_PREFIX}
  PREFIX            ${CMAKE_BINARY_DIR}/extern/hypre-cmake
  SOURCE_SUBDIR     src
  UPDATE_COMMAND    ""
  PATCH_COMMAND
    git reset --hard &&
    git clean -fd &&
    git apply "${HYPRE_PATCH_FILES}"
  CMAKE_ARGS        ${HYPRE_OPTIONS}
  TEST_COMMAND      ""
)
