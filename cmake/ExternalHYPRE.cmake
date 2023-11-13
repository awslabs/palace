# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Build HYPRE
#

# Force build order
set(HYPRE_DEPENDENCIES)

# Hypre does not add OpenMP flags
set(HYPRE_CFLAGS "${CMAKE_C_FLAGS}")
set(HYPRE_CXXFLAGS "${CMAKE_CXX_FLAGS}")
set(HYPRE_LDFLAGS "${CMAKE_EXE_LINKER_FLAGS}")
if(PALACE_WITH_OPENMP)
  find_package(OpenMP REQUIRED)
  set(HYPRE_CFLAGS "${OpenMP_C_FLAGS} ${HYPRE_CFLAGS}")
  set(HYPRE_CXXFLAGS "${OpenMP_CXX_FLAGS} ${HYPRE_CXXFLAGS}")
  string(REPLACE ";" " " HYPRE_OPENMP_LIBRARIES "${OpenMP_C_LIBRARIES}")
  set(HYPRE_LDFLAGS "${HYPRE_OPENMP_LIBRARIES} ${HYPRE_LDFLAGS}")
endif()

# Silence some CUDA/HIP include file warnings
if(HYPRE_CFLAGS MATCHES "-pedantic" AND (PALACE_WITH_CUDA OR PALACE_WITH_HIP))
  string(REGEX REPLACE "-pedantic" "" HYPRE_CFLAGS ${HYPRE_CFLAGS})
  string(REGEX REPLACE "-pedantic" "" HYPRE_CXXFLAGS ${HYPRE_CXXFLAGS})
endif()
if(PALACE_WITH_CUDA)
  set(HYPRE_CFLAGS "${HYPRE_CFLAGS} -isystem ${CUDA_DIR}/include")
  set(HYPRE_CXXFLAGS "${HYPRE_CXXFLAGS} -isystem ${CUDA_DIR}/include")
endif()
if(PALACE_WITH_HIP)
  set(HYPRE_CFLAGS "${HYPRE_CFLAGS} -isystem ${ROCM_DIR}/include")
  set(HYPRE_CXXFLAGS "${HYPRE_CXXFLAGS} -isystem ${ROCM_DIR}/include")
endif()

# Need to manually specificy MPI flags for test program compilation/linkage during configure
if(NOT MPI_FOUND)
  message(FATAL_ERROR "MPI is not found when trying to build HYPRE")
endif()
if(NOT CMAKE_C_COMPILER STREQUAL MPI_C_COMPILER)
  foreach(INCLUDE_DIR IN LISTS MPI_C_INCLUDE_DIRS)
    set(HYPRE_CFLAGS "${HYPRE_CFLAGS} -I${INCLUDE_DIR}")
    set(HYPRE_CXXFLAGS "${HYPRE_CXXFLAGS} -I${INCLUDE_DIR}")
  endforeach()
  string(REPLACE ";" " " HYPRE_MPI_LIBRARIES "${MPI_C_LIBRARIES}")
  set(HYPRE_LDFLAGS "${HYPRE_LDFLAGS} ${HYPRE_MPI_LIBRARIES}")
endif()

# Use Autotools build instead of CMake for HIP support
set(HYPRE_OPTIONS
  "CC=${CMAKE_C_COMPILER}"
  "CFLAGS=${HYPRE_CFLAGS}"
  "CXX=${CMAKE_CXX_COMPILER}"
  "CXXFLAGS=${HYPRE_CXXFLAGS}"
  "FC="
  "LDFLAGS=${HYPRE_LDFLAGS}"
  "--prefix=${CMAKE_INSTALL_PREFIX}"
  "--disable-fortran"
  "--with-MPI"
)
if(CMAKE_BUILD_TYPE MATCHES "Debug|debug|DEBUG")
  list(APPEND HYPRE_OPTIONS "--enable-debug")
else()
  list(APPEND HYPRE_OPTIONS "--disable-debug")
endif()
if(BUILD_SHARED_LIBS)
  list(APPEND HYPRE_OPTIONS "--enable-shared")
else()
  list(APPEND HYPRE_OPTIONS "--disable-shared")
endif()
if(PALACE_WITH_64BIT_INT)
  list(APPEND HYPRE_OPTIONS
    "--enable-mixedint"
    "--disable-bigint"
  )
else()
  list(APPEND HYPRE_OPTIONS
    "--disable-mixedint"
    "--disable-bigint"
  )
endif()
if(PALACE_WITH_OPENMP)
  list(APPEND HYPRE_OPTIONS "--with-openmp")
endif()

# User might specify the MPI compiler wrappers directly, otherwise we need to supply MPI
# as found from the CMake module
if(NOT CMAKE_C_COMPILER STREQUAL MPI_C_COMPILER)
  set(HYPRE_MPI_LIBRARIES)
  set(HYPRE_MPI_LIBRARY_DIRS)
  foreach(LIB IN LISTS MPI_C_LIBRARIES)
    get_filename_component(LIB_NAME ${LIB} NAME_WE)
    get_filename_component(LIB_DIR ${LIB} DIRECTORY)
    string(REGEX REPLACE "^lib" "" LIB_NAME ${LIB_NAME})
    list(APPEND HYPRE_MPI_LIBRARIES ${LIB_NAME})
    list(APPEND HYPRE_MPI_LIBRARY_DIRS ${LIB_DIR})
  endforeach()
  list(REMOVE_DUPLICATES HYPRE_MPI_LIBRARIES)
  string(REPLACE ";" " " HYPRE_MPI_LIBRARIES "${HYPRE_MPI_LIBRARIES}")
  list(REMOVE_DUPLICATES HYPRE_MPI_LIBRARY_DIRS)
  string(REPLACE ";" " " HYPRE_MPI_LIBRARY_DIRS "${HYPRE_MPI_LIBRARY_DIRS}")
  string(REPLACE ";" " " HYPRE_MPI_INCLUDE_DIRS "${MPI_C_INCLUDE_DIRS}")
  list(APPEND HYPRE_OPTIONS
    "--with-MPI-libs=${HYPRE_MPI_LIBRARIES}"
    "--with-MPI-lib-dirs=${HYPRE_MPI_LIBRARY_DIRS}"
    "--with-MPI-include=${HYPRE_MPI_INCLUDE_DIRS}"
  )
endif()

# Configure BLAS/LAPACK
if(NOT "${BLAS_LAPACK_LIBRARIES}" STREQUAL "")
  string(REPLACE "$<SEMICOLON>" " " HYPRE_BLAS_LAPACK_LIBRARIES "${BLAS_LAPACK_LIBRARIES}")
  list(APPEND HYPRE_OPTIONS
    "--with-blas-lib=${HYPRE_BLAS_LAPACK_LIBRARIES}"
    "--with-lapack-lib=${HYPRE_BLAS_LAPACK_LIBRARIES}"
  )
endif()

# Configure GPU support
if(PALACE_WITH_CUDA OR PALACE_WITH_HIP)
  list(APPEND HYPRE_OPTIONS
    "--disable-unified-memory"
  )
  if(PALACE_WITH_GPU_AWARE_MPI)
    list(APPEND HYPRE_OPTIONS
      "--enable-gpu-aware-mpi"
    )
  else()
    list(APPEND HYPRE_OPTIONS
      "--disable-gpu-aware-mpi"
    )
  endif()
endif()
if(PALACE_WITH_CUDA)
  list(APPEND HYPRE_OPTIONS
    "--with-cuda"
    "--with-cuda-home=${CUDA_DIR}"
    "--enable-curand"
    "--enable-cusparse"
    "--enable-device-memory-pool"
  )
  if(NOT "${CMAKE_CUDA_ARCHITECTURES}" STREQUAL "")
    list(GET CMAKE_CUDA_ARCHITECTURES 0 HYPRE_CUDA_ARCH)
    list(APPEND HYPRE_OPTIONS
      "--with-gpu-arch=${HYPRE_CUDA_ARCH}"
    )
  endif()
endif()
if(PALACE_WITH_HIP)
  list(APPEND HYPRE_OPTIONS
    "--with-hip"
    "ROCM_PATH=${ROCM_DIR}"
    "--enable-rocrand"
    "--enable-rocsparse"
  )
  if(NOT "${CMAKE_HIP_ARCHITECTURES}" STREQUAL "")
    list(GET CMAKE_HIP_ARCHITECTURES 0 HYPRE_HIP_ARCH)
    list(APPEND HYPRE_OPTIONS
      "--with-gpu-arch=${HYPRE_HIP_ARCH}"
    )
  endif()
endif()

string(REPLACE ";" "; " HYPRE_OPTIONS_PRINT "${HYPRE_OPTIONS}")
message(STATUS "HYPRE_OPTIONS: ${HYPRE_OPTIONS_PRINT}")

include(ExternalProject)
ExternalProject_Add(hypre
  DEPENDS           ${HYPRE_DEPENDENCIES}
  GIT_REPOSITORY    ${EXTERN_HYPRE_URL}
  GIT_TAG           ${EXTERN_HYPRE_GIT_TAG}
  SOURCE_DIR        ${CMAKE_BINARY_DIR}/extern/hypre
  INSTALL_DIR       ${CMAKE_INSTALL_PREFIX}
  PREFIX            ${CMAKE_BINARY_DIR}/extern/hypre-cmake
  BUILD_IN_SOURCE   TRUE
  SOURCE_SUBDIR     src
  UPDATE_COMMAND    ""
  CONFIGURE_COMMAND ./configure ${HYPRE_OPTIONS}
  TEST_COMMAND      ""
)
