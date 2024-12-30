# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Build MFEM
#

# Force build order
if(PALACE_BUILD_EXTERNAL_DEPS)
  set(MFEM_DEPENDENCIES hypre metis)
  if(PALACE_WITH_MUMPS)
    list(APPEND MFEM_DEPENDENCIES mumps)
  endif()
  if(PALACE_WITH_STRUMPACK)
    list(APPEND MFEM_DEPENDENCIES strumpack)
  endif()
  if(PALACE_WITH_SUPERLU)
    list(APPEND MFEM_DEPENDENCIES superlu_dist)
  endif()
  if(PALACE_WITH_SUNDIALS)
    list(APPEND MFEM_DEPENDENCIES sundials)
  endif()
else()
  set(MFEM_DEPENDENCIES)
endif()
if(PALACE_WITH_GSLIB)
  list(APPEND MFEM_DEPENDENCIES gslib)
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

# Silence some CUDA/HIP include file warnings
if(PALACE_WITH_CUDA OR PALACE_WITH_HIP)
  set(MFEM_CUDA_FLAGS "${CMAKE_CUDA_FLAGS}")
  set(MFEM_HIP_FLAGS "${CMAKE_HIP_FLAGS}")
  if(MFEM_CXX_FLAGS MATCHES "-pedantic")
    string(REGEX REPLACE "-pedantic" "" MFEM_CXX_FLAGS ${MFEM_CXX_FLAGS})
  endif()
  if(MFEM_CUDA_FLAGS MATCHES "-pedantic")
    string(REGEX REPLACE "-pedantic" "" MFEM_CUDA_FLAGS ${MFEM_CUDA_FLAGS})
  endif()
  if(MFEM_HIP_FLAGS MATCHES "-pedantic")
    string(REGEX REPLACE "-pedantic" "" MFEM_HIP_FLAGS ${MFEM_HIP_FLAGS})
  endif()
  if(MFEM_CUDA_FLAGS MATCHES "-ccbin")
    # MFEM adds this via CMAKE_CUDA_HOST_COMPILER
    string(REGEX REPLACE "-ccbin ([^ ]+)" "" MFEM_CUDA_FLAGS ${MFEM_CUDA_FLAGS})
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

# Replace mfem abort calls with exceptions for testing, default off
set(PALACE_MFEM_USE_EXCEPTIONS CACHE NO "MFEM throw exceptsions instead of abort calls")

set(MFEM_OPTIONS ${PALACE_SUPERBUILD_DEFAULT_ARGS})
list(APPEND MFEM_OPTIONS
  "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
  "-DCMAKE_CXX_FLAGS=${MFEM_CXX_FLAGS}"
  "-DMFEM_USE_MPI=YES"
  "-DMFEM_USE_OPENMP=${PALACE_WITH_OPENMP}"
  "-DMFEM_THREAD_SAFE=${PALACE_WITH_OPENMP}"
  "-DMFEM_USE_SUPERLU=${PALACE_WITH_SUPERLU}"
  "-DMFEM_USE_STRUMPACK=${PALACE_WITH_STRUMPACK}"
  "-DMFEM_USE_MUMPS=${PALACE_WITH_MUMPS}"
  "-DMFEM_USE_ZLIB=${PALACE_MFEM_WITH_ZLIB}"
  "-DMFEM_USE_LIBUNWIND=${PALACE_MFEM_WITH_LIBUNWIND}"
  "-DMFEM_USE_METIS_5=YES"
  "-DMFEM_USE_CEED=NO"
  "-DMFEM_USE_SUNDIALS=${PALACE_WITH_SUNDIALS}"
  "-DMFEM_USE_EXCEPTIONS=${PALACE_MFEM_USE_EXCEPTIONS}"
)
if(PALACE_WITH_STRUMPACK OR PALACE_WITH_MUMPS)
  list(APPEND MFEM_OPTIONS
    "-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}"
    "-DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}"
  )
endif()

# Configure BLAS/LAPACK for dependencies
if(NOT "${BLAS_LAPACK_LIBRARIES}" STREQUAL "")
  list(APPEND MFEM_OPTIONS
    # "-DMFEM_USE_LAPACK=YES"
    "-DBLAS_LIBRARIES=${BLAS_LAPACK_LIBRARIES}"
    "-DLAPACK_LIBRARIES=${BLAS_LAPACK_LIBRARIES}"
  )
endif()

# Configure GPU support
if(PALACE_WITH_CUDA)
  list(APPEND MFEM_OPTIONS
    "-DMFEM_USE_CUDA=YES"
    "-DCMAKE_CUDA_COMPILER=${CMAKE_CUDA_COMPILER}"
    "-DCMAKE_CUDA_FLAGS=${MFEM_CUDA_FLAGS}"
  )
  if(NOT "${CMAKE_CUDA_ARCHITECTURES}" STREQUAL "")
    list(APPEND MFEM_OPTIONS
      "-DCMAKE_CUDA_ARCHITECTURES=${CMAKE_CUDA_ARCHITECTURES}"
      "-DCUDA_ARCH=${CMAKE_CUDA_ARCHITECTURES}"
    )
  endif()
else()
  list(APPEND MFEM_OPTIONS
    "-DMFEM_USE_CUDA=NO"
  )
endif()
if(PALACE_WITH_HIP)
  list(APPEND MFEM_OPTIONS
    "-DMFEM_USE_HIP=YES"
    "-DROCM_PATH=${ROCM_DIR}"
    "-DCMAKE_HIP_COMPILER=${CMAKE_HIP_COMPILER}"
    "-DCMAKE_HIP_FLAGS=${MFEM_HIP_FLAGS}"
  )
  if(NOT "${CMAKE_HIP_ARCHITECTURES}" STREQUAL "")
    list(APPEND MFEM_OPTIONS
      "-DCMAKE_HIP_ARCHITECTURES=${CMAKE_HIP_ARCHITECTURES}"
      "-DHIP_ARCH=${CMAKE_HIP_ARCHITECTURES}"
    )
  endif()
else()
  list(APPEND MFEM_OPTIONS
    "-DMFEM_USE_HIP=NO"
  )
endif()

# MFEM with GSLIB is always built internally
if(PALACE_WITH_GSLIB)
  list(APPEND MFEM_OPTIONS
    "-DMFEM_USE_GSLIB=YES"
    "-DGSLIB_DIR=${CMAKE_INSTALL_PREFIX}"
  )
endif()

# Configure the rest of MFEM's dependencies
if(PALACE_BUILD_EXTERNAL_DEPS)
  list(APPEND MFEM_OPTIONS
    "-DMETIS_LIBRARIES=${METIS_LIBRARIES}"
    "-DMETIS_INCLUDE_DIRS=${CMAKE_INSTALL_PREFIX}/include"
    "-DHYPRE_DIR=${CMAKE_INSTALL_PREFIX}"
    "-DHYPRE_REQUIRED_PACKAGES=LAPACK$<SEMICOLON>BLAS"
  )
  if(PALACE_WITH_SUPERLU OR PALACE_WITH_STRUMPACK)
    list(APPEND MFEM_OPTIONS
      "-DParMETIS_LIBRARIES=${PARMETIS_LIBRARIES}$<SEMICOLON>${METIS_LIBRARIES}"
      "-DParMETIS_INCLUDE_DIRS=${CMAKE_INSTALL_PREFIX}/include"
    )
  endif()

  # HYPRE is built with cusparse, curand (or HIP counterparts), and these are added to
  # HYPRE_LIBRARIES by the MFEM CMake build. However, this ignores the include directories
  # (for #include <cusparse.h>, for example), which we can add this way via CMake defining
  # CUDAToolkit_INCLUDE_DIRS (and the HIP counterpart).
  if(PALACE_WITH_CUDA)
    list(APPEND MFEM_OPTIONS
      "-DHYPRE_REQUIRED_PACKAGES=CUDAToolkit"
    )
  endif()
  if(PALACE_WITH_HIP)
    list(APPEND MFEM_OPTIONS
      "-DHYPRE_REQUIRED_PACKAGES=rocsparse"
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

  # Find cuBLAS and cuSOLVER (or ROCm counterparts)
  if(PALACE_WITH_SUPERLU OR PALACE_WITH_STRUMPACK)
    if(PALACE_WITH_CUDA)
      find_package(CUDAToolkit REQUIRED)
      get_target_property(SUPERLU_STRUMPACK_CUBLAS_LIBRARY CUDA::cublas LOCATION)
      get_target_property(SUPERLU_STRUMPACK_CUBLASLT_LIBRARY CUDA::cublasLt LOCATION)
      get_target_property(SUPERLU_STRUMPACK_CUSOLVER_LIBRARY CUDA::cusolver LOCATION)
      set(SUPERLU_STRUMPACK_CUDA_LIBRARIES
        ${SUPERLU_STRUMPACK_CUBLAS_LIBRARY}
        ${SUPERLU_STRUMPACK_CUBLASLT_LIBRARY}
        ${SUPERLU_STRUMPACK_CUSOLVER_LIBRARY}
      )
    endif()
    if(PALACE_WITH_HIP)
      find_package(hipblas REQUIRED)
      find_package(rocblas REQUIRED)
      find_package(rocsolver REQUIRED)
      get_target_property(SUPERLU_STRUMPACK_HIPBLAS_LIBRARY roc::hipblas LOCATION)
      get_target_property(SUPERLU_STRUMPACK_ROCBLAS_LIBRARY roc::rocblas LOCATION)
      get_target_property(SUPERLU_STRUMPACK_ROCSOLVER_LIBRARY roc::rocsolver LOCATION)
      set(SUPERLU_STRUMPACK_ROCM_LIBRARIES
        ${SUPERLU_STRUMPACK_ROCBLAS_LIBRARY}
        ${SUPERLU_STRUMPACK_HIPBLAS_LIBRARY}
        ${SUPERLU_STRUMPACK_ROCSOLVER_LIBRARY}
      )
    endif()
  endif()

  # Configure SuperLU_DIST
  if(PALACE_WITH_SUPERLU)
    set(SUPERLU_REQUIRED_PACKAGES "ParMETIS" "METIS" "LAPACK" "BLAS" "MPI")
    set(SUPERLU_REQUIRED_LIBRARIES)
    if(PALACE_WITH_OPENMP)
      list(APPEND SUPERLU_REQUIRED_PACKAGES "OpenMP")
    endif()
    if(PALACE_WITH_CUDA)
      list(APPEND SUPERLU_REQUIRED_PACKAGES "CUDAToolkit")
      list(APPEND SUPERLU_REQUIRED_LIBRARIES ${SUPERLU_STRUMPACK_CUDA_LIBRARIES})
    endif()
    if(PALACE_WITH_HIP)
      list(APPEND SUPERLU_REQUIRED_PACKAGES "hipblas$<SEMICOLON>rocblas")
      list(APPEND SUPERLU_REQUIRED_LIBRARIES ${SUPERLU_STRUMPACK_ROCM_LIBRARIES})
    endif()
    string(REPLACE ";" "$<SEMICOLON>" SUPERLU_REQUIRED_PACKAGES "${SUPERLU_REQUIRED_PACKAGES}")
    string(REPLACE ";" "$<SEMICOLON>" SUPERLU_REQUIRED_LIBRARIES "${SUPERLU_REQUIRED_LIBRARIES}")
    list(APPEND MFEM_OPTIONS
      "-DSuperLUDist_DIR=${CMAKE_INSTALL_PREFIX}"
      "-DSuperLUDist_REQUIRED_PACKAGES=${SUPERLU_REQUIRED_PACKAGES}"
    )
    if(NOT "${SUPERLU_REQUIRED_LIBRARIES}" STREQUAL "")
      list(APPEND MFEM_OPTIONS
        "-DSuperLUDist_REQUIRED_LIBRARIES=${SUPERLU_REQUIRED_LIBRARIES}"
      )
    endif()
  endif()

  # Configure STRUMPACK
  if(PALACE_WITH_STRUMPACK)
    set(STRUMPACK_REQUIRED_PACKAGES "ParMETIS" "METIS" "LAPACK" "BLAS" "MPI" "MPI_Fortran")
    set(STRUMPACK_REQUIRED_LIBRARIES ${SCALAPACK_LIBRARIES} ${STRUMPACK_MUMPS_GFORTRAN_LIBRARY})
    if(NOT "${STRUMPACK_EXTRA_LIBRARIES}" STREQUAL "")
      list(PREPEND STRUMPACK_REQUIRED_LIBRARIES ${STRUMPACK_EXTRA_LIBRARIES})
    endif()
    if(PALACE_WITH_OPENMP)
      list(APPEND STRUMPACK_REQUIRED_PACKAGES "OpenMP")
    endif()
    if(PALACE_WITH_CUDA)
      list(APPEND STRUMPACK_REQUIRED_PACKAGES "CUDAToolkit")
      list(APPEND STRUMPACK_REQUIRED_LIBRARIES ${SUPERLU_STRUMPACK_CUDA_LIBRARIES})
    endif()
    if(PALACE_WITH_HIP)
      list(APPEND STRUMPACK_REQUIRED_PACKAGES "hipblas$<SEMICOLON>rocblas")
      list(APPEND STRUMPACK_REQUIRED_LIBRARIES ${SUPERLU_STRUMPACK_ROCM_LIBRARIES})
    endif()
    string(REPLACE ";" "$<SEMICOLON>" STRUMPACK_REQUIRED_PACKAGES "${STRUMPACK_REQUIRED_PACKAGES}")
    string(REPLACE ";" "$<SEMICOLON>" STRUMPACK_REQUIRED_LIBRARIES "${STRUMPACK_REQUIRED_LIBRARIES}")
    list(APPEND MFEM_OPTIONS
      "-DSTRUMPACK_DIR=${CMAKE_INSTALL_PREFIX}"
      "-DSTRUMPACK_REQUIRED_PACKAGES=${STRUMPACK_REQUIRED_PACKAGES}"
      "-DSTRUMPACK_REQUIRED_LIBRARIES=${STRUMPACK_REQUIRED_LIBRARIES}"
    )
  endif()

  # Configure MUMPS
  if(PALACE_WITH_MUMPS)
    set(MUMPS_REQUIRED_PACKAGES "METIS" "LAPACK" "BLAS" "MPI" "MPI_Fortran" "Threads")
    if(PALACE_WITH_OPENMP)
      list(APPEND MUMPS_REQUIRED_PACKAGES "OpenMP")
    endif()
    string(REPLACE ";" "$<SEMICOLON>" MUMPS_REQUIRED_PACKAGES "${MUMPS_REQUIRED_PACKAGES}")
    list(APPEND MFEM_OPTIONS
      "-DMUMPS_DIR=${CMAKE_INSTALL_PREFIX}"
      "-DMUMPS_REQUIRED_PACKAGES=${MUMPS_REQUIRED_PACKAGES}"
      "-DMUMPS_REQUIRED_LIBRARIES=${SCALAPACK_LIBRARIES}$<SEMICOLON>${STRUMPACK_MUMPS_GFORTRAN_LIBRARY}"
    )
  endif()

  # Configure SUNDIALS
  if(PALACE_WITH_SUNDIALS)
    set(SUNDIALS_REQUIRED_PACKAGES "LAPACK" "BLAS" "MPI")
    if(PALACE_WITH_OPENMP)
      list(APPEND SUNDIALS_REQUIRED_PACKAGES "OpenMP")
    endif()
    if(PALACE_WITH_CUDA)
      list(APPEND SUNDIALS_REQUIRED_PACKAGES "CUDAToolkit")
      list(APPEND SUNDIALS_REQUIRED_LIBRARIES ${SUPERLU_STRUMPACK_CUDA_LIBRARIES})
    endif()
    string(REPLACE ";" "$<SEMICOLON>" SUNDIALS_REQUIRED_PACKAGES "${SUNDIALS_REQUIRED_PACKAGES}")
    string(REPLACE ";" "$<SEMICOLON>" SUNDIALS_REQUIRED_LIBRARIES "${SUNDIALS_REQUIRED_LIBRARIES}")
    list(APPEND MFEM_OPTIONS
      "-DSUNDIALS_DIR=${CMAKE_INSTALL_PREFIX}"
      "-DSUNDIALS_REQUIRED_PACKAGES=${SUNDIALS_REQUIRED_PACKAGES}"
    )
    if(NOT "${SUNDIALS_REQUIRED_LIBRARIES}" STREQUAL "")
      list(APPEND MFEM_OPTIONS
        "-DSUNDIALS_REQUIRED_LIBRARIES=${SUNDIALS_REQUIRED_LIBRARIES}"
      )
    endif()
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
    "SuperLUDist"
    "STRUMPACK"
    "MUMPS"
    "SUNDIALS"
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
  "${CMAKE_SOURCE_DIR}/extern/patch/mfem/patch_mfem_device_fixes.diff"
  "${CMAKE_SOURCE_DIR}/extern/patch/mfem/patch_mesh_vis_dev.diff"
  "${CMAKE_SOURCE_DIR}/extern/patch/mfem/patch_mesh_prism_vtu_fix.diff"
  "${CMAKE_SOURCE_DIR}/extern/patch/mfem/patch_par_tet_mesh_fix_dev.diff"
  "${CMAKE_SOURCE_DIR}/extern/patch/mfem/patch_gmsh_parser_performance.diff"
  "${CMAKE_SOURCE_DIR}/extern/patch/mfem/patch_gmsh_reader_periodic_bugfix.diff"
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
    git reset --hard &&
    git clean -fd &&
    git apply "${MFEM_PATCH_FILES}"
  CONFIGURE_COMMAND ${CMAKE_COMMAND} <SOURCE_DIR> "${MFEM_OPTIONS}"
  TEST_COMMAND      ${CMAKE_MAKE_PROGRAM} ex1 ex1p
)
