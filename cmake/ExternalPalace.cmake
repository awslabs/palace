# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Build Palace
#

# The superbuild owns the exact dependency source revisions. Target existence is the
# authority for which pins are reported, avoiding a second copy of the feature-selection
# predicates in extern/CMakeLists.txt.
include(BuildDependencyMetadata)
if(PALACE_BUILD_EXTERNAL_DEPS)
  # Force the same direct build order as before, but let configured targets determine which
  # optional eigensolver is present.
  set(PALACE_DEPENDENCIES)
  foreach(_palace_dependency IN ITEMS
      mfem libCEED json json_schema_validator fmt eigen scn slepc arpack-ng)
    if(TARGET "${_palace_dependency}")
      list(APPEND PALACE_DEPENDENCIES "${_palace_dependency}")
    endif()
  endforeach()

  # Keep these labels aligned with the Spack helper where the same dependency is present.
  # Target/spec presence is authoritative in each producer.
  set(_palace_build_dependencies)
  palace_append_external_build_dependency(_palace_build_dependencies strumpack STRUMPACK
    "${EXTERN_STRUMPACK_GIT_TAG}")
  palace_append_external_build_dependency(_palace_build_dependencies arpack-ng arpack_ng
    "${EXTERN_ARPACK_GIT_TAG}")
  palace_append_external_build_dependency(_palace_build_dependencies butterflypack
    butterflypack "${EXTERN_BUTTERFLYPACK_GIT_TAG}")
  palace_append_external_build_dependency(_palace_build_dependencies eigen eigen
    "${EXTERN_EIGEN_VERSION}")
  palace_append_external_build_dependency(_palace_build_dependencies fmt fmt
    "${EXTERN_FMT_VERSION}")
  palace_append_external_build_dependency(_palace_build_dependencies gslib gslib
    "${EXTERN_GSLIB_GIT_TAG}")
  palace_append_external_build_dependency(_palace_build_dependencies hypre hypre
    "${EXTERN_HYPRE_GIT_TAG}")
  palace_append_external_build_dependency(_palace_build_dependencies json json
    "${EXTERN_JSON_VERSION}")
  palace_append_external_build_dependency(_palace_build_dependencies json_schema_validator
    json_schema_validator "${EXTERN_JSON_SCHEMA_VALIDATOR_VERSION}")
  palace_append_external_build_dependency(_palace_build_dependencies libCEED libCEED
    "${EXTERN_LIBCEED_GIT_TAG}")
  palace_append_external_build_dependency(_palace_build_dependencies libxsmm libxsmm
    "${EXTERN_LIBXSMM_GIT_TAG}")
  palace_append_external_build_dependency(_palace_build_dependencies magma magma
    "${EXTERN_MAGMA_GIT_TAG}")
  palace_append_external_build_dependency(_palace_build_dependencies metis metis
    "${EXTERN_METIS_GIT_TAG}")
  palace_append_external_build_dependency(_palace_build_dependencies mfem mfem
    "${EXTERN_MFEM_GIT_TAG}")
  palace_append_external_build_dependency(_palace_build_dependencies mumps mumps
    "${EXTERN_MUMPS_GIT_TAG}")
  palace_append_external_build_dependency(_palace_build_dependencies parmetis parmetis
    "${EXTERN_PARMETIS_GIT_TAG}")
  palace_append_external_build_dependency(_palace_build_dependencies petsc petsc
    "${EXTERN_PETSC_GIT_TAG}")
  palace_append_external_build_dependency(_palace_build_dependencies scalapack scalapack
    "${EXTERN_SCALAPACK_GIT_TAG}")
  palace_append_external_build_dependency(_palace_build_dependencies scn scn
    "${EXTERN_SCN_VERSION}")
  palace_append_external_build_dependency(_palace_build_dependencies slepc slepc
    "${EXTERN_SLEPC_GIT_TAG}")
  palace_append_external_build_dependency(_palace_build_dependencies sundials sundials
    "${EXTERN_SUNDIALS_GIT_TAG}")
  palace_append_external_build_dependency(_palace_build_dependencies superlu_dist
    superlu_dist "${EXTERN_SUPERLU_GIT_TAG}")
  palace_append_external_build_dependency(_palace_build_dependencies umpire umpire
    "${EXTERN_UMPIRE_GIT_TAG}")
  palace_append_external_build_dependency(_palace_build_dependencies zfp zfp
    "${EXTERN_ZFP_GIT_TAG}")

  list(SORT _palace_build_dependencies)
  list(JOIN _palace_build_dependencies "|" PALACE_BUILD_DEPENDENCIES)
endif()

set(PALACE_OPTIONS ${PALACE_SUPERBUILD_DEFAULT_ARGS})
list(APPEND PALACE_OPTIONS
  "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
  "-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
  "-DPALACE_WITH_OPENMP=${PALACE_WITH_OPENMP}"
  "-DPALACE_WITH_SLEPC=${PALACE_WITH_SLEPC}"
  "-DPALACE_WITH_ARPACK=${PALACE_WITH_ARPACK}"
  "-DPALACE_WITH_SUNDIALS=${PALACE_WITH_SUNDIALS}"
  "-DPALACE_WITH_STRUMPACK=${PALACE_WITH_STRUMPACK}"
  "-DPALACE_WITH_SUPERLU=${PALACE_WITH_SUPERLU}"
  "-DPALACE_WITH_MUMPS=${PALACE_WITH_MUMPS}"
  "-DPALACE_WITH_GSLIB=${PALACE_WITH_GSLIB}"
  "-DPALACE_BUILD_DEPENDENCIES=${PALACE_BUILD_DEPENDENCIES}"
  "-DANALYZE_SOURCES_CLANG_TIDY=${ANALYZE_SOURCES_CLANG_TIDY}"
  "-DANALYZE_SOURCES_CPPCHECK=${ANALYZE_SOURCES_CPPCHECK}"
  "-DPALACE_BUILD_EXTERNAL_DEPS=${PALACE_BUILD_EXTERNAL_DEPS}" # For Catch2
  "-DPALACE_BUILD_WITH_COVERAGE=${PALACE_BUILD_WITH_COVERAGE}"
  "-DPALACE_BUILD_WITH_SANITIZERS=${PALACE_BUILD_WITH_SANITIZERS}"
)
if(NOT "${MFEM_DIR}" STREQUAL "")
  list(APPEND PALACE_OPTIONS "-DMFEM_DIR=${MFEM_DIR}")
elseif(PALACE_BUILD_EXTERNAL_DEPS)
  list(APPEND PALACE_OPTIONS "-DMFEM_DIR=${CMAKE_INSTALL_PREFIX}")
endif()
if(NOT "${MUMPS_DIR}" STREQUAL "")
  list(APPEND PALACE_OPTIONS "-DMUMPS_DIR=${MUMPS_DIR}")
elseif(PALACE_BUILD_EXTERNAL_DEPS AND PALACE_WITH_MUMPS)
  list(APPEND PALACE_OPTIONS "-DMUMPS_DIR=${CMAKE_INSTALL_PREFIX}")
endif()
if(NOT "${STRUMPACK_DIR}" STREQUAL "")
  list(APPEND PALACE_OPTIONS "-DSTRUMPACK_DIR=${STRUMPACK_DIR}")
elseif(PALACE_BUILD_EXTERNAL_DEPS AND PALACE_WITH_STRUMPACK)
  list(APPEND PALACE_OPTIONS "-DSTRUMPACK_DIR=${CMAKE_INSTALL_PREFIX}")
endif()
if(NOT "${SUPERLU_DIST_DIR}" STREQUAL "")
  list(APPEND PALACE_OPTIONS "-DSUPERLU_DIST_DIR=${SUPERLU_DIST_DIR}")
elseif(PALACE_BUILD_EXTERNAL_DEPS AND PALACE_WITH_SUPERLU)
  list(APPEND PALACE_OPTIONS "-DSUPERLU_DIST_DIR=${CMAKE_INSTALL_PREFIX}")
endif()
if(NOT "${METIS_DIR}" STREQUAL "")
  list(APPEND PALACE_OPTIONS "-DMETIS_DIR=${METIS_DIR}")
elseif(PALACE_BUILD_EXTERNAL_DEPS)
  list(APPEND PALACE_OPTIONS "-DMETIS_DIR=${CMAKE_INSTALL_PREFIX}")
endif()
if(NOT "${PARMETIS_DIR}" STREQUAL "")
  list(APPEND PALACE_OPTIONS "-DPARMETIS_DIR=${PARMETIS_DIR}")
elseif(PALACE_BUILD_EXTERNAL_DEPS)
  list(APPEND PALACE_OPTIONS "-DPARMETIS_DIR=${CMAKE_INSTALL_PREFIX}")
endif()
if(NOT "${HYPRE_DIR}" STREQUAL "")
  list(APPEND PALACE_OPTIONS "-DHYPRE_DIR=${HYPRE_DIR}")
elseif(PALACE_BUILD_EXTERNAL_DEPS)
  list(APPEND PALACE_OPTIONS "-DHYPRE_DIR=${CMAKE_INSTALL_PREFIX}")
endif()
if(PALACE_BUILD_EXTERNAL_DEPS AND PALACE_WITH_SUNDIALS)
  list(APPEND PALACE_OPTIONS "-DSUNDIALS_DIR=${CMAKE_INSTALL_PREFIX}")
endif()
if(PALACE_BUILD_EXTERNAL_DEPS AND PALACE_WITH_GSLIB)
  list(APPEND PALACE_OPTIONS "-DGSLIB_DIR=${CMAKE_INSTALL_PREFIX}")
endif()
if(PALACE_WITH_ARPACK)
  list(APPEND PALACE_OPTIONS
    "-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}"
    "-DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}"
  )
endif()

# Configure LAPACK dependency
if(NOT "${BLAS_LAPACK_LIBRARIES}" STREQUAL "")
  list(APPEND PALACE_OPTIONS
    "-DBLAS_LIBRARIES=${BLAS_LAPACK_LIBRARIES}"
    "-DLAPACK_LIBRARIES=${BLAS_LAPACK_LIBRARIES}"
  )
endif()

# Configure GPU support
if(PALACE_WITH_CUDA)
  list(APPEND PALACE_OPTIONS
    "-DPALACE_WITH_CUDA=ON"
    "-DPALACE_WITH_GPU_AWARE_MPI=${PALACE_WITH_GPU_AWARE_MPI}"
    "-DCMAKE_CUDA_COMPILER=${CMAKE_CUDA_COMPILER}"
    "-DCMAKE_CUDA_FLAGS=${CMAKE_CUDA_FLAGS}"
  )
  if(NOT "${CMAKE_CUDA_ARCHITECTURES}" STREQUAL "")
    list(APPEND PALACE_OPTIONS
      "-DCMAKE_CUDA_ARCHITECTURES=${CMAKE_CUDA_ARCHITECTURES}"
    )
  endif()
else()
  list(APPEND PALACE_OPTIONS
    "-DPALACE_WITH_CUDA=OFF"
  )
endif()
if(PALACE_WITH_HIP)
  list(APPEND PALACE_OPTIONS
    "-DPALACE_WITH_HIP=ON"
    "-DPALACE_WITH_GPU_AWARE_MPI=${PALACE_WITH_GPU_AWARE_MPI}"
    "-DCMAKE_HIP_COMPILER=${CMAKE_HIP_COMPILER}"
    "-DCMAKE_HIP_FLAGS=${CMAKE_HIP_FLAGS}"
  )
  if(NOT "${CMAKE_HIP_ARCHITECTURES}" STREQUAL "")
    list(APPEND PALACE_OPTIONS
      "-DCMAKE_HIP_ARCHITECTURES=${CMAKE_HIP_ARCHITECTURES}"
    )
  endif()
else()
  list(APPEND PALACE_OPTIONS
    "-DPALACE_WITH_HIP=OFF"
  )
endif()

string(REPLACE ";" "; " PALACE_OPTIONS_PRINT "${PALACE_OPTIONS}")
message(STATUS "PALACE_OPTIONS: ${PALACE_OPTIONS_PRINT}")

include(ExternalProject)
ExternalProject_Add(palace
  DEPENDS           ${PALACE_DEPENDENCIES}
  SOURCE_DIR        ${CMAKE_SOURCE_DIR}/palace
  BINARY_DIR        ${CMAKE_BINARY_DIR}/palace-build
  INSTALL_DIR       ${CMAKE_INSTALL_PREFIX}
  PREFIX            ${CMAKE_BINARY_DIR}/palace-cmake
  BUILD_ALWAYS      TRUE
  DOWNLOAD_COMMAND  ""
  CONFIGURE_COMMAND ${CMAKE_COMMAND} <SOURCE_DIR> "${PALACE_OPTIONS}"
  TEST_COMMAND      ""
)

# Add target for Palace unit tests
ExternalProject_Add_Step(palace tests
  COMMAND           ${CMAKE_MAKE_PROGRAM} unit-tests
  COMMAND           ${CMAKE_COMMAND} --install <BINARY_DIR>/test/unit --prefix ${CMAKE_INSTALL_PREFIX}
  DEPENDEES         install
  DEPENDERS         ""
  COMMENT           "Building and installing unit tests for 'palace'"
  WORKING_DIRECTORY <BINARY_DIR>
  EXCLUDE_FROM_MAIN TRUE
)
ExternalProject_Add_StepTargets(palace tests)
