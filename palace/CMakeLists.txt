# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# CMake configuration for the main Palace application
#

# CMake 3.21 was released in Jul. 2021 (required for HIP support)
cmake_minimum_required(VERSION 3.21)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Prohibit in-source builds
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR "In-source builds are prohibited")
endif()

# C++17 required for std::filesystem, among others
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

# Initialize the project
project(palace LANGUAGES CXX VERSION 0.15.0)

# Define build settings and defaults
set(PALACE_WITH_OPENMP OFF CACHE BOOL "Use OpenMP for shared-memory parallelism")
set(PALACE_WITH_CUDA OFF CACHE BOOL "Use CUDA for NVIDIA GPU support")
set(PALACE_WITH_HIP OFF CACHE BOOL "Use HIP for AMD or NVIDIA GPU support")

set(PALACE_WITH_SUPERLU ON CACHE BOOL "Build with SuperLU_DIST sparse direct solver")
set(PALACE_WITH_STRUMPACK OFF CACHE BOOL "Build with STRUMPACK sparse direct solver")
set(PALACE_WITH_MUMPS OFF CACHE BOOL "Build with MUMPS sparse direct solver")
set(PALACE_WITH_SLEPC ON CACHE BOOL "Build with SLEPc eigenvalue solver")
set(PALACE_WITH_ARPACK OFF CACHE BOOL "Build with ARPACK eigenvalue solver")
set(PALACE_WITH_SUNDIALS ON CACHE BOOL "Build with SUNDIALS time integrators")
set(PALACE_WITH_GSLIB ON CACHE BOOL "Build with GSLIB for high-order field interpolation")

set(ANALYZE_SOURCES_CLANG_TIDY OFF CACHE BOOL "Run static analysis checks using clang-tidy")
set(ANALYZE_SOURCES_CPPCHECK OFF CACHE BOOL "Run static analysis checks using cppcheck")

set(PALACE_BUILD_WITH_COVERAGE OFF CACHE BOOL "Compile Palace with coverage flags (source-based for LLVM; gcov for GCC)")

# Help find third-party dependencies
set(MFEM_DIR "" CACHE STRING
  "Path to MFEM build or installation directory (not required if already on CMAKE_PREFIX_PATH)"
)
set(LIBCEED_DIR "" CACHE STRING
  "Path to libCEED build or installation directory (not required if already on CMAKE_PREFIX_PATH or PKG_CONFIG_PATH)"
)
set(PETSC_DIR "" CACHE STRING
  "Path to PETSc build or installation directory (not required if already on CMAKE_PREFIX_PATH or PKG_CONFIG_PATH)"
)
set(SLEPC_DIR "" CACHE STRING
  "Path to SLEPc build or installation directory (not required if already on CMAKE_PREFIX_PATH or PKG_CONFIG_PATH)"
)
set(ARPACK_DIR "" CACHE STRING
  "Path to ARPACK build or installation directory (not required if already on CMAKE_PREFIX_PATH)"
)
set(HYPRE_DIR "" CACHE STRING
  "Path to HYPRE build or installation directory (not required if already on CMAKE_PREFIX_PATH)"
)
set(MUMPS_DIR "" CACHE STRING
  "Path to MUMPS build or installation directory (not required if already on CMAKE_PREFIX_PATH)"
)
set(STRUMPACK_DIR "" CACHE STRING
  "Path to STRUMPACK build or installation directory (not required if already on CMAKE_PREFIX_PATH)"
)
set(SUPERLU_DIST_DIR "" CACHE STRING
  "Path to SuperLU_DIST build or installation directory (not required if already on CMAKE_PREFIX_PATH)"
)
set(METIS_DIR "" CACHE STRING
  "Path to METIS build or installation directory (not required if already on CMAKE_PREFIX_PATH)"
)
set(PARMETIS_DIR "" CACHE STRING
  "Path to ParMETIS build or installation directory (not required if already on CMAKE_PREFIX_PATH)"
)
set(SUNDIALS_DIR "" CACHE STRING
  "Path to SUNDIALS build or installation directory (not required if already on CMAKE_PREFIX_PATH)"
)
set(GSLIB_DIR "" CACHE STRING
  "Path to GSLIB build or installation directory (not required if already on CMAKE_PREFIX_PATH)"
)
set(SCALAPACK_DIR "" CACHE STRING
  "Path to ScaLAPACK build or installation directory (not required if already on CMAKE_PREFIX_PATH)"
)

# Enable Fortran if required
if(PALACE_WITH_ARPACK OR PALACE_WITH_STRUMPACK OR PALACE_WITH_MUMPS)
  enable_language(Fortran)
endif()

# Enable CUDA/HIP if required
if(PALACE_WITH_CUDA AND PALACE_WITH_HIP)
  message(FATAL_ERROR "PALACE_WITH_CUDA is not compatible with PALACE_WITH_HIP")
endif()
if(PALACE_WITH_CUDA)
  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-extended-lambda")
elseif(PALACE_WITH_HIP)
  enable_language(HIP)
endif()

# Set a default build type if none was provided
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting CMAKE_BUILD_TYPE to 'Release' as none was specified")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Specifies the build type ('Debug' or 'Release', for example)" FORCE
  )
endif()

# Centralize Error handling
add_library(common_warnings INTERFACE)
target_compile_options(
  common_warnings
  INTERFACE -Wall
            -Wpedantic
            -Wno-deprecated-declarations
)

if(PALACE_BUILD_WITH_COVERAGE)
  add_library(coverage_flags INTERFACE)
  target_compile_options(
    coverage_flags
    INTERFACE
        $<$<CXX_COMPILER_ID:AppleClang,Clang,IntelLLVM>:-fprofile-instr-generate>
        $<$<CXX_COMPILER_ID:AppleClang,Clang,IntelLLVM>:-fcoverage-mapping>
        $<$<CXX_COMPILER_ID:GNU>:--coverage>
        $<$<CXX_COMPILER_ID:GNU>:-fprofile-abs-path>
  )
  target_link_options(
    coverage_flags
    INTERFACE
        $<$<CXX_COMPILER_ID:AppleClang,Clang,IntelLLVM>:-fprofile-instr-generate>
        $<$<CXX_COMPILER_ID:AppleClang,Clang,IntelLLVM>:-fcoverage-mapping>
        $<$<CXX_COMPILER_ID:GNU>:--coverage>
  )
endif()

# Add extra CMake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Find MPI
find_package(MPI REQUIRED)
if(PALACE_WITH_ARPACK OR PALACE_WITH_STRUMPACK OR PALACE_WITH_MUMPS)
  find_package(MPI REQUIRED COMPONENTS Fortran)
endif()

# Find OpenMP
if(PALACE_WITH_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

# Find LAPACK
find_package(LAPACK REQUIRED)

# Find Umpire and Camp (needed when CUDA/HIP is enabled)
if(PALACE_WITH_CUDA OR PALACE_WITH_HIP)
  find_package(umpire REQUIRED CONFIG)
  find_package(camp REQUIRED CONFIG)
endif()

# Find nlohmann/json
find_package(nlohmann_json REQUIRED CONFIG)
message(STATUS "Found nlohmann/json: ${nlohmann_json_VERSION} in ${nlohmann_json_DIR}")

# Find nlohmann/json-schema-validator
find_package(nlohmann_json_schema_validator REQUIRED CONFIG)
message(STATUS "Found nlohmann_json_schema_validator: ${nlohmann_json_schema_validator_VERSION} in ${nlohmann_json_schema_validator_DIR}")

# Find fmt
find_package(fmt REQUIRED CONFIG)
message(STATUS "Found fmt: ${fmt_VERSION} in ${fmt_DIR}")

# Find scn
find_package(scn REQUIRED CONFIG)
message(STATUS "Found scn: ${scn_VERSION} in ${scn_DIR}")

# Find Eigen
find_package(Eigen3 REQUIRED CONFIG)
message(STATUS "Found Eigen: ${Eigen3_VERSION} in ${Eigen3_DIR}")

# Find MFEM (recent MFEM targets link to target Threads::Threads, for some reason)
find_package(Threads REQUIRED)

# MFEM in Spack is currently using autotools, and this find_package will fail
# (hence the QUIET), but this is not a problem.
if(NOT "${MFEM_DIR}" STREQUAL "")
  # MFEM_DIR provided, assume MFEM is available (e.g., from Spack with autotools)
  set(MFEM_ROOT ${MFEM_DIR})
  find_library(MFEM_LIBRARY NAMES mfem HINTS ${MFEM_DIR}/lib REQUIRED)
  message(STATUS "Found MFEM: ${MFEM_LIBRARY}")
else()
  # Try to find MFEM via CMake config
  find_package(MFEM QUIET CONFIG)
  if(MFEM_FOUND)
    message(STATUS "Found MFEM: ${MFEM_VERSION} in ${MFEM_DIR}")
    if(NOT MFEM_USE_MPI)
      message(FATAL_ERROR "Build requires MFEM with MPI support")
    endif()
    set(MFEM_LIBRARY mfem)
  else()
    message(FATAL_ERROR "MFEM not found. Set MFEM_DIR or ensure MFEM is on CMAKE_PREFIX_PATH")
  endif()
endif()


# Find libCEED
include(PkgConfigHelpers)
set(LIBCEED_TEST_DEPS)
if(PALACE_WITH_OPENMP)
  list(APPEND LIBCEED_TEST_DEPS OpenMP::OpenMP_CXX)
endif()
if(PALACE_WITH_CUDA)
  if(BUILD_SHARED_LIBS)
    list(APPEND LIBCEED_TEST_DEPS CUDA::cudart CUDA::cublas)
  else()
    list(APPEND LIBCEED_TEST_DEPS CUDA::cudart_static CUDA::cublas_static)
  endif()
elseif(PALACE_WITH_HIP)
  list(APPEND LIBCEED_TEST_DEPS hip::host)
endif()
find_libceed_pkgconfig("${LIBCEED_TEST_DEPS}" LIBCEED_TARGET)
if("${LIBCEED_TARGET}" STREQUAL "")
  message(FATAL_ERROR "libCEED could not be found, be sure to set LIBCEED_DIR")
endif()

# Find PETSc and SLEPc
if(PALACE_WITH_SLEPC)
  set(PETSC_TEST_DEPS MPI::MPI_CXX)
  if(PALACE_WITH_OPENMP)
    list(APPEND PETSC_TEST_DEPS OpenMP::OpenMP_CXX)
  endif()
  find_petsc_pkgconfig("${PETSC_TEST_DEPS}" PETSC_TARGET)
  if("${PETSC_TARGET}" STREQUAL "")
    message(FATAL_ERROR "PETSc could not be found, be sure to set PETSC_DIR")
  endif()
  find_slepc_pkgconfig("${PETSC_TARGET};${PETSC_TEST_DEPS}" SLEPC_TARGET)
  if("${SLEPC_TARGET}" STREQUAL "")
    message(FATAL_ERROR "SLEPc could not be found, be sure to set SLEPC_DIR")
  endif()
elseif(NOT PALACE_WITH_ARPACK)
  message(FATAL_ERROR "Build requires at least one of ARPACK or SLEPc dependencies")
endif()

# Find ARPACK
if(PALACE_WITH_ARPACK)
  if(NOT "${ARPACK_DIR}" STREQUAL "")
    set(arpackng_ROOT ${ARPACK_DIR})
  endif()
  find_package(arpackng REQUIRED CONFIG)
  message(STATUS "Found ARPACK: ${arpackng_VERSION} in ${arpackng_DIR}")
elseif(NOT PALACE_WITH_SLEPC)
  message(FATAL_ERROR "Build requires at least one of ARPACK or SLEPc dependencies")
endif()

# Optionally configure static analysis
include(StaticAnalysisHelpers)
if(ANALYZE_SOURCES_CLANG_TIDY)
  configure_clang_tidy()
else()
  message(STATUS "Static analysis with clang-tidy not requested")
endif()
if(ANALYZE_SOURCES_CPPCHECK)
  configure_cppcheck()
else()
  message(STATUS "Static analysis with cppcheck not requested")
endif()

# Add library target
set(LIB_TARGET_NAME libpalace)
add_library(${LIB_TARGET_NAME} "")
target_include_directories(${LIB_TARGET_NAME} PUBLIC ${CMAKE_SOURCE_DIR})
target_include_directories(${LIB_TARGET_NAME} PUBLIC ${CMAKE_BINARY_DIR}/generated)

# Add source files
add_subdirectory(drivers)
add_subdirectory(fem)
add_subdirectory(linalg)
add_subdirectory(models)
add_subdirectory(utils)

# Add executable target
set(TARGET_NAME palace)
add_executable(${TARGET_NAME} ${CMAKE_SOURCE_DIR}/main.cpp)
target_link_libraries(${TARGET_NAME} PRIVATE ${LIB_TARGET_NAME})

# Add binary extension for build architecture
if(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm")
  set(TARGET_EXTENSION "arm64")
else()
  set(TARGET_EXTENSION "x86_64")
endif()
set_target_properties(${TARGET_NAME}
  PROPERTIES
  OUTPUT_NAME "${TARGET_NAME}-${TARGET_EXTENSION}"
  SUFFIX ".bin"
)
set_target_properties(${LIB_TARGET_NAME}
  PROPERTIES OUTPUT_NAME "${TARGET_NAME}"
)

# Handle device source code
#
# TARGET_SOURCES_DEVICE defines a list of files that contain device code.
# We manually specify this list because we cannot use compilers like nvcc for
# everything, since they tend to be slower and do not support everything we need.
if(NOT "${TARGET_SOURCES_DEVICE}" STREQUAL "" AND (PALACE_WITH_CUDA OR PALACE_WITH_HIP))
  if(PALACE_WITH_CUDA)
    set(LANGUAGE_PROPERTY CUDA)
  elseif(PALACE_WITH_HIP)
    set(LANGUAGE_PROPERTY HIP)
  endif()
  set(COMPILE_OPTIONS_PROPERTY "-Wno-pedantic")
  if(PALACE_WITH_OPENMP)
    set(COMPILE_OPTIONS_PROPERTY "${COMPILE_OPTIONS_PROPERTY} ${OpenMP_CXX_FLAGS}")
  endif()
  set_property(
    SOURCE ${TARGET_SOURCES_DEVICE}
    PROPERTY LANGUAGE ${LANGUAGE_PROPERTY}
  )
  set_property(
    SOURCE ${TARGET_SOURCES_DEVICE}
    APPEND PROPERTY COMPILE_OPTIONS "${COMPILE_OPTIONS_PROPERTY}"
  )
endif()

if(PALACE_WITH_CUDA)
  if (BUILD_SHARED_LIBS)
    target_link_libraries(${LIB_TARGET_NAME} PUBLIC CUDA::cudart CUDA::cublas)
  else()
    target_link_libraries(${LIB_TARGET_NAME} PUBLIC CUDA::cudart_static CUDA::cublas_static)
  endif()
endif()

# Add JIT source file path definition for libCEED
set_property(
  SOURCE ${CMAKE_SOURCE_DIR}/main.cpp
  APPEND PROPERTY COMPILE_DEFINITIONS "PALACE_LIBCEED_JIT_SOURCE;PALACE_LIBCEED_JIT_SOURCE_DIR=\"${CMAKE_INSTALL_PREFIX}/include/palace/\""
)

# Generate embedded schema header from JSON schema files
include(${CMAKE_SOURCE_DIR}/../cmake/EmbedSchema.cmake)

# Add Git revision information (forces reconfigure when Git status changes)
include(GetGitDescription)
git_describe(GIT_COMMIT_ID)
message(STATUS "Git string: ${GIT_COMMIT_ID}")
if(NOT GIT_COMMIT_ID MATCHES "NOTFOUND")
  set_property(
    SOURCE ${CMAKE_SOURCE_DIR}/main.cpp
    APPEND PROPERTY COMPILE_DEFINITIONS "PALACE_GIT_COMMIT;PALACE_GIT_COMMIT_ID=\"${GIT_COMMIT_ID}\""
  )
endif()

# Check C++ compiler support for constexpr std::sqrt and std::filesystem
include(CheckCompilerFeatureSupport)
if(NOT DEFINED STD_FS_LIBRARIES_CACHE)
  check_std_fs_support(STD_FS_SUPPORT STD_FS_LIBRARIES)
  if(NOT STD_FS_SUPPORT)
    message(FATAL_ERROR "Could not compile a C++ program using std::filesystem")
  endif()
  set(STD_FS_LIBRARIES_CACHE ${STD_FS_LIBRARIES} CACHE INTERNAL "")
endif()
if(NOT "${STD_FS_LIBRARIES_CACHE}" STREQUAL "")
  target_link_libraries(${LIB_TARGET_NAME}
    PUBLIC ${STD_FS_LIBRARIES_CACHE}
  )
endif()

# Link with third-party dependencies
if(PALACE_WITH_SLEPC)
  target_link_libraries(${LIB_TARGET_NAME}
    PUBLIC ${SLEPC_TARGET} ${PETSC_TARGET}
  )
  target_compile_definitions(${LIB_TARGET_NAME}
    PUBLIC PALACE_WITH_SLEPC
  )
endif()
if(PALACE_WITH_ARPACK)
  target_link_libraries(${LIB_TARGET_NAME}
    PUBLIC PARPACK::PARPACK ARPACK::ARPACK ${MPI_Fortran_LIBRARIES}
  )
  target_compile_definitions(${LIB_TARGET_NAME}
    PUBLIC PALACE_WITH_ARPACK
  )
endif()
if(PALACE_WITH_OPENMP)
  target_link_libraries(${LIB_TARGET_NAME}
    PUBLIC OpenMP::OpenMP_CXX
  )
endif()
target_link_libraries(${LIB_TARGET_NAME}
  PUBLIC ${MFEM_LIBRARY} ${LIBCEED_TARGET} nlohmann_json::nlohmann_json
         nlohmann_json_schema_validator::validator fmt::fmt scn::scn
         Eigen3::Eigen LAPACK::LAPACK MPI::MPI_CXX
)

# Add MFEM include directory for Spack builds
if(NOT "${MFEM_DIR}" STREQUAL "")
  target_include_directories(${LIB_TARGET_NAME} PUBLIC ${MFEM_DIR}/include)
endif()
if(PALACE_WITH_CUDA OR PALACE_WITH_HIP)
  target_link_libraries(${LIB_TARGET_NAME}
    PUBLIC umpire camp
  )
endif()

# Set package roots from *_DIR variables
if(NOT "${HYPRE_DIR}" STREQUAL "")
  set(HYPRE_ROOT ${HYPRE_DIR})
endif()
if(NOT "${STRUMPACK_DIR}" STREQUAL "")
  set(strumpack_ROOT ${STRUMPACK_DIR})
endif()
if(NOT "${METIS_DIR}" STREQUAL "")
  set(METIS_ROOT ${METIS_DIR})
endif()
if(NOT "${PARMETIS_DIR}" STREQUAL "")
  set(ParMETIS_ROOT ${PARMETIS_DIR})
endif()
if(NOT "${SUPERLU_DIST_DIR}" STREQUAL "")
  set(superlu_dist_ROOT ${SUPERLU_DIST_DIR})
endif()
if(NOT "${MUMPS_DIR}" STREQUAL "")
  set(MUMPS_ROOT ${MUMPS_DIR})
endif()
if(NOT "${SUNDIALS_DIR}" STREQUAL "")
  set(SUNDIALS_ROOT ${SUNDIALS_DIR})
endif()
if(NOT "${GSLIB_DIR}" STREQUAL "")
  set(GSLIB_ROOT ${GSLIB_DIR})
endif()
if(NOT "${SCALAPACK_DIR}" STREQUAL "")
  set(SCALAPACK_ROOT ${SCALAPACK_DIR})
endif()

# Find METIS and ParMETIS early (needed for MUMPS/STRUMPACK/SuperLU CMake configs)
# but link them later for correct static library ordering
find_package(METIS QUIET)
if(NOT TARGET METIS::METIS)
  find_library(METIS_LIBRARY NAMES metis HINTS ${METIS_ROOT}/lib REQUIRED)
  find_path(METIS_INCLUDE_DIR NAMES metis.h HINTS ${METIS_ROOT}/include REQUIRED)
  add_library(METIS::METIS UNKNOWN IMPORTED)
  set_target_properties(METIS::METIS PROPERTIES
    IMPORTED_LOCATION "${METIS_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${METIS_INCLUDE_DIR}")
  message(STATUS "Found METIS: ${METIS_LIBRARY}")
else()
  message(STATUS "Found METIS: ${METIS_VERSION}")
endif()

find_package(ParMETIS QUIET)
if(NOT TARGET ParMETIS::ParMETIS)
  find_library(PARMETIS_LIBRARY NAMES parmetis HINTS ${ParMETIS_ROOT}/lib REQUIRED)
  find_path(PARMETIS_INCLUDE_DIR NAMES parmetis.h HINTS ${ParMETIS_ROOT}/include REQUIRED)
  add_library(ParMETIS::ParMETIS UNKNOWN IMPORTED)
  set_target_properties(ParMETIS::ParMETIS PROPERTIES
    IMPORTED_LOCATION "${PARMETIS_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${PARMETIS_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "METIS::METIS")
  message(STATUS "Found ParMETIS: ${PARMETIS_LIBRARY}")
else()
  message(STATUS "Found ParMETIS: ${ParMETIS_VERSION}")
endif()

# Find ScaLAPACK early (needed for MUMPS/STRUMPACK CMake configs)
if(PALACE_WITH_STRUMPACK OR PALACE_WITH_MUMPS)
  find_package(scalapack QUIET)
  if(NOT TARGET SCALAPACK::SCALAPACK AND NOT TARGET scalapack)
    find_library(SCALAPACK_LIBRARY NAMES scalapack mkl_scalapack_lp64 mkl_scalapack_ilp64
      HINTS ${SCALAPACK_ROOT}/lib ${CMAKE_INSTALL_PREFIX}/lib REQUIRED)
    add_library(SCALAPACK::SCALAPACK UNKNOWN IMPORTED)
    set_target_properties(SCALAPACK::SCALAPACK PROPERTIES
      IMPORTED_LOCATION "${SCALAPACK_LIBRARY}")
    message(STATUS "Found ScaLAPACK: ${SCALAPACK_LIBRARY}")
  elseif(TARGET scalapack AND NOT TARGET SCALAPACK::SCALAPACK)
    add_library(SCALAPACK::SCALAPACK ALIAS scalapack)
    message(STATUS "Found ScaLAPACK: scalapack target")
  else()
    message(STATUS "Found ScaLAPACK: SCALAPACK::SCALAPACK target")
  endif()
endif()

# Link ZLIB
find_package(ZLIB)
if(ZLIB_FOUND)
  message(STATUS "Found ZLIB: ${ZLIB_VERSION_STRING}")
  target_link_libraries(${LIB_TARGET_NAME} PUBLIC ZLIB::ZLIB)
endif()

# Link HYPRE
find_package(HYPRE QUIET)
if(TARGET HYPRE::HYPRE)
  message(STATUS "Found HYPRE: ${HYPRE_VERSION}")
  target_link_libraries(${LIB_TARGET_NAME} PUBLIC HYPRE::HYPRE)
else()
  find_library(HYPRE_LIBRARY NAMES HYPRE HINTS ${HYPRE_ROOT}/lib REQUIRED)
  find_path(HYPRE_INCLUDE_DIR NAMES HYPRE.h HINTS ${HYPRE_ROOT}/include REQUIRED)
  message(STATUS "Found HYPRE: ${HYPRE_LIBRARY}")
  target_link_libraries(${LIB_TARGET_NAME} PUBLIC ${HYPRE_LIBRARY})
  target_include_directories(${LIB_TARGET_NAME} PUBLIC ${HYPRE_INCLUDE_DIR})
endif()

# Link SuperLU_DIST
if(PALACE_WITH_SUPERLU)
  find_package(superlu_dist QUIET)
  if(TARGET superlu_dist::superlu_dist)
    message(STATUS "Found SuperLU_DIST: ${superlu_dist_VERSION}")
    target_link_libraries(${LIB_TARGET_NAME} PUBLIC superlu_dist::superlu_dist)
  else()
    find_library(SUPERLU_DIST_LIBRARY NAMES superlu_dist HINTS ${superlu_dist_ROOT}/lib REQUIRED)
    find_path(SUPERLU_DIST_INCLUDE_DIR NAMES superlu_defs.h HINTS ${superlu_dist_ROOT}/include REQUIRED)
    message(STATUS "Found SuperLU_DIST: ${SUPERLU_DIST_LIBRARY}")
    target_link_libraries(${LIB_TARGET_NAME} PUBLIC ${SUPERLU_DIST_LIBRARY})
    target_include_directories(${LIB_TARGET_NAME} PUBLIC ${SUPERLU_DIST_INCLUDE_DIR})
  endif()
endif()

# Link STRUMPACK
if(PALACE_WITH_STRUMPACK)
  find_package(strumpack QUIET)
  if(TARGET strumpack::strumpack)
    message(STATUS "Found STRUMPACK: ${strumpack_VERSION}")
    target_link_libraries(${LIB_TARGET_NAME} PUBLIC strumpack::strumpack)
  else()
    find_library(STRUMPACK_LIBRARY NAMES strumpack HINTS ${strumpack_ROOT}/lib REQUIRED)
    find_path(STRUMPACK_INCLUDE_DIR NAMES StrumpackSparseSolverMPIDist.hpp HINTS ${strumpack_ROOT}/include REQUIRED)
    message(STATUS "Found STRUMPACK: ${STRUMPACK_LIBRARY}")
    target_link_libraries(${LIB_TARGET_NAME} PUBLIC ${STRUMPACK_LIBRARY})
    target_include_directories(${LIB_TARGET_NAME} PUBLIC ${STRUMPACK_INCLUDE_DIR})
  endif()

  # STRUMPACK dependencies (ScaLAPACK target created earlier)
  target_link_libraries(${LIB_TARGET_NAME} PUBLIC SCALAPACK::SCALAPACK ${MPI_Fortran_LIBRARIES})

  find_package(ButterflyPACK QUIET)
  if(ButterflyPACK_FOUND)
    message(STATUS "Found ButterflyPACK")
    target_link_libraries(${LIB_TARGET_NAME} PUBLIC
      ButterflyPACK::dbutterflypack ButterflyPACK::zbutterflypack
      ButterflyPACK::sbutterflypack ButterflyPACK::cbutterflypack)
  endif()

  foreach(pkg zfp slate lapackpp blaspp ptscotch)
    find_package(${pkg} QUIET)
    if(${pkg}_FOUND)
      message(STATUS "Found ${pkg}")
      target_link_libraries(${LIB_TARGET_NAME} PUBLIC ${pkg})
    endif()
  endforeach()
endif()

# Link MUMPS
if(PALACE_WITH_MUMPS)
  find_package(MUMPS QUIET)
  if(TARGET MUMPS::MUMPS)
    message(STATUS "Found MUMPS: ${MUMPS_VERSION}")
    target_link_libraries(${LIB_TARGET_NAME} PUBLIC MUMPS::MUMPS)
  else()
    find_library(DMUMPS_LIBRARY NAMES dmumps HINTS ${MUMPS_ROOT}/lib REQUIRED)
    find_library(MUMPS_COMMON_LIBRARY NAMES mumps_common HINTS ${MUMPS_ROOT}/lib REQUIRED)
    find_library(PORD_LIBRARY NAMES pord HINTS ${MUMPS_ROOT}/lib REQUIRED)
    find_path(MUMPS_INCLUDE_DIR NAMES dmumps_c.h HINTS ${MUMPS_ROOT}/include REQUIRED)
    message(STATUS "Found MUMPS: ${DMUMPS_LIBRARY}")
    target_link_libraries(${LIB_TARGET_NAME} PUBLIC
      ${DMUMPS_LIBRARY} ${MUMPS_COMMON_LIBRARY} ${PORD_LIBRARY}
      SCALAPACK::SCALAPACK ${MPI_Fortran_LIBRARIES})
    target_include_directories(${LIB_TARGET_NAME} PUBLIC ${MUMPS_INCLUDE_DIR})
  endif()
endif()

# Link METIS and ParMETIS (must come after SuperLU_DIST, STRUMPACK, MUMPS for static linking)
target_link_libraries(${LIB_TARGET_NAME} PUBLIC ParMETIS::ParMETIS METIS::METIS)

# Link SUNDIALS
if(PALACE_WITH_SUNDIALS)
  # SUNDIALS needs MPI_C, but MPI is not a valid target because we only
  # specify CXX as language for palace.
  if(NOT TARGET MPI::MPI_C AND TARGET MPI::MPI_CXX)
    add_library(MPI::MPI_C ALIAS MPI::MPI_CXX)
  endif()
  find_package(SUNDIALS REQUIRED)
  message(STATUS "Found SUNDIALS: ${SUNDIALS_VERSION}")
  # Check if we are building static or shared libraries
  if(TARGET SUNDIALS::arkode)
    set(_sundials_suffix "")
  else()
    set(_sundials_suffix "_static")
  endif()
  target_link_libraries(${LIB_TARGET_NAME} PUBLIC
    SUNDIALS::arkode${_sundials_suffix} SUNDIALS::cvodes${_sundials_suffix}
    SUNDIALS::kinsol${_sundials_suffix} SUNDIALS::nvecparallel${_sundials_suffix}
    SUNDIALS::nvecmpiplusx${_sundials_suffix})
  if(PALACE_WITH_CUDA)
    target_link_libraries(${LIB_TARGET_NAME} PUBLIC SUNDIALS::nveccuda${_sundials_suffix})
  endif()
endif()

# Link GSLIB
if(PALACE_WITH_GSLIB)
  find_library(GSLIB_LIBRARY NAMES gs gslib HINTS ${GSLIB_ROOT}/lib REQUIRED)
  message(STATUS "Found GSLIB: ${GSLIB_LIBRARY}")
  target_link_libraries(${LIB_TARGET_NAME} PUBLIC ${GSLIB_LIBRARY})
endif()

target_link_libraries(
  ${LIB_TARGET_NAME}
  PRIVATE $<BUILD_INTERFACE:$<$<CONFIG:Debug,RelWithDebInfo>:common_warnings>>
          $<$<BOOL:${PALACE_BUILD_WITH_COVERAGE}>:coverage_flags>)
if(PALACE_WITH_GSLIB)
  target_compile_definitions(${LIB_TARGET_NAME}
    PUBLIC PALACE_WITH_GSLIB
  )
endif()
if(PALACE_WITH_GPU_AWARE_MPI)
  target_compile_definitions(${LIB_TARGET_NAME}
    PUBLIC PALACE_WITH_GPU_AWARE_MPI
  )
endif()


# Install target and helper scripts
install(
  TARGETS ${TARGET_NAME} ${LIB_TARGET_NAME}
  RUNTIME DESTINATION bin
)
install(
  FILES ${CMAKE_SOURCE_DIR}/../scripts/palace ${CMAKE_SOURCE_DIR}/../scripts/validate-config
  DESTINATION bin
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE
              WORLD_READ WORLD_EXECUTE
)
install(
  DIRECTORY ${CMAKE_SOURCE_DIR}/../scripts/schema
  DESTINATION bin
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ
)
install(
  DIRECTORY ${CMAKE_SOURCE_DIR}/fem/qfunctions
  DESTINATION include/palace
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ
)

# Add tests (disabled by default)
add_subdirectory(../test/unit ${CMAKE_BINARY_DIR}/test/unit EXCLUDE_FROM_ALL)

# Status messages for build settings
message(STATUS "CMake build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Building for architecture: ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "Summary of extra compiler flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "Installation directory: ${CMAKE_INSTALL_PREFIX}/bin")
