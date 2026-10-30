# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# CMake configuration for the main Palace application
#

# CMake 3.27 was released in Jul. 2023 (required for AOCL in FindBLAS/FindLAPACK)
cmake_minimum_required(VERSION 3.27)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Prohibit in-source builds
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR "In-source builds are prohibited")
endif()

# C++ standard. Override via -DCMAKE_CXX_STANDARD=<VALUE>
if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 20)
endif()
if(CMAKE_CXX_STANDARD LESS 20)
  message(FATAL_ERROR "Palace requires C++20 or later (got ${CMAKE_CXX_STANDARD})")
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Initialize the project
project(palace LANGUAGES CXX VERSION 0.18.1)

# Define build settings and defaults. Palace uses OpenMP, CUDA, or HIP when MFEM does.
set(PALACE_WITH_SLEPC ON CACHE BOOL "Build with SLEPc eigenvalue solver")
set(PALACE_WITH_ARPACK OFF CACHE BOOL "Build with ARPACK eigenvalue solver")

set(ANALYZE_SOURCES_CLANG_TIDY OFF CACHE BOOL "Run static analysis checks using clang-tidy")
set(ANALYZE_SOURCES_CPPCHECK OFF CACHE BOOL "Run static analysis checks using cppcheck")

set(PALACE_BUILD_WITH_COVERAGE OFF CACHE BOOL "Compile Palace with coverage flags (source-based for LLVM; gcov for GCC)")
set(PALACE_BUILD_WITH_SANITIZERS OFF CACHE BOOL "Compile Palace with sanitizers (AddressSanitizer + UndefinedBehaviorSanitizer)")
set(PALACE_WITH_ENZYME OFF CACHE BOOL "Build Enzyme support for Palace unit tests")

# Help find third-party dependencies
set(MFEM_DIR "" CACHE STRING
  "Path to MFEM installation or CMake package directory (not required if already on CMAKE_PREFIX_PATH)"
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
set(ENZYME_DIR "" CACHE STRING
  "Path to Enzyme build or installation directory (not required if already on CMAKE_PREFIX_PATH)"
)

# Find MFEM. The mfem target also links HYPRE, METIS, the sparse direct solvers, SUNDIALS,
# GSLIB, Umpire, and the other libraries MFEM was built with. Recent MFEM targets link to
# target Threads::Threads, for some reason.
find_package(Threads REQUIRED)
find_package(MFEM REQUIRED CONFIG HINTS "${MFEM_DIR}")
message(STATUS "Found MFEM: ${MFEM_VERSION} in ${MFEM_DIR}")
if(NOT MFEM_USE_MPI)
  message(FATAL_ERROR "Build requires MFEM with MPI support")
endif()

# Enable Fortran if required
if(PALACE_WITH_ARPACK)
  enable_language(Fortran)
endif()

# Enable CUDA/HIP if required
if(MFEM_USE_CUDA)
  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)
  # C++20 device code: nvcc supports -std=c++20 from 12.0
  if(CMAKE_CUDA_COMPILER_ID STREQUAL "NVIDIA" AND
     CMAKE_CUDA_COMPILER_VERSION VERSION_LESS 12.0)
    message(FATAL_ERROR
      "CUDA builds require CUDA 12.0 or later (found ${CMAKE_CUDA_COMPILER_VERSION})")
  endif()
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-extended-lambda")
elseif(MFEM_USE_HIP)
  enable_language(HIP)
endif()

# Device translation units share headers with host code, so compile them with the same C++
# standard (otherwise nvcc/hipcc use their own default)
set(CMAKE_CUDA_STANDARD ${CMAKE_CXX_STANDARD})
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_EXTENSIONS OFF)
set(CMAKE_HIP_STANDARD ${CMAKE_CXX_STANDARD})
set(CMAKE_HIP_STANDARD_REQUIRED ON)
set(CMAKE_HIP_EXTENSIONS OFF)

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
        $<$<CXX_COMPILER_ID:GNU>:-fprofile-dir=${CMAKE_BINARY_DIR}/coverage.%p>
  )
  target_link_options(
    coverage_flags
    INTERFACE
        $<$<CXX_COMPILER_ID:AppleClang,Clang,IntelLLVM>:-fprofile-instr-generate>
        $<$<CXX_COMPILER_ID:AppleClang,Clang,IntelLLVM>:-fcoverage-mapping>
        $<$<CXX_COMPILER_ID:GNU>:--coverage>
  )
endif()

if(PALACE_BUILD_WITH_SANITIZERS)
  if(CMAKE_SYSTEM_NAME STREQUAL "Darwin" AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64"
     AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(FATAL_ERROR
      "PALACE_BUILD_WITH_SANITIZERS is not compatible with GCC on macOS Apple Silicon. "
      "GCC does not ship libasan on this platform. Use Clang instead "
      "(e.g. -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++)."
    )
  endif()
  add_library(sanitizer_flags INTERFACE)
  target_compile_options(sanitizer_flags INTERFACE -fsanitize=address,undefined -fno-omit-frame-pointer)
  target_link_options(sanitizer_flags INTERFACE -fsanitize=address,undefined)
endif()

# Add extra CMake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Find MPI
find_package(MPI REQUIRED)
if(PALACE_WITH_ARPACK)
  find_package(MPI REQUIRED COMPONENTS Fortran)
endif()

# Find OpenMP
if(MFEM_USE_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

# Find LAPACK
find_package(LAPACK REQUIRED)

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

# Find Enzyme for unit tests
if(PALACE_WITH_ENZYME)
  if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    message(FATAL_ERROR
      "PALACE_WITH_ENZYME requires upstream Clang built against the same LLVM as Enzyme "
      "(got ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION})"
    )
  endif()
  find_package(Enzyme REQUIRED CONFIG HINTS "${ENZYME_DIR}")
  set(ENZYME_LLVM_VERSION
    "${Enzyme_LLVM_VERSION_MAJOR}.${Enzyme_LLVM_VERSION_MINOR}.${Enzyme_LLVM_VERSION_PATCH}"
  )
  if(NOT TARGET ClangEnzymeFlags OR
     NOT TARGET ClangEnzyme-${Enzyme_LLVM_VERSION_MAJOR})
    message(FATAL_ERROR
      "The Enzyme package does not provide its Clang compiler plugin targets. "
      "Rebuild Enzyme with the matching Clang development package."
    )
  endif()
  if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL ENZYME_LLVM_VERSION)
    message(FATAL_ERROR
      "Enzyme requires matching Clang and LLVM versions: compiler is "
      "${CMAKE_CXX_COMPILER_VERSION}, Enzyme uses LLVM ${ENZYME_LLVM_VERSION}"
    )
  endif()
  message(STATUS "Found Enzyme for LLVM ${ENZYME_LLVM_VERSION} in ${Enzyme_DIR}")
endif()

# Find libCEED
include(PkgConfigHelpers)
set(GPU_TEST_DEPS)
if(MFEM_USE_CUDA)
  if(BUILD_SHARED_LIBS)
    list(APPEND GPU_TEST_DEPS CUDA::cudart CUDA::cublas CUDA::cuda_driver)
  else()
    list(APPEND GPU_TEST_DEPS CUDA::cudart_static CUDA::cublas_static CUDA::cuda_driver)
  endif()
elseif(MFEM_USE_HIP)
  list(APPEND GPU_TEST_DEPS hip::host)
endif()
set(LIBCEED_TEST_DEPS ${GPU_TEST_DEPS})
if(MFEM_USE_OPENMP)
  list(APPEND LIBCEED_TEST_DEPS OpenMP::OpenMP_CXX)
endif()
find_libceed_pkgconfig("${LIBCEED_TEST_DEPS}" LIBCEED_TARGET)
if("${LIBCEED_TARGET}" STREQUAL "")
  message(FATAL_ERROR "libCEED could not be found, be sure to set LIBCEED_DIR")
endif()

# Find PETSc and SLEPc
if(PALACE_WITH_SLEPC)
  set(PETSC_TEST_DEPS MPI::MPI_CXX ${GPU_TEST_DEPS})
  if(MFEM_USE_OPENMP)
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

target_sources(${LIB_TARGET_NAME}
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/driver.cpp
)

# Add executable target
set(TARGET_NAME palace)
add_executable(${TARGET_NAME} ${CMAKE_SOURCE_DIR}/main.cpp)
target_link_libraries(${TARGET_NAME}
  PRIVATE ${LIB_TARGET_NAME}
          $<$<BOOL:${PALACE_BUILD_WITH_SANITIZERS}>:sanitizer_flags>)

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
if(NOT "${TARGET_SOURCES_DEVICE}" STREQUAL "" AND (MFEM_USE_CUDA OR MFEM_USE_HIP))
  if(MFEM_USE_CUDA)
    set(LANGUAGE_PROPERTY CUDA)
  elseif(MFEM_USE_HIP)
    set(LANGUAGE_PROPERTY HIP)
  endif()
  set(COMPILE_OPTIONS_PROPERTY "-Wno-pedantic")
  if(MFEM_USE_OPENMP)
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

if(MFEM_USE_CUDA)
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

# Add Git revision information (forces reconfigure when Git status changes). Builds without
# usable Git metadata (e.g. Spack) can supply it via PALACE_GIT_COMMIT_ID, as a CMake variable
# or, failing that, an environment variable.
set(PALACE_GIT_COMMIT_ID "" CACHE STRING "Version string to embed instead of git describe")
if(PALACE_GIT_COMMIT_ID)
  set(GIT_COMMIT_ID "${PALACE_GIT_COMMIT_ID}")
elseif(NOT "$ENV{PALACE_GIT_COMMIT_ID}" STREQUAL "")
  set(GIT_COMMIT_ID "$ENV{PALACE_GIT_COMMIT_ID}")
else()
  include(GetGitDescription)
  git_describe(GIT_COMMIT_ID)
endif()
message(STATUS "Git string: ${GIT_COMMIT_ID}")
if(NOT GIT_COMMIT_ID MATCHES "NOTFOUND")
  set_property(
    SOURCE ${CMAKE_SOURCE_DIR}/main.cpp
    APPEND PROPERTY COMPILE_DEFINITIONS "PALACE_GIT_COMMIT;PALACE_GIT_COMMIT_ID=\"${GIT_COMMIT_ID}\""
  )
endif()

# Identifier of the producing build (e.g. the Spack DAG hash), printed by --version
set(PALACE_BUILD_ID "" CACHE STRING "Build identifier to embed (e.g. Spack DAG hash)")
if(PALACE_BUILD_ID)
  set_property(
    SOURCE ${CMAKE_SOURCE_DIR}/main.cpp
    APPEND PROPERTY COMPILE_DEFINITIONS "PALACE_BUILD_ID=\"${PALACE_BUILD_ID}\""
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
if(MFEM_USE_OPENMP)
  target_link_libraries(${LIB_TARGET_NAME}
    PUBLIC OpenMP::OpenMP_CXX
  )
endif()
target_link_libraries(${LIB_TARGET_NAME}
  PUBLIC mfem ${LIBCEED_TARGET} nlohmann_json::nlohmann_json
         nlohmann_json_schema_validator::validator fmt::fmt scn::scn
         Eigen3::Eigen LAPACK::LAPACK MPI::MPI_CXX
)
if(MFEM_USE_CUDA)
  # libCEED uses the CUDA driver API. FindCUDAToolkit resolves this target to
  # the toolkit stub at build time without adding the stub directory to RPATH;
  # the NVIDIA driver provides libcuda.so.1 at runtime.
  target_link_libraries(${LIB_TARGET_NAME} PUBLIC CUDA::cuda_driver)
endif()

# Optional STRUMPACK dependencies, linked for their version queries
if(MFEM_USE_STRUMPACK)
  find_package(ButterflyPACK QUIET)
  foreach(pkg zfp slate lapackpp blaspp)
    find_package(${pkg} QUIET)
    if(${pkg}_FOUND)
      message(STATUS "Found ${pkg}")
      if(TARGET ${pkg}::${pkg})
        target_link_libraries(${LIB_TARGET_NAME} PUBLIC ${pkg}::${pkg})
      else()
        target_link_libraries(${LIB_TARGET_NAME} PUBLIC ${pkg})
      endif()
    endif()
  endforeach()
endif()

# Only utils/gslibversion.cpp includes the GSLIB header, for the version query
if(MFEM_USE_GSLIB)
  target_sources(${LIB_TARGET_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/utils/gslibversion.cpp)
endif()

target_link_libraries(
  ${LIB_TARGET_NAME}
  PRIVATE $<BUILD_INTERFACE:$<$<CONFIG:Debug,RelWithDebInfo>:common_warnings>>
          $<$<BOOL:${PALACE_BUILD_WITH_COVERAGE}>:coverage_flags>
          $<$<BOOL:${PALACE_BUILD_WITH_SANITIZERS}>:sanitizer_flags>)

# Tell utils/versions.cpp which optional dependencies are linked. Versions of dependencies
# without a header or runtime query that Palace can reach are taken from CMake.
set(VERSION_DEFINITIONS)
foreach(pkg zfp slate lapackpp blaspp)
  if(MFEM_USE_STRUMPACK AND ${pkg}_FOUND)
    string(TOUPPER ${pkg} _pkg)
    list(APPEND VERSION_DEFINITIONS PALACE_WITH_${_pkg})
  endif()
endforeach()
# The BLAS/LAPACK vendor selects which vendor-specific version query to compile in. Match
# BLA_VENDOR and library file names only, since directory names (e.g. hashes) are arbitrary.
set(_blas_lapack "${BLA_VENDOR}")
foreach(_lib IN LISTS LAPACK_LIBRARIES BLAS_LIBRARIES)
  get_filename_component(_lib "${_lib}" NAME)
  list(APPEND _blas_lapack "${_lib}")
endforeach()
string(TOLOWER "${_blas_lapack}" _blas_lapack)
if(_blas_lapack MATCHES "(^|;)(lib)?openblas")
  list(APPEND VERSION_DEFINITIONS PALACE_BLAS_OPENBLAS)
  # ILP64 OpenBLAS built with SYMBOLSUFFIX=64_ is named libopenblas64_
  if(_blas_lapack MATCHES "(^|;)(lib)?openblas[op]?64_")
    list(APPEND VERSION_DEFINITIONS PALACE_BLAS_OPENBLAS64)
  endif()
elseif(_blas_lapack MATCHES "(^|;)((lib)?mkl_|intel10)")
  list(APPEND VERSION_DEFINITIONS PALACE_BLAS_MKL)
elseif(_blas_lapack MATCHES "(^|;)((lib)?armpl|arm)")
  list(APPEND VERSION_DEFINITIONS PALACE_BLAS_ARMPL)
elseif(_blas_lapack MATCHES "(^|;)((lib)?blis|aocl)")
  list(APPEND VERSION_DEFINITIONS PALACE_BLAS_AOCL)
endif()
if(nlohmann_json_schema_validator_VERSION)
  list(APPEND VERSION_DEFINITIONS
    "PALACE_JSON_SCHEMA_VALIDATOR_VERSION=\"${nlohmann_json_schema_validator_VERSION}\"")
endif()
if(PALACE_WITH_ARPACK AND arpackng_VERSION)
  list(APPEND VERSION_DEFINITIONS "PALACE_ARPACK_VERSION=\"${arpackng_VERSION}\"")
endif()
if(MFEM_USE_STRUMPACK AND ButterflyPACK_FOUND AND ButterflyPACK_VERSION)
  list(APPEND VERSION_DEFINITIONS
    "PALACE_BUTTERFLYPACK_VERSION=\"${ButterflyPACK_VERSION}\"")
endif()
# libxsmm and MAGMA are only linked by libCEED, so check its pkg-config dependencies. Query
# pkg-config directly, since pkg_check_modules caches its results across reconfigures.
if(libCEED_STATIC_LIBRARIES MATCHES "(^|;)(:lib)?xsmm(\\.a)?(;|$)")
  execute_process(
    COMMAND ${PKG_CONFIG_EXECUTABLE} --modversion libxsmm
    OUTPUT_VARIABLE _libxsmm_version OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _result ERROR_QUIET
  )
  if(_result EQUAL 0 AND _libxsmm_version)
    list(APPEND VERSION_DEFINITIONS "PALACE_LIBXSMM_VERSION=\"${_libxsmm_version}\"")
  endif()
endif()
if(libCEED_STATIC_LIBRARIES MATCHES "(^|;)(:lib)?magma(\\.a)?(;|$)")
  # magma.pc carries no version, so read it from the headers
  pkg_get_variable(_magma_includedir magma includedir)
  find_file(MAGMA_TYPES_HEADER magma_types.h HINTS ${_magma_includedir})
  if(MAGMA_TYPES_HEADER)
    set(_magma_version)
    foreach(_part MAJOR MINOR MICRO)
      file(STRINGS ${MAGMA_TYPES_HEADER} _define
        REGEX "^#define MAGMA_VERSION_${_part}[ \t]+[0-9]+")
      if(NOT _define MATCHES "^#define MAGMA_VERSION_${_part}[ \t]+([0-9]+)")
        set(_magma_version)
        break()
      endif()
      list(APPEND _magma_version ${CMAKE_MATCH_1})
    endforeach()
    if(_magma_version)
      list(JOIN _magma_version "." _magma_version)
      list(APPEND VERSION_DEFINITIONS "PALACE_MAGMA_VERSION=\"${_magma_version}\"")
    endif()
  endif()
endif()
set_property(
  SOURCE ${CMAKE_SOURCE_DIR}/utils/versions.cpp
  APPEND PROPERTY COMPILE_DEFINITIONS "${VERSION_DEFINITIONS}"
)


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
enable_testing()
add_subdirectory(../test/unit ${CMAKE_BINARY_DIR}/test/unit EXCLUDE_FROM_ALL)

# Status messages for build settings
message(STATUS "CMake build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Building for architecture: ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "Summary of extra compiler flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "Installation directory: ${CMAKE_INSTALL_PREFIX}/bin")
