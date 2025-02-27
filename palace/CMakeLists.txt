# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# CMake configuration for the main Palace application
#

# CMake 3.21 was released in Jul. 2021 (required for HIP support)
cmake_minimum_required(VERSION 3.21)

# Prohibit in-source builds
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR "In-source builds are prohibited")
endif()

# C++17 required for std::filesystem, among others
# C++20 for ranges, lambdas
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

# Initialize the project
project(palace LANGUAGES CXX VERSION 0.13.0)

# Define build settings and defaults
set(PALACE_WITH_OPENMP OFF CACHE BOOL "Use OpenMP for shared-memory parallelism")
set(PALACE_WITH_CUDA OFF CACHE BOOL "Use CUDA for NVIDIA GPU support")
set(PALACE_WITH_HIP OFF CACHE BOOL "Use HIP for AMD or NVIDIA GPU support")

set(PALACE_WITH_SLEPC ON CACHE BOOL "Build with SLEPc eigenvalue solver")
set(PALACE_WITH_ARPACK OFF CACHE BOOL "Build with ARPACK eigenvalue solver")

set(ANALYZE_SOURCES_CLANG_TIDY OFF CACHE BOOL "Run static analysis checks using clang-tidy")
set(ANALYZE_SOURCES_CPPCHECK OFF CACHE BOOL "Run static analysis checks using cppcheck")

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

# Enable Fortran if required
if(PALACE_WITH_ARPACK)
  enable_language(Fortran)
endif()

# Enable CUDA/HIP if required
if(PALACE_WITH_CUDA AND PALACE_WITH_HIP)
  message(FATAL_ERROR "PALACE_WITH_CUDA is not compatible with PALACE_WITH_HIP")
endif()
if(PALACE_WITH_CUDA)
  enable_language(CUDA)
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

# Add extra CMake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Find MPI
find_package(MPI REQUIRED)

# Find OpenMP
if(PALACE_WITH_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

# Find LAPACK
find_package(LAPACK REQUIRED)

# Find nlohmann/json
find_package(nlohmann_json REQUIRED CONFIG)
message(STATUS "Found nlohmann/json: ${nlohmann_json_VERSION} in ${nlohmann_json_DIR}")

# Find fmt
find_package(fmt REQUIRED CONFIG)
message(STATUS "Found fmt: ${fmt_VERSION} in ${fmt_DIR}")

# Find Eigen
find_package(Eigen3 REQUIRED CONFIG)
message(STATUS "Found Eigen: ${Eigen3_VERSION} in ${Eigen3_DIR}")

# Find MFEM (recent MFEM targets link to target Threads::Threads, for some reason)
if(NOT "${MFEM_DIR}" STREQUAL "")
  set(MFEM_ROOT ${MFEM_DIR})
endif()
find_package(Threads REQUIRED)
find_package(MFEM REQUIRED CONFIG)
message(STATUS "Found MFEM: ${MFEM_VERSION} in ${MFEM_DIR}")
if(NOT MFEM_USE_MPI)
  message(FATAL_ERROR "Build requires MFEM with MPI support")
endif()
# if(MFEM_CXX_FLAGS)
#   # Pull compiler flags from MFEM for OpenMP and optimizations
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MFEM_CXX_FLAGS}")
# endif()

# Find libCEED
include(PkgConfigHelpers)
set(LIBCEED_TEST_DEPS)
if(PALACE_WITH_OPENMP)
  list(APPEND LIBCEED_TEST_DEPS OpenMP::OpenMP_CXX)
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

# Resolve linker error with CUDA and static linkage
if(PALACE_WITH_CUDA)
  set_target_properties(${TARGET_NAME} PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
endif()

# Add JIT source file path definition for libCEED
set_property(
  SOURCE ${CMAKE_SOURCE_DIR}/main.cpp
  APPEND PROPERTY COMPILE_DEFINITIONS "PALACE_LIBCEED_JIT_SOURCE;PALACE_LIBCEED_JIT_SOURCE_DIR=\"${CMAKE_INSTALL_PREFIX}/include/palace/\""
)

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
if(NOT DEFINED CONSTEXPR_SQRT_SUPPORT_CACHE)
  check_constexpr_sqrt_support(CONSTEXPR_SQRT_SUPPORT)
  set(CONSTEXPR_SQRT_SUPPORT_CACHE ${CONSTEXPR_SQRT_SUPPORT} CACHE INTERNAL "")
endif()
if(CONSTEXPR_SQRT_SUPPORT_CACHE)
  target_compile_definitions(${LIB_TARGET_NAME}
    PUBLIC PALACE_WITH_CONSTEXPR_SQRT
  )
endif()
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
  PUBLIC mfem ${LIBCEED_TARGET} nlohmann_json::nlohmann_json fmt::fmt
         Eigen3::Eigen LAPACK::LAPACK MPI::MPI_CXX
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
add_subdirectory(../test/unit ${CMAKE_BINARY_DIR}/test/unit EXCLUDE_FROM_ALL)

# Status messages for build settings
message(STATUS "CMake build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Building for architecture: ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "Summary of extra compiler flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "Installation directory: ${CMAKE_INSTALL_PREFIX}/bin")
