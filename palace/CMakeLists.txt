# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# CMake configuration for the Palace application
#

# CMake 3.13 was released in Nov. 2018
cmake_minimum_required(VERSION 3.13)

# Prohibit in-source builds
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are prohibited")
endif()

# C++17 required for std::filesystem, among others
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

# Initialize the project
project(palace LANGUAGES CXX VERSION 0.11.0)

# Define build settings and defaults
set(PALACE_WITH_OPENMP OFF CACHE BOOL "Use OpenMP")
set(PALACE_WITH_GSLIB ON CACHE BOOL "Build with GSLIB library for high-order field interpolation")
set(PALACE_WITH_SUPERLU ON CACHE BOOL "Build with SuperLU_DIST sparse direct solver")
set(PALACE_WITH_STRUMPACK OFF CACHE BOOL "Build with STRUMPACK sparse direct solver")
set(PALACE_WITH_MUMPS OFF CACHE BOOL "Build with MUMPS sparse direct solver")
set(PALACE_WITH_SLEPC ON CACHE BOOL "Build with SLEPc eigenvalue solver")
set(PALACE_WITH_ARPACK OFF CACHE BOOL "Build with ARPACK eigenvalue solver")
set(PALACE_WITH_INTERNAL_JSON ON CACHE BOOL "Use internal nlohmann/json library")
set(PALACE_WITH_INTERNAL_FMT ON CACHE BOOL "Use internal fmt library")
set(PALACE_WITH_INTERNAL_MFEM ON CACHE BOOL "Use internal MFEM library")
set(ANALYZE_SOURCES_CLANG_TIDY OFF CACHE BOOL "Run static analysis checks using clang-tidy")
set(ANALYZE_SOURCES_CPPCHECK OFF CACHE BOOL "Run static analysis checks using cppcheck")

# Help find third-party dependencies
set(PETSC_DIR "" CACHE STRING "Path to PETSc build or installation directory (not required if already on PKG_CONFIG_PATH)")
set(SLEPC_DIR "" CACHE STRING "Path to SLEPc build or installation directory (not required if already on PKG_CONFIG_PATH)")
set(ARPACK_DIR "" CACHE STRING "Path to ARPACK build or installation directory (not required if already on CMAKE_PREFIX_PATH)")

# Help find dependencies for the internal MFEM build
# If we trust MFEM's Find<PACKAGE>.cmake module, we can just set <PACKAGE>_DIR and, if
# needed, <PACKAGE>_REQUIRED_PACKAGES. The extra <PACKAGE>_REQUIRED_LIBRARIES can be used
# to add any additional dependency libraries.
set(PALACE_MFEM_DEPS
  "METIS"
  "HYPRE"
  "SuperLUDist"
  "STRUMPACK"
  "MUMPS"
  "GSLIB"
)
foreach(DEP IN LISTS PALACE_MFEM_DEPS)
  set(${DEP}_DIR "" CACHE STRING "Path to ${DEP} build or installation directory")
  set(${DEP}_REQUIRED_PACKAGES "" CACHE STRING "List of additional required packages for ${DEP}")
  set(${DEP}_REQUIRED_LIBRARIES "" CACHE STRING "List of additional required libraries for ${DEP}")
  # set(${DEP}_LIBRARIES "" CACHE STRING "List of library files for ${DEP}")
  # set(${DEP}_INCLUDE_DIRS "" CACHE STRING "Path to ${DEP} include directories")
endforeach()

# Enable Fortran if required
if(PALACE_WITH_STRUMPACK OR PALACE_WITH_MUMPS OR PALACE_WITH_ARPACK)
  enable_language(Fortran)
endif()

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose CMake build type" FORCE)
endif()

# Set a default for dependency library builds if none was provided
if(NOT DEFINED BUILD_SHARED_LIBS)
  message(STATUS "Setting BUILD_SHARED_LIBS to 'OFF' as it was not specified")
  set(BUILD_SHARED_LIBS OFF CACHE BOOL "Global flag to cause add_library() to create shared libraries if ON")
endif()

# Add extra CMake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../cmake")

# Add third-party dependency builds
add_subdirectory(deps)
if(NOT PALACE_WITH_INTERNAL_JSON)
  find_package(nlohmann_json REQUIRED CONFIG)
endif()
if(NOT PALACE_WITH_INTERNAL_FMT)
  find_package(fmt REQUIRED CONFIG)
endif()
if(NOT PALACE_WITH_INTERNAL_MFEM)
  message(FATAL_ERROR "Palace only supports builds with internal MFEM")
endif()

# Find MPI
find_package(MPI REQUIRED)

# Find OpenMP
if(PALACE_WITH_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

# Find PETSc and SLEPc
include(PkgConfigHelpers)
set(PETSC_TEST_DEPS MPI::MPI_CXX)
if(PALACE_WITH_OPENMP)
  list(APPEND PETSC_TEST_DEPS OpenMP::OpenMP_CXX)
endif()
find_petsc_pkgconfig("${PETSC_TEST_DEPS}" PETSC_TARGET)
if("${PETSC_TARGET}" STREQUAL "")
  message(FATAL_ERROR "PETSc could not be found, be sure to set PETSC_DIR")
endif()
if(PALACE_WITH_SLEPC)
  find_slepc_pkgconfig("${PETSC_TARGET};${PETSC_TEST_DEPS}" SLEPC_TARGET)
  if("${SLEPC_TARGET}" STREQUAL "")
    message(FATAL_ERROR "SLEPc could not be found, be sure to set SLEPC_DIR")
  endif()
endif()

# Find ARPACK
if(PALACE_WITH_ARPACK)
  set(arpackng_ROOT ${ARPACK_DIR})
  find_package(arpackng REQUIRED CONFIG)
  message(STATUS "Found ARPACK: ${arpackng_VERSION}")
elseif(NOT PALACE_WITH_SLEPC)
  message(FATAL_ERROR "Build requires at least one of ARPACK or SLEPc dependencies")
endif()

# Add executable target
set(TARGET_NAME palace)
add_executable(${TARGET_NAME} "")

# Add source files
target_sources(${TARGET_NAME}
  PRIVATE ${CMAKE_SOURCE_DIR}/main.cpp
)
target_include_directories(${TARGET_NAME}
  PRIVATE ${CMAKE_SOURCE_DIR}
)
add_subdirectory(drivers)
add_subdirectory(fem)
add_subdirectory(linalg)
add_subdirectory(utils)

# Add binary extension for build architecture
if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "aarch64|arm")
  set(TARGET_EXTENSION "arm64")
else()
  set(TARGET_EXTENSION "x86_64")
endif()
set_target_properties(${TARGET_NAME} PROPERTIES
  OUTPUT_NAME "${TARGET_NAME}-${TARGET_EXTENSION}"
  SUFFIX ".bin"
)

# Add Git revision information (forces reconfigure when Git status changes)
include(GetGitDescription)
git_describe(GIT_COMMIT_ID)
message(STATUS "Git string: ${GIT_COMMIT_ID}")
if(NOT "${GIT_COMMIT_ID}" MATCHES "NOTFOUND")
  file(WRITE ${CMAKE_BINARY_DIR}/config/gitversion.hpp
"#pragma once
const char *GetGitCommit()
{
  static const char *commit = \"@GIT_COMMIT_ID@\";
  return commit;
}
")
  configure_file(${CMAKE_BINARY_DIR}/config/gitversion.hpp ${CMAKE_BINARY_DIR}/config @ONLY)
  target_include_directories(${TARGET_NAME}
    PRIVATE ${CMAKE_BINARY_DIR}/config
  )
  target_compile_definitions(${TARGET_NAME}
    PRIVATE PALACE_GIT_COMMIT
  )
endif()

# Check C++ compiler support for constexpr std::sqrt and std::filesystem
include(CheckCompilerFeatureSupport)
if(NOT DEFINED CONSTEXPR_SQRT_SUPPORT_CACHE)
  check_constexpr_sqrt_support(CONSTEXPR_SQRT_SUPPORT)
  set(CONSTEXPR_SQRT_SUPPORT_CACHE ${CONSTEXPR_SQRT_SUPPORT} CACHE INTERNAL "")
endif()
if(CONSTEXPR_SQRT_SUPPORT_CACHE)
  target_compile_definitions(${TARGET_NAME}
    PRIVATE PALACE_WITH_CONSTEXPR_SQRT
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
  target_link_libraries(${TARGET_NAME}
    PRIVATE ${STD_FS_LIBRARIES_CACHE}
  )
endif()

# Link with third-party dependencies
if(PALACE_WITH_SLEPC)
  target_link_libraries(${TARGET_NAME}
    PRIVATE ${SLEPC_TARGET}
  )
  target_compile_definitions(${TARGET_NAME}
    PRIVATE PALACE_WITH_SLEPC
  )
endif()
if(PALACE_WITH_ARPACK)
  target_link_libraries(${TARGET_NAME}
    PRIVATE PARPACK::PARPACK ARPACK::ARPACK MPI::MPI_Fortran
  )
  target_compile_definitions(${TARGET_NAME}
    PRIVATE PALACE_WITH_ARPACK
  )
endif()
target_link_libraries(${TARGET_NAME}
  PRIVATE ${PETSC_TARGET} mfem nlohmann_json::nlohmann_json fmt::fmt MPI::MPI_CXX
)

# Install target and helper scripts
install(TARGETS ${TARGET_NAME}
  RUNTIME DESTINATION bin
)
install(FILES
  ${CMAKE_SOURCE_DIR}/../scripts/palace ${CMAKE_SOURCE_DIR}/../scripts/validate-config
  DESTINATION bin
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)
install(DIRECTORY
  ${CMAKE_SOURCE_DIR}/../scripts/schema
  DESTINATION bin
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ
)

# Optionally configure static analysis
include(StaticAnalysisHelpers)
if(ANALYZE_SOURCES_CLANG_TIDY)
  configure_clang_tidy(${TARGET_NAME})
else()
  message(STATUS "Static analysis with clang-tidy not requested")
endif()
if(ANALYZE_SOURCES_CPPCHECK)
  configure_cppcheck(${TARGET_NAME})
else()
  message(STATUS "Static analysis with cppcheck not requested")
endif()

# Status messages for build settings
message(STATUS "CMake build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Building for architecture: ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "Summary of extra compiler flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "Installation directory: ${CMAKE_INSTALL_PREFIX}/bin")
