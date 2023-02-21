# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Add dependencies of Palace to the build
#

# Git tags and version numbers for internal builds
set(PALACE_INTERNAL_JSON_URL
  "https://github.com/nlohmann/json/releases/download/v3.11.2/json.tar.xz" CACHE STRING
  "URL for internal nlohmann/json build"
)
set(PALACE_INTERNAL_FMT_URL
  "https://github.com/fmtlib/fmt/archive/refs/tags/9.1.0.tar.gz" CACHE STRING
  "URL for internal fmt build"
)
set(PALACE_INTERNAL_MFEM_GIT_TAG
  "ea5c4664864cb84c0a453f527d6befec05a8b8a9" CACHE STRING  # master @ 02/24/2023
  "Git tag for internal MFEM build"
)

# Avoid DOWNLOAD_EXTRACT_TIMESTAMP warning
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
  cmake_policy(SET CMP0135 NEW)
endif()

# Get JSON library
if(PALACE_WITH_INTERNAL_JSON)
  message(STATUS "================= Configuring nlohmann/json dependency =================")

  # Need to manually turn on install step for JSON (for shared library builds)
  set(JSON_Install ON CACHE BOOL "")

  include(FetchContent)
  FetchContent_Declare(json
    URL ${PALACE_INTERNAL_JSON_URL}
  )
  FetchContent_MakeAvailable(json)
endif()

# Get fmt library
if(PALACE_WITH_INTERNAL_FMT)
  message(STATUS "====================== Configuring fmt dependency ======================")

  # Need to manually turn on install step for fmt (for shared library builds)
  set(FMT_INSTALL ON CACHE BOOL "")

  include(FetchContent)
  FetchContent_Declare(fmt
    URL ${PALACE_INTERNAL_FMT_URL}
  )
  FetchContent_MakeAvailable(fmt)
endif()

# Get MFEM library
if(PALACE_WITH_INTERNAL_MFEM)
  message(STATUS "===================== Configuring MFEM dependency ======================")

  # Configure the MFEM build by setting variables recognized by the MFEM CMake build
  # Paths, dependency packages, and dependency libraries should be set by cache
  set(MFEM_USE_MPI YES CACHE BOOL "")
  set(MFEM_USE_METIS_5 YES CACHE BOOL "")
  set(MFEM_USE_OPENMP ${PALACE_WITH_OPENMP} CACHE BOOL "")
  set(MFEM_USE_GSLIB ${PALACE_WITH_GSLIB} CACHE BOOL "")
  set(MFEM_USE_SUPERLU ${PALACE_WITH_SUPERLU} CACHE BOOL "")
  set(MFEM_USE_STRUMPACK ${PALACE_WITH_STRUMPACK} CACHE BOOL "")
  set(MFEM_USE_MUMPS ${PALACE_WITH_MUMPS} CACHE BOOL "")

  # Silence #pragma omp warnings when not building with OpenMP
  if(PALACE_WITH_STRUMPACK AND NOT PALACE_WITH_OPENMP)
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag(-Wno-unknown-pragmas SUPPORTS_NOPRAGMA_WARNING)
    if(SUPPORTS_NOPRAGMA_WARNING)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
    endif()
  endif()

  # Find optional MFEM dependencies with CMake because once passed to MFEM, they are
  # required
  find_package(ZLIB)
  if(ZLIB_FOUND)
    message(STATUS "Building MFEM with zlib support for binary output compression")
    set(MFEM_USE_ZLIB YES CACHE BOOL "")
  else()
    message(STATUS "Building without zlib and binary output compression support")
    set(MFEM_USE_ZLIB NO CACHE BOOL "")
  endif()
  set(MFEM_USE_LIBUNWIND NO CACHE BOOL "")

  # A number of patches to MFEM for our uses
  set(PALACE_MFEM_PATCH_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/patch/mfem/patch_integ.diff"
    "${CMAKE_CURRENT_SOURCE_DIR}/patch/mfem/patch_mesh.diff"
    "${CMAKE_CURRENT_SOURCE_DIR}/patch/mfem/patch_part.diff"
    "${CMAKE_CURRENT_SOURCE_DIR}/patch/mfem/patch_direct.diff"
  )

  include(FetchContent)
  FetchContent_Declare(mfem
    GIT_REPOSITORY https://github.com/mfem/mfem.git
    GIT_TAG        ${PALACE_INTERNAL_MFEM_GIT_TAG}
    PATCH_COMMAND
      git reset --hard && git clean -fd && git apply "${PALACE_MFEM_PATCH_FILES}"
  )
  FetchContent_MakeAvailable(mfem)
endif()

message(STATUS "================== Finished dependency configuration ===================")
