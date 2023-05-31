# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Build Scotch and PTScotch to provide METIS and ParMETIS API
#

# Force build order
set(SCOTCH_DEPENDENCIES)

set(SCOTCH_OPTIONS ${PALACE_SUPERBUILD_DEFAULT_ARGS})
list(APPEND SCOTCH_OPTIONS
  "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
  "-DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}"
  "-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}"
  "-DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}"
  "-DBUILD_PTSCOTCH=ON"
  "-DBUILD_LIBESMUMPS=ON"
  "-DBUILD_LIBSCOTCHMETIS=ON"
  "-DINSTALL_METIS_HEADERS=ON"
)

# Configure 64-bit indices
if(PALACE_WITH_64BIT_INT)
  list(APPEND SCOTCH_OPTIONS
    "-DINTSIZE=64"
  )
endif()

# Configure multithreading
if(PALACE_WITH_OPENMP)
  list(APPEND SCOTCH_OPTIONS
    "-DTHREADS=ON"
    "-DMPI_THREAD_MULTIPLE=ON"
  )
else()
  list(APPEND SCOTCH_OPTIONS
    "-DTHREADS=OFF"
    "-DMPI_THREAD_MULTIPLE=OFF"
  )
endif()

string(REPLACE ";" "; " SCOTCH_OPTIONS_PRINT "${SCOTCH_OPTIONS}")
message(STATUS "SCOTCH_OPTIONS: ${SCOTCH_OPTIONS_PRINT}")

# Some build fixes
set(SCOTCH_PATCH_FILES
  "${CMAKE_CURRENT_SOURCE_DIR}/patch/scotch/patch_build.diff"
)

include(ExternalProject)
ExternalProject_Add(scotch
  DEPENDS           ${SCOTCH_DEPENDENCIES}
  GIT_REPOSITORY    ${CMAKE_CURRENT_SOURCE_DIR}/scotch
  GIT_TAG           ${EXTERN_SCOTCH_GIT_TAG}
  SOURCE_DIR        ${CMAKE_CURRENT_BINARY_DIR}/scotch
  BINARY_DIR        ${CMAKE_CURRENT_BINARY_DIR}/scotch-build
  INSTALL_DIR       ${CMAKE_INSTALL_PREFIX}
  PREFIX            ${CMAKE_CURRENT_BINARY_DIR}/scotch-cmake
  UPDATE_COMMAND    ""
  PATCH_COMMAND     git apply "${SCOTCH_PATCH_FILES}"
  CONFIGURE_COMMAND cmake <SOURCE_DIR> "${SCOTCH_OPTIONS}"
  TEST_COMMAND      ""
)
list(APPEND METIS_PARMETIS_DEPENDENCIES gklib)

include(GNUInstallDirs)
if(BUILD_SHARED_LIBS)
  set(_SCOTCH_LIB_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
else()
  set(_SCOTCH_LIB_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()
set(_SCOTCH_LIBRARIES "")
set(_PTSCOTCH_LIBRARIES "")
foreach(LIB scotcherrexit scotcherr scotch scotchmetisv5 esmumps)
  set(_SCOTCH_LIBRARIES ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/lib${LIB}${_SCOTCH_LIB_SUFFIX}$<SEMICOLON>${_SCOTCH_LIBRARIES})
endforeach()
foreach(LIB ptscotcherrexit ptscotcherr ptscotch ptscotchparmetisv3 ptesmumps)
  set(_PTSCOTCH_LIBRARIES ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/lib${LIB}${_SCOTCH_LIB_SUFFIX}$<SEMICOLON>${_PTSCOTCH_LIBRARIES})
endforeach()
set(METIS_LIBRARIES ${_SCOTCH_LIBRARIES} CACHE STRING "List of library files for METIS")
set(PARMETIS_LIBRARIES ${_PTSCOTCH_LIBRARIES} CACHE STRING "List of library files for ParMETIS")
set(METIS_INCLUDE_DIRS ${CMAKE_INSTALL_PREFIX}/include CACHE STRING "Path to METIS include directories")
