# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Add source files and subdirectories.
#

target_sources(${LIB_TARGET_NAME}
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/amg.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ams.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/blockprecond.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/arpack.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/chebyshev.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/conformingprolongation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cudss.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/deflation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/densematrix.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/distrelaxation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/divfree.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/entityjacobi.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/errorestimator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/floquetcorrection.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/fusedvector.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/gmg.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/gmresprobe.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/hcurl.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/hypre.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/jacobi.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ksp.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/iterative.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/lowrank.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/mumps.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/nleps.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/operator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/rap.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/slepc.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/sliverpatch.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solverstats.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solver.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/spectralnorm.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/strumpack.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/superlu.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/vector.cpp
)

# The CSR row sums round every product, to match the arithmetic of hypre's host products.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|IntelLLVM")
  set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/fusedvector.cpp
    TARGET_DIRECTORY ${LIB_TARGET_NAME}
    PROPERTIES COMPILE_OPTIONS "-ffp-contract=off"
  )
endif()

# Handle device source code
set(TARGET_SOURCES_DEVICE
  ${TARGET_SOURCES_DEVICE}
  ${CMAKE_CURRENT_SOURCE_DIR}/blockprecond.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/chebyshev.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cudss.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/densematrix.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/jacobi.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/lowrank.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/operator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/vector.cpp
  PARENT_SCOPE
)
