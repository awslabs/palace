# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Add source files and subdirectories.
#

target_sources(${TARGET_NAME}
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/amg.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ams.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/arpack.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/chebyshev.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/curlcurl.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/distrelaxation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/divfree.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/feast.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/errorestimator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/gmg.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/hypre.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ksp.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/mumps.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/pc.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/petsc.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/slepc.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/strumpack.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/superlu.cpp
)
