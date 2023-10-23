# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Add source files and subdirectories.
#

target_sources(${LIB_TARGET_NAME}
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/basesolver.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/drivensolver.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/eigensolver.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/electrostaticsolver.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/magnetostaticsolver.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/transientsolver.cpp
)
