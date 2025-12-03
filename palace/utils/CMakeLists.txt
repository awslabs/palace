# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Add source files and subdirectories.
#

target_sources(${LIB_TARGET_NAME}
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/configfile.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/device.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/dorfler.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/geodata.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/geodata_impl.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/iodata.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/meshio.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/omp.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/tablecsv.cpp
)
