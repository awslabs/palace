# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Add source files and subdirectories.
#

target_sources(${LIB_TARGET_NAME}
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/curlcurloperator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/domainpostoperator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/farfieldboundaryoperator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/laplaceoperator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/lumpedportoperator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/materialoperator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/postoperator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/portexcitationhelper.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/romoperator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/spaceoperator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/surfaceconductivityoperator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/surfacecurrentoperator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/surfaceimpedanceoperator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/surfacepostoperator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/timeoperator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/waveportoperator.cpp
)
