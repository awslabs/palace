# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Add source files and subdirectories.
#

target_sources(${TARGET_NAME}
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/bilinearform.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/coefficient.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/integrator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/interpolator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/lumpedelement.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/integ/curlcurl.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/integ/curlcurlmass.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/integ/diffusion.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/integ/diffusionmass.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/integ/divdiv.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/integ/divdivmass.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/integ/mass.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/integ/mixedveccurl.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/integ/mixedvecgrad.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/integ/vecfemass.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/libceed/basis.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/libceed/operator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/libceed/restriction.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/libceed/utils.cpp
)
