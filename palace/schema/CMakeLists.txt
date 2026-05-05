# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Schema emission helper under utils/; they depend on reflect-cpp (pulled
# in via FetchContent from the parent palace/CMakeLists.txt).
#
# Two targets:
#   1. palace_schema_utils — static lib built from utils/schema.cpp.
#      Exposes the utils/ headers via its public include directory so the
#      emit binary (and any future consumer) can include them as
#      "schema/utils/<file>.hpp" from the palace source root.
#   2. palace_emit_schema — console helper that prints the generated schema
#      for `palace::schema::PalaceConfig` to stdout. Invoked at build time
#      by cmake/EmbedSchema.cmake.
#
# Neither target is linked into libpalace in Phase 1; reflect-cpp stays out
# of the runtime binary.

add_library(palace_schema_utils STATIC utils/schema.cpp)
# PUBLIC includes the palace/ source root so consumers can write
# `#include "schema/utils/generator.hpp"`.
target_include_directories(palace_schema_utils PUBLIC ${CMAKE_SOURCE_DIR})
target_link_libraries(palace_schema_utils PUBLIC reflectcpp)
target_compile_features(palace_schema_utils PUBLIC cxx_std_20)

add_executable(palace_emit_schema emit_schema.cpp)
target_link_libraries(palace_emit_schema PRIVATE palace_schema_utils)
target_compile_features(palace_emit_schema PRIVATE cxx_std_20)
target_compile_definitions(palace_emit_schema PRIVATE
    PALACE_VERSION="${PROJECT_VERSION}"
)
