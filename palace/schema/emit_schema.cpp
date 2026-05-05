// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Standalone helper that prints the Palace JSON schema (draft 2020-12) to
// stdout. Invoked from CMake at build time; the output is captured into
// build/generated/schema/config-schema.json and then embedded into libpalace
// via the existing embed_schema.cmake helper.

#include <iostream>

#include "schema/utils/generator.hpp"
#include "schema/types/config.hpp"

#ifndef PALACE_VERSION
#  define PALACE_VERSION "unknown"
#endif

int main()
{
    std::cout << palace::schema::utils::schema<palace::schema::PalaceConfig>(
        {.emit_defaults = true, .version = PALACE_VERSION});
    return 0;
}
