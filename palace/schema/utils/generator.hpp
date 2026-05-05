// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_SCHEMA_UTILS_GENERATOR_HPP
#define PALACE_SCHEMA_UTILS_GENERATOR_HPP

#include <string>

#include "annotations.hpp"

namespace palace::schema::utils {

// Emit a JSON Schema (draft 2020-12) string for `T`.
//
// With the default `SchemaOptions{}` (`emit_defaults = true`), the schema is
// produced by reflect-cpp and then augmented with `"default": ...` entries
// harvested from `T{}` for every field that has a C++ default; fields with
// defaults are also pruned from the parent's `required` list.
//
// With `SchemaOptions{.emit_defaults = false}`, neither pass runs — the
// output carries reflect-cpp's native `required` list (every non-optional
// field is required) and no `default` keys.
template <class T>
std::string schema(SchemaOptions opts = {});

}  // namespace palace::schema::utils

#include "schema_impl.hpp"

#endif  // PALACE_SCHEMA_UTILS_SCHEMA_HPP
