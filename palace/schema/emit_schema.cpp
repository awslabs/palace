// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Standalone helper that prints the Palace JSON schema (draft 2020-12) to
// stdout. Invoked from CMake at build time; the output is captured into
// build/generated/schema/config-schema.json and then embedded into libpalace
// via the existing embed_schema.cmake helper.

#include <iostream>
#include <string>
#include <variant>

#include <rfl/Generic.hpp>
#include <rfl/json.hpp>

#include "schema/types/config.hpp"
#include "schema/utils/generator.hpp"
#include "schema/version.hpp"

namespace
{

// Root-level cross-field rule: when `Problem.Type` picks a specific
// simulation mode, the matching `Solver.<Mode>` sub-block is required. This
// is schema-level cross-field logic that cannot be expressed from the type
// alone; the fragment is copied verbatim from the pre-rewrite
// `scripts/schema/config-schema.json` so runtime validators (including
// Julia's `JSONSchema.jl`) enforce the same contract as Palace's drivers.
constexpr const char *kRootConditional = R"([
    {
      "if": {
        "properties": {
          "Problem": {
            "properties": { "Type": { "const": "Driven" } },
            "required": ["Type"]
          }
        }
      },
      "then": {
        "properties": {
          "Solver": { "required": ["Driven"] }
        }
      }
    },
    {
      "if": {
        "properties": {
          "Problem": {
            "properties": { "Type": { "const": "Eigenmode" } },
            "required": ["Type"]
          }
        }
      },
      "then": {
        "properties": {
          "Solver": { "required": ["Eigenmode"] }
        }
      }
    },
    {
      "if": {
        "properties": {
          "Problem": {
            "properties": { "Type": { "const": "Transient" } },
            "required": ["Type"]
          }
        }
      },
      "then": {
        "properties": {
          "Solver": { "required": ["Transient"] }
        }
      }
    },
    {
      "if": {
        "properties": {
          "Problem": {
            "properties": { "Type": { "const": "BoundaryMode" } },
            "required": ["Type"]
          }
        }
      },
      "then": {
        "properties": {
          "Solver": { "required": ["BoundaryMode"] }
        }
      }
    }
  ])";

// Parse the conditional-block JSON array and splice it into the root
// schema's `allOf` key. Any pre-existing `allOf` on the root is replaced —
// reflect-cpp does not emit one, and this is the single intended producer.
std::string inject_root_allof(std::string schema_json)
{
  auto root_r = rfl::json::read<rfl::Generic>(schema_json);
  if (!root_r)
  {
    return schema_json;
  }
  auto &root_var = root_r->value();
  if (!std::holds_alternative<rfl::Generic::Object>(root_var))
  {
    return schema_json;
  }
  auto &root_obj = std::get<rfl::Generic::Object>(root_var);

  auto cond_r = rfl::json::read<rfl::Generic>(kRootConditional);
  if (!cond_r)
  {
    return schema_json;
  }
  root_obj["allOf"] = std::move(*cond_r);
  return rfl::json::write(*root_r, rfl::json::pretty);
}

}  // namespace

int main()
{
  auto s = palace::schema::utils::schema<palace::schema::PalaceConfiguration>(
      {.emit_defaults = true,
       .version = std::string(palace::schema::schema_version),
       .defs_prefix = "palace__schema__"});
  s = inject_root_allof(std::move(s));
  std::cout << s;
  return 0;
}
