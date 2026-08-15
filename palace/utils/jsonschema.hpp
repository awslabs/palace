// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_JSONSCHEMA_HPP
#define PALACE_UTILS_JSONSCHEMA_HPP

#include <string>
#include <nlohmann/json_fwd.hpp>

namespace palace
{

// Return the schema version extracted from the "$id" URN (e.g. "1-0-0").
std::string GetSchemaVersion();

// Rewrite enum-valued strings in the config to the canonical spelling declared in the
// schema, matching case-insensitively (e.g. "superlu" -> "SuperLU"). Values that already
// match exactly, or that do not match any allowed spelling, are left untouched so that
// schema validation still reports genuinely invalid values. Run this before ValidateConfig,
// which matches enum values case-sensitively.
void CanonicalizeEnumCase(nlohmann::json &config);

// Validate a full JSON config against the embedded root schema.
// Returns empty string on success, error message on failure.
std::string ValidateConfig(const nlohmann::json &config);

// Validate a JSON fragment against a named sub-schema (e.g., "LumpedPort", "Materials").
// Searches depth-first through schema properties to find the matching key.
// Returns empty string on success, error message on failure.
std::string ValidateConfig(const nlohmann::json &config, const std::string &schema_key);

}  // namespace palace

#endif  // PALACE_UTILS_JSONSCHEMA_HPP
