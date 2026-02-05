// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_JSONSCHEMA_HPP
#define PALACE_UTILS_JSONSCHEMA_HPP

#include <string>
#include <nlohmann/json_fwd.hpp>

namespace palace
{

// Validate a full JSON config against the embedded root schema.
// Returns empty string on success, error message on failure.
std::string ValidateConfig(const nlohmann::json &config);

// Validate a JSON fragment against a named sub-schema (e.g., "LumpedPort", "Materials").
// Searches depth-first through schema properties to find the matching key.
// Returns empty string on success, error message on failure.
std::string ValidateConfig(const nlohmann::json &config, const std::string &schema_key);

}  // namespace palace

#endif  // PALACE_UTILS_JSONSCHEMA_HPP
