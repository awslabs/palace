// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_SCHEMA_VERSION_HPP
#define PALACE_SCHEMA_VERSION_HPP

#include <string_view>

namespace palace::schema {

// SchemaVer version of the Palace configuration schema.
// See https://docs.snowplow.io/docs/api-reference/iglu/common-architecture/schemaver/
//
// Format: MODEL-REVISION-ADDITION (three non-negative integers, hyphen-
// separated, starts at 1-0-0). This tracks *schema compatibility*, not
// Palace's release version — unlike `PROJECT_VERSION`, the schema version
// only moves when the shape of the config changes.
//
// Bump guidance:
//
//   MODEL    — incompatible change. A consumer holding the old version
//              can no longer validate a document produced under the new
//              version, or vice versa. Examples: renamed field, removed
//              enum value, tightened validator, new required field
//              without a default.
//
//   REVISION — compatible with the previous MODEL but introduces a
//              change a consumer should take note of. Examples: new
//              enumerator in an enum that's also emitted on read, new
//              optional field whose presence affects behavior, loosened
//              validator that accepts more inputs.
//
//   ADDITION — purely additive, zero-impact on existing consumers.
//              Examples: new description text, new
//              `x-palace-advanced` flag, new optional field with an
//              unchanged default, new oneOf arm whose absence is still
//              valid input.
//
// Palace starts this schema at 1-0-0. Subsequent PRs that touch
// `palace/schema/types/` or the post-emit passes must update this
// constant in the same commit.
inline constexpr std::string_view schema_version = "1-0-0";

}  // namespace palace::schema

#endif  // PALACE_SCHEMA_VERSION_HPP
