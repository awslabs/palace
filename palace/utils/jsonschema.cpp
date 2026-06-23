// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "jsonschema.hpp"

#include <algorithm>
#include <sstream>
#include <nlohmann/json-schema.hpp>
#include "communication.hpp"
#include "embedded_schema.hpp"

namespace palace
{

using json = nlohmann::json;
using json_validator = nlohmann::json_schema::json_validator;
using error_handler = nlohmann::json_schema::error_handler;

// Root schema entry point. The schema is generated at build time by
// palace_emit_schema (see cmake/EmbedSchema.cmake); it is a single draft
// 2020-12 document with all sub-section shapes under `$defs`.
constexpr const char *root_schema_file = "config-schema.json";

namespace
{

// Loader used by nlohmann_json_schema_validator's root_schema setup when it
// encounters a `$ref`. Since the generated schema is self-contained (all refs
// are in-document `#/$defs/...`), the loader only ever sees the root URI and
// returns the single embedded schema.
class EmbeddedSchemaLoader
{
public:
  json operator()(const nlohmann::json_uri &, json &schema)
  {
    const auto &schema_map = schema::GetSchemaMap();
    auto it = schema_map.find(root_schema_file);
    if (it == schema_map.end())
    {
      throw std::runtime_error("Root schema not found in embedded schemas");
    }
    schema = json::parse(it->second);
    return schema;
  }
};

// Search for a schema property by key name within the single embedded schema
// document, collecting all matches. Recurses through `properties`, resolves
// in-document `#/$defs/...` refs, and follows `items` for arrays. External-file
// $ref support was removed when the schema collapsed to a single generated
// file. The depth parameter guards against pathological self-referential $defs
// cycles.
void FindAllSchemasByKey(const json &schema, const std::string &key, const json &root_defs,
                         std::vector<json> &results, int depth = 0)
{
  constexpr int kMaxDepth = 32;
  if (depth > kMaxDepth)
  {
    // Hitting the cap means a self-referential $defs cycle in the schema, which
    // is a developer error rather than bad user input. Warn so it surfaces
    // instead of silently truncating the search.
    Mpi::Warning("Schema search for key '{}' exceeded max depth {}; check for a "
                 "self-referential $defs cycle\n",
                 key, kMaxDepth);
    return;
  }
  if (!schema.is_object())
  {
    return;
  }

  // Track $defs from this level (carried down so nested $ref can resolve).
  const json &defs = schema.contains("$defs") ? schema["$defs"] : root_defs;

  // The generated schema root is a $ref into $defs, and nested property bodies
  // can also be pure references. Follow local refs before looking for child
  // properties so named-fragment validation works with the single-file schema.
  if (auto ref_it = schema.find("$ref"); ref_it != schema.end() && ref_it->is_string())
  {
    const auto ref_raw = ref_it->get<std::string>();
    if (ref_raw.rfind("#/$defs/", 0) == 0)
    {
      const auto def_name = ref_raw.substr(8);  // strlen("#/$defs/")
      if (defs.is_object() && defs.contains(def_name))
      {
        FindAllSchemasByKey(defs[def_name], key, defs, results, depth + 1);
      }
    }
  }

  // Check properties at this level.
  auto props_it = schema.find("properties");
  if (props_it != schema.end() && props_it->contains(key))
  {
    json result = (*props_it)[key];
    // For array items, use the items schema.
    if (result.contains("items"))
    {
      result = result["items"];
    }
    // Attach $defs so in-document $ref can resolve.
    if (!defs.is_null() && !defs.empty())
    {
      result["$defs"] = defs;
    }
    results.push_back(result);
  }

  // Recurse into properties, following in-document $ref along the way.
  if (props_it != schema.end())
  {
    for (const auto &[_k, v] : props_it->items())
    {
      if (v.contains("$ref"))
      {
        // Internal $defs reference: resolve and recurse so nested-via-$ref
        // schemas (e.g. BoundaryPostprocessing.FarField) are still findable.
        const auto ref_raw = v["$ref"].get<std::string>();
        if (ref_raw.rfind("#/$defs/", 0) == 0)
        {
          const auto def_name = ref_raw.substr(8);  // strlen("#/$defs/")
          if (defs.is_object() && defs.contains(def_name))
          {
            FindAllSchemasByKey(defs[def_name], key, defs, results, depth + 1);
          }
        }
        // Non-local $ref shapes are unreachable in the generated schema.
      }
      else
      {
        FindAllSchemasByKey(v, key, defs, results, depth + 1);
      }
    }
  }

  // Check items for arrays.
  if (auto items_it = schema.find("items"); items_it != schema.end())
  {
    FindAllSchemasByKey(*items_it, key, defs, results, depth + 1);
  }
}

// Search for a schema property by key name, checking each level before recursing.
// Returns the schema for that key (with $defs preserved), or null if not found.
// Warns and returns null if the key is ambiguous (found in multiple schema files).
json FindSchemaByKey(const json &schema, const std::string &key,
                     const json &root_defs = json())
{
  std::vector<json> results;
  FindAllSchemasByKey(schema, key, root_defs, results);

  if (results.empty())
  {
    return json();
  }
  if (results.size() > 1)
  {
    Mpi::Warning("Ambiguous schema key '{}' found {} times; use a more specific path\n",
                 key, results.size());
    return json();
  }
  return results[0];
}

// Resolve a $ref in a schema node. Returns the resolved schema, or null on failure.
// The defs parameter provides $defs context for internal references.
json ResolveRef(const json &node, const json &defs)
{
  if (!node.contains("$ref"))
  {
    return node;
  }
  std::string ref = node["$ref"].get<std::string>();
  if (ref.empty())
  {
    return json();  // Invalid empty $ref.
  }
  // The generated schema only emits in-document `#/$defs/...` refs; external-
  // file refs disappeared when scripts/schema/config/*.json was replaced by a
  // single build-time-generated config-schema.json.
  if (ref.rfind("#/$defs/", 0) == 0)
  {
    std::string def_name = ref.substr(8);
    if (defs.contains(def_name))
    {
      return defs[def_name];
    }
  }
  return json();
}

// Find enum values in schema by following a JSON pointer path.
json FindEnumInSchema(const json &schema, const std::string &ptr)
{
  if (ptr.empty() || ptr == "/")
  {
    return json();
  }

  // Parse path into tokens (skip numeric indices for arrays).
  std::vector<std::string> tokens;
  std::size_t pos = 0;
  while (pos < ptr.size())
  {
    if (ptr[pos] == '/')
    {
      pos++;
      continue;
    }
    std::size_t next = ptr.find('/', pos);
    std::string token = ptr.substr(pos, next - pos);
    // Skip array indices.
    if (token.empty() || !std::all_of(token.begin(), token.end(), ::isdigit))
    {
      tokens.push_back(token);
    }
    pos = next;
  }

  // Get $defs from root schema for resolving internal refs.
  json defs = schema.value("$defs", json::object());

  // Navigate schema following the path.
  json current = schema;
  for (const auto &token : tokens)
  {
    // Resolve $ref chain.
    while (current.contains("$ref"))
    {
      json resolved = ResolveRef(current, defs);
      if (resolved.is_null())
      {
        return json();
      }
      // Pick up $defs from resolved schema.
      if (resolved.contains("$defs"))
      {
        defs = resolved["$defs"];
      }
      current = resolved;
    }

    // Look in properties.
    if (current.contains("properties") && current["properties"].contains(token))
    {
      current = current["properties"][token];
      // Follow into array items.
      if (current.contains("items"))
      {
        current = current["items"];
      }
    }
    else
    {
      return json();
    }
  }

  // Final $ref resolution.
  while (current.contains("$ref"))
  {
    json resolved = ResolveRef(current, defs);
    if (resolved.is_null())
    {
      break;
    }
    current = resolved;
  }

  // Extract enum values (handle anyOf).
  if (current.contains("enum"))
  {
    return current["enum"];
  }
  if (current.contains("anyOf"))
  {
    for (const auto &opt : current["anyOf"])
    {
      if (opt.contains("enum"))
      {
        return opt["enum"];
      }
    }
  }
  return json();
}

std::string ValidateBoundaryMutualExclusion(const json &boundaries)
{
  if (!boundaries.is_object())
  {
    return "";
  }
  if (boundaries.contains("PEC") && boundaries.contains("Ground"))
  {
    return "At [\"Boundaries\"]: properties 'PEC' and 'Ground' are mutually exclusive\n";
  }
  if (boundaries.contains("PMC") && boundaries.contains("ZeroCharge"))
  {
    return "At [\"Boundaries\"]: properties 'PMC' and 'ZeroCharge' are mutually "
           "exclusive\n";
  }
  return "";
}

}  // namespace

// Custom error handler that formats errors with documentation-style paths.
class SchemaErrorHandler : public error_handler
{
  std::ostringstream errors;
  bool has_error = false;
  const json *schema;

  // Convert JSON pointer "/Boundaries/LumpedPort/0" to ["Boundaries"]["LumpedPort"][0]
  static std::string FormatPath(const std::string &ptr)
  {
    if (ptr.empty() || ptr == "/")
    {
      return "config";
    }
    std::ostringstream result;
    std::size_t pos = 0;
    while (pos < ptr.size())
    {
      if (ptr[pos] == '/')
      {
        pos++;
        continue;
      }
      std::size_t next = ptr.find('/', pos);
      std::string token = ptr.substr(pos, next - pos);
      // Check if token is a number (array index).
      bool is_index = !token.empty() && std::all_of(token.begin(), token.end(), ::isdigit);
      if (is_index)
      {
        result << "[" << token << "]";
      }
      else
      {
        result << "[\"" << token << "\"]";
      }
      pos = next;
    }
    return result.str();
  }

public:
  explicit SchemaErrorHandler(const json *s = nullptr) : schema(s) {}

  void error(const json::json_pointer &ptr, const json &instance,
             const std::string &message) override
  {
    errors << "At " << FormatPath(ptr.to_string()) << ": " << message;
    // Enhance type mismatch errors with actual type. These message strings are
    // implementation details of json-schema-validator 2.4.0; update if upgrading.
    // Use find() so the enhancement also fires for oneOf/anyOf-wrapped messages
    // like "[combination: oneOf / case#0] unexpected instance type".
    if (message.find("unexpected instance type") != std::string::npos)
    {
      errors << " (got " << instance.type_name() << ")";
    }
    // Enhance enum errors with valid values.
    else if (schema &&
             message.find("instance not found in required enum") != std::string::npos)
    {
      json enum_values = FindEnumInSchema(*schema, ptr.to_string());
      if (!enum_values.is_null() && enum_values.is_array() && !enum_values.empty())
      {
        errors << "; valid values: ";
        for (std::size_t i = 0; i < enum_values.size(); i++)
        {
          if (i > 0)
          {
            errors << ", ";
          }
          errors << enum_values[i].dump();
        }
      }
    }
    errors << "\n";
    has_error = true;
  }

  explicit operator bool() const { return has_error; }
  std::string get_errors() const { return errors.str(); }
};

std::string ValidateConfig(const nlohmann::json &config)
{
  if (auto boundaries = config.find("Boundaries"); boundaries != config.end())
  {
    auto err = ValidateBoundaryMutualExclusion(*boundaries);
    if (!err.empty())
    {
      return err;
    }
  }

  const auto &schema_map = schema::GetSchemaMap();
  auto it = schema_map.find(root_schema_file);
  if (it == schema_map.end())
  {
    return "Root schema not found in embedded schemas";
  }

  json schema;
  try
  {
    schema = json::parse(it->second);
  }
  catch (json::parse_error &e)
  {
    return std::string("Failed to parse schema: ") + e.what();
  }

  json_validator validator{EmbeddedSchemaLoader()};
  try
  {
    validator.set_root_schema(schema);
    SchemaErrorHandler handler(&schema);
    validator.validate(config, handler);
    if (handler)
    {
      return handler.get_errors();
    }
    return "";
  }
  catch (std::exception &e)
  {
    return e.what();
  }
}

std::string ValidateConfig(const nlohmann::json &config, const std::string &schema_key)
{
  if (schema_key == "Boundaries")
  {
    auto err = ValidateBoundaryMutualExclusion(config);
    if (!err.empty())
    {
      return err;
    }
  }

  const auto &schema_map = schema::GetSchemaMap();
  auto it = schema_map.find(root_schema_file);
  if (it == schema_map.end())
  {
    return "Root schema not found in embedded schemas";
  }

  json root_schema;
  try
  {
    root_schema = json::parse(it->second);
  }
  catch (json::parse_error &e)
  {
    return std::string("Failed to parse schema: ") + e.what();
  }

  // Find the sub-schema by key.
  json sub_schema = FindSchemaByKey(root_schema, schema_key);
  if (sub_schema.is_null())
  {
    return "Schema key not found: " + schema_key;
  }

  json_validator validator{EmbeddedSchemaLoader()};
  try
  {
    validator.set_root_schema(sub_schema);
    SchemaErrorHandler handler(&sub_schema);
    validator.validate(config, handler);
    if (handler)
    {
      return handler.get_errors();
    }
    return "";
  }
  catch (std::exception &e)
  {
    return e.what();
  }
}

}  // namespace palace
