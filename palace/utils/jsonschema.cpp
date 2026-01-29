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

// Root schema entry point. Schemas use JSON Schema draft-07.
constexpr const char *root_schema_file = "config-schema.json";

namespace
{

// Strip leading path prefixes (/, ./) to normalize schema references.
std::string StripPathPrefix(std::string path)
{
  while (!path.empty() && (path[0] == '/' || path[0] == '.'))
  {
    path = path.substr(1);
  }
  return path;
}

// Loader for resolving $ref in schemas using embedded schema strings.
class EmbeddedSchemaLoader
{
public:
  json operator()(const nlohmann::json_uri &uri, json &schema)
  {
    std::string path = StripPathPrefix(uri.path());
    const auto &schema_map = schema::GetSchemaMap();
    auto it = schema_map.find(path);
    if (it == schema_map.end())
    {
      throw std::runtime_error("Schema not found: " + path);
    }
    schema = json::parse(it->second);
    return schema;
  }
};

// Search for a schema property by key name, collecting all matches.
// Helper for FindSchemaByKey - populates results vector.
void FindAllSchemasByKey(const json &schema, const std::string &key, const json &root_defs,
                         std::vector<json> &results)
{
  if (!schema.is_object())
  {
    return;
  }

  // Track $defs from this level.
  json defs = root_defs;
  if (schema.contains("$defs"))
  {
    defs = schema["$defs"];
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
    // Attach $defs so $ref can resolve.
    if (!defs.is_null() && !defs.empty())
    {
      result["$defs"] = defs;
    }
    results.push_back(result);
  }

  // Recurse into properties.
  if (props_it != schema.end())
  {
    for (const auto &[k, v] : props_it->items())
    {
      // Handle $ref.
      if (v.contains("$ref"))
      {
        std::string ref_raw = v["$ref"].get<std::string>();
        // Skip internal $defs references (handled elsewhere).
        if (ref_raw.find("#/$defs/") == 0)
        {
          continue;
        }
        std::string ref = StripPathPrefix(ref_raw);
        const auto &schema_map = schema::GetSchemaMap();
        if (auto it = schema_map.find(ref); it != schema_map.end())
        {
          FindAllSchemasByKey(json::parse(it->second), key, defs, results);
        }
        else
        {
          Mpi::Warning("Could not resolve schema $ref: {}\n", ref);
        }
      }
      else
      {
        FindAllSchemasByKey(v, key, defs, results);
      }
    }
  }

  // Check items for arrays.
  if (auto items_it = schema.find("items"); items_it != schema.end())
  {
    FindAllSchemasByKey(*items_it, key, defs, results);
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

// Find enum values in schema by following a JSON pointer path.
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
  if (ref[0] == '.')
  {
    const auto &schema_map = schema::GetSchemaMap();
    if (auto it = schema_map.find(StripPathPrefix(ref)); it != schema_map.end())
    {
      return json::parse(it->second);
    }
  }
  else if (ref.substr(0, 8) == "#/$defs/")
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
    // Enhance type mismatch errors with actual type.
    if (message == "unexpected instance type")
    {
      errors << " (got " << instance.type_name() << ")";
    }
    // Enhance enum errors with valid values.
    else if (schema && message == "instance not found in required enum")
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

  operator bool() const { return has_error; }
  std::string get_errors() const { return errors.str(); }
};

std::string ValidateConfig(const nlohmann::json &config)
{
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
