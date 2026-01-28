// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "jsonschema.hpp"

#include <fstream>
#include <nlohmann/json-schema.hpp>

namespace palace
{

using json = nlohmann::json;
using json_validator = nlohmann::json_schema::json_validator;

namespace
{

// Loader for resolving $ref in schemas.
class SchemaLoader
{
  std::string base_dir;

public:
  SchemaLoader(const std::string &dir) : base_dir(dir) {}

  json operator()(const nlohmann::json_uri &uri, json &schema)
  {
    std::string path = uri.path();
    if (path.substr(0, 2) == "./")
    {
      path = path.substr(2);
    }
    std::string full_path = base_dir + "/" + path;
    std::ifstream f(full_path);
    if (!f.is_open())
    {
      throw std::runtime_error("Failed to open schema: " + full_path);
    }
    schema = json::parse(f);
    return schema;
  }
};

// Depth-first search for a schema property by key name.
// Returns the schema for that key (with $defs preserved), or null if not found.
json FindSchemaByKey(const json &schema, const std::string &key,
                     const std::string &schema_dir, const json &root_defs = json())
{
  if (!schema.is_object())
  {
    return json();
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
    return result;
  }

  // Recurse into properties.
  if (props_it != schema.end())
  {
    for (const auto &[k, v] : props_it->items())
    {
      // Handle $ref.
      if (v.contains("$ref"))
      {
        std::string ref = v["$ref"].get<std::string>();
        if (ref.substr(0, 2) == "./")
        {
          ref = ref.substr(2);
        }
        std::string path = schema_dir + "/" + ref;
        std::ifstream f(path);
        if (f.is_open())
        {
          json resolved = json::parse(f);
          auto result = FindSchemaByKey(resolved, key, schema_dir, defs);
          if (!result.is_null())
          {
            return result;
          }
        }
      }
      else
      {
        auto result = FindSchemaByKey(v, key, schema_dir, defs);
        if (!result.is_null())
        {
          return result;
        }
      }
    }
  }

  // Check items for arrays.
  auto items_it = schema.find("items");
  if (items_it != schema.end())
  {
    auto result = FindSchemaByKey(*items_it, key, schema_dir, defs);
    if (!result.is_null())
    {
      return result;
    }
  }

  return json();
}

}  // namespace

std::string ValidateConfig(const json &config, const std::string &schema_dir)
{
  std::string schema_path = schema_dir + "/config-schema.json";
  std::ifstream f(schema_path);
  if (!f.is_open())
  {
    return "Failed to open schema file: " + schema_path;
  }

  json schema;
  try
  {
    schema = json::parse(f);
  }
  catch (json::parse_error &e)
  {
    return std::string("Failed to parse schema: ") + e.what();
  }

  json_validator validator{SchemaLoader(schema_dir)};
  try
  {
    validator.set_root_schema(schema);
    validator.validate(config);
    return "";
  }
  catch (std::exception &e)
  {
    return e.what();
  }
}

std::string ValidateConfig(const json &config, const std::string &schema_dir,
                           const std::string &schema_key)
{
  std::string schema_path = schema_dir + "/config-schema.json";
  std::ifstream f(schema_path);
  if (!f.is_open())
  {
    return "Failed to open schema file: " + schema_path;
  }

  json root_schema;
  try
  {
    root_schema = json::parse(f);
  }
  catch (json::parse_error &e)
  {
    return std::string("Failed to parse schema: ") + e.what();
  }

  // Find the sub-schema by key.
  json sub_schema = FindSchemaByKey(root_schema, schema_key, schema_dir);
  if (sub_schema.is_null())
  {
    return "Schema key not found: " + schema_key;
  }

  json_validator validator{SchemaLoader(schema_dir)};
  try
  {
    validator.set_root_schema(sub_schema);
    validator.validate(config);
    return "";
  }
  catch (std::exception &e)
  {
    return e.what();
  }
}

}  // namespace palace
