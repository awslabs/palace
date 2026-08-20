// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "jsonschema.hpp"

#include <algorithm>
#include <cctype>
#include <sstream>
#include <string_view>
#include <unordered_set>
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

// Search for a schema property by key name, collecting all matches.
// Helper for FindSchemaByKey - populates results vector. The depth parameter
// guards against pathological self-referential $defs cycles.
void FindAllSchemasByKey(const json &schema, const std::string &key, const json &root_defs,
                         std::vector<json> &results, int depth = 0)
{
  constexpr int max_depth = 32;
  if (depth > max_depth)
  {
    // Hitting the cap means a self-referential $defs cycle in the schema, which
    // is a developer error rather than bad user input. Warn so it surfaces
    // instead of silently truncating the search.
    Mpi::Warning("Schema search for key '{}' exceeded max depth {}; check for a "
                 "self-referential $defs cycle\n",
                 key, max_depth);
    return;
  }
  if (!schema.is_object())
  {
    return;
  }

  // Track $defs from this level.
  const json &defs = schema.contains("$defs") ? schema["$defs"] : root_defs;

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
      // Handle $ref: resolve internal #/$defs/ references.
      if (v.contains("$ref"))
      {
        std::string ref_raw = v["$ref"].get<std::string>();
        if (ref_raw.find("#/$defs/") == 0)
        {
          std::string def_name = ref_raw.substr(8);
          if (!defs.is_null() && defs.contains(def_name))
          {
            FindAllSchemasByKey(defs[def_name], key, defs, results, depth + 1);
          }
        }
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
// Warns and returns null if the key is ambiguous (found in multiple schema scopes).
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
  if (ref.substr(0, 8) == "#/$defs/")
  {
    std::string def_name = ref.substr(8);
    if (defs.contains(def_name))
    {
      return defs[def_name];
    }
  }
  return json();
}

void AppendUnique(json &values, const json &candidate)
{
  if (std::none_of(values.begin(), values.end(),
                   [&candidate](const json &value) { return value == candidate; }))
  {
    values.push_back(candidate);
  }
}

// Collect enum-like values at an exact instance path. Composition branches are searched
// with the same unconsumed path, so a property is never matched merely because it has the
// same name elsewhere in the schema.
void CollectEnumValues(const json &schema, const json &defs,
                       const std::vector<std::string> &tokens, std::size_t token_index,
                       json &values, int depth = 0)
{
  constexpr int max_depth = 32;
  if (!schema.is_object() || depth > max_depth)
  {
    return;
  }

  const json &local_defs = schema.contains("$defs") ? schema["$defs"] : defs;
  if (schema.contains("$ref"))
  {
    json resolved = ResolveRef(schema, local_defs);
    if (!resolved.is_null())
    {
      CollectEnumValues(resolved, local_defs, tokens, token_index, values, depth + 1);
    }
    return;
  }

  if (token_index == tokens.size())
  {
    if (auto enum_it = schema.find("enum"); enum_it != schema.end() && enum_it->is_array())
    {
      for (const auto &value : *enum_it)
      {
        AppendUnique(values, value);
      }
    }
    if (auto const_it = schema.find("const"); const_it != schema.end())
    {
      AppendUnique(values, *const_it);
    }
  }
  else
  {
    const std::string &token = tokens[token_index];
    bool is_index =
        !token.empty() && std::all_of(token.begin(), token.end(),
                                      [](unsigned char c) { return std::isdigit(c); });
    if (is_index)
    {
      if (auto items_it = schema.find("items"); items_it != schema.end())
      {
        CollectEnumValues(*items_it, local_defs, tokens, token_index + 1, values,
                          depth + 1);
      }
    }
    else if (auto properties_it = schema.find("properties");
             properties_it != schema.end() && properties_it->contains(token))
    {
      CollectEnumValues((*properties_it)[token], local_defs, tokens, token_index + 1,
                        values, depth + 1);
    }
  }

  // Composition can occur either at the enum itself or at an object/array containing the
  // property. Preserve the current path position while exploring each possible branch.
  for (const char *keyword : {"oneOf", "anyOf", "allOf"})
  {
    if (auto branches_it = schema.find(keyword);
        branches_it != schema.end() && branches_it->is_array())
    {
      for (const auto &branch : *branches_it)
      {
        CollectEnumValues(branch, local_defs, tokens, token_index, values, depth + 1);
      }
    }
  }
  // A conditional predicate does not itself constrain the instance; either result branch
  // can contribute a spelling. Schema validation remains responsible for applicability.
  for (const char *keyword : {"then", "else"})
  {
    if (auto branch_it = schema.find(keyword); branch_it != schema.end())
    {
      CollectEnumValues(*branch_it, local_defs, tokens, token_index, values, depth + 1);
    }
  }
}

// Find enum values in schema by following an exact JSON pointer path.
json FindEnumInSchema(const json &schema, const std::string &ptr)
{
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
    if (!token.empty())
    {
      tokens.push_back(std::move(token));
    }
    pos = next;
  }

  json values = json::array();
  CollectEnumValues(schema, schema.value("$defs", json::object()), tokens, 0, values);
  return values.empty() ? json() : values;
}

// Find an allowed string which differs from the input only by ASCII case.
const json *FindCaseInsensitiveMatch(const json &allowed, std::string_view input)
{
  for (const auto &candidate : allowed)
  {
    if (!candidate.is_string())
    {
      continue;
    }
    const auto &value = candidate.get_ref<const std::string &>();
    bool matches = input != value && input.size() == value.size() &&
                   std::equal(input.begin(), input.end(), value.begin(),
                              [](unsigned char a, unsigned char b)
                              { return std::tolower(a) == std::tolower(b); });
    if (matches)
    {
      return &candidate;
    }
  }
  return nullptr;
}

}  // namespace

// Custom error handler that formats errors with documentation-style paths.
class SchemaErrorHandler : public error_handler
{
  std::ostringstream errors;
  bool has_error = false;
  const json *schema;
  std::unordered_set<std::string> reported_enum_paths;

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
    const std::string path = ptr.to_string();
    // These message strings are implementation details of json-schema-validator 2.4.0;
    // update them when upgrading the dependency.
    const bool is_enum_failure =
        message.find("instance not found in required enum") != std::string::npos ||
        message.find("instance not const") != std::string::npos ||
        message.find("no subschema has succeeded") != std::string::npos;
    json enum_values;
    if ((schema != nullptr) && is_enum_failure)
    {
      enum_values = FindEnumInSchema(*schema, path);
      if (!enum_values.is_null() && enum_values.is_array() && !enum_values.empty() &&
          !reported_enum_paths.insert(path).second)
      {
        // Logical combinations can report the same enum failure once per branch. The first
        // message already contains the union of values at this exact instance path.
        return;
      }
    }

    errors << "At " << FormatPath(path) << ": " << message;
    // Enhance type mismatch errors with actual type. Use find() so the enhancement also
    // fires for oneOf/anyOf-wrapped messages
    // like "[combination: oneOf / case#0] unexpected instance type".
    if (message.find("unexpected instance type") != std::string::npos)
    {
      errors << " (got " << instance.type_name() << ")";
    }
    else if (!enum_values.is_null() && enum_values.is_array() && !enum_values.empty())
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
      if (instance.is_string())
      {
        const auto &input = instance.get_ref<const std::string &>();
        if (const json *suggestion = FindCaseInsensitiveMatch(enum_values, input))
        {
          errors << ". Did you mean " << suggestion->dump() << "?";
        }
      }
    }
    errors << "\n";
    has_error = true;
  }

  explicit operator bool() const { return has_error; }
  std::string get_errors() const { return errors.str(); }
};

std::string GetSchemaVersion()
{
  const auto &schema_map = schema::GetSchemaMap();
  auto it = schema_map.find(root_schema_file);
  if (it != schema_map.end())
  {
    const json schema = json::parse(it->second);
    if (schema.contains("$id"))
    {
      const std::string &id = schema["$id"];
      constexpr std::string_view prefix = "urn:palace:schema:";
      if (id.substr(0, prefix.size()) == prefix)
      {
        return id.substr(prefix.size());
      }
    }
  }
  return "unknown";
}

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

  // No schema loader is supplied: the single-file schema uses only internal
  // ("#/$defs/...") references. An external $ref would make the validator throw a
  // descriptive error, caught and returned below.
  json_validator validator;
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

  // No schema loader is supplied: see the note in the single-argument overload above.
  json_validator validator;
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
