// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

#include <rfl.hpp>
#include <rfl/Generic.hpp>
#include <rfl/json.hpp>

// Needed for the `EnumDescEntry` / `FieldFlag` struct definitions, which
// are shared with the compile-time type-graph walker in the template
// header.
#include "schema_impl.hpp"

namespace palace::schema::utils::detail
{

namespace
{

// Walk the `properties` map of one `$defs` entry in lock-step with the
// parsed JSON of a default-constructed instance of that entry's type. For
// each property:
//   * If the defaults JSON has a matching key, write it as `"default"` at
//     the property body.
//   * If the defaults JSON does not have the key, treat the field as an
//     empty std::optional and inject `"default": null` — reflect-cpp's JSON
//     writer drops null optionals from the serialised output, so absence is
//     the only signal we have for that case, and `T{}` must be constructible
//     (so non-optional missing keys cannot happen at this point).
//
// We deliberately do NOT recurse into nested `$ref` / `items.$ref` / inline
// `properties` here. Each reachable aggregate type has its own entry in the
// outer defaults map, and that entry gets visited on its own pass — so
// recursion would just overwrite the same slots we're about to fill from the
// outer loop. This is the core architectural difference vs. the old
// $ref-rooted walk: visiting every `$defs` entry directly uniformly covers
// `items.$ref` (arrays of structs) and TaggedUnion arms, neither of which
// the old walk could reach.
void inject_properties(rfl::Generic::Object &schema_entry,
                       const rfl::Generic::Object &defaults_map)
{
  if (schema_entry.count("properties") == 0)
    return;
  auto &props_var = schema_entry.at("properties").value();
  if (!std::holds_alternative<rfl::Generic::Object>(props_var))
    return;
  auto &props = std::get<rfl::Generic::Object>(props_var);

  for (auto &[name, prop_generic] : props)
  {
    if (!std::holds_alternative<rfl::Generic::Object>(prop_generic.value()))
    {
      continue;
    }
    auto &prop_obj = std::get<rfl::Generic::Object>(prop_generic.value());

    // Properties carrying a `$ref` inherit their default from the
    // referenced `$defs` entry; stamping a `"default"` here would
    // duplicate that information and bloat the schema. Skip them.
    if (prop_obj.count("$ref") != 0)
    {
      continue;
    }

    if (defaults_map.count(name) == 0)
    {
      prop_obj["default"] = rfl::Generic(std::nullopt);
      continue;
    }
    // operator[] returns an existing slot or creates a new one, so we
    // overwrite rather than duplicate if the pass runs twice on the same
    // schema.
    prop_obj["default"] = defaults_map.at(name);
  }
}

// Walk the entire schema tree and add `"additionalProperties": false` to every
// object that declares `properties`. reflect-cpp does not emit this keyword,
// so external validators (e.g. Python jsonschema) would otherwise accept
// unknown keys even though our loader rejects them via rfl::NoExtraFields.
// Keeping the schema aligned with the runtime behavior avoids that drift.
void inject_additional_properties_false_walk(rfl::Generic &node)
{
  if (std::holds_alternative<rfl::Generic::Object>(node.value()))
  {
    auto &obj = std::get<rfl::Generic::Object>(node.value());
    if (obj.count("properties") != 0 && obj.count("additionalProperties") == 0)
    {
      obj["additionalProperties"] = rfl::Generic(false);
    }
    for (auto &[_, v] : obj)
    {
      inject_additional_properties_false_walk(v);
    }
  }
  else if (std::holds_alternative<rfl::Generic::Array>(node.value()))
  {
    auto &arr = std::get<rfl::Generic::Array>(node.value());
    for (auto &v : arr)
    {
      inject_additional_properties_false_walk(v);
    }
  }
}

// Collapse a multi-constraint `allOf` block produced by rfl's schema
// emitter into a single flat body. Each rfl::Validator contributes its
// own `{type, <keyword>}` object; a Closed<int, 1, 2> (two validators)
// comes out as `allOf: [{type:"integer", minimum:1}, {type:"integer",
// maximum:2}]`. Callers prefer `{type:"integer", minimum:1, maximum:2}`
// — same meaning, simpler schema.
//
// Preconditions for merging: every arm of `allOf` must be an object,
// contain only keys from the numeric-constraint whitelist, and agree on
// `type`. When any arm violates those, the `allOf` is left as-is (the
// emitter uses `allOf` for legitimately compound schemas too, and we
// only want to flatten the Validator-induced shape).
bool is_numeric_constraint_key(std::string_view k)
{
  return k == "type" || k == "minimum" || k == "maximum" || k == "exclusiveMinimum" ||
         k == "exclusiveMaximum" || k == "multipleOf";
}

bool arm_is_mergeable(const rfl::Generic &arm,
                      const std::optional<std::string> &expected_type)
{
  if (!std::holds_alternative<rfl::Generic::Object>(arm.value()))
    return false;
  const auto &obj = std::get<rfl::Generic::Object>(arm.value());
  for (const auto &entry : obj)
  {
    if (!is_numeric_constraint_key(entry.first))
      return false;
  }
  if (expected_type && obj.count("type") != 0)
  {
    const auto &t = obj.at("type").value();
    if (!std::holds_alternative<std::string>(t))
      return false;
    if (std::get<std::string>(t) != *expected_type)
      return false;
  }
  return true;
}

void flatten_validator_allof_walk(rfl::Generic &node)
{
  if (std::holds_alternative<rfl::Generic::Object>(node.value()))
  {
    auto &obj = std::get<rfl::Generic::Object>(node.value());
    if (obj.count("allOf") != 0 &&
        std::holds_alternative<rfl::Generic::Array>(obj.at("allOf").value()))
    {
      auto &arr = std::get<rfl::Generic::Array>(obj.at("allOf").value());
      // Determine the expected `type` across arms — use the first arm's
      // type if any arm declares one.
      std::optional<std::string> expected_type;
      for (const auto &arm : arr)
      {
        if (!std::holds_alternative<rfl::Generic::Object>(arm.value()))
        {
          break;
        }
        const auto &a = std::get<rfl::Generic::Object>(arm.value());
        if (a.count("type") != 0 &&
            std::holds_alternative<std::string>(a.at("type").value()))
        {
          expected_type = std::get<std::string>(a.at("type").value());
          break;
        }
      }
      bool all_mergeable = !arr.empty();
      for (const auto &arm : arr)
      {
        if (!arm_is_mergeable(arm, expected_type))
        {
          all_mergeable = false;
          break;
        }
      }
      if (all_mergeable)
      {
        // Collect the union of constraint keys. `allOf` arms are
        // guaranteed by rfl's emitter to touch disjoint constraint
        // keys, so naive overwrite is safe; but we defensively skip
        // duplicates to preserve determinism if that ever changes.
        rfl::Generic::Object merged;
        for (auto &arm : arr)
        {
          auto &a = std::get<rfl::Generic::Object>(arm.value());
          for (auto &kv : a)
          {
            if (merged.count(kv.first) != 0)
              continue;
            merged.insert(std::move(kv.first), std::move(kv.second));
          }
        }
        // Splice merged keys into the parent object in place of
        // `allOf`. Preserve surrounding keys (description, default,
        // x-palace-*) by iterating the parent and rebuilding with
        // constraints slotted where `allOf` sat.
        rfl::Generic::Object rebuilt;
        for (auto &entry : obj)
        {
          if (entry.first == "allOf")
          {
            for (auto &mkv : merged)
            {
              if (rebuilt.count(mkv.first) != 0)
                continue;
              rebuilt.insert(std::move(mkv.first), std::move(mkv.second));
            }
            continue;
          }
          rebuilt.insert(std::move(entry.first), std::move(entry.second));
        }
        obj = std::move(rebuilt);
      }
    }
    for (auto &[_, v] : obj)
    {
      flatten_validator_allof_walk(v);
    }
  }
  else if (std::holds_alternative<rfl::Generic::Array>(node.value()))
  {
    auto &arr = std::get<rfl::Generic::Array>(node.value());
    for (auto &v : arr)
    {
      flatten_validator_allof_walk(v);
    }
  }
}

}  // namespace

std::string flatten_validator_allof(std::string schema_json)
{
  auto schema_gen_r = rfl::json::read<rfl::Generic>(schema_json);
  if (!schema_gen_r)
    return schema_json;
  flatten_validator_allof_walk(*schema_gen_r);
  return rfl::json::write(*schema_gen_r, rfl::json::pretty);
}

std::string
inject_defaults(std::string schema_json,
                const std::unordered_map<std::string, std::string> &defaults_map)
{
  auto schema_gen_r = rfl::json::read<rfl::Generic>(schema_json);
  if (!schema_gen_r)
  {
    return schema_json;
  }

  auto &schema_var = schema_gen_r->value();
  if (!std::holds_alternative<rfl::Generic::Object>(schema_var))
  {
    return schema_json;
  }
  auto &root_obj = std::get<rfl::Generic::Object>(schema_var);

  // Iterate every `$defs` entry and fill in defaults from the pre-collected
  // map. Entries whose $defs name isn't in the map (types we couldn't reach
  // via the compile-time type-graph walk) are left alone. Keeping the loop
  // uniform — no special case for root $ref vs TaggedUnion anyOf — is what
  // closes the `items.$ref` / TaggedUnion gaps.
  if (root_obj.count("$defs") != 0)
  {
    auto &defs_var = root_obj.at("$defs").value();
    if (std::holds_alternative<rfl::Generic::Object>(defs_var))
    {
      auto &defs = std::get<rfl::Generic::Object>(defs_var);
      for (auto &[def_name, def_generic] : defs)
      {
        const auto it = defaults_map.find(def_name);
        if (it == defaults_map.end())
          continue;
        if (!std::holds_alternative<rfl::Generic::Object>(def_generic.value()))
        {
          continue;
        }
        auto &def_obj = std::get<rfl::Generic::Object>(def_generic.value());

        auto parsed = rfl::json::read<rfl::Generic>(it->second);
        if (!parsed)
          continue;
        if (!std::holds_alternative<rfl::Generic::Object>(parsed->value()))
        {
          continue;
        }
        const auto &parsed_map = std::get<rfl::Generic::Object>(parsed->value());

        inject_properties(def_obj, parsed_map);
      }
    }
  }

  return rfl::json::write(*schema_gen_r, rfl::json::pretty);
}

// Standalone pass equivalent of the old coupled post-step inside
// `inject_defaults`. Decoupled so schemas emitted with `emit_defaults=false`
// can still pick up the `additionalProperties: false` hardening.
std::string inject_additional_properties_false(std::string schema_json)
{
  auto schema_gen_r = rfl::json::read<rfl::Generic>(schema_json);
  if (!schema_gen_r)
    return schema_json;
  inject_additional_properties_false_walk(*schema_gen_r);
  return rfl::json::write(*schema_gen_r, rfl::json::pretty);
}

// Inject arm-level discriminator descriptions into a TaggedUnion schema. For
// each entry in `$defs`, look at `properties[tag_name]`. If that property has
// a single-value `enum` whose value matches one of the `(discriminator, desc)`
// pairs, write the description in as `"description"` on the property body.
// Existing descriptions are preserved so a future direct annotation on the
// literal (should rfl ever allow it) would win over the side-channel one.
std::string
inject_tag_descriptions(std::string schema_json, std::string tag_name,
                        std::vector<std::pair<std::string, std::string>> arm_descs)
{
  auto schema_gen_r = rfl::json::read<rfl::Generic>(schema_json);
  if (!schema_gen_r)
  {
    return schema_json;
  }

  auto &schema_var = schema_gen_r->value();
  if (!std::holds_alternative<rfl::Generic::Object>(schema_var))
  {
    return schema_json;
  }
  auto &root_obj = std::get<rfl::Generic::Object>(schema_var);

  if (root_obj.count("$defs") == 0)
  {
    return schema_json;
  }
  auto &defs_var = root_obj.at("$defs").value();
  if (!std::holds_alternative<rfl::Generic::Object>(defs_var))
  {
    return schema_json;
  }
  auto &defs = std::get<rfl::Generic::Object>(defs_var);

  for (auto &[_def_name, def_g] : defs)
  {
    if (!std::holds_alternative<rfl::Generic::Object>(def_g.value()))
    {
      continue;
    }
    auto &def_obj = std::get<rfl::Generic::Object>(def_g.value());
    if (def_obj.count("properties") == 0)
      continue;
    auto &props_var = def_obj.at("properties").value();
    if (!std::holds_alternative<rfl::Generic::Object>(props_var))
      continue;
    auto &props = std::get<rfl::Generic::Object>(props_var);

    if (props.count(tag_name) == 0)
      continue;
    auto &tag_prop_var = props.at(tag_name).value();
    if (!std::holds_alternative<rfl::Generic::Object>(tag_prop_var))
      continue;
    auto &tag_prop = std::get<rfl::Generic::Object>(tag_prop_var);

    // Require the single-element enum shape rfl emits for
    // rfl::Literal<"Tag"> so we don't accidentally label a non-literal
    // field named `Type`.
    if (tag_prop.count("enum") == 0)
      continue;
    auto &enum_var = tag_prop.at("enum").value();
    if (!std::holds_alternative<rfl::Generic::Array>(enum_var))
      continue;
    auto &enum_arr = std::get<rfl::Generic::Array>(enum_var);
    if (enum_arr.size() != 1)
      continue;
    if (!std::holds_alternative<std::string>(enum_arr[0].value()))
      continue;
    const auto &discriminator = std::get<std::string>(enum_arr[0].value());

    for (const auto &[disc, desc] : arm_descs)
    {
      if (disc != discriminator)
        continue;
      if (tag_prop.count("description") != 0)
        break;
      tag_prop["description"] = rfl::Generic(desc);
      break;
    }
  }

  return rfl::json::write(*schema_gen_r, rfl::json::pretty);
}

// Replace each matching `$defs` entry's `required` array with the explicit
// list supplied in `required_map`. An entry whose list is empty drops the
// `required` key entirely (rather than leaving `"required": []`, which is
// noise). Entries absent from the map are left alone — this only happens
// for `$defs` entries we couldn't reach via the compile-time walker.
std::string
set_required(std::string schema_json,
             const std::unordered_map<std::string, std::vector<std::string>> &required_map)
{
  if (required_map.empty())
    return schema_json;
  auto schema_gen_r = rfl::json::read<rfl::Generic>(schema_json);
  if (!schema_gen_r)
    return schema_json;

  auto &schema_var = schema_gen_r->value();
  if (!std::holds_alternative<rfl::Generic::Object>(schema_var))
  {
    return schema_json;
  }
  auto &root_obj = std::get<rfl::Generic::Object>(schema_var);

  if (root_obj.count("$defs") == 0)
    return schema_json;
  auto &defs_var = root_obj.at("$defs").value();
  if (!std::holds_alternative<rfl::Generic::Object>(defs_var))
  {
    return schema_json;
  }
  auto &defs = std::get<rfl::Generic::Object>(defs_var);

  for (auto &[def_name, def_g] : defs)
  {
    const auto it = required_map.find(def_name);
    if (it == required_map.end())
      continue;
    if (!std::holds_alternative<rfl::Generic::Object>(def_g.value()))
    {
      continue;
    }
    auto &def_obj = std::get<rfl::Generic::Object>(def_g.value());

    if (it->second.empty())
    {
      // Drop the key entirely. rfl::Object has no erase; rebuild
      // without the `required` slot.
      rfl::Generic::Object filtered;
      for (auto &entry : def_obj)
      {
        if (entry.first == "required")
          continue;
        filtered.insert(std::move(entry.first), std::move(entry.second));
      }
      def_obj = std::move(filtered);
      continue;
    }

    rfl::Generic::Array arr;
    for (const auto &n : it->second)
    {
      arr.emplace_back(rfl::Generic(n));
    }
    // operator[] overwrites the existing slot, preserving key order
    // (and inserts if missing — happens when reflect-cpp emitted no
    // required array for this entry, e.g. all fields were optional).
    def_obj["required"] = rfl::Generic(std::move(arr));
  }

  return rfl::json::write(*schema_gen_r, rfl::json::pretty);
}

// --- Schema-driven required-key check (loader pre-flight) -------------------
//
// Handles the subset of JSON Schema keywords palace::schema::utils emits: $ref, properties,
// required, items (array), anyOf (nullable T-or-null), oneOf (discriminated
// union). Deliberately minimal — we're not reimplementing a full validator;
// rfl::json::read handles type/range/enum enforcement. The purpose is to
// catch the "DefaultIfMissing silently filled a required field" case before
// rfl sees it.
namespace
{

const rfl::Generic *resolve_ref(const rfl::Generic::Object &root_obj, std::string_view ref)
{
  constexpr std::string_view prefix = "#/$defs/";
  if (ref.size() <= prefix.size() || ref.substr(0, prefix.size()) != prefix)
  {
    return nullptr;
  }
  auto name = std::string(ref.substr(prefix.size()));
  if (root_obj.count("$defs") == 0)
    return nullptr;
  const auto &defs_var = root_obj.at("$defs").value();
  if (!std::holds_alternative<rfl::Generic::Object>(defs_var))
    return nullptr;
  const auto &defs = std::get<rfl::Generic::Object>(defs_var);
  if (defs.count(name) == 0)
    return nullptr;
  return &defs.at(name);
}

bool is_null_schema(const rfl::Generic &node)
{
  if (!std::holds_alternative<rfl::Generic::Object>(node.value()))
    return false;
  const auto &obj = std::get<rfl::Generic::Object>(node.value());
  if (obj.count("type") == 0)
    return false;
  const auto &t_var = obj.at("type").value();
  if (!std::holds_alternative<std::string>(t_var))
    return false;
  return std::get<std::string>(t_var) == "null";
}

// Extract the tag value, if any, from a schema arm's `properties.<*>.enum`.
// Returns (tag_name, tag_value) for the single property whose body carries a
// one-element enum (the discriminator convention reflect-cpp uses for
// rfl::Literal). Returns nullopt if no such property exists.
std::optional<std::pair<std::string, std::string>>
arm_discriminator(const rfl::Generic::Object &root_obj, const rfl::Generic &arm_schema)
{
  // Resolve $ref if present.
  const rfl::Generic *resolved = &arm_schema;
  if (std::holds_alternative<rfl::Generic::Object>(arm_schema.value()))
  {
    const auto &obj = std::get<rfl::Generic::Object>(arm_schema.value());
    if (obj.count("$ref") != 0)
    {
      const auto &ref_var = obj.at("$ref").value();
      if (std::holds_alternative<std::string>(ref_var))
      {
        const auto *r = resolve_ref(root_obj, std::get<std::string>(ref_var));
        if (r)
          resolved = r;
      }
    }
  }
  if (!std::holds_alternative<rfl::Generic::Object>(resolved->value()))
  {
    return std::nullopt;
  }
  const auto &obj = std::get<rfl::Generic::Object>(resolved->value());
  if (obj.count("properties") == 0)
    return std::nullopt;
  const auto &props_var = obj.at("properties").value();
  if (!std::holds_alternative<rfl::Generic::Object>(props_var))
    return std::nullopt;
  const auto &props = std::get<rfl::Generic::Object>(props_var);
  for (const auto &[name, body] : props)
  {
    if (!std::holds_alternative<rfl::Generic::Object>(body.value()))
      continue;
    const auto &body_obj = std::get<rfl::Generic::Object>(body.value());
    if (body_obj.count("enum") == 0)
      continue;
    const auto &enum_var = body_obj.at("enum").value();
    if (!std::holds_alternative<rfl::Generic::Array>(enum_var))
      continue;
    const auto &enum_arr = std::get<rfl::Generic::Array>(enum_var);
    if (enum_arr.size() != 1)
      continue;
    if (!std::holds_alternative<std::string>(enum_arr[0].value()))
      continue;
    return std::make_pair(name, std::get<std::string>(enum_arr[0].value()));
  }
  return std::nullopt;
}

void check_node(const rfl::Generic::Object &root_obj, const rfl::Generic &schema_node,
                const rfl::Generic &input_node, std::string path,
                std::vector<std::string> &out)
{
  if (!std::holds_alternative<rfl::Generic::Object>(schema_node.value()))
    return;
  const auto &schema_obj = std::get<rfl::Generic::Object>(schema_node.value());

  // $ref: resolve once, then recurse with the referenced entry.
  if (schema_obj.count("$ref") != 0)
  {
    const auto &ref_var = schema_obj.at("$ref").value();
    if (!std::holds_alternative<std::string>(ref_var))
      return;
    const auto *referenced = resolve_ref(root_obj, std::get<std::string>(ref_var));
    if (referenced)
    {
      check_node(root_obj, *referenced, input_node, path, out);
    }
    return;
  }

  // anyOf: palace::schema::utils emits this only for `std::optional<T>` (= [T-schema,
  // null]). If input is null, the null branch matches — nothing to check. Else recurse into
  // the first non-null branch.
  if (schema_obj.count("anyOf") != 0)
  {
    if (std::holds_alternative<std::nullopt_t>(input_node.value()))
      return;
    const auto &arr_var = schema_obj.at("anyOf").value();
    if (!std::holds_alternative<rfl::Generic::Array>(arr_var))
      return;
    const auto &arr = std::get<rfl::Generic::Array>(arr_var);
    for (const auto &branch : arr)
    {
      if (is_null_schema(branch))
        continue;
      check_node(root_obj, branch, input_node, path, out);
      return;
    }
    return;
  }

  // oneOf: discriminated union. Find the arm whose discriminator value
  // matches input's tag field. If none matches, skip (rfl will reject at
  // parse time with a proper tagged-union diagnostic).
  if (schema_obj.count("oneOf") != 0)
  {
    if (!std::holds_alternative<rfl::Generic::Object>(input_node.value()))
      return;
    const auto &input_obj = std::get<rfl::Generic::Object>(input_node.value());
    const auto &arr_var = schema_obj.at("oneOf").value();
    if (!std::holds_alternative<rfl::Generic::Array>(arr_var))
      return;
    const auto &arr = std::get<rfl::Generic::Array>(arr_var);
    for (const auto &branch : arr)
    {
      auto disc = arm_discriminator(root_obj, branch);
      if (!disc)
        continue;
      if (input_obj.count(disc->first) == 0)
        continue;
      const auto &input_tag_var = input_obj.at(disc->first).value();
      if (!std::holds_alternative<std::string>(input_tag_var))
        continue;
      if (std::get<std::string>(input_tag_var) != disc->second)
        continue;
      check_node(root_obj, branch, input_node, path, out);
      return;
    }
    return;
  }

  // items: each array element gets recursed against the items schema.
  if (schema_obj.count("items") != 0)
  {
    if (!std::holds_alternative<rfl::Generic::Array>(input_node.value()))
      return;
    const auto &input_arr = std::get<rfl::Generic::Array>(input_node.value());
    const auto &items = schema_obj.at("items");
    for (std::size_t i = 0; i < input_arr.size(); ++i)
    {
      check_node(root_obj, items, input_arr[i], path + "[" + std::to_string(i) + "]", out);
    }
    return;
  }

  // properties + required: the object case.
  if (schema_obj.count("properties") == 0)
    return;
  if (!std::holds_alternative<rfl::Generic::Object>(input_node.value()))
    return;
  const auto &input_obj = std::get<rfl::Generic::Object>(input_node.value());

  // Missing-required check.
  if (schema_obj.count("required") != 0)
  {
    const auto &req_var = schema_obj.at("required").value();
    if (std::holds_alternative<rfl::Generic::Array>(req_var))
    {
      const auto &req_arr = std::get<rfl::Generic::Array>(req_var);
      for (const auto &r : req_arr)
      {
        if (!std::holds_alternative<std::string>(r.value()))
          continue;
        const auto &name = std::get<std::string>(r.value());
        if (input_obj.count(name) == 0)
        {
          out.push_back(path.empty() ? name : path + "." + name);
        }
      }
    }
  }

  // Recurse into each present property that the schema knows about.
  const auto &props_var = schema_obj.at("properties").value();
  if (!std::holds_alternative<rfl::Generic::Object>(props_var))
    return;
  const auto &props = std::get<rfl::Generic::Object>(props_var);
  for (const auto &[name, input_val] : input_obj)
  {
    if (props.count(name) == 0)
      continue;  // extras — NoExtraFields rejects
    const auto &prop_schema = props.at(name);
    std::string child_path = path.empty() ? name : path + "." + name;
    check_node(root_obj, prop_schema, input_val, std::move(child_path), out);
  }
}

}  // namespace

std::vector<std::string> check_required_keys(const std::string &schema_json,
                                             const std::string &input_json)
{
  auto schema_r = rfl::json::read<rfl::Generic>(schema_json);
  auto input_r = rfl::json::read<rfl::Generic>(input_json);
  if (!schema_r || !input_r)
    return {};
  if (!std::holds_alternative<rfl::Generic::Object>(schema_r->value()))
    return {};
  const auto &root_obj = std::get<rfl::Generic::Object>(schema_r->value());

  std::vector<std::string> out;
  // The root usually carries a `$ref` to a $defs entry; resolve and recurse.
  if (root_obj.count("$ref") != 0)
  {
    const auto &ref_var = root_obj.at("$ref").value();
    if (std::holds_alternative<std::string>(ref_var))
    {
      const auto *referenced = resolve_ref(root_obj, std::get<std::string>(ref_var));
      if (referenced)
      {
        check_node(root_obj, *referenced, *input_r, "", out);
      }
    }
  }
  else
  {
    check_node(root_obj, *schema_r, *input_r, "", out);
  }
  return out;
}

// Rename the root's combiner keyword (e.g. `anyOf` -> `oneOf`). No-op if the
// root doesn't carry the `from` key, so this is safe to call unconditionally
// for types that opt in via `palace::schema::utils::schema_composition<T>`. Preserves the
// surrounding keys (`$schema`, `$defs`, etc.) and the array value unchanged.
std::string rewrite_root_composition(std::string schema_json, std::string from,
                                     std::string to)
{
  auto schema_gen_r = rfl::json::read<rfl::Generic>(schema_json);
  if (!schema_gen_r)
    return schema_json;

  auto &schema_var = schema_gen_r->value();
  if (!std::holds_alternative<rfl::Generic::Object>(schema_var))
    return schema_json;
  auto &root_obj = std::get<rfl::Generic::Object>(schema_var);

  if (root_obj.count(from) == 0)
    return schema_json;

  // rfl::Object preserves insertion order and has no erase. Rename by
  // scanning the underlying pair list and mutating the key in place.
  for (auto &entry : root_obj)
  {
    if (entry.first == from)
    {
      entry.first = to;
      break;
    }
  }

  return rfl::json::write(*schema_gen_r, rfl::json::pretty);
}

// Rewrite an inline-enum property body in place to PR-716's per-value
// oneOf form. Expects `prop_obj` to carry `{"type": "string", "enum":
// [...string-array...]}`. The `pairs` vector supplies (value, description)
// entries; any enum values missing from `pairs` are still emitted as a
// `{"const": value}` arm without a `description`, keeping the schema
// shape-equivalent. `type` and `enum` are dropped; surrounding keys
// (`description`, `default`) are preserved.
namespace
{

bool rewrite_enum_body(rfl::Generic::Object &prop_obj,
                       const std::vector<std::pair<std::string, std::string>> &pairs)
{
  if (prop_obj.count("type") == 0 || prop_obj.count("enum") == 0)
  {
    return false;
  }
  const auto &type_var = prop_obj.at("type").value();
  if (!std::holds_alternative<std::string>(type_var))
    return false;
  if (std::get<std::string>(type_var) != "string")
    return false;

  const auto &enum_var = prop_obj.at("enum").value();
  if (!std::holds_alternative<rfl::Generic::Array>(enum_var))
    return false;
  const auto &enum_arr = std::get<rfl::Generic::Array>(enum_var);

  rfl::Generic::Array one_of;
  for (const auto &v : enum_arr)
  {
    if (!std::holds_alternative<std::string>(v.value()))
      return false;
    const auto &val = std::get<std::string>(v.value());
    rfl::Generic::Object arm;
    arm["const"] = rfl::Generic(val);
    for (const auto &[p_val, p_desc] : pairs)
    {
      if (p_val != val)
        continue;
      // Empty `p_desc` is the caller's opt-out: PALACE_SCHEMA_ENUM
      // forces every enumerator to appear in the description array
      // so the enum body and the descs stay synchronized, and `""`
      // marks "no user-facing description" — skip the key to keep
      // the arm as a bare `{"const": ...}`.
      if (!p_desc.empty())
      {
        arm["description"] = rfl::Generic(p_desc);
      }
      break;
    }
    one_of.emplace_back(std::move(arm));
  }

  // rfl::Object (ordered map) — rename `enum` to `oneOf` in place, drop
  // `type`. Rename first so the new key lands in the same slot the old
  // `enum` occupied, which keeps key ordering stable across rebuilds.
  for (auto &entry : prop_obj)
  {
    if (entry.first == "enum")
    {
      entry.first = "oneOf";
      entry.second = rfl::Generic(std::move(one_of));
      break;
    }
  }
  // Remove `type` by rebuilding the object without it.
  rfl::Generic::Object filtered;
  for (auto &entry : prop_obj)
  {
    if (entry.first == "type")
      continue;
    filtered.insert(std::move(entry.first), std::move(entry.second));
  }
  prop_obj = std::move(filtered);
  return true;
}

}  // namespace

std::string inject_enum_descriptions(std::string schema_json,
                                     const std::vector<EnumDescEntry> &entries)
{
  if (entries.empty())
    return schema_json;
  auto schema_gen_r = rfl::json::read<rfl::Generic>(schema_json);
  if (!schema_gen_r)
    return schema_json;

  auto &schema_var = schema_gen_r->value();
  if (!std::holds_alternative<rfl::Generic::Object>(schema_var))
    return schema_json;
  auto &root_obj = std::get<rfl::Generic::Object>(schema_var);

  if (root_obj.count("$defs") == 0)
    return schema_json;
  auto &defs_var = root_obj.at("$defs").value();
  if (!std::holds_alternative<rfl::Generic::Object>(defs_var))
    return schema_json;
  auto &defs = std::get<rfl::Generic::Object>(defs_var);

  for (const auto &entry : entries)
  {
    if (defs.count(entry.struct_name) == 0)
      continue;
    auto &def_var = defs.at(entry.struct_name).value();
    if (!std::holds_alternative<rfl::Generic::Object>(def_var))
      continue;
    auto &def_obj = std::get<rfl::Generic::Object>(def_var);
    if (def_obj.count("properties") == 0)
      continue;
    auto &props_var = def_obj.at("properties").value();
    if (!std::holds_alternative<rfl::Generic::Object>(props_var))
      continue;
    auto &props = std::get<rfl::Generic::Object>(props_var);
    if (props.count(entry.field_name) == 0)
      continue;
    auto &prop_var = props.at(entry.field_name).value();
    if (!std::holds_alternative<rfl::Generic::Object>(prop_var))
      continue;
    auto &prop_obj = std::get<rfl::Generic::Object>(prop_var);
    rewrite_enum_body(prop_obj, entry.pairs);
  }

  return rfl::json::write(*schema_gen_r, rfl::json::pretty);
}

// Splice a `"oneOf": [{"required": [...]}, ...]` block into each listed
// `$defs` entry. This coexists with the entry's ordinary `required`
// array — the two apply independently: every item in `required` is
// universally required, while `oneOf` demands exactly one arm be
// satisfied. Typical use is "one of (Attributes | Elements)" on
// LumpedPort / SurfaceCurrent where the struct body can be described
// two ways, and PR 716's schema enforces the choice at validation time.
std::string inject_oneof_required(std::string schema_json,
                                  const std::vector<OneOfRequired> &entries)
{
  if (entries.empty())
    return schema_json;
  auto schema_gen_r = rfl::json::read<rfl::Generic>(schema_json);
  if (!schema_gen_r)
    return schema_json;

  auto &schema_var = schema_gen_r->value();
  if (!std::holds_alternative<rfl::Generic::Object>(schema_var))
    return schema_json;
  auto &root_obj = std::get<rfl::Generic::Object>(schema_var);

  if (root_obj.count("$defs") == 0)
    return schema_json;
  auto &defs_var = root_obj.at("$defs").value();
  if (!std::holds_alternative<rfl::Generic::Object>(defs_var))
    return schema_json;
  auto &defs = std::get<rfl::Generic::Object>(defs_var);

  for (const auto &e : entries)
  {
    if (defs.count(e.struct_name) == 0)
      continue;
    auto &def_var = defs.at(e.struct_name).value();
    if (!std::holds_alternative<rfl::Generic::Object>(def_var))
      continue;
    auto &def_obj = std::get<rfl::Generic::Object>(def_var);

    rfl::Generic::Array one_of;
    for (const auto &arm : e.arms)
    {
      rfl::Generic::Object arm_obj;
      rfl::Generic::Array req_arr;
      for (const auto &name : arm)
      {
        req_arr.emplace_back(rfl::Generic(name));
      }
      arm_obj["required"] = rfl::Generic(std::move(req_arr));
      one_of.emplace_back(std::move(arm_obj));
    }
    def_obj["oneOf"] = rfl::Generic(std::move(one_of));
  }

  return rfl::json::write(*schema_gen_r, rfl::json::pretty);
}

// Stamp `"x-palace-<flag>": true` on each listed property. `flag` is
// typically `"advanced"` or `"deprecated"`; any other value rides along
// under the `x-palace-` prefix (JSON Schema's reserved vendor-extension
// namespace, so validators ignore them).
std::string inject_custom_keywords(std::string schema_json,
                                   const std::vector<FieldFlag> &flags)
{
  if (flags.empty())
    return schema_json;
  auto schema_gen_r = rfl::json::read<rfl::Generic>(schema_json);
  if (!schema_gen_r)
    return schema_json;

  auto &schema_var = schema_gen_r->value();
  if (!std::holds_alternative<rfl::Generic::Object>(schema_var))
    return schema_json;
  auto &root_obj = std::get<rfl::Generic::Object>(schema_var);

  if (root_obj.count("$defs") == 0)
    return schema_json;
  auto &defs_var = root_obj.at("$defs").value();
  if (!std::holds_alternative<rfl::Generic::Object>(defs_var))
    return schema_json;
  auto &defs = std::get<rfl::Generic::Object>(defs_var);

  for (const auto &f : flags)
  {
    if (defs.count(f.struct_name) == 0)
      continue;
    auto &def_var = defs.at(f.struct_name).value();
    if (!std::holds_alternative<rfl::Generic::Object>(def_var))
      continue;
    auto &def_obj = std::get<rfl::Generic::Object>(def_var);
    if (def_obj.count("properties") == 0)
      continue;
    auto &props_var = def_obj.at("properties").value();
    if (!std::holds_alternative<rfl::Generic::Object>(props_var))
      continue;
    auto &props = std::get<rfl::Generic::Object>(props_var);
    if (props.count(f.field_name) == 0)
      continue;
    auto &prop_var = props.at(f.field_name).value();
    if (!std::holds_alternative<rfl::Generic::Object>(prop_var))
      continue;
    auto &prop_obj = std::get<rfl::Generic::Object>(prop_var);
    prop_obj["x-palace-" + f.flag] = rfl::Generic(true);
  }

  return rfl::json::write(*schema_gen_r, rfl::json::pretty);
}

// Stamp `"version": <value>` on the root object. JSON Schema allows custom
// root keywords, so this is a lightweight schema-level annotation (it does
// not affect instance validation).
std::string inject_root_version(std::string schema_json, std::string version)
{
  auto schema_gen_r = rfl::json::read<rfl::Generic>(schema_json);
  if (!schema_gen_r)
    return schema_json;

  auto &schema_var = schema_gen_r->value();
  if (!std::holds_alternative<rfl::Generic::Object>(schema_var))
    return schema_json;
  auto &root_obj = std::get<rfl::Generic::Object>(schema_var);

  root_obj["version"] = rfl::Generic(std::move(version));

  return rfl::json::write(*schema_gen_r, rfl::json::pretty);
}

// Recursively walk a JSON tree and rewrite every `$ref` string that begins
// with `"#/$defs/<prefix>"` by removing `<prefix>`. Leaves all other
// strings alone. Used to keep `$ref` targets in lock-step with the
// `$defs` key renames done by `strip_defs_prefix`.
namespace
{
void rewrite_refs(rfl::Generic &node, const std::string &prefix)
{
  const std::string needle = "#/$defs/" + prefix;
  if (std::holds_alternative<rfl::Generic::Object>(node.value()))
  {
    auto &obj = std::get<rfl::Generic::Object>(node.value());
    for (auto &entry : obj)
    {
      if (entry.first == "$ref" &&
          std::holds_alternative<std::string>(entry.second.value()))
      {
        auto &s = std::get<std::string>(entry.second.value());
        if (s.rfind(needle, 0) == 0)
        {
          s = "#/$defs/" + s.substr(needle.size());
        }
      }
      rewrite_refs(entry.second, prefix);
    }
  }
  else if (std::holds_alternative<rfl::Generic::Array>(node.value()))
  {
    for (auto &v : std::get<rfl::Generic::Array>(node.value()))
    {
      rewrite_refs(v, prefix);
    }
  }
}
}  // namespace

// Strip `prefix` from every `$defs` entry name and matching `$ref` string.
// reflect-cpp emits `$defs` keys from `rfl::parsing::make_type_name<T>()`,
// which renders `palace::schema::Foo` as `palace__schema__Foo`. Callers
// pass `"palace__schema__"` to collapse those to plain `Foo`. Entries not
// starting with the prefix, and refs pointing elsewhere, are untouched.
std::string strip_defs_prefix(std::string schema_json, std::string prefix)
{
  if (prefix.empty())
    return schema_json;
  auto schema_gen_r = rfl::json::read<rfl::Generic>(schema_json);
  if (!schema_gen_r)
    return schema_json;

  auto &schema_var = schema_gen_r->value();
  if (!std::holds_alternative<rfl::Generic::Object>(schema_var))
    return schema_json;
  auto &root_obj = std::get<rfl::Generic::Object>(schema_var);

  if (root_obj.count("$defs") != 0)
  {
    auto &defs_var = root_obj.at("$defs").value();
    if (std::holds_alternative<rfl::Generic::Object>(defs_var))
    {
      auto &defs = std::get<rfl::Generic::Object>(defs_var);
      // rfl::Object preserves insertion order and exposes pair access via
      // its iterator — rename keys in place so the relative order of
      // definitions is preserved across the rewrite.
      for (auto &entry : defs)
      {
        if (entry.first.rfind(prefix, 0) == 0)
        {
          entry.first = entry.first.substr(prefix.size());
        }
      }
    }
  }

  rewrite_refs(*schema_gen_r, prefix);

  return rfl::json::write(*schema_gen_r, rfl::json::pretty);
}

}  // namespace palace::schema::utils::detail
