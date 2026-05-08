// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_SCHEMA_UTILS_SCHEMA_IMPL_HPP
#define PALACE_SCHEMA_UTILS_SCHEMA_IMPL_HPP

#include <array>
#include <optional>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <rfl.hpp>
#include <rfl/DefaultVal.hpp>
#include <rfl/Description.hpp>
#include <rfl/Rename.hpp>
#include <rfl/TaggedUnion.hpp>
#include <rfl/Tuple.hpp>
#include <rfl/Validator.hpp>
#include <rfl/json.hpp>
#include <rfl/named_tuple_t.hpp>
#include <rfl/internal/StringLiteral.hpp>
#include <rfl/internal/tag_t.hpp>
#include <rfl/parsing/make_type_name.hpp>

#include "annotations.hpp"

namespace palace::schema::utils::detail
{

// Defined in src/schema.cpp. Takes the raw reflect-cpp schema string plus a
// map of `$defs`-key-name -> serialized default JSON for that entry's type.
// For each `$defs` entry whose key appears in the map, walks the entry's
// `properties` in lock-step with the parsed defaults and injects `"default"`
// at each matching leaf. Entries not in the map (or properties we can't match)
// are left alone except for the optional-field `"default": null` fallback.
std::string
inject_defaults(std::string schema_json,
                const std::unordered_map<std::string, std::string> &defaults_map);

// Defined in src/schema.cpp. Walks the whole schema tree and stamps
// `"additionalProperties": false` onto every object that declares
// `properties` and doesn't already set the key. reflect-cpp omits this
// keyword entirely, so external validators (Python jsonschema, etc.) would
// otherwise accept unknown keys. Independent of the default-injection pass so
// it can run in schemas emitted without defaults too.
std::string inject_additional_properties_false(std::string schema_json);

// Defined in src/schema.cpp. Collapses multi-constraint `allOf` blocks
// produced by chained rfl::Validator<T, C1, C2, ...> emissions into a
// single flat body. E.g. `Closed<int, 1, 2>` originally serialises as
// `{allOf: [{type:"integer", minimum:1}, {type:"integer", maximum:2}]}`
// and becomes `{type:"integer", minimum:1, maximum:2}`. Arms that mix
// non-constraint keywords or disagree on `type` are left untouched — the
// pass only flattens the specific shape rfl's Validator emitter
// produces.
std::string flatten_validator_allof(std::string schema_json);

// Defined in src/schema.cpp. Walks `$defs` and replaces each entry's
// `required` array with the explicit list supplied in `required_map`. An
// entry absent from the map (no fields marked required) gets `required`
// dropped entirely. This implements the "required ⟺ explicitly marked
// via PALACE_SCHEMA_DESC_REQUIRED or PALACE_SCHEMA_TAG" rule, inverting
// reflect-cpp's default-everything-required behaviour.
std::string
set_required(std::string schema_json,
             const std::unordered_map<std::string, std::vector<std::string>> &required_map);

// Defined in src/schema.cpp. Walks a palace::schema::utils-emitted schema alongside an
// input JSON document, returning dotted paths of every key listed in a `required` array
// that is absent from its corresponding object. Handles $ref navigation, nullable `anyOf:
// [T, null]`, discriminated `oneOf` via discriminator lookup, and array `items` recursion.
// Other keywords (type, minimum, enum, etc.) are deferred to rfl::json::read which enforces
// them at parse time.
//
// Returns an empty vector when the input has every required key at every
// level. This is the bridge that lets `palace::schema::utils::load<T>` reject JSON the
// loader's `rfl::DefaultIfMissing` would otherwise silently fill.
std::vector<std::string> check_required_keys(const std::string &schema_json,
                                             const std::string &input_json);

// Defined in src/schema.cpp. Walks `$defs` and, for each def whose
// `properties[tag_name]` single-element `enum` matches one of `arm_descs`'
// discriminator values, writes a `"description"` onto that property body.
// Existing descriptions are preserved (the injection only fills empty slots).
std::string
inject_tag_descriptions(std::string schema_json, std::string tag_name,
                        std::vector<std::pair<std::string, std::string>> arm_descs);

// Defined in src/schema.cpp. Renames the root's combiner keyword (e.g.
// `anyOf` to `oneOf`) when the user opts in via
// `palace::schema::utils::schema_composition<T>`. No-op if the root does not contain the
// `from` key.
std::string rewrite_root_composition(std::string schema_json, std::string from,
                                     std::string to);

// Defined in src/schema.cpp. Writes `"$id": "urn:palace:schema:<version>"`
// onto the root of the schema, placed immediately after `$schema`.
// Overwrites any existing `$id` so callers can re-run the pass safely.
// Skipped by the caller when `version` is empty.
std::string inject_root_id(std::string schema_json, std::string version);

// Defined in src/schema.cpp. Strips `prefix` from every `$defs` entry name
// and rewrites any `$ref` string pointing at those entries to match. Other
// strings in the document are untouched — the rename operates only on the
// `$defs` map's keys and `$ref` string values. No-op when `prefix` is
// empty or when no entries start with it.
std::string strip_defs_prefix(std::string schema_json, std::string prefix);

// Defined in src/schema.cpp. For each `(struct_name, field_name, pairs)`
// triple, locates `$defs[struct_name]/properties[field_name]`. If the body
// matches reflect-cpp's inline-enum shape (`{"type": "string", "enum":
// [...]}`), it is rewritten to PR-716's `{"oneOf": [{"const": v,
// "description": d}, ...]}` form using `pairs` as the (value, description)
// map. Outer `description` / `default` / other keys are preserved; `type`
// and `enum` are removed since `oneOf` subsumes them.
struct EnumDescEntry
{
  std::string struct_name;
  std::string field_name;
  std::vector<std::pair<std::string, std::string>> pairs;
};
std::string inject_enum_descriptions(std::string schema_json,
                                     const std::vector<EnumDescEntry> &entries);

// Defined in src/schema.cpp. For each `(struct_name, field_name, flag)`
// triple, stamps `"x-palace-<flag>": true` onto the matching property body
// under `$defs[struct_name]/properties[field_name]`. Flag strings are
// `"advanced"` / `"deprecated"` in practice; any other string is accepted
// and written through verbatim.
struct FieldFlag
{
  std::string struct_name;
  std::string field_name;
  std::string flag;
};
std::string inject_custom_keywords(std::string schema_json,
                                   const std::vector<FieldFlag> &flags);

// Defined in src/schema.cpp. For each entry, splices
// `{"oneOf": [{"required": arms[0]}, {"required": arms[1]}, ...]}` into
// `$defs[struct_name]`. Used to express "this struct must satisfy exactly
// one of these required-field alternatives" (e.g. LumpedPort: either
// `Attributes` or `Elements`). Coexists with the plain `required` list
// — the `oneOf` is added alongside, not replacing.
struct OneOfRequired
{
  std::string struct_name;
  std::vector<std::vector<std::string>> arms;
};
std::string inject_oneof_required(std::string schema_json,
                                  const std::vector<OneOfRequired> &entries);

// Defined in src/schema.cpp. For each entry, hoists
// `$defs[struct_name]/properties[field_name]` into a shared
// `$defs[alias_name]` entry (keeping only the schema-describing keys —
// `description` / `default` / `title` / `examples` stay at the field site)
// and replaces the field body with `{"$ref": "#/$defs/<alias>"}`. The first
// occurrence of each alias supplies the canonical body; later occurrences
// are just rewritten to `$ref`. Runs after every body-mutating pass
// (defaults, required, enum descriptions, custom keywords) so the hoisted
// body reflects the final property shape.
struct FieldAlias
{
  std::string struct_name;
  std::string field_name;
  std::string alias_name;
};
std::string inject_field_aliases(std::string schema_json,
                                 const std::vector<FieldAlias> &entries);

// Defined in src/schema.cpp. Counterpart to `inject_field_aliases` for
// `rfl::Variant<...>` fields whose alternatives expand inline as
// `anyOf: [{<arm0 body>}, {<arm1 body>}, ...]`. For each entry, locates
// `$defs[struct_name]/properties[field_name]/anyOf[arm_index]`, moves the
// arm's body into `$defs[alias_name]` (first occurrence wins), and replaces
// the arm with `{"$ref": "#/$defs/<alias_name>"}`. Sibling keys on the
// outer property (description, default) stay where they are — only the
// arm body is rewritten.
struct VariantArmAlias
{
  std::string struct_name;
  std::string field_name;
  std::size_t arm_index;
  std::string alias_name;
};
std::string inject_variant_arm_aliases(std::string schema_json,
                                       const std::vector<VariantArmAlias> &entries);

// Defined in src/schema.cpp. Renames the combiner keyword on a specific
// nested property body. Each entry locates
// `$defs[struct_name]/properties[field_name]` (descending through `items`
// for array-of-variant fields) and replaces the `from` key with `to` if
// present. Used to flip the inline `rfl::Variant` composition from rfl's
// default `anyOf` to PR-716's `oneOf`. The whole-schema variant of this
// pass (`rewrite_root_composition`) only touches the root, so nested
// variants need a separate field-targeted rewrite.
struct FieldCompositionRewrite
{
  std::string struct_name;
  std::string field_name;
  std::string from;
  std::string to;
};
std::string inject_field_composition(std::string schema_json,
                                     const std::vector<FieldCompositionRewrite> &entries);

// Compile-time dispatch helper: pattern-matches a `rfl::TaggedUnion<tag,
// Arms...>` via a pointer parameter so we can fold over the arm pack and pull
// each arm's discriminator literal out of `rfl::internal::tag_t`. Only arms
// that provide a `tag_description` static member (via CJS_TAG) contribute a
// pair; if no arm has one, we skip the schema round-trip entirely.
template <rfl::internal::StringLiteral tag, class... Arms>
std::string dispatch_tag_descriptions(std::string s, rfl::TaggedUnion<tag, Arms...> *)
{
  std::vector<std::pair<std::string, std::string>> pairs;
  auto push = [&]<class Arm>()
  {
    if constexpr (has_tag_description<Arm>)
    {
      // rfl::internal::tag_t<discriminator, Arm> resolves to the
      // rfl::Literal<"Tag"> type that `make_tag` would return for this
      // arm, which in our case is the default-constructible single-
      // field Literal for the arm's discriminator. `.name()` yields the
      // tag string itself ("Point", "Linear", "Log").
      using Lit = rfl::internal::tag_t<tag, Arm>;
      pairs.emplace_back(Lit{}.name(), std::string(Arm::tag_description));
    }
  };
  (push.template operator()<Arms>(), ...);
  if (pairs.empty())
    return s;
  return inject_tag_descriptions(std::move(s), tag.str(), std::move(pairs));
}

template <class T>
std::string inject_tag_descriptions_for(std::string s)
{
  return dispatch_tag_descriptions(std::move(s), static_cast<T *>(nullptr));
}

// --- Type-graph walk to collect per-$defs-entry defaults ---------------------
//
// reflect-cpp emits one `$defs` entry per distinct aggregate struct type
// reachable from T. To inject `"default": ...` at every `$defs` entry's
// properties — not just the root — we walk the reachable type graph at compile
// time, and for each struct type record a pair of:
//
//   (rfl::parsing::make_type_name<Type>(),  rfl::json::write(Type{}))
//
// i.e. the string reflect-cpp uses as the `$defs` key, and the JSON dump of
// the default-constructed instance. The runtime pass in src/schema.cpp then
// iterates `$defs` and, for each entry whose name is in the map, parses the
// matching JSON and walks the entry's `properties` in lock-step.
//
// Walking rules:
//   * Strip annotation wrappers first: rfl::Description<text, U>,
//     rfl::Rename<name, U>, rfl::Validator<U, ...>, rfl::DefaultVal<U>.
//   * Primitives / enums / strings / rfl::Literal → stop.
//   * std::vector<U> / std::array<U, N> / std::optional<U> → recurse into U.
//   * rfl::TaggedUnion<tag, Arms...> → recurse into each Arm (the union itself
//     never becomes a `$defs` entry — only the arms do).
//   * Any other class/struct type → treat as an rfl-reflectable aggregate:
//     add `(make_type_name<T>(), json::write(T{}))` to the map, then recurse
//     into each field's stripped type. Requires `T{}` to be default
//     constructible, which is already a precondition of the existing pass.
//
// Cycles are broken by a seen-set keyed on the $defs name: self-referential
// or mutually-recursive types get visited exactly once.

// Strip annotation wrappers recursively: peel Description/Rename/Validator/
// DefaultVal layers until the underlying type surfaces. The wrappers nest in
// arbitrary orders (Validator inside Description is the common case, see
// Problem::Verbose), so the recursion has to continue until the type
// stops matching any wrapper specialization.
template <class T>
struct strip_ann
{
  using type = T;
};

template <rfl::internal::StringLiteral name, class T>
struct strip_ann<rfl::Description<name, T>>
{
  using type = typename strip_ann<T>::type;
};

// Peel Palace-specific Description subtypes the same way — the walker
// keeps the outer type for trait-based `is_desc_required_v` / `is_desc_advanced_v`
// detection, but everywhere `Clean` is computed we want the underlying
// value type.
template <rfl::internal::StringLiteral name, class T>
struct strip_ann<DescRequired<name, T>>
{
  using type = typename strip_ann<T>::type;
};

template <rfl::internal::StringLiteral name, class T>
struct strip_ann<DescAdvanced<name, T>>
{
  using type = typename strip_ann<T>::type;
};

template <rfl::internal::StringLiteral name, class T>
struct strip_ann<DescDeprecated<name, T>>
{
  using type = typename strip_ann<T>::type;
};

template <rfl::internal::StringLiteral name, class T>
struct strip_ann<rfl::Rename<name, T>>
{
  using type = typename strip_ann<T>::type;
};

template <class T, class V, class... Vs>
struct strip_ann<rfl::Validator<T, V, Vs...>>
{
  using type = typename strip_ann<T>::type;
};

template <class T>
struct strip_ann<rfl::DefaultVal<T>>
{
  using type = typename strip_ann<T>::type;
};

template <class T>
using strip_ann_t = typename strip_ann<std::remove_cvref_t<T>>::type;

// Shape traits for container/union types we unwrap without treating as a
// `$defs` entry in their own right. These match the shapes reflect-cpp emits
// inline (arrays, anyOf-wrapped optionals, TaggedUnion `anyOf`/`oneOf`).
template <class>
struct is_vector : std::false_type
{
};
template <class U, class A>
struct is_vector<std::vector<U, A>> : std::true_type
{
};

template <class>
struct is_std_array : std::false_type
{
};
template <class U, std::size_t N>
struct is_std_array<std::array<U, N>> : std::true_type
{
};

template <class>
struct is_optional : std::false_type
{
};
template <class U>
struct is_optional<std::optional<U>> : std::true_type
{
};

template <class>
struct is_rfl_tagged_union : std::false_type
{
};
template <rfl::internal::StringLiteral tag, class... Arms>
struct is_rfl_tagged_union<rfl::TaggedUnion<tag, Arms...>> : std::true_type
{
};

// Peel one layer of `std::vector<U>` / `std::array<U, N>` to get the
// element type. Lazy: when T has no inner element type, `type` resolves
// to T itself, so callers can probe the wrapped variant without
// force-instantiating `T::value_type` on non-container types (which is
// what `std::conditional_t` would do).
template <class T>
struct container_inner
{
  using type = T;
};
template <class U, class A>
struct container_inner<std::vector<U, A>>
{
  using type = U;
};
template <class U, std::size_t N>
struct container_inner<std::array<U, N>>
{
  using type = U;
};
template <class T>
using container_inner_t = typename container_inner<T>::type;

// rfl::Variant<Alts...> is emitted inline by reflect-cpp as `anyOf: [...]`;
// it never becomes a `$defs` entry of its own. The walker must recurse into
// each alternative to collect defaults/required/etc. from reachable structs,
// but must not treat the Variant itself as a reflectable aggregate (it is
// not default-schema-emittable as a struct, and default-constructing one
// would invoke the first alternative's default ctor, which for a validated
// string would validate an empty value).
template <class>
struct is_rfl_variant : std::false_type
{
};
template <class... Alts>
struct is_rfl_variant<rfl::Variant<Alts...>> : std::true_type
{
};

template <class>
struct is_rfl_literal : std::false_type
{
};
template <rfl::internal::StringLiteral... names>
struct is_rfl_literal<rfl::Literal<names...>> : std::true_type
{
};

// Accumulator for the type-graph walk. Plain struct so the recursion can pass
// it by reference without per-call template-parameter noise.
//
// `required_map` records the names of fields explicitly marked required
// (via `PALACE_SCHEMA_DESC_REQUIRED` or `PALACE_SCHEMA_TAG`). Keyed by
// `$defs` name; the runtime `set_required` pass uses this to REPLACE
// reflect-cpp's native required array — anything not listed becomes
// optional.
struct DefaultsAccum
{
  std::unordered_map<std::string, std::string> entries;
  std::unordered_map<std::string, std::vector<std::string>> required_map;
  // Per-enum-value descriptions and per-field custom-keyword flags are
  // harvested while walking the same type graph so consumers can run a
  // single pass. `enum_descs` is flattened into one entry per
  // (struct, enum-field) pair; `field_flags` is one entry per flagged
  // field.
  std::vector<EnumDescEntry> enum_descs;
  std::vector<FieldFlag> field_flags;
  // Cross-field `oneOf` required alternatives — one entry per struct
  // type that specializes `schema_oneof_required<T>`.
  std::vector<OneOfRequired> oneof_required;
  // Named-alias hoist sites — one entry per field whose cleaned type
  // specializes `schema_alias_name<T>`.
  std::vector<FieldAlias> field_aliases;
  // Per-arm alias hoist sites for `rfl::Variant<...>` fields: one entry
  // per (struct, field, arm-index) whose cleaned arm type specializes
  // `schema_alias_name<T>`.
  std::vector<VariantArmAlias> variant_arm_aliases;
  // Field-level composition rewrites for `rfl::Variant<...>` fields whose
  // variant type opts into `Compose::OneOf` via `schema_composition<U>`.
  // The post-emit pass flips the inline `anyOf` to `oneOf` at the
  // matching field site (descending through `items` for array fields).
  std::vector<FieldCompositionRewrite> field_compositions;
  // One entry per `rfl::TaggedUnion<tag, Arms...>` reached during the
  // type-graph walk. Each entry carries the discriminator field name
  // and the (tag-value, description) pairs harvested from arms whose
  // type provides `tag_description`. Used by the runtime
  // `inject_tag_descriptions` pass to rewrite every arm's discriminator
  // property body from `{type:string, enum:[X], default:X}` into
  // `{const:X, description:...}` regardless of nesting depth.
  struct TaggedUnionDispatch
  {
    std::string tag_name;
    std::vector<std::pair<std::string, std::string>> arm_descs;
  };
  std::vector<TaggedUnionDispatch> tagged_unions;
  std::unordered_set<std::string> seen;
};

// Forward declaration — mutually recursive with walk_struct_fields.
template <class T>
void visit_type(DefaultsAccum &acc);

// Walk the field-value types of a reflect-cpp-reflectable struct. `Values` is
// `rfl::Tuple<T1, T2, ...>` where each Ti is the struct's declared member
// type (possibly still wrapped in Description/Validator/etc.). We strip the
// wrappers and recurse. This helper exists as a struct so it can be partially
// specialized on rfl::Tuple<...>; function templates cannot do that.
template <class Tup>
struct walk_tuple;

template <class... Ts>
struct walk_tuple<rfl::Tuple<Ts...>>
{
  static void apply(DefaultsAccum &acc) { (visit_type<strip_ann_t<Ts>>(acc), ...); }
};

template <class T>
void walk_struct_fields(DefaultsAccum &acc)
{
  using Values = typename rfl::named_tuple_t<T>::Values;
  walk_tuple<Values>::apply(acc);
}

// Per-field walker that collects three kinds of annotation in one pass:
//
//   * Enum-description entries: for every `rfl::Field<name, Type>` whose
//     cleaned value type is an enum E with non-empty
//     `enum_descriptions<E>`, emit an EnumDescEntry keyed on U's `$defs`
//     name and the field name.
//
//   * Custom-keyword flags: for every field declared via
//     `PALACE_SCHEMA_DESC_ADVANCED` / `PALACE_SCHEMA_DESC_DEPRECATED`, the
//     field's declared type is `DescAdvanced<>` / `DescDeprecated<>` —
//     detected via `is_desc_advanced_v` / `is_desc_deprecated_v` traits.
//
//   * Required-field marks: a field declared via
//     `PALACE_SCHEMA_DESC_REQUIRED` has declared type `DescRequired<>`,
//     detected via `is_desc_required_v`. Tagged-union discriminators
//     declared via `PALACE_SCHEMA_TAG` have type `rfl::Literal<Value>`,
//     and the macro emits a hidden-friend `palace_schema_field_required`
//     keyed on `const rfl::Literal<Value>*` — probed via ADL using the
//     cleaned field type. The Literal type is unique per arm (different
//     `Value`), so multiple arms' friends never collide at namespace
//     scope.
template <class U, class FieldsTup>
struct field_annotations_walker;

template <class U, class... Fs>
struct field_annotations_walker<U, rfl::Tuple<Fs...>>
{
  static void apply(const std::string &struct_name, std::vector<EnumDescEntry> &enum_out,
                    std::vector<FieldFlag> &flag_out,
                    std::vector<std::string> &required_out,
                    std::vector<FieldAlias> &alias_out,
                    std::vector<VariantArmAlias> &variant_arm_alias_out,
                    std::vector<FieldCompositionRewrite> &field_composition_out)
  {
    (check_one<Fs>(struct_name, enum_out, flag_out, required_out, alias_out,
                   variant_arm_alias_out, field_composition_out),
     ...);
  }

  template <class F>
  static void
  check_one(const std::string &struct_name, std::vector<EnumDescEntry> &enum_out,
            std::vector<FieldFlag> &flag_out, std::vector<std::string> &required_out,
            std::vector<FieldAlias> &alias_out,
            std::vector<VariantArmAlias> &variant_arm_alias_out,
            std::vector<FieldCompositionRewrite> &field_composition_out)
  {
    check_one_impl(struct_name, enum_out, flag_out, required_out, alias_out,
                   variant_arm_alias_out, field_composition_out, static_cast<F *>(nullptr));
  }

  template <rfl::internal::StringLiteral name, class Type>
  static void
  check_one_impl(const std::string &struct_name, std::vector<EnumDescEntry> &enum_out,
                 std::vector<FieldFlag> &flag_out, std::vector<std::string> &required_out,
                 std::vector<FieldAlias> &alias_out,
                 std::vector<VariantArmAlias> &variant_arm_alias_out,
                 std::vector<FieldCompositionRewrite> &field_composition_out,
                 rfl::Field<name, Type> *)
  {
    using Raw = std::remove_cvref_t<Type>;
    using Clean = strip_ann_t<Raw>;
    // The enum-description rewrite handles direct enum properties,
    // which covers every PR-716 per-value-description site.
    if constexpr (std::is_enum_v<Clean>)
    {
      if constexpr (enum_descriptions<Clean>::value.size() != 0)
      {
        EnumDescEntry entry;
        entry.struct_name = struct_name;
        entry.field_name.assign(std::string_view{name.str()});
        for (auto [v, d] : enum_descriptions<Clean>::value)
        {
          entry.pairs.emplace_back(std::string(v), std::string(d));
        }
        enum_out.push_back(std::move(entry));
      }
    }
    // Advanced / deprecated flags: declared type IS the wrapper subclass,
    // so detection is a pure type-trait check on the un-stripped type.
    if constexpr (is_desc_advanced<Raw>::value)
    {
      FieldFlag f;
      f.struct_name = struct_name;
      f.field_name.assign(std::string_view{name.str()});
      f.flag = "advanced";
      flag_out.push_back(std::move(f));
    }
    if constexpr (is_desc_deprecated<Raw>::value)
    {
      FieldFlag f;
      f.struct_name = struct_name;
      f.field_name.assign(std::string_view{name.str()});
      f.flag = "deprecated";
      flag_out.push_back(std::move(f));
    }
    // Named-alias hoist: record a FieldAlias entry for every field whose
    // cleaned type T has a non-empty `schema_alias_name<T>` specialization.
    // The post-emit pass uses this to promote inline-emitted bodies into
    // shared `$defs` entries — e.g. `std::array<double, 3>` becomes
    // `$defs/Vector3`.
    if constexpr (!schema_alias_name<Clean>::value.empty())
    {
      FieldAlias a;
      a.struct_name = struct_name;
      a.field_name.assign(std::string_view{name.str()});
      a.alias_name.assign(schema_alias_name<Clean>::value);
      alias_out.push_back(std::move(a));
    }
    // Per-arm alias / composition hoist: when the cleaned field type
    // *contains* a variant (directly or through a `std::vector` /
    // `std::array` wrapper, since reflect-cpp emits the variant inline at
    // every container element), record per-arm aliases and an optional
    // `anyOf → oneOf` field-composition rewrite for the inline body.
    //
    // Used to wire:
    //   * `Direction` variants → shared `$defs/PortDirection` /
    //     `$defs/DipoleDirection` / `$defs/Vector3` entries.
    //   * `LumpedPort` / `SurfaceCurrent` array elements → arm-specific
    //     `$defs/LumpedPortAttributes` etc., with `oneOf` semantics.
    {
      // Strip one layer of std::vector<>/std::array<> if present so the
      // inner type can be inspected as a variant. `std::conditional_t`
      // would force-instantiate `Clean::value_type` even on the false
      // branch, so dispatch through a helper trait (`container_inner_t`)
      // defined elsewhere in this header that resolves to `Clean` when
      // it isn't a container.
      using Inner = strip_ann_t<container_inner_t<Clean>>;
      if constexpr (is_rfl_variant<Inner>::value)
      {
        []<class... Alts>(const std::string &sn, std::string_view fn,
                          std::vector<VariantArmAlias> &out, rfl::Variant<Alts...> *)
        {
          std::size_t idx = 0;
          auto emit_one = [&]<class Alt>()
          {
            using AltClean = strip_ann_t<Alt>;
            if constexpr (!schema_alias_name<AltClean>::value.empty())
            {
              VariantArmAlias v;
              v.struct_name = sn;
              v.field_name.assign(fn);
              v.arm_index = idx;
              v.alias_name.assign(schema_alias_name<AltClean>::value);
              out.push_back(std::move(v));
            }
            ++idx;
          };
          (emit_one.template operator()<Alts>(), ...);
        }(struct_name, std::string_view{name.str()}, variant_arm_alias_out,
          static_cast<Inner *>(nullptr));
        // Field-level composition: if the variant alias opted into
        // `Compose::OneOf`, queue an `anyOf → oneOf` rewrite at this
        // field site.
        if constexpr (schema_composition<Inner>::value == Compose::OneOf)
        {
          FieldCompositionRewrite c;
          c.struct_name = struct_name;
          c.field_name.assign(std::string_view{name.str()});
          c.from = "anyOf";
          c.to = "oneOf";
          field_composition_out.push_back(std::move(c));
        }
      }
    }
    // Required: either the field's declared type is DescRequired<> (via
    // PALACE_SCHEMA_DESC_REQUIRED), or it is an `rfl::Literal<...>` (a
    // tagged-union discriminator declared via PALACE_SCHEMA_TAG — the
    // discriminator is always required for the union to parse).
    if constexpr (is_desc_required<Raw>::value || is_rfl_literal<Clean>::value)
    {
      required_out.emplace_back(std::string_view{name.str()});
    }
  }
};

// Visit one type-graph node. Dispatches to the right recursion (containers,
// optionals, TaggedUnion arms, or struct fields) and records a `$defs` entry
// when the type is itself a reflectable aggregate.
template <class T>
void visit_type(DefaultsAccum &acc)
{
  using U = std::remove_cvref_t<T>;

  if constexpr (is_vector<U>::value)
  {
    visit_type<strip_ann_t<typename U::value_type>>(acc);
  }
  else if constexpr (is_std_array<U>::value)
  {
    visit_type<strip_ann_t<typename U::value_type>>(acc);
  }
  else if constexpr (is_optional<U>::value)
  {
    visit_type<strip_ann_t<typename U::value_type>>(acc);
  }
  else if constexpr (is_rfl_tagged_union<U>::value)
  {
    // Peel the union at the type level so the arm pack becomes visible.
    // Each arm gets its own `$defs` entry; the union itself does not.
    // While here, harvest the discriminator name and per-arm tag
    // descriptions so the runtime pass can rewrite each arm's
    // discriminator property body into a `{const, description}` pair.
    []<rfl::internal::StringLiteral tag, class... Arms>(DefaultsAccum &a,
                                                        rfl::TaggedUnion<tag, Arms...> *)
    {
      DefaultsAccum::TaggedUnionDispatch d;
      d.tag_name = tag.str();
      auto push_desc = [&]<class Arm>()
      {
        if constexpr (has_tag_description<Arm>)
        {
          using Lit = rfl::internal::tag_t<tag, Arm>;
          d.arm_descs.emplace_back(Lit{}.name(), std::string(Arm::tag_description));
        }
      };
      (push_desc.template operator()<Arms>(), ...);
      a.tagged_unions.push_back(std::move(d));
      (visit_type<Arms>(a), ...);
    }(acc, static_cast<U *>(nullptr));
  }
  else if constexpr (is_rfl_variant<U>::value)
  {
    // rfl::Variant alternatives are emitted inline as `anyOf` — no `$defs`
    // entry for the variant itself, but each alternative may still reach
    // further reflectable structs that need to be walked.
    []<class... Alts>(DefaultsAccum &a, rfl::Variant<Alts...> *)
    { (visit_type<strip_ann_t<Alts>>(a), ...); }(acc, static_cast<U *>(nullptr));
  }
  else if constexpr (is_rfl_literal<U>::value)
  {
    // rfl::Literal is a primitive-ish at the schema level (string + enum).
  }
  else if constexpr (std::is_arithmetic_v<U> || std::is_enum_v<U> ||
                     std::is_same_v<U, std::string>)
  {
    // Plain leaves — nothing to add.
  }
  else if constexpr (std::is_class_v<U> && std::is_default_constructible_v<U>)
  {
    // Treat as a reflectable aggregate. `make_type_name` produces the
    // same `<ns>__<Name>` string reflect-cpp uses in $defs keys, so the
    // runtime pass can look it up directly.
    auto name = rfl::parsing::make_type_name<U>();
    if (acc.seen.count(name) != 0)
      return;
    acc.seen.insert(name);
    acc.entries.emplace(name, rfl::json::write(U{}));
    using Fields = typename rfl::named_tuple_t<U>::Fields;
    std::vector<std::string> required;
    field_annotations_walker<U, Fields>::apply(
        name, acc.enum_descs, acc.field_flags, required, acc.field_aliases,
        acc.variant_arm_aliases, acc.field_compositions);
    // Always record an entry (even if empty) so the `set_required` pass
    // knows to replace reflect-cpp's native required array — the empty
    // list signals "nothing required", which collapses to an entry with
    // no `required` key.
    acc.required_map.emplace(name, std::move(required));
    // Pick up a `schema_oneof_required<U>` specialization if the user
    // provided one. The primary template returns an empty vector and
    // produces no entry.
    {
      auto arms = schema_oneof_required<U>::value();
      if (!arms.empty())
      {
        OneOfRequired one;
        one.struct_name = name;
        one.arms = std::move(arms);
        acc.oneof_required.push_back(std::move(one));
      }
    }
    walk_struct_fields<U>(acc);
  }
  // Anything else (function pointers, unsupported types) silently falls
  // through — the runtime pass will just not inject defaults for them.
}

template <class T>
DefaultsAccum collect_defaults_accum()
{
  DefaultsAccum acc;
  visit_type<T>(acc);
  return acc;
}

}  // namespace palace::schema::utils::detail

namespace palace::schema::utils
{

template <class T>
std::string schema(SchemaOptions opts)
{
  // `required` rule: a field is required iff it is explicitly marked via
  // `PALACE_SCHEMA_DESC_REQUIRED` (or, for tagged-union discriminators,
  // `PALACE_SCHEMA_TAG`). Anything declared with plain `PALACE_SCHEMA_DESC`
  // — including non-optional scalars with in-class initializers — is
  // optional. This inverts reflect-cpp's default (every non-`std::optional`
  // field is required) so authors don't have to fight a heuristic that
  // can't distinguish `bool x = false;` from `bool x;`.
  auto s = rfl::json::to_schema<T>(rfl::json::pretty);
  // Flatten Validator-induced `allOf` blocks first so every later pass
  // — defaults, required, enum descriptions, custom keywords — sees
  // property bodies in their final shape. Keeps the post-emit pipeline
  // uniform instead of having each pass re-derive the collapse rule.
  s = detail::flatten_validator_allof(std::move(s));
  // One type-graph walk produces the inputs for every per-$defs-entry
  // pass: default JSON snapshots, required lists, enum descriptions, and
  // custom-keyword flags. Walk unconditionally — `emit_defaults=false`
  // still needs the enum + flag passes to match PR-716 schema shape.
  auto accum = detail::collect_defaults_accum<T>();
  if (opts.emit_defaults)
  {
    s = detail::inject_defaults(std::move(s), accum.entries);
  }
  s = detail::set_required(std::move(s), accum.required_map);
  s = detail::inject_additional_properties_false(std::move(s));
  if (!accum.enum_descs.empty())
  {
    s = detail::inject_enum_descriptions(std::move(s), accum.enum_descs);
  }
  if (!accum.field_flags.empty())
  {
    s = detail::inject_custom_keywords(std::move(s), accum.field_flags);
  }
  if (!accum.oneof_required.empty())
  {
    s = detail::inject_oneof_required(std::move(s), accum.oneof_required);
  }
  // Variant-arm aliasing must run before the whole-field alias pass: a
  // `Direction` field uses `rfl::Variant<Label, Vector3>`, where the
  // numeric arm shares the `$defs/Vector3` alias with plain
  // `std::array<double, 3>` field sites. Hoisting the arm first installs
  // the canonical body; the later field-alias pass for plain `Vector3`
  // fields then sees the entry already populated and only writes a
  // `$ref`.
  if (!accum.variant_arm_aliases.empty())
  {
    s = detail::inject_variant_arm_aliases(std::move(s), accum.variant_arm_aliases);
  }
  // Field-level composition rewrite (`anyOf` → `oneOf`) for variant
  // fields opted in via `schema_composition<Variant<...>>`. Runs after
  // the arm-alias pass — at this point each arm is already `{$ref:...}`,
  // so the only edit at the field site is the combiner keyword.
  if (!accum.field_compositions.empty())
  {
    s = detail::inject_field_composition(std::move(s), accum.field_compositions);
  }
  // Runs after every body-mutating pass so the hoisted `$defs` body
  // already reflects defaults/required/enum-description/custom-keyword
  // rewrites.
  if (!accum.field_aliases.empty())
  {
    s = detail::inject_field_aliases(std::move(s), accum.field_aliases);
  }
  // Run for every TaggedUnion encountered during the type-graph walk,
  // not only when the root T is itself a tagged union. Nested unions
  // (e.g. `Sample` under `Driven.Samples`) need the same body rewrite.
  for (const auto &d : accum.tagged_unions)
  {
    if (d.arm_descs.empty())
      continue;
    s = detail::inject_tag_descriptions(std::move(s), d.tag_name, d.arm_descs);
  }
  if constexpr (schema_composition<T>::value == Compose::OneOf)
  {
    s = detail::rewrite_root_composition(std::move(s), "anyOf", "oneOf");
  }
  if (!opts.version.empty())
  {
    s = detail::inject_root_id(std::move(s), std::move(opts.version));
  }
  // Runs last so earlier passes can key on reflect-cpp's native
  // fully-qualified `$defs` names (e.g. `palace__schema__ProblemData`)
  // without having to account for the rename.
  if (!opts.defs_prefix.empty())
  {
    s = detail::strip_defs_prefix(std::move(s), std::move(opts.defs_prefix));
  }
  return s;
}

}  // namespace palace::schema::utils

#endif  // PALACE_SCHEMA_UTILS_SCHEMA_IMPL_HPP
