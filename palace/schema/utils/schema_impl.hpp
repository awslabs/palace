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
#include <rfl/internal/StringLiteral.hpp>
#include <rfl/internal/tag_t.hpp>
#include <rfl/json.hpp>
#include <rfl/named_tuple_t.hpp>
#include <rfl/parsing/make_type_name.hpp>

#include "annotations.hpp"

namespace palace::schema::utils::detail {

// Defined in src/schema.cpp. Takes the raw reflect-cpp schema string plus a
// map of `$defs`-key-name -> serialized default JSON for that entry's type.
// For each `$defs` entry whose key appears in the map, walks the entry's
// `properties` in lock-step with the parsed defaults and injects `"default"`
// at each matching leaf. Entries not in the map (or properties we can't match)
// are left alone except for the optional-field `"default": null` fallback.
std::string inject_defaults(
    std::string schema_json,
    const std::unordered_map<std::string, std::string>& defaults_map);

// Defined in src/schema.cpp. Walks the whole schema tree and stamps
// `"additionalProperties": false` onto every object that declares
// `properties` and doesn't already set the key. reflect-cpp omits this
// keyword entirely, so external validators (Python jsonschema, etc.) would
// otherwise accept unknown keys. Independent of the default-injection pass so
// it can run in schemas emitted without defaults too.
std::string inject_additional_properties_false(std::string schema_json);

// Defined in src/schema.cpp. Walks `$defs` and, for each entry whose name
// appears in `prune_map`, removes the listed field names from the entry's
// `required` array. This implements the "required ⟺ not optional AND no
// meaningful default" rule — fields with a C++ default that overrides the
// field type's zero-init baseline get dropped from `required`, matching the
// loader's silent-default behavior for those fields.
std::string prune_required(
    std::string schema_json,
    const std::unordered_map<std::string, std::vector<std::string>>& prune_map);

// Defined in src/schema.cpp. Walks a palace::schema::utils-emitted schema alongside an input
// JSON document, returning dotted paths of every key listed in a `required`
// array that is absent from its corresponding object. Handles $ref navigation,
// nullable `anyOf: [T, null]`, discriminated `oneOf` via discriminator lookup,
// and array `items` recursion. Other keywords (type, minimum, enum, etc.) are
// deferred to rfl::json::read which enforces them at parse time.
//
// Returns an empty vector when the input has every required key at every
// level. This is the bridge that lets `palace::schema::utils::load<T>` reject JSON the loader's
// `rfl::DefaultIfMissing` would otherwise silently fill.
std::vector<std::string> check_required_keys(const std::string& schema_json,
                                              const std::string& input_json);

// Defined in src/schema.cpp. Walks `$defs` and, for each def whose
// `properties[tag_name]` single-element `enum` matches one of `arm_descs`'
// discriminator values, writes a `"description"` onto that property body.
// Existing descriptions are preserved (the injection only fills empty slots).
std::string inject_tag_descriptions(std::string schema_json,
                                    std::string tag_name,
                                    std::vector<std::pair<std::string, std::string>> arm_descs);

// Defined in src/schema.cpp. Renames the root's combiner keyword (e.g.
// `anyOf` to `oneOf`) when the user opts in via `palace::schema::utils::schema_composition<T>`.
// No-op if the root does not contain the `from` key.
std::string rewrite_root_composition(std::string schema_json,
                                     std::string from,
                                     std::string to);

// Defined in src/schema.cpp. Writes `"version": <value>` onto the root of
// the schema. Overwrites any existing value so callers can re-run the pass
// safely. Skipped by the caller when `version` is empty.
std::string inject_root_version(std::string schema_json, std::string version);

// Defined in src/schema.cpp. For each `(struct_name, field_name, pairs)`
// triple, locates `$defs[struct_name]/properties[field_name]`. If the body
// matches reflect-cpp's inline-enum shape (`{"type": "string", "enum":
// [...]}`), it is rewritten to PR-716's `{"oneOf": [{"const": v,
// "description": d}, ...]}` form using `pairs` as the (value, description)
// map. Outer `description` / `default` / other keys are preserved; `type`
// and `enum` are removed since `oneOf` subsumes them.
struct EnumDescEntry {
    std::string struct_name;
    std::string field_name;
    std::vector<std::pair<std::string, std::string>> pairs;
};
std::string inject_enum_descriptions(std::string schema_json,
                                     const std::vector<EnumDescEntry>& entries);

// Defined in src/schema.cpp. For each `(struct_name, field_name, flag)`
// triple, stamps `"x-palace-<flag>": true` onto the matching property body
// under `$defs[struct_name]/properties[field_name]`. Flag strings are
// `"advanced"` / `"deprecated"` in practice; any other string is accepted
// and written through verbatim.
struct FieldFlag {
    std::string struct_name;
    std::string field_name;
    std::string flag;
};
std::string inject_custom_keywords(std::string schema_json,
                                   const std::vector<FieldFlag>& flags);

// Compile-time dispatch helper: pattern-matches a `rfl::TaggedUnion<tag,
// Arms...>` via a pointer parameter so we can fold over the arm pack and pull
// each arm's discriminator literal out of `rfl::internal::tag_t`. Only arms
// that provide a `tag_description` static member (via CJS_TAG) contribute a
// pair; if no arm has one, we skip the schema round-trip entirely.
template <rfl::internal::StringLiteral tag, class... Arms>
std::string dispatch_tag_descriptions(std::string s,
                                      rfl::TaggedUnion<tag, Arms...>*) {
    std::vector<std::pair<std::string, std::string>> pairs;
    auto push = [&]<class Arm>() {
        if constexpr (has_tag_description<Arm>) {
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
    if (pairs.empty()) return s;
    return inject_tag_descriptions(std::move(s), tag.str(), std::move(pairs));
}

template <class T>
std::string inject_tag_descriptions_for(std::string s) {
    return dispatch_tag_descriptions(std::move(s), static_cast<T*>(nullptr));
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
// ProblemData::Verbose), so the recursion has to continue until the type
// stops matching any wrapper specialization.
template <class T>
struct strip_ann { using type = T; };

template <rfl::internal::StringLiteral name, class T>
struct strip_ann<rfl::Description<name, T>> {
    using type = typename strip_ann<T>::type;
};

template <rfl::internal::StringLiteral name, class T>
struct strip_ann<rfl::Rename<name, T>> {
    using type = typename strip_ann<T>::type;
};

template <class T, class V, class... Vs>
struct strip_ann<rfl::Validator<T, V, Vs...>> {
    using type = typename strip_ann<T>::type;
};

template <class T>
struct strip_ann<rfl::DefaultVal<T>> {
    using type = typename strip_ann<T>::type;
};

template <class T>
using strip_ann_t = typename strip_ann<std::remove_cvref_t<T>>::type;

// Shape traits for container/union types we unwrap without treating as a
// `$defs` entry in their own right. These match the shapes reflect-cpp emits
// inline (arrays, anyOf-wrapped optionals, TaggedUnion `anyOf`/`oneOf`).
template <class>       struct is_vector                 : std::false_type {};
template <class U, class A> struct is_vector<std::vector<U, A>> : std::true_type {};

template <class>                struct is_std_array     : std::false_type {};
template <class U, std::size_t N>
struct is_std_array<std::array<U, N>>                   : std::true_type {};

template <class>       struct is_optional               : std::false_type {};
template <class U>     struct is_optional<std::optional<U>> : std::true_type {};

template <class>       struct is_rfl_tagged_union       : std::false_type {};
template <rfl::internal::StringLiteral tag, class... Arms>
struct is_rfl_tagged_union<rfl::TaggedUnion<tag, Arms...>> : std::true_type {};

template <class>       struct is_rfl_literal            : std::false_type {};
template <rfl::internal::StringLiteral... names>
struct is_rfl_literal<rfl::Literal<names...>>           : std::true_type {};

// "Scalar" for the prune rule: a leaf type whose JSON representation is a
// single primitive token — arithmetic, enum, string, or rfl::Literal (a
// one-of-these-strings tag). Scalars are the only kind of field that stays
// in `required` under the prune rule; aggregates and containers are always
// pruneable because a parent-level `= {}` defers to the child's own
// defaults, which is semantically optional.
template <class T>
inline constexpr bool is_scalar_v =
    std::is_arithmetic_v<T> || std::is_enum_v<T> ||
    std::is_same_v<T, std::string> || is_rfl_literal<T>::value;

// Accumulator for the type-graph walk. Plain struct so the recursion can pass
// it by reference without per-call template-parameter noise.
//
// `prune_map` records the names of fields whose C++ default differs from the
// field type's zero-init baseline. Those fields have a "meaningful" default
// and get dropped from the schema's `required` list (and become defaultable
// in the loader). Keyed by `$defs` name so the runtime pass can look up per
// aggregate type.
struct DefaultsAccum {
    std::unordered_map<std::string, std::string> entries;
    std::unordered_map<std::string, std::vector<std::string>> prune_map;
    // Per-enum-value descriptions and per-field custom-keyword flags are
    // harvested while walking the same type graph so consumers can run a
    // single pass. `enum_descs` is flattened into one entry per
    // (struct, enum-field) pair; `field_flags` is one entry per flagged
    // field.
    std::vector<EnumDescEntry> enum_descs;
    std::vector<FieldFlag> field_flags;
    std::unordered_set<std::string> seen;
};

// Forward declaration — mutually recursive with walk_struct_fields.
template <class T>
void visit_type(DefaultsAccum& acc);

// Walk the field-value types of a reflect-cpp-reflectable struct. `Values` is
// `rfl::Tuple<T1, T2, ...>` where each Ti is the struct's declared member
// type (possibly still wrapped in Description/Validator/etc.). We strip the
// wrappers and recurse. This helper exists as a struct so it can be partially
// specialized on rfl::Tuple<...>; function templates cannot do that.
template <class Tup>
struct walk_tuple;

template <class... Ts>
struct walk_tuple<rfl::Tuple<Ts...>> {
    static void apply(DefaultsAccum& acc) {
        (visit_type<strip_ann_t<Ts>>(acc), ...);
    }
};

template <class T>
void walk_struct_fields(DefaultsAccum& acc) {
    using Values = typename rfl::named_tuple_t<T>::Values;
    walk_tuple<Values>::apply(acc);
}

// Per-field prune detector. For each `rfl::Field<name, Type>` in U's named
// tuple, strip annotations to the clean value type, then compare the actual
// value inside U{} (as serialised via rfl::json::write(U{})) against the
// baseline `rfl::json::write(Clean{})`. Mismatch means "user overrode the
// type's zero-init default" → the field is optional with that default.
//
// std::optional fields are skipped (they're already optional structurally —
// reflect-cpp emits them as anyOf[T, null] and never lists them in required).
// Non-reachable or unserialisable fields are left alone (the prune list stays
// conservative — fewer prunes = more-required schema, which is the safer
// direction).
template <class T, class FieldsTup>
struct prune_fields_walker;

template <class T, class... Fs>
struct prune_fields_walker<T, rfl::Tuple<Fs...>> {
    static void apply(const rfl::Generic::Object& u_obj,
                      std::vector<std::string>& out) {
        (check_one<Fs>(u_obj, out), ...);
    }

    template <class F>
    static void check_one(const rfl::Generic::Object& u_obj,
                          std::vector<std::string>& out) {
        check_one_impl(u_obj, out, static_cast<F*>(nullptr));
    }

    template <rfl::internal::StringLiteral name, class Type>
    static void check_one_impl(const rfl::Generic::Object& u_obj,
                               std::vector<std::string>& out,
                               rfl::Field<name, Type>*) {
        using Clean = strip_ann_t<std::remove_cvref_t<Type>>;
        if constexpr (is_optional<Clean>::value) {
            return;  // optional fields are never in `required` to begin with
        } else if constexpr (!is_scalar_v<Clean>) {
            // Aggregates / vectors / arrays / tagged unions: always
            // pruneable. A parent writing `Field = {}` (or any other form)
            // defers to the child's own defaults, which is the idiomatic C++
            // way to say "this sub-section is optional, use its defaults if
            // absent". Schema-side this translates to an entry with a
            // documented `"default"` and no seat in `required`.
            out.push_back(std::string{std::string_view{name.str()}});
        } else {
            // Scalars: stay required iff the parent didn't override the
            // type's zero-init (e.g. `int Index;` or `int Index = 0;` stays
            // required; `int Verbose = 1;` gets pruned with default 1).
            std::string field_name{std::string_view{name.str()}};
            if (u_obj.count(field_name) == 0) return;
            auto actual_str = rfl::json::write(u_obj.at(field_name));
            auto baseline_str = rfl::json::write(Clean{});
            if (actual_str != baseline_str) {
                out.push_back(std::move(field_name));
            }
        }
    }
};

// Per-field walker that collects two kinds of annotation in one pass:
//
//   * Enum-description entries: for every `rfl::Field<name, Type>` whose
//     cleaned value type is an enum E with non-empty
//     `enum_descriptions<E>`, emit an EnumDescEntry keyed on U's `$defs`
//     name and the field name. Enums without a specialization stay flat.
//
//   * Custom-keyword flags: for every field declared via
//     `PALACE_SCHEMA_DESC_ADVANCED` / `PALACE_SCHEMA_DESC_DEPRECATED`, the
//     macro injects a hidden-friend overload of `palace_schema_field_flag`
//     keyed on `FieldFlagTag<"FieldName">`. We discover it via ADL by
//     passing a `const U*` as the second argument (hidden friends are
//     found through their enclosing class's associated namespace/type
//     set). The `requires` check falls through silently for fields
//     without a flag.
template <class U, class FieldsTup>
struct field_annotations_walker;

template <class U, class... Fs>
struct field_annotations_walker<U, rfl::Tuple<Fs...>> {
    static void apply(const std::string& struct_name,
                      std::vector<EnumDescEntry>& enum_out,
                      std::vector<FieldFlag>& flag_out) {
        (check_one<Fs>(struct_name, enum_out, flag_out), ...);
    }

    template <class F>
    static void check_one(const std::string& struct_name,
                          std::vector<EnumDescEntry>& enum_out,
                          std::vector<FieldFlag>& flag_out) {
        check_one_impl(struct_name, enum_out, flag_out,
                       static_cast<F*>(nullptr));
    }

    template <rfl::internal::StringLiteral name, class Type>
    static void check_one_impl(const std::string& struct_name,
                               std::vector<EnumDescEntry>& enum_out,
                               std::vector<FieldFlag>& flag_out,
                               rfl::Field<name, Type>*) {
        using Clean = strip_ann_t<std::remove_cvref_t<Type>>;
        // The enum-description rewrite handles direct enum properties,
        // which covers every PR-716 per-value-description site. Enums
        // inside optionals, vectors, or tagged unions fall back to their
        // flat shape.
        if constexpr (std::is_enum_v<Clean>) {
            if constexpr (enum_descriptions<Clean>::value.size() != 0) {
                EnumDescEntry entry;
                entry.struct_name = struct_name;
                entry.field_name.assign(std::string_view{name.str()});
                for (auto [v, d] : enum_descriptions<Clean>::value) {
                    entry.pairs.emplace_back(std::string(v), std::string(d));
                }
                enum_out.push_back(std::move(entry));
            }
        }
        // ADL-probe for a PALACE_SCHEMA_DESC_ADVANCED /
        // PALACE_SCHEMA_DESC_DEPRECATED hidden friend on U. The friend's
        // first parameter type `FieldFlagTag<name>` discriminates per
        // field; the second `const auto*` accepts the `const U*` we pass
        // so U's hidden friends are looked up via ADL.
        if constexpr (requires {
            palace_schema_field_flag(FieldFlagTag<name>{},
                                     static_cast<const U*>(nullptr));
        }) {
            constexpr auto flag_sv = palace_schema_field_flag(
                FieldFlagTag<name>{}, static_cast<const U*>(nullptr));
            FieldFlag f;
            f.struct_name = struct_name;
            f.field_name.assign(std::string_view{name.str()});
            f.flag.assign(flag_sv);
            flag_out.push_back(std::move(f));
        }
    }
};

// Compute the prune list for one aggregate U. Serialises U{} once, then
// walks its compile-time field list against the parsed object.
template <class U>
std::vector<std::string> collect_prune_list_for() {
    auto u_json = rfl::json::write(U{});
    auto u_parsed = rfl::json::read<rfl::Generic>(u_json);
    if (!u_parsed) return {};
    if (!std::holds_alternative<rfl::Generic::Object>(u_parsed->value())) {
        return {};
    }
    const auto& u_obj = std::get<rfl::Generic::Object>(u_parsed->value());
    std::vector<std::string> out;
    using Fields = typename rfl::named_tuple_t<U>::Fields;
    prune_fields_walker<U, Fields>::apply(u_obj, out);
    return out;
}

// Visit one type-graph node. Dispatches to the right recursion (containers,
// optionals, TaggedUnion arms, or struct fields) and records a `$defs` entry
// when the type is itself a reflectable aggregate.
template <class T>
void visit_type(DefaultsAccum& acc) {
    using U = std::remove_cvref_t<T>;

    if constexpr (is_vector<U>::value) {
        visit_type<strip_ann_t<typename U::value_type>>(acc);
    } else if constexpr (is_std_array<U>::value) {
        visit_type<strip_ann_t<typename U::value_type>>(acc);
    } else if constexpr (is_optional<U>::value) {
        visit_type<strip_ann_t<typename U::value_type>>(acc);
    } else if constexpr (is_rfl_tagged_union<U>::value) {
        // Peel the union at the type level so the arm pack becomes visible.
        // Each arm gets its own `$defs` entry; the union itself does not.
        []<rfl::internal::StringLiteral tag, class... Arms>(
            DefaultsAccum& a, rfl::TaggedUnion<tag, Arms...>*) {
            (visit_type<Arms>(a), ...);
        }(acc, static_cast<U*>(nullptr));
    } else if constexpr (is_rfl_literal<U>::value) {
        // rfl::Literal is a primitive-ish at the schema level (string + enum).
    } else if constexpr (std::is_arithmetic_v<U> || std::is_enum_v<U> ||
                         std::is_same_v<U, std::string>) {
        // Plain leaves — nothing to add.
    } else if constexpr (std::is_class_v<U> && std::is_default_constructible_v<U>) {
        // Treat as a reflectable aggregate. `make_type_name` produces the
        // same `<ns>__<Name>` string reflect-cpp uses in $defs keys, so the
        // runtime pass can look it up directly.
        auto name = rfl::parsing::make_type_name<U>();
        if (acc.seen.count(name) != 0) return;
        acc.seen.insert(name);
        acc.entries.emplace(name, rfl::json::write(U{}));
        auto prune_list = collect_prune_list_for<U>();
        if (!prune_list.empty()) {
            acc.prune_map.emplace(name, std::move(prune_list));
        }
        using Fields = typename rfl::named_tuple_t<U>::Fields;
        field_annotations_walker<U, Fields>::apply(
            name, acc.enum_descs, acc.field_flags);
        walk_struct_fields<U>(acc);
    }
    // Anything else (function pointers, unsupported types) silently falls
    // through — the runtime pass will just not inject defaults for them.
}

template <class T>
DefaultsAccum collect_defaults_accum() {
    DefaultsAccum acc;
    visit_type<T>(acc);
    return acc;
}

}  // namespace palace::schema::utils::detail

namespace palace::schema::utils {

template <class T>
std::string schema(SchemaOptions opts) {
    // `required` rule: a field is required iff ALL of
    //   (1) its type is not `std::optional<T>`,
    //   (2) its clean type is a scalar primitive (arithmetic / enum /
    //       std::string / rfl::Literal), and
    //   (3) `T{}.field` equals `FieldType{}` (the field did not override its
    //       type's zero-init).
    // Anything that's not a scalar — nested aggregates, vectors, arrays,
    // tagged unions — is always pruneable: a parent-level `Field = {}`
    // semantically defers to the child's own defaults, which matches "this
    // sub-section is optional" in the JSON-config sense. The child type is
    // responsible for marking its own scalar fields as required.
    // reflect-cpp natively emits every non-optional field in `required`;
    // the `prune_required` pass below drops everything that shouldn't be
    // there per the rule. This keeps the schema in sync with the loader.
    auto s = rfl::json::to_schema<T>(rfl::json::pretty);
    // One type-graph walk produces the inputs for every per-$defs-entry
    // pass: default JSON snapshots, prune lists, enum descriptions, and
    // custom-keyword flags. Walk unconditionally — `emit_defaults=false`
    // still needs the enum + flag passes to match PR-716 schema shape.
    auto accum = detail::collect_defaults_accum<T>();
    if (opts.emit_defaults) {
        s = detail::inject_defaults(std::move(s), accum.entries);
        s = detail::prune_required(std::move(s), accum.prune_map);
    }
    s = detail::inject_additional_properties_false(std::move(s));
    if (!accum.enum_descs.empty()) {
        s = detail::inject_enum_descriptions(std::move(s), accum.enum_descs);
    }
    if (!accum.field_flags.empty()) {
        s = detail::inject_custom_keywords(std::move(s), accum.field_flags);
    }
    if constexpr (is_tagged_union<T>::value) {
        s = detail::inject_tag_descriptions_for<T>(std::move(s));
    }
    if constexpr (schema_composition<T>::value == Compose::OneOf) {
        s = detail::rewrite_root_composition(std::move(s), "anyOf", "oneOf");
    }
    if (!opts.version.empty()) {
        s = detail::inject_root_version(std::move(s), std::move(opts.version));
    }
    return s;
}

}  // namespace palace::schema::utils

#endif  // PALACE_SCHEMA_UTILS_SCHEMA_IMPL_HPP
