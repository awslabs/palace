// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_SCHEMA_UTILS_ANNOTATIONS_HPP
#define PALACE_SCHEMA_UTILS_ANNOTATIONS_HPP

#include <concepts>
#include <cstddef>
#include <string>
#include <string_view>
#include <type_traits>

#include <rfl.hpp>

namespace palace::schema::utils {

// --- Schema emission options -------------------------------------------------
//
// `emit_defaults` controls whether `palace::schema::utils::schema<T>()`
// injects a `"default"` onto each property that has a C++ in-class
// initializer. Two modes:
//
//   true  (the default) — every property with a C++ initializer gets a
//         `"default"` in the schema. Useful for documentation and for
//         renderers that surface defaults in the generated docs.
//
//   false — no `"default"` keys are emitted. Matches Palace's legacy
//         hand-maintained schemas, which carried only the strict contract.
//
// The `required` list is unaffected by this flag — it always contains every
// non-`std::optional` field (matching reflect-cpp's native inference).
// Other passes (additionalProperties:false, TaggedUnion tag descriptions,
// schema_composition rewrite) also run regardless.
//
// `version`, if non-empty, is emitted as a `"version"` key on the root
// object of the schema.
struct SchemaOptions {
    bool emit_defaults = true;
    std::string version;
};

// Convenience alias for rfl::Description. Mirrors the StringLiteral<N> class-
// type non-type template parameter so callers can write
// `palace::schema::utils::Desc<"text", T>` directly.
template <rfl::internal::StringLiteral text, class T>
using Desc = rfl::Description<text, T>;

// --- TaggedUnion arm-discriminator descriptions -----------------------------
//
// rfl::TaggedUnion's arm types require a bare `rfl::Literal<"Tag">` member as
// the discriminator — any wrapper (Description, Validator, DefaultVal) around
// it breaks the dispatch. Descriptions travel through a
// `static constexpr std::string_view tag_description` side channel on the arm
// struct. Static members are not visible to reflect-cpp's aggregate-init
// field probe, so this is a clean channel that does not disturb
// (de)serialization.
template <class T>
concept has_tag_description = requires {
    { T::tag_description } -> std::convertible_to<std::string_view>;
};

// Trait matching `rfl::TaggedUnion<tag, Arms...>` so `schema<T>()` can
// dispatch the tag-description injection pass only for tagged unions.
template <class T>
struct is_tagged_union : std::false_type {};

template <rfl::internal::StringLiteral tag, class... Arms>
struct is_tagged_union<rfl::TaggedUnion<tag, Arms...>> : std::true_type {
    static constexpr auto discriminator = tag;
};

// --- Schema composition keyword control ------------------------------------
//
// JSON Schema offers several "combiner" keywords for schemas made of multiple
// sub-schemas: oneOf (exactly one matches), anyOf (at least one matches),
// allOf (all match simultaneously). reflect-cpp emits `anyOf` for every
// rfl::TaggedUnion — correct at runtime, but Palace's convention is `oneOf`
// for discriminated unions. Specialize `schema_composition<T>` for the
// user-side alias to rewrite the combiner keyword after emission.
enum class Compose { AnyOf, OneOf };

template <class T>
struct schema_composition {
    static constexpr auto value = Compose::AnyOf;
};

// --- Validator sugar -------------------------------------------------------
//
// Shortcuts for the single- and two-bound numeric-range validators. Interval
// notation (Closed / Open / LeftOpen / RightOpen) keeps the call site legible.

// value >= lo
template <class T, auto lo>
using Min = rfl::Validator<T, rfl::Minimum<lo>>;

// value <= hi
template <class T, auto hi>
using Max = rfl::Validator<T, rfl::Maximum<hi>>;

// value > lo  (exclusive minimum)
template <class T, auto lo>
using XMin = rfl::Validator<T, rfl::ExclusiveMinimum<lo>>;

// value < hi  (exclusive maximum)
template <class T, auto hi>
using XMax = rfl::Validator<T, rfl::ExclusiveMaximum<hi>>;

// [lo, hi] — inclusive on both ends
template <class T, auto lo, auto hi>
using Closed = rfl::Validator<T, rfl::Minimum<lo>, rfl::Maximum<hi>>;

// (lo, hi) — exclusive on both ends
template <class T, auto lo, auto hi>
using Open = rfl::Validator<T, rfl::ExclusiveMinimum<lo>, rfl::ExclusiveMaximum<hi>>;

// (lo, hi] — exclusive minimum, inclusive maximum
template <class T, auto lo, auto hi>
using LeftOpen = rfl::Validator<T, rfl::ExclusiveMinimum<lo>, rfl::Maximum<hi>>;

// [lo, hi) — inclusive minimum, exclusive maximum
template <class T, auto lo, auto hi>
using RightOpen = rfl::Validator<T, rfl::Minimum<lo>, rfl::ExclusiveMaximum<hi>>;

}  // namespace palace::schema::utils

// Attach a description to a TaggedUnion arm's discriminator literal. Call-site
// form:
//   PALACE_SCHEMA_TAG(Type, "Point", "Explicit frequency list.");
// The trailing `;` terminates the rfl::Literal member declaration; the inner
// `;` separates the two member declarations.
#define PALACE_SCHEMA_TAG(FieldName, Value, Desc)                              \
    static constexpr ::std::string_view tag_description = (Desc);              \
    ::rfl::Literal<Value> FieldName

// Declare a described member. Call-site form:
//   PALACE_SCHEMA_DESC(Verbose, "Level of printing",
//                      palace::schema::utils::Min<int, 0>) = 1;
// Type is the trailing variadic so commas inside template arguments pass
// through the preprocessor without extra parentheses.
#define PALACE_SCHEMA_DESC(FieldName, Description, .../*Type*/)                \
    ::palace::schema::utils::Desc<Description, __VA_ARGS__> FieldName

#endif  // PALACE_SCHEMA_UTILS_ANNOTATIONS_HPP
