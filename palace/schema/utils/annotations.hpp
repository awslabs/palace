// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_SCHEMA_UTILS_ANNOTATIONS_HPP
#define PALACE_SCHEMA_UTILS_ANNOTATIONS_HPP

#include <array>
#include <concepts>
#include <cstddef>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

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

// --- Per-enum-value descriptions -------------------------------------------
//
// reflect-cpp emits `enum class E` as an inline `{"type": "string", "enum":
// [...]}` body on every property of type E. Palace's docs want per-value
// descriptions in the PR-716 form: `{"oneOf": [{"const": "Foo",
// "description": "..."}, ...]}`. Specialize `enum_descriptions<E>` with a
// (value, description) array and the `inject_enum_descriptions` post-emit
// pass rewrites every matching property body.
//
// Empty specializations (the primary template) leave the flat enum alone.
template <class E>
struct enum_descriptions {
    static constexpr auto value =
        std::array<std::pair<std::string_view, std::string_view>, 0>{};
};

// --- Per-field custom-keyword flags ----------------------------------------
//
// JSON Schema reserves the `x-*` prefix for vendor extensions. Palace uses
// `x-palace-advanced` / `x-palace-deprecated` to flag fine-tuning knobs and
// legacy keys — the Julia doc generator renders a matching badge. Flags are
// attached at the field declaration via the `PALACE_SCHEMA_DESC_ADVANCED` /
// `PALACE_SCHEMA_DESC_DEPRECATED` macros (below), which expand to the
// annotated field plus a hidden-friend function keyed on a
// `FieldFlagTag<"FieldName">` tag type. The walker in schema_impl.hpp
// discovers the friend via ADL on the struct type and collects
// `(struct_name, field_name, flag)` triples for the post-emit
// `inject_custom_keywords` pass.
template <rfl::internal::StringLiteral FieldName>
struct FieldFlagTag {};

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

// Declare a described member that also carries an `x-palace-advanced`
// flag. Expands to the same `Desc<>` field `PALACE_SCHEMA_DESC` produces
// plus a hidden-friend overload selected by `FieldFlagTag<"FieldName">`;
// the walker in schema_impl.hpp finds that friend via ADL on the enclosing
// struct and records a `(struct, field, "advanced")` triple for
// `inject_custom_keywords`. The second parameter (`_palace_struct_tag_t*`)
// exists purely to bring the enclosing class into the ADL lookup so hidden
// friends are considered. Because the static-ness of the friend is
// invisible to reflect-cpp's aggregate-init field probe, this has no
// effect on (de)serialisation.
#define PALACE_SCHEMA_DESC_ADVANCED(FieldName, Description, .../*Type*/)       \
    friend constexpr ::std::string_view palace_schema_field_flag(              \
        ::palace::schema::utils::FieldFlagTag<#FieldName>,                     \
        const auto*) { return "advanced"; }                                    \
    ::palace::schema::utils::Desc<Description, __VA_ARGS__> FieldName

// As above, but stamps `x-palace-deprecated`. Intended for legacy keys the
// loader still accepts but docs should label as discouraged. The field
// declaration comes last so a call-site `= default;` initializer attaches
// to it, mirroring `PALACE_SCHEMA_DESC`.
#define PALACE_SCHEMA_DESC_DEPRECATED(FieldName, Description, .../*Type*/)     \
    friend constexpr ::std::string_view palace_schema_field_flag(              \
        ::palace::schema::utils::FieldFlagTag<#FieldName>,                     \
        const auto*) { return "deprecated"; }                                  \
    ::palace::schema::utils::Desc<Description, __VA_ARGS__> FieldName

// --- PALACE_SCHEMA_ENUM ----------------------------------------------------
//
// Declare an `enum class` together with its `enum_descriptions` spec in one
// place. Each entry is a `(Enumerator, "description")` tuple. Expansion
// produces the enum body and the trait specialization in the calling
// namespace; the specialization uses the fully-qualified name so it is
// legal as long as the call site is inside or enclosing `palace::schema`.
//
// Call site (inside `namespace palace::schema { ... }`):
//
//   PALACE_SCHEMA_ENUM(ProblemType,
//       (Eigenmode, "Perform an undamped or damped eigenfrequency analysis."),
//       (Driven,    "Perform a frequency-domain driven simulation."));
//
// Entries may supply an empty description (`""`) — the post-emit pass
// treats empty descriptions as "no description", so the generated schema
// still emits a bare `{"const": ...}` arm for those values. That's the
// knob to use for enums where only some enumerators have user-visible
// docs (see `KrylovSolver::MINRES` etc.).
//
// Implementation: small FOR_EACH-over-variadic, capped at 10 entries
// (Palace's largest enum has 8). Extend by adding PAL_SCHEMA_FE_N macros
// and padding PAL_SCHEMA_PICK_FE if a larger enum ever lands.

#define PAL_SCHEMA_EXPAND(x) x
#define PAL_SCHEMA_STR_IMPL(x) #x
#define PAL_SCHEMA_STR(x) PAL_SCHEMA_STR_IMPL(x)

// (Name, "desc") → Name
#define PAL_SCHEMA_TUPLE_NAME_IMPL(name, desc) name
#define PAL_SCHEMA_TUPLE_NAME(t) PAL_SCHEMA_EXPAND(PAL_SCHEMA_TUPLE_NAME_IMPL t)
// (Name, "desc") → "desc"
#define PAL_SCHEMA_TUPLE_DESC_IMPL(name, desc) desc
#define PAL_SCHEMA_TUPLE_DESC(t) PAL_SCHEMA_EXPAND(PAL_SCHEMA_TUPLE_DESC_IMPL t)

#define PAL_SCHEMA_FE_1(F, x)       F(x)
#define PAL_SCHEMA_FE_2(F, x, ...)  F(x) PAL_SCHEMA_FE_1(F, __VA_ARGS__)
#define PAL_SCHEMA_FE_3(F, x, ...)  F(x) PAL_SCHEMA_FE_2(F, __VA_ARGS__)
#define PAL_SCHEMA_FE_4(F, x, ...)  F(x) PAL_SCHEMA_FE_3(F, __VA_ARGS__)
#define PAL_SCHEMA_FE_5(F, x, ...)  F(x) PAL_SCHEMA_FE_4(F, __VA_ARGS__)
#define PAL_SCHEMA_FE_6(F, x, ...)  F(x) PAL_SCHEMA_FE_5(F, __VA_ARGS__)
#define PAL_SCHEMA_FE_7(F, x, ...)  F(x) PAL_SCHEMA_FE_6(F, __VA_ARGS__)
#define PAL_SCHEMA_FE_8(F, x, ...)  F(x) PAL_SCHEMA_FE_7(F, __VA_ARGS__)
#define PAL_SCHEMA_FE_9(F, x, ...)  F(x) PAL_SCHEMA_FE_8(F, __VA_ARGS__)
#define PAL_SCHEMA_FE_10(F, x, ...) F(x) PAL_SCHEMA_FE_9(F, __VA_ARGS__)

#define PAL_SCHEMA_PICK_FE(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, N, ...) N
#define PAL_SCHEMA_FOR_EACH(F, ...)                                            \
    PAL_SCHEMA_EXPAND(PAL_SCHEMA_PICK_FE(__VA_ARGS__,                          \
        PAL_SCHEMA_FE_10, PAL_SCHEMA_FE_9,  PAL_SCHEMA_FE_8,                   \
        PAL_SCHEMA_FE_7,  PAL_SCHEMA_FE_6,  PAL_SCHEMA_FE_5,                   \
        PAL_SCHEMA_FE_4,  PAL_SCHEMA_FE_3,  PAL_SCHEMA_FE_2,                   \
        PAL_SCHEMA_FE_1)(F, __VA_ARGS__))

#define PAL_SCHEMA_ENUM_EMIT_NAME(t) PAL_SCHEMA_TUPLE_NAME(t),
#define PAL_SCHEMA_ENUM_EMIT_PAIR(t)                                           \
    ::std::pair<::std::string_view, ::std::string_view>{                       \
        PAL_SCHEMA_STR(PAL_SCHEMA_TUPLE_NAME(t)),                              \
        PAL_SCHEMA_TUPLE_DESC(t)},

#define PALACE_SCHEMA_ENUM(EnumName, ...)                                      \
    enum class EnumName {                                                      \
        PAL_SCHEMA_FOR_EACH(PAL_SCHEMA_ENUM_EMIT_NAME, __VA_ARGS__)            \
    };                                                                         \
    template <>                                                                \
    struct ::palace::schema::utils::enum_descriptions<EnumName> {              \
        static constexpr auto value = ::std::array{                            \
            PAL_SCHEMA_FOR_EACH(PAL_SCHEMA_ENUM_EMIT_PAIR, __VA_ARGS__)        \
        };                                                                     \
    }

#endif  // PALACE_SCHEMA_UTILS_ANNOTATIONS_HPP
