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
#include <vector>

#include <rfl.hpp>
#include <rfl/Description.hpp>
#include <rfl/internal/StringLiteral.hpp>
#include <rfl/internal/is_description.hpp>

namespace palace::schema::utils
{

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
// `required` is not controlled by this flag. It is driven entirely by
// explicit markers (`PALACE_SCHEMA_DESC_REQUIRED`, `PALACE_SCHEMA_TAG`);
// anything else is optional and the loader substitutes the C++ default.
//
// `version`, if non-empty, is emitted as a `"version"` key on the root
// object of the schema.
//
// `defs_prefix`, if non-empty, is stripped from every `$defs` entry name
// and every `$ref` string in the final schema. reflect-cpp derives each
// `$defs` key from `rfl::parsing::make_type_name<T>()`, which renders
// `palace::schema::Foo` as `palace__schema__Foo`. Passing
// `"palace__schema__"` here collapses those to just `Foo`, matching the
// hand-written-schema conventions in `scripts/schema/`.
struct SchemaOptions
{
  bool emit_defaults = true;
  std::string version;
  std::string defs_prefix;
};

// Convenience alias for rfl::Description. Mirrors the StringLiteral<N> class-
// type non-type template parameter so callers can write
// `palace::schema::utils::Desc<"text", T>` directly.
template <rfl::internal::StringLiteral text, class T>
using Desc = rfl::Description<text, T>;

// --- Annotated description subtypes ----------------------------------------
//
// Three thin subclasses of `rfl::Description<text, T>` carry Palace-specific
// per-field metadata entirely at the type level. Because they inherit from
// `Description`, they preserve its `ReflectionType`, `Content`, and
// `reflection()` members — reflect-cpp's JSON reader/writer route through
// them transparently (via `has_reflection_type_v`), and a companion
// `is_description` specialization below makes the schema emitter also treat
// them as plain descriptions. The walker in `schema_impl.hpp` detects the
// subtype on each field via a type trait and emits:
//
//   DescRequired<text, T>   → property listed in its parent's `required`
//   DescAdvanced<text, T>   → `"x-palace-advanced": true` on the property
//   DescDeprecated<text, T> → `"x-palace-deprecated": true` on the property
//
// No ADL tricks, no class-scope boilerplate: the marker is a property of
// the field's declared type, which is all the walker needs.
template <rfl::internal::StringLiteral text, class T>
struct DescRequired : rfl::Description<text, T>
{
  using rfl::Description<text, T>::Description;
  using rfl::Description<text, T>::operator=;
};

template <rfl::internal::StringLiteral text, class T>
struct DescAdvanced : rfl::Description<text, T>
{
  using rfl::Description<text, T>::Description;
  using rfl::Description<text, T>::operator=;
};

template <rfl::internal::StringLiteral text, class T>
struct DescDeprecated : rfl::Description<text, T>
{
  using rfl::Description<text, T>::Description;
  using rfl::Description<text, T>::operator=;
};

// --- Type traits used by the walker ----------------------------------------

template <class>
struct is_desc_required : std::false_type
{
};
template <rfl::internal::StringLiteral text, class T>
struct is_desc_required<DescRequired<text, T>> : std::true_type
{
};

template <class>
struct is_desc_advanced : std::false_type
{
};
template <rfl::internal::StringLiteral text, class T>
struct is_desc_advanced<DescAdvanced<text, T>> : std::true_type
{
};

template <class>
struct is_desc_deprecated : std::false_type
{
};
template <rfl::internal::StringLiteral text, class T>
struct is_desc_deprecated<DescDeprecated<text, T>> : std::true_type
{
};

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
struct is_tagged_union : std::false_type
{
};

template <rfl::internal::StringLiteral tag, class... Arms>
struct is_tagged_union<rfl::TaggedUnion<tag, Arms...>> : std::true_type
{
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
enum class Compose
{
  AnyOf,
  OneOf
};

template <class T>
struct schema_composition
{
  static constexpr auto value = Compose::AnyOf;
};

// --- Cross-field `oneOf` required alternatives ----------------------------
//
// Some structs satisfy one of several mutually-exclusive "shapes" — e.g.
// a LumpedPort item must provide either `Attributes` (single-element
// port) or `Elements` (multi-element port) but not both. That is not
// representable by a simple `required` list; JSON Schema expresses it as
// `oneOf: [{required:[...A]}, {required:[...B]}]`. Specialize
// `schema_oneof_required<T>` for any struct that needs this shape; the
// post-emit `inject_oneof_required` pass splices the resulting `oneOf`
// into `$defs[type_name(T)]`, alongside the ordinary `required` array.
//
// Returning the value as a `std::vector<std::vector<std::string>>` from a
// static member function keeps the API flexible (arms of different
// lengths) while avoiding constexpr-array-of-arrays gymnastics.
template <class T>
struct schema_oneof_required
{
  static std::vector<std::vector<std::string>> value() { return {}; }
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
struct enum_descriptions
{
  static constexpr auto value =
      std::array<std::pair<std::string_view, std::string_view>, 0>{};
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

// Teach reflect-cpp's schema emitter to treat our Description subtypes as
// plain descriptions. `is_description_v` drives the `make_description`
// branch in `rfl::parsing::Parser_default::to_schema`, which uses
// `U::Content` and `U::Type` — both inherited from `rfl::Description` by
// our subtypes — to render the property body. Without these
// specializations, the emitter would fall through to `has_reflection_type_v`
// and register each subtype as a separate `$defs` entry.
namespace rfl::internal
{
template <StringLiteral text, class T>
class is_description<::palace::schema::utils::DescRequired<text, T>> : public std::true_type
{
};
template <StringLiteral text, class T>
class is_description<::palace::schema::utils::DescAdvanced<text, T>> : public std::true_type
{
};
template <StringLiteral text, class T>
class is_description<::palace::schema::utils::DescDeprecated<text, T>>
  : public std::true_type
{
};
}  // namespace rfl::internal

// Attach a description to a TaggedUnion arm's discriminator literal.
// Call-site form:
//   PALACE_SCHEMA_TAG(Type, "Point", "Explicit frequency list.");
// The description travels through a `tag_description` side channel (see
// the `has_tag_description` concept above). The walker marks any field
// of type `rfl::Literal<...>` as required, so the discriminator is
// automatically part of the arm's `required` array — no extra machinery
// here.
#define PALACE_SCHEMA_TAG(FieldName, Value, Desc)               \
  static constexpr ::std::string_view tag_description = (Desc); \
  ::rfl::Literal<Value> FieldName

// Declare a described member. Call-site form:
//   PALACE_SCHEMA_DESC(Verbose, "Level of printing",
//                      palace::schema::utils::Min<int, 0>) = 1;
// Type is the trailing variadic so commas inside template arguments pass
// through the preprocessor without extra parentheses.
#define PALACE_SCHEMA_DESC(FieldName, Description, ... /*Type*/) \
  ::palace::schema::utils::Desc<Description, __VA_ARGS__> FieldName

// Declare a described member carrying the `x-palace-advanced` flag. The
// walker detects the `DescAdvanced<>` wrapper type at field-probe time and
// stamps `"x-palace-advanced": true` on the emitted property body.
#define PALACE_SCHEMA_DESC_ADVANCED(FieldName, Description, ... /*Type*/) \
  ::palace::schema::utils::DescAdvanced<Description, __VA_ARGS__> FieldName

// Declare a described member carrying the `x-palace-deprecated` flag.
// Intended for legacy keys the loader still accepts but docs should label
// as discouraged.
#define PALACE_SCHEMA_DESC_DEPRECATED(FieldName, Description, ... /*Type*/) \
  ::palace::schema::utils::DescDeprecated<Description, __VA_ARGS__> FieldName

// Declare a described member and mark it as required in the emitted
// schema. The walker replaces reflect-cpp's native `required` array with
// only the fields whose declared type is `DescRequired<>` — anything
// declared via plain `PALACE_SCHEMA_DESC` is optional and the loader
// substitutes the C++ initializer when the key is absent.
#define PALACE_SCHEMA_DESC_REQUIRED(FieldName, Description, ... /*Type*/) \
  ::palace::schema::utils::DescRequired<Description, __VA_ARGS__> FieldName

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

#define PAL_SCHEMA_FE_1(F, x) F(x)
#define PAL_SCHEMA_FE_2(F, x, ...) F(x) PAL_SCHEMA_FE_1(F, __VA_ARGS__)
#define PAL_SCHEMA_FE_3(F, x, ...) F(x) PAL_SCHEMA_FE_2(F, __VA_ARGS__)
#define PAL_SCHEMA_FE_4(F, x, ...) F(x) PAL_SCHEMA_FE_3(F, __VA_ARGS__)
#define PAL_SCHEMA_FE_5(F, x, ...) F(x) PAL_SCHEMA_FE_4(F, __VA_ARGS__)
#define PAL_SCHEMA_FE_6(F, x, ...) F(x) PAL_SCHEMA_FE_5(F, __VA_ARGS__)
#define PAL_SCHEMA_FE_7(F, x, ...) F(x) PAL_SCHEMA_FE_6(F, __VA_ARGS__)
#define PAL_SCHEMA_FE_8(F, x, ...) F(x) PAL_SCHEMA_FE_7(F, __VA_ARGS__)
#define PAL_SCHEMA_FE_9(F, x, ...) F(x) PAL_SCHEMA_FE_8(F, __VA_ARGS__)
#define PAL_SCHEMA_FE_10(F, x, ...) F(x) PAL_SCHEMA_FE_9(F, __VA_ARGS__)

#define PAL_SCHEMA_PICK_FE(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, N, ...) N
#define PAL_SCHEMA_FOR_EACH(F, ...)                                                       \
  PAL_SCHEMA_EXPAND(PAL_SCHEMA_PICK_FE(__VA_ARGS__, PAL_SCHEMA_FE_10, PAL_SCHEMA_FE_9,    \
                                       PAL_SCHEMA_FE_8, PAL_SCHEMA_FE_7, PAL_SCHEMA_FE_6, \
                                       PAL_SCHEMA_FE_5, PAL_SCHEMA_FE_4, PAL_SCHEMA_FE_3, \
                                       PAL_SCHEMA_FE_2, PAL_SCHEMA_FE_1)(F, __VA_ARGS__))

#define PAL_SCHEMA_ENUM_EMIT_NAME(t) PAL_SCHEMA_TUPLE_NAME(t),
#define PAL_SCHEMA_ENUM_EMIT_PAIR(t)                   \
  ::std::pair<::std::string_view, ::std::string_view>{ \
      PAL_SCHEMA_STR(PAL_SCHEMA_TUPLE_NAME(t)), PAL_SCHEMA_TUPLE_DESC(t)},

#define PALACE_SCHEMA_ENUM(EnumName, ...)                                          \
  enum class EnumName                                                              \
  {                                                                                \
    PAL_SCHEMA_FOR_EACH(PAL_SCHEMA_ENUM_EMIT_NAME, __VA_ARGS__)                    \
  };                                                                               \
  template <>                                                                      \
  struct ::palace::schema::utils::enum_descriptions<EnumName>                      \
  {                                                                                \
    static constexpr auto value =                                                  \
        ::std::array{PAL_SCHEMA_FOR_EACH(PAL_SCHEMA_ENUM_EMIT_PAIR, __VA_ARGS__)}; \
  }

#endif  // PALACE_SCHEMA_UTILS_ANNOTATIONS_HPP
