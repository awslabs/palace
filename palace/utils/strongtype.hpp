// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_STRONGTYPE_HPP
#define PALACE_UTILS_STRONGTYPE_HPP

#include <fmt/format.h>
#include <nlohmann/json.hpp>

template <typename Tag, typename T = std::size_t>
class StrongT
{
private:
  T value;

public:
  // Explicit initialization only
  explicit constexpr StrongT(T v) noexcept : value(v) {}

  // Explicit conversion only
  explicit constexpr operator T() const noexcept { return value; }

  [[nodiscard]] constexpr T get() const noexcept { return value; }

  // Comparison boilerplate
  constexpr bool operator==(const StrongT &other) const noexcept
  {
    return value == other.value;
  }
  constexpr bool operator!=(const StrongT &other) const noexcept
  {
    return value != other.value;
  }
  constexpr bool operator<(const StrongT &other) const noexcept
  {
    return value < other.value;
  }
  constexpr bool operator<=(const StrongT &other) const noexcept
  {
    return value <= other.value;
  }
  constexpr bool operator>(const StrongT &other) const noexcept
  {
    return value > other.value;
  }
  constexpr bool operator>=(const StrongT &other) const noexcept
  {
    return value >= other.value;
  }
};

// fmt for StrongT using default of underlying T
template <typename Tag, typename T>
struct fmt::formatter<StrongT<Tag, T>> : fmt::formatter<T>
{
  template <typename FormatContext>
  auto format(const StrongT<Tag, T> &v, FormatContext &ctx)
  {
    return fmt::formatter<T>::format(v.get(), ctx);
  }
};

// JSON serialization / deserialization for StrongT using underlying T
// Specialize full adl_serializer since we don't have a default constructor
namespace nlohmann
{
template <typename Tag, typename T>
struct adl_serializer<StrongT<Tag, T>>
{
  static StrongT<Tag, T> from_json(const json &j)
  {
    return StrongT<Tag, T>{j.template get<T>()};
  }

  static void to_json(json &j, const StrongT<Tag, T> &t) { j = t.get(); }
};
}  // namespace nlohmann

struct ExcitationIdxTag
{
};
using ExcitationIdx = StrongT<ExcitationIdxTag, std::size_t>;

#endif  // PALACE_UTILS_STRONGTYPE_HPP
