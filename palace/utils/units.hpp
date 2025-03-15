// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_UNITS_HPP
#define PALACE_UTILS_UNITS_HPP

#include <algorithm>
#include <array>
#include <vector>

#include "utils/constants.hpp"

namespace palace
{

class Units
{
  // Mesh unit length from config["Model"]["L0"] in [m]. Sets length of mesh coordinate.
  double L0_m;

  // Characteristic reference length [m] and time [ns] for (non)dimensionalization.
  // Note: In input, config["Model"]["Lc"] is measured in units of L0; here Lc_m is in [m]
  double Lc_m;   // [m]
  double tc_ns;  // [ns]

  // Characteristic reference magnetic field strength Hc² = 1 / (Zc * Lc²) A/m (with Ec =
  // Hc Zc). Yields Pc = Hc² Zc Lc² = 1 W. Precompute.
  double Hc;

public:
  Units(double L0_, double Lc_)
    : L0_m(L0_), Lc_m(Lc_), tc_ns(1.0e9 * Lc_m / electromagnetics::c0_),
      Hc(1.0 / std::sqrt(electromagnetics::Z0_ * Lc_m * Lc_m)) {};

  // Return the mesh scaling factor in units L0 x [m] for mesh IO.
  double GetMeshLengthRelativeScale() const { return Lc_m / L0_m; }

  // Redimensionalize values for output. Outputs which depend on the fields assume a
  // characteristic reference magnetic field strength Hc such that Pc = 1 W, where Pc is the
  // characteristic reference power.
  enum class ValueType : std::uint8_t
  {
    TIME,          // [ns]
    FREQUENCY,     // [GHz]
    LENGTH,        // [m]
    IMPEDANCE,     // [Ω]
    INDUCTANCE,    // [H] = [Ωs]
    CAPACITANCE,   // [F] = [s/Ω]
    CONDUCTIVITY,  // [S/m]
    VOLTAGE,       // [V]
    CURRENT,       // [A]
    POWER,         // [W]
    ENERGY,        // [J]
    FIELD_E,       // [V/m]
    FIELD_D,       // [C/m²] = [A⋅s/m²]
    FIELD_H,       // [A/m]
    FIELD_B        // [Wb/m²] = [V⋅s/m²]
  };

  // Helper class to essentially allow static_assert(false) in constexpr if branch.
  // Obsolsee in current compiles, but can only remove in C++23.
  template <ValueType T>
  struct always_false : std::false_type
  {
  };

  template <ValueType unit>
  double GetScaleFactor() const
  {
    if constexpr (unit == ValueType::TIME)
    {
      return tc_ns;  // [ns]
    }
    else if constexpr (unit == ValueType::FREQUENCY)
    {
      return 1.0 / (2.0 * M_PI * tc_ns);  // [GHz/rad]
    }
    else if constexpr (unit == ValueType::LENGTH)
    {
      return Lc_m;  // [m]
    }
    else if constexpr (unit == ValueType::IMPEDANCE)
    {
      return electromagnetics::Z0_;  // [Ω]
    }
    else if constexpr (unit == ValueType::INDUCTANCE)
    {
      return electromagnetics::mu0_ * Lc_m;  // [H]
    }
    else if constexpr (unit == ValueType::CAPACITANCE)
    {
      return electromagnetics::epsilon0_ * Lc_m;  // [F]
    }
    else if constexpr (unit == ValueType::CONDUCTIVITY)
    {
      return 1.0 / (electromagnetics::Z0_ * Lc_m);  // [S/m]
    }
    else if constexpr (unit == ValueType::VOLTAGE)
    {
      return Hc * electromagnetics::Z0_ * Lc_m;  // [V]
    }
    else if constexpr (unit == ValueType::CURRENT)
    {
      return Hc * Lc_m;  // [A]
    }
    else if constexpr (unit == ValueType::POWER)
    {
      return Hc * Hc * electromagnetics::Z0_ * Lc_m * Lc_m;  // [W]
    }
    else if constexpr (unit == ValueType::ENERGY)
    {
      return Hc * Hc * electromagnetics::Z0_ * Lc_m * Lc_m * tc_ns;  // [J]
    }
    else if constexpr (unit == ValueType::FIELD_E)
    {
      return Hc * electromagnetics::Z0_;  // [V/m]
    }
    else if constexpr (unit == ValueType::FIELD_D)
    {
      return electromagnetics::epsilon0_ * Hc * electromagnetics::Z0_;  // [C/m²]
    }
    else if constexpr (unit == ValueType::FIELD_H)
    {
      return Hc;  // [A/m]
    }
    else if constexpr (unit == ValueType::FIELD_B)
    {
      return electromagnetics::mu0_ * Hc;  // [Wb/m²]
    }
    else
    {
      static_assert(always_false<unit>::value, "ValueType unkown");
    }
  }

  template <ValueType unit, typename T>
  auto Dimensionalize(T value) const
  {
    return value * GetScaleFactor<unit>();
  }

  template <ValueType unit, typename T, std::size_t N>
  auto Dimensionalize(const std::array<T, N> &value) const
  {
    auto out = value;
    std::transform(out.begin(), out.end(), out.begin(),
                   [this](T v) { return Dimensionalize<unit>(v); });
    return out;
  }

  template <ValueType unit, typename T>
  auto Dimensionalize(const std::vector<T> &value) const
  {
    auto out = value;
    std::transform(out.begin(), out.end(), out.begin(),
                   [this](T v) { return Dimensionalize<unit>(v); });
    return out;
  }

  template <ValueType unit, typename T>
  auto NonDimensionalize(T value) const
  {
    return value / GetScaleFactor<unit>();
  }

  template <ValueType unit, typename T, std::size_t N>
  auto NonDimensionalize(const std::array<T, N> &value) const
  {
    auto out = value;
    std::transform(out.begin(), out.end(), out.begin(),
                   [this](T v) { return NonDimensionalize<unit>(v); });
    return out;
  }

  template <ValueType unit, typename T>
  auto NonDimensionalize(const std::vector<T> &value) const
  {
    auto out = value;
    std::transform(out.begin(), out.end(), out.begin(),
                   [this](T v) { return NonDimensionalize<unit>(v); });
    return out;
  }
};

}  // namespace palace

#endif  // PALACE_UTILS_UNITS_HPP
