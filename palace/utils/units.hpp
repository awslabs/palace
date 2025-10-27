// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_UNITS_HPP
#define PALACE_UTILS_UNITS_HPP

#include <algorithm>
#include <array>
#include <vector>
#include "fem/gridfunction.hpp"
#include "linalg/vector.hpp"
#include "utils/constants.hpp"

namespace palace
{

class Units
{
  // Mesh unit length from config["Model"]["L0"] in [m]. Sets length of mesh coordinate.
  double L0_m;

  // Characteristic reference length [m] and time [ns] for (non)dimensionalization.
  // Note: In input, config["Model"]["Lc"] is measured in units of L0; here Lc_m is in [m].
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
  // static_assert(false) fixed as defect report for C++23 (P2593).
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
      return 1.0 / tc_ns;  // [GHz]
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
      static_assert(always_false<unit>::value, "ValueType unknown");
    }
  }

  // Return a copy of value scaled by given unit.
  template <ValueType unit, typename T>
  auto Dimensionalize(T value) const
  {
    return value * GetScaleFactor<unit>();
  }

  // Return a copy of array scaled by given unit.
  template <ValueType unit, typename T, std::size_t N>
  auto Dimensionalize(std::array<T, N> value) const
  {
    std::transform(value.begin(), value.end(), value.begin(),
                   [this](T v) { return Dimensionalize<unit>(v); });
    return value;
  }

  // Return a copy of vector scaled by given unit.
  template <ValueType unit, typename T>
  auto Dimensionalize(std::vector<T> value) const
  {
    std::transform(value.begin(), value.end(), value.begin(),
                   [this](T v) { return Dimensionalize<unit>(v); });
    return value;
  }

  // Return a copy Vector or ComplexVector scaled by given unit.
  template <ValueType unit, typename T>
  auto Dimensionalize(T value) const
      -> std::enable_if_t<std::is_same_v<T, Vector> || std::is_same_v<T, ComplexVector>, T>
  {
    value *= GetScaleFactor<unit>();
    return value;
  }

  // Return a copy of the gridfunction scaled by given unit.
  // Neighbor face data is also scaled.
  template <ValueType unit>
  auto Dimensionalize(mfem::ParGridFunction value) const
  {
    value *= GetScaleFactor<unit>();
    value.FaceNbrData() *= GetScaleFactor<unit>();
    return value;
  }

  // Return a copy of the gridfunction scaled by given unit.
  // Neighbor face data is also scaled.
  template <ValueType unit>
  auto Dimensionalize(GridFunction value) const
  {
    value *= GetScaleFactor<unit>();
    value.Real().FaceNbrData() *= GetScaleFactor<unit>();
    if (value.HasImag())
    {
      value.Imag().FaceNbrData() *= GetScaleFactor<unit>();
    }
    return value;
  }

  // Scale value in-place by given unit.
  template <ValueType unit, typename T>
  void DimensionalizeInPlace(T &value) const
  {
    value *= GetScaleFactor<unit>();
  }

  // Scale gridfunction and its neighbor face data in-place by given unit.
  template <ValueType unit>
  void DimensionalizeInPlace(mfem::ParGridFunction &value) const
  {
    value *= GetScaleFactor<unit>();
    value.FaceNbrData() *= GetScaleFactor<unit>();
  }

  // Scale complex gridfunction and its neighbor face data in-place by given unit.
  template <ValueType unit>
  auto Dimensionalize(GridFunction &value) const
  {
    value *= GetScaleFactor<unit>();
    value.Real().FaceNbrData() *= GetScaleFactor<unit>();
    if (value.HasImag())
    {
      value.Imag().FaceNbrData() *= GetScaleFactor<unit>();
    }
  }

  // Return a copy of value scaled by inverse of given unit.
  template <ValueType unit, typename T>
  auto Nondimensionalize(T value) const
  {
    return value / GetScaleFactor<unit>();
  }

  // Return a copy of array scaled by inverse of given unit.
  template <ValueType unit, typename T, std::size_t N>
  auto Nondimensionalize(std::array<T, N> value) const
  {
    std::transform(value.begin(), value.end(), value.begin(),
                   [this](T v) { return Nondimensionalize<unit>(v); });
    return value;
  }

  // Return a copy of array scaled by inverse of given unit.
  template <ValueType unit, typename T>
  auto Nondimensionalize(std::vector<T> value) const
  {
    std::transform(value.begin(), value.end(), value.begin(),
                   [this](T v) { return Nondimensionalize<unit>(v); });
    return value;
  }

  // Return a copy of Vector of ComplexVector scaled by inverse of given unit.
  template <ValueType unit, typename T>
  auto Nondimensionalize(T value) const
      -> std::enable_if_t<std::is_same_v<T, Vector> || std::is_same_v<T, ComplexVector>, T>
  {
    value *= (1.0 / GetScaleFactor<unit>());
    return value;
  }

  // Return a copy of the gridfunction scaled by given unit.
  // Neighbor face data is also scaled.
  template <ValueType unit>
  auto Nondimensionalize(mfem::ParGridFunction value) const
  {
    value *= (1.0 / GetScaleFactor<unit>());
    value.FaceNbrData() *= (1.0 / GetScaleFactor<unit>());
    return value;
  }

  // Return a copy of the gridfunction scaled by given unit.
  // Neighbor face data is also scaled.
  template <ValueType unit>
  auto Nondimensionalize(GridFunction value) const
  {
    value *= (1.0 / GetScaleFactor<unit>());
    value.Real().FaceNbrData() *= (1.0 / GetScaleFactor<unit>());
    if (value.HasImag())
    {
      value.Imag().FaceNbrData() *= (1.0 / GetScaleFactor<unit>());
    }
    return value;
  }

  // Scale value in-place by inverse of given unit.
  template <ValueType unit, typename T>
  void NondimensionalizeInPlace(T &value) const
  {
    value *= (1.0 / GetScaleFactor<unit>());
  }

  // Scale gridfunction and its neighbor face data in-place by inverse of given unit.
  template <ValueType unit>
  void NondimensionalizeInPlace(mfem::ParGridFunction &value) const
  {
    value *= (1.0 / GetScaleFactor<unit>());
    value.FaceNbrData() *= (1.0 / GetScaleFactor<unit>());
  }

  // Scale complex gridfunction and its neighbor face data in-place by inverse of given
  // unit.
  template <ValueType unit>
  auto Nondimensionalize(GridFunction &value) const
  {
    value *= (1.0 / GetScaleFactor<unit>());
    value.Real().FaceNbrData() *= (1.0 / GetScaleFactor<unit>());
    if (value.HasImag())
    {
      value.Imag().FaceNbrData() *= (1.0 / GetScaleFactor<unit>());
    }
  }
};

}  // namespace palace

#endif  // PALACE_UTILS_UNITS_HPP
