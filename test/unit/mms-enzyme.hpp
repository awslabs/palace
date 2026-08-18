// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_TEST_MMS_ENZYME_HPP
#define PALACE_TEST_MMS_ENZYME_HPP

#include <array>
#include <type_traits>

#include <mfem.hpp>

template <typename ReturnType, typename... Args>
ReturnType __enzyme_autodiff(Args...);

template <typename ReturnType, typename... Args>
ReturnType __enzyme_fwddiff(Args...);

extern int enzyme_dup;
extern int enzyme_dupnoneed;

namespace palace::testing
{

inline constexpr int kMmsDimension = 3;
using MmsCoordinates = std::array<double, kMmsDimension>;
using MmsHessian = std::array<MmsCoordinates, kMmsDimension>;

namespace internal
{

template <auto Potential>
void EvaluatePotential(const MmsCoordinates &x, double &value)
{
  value = Potential(x);
}

template <auto Potential>
void EvaluateGradient(const MmsCoordinates &x, MmsCoordinates &gradient)
{
  gradient.fill(0.0);
  double value = 0.0;
  double seed = 1.0;
  __enzyme_autodiff<void>(EvaluatePotential<Potential>, enzyme_dup, &x, &gradient,
                          enzyme_dupnoneed, &value, &seed);
}

}  // namespace internal

// Test-only manufactured scalar field generated from a user-supplied potential. Enzyme
// differentiates fixed-size coordinate storage; MFEM vectors are adapted at the boundary.
template <auto Potential>
class EnzymeMmsScalar
{
  static_assert(std::is_invocable_r_v<double, decltype(Potential), const MmsCoordinates &>,
                "The manufactured potential must be double(const MmsCoordinates &)");

public:
  explicit EnzymeMmsScalar(const MmsCoordinates &epsilon) : epsilon(epsilon) {}

  const MmsCoordinates &Permittivity() const { return epsilon; }

  double Value(const mfem::Vector &x) const { return Potential(ToCoordinates(x)); }

  MmsCoordinates Gradient(const mfem::Vector &x) const
  {
    MmsCoordinates gradient{};
    internal::EvaluateGradient<Potential>(ToCoordinates(x), gradient);
    return gradient;
  }

  MmsHessian Hessian(const mfem::Vector &x) const
  {
    const auto coordinates = ToCoordinates(x);
    MmsHessian hessian{};
    for (int j = 0; j < kMmsDimension; j++)
    {
      MmsCoordinates direction{};
      direction[j] = 1.0;
      MmsCoordinates gradient{};
      MmsCoordinates hessian_column{};
      __enzyme_fwddiff<void>(internal::EvaluateGradient<Potential>, enzyme_dup,
                             &coordinates, &direction, enzyme_dupnoneed, &gradient,
                             &hessian_column);
      for (int i = 0; i < kMmsDimension; i++)
      {
        hessian[i][j] = hessian_column[i];
      }
    }
    return hessian;
  }

  void ElectricField(const mfem::Vector &x, mfem::Vector &electric) const
  {
    const auto gradient = Gradient(x);
    for (int i = 0; i < kMmsDimension; i++)
    {
      electric[i] = -gradient[i];
    }
  }

  double ChargeDensity(const mfem::Vector &x) const
  {
    const auto hessian = Hessian(x);
    double rho = 0.0;
    for (int i = 0; i < kMmsDimension; i++)
    {
      rho -= epsilon[i] * hessian[i][i];
    }
    return rho;
  }

  double NormalFlux(const mfem::Vector &x, const MmsCoordinates &normal) const
  {
    const auto gradient = Gradient(x);
    double flux = 0.0;
    for (int i = 0; i < kMmsDimension; i++)
    {
      flux += normal[i] * epsilon[i] * gradient[i];
    }
    return flux;
  }

private:
  static MmsCoordinates ToCoordinates(const mfem::Vector &x)
  {
    MFEM_VERIFY(x.Size() == kMmsDimension, "Expected three-dimensional MMS coordinates");
    MmsCoordinates coordinates{};
    for (int i = 0; i < kMmsDimension; i++)
    {
      coordinates[i] = x[i];
    }
    return coordinates;
  }

  MmsCoordinates epsilon;
};

}  // namespace palace::testing

#endif  // PALACE_TEST_MMS_ENZYME_HPP
