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
using MmsJacobian = std::array<MmsCoordinates, kMmsDimension>;

namespace internal
{

inline MmsCoordinates ToCoordinates(const mfem::Vector &x)
{
  MFEM_VERIFY(x.Size() == kMmsDimension, "Expected three-dimensional MMS coordinates");
  MmsCoordinates coordinates{};
  for (int i = 0; i < kMmsDimension; i++)
  {
    coordinates[i] = x[i];
  }
  return coordinates;
}

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

template <auto Potential>
void EvaluateVectorPotential(const MmsCoordinates &x, MmsCoordinates &value)
{
  value = Potential(x);
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

  double Value(const mfem::Vector &x) const
  {
    return Potential(internal::ToCoordinates(x));
  }

  MmsCoordinates Gradient(const mfem::Vector &x) const
  {
    MmsCoordinates gradient{};
    internal::EvaluateGradient<Potential>(internal::ToCoordinates(x), gradient);
    return gradient;
  }

  MmsHessian Hessian(const mfem::Vector &x) const
  {
    const auto coordinates = internal::ToCoordinates(x);
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
  MmsCoordinates epsilon;
};

// Test-only manufactured vector field generated from a user-supplied vector potential A.
// Only first derivatives are needed: the Jacobian J_ij = ∂A_i/∂x_j comes from one
// forward-mode pass per coordinate direction (J e_j is column j), and B = ∇×A is read off
// its antisymmetric part.
template <auto Potential>
class EnzymeMmsVector
{
  static_assert(
      std::is_invocable_r_v<MmsCoordinates, decltype(Potential), const MmsCoordinates &>,
      "The manufactured vector potential must be MmsCoordinates(const MmsCoordinates &)");

public:
  explicit EnzymeMmsVector(const MmsCoordinates &mu) : mu(mu) {}

  const MmsCoordinates &Permeability() const { return mu; }

  void Value(const mfem::Vector &x, mfem::Vector &potential) const
  {
    const auto value = Potential(internal::ToCoordinates(x));
    for (int i = 0; i < kMmsDimension; i++)
    {
      potential[i] = value[i];
    }
  }

  MmsJacobian Jacobian(const mfem::Vector &x) const
  {
    const auto coordinates = internal::ToCoordinates(x);
    MmsJacobian jacobian{};
    for (int j = 0; j < kMmsDimension; j++)
    {
      MmsCoordinates direction{};
      direction[j] = 1.0;
      MmsCoordinates value{};
      MmsCoordinates jacobian_column{};
      __enzyme_fwddiff<void>(internal::EvaluateVectorPotential<Potential>, enzyme_dup,
                             &coordinates, &direction, enzyme_dupnoneed, &value,
                             &jacobian_column);
      for (int i = 0; i < kMmsDimension; i++)
      {
        jacobian[i][j] = jacobian_column[i];
      }
    }
    return jacobian;
  }

  // B = ∇×A = (∂A_z/∂y - ∂A_y/∂z, ∂A_x/∂z - ∂A_z/∂x, ∂A_y/∂x - ∂A_x/∂y).
  void MagneticFluxDensity(const mfem::Vector &x, mfem::Vector &flux_density) const
  {
    const auto J = Jacobian(x);
    flux_density[0] = J[2][1] - J[1][2];
    flux_density[1] = J[0][2] - J[2][0];
    flux_density[2] = J[1][0] - J[0][1];
  }

private:
  MmsCoordinates mu;
};

}  // namespace palace::testing

#endif  // PALACE_TEST_MMS_ENZYME_HPP
