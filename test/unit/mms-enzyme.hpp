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
extern int enzyme_const;

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

// Jacobian J_ij = ∂F_i/∂x_j of a vector field by one forward-mode pass per coordinate
// direction (J e_j is column j).
template <auto Potential>
MmsJacobian EvaluateJacobian(const MmsCoordinates &x)
{
  MmsJacobian jacobian{};
  for (int j = 0; j < kMmsDimension; j++)
  {
    MmsCoordinates direction{};
    direction[j] = 1.0;
    MmsCoordinates value{};
    MmsCoordinates jacobian_column{};
    __enzyme_fwddiff<void>(EvaluateVectorPotential<Potential>, enzyme_dup, &x, &direction,
                           enzyme_dupnoneed, &value, &jacobian_column);
    for (int i = 0; i < kMmsDimension; i++)
    {
      jacobian[i][j] = jacobian_column[i];
    }
  }
  return jacobian;
}

// ∇×F = (∂F_z/∂y − ∂F_y/∂z, ∂F_x/∂z − ∂F_z/∂x, ∂F_y/∂x − ∂F_x/∂y) from the Jacobian.
inline MmsCoordinates CurlFromJacobian(const MmsJacobian &J)
{
  return {J[2][1] - J[1][2], J[0][2] - J[2][0], J[1][0] - J[0][1]};
}

// H = μ⁻¹∇×F for a diagonal μ, written as a function of x alone (μ is a parameter) so
// that a second forward-mode pass can differentiate it. The inner derivative is itself an
// Enzyme forward-mode pass, so the outer pass is forward-over-forward.
template <auto Potential>
void EvaluateScaledCurl(const MmsCoordinates &x, const MmsCoordinates &mu,
                        MmsCoordinates &scaled_curl)
{
  const auto curl = CurlFromJacobian(EvaluateJacobian<Potential>(x));
  for (int i = 0; i < kMmsDimension; i++)
  {
    scaled_curl[i] = curl[i] / mu[i];
  }
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

// Test-only manufactured vector field F generated from a user-supplied expression: the
// vector potential A of the magnetostatic test or the electric field E of the driven test.
// The Jacobian J_ij = ∂F_i/∂x_j comes from one forward-mode pass per coordinate direction,
// and ∇×F is read off its antisymmetric part. The strong-form curl-curl term ∇×(μ⁻¹∇×F)
// differentiates H = μ⁻¹∇×F once more, as a whole, rather than assembling it from the 27
// second derivatives ∂²F_i/∂x_j∂x_k: only EvaluateScaledCurl knows how μ enters, so a full
// (for example rotated) permeability tensor would change that one function.
template <auto Potential>
class EnzymeMmsVector
{
  static_assert(
      std::is_invocable_r_v<MmsCoordinates, decltype(Potential), const MmsCoordinates &>,
      "The manufactured vector field must be MmsCoordinates(const MmsCoordinates &)");

public:
  explicit EnzymeMmsVector(const MmsCoordinates &mu) : mu(mu) {}

  const MmsCoordinates &Permeability() const { return mu; }

  void Value(const mfem::Vector &x, mfem::Vector &value) const
  {
    const auto field = Potential(internal::ToCoordinates(x));
    for (int i = 0; i < kMmsDimension; i++)
    {
      value[i] = field[i];
    }
  }

  MmsJacobian Jacobian(const mfem::Vector &x) const
  {
    return internal::EvaluateJacobian<Potential>(internal::ToCoordinates(x));
  }

  // ∇×F.
  void Curl(const mfem::Vector &x, mfem::Vector &curl) const
  {
    const auto value = internal::CurlFromJacobian(Jacobian(x));
    for (int i = 0; i < kMmsDimension; i++)
    {
      curl[i] = value[i];
    }
  }

  // The magnetostatic name for ∇×A.
  void MagneticFluxDensity(const mfem::Vector &x, mfem::Vector &flux_density) const
  {
    Curl(x, flux_density);
  }

  // ∇×(μ⁻¹∇×F): the Jacobian of H = μ⁻¹∇×F by forward-over-forward differentiation (μ
  // is passed as a constant, not differentiated), then its curl.
  void CurlOfScaledCurl(const mfem::Vector &x, mfem::Vector &curl_curl) const
  {
    const auto coordinates = internal::ToCoordinates(x);
    MmsJacobian scaled_curl_jacobian{};
    for (int j = 0; j < kMmsDimension; j++)
    {
      MmsCoordinates direction{};
      direction[j] = 1.0;
      MmsCoordinates scaled_curl{};
      MmsCoordinates column{};
      __enzyme_fwddiff<void>(internal::EvaluateScaledCurl<Potential>, enzyme_dup,
                             &coordinates, &direction, enzyme_const, &mu, enzyme_dupnoneed,
                             &scaled_curl, &column);
      for (int i = 0; i < kMmsDimension; i++)
      {
        scaled_curl_jacobian[i][j] = column[i];
      }
    }
    const auto value = internal::CurlFromJacobian(scaled_curl_jacobian);
    for (int i = 0; i < kMmsDimension; i++)
    {
      curl_curl[i] = value[i];
    }
  }

private:
  MmsCoordinates mu;
};

}  // namespace palace::testing

#endif  // PALACE_TEST_MMS_ENZYME_HPP
