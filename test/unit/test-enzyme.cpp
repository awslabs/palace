// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <vector>

#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace Catch::Matchers;

template <typename ReturnType, typename... Args>
ReturnType __enzyme_autodiff(Args...);

namespace
{

template <typename VectorType>
void SquaredNorm(const VectorType &x, double &value)
{
  for (int i = 0; i < 3; i++)
  {
    value += x[i] * x[i];
  }
}

template <typename VectorType>
void CheckSquaredNormGradient()
{
  VectorType x(3);
  VectorType gradient(3);
  const std::array<double, 3> values{2.0, -3.0, 0.5};
  for (int i = 0; i < 3; i++)
  {
    x[i] = values[i];
    gradient[i] = 0.0;
  }

  double value = 0.0;
  double seed = 1.0;
  __enzyme_autodiff<void>(SquaredNorm<VectorType>, &x, &gradient, &value, &seed);

  for (int i = 0; i < 3; i++)
  {
    CHECK_THAT(gradient[i], WithinAbs(2.0 * values[i], 1.0e-14));
  }
}

}  // namespace

TEST_CASE("Enzyme differentiates MFEM-compatible vectors", "[enzyme][Serial]")
{
  SECTION("mfem::Vector")
  {
    CheckSquaredNormGradient<mfem::Vector>();
  }
  SECTION("std::vector")
  {
    CheckSquaredNormGradient<std::vector<double>>();
  }
}
