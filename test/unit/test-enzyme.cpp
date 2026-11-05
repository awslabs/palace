// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "mms-enzyme.hpp"

using namespace palace::testing;
using namespace Catch::Matchers;

namespace
{

double QuadraticPotential(const MmsCoordinates &x)
{
  return x[0] * x[0] + 2.0 * x[1] * x[1] + 3.0 * x[2] * x[2] + 4.0 * x[0] * x[1] +
         6.0 * x[0] * x[2] + 5.0 * x[1] * x[2];
}

}  // namespace

TEST_CASE("Enzyme generates gradient and Hessian from a scalar potential",
          "[enzyme][Serial]")
{
  const EnzymeMmsScalar<QuadraticPotential> mms({1.0, 1.0, 1.0});
  mfem::Vector x({2.0, -3.0, 0.5});
  const MmsCoordinates expected_gradient = {2.0 * x[0] + 4.0 * x[1] + 6.0 * x[2],
                                            4.0 * x[1] + 4.0 * x[0] + 5.0 * x[2],
                                            6.0 * x[2] + 6.0 * x[0] + 5.0 * x[1]};
  constexpr MmsHessian expected_hessian = {MmsCoordinates{2.0, 4.0, 6.0},
                                           MmsCoordinates{4.0, 4.0, 5.0},
                                           MmsCoordinates{6.0, 5.0, 6.0}};

  const auto gradient = mms.Gradient(x);
  const auto hessian = mms.Hessian(x);
  for (int i = 0; i < kMmsDimension; i++)
  {
    CHECK_THAT(gradient[i], WithinAbs(expected_gradient[i], 1.0e-14));
    for (int j = 0; j < kMmsDimension; j++)
    {
      CHECK_THAT(hessian[i][j], WithinAbs(expected_hessian[i][j], 1.0e-14));
    }
  }
}
