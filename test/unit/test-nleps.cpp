// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <complex>
#include <numbers>
#include <catch2/catch_test_macros.hpp>
#include "linalg/eps.hpp"

using namespace palace;

// The Newton step on the pole-cleared d(λ) T(λ) uses L = d′/d, and the line search compares
// residuals scaled by |d|.
TEST_CASE("Pole-cleared Newton terms", "[nleps][Serial]")
{
  const std::complex<double> p(-0.1, 1.5);
  const nleps::PoleTerms poles{{p, std::conj(p)}};
  auto d = [&](std::complex<double> l) { return (l - p) * (l - std::conj(p)); };

  const std::complex<double> l(-0.02, 1.7), s(0.0, 2.0);
  const double h = 1.0e-6;
  const std::complex<double> dd = (d(l + h) - d(l - h)) / (2.0 * h);
  CHECK(std::abs(poles.LogDerivative(l) - dd / d(l)) <= 1.0e-8 * std::abs(dd / d(l)));
  CHECK(std::abs(poles.RelativeScale(l, s) - std::abs(d(l) / d(s))) <= 1.0e-14);

  CHECK(poles.IsAtPole(p));
  CHECK_FALSE(poles.IsAtPole(l));

  // A sample node on a pole is moved off it: here the middle node of five.
  {
    const nleps::PoleTerms undamped{{std::complex<double>(0.0, 1.0)}};
    const std::complex<double> c(0.0, 1.0), r(0.0, 0.5);
    for (int j = 0; j < 5; j++)
    {
      CHECK_FALSE(undamped.IsAtPole(c + r * undamped.SampleNode(j, 5, c, r), 1.0e-6));
    }
    CHECK(undamped.IsAtPole(c + r * std::cos(std::numbers::pi * 5 / 10.0), 1.0e-6));
  }

  // Pole proximity is relative, so it does not depend on the frequency scale.
  for (double scale : {1.0e-7, 1.0, 1.0e7})
  {
    CAPTURE(scale);
    const nleps::PoleTerms scaled{{scale * p}};
    CHECK(scaled.IsAtPole(scale * p * (1.0 + 1.0e-9)));
    CHECK_FALSE(scaled.IsAtPole(scale * p * 0.95));
  }
  CHECK(poles.Separates(std::complex<double>(0.0, 1.0), s));
  CHECK_FALSE(poles.Separates(l, s));
}
