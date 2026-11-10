// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "test-mms-magnetostatic.hpp"

#include <cmath>

#include <catch2/catch_test_macros.hpp>

namespace palace::testing::mms::magnetostatic
{
namespace
{

// These hand-derived quantities remain independent of Enzyme and serve as the analytic
// solver baseline and pointwise oracle for the generated backend.

// A = (a(y)a(z), a(z)a(x), a(x)a(y)) with a(t) = t(1 - t), so each component vanishes on
// the four faces where it is tangential. Its curl B is reproduced exactly (see
// CheckPolynomialExactness) at order 2 on hexahedra, where B_x ∈ Q_{2,1,1} lies in the
// order-1 Raviart-Thomas space, and at order 4 on tetrahedra, where B has total degree 3.
double Bubble(double t)
{
  return t * (1.0 - t);
}

double BubbleDerivative(double t)
{
  return 1.0 - 2.0 * t;
}

void PolynomialPotential(const mfem::Vector &x, mfem::Vector &A)
{
  A(0) = Bubble(x(1)) * Bubble(x(2));
  A(1) = Bubble(x(2)) * Bubble(x(0));
  A(2) = Bubble(x(0)) * Bubble(x(1));
}

// ∇×A = (a(x)(a'(y) - a'(z)), a(y)(a'(z) - a'(x)), a(z)(a'(x) - a'(y))).
void PolynomialFluxDensity(const mfem::Vector &x, mfem::Vector &B)
{
  B(0) = Bubble(x(0)) * (BubbleDerivative(x(1)) - BubbleDerivative(x(2)));
  B(1) = Bubble(x(1)) * (BubbleDerivative(x(2)) - BubbleDerivative(x(0)));
  B(2) = Bubble(x(2)) * (BubbleDerivative(x(0)) - BubbleDerivative(x(1)));
}

// A = (sin πy sin πz, sin πz sin πx, sin πx sin πy): smooth, tangentially zero on every
// face of the unit cube, and divergence-free.
void TrigPotential(const mfem::Vector &x, mfem::Vector &A)
{
  const double sx = std::sin(M_PI * x(0)), sy = std::sin(M_PI * x(1)),
               sz = std::sin(M_PI * x(2));
  A(0) = sy * sz;
  A(1) = sz * sx;
  A(2) = sx * sy;
}

// ∇×A = π (sin πx (cos πy - cos πz), sin πy (cos πz - cos πx), sin πz (cos πx - cos πy)).
void TrigFluxDensity(const mfem::Vector &x, mfem::Vector &B)
{
  const double sx = std::sin(M_PI * x(0)), sy = std::sin(M_PI * x(1)),
               sz = std::sin(M_PI * x(2));
  const double cx = std::cos(M_PI * x(0)), cy = std::cos(M_PI * x(1)),
               cz = std::cos(M_PI * x(2));
  B(0) = M_PI * sx * (cy - cz);
  B(1) = M_PI * sy * (cz - cx);
  B(2) = M_PI * sz * (cx - cy);
}

}  // namespace

namespace analytic
{

const MmsCase &PolynomialCase()
{
  static const MmsCase mms{PolynomialPotential, PolynomialFluxDensity, kAnisotropicMuR};
  return mms;
}

const MmsCase &TrigCase()
{
  static const MmsCase mms{TrigPotential, TrigFluxDensity, kAnisotropicMuR};
  return mms;
}

}  // namespace analytic
}  // namespace palace::testing::mms::magnetostatic

using namespace palace::testing::mms::magnetostatic;

TEST_CASE("Magnetostatic MMS analytic is exact for a polynomial in the FE space",
          "[magnetostatic][mms][analytic][Serial][Parallel]")
{
  CheckPolynomialExactness(analytic::PolynomialCase());
}

TEST_CASE("Magnetostatic MMS analytic converges at the optimal rate",
          "[magnetostatic][mms][analytic][Serial][Parallel]")
{
  CheckOptimalConvergence(analytic::TrigCase());
}
