// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "test-mms-driven.hpp"
#include "test-mms-magnetostatic.hpp"

#include <cmath>

#include <catch2/catch_test_macros.hpp>

namespace palace::testing::mms::driven
{
namespace analytic
{

// The driven test reuses the hand-derived vector fields of the magnetostatic analytic
// backend (test-mms-magnetostatic-analytic.cpp), which satisfy n×E = 0 on every face of
// the unit cube, together with their curls B = ∇×E. The strong-form source additionally
// needs ∇×(μ⁻¹B), derived below for the diagonal μ_r: with H_i = B_i/μ_i,
//   ∇×H = (∂_y H_z − ∂_z H_y, ∂_z H_x − ∂_x H_z, ∂_x H_y − ∂_y H_x).
//
// The polynomial is the bubble field (a(y)a(z), a(z)a(x), a(x)a(y)), a(t) = t(1 − t),
// unchanged. The smooth case adds a gradient to the trigonometric field:
//   E = (sin πy sin πz, sin πz sin πx, sin πx sin πy) + ∇φ,   φ = sin πx sin πy sin πz / π.
// The trigonometric field alone has an x-component independent of x (and likewise for y,
// z), and the hexahedral ND_p x-component, Q_{p-1,p,p}, has its reduced degree exactly in
// x. Such a field never probes the reduced direction, so its hexahedral field error
// converges one order too fast, at h^(p+1). ∇φ depends on every coordinate, vanishes
// tangentially on ∂Ω because φ = 0 there, and has zero curl, so ∇×E and ∇×(μ⁻¹∇×E) are
// unchanged. It also gives E a gradient component, which the mass term sees but the
// curl-curl term (and so the magnetostatic test) does not.

namespace
{

double Bubble(double t)
{
  return t * (1.0 - t);
}

// B = (a(x)(a'(y) − a'(z)), a(y)(a'(z) − a'(x)), a(z)(a'(x) − a'(y))) with a'' = −2, so
// for example ∂_y B_z = −a(z)a''(y) = 2a(z) and ∂_z B_y = a(y)a''(z) = −2a(y), giving
//   ∇×(μ⁻¹B) = 2 (a(y)/μ_y + a(z)/μ_z, a(z)/μ_z + a(x)/μ_x, a(x)/μ_x + a(y)/μ_y).
// For μ = 1 this is −ΔE, as expected for the divergence-free bubble field.
void PolynomialCurlCurl(const mfem::Vector &x, mfem::Vector &f)
{
  const auto &mu = kAnisotropicMuR;
  const double ax = Bubble(x(0)) / mu[0], ay = Bubble(x(1)) / mu[1],
               az = Bubble(x(2)) / mu[2];
  f(0) = 2.0 * (ay + az);
  f(1) = 2.0 * (az + ax);
  f(2) = 2.0 * (ax + ay);
}

// B = π (sin πx (cos πy − cos πz), sin πy (cos πz − cos πx), sin πz (cos πx − cos πy)), so
// for example ∂_y B_z = π² sin πy sin πz and ∂_z B_y = −π² sin πy sin πz, giving
//   ∇×(μ⁻¹B) = π² (sin πy sin πz (1/μ_y + 1/μ_z), sin πz sin πx (1/μ_z + 1/μ_x),
//                  sin πx sin πy (1/μ_x + 1/μ_y)).
// For μ = 1 this is 2π² times the trigonometric field, i.e. −ΔE for a divergence-free E.
void TrigCurlCurl(const mfem::Vector &x, mfem::Vector &f)
{
  const auto &mu = kAnisotropicMuR;
  const double sx = std::sin(M_PI * x(0)), sy = std::sin(M_PI * x(1)),
               sz = std::sin(M_PI * x(2));
  constexpr double pi2 = M_PI * M_PI;
  f(0) = pi2 * sy * sz * (1.0 / mu[1] + 1.0 / mu[2]);
  f(1) = pi2 * sz * sx * (1.0 / mu[2] + 1.0 / mu[0]);
  f(2) = pi2 * sx * sy * (1.0 / mu[0] + 1.0 / mu[1]);
}

void TrigElectricField(const mfem::Vector &x, mfem::Vector &E)
{
  magnetostatic::analytic::TrigCase().potential(x, E);
  const double sx = std::sin(M_PI * x(0)), sy = std::sin(M_PI * x(1)),
               sz = std::sin(M_PI * x(2));
  const double cx = std::cos(M_PI * x(0)), cy = std::cos(M_PI * x(1)),
               cz = std::cos(M_PI * x(2));
  E(0) += cx * sy * sz;
  E(1) += sx * cy * sz;
  E(2) += sx * sy * cz;
}

}  // namespace

const MmsCase &PolynomialCase()
{
  const auto &field = magnetostatic::analytic::PolynomialCase();
  static const MmsCase mms{field.potential, field.magnetic_flux_density, PolynomialCurlCurl,
                           kAnisotropicMuR, kAnisotropicEpsilonR,        kOmega};
  return mms;
}

const MmsCase &TrigCase()
{
  static const MmsCase mms{
      TrigElectricField,    magnetostatic::analytic::TrigCase().magnetic_flux_density,
      TrigCurlCurl,         kAnisotropicMuR,
      kAnisotropicEpsilonR, kOmega};
  return mms;
}

}  // namespace analytic
}  // namespace palace::testing::mms::driven

using namespace palace::testing::mms::driven;

TEST_CASE("Driven MMS analytic is exact for a polynomial in the FE space",
          "[driven][mms][analytic][Serial][Parallel]")
{
  CheckPolynomialExactness(analytic::PolynomialCase());
}

TEST_CASE("Driven MMS analytic converges at the optimal rate",
          "[driven][mms][analytic][Serial][Parallel]")
{
  CheckOptimalConvergence(analytic::TrigCase());
}
