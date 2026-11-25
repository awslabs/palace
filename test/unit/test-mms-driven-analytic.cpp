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

// Neumann-face fields. x = 1 is natural, so E need only vanish tangentially on the five PEC
// faces: n×E = 0 on y = 0, 1 needs E_x = E_z = 0, on z = 0, 1 needs E_x = E_y = 0, and on
// x = 0 needs E_y = E_z = 0, while on x = 1 the trace (E_y, E_z) stays free. The x-factor
// of E_y and E_z therefore only has to vanish at x = 0 (x for the polynomial, sin(πx/2) for
// the smooth field), which also makes the Neumann data n×(μ⁻¹∇×E) on x = 1 nonzero.
//
// Polynomial: E = (a(y)a(z), x a(z), x a(y)), with a(t) = t(1 − t), a' = 1 − 2t, a'' = −2.
// Its degrees match the PEC bubble field (each component in ND_2 on hexahedra, total
// degree 4 on tetrahedra), so the same orders are exact. Its trace on x = 1 is
// (E_y, E_z) = (a(z), a(y)).
void NeumannPolynomialField(const mfem::Vector &x, mfem::Vector &E)
{
  E(0) = Bubble(x(1)) * Bubble(x(2));
  E(1) = x(0) * Bubble(x(2));
  E(2) = x(0) * Bubble(x(1));
}

double BubbleDerivative(double t)
{
  return 1.0 - 2.0 * t;
}

// ∇×E = (∂_y E_z − ∂_z E_y, ∂_z E_x − ∂_x E_z, ∂_x E_y − ∂_y E_x)
//     = (x (a'(y) − a'(z)), a(y) (a'(z) − 1), a(z) (1 − a'(y))).
void NeumannPolynomialCurl(const mfem::Vector &x, mfem::Vector &B)
{
  const double ay = Bubble(x(1)), az = Bubble(x(2));
  const double dy = BubbleDerivative(x(1)), dz = BubbleDerivative(x(2));
  B(0) = x(0) * (dy - dz);
  B(1) = ay * (dz - 1.0);
  B(2) = az * (1.0 - dy);
}

// With H = μ⁻¹B and a'' = −2: ∂_y H_z = 2a(z)/μ_z and ∂_z H_y = −2a(y)/μ_y; ∂_z H_x =
// 2x/μ_x and ∂_x H_z = 0; ∂_x H_y = 0 and ∂_y H_x = −2x/μ_x. So
//   ∇×(μ⁻¹B) = (2 (a(y)/μ_y + a(z)/μ_z), 2x/μ_x, 2x/μ_x).
// For μ = 1 this is −ΔE, as it must be for the divergence-free E.
void NeumannPolynomialCurlCurl(const mfem::Vector &x, mfem::Vector &f)
{
  const auto &mu = kAnisotropicMuR;
  f(0) = 2.0 * (Bubble(x(1)) / mu[1] + Bubble(x(2)) / mu[2]);
  f(1) = 2.0 * x(0) / mu[0];
  f(2) = 2.0 * x(0) / mu[0];
}

// Smooth: the PEC construction with sin πx replaced by s = sin(πx/2), which vanishes at
// x = 0 but not at x = 1, plus the gradient of φ = s sin πy sin πz / π (zero on the PEC
// faces, so its tangential gradient vanishes there; nonzero on x = 1). As for the PEC
// field, ∇φ gives every component a dependence on its own coordinate:
//   E = (sin πy sin πz, s sin πz, s sin πy) + ∇φ,
//   ∇φ = (½ cos(πx/2) sin πy sin πz, s cos πy sin πz, s sin πy cos πz).
void NeumannTrigField(const mfem::Vector &x, mfem::Vector &E)
{
  const double s = std::sin(0.5 * M_PI * x(0)), c = std::cos(0.5 * M_PI * x(0));
  const double sy = std::sin(M_PI * x(1)), sz = std::sin(M_PI * x(2));
  const double cy = std::cos(M_PI * x(1)), cz = std::cos(M_PI * x(2));
  E(0) = sy * sz + 0.5 * c * sy * sz;
  E(1) = s * sz + s * cy * sz;
  E(2) = s * sy + s * sy * cz;
}

// ∇φ has zero curl, so with s' = (π/2) cos(πx/2) the curl of the first part is
//   ∇×E = (π s (cos πy − cos πz), sin πy (π cos πz − s'), sin πz (s' − π cos πy)).
void NeumannTrigCurl(const mfem::Vector &x, mfem::Vector &B)
{
  const double s = std::sin(0.5 * M_PI * x(0)),
               ds = 0.5 * M_PI * std::cos(0.5 * M_PI * x(0));
  const double sy = std::sin(M_PI * x(1)), sz = std::sin(M_PI * x(2));
  const double cy = std::cos(M_PI * x(1)), cz = std::cos(M_PI * x(2));
  B(0) = M_PI * s * (cy - cz);
  B(1) = sy * (M_PI * cz - ds);
  B(2) = sz * (ds - M_PI * cy);
}

// With H = μ⁻¹B and s'' = −(π²/4) s: ∂_y B_z = π² sin πy sin πz and ∂_z B_y = −π² sin πy
// sin πz; ∂_z B_x = π² s sin πz and ∂_x B_z = −(π²/4) s sin πz; ∂_x B_y = (π²/4) s sin πy
// and ∂_y B_x = −π² s sin πy. So
//   ∇×(μ⁻¹B) = π² (sin πy sin πz (1/μ_y + 1/μ_z), s sin πz (1/μ_x + 1/(4μ_z)),
//                  s sin πy (1/μ_x + 1/(4μ_y))).
// For μ = 1 the y-component is (5π²/4) s sin πz = −Δ(s sin πz), as expected.
void NeumannTrigCurlCurl(const mfem::Vector &x, mfem::Vector &f)
{
  const auto &mu = kAnisotropicMuR;
  const double s = std::sin(0.5 * M_PI * x(0));
  const double sy = std::sin(M_PI * x(1)), sz = std::sin(M_PI * x(2));
  constexpr double pi2 = M_PI * M_PI;
  f(0) = pi2 * sy * sz * (1.0 / mu[1] + 1.0 / mu[2]);
  f(1) = pi2 * s * sz * (1.0 / mu[0] + 0.25 / mu[2]);
  f(2) = pi2 * s * sy * (1.0 / mu[0] + 0.25 / mu[1]);
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

const MmsCase &NeumannPolynomialCase()
{
  static const MmsCase mms{NeumannPolynomialField,    NeumannPolynomialCurl,
                           NeumannPolynomialCurlCurl, kAnisotropicMuR,
                           kAnisotropicEpsilonR,      kOmega};
  return mms;
}

const MmsCase &NeumannTrigCase()
{
  static const MmsCase mms{NeumannTrigField, NeumannTrigCurl,      NeumannTrigCurlCurl,
                           kAnisotropicMuR,  kAnisotropicEpsilonR, kOmega};
  return mms;
}

}  // namespace analytic
}  // namespace palace::testing::mms::driven

using namespace palace::testing::mms::driven;

TEST_CASE("Driven MMS analytic mass projection is exact for a polynomial in the FE space",
          "[driven][mms][analytic][Serial][Parallel]")
{
  CheckMassProjectionExactness(analytic::PolynomialCase());
}

TEST_CASE("Driven MMS analytic mass projection converges at the optimal rate",
          "[driven][mms][analytic][Serial][Parallel]")
{
  CheckMassProjectionConvergence(analytic::TrigCase());
}

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

TEST_CASE(
    "Driven MMS analytic with a Neumann face is exact for a polynomial in the FE space",
    "[driven][mms][analytic][Serial][Parallel]")
{
  CheckNeumannExactness(analytic::NeumannPolynomialCase());
}

TEST_CASE("Driven MMS analytic with a Neumann face converges at the optimal rate",
          "[driven][mms][analytic][Serial][Parallel]")
{
  CheckNeumannConvergence(analytic::NeumannTrigCase());
}
