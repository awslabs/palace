// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "mms-enzyme.hpp"
#include "test-mms-driven.hpp"

#include <array>
#include <cmath>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace palace::testing;
using namespace Catch::Matchers;

namespace palace::testing::mms::driven
{
namespace
{

void CheckPointwiseAgreement(const MmsCase &enzyme, const MmsCase &analytic,
                             const std::vector<std::array<double, 3>> &points)
{
  // μ_r enters curl_curl; ε_r and ω enter the full source assembled by the solve. All
  // three must match for the cases to describe the same problem.
  REQUIRE(enzyme.mu_r == analytic.mu_r);
  REQUIRE(enzyme.epsilon_r == analytic.epsilon_r);
  REQUIRE(enzyme.omega == analytic.omega);
  for (const auto &point : points)
  {
    mfem::Vector x(3);
    for (int i = 0; i < 3; i++)
    {
      x[i] = point[i];
    }
    mfem::Vector enzyme_field(3), analytic_field(3);
    mfem::Vector enzyme_curl(3), analytic_curl(3);
    mfem::Vector enzyme_curl_curl(3), analytic_curl_curl(3);
    enzyme.electric_field(x, enzyme_field);
    analytic.electric_field(x, analytic_field);
    enzyme.curl(x, enzyme_curl);
    analytic.curl(x, analytic_curl);
    enzyme.curl_curl(x, enzyme_curl_curl);
    analytic.curl_curl(x, analytic_curl_curl);
    for (int i = 0; i < 3; i++)
    {
      CHECK_THAT(enzyme_field[i], WithinAbs(analytic_field[i], 1.0e-14));
      CHECK_THAT(enzyme_curl[i], WithinAbs(analytic_curl[i], 5.0e-13));
      // Second derivatives of the trigonometric field carry a factor π² ≈ 10.
      CHECK_THAT(enzyme_curl_curl[i], WithinAbs(analytic_curl_curl[i], 5.0e-12));
    }
  }
}

}  // namespace

// Enzyme backend: the manufactured fields and the cases generated from them, the
// counterparts of analytic::PolynomialCase(), analytic::TrigCase(), and the Neumann-face
// analytic::NeumannPolynomialCase() and analytic::NeumannTrigCase(). File-local, since
// only this file (compiled only with Enzyme) uses them.
namespace enzyme
{
namespace
{

// These electric fields are the only problem-specific mathematical expressions in the
// Enzyme backend. ∇×E and the strong-form term ∇×(μ⁻¹∇×E) are generated from them. They
// match the analytic backend: the bubble polynomial, and the trigonometric field plus the
// gradient of φ = sin πx sin πy sin πz / π (see test-mms-driven-analytic.cpp).
MmsCoordinates PolynomialField(const MmsCoordinates &x)
{
  auto bubble = [](double t) { return t * (1.0 - t); };
  return {bubble(x[1]) * bubble(x[2]), bubble(x[2]) * bubble(x[0]),
          bubble(x[0]) * bubble(x[1])};
}

MmsCoordinates TrigField(const MmsCoordinates &x)
{
  const double sx = std::sin(M_PI * x[0]), sy = std::sin(M_PI * x[1]),
               sz = std::sin(M_PI * x[2]);
  const double cx = std::cos(M_PI * x[0]), cy = std::cos(M_PI * x[1]),
               cz = std::cos(M_PI * x[2]);
  return {sy * sz + cx * sy * sz, sz * sx + sx * cy * sz, sx * sy + sx * sy * cz};
}

// Neumann-face fields, which vanish tangentially on every face except x = 1 (see
// test-mms-driven-analytic.cpp): E = (a(y)a(z), x a(z), x a(y)), and the smooth field with
// s = sin(πx/2) in place of sin πx plus the gradient of φ = s sin πy sin πz / π.
MmsCoordinates NeumannPolynomialField(const MmsCoordinates &x)
{
  auto bubble = [](double t) { return t * (1.0 - t); };
  return {bubble(x[1]) * bubble(x[2]), x[0] * bubble(x[2]), x[0] * bubble(x[1])};
}

MmsCoordinates NeumannTrigField(const MmsCoordinates &x)
{
  const double s = std::sin(0.5 * M_PI * x[0]), c = std::cos(0.5 * M_PI * x[0]);
  const double sy = std::sin(M_PI * x[1]), sz = std::sin(M_PI * x[2]);
  const double cy = std::cos(M_PI * x[1]), cz = std::cos(M_PI * x[2]);
  return {sy * sz + 0.5 * c * sy * sz, s * sz + s * cy * sz, s * sy + s * sy * cz};
}

const EnzymeMmsVector<PolynomialField> kPolynomialMms(kAnisotropicMuR);
const EnzymeMmsVector<TrigField> kTrigMms(kAnisotropicMuR);
const EnzymeMmsVector<NeumannPolynomialField> kNeumannPolynomialMms(kAnisotropicMuR);
const EnzymeMmsVector<NeumannTrigField> kNeumannTrigMms(kAnisotropicMuR);

template <typename EnzymeMms>
MmsCase MakeCase(const EnzymeMms &mms)
{
  const auto *mms_ptr = &mms;
  return {
      [mms_ptr](const mfem::Vector &x, mfem::Vector &field) { mms_ptr->Value(x, field); },
      [mms_ptr](const mfem::Vector &x, mfem::Vector &curl) { mms_ptr->Curl(x, curl); },
      [mms_ptr](const mfem::Vector &x, mfem::Vector &curl_curl)
      { mms_ptr->CurlOfScaledCurl(x, curl_curl); },
      mms.Permeability(),
      kAnisotropicEpsilonR,
      kOmega,
  };
}

const MmsCase &PolynomialCase()
{
  static const MmsCase mms = MakeCase(kPolynomialMms);
  return mms;
}

const MmsCase &TrigCase()
{
  static const MmsCase mms = MakeCase(kTrigMms);
  return mms;
}

const MmsCase &NeumannPolynomialCase()
{
  static const MmsCase mms = MakeCase(kNeumannPolynomialMms);
  return mms;
}

const MmsCase &NeumannTrigCase()
{
  static const MmsCase mms = MakeCase(kNeumannTrigMms);
  return mms;
}

}  // namespace
}  // namespace enzyme
}  // namespace palace::testing::mms::driven

using namespace palace::testing::mms::driven;

TEST_CASE("Driven MMS Enzyme quantities match the analytic cases",
          "[driven][mms][enzyme][Serial]")
{
  // Independent check of the hand-derived strong-form source: Enzyme differentiates only
  // the fields above, twice for ∇×(μ⁻¹∇×E).
  const std::vector<std::array<double, 3>> points = {
      {0.17, 0.23, 0.31}, {0.38, 0.41, 0.62}, {0.73, 0.14, 0.87}};
  CheckPointwiseAgreement(enzyme::PolynomialCase(), analytic::PolynomialCase(), points);
  CheckPointwiseAgreement(enzyme::TrigCase(), analytic::TrigCase(), points);

  // The Neumann-face fields (also used by the impedance problem), including a point on
  // the face x = 1, where the curl supplies the surface data n×(μ⁻¹∇×E). The end-to-end
  // solves of both problems follow below.
  auto neumann_points = points;
  neumann_points.push_back({1.0, 0.27, 0.58});
  CheckPointwiseAgreement(enzyme::NeumannPolynomialCase(),
                          analytic::NeumannPolynomialCase(), neumann_points);
  CheckPointwiseAgreement(enzyme::NeumannTrigCase(), analytic::NeumannTrigCase(),
                          neumann_points);
}

TEST_CASE("Driven MMS Enzyme reproduces the analytic errors",
          "[driven][mms][enzyme][Serial][Parallel]")
{
  // Run the full solve with Enzyme-generated and hand-derived data; errors must agree. The
  // smooth field keeps the Enzyme derivatives from being exact, unlike a low-degree
  // polynomial.
  CheckBackendAgreement(enzyme::TrigCase(), analytic::TrigCase());
}

TEST_CASE("Driven MMS Enzyme with a Neumann face reproduces the analytic errors",
          "[driven][mms][enzyme][Serial][Parallel]")
{
  // As above for the Neumann-face problem, which also takes the Enzyme curl into the
  // surface data on the natural face.
  CheckNeumannBackendAgreement(enzyme::NeumannTrigCase(), analytic::NeumannTrigCase());
}

TEST_CASE("Driven MMS Enzyme with an impedance face reproduces the analytic errors",
          "[driven][mms][enzyme][Serial][Parallel]")
{
  // As above for the impedance problem, one end-to-end solve per problem. It adds no
  // derivative coverage beyond the Neumann-face solve: the Enzyme curl enters only the
  // n×H part of the surface data, as there, and the γ terms use field values, which need
  // no derivative.
  CheckImpedanceBackendAgreement(enzyme::NeumannTrigCase(), analytic::NeumannTrigCase());
}
