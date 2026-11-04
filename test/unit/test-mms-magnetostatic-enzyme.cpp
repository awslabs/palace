// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "mms-enzyme.hpp"
#include "test-mms-magnetostatic.hpp"

#include <cmath>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace palace::testing;
using namespace Catch::Matchers;

namespace palace::testing::mms::magnetostatic
{
namespace
{

void CheckPointwiseAgreement(const MmsCase &enzyme, const MmsCase &analytic,
                             const std::vector<std::array<double, 3>> &points)
{
  REQUIRE(enzyme.mu_r == analytic.mu_r);
  for (const auto &point : points)
  {
    mfem::Vector x(3);
    for (int i = 0; i < 3; i++)
    {
      x[i] = point[i];
    }
    mfem::Vector enzyme_potential(3), analytic_potential(3);
    mfem::Vector enzyme_flux_density(3), analytic_flux_density(3);
    enzyme.potential(x, enzyme_potential);
    analytic.potential(x, analytic_potential);
    enzyme.magnetic_flux_density(x, enzyme_flux_density);
    analytic.magnetic_flux_density(x, analytic_flux_density);
    for (int i = 0; i < 3; i++)
    {
      CHECK_THAT(enzyme_potential[i], WithinAbs(analytic_potential[i], 1.0e-14));
      CHECK_THAT(enzyme_flux_density[i], WithinAbs(analytic_flux_density[i], 5.0e-13));
    }
  }
}

// These vector potentials are the only problem-specific mathematical expressions in the
// Enzyme backend. The flux density B = ∇×A is generated from them.
MmsCoordinates PolynomialPotential(const MmsCoordinates &x)
{
  auto bubble = [](double t) { return t * (1.0 - t); };
  return {bubble(x[1]) * bubble(x[2]), bubble(x[2]) * bubble(x[0]),
          bubble(x[0]) * bubble(x[1])};
}

MmsCoordinates TrigPotential(const MmsCoordinates &x)
{
  const double sx = std::sin(M_PI * x[0]), sy = std::sin(M_PI * x[1]),
               sz = std::sin(M_PI * x[2]);
  return {sy * sz, sz * sx, sx * sy};
}

const EnzymeMmsVector<PolynomialPotential> kPolynomialMms(kAnisotropicMuR);
const EnzymeMmsVector<TrigPotential> kTrigMms(kAnisotropicMuR);

template <typename EnzymeMms>
MmsCase MakeCase(const EnzymeMms &mms)
{
  const auto *mms_ptr = &mms;
  return {
      [mms_ptr](const mfem::Vector &x, mfem::Vector &potential)
      { mms_ptr->Value(x, potential); },
      [mms_ptr](const mfem::Vector &x, mfem::Vector &flux_density)
      { mms_ptr->MagneticFluxDensity(x, flux_density); },
      mms.Permeability(),
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

}  // namespace
}  // namespace palace::testing::mms::magnetostatic

using namespace palace::testing::mms::magnetostatic;

TEST_CASE("Magnetostatic MMS Enzyme quantities match the analytic cases",
          "[magnetostatic][mms][enzyme][Serial]")
{
  const std::vector<std::array<double, 3>> points = {
      {0.17, 0.23, 0.31}, {0.38, 0.41, 0.62}, {0.73, 0.14, 0.87}};
  CheckPointwiseAgreement(PolynomialCase(), analytic::PolynomialCase(), points);
  CheckPointwiseAgreement(TrigCase(), analytic::TrigCase(), points);
}

TEST_CASE("Magnetostatic MMS Enzyme reproduces the analytic errors",
          "[magnetostatic][mms][enzyme][Serial][Parallel]")
{
  // Run the full solve with Enzyme-generated and hand-derived data; errors must agree. A
  // trigonometric potential keeps the Enzyme derivatives from being exact, unlike a
  // low-degree polynomial.
  CheckBackendAgreement(TrigCase(), analytic::TrigCase());
}
