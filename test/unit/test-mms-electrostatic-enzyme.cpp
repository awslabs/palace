// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "mms-enzyme.hpp"
#include "test-mms-electrostatic.hpp"

#include <cmath>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace palace::testing;
using namespace Catch::Matchers;

namespace palace::testing::mms::electrostatic
{
namespace
{

// These potentials are the only problem-specific mathematical expressions in the Enzyme
// backend. Electric field, charge density, and Neumann flux are generated from them.
double SinPotential(const MmsCoordinates &x)
{
  return std::sin(M_PI * x[0]) * std::sin(2.0 * M_PI * x[1]) * std::sin(M_PI * x[2]);
}

double CosPotential(const MmsCoordinates &x)
{
  return std::cos(M_PI * x[0]) * std::cos(2.0 * M_PI * x[1]) * std::cos(M_PI * x[2]);
}

double PolynomialPotential(const MmsCoordinates &x)
{
  return x[0] * x[0] + 2.0 * x[1] * x[1] + 3.0 * x[2] * x[2];
}

const EnzymeMmsScalar<SinPotential> kSinMms(kIsotropicEpsilonR);
const EnzymeMmsScalar<CosPotential> kCosMms(kIsotropicEpsilonR);
const EnzymeMmsScalar<PolynomialPotential> kPolynomialMms(kAnisotropicEpsilonR);
const EnzymeMmsScalar<PolynomialPotential> kCurvedPolynomialMms(kIsotropicEpsilonR);

template <typename EnzymeMms>
MmsCase MakeCase(const EnzymeMms &mms, bool nonzero_dirichlet,
                 ScalarFunction neumann_flux = {})
{
  const auto *mms_ptr = &mms;
  return {[mms_ptr](const mfem::Vector &x) { return mms_ptr->Value(x); },
          [mms_ptr](const mfem::Vector &x, mfem::Vector &electric)
          { mms_ptr->ElectricField(x, electric); },
          [mms_ptr](const mfem::Vector &x) { return mms_ptr->ChargeDensity(x); },
          mms.Permittivity(),
          nonzero_dirichlet,
          neumann_flux};
}

const MmsCase &HomogeneousCase()
{
  static const MmsCase mms = MakeCase(kSinMms, false);
  return mms;
}

const MmsCase &NonHomogeneousCase()
{
  static const MmsCase mms = MakeCase(kCosMms, true);
  return mms;
}

const MmsCase &PolynomialCase()
{
  static const MmsCase mms = MakeCase(kPolynomialMms, true,
                                      [](const mfem::Vector &x)
                                      {
                                        constexpr MmsCoordinates normal = {1.0, 0.0, 0.0};
                                        return kPolynomialMms.NormalFlux(x, normal);
                                      });
  return mms;
}

const MmsCase &CurvedCylinderCase()
{
  static const MmsCase mms = MakeCase(kCurvedPolynomialMms, true);
  return mms;
}

}  // namespace
}  // namespace palace::testing::mms::electrostatic

using namespace palace::testing::mms::electrostatic;

TEST_CASE("Electrostatic MMS Enzyme quantities match the trigonometric analytic cases",
          "[mms][electrostatic][enzyme][Serial]")
{
  const std::vector<std::array<double, 3>> points = {
      {0.17, 0.23, 0.31}, {0.38, 0.41, 0.62}, {0.73, 0.14, 0.87}};
  CheckPointwiseAgreement(HomogeneousCase(), analytic::HomogeneousCase(), points);
  CheckPointwiseAgreement(NonHomogeneousCase(), analytic::NonHomogeneousCase(), points);
}

TEST_CASE("Electrostatic MMS Enzyme quantities match the polynomial analytic case",
          "[mms][electrostatic][enzyme][Serial]")
{
  const std::vector<std::array<double, 3>> points = {
      {1.0, 0.2, 0.3}, {1.0, 0.5, 0.7}, {1.0, 0.9, 0.1}};
  const auto &analytic_polynomial = analytic::PolynomialCase();
  CheckPointwiseAgreement(PolynomialCase(), analytic_polynomial, points,
                          /*check_neumann=*/true);

  constexpr MmsCoordinates normal = {2.0 / 3.0, -1.0 / 3.0, 2.0 / 3.0};
  for (const auto &point : points)
  {
    mfem::Vector x(3);
    for (int i = 0; i < kMmsDimension; i++)
    {
      x[i] = point[i];
    }
    mfem::Vector analytic_electric(kMmsDimension);
    analytic_polynomial.electric_field(x, analytic_electric);
    double expected_flux = 0.0;
    for (int i = 0; i < kMmsDimension; i++)
    {
      expected_flux -= normal[i] * analytic_polynomial.epsilon_r[i] * analytic_electric[i];
    }
    CHECK_THAT(kPolynomialMms.NormalFlux(x, normal), WithinAbs(expected_flux, 5.0e-13));
  }
}

TEST_CASE("Electrostatic MMS Enzyme solves to small potential and field errors",
          "[mms][electrostatic][enzyme][Serial][Parallel]")
{
  CheckSmallErrors(HomogeneousCase(), NonHomogeneousCase());
}

TEST_CASE("Electrostatic MMS Enzyme is exact for a polynomial in the FE space",
          "[mms][electrostatic][enzyme][Serial][Parallel]")
{
  CheckPolynomialExactness(PolynomialCase());
}

TEST_CASE("Electrostatic MMS Enzyme handles a Neumann boundary",
          "[mms][electrostatic][enzyme][Serial][Parallel]")
{
  CheckNeumannBoundary(PolynomialCase());
}

TEST_CASE("Electrostatic MMS Enzyme converges at the optimal rate",
          "[mms][electrostatic][enzyme][Serial][Parallel]")
{
  CheckOptimalConvergence(NonHomogeneousCase());
}

TEST_CASE("Electrostatic MMS Enzyme converges on curved tetrahedra",
          "[mms][electrostatic][enzyme][Serial][Parallel]")
{
  CheckCurvedConvergence(CurvedCylinderCase());
}
