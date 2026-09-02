// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "test-mms-electrostatic.hpp"

#include <cmath>

#include <catch2/catch_test_macros.hpp>

namespace palace::testing::mms::electrostatic
{
namespace
{

// These hand-derived quantities remain independent of Enzyme and serve as the analytic
// solver baseline and pointwise oracle for the generated backend.
double SinPotential(const mfem::Vector &x)
{
  return std::sin(M_PI * x[0]) * std::sin(2.0 * M_PI * x[1]) * std::sin(M_PI * x[2]);
}

void SinElectricField(const mfem::Vector &x, mfem::Vector &electric)
{
  electric[0] =
      -M_PI * std::cos(M_PI * x[0]) * std::sin(2.0 * M_PI * x[1]) * std::sin(M_PI * x[2]);
  electric[1] = -2.0 * M_PI * std::sin(M_PI * x[0]) * std::cos(2.0 * M_PI * x[1]) *
                std::sin(M_PI * x[2]);
  electric[2] =
      -M_PI * std::sin(M_PI * x[0]) * std::sin(2.0 * M_PI * x[1]) * std::cos(M_PI * x[2]);
}

double SinChargeDensity(const mfem::Vector &x)
{
  return 6.0 * M_PI * M_PI * SinPotential(x);
}

double CosPotential(const mfem::Vector &x)
{
  return std::cos(M_PI * x[0]) * std::cos(2.0 * M_PI * x[1]) * std::cos(M_PI * x[2]);
}

void CosElectricField(const mfem::Vector &x, mfem::Vector &electric)
{
  electric[0] =
      M_PI * std::sin(M_PI * x[0]) * std::cos(2.0 * M_PI * x[1]) * std::cos(M_PI * x[2]);
  electric[1] = 2.0 * M_PI * std::cos(M_PI * x[0]) * std::sin(2.0 * M_PI * x[1]) *
                std::cos(M_PI * x[2]);
  electric[2] =
      M_PI * std::cos(M_PI * x[0]) * std::cos(2.0 * M_PI * x[1]) * std::sin(M_PI * x[2]);
}

double CosChargeDensity(const mfem::Vector &x)
{
  return 6.0 * M_PI * M_PI * CosPotential(x);
}

double PolynomialPotential(const mfem::Vector &x)
{
  return x[0] * x[0] + 2.0 * x[1] * x[1] + 3.0 * x[2] * x[2];
}

void PolynomialElectricField(const mfem::Vector &x, mfem::Vector &electric)
{
  electric[0] = -2.0 * x[0];
  electric[1] = -4.0 * x[1];
  electric[2] = -6.0 * x[2];
}

double PolynomialChargeDensity(const mfem::Vector &)
{
  return -2.0 * (kAnisotropicEpsilonR[0] + 2.0 * kAnisotropicEpsilonR[1] +
                 3.0 * kAnisotropicEpsilonR[2]);
}

double CurvedPolynomialChargeDensity(const mfem::Vector &)
{
  return -12.0;
}

double PolynomialNeumannFlux(const mfem::Vector &x)
{
  return 2.0 * kAnisotropicEpsilonR[0] * x[0];
}

}  // namespace

namespace analytic
{

const MmsCase &HomogeneousCase()
{
  static const MmsCase mms{
      SinPotential, SinElectricField, SinChargeDensity, kIsotropicEpsilonR, false, {}};
  return mms;
}

const MmsCase &NonHomogeneousCase()
{
  static const MmsCase mms{
      CosPotential, CosElectricField, CosChargeDensity, kIsotropicEpsilonR, true, {}};
  return mms;
}

const MmsCase &PolynomialCase()
{
  static const MmsCase mms{PolynomialPotential,
                           PolynomialElectricField,
                           PolynomialChargeDensity,
                           kAnisotropicEpsilonR,
                           true,
                           PolynomialNeumannFlux};
  return mms;
}

const MmsCase &CurvedCylinderCase()
{
  static const MmsCase mms{PolynomialPotential,
                           PolynomialElectricField,
                           CurvedPolynomialChargeDensity,
                           kIsotropicEpsilonR,
                           true,
                           {}};
  return mms;
}

}  // namespace analytic
}  // namespace palace::testing::mms::electrostatic

using namespace palace::testing::mms::electrostatic;

TEST_CASE("Electrostatic MMS analytic solves to small potential and field errors",
          "[electrostatic][mms][analytic][Serial][Parallel]")
{
  CheckSmallErrors(analytic::HomogeneousCase(), analytic::NonHomogeneousCase());
}

TEST_CASE("Electrostatic MMS analytic is exact for a polynomial in the FE space",
          "[electrostatic][mms][analytic][Serial][Parallel]")
{
  CheckPolynomialExactness(analytic::PolynomialCase());
}

TEST_CASE("Electrostatic MMS analytic handles a Neumann boundary",
          "[electrostatic][mms][analytic][Serial][Parallel]")
{
  CheckNeumannBoundary(analytic::PolynomialCase());
}

TEST_CASE("Electrostatic MMS analytic converges at the optimal rate",
          "[electrostatic][mms][analytic][Serial][Parallel]")
{
  CheckOptimalConvergence(analytic::NonHomogeneousCase());
}

TEST_CASE("Electrostatic MMS analytic converges on curved tetrahedra",
          "[electrostatic][mms][analytic][Serial][Parallel]")
{
  CheckCurvedConvergence(analytic::CurvedCylinderCase());
}
