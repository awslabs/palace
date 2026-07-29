// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <complex>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include "drivers/eigensolver.hpp"

using palace::BuildComplexSumOperator;
using palace::ComplexWrapperOperator;
using palace::EigenSolverBackend;
using palace::EigenvalueSolver;
using palace::NonlinearEigenSolver;
using palace::internal::BuildResponseCorrectedMass;
using palace::internal::EigenvalueToAngularFrequency;
using palace::internal::GetEigenvalueTarget;
using palace::internal::MaximumWeightModeAssignment;

TEST_CASE("Eigenmode response-corrected mass", "[eigensolver][surface-response][Serial]")
{
  mfem::DenseMatrix mass_real(2), mass_imag(2), response(2);
  mass_real = 0.0;
  mass_imag = 0.0;
  response = 0.0;
  mass_real(0, 0) = 2.0;
  mass_real(1, 1) = 3.0;
  mass_imag(0, 1) = 0.5;
  mass_imag(1, 0) = 0.5;
  response(0, 0) = -0.25;
  response(1, 1) = 0.75;

  ComplexWrapperOperator mass(&mass_real, &mass_imag);
  auto corrected = BuildResponseCorrectedMass(mass, response);
  REQUIRE(corrected.real);
  REQUIRE(corrected.op);
  CHECK(corrected.op->Real() == corrected.real.get());
  CHECK(corrected.op->Imag() == &mass_imag);

  mfem::Vector x({2.0, -1.0}), y(2);
  corrected.op->Real()->Mult(x, y);
  CHECK(y[0] == Catch::Approx(3.5));
  CHECK(y[1] == Catch::Approx(-3.75));

  auto sum = BuildComplexSumOperator({{2.0, corrected.op.get()}});
  palace::ComplexVector xc(2), yc(2);
  xc.Real() = x;
  xc.Imag() = 0.0;
  sum->Mult(xc, yc);
  CHECK(yc.Real()[0] == Catch::Approx(7.0));
  CHECK(yc.Real()[1] == Catch::Approx(-7.5));
  CHECK(yc.Imag()[0] == Catch::Approx(-1.0));
  CHECK(yc.Imag()[1] == Catch::Approx(2.0));
}

TEST_CASE("Eigenmode spectral conventions", "[eigensolver][surface-response][Serial]")
{
  using Which = EigenvalueSolver::WhichType;

  SECTION("Linear generalized eigenproblem")
  {
    auto target = GetEigenvalueTarget(3.0, false, EigenSolverBackend::SLEPC,
                                      NonlinearEigenSolver::HYBRID);
    CHECK(target.shift == std::complex<double>(9.0, 0.0));
    CHECK(target.which == Which::TARGET_REAL);

    target = GetEigenvalueTarget(3.0, false, EigenSolverBackend::ARPACK,
                                 NonlinearEigenSolver::HYBRID);
    CHECK(target.which == Which::LARGEST_REAL);
    CHECK(EigenvalueToAngularFrequency({9.0, 0.0}, false) ==
          std::complex<double>(3.0, 0.0));
  }

  SECTION("Polynomial and nonlinear eigenproblems")
  {
    auto target = GetEigenvalueTarget(3.0, true, EigenSolverBackend::SLEPC,
                                      NonlinearEigenSolver::HYBRID);
    CHECK(target.shift == std::complex<double>(0.0, 3.0));
    CHECK(target.which == Which::TARGET_IMAGINARY);

    target = GetEigenvalueTarget(3.0, true, EigenSolverBackend::SLEPC,
                                 NonlinearEigenSolver::SLP);
    CHECK(target.which == Which::TARGET_MAGNITUDE);

    target = GetEigenvalueTarget(3.0, true, EigenSolverBackend::ARPACK,
                                 NonlinearEigenSolver::HYBRID);
    CHECK(target.which == Which::SMALLEST_IMAGINARY);

    const auto omega = EigenvalueToAngularFrequency({-0.25, 3.0}, true);
    CHECK(omega.real() == Catch::Approx(3.0));
    CHECK(omega.imag() == Catch::Approx(0.25));
  }
}

TEST_CASE("Eigenmode maximum-weight assignment", "[eigensolver][surface-response][Serial]")
{
  SECTION("Empty")
  {
    CHECK(MaximumWeightModeAssignment({}).empty());
  }

  SECTION("Global optimum differs from row-greedy matching")
  {
    const auto assignment = MaximumWeightModeAssignment({{0.90, 0.80}, {0.85, 0.10}});
    REQUIRE(assignment.size() == 2);
    CHECK(assignment[0] == 1);
    CHECK(assignment[1] == 0);
  }

  SECTION("Rectangular")
  {
    const auto assignment =
        MaximumWeightModeAssignment({{0.20, 0.95, 0.10}, {0.80, 0.10, 0.70}});
    REQUIRE(assignment.size() == 2);
    CHECK(assignment[0] == 1);
    CHECK(assignment[1] == 0);
  }

  SECTION("Invalid")
  {
    CHECK_THROWS(MaximumWeightModeAssignment({{1.0}, {0.5}}));
    CHECK_THROWS(MaximumWeightModeAssignment({{1.0, 0.0}, {0.5}}));
    CHECK_THROWS(MaximumWeightModeAssignment({{1.0, -0.1}}));
  }
}
