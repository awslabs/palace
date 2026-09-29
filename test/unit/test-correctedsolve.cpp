// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <memory>
#include <vector>

#include <mfem.hpp>
#include "drivers/electrostaticsolver.hpp"
#include "linalg/iterative.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

namespace palace
{
using namespace Catch::Matchers;

namespace
{

// A diagonal "thin-metal stiffness" on one rank and its exact (Jacobi) preconditioner: the
// fail-closed rule is a property of the corrected PCG iteration, not of the distribution.
class DiagonalStiffness : public Operator
{
private:
  std::vector<double> diagonal;

public:
  explicit DiagonalStiffness(std::vector<double> d)
    : Operator(static_cast<int>(d.size())), diagonal(std::move(d))
  {
  }
  void Mult(const Vector &x, Vector &y) const override
  {
    for (int i = 0; i < x.Size(); i++)
    {
      y(i) = diagonal[i] * x(i);
    }
  }
  const std::vector<double> &Diagonal() const { return diagonal; }
};

class DiagonalInverse : public Solver<Operator>
{
private:
  std::vector<double> diagonal;

public:
  void SetOperator(const Operator &op) override
  {
    diagonal = dynamic_cast<const DiagonalStiffness &>(op).Diagonal();
    this->height = op.Height();
    this->width = op.Width();
  }
  void Mult(const Vector &x, Vector &y) const override
  {
    for (int i = 0; i < x.Size(); i++)
    {
      y(i) = x(i) / diagonal[i];
    }
  }
};

// A synthetic response correction Pᵀ D P of rank two: P picks two "trace" coefficients
// (dofs 0 and 1) and D is a 2 x 2 symmetric defect. D = diag(-d, -d) with d below the
// stiffness keeps K + Pᵀ D P positive definite; d above it makes the corrected operator
// indefinite while K stays the preconditioner (the hexagon's situation).
class SyntheticCorrection : public Operator
{
private:
  double defect;

public:
  SyntheticCorrection(int size, double defect) : Operator(size), defect(defect) {}
  void Mult(const Vector &x, Vector &y) const override
  {
    y = 0.0;
    y(0) = -defect * x(0);
    y(1) = -defect * x(1);
  }
};

}  // namespace

// The self-consistent corrected solve (drivers/electrostaticsolver.hpp SolveCorrectedField)
// must reject a corrected operator that is not positive definite: PCG on such an operator
// meets a search direction with (Ap, p) <= 0 and otherwise "converges" to a meaningless
// field (sc-closure diagnostics 2026-09-29, the 120-degree corner model at p5: five such
// directions, corrected domain correction 300x the fixed-trace one). A positive definite
// correction is accepted, its solve is byte-identical to the plain PCG solve, and the raw
// solver settings are restored.
TEST_CASE("SolveCorrectedField fail-closed on an indefinite correction",
          "[cgsolver][Serial]")
{
  const std::vector<double> stiffness = {1.0, 2.0, 3.0, 5.0, 8.0, 13.0};
  const int n = static_cast<int>(stiffness.size());
  DiagonalStiffness K(stiffness);
  auto MakeKsp = [&]()
  {
    auto cg = std::make_unique<CgSolver<Operator>>(MPI_COMM_SELF, 0);
    cg->SetRelTol(1.0e-10);
    cg->SetMaxIter(100);
    KspSolver ksp(std::move(cg), std::make_unique<DiagonalInverse>());
    ksp.SetOperators(K, K);
    return ksp;
  };
  Vector rhs(n);
  for (int i = 0; i < n; i++)
  {
    rhs(i) = 1.0 + 0.5 * i;
  }

  SECTION("Positive definite correction: accepted, exact solution, settings restored")
  {
    auto ksp = MakeKsp();
    SyntheticCorrection C(n, 0.5);  // K + C = diag(0.5, 1.5, 3, 5, 8, 13) > 0.
    SumOperator corrected(K, C);
    Vector x(n);
    x = 0.0;
    const auto record = SolveCorrectedField(ksp, K, corrected, 1.0e-8, false, rhs, x);
    CHECK(record.converged);
    CHECK(record.accepted);
    CHECK(record.negative_curvature_count == 0);
    CHECK(record.first_negative_curvature_iteration == -1);
    CHECK(record.iterations > 0);
    CHECK(record.ritz_min > 0.0);
    CHECK(record.ritz_max >= record.ritz_min);
    // Ritz values of K⁻¹ (K + C) = diag(0.5, 0.75, 1, 1, 1, 1) lie in [0.5, 1].
    CHECK_THAT(record.ritz_min, WithinAbs(0.5, 1.0e-6));
    CHECK_THAT(record.ritz_max, WithinAbs(1.0, 1.0e-6));
    CHECK_THAT(x(0), WithinRel(rhs(0) / 0.5, 1.0e-7));
    CHECK_THAT(x(1), WithinRel(rhs(1) / 1.5, 1.0e-7));
    CHECK_THAT(x(5), WithinRel(rhs(5) / 13.0, 1.0e-7));
    // The raw solver settings are restored: tolerance, operator, no CG history recorded.
    CHECK_THAT(ksp.GetRelTol(), WithinAbs(1.0e-10, 0.0));
    Vector y(n);
    y = 0.0;
    ksp.Mult(rhs, y);
    CHECK_THAT(y(0), WithinRel(rhs(0) / stiffness[0], 1.0e-9));
    CHECK(ksp.GetCgAlphaHistory().empty());
  }

  SECTION("Indefinite correction: PCG converges but the field is rejected")
  {
    auto ksp = MakeKsp();
    SyntheticCorrection C(n,
                          3.0);  // K + C = diag(-2, -1, 3, 5, 8, 13): two negative modes.
    SumOperator corrected(K, C);
    Vector x(n);
    x = 0.0;
    const auto record = SolveCorrectedField(ksp, K, corrected, 1.0e-8, false, rhs, x);
    // Preconditioned CG on a diagonal operator with three distinct preconditioned
    // eigenvalues (-2, -0.5, 1) terminates in three steps and reports convergence.
    CHECK(record.converged);
    CHECK(record.negative_curvature_count >= 1);
    CHECK(record.first_negative_curvature_iteration >= 0);
    CHECK(record.first_negative_curvature_iteration < record.iterations);
    CHECK(record.ritz_min < 0.0);
    CHECK_FALSE(record.accepted);
    CHECK_THAT(ksp.GetRelTol(), WithinAbs(1.0e-10, 0.0));
  }
}

}  // namespace palace
