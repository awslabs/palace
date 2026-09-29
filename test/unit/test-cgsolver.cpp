// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <vector>

#include <mfem.hpp>
#include "linalg/iterative.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

namespace palace
{
using namespace Catch::Matchers;

namespace
{

// A diagonal operator on one rank (the CG record is a property of the iteration, not of the
// distribution).
class DiagonalTestOperator : public Operator
{
private:
  std::vector<double> diagonal;

public:
  explicit DiagonalTestOperator(std::vector<double> d)
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
};

}  // namespace

// The self-consistent response-corrected electrostatic solve runs PCG on K + Pᵀ D P with
// AMG(K) as the preconditioner. Palace's CgSolver guards (Ap, p) > 0 by MFEM_ASSERT only,
// so a Release binary accepts an indefinite operator silently (sc-closure diagnostics
// 2026-09-29: the 120-degree corner model at p5). The recorded CG coefficient history must
// expose the negative curvature and the Ritz values of the preconditioned operator.
TEST_CASE("CgSolver curvature record", "[cgsolver][Serial]")
{
  // Unpreconditioned CG on a diagonal operator: the Ritz values converge to the
  // eigenvalues.
  auto Solve = [](std::vector<double> diagonal, bool record)
  {
    DiagonalTestOperator A(std::move(diagonal));
    CgSolver<Operator> cg(MPI_COMM_SELF, 0);
    cg.SetOperator(A);
    cg.SetRelTol(1.0e-12);
    cg.SetMaxIter(50);
    cg.EnableCgHistory(record);
    Vector b(A.Height()), x(A.Height());
    for (int i = 0; i < b.Size(); i++)
    {
      b(i) = 1.0 + 0.1 * i;
    }
    x = 0.0;
    cg.Mult(b, x);
    return std::tuple{cg.GetCgAlphaHistory(), cg.GetCgBetaRatioHistory(),
                      cg.GetCgNegativeCurvatureCount(), cg.GetConverged(), x};
  };

  SECTION("Positive definite: no negative curvature, positive alphas")
  {
    const auto [alpha, beta_ratio, negative, converged, x] =
        Solve({1.0, 2.0, 4.0, 8.0}, true);
    CHECK(converged);
    CHECK(negative == 0);
    REQUIRE(alpha.size() == beta_ratio.size());
    CHECK(alpha.size() >= 4);  // Four distinct eigenvalues: CG needs four steps.
    for (const double a : alpha)
    {
      CHECK(a > 0.0);
    }
    CHECK(beta_ratio[0] == 0.0);
    for (std::size_t k = 1; k < beta_ratio.size(); k++)
    {
      CHECK(beta_ratio[k] >= 0.0);
    }
    // The solution is the exact one.
    for (int i = 0; i < x.Size(); i++)
    {
      CHECK_THAT(x(i), WithinRel((1.0 + 0.1 * i) / std::pow(2.0, i), 1.0e-8));
    }
  }

  SECTION("Indefinite: the record shows a search direction with (Ap, p) <= 0")
  {
    // One negative eigenvalue: CG on this operator meets a search direction of negative
    // curvature (in exact arithmetic within the first four steps, the Krylov space reaching
    // the negative eigenvector).
    const auto [alpha, beta_ratio, negative, converged, x] =
        Solve({1.0, 2.0, -0.5, 8.0}, true);
    (void)converged;
    (void)x;
    CHECK(negative >= 1);
    bool negative_alpha = false;
    for (const double a : alpha)
    {
      negative_alpha = negative_alpha || a < 0.0;
    }
    CHECK(negative_alpha);
  }

  SECTION("Recording off: nothing is stored")
  {
    const auto [alpha, beta_ratio, negative, converged, x] =
        Solve({1.0, 2.0, -0.5, 8.0}, false);
    (void)converged;
    (void)x;
    CHECK(alpha.empty());
    CHECK(beta_ratio.empty());
    CHECK(negative == 0);
  }
}

}  // namespace palace
