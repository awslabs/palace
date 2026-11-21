// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "linalg/initialguess.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "test-linalg-helpers.hpp"
#include "utils/communication.hpp"

namespace palace
{

using namespace Catch::Matchers;
using test::SetRandomTestVector;
using test::TridiagonalTestOperator;
using test::UseDevice;

// The minimal-residual combination of a few vectors, used for the initial guesses of
// frequency sweeps.
TEST_CASE("MinimalResidualCombination", "[initialguess][Serial][Parallel][GPU]")
{
  constexpr int n = 200;
  MPI_Comm comm = Mpi::World();
  TridiagonalTestOperator A(n);

  auto Check =
      [&]<typename OperType, typename VecType, typename ScalarType>(
          const OperType &Aop, const VecType &, ScalarType c1, ScalarType c2, ScalarType c3)
  {
    VecType X1(n), X2(n), X3(n), x(n), y(n), b(n), xs(n), zero(n);
    UseDevice(X1, X2, X3, x, y, b, xs, zero);
    SetRandomTestVector(X1, 1);
    SetRandomTestVector(X2, 2);
    SetRandomTestVector(X3, 3);
    zero = 0.0;
    linalg::SolutionHistory<OperType> history(4, n);
    auto Combine = [&](std::vector<const VecType *> X, const VecType &rhs, VecType &z)
    {
      // Store the vectors with the first one as the most recent.
      history.Clear();
      for (auto it = X.rbegin(); it != X.rend(); ++it)
      {
        history.Push(**it);
      }
      history.MinimalResidualCombination(comm, Aop, rhs, z);
    };
    auto Residual = [&](const VecType &rhs, const VecType &z)
    {
      VecType r(n);
      UseDevice(r);
      Aop.Mult(z, r);
      linalg::AXPBY(1.0, rhs, -1.0, r);
      return linalg::Norml2(comm, r);
    };
    auto RelDiff = [&](const VecType &u, const VecType &v)
    {
      VecType d(n);
      UseDevice(d);
      d = u;
      d -= v;
      return linalg::Norml2(comm, d) / linalg::Norml2(comm, v);
    };

    // A solution in the span is recovered.
    {
      INFO("Exact solution in the span");
      xs = X1;
      xs *= c1;
      xs.Add(c2, X2);
      xs.Add(c3, X3);
      Aop.Mult(xs, b);
      Combine({&X1, &X2, &X3}, b, x);
      CHECK(RelDiff(x, xs) <= 1.0e-10);
      CHECK(Residual(b, x) <= 1.0e-12 * linalg::Norml2(comm, b));
    }

    // Otherwise, the residual is never larger than with the first vector (the previous
    // solution) scaled or not, and decreases with more vectors; the combination is a
    // minimizer.
    {
      INFO("Right-hand side not in the span");
      SetRandomTestVector(b, 4);
      VecType u(n);
      UseDevice(u);
      Aop.Mult(X1, u);
      y = X1;
      y *= linalg::Dot(comm, b, u) / linalg::Dot(comm, u, u);
      const double res_1 = Residual(b, y);
      CHECK(Residual(b, X1) >= res_1);
      Combine({&X1}, b, x);
      CHECK_THAT(Residual(b, x), WithinRel(res_1, 1.0e-12));
      Combine({&X1, &X2}, b, x);
      const double res_2 = Residual(b, x);
      CHECK(res_2 <= res_1 * (1.0 + 1.0e-12));
      Combine({&X1, &X2, &X3}, b, x);
      const double res_3 = Residual(b, x);
      CHECK(res_3 <= res_2 * (1.0 + 1.0e-12));
      CHECK(res_3 <= linalg::Norml2(comm, b));
      for (const VecType *v : {&X1, &X2, &X3})
      {
        for (double t : {1.0e-3, -1.0e-3})
        {
          y = x;
          y.Add(t, *v);
          CHECK(Residual(b, y) >= res_3);
        }
      }

      // Numerically dependent vectors, and a zero first vector, are dropped.
      VecType X4(n);
      UseDevice(X4);
      X4 = X1;
      X4 *= 2.0;
      X4 -= X2;
      Combine({&X1, &X2}, b, y);
      Combine({&X1, &X2, &X4, &X1}, b, x);
      CHECK(RelDiff(x, y) <= 1.0e-8);
      Combine({&zero, &X1, &X2}, b, x);
      CHECK(RelDiff(x, y) <= 1.0e-8);

      // Only zero vectors give a zero combination.
      Combine({&zero, &zero}, b, x);
      CHECK(linalg::Norml2(comm, x) == 0.0);
    }

    // Samples of a smooth sweep, the most recent first, are nearly dependent, and so are
    // their images: a solution in their span (here the next sample of a cubic in t) is
    // still recovered to round-off, which needs the reorthogonalization.
    {
      INFO("Nearly dependent vectors");
      constexpr double delta = 3.0e-3;
      std::array<VecType, 4> Y, S;
      for (int k = 0; k < 4; k++)
      {
        Y[k].SetSize(n);
        UseDevice(Y[k]);
        SetRandomTestVector(Y[k], 5 + k);
      }
      auto Sample = [&](double t, VecType &z)
      {
        z.SetSize(n);
        UseDevice(z);
        z = Y[3];
        for (int k = 2; k >= 0; k--)
        {
          z *= t;
          z += Y[k];
        }
      };
      for (int j = 0; j < 4; j++)
      {
        Sample((4 - j) * delta, S[j]);
      }
      Sample(5.0 * delta, xs);
      Aop.Mult(xs, b);
      Combine({&S[0], &S[1], &S[2], &S[3]}, b, x);
      CHECK(Residual(b, x) <= 1.0e-12 * linalg::Norml2(comm, b));
      CHECK(RelDiff(x, xs) <= 1.0e-10);
    }

    // With a nearly dependent second vector and a large right-hand side, the coefficients
    // overflow, and so does the residual of the combination: the result is the scaled first
    // vector.
    {
      INFO("Overflowing coefficients");
      VecType Xd(n), u(n);
      UseDevice(Xd, u);
      Xd = X1;
      Xd.Add(1.0e-9, X2);
      xs = X1;
      xs += X2;
      xs *= 1.0e300;
      Aop.Mult(xs, b);
      Combine({&X1, &Xd}, b, x);
      Aop.Mult(X1, u);
      y = X1;
      y *= 1.0e-300 * linalg::Dot(comm, b, u) / linalg::Dot(comm, u, u);
      x *= 1.0e-300;
      CHECK(RelDiff(x, y) <= 1.0e-12);
    }
  };

  SECTION("Real")
  {
    Check(static_cast<const Operator &>(A), Vector(), 0.7, -1.3, 0.4);
  }

  SECTION("Complex")
  {
    TridiagonalTestOperator S(n, 0.5, 0.25, 3);
    ComplexWrapperOperator Ac(&A, &S);
    using C = std::complex<double>;
    Check(static_cast<const ComplexOperator &>(Ac), ComplexVector(), C(0.7, 0.2),
          C(-1.3, 0.5), C(0.1, -0.4));
  }
}

// The history of previous solutions used for the initial guesses of frequency sweeps.
TEST_CASE("SolutionHistory", "[initialguess][Serial][Parallel][GPU]")
{
  // The most recent solutions first, up to the capacity, until cleared.
  constexpr int n = 10;
  linalg::SolutionHistory<ComplexOperator> history(3, n);
  CHECK(history.Vectors().empty());
  ComplexVector x(n);
  UseDevice(x);
  auto First = [](const ComplexVector *v) { return v->Real().HostRead()[0]; };
  for (int k = 1; k <= 5; k++)
  {
    x = double(k);
    history.Push(x);
    CHECK(history.Size() == std::min<std::size_t>(k, 3));
  }
  auto X = history.Vectors();
  REQUIRE(X.size() == 3);
  CHECK(First(X[0]) == 5.0);
  CHECK(First(X[1]) == 4.0);
  CHECK(First(X[2]) == 3.0);
  x = 6.0;
  CHECK(First(X[0]) == 5.0);  // A copy.
  history.Clear();
  CHECK(history.Vectors().empty());
  history.Push(x);
  X = history.Vectors();
  REQUIRE(X.size() == 1);
  CHECK(First(X[0]) == 6.0);
  history.Release();
  history.Push(x);
  CHECK(history.Vectors().empty());
  linalg::SolutionHistory<ComplexOperator> none(0, n);
  none.Push(x);
  CHECK(none.Vectors().empty());
}

}  // namespace palace
