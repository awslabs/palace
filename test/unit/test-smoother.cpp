// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "linalg/chebyshev.hpp"
#include "linalg/jacobi.hpp"
#include "utils/communication.hpp"

namespace palace
{

using namespace Catch::Matchers;

namespace
{

// Symmetric 2×2 operator [100, 10; 10, 2], scaled by sign (-1 gives a negative definite
// operator with a negative diagonal), optionally followed by an eliminated essential DOF
// with unit diagonal (decoupled from the other DOFs, as for Operator::DIAG_ONE).
class TestOperator : public Operator
{
private:
  const double sign;
  const bool eliminated_dof;

public:
  TestOperator(double sign = 1.0, bool eliminated_dof = false)
    : Operator(eliminated_dof ? 3 : 2), sign(sign), eliminated_dof(eliminated_dof)
  {
  }

  void AssembleDiagonal(Vector &diag) const override
  {
    diag.SetSize(height);
    diag[0] = sign * 100.0;
    diag[1] = sign * 2.0;
    if (eliminated_dof)
    {
      diag[2] = 1.0;
    }
  }

  void Mult(const Vector &x, Vector &y) const override
  {
    y.SetSize(height);
    const double x0 = x[0], x1 = x[1];
    y[0] = sign * (100.0 * x0 + 10.0 * x1);
    y[1] = sign * (10.0 * x0 + 2.0 * x1);
    if (eliminated_dof)
    {
      y[2] = x[2];
    }
  }

  void MultTranspose(const Vector &x, Vector &y) const override { Mult(x, y); }
};

class DiagonalTestOperator : public Operator
{
private:
  const bool has_zero_entry;

public:
  DiagonalTestOperator(bool has_zero_entry) : Operator(2), has_zero_entry(has_zero_entry) {}

  void AssembleDiagonal(Vector &diag) const override
  {
    diag.SetSize(2);
    diag[0] = 1.0;
    diag[1] = has_zero_entry ? 0.0 : 1.0;
  }

  void Mult(const Vector &x, Vector &y) const override
  {
    y.SetSize(2);
    y[0] = x[0];
    y[1] = has_zero_entry ? 0.0 : x[1];
  }

  void MultTranspose(const Vector &x, Vector &y) const override { Mult(x, y); }
};

// Apply the smoother to x = D [1, ..., 1]ᵀ for the diagonal D of A, so that D⁻¹x is a
// vector of ones: a single application of order one (Chebyshev) or a damped Jacobi step
// gives y = s D⁻¹x for the expected scaling s of the smoother.
template <typename Smoother>
void Check(Smoother &smoother, const Operator &A, double expected)
{
  smoother.SetOperator(A);
  Vector x(A.Height()), y(A.Height());
  A.AssembleDiagonal(x);
  smoother.Mult(x, y);
  for (int i = 0; i < y.Size(); i++)
  {
    CHECK_THAT(y[i], WithinRel(expected, 2.0e-3));
  }
}

template <typename Smoother>
void Check(Smoother &smoother, const ComplexOperator &A, double expected)
{
  smoother.SetOperator(A);
  ComplexVector x(A.Height()), y(A.Height());
  x = 0.0;
  A.Real()->AssembleDiagonal(x.Real());
  smoother.Mult(x, y);
  for (int i = 0; i < y.Size(); i++)
  {
    CHECK_THAT(y.Real()[i], WithinRel(expected, 2.0e-3));
  }
  CHECK(y.Imag().Normlinf() == 0.0);
}

void CheckSmoothers(Operator &A, double lambda_max)
{
  ComplexWrapperOperator Ac(&A, nullptr);
  ChebyshevSmoother<Operator> chebyshev(Mpi::World(), 1, 1, 1.0);
  Check(chebyshev, A, 4.0 / (3.0 * lambda_max));
  ChebyshevSmoother<ComplexOperator> chebyshev_c(Mpi::World(), 1, 1, 1.0);
  Check(chebyshev_c, Ac, 4.0 / (3.0 * lambda_max));
  JacobiSmoother<Operator> jacobi(Mpi::World(), 0.0);
  Check(jacobi, A, 2.0 / lambda_max);
  JacobiSmoother<ComplexOperator> jacobi_c(Mpi::World(), 0.0);
  Check(jacobi_c, Ac, 2.0 / lambda_max);
}

}  // namespace

TEST_CASE("Smoother estimates use the Hermitian Jacobi similarity",
          "[smoother][Serial][Parallel]")
{
  // D⁻¹A is not symmetric, while D⁻¹ᐟ²AD⁻¹ᐟ² has eigenvalues 1 ± 1/√2.
  TestOperator A;
  CheckSmoothers(A, 1.0 + 1.0 / std::sqrt(2.0));
}

TEST_CASE("Smoother estimates support negative definite operators",
          "[smoother][Serial][Parallel]")
{
  // Negating A and its diagonal D leaves D⁻¹A, and so the spectral estimate, unchanged (the
  // Hermitian similarity scales by |D|⁻¹ᐟ²). An eliminated essential DOF with unit diagonal
  // adds the eigenvalue 1 to D⁻¹A without changing its maximum, also when the diagonal then
  // has mixed signs. The smoothed values change sign with D⁻¹, which the check absorbs by
  // applying the smoother to D [1, ..., 1]ᵀ.
  const double lambda_max = 1.0 + 1.0 / std::sqrt(2.0);
  TestOperator A(-1.0), A_ess(-1.0, true), A_spd_ess(1.0, true);
  CheckSmoothers(A, lambda_max);
  CheckSmoothers(A_ess, lambda_max);
  CheckSmoothers(A_spd_ess, lambda_max);
}

TEST_CASE("Smoother estimates reject singular Jacobi diagonals",
          "[smoother][Serial][Parallel]")
{
  // Only one rank has an invalid entry so the parallel test exercises the collective
  // validation before entering a collective eigensolver.
  DiagonalTestOperator A(Mpi::Root(Mpi::World()));
  ComplexWrapperOperator Ac(&A, nullptr);

  ChebyshevSmoother<Operator> chebyshev(Mpi::World(), 1, 1, 1.0);
  CHECK_THROWS_WITH(chebyshev.SetOperator(A),
                    ContainsSubstring("finite, nonzero operator diagonal"));
  ChebyshevSmoother<ComplexOperator> chebyshev_c(Mpi::World(), 1, 1, 1.0);
  CHECK_THROWS_WITH(chebyshev_c.SetOperator(Ac),
                    ContainsSubstring("finite, nonzero real operator diagonal"));
  JacobiSmoother<Operator> jacobi(Mpi::World(), 0.0);
  CHECK_THROWS_WITH(jacobi.SetOperator(A),
                    ContainsSubstring("finite, nonzero operator diagonal"));
  JacobiSmoother<ComplexOperator> jacobi_c(Mpi::World(), 0.0);
  CHECK_THROWS_WITH(jacobi_c.SetOperator(Ac),
                    ContainsSubstring("finite, nonzero real operator diagonal"));
}

TEST_CASE("Automatic Jacobi rejects an invalid damping denominator",
          "[smoother][Serial][Parallel]")
{
  TestOperator A;
  JacobiSmoother<Operator> jacobi(Mpi::World(), 0.0, 0.0);
  CHECK_THROWS_WITH(jacobi.SetOperator(A),
                    ContainsSubstring("finite, strictly positive damping denominator"));
}

}  // namespace palace
