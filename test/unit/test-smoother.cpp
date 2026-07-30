// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "linalg/chebyshev.hpp"
#include "linalg/jacobi.hpp"
#include "utils/communication.hpp"

namespace palace
{

using namespace Catch::Matchers;

namespace
{

class TestOperator : public Operator
{
public:
  TestOperator() : Operator(2) {}

  void AssembleDiagonal(Vector &diag) const override
  {
    diag.SetSize(2);
    diag[0] = 100.0;
    diag[1] = 2.0;
  }

  void Mult(const Vector &x, Vector &y) const override
  {
    y.SetSize(2);
    const double x0 = x[0], x1 = x[1];
    y[0] = 100.0 * x0 + 10.0 * x1;
    y[1] = 10.0 * x0 + 2.0 * x1;
  }

  void MultTranspose(const Vector &x, Vector &y) const override { Mult(x, y); }
};

template <typename Smoother>
void Check(Smoother &smoother, const Operator &A, double expected)
{
  smoother.SetOperator(A);
  Vector x(2), y(2);
  x[0] = 100.0;
  x[1] = 2.0;
  smoother.Mult(x, y);
  CHECK_THAT(y[0], WithinRel(expected, 2.0e-3));
  CHECK_THAT(y[1], WithinRel(expected, 2.0e-3));
}

template <typename Smoother>
void Check(Smoother &smoother, const ComplexOperator &A, double expected)
{
  smoother.SetOperator(A);
  ComplexVector x(2), y(2);
  x = 0.0;
  x.Real()[0] = 100.0;
  x.Real()[1] = 2.0;
  smoother.Mult(x, y);
  CHECK_THAT(y.Real()[0], WithinRel(expected, 2.0e-3));
  CHECK_THAT(y.Real()[1], WithinRel(expected, 2.0e-3));
  CHECK(y.Imag().Normlinf() == 0.0);
}

}  // namespace

TEST_CASE("Smoother estimates use the Hermitian Jacobi similarity",
          "[smoother][Serial][Parallel]")
{
  // D⁻¹A is not symmetric, while D⁻¹ᐟ²AD⁻¹ᐟ² has eigenvalues 1 ± 1/√2.
  TestOperator A;
  ComplexWrapperOperator Ac(&A, nullptr);
  const double lambda_max = 1.0 + 1.0 / std::sqrt(2.0);

  ChebyshevSmoother<Operator> chebyshev(Mpi::World(), 1, 1, 1.0);
  Check(chebyshev, A, 4.0 / (3.0 * lambda_max));
  ChebyshevSmoother<ComplexOperator> chebyshev_c(Mpi::World(), 1, 1, 1.0);
  Check(chebyshev_c, Ac, 4.0 / (3.0 * lambda_max));
  JacobiSmoother<Operator> jacobi(Mpi::World(), 0.0);
  Check(jacobi, A, 2.0 / lambda_max);
  JacobiSmoother<ComplexOperator> jacobi_c(Mpi::World(), 0.0);
  Check(jacobi_c, Ac, 2.0 / lambda_max);
}

}  // namespace palace
