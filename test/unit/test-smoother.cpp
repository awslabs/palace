// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <memory>
#include <numbers>
#include <type_traits>
#include <utility>
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

// The tridiagonal [-1, 2 + shift i / n, -1] stencil on the n local rows of each process,
// counting its applications. The Jacobi-scaled spectrum is dense near its maximum, so an
// eigenvalue estimate from a random start needs many applications, and the eigenvectors
// change with the shift.
class CountingLaplacianOperator : public Operator
{
private:
  const double shift;

public:
  mutable int applications = 0;

  CountingLaplacianOperator(int n, double shift = 0.0) : Operator(n), shift(shift) {}

  double Diagonal(int i) const { return 2.0 + shift * i / height; }

  void AssembleDiagonal(Vector &diag) const override
  {
    diag.SetSize(height);
    for (int i = 0; i < height; i++)
    {
      diag(i) = Diagonal(i);
    }
  }

  void Mult(const Vector &x, Vector &y) const override
  {
    applications++;
    y.SetSize(height);
    const auto *X = x.HostRead();
    auto *Y = y.HostWrite();
    for (int i = 0; i < height; i++)
    {
      Y[i] = Diagonal(i) * X[i] - ((i > 0) ? X[i - 1] : 0.0) -
             ((i < height - 1) ? X[i + 1] : 0.0);
    }
  }

  void MultTranspose(const Vector &x, Vector &y) const override { Mult(x, y); }
};

// Set an operator, a nearby one with other eigenvectors, and the first again: each estimate
// after the first starts from the dominant vector of the previous one, so it needs fewer
// applications of the operator than the estimate of a new smoother from a random start,
// and gives the same smoother up to the accuracy of the estimates (the spectrum is
// clustered at its maximum).
template <typename OperType, typename MakeSmoother>
void CheckWarmStart(MakeSmoother make, const OperType &op, const OperType &op_near,
                    const CountingLaplacianOperator &A,
                    const CountingLaplacianOperator &A_near)
{
  using VecType =
      std::conditional_t<std::is_same_v<OperType, ComplexOperator>, ComplexVector, Vector>;
  auto smoother = make();
  smoother->SetOperator(op);
  for (const auto &[o, counter] : {std::pair{&op_near, &A_near}, std::pair{&op, &A}})
  {
    counter->applications = 0;
    smoother->SetOperator(*o);
    const int warm = counter->applications;
    auto reference = make();
    counter->applications = 0;
    reference->SetOperator(*o);
    const int cold = counter->applications;
    CAPTURE(cold, warm);
    CHECK(warm < cold);

    VecType x(o->Width()), y(o->Height()), y_ref(o->Height());
    x = 1.0;
    smoother->Mult(x, y);
    reference->Mult(x, y_ref);
    y -= y_ref;
    CHECK(linalg::Norml2(Mpi::World(), y) <= 1.0e-2 * linalg::Norml2(Mpi::World(), y_ref));
  }
}

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
  CheckSmoothers(A, 1.0 + 1.0 / std::numbers::sqrt2);
}

TEST_CASE("Smoother estimates support negative definite operators",
          "[smoother][Serial][Parallel]")
{
  // Negating A and its diagonal D leaves D⁻¹A, and so the spectral estimate, unchanged (the
  // Hermitian similarity scales by |D|⁻¹ᐟ²). An eliminated essential DOF with unit diagonal
  // adds the eigenvalue 1 to D⁻¹A without changing its maximum, also when the diagonal then
  // has mixed signs. The smoothed values change sign with D⁻¹, which the check absorbs by
  // applying the smoother to D [1, ..., 1]ᵀ.
  const double lambda_max = 1.0 + 1.0 / std::numbers::sqrt2;
  TestOperator A(-1.0), A_ess(-1.0, true), A_spd_ess(1.0, true);
  CheckSmoothers(A, lambda_max);
  CheckSmoothers(A_ess, lambda_max);
  CheckSmoothers(A_spd_ess, lambda_max);
}

TEST_CASE("Spectral norm estimates", "[smoother][Serial][Parallel]")
{
  // The eigenvalues of [100, 10; 10, 2] are 51 ± √2501. A Krylov space of dimension at most
  // three captures the spectrum of the block diagonal operator, so the estimates are exact
  // up to round-off. The returned vector is the dominant (right singular) vector.
  const double lambda_max = 51.0 + std::sqrt(2501.0);
  const auto comm = Mpi::World();
  TestOperator A, A_neg(-1.0), A_ind(-1.0, true);
  for (const auto *op : {&A, &A_neg, &A_ind})
  {
    ComplexWrapperOperator Ar(const_cast<TestOperator *>(op), nullptr),
        Ac(const_cast<TestOperator *>(op), const_cast<TestOperator *>(op));
    for (bool sym : {true, false})
    {
      Vector x0;
      CHECK_THAT(linalg::SpectralNorm(comm, *op, sym, &x0), WithinRel(lambda_max, 1.0e-10));
      Vector y(x0.Size());
      op->Mult(x0, y);
      CHECK_THAT(linalg::Norml2(comm, y), WithinRel(lambda_max, 1.0e-8));
      CHECK_THAT(linalg::SpectralNorm(comm, Ar, sym), WithinRel(lambda_max, 1.0e-10));
    }
    ComplexVector x0;
    CHECK_THAT(linalg::SpectralNorm(comm, Ac, false, &x0),
               WithinRel(std::numbers::sqrt2 * lambda_max, 1.0e-10));
    ComplexVector y(x0.Size());
    Ac.Mult(x0, y);
    CHECK_THAT(linalg::Norml2(comm, y),
               WithinRel(std::numbers::sqrt2 * lambda_max, 1.0e-8));
  }
}

TEST_CASE("Spectral norm estimates of clustered and indefinite spectra",
          "[smoother][Serial][Parallel]")
{
  const auto comm = Mpi::World();
  {
    // The eigenvalues 2 - 2 cos(kπ / (n + 1)) of the tridiagonal [-1, 2, -1] cluster at the
    // maximum, so the iteration restarts before it converges, and the estimate (the Ritz
    // value plus its residual norm) is only accurate to a fraction of a percent.
    constexpr int n = 256;
    CountingLaplacianOperator A(n);
    const double lambda = linalg::SpectralNorm(comm, A, true);
    CHECK(A.applications > 8);
    CHECK_THAT(lambda, WithinRel(2.0 + 2.0 * std::cos(std::numbers::pi / (n + 1)), 1.0e-2));
  }
  {
    // An isolated eigenvalue 0.99 converges first, while the cluster in [-1, -0.9] still
    // moves towards the eigenvalue of largest magnitude, -1.
    constexpr int n = 40;
    Vector d(n);
    d(0) = 0.99;
    for (int i = 1; i < n; i++)
    {
      d(i) = -1.0 + 0.1 * (i - 1) / (n - 2);
    }
    DiagonalOperator D(d);
    CHECK_THAT(linalg::SpectralNorm(comm, D, true), WithinRel(1.0, 5.0e-3));
  }
}

TEST_CASE("Repeated spectral norm estimates", "[smoother][Serial][Parallel]")
{
  const auto comm = Mpi::World();
  {
    // Estimates which start from an exact eigenvector of the operator.
    const double lambda_max = 51.0 + std::sqrt(2501.0);
    TestOperator A;
    ComplexWrapperOperator Ac(&A, &A);
    Vector x0;
    ComplexVector x0_c;
    for (int i = 0; i < 3; i++)
    {
      CHECK_THAT(linalg::SpectralNorm(comm, A, true, &x0), WithinRel(lambda_max, 1.0e-10));
      CHECK_THAT(linalg::SpectralNorm(comm, Ac, false, &x0_c),
                 WithinRel(std::numbers::sqrt2 * lambda_max, 1.0e-10));
    }
    // The Jacobi-scaled operator of the smoother has the eigenvalues 1 ± 1/√2.
    ChebyshevSmoother<Operator> chebyshev(comm, 1, 1, 1.0);
    for (int i = 0; i < 3; i++)
    {
      Check(chebyshev, A, 4.0 / (3.0 * (1.0 + 1.0 / std::numbers::sqrt2)));
    }
  }
  {
    // The dominant eigenvector of the previous operator has the eigenvalue 0.9, while the
    // eigenvalue 1.1 of an orthogonal eigenvector is now the largest. The small random
    // component of each starting vector, iterated a few times, lets the estimates find it
    // again after a few of them.
    constexpr int n = 4000;
    Vector d(n), x0(n);
    for (int i = 0; i < n; i++)
    {
      d(i) = 0.8 * i / n;
    }
    d(0) = 0.9;
    d(1) = 1.1;
    x0 = 0.0;
    x0(0) = 1.0;
    DiagonalOperator D(d);
    double lambda = 0.0;
    int estimates = 0;
    while (estimates < 20 && lambda < 1.1 * (1.0 - 1.0e-3))
    {
      lambda = linalg::SpectralNorm(comm, D, true, &x0);
      estimates++;
    }
    CAPTURE(estimates);
    CHECK(estimates <= 6);
    CHECK_THAT(lambda, WithinRel(1.1, 1.0e-2));
  }
  {
    // A zero starting vector gives a cold start.
    TestOperator A;
    Vector x0(A.Height());
    x0 = 0.0;
    CHECK_THAT(linalg::SpectralNorm(comm, A, true, &x0),
               WithinRel(51.0 + std::sqrt(2501.0), 1.0e-10));
  }
}

TEST_CASE("Chebyshev estimates start from the previous dominant vector",
          "[smoother][Serial][Parallel]")
{
  constexpr int n = 32;
  CountingLaplacianOperator A(n), A_near(n, 0.2);
  // A real operator takes the Hermitian estimate in real arithmetic, (1 + i) A the general
  // complex one.
  ComplexWrapperOperator Ar(&A, nullptr), Ac(&A, &A), Ar_near(&A_near, nullptr),
      Ac_near(&A_near, &A_near);
  const auto comm = Mpi::World();
  CheckWarmStart(
      [&]() { return std::make_unique<ChebyshevSmoother<Operator>>(comm, 1, 2, 1.0); },
      static_cast<const Operator &>(A), static_cast<const Operator &>(A_near), A, A_near);
  {
    // A complex operator with no imaginary part takes the same real estimate as the real
    // operator.
    ChebyshevSmoother<Operator> smoother(comm, 1, 2, 1.0);
    ChebyshevSmoother<ComplexOperator> smoother_c(comm, 1, 2, 1.0);
    A.applications = 0;
    smoother.SetOperator(A);
    const int real = A.applications;
    A.applications = 0;
    smoother_c.SetOperator(Ar);
    CHECK(A.applications == real);
  }
  for (const auto &[op, op_near] : {std::pair{&Ar, &Ar_near}, std::pair{&Ac, &Ac_near}})
  {
    const ComplexOperator &o = *op, &o_near = *op_near;
    CheckWarmStart(
        [&]()
        { return std::make_unique<ChebyshevSmoother<ComplexOperator>>(comm, 1, 2, 1.0); },
        o, o_near, A, A_near);
    CheckWarmStart(
        [&]()
        {
          return std::make_unique<ChebyshevSmoother1stKind<ComplexOperator>>(comm, 1, 2,
                                                                             1.0, 0.0);
        },
        o, o_near, A, A_near);
  }
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
