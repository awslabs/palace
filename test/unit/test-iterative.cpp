// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <type_traits>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "linalg/iterative.hpp"
#include "linalg/jacobi.hpp"
#include "linalg/operator.hpp"
#include "utils/communication.hpp"
#include "utils/device.hpp"
#include "utils/timer.hpp"

namespace palace
{

using namespace Catch::Matchers;

namespace
{

class SpdTestOperator : public Operator
{
public:
  SpdTestOperator() : Operator(2) {}

  // JacobiSmoother::SetOperator calls this, so it is required for the FGMRES case.
  void AssembleDiagonal(Vector &diag) const override
  {
    diag.SetSize(2);
    diag[0] = 4.0;
    diag[1] = 3.0;
  }

  void Mult(const Vector &x, Vector &y) const override
  {
    y.SetSize(2);
    y[0] = 4.0 * x[0] + x[1];
    y[1] = x[0] + 3.0 * x[1];
  }

  void MultTranspose(const Vector &x, Vector &y) const override { Mult(x, y); }
};

// The symmetric stencil [s, d_i, s] on the local rows of each process, with
// d_i = d (4 + i % 7).
class TridiagonalTestOperator : public Operator
{
  const double s, d;

public:
  TridiagonalTestOperator(int n, double s = -1.0, double d = 1.0) : Operator(n), s(s), d(d)
  {
  }

  void AssembleDiagonal(Vector &diag) const override
  {
    diag.SetSize(height);
    auto *D = diag.HostWrite();
    for (int i = 0; i < height; i++)
    {
      D[i] = d * (4.0 + i % 7);
    }
  }

  void Mult(const Vector &x, Vector &y) const override
  {
    y.SetSize(height);
    const auto *X = x.HostRead();
    auto *Y = y.HostWrite();
    for (int i = 0; i < height; i++)
    {
      Y[i] = d * (4.0 + i % 7) * X[i] +
             s * (((i > 0) ? X[i - 1] : 0.0) + ((i < height - 1) ? X[i + 1] : 0.0));
    }
  }

  void MultTranspose(const Vector &x, Vector &y) const override { Mult(x, y); }
};

// A Jacobi smoother which counts its applications. CgSolver does not fuse a subclass into
// its own updates, since it may change the application.
template <typename OperType>
class CountingJacobiSmoother : public JacobiSmoother<OperType>
{
  using VecType = typename Solver<OperType>::VecType;

public:
  mutable int applications = 0;
  using JacobiSmoother<OperType>::JacobiSmoother;
  void Mult(const VecType &x, VecType &y) const override
  {
    applications++;
    JacobiSmoother<OperType>::Mult(x, y);
  }
};

// Run a fixed number of CG iterations with a Jacobi preconditioner and return the number of
// preconditioner applications in the timing report.
template <typename OperType, typename VecType>
int SolveJacobiPcg(const std::type_identity_t<OperType> &op,
                   JacobiSmoother<OperType> &jacobi, const VecType &b, VecType &x,
                   int max_it)
{
  jacobi.SetOperator(op);
  CgSolver<OperType> ksp(Mpi::World(), 0);
  ksp.SetOperator(op);
  ksp.SetPreconditioner(jacobi);
  ksp.SetMaxIter(max_it);
  ksp.EnableTimer();
  x = 0.0;
  const int count = BlockTimer::GlobalTimer().Counts(Timer::KSP_PRECONDITIONER);
  ksp.Mult(b, x);
  CHECK(ksp.GetNumIterations() == max_it);
  return BlockTimer::GlobalTimer().Counts(Timer::KSP_PRECONDITIONER) - count;
}

}  // namespace

// A zero residual with eps = 0 used to fall through the strict convergence test and
// divide by zero on the first iteration, returning a NaN solution vector.
TEST_CASE("CgSolver zero residual", "[iterative][Serial][Parallel]")
{
  SpdTestOperator A;
  CgSolver<Operator> ksp(Mpi::World(), 0);
  ksp.SetOperator(A);
  ksp.SetRelTol(1.0e-10);
  ksp.SetMaxIter(10);

  SECTION("Zero right-hand side")
  {
    Vector b(2), x(2);
    b = 0.0;
    x = 0.0;
    ksp.Mult(b, x);
    CHECK(ksp.GetNumIterations() == 0);
    CHECK(ksp.GetConverged());
    CHECK_THAT(ksp.GetFinalRes(), WithinAbs(0.0, 1.0e-14));
    CHECK(std::isfinite(x[0]));
    CHECK(std::isfinite(x[1]));
  }

  SECTION("Exact initial guess")
  {
    // Regression coverage: eps > 0 here, so this passed before the zero-residual fix too.
    Vector b(2), x(2);
    b = 0.0;
    b[0] = 2.0;  // b = A x for x = (1, -2)
    b[1] = -5.0;
    x = 0.0;
    x[0] = 1.0;
    x[1] = -2.0;
    ksp.SetInitialGuess(true);
    ksp.Mult(b, x);
    CHECK(ksp.GetNumIterations() == 0);
    CHECK(ksp.GetConverged());
    CHECK_THAT(ksp.GetFinalRes(), WithinAbs(0.0, 1.0e-14));
    CHECK_THAT(x[0], WithinAbs(1.0, 1.0e-14));
    CHECK_THAT(x[1], WithinAbs(-2.0, 1.0e-14));
  }
}

TEST_CASE("GmresSolver zero residual", "[iterative][Serial][Parallel]")
{
  SpdTestOperator A;
  GmresSolver<Operator> ksp(Mpi::World(), 0);
  ksp.SetOperator(A);
  ksp.SetRelTol(1.0e-10);
  ksp.SetMaxIter(10);

  SECTION("Zero right-hand side")
  {
    Vector b(2), x(2);
    b = 0.0;
    x = 0.0;
    ksp.Mult(b, x);
    CHECK(ksp.GetNumIterations() == 0);
    CHECK(ksp.GetConverged());
    CHECK_THAT(ksp.GetFinalRes(), WithinAbs(0.0, 1.0e-14));
    CHECK(std::isfinite(x[0]));
    CHECK(std::isfinite(x[1]));
  }

  SECTION("Exact initial guess")
  {
    // Regression coverage: eps > 0 here, so this passed before the zero-residual fix too.
    Vector b(2), x(2);
    b = 0.0;
    b[0] = 2.0;
    b[1] = -5.0;
    x = 0.0;
    x[0] = 1.0;
    x[1] = -2.0;
    ksp.SetInitialGuess(true);
    ksp.Mult(b, x);
    CHECK(ksp.GetNumIterations() == 0);
    CHECK(ksp.GetConverged());
    CHECK_THAT(ksp.GetFinalRes(), WithinAbs(0.0, 1.0e-14));
    CHECK_THAT(x[0], WithinAbs(1.0, 1.0e-14));
    CHECK_THAT(x[1], WithinAbs(-2.0, 1.0e-14));
  }
}

TEST_CASE("FgmresSolver zero residual", "[iterative][Serial][Parallel]")
{
  SpdTestOperator A;
  JacobiSmoother<Operator> jacobi(Mpi::World());
  jacobi.SetOperator(A);
  FgmresSolver<Operator> ksp(Mpi::World(), 0);
  ksp.SetOperator(A);
  ksp.SetPreconditioner(jacobi);
  ksp.SetRelTol(1.0e-10);
  ksp.SetMaxIter(10);

  SECTION("Zero right-hand side")
  {
    Vector b(2), x(2);
    b = 0.0;
    x = 0.0;
    ksp.Mult(b, x);
    CHECK(ksp.GetNumIterations() == 0);
    CHECK(ksp.GetConverged());
    CHECK_THAT(ksp.GetFinalRes(), WithinAbs(0.0, 1.0e-14));
    CHECK(std::isfinite(x[0]));
    CHECK(std::isfinite(x[1]));
  }
}

// Without a preconditioner, the initial-guess path took the square root of the right-hand
// side norm, so the relative tolerance was measured against sqrt(||b||) instead of ||b||.
TEST_CASE("CgSolver initial residual with initial guess", "[iterative][Serial][Parallel]")
{
  SpdTestOperator A;
  CgSolver<Operator> ksp(Mpi::World(), 0);
  ksp.SetOperator(A);
  ksp.SetRelTol(1.0e-10);
  ksp.SetMaxIter(10);
  ksp.SetInitialGuess(true);

  Vector b(2), x(2);
  b[0] = 3.0;  // ||b|| = 5 on each rank
  b[1] = 4.0;
  x = 0.0;
  ksp.Mult(b, x);
  const double b_norm = 5.0 * std::sqrt(Mpi::Size(Mpi::World()));
  CHECK_THAT(ksp.GetInitialRes(), WithinRel(b_norm, 1.0e-14));
}

TEST_CASE("CgSolver zero residual (complex)", "[iterative][Serial][Parallel]")
{
  SpdTestOperator A;
  ComplexWrapperOperator Ac(&A, nullptr);
  CgSolver<ComplexOperator> ksp(Mpi::World(), 0);
  ksp.SetOperator(Ac);
  ksp.SetRelTol(1.0e-10);
  ksp.SetMaxIter(10);

  SECTION("Zero right-hand side")
  {
    ComplexVector b(2), x(2);
    b = 0.0;
    x = 0.0;
    ksp.Mult(b, x);
    CHECK(ksp.GetNumIterations() == 0);
    CHECK(ksp.GetConverged());
    CHECK_THAT(ksp.GetFinalRes(), WithinAbs(0.0, 1.0e-14));
    CHECK(std::isfinite(x.Real()[0]));
    CHECK(std::isfinite(x.Real()[1]));
  }
}

// On a single-threaded host, CG applies a Jacobi preconditioner with a real diagonal in the
// same pass as its residual update and the following inner product. Compare with the
// separate operations, which it uses for a subclass of JacobiSmoother.
TEST_CASE("CgSolver fused Jacobi update", "[iterative][Serial][Parallel]")
{
  constexpr int n = 101, max_it = 10;
  MPI_Comm comm = Mpi::World();
  TridiagonalTestOperator A(n);
  // The fused update replaces all but the first application of the preconditioner.
  const int fused_applications = utils::IsSerialHost() ? 1 : max_it + 1;

  SECTION("Real")
  {
    Vector b(n), x(n), y(n);
    for (int i = 0; i < n; i++)
    {
      b[i] = std::sin(1.0 + i + Mpi::Rank(comm));
    }
    JacobiSmoother<Operator> jacobi(comm);
    CountingJacobiSmoother<Operator> counting(comm);
    CHECK(SolveJacobiPcg(A, jacobi, b, x, max_it) == fused_applications);
    CHECK(SolveJacobiPcg(A, counting, b, y, max_it) == max_it + 1);
    CHECK(counting.applications == max_it + 1);
    // The same operations, up to how the compiler contracts and vectorizes the loops.
    y -= x;
    CHECK(linalg::Norml2(comm, y) <= 1.0e-12 * linalg::Norml2(comm, x));
  }

  SECTION("Complex")
  {
    // The imaginary part has a zero diagonal, so the Jacobi preconditioner is real while
    // the operator is not Hermitian and the CG step sizes are complex.
    TridiagonalTestOperator Ai(n, 0.5, 0.0);
    ComplexWrapperOperator Ac(&A, &Ai);
    ComplexVector b(n), x(n), y(n);
    for (int i = 0; i < n; i++)
    {
      b.Real()[i] = std::sin(1.0 + i + Mpi::Rank(comm));
      b.Imag()[i] = std::cos(2.0 * i);
    }
    JacobiSmoother<ComplexOperator> jacobi(comm);
    CountingJacobiSmoother<ComplexOperator> counting(comm);
    CHECK(SolveJacobiPcg(Ac, jacobi, b, x, max_it) == fused_applications);
    CHECK(SolveJacobiPcg(Ac, counting, b, y, max_it) == max_it + 1);
    // The fused update takes the imaginary part of (z, r) = rᴴ D r as zero, the separate
    // operations compute it up to rounding.
    y -= x;
    CHECK(linalg::Norml2(comm, y) <= 1.0e-12 * linalg::Norml2(comm, x));
  }

  SECTION("Complex diagonal")
  {
    // The diagonal is complex on the first process only, and all processes use the separate
    // operations.
    TridiagonalTestOperator Ai(n, 0.5, Mpi::Root(comm) ? 0.25 : 0.0);
    ComplexWrapperOperator Ac(&A, &Ai);
    ComplexVector b(n), x(n), y(n);
    for (int i = 0; i < n; i++)
    {
      b.Real()[i] = std::sin(1.0 + i + Mpi::Rank(comm));
      b.Imag()[i] = std::cos(2.0 * i);
    }
    JacobiSmoother<ComplexOperator> jacobi(comm);
    CountingJacobiSmoother<ComplexOperator> counting(comm);
    CHECK(SolveJacobiPcg(Ac, jacobi, b, x, max_it) == max_it + 1);
    CHECK(SolveJacobiPcg(Ac, counting, b, y, max_it) == max_it + 1);
    y -= x;
    CHECK(linalg::Norml2(comm, y) == 0.0);
  }
}

}  // namespace palace
