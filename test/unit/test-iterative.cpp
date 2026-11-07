// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
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
// d_i = d (4 + i % m).
class TridiagonalTestOperator : public Operator
{
  const double s, d;
  const int m;

public:
  TridiagonalTestOperator(int n, double s = -1.0, double d = 1.0, int m = 7)
    : Operator(n), s(s), d(d), m(m)
  {
  }

  void AssembleDiagonal(Vector &diag) const override
  {
    diag.SetSize(height);
    auto *D = diag.HostWrite();
    for (int i = 0; i < height; i++)
    {
      D[i] = d * (4.0 + i % m);
    }
  }

  void Mult(const Vector &x, Vector &y) const override
  {
    y.SetSize(height);
    const auto *X = x.HostRead();
    auto *Y = y.HostWrite();
    for (int i = 0; i < height; i++)
    {
      Y[i] = d * (4.0 + i % m) * X[i] +
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

// The skew-symmetric stencil [-s, 0, s] on the local rows of each process.
class SkewTridiagonalTestOperator : public Operator
{
  const double s;

public:
  SkewTridiagonalTestOperator(int n, double s) : Operator(n), s(s) {}

  void AssembleDiagonal(Vector &diag) const override
  {
    diag.SetSize(height);
    diag = 0.0;
  }

  void Mult(const Vector &x, Vector &y) const override
  {
    y.SetSize(height);
    const auto *X = x.HostRead();
    auto *Y = y.HostWrite();
    for (int i = 0; i < height; i++)
    {
      Y[i] = s * (((i < height - 1) ? X[i + 1] : 0.0) - ((i > 0) ? X[i - 1] : 0.0));
    }
  }
};

// The symmetric positive definite preconditioner [1/4, 1, 1/4] on the local rows of each
// process (on the real and imaginary parts of complex vectors), which is not a Jacobi
// smoother.
template <typename OperType>
class TridiagonalTestSolver : public Solver<OperType>
{
  using VecType = typename Solver<OperType>::VecType;

  static void Apply(const Vector &x, Vector &y)
  {
    const int n = x.Size();
    const auto *X = x.HostRead();
    auto *Y = y.HostWrite();
    for (int i = 0; i < n; i++)
    {
      Y[i] = X[i] + 0.25 * (((i > 0) ? X[i - 1] : 0.0) + ((i < n - 1) ? X[i + 1] : 0.0));
    }
  }

public:
  void SetOperator(const OperType &op) override {}

  void Mult(const VecType &x, VecType &y) const override
  {
    if constexpr (std::is_same_v<VecType, ComplexVector>)
    {
      Apply(x.Real(), y.Real());
      Apply(x.Imag(), y.Imag());
    }
    else
    {
      Apply(x, y);
    }
  }
};

// An operator which applies another one for its first applications and then returns zero,
// so that (Ap, p) = 0 in a later CG iteration.
template <typename OperType>
class BreakdownTestOperator : public OperType
{
  using VecType =
      std::conditional_t<std::is_same_v<OperType, ComplexOperator>, ComplexVector, Vector>;
  const OperType &A;
  const int applications;
  mutable int count = 0;

public:
  BreakdownTestOperator(const OperType &A, int applications)
    : OperType(A.Height()), A(A), applications(applications)
  {
  }

  void Mult(const VecType &x, VecType &y) const override
  {
    if (count++ < applications)
    {
      A.Mult(x, y);
    }
    else
    {
      y.SetSize(x.Size());
      y = 0.0;
    }
  }
};

// Preconditioned CG with separate vector operations and x updated with each search
// direction in its iteration, which stops as CgSolver does. Returns the number of
// iterations.
template <typename OperType, typename VecType>
int ReferencePcg(const OperType &A, const Solver<OperType> *B, const VecType &b, VecType &x,
                 bool initial_guess, int max_it, double rel_tol)
{
  MPI_Comm comm = Mpi::World();
  const int n = b.Size();
  VecType r(n), z(n), p(n), w(n);
  auto ApplyB = [B](const VecType &u, VecType &v)
  {
    if (B)
    {
      B->Mult(u, v);
    }
    else
    {
      v = u;
    }
  };
  if (initial_guess)
  {
    A.Mult(x, r);
    linalg::AXPBY(1.0, b, -1.0, r);
  }
  else
  {
    r = b;
    x = 0.0;
  }
  ApplyB(r, z);
  auto beta = linalg::Dot(comm, z, r), beta_prev = beta;
  double initial_res = std::sqrt(std::abs(beta));
  if (initial_guess)
  {
    ApplyB(b, w);
    initial_res = std::sqrt(std::abs(linalg::Dot(comm, w, b)));
  }
  int it = 0;
  for (; it < max_it && std::sqrt(std::abs(beta)) > rel_tol * initial_res; it++)
  {
    if (it == 0)
    {
      p = z;
    }
    else
    {
      p *= beta / beta_prev;
      p += z;
    }
    A.Mult(p, w);
    const auto alpha = beta / linalg::Dot(comm, w, p);
    x.Add(alpha, p);
    r.Add(-alpha, w);
    ApplyB(r, z);
    beta_prev = beta;
    beta = linalg::Dot(comm, z, r);
  }
  return it;
}

// Compare CgSolver with ReferencePcg on each way out of the iterations, with and without an
// initial guess x0. One solver is reused for all solves, so most start with the search
// direction of the previous one in its workspace.
template <typename OperType, typename VecType>
void CheckCgExitPaths(const OperType &A, const Solver<OperType> *B, const VecType &b,
                      const VecType &x0, bool converges)
{
  MPI_Comm comm = Mpi::World();
  const int n = b.Size();
  CgSolver<OperType> ksp(comm, 0);
  ksp.SetOperator(A);
  if (B)
  {
    ksp.SetPreconditioner(*B);
  }
  auto Solve = [&](const OperType &op, const VecType &rhs, bool initial_guess, int max_it,
                   double rel_tol, VecType &x)
  {
    ksp.SetOperator(op);
    ksp.SetInitialGuess(initial_guess);
    ksp.SetMaxIter(max_it);
    ksp.SetRelTol(rel_tol);
    x = x0;
    ksp.Mult(rhs, x);
    return ksp.GetNumIterations();
  };
  auto Reference =
      [&](const VecType &rhs, bool initial_guess, int max_it, double rel_tol, VecType &x)
  {
    x = x0;
    return ReferencePcg(A, B, rhs, x, initial_guess, max_it, rel_tol);
  };
  VecType x(n), y(n), zero(n);
  zero = 0.0;
  auto RelDiff = [&](const VecType &u, const VecType &v)
  {
    VecType d(n);
    d = u;
    d -= v;
    return linalg::Norml2(comm, d) / linalg::Norml2(comm, v);
  };

  for (bool initial_guess : {false, true})
  {
    CAPTURE(initial_guess);
    for (int max_it : {1, 2, 10})
    {
      CAPTURE(max_it);
      CHECK(Solve(A, b, initial_guess, max_it, 0.0, x) == max_it);
      CHECK(Reference(b, initial_guess, max_it, 0.0, y) == max_it);
      CHECK(RelDiff(x, y) <= 1.0e-12);
    }
    if (converges)
    {
      const int it = Solve(A, b, initial_guess, 1000, 1.0e-10, x);
      CHECK(it < 1000);
      CHECK(Reference(b, initial_guess, 1000, 1.0e-10, y) == it);
      CHECK(RelDiff(x, y) <= 1.0e-12);
      A.Mult(x, y);
      y -= b;
      CHECK(linalg::Norml2(comm, y) <= 1.0e-8 * linalg::Norml2(comm, b));
    }

    // Without iterations, x is the initial guess or zero.
    CHECK(Solve(A, b, initial_guess, 0, 0.0, x) == 0);
    x -= initial_guess ? x0 : zero;
    CHECK(linalg::Norml2(comm, x) == 0.0);

    // With (Ap, p) = 0 in the fourth iteration, x has been updated with the first three
    // search directions.
    BreakdownTestOperator<OperType> Ab(A, initial_guess ? 4 : 3);
    CHECK_THROWS(Solve(Ab, b, initial_guess, 10, 0.0, x));
    Reference(b, initial_guess, 3, 0.0, y);
    CHECK(RelDiff(x, y) <= 1.0e-12);
  }

  // A zero right-hand side or an exact initial guess converge without iterations.
  CHECK(Solve(A, zero, false, 10, 1.0e-10, x) == 0);
  CHECK(linalg::Norml2(comm, x) == 0.0);
  A.Mult(x0, y);
  CHECK(Solve(A, y, true, 10, 1.0e-10, x) == 0);
  x -= x0;
  CHECK(linalg::Norml2(comm, x) == 0.0);
}

template <typename OperType, typename VecType>
void CheckCgExitPaths(
    const OperType &A,
    const std::vector<std::pair<std::string, const Solver<OperType> *>> &pcs,
    const VecType &b, const VecType &x0, bool converges)
{
  for (const auto &[name, B] : pcs)
  {
    INFO("Preconditioner: " << name);
    CheckCgExitPaths(A, B, b, x0, converges);
  }
}

// Solve with CgSolver and ReferencePcg to a relative tolerance of 10⁻⁸ and compare.
template <typename OperType, typename VecType>
void CheckCgSolve(const OperType &A, const Solver<OperType> *B, const VecType &b,
                  const VecType &x0, bool initial_guess, int max_it)
{
  MPI_Comm comm = Mpi::World();
  CgSolver<OperType> ksp(comm, 0);
  ksp.SetOperator(A);
  if (B)
  {
    ksp.SetPreconditioner(*B);
  }
  ksp.SetInitialGuess(initial_guess);
  ksp.SetMaxIter(max_it);
  ksp.SetRelTol(1.0e-8);
  VecType x(b.Size()), y(b.Size());
  x = x0;
  y = x0;
  ksp.Mult(b, x);
  CHECK(ksp.GetNumIterations() == ReferencePcg(A, B, b, y, initial_guess, max_it, 1.0e-8));
  y -= x;
  CHECK(linalg::Norml2(comm, y) <= 1.0e-12 * std::max(linalg::Norml2(comm, x), 1.0e-300));
}

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
  constexpr int n = 1001, max_it = 10;
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

// CG updates x with each search direction in the next iteration, and with the last one
// after the iterations. Check each way out of the iterations for real and complex
// operators, with no preconditioner, Jacobi preconditioners on the fused and separate
// paths, and another preconditioner.
TEST_CASE("CgSolver exit paths", "[iterative][Serial][Parallel]")
{
  constexpr int n = 1001;
  MPI_Comm comm = Mpi::World();
  TridiagonalTestOperator A(n);

  SECTION("Real")
  {
    Vector b(n), x0(n);
    for (int i = 0; i < n; i++)
    {
      b[i] = std::sin(1.0 + i + Mpi::Rank(comm));
      x0[i] = std::cos(0.5 * i);
    }
    JacobiSmoother<Operator> jacobi(comm);
    CountingJacobiSmoother<Operator> counting(comm);
    TridiagonalTestSolver<Operator> tridiagonal;
    jacobi.SetOperator(A);
    counting.SetOperator(A);
    CheckCgExitPaths<Operator>(A,
                               {{"none", nullptr},
                                {"Jacobi", &jacobi},
                                {"Jacobi subclass", &counting},
                                {"tridiagonal", &tridiagonal}},
                               b, x0, true);
  }

  SECTION("Complex")
  {
    ComplexVector b(n), x0(n);
    for (int i = 0; i < n; i++)
    {
      b.Real()[i] = std::sin(1.0 + i + Mpi::Rank(comm));
      b.Imag()[i] = std::cos(2.0 * i);
      x0.Real()[i] = std::cos(0.5 * i);
      x0.Imag()[i] = std::sin(0.3 * i);
    }
    auto Check = [&](const ComplexOperator &Ac, bool converges)
    {
      JacobiSmoother<ComplexOperator> jacobi(comm);
      CountingJacobiSmoother<ComplexOperator> counting(comm);
      TridiagonalTestSolver<ComplexOperator> tridiagonal;
      jacobi.SetOperator(Ac);
      counting.SetOperator(Ac);
      CheckCgExitPaths<ComplexOperator>(Ac,
                                        {{"none", nullptr},
                                         {"Jacobi", &jacobi},
                                         {"Jacobi subclass", &counting},
                                         {"tridiagonal", &tridiagonal}},
                                        b, x0, converges);
    };

    // A Hermitian operator, for which the CG step sizes are real up to rounding.
    SkewTridiagonalTestOperator K(n, 0.5);
    ComplexWrapperOperator Ah(&A, &K);
    {
      INFO("Hermitian");
      Check(Ah, true);
    }

    // Complex symmetric operators, for which they are complex, with a real diagonal (a real
    // Jacobi preconditioner) or a complex one which is not a multiple of a real one, so
    // that the inner products (z, r) and their ratios are also complex.
    TridiagonalTestOperator S0(n, 0.5, 0.0), S1(n, 0.5, 0.25, 3);
    ComplexWrapperOperator As0(&A, &S0), As1(&A, &S1);
    {
      INFO("Complex symmetric, real diagonal");
      Check(As0, false);
    }
    {
      INFO("Complex symmetric, complex diagonal");
      Check(As1, false);
    }
  }
}

// Local sizes around the chunk size of the fused Jacobi path, and an empty first process.
TEST_CASE("CgSolver local sizes", "[iterative][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  for (int n_root : {0, 1, 7, 255, 256, 257, 513})
  {
    CAPTURE(n_root);
    const int n = Mpi::Root(comm) ? n_root : 300;
    TridiagonalTestOperator A(n);
    SkewTridiagonalTestOperator K(n, 0.5);
    ComplexWrapperOperator Ah(&A, &K);
    Vector b(n), x0(n);
    ComplexVector bc(n), xc0(n);
    for (int i = 0; i < n; i++)
    {
      b[i] = bc.Real()[i] = std::sin(1.0 + i + Mpi::Rank(comm));
      bc.Imag()[i] = std::cos(2.0 * i);
      x0[i] = xc0.Real()[i] = std::cos(0.5 * i);
      xc0.Imag()[i] = std::sin(0.3 * i);
    }
    JacobiSmoother<Operator> jacobi(comm);
    JacobiSmoother<ComplexOperator> cjacobi(comm);
    TridiagonalTestSolver<Operator> tridiagonal;
    TridiagonalTestSolver<ComplexOperator> ctridiagonal;
    jacobi.SetOperator(A);
    cjacobi.SetOperator(Ah);
    for (bool initial_guess : {false, true})
    {
      for (int max_it : {1, 3, 50})
      {
        CAPTURE(initial_guess, max_it);
        for (const auto *B :
             std::array<const Solver<Operator> *, 3>{nullptr, &jacobi, &tridiagonal})
        {
          CheckCgSolve<Operator>(A, B, b, x0, initial_guess, max_it);
        }
        for (const auto *B : std::array<const Solver<ComplexOperator> *, 3>{
                 nullptr, &cjacobi, &ctridiagonal})
        {
          CheckCgSolve<ComplexOperator>(Ah, B, bc, xc0, initial_guess, max_it);
        }
      }
    }
  }
}

}  // namespace palace
