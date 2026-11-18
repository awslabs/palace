// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "linalg/iterative.hpp"
#include "linalg/jacobi.hpp"
#include "linalg/operator.hpp"
#include "test-linalg-helpers.hpp"
#include "utils/communication.hpp"
#include "utils/device.hpp"
#include "utils/timer.hpp"

namespace palace
{

using namespace Catch::Matchers;
using test::SetRandomTestVector;
using test::TridiagonalTestOperator;
using test::UseDevice;

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
    auto *D = diag.HostWrite();
    D[0] = 4.0;
    D[1] = 3.0;
  }

  void Mult(const Vector &x, Vector &y) const override
  {
    y.SetSize(2);
    const auto *X = x.HostRead();
    auto *Y = y.HostWrite();
    Y[0] = 4.0 * X[0] + X[1];
    Y[1] = X[0] + 3.0 * X[1];
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

// GmresSolver or FgmresSolver with access to its Krylov basis.
template <typename KspType>
class KrylovBasisTestSolver : public KspType
{
public:
  using KspType::KspType;
  const auto &GetBasis() const { return this->V; }
};

// Solve without restarts to a relative tolerance of 10⁻¹⁰ from a zero initial guess and
// check the residual and the orthonormality of the Krylov basis of all iterations (left
// unchanged by the solve, without right preconditioning in GMRES). The basis is orthonormal
// to round-off with CGS2, while MGS and CGS lose orthogonality as the residual decreases,
// to about 10⁻⁵ here.
template <typename KspType, typename OperType, typename VecType>
void CheckGmresSolve(const std::type_identity_t<OperType> &A, const Solver<OperType> *B,
                     const VecType &b, Orthogonalization orthog)
{
  MPI_Comm comm = Mpi::World();
  KrylovBasisTestSolver<KspType> ksp(comm, 0);
  ksp.SetOperator(A);
  if (B)
  {
    ksp.SetPreconditioner(*B);
  }
  ksp.SetOrthogonalization(orthog);
  ksp.SetMaxIter(100);
  ksp.SetRelTol(1.0e-10);
  VecType x(b.Size()), r(b.Size());
  x = 0.0;
  ksp.Mult(b, x);
  CHECK(ksp.GetConverged());
  CHECK(ksp.GetNumIterations() < 50);
  A.Mult(x, r);
  r -= b;
  CHECK(linalg::Norml2(comm, r) <= 1.0e-9 * linalg::Norml2(comm, b));
  const auto &V = ksp.GetBasis();
  double loss = 0.0;
  for (int i = 0; i <= ksp.GetNumIterations(); i++)
  {
    for (int j = 0; j <= i; j++)
    {
      loss =
          std::max(loss, std::abs(linalg::Dot(comm, V[i], V[j]) - ((i == j) ? 1.0 : 0.0)));
    }
  }
  CHECK(loss <= ((orthog == Orthogonalization::CGS2) ? 1.0e-13 : 1.0e-3));
}

// An operator which counts its applications.
template <typename OperType>
class CountingTestOperator : public OperType
{
  using VecType =
      std::conditional_t<std::is_same_v<OperType, ComplexOperator>, ComplexVector, Vector>;
  const OperType &A;

public:
  mutable int applications = 0;
  CountingTestOperator(const OperType &A) : OperType(A.Height()), A(A) {}

  void Mult(const VecType &x, VecType &y) const override
  {
    applications++;
    A.Mult(x, y);
  }
};

// A Jacobi smoother with its result rounded to single precision, a model of a
// mixed-precision preconditioner: its error is relative to its argument, and it is not
// linear.
template <typename OperType>
class RoundingJacobiSmoother : public JacobiSmoother<OperType>
{
  using VecType = typename Solver<OperType>::VecType;

  static void Round(Vector &y)
  {
    auto *Y = y.HostReadWrite();
    for (int i = 0; i < y.Size(); i++)
    {
      Y[i] = static_cast<float>(Y[i]);
    }
  }

public:
  using JacobiSmoother<OperType>::JacobiSmoother;

  void Mult(const VecType &x, VecType &y) const override
  {
    JacobiSmoother<OperType>::Mult(x, y);
    if constexpr (std::is_same_v<VecType, ComplexVector>)
    {
      Round(y.Real());
      Round(y.Imag());
    }
    else
    {
      Round(y);
    }
  }
};

// A GMRES (with the given preconditioner side) or FGMRES solver and its preconditioner.
template <typename OperType>
struct GmresTestCase
{
  std::string name;
  const Solver<OperType> *B;
  PreconditionerSide side;
  bool flexible;

  bool Left() const { return B && !flexible && side == PreconditionerSide::LEFT; }
};

// The statistics of a solve, with the numbers of applications of the operator and the
// preconditioner.
struct GmresTestStats
{
  int it, operator_applications, preconditioner_applications;
  bool converged;
  double initial_res, final_res;
};

template <typename OperType, typename VecType>
GmresTestStats SolveGmres(const GmresTestCase<OperType> &c, const OperType &A,
                          const VecType &b, VecType &x, bool initial_guess, double rel_tol,
                          int max_it, int restart_dim = -1)
{
  MPI_Comm comm = Mpi::World();
  std::unique_ptr<GmresSolver<OperType>> ksp;
  if (c.flexible)
  {
    ksp = std::make_unique<FgmresSolver<OperType>>(comm, 0);
  }
  else
  {
    ksp = std::make_unique<GmresSolver<OperType>>(comm, 0);
    ksp->SetPreconditionerSide(c.side);
  }
  CountingTestOperator<OperType> Ac(A);
  ksp->SetOperator(Ac);
  if (c.B)
  {
    ksp->SetPreconditioner(*c.B);
  }
  ksp->SetInitialGuess(initial_guess);
  ksp->SetRelTol(rel_tol);
  ksp->SetMaxIter(max_it);
  ksp->SetRestartDim(restart_dim);
  ksp->EnableTimer();
  const int count = BlockTimer::GlobalTimer().Counts(Timer::KSP_PRECONDITIONER);
  ksp->Mult(b, x);
  return {ksp->GetNumIterations(),
          Ac.applications,
          BlockTimer::GlobalTimer().Counts(Timer::KSP_PRECONDITIONER) - count,
          ksp->GetConverged(),
          ksp->GetInitialRes(),
          ksp->GetFinalRes()};
}

// The vectors u and w of the norm of the GMRES stopping test, ||w - u|| for the guess x:
// u = B A x and w = B b with left preconditioning, else u = A x and w = b.
template <typename OperType, typename VecType>
void GmresResidualTerms(const GmresTestCase<OperType> &c, const OperType &A,
                        const VecType &b, const VecType &x, VecType &u, VecType &w)
{
  u.SetSize(b.Size());
  w.SetSize(b.Size());
  UseDevice(u, w);
  A.Mult(x, u);
  w = b;
  if (c.Left())
  {
    VecType t(b.Size());
    UseDevice(t);
    c.B->Mult(u, t);
    u = t;
    c.B->Mult(b, w);
  }
}

// The residual of x in the norm of the GMRES stopping test.
template <typename OperType, typename VecType>
double GmresResidual(const GmresTestCase<OperType> &c, const OperType &A, const VecType &b,
                     const VecType &x)
{
  VecType u, w;
  GmresResidualTerms(c, A, b, x, u, w);
  w -= u;
  return linalg::Norml2(Mpi::World(), w);
}

// Check the scaled initial guess of GMRES and FGMRES, with no preconditioner, a Jacobi
// preconditioner, and another one, for the system A x = b with the exact solution xs and a
// scalar c ≠ 1.
template <typename OperType, typename VecType, typename ScalarType>
void CheckGmresInitialGuess(const OperType &A, const Solver<OperType> &jacobi,
                            const Solver<OperType> &tridiagonal, const VecType &xs,
                            ScalarType c)
{
  MPI_Comm comm = Mpi::World();
  const int n = xs.Size();
  VecType b(n), x(n), y(n), x0(n), noise(n);
  UseDevice(b, x, y, x0, noise);
  A.Mult(xs, b);
  auto RelDiff = [&](const VecType &u, const VecType &v)
  {
    VecType d(n);
    UseDevice(d);
    d = u;
    d -= v;
    return linalg::Norml2(comm, d) / linalg::Norml2(comm, v);
  };

  constexpr auto left = PreconditionerSide::LEFT, right = PreconditionerSide::RIGHT;
  const std::vector<GmresTestCase<OperType>> tests = {
      {"GMRES, none", nullptr, left, false},
      {"GMRES, left Jacobi", &jacobi, left, false},
      {"GMRES, right Jacobi", &jacobi, right, false},
      {"GMRES, left tridiagonal", &tridiagonal, left, false},
      {"GMRES, right tridiagonal", &tridiagonal, right, false},
      {"FGMRES, Jacobi", &jacobi, right, true},
      {"FGMRES, tridiagonal", &tridiagonal, right, true}};
  for (const auto &test : tests)
  {
    INFO("Solver: " << test.name);

    // Without an initial guess, for the reference norm of the stopping test (the
    // residual of a zero initial guess) and the numbers of applications.
    x = 1.0;
    const auto zero = SolveGmres(test, A, b, x, false, 1.0e-8, 1000, 5);
    CHECK(zero.converged);

    // Solve for rhs from the initial guess x, which gives the iterations, final residual,
    // and solution of the solve without an initial guess (returned first).
    auto SolveAsWithoutGuess = [&](const VecType &rhs)
    {
      y = 1.0;
      const auto ref = SolveGmres(test, A, rhs, y, false, 1.0e-8, 1000, 5);
      const auto stats = SolveGmres(test, A, rhs, x, true, 1.0e-8, 1000, 5);
      CHECK(stats.converged);
      CHECK(stats.it == ref.it);
      CHECK(stats.final_res == ref.final_res);
      y -= x;
      CHECK(linalg::Norml2(comm, y) == 0.0);
      return std::pair(ref, stats);
    };

    // A multiple of the exact solution is scaled to it, without iterations, applying A
    // once and B twice with left preconditioning (to b - A x and b), else not at all; with
    // left preconditioning and |α| > 2 (here α = 4), the residual is computed again, with
    // one more application of each.
    for (const ScalarType s : {c, ScalarType(0.25)})
    {
      INFO("Multiple " << s << " of the exact solution");
      const bool again = test.Left() && std::abs(s) < 0.5;
      x = xs;
      x *= s;
      const auto stats = SolveGmres(test, A, b, x, true, 1.0e-10, 1000, 5);
      CHECK(stats.it == 0);
      CHECK(stats.converged);
      CHECK(stats.operator_applications == (again ? 2 : 1));
      CHECK(stats.preconditioner_applications == (test.Left() ? (again ? 3 : 2) : 0));
      CHECK_THAT(stats.initial_res, WithinRel(zero.initial_res, 1.0e-12));
      CHECK(stats.final_res <= 1.0e-12 * stats.initial_res);
      CHECK(RelDiff(x, xs) <= 1.0e-12);
    }

    // For a guess much smaller than the solution (|α| ~ 10⁶), the residual of the scaled
    // guess is formed from A x without left preconditioning, or computed again with it, so
    // its rounding error is not amplified by |α|: the solver returns its true residual.
    {
      INFO("Small guess");
      SetRandomTestVector(noise, 14);
      x = xs;
      x.Add(1.0e-4, noise);
      x *= 1.0e-6;
      const auto stats = SolveGmres(test, A, b, x, true, 1.0, 1000, 5);
      CHECK(stats.it == 0);
      CHECK_THAT(stats.final_res, WithinRel(GmresResidual(test, A, b, x), 1.0e-8));
    }

    // A zero initial guess takes the path without an initial guess, with the same numbers
    // of applications.
    {
      INFO("Zero initial guess");
      x = 0.0;
      const auto [ref, stats] = SolveAsWithoutGuess(b);
      CHECK(stats.operator_applications == ref.operator_applications);
      CHECK(stats.preconditioner_applications == ref.preconditioner_applications);
      CHECK(stats.initial_res == ref.initial_res);
    }

    // A guess orthogonal to the right-hand side in the norm of the stopping test is scaled
    // to zero (α = 0 exactly: the supports of u and w are disjoint), which gives the
    // results without an initial guess, with one more application of A and, with left
    // preconditioning, of B.
    {
      INFO("Guess orthogonal in the stopping norm");
      VecType b2(n);
      UseDevice(b2);
      SetRandomTestVector(b2, 10);
      SetRandomTestVector(x, 11);
      linalg::SetSubVector(b2, 0, n / 2, 0.0);
      linalg::SetSubVector(x, n / 4, n, 0.0);
      const auto [ref, stats] = SolveAsWithoutGuess(b2);
      CHECK(stats.operator_applications == ref.operator_applications + 1);
      CHECK(stats.preconditioner_applications ==
            ref.preconditioner_applications + (test.Left() ? 1 : 0));
      CHECK_THAT(stats.initial_res, WithinRel(ref.initial_res, 1.0e-14));
    }

    // With a zero right-hand side, a nonzero guess is scaled to the zero solution.
    {
      INFO("Zero right-hand side");
      VecType zero_rhs(n);
      UseDevice(zero_rhs);
      zero_rhs = 0.0;
      SetRandomTestVector(x, 12);
      const auto stats = SolveGmres(test, A, zero_rhs, x, true, 1.0e-8, 1000, 5);
      CHECK(stats.it == 0);
      CHECK(stats.converged);
      CHECK(stats.initial_res == 0.0);
      CHECK(stats.final_res == 0.0);
      CHECK(linalg::Norml2(comm, x) == 0.0);
    }

    // A guess so large that the inner products of its scaling overflow is replaced by a
    // zero guess (α = 0, not a non-finite α): the results without an initial guess.
    {
      INFO("Overflowing guess");
      VecType b2(n);
      UseDevice(b2);
      b2 = b;
      b2 *= 1.0e10;
      SetRandomTestVector(x, 13);
      x *= 1.0e300;
      SolveAsWithoutGuess(b2);
    }

    // Other initial guesses x0 = a xs + e noise: random ones, a large one (worse than zero
    // unscaled), and perturbations of multiples of the exact solution. With a relative
    // tolerance of 1, the solver returns the scaled guess without iterations, since its
    // residual is smaller than that of a zero guess. The scaling factor minimizes the
    // residual.
    const std::array<std::pair<ScalarType, double>, 6> guesses = {
        {{0.0, 1.0}, {0.0, 1.0}, {0.0, 100.0}, {1.0, 0.01}, {-3.0, 1.0}, {c, 0.1}}};
    for (int k = 0; k < 6; k++)
    {
      CAPTURE(k);
      SetRandomTestVector(noise, k);
      x0 = xs;
      x0 *= guesses[k].first;
      x0.Add(guesses[k].second, noise);
      VecType u, w;
      GmresResidualTerms(test, A, b, x0, u, w);
      const auto alpha = linalg::Dot(comm, w, u) / linalg::Dot(comm, u, u);
      const double res_zero = linalg::Norml2(comm, w);
      const double res_unscaled = GmresResidual(test, A, b, x0);
      if (k == 2)
      {
        CHECK(res_unscaled > 10.0 * res_zero);
      }

      x = x0;
      const auto stats = SolveGmres(test, A, b, x, true, 1.0, 1000, 5);
      CHECK(stats.it == 0);
      CHECK_THAT(stats.initial_res, WithinRel(zero.initial_res, 1.0e-12));
      CHECK_THAT(stats.initial_res, WithinRel(res_zero, 1.0e-12));
      CHECK(stats.final_res <= stats.initial_res);
      CHECK(stats.final_res <= res_unscaled);
      CHECK_THAT(stats.final_res, WithinRel(GmresResidual(test, A, b, x), 1.0e-10));
      y = x0;
      y *= alpha;
      CHECK(RelDiff(x, y) <= 1.0e-12);
      std::vector<ScalarType> factors = {0.99, 1.01};
      if constexpr (std::is_same_v<VecType, ComplexVector>)
      {
        factors.insert(factors.end(), {ScalarType(1.0, 0.01), ScalarType(1.0, -0.01)});
      }
      for (ScalarType t : factors)
      {
        y = x;
        y *= t;
        CHECK(stats.final_res <= GmresResidual(test, A, b, y));
      }

      // The solve from the guess converges, also with restarts.
      x = x0;
      const auto solve = SolveGmres(test, A, b, x, true, 1.0e-8, 1000, 5);
      CHECK(solve.converged);
      CHECK(GmresResidual(test, A, b, x) <= 1.01e-8 * res_zero);
    }
  }
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
    const auto *X = x.HostRead();  // The solver may have used x on the device
    CHECK_THAT(X[0], WithinAbs(1.0, 1.0e-14));
    CHECK_THAT(X[1], WithinAbs(-2.0, 1.0e-14));
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

// GMRES and FGMRES with each orthogonalization on real and complex nonsymmetric systems.
TEST_CASE("GmresSolver orthogonalization", "[iterative][Serial][Parallel]")
{
  constexpr int n = 200;
  MPI_Comm comm = Mpi::World();
  const auto orthog =
      GENERATE(Orthogonalization::MGS, Orthogonalization::CGS, Orthogonalization::CGS2);
  INFO("Orthogonalization: " << static_cast<int>(orthog));
  TridiagonalTestOperator T(n);
  SkewTridiagonalTestOperator K(n, 0.5);
  SumOperator A(T, K);

  SECTION("Real")
  {
    Vector b(n);
    for (int i = 0; i < n; i++)
    {
      b[i] = std::sin(1.0 + i + Mpi::Rank(comm));
    }
    TridiagonalTestSolver<Operator> pc;
    CheckGmresSolve<GmresSolver<Operator>, Operator>(A, nullptr, b, orthog);
    CheckGmresSolve<FgmresSolver<Operator>, Operator>(A, &pc, b, orthog);
  }

  SECTION("Complex")
  {
    // The imaginary part is symmetric, so the operator is neither symmetric nor Hermitian.
    TridiagonalTestOperator Ai(n, 0.5, 0.25, 3);
    ComplexWrapperOperator Ac(&A, &Ai);
    ComplexVector b(n);
    for (int i = 0; i < n; i++)
    {
      b.Real()[i] = std::sin(1.0 + i + Mpi::Rank(comm));
      b.Imag()[i] = std::cos(2.0 * i);
    }
    TridiagonalTestSolver<ComplexOperator> pc;
    CheckGmresSolve<GmresSolver<ComplexOperator>, ComplexOperator>(Ac, nullptr, b, orthog);
    CheckGmresSolve<FgmresSolver<ComplexOperator>, ComplexOperator>(Ac, &pc, b, orthog);
  }
}

// GMRES and FGMRES scale a nonzero initial guess to minimize the initial residual in the
// norm of their stopping test, so it is never worse than a zero initial guess, and take
// the path without an initial guess for a zero one.
TEST_CASE("GmresSolver scaled initial guess", "[iterative][Serial][Parallel][GPU]")
{
  constexpr int n = 200;
  MPI_Comm comm = Mpi::World();
  TridiagonalTestOperator A(n);

  SECTION("Real")
  {
    Vector xs(n);
    for (int i = 0; i < n; i++)
    {
      xs[i] = std::cos(0.5 * i) + 0.1 * Mpi::Rank(comm);
    }
    UseDevice(xs);
    JacobiSmoother<Operator> jacobi(comm);
    TridiagonalTestSolver<Operator> tridiagonal;
    jacobi.SetOperator(A);
    CheckGmresInitialGuess<Operator>(A, jacobi, tridiagonal, xs, -2.5);
  }

  SECTION("Complex")
  {
    // A complex symmetric operator with a complex diagonal.
    TridiagonalTestOperator S(n, 0.5, 0.25, 3);
    ComplexWrapperOperator Ac(&A, &S);
    ComplexVector xs(n);
    for (int i = 0; i < n; i++)
    {
      xs.Real()[i] = std::cos(0.5 * i) + 0.1 * Mpi::Rank(comm);
      xs.Imag()[i] = std::sin(0.3 * i);
    }
    UseDevice(xs);
    JacobiSmoother<ComplexOperator> jacobi(comm);
    TridiagonalTestSolver<ComplexOperator> tridiagonal;
    jacobi.SetOperator(Ac);
    CheckGmresInitialGuess<ComplexOperator>(Ac, jacobi, tridiagonal, xs,
                                            std::complex<double>(0.3, -1.7));
  }
}

// With a preconditioner whose error is relative to its argument, a good initial guess
// converges to a residual below that error relative to the right-hand side, also when it
// is scaled: with left preconditioning, the scaled initial residual is formed from the
// residual of the guess (α ≈ 1), or computed again (|α| > 2, a guess much smaller than the
// solution), not from B A x and B b, whose errors are relative to B b.
TEST_CASE("GmresSolver scaled initial guess with an inexact preconditioner",
          "[iterative][Serial][Parallel][GPU]")
{
  constexpr int n = 200;
  constexpr double rel_tol = 1.0e-10;
  MPI_Comm comm = Mpi::World();
  TridiagonalTestOperator A(n);

  // The residual is checked in the norm of the stopping test with the exact Jacobi
  // preconditioner.
  auto Check =
      [&]<typename OperType, typename VecType>(const OperType &Aop, const VecType &xs)
  {
    JacobiSmoother<OperType> jacobi(comm);
    RoundingJacobiSmoother<OperType> rounding(comm);
    jacobi.SetOperator(Aop);
    rounding.SetOperator(Aop);
    constexpr auto left = PreconditionerSide::LEFT, right = PreconditionerSide::RIGHT;
    for (const GmresTestCase<OperType> &test :
         {GmresTestCase<OperType>{"GMRES, left", &rounding, left, false},
          GmresTestCase<OperType>{"GMRES, right", &rounding, right, false},
          GmresTestCase<OperType>{"FGMRES", &rounding, right, true}})
    {
      INFO("Solver: " << test.name);
      const GmresTestCase<OperType> exact{test.name, &jacobi, test.side, test.flexible};
      VecType b(n), x(n), noise(n);
      UseDevice(b, x, noise);
      Aop.Mult(xs, b);
      SetRandomTestVector(noise, 1);
      for (double scale : {1.0, 1.0001, 0.2, 0.02})
      {
        CAPTURE(scale);
        x = xs;
        x.Add(1.0e-6, noise);
        x *= scale;
        const auto stats = SolveGmres(test, Aop, b, x, true, rel_tol, 100);
        CHECK(stats.converged);
        VecType zero(n);
        UseDevice(zero);
        zero = 0.0;
        CHECK(GmresResidual(exact, Aop, b, x) <=
              2.0 * rel_tol * GmresResidual(exact, Aop, b, zero));
      }
    }
  };

  SECTION("Real")
  {
    Vector xs(n);
    for (int i = 0; i < n; i++)
    {
      xs[i] = std::cos(0.5 * i) + 0.1 * Mpi::Rank(comm);
    }
    UseDevice(xs);
    Check(static_cast<const Operator &>(A), xs);
  }

  SECTION("Complex")
  {
    TridiagonalTestOperator S(n, 0.5, 0.0);
    ComplexWrapperOperator Ac(&A, &S);
    ComplexVector xs(n);
    for (int i = 0; i < n; i++)
    {
      xs.Real()[i] = std::cos(0.5 * i) + 0.1 * Mpi::Rank(comm);
      xs.Imag()[i] = std::sin(0.3 * i);
    }
    UseDevice(xs);
    Check(static_cast<const ComplexOperator &>(Ac), xs);
  }
}

}  // namespace palace
