// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "linalg/deflation.hpp"
#include "linalg/iterative.hpp"
#include "linalg/jacobi.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
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

// A diagonal preconditioner which amplifies the first local component by 1e20: GMRES with
// left preconditioning reduces the preconditioned residual by 1e-20 in one iteration while
// the true residual of the identity operator is not reduced (false convergence).
class AmplifyingTestSolver : public Solver<Operator>
{
public:
  void SetOperator(const Operator &op) override {}

  void Mult(const Vector &x, Vector &y) const override
  {
    y.SetSize(x.Size());
    const auto *X = x.HostRead();
    auto *Y = y.HostWrite();
    for (int i = 0; i < x.Size(); i++)
    {
      Y[i] = (i == 0) ? 1.0e20 * X[i] : X[i];
    }
  }
};

class IdentityTestOperator : public Operator
{
public:
  IdentityTestOperator(int n) : Operator(n) {}

  void Mult(const Vector &x, Vector &y) const override { y = x; }
};

// The Jacobi preconditioner of TridiagonalTestOperator (with its default stencil), with the
// local component k amplified by amp: for large amp, GMRES with left preconditioning
// reduces the preconditioned residual without reducing the true residual.
class AmplifiedJacobiTestSolver : public Solver<Operator>
{
  const double amp;
  const int k;

public:
  AmplifiedJacobiTestSolver(double amp, int k) : amp(amp), k(k) {}

  void SetOperator(const Operator &op) override {}

  void Mult(const Vector &x, Vector &y) const override
  {
    y.SetSize(x.Size());
    const auto *X = x.HostRead();
    auto *Y = y.HostWrite();
    for (int i = 0; i < x.Size(); i++)
    {
      Y[i] = ((i == k) ? amp : 1.0) * X[i] / (4.0 + i % 7);
    }
  }
};

// A diagonal operator (or preconditioner) with the given local diagonal on each process.
class DiagonalTestOperator : public Solver<Operator>
{
  const std::vector<double> d;

public:
  DiagonalTestOperator(std::vector<double> d) : d(std::move(d))
  {
    height = width = static_cast<int>(this->d.size());
  }

  void SetOperator(const Operator &op) override {}

  void Mult(const Vector &x, Vector &y) const override
  {
    y.SetSize(x.Size());
    const auto *X = x.HostRead();
    auto *Y = y.HostWrite();
    for (int i = 0; i < x.Size(); i++)
    {
      Y[i] = d[i] * X[i];
    }
  }
};

// The complex diagonal operator with entries 1 + 0.6 sin(g) + 0.1 i cos(g) on the global
// rows g, except four entries close to zero (and one of them negative) on the first rows:
// with the identity preconditioner, a cluster around 1 and four outliers which slow down
// GMRES.
class DeflationTestOperator : public ComplexOperator
{
  std::vector<std::complex<double>> d;

public:
  DeflationTestOperator(MPI_Comm comm, int n) : ComplexOperator(n), d(n)
  {
    int offset = 0;
    MPI_Exscan(&n, &offset, 1, MPI_INT, MPI_SUM, comm);
    if (Mpi::Root(comm))
    {
      offset = 0;
    }
    const std::complex<double> outliers[4] = {
        {1.0e-3, 0.0}, {2.0e-3, 1.0e-4}, {-1.5e-3, 0.0}, {1.0e-3, 1.0e-3}};
    for (int i = 0; i < n; i++)
    {
      const int g = offset + i;
      d[i] = (g < 4) ? outliers[g]
                     : std::complex<double>(1.0 + 0.6 * std::sin(g), 0.1 * std::cos(g));
    }
  }

  void Mult(const ComplexVector &x, ComplexVector &y) const override
  {
    const auto *XR = x.Real().HostRead(), *XI = x.Imag().HostRead();
    auto *YR = y.Real().HostWrite(), *YI = y.Imag().HostWrite();
    for (std::size_t i = 0; i < d.size(); i++)
    {
      const std::complex<double> v = d[i] * std::complex<double>(XR[i], XI[i]);
      YR[i] = v.real();
      YI[i] = v.imag();
    }
  }
};

class IdentityComplexTestSolver : public Solver<ComplexOperator>
{
public:
  void SetOperator(const ComplexOperator &op) override {}

  void Mult(const ComplexVector &x, ComplexVector &y) const override { y = x; }
};

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

// The solves after the first one with the same operator are deflated with harmonic Ritz
// vectors from the Krylov subspace of the first solve: the outliers of the preconditioned
// operator no longer slow down GMRES, and the solution keeps the requested accuracy.
TEST_CASE("ComplexKspSolver deflation", "[iterative][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  constexpr int n = 50;
  constexpr double tol = 1.0e-10;
  const auto orthog = GENERATE(Orthogonalization::MGS, Orthogonalization::CGS2);
  DeflationTestOperator A(comm, n);
  auto MakeGmres = [&]()
  {
    auto gmres = std::make_unique<GmresSolver<ComplexOperator>>(comm, 0);
    gmres->SetRelTol(tol);
    gmres->SetMaxIter(200);
    gmres->SetRestartDim(200);
    gmres->SetOrthogonalization(orthog);
    return gmres;
  };
  auto gmres = MakeGmres();
  auto *ref_gmres = gmres.get();
  ComplexKspSolver ksp(std::move(gmres), std::make_unique<IdentityComplexTestSolver>());
  ksp.SetOperators(A, A);
  ksp.SetDeflation(4);
  CHECK(ksp.GetDeflationSize() == 0);

  ComplexVector b1(n), b2(n), x1(n), x2(n), r(n);
  for (auto *v : {&b1, &b2, &x1, &x2, &r})
  {
    v->UseDevice(true);
  }
  linalg::SetRandom(comm, b1, 1);
  linalg::SetRandom(comm, b2, 2);
  auto RelRes = [&](const ComplexVector &b, const ComplexVector &x)
  {
    A.Mult(x, r);
    linalg::AXPBY(1.0, b, -1.0, r);
    return linalg::Norml2(comm, r) / linalg::Norml2(comm, b);
  };

  // The first solve provides the deflation basis.
  ksp.Mult(b1, x1);
  const int its1 = ksp.NumTotalMultIterations();
  CHECK(ksp.GetConverged());
  CHECK(ksp.GetDeflationSize() == 4);
  CHECK(RelRes(b1, x1) < 1.0e3 * tol);

  // Reference solve without deflation.
  int its_ref;
  {
    auto ref = MakeGmres();
    IdentityComplexTestSolver B;
    ref->SetOperator(A);
    ref->SetPreconditioner(B);
    ComplexVector x(n);
    x.UseDevice(true);
    ref->Mult(b2, x);
    its_ref = ref->GetNumIterations();
  }
  ksp.Mult(b2, x2);
  const int its2 = ksp.NumTotalMultIterations() - its1;
  CHECK(ksp.GetConverged());
  CHECK(RelRes(b2, x2) < 1.0e3 * tol);
  CHECK(its2 + 10 <= its_ref);
  CHECK(ref_gmres->GetNumIterations() == its2);

  // New operators discard the basis, which the next solve computes again.
  ksp.SetOperators(A, A);
  CHECK(ksp.GetDeflationSize() == 0);
  ksp.Mult(b2, x2);
  CHECK(ksp.GetDeflationSize() == 4);
  CHECK(ksp.NumTotalMultIterations() - its1 - its2 == its_ref);
}

// An initial guess whose true residual is at least half that of zero (the bound of the
// adaptive sweep HDM solves) is discarded: the solve is then the solve from zero, bitwise.
// A better guess is kept, and its solve is not changed by the check.
TEST_CASE("GmresSolver initial guess check", "[iterative][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  constexpr int n = 50;
  DeflationTestOperator A(comm, n);
  IdentityComplexTestSolver B;
  auto Solve = [&](const ComplexVector &b, ComplexVector &x, bool guess, double tol)
  {
    GmresSolver<ComplexOperator> gmres(comm, 0);
    gmres.SetRelTol(1.0e-10);
    gmres.SetMaxIter(200);
    gmres.SetRestartDim(200);
    gmres.SetOperator(A);
    gmres.SetPreconditioner(B);
    gmres.SetInitialGuess(guess);
    gmres.SetInitialGuessCheck(tol);
    gmres.Mult(b, x);
    CHECK(gmres.GetConverged());
    return gmres.GetNumIterations();
  };
  ComplexVector b(n), x0(n), xz(n), x1(n), x2(n);
  for (auto *v : {&b, &x0, &xz, &x1, &x2})
  {
    v->UseDevice(true);
  }
  linalg::SetRandom(comm, b, 1);
  xz = 0.0;
  const int its_zero = Solve(b, xz, false, 0.0);
  auto Diff = [&](const ComplexVector &x, const ComplexVector &y)
  {
    x2 = x;
    x2 -= y;
    return linalg::Norml2(comm, x2);
  };

  SECTION("Guess no better than zero")
  {
    linalg::SetRandom(comm, x0, 2);
    x0 *= 10.0;
    x1 = x0;
    CHECK(Solve(b, x1, true, 0.5) == its_zero);
    CHECK(Diff(x1, xz) == 0.0);
  }

  SECTION("Guess with 0.7 times the residual of zero")
  {
    x0 = xz;
    x0 *= 0.3;
    x1 = x0;
    CHECK(Solve(b, x1, true, 0.5) == its_zero);
    CHECK(Diff(x1, xz) == 0.0);
  }

  SECTION("Guess better than zero")
  {
    x0 = xz;
    x0 *= 0.9;
    x1 = x0;
    const int its1 = Solve(b, x1, true, 0.5);
    ComplexVector x3(n);
    x3.UseDevice(true);
    x3 = x0;
    const int its2 = Solve(b, x3, true, 0.0);
    CHECK(its1 == its2);
    CHECK(Diff(x1, x3) == 0.0);
    CHECK(its1 <= its_zero);
  }
}

// The harmonic Ritz analysis of a solve predicts the iterations of GMRES with and without
// the outliers of the preconditioned operator, and the automatic deflation of the solves
// after the first one never needs more iterations than the solve without deflation.
TEST_CASE("ComplexKspSolver automatic deflation", "[iterative][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  constexpr int n = 50;
  constexpr double tol = 1.0e-10;
  DeflationTestOperator A(comm, n);
  IdentityComplexTestSolver B;
  auto MakeGmres = [&]()
  {
    auto gmres = std::make_unique<GmresSolver<ComplexOperator>>(comm, 0);
    gmres->SetRelTol(tol);
    gmres->SetMaxIter(200);
    gmres->SetRestartDim(200);
    gmres->SetOrthogonalization(Orthogonalization::CGS2);
    return gmres;
  };
  ComplexVector b1(n), b2(n), x(n), r(n);
  for (auto *v : {&b1, &b2, &x, &r})
  {
    v->UseDevice(true);
  }
  linalg::SetRandom(comm, b1, 1);
  linalg::SetRandom(comm, b2, 2);

  // Model: GMRES for the projected operator reproduces the iterations of the solve, and
  // predicts many fewer with the four outliers deflated.
  {
    auto gmres = MakeGmres();
    gmres->SetOperator(A);
    gmres->SetPreconditioner(B);
    x = 0.0;
    gmres->Mult(b1, x);
    const int its = gmres->GetNumIterations();
    const std::vector<ComplexVector> *V;
    const std::complex<double> *H;
    int ld;
    const int m = gmres->GetArnoldiRelation(V, H, ld);
    REQUIRE(m >= its - 2);
    DeflationAnalysis analysis;
    analysis.Compute(comm, H, ld, m, 8, {4, 8}, tol);
    CHECK(std::abs(analysis.its_plain - its) <= 2);
    REQUIRE(analysis.its_deflated.size() == 2);
    CHECK(analysis.its_deflated[0] + 10 <= analysis.its_plain);
    CHECK(std::abs(analysis.theta[0]) < 0.01);
    CHECK(std::abs(analysis.theta[3]) < 0.01);
    CHECK(std::abs(analysis.theta[4]) > 0.1);

    // A basis of 8 vectors of which 4 are active deflates the outliers, and activating
    // fewer vectors after more gives the same preconditioner.
    DeflatedPreconditioner defl(comm);
    REQUIRE(defl.Build(A, B, *V, analysis, 8) == 8);
    CHECK(defl.Size() == 0);
    CHECK(defl.Prefix(4) == 4);
    auto Solve = [&](const Solver<ComplexOperator> &pc)
    {
      auto solver = MakeGmres();
      solver->SetOperator(A);
      solver->SetPreconditioner(pc);
      x = 0.0;
      solver->Mult(b2, x);
      CHECK(solver->GetConverged());
      return solver->GetNumIterations();
    };
    const int its_plain = Solve(B);
    REQUIRE(defl.Activate(4) == 4);
    const int its4 = Solve(defl);
    CHECK(its4 + 10 <= its_plain);
    REQUIRE(defl.Activate(8) == 8);
    Solve(defl);
    REQUIRE(defl.Activate(4) == 4);
    CHECK(Solve(defl) == its4);
  }

  // Automatic deflation: the second solve is deflated unless the cost model rejects it, and
  // never takes more iterations than without deflation.
  int its_ref;
  {
    auto ref = MakeGmres();
    ref->SetOperator(A);
    ref->SetPreconditioner(B);
    x = 0.0;
    ref->Mult(b2, x);
    its_ref = ref->GetNumIterations();
  }
  ComplexKspSolver ksp(MakeGmres(), std::make_unique<IdentityComplexTestSolver>());
  ksp.SetOperators(A, A);
  ksp.SetDeflation(-1);
  ksp.Mult(b1, x);
  const int its1 = ksp.NumTotalMultIterations();
  ksp.Mult(b2, x);
  const int its2 = ksp.NumTotalMultIterations() - its1;
  CHECK(ksp.GetConverged());
  A.Mult(x, r);
  linalg::AXPBY(1.0, b2, -1.0, r);
  CHECK(linalg::Norml2(comm, r) / linalg::Norml2(comm, b2) < 1.0e3 * tol);
  CHECK(its2 <= its_ref);
  CHECK((ksp.GetDeflationSize() == 0 || its2 + 10 <= its_ref));
}

// The automatic deflation tries larger sets of vectors while the deflation pays, and goes
// back (or stops) when a set does not pay; a set which did not pay is not tried again.
TEST_CASE("Automatic deflation levels", "[iterative][Serial]")
{
  using Action = DeflationLevels::Action;
  DeflationLevels levels;
  levels.Reset({4});
  CHECK(levels.Size() == 4);
  CHECK(levels.Verify(true) == Action::KEEP);
  CHECK(levels.Verify(false) == Action::STOP);

  levels.Reset({4, 4, 10, 20});
  CHECK(levels.Verify(true) == Action::GROW);
  CHECK(levels.Size() == 10);
  CHECK(levels.PreviousSize() == 4);
  CHECK(levels.Trying());
  CHECK(levels.Compare(true) == Action::GROW);
  CHECK(levels.Size() == 20);
  CHECK(levels.Compare(false) == Action::BACK);
  CHECK(levels.Size() == 10);
  CHECK(!levels.Trying());
  CHECK(levels.Verify(true) == Action::KEEP);
  CHECK(levels.Verify(false) == Action::BACK);
  CHECK(levels.Size() == 4);
  CHECK(levels.Verify(true) == Action::KEEP);
  CHECK(levels.Verify(false) == Action::STOP);
}

TEST_CASE("GmresSolver false convergence", "[iterative][Serial][Parallel]")
{
  IdentityTestOperator A(3);
  AmplifyingTestSolver B;
  GmresSolver<Operator> ksp(Mpi::World(), 0);
  ksp.SetOperator(A);
  ksp.SetPreconditioner(B);
  ksp.SetRelTol(1.0e-12);
  ksp.SetMaxIter(10);
  Vector b(3), x(3), r(3);
  b = 1.0;
  x = 0.0;

  SECTION("Unchecked")
  {
    // Without the check (the default), the garbage solution is reported as converged.
    ksp.Mult(b, x);
    CHECK(ksp.GetNumIterations() == 1);
    CHECK(ksp.GetConverged());
    CHECK(!ksp.GetFalseConverged());
    A.Mult(x, r);
    r -= b;
    CHECK(linalg::Norml2(Mpi::World(), r) > 0.5 * linalg::Norml2(Mpi::World(), b));
  }

  SECTION("Checked")
  {
    // The check fails after the first iteration. The continuation with right
    // preconditioning reduces the true residual only far from the iterate (it changes the
    // amplified component), so the solve is falsely converged and returns the iterate.
    ksp.SetTrueResidualCheck(0.5);
    ksp.Mult(b, x);
    CHECK(ksp.GetNumIterations() > 1);
    CHECK(!ksp.GetConverged());
    CHECK(ksp.GetFalseConverged());
    A.Mult(x, r);
    r -= b;
    CHECK(linalg::Norml2(Mpi::World(), r) > 0.5 * linalg::Norml2(Mpi::World(), b));
  }

  SECTION("Checked with a singular preconditioner")
  {
    // A preconditioner which annihilates a component of the residual: the right-
    // preconditioned continuation cannot reduce the true residual and the solve is falsely
    // converged.
    DiagonalTestOperator B0({1.0, 1.0, 0.0});
    ksp.SetPreconditioner(B0);
    ksp.SetTrueResidualCheck(0.5);
    ksp.Mult(b, x);
    CHECK(!ksp.GetConverged());
    CHECK(ksp.GetFalseConverged());
  }

  SECTION("Stiff error of the initial guess")
  {
    // A soft (resonant) component amplified 10⁶ times and a stiff one, with a good
    // preconditioner: the initial guess has an error in the stiff component which the
    // preconditioned residual does not see (1.5e-12 of norm(B b) after the first iteration)
    // but whose true residual is about norm(b). GMRES with left preconditioning converges
    // in one iteration, the check fails and the solve from zero passes it.
    DiagonalTestOperator A2({1.0e-6, 1.0, 1.0, 1.0e6});
    DiagonalTestOperator B2({1.1e6, 0.95, 1.02, 0.9e-6});
    GmresSolver<Operator> ksp2(Mpi::World(), 0);
    ksp2.SetOperator(A2);
    ksp2.SetPreconditioner(B2);
    ksp2.SetRelTol(1.0e-10);
    ksp2.SetMaxIter(10);
    ksp2.SetInitialGuess(true);
    Vector b2(4), x2(4), r2(4);
    b2 = 1.0;
    b2[3] = 0.0;
    x2[0] = 1.0e6 + 1.0e3;  // Exact solution (10⁶, 1, 1, 0) with soft and stiff errors
    x2[1] = 1.0;
    x2[2] = 1.0;
    x2[3] = 1.0e-5;
    Vector x3(x2);
    ksp2.Mult(b2, x2);
    CHECK(ksp2.GetNumIterations() == 1);
    CHECK(ksp2.GetConverged());
    A2.Mult(x2, r2);
    r2 -= b2;
    CHECK(linalg::Norml2(Mpi::World(), r2) > 0.5 * linalg::Norml2(Mpi::World(), b2));

    ksp2.SetTrueResidualCheck(0.5);
    ksp2.Mult(b2, x3);
    CHECK(ksp2.GetNumIterations() > 1);
    CHECK(ksp2.GetConverged());
    CHECK(!ksp2.GetFalseConverged());
    A2.Mult(x3, r2);
    r2 -= b2;
    CHECK(linalg::Norml2(Mpi::World(), r2) <= 1.0e-6 * linalg::Norml2(Mpi::World(), b2));
    CHECK_THAT(x3[0], WithinRel(1.0e6, 1.0e-12));
    CHECK_THAT(x3[3], WithinAbs(0.0, 1.0e-12));
    const std::vector<Vector> *V;
    const double *H;
    int ld;
    CHECK(ksp2.GetArnoldiRelation(V, H, ld) > 0);
  }

  SECTION("Accurate solution with a stiff residual")
  {
    // From zero, with a stiff component of the right-hand side whose residual the
    // preconditioned residual does not see: GMRES with left preconditioning converges in
    // one iteration to a solution accurate to 3e-10 relative, with a true residual of 0.8
    // norm(b). One iteration with right preconditioning removes it, changing the solution
    // by 3e-10 relative: the solution of the left-preconditioned solve is kept.
    DiagonalTestOperator A4({1.0e-6, 1.0, 1.0, 1.0e6});
    DiagonalTestOperator B4({1.0e6, 1.0, 1.0, 2.0e-7});
    GmresSolver<Operator> ksp4(Mpi::World(), 0);
    ksp4.SetOperator(A4);
    ksp4.SetPreconditioner(B4);
    ksp4.SetRelTol(1.0e-10);
    ksp4.SetMaxIter(10);
    ksp4.SetTrueResidualCheck(0.5);
    Vector b4(4), x4(4), r4(4);
    b4 = 1.0;
    b4[3] = 400.0;  // Exact solution (10⁶, 1, 1, 4e-4)
    x4 = 0.0;
    ksp4.Mult(b4, x4);
    CHECK(ksp4.GetNumIterations() == 2);
    CHECK(ksp4.GetConverged());
    CHECK(!ksp4.GetFalseConverged());
    CHECK_THAT(x4[0], WithinRel(1.0e6, 1.0e-12));
    CHECK_THAT(x4[3], WithinAbs(8.0e-5, 1.0e-12));
    A4.Mult(x4, r4);
    r4 -= b4;
    CHECK(linalg::Norml2(Mpi::World(), r4) > 0.5 * linalg::Norml2(Mpi::World(), b4));
  }

  SECTION("Amplified Jacobi preconditioner")
  {
    // A well-conditioned system with a Jacobi preconditioner which amplifies one component:
    // GMRES with left preconditioning converges at once without reducing the true residual.
    // The continuation with right preconditioning reduces it by moving the solution by
    // about its norm, so the solve is falsely converged; without amplification, the check
    // passes.
    const int n = 40;
    TridiagonalTestOperator A5(n);
    Vector xe(n), b5(n), x5(n);
    {
      auto *XE = xe.HostWrite();
      for (int i = 0; i < n; i++)
      {
        XE[i] = 1.0 + 0.1 * i + 0.01 * (Mpi::Rank(Mpi::World()) + 1);
      }
    }
    A5.Mult(xe, b5);
    for (double amp : {1.0, 1.0e14, 1.0e20})
    {
      AmplifiedJacobiTestSolver B5(amp, 3);
      GmresSolver<Operator> ksp5(Mpi::World(), 0);
      ksp5.SetOperator(A5);
      ksp5.SetPreconditioner(B5);
      ksp5.SetRelTol(1.0e-12);
      ksp5.SetMaxIter(200);
      ksp5.SetStagnation(100, 1.0e-3);
      ksp5.SetTrueResidualCheck(0.5);
      x5 = 0.0;
      ksp5.Mult(b5, x5);
      if (amp == 1.0)
      {
        CHECK(ksp5.GetConverged());
        CHECK(!ksp5.GetFalseConverged());
        x5 -= xe;
        CHECK(linalg::Norml2(Mpi::World(), x5) <=
              1.0e-9 * linalg::Norml2(Mpi::World(), xe));
      }
      else
      {
        CHECK(!ksp5.GetConverged());
        CHECK(ksp5.GetFalseConverged());
      }
    }
  }

  SECTION("Arnoldi relation of a continued solve")
  {
    // A left-preconditioned solve provides the Arnoldi relation of B A (for the deflation
    // of later solves), its right-preconditioned continuation after a failed check none: it
    // would be the one of A B, whose Krylov vectors do not approximate eigenvectors of B A.
    // The stiff component of the right-hand side (B A = 0.02, apart from the others, which
    // GMRES reduces little) is invisible in the preconditioned residual, the others need
    // several iterations.
    DiagonalTestOperator A3({1.0e-6, 1.0, 2.0, 3.0, 4.0, 5.0, 1.0e6});
    DiagonalTestOperator B3({1.1e6, 1.0, 0.55, 0.3, 0.3, 0.16, 2.0e-8});
    GmresSolver<Operator> ksp3(Mpi::World(), 0);
    ksp3.SetOperator(A3);
    ksp3.SetPreconditioner(B3);
    ksp3.SetRelTol(1.0e-10);
    ksp3.SetMaxIter(20);
    Vector b3(7), x3(7);
    b3 = 1.0;
    b3[6] = 400.0;
    const std::vector<Vector> *V;
    const double *H;
    int ld;
    x3 = 0.0;
    ksp3.Mult(b3, x3);
    const int its = ksp3.GetNumIterations();
    CHECK(its >= 3);
    CHECK(ksp3.GetArnoldiRelation(V, H, ld) > 0);

    ksp3.SetTrueResidualCheck(0.5);
    x3 = 0.0;
    ksp3.Mult(b3, x3);
    CHECK(ksp3.GetNumIterations() > its);
    CHECK(ksp3.GetConverged());
    CHECK(!ksp3.GetFalseConverged());
    CHECK(ksp3.GetArnoldiRelation(V, H, ld) == 0);
  }

  SECTION("Stable preconditioner")
  {
    // The check passes for a converged solve with a stable preconditioner.
    SpdTestOperator A2;
    TridiagonalTestSolver<Operator> B2;
    GmresSolver<Operator> ksp2(Mpi::World(), 0);
    ksp2.SetOperator(A2);
    ksp2.SetPreconditioner(B2);
    ksp2.SetRelTol(1.0e-12);
    ksp2.SetMaxIter(10);
    ksp2.SetTrueResidualCheck(0.5);
    Vector b2(2), x2(2);
    b2 = 1.0;
    x2 = 0.0;
    ksp2.Mult(b2, x2);
    CHECK(ksp2.GetConverged());
    CHECK(!ksp2.GetFalseConverged());
  }
}

TEST_CASE("KspSolver allowed false convergence", "[iterative][Serial][Parallel]")
{
  // An allowed false convergence returns the invalid solution to the caller, which checks
  // GetFalseConverged (the HDM solve of an adaptive sweep then solves again from zero). The
  // preconditioner annihilates a component of the residual, which the right-preconditioned
  // continuation of the solve cannot reduce.
  IdentityTestOperator A(3);
  auto gmres = std::make_unique<GmresSolver<Operator>>(Mpi::World(), 0);
  gmres->SetRelTol(1.0e-12);
  gmres->SetMaxIter(10);
  gmres->SetTrueResidualCheck(0.5);
  KspSolver ksp(std::move(gmres),
                std::make_unique<DiagonalTestOperator>(std::vector<double>{1.0, 1.0, 0.0}),
                "Test");
  ksp.SetOperators(A, A);
  ksp.SetAllowFalseConvergence(true);
  Vector b(3), x(3);
  b = 1.0;
  x = 0.0;
  ksp.Mult(b, x, false);
  CHECK(ksp.GetFalseConverged());
  CHECK(!ksp.GetConverged());

  // A stable solve after it is not affected.
  x = 0.0;
  ksp.SetAllowFalseConvergence(false);
  auto gmres2 = std::make_unique<GmresSolver<Operator>>(Mpi::World(), 0);
  gmres2->SetRelTol(1.0e-12);
  gmres2->SetMaxIter(10);
  gmres2->SetTrueResidualCheck(0.5);
  KspSolver ksp2(std::move(gmres2), nullptr, "Test");
  ksp2.SetOperators(A, A);
  ksp2.Mult(b, x);
  CHECK(ksp2.GetConverged());
  CHECK(!ksp2.GetFalseConverged());
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

// With the Jacobi stopping norm, CG stops at the first iterate whose residual r satisfies
// √(D⁻¹ r, r) <= tol √(D⁻¹ b, b), whatever the preconditioner (here none, or Jacobi), with
// or without an initial guess.
TEST_CASE("CgSolver Jacobi stopping norm", "[iterative][Serial][Parallel]")
{
  constexpr int n = 100, max_it = 1000;
  constexpr double tol = 1.0e-6;
  MPI_Comm comm = Mpi::World();
  TridiagonalTestOperator A(n);
  ComplexWrapperOperator Ac(&A, nullptr);
  Vector dinv;
  A.AssembleDiagonal(dinv);
  dinv.Reciprocal();

  auto Check =
      [&]<OperatorType OperType, typename VecType>(const OperType &op, const VecType &b)
  {
    auto JacobiNorm = [&](const VecType &v)
    {
      double s = 0.0;
      for (int i = 0; i < n; i++)
      {
        if constexpr (std::is_same_v<VecType, ComplexVector>)
        {
          s += (v.Real()[i] * v.Real()[i] + v.Imag()[i] * v.Imag()[i]) * dinv[i];
        }
        else
        {
          s += v[i] * v[i] * dinv[i];
        }
      }
      Mpi::GlobalSum(1, &s, comm);
      return std::sqrt(s);
    };
    auto ResidualNorm = [&](const VecType &x)
    {
      VecType r(n);
      op.Mult(x, r);
      r -= b;
      return JacobiNorm(r);
    };
    struct Result
    {
      int its;
      double initial, final;
      bool converged;
    };
    auto Solve = [&](const Solver<OperType> *B, bool jacobi_norm, int its, double rel_tol,
                     VecType &x, bool guess = false)
    {
      CgSolver<OperType> ksp(comm, 0);
      ksp.SetOperator(op);
      if (B)
      {
        ksp.SetPreconditioner(*B);
      }
      ksp.SetJacobiStoppingNorm(jacobi_norm ? &dinv : nullptr);
      ksp.SetRelTol(rel_tol);
      ksp.SetAbsTol(0.0);
      ksp.SetMaxIter(its);
      ksp.SetInitialGuess(guess);
      ksp.Mult(b, x);
      return Result{ksp.GetNumIterations(), ksp.GetInitialRes(), ksp.GetFinalRes(),
                    ksp.GetConverged()};
    };
    const double b_norm = JacobiNorm(b);
    VecType x(n), y(n);
    x = 0.0;
    const auto res = Solve(nullptr, true, max_it, tol, x);
    CHECK(res.converged);
    CHECK_THAT(res.initial, WithinRel(b_norm, 1.0e-12));
    CHECK(res.final <= tol * b_norm);
    CHECK_THAT(res.final, WithinRel(ResidualNorm(x), 1.0e-6));

    // The previous iterate does not satisfy the stopping test.
    REQUIRE(res.its > 1);
    y = 0.0;
    Solve(nullptr, true, res.its - 1, 0.0, y);
    CHECK(ResidualNorm(y) > tol * b_norm);

    // The default norm differs.
    y = 0.0;
    CHECK(std::abs(Solve(nullptr, false, max_it, tol, y).initial - b_norm) > 0.1 * b_norm);

    // With an initial guess, the tolerance stays relative to the norm of b.
    for (int i = 0; i < n; i++)
    {
      if constexpr (std::is_same_v<VecType, ComplexVector>)
      {
        y.Real()[i] = 0.1 * std::cos(i);
        y.Imag()[i] = -0.2 * std::sin(i);
      }
      else
      {
        y[i] = 0.1 * std::cos(i);
      }
    }
    const double r0_norm = ResidualNorm(y);
    REQUIRE(std::abs(r0_norm - b_norm) > 0.01 * b_norm);
    const auto res_guess = Solve(nullptr, true, max_it, tol, y, true);
    CHECK(res_guess.converged);
    CHECK_THAT(res_guess.initial, WithinRel(b_norm, 1.0e-12));
    CHECK(res_guess.final <= tol * b_norm);
    CHECK_THAT(res_guess.final, WithinRel(ResidualNorm(y), 1.0e-6));

    // With Jacobi preconditioning, the stopping test does not change.
    JacobiSmoother<OperType> jacobi(comm);
    jacobi.SetOperator(op);
    x = 0.0;
    y = 0.0;
    const auto res_jacobi = Solve(&jacobi, true, max_it, tol, x);
    const auto res_ref = Solve(&jacobi, false, max_it, tol, y);
    CHECK(res_jacobi.its == res_ref.its);
    CHECK_THAT(res_jacobi.initial, WithinRel(res_ref.initial, 1.0e-12));
    CHECK_THAT(res_jacobi.final, WithinRel(res_ref.final, 1.0e-10));
  };

  // The same right-hand side on all processes, so that a missing reduction of the norm
  // fails the checks rather than letting the processes disagree on convergence.
  SECTION("Real")
  {
    Vector b(n);
    for (int i = 0; i < n; i++)
    {
      b[i] = std::sin(1.0 + i);
    }
    Check(static_cast<const Operator &>(A), b);
  }

  SECTION("Complex")
  {
    ComplexVector b(n);
    for (int i = 0; i < n; i++)
    {
      b.Real()[i] = std::sin(1.0 + i);
      b.Imag()[i] = std::cos(2.0 * i);
    }
    Check(static_cast<const ComplexOperator &>(Ac), b);
  }
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

// The upper bidiagonal stencil [a + r cos(θᵢ + t), s] on the local rows of each process,
// with θᵢ = 2π frac(i φ) for the golden ratio φ: with the imaginary part r sin(θᵢ) (a = 0,
// s = 0, t = -π/2), the eigenvalues fill the circle of radius r around a, and GMRES
// converges at the rate r.
class BidiagonalTestOperator : public Operator
{
  const double a, r, s, t;

public:
  BidiagonalTestOperator(int n, double a, double r, double s, double t = 0.0)
    : Operator(n), a(a), r(r), s(s), t(t)
  {
  }

  void AssembleDiagonal(Vector &diag) const override
  {
    diag.SetSize(height);
    auto *D = diag.HostWrite();
    for (int i = 0; i < height; i++)
    {
      D[i] = a + r * std::cos(2.0 * M_PI * std::fmod(0.6180339887498949 * i, 1.0) + t);
    }
  }

  void Mult(const Vector &x, Vector &y) const override
  {
    y.SetSize(height);
    const auto *X = x.HostRead();
    auto *Y = y.HostWrite();
    for (int i = 0; i < height; i++)
    {
      const double d =
          a + r * std::cos(2.0 * M_PI * std::fmod(0.6180339887498949 * i, 1.0) + t);
      Y[i] = d * X[i] + ((i < height - 1) ? s * X[i + 1] : 0.0);
    }
  }
};

// Solve with MGS, CGS and CGS2 (delayed reorthogonalization, one reduction per iteration):
// CGS2 needs the same number of iterations as MGS and reaches the same solution.
template <typename OperType>
void CheckGmresOrthogonalization(const OperType &A, const Solver<OperType> *B,
                                 bool flexible, PreconditionerSide side, int restart,
                                 double tol)
{
  using VecType =
      std::conditional_t<std::is_same_v<OperType, ComplexOperator>, ComplexVector, Vector>;
  MPI_Comm comm = Mpi::World();
  const int n = A.Height();
  VecType b(n), r(n), Br(n), Bb(n);
  b.UseDevice(true);
  linalg::SetRandom(comm, b, 7);
  const std::array<Orthogonalization, 3> types = {
      Orthogonalization::MGS, Orthogonalization::CGS, Orthogonalization::CGS2};
  std::array<int, 3> its = {};
  std::array<VecType, 3> xs;
  for (int k = 0; k < 3; k++)
  {
    CAPTURE(k, flexible, restart);
    std::unique_ptr<GmresSolver<OperType>> ksp;
    if (flexible)
    {
      ksp = std::make_unique<FgmresSolver<OperType>>(comm, 0);
    }
    else
    {
      ksp = std::make_unique<GmresSolver<OperType>>(comm, 0);
      ksp->SetPreconditionerSide(side);
    }
    ksp->SetOperator(A);
    if (B)
    {
      ksp->SetPreconditioner(*B);
    }
    ksp->SetOrthogonalization(types[k]);
    ksp->SetRelTol(tol);
    ksp->SetMaxIter(1000);
    if (restart > 0)
    {
      ksp->SetRestartDim(restart);
    }
    xs[k].SetSize(n);
    xs[k].UseDevice(true);
    xs[k] = 0.0;
    ksp->Mult(b, xs[k]);
    CHECK(ksp->GetConverged());
    its[k] = ksp->GetNumIterations();

    // The true (preconditioned for left preconditioning) residual.
    A.Mult(xs[k], r);
    r -= b;
    const bool left = (B && !flexible && side == PreconditionerSide::LEFT);
    if (left)
    {
      B->Mult(r, Br);
      B->Mult(b, Bb);
    }
    const double res =
        linalg::Norml2(comm, left ? Br : r) / linalg::Norml2(comm, left ? Bb : b);
    if (types[k] != Orthogonalization::CGS)  // CGS loses orthogonality at these tolerances
    {
      CHECK(res <= 10.0 * tol);
    }
  }
  CAPTURE(its[0], its[1], its[2]);
  CHECK(its[0] > 5);
  CHECK(std::abs(its[2] - its[0]) <= 1);
  VecType d(n);
  d = xs[2];
  d -= xs[0];
  CHECK(linalg::Norml2(comm, d) <= 1.0e-6 * linalg::Norml2(comm, xs[0]));
}

TEST_CASE("GmresSolver orthogonalization", "[iterative][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const int n = 400 + 37 * Mpi::Rank(comm);
  BidiagonalTestOperator Ar(n, 1.0, 0.8, 0.3), Ai(n, 0.0, 0.8, 0.0, -0.5 * M_PI);
  ComplexWrapperOperator A(&Ar, &Ai);
  JacobiSmoother<Operator> jacobi(comm);
  jacobi.SetOperator(Ar);
  TridiagonalTestSolver<ComplexOperator> ctridiagonal;
  for (double tol : {1.0e-10, 1.0e-12})
  {
    CAPTURE(tol);
    CheckGmresOrthogonalization<Operator>(Ar, nullptr, false, PreconditionerSide::LEFT, 0,
                                          tol);
    CheckGmresOrthogonalization<Operator>(Ar, &jacobi, false, PreconditionerSide::LEFT, 0,
                                          tol);
    CheckGmresOrthogonalization<ComplexOperator>(A, nullptr, false,
                                                 PreconditionerSide::LEFT, 0, tol);
    CheckGmresOrthogonalization<ComplexOperator>(A, &ctridiagonal, false,
                                                 PreconditionerSide::LEFT, 0, tol);
    CheckGmresOrthogonalization<ComplexOperator>(A, &ctridiagonal, false,
                                                 PreconditionerSide::RIGHT, 0, tol);
    CheckGmresOrthogonalization<ComplexOperator>(A, &ctridiagonal, true,
                                                 PreconditionerSide::RIGHT, 0, tol);
    CheckGmresOrthogonalization<ComplexOperator>(A, &ctridiagonal, false,
                                                 PreconditionerSide::LEFT, 25, tol);
    CheckGmresOrthogonalization<ComplexOperator>(A, &ctridiagonal, true,
                                                 PreconditionerSide::RIGHT, 25, tol);
  }
}

}  // namespace palace
