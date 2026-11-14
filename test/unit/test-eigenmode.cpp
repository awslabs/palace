// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <complex>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include "linalg/eps.hpp"
#include "linalg/iterative.hpp"
#include "linalg/jacobi.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/slepc.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

namespace palace
{

#if defined(PALACE_WITH_SLEPC)

namespace
{

// The symmetric tridiagonal stencil [-1, 2 + i / n, -1] minus s times the identity (real,
// as a complex operator) on a single process: distinct eigenvalues in (0, 5).
class TridiagonalEigenTestOperator : public ComplexOperator
{
private:
  const double s;

public:
  TridiagonalEigenTestOperator(int n, double s = 0.0) : ComplexOperator(n), s(s) {}

  bool IsReal() const override { return true; }

  void AssembleDiagonal(ComplexVector &diag) const override
  {
    diag.SetSize(height);
    diag = 0.0;
    auto *D = diag.Real().HostWrite();
    for (int i = 0; i < height; i++)
    {
      D[i] = 2.0 + static_cast<double>(i) / height - s;
    }
  }

  void Mult(const ComplexVector &x, ComplexVector &y) const override
  {
    y.SetSize(height);
    const auto *XR = x.Real().HostRead();
    const auto *XI = x.Imag().HostRead();
    auto *YR = y.Real().HostWrite();
    auto *YI = y.Imag().HostWrite();
    for (int i = 0; i < height; i++)
    {
      const double d = 2.0 + static_cast<double>(i) / height - s;
      YR[i] =
          d * XR[i] - ((i > 0) ? XR[i - 1] : 0.0) - ((i < height - 1) ? XR[i + 1] : 0.0);
      YI[i] =
          d * XI[i] - ((i > 0) ? XI[i - 1] : 0.0) - ((i < height - 1) ? XI[i + 1] : 0.0);
    }
  }

  void AddMult(const ComplexVector &x, ComplexVector &y,
               const std::complex<double> a = 1.0) const override
  {
    ComplexVector t;
    Mult(x, t);
    y.Add(a, t);
  }
};

class IdentityEigenTestOperator : public ComplexOperator
{
public:
  IdentityEigenTestOperator(int n) : ComplexOperator(n) {}

  bool IsReal() const override { return true; }

  void Mult(const ComplexVector &x, ComplexVector &y) const override { y = x; }

  void AddMult(const ComplexVector &x, ComplexVector &y,
               const std::complex<double> a = 1.0) const override
  {
    y.Add(a, x);
  }
};

struct TridiagonalEigenTestResult
{
  std::vector<std::complex<double>> lambda;
  std::vector<ComplexVector> x;
  int its;  // Number of restarts
};

// Shift-and-invert Krylov-Schur solve for the nev eigenvalues of the tridiagonal operator
// closest to the shift s, with GMRES and Jacobi preconditioning for the linear solves, from
// the start vector v0 (a warm start if warm).
TridiagonalEigenTestResult SolveTridiagonalEigenTest(int n, double s, int nev, int ncv,
                                                     double tol, double lin_tol,
                                                     const ComplexVector &v0, bool warm)
{
  MPI_Comm comm = Mpi::World();
  TridiagonalEigenTestOperator K(n), A(n, s);
  IdentityEigenTestOperator M(n);
  auto gmres = std::make_unique<GmresSolver<ComplexOperator>>(comm, 0);
  gmres->SetRelTol(lin_tol);
  gmres->SetMaxIter(10 * n);
  gmres->SetRestartDim(n);
  ComplexKspSolver ksp(std::move(gmres),
                       std::make_unique<JacobiSmoother<ComplexOperator>>(comm));
  ksp.SetOperators(A, A);
  slepc::SlepcEPSSolver eps(comm, 0);
  eps.SetType(slepc::SlepcEigenvalueSolver::Type::KRYLOVSCHUR);
  eps.SetProblemType(slepc::SlepcEigenvalueSolver::ProblemType::GEN_NON_HERMITIAN);
  eps.SetOperators(K, M, EigenvalueSolver::ScaleType::NONE);
  eps.SetNumModes(nev, ncv);
  eps.SetTol(tol);
  eps.SetWhichEigenpairs(EigenvalueSolver::WhichType::TARGET_REAL);
  eps.SetLinearSolver(ksp);
  eps.SetShiftInvert(s);
  eps.SetInitialSpace(v0);
  eps.SetWarmStart(warm);
  const int num_conv = eps.Solve();
  REQUIRE(num_conv >= nev);
  TridiagonalEigenTestResult result;
  result.its = eps.GetNumIterations();
  for (int i = 0; i < nev; i++)
  {
    result.lambda.push_back(eps.GetEigenvalue(i));
    result.x.emplace_back(n);
    eps.GetEigenvector(i, result.x.back());
  }
  return result;
}

}  // namespace

TEST_CASE("Eigenmode stable Ritz value convergence", "[eigenmode][Serial]")
{
  // A warm-started shift-and-invert solve accepts an eigenpair only once its Ritz value is
  // the same eigenpair as at the previous restart. A start vector which already combines
  // the wanted eigenvectors then needs one more restart than its first Krylov cycle, which
  // converges on the error estimates alone without the warm start.
  constexpr int n = 120, nev = 2, ncv = 8;
  const double s = 1.5, tol = 1.0e-10, lin_tol = 1.0e-12;
  ComplexVector v0(n);
  linalg::SetRandom(Mpi::World(), v0);
  const auto cold = SolveTridiagonalEigenTest(n, s, nev, ncv, tol, lin_tol, v0, false);
  CHECK(cold.its > 1);

  ComplexVector v(n), r(n);
  linalg::SetRandom(Mpi::World(), r, 1);
  v = cold.x[0];
  v.Add(1.0, cold.x[1]);
  v.Add(1.0e-12 / std::sqrt(static_cast<double>(n)), r);
  const auto first = SolveTridiagonalEigenTest(n, s, nev, ncv, tol, lin_tol, v, false);
  const auto warm = SolveTridiagonalEigenTest(n, s, nev, ncv, tol, lin_tol, v, true);
  CHECK(first.its == 1);
  CHECK(warm.its == 2);
  for (int i = 0; i < nev; i++)
  {
    CHECK(std::abs(first.lambda[i] - cold.lambda[i]) <= 1.0e-8 * std::abs(cold.lambda[i]));
    CHECK(std::abs(warm.lambda[i] - cold.lambda[i]) <= 1.0e-8 * std::abs(cold.lambda[i]));
  }

  // Inexact linear solves, for a shift below the spectrum (GMRES needs the full dimension
  // for any tolerance at the interior shift): the Ritz values keep moving by about the
  // error of the linear solves, but the warm-started solve still stops, after at most a few
  // more restarts.
  const auto exact = SolveTridiagonalEigenTest(n, 0.0, nev, ncv, tol, lin_tol, v0, false);
  const auto inexact = SolveTridiagonalEigenTest(n, 0.0, nev, ncv, tol, 1.0e-6, v0, true);
  CHECK(inexact.its <= exact.its + 3);
  for (int i = 0; i < nev; i++)
  {
    CHECK(std::abs(inexact.lambda[i] - exact.lambda[i]) <=
          1.0e-4 * std::abs(exact.lambda[i]));
  }
}

#endif

}  // namespace palace
