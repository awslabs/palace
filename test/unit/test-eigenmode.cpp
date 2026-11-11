// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <complex>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include "fem/fieldtransfer.hpp"
#include "linalg/eps.hpp"
#include "linalg/iterative.hpp"
#include "linalg/jacobi.hpp"
#include "linalg/ksp.hpp"
#include "linalg/slepc.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

#if defined(PALACE_WITH_SLEPC)
#include <slepc.h>
#endif

using namespace palace;

TEST_CASE("Eigenmode linear solver relaxation", "[eigenmode][Serial]")
{
  // The tolerance grows as the error estimates of the wanted eigenpairs decrease, between
  // the linear solver tolerance and 1e-4.
  const double eig_tol = 1.0e-11, lin_tol = 1.0e-12, c = 1.0e-3;
  CHECK(GetRelaxedLinearSolverTol(c, eig_tol, lin_tol, 1.0) == lin_tol);
  CHECK(GetRelaxedLinearSolverTol(c, eig_tol, lin_tol, 1.0e-3) == Catch::Approx(1.0e-11));
  CHECK(GetRelaxedLinearSolverTol(c, eig_tol, lin_tol, 1.0e-8) == Catch::Approx(1.0e-6));
  CHECK(GetRelaxedLinearSolverTol(c, eig_tol, lin_tol, 1.0e-12) == 1.0e-4);

  // No relaxation without wanted unconverged eigenpairs or with c = 0, and never below the
  // linear solver tolerance.
  CHECK(GetRelaxedLinearSolverTol(c, eig_tol, lin_tol, 0.0) == lin_tol);
  CHECK(GetRelaxedLinearSolverTol(0.0, eig_tol, lin_tol, 1.0e-8) == lin_tol);
  CHECK(GetRelaxedLinearSolverTol(c, eig_tol, 1.0e-3, 1.0e-12) == 1.0e-3);

  // Less relaxation for wanted eigenvalues far in the interior of the transformed spectrum,
  // whose error estimates amplify the error of the linear solves by the ratio, and for
  // eigenvalues close to others (relative separation gap): c <= 0.01 gap / ratio.
  CHECK(GetRelaxedLinearSolverTol(c, eig_tol, lin_tol, 1.0e-8, 10.0) ==
        Catch::Approx(1.0e-6));
  CHECK(GetRelaxedLinearSolverTol(c, eig_tol, lin_tol, 1.0e-8, 100.0) ==
        Catch::Approx(1.0e-7));
  CHECK(GetRelaxedLinearSolverTol(c, eig_tol, lin_tol, 1.0e-8, 1.0, 1.0e-2) ==
        Catch::Approx(1.0e-7));
  CHECK(GetRelaxedLinearSolverTol(c, eig_tol, lin_tol, 1.0e-8, 20.0, 1.0e-3) ==
        Catch::Approx(5.0e-10));
  CHECK(GetRelaxedLinearSolverTol(c, eig_tol, lin_tol, 1.0e-8, 1.0e4, 1.0e-4) == lin_tol);
}

#if defined(PALACE_WITH_SLEPC)

namespace
{

// The symmetric tridiagonal stencil [-1, 2 + i / n, -1] minus s times the identity (real,
// as a complex operator) on a single process: distinct eigenvalues in (0, 5).
class TridiagonalEigenTestOperator : public ComplexOperator
{
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
  int its;         // Number of restarts
  int num_solves;  // Number of linear solves
  int num_its;     // Total number of linear solver iterations
};

// Shift-and-invert Krylov-Schur solve for the nev eigenvalues of the tridiagonal operator
// closest to the target s (or t > 0, with the shift s), with GMRES and Jacobi
// preconditioning for the linear solves, from the start vector v0 (a warm start if warm),
// with the linear solver relaxation constant.
TridiagonalEigenTestResult SolveTridiagonalEigenTest(int n, double s, int nev, int ncv,
                                                     double tol, double lin_tol,
                                                     const ComplexVector &v0, bool warm,
                                                     double relaxation = 0.0,
                                                     double t = 0.0)
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
  if (t > 0.0 && t != s)
  {
    eps.SetTarget(t);
  }
  eps.SetInitialSpace(v0);
  eps.SetWarmStart(warm);
  eps.SetLinearSolverRelaxation(relaxation);
  const int num_conv = eps.Solve();
  REQUIRE(num_conv >= nev);
  TridiagonalEigenTestResult result;
  result.its = eps.GetNumIterations();
  result.num_solves = ksp.NumTotalMult();
  result.num_its = ksp.NumTotalMultIterations();
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
  // The shift-and-invert solve accepts an eigenpair only once its Ritz value is the same
  // eigenpair as at the previous restart. A warm start which already combines the wanted
  // eigenvectors then needs one more restart than its first Krylov cycle, while the first
  // cycle of a cold solve converges on the error estimates alone.
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
  // for any tolerance at the interior shift): the solve still stops, after at most a few
  // more restarts.
  const auto exact = SolveTridiagonalEigenTest(n, 0.0, nev, ncv, tol, lin_tol, v0, false);
  const auto inexact = SolveTridiagonalEigenTest(n, 0.0, nev, ncv, tol, 1.0e-6, v0, false);
  CHECK(inexact.its <= exact.its + 3);
  for (int i = 0; i < nev; i++)
  {
    CHECK(std::abs(inexact.lambda[i] - exact.lambda[i]) <=
          1.0e-4 * std::abs(exact.lambda[i]));
  }
}

TEST_CASE("Eigenmode relaxed linear solves", "[eigenmode][Serial]")
{
  // With the relaxation of the linear solver tolerance, the eigenvalues agree with those of
  // strict linear solves to the eigenvalue solver tolerance, for fewer linear solver
  // iterations per linear solve (shift below the spectrum, see above).
  constexpr int n = 200, nev = 3, ncv = 10;
  const double s = 0.0, tol = 1.0e-10, lin_tol = 1.0e-12;
  ComplexVector v0(n);
  linalg::SetRandom(Mpi::World(), v0);
  const auto strict = SolveTridiagonalEigenTest(n, s, nev, ncv, tol, lin_tol, v0, false);
  const auto relaxed =
      SolveTridiagonalEigenTest(n, s, nev, ncv, tol, lin_tol, v0, false, 1.0e-3);
  CHECK(strict.its > 1);
  CHECK(relaxed.num_its * strict.num_solves < strict.num_its * relaxed.num_solves);
  for (int i = 0; i < nev; i++)
  {
    CHECK(std::abs(relaxed.lambda[i] - strict.lambda[i]) <=
          10.0 * tol * std::abs(strict.lambda[i]));
  }
}

TEST_CASE("Eigenmode shift apart from the target", "[eigenmode][Serial]")
{
  // The nev eigenvalues above the target t nearest to it, with the shift of the spectral
  // transformation at the target or in the middle of the target and the farthest of them.
  constexpr int n = 120, nev = 4, ncv = 12;
  const double t = 1.5, tol = 1.0e-10, lin_tol = 1.0e-12;
  ComplexVector v0(n);
  linalg::SetRandom(Mpi::World(), v0);
  const auto ref = SolveTridiagonalEigenTest(n, t, nev, ncv, tol, lin_tol, v0, false);
  const double s = 0.5 * (t + ref.lambda[nev - 1].real());
  const auto shifted =
      SolveTridiagonalEigenTest(n, s, nev, ncv, tol, lin_tol, v0, false, 0.0, t);
  for (int i = 0; i < nev; i++)
  {
    CHECK(ref.lambda[i].real() >= t);
    CHECK(std::abs(shifted.lambda[i] - ref.lambda[i]) <= 1.0e-8 * std::abs(ref.lambda[i]));
  }
}

#endif

TEST_CASE("Eigenmode shift from eigenvalue estimates", "[eigenmode][Serial]")
{
  // The shift is τ + f (ωₙ - τ), f <= 1/2, unless it is within a quarter of the interval
  // between consecutive estimates (or τ) around it from one of them.
  auto Mode = [](double w, double q)
  { return std::complex<double>(-0.5 * w / q, w); };  // λ = iω, Q ≈ ω / (-2 Re{λ})

  // Wanted modes ω = 4.15, 4.35, ..., 6.15 above the target 4 (Q = 1e5), then 6.35.
  const double t = 4.0;
  constexpr int n = 11;
  std::vector<std::complex<double>> lambda;
  for (int k = 0; k <= n; k++)
  {
    lambda.push_back(Mode(4.15 + 0.2 * k, 1.0e5));
  }
  CHECK(GetEigenvalueShift(t, lambda, n) == Catch::Approx(5.075));
  CHECK(GetEigenvalueShift(t, lambda, n, 1.0) == Catch::Approx(5.075));
  CHECK(GetEigenvalueShift(t, lambda, n, 0.3) == Catch::Approx(4.645));
  CHECK(GetEigenvalueShift(t, lambda, n, 0.36) == Catch::Approx(4.65));  // Near 4.75
  CHECK(GetEigenvalueShift(t, lambda, n, 0.03) == Catch::Approx(4.0645));
  CHECK(GetEigenvalueShift(t, lambda, n, 0.01) == t);  // Near the target
  CHECK(GetEigenvalueShift(t, lambda, n, 0.0) == t);

  // A strongly damped estimate does not move the shift.
  auto damped = lambda;
  damped.push_back(Mode(5.08, 1.0));
  CHECK(GetEigenvalueShift(t, damped, n) == Catch::Approx(5.075));

  // Sparse spectra: the shift stays clear of the eigenvalues; a strongly damped wanted
  // mode counts for ωₙ.
  const std::vector<std::complex<double>> sparse = {Mode(5.2, 1.0e6), Mode(15.5, 1.0e6),
                                                    Mode(25.7, 1.0e6), Mode(31.8, 0.7),
                                                    Mode(35.9, 1.0e6)};
  CHECK(GetEigenvalueShift(2.0, sparse, 3) == Catch::Approx(10.35));  // Near 15.5 (above)
  CHECK(GetEigenvalueShift(2.0, sparse, 4) == Catch::Approx(10.35));  // Near 15.5 (below)
  CHECK(GetEigenvalueShift(3.0, {Mode(5.36, 1.0e6), Mode(15.96, 1.0e6)}, 2) ==
        Catch::Approx(9.48));
  const std::vector<std::complex<double>> degenerate = {Mode(4.5, 1.0e6), Mode(4.5, 1.0e6),
                                                        Mode(5.0, 1.0e6), Mode(6.0, 1.0e6)};
  CHECK(GetEigenvalueShift(t, degenerate, 4, 0.3) == Catch::Approx(4.25));
  const std::vector<std::complex<double>> cluster = {Mode(3.0, 1.0e6), Mode(3.05, 1.0e6),
                                                     Mode(4.0, 1.0e6), Mode(5.0, 1.0e6)};
  CHECK(GetEigenvalueShift(2.0, cluster, 4, 0.34) ==
        Catch::Approx(2.5));  // Not in the pair
  CHECK(GetEigenvalueShift(3.0, {Mode(5.36, 1.0e6), Mode(15.96, 1.0e6)}, 2, 0.3) == 3.0);

  // Without estimates of every wanted mode, the shift is the target.
  CHECK(GetEigenvalueShift(
            t, std::vector<std::complex<double>>(lambda.begin(), lambda.begin() + n - 1),
            n) == t);

  // Check of a shifted solve: n converged and the shift at most in the middle of the target
  // and the farthest requested eigenvalue.
  CHECK(CheckEigenvalueShift(t, 5.075, lambda, n));
  CHECK(!CheckEigenvalueShift(t, 5.08, lambda, n));
  CHECK(!CheckEigenvalueShift(t, 5.0, lambda, n + 2));
  auto entered = lambda;
  entered.insert(entered.begin() + 3, Mode(4.7, 1.0e5));  // ωₙ = 5.95
  CHECK(CheckEigenvalueShift(t, 4.975, entered, n));
  CHECK(!CheckEigenvalueShift(t, 5.0, entered, n));
}

#if defined(PALACE_WITH_SLEPC)

TEST_CASE("Eigenmode quality factor region", "[eigenmode][Serial]")
{
  // Region for eigenvalues λ = iω above the target λ = i s, with Q = |ω| / (2 Im{ω}).
  const double s = 0.5, q_min = 0.54;
  RG rg;
  PalacePetscCall(RGCreate(PETSC_COMM_SELF, &rg));
  slepc::ConfigureQualityFactorRegion(rg, s, q_min);
  auto Inside = [&rg](std::complex<double> omega)
  {
    PetscScalar lambda = std::complex<double>(0.0, 1.0) * omega;
    PetscInt inside;
    PalacePetscCall(RGCheckInside(rg, 1, &lambda, nullptr, &inside));
    return inside > 0;
  };
  auto Damped = [](double re, double q)
  { return std::complex<double>(re, re / std::sqrt(4.0 * q * q - 1.0)); };

  // Damped modes above the target.
  CHECK(Inside({1.0, 0.0}));
  CHECK(Inside(Damped(1.0, 1.0e6)));
  CHECK(Inside(Damped(1.0, 0.6)));
  CHECK(Inside(Damped(1.0, 0.55)));
  CHECK(Inside(Damped(100.0, 0.55)));
  CHECK(!Inside(Damped(1.0, 0.53)));
  CHECK(!Inside(Damped(1.0, 0.501)));
  CHECK(!Inside(Damped(100.0, 0.53)));

  // Modes below the target.
  CHECK(!Inside({0.4, 0.0}));
  CHECK(!Inside(Damped(0.4, 10.0)));

  // Growing modes are kept for |Q| >= 10 only.
  CHECK(Inside(std::conj(Damped(1.0, 1.0e6))));
  CHECK(Inside(std::conj(Damped(1.0, 10.5))));
  CHECK(!Inside(std::conj(Damped(1.0, 9.5))));
  CHECK(!Inside(std::conj(Damped(1.0, 1.0))));
  PalacePetscCall(RGDestroy(&rg));
}

#endif

TEST_CASE("Vector field transfer through refinement and rebalancing",
          "[eigenmode][Serial][Parallel]")
{
  // H(curl) fields of a space of order p, linear in the coordinates, are kept exactly
  // through nonconforming refinement and rebalancing of their mesh. MFEM's transfer of
  // H(curl) grid functions through ParMesh::Rebalance does not support this.
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  const int order = GENERATE(1, 2);
  MPI_Comm comm = Mpi::World();
  auto smesh = mfem::Mesh::MakeCartesian3D(3, 3, 3, type, 1.0, 1.0, 1.0);
  smesh.EnsureNCMesh(true);
  mfem::ParMesh pmesh(comm, smesh);
  mfem::ND_FECollection nd_fec(order, pmesh.Dimension());

  // Two fields: a constant one, exact for every order, and a linear one for order 2.
  auto Field = [order](const mfem::Vector &x, int k, mfem::Vector &E)
  {
    E.SetSize(3);
    E(0) = 1.0 + ((order > 1) ? 2.0 * x(1) - x(2) : 0.0);
    E(1) = -0.5 + ((order > 1) ? 3.0 * x(0) + x(2) : 0.0);
    E(2) = 2.0 + ((order > 1) ? x(0) - 4.0 * x(1) : 0.0);
    if (k == 1)
    {
      E *= -0.25;
      E(0) += 1.0;
    }
  };
  auto TrueDofs = [&](mfem::ParFiniteElementSpace &fespace, int k)
  {
    mfem::VectorFunctionCoefficient coeff(
        3, [&Field, k](const mfem::Vector &x, mfem::Vector &E) { Field(x, k, E); });
    mfem::ParGridFunction gf(&fespace);
    gf.ProjectCoefficient(coeff);
    Vector tv(fespace.GetTrueVSize());
    gf.GetTrueDofs(tv);
    return tv;
  };

  VectorFieldTransfer transfer(pmesh, order);
  {
    mfem::ParFiniteElementSpace nd_fespace(&pmesh, &nd_fec);
    for (int k = 0; k < 2; k++)
    {
      transfer.AddField(nd_fespace, TrueDofs(nd_fespace, k));
    }
  }
  CHECK(transfer.Size() == 2u);
  CHECK(transfer.IsCurrent(pmesh));

  // Refine the elements in a corner, twice, then rebalance.
  for (int it = 0; it < 2; it++)
  {
    mfem::Array<int> marked;
    for (int e = 0; e < pmesh.GetNE(); e++)
    {
      mfem::Vector c;
      pmesh.GetElementCenter(e, c);
      if (c(0) < 0.5 && c(1) < 0.5)
      {
        marked.Append(e);
      }
    }
    pmesh.GeneralRefinement(marked, 1, 0);
    CHECK(!transfer.IsCurrent(pmesh));
    transfer.Update();
    CHECK(transfer.IsCurrent(pmesh));
  }
  pmesh.Rebalance();
  transfer.Update();
  CHECK(transfer.IsCurrent(pmesh));

  mfem::ParFiniteElementSpace nd_fespace(&pmesh, &nd_fec);
  for (int k = 0; k < 2; k++)
  {
    const auto ref = TrueDofs(nd_fespace, k);
    Vector tv;
    transfer.GetField(k, nd_fespace, tv);
    REQUIRE(tv.Size() == ref.Size());
    tv -= ref;
    const double err = linalg::Norml2(comm, tv), norm = linalg::Norml2(comm, ref);
    CHECK(err <= 1.0e-12 * norm);
  }
}
