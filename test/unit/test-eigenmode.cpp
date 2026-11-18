// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <complex>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include "fem/fieldtransfer.hpp"
#include "linalg/eps.hpp"
#include "linalg/iterative.hpp"
#include "linalg/jacobi.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/slepc.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

#if defined(PALACE_WITH_SLEPC)
#include <slepc.h>
#endif

namespace palace
{

TEST_CASE("Vector field transfer through refinement and rebalancing",
          "[eigenmode][Serial][Parallel][GPU]")
{
  // H(curl) fields of a space of order p, linear in the coordinates, are kept exactly
  // through nonconforming refinement and rebalancing of their mesh.
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
    tv.UseDevice(true);
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
    tv.UseDevice(true);
    transfer.GetField(k, nd_fespace, tv);
    REQUIRE(tv.Size() == ref.Size());
    tv -= ref;
    const double err = linalg::Norml2(comm, tv), norm = linalg::Norml2(comm, ref);
    CHECK(err <= 1.0e-12 * norm);
  }
}

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

// GMRES which records the relative tolerance of each solve.
class TolRecordingGmresSolver : public GmresSolver<ComplexOperator>
{
public:
  mutable std::vector<double> tols;

  using GmresSolver<ComplexOperator>::GmresSolver;

  void Mult(const ComplexVector &b, ComplexVector &x) const override
  {
    tols.push_back(GetRelTol());
    GmresSolver<ComplexOperator>::Mult(b, x);
  }
};

struct TridiagonalEigenTestResult
{
  std::vector<std::complex<double>> lambda;
  std::vector<ComplexVector> x;
  int its;                       // Number of restarts
  int lin_its;                   // Number of linear solver iterations
  std::vector<double> lin_tols;  // Tolerance of each linear solve
};

// Shift-and-invert Krylov-Schur solve for the nev eigenvalues of the tridiagonal operator
// closest to the shift s, with GMRES and Jacobi preconditioning for the linear solves, from
// the start vector v0, and with the linear solver tolerance first_tol for the first
// num_first linear solves (the first Krylov cycle if negative) if positive.
TridiagonalEigenTestResult SolveTridiagonalEigenTest(int n, double s, int nev, int ncv,
                                                     double tol, double lin_tol,
                                                     const ComplexVector &v0,
                                                     double first_tol = 0.0,
                                                     int num_first = -1)
{
  MPI_Comm comm = Mpi::World();
  TridiagonalEigenTestOperator K(n), A(n, s);
  IdentityEigenTestOperator M(n);
  auto gmres = std::make_unique<TolRecordingGmresSolver>(comm, 0);
  const auto *recorder = gmres.get();
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
  eps.SetFirstCycleLinearSolverTol(first_tol, num_first);
  const int num_conv = eps.Solve();
  REQUIRE(num_conv >= nev);
  TridiagonalEigenTestResult result;
  result.its = eps.GetNumIterations();
  result.lin_its = ksp.NumTotalMultIterations();
  result.lin_tols = recorder->tols;
  for (int i = 0; i < nev; i++)
  {
    result.lambda.push_back(eps.GetEigenvalue(i));
    result.x.emplace_back(n);
    eps.GetEigenvector(i, result.x.back());
  }
  return result;
}

}  // namespace

TEST_CASE("Eigenmode tight first linear solves", "[eigenmode][Serial]")
{
  // The tighter tolerance applies to the given number of first linear solves, or to the
  // first Krylov cycle (ncv solves) if negative.
  constexpr int n = 120, nev = 2, ncv = 8, num_first = 3;
  const double tol = 1.0e-8, lin_tol = 1.0e-10, first_tol = 1.0e-12;
  ComplexVector v0(n);
  linalg::SetRandom(Mpi::World(), v0);
  for (int num : {num_first, -1})
  {
    const auto result =
        SolveTridiagonalEigenTest(n, 0.0, nev, ncv, tol, lin_tol, v0, first_tol, num);
    const std::size_t num_tight = (num > 0) ? num : ncv;
    REQUIRE(result.its > 1);
    REQUIRE(result.lin_tols.size() > num_tight);
    for (std::size_t i = 0; i < result.lin_tols.size(); i++)
    {
      CHECK(result.lin_tols[i] == ((i < num_tight) ? first_tol : lin_tol));
    }
  }
}

TEST_CASE("Eigenmode warm start", "[eigenmode][Serial]")
{
  // From a start vector combining the wanted eigenvectors, the eigenvalues are only as
  // accurate as the linear solves of the first Krylov cycle, whose error the error estimate
  // of the eigenvalue solver does not see: a smaller tolerance for that cycle gives the
  // accuracy of exact solves. The shift is below the spectrum, so that GMRES converges for
  // any tolerance.
  constexpr int n = 120, nev = 2, ncv = 8;
  const double s = 0.0, tol = 1.0e-10, exact_tol = 1.0e-12, lin_tol = 1.0e-4;
  ComplexVector v0(n);
  linalg::SetRandom(Mpi::World(), v0);
  const auto exact = SolveTridiagonalEigenTest(n, s, nev, ncv, tol, exact_tol, v0);
  CHECK(exact.its > 1);

  ComplexVector v(n), r(n);
  linalg::SetRandom(Mpi::World(), r, 1);
  v = exact.x[0];
  v.Add(1.0, exact.x[1]);
  v.Add(1.0e-12 / std::sqrt(static_cast<double>(n)), r);
  const auto warm = SolveTridiagonalEigenTest(n, s, nev, ncv, tol, lin_tol, v);
  const auto warm_tight =
      SolveTridiagonalEigenTest(n, s, nev, ncv, tol, lin_tol, v, exact_tol);
  CHECK(warm_tight.its == 1);
  double err = 0.0, err_tight = 0.0;
  for (int i = 0; i < nev; i++)
  {
    const double a = std::abs(exact.lambda[i]);
    err = std::max(err, std::abs(warm.lambda[i] - exact.lambda[i]) / a);
    err_tight = std::max(err_tight, std::abs(warm_tight.lambda[i] - exact.lambda[i]) / a);
  }
  CHECK(err > 100.0 * tol);
  CHECK(err_tight <= tol);

  // The tolerance is restored after the first cycle: from a random vector, the solve takes
  // more linear solver iterations than without the smaller tolerance, and fewer than with
  // it for all cycles.
  const auto cold = SolveTridiagonalEigenTest(n, s, nev, ncv, tol, lin_tol, v0);
  const auto cold_tight =
      SolveTridiagonalEigenTest(n, s, nev, ncv, tol, lin_tol, v0, exact_tol);
  CHECK(cold_tight.lin_its > cold.lin_its);
  CHECK(cold_tight.lin_its < exact.lin_its);
}

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

}  // namespace palace
