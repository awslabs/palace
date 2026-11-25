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
#include "drivers/eigensolver.hpp"
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

TEST_CASE("Eigenmode warm-start gate", "[eigenmode][Serial]")
{
  // Records the solves of an adaptive run, cold (false) or warm-started (true) as the gate
  // decided, with their numbers of Krylov cycles, and returns the decision before each
  // solve and after the last.
  auto Run = [](bool q_region, const std::vector<std::pair<bool, int>> &solves)
  {
    WarmStartGate gate;
    std::vector<bool> next = {gate.Next(q_region)};
    for (const auto &[warm, its] : solves)
    {
      REQUIRE(warm == next.back());
      gate.Update(warm, its);
      next.push_back(gate.Next(q_region));
    }
    return next;
  };
  using V = std::vector<bool>;

  SECTION("Without the quality factor region")
  {
    // Only the cold solves of three cycles or more are followed by a warm start, whatever
    // the warm-started solves take.
    CHECK(Run(false, {{false, 2}, {false, 1}, {false, 3}, {true, 2}, {true, 1}}) ==
          V{false, false, false, true, true, true});
  }

  SECTION("With the quality factor region")
  {
    // Cold solves of two cycles are also followed by a warm start, as long as each
    // warm-started solve after one converges in one cycle.
    CHECK(Run(true, {{false, 1}, {false, 2}, {true, 1}, {true, 1}}) ==
          V{false, false, true, true, true});

    // After the first which does not, only cold solves of three cycles or more are.
    CHECK(
        Run(true, {{false, 2}, {true, 2}, {false, 2}, {false, 3}, {true, 2}, {true, 1}}) ==
        V{false, true, false, false, true, true, true});
  }
}

TEST_CASE("Eigenmode subspace growth", "[eigenmode][Serial]")
{
  // The configuration resolves an unset dimension to max(2 N, N + 15): 18 for 3 modes, and
  // 40 for 20.
  constexpr long long int size = 1000000;
  SubspaceGrowth growth;
  CHECK(growth.GetMaxSize(3, 18, size) == 18);

  // Without the quality factor region, or after solves of at most 10 cycles, the dimension
  // stays the configured one.
  growth.Update(false, 50);
  growth.Update(true, 10);
  CHECK(growth.GetMaxSize(3, 18, size) == 18);

  // With the region, after a solve of more than 10 cycles, the later solves use 6 more
  // vectors than the default, also after solves which converge fast.
  growth.Update(true, 11);
  CHECK(growth.GetMaxSize(3, 18, size) == 24);
  growth.Update(true, 1);
  CHECK(growth.GetMaxSize(3, 18, size) == 24);
  CHECK(growth.GetMaxSize(20, 40, size) == 46);

  // Another configured dimension is kept, and the dimension is at most the size of the
  // problem.
  CHECK(growth.GetMaxSize(3, 30, size) == 30);
  CHECK(growth.GetMaxSize(3, 18, 20) == 20);
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

// Diagonal (complex) operator on a single process.
class DiagonalEigenTestOperator : public ComplexOperator
{
private:
  const std::vector<std::complex<double>> d;

public:
  DiagonalEigenTestOperator(const std::vector<std::complex<double>> &d)
    : ComplexOperator(static_cast<int>(d.size())), d(d)
  {
  }

  bool IsReal() const override
  {
    return std::ranges::all_of(d, [](std::complex<double> x) { return x.imag() == 0.0; });
  }

  void AssembleDiagonal(ComplexVector &diag) const override
  {
    diag.SetSize(height);
    auto *DR = diag.Real().HostWrite();
    auto *DI = diag.Imag().HostWrite();
    for (int i = 0; i < height; i++)
    {
      DR[i] = d[i].real();
      DI[i] = d[i].imag();
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
      const std::complex<double> v = d[i] * std::complex<double>(XR[i], XI[i]);
      YR[i] = v.real();
      YI[i] = v.imag();
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

struct TridiagonalEigenTestResult
{
  std::vector<std::complex<double>> lambda;
  std::vector<ComplexVector> x;
  int its;      // Number of restarts
  int lin_its;  // Number of linear solver iterations
};

// Shift-and-invert Krylov-Schur solve for the nev eigenvalues of the tridiagonal operator
// closest to the shift s, with GMRES and Jacobi preconditioning for the linear solves, from
// the start vector v0, and with the linear solver tolerance first_tol for the first Krylov
// cycle if positive.
TridiagonalEigenTestResult SolveTridiagonalEigenTest(int n, double s, int nev, int ncv,
                                                     double tol, double lin_tol,
                                                     const ComplexVector &v0,
                                                     double first_tol = 0.0)
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
  eps.SetFirstCycleLinearSolverTol(first_tol);
  const int num_conv = eps.Solve();
  REQUIRE(num_conv >= nev);
  TridiagonalEigenTestResult result;
  result.its = eps.GetNumIterations();
  result.lin_its = ksp.NumTotalMultIterations();
  for (int i = 0; i < nev; i++)
  {
    result.lambda.push_back(eps.GetEigenvalue(i));
    result.x.emplace_back(n);
    eps.GetEigenvector(i, result.x.back());
  }
  return result;
}

// Quadratic problem K + λ C + λ² M with diagonal K and C and M = I, scaled, with
// eigenvalues λ = iω, Re{ω} = ω_r and Q = |ω| / (2 Im{ω}) for c = ω_r / (Q² - 1/4)^½ and
// k = Q² c², for the target λ = i. The mode closest to the target in Re{ω} has Q = 0.52,
// the next one Q = 100.
struct DampedEigenTestProblem
{
  static constexpr int n = 20;
  const std::complex<double> lambda_s = {0.0, 1.0};  // λ = iσ
  std::vector<double> omega_r = {1.05, 1.3}, q = {0.52, 100.0};
  std::vector<std::complex<double>> k, c, m, a;

  DampedEigenTestProblem() : k(n), c(n), m(n, 1.0), a(n)
  {
    for (int i = 2; i < n; i++)
    {
      omega_r.push_back(1.5 + 0.5 * i);
      q.push_back(50.0);
    }
    for (int i = 0; i < n; i++)
    {
      c[i] = omega_r[i] / std::sqrt(q[i] * q[i] - 0.25);
      k[i] = q[i] * q[i] * c[i] * c[i];
      a[i] = k[i] + lambda_s * c[i] + lambda_s * lambda_s * m[i];
    }
  }

  // Eigenvalue ω = λ / i of mode i.
  std::complex<double> Omega(int i) const { return {omega_r[i], 0.5 * c[i].real()}; }

  struct Result
  {
    std::complex<double> omega;  // First eigenvalue ω = λ / i
    int its;                     // Number of restarts
    double errest;               // SLEPc's relative error estimate of the first eigenpair
  };

  // Solve for one mode with SLEPc's Krylov-Schur solver for the linearization, for the
  // minimum quality factor q_min and tolerance tol, with ncv Krylov vectors, from a random
  // start vector or, if start >= 0, from the eigenvector of mode start plus noise times a
  // random vector (a warm start if warm).
  Result Solve(double q_min, int start = -1, bool warm = false, double tol = 1.0e-10,
               double noise = 1.0e-12, int ncv = 16) const
  {
    MPI_Comm comm = Mpi::World();
    DiagonalEigenTestOperator K(k), C(c), M(m), A(a);
    auto gmres = std::make_unique<GmresSolver<ComplexOperator>>(comm, 0);
    gmres->SetRelTol(1.0e-12);
    gmres->SetMaxIter(100);
    ComplexKspSolver ksp(std::move(gmres),
                         std::make_unique<JacobiSmoother<ComplexOperator>>(comm));
    ksp.SetOperators(A, A);
    slepc::SlepcPEPLinearSolver eps(comm, 0);
    eps.SetType(slepc::SlepcEigenvalueSolver::Type::KRYLOVSCHUR);
    eps.SetProblemType(slepc::SlepcEigenvalueSolver::ProblemType::GEN_NON_HERMITIAN);
    eps.SetOperators(K, C, M, EigenvalueSolver::ScaleType::NORM_2);
    eps.SetNumModes(1, ncv);
    eps.SetTol(tol);
    eps.SetLinearSolver(ksp);
    eps.SetShiftInvert(lambda_s);
    eps.SetWhichEigenpairs(EigenvalueSolver::WhichType::TARGET_IMAGINARY);
    eps.SetMinQualityFactor(q_min);
    if (start >= 0)
    {
      ComplexVector v(n), v_lambda(n), r(n);
      linalg::SetRandom(comm, r);
      v = 0.0;
      v.Real().HostReadWrite()[start] = 1.0;
      v.Add(noise / std::sqrt(static_cast<double>(n)), r);
      v_lambda = v;
      v_lambda *= std::complex<double>(-0.5 * c[start].real(), omega_r[start]);
      eps.SetInitialSpaceBlocks(v, v_lambda);
      eps.SetWarmStart(warm);
    }
    REQUIRE(eps.Solve() >= 1);
    PetscReal errest;
    PalacePetscCall(EPSGetErrorEstimate(eps, 0, &errest));
    return {eps.GetEigenvalue(0) / std::complex<double>(0.0, 1.0), eps.GetNumIterations(),
            errest};
  }
};

}  // namespace

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

  // Growing modes with the same bound on |Q|.
  CHECK(Inside(std::conj(Damped(1.0, 1.0e6))));
  CHECK(Inside(std::conj(Damped(1.0, 0.55))));
  CHECK(!Inside(std::conj(Damped(1.0, 0.53))));
  CHECK(!Inside(std::conj(Damped(100.0, 0.53))));
  PalacePetscCall(RGDestroy(&rg));
}

TEST_CASE("Eigenmode quality factor region in a damped solve", "[eigenmode][Serial]")
{
  // The solve for one mode returns the second mode of the damped test problem, with
  // Q = 100, with MinQ = 0.54, and the first one, with Q = 0.52, without.
  const DampedEigenTestProblem p;
  for (double q_min : {0.0, 0.54})
  {
    const int i = (q_min > 0.5) ? 1 : 0;
    const auto result = p.Solve(q_min);
    CHECK(std::abs(result.omega - p.Omega(i)) <= 1.0e-8 * p.omega_r[i]);
  }
}

TEST_CASE("Eigenmode first restart of a warm start in the quality factor region",
          "[eigenmode][Serial]")
{
  // With the quality factor region, a warm-started solve accepts an eigenpair at its first
  // restart only if its error estimate is below a tenth of the tolerance. The start vector
  // is the wanted eigenvector plus noise, which sets the error estimate at the first
  // restart, measured by a cold solve which converges there. The warm-started solve takes
  // one more restart when this estimate is 0.3 times the tolerance, but not when it is 0.05
  // times the tolerance, nor without the region.
  const DampedEigenTestProblem p;
  constexpr double noise = 1.0e-4;
  constexpr int ncv = 6;
  for (double q_min : {0.54, 0.0})
  {
    const int i = (q_min > 0.5) ? 1 : 0;
    const auto first = p.Solve(q_min, i, false, 1.0e-3, noise, ncv);
    REQUIRE(first.its == 1);
    REQUIRE(first.errest > 1.0e-12);
    for (double ratio : {0.3, 0.05})
    {
      const double tol = first.errest / ratio;
      const auto cold = p.Solve(q_min, i, false, tol, noise, ncv);
      const auto warm = p.Solve(q_min, i, true, tol, noise, ncv);
      CHECK(cold.its == 1);
      CHECK(warm.its == ((q_min > 0.5 && ratio > 0.1) ? 2 : 1));
      CHECK(std::abs(warm.omega - p.Omega(i)) <= 1.0e-6 * p.omega_r[i]);
    }
  }

  // The margin applies to the first restart only: a solve which accepts the eigenpair at
  // its second restart with an estimate of 0.3 times the tolerance takes two restarts
  // whether warm-started or not.
  const auto first = p.Solve(0.54, 1, false, 1.0e-3, 100.0 * noise, ncv / 2);
  REQUIRE(first.its == 1);
  const auto second = p.Solve(0.54, 1, false, 0.5 * first.errest, 100.0 * noise, ncv / 2);
  REQUIRE(second.its == 2);
  const double tol = second.errest / 0.3;
  REQUIRE(tol < first.errest);
  CHECK(p.Solve(0.54, 1, true, tol, 100.0 * noise, ncv / 2).its == 2);
}

#endif

}  // namespace palace
