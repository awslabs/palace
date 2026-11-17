// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// [experiment] Review of PR 1203: accuracy and cost of warm-started shift-and-invert
// Krylov-Schur solves with inexact linear solves, on a model problem.

#include <cmath>
#include <complex>
#include <cstdio>
#include <numbers>
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

#if defined(PALACE_WITH_SLEPC)

namespace palace
{

namespace
{

class TridiagOp : public ComplexOperator
{
private:
  const double s;

public:
  TridiagOp(int n, double s = 0.0) : ComplexOperator(n), s(s) {}
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
      YR[i] = d * XR[i] - ((i > 0) ? XR[i - 1] : 0.0) - ((i < height - 1) ? XR[i + 1] : 0.0);
      YI[i] = d * XI[i] - ((i > 0) ? XI[i - 1] : 0.0) - ((i < height - 1) ? XI[i + 1] : 0.0);
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

class IdOp : public ComplexOperator
{
public:
  IdOp(int n) : ComplexOperator(n) {}
  bool IsReal() const override { return true; }
  void Mult(const ComplexVector &x, ComplexVector &y) const override { y = x; }
  void AddMult(const ComplexVector &x, ComplexVector &y,
               const std::complex<double> a = 1.0) const override
  {
    y.Add(a, x);
  }
};

struct Result
{
  std::vector<std::complex<double>> lambda;
  std::vector<ComplexVector> x;
  int its, nlin, nlinits;
};

Result Solve(int n, double s, int nev, int ncv, double tol, double lin_tol,
             const ComplexVector &v0, bool settled, double first_tol)
{
  MPI_Comm comm = Mpi::World();
  TridiagOp K(n), A(n, s);
  IdOp M(n);
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
  eps.SetWarmStart(settled);
  if (first_tol > 0.0)
  {
    eps.SetFirstCycleLinearSolverTol(first_tol);
  }
  const int num_conv = eps.Solve();
  REQUIRE(num_conv >= nev);
  Result r;
  r.its = eps.GetNumIterations();
  r.nlin = ksp.NumTotalMult();
  r.nlinits = ksp.NumTotalMultIterations();
  for (int i = 0; i < nev; i++)
  {
    r.lambda.push_back(eps.GetEigenvalue(i));
    r.x.emplace_back(n);
    eps.GetEigenvector(i, r.x.back());
  }
  return r;
}

}  // namespace

TEST_CASE("Warm start experiment", "[exp][Serial]")
{
  const int n = 400;
  const double s = 0.0;
  for (int nev : {2, 6})
  {
    const int ncv = std::max(2 * nev, nev + 15);
    ComplexVector r0(n), noise(n);
    linalg::SetRandom(Mpi::World(), r0);
    const auto ref = Solve(n, s, nev, ncv, 1.0e-13, 1.0e-14, r0, false, 0.0);
    for (double delta : {1.0e-2, 1.0e-4})
    {
      // Warm start: exact eigenvectors perturbed by delta (relative), unit weights with
      // golden-ratio phases, plus 1e-3 of a random vector.
      ComplexVector v(n);
      v = 0.0;
      for (int i = 0; i < nev; i++)
      {
        ComplexVector xi(ref.x[i]);
        linalg::SetRandom(Mpi::World(), noise, 100 + i);
        const double nx = linalg::Norml2(Mpi::World(), xi);
        const double nn = linalg::Norml2(Mpi::World(), noise);
        xi.Add(delta * nx / nn, noise);
        const double nxi = linalg::Norml2(Mpi::World(), xi);
        const double phase =
            2.0 * std::numbers::pi * std::fmod(0.5 * (std::sqrt(5.0) - 1.0) * i, 1.0);
        v.AXPY(std::polar(1.0 / nxi, phase), xi);
      }
      const double nv = linalg::Norml2(Mpi::World(), v);
      const double nr = linalg::Norml2(Mpi::World(), r0);
      v.Add(1.0e-3 * nv / nr, r0);
      for (double tol : {1.0e-6, 1.0e-8, 1.0e-10})
      {
        for (double lin_tol : {1.0e-6, 1.0e-8, 1.0e-10})
        {
          if (lin_tol > 100.0 * tol && tol < 1.0e-8)
          {
            continue;
          }
          struct Var
          {
            const char *name;
            bool warm, settled;
            double first_tol;
          };
          const double ft = std::max(1.0e-2 * lin_tol, std::min(lin_tol, 1.0e-12));
          const Var vars[] = {{"cold", false, false, 0.0},
                              {"warm", true, false, 0.0},
                              {"warm+settled", true, true, 0.0},
                              {"warm+tight", true, false, ft},
                              {"warm+settled+tight(PR)", true, true, ft}};
          for (const auto &var : vars)
          {
            const auto res =
                Solve(n, s, nev, ncv, tol, lin_tol, var.warm ? v : r0, var.settled,
                      var.first_tol);
            double err = 0.0;
            for (int i = 0; i < nev; i++)
            {
              // Match to the nearest reference eigenvalue.
              double e = mfem::infinity();
              for (int j = 0; j < nev; j++)
              {
                e = std::min(e, std::abs(res.lambda[i] - ref.lambda[j]) /
                                    std::abs(ref.lambda[j]));
              }
              err = std::max(err, e);
            }
            std::printf("nev=%d delta=%.0e tol=%.0e lin_tol=%.0e %-24s its=%2d "
                        "solves=%3d gmres=%5d max|dl|/|l|=%.2e (%.2f tol)\n",
                        nev, delta, tol, lin_tol, var.name, res.its, res.nlin,
                        res.nlinits, err, err / tol);
          }
        }
      }
    }
  }
}

}  // namespace palace

#endif
