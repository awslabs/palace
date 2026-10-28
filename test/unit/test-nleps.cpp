// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <complex>
#include <memory>
#include <numbers>
#include <utility>
#include <vector>
#include <Eigen/Dense>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
#include "linalg/eps.hpp"
#include "linalg/iterative.hpp"
#include "linalg/ksp.hpp"
#include "linalg/nleps.hpp"
#include "linalg/rap.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

using namespace palace;

// The Newton step on the pole-cleared d(λ) T(λ) uses L = d′/d, and the line search compares
// residuals scaled by |d|.
TEST_CASE("Pole-cleared Newton terms", "[nleps][Serial]")
{
  const std::complex<double> p(-0.1, 1.5);
  const nleps::PoleTerms poles{{p, std::conj(p)}};
  auto d = [&](std::complex<double> l) { return (l - p) * (l - std::conj(p)); };

  const std::complex<double> l(-0.02, 1.7), s(0.0, 2.0);
  const double h = 1.0e-6;
  const std::complex<double> dd = (d(l + h) - d(l - h)) / (2.0 * h);
  CHECK(std::abs(poles.LogDerivative(l) - dd / d(l)) <= 1.0e-8 * std::abs(dd / d(l)));
  CHECK(std::abs(poles.RelativeScale(l, s) - std::abs(d(l) / d(s))) <= 1.0e-14);

  CHECK(poles.IsAtPole(p));
  CHECK_FALSE(poles.IsAtPole(l));

  // A sample node on a pole is moved off it: here the middle node of five.
  {
    const nleps::PoleTerms undamped{{std::complex<double>(0.0, 1.0)}};
    const std::complex<double> c(0.0, 1.0), r(0.0, 0.5);
    for (int j = 0; j < 5; j++)
    {
      CHECK_FALSE(undamped.IsAtPole(c + r * undamped.SampleNode(j, 5, c, r), 1.0e-6));
    }
    CHECK(undamped.IsAtPole(c + r * std::cos(std::numbers::pi * 5 / 10.0), 1.0e-6));
  }

  // Pole proximity is relative, so it does not depend on the frequency scale.
  for (double scale : {1.0e-7, 1.0, 1.0e7})
  {
    CAPTURE(scale);
    const nleps::PoleTerms scaled{{scale * p}};
    CHECK(scaled.IsAtPole(scale * p * (1.0 + 1.0e-9)));
    CHECK_FALSE(scaled.IsAtPole(scale * p * 0.95));
  }
  CHECK(poles.Separates(std::complex<double>(0.0, 1.0), s));
  CHECK_FALSE(poles.Separates(l, s));
}

namespace
{

// The nonlinear eigenvalue problem T(λ) x = (K + λ² (1 + χ(λ)) I) x = 0 with a diagonal
// K = diag(k_j), one degree of freedom per element of an L2 space of order 0 (so the
// operators are the actual parallel operators of the solver), and the damped Lorentz
// susceptibility χ(λ) = Δε ω₀² / (λ² + γ λ + ω₀²) of every dof, with poles at the roots
// of λ² + γ λ + ω₀². Each k_j has the exact roots of the quartic
//   λ⁴ + γ λ³ + (ω₀² (1 + Δε) + k_j) λ² + k_j γ λ + k_j ω₀² = 0,
// one on each side of the pole for each mode, as for a mode shape in a Lorentz medium.
constexpr double w0 = 1.0, delta_eps = 0.5, gamma_l = 0.05;

std::complex<double> Chi(std::complex<double> l)
{
  return delta_eps * w0 * w0 / (l * l + gamma_l * l + w0 * w0);
}

class DiagonalOperator : public Operator
{
  Vector d;

public:
  explicit DiagonalOperator(Vector &&d) : Operator(d.Size()), d(std::move(d)) {}
  void Mult(const Vector &x, Vector &y) const override
  {
    const auto *hd = d.HostRead();
    const auto *hx = x.HostRead();
    auto *hy = y.HostWrite();
    for (int i = 0; i < d.Size(); i++)
    {
      hy[i] = hd[i] * hx[i];
    }
  }
  void MultTranspose(const Vector &x, Vector &y) const override { Mult(x, y); }
};

// T(λ) as the diagonal of the (lagged) preconditioner, and its exact inverse times a scale.
class DiagonalComplexOperator : public ComplexOperator
{
public:
  std::vector<std::complex<double>> d;
  explicit DiagonalComplexOperator(std::vector<std::complex<double>> &&d)
    : ComplexOperator(static_cast<int>(d.size())), d(std::move(d))
  {
  }
  void Mult(const ComplexVector &x, ComplexVector &y) const override
  {
    const auto *xr = x.Real().HostRead();
    const auto *xi = x.Imag().HostRead();
    auto *yr = y.Real().HostWrite();
    auto *yi = y.Imag().HostWrite();
    for (std::size_t i = 0; i < d.size(); i++)
    {
      const auto v = d[i] * std::complex<double>(xr[i], xi[i]);
      yr[i] = v.real();
      yi[i] = v.imag();
    }
  }
};

class ScaledDiagonalInverse : public Solver<ComplexOperator>
{
  double scale;
  std::vector<std::complex<double>> inv;

public:
  explicit ScaledDiagonalInverse(double scale) : scale(scale) {}
  void SetOperator(const ComplexOperator &op) override
  {
    const auto *diag = dynamic_cast<const DiagonalComplexOperator *>(&op);
    MFEM_VERIFY(diag, "Unexpected preconditioner operator!");
    height = width = op.Height();
    inv.resize(diag->d.size());
    for (std::size_t i = 0; i < inv.size(); i++)
    {
      inv[i] = scale / diag->d[i];
    }
  }
  void Mult(const ComplexVector &x, ComplexVector &y) const override
  {
    DiagonalComplexOperator(std::vector<std::complex<double>>(inv)).Mult(x, y);
  }
};

// The seed eigenpairs of the linear eigensolver which precedes Quasi-Newton (the
// eigenmode driver passes the solved seed eigensolver; only its results are used).
class SeedEigenpairs : public EigenvalueSolver
{
  std::vector<std::complex<double>> lambda;
  std::vector<ComplexVector> x;

public:
  SeedEigenpairs(std::vector<std::complex<double>> &&lambda, std::vector<ComplexVector> &&x)
    : lambda(std::move(lambda)), x(std::move(x))
  {
  }
  void SetLinearSolver(ComplexKspSolver &) override {}
  void SetDivFreeProjector(const DivFreeSolver<ComplexVector> &) override {}
  void SetBMat(const Operator &) override {}
  double GetScalingGamma() const override { return 1.0; }
  double GetScalingDelta() const override { return 1.0; }
  void SetNumModes(int, int) override {}
  void SetTol(double) override {}
  void SetMaxIter(int) override {}
  void SetWhichEigenpairs(WhichType) override {}
  void SetShiftInvert(std::complex<double>, bool) override {}
  void SetInitialSpace(const ComplexVector &) override {}
  int Solve() override { return static_cast<int>(lambda.size()); }
  std::complex<double> GetEigenvalue(int i) const override { return lambda.at(i); }
  void GetEigenvector(int i, ComplexVector &y) const override { y = x.at(i); }
  double GetError(int, ErrorType) const override { return 0.0; }
  void RescaleEigenvectors(int) override {}
};

struct DiagonalLorentzProblem
{
  std::vector<double> k;  // Global, by element
  std::unique_ptr<Mesh> mesh;
  mfem::L2_FECollection fec;
  std::unique_ptr<FiniteElementSpace> fespace;
  std::vector<long long> gid;  // Global index of each local dof
  std::unique_ptr<ComplexParOperator> K, M;

  explicit DiagonalLorentzProblem(std::vector<double> k_) : k(std::move(k_)), fec(0, 2)
  {
    auto smesh = mfem::Mesh::MakeCartesian2D(static_cast<int>(k.size()), 1,
                                             mfem::Element::QUADRILATERAL, false,
                                             static_cast<double>(k.size()), 1.0);
    mesh = std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(Mpi::World(), smesh));
    fespace = std::make_unique<FiniteElementSpace>(*mesh, &fec);
    const int n = mesh->Get().GetNE();
    REQUIRE(fespace->GetVSize() == n);
    Vector kd(n), md(n);
    for (int i = 0; i < n; i++)
    {
      gid.push_back(mesh->Get().GetGlobalElementNum(i));
      kd(i) = k.at(gid.back());
      md(i) = 1.0;
    }
    K = std::make_unique<ComplexParOperator>(
        std::make_unique<DiagonalOperator>(std::move(kd)), nullptr, *fespace);
    M = std::make_unique<ComplexParOperator>(
        std::make_unique<DiagonalOperator>(std::move(md)), nullptr, *fespace);
  }

  int LocalSize() const { return static_cast<int>(gid.size()); }

  // A2(λ) = λ² χ(λ) I.
  std::unique_ptr<ComplexOperator> A2(std::complex<double> l) const
  {
    const auto a = l * l * Chi(l);
    Vector ar(LocalSize()), ai(LocalSize());
    ar = a.real();
    ai = a.imag();
    return std::make_unique<ComplexParOperator>(
        std::make_unique<DiagonalOperator>(std::move(ar)),
        std::make_unique<DiagonalOperator>(std::move(ai)), *fespace);
  }

  // The preconditioner matrix a0 K + a2 M + A2(λ) at λ = i ω.
  std::unique_ptr<ComplexOperator> Preconditioner(std::complex<double> a0,
                                                  std::complex<double> a2,
                                                  std::complex<double> omega) const
  {
    const auto l = std::complex<double>(0.0, 1.0) * omega;
    std::vector<std::complex<double>> d(LocalSize());
    for (int i = 0; i < LocalSize(); i++)
    {
      d[i] = a0 * k.at(gid[i]) + a2 + l * l * Chi(l);
    }
    return std::make_unique<DiagonalComplexOperator>(std::move(d));
  }

  // The unit vector of global dof j, scaled by s.
  ComplexVector Unit(long long j, std::complex<double> s = 1.0) const
  {
    ComplexVector x(LocalSize());
    x = 0.0;
    for (int i = 0; i < LocalSize(); i++)
    {
      if (gid[i] == j)
      {
        x.Real()(i) = s.real();
        x.Imag()(i) = s.imag();
      }
    }
    return x;
  }

  // ||T(λ) x|| / (||x|| (max k + |λ|² |1 + χ(λ)|)), computed directly.
  double BackwardError(std::complex<double> l, const ComplexVector &x) const
  {
    double rr = 0.0, xx = 0.0;
    const auto *xr = x.Real().HostRead();
    const auto *xi = x.Imag().HostRead();
    for (int i = 0; i < LocalSize(); i++)
    {
      const std::complex<double> xv(xr[i], xi[i]);
      rr += std::norm((k.at(gid[i]) + l * l * (1.0 + Chi(l))) * xv);
      xx += std::norm(xv);
    }
    Mpi::GlobalSum(1, &rr, Mpi::World());
    Mpi::GlobalSum(1, &xx, Mpi::World());
    const double scale =
        *std::ranges::max_element(k) + std::norm(l) * std::abs(1.0 + Chi(l));
    return std::sqrt(rr / xx) / scale;
  }

  // The exact roots of dof j, with Re λ ≤ 0 and Im λ > 0 (ω = λ / i with Re ω > 0).
  std::vector<std::complex<double>> Roots(double kj) const
  {
    Eigen::Matrix4cd companion = Eigen::Matrix4cd::Zero();
    companion(0, 0) = -gamma_l;
    companion(0, 1) = -(w0 * w0 * (1.0 + delta_eps) + kj);
    companion(0, 2) = -kj * gamma_l;
    companion(0, 3) = -kj * w0 * w0;
    companion(1, 0) = companion(2, 1) = companion(3, 2) = 1.0;
    const Eigen::ComplexEigenSolver<Eigen::Matrix4cd> roots(companion);
    REQUIRE(roots.info() == Eigen::Success);
    std::vector<std::complex<double>> out;
    for (const auto l : roots.eigenvalues())
    {
      if (l.imag() > 0.0)
      {
        out.push_back(l);
      }
    }
    return out;
  }
};

struct QuasiNewtonResult
{
  std::vector<std::complex<double>> lambda;
  std::vector<ComplexVector> x;
  std::vector<double> backward_error;
};

// Quasi-Newton as the eigenmode driver configures it, from the seeds θ_j with the unit
// eigenvectors of their dofs, at the target ω = target (σ = i target).
QuasiNewtonResult
SolveQuasiNewton(const DiagonalLorentzProblem &p,
                 const std::vector<std::pair<std::complex<double>, int>> &seeds, int nev,
                 double target, double pc_scale, double tol)
{
  std::vector<std::complex<double>> theta;
  std::vector<ComplexVector> vecs;
  for (const auto &[t, j] : seeds)
  {
    theta.push_back(t);
    vecs.push_back(p.Unit(j));
  }
  const int num_seeds = static_cast<int>(seeds.size());
  QuasiNewtonSolver qn(Mpi::World(),
                       std::make_unique<SeedEigenpairs>(std::move(theta), std::move(vecs)),
                       num_seeds, 0, true);
  qn.SetTol(tol);
  qn.SetMaxIter(100);
  qn.SetOperators(*p.K, *p.M, EigenvalueSolver::ScaleType::NONE);
  qn.SetExtraSystemMatrix([&p](std::complex<double> l) { return p.A2(l); });
  qn.SetExtraSystemPoles(
      {std::complex<double>(-0.5 * gamma_l, std::sqrt(w0 * w0 - 0.25 * gamma_l * gamma_l)),
       std::complex<double>(-0.5 * gamma_l,
                            -std::sqrt(w0 * w0 - 0.25 * gamma_l * gamma_l))});
  qn.SetPreconditionerUpdate([&p](std::complex<double> a0, std::complex<double>,
                                  std::complex<double> a2, std::complex<double> a3)
                             { return p.Preconditioner(a0, a2, a3); });
  qn.SetNumModes(nev, nev);
  qn.SetPreconditionerLag(10, 1.0e-4);
  qn.SetMaxRestart(2);
  auto gmres = std::make_unique<GmresSolver<ComplexOperator>>(Mpi::World(), 0);
  gmres->SetRelTol(1.0e-12);
  gmres->SetMaxIter(100);
  ComplexKspSolver ksp(std::move(gmres), std::make_unique<ScaledDiagonalInverse>(pc_scale));
  qn.SetLinearSolver(ksp);
  qn.SetShiftInvert(std::complex<double>(0.0, target));
  const int nconv = qn.Solve();
  QuasiNewtonResult out;
  for (int i = 0; i < nconv; i++)
  {
    out.lambda.push_back(qn.GetEigenvalue(i));
    auto &x = out.x.emplace_back(p.LocalSize());
    qn.GetEigenvector(i, x);
    out.backward_error.push_back(qn.GetError(i, EigenvalueSolver::ErrorType::BACKWARD));
  }
  return out;
}

// The roots of every dof between the target and the pole, sorted by frequency.
std::vector<std::complex<double>> TargetSideRoots(const DiagonalLorentzProblem &p,
                                                  double target)
{
  std::vector<std::complex<double>> out;
  for (const double kj : p.k)
  {
    for (const auto l : p.Roots(kj))
    {
      if (l.imag() > target && l.imag() < w0)
      {
        out.push_back(l);
      }
    }
  }
  std::ranges::sort(out, [](auto a, auto b) { return a.imag() < b.imag(); });
  return out;
}

void CheckEigenpairs(const DiagonalLorentzProblem &p, const QuasiNewtonResult &result,
                     const std::vector<std::complex<double>> &expected, double tol)
{
  REQUIRE(result.lambda.size() == expected.size());
  for (std::size_t i = 0; i < expected.size(); i++)
  {
    CAPTURE(i, result.lambda[i], expected[i], result.backward_error[i]);
    CHECK(std::abs(result.lambda[i] - expected[i]) <= 1.0e-8 * std::abs(expected[i]));
    CHECK(result.backward_error[i] <= tol);
    // The backward error of the returned pair, computed independently of the solver.
    CHECK(p.BackwardError(result.lambda[i], result.x[i]) <= 10.0 * tol);
  }
}

}  // namespace

// Each seed lies above the pole, where the polynomial seed of a mode shape in a Lorentz
// medium falls, and the target below it: the requested modes are the roots between the
// target and the pole, found from the target's side of the nonlinear Rayleigh functional,
// not the roots of the seeds' side.
TEST_CASE("Quasi-Newton finds the roots on the target's side of a pole",
          "[nleps][Serial][Parallel]")
{
  const DiagonalLorentzProblem p({1.5 * 1.5, 1.8 * 1.8, 2.1 * 2.1});
  constexpr double target = 0.3, tol = 1.0e-10;
  std::vector<std::pair<std::complex<double>, int>> seeds;
  for (int j = 0; j < 3; j++)
  {
    seeds.emplace_back(std::complex<double>(0.0, std::sqrt(p.k[j])), j);
  }
  const auto expected = TargetSideRoots(p, target);
  REQUIRE(expected.size() == 3);

  SECTION("Exact preconditioner")
  {
    CheckEigenpairs(p, SolveQuasiNewton(p, seeds, 3, target, 1.0, tol), expected, tol);
  }

  SECTION("Scaled preconditioner")
  {
    // The preconditioner scaled by 1e-6 is as good a preconditioner, but the
    // preconditioned residual of the Newton correction is 1e-6 times smaller. With an
    // absolute tolerance on it, as the inner solves had with 1e-12, GMRES returns a zero
    // correction above the eigensolver tolerance and the iteration stalls.
    CheckEigenpairs(p, SolveQuasiNewton(p, seeds, 3, target, 1.0e-6, tol), expected, tol);
  }
}

// A degenerate pair: after the first mode of the pair is locked, the deflated iterate of
// the second can collapse into the span of the first, with a small residual for the unit
// extended vector but not for the stored eigenvector. Both modes must be returned as
// eigenpairs with independent, actual backward errors below the tolerance.
namespace
{

QuasiNewtonResult SolveDegeneratePair(const DiagonalLorentzProblem &p, double target,
                                      double tol)
{
  std::vector<std::pair<std::complex<double>, int>> seeds;
  for (int j = 0; j < 3; j++)
  {
    seeds.emplace_back(std::complex<double>(0.0, std::sqrt(p.k[j])), j);
  }
  return SolveQuasiNewton(p, seeds, 3, target, 1.0, tol);
}

}  // namespace

TEST_CASE("Quasi-Newton deflation of a degenerate pair", "[nleps][Serial][Parallel]")
{
  const DiagonalLorentzProblem p({1.5 * 1.5, 1.5 * 1.5, 2.1 * 2.1});
  constexpr double target = 0.3, tol = 1.0e-10;
  const auto expected = TargetSideRoots(p, target);
  REQUIRE(expected.size() == 3);
  CheckEigenpairs(p, SolveDegeneratePair(p, target, tol), expected, tol);
}

// Known limitation: both modes of an exactly degenerate pair are returned with the same
// eigenvector (normalized overlap 1), although their eigenvalues and backward errors are
// correct. The eigenvectors are extracted from the invariant pair (X, H) through the
// eigenvectors of H, as on main; a repeated eigenvalue of H with the deflation coupling
// would give a single eigenvector, which is the likely cause. [!shouldfail] reports this
// test as passing while it fails, and as failing once the limitation is fixed.
TEST_CASE("Quasi-Newton eigenvectors of a degenerate pair are independent",
          "[nleps][Serial][Parallel][!shouldfail]")
{
  const DiagonalLorentzProblem p({1.5 * 1.5, 1.5 * 1.5, 2.1 * 2.1});
  constexpr double target = 0.3, tol = 1.0e-10;
  const auto result = SolveDegeneratePair(p, target, tol);
  REQUIRE(result.x.size() == 3);
  const auto overlap = std::abs(linalg::Dot(Mpi::World(), result.x[0], result.x[1])) /
                       (linalg::Norml2(Mpi::World(), result.x[0]) *
                        linalg::Norml2(Mpi::World(), result.x[1]));
  CHECK(overlap < 1.0 - 1.0e-6);
}
