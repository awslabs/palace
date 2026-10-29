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
#if defined(PALACE_WITH_SLEPC)
#include "linalg/slepc.hpp"
#endif
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

using namespace palace;

// The Newton step on the pole-cleared d(λ) T(λ) uses L = d′/d, and the line search compares
// residuals scaled by |d|.
// Decision coverage of the locking criterion of Quasi-Newton (not an end-to-end trajectory
// in which the iterate collapses): the extended residual is compared with the tolerance
// relative to the norm of the physical column, which an absolute threshold does not do.
TEST_CASE("Quasi-Newton locks a column relative to its norm", "[nleps][Serial]")
{
  constexpr double rtol = 1.0e-8;
  CHECK_FALSE(nleps::AcceptsInvariantColumn(1.0e-10, 1.0e-6, rtol));  // 1e-10 < rtol
  CHECK(nleps::AcceptsInvariantColumn(1.0e-16, 1.0e-6, rtol));
  CHECK_FALSE(nleps::AcceptsInvariantColumn(rtol * 1.0e-6, 1.0e-6, rtol));  // Strict
  CHECK(nleps::AcceptsInvariantColumn(0.5e-8, 1.0, rtol));
  CHECK_FALSE(nleps::AcceptsInvariantColumn(2.0e-8, 1.0, rtol));
  for (double c : {1.0e-12, 1.0e-3, 1.0e3})
  {
    // The decision depends on res / ||v|| only.
    CAPTURE(c);
    CHECK_FALSE(nleps::AcceptsInvariantColumn(c * 1.0e-10, c * 1.0e-6, rtol));
    CHECK(nleps::AcceptsInvariantColumn(c * 1.0e-16, c * 1.0e-6, rtol));
  }
  CHECK_FALSE(nleps::AcceptsInvariantColumn(0.0, 0.0, rtol));  // A zero column
}

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

// The diagonal plus c at row i0, column i1 (local indices; no coupling if either is -1).
class CoupledDiagonalOperator : public Operator
{
  Vector d;
  int i0, i1;
  double c;

public:
  CoupledDiagonalOperator(Vector &&d, int i0, int i1, double c)
    : Operator(d.Size()), d(std::move(d)), i0(i0), i1(i1), c(c)
  {
  }
  void Mult(const Vector &x, Vector &y) const override
  {
    const auto *hd = d.HostRead();
    const auto *hx = x.HostRead();
    auto *hy = y.HostWrite();
    for (int i = 0; i < d.Size(); i++)
    {
      hy[i] = hd[i] * hx[i];
    }
    if (i0 >= 0 && i1 >= 0)
    {
      hy[i0] += c * hx[i1];
    }
  }
};

struct DiagonalLorentzProblem
{
  std::vector<double> k;  // Global, by element
  double coupling;        // K(0, 1), one process only

  std::unique_ptr<Mesh> mesh;
  mfem::L2_FECollection fec;
  std::unique_ptr<FiniteElementSpace> fespace;
  std::vector<long long> gid;  // Global index of each local dof
  std::unique_ptr<ComplexParOperator> K, M;

  explicit DiagonalLorentzProblem(std::vector<double> k_, double coupling = 0.0)
    : k(std::move(k_)), coupling(coupling), fec(0, 2)
  {
    REQUIRE((coupling == 0.0 || Mpi::Size(Mpi::World()) == 1));
    auto smesh = mfem::Mesh::MakeCartesian2D(static_cast<int>(k.size()), 1,
                                             mfem::Element::QUADRILATERAL, false,
                                             static_cast<double>(k.size()), 1.0);
    mesh = std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(Mpi::World(), smesh));
    fespace = std::make_unique<FiniteElementSpace>(*mesh, &fec);
    const int n = mesh->Get().GetNE();
    REQUIRE(fespace->GetVSize() == n);
    // The global number of local element i is the offset of this process plus i. The
    // offset comes from a collective scan on every process: ParMesh::GetGlobalElementNum
    // computes it on first use, which a process without elements never makes, so with an
    // empty process the scan was unmatched and later collectives returned wrong numbers.
    long long num_local = n, offset = 0;
    MPI_Exscan(&num_local, &offset, 1, MPI_LONG_LONG, MPI_SUM, Mpi::World());
    if (Mpi::Rank(Mpi::World()) == 0)
    {
      offset = 0;  // MPI_Exscan leaves the result on the first process undefined
    }
    Vector kd(n), md(n);
    for (int i = 0; i < n; i++)
    {
      gid.push_back(offset + i);
      kd(i) = k.at(gid.back());
      md(i) = 1.0;
    }
    K = std::make_unique<ComplexParOperator>(
        std::make_unique<CoupledDiagonalOperator>(std::move(kd), coupling != 0.0 ? 0 : -1,
                                                  coupling != 0.0 ? 1 : -1, coupling),
        nullptr, *fespace);
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
      std::complex<double> tx = (k.at(gid[i]) + l * l * (1.0 + Chi(l))) * xv;
      if (coupling != 0.0 && gid[i] == 0)
      {
        tx += coupling * std::complex<double>(xr[1], xi[1]);  // One process: local = global
      }
      rr += std::norm(tx);
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
                 double target, double pc_scale, double tol, double mix = 0.0)
{
  // Each seed vector is the unit vector of its dof plus mix times those of the others.
  std::vector<std::complex<double>> theta;
  std::vector<ComplexVector> vecs;
  for (const auto &[t, j] : seeds)
  {
    theta.push_back(t);
    auto &x = vecs.emplace_back(p.Unit(j));
    for (std::size_t m = 0; m < p.k.size(); m++)
    {
      if (static_cast<int>(m) != j)
      {
        x.Add(mix, p.Unit(static_cast<long long>(m)));
      }
    }
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
    // The same with the preconditioner scaled by 1e-6. (Started from the target's side of
    // the pole, the Rayleigh roots are exact here, so this does not need small Newton
    // corrections; see the next test.)
    CheckEigenpairs(p, SolveQuasiNewton(p, seeds, 3, target, 1.0e-6, tol), expected, tol);
  }
}

// Seeds 2% below the roots on the target's side of the pole, with 5% of the other dofs in
// their vectors, so that Quasi-Newton takes several steps with vector corrections, with a
// preconditioner scaled by 1e-10: as good a preconditioner, but the preconditioned residual
// of the Newton correction is 1e-10 times smaller. With an absolute tolerance on it, as the
// inner solves had with 1e-12, GMRES returns a zero correction before the eigenpair
// residual reaches the tolerance, and the iteration stalls.
TEST_CASE("Quasi-Newton inner solves have no absolute tolerance",
          "[nleps][Serial][Parallel]")
{
  const DiagonalLorentzProblem p({1.5 * 1.5, 1.8 * 1.8, 2.1 * 2.1});
  constexpr double target = 0.3, tol = 1.0e-10;
  const auto expected = TargetSideRoots(p, target);
  REQUIRE(expected.size() == 3);
  // Root j belongs to dof j: the roots on the target's side increase with k.
  std::vector<std::pair<std::complex<double>, int>> seeds;
  for (int j = 0; j < 3; j++)
  {
    seeds.emplace_back(expected[j] * 0.98, j);
  }
  CheckEigenpairs(p, SolveQuasiNewton(p, seeds, 3, target, 1.0e-10, tol, 0.05), expected,
                  tol);
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

// Both modes of an exactly degenerate pair, each an eigenvector, must be returned as two
// independent fields with their own residuals. The second mode was locked with the same
// eigenvalue and an undetermined coupling to the first, a false Jordan block of H, and
// both were returned with the same field (normalized overlap 1).
TEST_CASE("Quasi-Newton eigenvectors of a degenerate pair are independent",
          "[nleps][Serial][Parallel]")
{
  const DiagonalLorentzProblem p({1.5 * 1.5, 1.5 * 1.5, 2.1 * 2.1});
  constexpr double target = 0.3, tol = 1.0e-10;
  const auto result = SolveDegeneratePair(p, target, tol);
  CheckEigenpairs(p, result, TargetSideRoots(p, target), tol);
  REQUIRE(result.x.size() == 3);
  const auto overlap = std::abs(linalg::Dot(Mpi::World(), result.x[0], result.x[1])) /
                       (linalg::Norml2(Mpi::World(), result.x[0]) *
                        linalg::Norml2(Mpi::World(), result.x[1]));
  CHECK(overlap < 1.0 - 1.0e-6);
}

namespace
{

double Overlap(const ComplexVector &x, const ComplexVector &y)
{
  return std::abs(linalg::Dot(Mpi::World(), x, y)) /
         (linalg::Norml2(Mpi::World(), x) * linalg::Norml2(Mpi::World(), y));
}

}  // namespace

// T(λ) = K + λ² (1 + χ(λ)) I with K upper triangular, K(0, 1) = 0.4: the eigenvector of
// the roots of the second dof, [0.4 / (k₁ - k₀), 1, 0], is not orthogonal to that of the
// first, e₀, while the invariant basis X is. Its field is X y for the eigenvector y of H.
TEST_CASE("Quasi-Newton eigenvectors of a nonnormal problem", "[nleps][Serial]")
{
  const DiagonalLorentzProblem p({1.5 * 1.5, 1.8 * 1.8, 2.1 * 2.1}, 0.4);
  constexpr double target = 0.3, tol = 1.0e-10;
  const auto result = SolveDegeneratePair(p, target, tol);
  CheckEigenpairs(p, result, TargetSideRoots(p, target), tol);
  REQUIRE(result.x.size() == 3);
  CAPTURE(Overlap(result.x[1], p.Unit(0)));
  CHECK(Overlap(result.x[1], p.Unit(0)) > 0.1);
}

// A Jordan block: with k₀ = k₁ and K(0, 1) = 0.3, T(λ) = [[0, 0.3], [0, 0]] at the root of
// the first two dofs, which has one eigenvector e₀ (algebraic multiplicity two). Every
// field returned at that root must be e₀, with its residual: there may be more than one
// such (parallel) field, but no independent one. This does not exercise the warning of an
// unrecovered field.
TEST_CASE("Quasi-Newton eigenvectors of a Jordan block", "[nleps][Serial]")
{
  const DiagonalLorentzProblem p({1.5 * 1.5, 1.5 * 1.5, 2.1 * 2.1}, 0.3);
  constexpr double target = 0.3, tol = 1.0e-10;
  const auto result = SolveDegeneratePair(p, target, tol);
  const auto root = TargetSideRoots(p, target).front();
  int found = 0;
  for (std::size_t i = 0; i < result.lambda.size(); i++)
  {
    CAPTURE(i, result.lambda[i], result.backward_error[i]);
    CHECK(p.BackwardError(result.lambda[i], result.x[i]) <= 10.0 * tol);
    if (std::abs(result.lambda[i] - root) <= 1.0e-6 * std::abs(root))
    {
      found++;
      CHECK(Overlap(result.x[i], p.Unit(0)) > 1.0 - 1.0e-6);
    }
  }
  CHECK(found >= 1);
}

// The seed solver of the Hybrid method, solved from a second start: its eigenpairs are
// those of both solves (here each the root ±1.5i of 2.25 + λ², the one nearest the shift
// first), and it supports one Solve of one problem. One dof: with two processes one has no
// local part of the problem and must still take part in the second solve.
TEST_CASE("Seed solver with a second start", "[nleps][Serial][Parallel]")
{
#if !defined(PALACE_WITH_SLEPC)
  SKIP("Requires SLEPc");
#else
  const DiagonalLorentzProblem p({2.25});
  Vector zero(p.LocalSize());
  zero = 0.0;
  ComplexParOperator C(
      std::make_unique<CoupledDiagonalOperator>(std::move(zero), -1, -1, 0.0), nullptr,
      *p.fespace);
  const std::complex<double> sigma(0.0, 1.4);
  DiagonalComplexOperator P(
      std::vector<std::complex<double>>(p.LocalSize(), 2.25 + sigma * sigma));
  auto gmres = std::make_unique<GmresSolver<ComplexOperator>>(Mpi::World(), 0);
  gmres->SetRelTol(1.0e-12);
  gmres->SetMaxIter(100);
  ComplexKspSolver ksp(std::move(gmres), std::make_unique<ScaledDiagonalInverse>(1.0));
  ksp.SetOperators(P, P);
  auto Configure = [&](slepc::SlepcPEPLinearSolver &eps, const DiagonalLorentzProblem &q)
  {
    eps.SetType(slepc::SlepcEigenvalueSolver::Type::KRYLOVSCHUR);
    eps.SetProblemType(slepc::SlepcEigenvalueSolver::ProblemType::GEN_NON_HERMITIAN);
    eps.SetOperators(*q.K, C, *q.M, EigenvalueSolver::ScaleType::NONE);
    eps.SetLinearSolver(ksp);
    eps.SetShiftInvert(sigma);
    eps.SetWhichEigenpairs(EigenvalueSolver::WhichType::TARGET_MAGNITUDE);
    eps.SetNumModes(1, 2);
    eps.SetTol(1.0e-10);
    eps.SetMaxIter(100);
    eps.SetInitialSpace(q.Unit(0));
  };
  // The second start is set before the problem.
  slepc::SlepcPEPLinearSolver eps(Mpi::World(), 0);
  eps.SetSecondInitialSpace(p.Unit(0, {0.3, 0.7}));
  Configure(eps, p);
  const int num = eps.Solve();
  REQUIRE(num >= 2);
  for (int i = 0; i < num; i++)
  {
    const auto l = eps.GetEigenvalue(i);
    CAPTURE(i, l, eps.GetError(i, EigenvalueSolver::ErrorType::ABSOLUTE));
    CHECK(std::abs(l * l + 2.25) <= 1.0e-8);
    CHECK(eps.GetError(i, EigenvalueSolver::ErrorType::ABSOLUTE) <= 1.0e-8);
  }
  const auto first = eps.GetEigenvalue(0);
  CHECK(std::abs(first - std::complex<double>(0.0, 1.5)) <= 1.0e-8);
  // Another problem (of the same or of another size) or Solve is an error, and leaves the
  // results.
  CHECK_THROWS(eps.SetOperators(*p.K, C, *p.M, EigenvalueSolver::ScaleType::NONE));
  const DiagonalLorentzProblem q({2.25, 4.0});
  CHECK_THROWS(eps.SetOperators(*q.K, *q.M, *q.M, EigenvalueSolver::ScaleType::NONE));
  CHECK_THROWS(eps.Solve());
  CHECK(eps.GetEigenvalue(0) == first);

  // A second start set after the problem is an error, also after a solve of it.
  slepc::SlepcPEPLinearSolver late(Mpi::World(), 0);
  Configure(late, p);
  CHECK_THROWS(late.SetSecondInitialSpace(p.Unit(0, {0.3, 0.7})));
  REQUIRE(late.Solve() >= 1);
  CHECK_THROWS(late.SetSecondInitialSpace(p.Unit(0, {0.3, 0.7})));
  CHECK(std::abs(late.GetEigenvalue(0) - std::complex<double>(0.0, 1.5)) <= 1.0e-8);

  // A second start of another size than the problem.
  slepc::SlepcPEPLinearSolver wrong(Mpi::World(), 0);
  wrong.SetSecondInitialSpace(ComplexVector(p.LocalSize() + 1));
  Configure(wrong, p);
  CHECK_THROWS(wrong.Solve());
#endif
}
