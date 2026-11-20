// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <complex>
#include <iterator>
#include <memory>
#include <regex>
#include <string>
#include <utility>
#include <vector>
#include <Eigen/Dense>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
#include "linalg/eps.hpp"
#if defined(PALACE_WITH_ARPACK)
#include "linalg/arpack.hpp"
#endif
#include "linalg/iterative.hpp"
#include "linalg/ksp.hpp"
#include "linalg/nleps.hpp"
#include "linalg/rap.hpp"
#if defined(PALACE_WITH_SLEPC)
#include "linalg/slepc.hpp"
#endif
#include "linalg/vector.hpp"
#include "test-helpers.hpp"
#include "utils/communication.hpp"

using namespace palace;

// The SLP Jacobian of the pole-cleared d(λ) T(λ) uses L = d′/d.
TEST_CASE("Pole-cleared Newton terms", "[nleps][Serial]")
{
  const std::complex<double> p(-0.1, 1.5);
  const nleps::PoleTerms poles{{p, std::conj(p)}};
  auto d = [&](std::complex<double> l) { return (l - p) * (l - std::conj(p)); };

  const std::complex<double> l(-0.02, 1.7);
  const double h = 1.0e-6;
  const std::complex<double> dd = (d(l + h) - d(l - h)) / (2.0 * h);
  CHECK(std::abs(poles.LogDerivative(l) - dd / d(l)) <= 1.0e-8 * std::abs(dd / d(l)));
  CHECK_FALSE(poles.Empty());
  CHECK(nleps::PoleTerms{}.Empty());
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
  void AssembleDiagonal(Vector &diag) const override { diag = d; }
};

// T(λ) as a diagonal matrix, and its exact inverse times a scale.
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

// Append decoupled dofs with the diagonal entry pad to the global vector v of a problem of
// v.size() dofs, so that it has more dofs than processes: MFEM's partitioning then leaves
// no process without a dof, which the eigenvalue solvers require on every process (an
// empty process fails a verification, and the others wait in a collective). With
// min_local > 1 it has twice min_local dofs per process, for a partitioning which is not
// even: PARPACK requires its number of basis vectors ≤ the local size on every process.
// The tests keep these dofs out of their solves: zero in a start vector, and far from the
// target.
void AddDecoupledDofs(std::vector<double> &v, double pad, int min_local = 1)
{
  const auto size = static_cast<std::size_t>(Mpi::Size(Mpi::World()));
  v.resize(std::max(v.size(), (min_local > 1) ? 2 * min_local * size : size + 1), pad);
}

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
        std::make_unique<DiagonalOperator>(std::move(kd)), nullptr, *fespace);
    M = std::make_unique<ComplexParOperator>(
        std::make_unique<DiagonalOperator>(std::move(md)), nullptr, *fespace);
  }

  int LocalSize() const { return static_cast<int>(gid.size()); }

  // The smallest number of dofs of a process.
  int MinLocalSize() const
  {
    int n = LocalSize();
    Mpi::GlobalMin(1, &n, Mpi::World());
    return n;
  }

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

// T(λ) = K + λ² M + λ² χ(λ) B on the dofs of a DiagonalLorentzProblem, with B = diag(b) the
// indicator of the dofs which have the Lorentz medium. λ² χ(λ) B is a rational term with
// the two poles of χ and its residues there.
struct RationalDiagonalProblem
{
  DiagonalLorentzProblem p;
  std::vector<double> b;  // Global, by element
  std::unique_ptr<ComplexParOperator> B;
  std::vector<std::complex<double>> poles, coeffs;

  RationalDiagonalProblem(std::vector<double> k, std::vector<double> b_)
    : p(std::move(k)), b(std::move(b_))
  {
    Vector bd(p.LocalSize());
    for (int i = 0; i < p.LocalSize(); i++)
    {
      bd(i) = b.at(p.gid[i]);
    }
    B = std::make_unique<ComplexParOperator>(
        std::make_unique<DiagonalOperator>(std::move(bd)), nullptr, *p.fespace);
    const std::complex<double> pole(-0.5 * gamma_l,
                                    std::sqrt(w0 * w0 - 0.25 * gamma_l * gamma_l));
    poles = {pole, std::conj(pole)};
    coeffs = {delta_eps * w0 * w0 / (poles[0] - poles[1]),
              delta_eps * w0 * w0 / (poles[1] - poles[0])};
  }

  // The diagonal entry of T(λ) of global dof j, and T(λ).
  std::complex<double> T(long long j, std::complex<double> l) const
  {
    return p.k.at(j) + l * l * (1.0 + b.at(j) * Chi(l));
  }
  std::unique_ptr<DiagonalComplexOperator> Matrix(std::complex<double> l) const
  {
    std::vector<std::complex<double>> d(p.LocalSize());
    for (int i = 0; i < p.LocalSize(); i++)
    {
      d[i] = T(p.gid[i], l);
    }
    return std::make_unique<DiagonalComplexOperator>(std::move(d));
  }

  // The eigenvalues of global dof j with Im λ > 0: two with the medium, one without.
  std::vector<std::complex<double>> Roots(long long j) const
  {
    if (b.at(j) != 0.0)
    {
      return p.Roots(p.k.at(j));
    }
    return {std::complex<double>(0.0, std::sqrt(p.k.at(j)))};
  }

  // ||T(λ) x|| / ||x||, computed directly.
  double Residual(std::complex<double> l, const ComplexVector &x) const
  {
    double rr = 0.0, xx = 0.0;
    const auto *xr = x.Real().HostRead();
    const auto *xi = x.Imag().HostRead();
    for (int i = 0; i < p.LocalSize(); i++)
    {
      const std::complex<double> xv(xr[i], xi[i]);
      rr += std::norm(T(p.gid[i], l) * xv);
      xx += std::norm(xv);
    }
    Mpi::GlobalSum(1, &rr, Mpi::World());
    Mpi::GlobalSum(1, &xx, Mpi::World());
    return std::sqrt(rr / xx);
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

}  // namespace

// The eigensolvers which linearize the quadratic problem solve a problem with rational
// terms as a linear one: the eigenvalues nearest the shift are the exact roots, on both
// sides of the poles (|p±| = 1), with one of a dof without the medium among them.
TEST_CASE("Linearized eigensolvers solve rational terms exactly",
          "[nleps][Serial][Parallel]")
{
  std::vector<double> k = {1.5,  2.5,  4.0,  6.5,  3.0,  12.0,
                           16.0, 20.0, 25.0, 30.0, 36.0, 42.0},
                      b = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  // The vectors of ARPACK have 2 + 2 blocks (the poles) of the dofs.
  constexpr int nev = 9, ncv = 24, min_local = ncv / 4;
  const auto num_coupled = static_cast<long long>(k.size());
  AddDecoupledDofs(k, 1.0e4, min_local);
  AddDecoupledDofs(b, 0.0, min_local);
  const RationalDiagonalProblem q(k, b);
  const std::complex<double> sigma(0.0, 1.6);
  auto Nearer = [sigma](auto a, auto b)
  { return std::abs(a - sigma) < std::abs(b - sigma); };
  std::vector<std::complex<double>> expected;
  for (std::size_t j = 0; j < q.b.size(); j++)
  {
    std::ranges::copy(q.Roots(j), std::back_inserter(expected));
  }
  std::ranges::sort(expected, Nearer);
  REQUIRE(expected.size() > nev);
  // The requested eigenvalues are separated from the next one, and from the eigenvalues
  // with Im λ < 0, which are farther from the shift than |σ|.
  REQUIRE(std::abs(expected[nev] - sigma) > 1.1 * std::abs(expected[nev - 1] - sigma));
  REQUIRE(std::abs(expected[nev] - sigma) < std::abs(sigma));

  auto P = q.Matrix(sigma);
  auto gmres = std::make_unique<GmresSolver<ComplexOperator>>(Mpi::World(), 0);
  gmres->SetRelTol(1.0e-14);
  gmres->SetMaxIter(100);
  ComplexKspSolver ksp(std::move(gmres), std::make_unique<ScaledDiagonalInverse>(1.0));
  ksp.SetOperators(*P, *P);
  ComplexVector start(q.p.LocalSize());
  {
    auto *vr = start.Real().HostWrite();
    auto *vi = start.Imag().HostWrite();
    for (int i = 0; i < q.p.LocalSize(); i++)
    {
      const bool coupled = (q.p.gid[i] < num_coupled);
      vr[i] = coupled ? 1.0 + 0.1 * q.p.gid[i] : 0.0;
      vi[i] = coupled ? 0.3 - 0.05 * q.p.gid[i] : 0.0;
    }
  }
  const auto scale =
      GENERATE(EigenvalueSolver::ScaleType::NONE, EigenvalueSolver::ScaleType::NORM_2);
  auto Check = [&](EigenvalueSolver &eigen)
  {
    eigen.AddPoleTerms(*q.B, q.poles, q.coeffs);
    eigen.SetOperators(*q.p.K, *q.p.M, scale);
    eigen.SetLinearSolver(ksp);
    eigen.SetShiftInvert(sigma);
    eigen.SetWhichEigenpairs(EigenvalueSolver::WhichType::TARGET_MAGNITUDE);
    eigen.SetNumModes(nev, ncv);
    eigen.SetTol(1.0e-10);
    eigen.SetMaxIter(200);
    eigen.SetInitialSpace(start);
    const int num = eigen.Solve();
    REQUIRE(num >= nev);
    std::vector<std::complex<double>> found;
    ComplexVector x(q.p.LocalSize());
    for (int i = 0; i < nev; i++)
    {
      const auto l = eigen.GetEigenvalue(i);
      found.push_back(l);
      eigen.GetEigenvector(i, x);
      CAPTURE(i, l);
      CHECK(q.Residual(l, x) <= 1.0e-8);
      CHECK(eigen.GetError(i, EigenvalueSolver::ErrorType::BACKWARD) <= 1.0e-8);
    }
    std::ranges::sort(found, Nearer);
    for (int i = 0; i < nev; i++)
    {
      CAPTURE(i, found[i], expected[i]);
      CHECK(std::abs(found[i] - expected[i]) <= 1.0e-8);
    }
  };
#if defined(PALACE_WITH_SLEPC)
  SECTION("SLEPc")
  {
    slepc::SlepcPEPLinearSolver eigen(Mpi::World(), 0);
    eigen.SetType(slepc::SlepcEigenvalueSolver::Type::KRYLOVSCHUR);
    eigen.SetProblemType(slepc::SlepcEigenvalueSolver::ProblemType::GEN_NON_HERMITIAN);
    Check(eigen);
  }
#endif
#if defined(PALACE_WITH_ARPACK)
  SECTION("ARPACK")
  {
    REQUIRE(q.p.MinLocalSize() >= min_local);
    arpack::ArpackPEPSolver eigen(Mpi::World(), 0);
    Check(eigen);
  }
#endif
}

// ARPACK reaching its iteration limit with fewer eigenvalues converged than requested
// (pznaupd info = 1) is not an error: it returns the number converged and those
// eigenpairs, as SLEPc does, and the caller decides whether that is enough. With enough
// iterations the same solve converges all of them.
TEST_CASE("ARPACK returns the eigenpairs converged within its iteration limit",
          "[nleps][Serial][Parallel]")
{
#if !defined(PALACE_WITH_ARPACK)
  SKIP("Requires the ARPACK eigensolver");
#else
  std::vector<double> k(60);
  for (std::size_t j = 0; j < k.size(); j++)
  {
    k[j] = 1.0 + j + 0.37 * std::sin(1.0 + j);  // Distinct
  }
  constexpr int nev = 8, ncv = nev + 2;
  const auto num_coupled = static_cast<long long>(k.size());
  AddDecoupledDofs(k, 1.0e4, ncv);
  const DiagonalLorentzProblem p(k);
  REQUIRE(p.MinLocalSize() >= ncv);
  const double sigma = 10.3;
  std::vector<std::complex<double>> d(p.LocalSize());
  ComplexVector start(p.LocalSize());
  for (int i = 0; i < p.LocalSize(); i++)
  {
    const bool coupled = (p.gid[i] < num_coupled);
    d[i] = p.k.at(p.gid[i]) - sigma;
    start.Real()(i) = coupled ? 1.0 + 0.1 * p.gid[i] : 0.0;
    start.Imag()(i) = coupled ? 0.3 - 0.05 * p.gid[i] : 0.0;
  }
  const DiagonalComplexOperator A(std::move(d));
  auto gmres = std::make_unique<GmresSolver<ComplexOperator>>(Mpi::World(), 0);
  gmres->SetRelTol(1.0e-14);
  gmres->SetMaxIter(100);
  ComplexKspSolver ksp(std::move(gmres), std::make_unique<ScaledDiagonalInverse>(1.0));
  ksp.SetOperators(A, A);
  auto Solve = [&](int max_it)
  {
    arpack::ArpackEPSSolver eigen(Mpi::World(), 0);
    eigen.SetOperators(*p.K, *p.M, EigenvalueSolver::ScaleType::NONE);
    eigen.SetLinearSolver(ksp);
    eigen.SetShiftInvert(sigma);
    eigen.SetWhichEigenpairs(EigenvalueSolver::WhichType::LARGEST_MAGNITUDE);
    eigen.SetNumModes(nev, ncv);
    eigen.SetTol(1.0e-10);
    eigen.SetMaxIter(max_it);
    eigen.SetInitialSpace(start);
    const int num = eigen.Solve();
    for (int i = 0; i < std::min(num, nev); i++)
    {
      const auto l = eigen.GetEigenvalue(i);
      CAPTURE(max_it, i, l);
      CHECK(std::ranges::any_of(k, [l](double kj)
                                { return std::abs(l - kj) <= 1.0e-8 * kj; }));
      CHECK(eigen.GetError(i, EigenvalueSolver::ErrorType::BACKWARD) <= 1.0e-8);
    }
    return num;
  };
  int num = 0;
  CHECK_NOTHROW(num = Solve(1));
  CHECK(num < nev);
  CHECK(Solve(300) >= nev);
#endif
}

#if defined(PALACE_WITH_ARPACK)
// ARPACK with the permutation of its eigenpairs set directly: the target-side compaction of
// a TARGET_IMAGINARY solve keeps the eigenpairs at some of its raw indices, which a solve
// cannot be made to produce in a given order.
class ArpackEPSSolverKept : public arpack::ArpackEPSSolver
{
public:
  using arpack::ArpackEPSSolver::ArpackEPSSolver;
  void Keep(const std::vector<int> &raw) { std::ranges::copy(raw, perm.get()); }
};
#endif

// After a target-side compaction which discards a raw index before a kept one (raw index 0
// of 4 here), the kept eigenpair i is at a raw index perm[i] which can be past the number
// kept. RescaleEigenvectors (as the driver calls it after SetBMat) must scale and measure
// each kept eigenvector at its raw index.
TEST_CASE("ARPACK rescales the eigenpairs kept at any raw index",
          "[nleps][Serial][Parallel]")
{
#if !defined(PALACE_WITH_ARPACK)
  SKIP("Requires the ARPACK eigensolver");
#else
  std::vector<double> k(60);
  for (std::size_t j = 0; j < k.size(); j++)
  {
    k[j] = 1.0 + j + 0.37 * std::sin(1.0 + j);  // Distinct
  }
  constexpr int nev = 4, ncv = 2 * nev + 2;
  AddDecoupledDofs(k, 1.0e4, ncv);  // Far from the target, in the random start
  const DiagonalLorentzProblem p(k);
  REQUIRE(p.MinLocalSize() >= ncv);
  const double sigma = 10.3;
  std::vector<std::complex<double>> d(p.LocalSize());
  Vector b(p.LocalSize());
  b.UseDevice(true);
  for (int i = 0; i < p.LocalSize(); i++)
  {
    d[i] = p.k.at(p.gid[i]) - sigma;
    b(i) = 1.0 + 0.5 * p.gid[i];  // A B norm which is not that of ARPACK's basis
  }
  const DiagonalComplexOperator A(std::move(d));
  const palace::DiagonalOperator B(b);
  auto gmres = std::make_unique<GmresSolver<ComplexOperator>>(Mpi::World(), 0);
  gmres->SetRelTol(1.0e-14);
  gmres->SetMaxIter(100);
  ComplexKspSolver ksp(std::move(gmres), std::make_unique<ScaledDiagonalInverse>(1.0));
  ksp.SetOperators(A, A);
  ArpackEPSSolverKept eigen(Mpi::World(), 0);
  eigen.SetOperators(*p.K, *p.M, EigenvalueSolver::ScaleType::NONE);
  eigen.SetLinearSolver(ksp);
  eigen.SetShiftInvert(sigma);
  eigen.SetWhichEigenpairs(EigenvalueSolver::WhichType::LARGEST_MAGNITUDE);
  eigen.SetNumModes(nev, ncv);
  eigen.SetTol(1.0e-10);
  eigen.SetMaxIter(300);
  REQUIRE(eigen.Solve() >= nev);
  eigen.SetBMat(B);
  const std::vector<int> kept = {1, 2, 3};
  eigen.Keep(kept);
  eigen.RescaleEigenvectors(static_cast<int>(kept.size()));
  ComplexVector x(p.LocalSize()), Bx(p.LocalSize());
  x.UseDevice(true);
  Bx.UseDevice(true);
  for (std::size_t i = 0; i < kept.size(); i++)
  {
    const auto l = eigen.GetEigenvalue(static_cast<int>(i));
    CAPTURE(i, l);
    CHECK(
        std::ranges::any_of(k, [l](double kj) { return std::abs(l - kj) <= 1.0e-8 * kj; }));
    eigen.GetEigenvector(static_cast<int>(i), x);
    CHECK(std::abs(linalg::Norml2(Mpi::World(), x, B, Bx) - 1.0) <= 1.0e-12);
    CHECK(eigen.GetError(static_cast<int>(i), EigenvalueSolver::ErrorType::BACKWARD) <=
          1.0e-8);
  }
#endif
}

// Quasi-Newton corrects its iterate with the inverse of T(σ) at σ = λ₀, its eigenvalue
// estimate at the last preconditioner update. Below PreconditionerLagTol the update was
// skipped, and the iteration then converges linearly with the stale operator. Here the
// seed is the eigenvector of dof 0 plus 1e-5 of the others, and its eigenvalue estimate λ₀
// is off the root λ* of dof 0 by Δ / t₀′(λ*), Δ = 0.4. The eigenvalue converges within the
// first lag window, where the residual falls below 1e-4. The component of dof 1, whose
// t₁(λ*) = δk = Δ / 19 is small, then contracts by Δ / (δk + Δ) = 0.95 per step: with the
// operator of λ₀ it needs about 150 iterations to the tolerance 1e-10, more than the 100
// of one attempt. The update at the check where that rate is too slow converges the
// eigenpair in about 22. On more processes than dofs, decoupled dofs, zero in the seed,
// keep the iteration that of the three dofs.
TEST_CASE("Quasi-Newton updates a stale preconditioner below PreconditionerLagTol",
          "[nleps][Serial][Parallel]")
{
  constexpr double k0 = 4.0, Delta = 0.4, rho = 0.95, tol = 1.0e-10;
  std::vector<double> k = {k0, k0 + Delta * (1.0 - rho) / rho, 9.0};
  const auto num_coupled = static_cast<long long>(k.size());
  AddDecoupledDofs(k, 1.0e4);
  const DiagonalLorentzProblem p(k);
  const auto roots = p.Roots(k0);
  const auto ls =
      *std::ranges::max_element(roots, [](auto a, auto b) { return a.imag() < b.imag(); });
  // t₀′(λ) = 2λ (1 + χ(λ)) + λ² χ′(λ).
  const auto den = ls * ls + gamma_l * ls + w0 * w0;
  const auto dt = 2.0 * ls * (1.0 + Chi(ls)) -
                  ls * ls * delta_eps * w0 * w0 * (2.0 * ls + gamma_l) / (den * den);
  const auto l0 = ls + Delta / dt;
  ComplexVector x0(p.LocalSize());
  x0 = 0.0;
  for (int i = 0; i < p.LocalSize(); i++)
  {
    x0.Real()(i) = (p.gid[i] == 0) ? 1.0 : (p.gid[i] < num_coupled) ? 1.0e-5 : 0.0;
  }
  std::vector<std::complex<double>> theta = {l0};
  std::vector<ComplexVector> vecs = {x0};
  QuasiNewtonSolver qn(Mpi::World(),
                       std::make_unique<SeedEigenpairs>(std::move(theta), std::move(vecs)),
                       1, 1, true);
  qn.SetTol(tol);
  qn.SetMaxIter(100);
  qn.SetOperators(*p.K, *p.M, EigenvalueSolver::ScaleType::NONE);
  qn.SetExtraSystemMatrix([&p](std::complex<double> l) { return p.A2(l); });
  qn.SetPreconditionerUpdate([&p](std::complex<double> a0, std::complex<double>,
                                  std::complex<double> a2, std::complex<double> a3)
                             { return p.Preconditioner(a0, a2, a3); });
  qn.SetNumModes(1, 1);
  qn.SetPreconditionerLag(10, 1.0e-4);
  qn.SetMaxRestart(0);
  auto gmres = std::make_unique<GmresSolver<ComplexOperator>>(Mpi::World(), 0);
  gmres->SetRelTol(1.0e-12);
  gmres->SetMaxIter(100);
  ComplexKspSolver ksp(std::move(gmres), std::make_unique<ScaledDiagonalInverse>(1.0));
  qn.SetLinearSolver(ksp);
  qn.SetShiftInvert(std::complex<double>(0.0, 1.5));
  int nconv = 0;
  const auto out = test::CaptureStdout([&] { nconv = qn.Solve(); });
  REQUIRE(nconv == 1);
  CAPTURE(ls, l0, qn.GetEigenvalue(0));
  CHECK(std::abs(qn.GetEigenvalue(0) - ls) <= 1.0e-8 * std::abs(ls));
  CHECK(qn.GetError(0, EigenvalueSolver::ErrorType::BACKWARD) <= tol);
  if (Mpi::Root(Mpi::World()))
  {
    // The first attempt, from the seed, converges.
    std::smatch m;
    REQUIRE(std::regex_search(
        out, m, std::regex(R"(Eigenvalue 0, Quasi-Newton (converged in (\d+)|did not))")));
    CAPTURE(m[0].str());
    REQUIRE(m[2].matched);
    CHECK(std::stoi(m[2].str()) < 40);
  }
}
