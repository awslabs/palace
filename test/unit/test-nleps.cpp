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

namespace
{

// Quasi-Newton on a DiagonalLorentzProblem from the initial guesses (eigenvalue estimate,
// global dof of the eigenvector), with num modes requested above the target σ = i s.
// T(λ) is diagonal, so the iteration from the eigenvector of a dof stays on that dof: it
// converges to a root of that dof, whatever the eigenvalue estimate.
struct QuasiNewtonResult
{
  std::vector<std::complex<double>> eig;
  std::vector<double> backward_error;
  std::vector<double> weight;  // Weight of the largest entry in each eigenvector
  std::vector<int> dof;        // Global dof of the largest entry in each eigenvector
  std::string out;
};

using GuessVector = std::vector<std::pair<long long, double>>;  // (global dof, value)

QuasiNewtonResult
RunQuasiNewtonV(const DiagonalLorentzProblem &p,
                const std::vector<std::pair<std::complex<double>, GuessVector>> &guesses,
                int num, double s, double tol, int ksp_its = 100,
                std::complex<double> pc_omega = 0.0, int max_its = 100)
{
  std::vector<std::complex<double>> theta;
  std::vector<ComplexVector> vecs;
  for (const auto &[l, entries] : guesses)
  {
    ComplexVector x0(p.LocalSize());
    x0 = 0.0;
    for (int i = 0; i < p.LocalSize(); i++)
    {
      for (const auto &[j, val] : entries)
      {
        if (p.gid[i] == j)
        {
          x0.Real()(i) = val;
        }
      }
    }
    theta.push_back(l);
    vecs.push_back(std::move(x0));
  }
  const auto num_guesses = static_cast<int>(theta.size());
  QuasiNewtonSolver qn(Mpi::World(),
                       std::make_unique<SeedEigenpairs>(std::move(theta), std::move(vecs)),
                       num_guesses, 1, true);
  qn.SetTol(tol);
  qn.SetMaxIter(max_its);
  qn.SetOperators(*p.K, *p.M, EigenvalueSolver::ScaleType::NONE);
  qn.SetExtraSystemMatrix([&p](std::complex<double> l) { return p.A2(l); });
  qn.SetPreconditionerUpdate(
      [&p, pc_omega](std::complex<double> a0, std::complex<double>, std::complex<double> a2,
                     std::complex<double> a3)
      {
        if (pc_omega == 0.0)
        {
          return p.Preconditioner(a0, a2, a3);
        }
        const auto l = std::complex<double>(0.0, 1.0) * pc_omega;
        return p.Preconditioner(1.0, l * l, pc_omega);
      });
  qn.SetNumModes(num, num);
  qn.SetPreconditionerLag(10, 1.0e-4);
  qn.SetMaxRestart(2);
  auto gmres = std::make_unique<GmresSolver<ComplexOperator>>(Mpi::World(), 0);
  gmres->SetRelTol(1.0e-12);
  gmres->SetMaxIter(ksp_its);
  ComplexKspSolver ksp(std::move(gmres), std::make_unique<ScaledDiagonalInverse>(1.0));
  qn.SetLinearSolver(ksp);
  qn.SetShiftInvert(std::complex<double>(0.0, s));
  QuasiNewtonResult r;
  int nconv = 0;
  r.out = test::CaptureStdout([&] { nconv = qn.Solve(); });
  ComplexVector x(p.LocalSize());
  for (int m = 0; m < nconv; m++)
  {
    r.eig.push_back(qn.GetEigenvalue(m));
    r.backward_error.push_back(qn.GetError(m, EigenvalueSolver::ErrorType::BACKWARD));
    // The dof of the eigenvector: its largest entry, and the fraction of ‖x‖² there.
    qn.GetEigenvector(m, x);
    double norm2 = 0.0, max2 = 0.0;
    int dof = -1;
    for (int i = 0; i < p.LocalSize(); i++)
    {
      const double a2 = std::norm(std::complex<double>(x.Real()(i), x.Imag()(i)));
      norm2 += a2;
      if (a2 > max2)
      {
        max2 = a2;
        dof = static_cast<int>(p.gid[i]);
      }
    }
    Mpi::GlobalSum(1, &norm2, Mpi::World());
    Mpi::GlobalMaxLoc(1, &max2, &dof, Mpi::World());
    r.weight.push_back(max2 / norm2);
    r.dof.push_back(dof);
  }
  return r;
}

// The diagonal entries of a DiagonalLorentzProblem whose dofs 1 to n have the entries kc,
// after a decoupled dof 0, all on the first process. The c vector of the Quasi-Newton
// eigenvalue update holds the same entries on every process (from srand(111)), and on macOS
// the first of them, 8.8e-4, is small: with an eigenvector on the first dof of a process,
// the iteration then depends on the platform.
std::vector<double> CoupledDofs(const std::vector<double> &kc)
{
  std::vector<double> k = {1.0e4};
  k.insert(k.end(), kc.begin(), kc.end());
  AddDecoupledDofs(k, 1.0e4, static_cast<int>(k.size()));
  return k;
}

// Check that no process has its first dof among dofs 1 to n.
void RequireNotFirstDofs(const DiagonalLorentzProblem &p, long long n)
{
  int first_coupled = (p.LocalSize() > 0 && p.gid[0] >= 1 && p.gid[0] <= n) ? 1 : 0;
  Mpi::GlobalMax(1, &first_coupled, Mpi::World());
  REQUIRE(first_coupled == 0);
}

QuasiNewtonResult
RunQuasiNewton(const DiagonalLorentzProblem &p,
               const std::vector<std::pair<std::complex<double>, long long>> &guesses,
               int num, double s, double tol)
{
  std::vector<std::pair<std::complex<double>, GuessVector>> v;
  for (const auto &[l, j] : guesses)
  {
    v.push_back({l, {{j, 1.0}}});
  }
  return RunQuasiNewtonV(p, v, num, s, tol);
}

// The root of dof k with the larger (upper) or smaller (lower) imaginary part.
std::complex<double> UpperRoot(const DiagonalLorentzProblem &p, double k)
{
  const auto roots = p.Roots(k);
  return *std::ranges::max_element(roots,
                                   [](auto a, auto b) { return a.imag() < b.imag(); });
}

std::complex<double> LowerRoot(const DiagonalLorentzProblem &p, double k)
{
  const auto roots = p.Roots(k);
  return *std::ranges::min_element(roots,
                                   [](auto a, auto b) { return a.imag() < b.imag(); });
}

}  // namespace

// The iteration from an initial guess can converge to an eigenvalue other than the one the
// guess approximates, farther from the target. Here the guess nearest the target 1.5i (at
// 2.16i, with the eigenvector of dof 2) converges to the upper root of dof 2 at 2.567i,
// and the next one (at 2.20i, with the eigenvector of dof 1) to the upper root of dof 1
// at 2.153i, nearer the target. The guesses nearer the target than the eigenvalues found
// are refined too, and the nearest eigenvalue is returned: stopping at the first converged
// eigenvalue would return the root of dof 2.
TEST_CASE("Quasi-Newton returns the eigenvalue nearest the target of those found",
          "[nleps][Serial][Parallel]")
{
  constexpr double tol = 1.0e-10;
  const DiagonalLorentzProblem p(CoupledDofs({4.0, 6.0, 9.0}));
  RequireNotFirstDofs(p, 3);
  const auto l0 = UpperRoot(p, 4.0), l1 = UpperRoot(p, 6.0);
  REQUIRE(l0.imag() < 2.16);
  REQUIRE(l1.imag() > 2.20);
  const auto r = RunQuasiNewton(
      p, {{std::complex<double>(0.0, 2.16), 2}, {std::complex<double>(0.0, 2.20), 1}}, 1,
      1.5, tol);
  CAPTURE(r.out);
  REQUIRE(r.eig.size() == 1);
  CAPTURE(r.eig[0], l0, l1);
  CHECK(std::abs(r.eig[0] - l0) <= 1.0e-8 * std::abs(l0));
  CHECK(r.backward_error[0] <= tol);
  CHECK(r.weight[0] >= 1.0 - 1.0e-12);
  if (Mpi::Root(Mpi::World()))
  {
    // The first eigenvalue found is the root of dof 2.
    std::smatch m;
    const bool converged = std::regex_search(
        r.out, m,
        std::regex(R"(Eigenvalue 0, Quasi-Newton converged in \d+ iterations )"
                   R"(\(([+-]?[0-9.]+e[+-][0-9]+)([+-][0-9.]+e[+-][0-9]+)i\))"));
    CHECK(converged);
    CHECK((converged && std::abs(std::stod(m[2].str()) - l1.imag()) < 1.0e-3));
  }
}

// The initial guesses left once the requested eigenvalues are found are refined if they are
// nearer the target than the eigenvalue found, in frequency and in the complex plane. Here
// the guess with the eigenvector of dof 2, the upper root of dof 2 at 2.567i, is the most
// accurate, and is refined first. The other guess, with the eigenvector of dof 1, has the
// frequency of the upper root of dof 1 at 2.153i, nearer the target 1.5i, and a damping
// (real part) larger by d. With d = 0 it is nearer the target in the complex plane too: it
// is refined, and the root of dof 1 is returned. With d = 1 it is farther, as the guess of
// a strongly damped mode, and is not refined.
TEST_CASE("Quasi-Newton refines the guesses nearer the target in the complex plane",
          "[nleps][Serial][Parallel]")
{
  constexpr double tol = 1.0e-10;
  const DiagonalLorentzProblem p(CoupledDofs({4.0, 6.0, 9.0}));
  RequireNotFirstDofs(p, 3);
  const auto l0 = UpperRoot(p, 4.0), l1 = UpperRoot(p, 6.0);
  const double d = GENERATE(0.0, 1.0);
  const std::complex<double> s(0.0, 1.5), guess0(l0.real() - d, l0.imag() + 1.0e-3);
  REQUIRE((std::abs(guess0 - s) < std::abs(l1 - s)) == (d == 0.0));
  const auto r = RunQuasiNewton(p, {{l1, 2}, {guess0, 1}}, 1, s.imag(), tol);
  CAPTURE(r.out, d);
  REQUIRE(r.eig.size() == 1);
  const auto expected = (d == 0.0) ? l0 : l1;
  CAPTURE(r.eig[0], expected);
  CHECK(std::abs(r.eig[0] - expected) <= 1.0e-8 * std::abs(expected));
  CHECK(r.backward_error[0] <= tol);
  CHECK(r.dof[0] == ((d == 0.0) ? 1 : 2));
  if (Mpi::Root(Mpi::World()))
  {
    const std::regex converged(R"(Quasi-Newton converged)");
    const auto num_converged =
        std::distance(std::sregex_iterator(r.out.begin(), r.out.end(), converged),
                      std::sregex_iterator());
    CHECK(num_converged == ((d == 0.0) ? 2 : 1));
  }
}

// A degenerate eigenvalue at the last position of those returned is returned with all its
// eigenpairs found. Here dofs 2 and 3 have the same roots, and their guesses, the most
// accurate, give the two requested eigenpairs. The guess of dof 1 is nearer the target, and
// its root is found next: the two eigenpairs of the root of dofs 2 and 3 are returned with
// it, three eigenpairs for the two requested. The eigenvectors of the degenerate eigenvalue
// are any two independent combinations of those of dofs 2 and 3.
TEST_CASE("Quasi-Newton returns all the eigenpairs found of the last eigenvalue",
          "[nleps][Serial][Parallel]")
{
  constexpr double tol = 1.0e-10;
  const DiagonalLorentzProblem p(CoupledDofs({4.0, 6.0, 6.0, 9.0}));
  RequireNotFirstDofs(p, 4);
  const auto l0 = UpperRoot(p, 4.0), l1 = UpperRoot(p, 6.0);
  const std::complex<double> offset(1.0e-6, 1.0e-6);
  const auto r = RunQuasiNewton(
      p, {{l1 + offset, 2}, {l1 + offset, 3}, {l0 + 1.0e3 * offset, 1}}, 2, 1.5, tol);
  CAPTURE(r.out);
  REQUIRE(r.eig.size() == 3);
  const std::vector<std::complex<double>> expected = {l0, l1, l1};
  for (std::size_t m = 0; m < expected.size(); m++)
  {
    CAPTURE(m, r.eig[m], expected[m]);
    CHECK(std::abs(r.eig[m] - expected[m]) <= 1.0e-8 * std::abs(expected[m]));
    CHECK(r.backward_error[m] <= tol);
  }
  CHECK(r.weight[0] >= 1.0 - 1.0e-12);
  CHECK(r.dof[0] == 1);
  CHECK(std::min(r.dof[1], r.dof[2]) == 2);
  CHECK(std::max(r.dof[1], r.dof[2]) == 3);
}

// The eigenvectors of distinct eigenvalues of a nonlinear eigenvalue problem can be
// linearly dependent: here every dof has two roots above the target 0.5i (one on each side
// of the Lorentz pole at |λ| = 1) with the same eigenvector. The lower roots of dofs 1 and
// 2 are found and deflated first; the eigenvector of each upper root is then in the span
// of the deflated eigenvectors. Such an eigenpair cannot be deflated with minimality index
// 1 (its part orthogonal to the deflated eigenvectors is rounding error), but it is an
// eigenpair of T(λ): all four are returned, converged, each on its dof. Deflating them
// corrupted the deflated problem: the fourth eigenpair was returned with a backward error
// above the tolerance.
TEST_CASE("Quasi-Newton finds eigenpairs with linearly dependent eigenvectors",
          "[nleps][Serial][Parallel]")
{
  constexpr double tol = 1.0e-10;
  const DiagonalLorentzProblem p(CoupledDofs({4.0, 9.0}));
  RequireNotFirstDofs(p, 2);
  const std::vector<std::pair<std::complex<double>, long long>> roots = {
      {LowerRoot(p, 4.0), 1},
      {LowerRoot(p, 9.0), 2},
      {UpperRoot(p, 4.0), 1},
      {UpperRoot(p, 9.0), 2}};
  // The initial guesses are the roots, off by 1e-4, with the eigenvectors of their dofs.
  std::vector<std::pair<std::complex<double>, long long>> guesses;
  for (const auto &[l, j] : roots)
  {
    guesses.emplace_back(l + std::complex<double>(1.0e-4, 1.0e-4), j);
  }
  const auto r = RunQuasiNewton(p, guesses, 4, 0.5, tol);
  CAPTURE(r.out);
  REQUIRE(r.eig.size() == 4);
  for (std::size_t m = 0; m < roots.size(); m++)
  {
    CAPTURE(m, r.eig[m], roots[m].first);
    CHECK(std::abs(r.eig[m] - roots[m].first) <= 1.0e-8 * std::abs(roots[m].first));
    CHECK(r.backward_error[m] <= tol);
    CHECK(r.weight[m] >= 1.0 - 1.0e-12);
  }
}

// An eigenpair which is not deflated can be found again, from another guess: it is then
// recognized, and returned once. Here, as in the test above, the upper roots of dofs 1 and
// 2 are not deflated, and a second guess of the upper root of dof 1 converges to it again,
// in each of its attempts.
TEST_CASE("Quasi-Newton recognizes an eigenpair found before", "[nleps][Serial][Parallel]")
{
  constexpr double tol = 1.0e-10;
  const DiagonalLorentzProblem p(CoupledDofs({4.0, 9.0}));
  RequireNotFirstDofs(p, 2);
  const std::vector<std::pair<std::complex<double>, long long>> roots = {
      {LowerRoot(p, 4.0), 1},
      {LowerRoot(p, 9.0), 2},
      {UpperRoot(p, 4.0), 1},
      {UpperRoot(p, 9.0), 2}};
  std::vector<std::pair<std::complex<double>, long long>> guesses;
  for (const auto &[l, j] : roots)
  {
    guesses.emplace_back(l + std::complex<double>(1.0e-4, 1.0e-4), j);
  }
  guesses.insert(guesses.begin() + 3, guesses[2]);
  const auto r = RunQuasiNewton(p, guesses, 4, 0.5, tol);
  CAPTURE(r.out);
  REQUIRE(r.eig.size() == 4);
  for (std::size_t m = 0; m < roots.size(); m++)
  {
    CAPTURE(m, r.eig[m], roots[m].first);
    CHECK(std::abs(r.eig[m] - roots[m].first) <= 1.0e-8 * std::abs(roots[m].first));
    CHECK(r.backward_error[m] <= tol);
  }
  if (Mpi::Root(Mpi::World()))
  {
    CHECK(r.out.find(", found before") != std::string::npos);
  }
}

// An attempt which stalls is ended, and the guess is refined again without the stall exit
// when it can change the result. Here the linear solves of the iteration are inexact: one
// GMRES iteration, with the diagonal preconditioner of T(λ) at the fixed λ = 1.8i. From a
// guess of the upper root of dof 1 with a small part on the decoupled dof 0, the residual
// then stays above the tolerance 1e-12 from the 10th iteration on: the attempt stalls after
// 31 iterations (two windows of 10 without a reduction of 10%).
TEST_CASE("Quasi-Newton ends the attempts which stall", "[nleps][Serial][Parallel]")
{
  constexpr double tol = 1.0e-12;
  const DiagonalLorentzProblem p(CoupledDofs({4.0, 6.0, 9.0}));
  RequireNotFirstDofs(p, 3);
  const auto l0 = UpperRoot(p, 4.0), l1 = UpperRoot(p, 6.0);
  const std::pair<std::complex<double>, GuessVector> stalling = {
      l0 + std::complex<double>(0.0, 1.0e-4), {{1, 1.0}, {0, 1.0e-8}}};
  // The iterations of the first attempt, which stalls, and the output after it.
  auto Stall = [](const std::string &out) -> std::pair<int, std::string>
  {
    std::smatch m;
    if (!std::regex_search(
            out, m,
            std::regex(R"(Quasi-Newton stalled after (\d+) iterations, restarting\.)")))
    {
      return {0, ""};
    }
    return {std::stoi(m[1].str()), m.suffix().str()};
  };
  SECTION("Fewer eigenvalues found than requested")
  {
    // With no eigenvalue found, the only guess is refined again without the stall exit, to
    // the 100 iterations of an attempt, before the random guess.
    const auto r = RunQuasiNewtonV(p, {stalling}, 1, 1.5, tol, 1, 1.8);
    CAPTURE(r.out);
    if (Mpi::Root(Mpi::World()))
    {
      const auto [its, after] = Stall(r.out);
      CHECK(its == 31);
      const auto again = after.find("Refining again the 1 initial guess which stalled");
      CHECK(again != std::string::npos);
      CHECK((again != std::string::npos &&
             after.find("Quasi-Newton did not converge in 100 iterations", again) !=
                 std::string::npos));
    }
  }
  SECTION("Guess nearer the target than the eigenvalue found")
  {
    // The next guess, of the upper root of dof 2, farther from the target, converges: the
    // guess which stalled is then refined again, in an attempt without the stall exit.
    const auto r =
        RunQuasiNewtonV(p, {stalling, {l1 + std::complex<double>(0.0, 1.0e-4), {{2, 1.0}}}},
                        1, 1.5, tol, 1, 1.8);
    CAPTURE(r.out);
    REQUIRE(r.eig.size() == 1);
    CHECK((std::abs(r.eig[0] - l1) <= 1.0e-8 * std::abs(l1) ||
           std::abs(r.eig[0] - l0) <= 1.0e-8 * std::abs(l0)));
    if (Mpi::Root(Mpi::World()))
    {
      const auto [its, after] = Stall(r.out);
      CHECK(its == 31);
      const auto found = after.find("Eigenvalue 0, Quasi-Newton converged in");
      CHECK(found != std::string::npos);
      CHECK((found != std::string::npos &&
             after.find("0 NLEPS (nconv=1, restart=0)", found) != std::string::npos));
      CHECK(after.find("stalled") == std::string::npos);
    }
  }
}
