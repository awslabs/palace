// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <complex>
#include <iterator>
#include <memory>
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
#include "linalg/rap.hpp"
#if defined(PALACE_WITH_SLEPC)
#include "linalg/slepc.hpp"
#endif
#include "linalg/vector.hpp"
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

}  // namespace

// The eigensolvers which linearize the quadratic problem solve a problem with rational
// terms as a linear one: the eigenvalues nearest the shift are the exact roots, on both
// sides of the poles (|p±| = 1), with one of a dof without the medium among them.
TEST_CASE("Linearized eigensolvers solve rational terms exactly",
          "[nleps][Serial][Parallel]")
{
  const RationalDiagonalProblem q(
      {1.5, 2.5, 4.0, 6.5, 3.0, 12.0, 16.0, 20.0, 25.0, 30.0, 36.0, 42.0},
      {1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0});
  const std::complex<double> sigma(0.0, 1.6);
  constexpr int nev = 9;
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
      vr[i] = 1.0 + 0.1 * q.p.gid[i];
      vi[i] = 0.3 - 0.05 * q.p.gid[i];
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
    eigen.SetNumModes(nev, 24);
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
  const DiagonalLorentzProblem p(k);
  const double sigma = 10.3;
  std::vector<std::complex<double>> d(p.LocalSize());
  ComplexVector start(p.LocalSize());
  for (int i = 0; i < p.LocalSize(); i++)
  {
    d[i] = p.k.at(p.gid[i]) - sigma;
    start.Real()(i) = 1.0 + 0.1 * p.gid[i];
    start.Imag()(i) = 0.3 - 0.05 * p.gid[i];
  }
  const DiagonalComplexOperator A(std::move(d));
  auto gmres = std::make_unique<GmresSolver<ComplexOperator>>(Mpi::World(), 0);
  gmres->SetRelTol(1.0e-14);
  gmres->SetMaxIter(100);
  ComplexKspSolver ksp(std::move(gmres), std::make_unique<ScaledDiagonalInverse>(1.0));
  ksp.SetOperators(A, A);
  constexpr int nev = 8;
  auto Solve = [&](int max_it)
  {
    arpack::ArpackEPSSolver eigen(Mpi::World(), 0);
    eigen.SetOperators(*p.K, *p.M, EigenvalueSolver::ScaleType::NONE);
    eigen.SetLinearSolver(ksp);
    eigen.SetShiftInvert(sigma);
    eigen.SetWhichEigenpairs(EigenvalueSolver::WhichType::LARGEST_MAGNITUDE);
    eigen.SetNumModes(nev, nev + 2);
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
