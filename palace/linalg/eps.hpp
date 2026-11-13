// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_EPS_HPP
#define PALACE_LINALG_EPS_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <vector>
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

namespace palace
{

template <typename VecType>
class DivFreeSolver;

namespace nleps
{

//
// Known poles p_k of A2(λ), as in dispersive materials. The SLP nonlinear eigensolver
// iterates on d(λ) T(λ), d(λ) = Π_k (λ - p_k): the eigenpairs are those of T, but the local
// linear model has no pole, so a step does not cross a nearby pole into the basin of the
// root on its far side. Since d is scalar, it cancels from the correction T⁻¹ (d T) and
// only the derivative changes: (d T)′ / d = T′ + L(λ) T with L = d′/d.
//
struct PoleTerms
{
  std::vector<std::complex<double>> poles;

  bool Empty() const { return poles.empty(); }

  // L(λ) = d′(λ) / d(λ).
  std::complex<double> LogDerivative(std::complex<double> l) const
  {
    std::complex<double> L = 0.0;
    for (const auto p : poles)
    {
      L += 1.0 / (l - p);
    }
    return L;
  }
};

}  // namespace nleps

//
// Rational terms f(λ) B of a quadratic eigenvalue problem, for the eigensolvers which
// linearize the problem, with a fixed matrix B and a scalar rational function
//   f(λ) = d₀ + d₁ λ + d₂ λ² + Σⱼ cⱼ / (λ - pⱼ) .
// The polynomial part adds to K, C, and M, and each pole pⱼ adds one auxiliary block
//   wⱼ = κ B x / ((λ - pⱼ) β) ,   (λ - pⱼ) wⱼ = κ B x / β ,
// with β the largest diagonal entry of B and κ = σ² / γ, for the shift σ and the eigenvalue
// scale γ of the solver, so that cⱼ B x / (λ - pⱼ) = cⱼ β wⱼ / κ (Su and Bai, Solving
// rational eigenvalue problems via linearization, SIAM J. Matrix Anal. Appl. 2011). With
// them, a problem whose only nonlinearity is rational is a linear eigenvalue problem.
//
// The scale κ does not change the eigenvalues, but it weighs the blocks against x in the
// Krylov basis, and so sets how the error of an inexact linear solve enters the modes. With
// κ = σ² / γ, a block has the size of λ² B x / (γ (λ - pⱼ) β) near the shift, as the
// auxiliary unknowns λ² x / (γ (λ - pⱼ)) of the earlier material pole terms c λ² / (λ - p).
// For a Lorentz cavity with the shift 1.6 % from a pole and the default linear solver
// tolerance (regression case cylinder_cavity_lorentz_pairs, 8 modes), the largest
// backward error is then 5.9e-8, against 3.2e-8 with that earlier form. With κ = σ - pⱼ,
// the right-hand side of the shifted solve has terms of the size 1 / (σ - pⱼ)², and the
// largest backward error is 5.4e-7, with modes outside the tolerance of the case; with
// κ = γ, it is 2.3e-7 (and a mode of the Lorentz-filled box of the unit tests is above its
// bound 1e-8).
//
// A block holds B x, not x: from a start with zero auxiliary blocks, the blocks of the
// Krylov vectors stay in the range of B. A block of x restricted to the support of B
// would give the linearization the eigenvalue pⱼ for every vector of the kernel of B on
// its support, with a zero x: the surface gradients for the boundary curl-curl matrix of
// a second-order absorbing boundary, for example.
//
class RationalTerms
{
private:
  // The matrices B (not owned) with their polynomial coefficients {d₀, d₁, d₂}, scales β,
  // and the range [begin, end) of their poles.
  struct Term
  {
    const ComplexOperator *B;
    std::array<std::complex<double>, 3> poly;
    double scale;
    std::size_t begin, end;
  };
  std::vector<Term> terms;
  std::vector<std::size_t> term_of_pole;
  std::vector<std::complex<double>> poles, residues;

  // Workspace vector.
  mutable ComplexVector z;

public:
  // Add a term f(λ) B, f(λ) = poly[0] + poly[1] λ + poly[2] λ² + Σⱼ c[j] / (λ - p[j]),
  // with distinct poles (on every process).
  void Add(MPI_Comm comm, const ComplexOperator &B,
           const std::array<std::complex<double>, 3> &poly,
           const std::vector<std::complex<double>> &p,
           const std::vector<std::complex<double>> &c)
  {
    MFEM_VERIFY(p.size() == c.size(),
                "A rational term of an eigenvalue problem needs one residue per pole!");
    for (std::size_t j = 0; j < p.size(); j++)
    {
      for (std::size_t i = 0; i < j; i++)
      {
        // Equal poles would give two auxiliary blocks for one unknown.
        MFEM_VERIFY(p[i] != p[j],
                    "A rational term of an eigenvalue problem has a repeated pole!");
      }
      term_of_pole.push_back(terms.size());
    }

    // The scale β of the blocks: the largest diagonal entry of B.
    ComplexVector diag(B.Height());
    diag.UseDevice(true);
    B.AssembleDiagonal(diag);
    double scale = 0.0;
    const auto *dr = diag.Real().HostRead();
    const auto *di = diag.Imag().HostRead();
    for (int i = 0; i < B.Height(); i++)
    {
      scale = std::max(scale, std::abs(std::complex<double>(dr[i], di[i])));
    }
    Mpi::GlobalMax(1, &scale, comm);
    terms.push_back(
        {&B, poly, (scale > 0.0) ? scale : 1.0, poles.size(), poles.size() + p.size()});
    poles.insert(poles.end(), p.begin(), p.end());
    residues.insert(residues.end(), c.begin(), c.end());
  }

  bool Empty() const { return terms.empty(); }
  std::size_t NumTerms() const { return terms.size(); }
  std::size_t NumPoles() const { return poles.size(); }
  std::complex<double> Pole(std::size_t k) const { return poles[k]; }
  std::complex<double> Residue(std::size_t k) const { return residues[k]; }

  // The polynomial coefficient dᵢ of term t.
  std::complex<double> Poly(std::size_t t, int i) const { return terms[t].poly[i]; }

  // The pole part Σⱼ cⱼ / (λ - pⱼ) of term t, and its value f(λ).
  std::complex<double> PoleValue(std::size_t t, std::complex<double> l) const
  {
    std::complex<double> f = 0.0;
    for (std::size_t k = terms[t].begin; k < terms[t].end; k++)
    {
      f += residues[k] / (l - poles[k]);
    }
    return f;
  }
  std::complex<double> Value(std::size_t t, std::complex<double> l) const
  {
    const auto &d = terms[t].poly;
    return d[0] + l * (d[1] + l * d[2]) + PoleValue(t, l);
  }

  // The linearization requires a shift-and-invert spectral transformation, with a shift
  // which is not a pole.
  void CheckShift(bool sinvert, std::complex<double> sigma) const
  {
    for (const auto p : poles)
    {
      MFEM_VERIFY(sinvert && sigma != p,
                  "An eigenvalue problem with rational terms requires shift-and-invert, "
                  "with a shift which is not a pole!");
    }
    MFEM_VERIFY(sinvert || terms.empty(),
                "An eigenvalue problem with rational terms requires shift-and-invert!");
  }

  // y += Σₜ a(t) Bₜ x, for a scalar a(t) of each term t.
  template <typename Coefficient>
  void AddMult(const ComplexVector &x, Coefficient &&a, ComplexVector &y) const
  {
    for (std::size_t t = 0; t < terms.size(); t++)
    {
      const std::complex<double> s = a(t);
      if (s != 0.0)
      {
        terms[t].B->AddMult(x, y, s);
      }
    }
  }

  // y += Σₖ b(k) β wₖ, for the auxiliary blocks w and a scalar b(k) of each pole k.
  template <typename Coefficient>
  void AddAuxiliary(const std::vector<ComplexVector> &w, Coefficient &&b,
                    ComplexVector &y) const
  {
    for (std::size_t k = 0; k < poles.size(); k++)
    {
      y.AXPY(b(k) * terms[term_of_pole[k]].scale, w[k]);
    }
  }

  // yₖ = a(k) B x / β + b(k) wₖ for every pole k, with one product with B for each term.
  template <typename CoefficientA, typename CoefficientB>
  void SetAuxiliary(const ComplexVector &x, CoefficientA &&a,
                    const std::vector<ComplexVector> &w, CoefficientB &&b,
                    std::vector<ComplexVector> &y) const
  {
    for (const auto &term : terms)
    {
      z.SetSize(x.Size());
      z.UseDevice(true);
      term.B->Mult(x, z);
      for (std::size_t k = term.begin; k < term.end; k++)
      {
        y[k].AXPBYPCZ(a(k) / term.scale, z, b(k), w[k], 0.0);
      }
    }
  }
};

//
// Pure abstract base class for solving generalized linear eigenvalue problems problems or
// quadratic polynomial eigenvalue problems.
//
class EigenvalueSolver
{
public:
  enum class ScaleType
  {
    NONE,
    NORM_2
  };

  enum class WhichType
  {
    LARGEST_MAGNITUDE,
    SMALLEST_MAGNITUDE,
    LARGEST_REAL,
    SMALLEST_REAL,
    LARGEST_IMAGINARY,
    SMALLEST_IMAGINARY,
    TARGET_MAGNITUDE,
    TARGET_REAL,
    TARGET_IMAGINARY
  };

  enum class ErrorType
  {
    ABSOLUTE,
    RELATIVE,
    BACKWARD
  };

public:
  EigenvalueSolver() = default;
  virtual ~EigenvalueSolver() = default;

  // Set operators for the generalized eigenvalue problem, quadratic polynomial
  // eigenvalue problem, or nonlinear eigenvalue problem.
  virtual void SetOperators(const ComplexOperator &K, const ComplexOperator &M,
                            ScaleType type)
  {
    MFEM_ABORT("SetOperators not defined!");
  }

  virtual void SetOperators(const ComplexOperator &K, const ComplexOperator &C,
                            const ComplexOperator &M, ScaleType type)
  {
    MFEM_ABORT("SetOperators not defined!");
  }

  virtual void SetOperators(const ComplexOperator &K, const ComplexOperator &M,
                            std::function<const ComplexOperator &(std::complex<double>)> A2,
                            ScaleType type)
  {
    MFEM_ABORT("SetOperators not defined!");
  }

  // Add the rational term f(λ) B to the quadratic eigenvalue problem, f(λ) = poly[0] +
  // poly[1] λ + poly[2] λ² + Σⱼ residues[j] / (λ - poles[j]) with distinct poles, for the
  // solvers which linearize it (see RationalTerms). Must be called on every process before
  // SetOperators. The linear solver (SetLinearSolver) must then solve with T(σ), which
  // includes every rational term at the shift σ.
  virtual void AddRationalTerm(const ComplexOperator &B,
                               const std::array<std::complex<double>, 3> &poly,
                               const std::vector<std::complex<double>> &poles,
                               const std::vector<std::complex<double>> &residues)
  {
    MFEM_ABORT("AddRationalTerm not defined!");
  }

  // Set the frequency-dependent A2 matrix function A2(λ) for the nonlinear eigensolvers.
  virtual void SetExtraSystemMatrix(
      std::function<std::unique_ptr<ComplexOperator>(std::complex<double>)>)
  {
    MFEM_ABORT("SetExtraSystemMatrix not defined!");
  }

  // Set the known poles of A2(λ) for the SLP nonlinear eigensolver (see nleps::PoleTerms).
  virtual void SetExtraSystemPoles(const std::vector<std::complex<double>> &)
  {
    MFEM_ABORT("SetExtraSystemPoles not defined!");
  }

  virtual void SetPreconditionerUpdate(std::function<std::unique_ptr<ComplexOperator>(
                                           std::complex<double>, std::complex<double>,
                                           std::complex<double>, std::complex<double>)>)
  {
    MFEM_ABORT("SetPreconditionerUpdate not defined!");
  }

  // For the linear generalized case, the linear solver should be configured to compute the
  // action of M⁻¹ (with no spectral transformation) or (K - σ M)⁻¹. For the quadratic
  // case, the linear solver should be configured to compute the action of M⁻¹ (with no
  // spectral transformation) or P(σ)⁻¹.
  virtual void SetLinearSolver(ComplexKspSolver &ksp) = 0;

  // Set the projection operator for enforcing the divergence-free constraint.
  virtual void SetDivFreeProjector(const DivFreeSolver<ComplexVector> &divfree) = 0;

  // Set optional B matrix used for weighted inner products. This must be set explicitly
  // even for generalized problems, otherwise the identity will be used.
  virtual void SetBMat(const Operator &B) = 0;

  // Get scaling factors used by the solver.
  virtual double GetScalingGamma() const = 0;
  virtual double GetScalingDelta() const = 0;

  // Set the number of required eigenmodes.
  virtual void SetNumModes(int num_eig, int num_vec = 0) = 0;

  // Set solver tolerance.
  virtual void SetTol(double tol) = 0;

  // Set maximum number of Arnoldi update iterations.
  virtual void SetMaxIter(int max_it) = 0;

  // Set target spectrum for the eigensolver. When a spectral transformation is used, this
  // applies to the spectrum of the shifted operator.
  virtual void SetWhichEigenpairs(WhichType type) = 0;

  // Set shift-and-invert spectral transformation.
  virtual void SetShiftInvert(std::complex<double> s, bool precond = false) = 0;

  // Set an initial vector for the solution subspace.
  virtual void SetInitialSpace(const ComplexVector &v) = 0;

  // Solve the eigenvalue problem. Returns the number of converged eigenvalues.
  virtual int Solve() = 0;

  // Get the corresponding eigenvalue.
  virtual std::complex<double> GetEigenvalue(int i) const = 0;

  // Get the corresponding eigenvector. Eigenvectors are normalized such that ||x||₂ = 1,
  // unless the B-matrix is set for weighted inner products.
  virtual void GetEigenvector(int i, ComplexVector &x) const = 0;

  // Get the corresponding eigenpair error.
  virtual double GetError(int i, ErrorType type) const = 0;

  // Re-normalize the given number of eigenvectors, for example if the matrix B for weighted
  // inner products has changed. This does not perform re-orthogonalization with respect to
  // the new matrix, only normalization.
  virtual void RescaleEigenvectors(int num_eig) = 0;
};

}  // namespace palace

#endif  // PALACE_LINALG_EPS_HPP
