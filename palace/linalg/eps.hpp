// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_EPS_HPP
#define PALACE_LINALG_EPS_HPP

#include <complex>
#include <vector>
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

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
// Rational terms (Σⱼ cⱼ λ² / (λ - pⱼ)) B of a quadratic eigenvalue problem, for the
// eigensolvers which linearize the problem: a term is exact with one auxiliary block
// λ² / (λ - pⱼ) x for each pole (Su and Bai, Solving rational eigenvalue problems via
// linearization, SIAM J. Matrix Anal. Appl. 2011), so that a problem whose only
// nonlinearity is rational is a linear eigenvalue problem. The auxiliary unknown of a pole
// exists only on the true dofs where its B has a nonzero row: a full-size block would give
// the linearization the eigenvalue pⱼ once for each true dof outside the support of B,
// with a zero x, which is not an eigenvalue of the problem.
//
class RationalTerms
{
private:
  // The matrices B (not owned) with the masks of their supports (1 on the true dofs where
  // B has a nonzero row, 0 elsewhere), and for each pole its term and its coefficient.
  std::vector<const ComplexOperator *> ops;
  std::vector<Vector> masks;
  std::vector<std::size_t> terms;
  std::vector<std::complex<double>> poles, coeffs;

  // Workspace vector.
  mutable ComplexVector z;

public:
  // Add a term. Its poles are nonzero (a zero pole is a term c λ B) and distinct.
  void Add(const ComplexOperator &B, const std::vector<std::complex<double>> &p,
           const std::vector<std::complex<double>> &c)
  {
    MFEM_VERIFY(!p.empty() && p.size() == c.size(),
                "A rational term of an eigenvalue problem needs one coefficient per pole!");
    for (std::size_t j = 0; j < p.size(); j++)
    {
      MFEM_VERIFY(p[j] != 0.0, "A rational term of an eigenvalue problem has a zero pole!");
      for (std::size_t i = 0; i < j; i++)
      {
        // Equal poles would give two auxiliary blocks for one unknown.
        MFEM_VERIFY(p[i] != p[j],
                    "A rational term of an eigenvalue problem has a repeated pole!");
      }
      terms.push_back(ops.size());
    }
    ops.push_back(&B);
    poles.insert(poles.end(), p.begin(), p.end());
    coeffs.insert(coeffs.end(), c.begin(), c.end());

    // The support of B from its diagonal: a mass-type matrix has a nonzero diagonal entry
    // on every true dof it acts on, and essential dofs are eliminated with a zero diagonal.
    ComplexVector diag(B.Height());
    B.AssembleDiagonal(diag);
    auto &mask = masks.emplace_back(B.Height());
    mask.UseDevice(true);
    const auto *dr = diag.Real().HostRead();
    const auto *di = diag.Imag().HostRead();
    auto *m = mask.HostWrite();
    for (int i = 0; i < B.Height(); i++)
    {
      m[i] = (dr[i] != 0.0 || di[i] != 0.0) ? 1.0 : 0.0;
    }
  }

  std::size_t NumPoles() const { return poles.size(); }
  std::complex<double> Pole(std::size_t k) const { return poles[k]; }

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
  }

  // y += Σₜ Bₜ (Σⱼ f(cₜⱼ, pₜⱼ) aₜⱼ) for the auxiliary blocks a.
  template <typename Function>
  void AddMult(const std::vector<ComplexVector> &a, Function &&f, ComplexVector &y) const
  {
    std::size_t k = 0;
    for (std::size_t t = 0; t < ops.size(); t++)
    {
      z.SetSize(ops[t]->Height());
      z.UseDevice(true);
      z = 0.0;
      for (; k < poles.size() && terms[k] == t; k++)
      {
        z.AXPY(f(coeffs[k], poles[k]), a[k]);
      }
      ops[t]->AddMult(z, y, std::complex<double>(1.0, 0.0));
    }
  }

  // y += Σₜ (Σⱼ f(cₜⱼ, pₜⱼ)) Bₜ x.
  template <typename Function>
  void AddMult(const ComplexVector &x, Function &&f, ComplexVector &y) const
  {
    std::size_t k = 0;
    for (std::size_t t = 0; t < ops.size(); t++)
    {
      std::complex<double> s = 0.0;
      for (; k < poles.size() && terms[k] == t; k++)
      {
        s += f(coeffs[k], poles[k]);
      }
      ops[t]->AddMult(x, y, s);
    }
  }

  // Restrict the auxiliary blocks to the supports of the B of their terms.
  void Mask(std::vector<ComplexVector> &a) const
  {
    for (std::size_t k = 0; k < poles.size(); k++)
    {
      a[k].Real() *= masks[terms[k]];
      a[k].Imag() *= masks[terms[k]];
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

  // Add the rational term (Σⱼ cⱼ λ² / (λ - pⱼ)) B to the quadratic eigenvalue problem, for
  // the solvers which linearize it (see RationalTerms). Must be called on every process
  // before SetOperators. The linear solver (SetLinearSolver) must then solve with T(σ),
  // which includes every rational term at the shift σ.
  virtual void AddPoleTerms(const ComplexOperator &B,
                            const std::vector<std::complex<double>> &poles,
                            const std::vector<std::complex<double>> &coeffs)
  {
    MFEM_ABORT("AddPoleTerms not defined!");
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
