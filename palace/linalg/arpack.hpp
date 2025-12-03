// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_ARPACK_HPP
#define PALACE_LINALG_ARPACK_HPP

#if defined(PALACE_WITH_ARPACK)

#include <complex>
#include <memory>
#include <mpi.h>
#include "linalg/eps.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

namespace arpack
{

//
// A wrapper for the ARPACK/PARPACK library for generalized linear eigenvalue problems or
// quadratic polynomial eigenvalue problems. Shift-and-invert spectral transformations are
// used to compute interior eigenvalues. Currently only implemented for complex scalar
// interface.
//
class ArpackEigenvalueSolver : public EigenvalueSolver
{
protected:
  // MPI communicator for PARPACK.
  MPI_Comm comm;

  // Control print level for debugging.
  int print;

  // Status variable for ARPACK.
  int info;

  // Number eigenvalues to be computed, subspace dimension, and problem size.
  int nev, ncv, n;

  // Relative eigenvalue error convergence tolerance for the solver.
  double rtol;

  // Maximum number of Arnoldi update iterations.
  int arpack_it;

  // Specifies which part of the spectrum to search for.
  EigenvalueSolver::WhichType which_type;

  // Variables for scaling, from Higham et al., IJNME 2008.
  double gamma, delta;

  // Parameters defining the spectral transformation.
  std::complex<double> sigma;
  bool sinvert;

  // Storage for computed eigenvalues.
  std::unique_ptr<std::complex<double>[]> eig;
  std::unique_ptr<int[]> perm;

  // Storage for Arnoldi basis vectors.
  std::unique_ptr<std::complex<double>[]> V;

  // Storage for computed residual norms and eigenvector scalings.
  std::unique_ptr<double[]> res, xscale;

  // On input used to define optional initial guess, on output stores final residual
  // vector.
  std::unique_ptr<std::complex<double>[]> r;

  // Reference to linear solver used for operator action for M⁻¹ (with no spectral
  // transformation) or (K - σ M)⁻¹ (generalized EVP with shift-and- invert) or P(σ)⁻¹
  // (polynomial with shift-and-invert) (not owned).
  const ComplexKspSolver *opInv;

  // Reference to solver for projecting an intermediate vector onto a divergence-free space
  // (not owned).
  const DivFreeSolver<ComplexVector> *opProj;

  // Reference to matrix used for weighted inner products (not owned). May be nullptr, in
  // which case identity is used.
  const Operator *opB;

  // Workspace vector for operator applications.
  mutable ComplexVector x1, y1, z1;

  // Perform the ARPACK RCI loop.
  int SolveInternal(int n, std::complex<double> *r, std::complex<double> *V,
                    std::complex<double> *eig, int *perm);

  // Helper routine for parameter checking.
  void CheckParameters() const;

  // Helper routines for ARPACK RCI.
  virtual void ApplyOp(const std::complex<double> *px, std::complex<double> *py) const = 0;
  virtual void ApplyOpB(const std::complex<double> *px, std::complex<double> *py) const = 0;

  // Helper routine for computing the eigenvector normalization.
  double GetEigenvectorNorm(const ComplexVector &x, ComplexVector &Bx) const;

  // Helper routine for computing the eigenpair residual.
  virtual double GetResidualNorm(std::complex<double> l, const ComplexVector &x,
                                 ComplexVector &r) const = 0;

  // Helper routine for computing the backward error.
  virtual double GetBackwardScaling(std::complex<double> l) const = 0;

  // Return problem type name.
  virtual const char *GetName() const = 0;

public:
  ArpackEigenvalueSolver(MPI_Comm comm, int print);

  //  For the linear generalized case, the linear solver should be configured to compute the
  //  action of M⁻¹ (with no spectral transformation) or (K - σ M)⁻¹. For the quadratic
  //  case, the linear solver should be configured to compute the action of M⁻¹ (with no
  //  spectral transformation) or P(σ)⁻¹.
  void SetLinearSolver(ComplexKspSolver &ksp) override;

  // Set the projection operator for enforcing the divergence-free constraint.
  void SetDivFreeProjector(const DivFreeSolver<ComplexVector> &divfree) override;

  // Set optional B matrix used for weighted inner products. This must be set explicitly
  // even for generalized problems, otherwise the identity will be used.
  void SetBMat(const Operator &B) override;

  // Get scaling factors used by the solver.
  double GetScalingGamma() const override { return gamma; }
  double GetScalingDelta() const override { return delta; }

  // Set the number of required eigenmodes.
  void SetNumModes(int num_eig, int num_vec = 0) override;

  // Set solver tolerance.
  void SetTol(double tol) override;

  // Set maximum number of Arnoldi update iterations.
  void SetMaxIter(int max_it) override;

  // Set target spectrum for the eigensolver. When a spectral transformation is used, this
  // applies to the spectrum of the shifted operator.
  void SetWhichEigenpairs(WhichType type) override;

  // Set shift-and-invert spectral transformation.
  void SetShiftInvert(std::complex<double> s, bool precond = false) override;

  // Set an initial vector for the solution subspace.
  void SetInitialSpace(const ComplexVector &v) override;

  // Solve the eigenvalue problem. Returns the number of converged eigenvalues.
  int Solve() override = 0;

  // Get the corresponding eigenvalue.
  std::complex<double> GetEigenvalue(int i) const override;

  // Get the corresponding eigenvector. Eigenvectors are normalized such that ||x||₂ = 1,
  // unless the B-matrix is set for weighted inner products.
  void GetEigenvector(int i, ComplexVector &x) const override;

  // Get the corresponding eigenpair error.
  double GetError(int i, ErrorType type) const override;

  // Re-normalize the given number of eigenvectors, for example if the matrix B for weighted
  // inner products has changed. This does not perform re-orthogonalization with respect to
  // the new matrix, only normalization.
  void RescaleEigenvectors(int num_eig) override;
};

// Generalized eigenvalue problem solver: K x = λ M x .
class ArpackEPSSolver : public ArpackEigenvalueSolver
{
private:
  // References to matrices defining the generalized eigenvalue problem (not owned).
  const ComplexOperator *opK, *opM;

  // Operator norms for scaling.
  mutable double normK, normM;

protected:
  void ApplyOp(const std::complex<double> *px, std::complex<double> *py) const override;
  void ApplyOpB(const std::complex<double> *px, std::complex<double> *py) const override;

  double GetResidualNorm(std::complex<double> l, const ComplexVector &x,
                         ComplexVector &r) const override;

  double GetBackwardScaling(std::complex<double> l) const override;

  const char *GetName() const override { return "EPS"; }

public:
  ArpackEPSSolver(MPI_Comm comm, int print);

  using ArpackEigenvalueSolver::SetOperators;
  void SetOperators(const ComplexOperator &K, const ComplexOperator &M,
                    ScaleType type) override;

  int Solve() override;
};

// Quadratic eigenvalue problem solver: P(λ) x = (K + λ C + λ² M) x = 0 .
class ArpackPEPSolver : public ArpackEigenvalueSolver
{
private:
  // References to matrices defining the quadratic polynomial eigenvalue problem
  // (not owned).
  const ComplexOperator *opK, *opC, *opM;

  // Operator norms for scaling.
  mutable double normK, normC, normM;

  // Workspace vectors for operator applications.
  mutable ComplexVector x2, y2;

protected:
  void ApplyOp(const std::complex<double> *px, std::complex<double> *py) const override;
  void ApplyOpB(const std::complex<double> *px, std::complex<double> *py) const override;

  double GetResidualNorm(std::complex<double> l, const ComplexVector &x,
                         ComplexVector &r) const override;

  double GetBackwardScaling(std::complex<double> l) const override;

  const char *GetName() const override { return "PEP"; }

public:
  ArpackPEPSolver(MPI_Comm comm, int print);

  using ArpackEigenvalueSolver::SetOperators;
  void SetOperators(const ComplexOperator &K, const ComplexOperator &C,
                    const ComplexOperator &M, ScaleType type) override;

  int Solve() override;
};

}  // namespace arpack

}  // namespace palace

#endif

#endif  // PALACE_LINALG_ARPACK_HPP
