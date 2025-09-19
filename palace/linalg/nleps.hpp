// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_NLEPS_HPP
#define PALACE_LINALG_NLEPS_HPP

#include <complex>
#include <memory>
#include <mpi.h>
#include "linalg/eps.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/spaceoperator.hpp"

namespace palace
{

//
// Abstract base class for nonlinear eigensolvers.
// Currently only implemented for complex scalar interface.
//
class NonLinearEigenvalueSolver : public EigenvalueSolver
{
protected:
  // MPI communicator.
  MPI_Comm comm;

  // Control print level for debugging.
  int print;

  // Number of eigenvalues to compute and problem size.
  int nev, n;

  // Relative eigenvalue error convergence tolerance for the solver.
  double rtol;

  // Maximum number of iterations.
  int nleps_it;

  // Specifies which part of the spectrum to search for.
  EigenvalueSolver::WhichType which_type;

  // Variables for scaling, from Higham et al., IJNME 2008.
  double gamma, delta;

  // Parameters defining the spectral transformation.
  std::complex<double> sigma;
  bool sinvert;

  // Storage for computed eigenvalues and eigenvectors.
  std::vector<std::complex<double>> eigenvalues;
  std::vector<ComplexVector> eigenvectors;
  std::unique_ptr<int[]> perm;

  // Storage for eigenpairs initial guesses.
  std::vector<std::complex<double>> init_eigenvalues;
  std::vector<ComplexVector> init_eigenvectors;

  // Storage for computed residual norms and eigenvector scalings.
  std::unique_ptr<double[]> res, xscale;

  // Reference to linear solver used for operator action of P(σ)⁻¹ (not owned).
  ComplexKspSolver *opInv;

  // Reference to solver for projecting an intermediate vector onto a divergence-free space
  // (not owned).
  const DivFreeSolver<ComplexVector> *opProj;

  // Reference to SpaceOperator to compute A2 matrix.
  SpaceOperator *space_op;

  // Reference to matrix used for weighted inner products (not owned). May be nullptr, in
  // which case identity is used.
  const Operator *opB;

  // Workspace vector for operator applications.
  mutable ComplexVector x1, y1;

  // Helper routine for computing the eigenvector normalization.
  double GetEigenvectorNorm(const ComplexVector &x, ComplexVector &Bx) const;

  // Helper routine for computing the eigenpair residual.
  virtual double GetResidualNorm(std::complex<double> l, const ComplexVector &x,
                                 ComplexVector &r) const = 0;

  // Helper routine for computing the backward error.
  virtual double GetBackwardScaling(std::complex<double> l) const = 0;

  // Get the associated MPI communicator.
  virtual MPI_Comm GetComm() const = 0;

  // Return problem type name.
  virtual const char *GetName() const = 0;

public:
  NonLinearEigenvalueSolver(MPI_Comm comm, int print);

  // The linear solver will be configured to compute the action of T(σ)⁻¹
  // where σ is the current eigenvalue estimate.
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

  // Set initial guess for the eigenpairs.
  void SetInitialGuess(const std::vector<std::complex<double>> &init_eig,
                       const std::vector<ComplexVector> &init_V,
                       const std::vector<double> &init_errors);

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

// Quasi-Newton nonlinear eigenvalue solver for (K + λ C + λ² M + A2(λ)) x = 0.
class QuasiNewtonSolver : public NonLinearEigenvalueSolver
{
private:
  // References to matrices defining the nonlinear eigenvalue problem
  // (not owned).
  const ComplexOperator *opK, *opC, *opM;

  // Operators used in the iterative linear solver.
  std::unique_ptr<ComplexOperator> opA2, opA, opP;

  // Operator norms for scaling.
  mutable double normK, normC, normM;

  // Update frequency of the preconditioner during Newton iterations.
  int preconditioner_lag;

  // Update tolerance of the preconditioner (no update below tol).
  double preconditioner_tol;

  // Maximum number of Newton attempts with the same initial guess.
  int max_restart;

protected:
  double GetResidualNorm(std::complex<double> l, const ComplexVector &x,
                         ComplexVector &r) const override;

  double GetBackwardScaling(std::complex<double> l) const override;

  const char *GetName() const override { return "QuasiNewton"; }

public:
  QuasiNewtonSolver(MPI_Comm comm, int print);

  using NonLinearEigenvalueSolver::SetOperators;
  void SetOperators(SpaceOperator &space_op, const ComplexOperator &K,
                    const ComplexOperator &M, ScaleType type) override;
  void SetOperators(SpaceOperator &space_op, const ComplexOperator &K,
                    const ComplexOperator &C, const ComplexOperator &M,
                    ScaleType type) override;

  // Set the update frequency of the preconditioner.
  void SetPreconditionerLag(int preconditioner_update_freq,
                            double preconditioner_update_tol);

  // Set the maximum number of restarts with the same initial guess.
  void SetMaxRestart(int max_num_restart);

  // Solve the nonlinear eigenvalue problem.
  int Solve() override;

  MPI_Comm GetComm() const override { return comm; }
};

//
// Interpolation operators to approximate the nonlinear A2 operator.
//
class Interpolation
{
public:
  Interpolation() = default;
  virtual ~Interpolation() = default;
  virtual void Interpolate(int order, const std::complex<double> sigma_min,
                           const std::complex<double> sigma_max) = 0;
  virtual void Mult(int order, const ComplexVector &x, ComplexVector &y) const = 0;
  virtual void AddMult(int order, const ComplexVector &x, ComplexVector &y,
                       std::complex<double> a = 1.0) const = 0;
};

// Newton polynomial interpolation to approximate the nonlinear A2 operator.
class NewtonInterpolationOperator : public Interpolation
{
private:
  // Reference to SpaceOperator to compute A2 matrix.
  SpaceOperator *space_op;

  // Number of points used in the interpolation.
  int num_points;

  // Interpolation points.
  std::vector<std::complex<double>> points;

  // Monomial basis coefficients.
  std::vector<std::vector<std::complex<double>>> coeffs;

  // Divided difference operators.
  std::vector<std::vector<std::unique_ptr<ComplexOperator>>> ops;

  // Workspace objects for solver application.
  mutable ComplexVector rhs;

public:
  NewtonInterpolationOperator(SpaceOperator &space_op);

  // Interpolate the A2 matrix between sigma_min and sigma_max with a Newton polynomial.
  void Interpolate(int order, const std::complex<double> sigma_min,
                   const std::complex<double> sigma_max);

  // Perform multiplication with interpolation operator of specified order.
  void Mult(int order, const ComplexVector &x, ComplexVector &y) const;
  void AddMult(int order, const ComplexVector &x, ComplexVector &y,
               std::complex<double> a = 1.0) const;
};

}  // namespace palace

#endif  // PALACE_LINALG_NLEPS_HPP
