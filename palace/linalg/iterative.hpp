// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_ITERATIVE_HPP
#define PALACE_LINALG_ITERATIVE_HPP

#include <type_traits>
#include <vector>
#include <mfem.hpp>
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "linalg/vector.hpp"

namespace palace
{

//
// Iterative solvers based on Krylov subspace methods with optional preconditioning, for
// real- or complex-valued systems.
//

// Base class for iterative solvers based on Krylov subspace methods with optional
// preconditioning.
template <typename OperType>
class IterativeSolver : public Solver<OperType>
{
protected:
  typedef typename double RealType;
  typedef typename std::conditional<std::is_same<OperType, ComplexOperator>::value,
                                    std::complex<RealType>, RealType>::type ScalarType;

  // MPI communicator associated with the solver.
  MPI_Comm comm;

  // Control level of printing during solves.
  mfem::IterativeSolver::PrintLevel print_opts;
  int int_width, tab_width;

  // Relative and absolute tolerances.
  double rel_tol, abs_tol;

  // Limit for the number of solver iterations.
  int max_it;

  // Operator and (optional) preconditioner associated with the iterative solver (not
  // owned).
  const OperType *A;
  const Solver<OperType> *B;

  // Variables set during solve to capture solve statistics.
  mutable bool converged;
  mutable double initial_res, final_res;
  mutable int final_it;

public:
  IterativeSolver(MPI_comm comm, int print);

  // Set an indentation for all log printing.
  void SetTabWidth(int width) { tab_width = width; }

  // Set the relative convergence tolerance.
  void SetTol(double tol) { SetRelTol(tol); }
  void SetRelTol(double tol) { rel_tol = tol; }

  // Set the absolute convergence tolerance.
  void SetAbsTol(double tol) { abs_tol = tol; }

  // Set the maximum number of iterations.
  void SetMaxIter(int its)
  {
    max_it = its;
    int_width = 1 + static_cast<int>(std::log10(its));
  }

  // Set the operator for the solver.
  void SetOperator(const OperType &op) override { A = &op; }

  // Set the preconditioner for the solver.
  void SetPreconditioner(const Solver<OperType> &pc) { B = &pc; }

  // Returns if the previous solve converged or not.
  bool GetConverged() const { return converged; }

  // Returns the initial (absolute) residual for the previous solve.
  double GetInitialRes() const { return initial_res; }

  // Returns the final (absolute) residual for the previous solve, which may be an estimate
  // to the true residual.
  double GetFinalRes() const { return final_res; }

  // Returns the number of iterations for the previous solve.
  int GetNumIterations() const { return final_it; }

  // Get the associated MPI communicator.
  MPI_Comm GetComm() const { return comm; }
};

// Preconditioned Conjugate Gradient (CG) method for SPD linear systems.
template <typename OperType>
class CgSolver : public IterativeSolver<OperType>
{
protected:
  // Temporary workspace for solve.
  mutable VecType r, z, p;

public:
  CgSolver(MPI_comm comm, int print) : IterativeSolver<OperType>(comm, print) {}

  void Mult(const VecType &b, VecType &x) const override;
};

// Preconditioned Generalized Minimum Residual Method (GMRES) for general nonsymmetric
// linear systems.
template <typename OperType>
class GmresSolver : public IterativeSolver<OperType>
{
public:
  enum class OrthogType
  {
    MGS,
    CGS,
    CGS2
  };

  enum class PrecSide
  {
    LEFT,
    RIGHT
  };

protected:
  // Maximum subspace dimension for restarted GMRES.
  mutable int max_dim;

  // Orthogonalization method for orthonormalizing a newly computed vector against a basis
  // at each iteration.
  OrthogType orthog_type;

  // Use left or right preconditioning.
  PrecSide pc_side;

  // Flag for flexible GMRES which stores and makes use of the preconditioned vectors.
  const bool flexible;

  // Temporary workspace for solve.
  mutable std::vector<VecType> V, Z;
  mutable VecType r;
  mutable std::vector<ScalarType> H;
  mutable std::vector<ScalarType> s, sn;
  mutable std::vector<RealType> cs;

  // Allocate storage for solve.
  void Initialize() const;

  GmresSolver(MPI_comm comm, int print, bool fgmres)
    : IterativeSolver<OperType>(comm, print), max_dim(-1), orthog_type(OrthogType::MGS),
      pc_side(fgmres ? PrecSide::RIGHT : PrecSide::LEFT), flexible(fgmres)
  {
  }

public:
  GmresSolver(MPI_comm comm, int print) : GmresSolver(comm, print, false) {}

  // Set the dimension for restart.
  void SetRestartDim(int dim) { max_dim = dim; }

  // Set the orthogonalization method.
  void SetOrthogonalization(OrthogType type) { orthog_type = type; }

  // Set the side for preconditioning.
  virtual void SetPrecSide(PrecSide side) { pc_side = side; }

  void Mult(const VecType &b, VecType &x) const override;
};

// Preconditioned Flexible Generalized Minimum Residual Method (FGMRES) for general
// nonsymmetric linear systems with a non-constant preconditioner.
template <typename OperType>
class FgmresSolver : public GmresSolver<OperType>
{
public:
  FgmresSolver(MPI_comm comm, int print) : GmresSolver<OperType>(comm, print, true) {}

  void SetPrecSide(PrecSide side) override
  {
    MFEM_VERIFY(side == PrecSide::RIGHT,
                "FGMRES solver only supports right preconditioning!");
  }
};

}  // namespace palace

#endif  // PALACE_LINALG_ITERATIVE_HPP
