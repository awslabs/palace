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
#include "utils/labels.hpp"

namespace palace
{

//
// Iterative solvers based on Krylov subspace methods with optional preconditioning, for
// real- or complex-valued systems.
//

// Base class for iterative solvers based on Krylov subspace methods with optional
// preconditioning.
template <OperatorType OperType>
class IterativeSolver : public Solver<OperType>
{
protected:
  using RealType = double;
  using ScalarType = std::conditional_t<std::is_same_v<OperType, ComplexOperator>,
                                        std::complex<RealType>, RealType>;

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
  mutable int final_it, final_restarts;

  // Set when the true residual contradicts a converged residual estimate.
  mutable bool false_converged = false;

  // Enable timer contribution for Timer::PRECONDITIONER.
  bool use_timer;

  // Performance tally ids: preconditioned operator application and orthogonalization.
  int perf_apply = -1, perf_orth = -1;

public:
  IterativeSolver(MPI_Comm comm, int print);

  // Set an indentation for all log printing.
  void SetTabWidth(int width) { tab_width = width; }

  // Get and set the control level of printing during solves.
  const mfem::IterativeSolver::PrintLevel &GetPrintLevel() const { return print_opts; }
  void SetPrintLevel(const mfem::IterativeSolver::PrintLevel &opts) { print_opts = opts; }

  // Set the relative convergence tolerance.
  void SetTol(double tol) { SetRelTol(tol); }
  void SetRelTol(double tol) { rel_tol = tol; }

  // Get the relative convergence tolerance.
  double GetRelTol() const { return rel_tol; }

  // Get the absolute convergence tolerance.
  double GetAbsTol() const { return abs_tol; }

  // Set the absolute convergence tolerance.
  void SetAbsTol(double tol) { abs_tol = tol; }

  // Get the maximum number of iterations.
  int GetMaxIter() const { return max_it; }

  // Set the maximum number of iterations.
  void SetMaxIter(int its)
  {
    max_it = its;
    int_width = 1 + static_cast<int>(std::log10(its));
  }

  // Set the operator for the solver.
  void SetOperator(const OperType &op) override
  {
    A = &op;
    this->height = op.Height();
    this->width = op.Width();
  }

  // Set the preconditioner for the solver.
  void SetPreconditioner(const Solver<OperType> &pc) { B = &pc; }

  void SetPerfTallies(const std::string &prefix, Timer::Index parent) override;

  // Returns if the previous solve converged or not.
  bool GetConverged() const { return converged && (rel_tol > 0.0 || abs_tol > 0.0); }

  // Returns if the previous solve converged on a residual estimate which the true residual
  // contradicts (the solution is invalid and GetConverged() is false).
  bool GetFalseConverged() const { return false_converged; }

  // Returns the initial (absolute) residual for the previous solve.
  double GetInitialRes() const { return initial_res; }

  // Returns the final (absolute) residual for the previous solve, which may be an estimate
  // to the true residual.
  double GetFinalRes() const { return final_res; }

  // Returns the number of iterations for the previous solve.
  int GetNumIterations() const { return final_it; }

  // Bytes held by the Krylov basis (zero for short-recurrence methods).
  virtual std::size_t GetBasisBytes() const { return 0; }
  // Returns the number of restarts for the previous solve (GMRES and FGMRES).
  int GetNumRestarts() const { return final_restarts; }

  // Get the associated MPI communicator.
  MPI_Comm GetComm() const { return comm; }

  // Activate preconditioner timing during solves.
  void EnableTimer() { use_timer = true; }
};

// Preconditioned Conjugate Gradient (CG) method for SPD linear systems.
template <OperatorType OperType>
class CgSolver : public IterativeSolver<OperType>
{
protected:
  using VecType = typename Solver<OperType>::VecType;
  using RealType = typename IterativeSolver<OperType>::RealType;
  using ScalarType = typename IterativeSolver<OperType>::ScalarType;

  using IterativeSolver<OperType>::comm;
  using IterativeSolver<OperType>::print_opts;
  using IterativeSolver<OperType>::int_width;
  using IterativeSolver<OperType>::tab_width;

  using IterativeSolver<OperType>::rel_tol;
  using IterativeSolver<OperType>::abs_tol;
  using IterativeSolver<OperType>::max_it;

  using IterativeSolver<OperType>::A;
  using IterativeSolver<OperType>::B;

  using IterativeSolver<OperType>::converged;
  using IterativeSolver<OperType>::initial_res;
  using IterativeSolver<OperType>::final_res;
  using IterativeSolver<OperType>::final_it;
  using IterativeSolver<OperType>::final_restarts;

  // Temporary workspace for solve.
  mutable VecType r, z, p;

  // Inverse diagonal of the operator for the Jacobi stopping norm, if any.
  const Vector *dinv_norm = nullptr;

public:
  CgSolver(MPI_Comm comm, int print) : IterativeSolver<OperType>(comm, print) {}

  // Stop on the residual norm of Jacobi preconditioning, √(D⁻¹ r, r) for the given inverse
  // diagonal D⁻¹ of the operator (real, on the host), instead of √(B r, r) for the
  // preconditioner B (nullptr): the accuracy at a given tolerance is then about the same as
  // with Jacobi preconditioning. It is computed in the same pass and reduction as (B r, r).
  void SetJacobiStoppingNorm(const Vector *dinv) { dinv_norm = dinv; }

  void Mult(const VecType &b, VecType &x) const override;
};

// Preconditioned Generalized Minimum Residual Method (GMRES) for general nonsymmetric
// linear systems.
template <OperatorType OperType>
class GmresSolver : public IterativeSolver<OperType>
{
protected:
  using VecType = typename Solver<OperType>::VecType;
  using RealType = typename IterativeSolver<OperType>::RealType;
  using ScalarType = typename IterativeSolver<OperType>::ScalarType;

  using IterativeSolver<OperType>::comm;
  using IterativeSolver<OperType>::print_opts;
  using IterativeSolver<OperType>::int_width;
  using IterativeSolver<OperType>::tab_width;

  using IterativeSolver<OperType>::rel_tol;
  using IterativeSolver<OperType>::abs_tol;
  using IterativeSolver<OperType>::max_it;

  using IterativeSolver<OperType>::A;
  using IterativeSolver<OperType>::B;

  using IterativeSolver<OperType>::converged;
  using IterativeSolver<OperType>::initial_res;
  using IterativeSolver<OperType>::final_res;
  using IterativeSolver<OperType>::final_it;
  using IterativeSolver<OperType>::final_restarts;

  // Maximum subspace dimension for restarted GMRES.
  mutable int max_dim;

  // Orthogonalization method for orthonormalizing a newly computed vector against a basis
  // at each iteration.
  Orthogonalization gs_orthog;

  // Use left or right preconditioning.
  PreconditionerSide pc_side;

  // Temporary workspace for solve.
  mutable std::vector<VecType> V;
  mutable VecType r;
  mutable std::vector<ScalarType> H;
  mutable std::vector<ScalarType> s, sn;
  mutable std::vector<RealType> cs;

  // The Hessenberg matrix before the Givens rotations, and the workspace of the CGS2
  // Arnoldi process, which reorthogonalizes each basis vector one iteration late (see
  // DelayedCGS2Step): the buffer of the single global reduction per iteration (the
  // reorthogonalization coefficients of the previous vector come first), the norm of the
  // previous vector after its reorthogonalization, and the entry of s before the latest
  // rotation.
  mutable std::vector<ScalarType> Hraw, dots;
  mutable RealType nu;
  mutable ScalarType s_prev;

  // Dimension of the Arnoldi relation of the last cycle of the most recent solve (see
  // GetArnoldiRelation).
  mutable int arnoldi_dim = 0;

  // Stagnation test: stop when the residual norm decreased by less than the relative
  // fraction stag_tol over the last stag_it iterations (stag_it = 0 disables the test).
  int stag_it;
  double stag_tol;
  mutable std::vector<RealType> res_hist;
  mutable bool stagnated;

  // Bound for the true residual of converged left-preconditioned solves, see
  // SetTrueResidualCheck.
  double true_res_tol;

  // Bound for the true residual of the initial guess, see SetInitialGuessCheck.
  double guess_tol = 0.0;

  // Record the residual norm of the current iteration and return whether the solve has
  // stagnated.
  bool Stagnated(RealType res) const
  {
    res_hist.push_back(res);
    const auto n = res_hist.size();
    stagnated = (stag_it > 0 && n > static_cast<std::size_t>(stag_it) &&
                 res > (1.0 - stag_tol) * res_hist[n - 1 - stag_it]);
    return stagnated;
  }

  // Allocate storage for solve.
  virtual void Initialize() const;
  virtual void Update(int j) const;

  // Orthogonalize V[j + 1] = M V[j] against V[0], ..., V[j], normalize it, and apply the
  // Givens rotations to the new column j of the Hessenberg matrix and to s. With CGS2, V[j]
  // and V[j + 1] are only finalized at the next iteration and column j - 1 is completed and
  // rotated again first (see DelayedCGS2Step).
  void Arnoldi(int j) const;

  // Arnoldi step of classical Gram-Schmidt with delayed reorthogonalization (DCGS2): one
  // global reduction per iteration and two passes over the basis instead of the 2 (j + 1)
  // passes and j + 2 reductions of MGS. See Swirydowicz et al., Numer. Linear Algebra Appl.
  // 28 (2021), e2343, and Bielich et al., Parallel Comput. 112 (2022), 102940. For j > 0,
  // on entry V[j] = u is the previous vector after one Gram-Schmidt pass, before its
  // reorthogonalization and normalization, and V[j + 1] = y = M u. The step
  // reorthogonalizes and normalizes V[j], completing column j - 1 of Hraw, computes column
  // j from the projections of y (using M Q = Q H for the reorthogonalization of u), with
  // the subdiagonal entry from the Pythagorean theorem, and leaves V[j + 1] after one
  // Gram-Schmidt pass. The convergence test uses this column, which only misses the
  // reorthogonalization of the next iteration.
  void DelayedCGS2Step(int j) const;

  // Generate and apply the Givens rotation of column i of H (after the previous rotations).
  void RotateColumn(int i) const;

public:
  GmresSolver(MPI_Comm comm, int print)
    : IterativeSolver<OperType>(comm, print), max_dim(-1),
      gs_orthog(Orthogonalization::MGS), pc_side(PreconditionerSide::LEFT), nu(1.0),
      s_prev(0.0), stag_it(0), stag_tol(0.0), stagnated(false), true_res_tol(0.0)
  {
  }

  // Set the dimension for restart.
  void SetRestartDim(int dim) { max_dim = dim; }

  // Stop when the residual norm decreased by less than the relative fraction tol over the
  // last its iterations (its = 0 disables the test).
  void SetStagnation(int its, double tol)
  {
    stag_it = its;
    stag_tol = tol;
  }

  // Whether the last solve stopped on the stagnation test.
  bool GetStagnated() const { return stagnated; }

  // With left preconditioning, the stopping test is on the preconditioned residual, which
  // an unstable preconditioner can reduce without reducing the true residual (false
  // convergence). A solve converged on the relative tolerance then requires norm(b - Ax) <=
  // max(tol, 1e7 rel_tol) norm(b), at the cost of one operator application (tol <= 0, the
  // default, disables the check). The bound allows for nearly singular (shift-and-invert)
  // systems, whose true residual can exceed the preconditioned one by orders of magnitude.
  // The preconditioned residual also weighs the components of the residual very unevenly
  // when the preconditioner scales them by many orders of magnitude (stiff rows of
  // multiscale meshes, solutions amplified by high-Q resonances), so an accurate solution
  // can fail the check, in particular from an initial guess with residual components which
  // the preconditioned residual does not see. A solve from a nonzero initial guess which
  // fails the check is solved again from zero. If the check still fails, the iterate is
  // continued with right preconditioning, which minimizes the true residual, until the true
  // residual is half the bound (for up to the maximum number of iterations more): if the
  // continuation satisfies the check within a relative distance of 1e-4 of the iterate, the
  // iterate is kept; otherwise (an unstable preconditioner moves it by about its norm) the
  // solve is falsely converged and returns the iterate.
  void SetTrueResidualCheck(double tol) { true_res_tol = tol; }

  // With an initial guess x, start from zero instead if norm(b - A x) >= tol norm(b): a
  // guess whose true residual is no smaller than that of zero can still have a small
  // preconditioned residual, while its residual components which the preconditioner does
  // not see (for example the stiff rows of multiscale meshes) remain in the solution of a
  // left-preconditioned solve and fail the true residual check. Costs two inner products;
  // the solves from a guess which is kept are unchanged. tol <= 0 (the default) disables
  // it.
  void SetInitialGuessCheck(double tol) { guess_tol = tol; }

  // Set the orthogonalization method.
  void SetOrthogonalization(Orthogonalization orthog) { gs_orthog = orthog; }

  // Set the side for preconditioning.
  virtual void SetPreconditionerSide(PreconditionerSide side) { pc_side = side; }
  PreconditionerSide GetPreconditionerSide() const { return pc_side; }

  // The Arnoldi relation B A V[0:m] = V[0:m+1] H (left preconditioning, or A B V[0:m] =
  // V[0:m+1] H with right preconditioning) of the last cycle of the most recent solve:
  // returns m, the number of completed Arnoldi steps (0 if none, or if the last cycle was
  // the right-preconditioned continuation of a left-preconditioned solve, see
  // SetTrueResidualCheck), and sets basis to the orthonormal basis vectors (the first m + 1
  // are used) and hbar to the (m + 1) x m upper Hessenberg matrix H, column-major with
  // leading dimension ld. Valid until the next solve.
  int GetArnoldiRelation(const std::vector<VecType> *&basis, const ScalarType *&hbar,
                         int &ld) const
  {
    basis = &V;
    hbar = Hraw.data();
    ld = max_dim + 1;
    return arnoldi_dim;
  }

  std::size_t GetBasisBytes() const override
  {
    std::size_t n = 0;
    for (const auto &v : V)
    {
      n += v.Size();
    }
    return n * sizeof(ScalarType);
  }

  void Mult(const VecType &b, VecType &x) const override;
};

// Preconditioned Flexible Generalized Minimum Residual Method (FGMRES) for general
// nonsymmetric linear systems with a non-constant preconditioner.
template <OperatorType OperType>
class FgmresSolver : public GmresSolver<OperType>
{
protected:
  using VecType = typename GmresSolver<OperType>::VecType;
  using RealType = typename GmresSolver<OperType>::RealType;
  using ScalarType = typename GmresSolver<OperType>::ScalarType;

  using GmresSolver<OperType>::comm;
  using GmresSolver<OperType>::print_opts;
  using GmresSolver<OperType>::int_width;
  using GmresSolver<OperType>::tab_width;

  using GmresSolver<OperType>::rel_tol;
  using GmresSolver<OperType>::abs_tol;
  using GmresSolver<OperType>::max_it;

  using GmresSolver<OperType>::A;
  using GmresSolver<OperType>::B;

  using GmresSolver<OperType>::converged;
  using GmresSolver<OperType>::initial_res;
  using GmresSolver<OperType>::final_res;
  using GmresSolver<OperType>::final_it;
  using GmresSolver<OperType>::final_restarts;

  using GmresSolver<OperType>::max_dim;
  using GmresSolver<OperType>::gs_orthog;
  using GmresSolver<OperType>::pc_side;
  using GmresSolver<OperType>::V;
  using GmresSolver<OperType>::H;
  using GmresSolver<OperType>::s;
  using GmresSolver<OperType>::sn;
  using GmresSolver<OperType>::cs;
  using GmresSolver<OperType>::dots;
  using GmresSolver<OperType>::nu;

  // Temporary workspace for solve.
  mutable std::vector<VecType> Z;

  // With CGS2, Z[j] is the preconditioned vector before the delayed reorthogonalization of
  // V[j], and the preconditioned basis is Z T with T upper triangular.
  mutable std::vector<ScalarType> T;

  // Allocate storage for solve.
  void Initialize() const override;
  void Update(int j) const override;

public:
  FgmresSolver(MPI_Comm comm, int print) : GmresSolver<OperType>(comm, print)
  {
    pc_side = PreconditionerSide::RIGHT;
  }

  void SetPreconditionerSide(const PreconditionerSide side) override
  {
    MFEM_VERIFY(side == PreconditionerSide::RIGHT,
                "FGMRES solver only supports right preconditioning!");
  }

  std::size_t GetBasisBytes() const override
  {
    std::size_t n = 0;
    for (const auto &z : Z)
    {
      n += z.Size();
    }
    return GmresSolver<OperType>::GetBasisBytes() + n * sizeof(ScalarType);
  }

  void Mult(const VecType &b, VecType &x) const override;
};

}  // namespace palace

#endif  // PALACE_LINALG_ITERATIVE_HPP
