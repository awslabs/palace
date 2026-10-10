// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "ksp.hpp"

#include <mfem.hpp>
#include "fem/fespace.hpp"
#include "linalg/amg.hpp"
#include "linalg/ams.hpp"
#include "linalg/cudss.hpp"
#include "linalg/gmg.hpp"
#include "linalg/jacobi.hpp"
#include "linalg/mumps.hpp"
#include "linalg/strumpack.hpp"
#include "linalg/superlu.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace
{

template <OperatorType OperType>
std::unique_ptr<IterativeSolver<OperType>>
ConfigureKrylovSolver(const config::LinearSolverData &linear, int verbose, MPI_Comm comm)
{
  // Create the solver.
  std::unique_ptr<IterativeSolver<OperType>> ksp;
  const auto type = linear.krylov_solver;
  const int print = verbose;
  switch (type)
  {
    case KrylovSolver::CG:
      ksp = std::make_unique<CgSolver<OperType>>(comm, print);
      break;
    case KrylovSolver::GMRES:
      {
        auto gmres = std::make_unique<GmresSolver<OperType>>(comm, print);
        gmres->SetRestartDim(linear.max_size);
        ksp = std::move(gmres);
      }
      break;
    case KrylovSolver::FGMRES:
      {
        auto fgmres = std::make_unique<FgmresSolver<OperType>>(comm, print);
        fgmres->SetRestartDim(linear.max_size);
        ksp = std::move(fgmres);
      }
      break;
    case KrylovSolver::MINRES:
    case KrylovSolver::BICGSTAB:
    case KrylovSolver::DEFAULT:
      MFEM_ABORT("Unexpected solver type for Krylov solver configuration!");
      break;
  }
  ksp->SetInitialGuess(linear.initial_guess);
  ksp->SetRelTol(linear.tol);
  ksp->SetAbsTol(linear.abs_tol);
  ksp->SetMaxIter(linear.max_it);

  // Configure preconditioning side (only for GMRES).
  if (linear.pc_side != PreconditionerSide::DEFAULT && type != KrylovSolver::GMRES)
  {
    Mpi::Warning(comm,
                 "Preconditioner side will be ignored for non-GMRES iterative solvers!\n");
  }
  else
  {
    if (type == KrylovSolver::GMRES || type == KrylovSolver::FGMRES)
    {
      auto *gmres = static_cast<GmresSolver<OperType> *>(ksp.get());
      switch (linear.pc_side)
      {
        case PreconditionerSide::LEFT:
          gmres->SetPreconditionerSide(PreconditionerSide::LEFT);
          break;
        case PreconditionerSide::RIGHT:
          gmres->SetPreconditionerSide(PreconditionerSide::RIGHT);
          break;
        case PreconditionerSide::DEFAULT:
          // Do nothing. Set in ctors.
          break;
      }
    }
  }

  // Configure orthogonalization method for GMRES/FMGRES.
  if (type == KrylovSolver::GMRES || type == KrylovSolver::FGMRES)
  {
    // Because FGMRES inherits from GMRES, this is OK.
    auto *gmres = static_cast<GmresSolver<OperType> *>(ksp.get());
    gmres->SetOrthogonalization(linear.gs_orthog);
  }

  // Configure timing for the primary linear solver.
  ksp->EnableTimer();

  return ksp;
}

template <typename OperType, typename T, typename... U>
auto MakeWrapperSolver(const config::LinearSolverData &linear, U &&...args)
{
  // Sparse direct solver types copy the input matrix, so there is no need to save the
  // parallel assembled operator.
  constexpr bool save_assembled = !(false ||
#if defined(MFEM_USE_SUPERLU)
                                    std::is_same_v<T, SuperLUSolver> ||
#endif
#if defined(MFEM_USE_STRUMPACK)
                                    std::is_same_v<T, StrumpackSolver> ||
                                    std::is_same_v<T, StrumpackMixedPrecisionSolver> ||
#endif
#if defined(MFEM_USE_MUMPS)
                                    std::is_same_v<T, MumpsSolver> ||
#endif
#if defined(MFEM_USE_CUDSS)
                                    std::is_same_v<T, CuDSSSolver> ||
#endif
                                    false);
  return std::make_unique<MfemWrapperSolver<OperType>>(
      std::make_unique<T>(std::forward<U>(args)...), save_assembled,
      linear.complex_coarse_solve, linear.drop_small_entries, linear.reorder_reuse);
}

// Sparse direct solver for the PML subdomain correction of the geometric multigrid
// preconditioner: the solver specified by the "Type" option if it is a sparse direct
// solver, otherwise the default sparse direct solver. The PML subdomain matrix is always
// factored in complex form, and is symmetric (but indefinite) unless the system matrix is
// not. The factorization is never compressed or in mixed precision, and the accuracy of
// SuperLU_DIST factorizations is checked (its static pivoting can be unstable for this
// indefinite matrix, in which case tiny pivot replacement and iterative refinement are
// used): approximate solves can lead the (left-preconditioned) Krylov solver to converge to
// an inaccurate solution.
std::unique_ptr<Solver<ComplexOperator>>
ConfigurePMLSubdomainSolver(const config::LinearSolverData &linear,
                            MatrixSymmetry pc_mat_sym, int verbose, MPI_Comm comm)
{
  auto type =
      (linear.type == LinearSolver::STRUMPACK_MP) ? LinearSolver::STRUMPACK : linear.type;
  if (type != LinearSolver::SUPERLU && type != LinearSolver::STRUMPACK &&
      type != LinearSolver::MUMPS && type != LinearSolver::CUDSS)
  {
#if defined(MFEM_USE_SUPERLU)
    type = LinearSolver::SUPERLU;
#elif defined(MFEM_USE_STRUMPACK)
    type = LinearSolver::STRUMPACK;
#elif defined(MFEM_USE_MUMPS)
    type = LinearSolver::MUMPS;
#else
    MFEM_ABORT("PML subdomain correction requires building with SuperLU_DIST, STRUMPACK, "
               "or MUMPS!");
#endif
  }
  auto linear_pml = linear;
  linear_pml.complex_coarse_solve = true;
  linear_pml.strumpack_compression_type = SparseCompression::NONE;
  const auto sym =
      (pc_mat_sym == MatrixSymmetry::UNSYMMETRIC) ? pc_mat_sym : MatrixSymmetry::SYMMETRIC;
  const int print = verbose - 1;
  std::unique_ptr<Solver<ComplexOperator>> pc;
  switch (type)
  {
    case LinearSolver::SUPERLU:
#if defined(MFEM_USE_SUPERLU)
      pc = MakeWrapperSolver<ComplexOperator, SuperLUSolver>(
          linear_pml, comm, linear.sym_factorization, linear.superlu_3d,
          linear.reorder_reuse, print, true);
#endif
      break;
    case LinearSolver::STRUMPACK:
#if defined(MFEM_USE_STRUMPACK)
      pc = MakeWrapperSolver<ComplexOperator, StrumpackSolver>(
          linear_pml, comm, linear.sym_factorization, linear_pml.strumpack_compression_type,
          linear.strumpack_lr_tol, linear.strumpack_butterfly_l,
          linear.strumpack_lossy_precision, linear.reorder_reuse, print);
#endif
      break;
    case LinearSolver::MUMPS:
#if defined(MFEM_USE_MUMPS)
      pc = MakeWrapperSolver<ComplexOperator, MumpsSolver>(linear_pml, comm, sym,
                                                           linear.sym_factorization, 0.0,
                                                           linear.reorder_reuse, print);
#endif
      break;
    case LinearSolver::CUDSS:
#if defined(MFEM_USE_CUDSS)
      pc = MakeWrapperSolver<ComplexOperator, CuDSSSolver>(
          linear_pml, comm, sym, linear.sym_factorization, linear.reorder_reuse, print);
#endif
      break;
    default:
      break;
  }
  MFEM_VERIFY(pc, "Unable to construct the sparse direct solver for the PML subdomain "
                  "correction!");
  return pc;
}

template <OperatorType OperType>
std::unique_ptr<Solver<OperType>>
ConfigurePreconditionerSolver(const config::LinearSolverData &linear,
                              MatrixSymmetry pc_mat_sym, int verbose, MPI_Comm comm,
                              FiniteElementSpaceHierarchy &fespaces,
                              FiniteElementSpaceHierarchy *aux_fespaces)
{
  // Create the real-valued solver first.
  std::unique_ptr<Solver<OperType>> pc;
  const auto type = linear.type;
  const int print = verbose - 1;
  const bool coarse_solver = fespaces.GetNumLevels() > 1;
  switch (type)
  {
    case LinearSolver::AMS:
      // Can either be the coarse solve for geometric multigrid or the solver at the finest
      // space (in which case fespaces.GetNumLevels() == 1).
      MFEM_VERIFY(aux_fespaces, "AMS solver relies on both primary space "
                                "and auxiliary spaces for construction!");
      pc = MakeWrapperSolver<OperType, HypreAmsSolver>(
          linear, fespaces.GetFESpaceAtLevel(0), aux_fespaces->GetFESpaceAtLevel(0),
          coarse_solver ? linear.ams_max_it : linear.mg_cycle_it, linear.mg_smooth_it,
          linear.ams_vector_interp, linear.ams_singular_op, linear.amg_agg_coarsen, print);
      break;
    case LinearSolver::BOOMER_AMG:
      pc = MakeWrapperSolver<OperType, BoomerAmgSolver>(
          linear, coarse_solver ? 1 : linear.mg_cycle_it, linear.mg_smooth_it,
          linear.amg_agg_coarsen, print);
      break;
    case LinearSolver::SUPERLU:
#if defined(MFEM_USE_SUPERLU)
      pc = MakeWrapperSolver<OperType, SuperLUSolver>(
          linear, comm, linear.sym_factorization, linear.superlu_3d, linear.reorder_reuse,
          print);
#endif
      break;
    case LinearSolver::STRUMPACK:
#if defined(MFEM_USE_STRUMPACK)
      pc = MakeWrapperSolver<OperType, StrumpackSolver>(
          linear, comm, linear.sym_factorization, linear.strumpack_compression_type,
          linear.strumpack_lr_tol, linear.strumpack_butterfly_l,
          linear.strumpack_lossy_precision, linear.reorder_reuse, print);
#endif
      break;
    case LinearSolver::STRUMPACK_MP:
#if defined(MFEM_USE_STRUMPACK)
      pc = MakeWrapperSolver<OperType, StrumpackMixedPrecisionSolver>(
          linear, comm, linear.sym_factorization, linear.strumpack_compression_type,
          linear.strumpack_lr_tol, linear.strumpack_butterfly_l,
          linear.strumpack_lossy_precision, linear.reorder_reuse, print);
#endif
      break;
    case LinearSolver::MUMPS:
#if defined(MFEM_USE_MUMPS)
      pc = MakeWrapperSolver<OperType, MumpsSolver>(
          linear, comm, pc_mat_sym, linear.sym_factorization,
          (linear.strumpack_compression_type == SparseCompression::BLR)
              ? linear.strumpack_lr_tol
              : 0.0,
          linear.reorder_reuse, print);
#endif
      break;
    case LinearSolver::JACOBI:
      pc = std::make_unique<JacobiSmoother<OperType>>(comm);
      break;
    case LinearSolver::CUDSS:
#if defined(MFEM_USE_CUDSS)
      pc = MakeWrapperSolver<OperType, CuDSSSolver>(
          linear, comm, pc_mat_sym, linear.sym_factorization, linear.reorder_reuse, print);
#endif
      break;
    case LinearSolver::DEFAULT:
      MFEM_ABORT("Unexpected solver type for preconditioner configuration!");
      break;
  }

  // Construct the actual solver, which has the right value type.
  if (fespaces.GetNumLevels() > 1)
  {
    // This will construct the multigrid hierarchy using pc as the coarse solver
    // (ownership of pc is transferred to the GeometricMultigridSolver). When a special
    // auxiliary space smoother for pre-/post-smoothing is not desired, the auxiliary
    // space is a nullptr here.
    auto gmg = [&]()
    {
      if (linear.mg_smooth_aux)
      {
        MFEM_VERIFY(aux_fespaces, "Multigrid with auxiliary space smoothers requires both "
                                  "primary space and auxiliary spaces for construction!");
        const auto G = fespaces.GetDiscreteInterpolators(*aux_fespaces);
        return std::make_unique<GeometricMultigridSolver<OperType>>(
            comm, std::move(pc), fespaces.GetProlongationOperators(), &G,
            linear.mg_cycle_it, linear.mg_smooth_it, linear.mg_smooth_order,
            linear.mg_smooth_sf_max, linear.mg_smooth_sf_min, linear.mg_smooth_cheby_4th);
      }
      else
      {
        return std::make_unique<GeometricMultigridSolver<OperType>>(
            comm, std::move(pc), fespaces.GetProlongationOperators(), nullptr,
            linear.mg_cycle_it, linear.mg_smooth_it, linear.mg_smooth_order,
            linear.mg_smooth_sf_max, linear.mg_smooth_sf_min, linear.mg_smooth_cheby_4th);
      }
    }();
    gmg->EnableTimer();  // Enable timing for primary geometric multigrid solver
    if constexpr (std::is_same_v<OperType, ComplexOperator>)
    {
      if (linear.pml_subdomain_solver == PMLSubdomainSolver::DIRECT)
      {
        gmg->SetPMLSubdomainSolver(
            ConfigurePMLSubdomainSolver(linear, pc_mat_sym, verbose, comm));
      }
    }
    return gmg;
  }
  else
  {
    return pc;
  }
}

}  // namespace

MatrixSymmetry GetPreconditionerMatrixSymmetry(const IoData &iodata)
{
  // The Floquet terms are not symmetric, also for a shifted preconditioner matrix.
  const auto &linear = iodata.solver.linear;
  if (iodata.boundaries.periodic.wave_vector != std::array<double, 3>{0.0, 0.0, 0.0})
  {
    return MatrixSymmetry::UNSYMMETRIC;
  }
  if (linear.pc_mat_shifted || iodata.problem.type == ProblemType::TRANSIENT ||
      iodata.problem.type == ProblemType::ELECTROSTATIC ||
      iodata.problem.type == ProblemType::MAGNETOSTATIC)
  {
    return MatrixSymmetry::SPD;
  }
  return MatrixSymmetry::SYMMETRIC;
}

template <OperatorType OperType>
BaseKspSolver<OperType>::BaseKspSolver(const config::LinearSolverData &linear,
                                       MatrixSymmetry pc_mat_sym, int verbose,
                                       FiniteElementSpaceHierarchy &fespaces,
                                       FiniteElementSpaceHierarchy *aux_fespaces)
  : BaseKspSolver(ConfigureKrylovSolver<OperType>(linear, verbose,
                                                  fespaces.GetFinestFESpace().GetComm()),
                  ConfigurePreconditionerSolver<OperType>(
                      linear, pc_mat_sym, verbose, fespaces.GetFinestFESpace().GetComm(),
                      fespaces, aux_fespaces))
{
  use_timer = true;
}

template <OperatorType OperType>
BaseKspSolver<OperType>::BaseKspSolver(const IoData &iodata,
                                       FiniteElementSpaceHierarchy &fespaces,
                                       FiniteElementSpaceHierarchy *aux_fespaces)
  : BaseKspSolver(iodata.solver.linear, GetPreconditionerMatrixSymmetry(iodata),
                  iodata.problem.verbose, fespaces, aux_fespaces)
{
}

template <OperatorType OperType>
BaseKspSolver<OperType>::BaseKspSolver(std::unique_ptr<IterativeSolver<OperType>> &&ksp,
                                       std::unique_ptr<Solver<OperType>> &&pc)
  : ksp(std::move(ksp)), pc(std::move(pc)), ksp_mult(0), ksp_mult_it(0), use_timer(false)
{
  if (this->pc)
  {
    this->ksp->SetPreconditioner(*this->pc);
  }
}

template <OperatorType OperType>
void BaseKspSolver<OperType>::SetOperators(const OperType &op, const OperType &pc_op)
{
  BlockTimer bt(Timer::KSP_SETUP, use_timer);
  ksp->SetOperator(op);
  if (pc)
  {
    const auto *mg_op = dynamic_cast<const BaseMultigridOperator<OperType> *>(&pc_op);
    const auto *mg_pc = dynamic_cast<const GeometricMultigridSolver<OperType> *>(pc.get());
    if (mg_op && !mg_pc)
    {
      pc->SetOperator(mg_op->GetFinestOperator());
    }
    else
    {
      pc->SetOperator(pc_op);
    }
  }
}

template <OperatorType OperType>
void BaseKspSolver<OperType>::Mult(const VecType &x, VecType &y) const
{
  BlockTimer bt(Timer::KSP, use_timer);
  ksp->Mult(x, y);
  if (!ksp->GetConverged())
  {
    Mpi::Warning(
        ksp->GetComm(),
        "Linear solver did not converge, norm(Ax-b)/norm(b) = {:.3e} (norm(b) = {:.3e})!\n",
        ksp->GetFinalRes() / ksp->GetInitialRes(), ksp->GetInitialRes());
  }
  ksp_mult++;
  ksp_mult_it += ksp->GetNumIterations();
}

template class BaseKspSolver<Operator>;
template class BaseKspSolver<ComplexOperator>;

}  // namespace palace
