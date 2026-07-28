// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "ksp.hpp"

#include <stdexcept>
#include <mfem.hpp>
#include "fem/fespace.hpp"
#include "fem/singularsystem.hpp"
#include "linalg/amg.hpp"
#include "linalg/ams.hpp"
#include "linalg/blockprecond.hpp"
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

template <typename OperType>
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
                                    std::is_same<T, SuperLUSolver>::value ||
#endif
#if defined(MFEM_USE_STRUMPACK)
                                    std::is_same<T, StrumpackSolver>::value ||
                                    std::is_same<T, StrumpackMixedPrecisionSolver>::value ||
#endif
#if defined(MFEM_USE_MUMPS)
                                    std::is_same<T, MumpsSolver>::value ||
#endif
#if defined(MFEM_USE_CUDSS)
                                    std::is_same<T, CuDSSSolver>::value ||
#endif
                                    false);
  return std::make_unique<MfemWrapperSolver<OperType>>(
      std::make_unique<T>(std::forward<U>(args)...), save_assembled,
      linear.complex_coarse_solve, linear.drop_small_entries, linear.reorder_reuse);
}

template <typename OperType>
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
  // Mirrors the prior derivation that used to be stored on LinearSolverData::pc_mat_sym.
  const auto &linear = iodata.solver.linear;
  if (linear.pc_mat_shifted || iodata.problem.type == ProblemType::TRANSIENT ||
      iodata.problem.type == ProblemType::ELECTROSTATIC ||
      iodata.problem.type == ProblemType::MAGNETOSTATIC)
  {
    return MatrixSymmetry::SPD;
  }
  if (iodata.boundaries.periodic.wave_vector == std::array<double, 3>{0.0, 0.0, 0.0})
  {
    return MatrixSymmetry::SYMMETRIC;
  }
  return MatrixSymmetry::UNSYMMETRIC;
}

std::unique_ptr<KspSolver>
MakeSingularPatchKspSolver(const IoData &iodata, FiniteElementSpaceHierarchy &fespaces,
                           const Operator &full_operator, const Operator &standard_operator,
                           const fem::singular::ParallelFeaturePatches &feature_patches)
{
  const auto &linear = iodata.solver.linear;
  const MPI_Comm comm = fespaces.GetFinestFESpace().GetComm();
  if (!UsesSingularPatchKspSolver(iodata))
  {
#if !defined(MFEM_USE_SUPERLU)
    if (linear.type == LinearSolver::BOOMER_AMG)
    {
      Mpi::Warning(comm, "Singular feature-patch preconditioning requires SuperLU. "
                         "Falling back to the "
                         "configured monolithic preconditioner!\n");
    }
#endif
    return std::make_unique<KspSolver>(iodata, fespaces);
  }

#if defined(MFEM_USE_SUPERLU)
  auto standard_pc =
      std::make_unique<BoomerAmgSolver>(linear.mg_cycle_it, linear.mg_smooth_it,
                                        linear.amg_agg_coarsen, iodata.problem.verbose - 1);
  std::vector<mfem::Array<int>> patch_dofs;
  std::vector<std::unique_ptr<mfem::Solver>> patch_solvers;
  std::vector<const Operator *> patch_operators;
  patch_dofs.reserve(feature_patches.patches.size());
  patch_solvers.reserve(feature_patches.patches.size());
  patch_operators.reserve(feature_patches.patches.size());
  for (const auto &patch : feature_patches.patches)
  {
    patch_dofs.push_back(patch.true_dofs);
    patch_solvers.push_back(
        std::make_unique<SuperLUSolver>(comm, linear.sym_factorization, linear.superlu_3d,
                                        linear.reorder_reuse, iodata.problem.verbose - 1));
    patch_operators.push_back(patch.matrix.get());
  }
  auto patch_pc = std::make_unique<AdditivePatchSolver>(
      full_operator.Height(), std::move(patch_dofs), std::move(patch_solvers));
  patch_pc->SetPatchOperators(patch_operators);
  auto pc = std::make_unique<SymmetricPatchSubspacePreconditioner>(
      standard_operator.Height(), std::move(standard_pc), std::move(patch_pc));
  pc->SetSubspaceOperators(full_operator, standard_operator);

  auto result = std::make_unique<KspSolver>(
      ConfigureKrylovSolver<Operator>(linear, iodata.problem.verbose, comm), std::move(pc));
  result->EnableTimer();
  return result;
#else
  throw std::logic_error("Unreachable singular feature-patch solver configuration!");
#endif
}

bool UsesSingularPatchKspSolver(const IoData &iodata)
{
#if defined(MFEM_USE_SUPERLU)
  return iodata.solver.linear.type == LinearSolver::BOOMER_AMG;
#else
  return false;
#endif
}

namespace
{

template <typename OperType>
std::unique_ptr<BaseKspSolver<OperType>>
MakeSingularDirectKspSolverImpl(const IoData &iodata, MPI_Comm comm)
{
#if defined(MFEM_USE_SUPERLU)
  auto linear = iodata.solver.linear;
  linear.krylov_solver = KrylovSolver::GMRES;
  linear.initial_guess = false;
  if (linear.sym_factorization == SymbolicFactorization::DEFAULT && Mpi::Size(comm) > 1)
  {
    // SuperLU_DIST's default METIS_AT_PLUS_A ordering can produce an inaccurate
    // distributed factorization for indefinite combined singular-element systems.
    // Prefer the distributed ordering, while preserving an explicit user choice.
    linear.sym_factorization = SymbolicFactorization::PARMETIS;
  }
  auto pc = MakeWrapperSolver<OperType, SuperLUSolver>(
      linear, comm, linear.sym_factorization, linear.superlu_3d, linear.reorder_reuse,
      iodata.problem.verbose - 1);
  auto result = std::make_unique<BaseKspSolver<OperType>>(
      ConfigureKrylovSolver<OperType>(linear, iodata.problem.verbose, comm), std::move(pc));
  result->EnableTimer();
  return result;
#else
  MFEM_ABORT("Combined singular-element systems require a Palace build with "
             "SuperLU_DIST!");
  return nullptr;
#endif
}

}  // namespace

std::unique_ptr<KspSolver> MakeSingularDirectKspSolver(const IoData &iodata, MPI_Comm comm)
{
  return MakeSingularDirectKspSolverImpl<Operator>(iodata, comm);
}

std::unique_ptr<ComplexKspSolver> MakeSingularComplexKspSolver(const IoData &iodata,
                                                               MPI_Comm comm)
{
  return MakeSingularDirectKspSolverImpl<ComplexOperator>(iodata, comm);
}

template <typename OperType>
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

template <typename OperType>
BaseKspSolver<OperType>::BaseKspSolver(const IoData &iodata,
                                       FiniteElementSpaceHierarchy &fespaces,
                                       FiniteElementSpaceHierarchy *aux_fespaces)
  : BaseKspSolver(iodata.solver.linear, GetPreconditionerMatrixSymmetry(iodata),
                  iodata.problem.verbose, fespaces, aux_fespaces)
{
}

template <typename OperType>
BaseKspSolver<OperType>::BaseKspSolver(std::unique_ptr<IterativeSolver<OperType>> &&ksp,
                                       std::unique_ptr<Solver<OperType>> &&pc)
  : ksp(std::move(ksp)), pc(std::move(pc)), ksp_mult(0), ksp_mult_it(0),
    ksp_mult_failures(0), use_timer(false)
{
  if (this->pc)
  {
    this->ksp->SetPreconditioner(*this->pc);
  }
}

template <typename OperType>
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

template <typename OperType>
void BaseKspSolver<OperType>::Mult(const VecType &x, VecType &y) const
{
  BlockTimer bt(Timer::KSP, use_timer);
  ksp->Mult(x, y);
  if (!ksp->GetConverged())
  {
    ksp_mult_failures++;
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
