// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "ksp.hpp"

#include <cmath>
#include <type_traits>
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
#include "utils/perfprobes.hpp"
#include "utils/perfstats.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace
{

// Forwards to a preconditioner, recording an inclusive tally.
template <OperatorType OperType>
class TalliedSolver : public Solver<OperType>
{
  using VecType = typename Solver<OperType>::VecType;
  const Solver<OperType> &pc;
  const int id;

public:
  TalliedSolver(const Solver<OperType> &pc, int id) : pc(pc), id(id) {}
  void SetOperator(const OperType &) override {}
  void Mult(const VecType &x, VecType &y) const override
  {
    perf::ScopedTally t(id);
    pc.Mult(x, y);
  }
};

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
    return gmg;
  }
  else
  {
    // Single-level preconditioner: report its setup (on first application) separately.
    if (auto *wrapper = dynamic_cast<MfemWrapperSolver<OperType> *>(pc.get()))
    {
      wrapper->EnableTimer(true);
    }
    return pc;
  }
}

// Growth rate of the error under the error propagator I - M A of the preconditioner M for
// A: the geometric mean of the growth over the second half of a short power iteration. It
// approaches the spectral radius ρ(I - M A) if I - M A is normal, and is a growth rate
// (between ρ and the norm) otherwise. Not finite if the iteration overflows or breaks down.
template <OperatorType OperType>
double ErrorPropagatorGrowth(const OperType &A, const Solver<OperType> &M, MPI_Comm comm)
{
  using VecType =
      std::conditional_t<std::is_same_v<OperType, ComplexOperator>, ComplexVector, Vector>;
  constexpr int its = 16;
  VecType x(A.Height()), r(A.Height()), z(A.Height());
  x.UseDevice(true);
  r.UseDevice(true);
  z.UseDevice(true);
  if constexpr (std::is_same_v<OperType, ComplexOperator>)
  {
    perf::SetRandom(x.Real(), comm, 2);
    perf::SetRandom(x.Imag(), comm, 3);
  }
  else
  {
    perf::SetRandom(x, comm, 2);
  }
  x *= 1.0 / linalg::Norml2(comm, x);
  double log_growth = 0.0;
  for (int k = 0; k < its; k++)
  {
    A.Mult(x, r);
    z = 0.0;
    M.Mult(r, z);
    x -= z;
    const double norm = linalg::Norml2(comm, x);
    if (norm == 0.0)
    {
      return 0.0;  // Exact preconditioner
    }
    if (!std::isfinite(norm))
    {
      return norm;
    }
    x *= 1.0 / norm;
    if (k >= its / 2)
    {
      log_growth += std::log(norm);
    }
  }
  return std::exp(log_growth / (its - its / 2));
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

template <OperatorType OperType>
BaseKspSolver<OperType>::BaseKspSolver(const config::LinearSolverData &linear,
                                       MatrixSymmetry pc_mat_sym, int verbose,
                                       FiniteElementSpaceHierarchy &fespaces,
                                       FiniteElementSpaceHierarchy *aux_fespaces)
  : BaseKspSolver(ConfigureKrylovSolver<OperType>(linear, verbose,
                                                  fespaces.GetFinestFESpace().GetComm()),
                  ConfigurePreconditionerSolver<OperType>(
                      linear, pc_mat_sym, verbose, fespaces.GetFinestFESpace().GetComm(),
                      fespaces, aux_fespaces),
                  perf::CurrentRole())
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
                                       std::unique_ptr<Solver<OperType>> &&pc,
                                       std::string_view role)
  : ksp(std::move(ksp)), pc(std::move(pc)), ksp_mult(0), ksp_mult_it(0), use_timer(false),
    role(role), perf_solve_tally(perf::RegisterTally(fmt::format("{}/LinearSolve", role))),
    perf_setup_tally(perf::RegisterTally(fmt::format("{}/LinearSolverSetup", role))),
    perf_basis(perf::RegisterCounter("Memory/KrylovBasis/" + this->role + "/Megabytes",
                                     perf::Reduce::STATS, perf::Accumulate::MAX)),
    perf_stats(this->role)
{
  MFEM_VERIFY(!this->role.empty(),
              "Linear solvers need a role for the performance report!");
  if (this->pc)
  {
    // Nested solvers run without exclusive timers: tally their preconditioner instead (not
    // on GPUs, where these points do not synchronize). The multigrid's own tallies, set
    // below, include its coarse solve (<role>/MG/Coarse).
    if (this->role != "Main" && !mfem::Device::Allows(mfem::Backend::DEVICE_MASK))
    {
      pc_tally = std::make_unique<TalliedSolver<OperType>>(
          *this->pc, perf::RegisterTally(this->role + "/Preconditioner"));
      this->ksp->SetPreconditioner(*pc_tally);
    }
    else
    {
      this->ksp->SetPreconditioner(*this->pc);
    }
  }
  this->ksp->SetPerfTallies(this->role + "/Krylov");
  if (this->pc)
  {
    this->pc->SetPerfTallies(this->role);
  }
}

template <OperatorType OperType>
void BaseKspSolver<OperType>::SetOperators(const OperType &op, const OperType &pc_op)
{
  BlockTimer bt(Timer::KSP_SETUP, use_timer);
  perf::ScopedTally st(perf_setup_tally);
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
    // Only the main solver is checked, as for the coarse solve (see MfemWrapperSolver).
    check_op = (check_pc && perf::Verification() && role == "Main") ? &op : nullptr;
  }
}

template <OperatorType OperType>
void BaseKspSolver<OperType>::Mult(const VecType &x, VecType &y) const
{
  BlockTimer bt(Timer::KSP, use_timer);
  {
    perf::ScopedTally st(perf_solve_tally);
    ksp->Mult(x, y);
  }
  if (!ksp->GetConverged())
  {
    Mpi::Warning(
        ksp->GetComm(),
        "Linear solver did not converge, norm(Ax-b)/norm(b) = {:.3e} (norm(b) = {:.3e})!\n",
        ksp->GetFinalRes() / ksp->GetInitialRes(), ksp->GetInitialRes());
  }
  ksp_mult++;
  ksp_mult_it += ksp->GetNumIterations();
  if (const auto bytes = ksp->GetBasisBytes())
  {
    perf::SetMax(perf_basis, static_cast<double>(bytes) / (1024.0 * 1024.0));
  }
  perf_stats.Record(ksp->GetNumIterations(), ksp->GetConverged(), ksp->GetNumRestarts(),
                    ksp->GetInitialRes(), ksp->GetFinalRes());
  if (check_op && ksp->GetNumIterations() > 0)
  {
    // After a solve which applied the preconditioner, so that a lazy coarse factorization
    // is not part of the check, and charged to the verification row and its own tally only.
    double growth;
    {
      BlockTimer bt(Timer::VERIFICATION);
      perf::ScopedTally t(perf::RegisterTally(role + "/Verification/ErrorPropagator"));
      perf::ScopedPause pause;
      growth = ErrorPropagatorGrowth(*check_op, *pc, ksp->GetComm());
    }
    const std::string base = "Verification/Preconditioner/" + role;
    perf::SetMax(perf::RegisterCounter(base + "/MaxErrorPropagatorGrowth",
                                       perf::Reduce::MAX, perf::Accumulate::MAX),
                 growth);
    perf::Add(
        perf::RegisterCounter(base + "/Checks", perf::Reduce::MAX, perf::Accumulate::SUM),
        1.0);
    check_op = nullptr;
    check_pc = false;
  }
}

template class BaseKspSolver<Operator>;
template class BaseKspSolver<ComplexOperator>;

}  // namespace palace
