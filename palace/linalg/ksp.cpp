// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "ksp.hpp"

#include <algorithm>
#include <cmath>
#include <optional>

#include <mfem.hpp>
#include "fem/fespace.hpp"
#include "linalg/amg.hpp"
#include "linalg/ams.hpp"
#include "linalg/cudss.hpp"
#include "linalg/deflation.hpp"
#include "linalg/gmg.hpp"
#include "linalg/jacobi.hpp"
#include "linalg/mumps.hpp"
#include "linalg/strumpack.hpp"
#include "linalg/superlu.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/iodata.hpp"
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
    gmres->SetStagnation(linear.stagnation_it, linear.stagnation_tol);
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

// Sparse direct solver for the preconditioner, or nullptr for other preconditioner types.
template <OperatorType OperType>
std::unique_ptr<Solver<OperType>> MakeDirectSolver(const config::LinearSolverData &linear,
                                                   MatrixSymmetry pc_mat_sym, int print,
                                                   MPI_Comm comm)
{
  switch (linear.type)
  {
    case LinearSolver::SUPERLU:
#if defined(MFEM_USE_SUPERLU)
      return MakeWrapperSolver<OperType, SuperLUSolver>(
          linear, comm, linear.sym_factorization, linear.superlu_3d, linear.reorder_reuse,
          print);
#endif
      break;
    case LinearSolver::STRUMPACK:
#if defined(MFEM_USE_STRUMPACK)
      return MakeWrapperSolver<OperType, StrumpackSolver>(
          linear, comm, linear.sym_factorization, linear.strumpack_compression_type,
          linear.strumpack_lr_tol, linear.strumpack_butterfly_l,
          linear.strumpack_lossy_precision, linear.reorder_reuse,
          linear.strumpack_selective_inversion > 0, print);
#endif
      break;
    case LinearSolver::STRUMPACK_MP:
#if defined(MFEM_USE_STRUMPACK)
      return MakeWrapperSolver<OperType, StrumpackMixedPrecisionSolver>(
          linear, comm, linear.sym_factorization, linear.strumpack_compression_type,
          linear.strumpack_lr_tol, linear.strumpack_butterfly_l,
          linear.strumpack_lossy_precision, linear.reorder_reuse,
          linear.strumpack_selective_inversion > 0, print);
#endif
      break;
    case LinearSolver::MUMPS:
#if defined(MFEM_USE_MUMPS)
      return MakeWrapperSolver<OperType, MumpsSolver>(
          linear, comm, pc_mat_sym, linear.sym_factorization,
          (linear.strumpack_compression_type == SparseCompression::BLR)
              ? linear.strumpack_lr_tol
              : 0.0,
          linear.reorder_reuse, print);
#endif
      break;
    case LinearSolver::CUDSS:
#if defined(MFEM_USE_CUDSS)
      return MakeWrapperSolver<OperType, CuDSSSolver>(
          linear, comm, pc_mat_sym, linear.sym_factorization, linear.reorder_reuse, print);
#endif
      break;
    default:
      break;
  }
  return nullptr;
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
          coarse_solver ? linear.ams_max_it : linear.mg_cycle_it,
          coarse_solver ? linear.ams_smooth_it : linear.mg_smooth_it,
          linear.ams_vector_interp, linear.ams_singular_op, linear.amg_agg_coarsen, print);
      break;
    case LinearSolver::BOOMER_AMG:
      pc = MakeWrapperSolver<OperType, BoomerAmgSolver>(
          linear, coarse_solver ? 1 : linear.mg_cycle_it, linear.mg_smooth_it,
          linear.amg_agg_coarsen, print);
      break;
    case LinearSolver::SUPERLU:
    case LinearSolver::STRUMPACK:
    case LinearSolver::STRUMPACK_MP:
    case LinearSolver::MUMPS:
    case LinearSolver::CUDSS:
      pc = MakeDirectSolver<OperType>(linear, pc_mat_sym, print, comm);
      break;
    case LinearSolver::JACOBI:
      pc = std::make_unique<JacobiSmoother<OperType>>(comm);
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
            linear.mg_smooth_sf_max, linear.mg_smooth_sf_min, linear.mg_smooth_cheby_4th,
            linear.mg_smooth_order_max, linear.mg_smooth_order_aux, linear.mg_aux_residual);
      }
      else
      {
        return std::make_unique<GeometricMultigridSolver<OperType>>(
            comm, std::move(pc), fespaces.GetProlongationOperators(), nullptr,
            linear.mg_cycle_it, linear.mg_smooth_it, linear.mg_smooth_order,
            linear.mg_smooth_sf_max, linear.mg_smooth_sf_min, linear.mg_smooth_cheby_4th,
            linear.mg_smooth_order_max, linear.mg_smooth_order_aux, linear.mg_aux_residual);
      }
    }();
    gmg->EnableTimer();  // Enable timing for primary geometric multigrid solver

    // Without an explicit choice, a sparse direct coarse solver of a complex operator
    // factorizes the real-valued approximation, and the multigrid solver switches to the
    // complex-valued matrix if a check at setup finds the approximation poor (see IoData
    // for the configurations which allow it).
    if constexpr (std::is_same_v<OperType, ComplexOperator>)
    {
      // Memory of the complex-valued coarse solver, over its first and later factorizations
      // and at their peaks, relative to the memory of the setup of the real-valued
      // approximation (see GeometricMultigridSolver), with margins over the measured
      // values. Only the solvers listed by ChoosesComplexCoarseSolve use it.
      double memory_ratio = 0.0;
      switch (type)
      {
        case LinearSolver::STRUMPACK:
        case LinearSolver::STRUMPACK_MP:
          memory_ratio = 6.0;
          break;
        case LinearSolver::SUPERLU:
        case LinearSolver::MUMPS:
          memory_ratio = 8.5;
          break;
        default:
          break;
      }
      if (linear.complex_coarse_check && !linear.complex_coarse_solve && memory_ratio > 0.0)
      {
        // As with an explicit complex-valued coarse solve: the auxiliary space smoothers
        // which follow the order of the primary space then use max(2 p, 4), and the
        // reordering reuse follows the same rule (see IoData).
        auto complex_linear = linear;
        complex_linear.complex_coarse_solve = true;
        complex_linear.reorder_reuse = linear.complex_coarse_reorder_reuse;
        const int p = fespaces.GetFinestFESpace().GetMaxElementOrder();
        gmg->SetComplexCoarseSolver(
            [complex_linear, pc_mat_sym, print, comm]()
            { return MakeDirectSolver<OperType>(complex_linear, pc_mat_sym, print, comm); },
            (linear.mg_smooth_order_aux > 0) ? 0 : std::max(2 * p, 4), memory_ratio);
      }
    }

    // By default, sliver patch solves for complex-valued operators only: for real-valued
    // ones, with their symmetric V-cycles, the patch correction does not speed up CG.
    const double sliver_angle =
        (linear.mg_sliver_angle < 0.0)
            ? (std::is_same_v<OperType, ComplexOperator> ? 2.0 : 0.0)
            : linear.mg_sliver_angle;
    if (sliver_angle > 0.0)
    {
      gmg->SetSliverPatches(sliver_angle, linear.mg_sliver_layers);
    }
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

// Exclusive timer row containing the untimed solves of a given role.
Timer::Index RoleParent(std::string_view role)
{
  if (role == "DivFree")
  {
    return Timer::DIV_FREE;
  }
  if (role == "Estimator")
  {
    return Timer::SOLVE_ESTIMATOR;
  }
  if (role == "Mass")
  {
    return Timer::TS;
  }
  return Timer::TOTAL;
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
                      fespaces, aux_fespaces))
{
  use_timer = true;

  // An inconsistent discretization can make a multigrid preconditioner unstable, and GMRES
  // with left preconditioning then reports convergence to an invalid solution: check the
  // true residual. Not with a sparse direct preconditioner, whose solves of nearly singular
  // (shift-and-invert) systems leave large true residuals.
  auto *gmres = dynamic_cast<GmresSolver<OperType> *>(ksp.get());
  if (gmres &&
      (dynamic_cast<const GeometricMultigridSolver<OperType> *>(pc.get()) ||
       linear.type == LinearSolver::AMS || linear.type == LinearSolver::BOOMER_AMG))
  {
    gmres->SetTrueResidualCheck(0.5);
  }
}

template <OperatorType OperType>
BaseKspSolver<OperType>::BaseKspSolver(const IoData &iodata,
                                       FiniteElementSpaceHierarchy &fespaces,
                                       FiniteElementSpaceHierarchy *aux_fespaces)
  : BaseKspSolver(iodata.solver.linear, GetPreconditionerMatrixSymmetry(iodata),
                  iodata.problem.verbose, fespaces, aux_fespaces)
{
  // The auxiliary space check of the single-precision smoother residuals is for driven and
  // transient solves: the shift-and-invert solves of eigenmode simulations converge with
  // the lagged residuals where it fails (see DistRelaxationSmoother).
  if (iodata.problem.type == ProblemType::EIGENMODE)
  {
    if (auto *gmg = dynamic_cast<GeometricMultigridSolver<OperType> *>(pc.get()))
    {
      gmg->SetAuxiliaryGradientCheck(false);
    }
  }
}

template <OperatorType OperType>
BaseKspSolver<OperType>::BaseKspSolver(std::unique_ptr<IterativeSolver<OperType>> &&ksp,
                                       std::unique_ptr<Solver<OperType>> &&pc,
                                       std::string_view role)
  : ksp(std::move(ksp)), pc(std::move(pc)), ksp_mult(0), ksp_mult_it(0), use_timer(false),
    role(role.empty() ? perf::CurrentRole() : role),
    perf_solve_tally(
        perf::RegisterTally(fmt::format("{}/LinearSolve", this->role), Timer::KSP)),
    perf_setup_tally(perf::RegisterTally(fmt::format("{}/LinearSolverSetup", this->role),
                                         Timer::KSP_SETUP)),
    perf_basis(perf::RegisterCounter("Memory/KrylovBasis/" + this->role + "/Megabytes",
                                     perf::Reduce::STATS, perf::Accumulate::MAX)),
    perf_stats(this->role)
{
  if (this->pc)
  {
    // Nested solvers run without exclusive timers: tally their preconditioner and coarse
    // solve instead (not on GPUs, where these points do not synchronize).
    if (this->role != "Main" && !mfem::Device::Allows(mfem::Backend::DEVICE_MASK))
    {
      const auto parent = perf::RoleParent(this->role);
      if (auto *gmg = dynamic_cast<GeometricMultigridSolver<OperType> *>(this->pc.get()))
      {
        gmg->SetCoarseTally(perf::RegisterTally(this->role + "/MG/Coarse", parent));
      }
      pc_tally = std::make_unique<TalliedSolver<OperType>>(
          *this->pc, perf::RegisterTally(this->role + "/Preconditioner", parent));
      this->ksp->SetPreconditioner(*pc_tally);
    }
    else
    {
      this->ksp->SetPreconditioner(*this->pc);
    }
  }
  const bool main = (this->role == "Main");
  const auto parent = main ? Timer::KSP : RoleParent(this->role);
  this->ksp->SetPerfTallies(this->role + "/Krylov", parent);
  if (this->pc)
  {
    this->pc->SetPerfTallies(this->role, main ? Timer::KSP_PRECONDITIONER : parent);
  }
}

template <OperatorType OperType>
BaseKspSolver<OperType>::~BaseKspSolver() = default;

template <OperatorType OperType>
const Solver<OperType> &BaseKspSolver<OperType>::GetBasePreconditioner() const
{
  return pc_tally ? *pc_tally : *pc;
}

template <OperatorType OperType>
void BaseKspSolver<OperType>::SetDeflation(int k)
{
  if constexpr (std::is_same_v<OperType, ComplexOperator>)
  {
    const auto *gmres = dynamic_cast<const GmresSolver<OperType> *>(ksp.get());
    if (k != 0 && (!pc || !gmres || dynamic_cast<const FgmresSolver<OperType> *>(gmres) ||
                   gmres->GetPreconditionerSide() != PreconditionerSide::LEFT))
    {
      Mpi::Warning(ksp->GetComm(), "Linear solver deflation is only available for GMRES "
                                   "with left preconditioning and is ignored!\n");
      k = 0;
    }
    defl_size = k;
    if (k != 0 && !defl)
    {
      defl = std::make_unique<DeflatedPreconditioner>(ksp->GetComm());
      defl_levels = std::make_unique<DeflationLevels>();
    }
    if (defl && pc)
    {
      defl->Clear();
      ksp->SetPreconditioner(GetBasePreconditioner());
    }
    defl_build = (k != 0);
    defl_off = defl_confirm = false;
  }
  else
  {
    MFEM_VERIFY(k == 0, "Linear solver deflation is only available for complex operators!");
  }
}

template <OperatorType OperType>
void BaseKspSolver<OperType>::SetInitialGuessCheck(double tol)
{
  if (auto *gmres = dynamic_cast<GmresSolver<OperType> *>(ksp.get()))
  {
    gmres->SetInitialGuessCheck(tol);
  }
}

template <OperatorType OperType>
int BaseKspSolver<OperType>::GetDeflationSize() const
{
  return defl ? defl->Size() : 0;
}

template <OperatorType OperType>
void BaseKspSolver<OperType>::BuildDeflation(double t_it) const
{
  if constexpr (std::is_same_v<OperType, ComplexOperator>)
  {
    // Fixed number of vectors, or the automatic choice among the candidates: the model
    // iteration counts with deflation, over those without, times the relative cost per
    // iteration of the deflation (its kernels measured on the Krylov basis, a global
    // reduction per application). Only outliers of B A (harmonic Ritz values with
    // |θ| < 0.25, B A being close to the identity elsewhere) are deflated at first:
    // deflating values of the bulk of the spectrum often barely shortens the solves (the
    // projected model overestimates their effect), while the deflated stopping norm,
    // stricter on the deflated directions, can lengthen them. A choice which saves less
    // than 5% is not taken. The basis also holds the vectors of the values below 0.5 and 1,
    // tried while the deflation pays (the Krylov subspace of the first solve, from a
    // generic right-hand side, approximates them better than those of the later solves).
    constexpr int kmax_auto = 64, min_its = 8;
    constexpr double margin = 0.95, outlier = 0.25, confident = 0.6, grow = 0.9;
    const auto &gmres = static_cast<const GmresSolver<OperType> &>(*ksp);
    const std::vector<VecType> *V;
    const std::complex<double> *H;
    int ld;
    const int m = gmres.GetArnoldiRelation(V, H, ld);
    if (m <= 0)
    {
      return;  // No Arnoldi relation (zero right-hand side): wait for the next solve.
    }
    defl_build = false;
    defl->Clear();
    ksp->SetPreconditioner(GetBasePreconditioner());
    const bool automatic = (defl_size < 0);
    if (m < min_its && automatic)
    {
      return;
    }
    MPI_Comm comm = ksp->GetComm();
    const int kmax = automatic ? std::min(kmax_auto, m / 2) : std::min(defl_size, m);
    std::vector<int> cand;
    if (automatic)
    {
      for (int k : {4, 8, 12, 16, 24, 32, 48, 64})
      {
        if (k <= kmax)
        {
          cand.push_back(k);
        }
      }
    }
    DeflationAnalysis analysis;
    analysis.Compute(comm, H, ld, m, kmax, cand, ksp->GetRelTol(),
                     automatic ? outlier : 0.0);
    const auto [t_vec, t_red] = MeasureDeflationCost(comm, *V, m);
    defl_cost_vec = t_vec / t_it;
    defl_cost_red = t_red / t_it;
    int k = kmax;
    double ratio = 1.0;
    if (automatic)
    {
      k = 0;
      double best = margin;
      for (std::size_t c = 0; c < analysis.k_candidates.size(); c++)
      {
        const int kc = analysis.k_candidates[c];
        const double r = static_cast<double>(analysis.its_deflated[c]) /
                         std::max(analysis.its_plain, 1) * DeflationCost(kc);
        if (kc <= analysis.num_outliers && r < best)
        {
          best = r;
          k = kc;
          ratio = static_cast<double>(analysis.its_deflated[c]) /
                  std::max(analysis.its_plain, 1);
        }
      }
    }
    if (k > 0)
    {
      // For the automatic choice, the basis holds the vectors of the values below 1, and
      // the nested sets are those below 0.25 (the choice), 0.5 and 1, each tried if the
      // model predicts at least 10% fewer iterations than with the previous one.
      std::vector<int> sets = {k};
      if (automatic)
      {
        auto Model = [&](int kc)
        {
          for (std::size_t c = 0; c < analysis.k_candidates.size(); c++)
          {
            if (analysis.k_candidates[c] == kc)
            {
              return analysis.its_deflated[c];
            }
          }
          return analysis.m + 1;
        };
        for (double f : {2.0, 4.0})
        {
          const int kc = analysis.Count(f * outlier);
          if (kc > sets.back() && Model(kc) <= grow * Model(sets.back()))
          {
            sets.push_back(kc);
          }
        }
      }
      const int kr =
          defl->Build(*defl_op, GetBasePreconditioner(), *V, analysis, sets.back());
      const int ka = defl->Activate(defl->Prefix(k));
      if (ka > 0)
      {
        ksp->SetPreconditioner(*defl);
        defl_confident = (ratio < confident);
        defl_solves = 0;
        defl_its.clear();
        defl_next_check = 8;
        if (automatic)
        {
          std::vector<int> sizes;
          for (int kc : sets)
          {
            sizes.push_back(defl->Prefix(kc));
          }
          defl_levels->Reset(sizes);
          Mpi::Print(
              comm,
              " Deflation: {:d} harmonic Ritz vectors (|theta| <= {:.3f}, {:d} below "
              "{:.2f}), predicted iterations {:.2f}x, cost {:.2f}x per iteration ({:d} "
              "more in reserve)\n",
              ka, std::abs(analysis.theta[k - 1]), analysis.num_outliers, outlier, ratio,
              DeflationCost(ka), kr - ka);
        }
        return;
      }
      defl->Clear();
    }
    if (automatic)
    {
      double r = 1.0;
      int kc = 0;
      for (std::size_t c = 0; c < analysis.k_candidates.size(); c++)
      {
        const double rc =
            static_cast<double>(analysis.its_deflated[c]) / std::max(analysis.its_plain, 1);
        if (analysis.k_candidates[c] <= analysis.num_outliers && rc < r)
        {
          r = rc;
          kc = analysis.k_candidates[c];
        }
      }
      Mpi::Print(comm,
                 " Deflation: none ({:d} harmonic Ritz values below {:.2f}; best predicted "
                 "iterations {:.2f}x with {:d} vectors, cost {:.2f}x per iteration)\n",
                 analysis.num_outliers, outlier, r, kc, DeflationCost(kc));
    }
  }
}

template <OperatorType OperType>
double BaseKspSolver<OperType>::DeflationCost(int k) const
{
  return 1.0 + (k > 0 ? k * defl_cost_vec + defl_cost_red : 0.0);
}

template <OperatorType OperType>
void BaseKspSolver<OperType>::CheckDeflation(const VecType &b, int its,
                                             DeflationCheckKind kind) const
{
  if constexpr (std::is_same_v<OperType, ComplexOperator>)
  {
    // A set of vectors pays when it saves at least 5% against no deflation, and a larger
    // set when it saves more than 3% against the previous one, at the cost of the deflation
    // per iteration. A verification compares the solve without deflation with the last
    // deflated ones. Their right-hand sides differ: the iterations without deflation vary
    // with the content of the right-hand side along the slow directions (by 20% or more),
    // those with deflation much less, and on a right-hand side with little of it the
    // deflation gains little. So when a verification finds no gain while the model
    // predicted a large one, the next right-hand side, after its deflated solve, is solved
    // without deflation up to the iterations at which the deflation pays, which decides
    // (confirmation). A comparison solves the same right-hand side with the previous set,
    // up to the iterations at which the current set pays.
    constexpr double margin = 0.95, again = 0.97;
    const int k = defl->Size();
    std::string msg;
    DeflationLevels::Action action;
    if (kind == DeflationCheckKind::VERIFICATION)
    {
      double mean = 0.0;
      for (int i : defl_its)
      {
        mean += i;
      }
      mean /= std::max<std::size_t>(defl_its.size(), 1);
      bool gain = (mean * DeflationCost(k) < margin * its);
      msg = fmt::format(
          " Deflation check: {:d} iterations without deflation, {:.1f} with {:d} "
          "vectors (cost {:.2f}x per iteration)",
          its, mean, k, DeflationCost(k));
      if (!gain && defl_confident)
      {
        defl_confirm = true;
        ksp->SetPreconditioner(*defl);
        Mpi::Print(ksp->GetComm(), "{}: confirming on the next solve\n", msg);
        return;
      }
      action = defl_levels->Verify(gain);
    }
    else if (kind == DeflationCheckKind::CONFIRMATION)
    {
      defl_confirm = false;
      const int limit = static_cast<int>(its * DeflationCost(k) / margin);
      const int its_ref = ReferenceSolve(b, 0, limit);
      const bool gain = (its_ref > limit);
      msg = fmt::format(" Deflation check: {:d} iterations with {:d} vectors (cost {:.2f}x "
                        "per iteration), {} without deflation",
                        its, k, DeflationCost(k),
                        gain ? fmt::format("over {:d}", limit) : std::to_string(its_ref));
      action = defl_levels->Verify(gain);
    }
    else
    {
      const int kp = defl_levels->PreviousSize();
      const int limit =
          static_cast<int>(its * DeflationCost(k) / (again * DeflationCost(kp)));
      const int its_ref = ReferenceSolve(b, kp, limit);
      const bool gain = (its_ref > limit);
      msg = fmt::format(" Deflation check: {:d} iterations with {:d} vectors (cost {:.2f}x "
                        "per iteration), {} with {:d} (cost {:.2f}x)",
                        its, k, DeflationCost(k),
                        gain ? fmt::format("over {:d}", limit) : std::to_string(its_ref),
                        kp, DeflationCost(kp));
      action = defl_levels->Compare(gain);
    }
    switch (action)
    {
      case DeflationLevels::Action::STOP:
        defl_off = true;
        defl->Clear();
        msg += ": no gain, disabled";
        break;
      case DeflationLevels::Action::BACK:
        defl->Activate(defl_levels->Size());
        msg += fmt::format(": no gain, back to {:d} vectors", defl->Size());
        break;
      case DeflationLevels::Action::GROW:
        if (defl->Activate(defl_levels->Size()) < defl_levels->Size())
        {
          defl_levels->Compare(false);
          defl->Activate(defl_levels->Size());
          msg += fmt::format(": singular for more vectors, kept {:d}", defl->Size());
        }
        else
        {
          msg += fmt::format(": trying {:d} vectors", defl->Size());
        }
        break;
      case DeflationLevels::Action::KEEP:
        break;
    }
    if (defl_off)
    {
      ksp->SetPreconditioner(GetBasePreconditioner());
    }
    else
    {
      ksp->SetPreconditioner(*defl);
    }
    defl_solves = 0;
    defl_its.clear();
    defl_next_check = 32;
    Mpi::Print(ksp->GetComm(), "{}\n", msg);
  }
}

template <OperatorType OperType>
int BaseKspSolver<OperType>::ReferenceSolve(const VecType &b, int k, int limit) const
{
  if constexpr (std::is_same_v<OperType, ComplexOperator>)
  {
    if (limit <= 0)
    {
      return limit + 1;
    }
    const int k0 = defl->Size();
    if (k > 0)
    {
      defl->Activate(k);
      ksp->SetPreconditioner(*defl);
    }
    else
    {
      ksp->SetPreconditioner(GetBasePreconditioner());
    }
    const int max_it = ksp->GetMaxIter();
    const auto print = ksp->GetPrintLevel();
    ksp->SetMaxIter(limit);
    ksp->SetPrintLevel(mfem::IterativeSolver::PrintLevel().None());
    VecType z(b.Size());
    z.UseDevice(true);
    z = 0.0;
    {
      perf::ScopedTally st(perf_solve_tally);
      ksp->Mult(b, z);
    }
    const int its = ksp->GetNumIterations();
    const bool converged = ksp->GetConverged();
    ksp->SetMaxIter(max_it);
    ksp->SetPrintLevel(print);
    defl->Activate(k0);
    ksp_mult_it += its;
    return converged ? its : limit + 1;
  }
  return 0;
}

template <OperatorType OperType>
void BaseKspSolver<OperType>::SetOperators(const OperType &op, const OperType &pc_op,
                                           bool reuse_coarse)
{
  BlockTimer bt(Timer::KSP_SETUP, use_timer);
  perf::ScopedTally st(perf_setup_tally);
  ksp->SetOperator(op);
  defl_op = &op;
  if (defl && pc)
  {
    defl->Clear();
    ksp->SetPreconditioner(GetBasePreconditioner());
    defl_build = (defl_size != 0);
    defl_off = defl_confirm = false;
  }
  if (pc)
  {
    const auto *mg_op = dynamic_cast<const BaseMultigridOperator<OperType> *>(&pc_op);
    auto *mg_pc = dynamic_cast<GeometricMultigridSolver<OperType> *>(pc.get());
    if (mg_op && !mg_pc)
    {
      pc->SetOperator(mg_op->GetFinestOperator());
    }
    else if (mg_pc)
    {
      mg_pc->SetOperator(pc_op, reuse_coarse);
    }
    else
    {
      pc->SetOperator(pc_op);
    }
  }
}

template <OperatorType OperType>
bool BaseKspSolver<OperType>::HasMultigridPreconditioner() const
{
  return dynamic_cast<const GeometricMultigridSolver<OperType> *>(pc.get()) != nullptr;
}

template <OperatorType OperType>
void BaseKspSolver<OperType>::Mult(const VecType &x, VecType &y,
                                   bool warn_unconverged) const
{
  BlockTimer bt(Timer::KSP, use_timer);
  // Automatic deflation: after some deflated solves, a solve without deflation verifies
  // that the deflation pays, and the solve after a larger set of vectors is activated
  // compares it with the previous set.
  std::optional<DeflationCheckKind> check;
  if constexpr (std::is_same_v<OperType, ComplexOperator>)
  {
    if (defl && defl->Size() > 0 && defl_size < 0)
    {
      if (defl_levels->Trying())
      {
        check = DeflationCheckKind::COMPARISON;
      }
      else if (defl_confirm)
      {
        check = DeflationCheckKind::CONFIRMATION;
      }
      else if (defl_solves >= defl_next_check)
      {
        check = DeflationCheckKind::VERIFICATION;
        ksp->SetPreconditioner(GetBasePreconditioner());
      }
    }
  }
  const double t_start = defl_build ? MPI_Wtime() : 0.0;
  {
    perf::ScopedTally st(perf_solve_tally);
    ksp->Mult(x, y);
  }
  const double t_solve = defl_build ? MPI_Wtime() - t_start : 0.0;
  const double rel_tol = ksp->GetRelTol();
  if (ksp->GetFalseConverged() && fallback_rel_tol > 0.0 && fallback_rel_tol < rel_tol)
  {
    Mpi::Warning(ksp->GetComm(),
                 "Repeating the linear solve with the relative tolerance {:.1e}!\n",
                 fallback_rel_tol);
    ksp_mult_it += ksp->GetNumIterations();
    ksp->SetRelTol(fallback_rel_tol);
    {
      perf::ScopedTally st(perf_solve_tally);
      ksp->Mult(x, y);
    }
    ksp->SetRelTol(rel_tol);
  }
  MFEM_VERIFY(
      allow_false_convergence || !ksp->GetFalseConverged(),
      "Linear solver false convergence: the preconditioned residual converged but the "
      "true residual did not, so the solution is invalid (unstable preconditioner)!");
  const int its = ksp->GetNumIterations(), restarts = ksp->GetNumRestarts();
  const bool converged = ksp->GetConverged();
  const double res0 = ksp->GetInitialRes(), res = ksp->GetFinalRes();
  false_converged = ksp->GetFalseConverged();
  if constexpr (std::is_same_v<OperType, ComplexOperator>)
  {
    if (defl_build && !defl_off)
    {
      // The first solve with these operators which iterates provides the deflation basis.
      BuildDeflation(t_solve / std::max(its, 1));
    }
    else if (check)
    {
      CheckDeflation(x, its, *check);
    }
    else if (defl && defl->Size() > 0)
    {
      defl_solves++;
      defl_its.push_back(its);
      if (defl_its.size() > 4)
      {
        defl_its.erase(defl_its.begin());
      }
    }
  }
  if (!converged && warn_unconverged)
  {
    Mpi::Warning(
        ksp->GetComm(),
        "Linear solver did not converge, norm(Ax-b)/norm(b) = {:.3e} (norm(b) = {:.3e})!\n",
        res / res0, res0);
  }
  ksp_mult++;
  ksp_mult_it += its;
  if (const auto bytes = ksp->GetBasisBytes() + (defl ? defl->GetBytes() : 0))
  {
    perf::SetMax(perf_basis, static_cast<double>(bytes) / (1024.0 * 1024.0));
  }
  perf_stats.Record(its, converged, restarts, res0, res);
}

template class BaseKspSolver<Operator>;
template class BaseKspSolver<ComplexOperator>;

void CoarseReuse::Solve(ComplexKspSolver &ksp, const ComplexOperator &A,
                        const ComplexOperator &P, double omega, const ComplexVector &b,
                        ComplexVector &x, bool initial_guess)
{
  // Only a multigrid preconditioner has a coarse level to keep: other preconditioners are
  // set up at every frequency, and their solves are not limited.
  const bool reuse = (tol > 0.0 && valid && ksp.HasMultigridPreconditioner() &&
                      std::abs(omega - omega_ref) <= tol * std::abs(omega_ref));
  ksp.SetOperators(A, P, reuse);
  const int it0 = ksp.NumTotalMultIterations();
  if (!reuse)
  {
    ksp.Mult(b, x);
    omega_ref = omega;
    its_ref = std::max(its_ref, ksp.NumTotalMultIterations() - it0);
    valid = true;
    return;
  }
  const int max_it = ksp.GetMaxIter();
  const int max_it_reuse =
      std::min(max_it, std::max(its_ref + 10, static_cast<int>((1.0 + growth) * its_ref)));
  ksp.SetMaxIter(max_it_reuse);
  ksp.Mult(b, x, false);
  ksp.SetMaxIter(max_it);
  if (!ksp.GetConverged())
  {
    Mpi::Print(" Reused coarse factorization did not converge in {:d} iterations, "
               "factoring at this frequency\n",
               max_it_reuse);
    ksp.SetOperators(A, P, false);
    ksp.SetInitialGuess(true);
    ksp.Mult(b, x);
    ksp.SetInitialGuess(initial_guess);
    omega_ref = omega;
  }
}

}  // namespace palace
