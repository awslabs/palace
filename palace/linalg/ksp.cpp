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
  if (linear.pc_side != PreconditionerSide::DEFAULT && type != KrylovSolver::GMRES &&
      type != KrylovSolver::FGMRES)
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

std::unique_ptr<Solver<Operator>>
ConfigureSingularElectrostaticCoarseSolver(const config::LinearSolverData &linear,
                                           MatrixSymmetry matrix_symmetry, int verbose,
                                           MPI_Comm comm)
{
  const int print = verbose - 1;
  switch (linear.type)
  {
    case LinearSolver::BOOMER_AMG:
      return MakeWrapperSolver<Operator, BoomerAmgSolver>(linear, 1, linear.mg_smooth_it,
                                                          linear.amg_agg_coarsen, print);
    case LinearSolver::SUPERLU:
#if defined(MFEM_USE_SUPERLU)
      return MakeWrapperSolver<Operator, SuperLUSolver>(
          linear, comm, linear.sym_factorization, linear.superlu_3d, linear.reorder_reuse,
          print);
#else
      MFEM_ABORT("Solver was not built with SuperLU_DIST support!");
      return {};
#endif
    case LinearSolver::STRUMPACK:
#if defined(MFEM_USE_STRUMPACK)
      return MakeWrapperSolver<Operator, StrumpackSolver>(
          linear, comm, linear.sym_factorization, linear.strumpack_compression_type,
          linear.strumpack_lr_tol, linear.strumpack_butterfly_l,
          linear.strumpack_lossy_precision, linear.reorder_reuse, print);
#else
      MFEM_ABORT("Solver was not built with STRUMPACK support!");
      return {};
#endif
    case LinearSolver::STRUMPACK_MP:
#if defined(MFEM_USE_STRUMPACK)
      return MakeWrapperSolver<Operator, StrumpackMixedPrecisionSolver>(
          linear, comm, linear.sym_factorization, linear.strumpack_compression_type,
          linear.strumpack_lr_tol, linear.strumpack_butterfly_l,
          linear.strumpack_lossy_precision, linear.reorder_reuse, print);
#else
      MFEM_ABORT("Solver was not built with STRUMPACK support!");
      return {};
#endif
    case LinearSolver::MUMPS:
#if defined(MFEM_USE_MUMPS)
      return MakeWrapperSolver<Operator, MumpsSolver>(
          linear, comm, matrix_symmetry, linear.sym_factorization,
          (linear.strumpack_compression_type == SparseCompression::BLR)
              ? linear.strumpack_lr_tol
              : 0.0,
          linear.reorder_reuse, print);
#else
      MFEM_ABORT("Solver was not built with MUMPS support!");
      return {};
#endif
    case LinearSolver::CUDSS:
#if defined(MFEM_USE_CUDSS)
      return MakeWrapperSolver<Operator, CuDSSSolver>(linear, comm, matrix_symmetry,
                                                      linear.sym_factorization,
                                                      linear.reorder_reuse, print);
#else
      MFEM_ABORT("Solver was not built with cuDSS support!");
      return {};
#endif
    case LinearSolver::JACOBI:
      return std::make_unique<JacobiSmoother<Operator>>(comm);
    case LinearSolver::AMS:
      MFEM_ABORT("Electrostatic singular multigrid does not support an AMS coarse solver!");
      return {};
    case LinearSolver::DEFAULT:
      MFEM_ABORT("Unexpected unresolved coarse solver type for singular electrostatics!");
      return {};
  }
  MFEM_ABORT("Unsupported coarse solver type for singular electrostatics!");
  return {};
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
                           const fem::singular::ParallelFeaturePatches &feature_patches,
                           const std::vector<const Operator *> &combined_prolongations,
                           const std::vector<mfem::Array<int>> &combined_essential_tdofs)
{
  const auto &linear = iodata.solver.linear;
  const MPI_Comm comm = fespaces.GetFinestFESpace().GetComm();
  const auto number_levels = fespaces.GetNumLevels();
  const auto *hierarchy = dynamic_cast<const MultigridOperator *>(&full_operator);
  MFEM_VERIFY(number_levels > 0 && combined_prolongations.size() + 1 == number_levels &&
                  combined_essential_tdofs.size() == number_levels &&
                  ((number_levels == 1 && !hierarchy) ||
                   (hierarchy && hierarchy->GetNumLevels() == number_levels)),
              "Singular electrostatic solver received an inconsistent multigrid "
              "hierarchy!");
  const Operator &coarse_operator =
      hierarchy ? hierarchy->GetOperatorAtLevel(0) : full_operator;

  std::unique_ptr<Solver<Operator>> coarse_solver;
  if (UsesSingularPatchKspSolver(iodata))
  {
#if defined(MFEM_USE_SUPERLU)
    auto standard_pc = std::make_unique<BoomerAmgSolver>(
        number_levels > 1 ? 1 : linear.mg_cycle_it, linear.mg_smooth_it,
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
      patch_solvers.push_back(std::make_unique<SuperLUSolver>(
          comm, linear.sym_factorization, linear.superlu_3d, linear.reorder_reuse,
          iodata.problem.verbose - 1));
      patch_operators.push_back(patch.matrix.get());
    }
    auto patch_pc = std::make_unique<AdditivePatchSolver>(
        coarse_operator.Height(), std::move(patch_dofs), std::move(patch_solvers));
    patch_pc->SetPatchOperators(patch_operators);
    auto patch_solver = std::make_unique<SymmetricPatchSubspacePreconditioner>(
        standard_operator.Height(), std::move(standard_pc), std::move(patch_pc));
    patch_solver->SetSubspaceOperators(coarse_operator, standard_operator);
    coarse_solver = std::move(patch_solver);
#else
    throw std::logic_error("Unreachable singular feature-patch solver configuration!");
#endif
  }
  else
  {
    coarse_solver = ConfigureSingularElectrostaticCoarseSolver(
        linear, GetPreconditionerMatrixSymmetry(iodata), iodata.problem.verbose, comm);
  }

  std::unique_ptr<Solver<Operator>> preconditioner;
  if (number_levels == 1)
  {
    preconditioner = std::move(coarse_solver);
  }
  else
  {
    MFEM_VERIFY(!linear.mg_smooth_aux,
                "Singular electrostatic multigrid does not use an auxiliary-space "
                "smoother!");
    std::vector<const mfem::Array<int> *> essential_tdofs(number_levels);
    for (std::size_t level = 0; level < number_levels; level++)
    {
      essential_tdofs[level] = &combined_essential_tdofs[level];
    }
    auto multigrid = std::make_unique<GeometricMultigridSolver<Operator>>(
        comm, std::move(coarse_solver), combined_prolongations, nullptr, linear.mg_cycle_it,
        linear.mg_smooth_it, linear.mg_smooth_order, linear.mg_smooth_sf_max,
        linear.mg_smooth_sf_min, linear.mg_smooth_cheby_4th, &essential_tdofs);
    multigrid->EnableTimer();
    preconditioner = std::move(multigrid);
  }

  auto result = std::make_unique<KspSolver>(
      ConfigureKrylovSolver<Operator>(linear, iodata.problem.verbose, comm),
      std::move(preconditioner));
  result->EnableTimer();
  return result;
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

namespace
{

class SingularAmsCoarseSolver : public Solver<ComplexOperator>
{
private:
  int standard_size;
  std::unique_ptr<Solver<ComplexOperator>> standard_solver, enrichment_solver;
  const ComplexOperator *full_operator = nullptr;
  std::unique_ptr<ComplexOperator> approximate_operator, standard_operator,
      enrichment_operator;
  mutable ComplexVector action, residual, standard_rhs, standard_correction, enrichment_rhs,
      enrichment_correction, full_enrichment_correction;

  static void ExtractSubVector(const ComplexVector &source, ComplexVector &destination,
                               int offset, int size)
  {
    destination.SetSize(size);
    const auto *source_real = source.Real().HostRead() + offset;
    const auto *source_imaginary = source.Imag().HostRead() + offset;
    auto *destination_real = destination.Real().HostWrite();
    auto *destination_imaginary = destination.Imag().HostWrite();
    for (int i = 0; i < size; i++)
    {
      destination_real[i] = source_real[i];
      destination_imaginary[i] = source_imaginary[i];
    }
  }

  static void AddSubVector(const ComplexVector &source, ComplexVector &destination,
                           int offset, std::complex<double> coefficient)
  {
    MFEM_VERIFY(offset >= 0 && offset + source.Size() <= destination.Size(),
                "Singular AMS correction received an inconsistent subvector!");
    const auto *source_real = source.Real().HostRead();
    const auto *source_imaginary = source.Imag().HostRead();
    auto *destination_real = destination.Real().HostReadWrite() + offset;
    auto *destination_imaginary = destination.Imag().HostReadWrite() + offset;
    const double coefficient_real = coefficient.real();
    const double coefficient_imaginary = coefficient.imag();
    for (int i = 0; i < source.Size(); i++)
    {
      const double value_real =
          coefficient_real * source_real[i] - coefficient_imaginary * source_imaginary[i];
      const double value_imaginary =
          coefficient_imaginary * source_real[i] + coefficient_real * source_imaginary[i];
      destination_real[i] += value_real;
      destination_imaginary[i] += value_imaginary;
    }
  }

  void AddStandardCorrection(const ComplexVector &source, ComplexVector &destination,
                             std::complex<double> coefficient = 1.0) const
  {
    ExtractSubVector(source, standard_rhs, 0, standard_size);
    standard_correction.SetSize(standard_size);
    standard_solver->Mult(standard_rhs, standard_correction);
    AddSubVector(standard_correction, destination, 0, coefficient);
  }

  void AddEnrichmentCorrection(const ComplexVector &source,
                               ComplexVector &destination) const
  {
    const int enrichment_size = source.Size() - standard_size;
    ExtractSubVector(source, enrichment_rhs, standard_size, enrichment_size);
    enrichment_correction.SetSize(enrichment_size);
    enrichment_solver->Mult(enrichment_rhs, enrichment_correction);
    AddSubVector(enrichment_correction, destination, standard_size, 1.0);
  }

public:
  SingularAmsCoarseSolver(int standard_size,
                          std::unique_ptr<Solver<ComplexOperator>> &&standard_solver,
                          std::unique_ptr<Solver<ComplexOperator>> &&enrichment_solver)
    : standard_size(standard_size), standard_solver(std::move(standard_solver)),
      enrichment_solver(std::move(enrichment_solver))
  {
    MFEM_VERIFY(this->standard_size > 0 && this->standard_solver && this->enrichment_solver,
                "Singular AMS coarse solver requires two nonempty subspace solvers!");
  }

  void SetOperator(const ComplexOperator &op) override
  {
    const auto *real_matrix = dynamic_cast<const mfem::HypreParMatrix *>(op.Real());
    const auto *imaginary_matrix = dynamic_cast<const mfem::HypreParMatrix *>(op.Imag());
    MFEM_VERIFY(real_matrix && (!op.Imag() || imaginary_matrix),
                "Singular AMS coarse solver requires compatible assembled real and "
                "imaginary operator parts!");
    approximate_operator.reset();
    if (imaginary_matrix)
    {
      auto approximate_matrix = std::unique_ptr<mfem::HypreParMatrix>(
          mfem::Add(1.0, *real_matrix, 1.0, *imaginary_matrix));
      approximate_operator =
          std::make_unique<ComplexWrapperOperator>(std::move(approximate_matrix), nullptr);
      full_operator = approximate_operator.get();
    }
    else
    {
      full_operator = &op;
    }
    const auto *matrix = dynamic_cast<const mfem::HypreParMatrix *>(full_operator->Real());
    MFEM_VERIFY(matrix && !full_operator->Imag() && standard_size < matrix->Height(),
                "Singular AMS coarse approximation is inconsistent!");

    mfem::Array<int> standard_indices(standard_size);
    for (int i = 0; i < standard_indices.Size(); i++)
    {
      standard_indices[i] = i;
    }
    const int enrichment_size = matrix->Height() - standard_size;
    mfem::Array<int> enrichment_indices(enrichment_size);
    for (int i = 0; i < enrichment_indices.Size(); i++)
    {
      enrichment_indices[i] = standard_size + i;
    }

#if MFEM_HYPRE_VERSION >= 21800
    auto standard_matrix =
        std::unique_ptr<Operator>(matrix->ExtractSubmatrix(standard_indices));
    auto enrichment_matrix =
        std::unique_ptr<Operator>(matrix->ExtractSubmatrix(enrichment_indices));
#else
    std::unique_ptr<Operator> standard_matrix, enrichment_matrix;
#endif
    MFEM_VERIFY(standard_matrix && enrichment_matrix,
                "Failed to extract singular AMS coarse diagonal blocks!");
    standard_operator =
        std::make_unique<ComplexWrapperOperator>(std::move(standard_matrix), nullptr);
    enrichment_operator =
        std::make_unique<ComplexWrapperOperator>(std::move(enrichment_matrix), nullptr);
    standard_solver->SetOperator(*standard_operator);
    enrichment_solver->SetOperator(*enrichment_operator);
    this->height = op.Height();
    this->width = op.Width();
  }

  void Mult(const ComplexVector &x, ComplexVector &y) const override
  {
    MFEM_VERIFY(full_operator && x.Size() == this->width,
                "Singular AMS coarse solver is not configured for this vector!");
    y.SetSize(this->height);
    y = 0.0;

    // Symmetric multiplicative correction:
    // B_s + (I - B_s A) B_e (I - A B_s).
    AddStandardCorrection(x, y);
    action.SetSize(this->height);
    residual.SetSize(this->height);
    full_operator->Mult(y, action);
    linalg::AXPBY(1.0, x, -1.0, action);
    residual = action;

    AddEnrichmentCorrection(residual, y);
    full_enrichment_correction.SetSize(this->height);
    full_enrichment_correction.UseDevice(y.UseDevice());
    full_enrichment_correction = 0.0;
    AddSubVector(enrichment_correction, full_enrichment_correction, standard_size, 1.0);
    full_operator->Mult(full_enrichment_correction, action);
    AddStandardCorrection(action, y, -1.0);
  }
};

}  // namespace

std::unique_ptr<Solver<ComplexOperator>>
MakeSingularComplexCoarseSolver(const config::LinearSolverData &linear, LinearSolver type,
                                MatrixSymmetry matrix_symmetry, int verbose, MPI_Comm comm,
                                int standard_size,
                                FiniteElementSpaceHierarchy *primary_fespaces,
                                FiniteElementSpaceHierarchy *auxiliary_fespaces)
{
  const int print = verbose - 1;
  switch (type)
  {
    case LinearSolver::AMS:
      {
        MFEM_VERIFY(primary_fespaces && auxiliary_fespaces && standard_size > 0 &&
                        !linear.complex_coarse_solve,
                    "Combined singular AMS requires standard H(curl) and H1 spaces, a "
                    "valid standard block size, and ComplexCoarseSolve = false!");
        auto standard_solver = MakeWrapperSolver<ComplexOperator, HypreAmsSolver>(
            linear, primary_fespaces->GetFESpaceAtLevel(0),
            auxiliary_fespaces->GetFESpaceAtLevel(0), linear.ams_max_it,
            linear.mg_smooth_it, linear.ams_vector_interp, linear.ams_singular_op,
            linear.amg_agg_coarsen, print);
        auto enrichment_solver = MakeWrapperSolver<ComplexOperator, BoomerAmgSolver>(
            linear, linear.ams_max_it, linear.mg_smooth_it, linear.amg_agg_coarsen, print);
        return std::make_unique<SingularAmsCoarseSolver>(
            standard_size, std::move(standard_solver), std::move(enrichment_solver));
      }
    case LinearSolver::BOOMER_AMG:
      return MakeWrapperSolver<ComplexOperator, BoomerAmgSolver>(
          linear, linear.ams_max_it, linear.mg_smooth_it, linear.amg_agg_coarsen, print);
    case LinearSolver::SUPERLU:
#if defined(MFEM_USE_SUPERLU)
      {
        auto coarse_linear = linear;
        if (coarse_linear.sym_factorization == SymbolicFactorization::DEFAULT &&
            Mpi::Size(comm) > 1)
        {
          coarse_linear.sym_factorization = SymbolicFactorization::PARMETIS;
        }
        return MakeWrapperSolver<ComplexOperator, SuperLUSolver>(
            coarse_linear, comm, coarse_linear.sym_factorization, coarse_linear.superlu_3d,
            coarse_linear.reorder_reuse, print);
      }
#else
      MFEM_ABORT("Solver was not built with SuperLU_DIST support!");
      return {};
#endif
    case LinearSolver::STRUMPACK:
#if defined(MFEM_USE_STRUMPACK)
      return MakeWrapperSolver<ComplexOperator, StrumpackSolver>(
          linear, comm, linear.sym_factorization, linear.strumpack_compression_type,
          linear.strumpack_lr_tol, linear.strumpack_butterfly_l,
          linear.strumpack_lossy_precision, linear.reorder_reuse, print);
#else
      MFEM_ABORT("Solver was not built with STRUMPACK support!");
      return {};
#endif
    case LinearSolver::STRUMPACK_MP:
#if defined(MFEM_USE_STRUMPACK)
      return MakeWrapperSolver<ComplexOperator, StrumpackMixedPrecisionSolver>(
          linear, comm, linear.sym_factorization, linear.strumpack_compression_type,
          linear.strumpack_lr_tol, linear.strumpack_butterfly_l,
          linear.strumpack_lossy_precision, linear.reorder_reuse, print);
#else
      MFEM_ABORT("Solver was not built with STRUMPACK support!");
      return {};
#endif
    case LinearSolver::MUMPS:
#if defined(MFEM_USE_MUMPS)
      return MakeWrapperSolver<ComplexOperator, MumpsSolver>(
          linear, comm, matrix_symmetry, linear.sym_factorization,
          (linear.strumpack_compression_type == SparseCompression::BLR)
              ? linear.strumpack_lr_tol
              : 0.0,
          linear.reorder_reuse, print);
#else
      MFEM_ABORT("Solver was not built with MUMPS support!");
      return {};
#endif
    case LinearSolver::CUDSS:
#if defined(MFEM_USE_CUDSS)
      return MakeWrapperSolver<ComplexOperator, CuDSSSolver>(linear, comm, matrix_symmetry,
                                                             linear.sym_factorization,
                                                             linear.reorder_reuse, print);
#else
      MFEM_ABORT("Solver was not built with cuDSS support!");
      return {};
#endif
    case LinearSolver::JACOBI:
      return std::make_unique<JacobiSmoother<ComplexOperator>>(comm);
    case LinearSolver::DEFAULT:
      MFEM_ABORT("Unexpected unresolved coarse solver type for singular Maxwell!");
      return {};
  }
  MFEM_ABORT("Unsupported coarse solver type for singular Maxwell!");
  return {};
}

std::unique_ptr<ComplexKspSolver> MakeSingularComplexKspSolver(
    const IoData &iodata, FiniteElementSpaceHierarchy &nd_fespaces,
    FiniteElementSpaceHierarchy &h1_fespaces,
    const std::vector<const Operator *> &combined_prolongations,
    const std::vector<const Operator *> &combined_gradients,
    const std::vector<mfem::Array<int>> &combined_nd_essential_tdofs,
    const std::vector<mfem::Array<int>> &combined_h1_essential_tdofs)
{
  const auto &linear = iodata.solver.linear;
  const auto number_levels = nd_fespaces.GetNumLevels();
  MFEM_VERIFY(number_levels > 0 && h1_fespaces.GetNumLevels() == number_levels &&
                  combined_prolongations.size() + 1 == number_levels &&
                  combined_gradients.size() == number_levels &&
                  combined_nd_essential_tdofs.size() == number_levels &&
                  combined_h1_essential_tdofs.size() == number_levels,
              "Combined singular Maxwell solver received an inconsistent hierarchy!");
  const MPI_Comm comm = nd_fespaces.GetFinestFESpace().GetComm();
  MFEM_VERIFY(linear.type != LinearSolver::AMS || !linear.complex_coarse_solve,
              "Singular Maxwell with AMS does not support ComplexCoarseSolve = true!");

  auto coarse_solver = MakeSingularComplexCoarseSolver(
      linear, linear.type, GetPreconditionerMatrixSymmetry(iodata), iodata.problem.verbose,
      comm, nd_fespaces.GetFESpaceAtLevel(0).GetTrueVSize(), &nd_fespaces, &h1_fespaces);

  std::unique_ptr<Solver<ComplexOperator>> preconditioner;
  if (number_levels == 1)
  {
    preconditioner = std::move(coarse_solver);
  }
  else
  {
    std::vector<const mfem::Array<int> *> primary_essential(number_levels),
        auxiliary_essential(number_levels);
    for (std::size_t level = 0; level < number_levels; level++)
    {
      primary_essential[level] = &combined_nd_essential_tdofs[level];
      auxiliary_essential[level] = &combined_h1_essential_tdofs[level];
    }
    const auto *gradients = linear.mg_smooth_aux ? &combined_gradients : nullptr;
    const auto *gradient_essential = linear.mg_smooth_aux ? &auxiliary_essential : nullptr;
    auto multigrid = std::make_unique<GeometricMultigridSolver<ComplexOperator>>(
        comm, std::move(coarse_solver), combined_prolongations, gradients,
        linear.mg_cycle_it, linear.mg_smooth_it, linear.mg_smooth_order,
        linear.mg_smooth_sf_max, linear.mg_smooth_sf_min, linear.mg_smooth_cheby_4th,
        &primary_essential, gradient_essential);
    multigrid->EnableTimer();
    preconditioner = std::move(multigrid);
  }
  auto result = std::make_unique<ComplexKspSolver>(
      ConfigureKrylovSolver<ComplexOperator>(linear, iodata.problem.verbose, comm),
      std::move(preconditioner));
  result->EnableTimer();
  return result;
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
