// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "ksp.hpp"

#include <mfem.hpp>
#include "fem/fespace.hpp"
#include "linalg/amg.hpp"
#include "linalg/ams.hpp"
#include "linalg/gmg.hpp"
#include "linalg/mumps.hpp"
#include "linalg/strumpack.hpp"
#include "linalg/superlu.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

namespace
{

template <typename OperType>
std::unique_ptr<IterativeSolver<OperType>> ConfigureKrylovSolver(MPI_Comm comm,
                                                                 const IoData &iodata)
{
  // Create the solver.
  std::unique_ptr<IterativeSolver<OperType>> ksp;
  const auto type = iodata.solver.linear.ksp_type;
  const int print = iodata.problem.verbose;
  switch (type)
  {
    case config::LinearSolverData::KspType::CG:
      ksp = std::make_unique<CgSolver<OperType>>(comm, print);
      break;
    case config::LinearSolverData::KspType::GMRES:
      {
        auto gmres = std::make_unique<GmresSolver<OperType>>(comm, print);
        gmres->SetRestartDim(iodata.solver.linear.max_size);
        ksp = std::move(gmres);
      }
      break;
    case config::LinearSolverData::KspType::FGMRES:
      {
        auto fgmres = std::make_unique<FgmresSolver<OperType>>(comm, print);
        fgmres->SetRestartDim(iodata.solver.linear.max_size);
        ksp = std::move(fgmres);
      }
      break;
    case config::LinearSolverData::KspType::MINRES:
    case config::LinearSolverData::KspType::BICGSTAB:
    case config::LinearSolverData::KspType::DEFAULT:
      MFEM_ABORT("Unexpected solver type for Krylov solver configuration!");
      break;
  }
  ksp->SetInitialGuess(iodata.solver.linear.initial_guess);
  ksp->SetRelTol(iodata.solver.linear.tol);
  ksp->SetMaxIter(iodata.solver.linear.max_it);

  // Configure preconditioning side (only for GMRES).
  if (iodata.solver.linear.pc_side_type != config::LinearSolverData::SideType::DEFAULT &&
      type != config::LinearSolverData::KspType::GMRES)
  {
    Mpi::Warning(comm,
                 "Preconditioner side will be ignored for non-GMRES iterative solvers!\n");
  }
  else
  {
    auto *gmres = static_cast<GmresSolver<OperType> *>(ksp.get());
    switch (iodata.solver.linear.pc_side_type)
    {
      case config::LinearSolverData::SideType::LEFT:
        gmres->SetPrecSide(GmresSolver<OperType>::PrecSide::LEFT);
        break;
      case config::LinearSolverData::SideType::RIGHT:
        gmres->SetPrecSide(GmresSolver<OperType>::PrecSide::RIGHT);
        break;
      case config::LinearSolverData::SideType::DEFAULT:
        // Do nothing
        break;
    }
  }

  // Configure orthogonalization method for GMRES/FMGRES.
  if (type == config::LinearSolverData::KspType::GMRES ||
      type == config::LinearSolverData::KspType::FGMRES)
  {
    // Because FGMRES inherits from GMRES, this is OK.
    auto *gmres = static_cast<GmresSolver<OperType> *>(ksp.get());
    switch (iodata.solver.linear.gs_orthog_type)
    {
      case config::LinearSolverData::OrthogType::MGS:
        gmres->SetOrthogonalization(GmresSolver<OperType>::OrthogType::MGS);
        break;
      case config::LinearSolverData::OrthogType::CGS:
        gmres->SetOrthogonalization(GmresSolver<OperType>::OrthogType::CGS);
        break;
      case config::LinearSolverData::OrthogType::CGS2:
        gmres->SetOrthogonalization(GmresSolver<OperType>::OrthogType::CGS2);
        break;
    }
  }

  return ksp;
}

template <typename OperType>
std::unique_ptr<Solver<OperType>>
ConfigurePreconditionerSolver(MPI_Comm comm, const IoData &iodata,
                              const FiniteElementSpaceHierarchy &fespaces,
                              const AuxiliaryFiniteElementSpaceHierarchy *aux_fespaces)
{
  // Create the real-valued solver first.
  std::unique_ptr<mfem::Solver> pc0;
  const auto type = iodata.solver.linear.type;
  const int print = iodata.problem.verbose - 1;
  switch (type)
  {
    case config::LinearSolverData::Type::AMS:
      // Can either be the coarse solve for geometric multigrid or the solver at the finest
      // space (in which case fespaces.GetNumLevels() == 1).
      MFEM_VERIFY(aux_fespaces, "AMS solver relies on both primary space "
                                "and auxiliary spaces for construction!");
      pc0 = std::make_unique<HypreAmsSolver>(iodata, fespaces.GetNumLevels() > 1,
                                             fespaces.GetFESpaceAtLevel(0),
                                             aux_fespaces->GetFESpaceAtLevel(0), print);
      break;
    case config::LinearSolverData::Type::BOOMER_AMG:
      pc0 = std::make_unique<BoomerAmgSolver>(iodata, fespaces.GetNumLevels() > 1, print);
      break;
    case config::LinearSolverData::Type::SUPERLU:
#if defined(MFEM_USE_SUPERLU)
      pc0 = std::make_unique<SuperLUSolver>(comm, iodata, print);
#else
      MFEM_ABORT("Solver was not built with SuperLU_DIST support, please choose a "
                 "different solver!");
#endif
      break;
    case config::LinearSolverData::Type::STRUMPACK:
#if defined(MFEM_USE_STRUMPACK)
      pc0 = std::make_unique<StrumpackSolver>(comm, iodata, print);
#else
      MFEM_ABORT("Solver was not built with STRUMPACK support, please choose a "
                 "different solver!");
#endif
      break;
    case config::LinearSolverData::Type::STRUMPACK_MP:
#if defined(MFEM_USE_STRUMPACK)
      pc0 = std::make_unique<StrumpackMixedPrecisionSolver>(comm, iodata, print);
#else
      MFEM_ABORT("Solver was not built with STRUMPACK support, please choose a "
                 "different solver!");
#endif
      break;
    case config::LinearSolverData::Type::MUMPS:
#if defined(MFEM_USE_MUMPS)
      pc0 = std::make_unique<MumpsSolver>(comm, iodata, print);
#else
      MFEM_ABORT(
          "Solver was not built with MUMPS support, please choose a different solver!");
#endif
      break;
    case config::LinearSolverData::Type::DEFAULT:
      MFEM_ABORT("Unexpected solver type for preconditioner configuration!");
      break;
  }

  // Construct the actual solver, which has the right value type.
  auto pc = std::make_unique<WrapperSolver<OperType>>(std::move(pc0));
  if (fespaces.GetNumLevels() > 1)
  {
    // This will construct the multigrid hierarchy using pc as the coarse solver
    // (ownership of pc is transferred to the GeometricMultigridSolver). When a special
    // auxiliary space smoother for pre-/post-smoothing is not desired, the auxiliary
    // space is a nullptr here.
    if (iodata.solver.linear.mg_smooth_aux)
    {
      MFEM_VERIFY(aux_fespaces, "Multigrid with auxiliary space smoothers requires both "
                                "primary space and auxiliary spaces for construction!");
      const auto G = aux_fespaces->GetDiscreteInterpolators();
      return std::make_unique<GeometricMultigridSolver<OperType>>(
          iodata, std::move(pc), fespaces.GetProlongationOperators(), &G);
    }
    else
    {
      return std::make_unique<GeometricMultigridSolver<OperType>>(
          iodata, std::move(pc), fespaces.GetProlongationOperators());
    }
  }
  else
  {
    return pc;
  }
}

}  // namespace

template <typename OperType>
BaseKspSolver<OperType>::BaseKspSolver(
    const IoData &iodata, const FiniteElementSpaceHierarchy &fespaces,
    const AuxiliaryFiniteElementSpaceHierarchy *aux_fespaces)
  : BaseKspSolver(
        ConfigureKrylovSolver<OperType>(fespaces.GetFinestFESpace().GetComm(), iodata),
        ConfigurePreconditionerSolver<OperType>(fespaces.GetFinestFESpace().GetComm(),
                                                iodata, fespaces, aux_fespaces))
{
}

template <typename OperType>
BaseKspSolver<OperType>::BaseKspSolver(std::unique_ptr<IterativeSolver<OperType>> &&ksp,
                                       std::unique_ptr<Solver<OperType>> &&pc)
  : ksp(std::move(ksp)), pc(std::move(pc)), ksp_mult(0), ksp_mult_it(0)
{
  this->ksp->SetPreconditioner(*this->pc);
}

template <typename OperType>
void BaseKspSolver<OperType>::SetOperators(const OperType &op, const OperType &pc_op)
{
  ksp->SetOperator(op);
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

template <typename OperType>
void BaseKspSolver<OperType>::Mult(const VecType &x, VecType &y) const
{
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
