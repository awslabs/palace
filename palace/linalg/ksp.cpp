// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "ksp.hpp"

#include "linalg/amg.hpp"
#include "linalg/ams.hpp"
#include "linalg/complex.hpp"
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

std::unique_ptr<mfem::IterativeSolver> ConfigureKrylovSolver(MPI_Comm comm,
                                                             const IoData &iodata)
{
  // Configure solver settings as needed based on inputs.
  config::LinearSolverData::KspType type = iodata.solver.linear.ksp_type;
  if (type == config::LinearSolverData::KspType::DEFAULT)
  {
    if (iodata.problem.type == config::ProblemData::Type::ELECTROSTATIC ||
        iodata.problem.type == config::ProblemData::Type::MAGNETOSTATIC ||
        iodata.problem.type == config::ProblemData::Type::TRANSIENT)
    {
      type = config::LinearSolverData::KspType::CG;
    }
    else
    {
      type = config::LinearSolverData::KspType::GMRES;
    }
  }
  mfem::IterativeSolver::PrintLevel print =
      mfem::IterativeSolver::PrintLevel().Warnings().Errors();
  if (iodata.problem.verbose > 0)
  {
    print.Summary();
    if (iodata.problem.verbose > 1)
    {
      print.Iterations();
      if (iodata.problem.verbose > 2)
      {
        print.All();
      }
    }
  }

  // XX TODO: We may want to replace the MFEM Krylov solvers with Hypre ones for performance
  //          (for examples, Hypre has a COGMRES solver which uses CGS (or CGS2) for
  //          orthogonalization). These will require some wrappers to allow operability with
  //          an mfem::Operator operator and mfem::Solver preconditioner.

  // Create the solver.
  std::unique_ptr<mfem::IterativeSolver> ksp;
  switch (type)
  {
    case config::LinearSolverData::KspType::CG:
      ksp = std::make_unique<mfem::CGSolver>(comm);
      break;
    case config::LinearSolverData::KspType::MINRES:
      ksp = std::make_unique<mfem::MINRESSolver>(comm);
      break;
    case config::LinearSolverData::KspType::GMRES:
      {
        auto gmres = std::make_unique<mfem::GMRESSolver>(comm);
        gmres->SetKDim(iodata.solver.linear.max_size);
        ksp = std::move(gmres);
      }
      break;
    case config::LinearSolverData::KspType::FGMRES:
      {
        auto fgmres = std::make_unique<mfem::FGMRESSolver>(comm);
        fgmres->SetKDim(iodata.solver.linear.max_size);
        ksp = std::move(fgmres);
      }
      break;
    case config::LinearSolverData::KspType::BICGSTAB:
      ksp = std::make_unique<mfem::BiCGSTABSolver>(comm);
      break;
    case config::LinearSolverData::KspType::DEFAULT:
    case config::LinearSolverData::KspType::INVALID:
      MFEM_ABORT("Unexpected solver type for Krylov solver configuration!");
      break;
  }
  ksp->iterative_mode = iodata.solver.linear.ksp_initial_guess;
  ksp->SetRelTol(iodata.solver.linear.tol);
  ksp->SetMaxIter(iodata.solver.linear.max_it);
  ksp->SetPrintLevel(print);
  return ksp;
}

std::unique_ptr<mfem::Solver>
ConfigurePreconditionerSolver(MPI_Comm comm, const IoData &iodata,
                              mfem::ParFiniteElementSpaceHierarchy &fespaces,
                              mfem::ParFiniteElementSpaceHierarchy *aux_fespaces)
{
  // Configure solver settings as needed based on inputs.
  config::LinearSolverData::Type type = iodata.solver.linear.type;
  if (type == config::LinearSolverData::Type::DEFAULT)
  {
    if (iodata.problem.type == config::ProblemData::Type::ELECTROSTATIC ||
        (iodata.problem.type == config::ProblemData::Type::TRANSIENT &&
         iodata.solver.transient.type == config::TransientSolverData::Type::CENTRAL_DIFF))
    {
      type = config::LinearSolverData::Type::BOOMER_AMG;
    }
    else if (iodata.problem.type == config::ProblemData::Type::MAGNETOSTATIC ||
             iodata.problem.type == config::ProblemData::Type::TRANSIENT)
    {
      type = config::LinearSolverData::Type::AMS;
    }
    else
    {
      // Prefer sparse direct solver for frequency domain problems if available.
#if defined(MFEM_USE_SUPERLU)
      type = config::LinearSolverData::Type::SUPERLU;
#elif defined(MFEM_USE_STRUMPACK)
      type = config::LinearSolverData::Type::STRUMPACK;
#elif defined(MFEM_USE_MUMPS)
      type = config::LinearSolverData::Type::MUMPS;
#else
      type = config::LinearSolverData::Type::AMS;
#endif
    }
  }
  int print = iodata.problem.verbose - 1;

  // Create the solver.
  std::unique_ptr<mfem::Solver> pc;
  switch (type)
  {
    case config::LinearSolverData::Type::AMS:
      // Can either be the coarse solve for geometric multigrid or the solver at the finest
      // space (in which case fespaces.GetNumLevels() == 1).
      MFEM_VERIFY(aux_fespaces, "AMS solver relies on both primary space "
                                "and auxiliary spaces for construction!");
      pc = std::make_unique<HypreAmsSolver>(iodata, fespaces.GetFESpaceAtLevel(0),
                                            aux_fespaces->GetFESpaceAtLevel(0), print);
      break;
    case config::LinearSolverData::Type::BOOMER_AMG:
      pc = std::make_unique<BoomerAmgSolver>(iodata, print);
      break;
    case config::LinearSolverData::Type::SUPERLU:
#if defined(MFEM_USE_SUPERLU)
      pc = std::make_unique<SuperLUSolver>(comm, iodata, print);
#else
      MFEM_ABORT("Solver was not built with SuperLU_DIST support, please choose a "
                 "different solver!");
#endif
      break;
    case config::LinearSolverData::Type::STRUMPACK:
#if defined(MFEM_USE_STRUMPACK)
      pc = std::make_unique<StrumpackSolver>(comm, iodata, print);
      break;
#endif
    case config::LinearSolverData::Type::STRUMPACK_MP:
#if defined(MFEM_USE_STRUMPACK) &&                                   \
    (STRUMPACK_VERSION_MAJOR >= 6 && STRUMPACK_VERSION_MINOR >= 3 && \
     STRUMPACK_VERSION_PATCH > 1)
      pc = std::make_unique<StrumpackMixedPrecisionSolver>(comm, iodata, print);
#else
      MFEM_ABORT("Solver was not built with STRUMPACK support or uses STRUMPACK older than "
                 "6.3.1 which does not include mixed-precision support, please choose a "
                 "different solver!");
#endif
      break;
    case config::LinearSolverData::Type::MUMPS:
#if defined(MFEM_USE_MUMPS)
      pc = std::make_unique<MumpsSolver>(comm, iodata, print);
#else
      MFEM_ABORT(
          "Solver was not built with MUMPS support, please choose a different solver!");
#endif
      break;
    default:
      MFEM_ABORT("Unexpected solver type for preconditioner configuration!");
      break;
  }
  if (iodata.solver.linear.mat_gmg)
  {
    // This will construct the multigrid hierarchy using pc as the coarse solver
    // (ownership of pc is transfered to the GeometricMultigridSolver). When a special
    // auxiliary space smoother for pre-/post-smoothing is not desired, the auxiliary
    // space is a nullptr here.
    return std::make_unique<GeometricMultigridSolver>(
        iodata, std::move(pc), fespaces,
        (iodata.problem.type != config::ProblemData::Type::MAGNETOSTATIC) ? aux_fespaces
                                                                          : nullptr);
  }
  else
  {
    return pc;
  }
}

class ComplexBlockDiagonalSolver : public mfem::Solver
{
private:
  std::unique_ptr<mfem::Solver> op_;

public:
  ComplexBlockDiagonalSolver(std::unique_ptr<mfem::Solver> &&op)
    : mfem::Solver(2 * op->Height(), 2 * op->Width()), op_(std::move(op))
  {
  }

  void SetOperator(const Operator &op) override {}

  void Mult(const Vector &x, Vector &y) const override
  {
    MFEM_ASSERT(x.Size() == 2 * op_->Width() && y.Size() == 2 * op_->Height(),
                "Incompatible dimensions for ComplexBlockDiagonalSolver::Mult!");
    mfem::Array<const Vector *> X(2);
    mfem::Array<Vector *> Y(2);
    Vector xr, xi, yr, yi;
    xr.MakeRef(const_cast<Vector &>(x), 0, op_->Width());
    xi.MakeRef(const_cast<Vector &>(x), op_->Width(), op_->Width());
    yr.MakeRef(y, 0, op_->Height());
    yi.MakeRef(y, op_->Height(), op_->Height());
    op_->ArrayMult(X, Y);
    yr.SyncAliasMemory(y);
    yi.SyncAliasMemory(y);
  }
};

}  // namespace

KspSolver::KspSolver(const IoData &iodata, mfem::ParFiniteElementSpaceHierarchy &fespaces,
                     mfem::ParFiniteElementSpaceHierarchy *aux_fespaces)
  : KspSolver(ConfigureKrylovSolver(fespaces.GetFinestFESpace().GetComm(), iodata),
              ConfigurePreconditionerSolver(fespaces.GetFinestFESpace().GetComm(), iodata,
                                            fespaces, aux_fespaces))
{
}

KspSolver::KspSolver(std::unique_ptr<mfem::IterativeSolver> &&ksp,
                     std::unique_ptr<mfem::Solver> &&pc)
  : mfem::Solver(), ksp_(std::move(ksp)), pc_(std::move(pc)), ksp_mult(0), ksp_mult_it(0)
{
  ksp_->SetPreconditioner(*pc_);
}

void KspSolver::SetOperator(const Operator &op, const Operator &pc_op)
{
  // Unset the preconditioner before so that IterativeSolver::SetOperator does not set the
  // preconditioner operator again.
  auto *gmg = dynamic_cast<GeometricMultigridSolver *>(pc_.get());
  if (gmg)
  {
    MFEM_ABORT("KspSolver with a GeometricMultigridSolver preconditioner must "
               "use the other signature for SetOperator!");
  }
  else
  {
    pc_->SetOperator(pc_op);
  }
  // ksp_->SetPreconditioner(nullptr);    //XX TODO WAITING MFEM PATCH
  ksp_->SetOperator(op);
  ksp_->SetPreconditioner(*pc_);
  height = op.Height();
  width = op.Width();
}

void KspSolver::SetOperator(const Operator &op,
                            const std::vector<std::unique_ptr<ParOperator>> &pc_ops,
                            const std::vector<std::unique_ptr<ParOperator>> *aux_pc_ops)
{
  // Unset the preconditioner before so that IterativeSolver::SetOperator does not set the
  // preconditioner operator again.
  auto *gmg = dynamic_cast<GeometricMultigridSolver *>(pc_.get());
  if (gmg)
  {
    gmg->SetOperator(pc_ops, aux_pc_ops);
  }
  else
  {
    MFEM_VERIFY(
        !aux_pc_ops,
        "Auxiliary space operators should not be specified for KspSolver::SetOperator "
        "unless the preconditioner is a GeometricMultigridSolver!");
    pc_->SetOperator(*pc_ops.back());
  }
  // ksp_->SetPreconditioner(nullptr);    //XX TODO WAITING MFEM PATCH
  ksp_->SetOperator(op);
  ksp_->SetPreconditioner(*pc_);
  height = op.Height();
  width = op.Width();
}

void KspSolver::Mult(const Vector &x, Vector &y) const
{
  ksp_->Mult(x, y);
  if (!ksp_->GetConverged())
  {
    Mpi::Warning(
        ksp_->GetComm(),
        "Linear solver did not converge, norm(Ax-b)/norm(b) = {:.3e} (norm(b) = {:.3e})!\n",
        ksp_->GetFinalRelNorm(), ksp_->GetInitialNorm());
  }
  ksp_mult++;
  ksp_mult_it += ksp_->GetNumIterations();
}

ComplexKspSolver::ComplexKspSolver(const IoData &iodata,
                                   mfem::ParFiniteElementSpaceHierarchy &fespaces,
                                   mfem::ParFiniteElementSpaceHierarchy *aux_fespaces)
  : KspSolver(ConfigureKrylovSolver(fespaces.GetFinestFESpace().GetComm(), iodata),
              std::make_unique<ComplexBlockDiagonalSolver>(ConfigurePreconditionerSolver(
                  fespaces.GetFinestFESpace().GetComm(), iodata, fespaces, aux_fespaces)))
{
}

ComplexKspSolver::ComplexKspSolver(std::unique_ptr<mfem::IterativeSolver> &&ksp,
                                   std::unique_ptr<mfem::Solver> &&pc)
  : KspSolver(std::move(ksp), std::make_unique<ComplexBlockDiagonalSolver>(std::move(pc)))
{
}

void ComplexKspSolver::SetOperator(const ComplexOperator &op, const Operator &pc_op)
{
  KspSolver::SetOperator(op, pc_op);  // XX TODO TEST THIS AT RUNTIME...
}

void ComplexKspSolver::SetOperator(
    const ComplexOperator &op, const std::vector<std::unique_ptr<ParOperator>> &pc_ops,
    const std::vector<std::unique_ptr<ParOperator>> *aux_pc_ops)
{
  KspSolver::SetOperator(op, pc_ops, aux_pc_ops);  // XX TODO TEST THIS AT RUNTIME...
}

void ComplexKspSolver::Mult(const ComplexVector &x, ComplexVector &y) const
{
  KspSolver::Mult(x, y);  // XX TODO TEST THIS AT RUNTIME...
  y.Sync();
}

}  // namespace palace
