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
    return std::make_unique<GeometricMultigridSolver>(iodata, std::move(pc), fespaces,
                                                      aux_fespaces);
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
  : mfem::Solver(), ksp_mult(0), ksp_mult_it(0)
{
  MFEM_VERIFY(fespaces.GetNumLevels() > 0,
              "Empty finite element space hierarchy linear solver setup!");
  MPI_Comm comm = fespaces.GetFESpaceAtLevel(0).GetComm();
  pc_ = ConfigurePreconditionerSolver(comm, iodata, fespaces, aux_fespaces);
  ksp_ = ConfigureKrylovSolver(comm, iodata);
  ksp_->SetPreconditioner(*pc_);
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
  : KspSolver()
{
  MFEM_VERIFY(fespaces.GetNumLevels() > 0,
              "Empty finite element space hierarchy linear solver setup!");
  MPI_Comm comm = fespaces.GetFESpaceAtLevel(0).GetComm();
  auto pcr = ConfigurePreconditionerSolver(comm, iodata, fespaces, aux_fespaces);
  pc_ = std::make_unique<ComplexBlockDiagonalSolver>(std::move(pcr));
  ksp_ = ConfigureKrylovSolver(comm, iodata);
  ksp_->SetPreconditioner(*pc_);
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

// XX TODO REMOVE

// KspSolver::KspSolver(MPI_Comm comm, const IoData &iodata, const std::string &prefix)
//   : clcustom(false), print(iodata.problem.verbose), print_opts(true), check_final(true),
//     solve(0)
// {
//   PalacePetscCall(KSPCreate(comm, &ksp));
//   PalacePetscCall(KSPSetOptionsPrefix(ksp, prefix.c_str()));
//   Configure(iodata);
//   ConfigureVerbose(print, prefix);
// }

// KspSolver::KspSolver(MPI_Comm comm, int print_lvl, const std::string &prefix)
//   : clcustom(false), print(print_lvl), print_opts(true), check_final(true), solve(0)
// {
//   PalacePetscCall(KSPCreate(comm, &ksp));
//   PalacePetscCall(KSPSetOptionsPrefix(ksp, prefix.c_str()));
//   ConfigureVerbose(print, prefix);
// }

// KspSolver::~KspSolver()
// {
//   MPI_Comm comm;
//   PalacePetscCall(PetscObjectGetComm(reinterpret_cast<PetscObject>(ksp), &comm));
//   PalacePetscCall(KSPDestroy(&ksp));
// }

// void KspSolver::Configure(const IoData &iodata)
// {
//   // Configure the Krylov solver. GMRES is the default solver for frequency domain
//   // problems.
//   switch (iodata.solver.linear.ksp_type)
//   {
//     case config::LinearSolverData::KspType::CG:
//       SetType(Type::CG);
//       break;
//     case config::LinearSolverData::KspType::CGSYM:
//       SetType(Type::CGSYM);
//       break;
//     case config::LinearSolverData::KspType::FCG:
//       SetType(Type::FCG);
//       break;
//     case config::LinearSolverData::KspType::MINRES:
//       SetType(Type::MINRES);
//       break;
//     case config::LinearSolverData::KspType::GMRES:
//     case config::LinearSolverData::KspType::DEFAULT:
//       SetType(Type::GMRES);
//       SetGMRESOptions(iodata.solver.linear.max_size, iodata.solver.linear.orthog_mgs,
//                       iodata.solver.linear.orthog_cgs2);
//       break;
//     case config::LinearSolverData::KspType::FGMRES:
//       SetType(Type::FGMRES);
//       SetGMRESOptions(iodata.solver.linear.max_size, iodata.solver.linear.orthog_mgs,
//                       iodata.solver.linear.orthog_cgs2);
//       break;
//     case config::LinearSolverData::KspType::BCGS:
//       SetType(Type::BCGS);
//       break;
//     case config::LinearSolverData::KspType::BCGSL:
//       SetType(Type::BCGSL);
//       break;
//     case config::LinearSolverData::KspType::FBCGS:
//       SetType(Type::FBCGS);
//       break;
//     case config::LinearSolverData::KspType::QMRCGS:
//       SetType(Type::QMRCGS);
//       break;
//     case config::LinearSolverData::KspType::TFQMR:
//       SetType(Type::TFQMR);
//       break;
//     default:
//       MFEM_ABORT("Unexpected type for KspSolver configuration!");
//       break;
//   }
//   SetTol(iodata.solver.linear.tol);
//   SetMaxIter(iodata.solver.linear.max_it);

//   // Reuse previous solution as guess for later solves if desired.
//   SetNonzeroInitialGuess(iodata.solver.linear.ksp_initial_guess);

//   // Optionally use left or right preconditioning (otherwise use PETSc default for the
//   given
//   // solver).
//   if (iodata.solver.linear.pc_side_type == config::LinearSolverData::SideType::LEFT)
//   {
//     PalacePetscCall(KSPSetPCSide(ksp, PC_LEFT));
//   }
//   else if (iodata.solver.linear.pc_side_type ==
//   config::LinearSolverData::SideType::RIGHT)
//   {
//     PalacePetscCall(KSPSetPCSide(ksp, PC_RIGHT));
//   }
// }

// void KspSolver::ConfigureVerbose(int print, const std::string &prefix)
// {
//   // Manage debugging output.
//   if (print > 0)
//   {
//     std::string opts = "-ksp_converged_reason";
//     if (print > 1)
//     {
//       opts.append(" -ksp_monitor");
//     }
//     if (print > 3)
//     {
//       opts.append(" -ksp_view");
//     }
//     if (prefix.length() > 0)
//     {
//       PetscOptionsPrefixPush(nullptr, prefix.c_str());
//     }
//     PetscOptionsInsertString(nullptr, opts.c_str());
//     if (prefix.length() > 0)
//     {
//       PetscOptionsPrefixPop(nullptr);
//     }
//   }
// }

// void KspSolver::SetType(KspSolver::Type type, bool piped)
// {
//   switch (type)
//   {
//     case Type::CG:
//       PalacePetscCall((piped) ? KSPSetType(ksp, KSPPIPECG) : KSPSetType(ksp, KSPCG));
//       PalacePetscCall(KSPCGSetType(ksp, KSP_CG_HERMITIAN));
//       break;
//     case Type::CGSYM:
//       PalacePetscCall((piped) ? KSPSetType(ksp, KSPPIPECG) : KSPSetType(ksp, KSPCG));
//       PalacePetscCall(KSPCGSetType(ksp, KSP_CG_SYMMETRIC));
//       break;
//     case Type::FCG:
//       PalacePetscCall(KSPSetType(ksp, KSPFCG));
//       break;
//     case Type::GMRES:
//       PalacePetscCall((piped) ? KSPSetType(ksp, KSPPGMRES) : KSPSetType(ksp, KSPGMRES));
//       break;
//     case Type::FGMRES:
//       PalacePetscCall((piped) ? KSPSetType(ksp, KSPPIPEFGMRES)
//                               : KSPSetType(ksp, KSPFGMRES));
//       break;
//     case Type::MINRES:
//       PalacePetscCall(KSPSetType(ksp, KSPMINRES));
//       break;
//     case Type::BCGS:
//       PalacePetscCall(KSPSetType(ksp, KSPBCGS));
//       break;
//     case Type::BCGSL:
//       PalacePetscCall(KSPSetType(ksp, KSPBCGSL));
//       PalacePetscCall(KSPBCGSLSetEll(ksp, 2));  // PETSc default
//       break;
//     case Type::FBCGS:
//       PalacePetscCall(KSPSetType(ksp, KSPFBCGS));
//       break;
//     case Type::QMRCGS:
//       PalacePetscCall(KSPSetType(ksp, KSPQMRCGS));
//       break;
//     case Type::TFQMR:
//       PalacePetscCall(KSPSetType(ksp, KSPTFQMR));
//       break;
//     case Type::CHOLESKY:
//       {
//         PC pc;
//         PalacePetscCall(KSPSetType(ksp, KSPPREONLY));
//         PalacePetscCall(KSPGetPC(ksp, &pc));
//         PalacePetscCall(PCSetType(pc, PCCHOLESKY));
//         SetCheckFinal(false);
//       }
//       break;
//     case Type::LU:
//       {
//         PC pc;
//         PalacePetscCall(KSPSetType(ksp, KSPPREONLY));
//         PalacePetscCall(KSPGetPC(ksp, &pc));
//         PalacePetscCall(PCSetType(pc, PCLU));
//         SetCheckFinal(false);
//       }
//       break;
//     default:
//       MFEM_ABORT("Unexpected type for KspSolver!");
//       break;
//   }
// }

// void KspSolver::SetTol(PetscReal tol)
// {
//   PalacePetscCall(KSPSetTolerances(ksp, tol, PETSC_DEFAULT, PETSC_DEFAULT,
//   PETSC_DEFAULT));
// }

// void KspSolver::SetAbsTol(PetscReal tol)
// {
//   PalacePetscCall(KSPSetTolerances(ksp, PETSC_DEFAULT, tol, PETSC_DEFAULT,
//   PETSC_DEFAULT));
// }

// void KspSolver::SetMaxIter(PetscInt maxits)
// {
//   PalacePetscCall(
//       KSPSetTolerances(ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, maxits));
// }

// void KspSolver::SetGMRESOptions(PetscInt maxsize, bool mgs, bool cgs2)
// {
//   PalacePetscCall(KSPGMRESSetRestart(ksp, maxsize));
//   if (mgs)
//   {
//     PalacePetscCall(
//         KSPGMRESSetOrthogonalization(ksp, KSPGMRESModifiedGramSchmidtOrthogonalization));
//   }
//   else if (cgs2)
//   {
//     PalacePetscCall(KSPGMRESSetCGSRefinementType(ksp, KSP_GMRES_CGS_REFINE_ALWAYS));
//   }
// }

// void KspSolver::SetTabLevel(PetscInt l)
// {
//   PalacePetscCall(PetscObjectSetTabLevel(reinterpret_cast<PetscObject>(ksp), l));
// }

// void KspSolver::SetNonzeroInitialGuess(bool guess)
// {
//   PalacePetscCall(KSPSetInitialGuessNonzero(ksp, guess ? PETSC_TRUE : PETSC_FALSE));
// }

// void KspSolver::SetOperator(const petsc::PetscParMatrix &A, bool copy_prefix)
// {
//   // If A is the same as before, PETSc will reuse things like symbolic factorizations
//   // automatically.
//   PalacePetscCall(KSPSetOperators(ksp, A, A));
//   if (copy_prefix)
//   {
//     // Set Mat prefix to be the same as KSP to enable setting command-line options.
//     const char *prefix;
//     PalacePetscCall(KSPGetOptionsPrefix(ksp, &prefix));
//     PalacePetscCall(MatSetOptionsPrefix(A, prefix));
//   }
// }

// void KspSolver::SetPreconditioner(const KspPreconditioner &op)
// {
//   // The PETSc shell preconditioner does not take ownership of the preconditioner object.
//   PC pc;
//   PalacePetscCall(KSPGetPC(ksp, &pc));
//   PalacePetscCall(PCSetType(pc, PCSHELL));
//   PalacePetscCall(PCShellSetContext(pc, (void *)&op));
//   PalacePetscCall(PCShellSetSetUp(pc, KspPreconditioner::PCSetUp));
//   PalacePetscCall(PCShellSetApply(pc, KspPreconditioner::PCApply));
//   PalacePetscCall(PCShellSetDestroy(pc, KspPreconditioner::PCDestroy));
// }

// void KspSolver::Customize() const
// {
//   if (!clcustom)
//   {
//     PalacePetscCall(KSPSetFromOptions(ksp));
//     if (print > 0 && print_opts)
//     {
//       PetscOptionsView(nullptr, PETSC_VIEWER_STDOUT_(GetComm()));
//       Mpi::Print(GetComm(), "\n");
//     }
//     clcustom = true;
//   }
// }

// void KspSolver::Mult(const petsc::PetscParVector &b, petsc::PetscParVector &x) const
// {
//   KSPConvergedReason reason;
//   PetscReal norm0 = 1.0, norm;
//   if (check_final)
//   {
//     norm0 = b.Norml2();
//   }
//   Customize();
//   PalacePetscCall(KSPSolve(ksp, b, x));
//   PalacePetscCall(KSPGetConvergedReason(ksp, &reason));
//   if (check_final && reason < 0)
//   {
//     Mat A;
//     Vec r;
//     PalacePetscCall(VecDuplicate(b, &r));
//     PalacePetscCall(KSPGetOperators(ksp, &A, nullptr));
//     PalacePetscCall(MatMult(A, x, r));
//     PalacePetscCall(VecAXPY(r, -1.0, b));
//     PalacePetscCall(VecNorm(r, NORM_2, &norm));
//     PalacePetscCall(VecDestroy(&r));
//     Mpi::Warning(GetComm(),
//                  "Linear solver did not converge, "
//                  "norm(Ax-b)/norm(b) = {:.3e} (norm(b) = {:.3e})!\n",
//                  norm / norm0, norm0);
//   }
//   solve++;
// }

// void KspSolver::Reset()
// {
//   PalacePetscCall(KSPReset(ksp));
// }

// PetscInt KspSolver::GetTotalNumMult() const
// {
//   return solve;
// }

// PetscInt KspSolver::GetNumIter() const
// {
//   PetscInt num_it;
//   PalacePetscCall(KSPGetIterationNumber(ksp, &num_it));
//   return num_it;
// }

// PetscInt KspSolver::GetTotalNumIter() const
// {
//   PetscInt num_it;
//   PalacePetscCall(KSPGetTotalIterations(ksp, &num_it));
//   return num_it;
// }

// MPI_Comm KspSolver::GetComm() const
// {
//   return ksp ? PetscObjectComm(reinterpret_cast<PetscObject>(ksp)) : MPI_COMM_NULL;
// }

// void KspSolver::SolveJacobi(const petsc::PetscParMatrix &A, const petsc::PetscParVector
// &b,
//                             petsc::PetscParVector &x, PetscInt sym, PetscReal tol,
//                             PetscInt max_it)
// {
//   MPI_Comm comm;
//   KSP ksp;
//   PC pc;
//   KSPConvergedReason reason;

//   comm = A.GetComm();
//   PalacePetscCall(KSPCreate(comm, &ksp));
//   PalacePetscCall(KSPSetOperators(ksp, A, A));
//   PalacePetscCall(KSPSetType(ksp, (sym == 1) ? KSPCG : KSPGMRES));
//   PalacePetscCall(KSPGetPC(ksp, &pc));
//   PalacePetscCall(PCSetType(pc, PCJACOBI));
//   PalacePetscCall(PCJacobiSetFixDiagonal(pc, PETSC_TRUE));
//   PalacePetscCall(KSPSetTolerances(ksp, tol, PETSC_DEFAULT, PETSC_DEFAULT, max_it));
//   // std::string opts = "-ksp_converged_reason -ksp_monitor";
//   // PetscOptionsInsertString(nullptr, opts.c_str());
//   // PalacePetscCall(KSPSetFromOptions(ksp));
//   x.SetZero();
//   PalacePetscCall(KSPSolve(ksp, b, x));
//   PalacePetscCall(KSPGetConvergedReason(ksp, &reason));
//   MFEM_VERIFY(reason > 0, "PETSc KSP did not converge!");
//   PalacePetscCall(KSPDestroy(&ksp));
// }

// void KspSolver::SolveDirect(const petsc::PetscParMatrix &A, const petsc::PetscParVector
// &b,
//                             petsc::PetscParVector &x, PetscInt sym)
// {
//   MPI_Comm comm;
//   KSP ksp;
//   PC pc;
//   KSPConvergedReason reason;

//   comm = A.GetComm();
//   PalacePetscCall(KSPCreate(comm, &ksp));
//   PalacePetscCall(KSPSetOperators(ksp, A, A));
//   PalacePetscCall(KSPSetType(ksp, KSPPREONLY));
//   PalacePetscCall(KSPGetPC(ksp, &pc));
// #if defined(PETSC_HAVE_MUMPS) || defined(PETSC_HAVE_SUPERLU_DIST)
//   PalacePetscCall(PCSetType(pc, (sym > 0) ? PCCHOLESKY : PCLU));
// #if defined(PETSC_HAVE_MUMPS)
//   PalacePetscCall(PCFactorSetMatSolverType(pc, MATSOLVERMUMPS));
// #elif defined(PETSC_HAVE_SUPERLU_DIST)
//   PalacePetscCall(PCFactorSetMatSolverType(pc, MATSOLVERSUPERLU_DIST));
// #endif
// #else
//   // Use PETSc default serial direct solver.
//   PalacePetscCall(PCSetType(pc, PCREDUNDANT));
//   PalacePetscCall(PCRedundantSetNumber(pc, Mpi::Size(comm)));
//   {
//     KSP ksp_in;
//     PC pc_in;
//     PalacePetscCall(PCRedundantGetKSP(pc, &ksp_in));
//     PalacePetscCall(KSPGetPC(ksp_in, &pc_in));
//     PalacePetscCall(PCSetType(pc_in, (sym > 0) ? PCCHOLESKY : PCLU));
//   }
// #endif
//   x.SetZero();
//   PalacePetscCall(KSPSolve(ksp, b, x));
//   PalacePetscCall(KSPGetConvergedReason(ksp, &reason));
//   MFEM_VERIFY(reason > 0, "PETSc KSP did not converge!");
//   PalacePetscCall(KSPDestroy(&ksp));
// }

}  // namespace palace
