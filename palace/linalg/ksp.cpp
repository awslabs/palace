// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "ksp.hpp"

#include <petsc.h>
#include <mfem.hpp>
#include "linalg/pc.hpp"
#include "linalg/petsc.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

KspSolver::KspSolver(MPI_Comm comm, const IoData &iodata, const std::string &prefix)
  : clcustom(false), print(iodata.problem.verbose), print_opts(true), check_final(true),
    solve(0)
{
  PalacePetscCall(KSPCreate(comm, &ksp));
  PalacePetscCall(KSPSetOptionsPrefix(ksp, prefix.c_str()));
  Configure(iodata);
  ConfigureVerbose(print, prefix);
}

KspSolver::KspSolver(MPI_Comm comm, int print_lvl, const std::string &prefix)
  : clcustom(false), print(print_lvl), print_opts(true), check_final(true), solve(0)
{
  PalacePetscCall(KSPCreate(comm, &ksp));
  PalacePetscCall(KSPSetOptionsPrefix(ksp, prefix.c_str()));
  ConfigureVerbose(print, prefix);
}

KspSolver::~KspSolver()
{
  MPI_Comm comm;
  PalacePetscCall(PetscObjectGetComm(reinterpret_cast<PetscObject>(ksp), &comm));
  PalacePetscCall(KSPDestroy(&ksp));
}

void KspSolver::Configure(const IoData &iodata)
{
  // Configure the Krylov solver. GMRES is the default solver for frequency domain
  // problems.
  switch (iodata.solver.linear.ksp_type)
  {
    case config::LinearSolverData::KspType::CG:
      SetType(Type::CG);
      break;
    case config::LinearSolverData::KspType::CGSYM:
      SetType(Type::CGSYM);
      break;
    case config::LinearSolverData::KspType::FCG:
      SetType(Type::FCG);
      break;
    case config::LinearSolverData::KspType::MINRES:
      SetType(Type::MINRES);
      break;
    case config::LinearSolverData::KspType::GMRES:
    case config::LinearSolverData::KspType::DEFAULT:
      SetType(Type::GMRES);
      SetGMRESOptions(iodata.solver.linear.max_size, iodata.solver.linear.orthog_mgs,
                      iodata.solver.linear.orthog_cgs2);
      break;
    case config::LinearSolverData::KspType::FGMRES:
      SetType(Type::FGMRES);
      SetGMRESOptions(iodata.solver.linear.max_size, iodata.solver.linear.orthog_mgs,
                      iodata.solver.linear.orthog_cgs2);
      break;
    case config::LinearSolverData::KspType::BCGS:
      SetType(Type::BCGS);
      break;
    case config::LinearSolverData::KspType::BCGSL:
      SetType(Type::BCGSL);
      break;
    case config::LinearSolverData::KspType::FBCGS:
      SetType(Type::FBCGS);
      break;
    case config::LinearSolverData::KspType::QMRCGS:
      SetType(Type::QMRCGS);
      break;
    case config::LinearSolverData::KspType::TFQMR:
      SetType(Type::TFQMR);
      break;
    case config::LinearSolverData::KspType::INVALID:
      MFEM_ABORT("Unexpected type for KspSolver configuration!");
      break;
  }
  SetTol(iodata.solver.linear.tol);
  SetMaxIter(iodata.solver.linear.max_it);

  // Reuse previous solution as guess for later solves if desired.
  SetNonzeroInitialGuess(iodata.solver.linear.ksp_initial_guess);

  // Optionally use left or right preconditioning (otherwise use PETSc default for the given
  // solver).
  if (iodata.solver.linear.pc_side_type == config::LinearSolverData::SideType::LEFT)
  {
    PalacePetscCall(KSPSetPCSide(ksp, PC_LEFT));
  }
  else if (iodata.solver.linear.pc_side_type == config::LinearSolverData::SideType::RIGHT)
  {
    PalacePetscCall(KSPSetPCSide(ksp, PC_RIGHT));
  }
}

void KspSolver::ConfigureVerbose(int print, const std::string &prefix)
{
  // Manage debugging output.
  if (print > 0)
  {
    std::string opts = "-ksp_converged_reason";
    if (print > 1)
    {
      opts.append(" -ksp_monitor");
    }
    if (print > 3)
    {
      opts.append(" -ksp_view");
    }
    if (prefix.length() > 0)
    {
      PetscOptionsPrefixPush(nullptr, prefix.c_str());
    }
    PetscOptionsInsertString(nullptr, opts.c_str());
    if (prefix.length() > 0)
    {
      PetscOptionsPrefixPop(nullptr);
    }
  }
}

void KspSolver::SetType(KspSolver::Type type, bool piped)
{
  switch (type)
  {
    case Type::CG:
      PalacePetscCall((piped) ? KSPSetType(ksp, KSPPIPECG) : KSPSetType(ksp, KSPCG));
      PalacePetscCall(KSPCGSetType(ksp, KSP_CG_HERMITIAN));
      break;
    case Type::CGSYM:
      PalacePetscCall((piped) ? KSPSetType(ksp, KSPPIPECG) : KSPSetType(ksp, KSPCG));
      PalacePetscCall(KSPCGSetType(ksp, KSP_CG_SYMMETRIC));
      break;
    case Type::FCG:
      PalacePetscCall(KSPSetType(ksp, KSPFCG));
      break;
    case Type::GMRES:
      PalacePetscCall((piped) ? KSPSetType(ksp, KSPPGMRES) : KSPSetType(ksp, KSPGMRES));
      break;
    case Type::FGMRES:
      PalacePetscCall((piped) ? KSPSetType(ksp, KSPPIPEFGMRES)
                              : KSPSetType(ksp, KSPFGMRES));
      break;
    case Type::MINRES:
      PalacePetscCall(KSPSetType(ksp, KSPMINRES));
      break;
    case Type::BCGS:
      PalacePetscCall(KSPSetType(ksp, KSPBCGS));
      break;
    case Type::BCGSL:
      PalacePetscCall(KSPSetType(ksp, KSPBCGSL));
      PalacePetscCall(KSPBCGSLSetEll(ksp, 2));  // PETSc default
      break;
    case Type::FBCGS:
      PalacePetscCall(KSPSetType(ksp, KSPFBCGS));
      break;
    case Type::QMRCGS:
      PalacePetscCall(KSPSetType(ksp, KSPQMRCGS));
      break;
    case Type::TFQMR:
      PalacePetscCall(KSPSetType(ksp, KSPTFQMR));
      break;
    case Type::CHOLESKY:
      {
        PC pc;
        PalacePetscCall(KSPSetType(ksp, KSPPREONLY));
        PalacePetscCall(KSPGetPC(ksp, &pc));
        PalacePetscCall(PCSetType(pc, PCCHOLESKY));
        SetCheckFinal(false);
      }
      break;
    case Type::LU:
      {
        PC pc;
        PalacePetscCall(KSPSetType(ksp, KSPPREONLY));
        PalacePetscCall(KSPGetPC(ksp, &pc));
        PalacePetscCall(PCSetType(pc, PCLU));
        SetCheckFinal(false);
      }
      break;
  }
}

void KspSolver::SetTol(PetscReal tol)
{
  PalacePetscCall(KSPSetTolerances(ksp, tol, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT));
}

void KspSolver::SetAbsTol(PetscReal tol)
{
  PalacePetscCall(KSPSetTolerances(ksp, PETSC_DEFAULT, tol, PETSC_DEFAULT, PETSC_DEFAULT));
}

void KspSolver::SetMaxIter(PetscInt maxits)
{
  PalacePetscCall(
      KSPSetTolerances(ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, maxits));
}

void KspSolver::SetGMRESOptions(PetscInt maxsize, bool mgs, bool cgs2)
{
  PalacePetscCall(KSPGMRESSetRestart(ksp, maxsize));
  if (mgs)
  {
    PalacePetscCall(
        KSPGMRESSetOrthogonalization(ksp, KSPGMRESModifiedGramSchmidtOrthogonalization));
  }
  else if (cgs2)
  {
    PalacePetscCall(KSPGMRESSetCGSRefinementType(ksp, KSP_GMRES_CGS_REFINE_ALWAYS));
  }
}

void KspSolver::SetTabLevel(PetscInt l)
{
  PalacePetscCall(PetscObjectSetTabLevel(reinterpret_cast<PetscObject>(ksp), l));
}

void KspSolver::SetNonzeroInitialGuess(bool guess)
{
  PalacePetscCall(KSPSetInitialGuessNonzero(ksp, guess ? PETSC_TRUE : PETSC_FALSE));
}

void KspSolver::SetOperator(const petsc::PetscParMatrix &A, bool copy_prefix)
{
  // If A is the same as before, PETSc will reuse things like symbolic factorizations
  // automatically.
  PalacePetscCall(KSPSetOperators(ksp, A, A));
  if (copy_prefix)
  {
    // Set Mat prefix to be the same as KSP to enable setting command-line options.
    const char *prefix;
    PalacePetscCall(KSPGetOptionsPrefix(ksp, &prefix));
    PalacePetscCall(MatSetOptionsPrefix(A, prefix));
  }
}

void KspSolver::SetPreconditioner(const KspPreconditioner &op)
{
  // The PETSc shell preconditioner does not take ownership of the preconditioner object.
  PC pc;
  PalacePetscCall(KSPGetPC(ksp, &pc));
  PalacePetscCall(PCSetType(pc, PCSHELL));
  PalacePetscCall(PCShellSetContext(pc, (void *)&op));
  PalacePetscCall(PCShellSetSetUp(pc, KspPreconditioner::PCSetUp));
  PalacePetscCall(PCShellSetApply(pc, KspPreconditioner::PCApply));
  PalacePetscCall(PCShellSetDestroy(pc, KspPreconditioner::PCDestroy));
}

void KspSolver::Customize() const
{
  if (!clcustom)
  {
    PalacePetscCall(KSPSetFromOptions(ksp));
    if (print > 0 && print_opts)
    {
      PetscOptionsView(nullptr, PETSC_VIEWER_STDOUT_(GetComm()));
      Mpi::Print(GetComm(), "\n");
    }
    clcustom = true;
  }
}

void KspSolver::Mult(const petsc::PetscParVector &b, petsc::PetscParVector &x) const
{
  KSPConvergedReason reason;
  PetscReal norm0 = 1.0, norm;
  if (check_final)
  {
    norm0 = b.Norml2();
  }
  Customize();
  PalacePetscCall(KSPSolve(ksp, b, x));
  PalacePetscCall(KSPGetConvergedReason(ksp, &reason));
  if (check_final && reason < 0)
  {
    Mat A;
    Vec r;
    PalacePetscCall(VecDuplicate(b, &r));
    PalacePetscCall(KSPGetOperators(ksp, &A, nullptr));
    PalacePetscCall(MatMult(A, x, r));
    PalacePetscCall(VecAXPY(r, -1.0, b));
    PalacePetscCall(VecNorm(r, NORM_2, &norm));
    PalacePetscCall(VecDestroy(&r));
    Mpi::Warning(GetComm(),
                 "Linear solver did not converge, "
                 "norm(Ax-b)/norm(b) = {:.3e} (norm(b) = {:.3e})!\n",
                 norm / norm0, norm0);
  }
  solve++;
}

void KspSolver::Reset()
{
  PalacePetscCall(KSPReset(ksp));
}

PetscInt KspSolver::GetTotalNumMult() const
{
  return solve;
}

PetscInt KspSolver::GetNumIter() const
{
  PetscInt num_it;
  PalacePetscCall(KSPGetIterationNumber(ksp, &num_it));
  return num_it;
}

PetscInt KspSolver::GetTotalNumIter() const
{
  PetscInt num_it;
  PalacePetscCall(KSPGetTotalIterations(ksp, &num_it));
  return num_it;
}

MPI_Comm KspSolver::GetComm() const
{
  return ksp ? PetscObjectComm(reinterpret_cast<PetscObject>(ksp)) : MPI_COMM_NULL;
}

void KspSolver::SolveJacobi(const petsc::PetscParMatrix &A, const petsc::PetscParVector &b,
                            petsc::PetscParVector &x, PetscInt sym, PetscReal tol,
                            PetscInt max_it)
{
  MPI_Comm comm;
  KSP ksp;
  PC pc;
  KSPConvergedReason reason;

  comm = A.GetComm();
  PalacePetscCall(KSPCreate(comm, &ksp));
  PalacePetscCall(KSPSetOperators(ksp, A, A));
  PalacePetscCall(KSPSetType(ksp, (sym == 1) ? KSPCG : KSPGMRES));
  PalacePetscCall(KSPGetPC(ksp, &pc));
  PalacePetscCall(PCSetType(pc, PCJACOBI));
  PalacePetscCall(PCJacobiSetFixDiagonal(pc, PETSC_TRUE));
  PalacePetscCall(KSPSetTolerances(ksp, tol, PETSC_DEFAULT, PETSC_DEFAULT, max_it));
  // std::string opts = "-ksp_converged_reason -ksp_monitor";
  // PetscOptionsInsertString(nullptr, opts.c_str());
  // PalacePetscCall(KSPSetFromOptions(ksp));
  x.SetZero();
  PalacePetscCall(KSPSolve(ksp, b, x));
  PalacePetscCall(KSPGetConvergedReason(ksp, &reason));
  MFEM_VERIFY(reason > 0, "PETSc KSP did not converge!");
  PalacePetscCall(KSPDestroy(&ksp));
}

void KspSolver::SolveDirect(const petsc::PetscParMatrix &A, const petsc::PetscParVector &b,
                            petsc::PetscParVector &x, PetscInt sym)
{
  MPI_Comm comm;
  KSP ksp;
  PC pc;
  KSPConvergedReason reason;

  comm = A.GetComm();
  PalacePetscCall(KSPCreate(comm, &ksp));
  PalacePetscCall(KSPSetOperators(ksp, A, A));
  PalacePetscCall(KSPSetType(ksp, KSPPREONLY));
  PalacePetscCall(KSPGetPC(ksp, &pc));
#if defined(PETSC_HAVE_MUMPS) || defined(PETSC_HAVE_SUPERLU_DIST)
  PalacePetscCall(PCSetType(pc, (sym > 0) ? PCCHOLESKY : PCLU));
#if defined(PETSC_HAVE_MUMPS)
  PalacePetscCall(PCFactorSetMatSolverType(pc, MATSOLVERMUMPS));
#elif defined(PETSC_HAVE_SUPERLU_DIST)
  PalacePetscCall(PCFactorSetMatSolverType(pc, MATSOLVERSUPERLU_DIST));
#endif
#else
  // Use PETSc default serial direct solver.
  PalacePetscCall(PCSetType(pc, PCREDUNDANT));
  PalacePetscCall(PCRedundantSetNumber(pc, Mpi::Size(comm)));
  {
    KSP ksp_in;
    PC pc_in;
    PalacePetscCall(PCRedundantGetKSP(pc, &ksp_in));
    PalacePetscCall(KSPGetPC(ksp_in, &pc_in));
    PalacePetscCall(PCSetType(pc_in, (sym > 0) ? PCCHOLESKY : PCLU));
  }
#endif
  x.SetZero();
  PalacePetscCall(KSPSolve(ksp, b, x));
  PalacePetscCall(KSPGetConvergedReason(ksp, &reason));
  MFEM_VERIFY(reason > 0, "PETSc KSP did not converge!");
  PalacePetscCall(KSPDestroy(&ksp));
}

}  // namespace palace
