// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "slepc.hpp"

#if defined(PALACE_WITH_SLEPC)

#include <petsc.h>
#include <slepc.h>
#include <mfem.hpp>
#include "linalg/divfree.hpp"
#include "linalg/ksp.hpp"
#include "utils/communication.hpp"

static PetscErrorCode __mat_apply_EPS_A(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_EPS_B(Mat, Vec, Vec);
static PetscErrorCode __pc_apply_EPS(PC, Vec, Vec);
static PetscErrorCode __mat_apply_PEPLinear_L0(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_PEPLinear_L1(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_PEPLinear_B(Mat, Vec, Vec);
static PetscErrorCode __pc_apply_PEPLinear(PC, Vec, Vec);
static PetscErrorCode __mat_apply_PEP_A0(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_PEP_A1(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_PEP_A2(Mat, Vec, Vec);
static PetscErrorCode __pc_apply_PEP(PC, Vec, Vec);

namespace palace::slepc
{

void Initialize(int &argc, char **&argv, const char rc_file[], const char help[])
{
  PalacePetscCall(SlepcInitialize(&argc, &argv, rc_file, help));
}

void Initialize()
{
  PalacePetscCall(SlepcInitializeNoArguments());
}

void Finalize()
{
  PalacePetscCall(SlepcFinalize());
}

PetscReal GetMaxSingularValue(const petsc::PetscParMatrix &A, PetscReal tol,
                              PetscInt maxits)
{
  // This method assumes the provided operator has the required operations for SLEPc's EPS
  // or SVD solvers, namely MATOP_MULT and MATOP_MULT_HERMITIAN_TRANSPOSE (if the matrix
  // is not Hermitian).
  PetscInt nconv;
  PetscReal sigma;
  if (A.GetHermitian())  // Returns true if symmetric and not PETSC_USE_COMPLEX
  {
    EPS eps;
    PetscScalar eig;
    PalacePetscCall(EPSCreate(A.GetComm(), &eps));
    PalacePetscCall(EPSSetOperators(eps, A, nullptr));
    PalacePetscCall(EPSSetProblemType(eps, EPS_HEP));
    PalacePetscCall(EPSSetWhichEigenpairs(eps, EPS_LARGEST_MAGNITUDE));
    PalacePetscCall(EPSSetDimensions(eps, 1, PETSC_DEFAULT, PETSC_DEFAULT));
    PalacePetscCall(EPSSetTolerances(eps, tol, maxits));
    PalacePetscCall(EPSSolve(eps));
    PalacePetscCall(EPSGetConverged(eps, &nconv));
    if (nconv < 1)
    {
      Mpi::Warning(A.GetComm(),
                   "SLEPc EPS solve did not converge for maximum singular value!\n");
      PalacePetscCall(EPSDestroy(&eps));
      return 0.0;
    }
    MFEM_VERIFY(nconv >= 1, " ");
    PalacePetscCall(EPSGetEigenvalue(eps, 0, &eig, nullptr));
    PalacePetscCall(EPSDestroy(&eps));
    MFEM_VERIFY(PetscImaginaryPart(eig) == 0.0,
                "Unexpected complex eigenvalue for Hermitian matrix (λ = " << eig << ")!");
    sigma = PetscAbsScalar(eig);
  }
  else
  {
    SVD svd;
    PalacePetscCall(SVDCreate(A.GetComm(), &svd));
    PalacePetscCall(SVDSetOperators(svd, A, nullptr));
    PalacePetscCall(SVDSetProblemType(svd, SVD_STANDARD));
    PalacePetscCall(SVDSetWhichSingularTriplets(svd, SVD_LARGEST));
    PalacePetscCall(SVDSetDimensions(svd, 1, PETSC_DEFAULT, PETSC_DEFAULT));
    PalacePetscCall(SVDSetTolerances(svd, tol, maxits));
    PalacePetscCall(SVDSolve(svd));
    PalacePetscCall(SVDGetConverged(svd, &nconv));
    if (nconv < 1)
    {
      Mpi::Warning(A.GetComm(),
                   "SLEPc SVD solve did not converge for maximum singular value!\n");
      PalacePetscCall(SVDDestroy(&svd));
      return 0.0;
    }
    MFEM_VERIFY(nconv >= 1, " ");
    PalacePetscCall(SVDGetSingularTriplet(svd, 0, &sigma, nullptr, nullptr));
    PalacePetscCall(SVDDestroy(&svd));
  }
  return sigma;
}

// Eigensolver base class methods

SlepcEigenSolver::SlepcEigenSolver(int print_lvl) : clcustom(false), print(print_lvl)
{
  sinvert = false;
  region = true;
  sigma = 0.0;
  gamma = delta = 1.0;

  res = nullptr;
  v0 = r0 = nullptr;
  opInv = nullptr;
  opProj = nullptr;
}

SlepcEigenSolver::~SlepcEigenSolver()
{
  delete[] res;
  delete v0;
  delete r0;
}

void SlepcEigenSolver::SetOperators(const petsc::PetscParMatrix &K,
                                    const petsc::PetscParMatrix &M,
                                    EigenSolverBase::ScaleType type)
{
  MFEM_ABORT("SetOperators not defined for base class SlepcEigenSolver!");
}

void SlepcEigenSolver::SetOperators(const petsc::PetscParMatrix &K,
                                    const petsc::PetscParMatrix &C,
                                    const petsc::PetscParMatrix &M,
                                    EigenSolverBase::ScaleType type)
{
  MFEM_ABORT("SetOperators not defined for base class SlepcEigenSolver!");
}

void SlepcEigenSolver::SetLinearSolver(const KspSolver &ksp)
{
  opInv = &ksp;
}

void SlepcEigenSolver::SetProjector(const DivFreeSolver &divfree)
{
  opProj = &divfree;
}

void SlepcEigenSolver::SetBMat(const petsc::PetscParMatrix &B)
{
  BV bv = GetBV();
  PalacePetscCall(BVSetMatrix(bv, B, PETSC_FALSE));
}

void SlepcEigenSolver::SetShiftInvert(double tr, double ti, bool precond)
{
  ST st = GetST();
  if (precond)
  {
    PalacePetscCall(STSetType(st, STPRECOND));
  }
  else
  {
    PalacePetscCall(STSetType(st, STSINVERT));
  }
  PalacePetscCall(STSetTransform(st, PETSC_TRUE));
  PalacePetscCall(STSetMatMode(st, ST_MATMODE_SHELL));
  sigma = tr + PETSC_i * ti;  // Wait until solve time to call EPS/PEPSetTarget
  sinvert = true;
}

void SlepcEigenSolver::SetOrthogonalization(bool mgs, bool cgs2)
{
  if (mgs || cgs2)
  {
    BVOrthogType type;
    BVOrthogRefineType refine;
    PetscReal eta;
    BVOrthogBlockType btype;
    BV bv = GetBV();
    if (mgs)
    {
      type = BV_ORTHOG_MGS;
      PalacePetscCall(BVGetOrthogonalization(bv, nullptr, &refine, &eta, &btype));
    }
    else  // cgs2
    {
      type = BV_ORTHOG_CGS;
      refine = BV_ORTHOG_REFINE_ALWAYS;
      eta = 1.0;
      PalacePetscCall(BVGetOrthogonalization(bv, nullptr, nullptr, nullptr, &btype));
    }
    PalacePetscCall(BVSetOrthogonalization(bv, type, refine, eta, btype));
  }
}

void SlepcEigenSolver::Customize()
{
  // Configure the KSP object for non-preconditioned spectral transformations.
  PetscBool precond;
  KSP ksp;
  ST st = GetST();
  PalacePetscCall(STGetKSP(st, &ksp));
  PalacePetscCall(
      PetscObjectTypeCompare(reinterpret_cast<PetscObject>(st), STPRECOND, &precond));
  if (!precond)
  {
    PalacePetscCall(KSPSetType(ksp, KSPPREONLY));
  }

  // Configure the region based on the given target if necessary.
  if (sinvert && region)
  {
    if (PetscImaginaryPart(sigma) == 0.0)
    {
      if (PetscRealPart(sigma) > 0.0)
      {
        SetRegion(PetscRealPart(sigma) / gamma, mfem::infinity(), -mfem::infinity(),
                  mfem::infinity());
      }
      else if (PetscRealPart(sigma) < 0.0)
      {
        SetRegion(-mfem::infinity(), PetscRealPart(sigma) / gamma, -mfem::infinity(),
                  mfem::infinity());
      }
    }
    else if (PetscRealPart(sigma) == 0.0)
    {
      if (PetscImaginaryPart(sigma) > 0.0)
      {
        SetRegion(-mfem::infinity(), mfem::infinity(), PetscImaginaryPart(sigma) / gamma,
                  mfem::infinity());
      }
      else if (PetscImaginaryPart(sigma) < 0.0)
      {
        SetRegion(-mfem::infinity(), mfem::infinity(), -mfem::infinity(),
                  PetscImaginaryPart(sigma) / gamma);
      }
    }
    else
    {
      MFEM_ABORT("Shift-and-invert with general complex eigenvalue target "
                 "is unsupported!");
    }
  }
}

void SlepcEigenSolver::SetPCShell(void *ctx, PetscErrorCode (*__pc_apply)(PC, Vec, Vec))
{
  // Configure linear solver for generalized problem or spectral transformation. This also
  // allows use of the divergence-free projector as a linear solve side-effect.
  KSP ksp;
  PC pc;
  ST st = GetST();
  PalacePetscCall(STGetKSP(st, &ksp));
  PalacePetscCall(KSPGetPC(ksp, &pc));

  // Configure the linear solver as a shell preconditioner.
  PalacePetscCall(PCSetType(pc, PCSHELL));
  PalacePetscCall(PCShellSetContext(pc, ctx));
  PalacePetscCall(PCShellSetApply(pc, __pc_apply));
}

void SlepcEigenSolver::SetRegion(PetscReal rminr, PetscReal rmaxr, PetscReal rmini,
                                 PetscReal rmaxi, bool complement)
{
  RG rg = GetRG();
  PalacePetscCall(RGSetType(rg, RGINTERVAL));
  PalacePetscCall(RGIntervalSetEndpoints(rg, rminr, rmaxr, rmini, rmaxi));
  if (complement)
  {
    PalacePetscCall(RGSetComplement(rg, PETSC_TRUE));
  }
}

void SlepcEigenSolver::GetBackTransform(PetscScalar eig, PetscReal &eigr,
                                        PetscReal &eigi) const
{
  eigr = gamma * PetscRealPart(eig);
  eigi = gamma * PetscImaginaryPart(eig);
}

void SlepcEigenSolver::GetError(int i, EigenSolverBase::ErrorType type, double &err) const
{
  PetscReal eigr, eigi;
  GetEigenvalue(i, eigr, eigi);
  PetscScalar eig = eigr + PETSC_i * eigi;
  if (res[i] <= 0.0)
  {
    GetEigenvector(i, *v0);
    GetResidual(eig, *v0, *r0);
    res[i] = r0->Norml2() / v0->Norml2();
  }
  switch (type)
  {
    case ErrorType::ABSOLUTE:
      err = res[i];
      break;
    case ErrorType::RELATIVE:
      err = res[i] / PetscAbsScalar(eig);
      break;
    case ErrorType::BACKWARD:
      err = res[i] / GetBackwardScaling(eig);
      break;
    default:
      MFEM_ABORT("Eigenpair error type not implemented!");
      break;
  }
}

// EPS specific methods

SlepcEPSSolverBase::SlepcEPSSolverBase(MPI_Comm comm, int print_lvl,
                                       const std::string &prefix)
  : SlepcEigenSolver(print_lvl)
{
  PalacePetscCall(EPSCreate(comm, &eps));
  PalacePetscCall(EPSSetOptionsPrefix(eps, prefix.c_str()));
  if (print > 0)
  {
    std::string opts = "-eps_monitor";
    if (print > 2)
    {
      opts.append(" -eps_view");
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
  A = B = nullptr;
}

SlepcEPSSolverBase::~SlepcEPSSolverBase()
{
  MPI_Comm comm;
  PalacePetscCall(PetscObjectGetComm(reinterpret_cast<PetscObject>(eps), &comm));
  PalacePetscCall(EPSDestroy(&eps));
  delete A;
  delete B;
}

void SlepcEPSSolverBase::SetNumModes(int numeig, int numvec)
{
  PalacePetscCall(
      EPSSetDimensions(eps, numeig, (numvec > 0) ? numvec : PETSC_DEFAULT, PETSC_DEFAULT));
}

void SlepcEPSSolverBase::SetTol(double tol)
{
  PalacePetscCall(EPSSetTolerances(eps, tol, PETSC_DEFAULT));
  PalacePetscCall(EPSSetConvergenceTest(eps, EPS_CONV_REL));
  // PalacePetscCall(EPSSetTrackAll(eps, PETSC_TRUE));
  // PalacePetscCall(EPSSetTrueResidual(eps, PETSC_TRUE));
}

void SlepcEPSSolverBase::SetMaxIter(int maxits)
{
  PalacePetscCall(
      EPSSetTolerances(eps, PETSC_DEFAULT, (maxits > 0) ? maxits : PETSC_DEFAULT));
}

void SlepcEPSSolverBase::SetWhichEigenpairs(EigenSolverBase::WhichType type)
{
  switch (type)
  {
    case WhichType::LARGEST_MAGNITUDE:
      PalacePetscCall(EPSSetWhichEigenpairs(eps, EPS_LARGEST_MAGNITUDE));
      region = false;
      break;
    case WhichType::SMALLEST_MAGNITUDE:
      PalacePetscCall(EPSSetWhichEigenpairs(eps, EPS_SMALLEST_MAGNITUDE));
      region = false;
      break;
    case WhichType::LARGEST_REAL:
      PalacePetscCall(EPSSetWhichEigenpairs(eps, EPS_LARGEST_REAL));
      break;
    case WhichType::SMALLEST_REAL:
      PalacePetscCall(EPSSetWhichEigenpairs(eps, EPS_SMALLEST_REAL));
      break;
    case WhichType::LARGEST_IMAGINARY:
      PalacePetscCall(EPSSetWhichEigenpairs(eps, EPS_LARGEST_IMAGINARY));
      break;
    case WhichType::SMALLEST_IMAGINARY:
      PalacePetscCall(EPSSetWhichEigenpairs(eps, EPS_SMALLEST_IMAGINARY));
      break;
    case WhichType::TARGET_MAGNITUDE:
      PalacePetscCall(EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE));
      region = false;
      break;
    case WhichType::TARGET_REAL:
      PalacePetscCall(EPSSetWhichEigenpairs(eps, EPS_TARGET_REAL));
      break;
    case WhichType::TARGET_IMAGINARY:
      PalacePetscCall(EPSSetWhichEigenpairs(eps, EPS_TARGET_IMAGINARY));
      break;
    default:
      MFEM_ABORT("Which eigenpair not implemented!");
      break;
  }
}

void SlepcEPSSolverBase::SetProblemType(SlepcEigenSolver::ProblemType type)
{
  switch (type)
  {
    case ProblemType::HERMITIAN:
      PalacePetscCall(EPSSetProblemType(eps, EPS_HEP));
      break;
    case ProblemType::NON_HERMITIAN:
      PalacePetscCall(EPSSetProblemType(eps, EPS_NHEP));
      break;
    case ProblemType::GEN_HERMITIAN:
      PalacePetscCall(EPSSetProblemType(eps, EPS_GHEP));
      break;
    case ProblemType::GEN_INDEFINITE:
      PalacePetscCall(EPSSetProblemType(eps, EPS_GHIEP));
      break;
    case ProblemType::GEN_NON_HERMITIAN:
      PalacePetscCall(EPSSetProblemType(eps, EPS_GNHEP));
      // PalacePetscCall(EPSSetProblemType(eps, EPS_PGNHEP));  // If B is SPD
      break;
    default:
      MFEM_ABORT("Problem type not implemented!");
      break;
  }
}

void SlepcEPSSolverBase::SetType(SlepcEigenSolver::Type type)
{
  switch (type)
  {
    case Type::KRYLOVSCHUR:
      PalacePetscCall(EPSSetType(eps, EPSKRYLOVSCHUR));
      break;
    case Type::ARPACK:
      PalacePetscCall(EPSSetType(eps, EPSARPACK));
      break;
    case Type::POWER:
      PalacePetscCall(EPSSetType(eps, EPSPOWER));
      break;
    case Type::SUBSPACE:
      PalacePetscCall(EPSSetType(eps, EPSSUBSPACE));
      break;
    case Type::JD:
      PalacePetscCall(EPSSetType(eps, EPSJD));
      region = false;
      break;
    default:
      MFEM_ABORT("Eigenvalue solver type not implemented!");
      break;
  }
}

void SlepcEPSSolverBase::SetInitialSpace(const petsc::PetscParVector &v)
{
  if (!v0)
  {
    v0 = new petsc::PetscParVector(v);
  }
  else
  {
    MFEM_VERIFY(v.GetSize() == v0->GetSize(),
                "Invalid modification of eigenvalue problem size!");
    v0->Copy(v);
  }
  Vec is[1];
  is[0] = *v0;
  PalacePetscCall(EPSSetInitialSpace(eps, 1, is));
}

void SlepcEPSSolverBase::Customize()
{
  SlepcEigenSolver::Customize();
  PalacePetscCall(EPSSetTarget(eps, sigma / gamma));
  if (!clcustom)
  {
    PalacePetscCall(EPSSetFromOptions(eps));
    // if (print > 0)  // These are printed by PETSc linear solver
    // {
    //   PetscOptionsView(nullptr, PETSC_VIEWER_STDOUT_(GetComm()));
    //   Mpi::Print(GetComm(), "\n");
    // }
    clcustom = true;
  }
}

int SlepcEPSSolverBase::Solve()
{
  MFEM_VERIFY(A && B && opInv, "Operators are not set for SlepcEPSSolverBase!");
  PetscInt numconv;
  Customize();
  PalacePetscCall(EPSSolve(eps));
  PalacePetscCall(EPSGetConverged(eps, &numconv));
  if (print > 0)
  {
    Mpi::Print(GetComm(), "\n");
    PalacePetscCall(EPSConvergedReasonView(eps, PETSC_VIEWER_STDOUT_(GetComm())));
    Mpi::Print(GetComm(),
               " Total number of linear systems solved: {:d}\n"
               " Total number of linear solver iterations: {:d}\n",
               opInv->GetTotalNumMult(), opInv->GetTotalNumIter());
  }
  delete[] res;
  res = new PetscReal[numconv];
  for (PetscInt i = 0; i < numconv; i++)
  {
    res[i] = -1.0;
  }
  return (int)numconv;
}

void SlepcEPSSolverBase::GetEigenvalue(int i, double &eigr, double &eigi) const
{
  PetscScalar eig;
  PalacePetscCall(EPSGetEigenvalue(eps, i, &eig, nullptr));
  GetBackTransform(eig, eigr, eigi);
}

void SlepcEPSSolverBase::GetEigenvector(int i, petsc::PetscParVector &v) const
{
  PalacePetscCall(EPSGetEigenvector(eps, i, v, nullptr));
}

BV SlepcEPSSolverBase::GetBV() const
{
  BV bv;
  PalacePetscCall(EPSGetBV(eps, &bv));
  return bv;
}

ST SlepcEPSSolverBase::GetST() const
{
  ST st;
  PalacePetscCall(EPSGetST(eps, &st));
  return st;
}

RG SlepcEPSSolverBase::GetRG() const
{
  RG rg;
  PalacePetscCall(EPSGetRG(eps, &rg));
  return rg;
}

MPI_Comm SlepcEPSSolverBase::GetComm() const
{
  return eps ? PetscObjectComm(reinterpret_cast<PetscObject>(eps)) : MPI_COMM_NULL;
}

SlepcEPSSolver::SlepcEPSSolver(MPI_Comm comm, int print_lvl, const std::string &prefix)
  : SlepcEPSSolverBase(comm, print_lvl, prefix)
{
  opK = opM = nullptr;
  normK = normM = 0.0;
}

void SlepcEPSSolver::SetOperators(const petsc::PetscParMatrix &K,
                                  const petsc::PetscParMatrix &M,
                                  EigenSolverBase::ScaleType type)
{
  // Construct shell matrices for the scaled operators which define the generalized
  // eigenvalue problem.
  bool first = (opK == nullptr);
  {
    Mat A_, B_;
    MPI_Comm comm = GetComm();
    PetscInt n = K.GetNumRows();
    PalacePetscCall(
        MatCreateShell(comm, n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &A_));
    PalacePetscCall(
        MatCreateShell(comm, n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &B_));
    PalacePetscCall(MatShellSetOperation(
        A_, MATOP_MULT,
        (void (*)()) static_cast<PetscErrorCode (*)(Mat, Vec, Vec)>(&__mat_apply_EPS_A)));
    PalacePetscCall(MatShellSetOperation(
        B_, MATOP_MULT,
        (void (*)()) static_cast<PetscErrorCode (*)(Mat, Vec, Vec)>(&__mat_apply_EPS_B)));
    delete A;
    delete B;
    A = new petsc::PetscParMatrix(A_, false);  // Inherits the PETSc Mat
    B = new petsc::PetscParMatrix(B_, false);
  }
  PalacePetscCall(EPSSetOperators(eps, *A, *B));
  opK = &K;
  opM = &M;
  if (first && type != ScaleType::NONE)
  {
    normK = opK->Norm2();
    normM = opM->Norm2();
    MFEM_VERIFY(normK > 0.0 && normM > 0.0, "Invalid matrix norms for EPS scaling!");
    gamma = normK / normM;  // Store γ² for linear problem
    delta = 2.0 / normK;
  }

  // Set up workspace.
  if (!v0)
  {
    v0 = new petsc::PetscParVector(K);
  }
  if (!r0)
  {
    r0 = new petsc::PetscParVector(K);
  }

  // Configure linear solver.
  if (first)
  {
    SetPCShell((void *)this, __pc_apply_EPS);
  }
}

void SlepcEPSSolver::GetResidual(PetscScalar eig, const petsc::PetscParVector &v,
                                 petsc::PetscParVector &r) const
{
  // r = (K - λ M) v for eigenvalue λ.
  opM->Mult(v, r);
  r.Scale(-eig);
  opK->MultAdd(v, r);
}

PetscReal SlepcEPSSolver::GetBackwardScaling(PetscScalar eig) const
{
  // Make sure not to use norms from scaling as this can be confusing if they are different.
  // Note that SLEPc typically uses ||.||∞, not the 2-norm.
  if (normK <= 0.0)
  {
    normK = opK->Norm2();
  }
  if (normM <= 0.0)
  {
    normM = opM->Norm2();
  }
  return normK + PetscAbsScalar(eig) * normM;
}

SlepcPEPLinearSolver::SlepcPEPLinearSolver(MPI_Comm comm, int print_lvl,
                                           const std::string &prefix)
  : SlepcEPSSolverBase(comm, print_lvl, prefix)
{
  opK = opC = opM = nullptr;
  normK = normC = normM = 0.0;
  B0 = nullptr;
  opB = nullptr;
  x1 = x2 = y1 = y2 = z = nullptr;
}

SlepcPEPLinearSolver::~SlepcPEPLinearSolver()
{
  delete B0;
  delete x1;
  delete x2;
  delete y1;
  delete y2;
  delete z;
}

void SlepcPEPLinearSolver::SetOperators(const petsc::PetscParMatrix &K,
                                        const petsc::PetscParMatrix &C,
                                        const petsc::PetscParMatrix &M,
                                        EigenSolverBase::ScaleType type)
{
  // Construct shell matrices for the scaled linearized operators which define the block 2x2
  // eigenvalue problem.
  bool first = (opK == nullptr);
  {
    Mat A_, B_;
    MPI_Comm comm = GetComm();
    PetscInt n = K.GetNumRows();
    PalacePetscCall(
        MatCreateShell(comm, 2 * n, 2 * n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &A_));
    PalacePetscCall(
        MatCreateShell(comm, 2 * n, 2 * n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &B_));
    PalacePetscCall(
        MatShellSetOperation(A_, MATOP_MULT,
                             (void (*)()) static_cast<PetscErrorCode (*)(Mat, Vec, Vec)>(
                                 &__mat_apply_PEPLinear_L0)));
    PalacePetscCall(
        MatShellSetOperation(B_, MATOP_MULT,
                             (void (*)()) static_cast<PetscErrorCode (*)(Mat, Vec, Vec)>(
                                 &__mat_apply_PEPLinear_L1)));
    delete A;
    delete B;
    A = new petsc::PetscParMatrix(A_, false);  // Inherits the PETSc Mat
    B = new petsc::PetscParMatrix(B_, false);
  }
  PalacePetscCall(EPSSetOperators(eps, *A, *B));
  opK = &K;
  opC = &C;
  opM = &M;
  if (first && type != ScaleType::NONE)
  {
    normK = opK->Norm2();
    normC = opC->Norm2();
    normM = opM->Norm2();
    MFEM_VERIFY(normK > 0.0 && normC > 0.0 && normM > 0.0,
                "Invalid matrix norms for PEP scaling!");
    gamma = std::sqrt(normK / normM);
    delta = 2.0 / (normK + gamma * normC);
  }

  // Set up workspace.
  if (!v0)
  {
    v0 = new petsc::PetscParVector(K);
  }
  if (!r0)
  {
    r0 = new petsc::PetscParVector(K);
  }
  if (!z)
  {
    z = new petsc::PetscParVector(*A);
  }
  if (!x1)
  {
    MPI_Comm comm = GetComm();
    PetscInt n = K.GetNumRows();
    delete x1;
    delete x2;
    delete y1;
    delete y2;
    x1 = new petsc::PetscParVector(comm, n, PETSC_DECIDE, nullptr);
    x2 = new petsc::PetscParVector(comm, n, PETSC_DECIDE, nullptr);
    y1 = new petsc::PetscParVector(comm, n, PETSC_DECIDE, nullptr);
    y2 = new petsc::PetscParVector(comm, n, PETSC_DECIDE, nullptr);
  }

  // Configure linear solver.
  if (first)
  {
    SetPCShell((void *)this, __pc_apply_PEPLinear);
  }
}

void SlepcPEPLinearSolver::SetBMat(const petsc::PetscParMatrix &B)
{
  // Construct an SPD linearized mass matrix for weighted inner products.
  Mat B0_;
  MPI_Comm comm = GetComm();
  PetscInt n = B.GetNumRows();
  PalacePetscCall(
      MatCreateShell(comm, 2 * n, 2 * n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &B0_));
  PalacePetscCall(
      MatShellSetOperation(B0_, MATOP_MULT,
                           (void (*)()) static_cast<PetscErrorCode (*)(Mat, Vec, Vec)>(
                               &__mat_apply_PEPLinear_B)));
  delete B0;
  B0 = new petsc::PetscParMatrix(B0_, false);  // Inherits the PETSc Mat
  opB = &B;
  SlepcEigenSolver::SetBMat(*B0);
}

void SlepcPEPLinearSolver::SetInitialSpace(const petsc::PetscParVector &v)
{
  if (!z)
  {
    z = new petsc::PetscParVector(v.GetComm(), 2 * v.GetSize(), PETSC_DECIDE);
  }
  else
  {
    MFEM_VERIFY(2 * v.GetSize() == z->GetSize(),
                "Invalid modification of eigenvalue problem size!");
  }
  {
    PetscScalar *pz = GetBlocks(*z, *x1, *x2);
    x1->Copy(v);
    x2->SetZero();
    RestoreBlocks(pz, *z, *x1, *x2);
  }
  Vec is[1];
  is[0] = *z;
  PalacePetscCall(EPSSetInitialSpace(eps, 1, is));
}

void SlepcPEPLinearSolver::GetEigenvector(int i, petsc::PetscParVector &v) const
{
  // Select the most accurate v for z = [v₁; v₂] from the linearized eigenvalue problem.
  PalacePetscCall(EPSGetEigenvector(eps, i, *z, nullptr));
  const PetscScalar *pz = GetBlocksRead(*z, *x1, *x2);
  {
    if (opB)
    {
      x1->Normalize(*opB, *r0);
    }
    else
    {
      x1->Normalize();
    }
    v.Copy(*x1);
  }
  RestoreBlocksRead(pz, *z, *x1, *x2);
}

void SlepcPEPLinearSolver::GetResidual(PetscScalar eig, const petsc::PetscParVector &v,
                                       petsc::PetscParVector &r) const
{
  // r = P(λ) v = (K + λ C + λ² M) v for eigenvalue λ.
  opM->Mult(v, r);
  r.Scale(eig);
  opC->MultAdd(v, r);
  r.Scale(eig);
  opK->MultAdd(v, r);
}

PetscReal SlepcPEPLinearSolver::GetBackwardScaling(PetscScalar eig) const
{
  // Make sure not to use norms from scaling as this can be confusing if they are different.
  // Note that SLEPc typically uses ||.||∞, not the 2-norm.
  if (normK <= 0.0)
  {
    normK = opK->Norm2();
  }
  if (normC <= 0.0)
  {
    normC = opC->Norm2();
  }
  if (normM <= 0.0)
  {
    normM = opM->Norm2();
  }
  PetscReal t = PetscAbsScalar(eig);
  return normK + t * normC + t * t * normM;
}

PetscScalar *SlepcPEPLinearSolver::GetBlocks(petsc::PetscParVector &v,
                                             petsc::PetscParVector &v1,
                                             petsc::PetscParVector &v2) const
{
  PetscInt n1 = v1.GetSize(), n2 = v2.GetSize();
  MFEM_VERIFY(n1 + n2 == v.GetSize(), "Unexpected size in PEP linearization!");
  PetscScalar *pv = v.GetArray();
  v1.PlaceArray(pv);
  v2.PlaceArray(pv + n1);
  return pv;
}

const PetscScalar *SlepcPEPLinearSolver::GetBlocksRead(const petsc::PetscParVector &v,
                                                       petsc::PetscParVector &v1,
                                                       petsc::PetscParVector &v2) const
{
  PetscInt n1 = v1.GetSize(), n2 = v2.GetSize();
  MFEM_VERIFY(n1 + n2 == v.GetSize(), "Unexpected size in PEP linearization!");
  const PetscScalar *pv = v.GetArrayRead();
  v1.PlaceArray(pv);
  v2.PlaceArray(pv + n1);
  return pv;
}

void SlepcPEPLinearSolver::RestoreBlocks(PetscScalar *pv, petsc::PetscParVector &v,
                                         petsc::PetscParVector &v1,
                                         petsc::PetscParVector &v2) const
{
  v1.ResetArray();
  v2.ResetArray();
  v.RestoreArray(pv);
}

void SlepcPEPLinearSolver::RestoreBlocksRead(const PetscScalar *pv,
                                             const petsc::PetscParVector &v,
                                             petsc::PetscParVector &v1,
                                             petsc::PetscParVector &v2) const
{
  v1.ResetArray();
  v2.ResetArray();
  v.RestoreArrayRead(pv);
}

// PEP specific methods

SlepcPEPSolverBase::SlepcPEPSolverBase(MPI_Comm comm, int print_lvl,
                                       const std::string &prefix)
  : SlepcEigenSolver(print_lvl)
{
  PalacePetscCall(PEPCreate(comm, &pep));
  PalacePetscCall(PEPSetOptionsPrefix(pep, prefix.c_str()));
  if (print > 0)
  {
    std::string opts = "-pep_monitor";
    if (print > 2)
    {
      opts.append(" -pep_view");
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
  A0 = A1 = A2 = nullptr;
}

SlepcPEPSolverBase::~SlepcPEPSolverBase()
{
  MPI_Comm comm;
  PalacePetscCall(PetscObjectGetComm(reinterpret_cast<PetscObject>(pep), &comm));
  PalacePetscCall(PEPDestroy(&pep));
  delete A0;
  delete A1;
  delete A2;
}

void SlepcPEPSolverBase::SetNumModes(int numeig, int numvec)
{
  PalacePetscCall(
      PEPSetDimensions(pep, numeig, (numvec > 0) ? numvec : PETSC_DEFAULT, PETSC_DEFAULT));
}

void SlepcPEPSolverBase::SetTol(double tol)
{
  PalacePetscCall(PEPSetTolerances(pep, tol, PETSC_DEFAULT));
  PalacePetscCall(PEPSetConvergenceTest(pep, PEP_CONV_REL));
  // PalacePetscCall(PEPSetTrackAll(pep, PETSC_TRUE));
}

void SlepcPEPSolverBase::SetMaxIter(int maxits)
{
  PalacePetscCall(
      PEPSetTolerances(pep, PETSC_DEFAULT, (maxits > 0) ? maxits : PETSC_DEFAULT));
}

void SlepcPEPSolverBase::SetWhichEigenpairs(EigenSolverBase::WhichType type)
{
  switch (type)
  {
    case WhichType::LARGEST_MAGNITUDE:
      PalacePetscCall(PEPSetWhichEigenpairs(pep, PEP_LARGEST_MAGNITUDE));
      region = false;
      break;
    case WhichType::SMALLEST_MAGNITUDE:
      PalacePetscCall(PEPSetWhichEigenpairs(pep, PEP_SMALLEST_MAGNITUDE));
      region = false;
      break;
    case WhichType::LARGEST_REAL:
      PalacePetscCall(PEPSetWhichEigenpairs(pep, PEP_LARGEST_REAL));
      break;
    case WhichType::SMALLEST_REAL:
      PalacePetscCall(PEPSetWhichEigenpairs(pep, PEP_SMALLEST_REAL));
      break;
    case WhichType::LARGEST_IMAGINARY:
      PalacePetscCall(PEPSetWhichEigenpairs(pep, PEP_LARGEST_IMAGINARY));
      break;
    case WhichType::SMALLEST_IMAGINARY:
      PalacePetscCall(PEPSetWhichEigenpairs(pep, PEP_SMALLEST_IMAGINARY));
      break;
    case WhichType::TARGET_MAGNITUDE:
      PalacePetscCall(PEPSetWhichEigenpairs(pep, PEP_TARGET_MAGNITUDE));
      region = false;
      break;
    case WhichType::TARGET_REAL:
      PalacePetscCall(PEPSetWhichEigenpairs(pep, PEP_TARGET_REAL));
      break;
    case WhichType::TARGET_IMAGINARY:
      PalacePetscCall(PEPSetWhichEigenpairs(pep, PEP_TARGET_IMAGINARY));
      break;
    default:
      MFEM_ABORT("Which eigenpair not implemented!");
      break;
  }
}

void SlepcPEPSolverBase::SetProblemType(SlepcEigenSolver::ProblemType type)
{
  switch (type)
  {
    case ProblemType::HERMITIAN:
    case ProblemType::GEN_HERMITIAN:
      PalacePetscCall(PEPSetProblemType(pep, PEP_HERMITIAN));
      break;
    case ProblemType::NON_HERMITIAN:
    case ProblemType::GEN_INDEFINITE:
    case ProblemType::GEN_NON_HERMITIAN:
      PalacePetscCall(PEPSetProblemType(pep, PEP_GENERAL));
      break;
    case ProblemType::HYPERBOLIC:
      PalacePetscCall(PEPSetProblemType(pep, PEP_HYPERBOLIC));
      break;
    case ProblemType::GYROSCOPIC:
      PalacePetscCall(PEPSetProblemType(pep, PEP_GYROSCOPIC));
      break;
    default:
      MFEM_ABORT("Problem type not implemented!");
      break;
  }
}

void SlepcPEPSolverBase::SetType(SlepcEigenSolver::Type type)
{
  switch (type)
  {
    case Type::TOAR:
      PalacePetscCall(PEPSetType(pep, PEPTOAR));
      break;
    case Type::STOAR:
      PalacePetscCall(PEPSetType(pep, PEPSTOAR));
      break;
    case Type::QARNOLDI:
      PalacePetscCall(PEPSetType(pep, PEPQARNOLDI));
      break;
    case Type::JD:
      PalacePetscCall(PEPSetType(pep, PEPJD));
      region = false;
      break;
    default:
      MFEM_ABORT("Eigenvalue solver type not implemented!");
      break;
  }
}

void SlepcPEPSolverBase::SetInitialSpace(const petsc::PetscParVector &v)
{
  if (!v0)
  {
    v0 = new petsc::PetscParVector(v);
  }
  else
  {
    MFEM_VERIFY(v.GetSize() == v0->GetSize(),
                "Invalid modification of eigenvalue problem size!");
    v0->Copy(v);
  }
  Vec is[1];
  is[0] = *v0;
  PalacePetscCall(PEPSetInitialSpace(pep, 1, is));
}

void SlepcPEPSolverBase::Customize()
{
  SlepcEigenSolver::Customize();
  PalacePetscCall(PEPSetTarget(pep, sigma / gamma));
  if (!clcustom)
  {
    PalacePetscCall(PEPSetFromOptions(pep));
    // if (print > 0)  // These are printed by PETSc linear solver
    // {
    //   PetscOptionsView(nullptr, PETSC_VIEWER_STDOUT_(GetComm()));
    //   Mpi::Print(GetComm(), "\n");
    // }
    clcustom = true;
  }
}

int SlepcPEPSolverBase::Solve()
{
  MFEM_VERIFY(A0 && A1 && A2 && opInv, "Operators are not set for SlepcPEPSolverBase!");
  PetscInt numconv;
  Customize();
  PalacePetscCall(PEPSolve(pep));
  PalacePetscCall(PEPGetConverged(pep, &numconv));
  if (print > 0)
  {
    Mpi::Print(GetComm(), "\n");
    PalacePetscCall(PEPConvergedReasonView(pep, PETSC_VIEWER_STDOUT_(GetComm())));
    Mpi::Print(GetComm(),
               " Total number of linear systems solved: {:d}\n"
               " Total number of linear solver iterations: {:d}\n",
               opInv->GetTotalNumMult(), opInv->GetTotalNumIter());
  }
  delete[] res;
  res = new PetscReal[numconv];
  for (PetscInt i = 0; i < numconv; i++)
  {
    res[i] = -1.0;
  }
  return (int)numconv;
}

void SlepcPEPSolverBase::GetEigenvalue(int i, double &eigr, double &eigi) const
{
  PetscScalar eig;
  PalacePetscCall(PEPGetEigenpair(pep, i, &eig, nullptr, nullptr, nullptr));
  GetBackTransform(eig, eigr, eigi);
}

void SlepcPEPSolverBase::GetEigenvector(int i, petsc::PetscParVector &v) const
{
  PalacePetscCall(PEPGetEigenpair(pep, i, nullptr, nullptr, v, nullptr));
}

BV SlepcPEPSolverBase::GetBV() const
{
  BV bv;
  PalacePetscCall(PEPGetBV(pep, &bv));
  return bv;
}

ST SlepcPEPSolverBase::GetST() const
{
  ST st;
  PalacePetscCall(PEPGetST(pep, &st));
  return st;
}

RG SlepcPEPSolverBase::GetRG() const
{
  RG rg;
  PalacePetscCall(PEPGetRG(pep, &rg));
  return rg;
}

MPI_Comm SlepcPEPSolverBase::GetComm() const
{
  return pep ? PetscObjectComm(reinterpret_cast<PetscObject>(pep)) : MPI_COMM_NULL;
}

SlepcPEPSolver::SlepcPEPSolver(MPI_Comm comm, int print_lvl, const std::string &prefix)
  : SlepcPEPSolverBase(comm, print_lvl, prefix)
{
  opK = opC = opM = nullptr;
  normK = normC = normM = 0.0;
}

void SlepcPEPSolver::SetOperators(const petsc::PetscParMatrix &K,
                                  const petsc::PetscParMatrix &C,
                                  const petsc::PetscParMatrix &M,
                                  EigenSolverBase::ScaleType type)
{
  // Construct shell matrices for the scaled operators which define the quadratic polynomial
  // eigenvalue problem.
  bool first = (opK == nullptr);
  {
    Mat A0_, A1_, A2_;
    MPI_Comm comm = GetComm();
    PetscInt n = K.GetNumRows();
    PalacePetscCall(
        MatCreateShell(comm, n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &A0_));
    PalacePetscCall(
        MatCreateShell(comm, n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &A1_));
    PalacePetscCall(
        MatCreateShell(comm, n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &A2_));
    PalacePetscCall(MatShellSetOperation(
        A0_, MATOP_MULT,
        (void (*)()) static_cast<PetscErrorCode (*)(Mat, Vec, Vec)>(&__mat_apply_PEP_A0)));
    PalacePetscCall(MatShellSetOperation(
        A1_, MATOP_MULT,
        (void (*)()) static_cast<PetscErrorCode (*)(Mat, Vec, Vec)>(&__mat_apply_PEP_A1)));
    PalacePetscCall(MatShellSetOperation(
        A2_, MATOP_MULT,
        (void (*)()) static_cast<PetscErrorCode (*)(Mat, Vec, Vec)>(&__mat_apply_PEP_A2)));
    delete A0;
    delete A1;
    delete A2;
    A0 = new petsc::PetscParMatrix(A0_, false);  // Inherits the PETSc Mat
    A1 = new petsc::PetscParMatrix(A1_, false);
    A2 = new petsc::PetscParMatrix(A2_, false);
  }
  Mat A[3] = {*A0, *A1, *A2};
  PalacePetscCall(PEPSetOperators(pep, 3, A));
  opK = &K;
  opC = &C;
  opM = &M;
  if (first && type != ScaleType::NONE)
  {
    normK = opK->Norm2();
    normC = opC->Norm2();
    normM = opM->Norm2();
    MFEM_VERIFY(normK > 0.0 && normC > 0.0 && normM > 0.0,
                "Invalid matrix norms for PEP scaling!");
    gamma = std::sqrt(normK / normM);
    delta = 2.0 / (normK + gamma * normC);
  }

  // Set up workspace.
  if (!v0)
  {
    v0 = new petsc::PetscParVector(K);
  }
  if (!r0)
  {
    r0 = new petsc::PetscParVector(K);
  }

  // Configure linear solver.
  if (first)
  {
    SetPCShell((void *)this, __pc_apply_PEP);
  }
}

void SlepcPEPSolver::GetResidual(PetscScalar eig, const petsc::PetscParVector &v,
                                 petsc::PetscParVector &r) const
{
  // r = P(λ) v = (K + λ C + λ² M) v for eigenvalue λ.
  opM->Mult(v, r);
  r.Scale(eig);
  opC->MultAdd(v, r);
  r.Scale(eig);
  opK->MultAdd(v, r);
}

PetscReal SlepcPEPSolver::GetBackwardScaling(PetscScalar eig) const
{
  // Make sure not to use norms from scaling as this can be confusing if they are different.
  // Note that SLEPc typically uses ||.||∞, not Frobenius.
  if (normK <= 0.0)
  {
    normK = opK->NormInf();
  }
  if (normC <= 0.0)
  {
    normC = opC->NormInf();
  }
  if (normM <= 0.0)
  {
    normM = opM->NormInf();
  }
  PetscReal t = PetscAbsScalar(eig);
  return normK + t * normC + t * t * normM;
}

}  // namespace palace::slepc

PetscErrorCode __mat_apply_EPS_A(Mat A, Vec x, Vec y)
{
  // Apply the operator: K (no transform) or M .
  palace::slepc::SlepcEPSSolver *slepc;
  palace::petsc::PetscParVector xx(x, true), yy(y, true);
  PetscFunctionBeginUser;

  PetscCall(MatShellGetContext(A, (void **)&slepc));
  MFEM_VERIFY(slepc, "Invalid PETSc shell matrix context for SLEPc!");
  {
    slepc->GetOpK()->Mult(xx, yy);
    yy.Scale(slepc->GetScalingDelta());
  }
  PetscFunctionReturn(0);
}

PetscErrorCode __mat_apply_EPS_B(Mat A, Vec x, Vec y)
{
  // Apply the operator: M (no transform) or (K - σ M) .
  palace::slepc::SlepcEPSSolver *slepc;
  palace::petsc::PetscParVector xx(x, true), yy(y, true);
  PetscFunctionBeginUser;

  PetscCall(MatShellGetContext(A, (void **)&slepc));
  MFEM_VERIFY(slepc, "Invalid PETSc shell matrix context for SLEPc!");
  {
    slepc->GetOpM()->Mult(xx, yy);
    yy.Scale(slepc->GetScalingDelta() * slepc->GetScalingGamma());
  }
  PetscFunctionReturn(0);
}

PetscErrorCode __pc_apply_EPS(PC pc, Vec x, Vec y)
{
  // Solve the linear system associated with the generalized eigenvalue problem: y = M⁻¹ x ,
  // or shift-and-invert spectral transformation: y =(K - σ M)⁻¹ x . Enforces the
  // divergence-free constraint using the supplied projector.
  palace::slepc::SlepcEPSSolver *slepc;
  palace::petsc::PetscParVector xx(x, true), yy(y, true);
  PetscFunctionBeginUser;

  PetscCall(PCShellGetContext(pc, (void **)&slepc));
  MFEM_VERIFY(slepc, "Invalid PETSc shell PC context for SLEPc!");
  slepc->GetKspSolver()->Mult(xx, yy);
  if (!slepc->IsShiftInvert())
  {
    yy.Scale(1.0 / (slepc->GetScalingDelta() * slepc->GetScalingGamma()));
  }
  else
  {
    yy.Scale(1.0 / slepc->GetScalingDelta());
  }

  // Debug
  // Mpi::Print(" Before projection: {:e}\n", yy.Norml2());

  if (slepc->GetDivFreeSolver())
  {
    slepc->GetDivFreeSolver()->Mult(yy);
  }

  // Debug
  // Mpi::Print(" After projection: {:e}\n", yy.Norml2());

  PetscFunctionReturn(0);
}

PetscErrorCode __mat_apply_PEPLinear_L0(Mat A, Vec x, Vec y)
{
  // Apply the linearized operator: L₀ (no transform) or L₁ . With:
  //               L₀ = [  0   I ]    L₁ = [ I  0 ]
  //                    [ -K  -C ] ,       [ 0  M ] .
  palace::slepc::SlepcPEPLinearSolver *slepc;
  palace::petsc::PetscParVector xx(x, true), yy(y, true);
  const PetscScalar *px;
  PetscScalar *py;
  PetscFunctionBeginUser;

  PetscCall(MatShellGetContext(A, (void **)&slepc));
  MFEM_VERIFY(slepc, "Invalid PETSc shell matrix context for SLEPc!");
  px = slepc->GetBlocksRead(xx, *slepc->GetX1(), *slepc->GetX2());
  py = slepc->GetBlocks(yy, *slepc->GetY1(), *slepc->GetY2());
  {
    slepc->GetY1()->Copy(*slepc->GetX2());
    slepc->GetOpC()->Mult(*slepc->GetX2(), *slepc->GetY2());
    slepc->GetY2()->Scale(slepc->GetScalingGamma());
    slepc->GetOpK()->MultAdd(*slepc->GetX1(), *slepc->GetY2());
    slepc->GetY2()->Scale(-slepc->GetScalingDelta());
  }
  slepc->RestoreBlocksRead(px, xx, *slepc->GetX1(), *slepc->GetX2());
  slepc->RestoreBlocks(py, yy, *slepc->GetY1(), *slepc->GetY2());
  PetscFunctionReturn(0);
}

PetscErrorCode __mat_apply_PEPLinear_L1(Mat A, Vec x, Vec y)
{
  // Apply the linearized operator: L₁ (no transform) or (L₀ - σ L₁) . With:
  //               L₀ = [  0   I ]    L₁ = [ I  0 ]
  //                    [ -K  -C ] ,       [ 0  M ] .
  palace::slepc::SlepcPEPLinearSolver *slepc;
  palace::petsc::PetscParVector xx(x, true), yy(y, true);
  const PetscScalar *px;
  PetscScalar *py;
  PetscFunctionBeginUser;

  PetscCall(MatShellGetContext(A, (void **)&slepc));
  MFEM_VERIFY(slepc, "Invalid PETSc shell matrix context for SLEPc!");
  px = slepc->GetBlocksRead(xx, *slepc->GetX1(), *slepc->GetX2());
  py = slepc->GetBlocks(yy, *slepc->GetY1(), *slepc->GetY2());
  {
    slepc->GetY1()->Copy(*slepc->GetX1());
    slepc->GetOpM()->Mult(*slepc->GetX2(), *slepc->GetY2());
    slepc->GetY2()->Scale(slepc->GetScalingDelta() * slepc->GetScalingGamma() *
                          slepc->GetScalingGamma());
  }
  slepc->RestoreBlocksRead(px, xx, *slepc->GetX1(), *slepc->GetX2());
  slepc->RestoreBlocks(py, yy, *slepc->GetY1(), *slepc->GetY2());
  PetscFunctionReturn(0);
}

PetscErrorCode __mat_apply_PEPLinear_B(Mat A, Vec x, Vec y)
{
  // Apply the linearized mass matrix L₁ using the supplied SPD mass matrix.
  palace::slepc::SlepcPEPLinearSolver *slepc;
  palace::petsc::PetscParVector xx(x, true), yy(y, true);
  const PetscScalar *px;
  PetscScalar *py;
  PetscFunctionBeginUser;

  PetscCall(MatShellGetContext(A, (void **)&slepc));
  MFEM_VERIFY(slepc, "Invalid PETSc shell matrix context for SLEPc!");
  px = slepc->GetBlocksRead(xx, *slepc->GetX1(), *slepc->GetX2());
  py = slepc->GetBlocks(yy, *slepc->GetY1(), *slepc->GetY2());
  {
    slepc->GetY1()->Copy(*slepc->GetX1());
    slepc->GetOpB()->Mult(*slepc->GetX2(), *slepc->GetY2());
    slepc->GetY2()->Scale(slepc->GetScalingDelta() * slepc->GetScalingGamma() *
                          slepc->GetScalingGamma());
  }
  slepc->RestoreBlocksRead(px, xx, *slepc->GetX1(), *slepc->GetX2());
  slepc->RestoreBlocks(py, yy, *slepc->GetY1(), *slepc->GetY2());
  PetscFunctionReturn(0);
}

PetscErrorCode __pc_apply_PEPLinear(PC pc, Vec x, Vec y)
{
  // Solve the linear system associated with the generalized eigenvalue problem after
  // linearization: y = L₁⁻¹ x , or with the shift-and-invert spectral transformation: y =
  // (L₀ - σ L₁)⁻¹ x . Enforces the divergence-free constraint using the supplied
  // projectors.
  palace::slepc::SlepcPEPLinearSolver *slepc;
  palace::petsc::PetscParVector xx(x, true), yy(y, true);
  const PetscScalar *px;
  PetscScalar *py;
  PetscScalar sigma;
  PetscReal gamma, delta;
  PetscFunctionBeginUser;

  PetscCall(PCShellGetContext(pc, (void **)&slepc));
  MFEM_VERIFY(slepc, "Invalid PETSc shell PC context for SLEPc!");
  sigma = slepc->GetTarget();
  gamma = slepc->GetScalingGamma();
  delta = slepc->GetScalingDelta();
  px = slepc->GetBlocksRead(xx, *slepc->GetX1(), *slepc->GetX2());
  py = slepc->GetBlocks(yy, *slepc->GetY1(), *slepc->GetY2());
  if (!slepc->IsShiftInvert())
  {
    slepc->GetKspSolver()->Mult(*slepc->GetX2(), *slepc->GetY2());
    slepc->GetY2()->Scale(1.0 / (delta * gamma * gamma));
    if (slepc->GetDivFreeSolver())
    {
      slepc->GetDivFreeSolver()->Mult(*slepc->GetY2());
    }
    slepc->GetY1()->Copy(*slepc->GetX1());
    if (slepc->GetDivFreeSolver())
    {
      slepc->GetDivFreeSolver()->Mult(*slepc->GetY1());
    }
  }
  else
  {
    slepc->GetY1()->AXPBY(-sigma / (delta * gamma), *slepc->GetX2(), 0.0);  // Temporarily
    slepc->GetOpK()->MultAdd(*slepc->GetX1(), *slepc->GetY1());
    slepc->GetKspSolver()->Mult(*slepc->GetY1(), *slepc->GetY2());

    // Debug
    // Mpi::Print(" Before projection: {:e}\n", slepc->GetY2()->Norml2());

    if (slepc->GetDivFreeSolver())
    {
      slepc->GetDivFreeSolver()->Mult(*slepc->GetY2());
    }

    // Debug
    // Mpi::Print(" After projection: {:e}\n", slepc->GetY2()->Norml2());

    slepc->GetY1()->AXPBYPCZ(gamma / sigma, *slepc->GetY2(), -gamma / sigma,
                             *slepc->GetX1(), 0.0);

    // Debug
    // Mpi::Print(" Before projection: {:e}\n", slepc->GetY1()->Norml2());

    if (slepc->GetDivFreeSolver())
    {
      slepc->GetDivFreeSolver()->Mult(*slepc->GetY1());
    }

    // Debug
    // Mpi::Print(" After projection: {:e}\n", slepc->GetY1()->Norml2());
  }
  slepc->RestoreBlocksRead(px, xx, *slepc->GetX1(), *slepc->GetX2());
  slepc->RestoreBlocks(py, yy, *slepc->GetY1(), *slepc->GetY2());
  PetscFunctionReturn(0);
}

PetscErrorCode __mat_apply_PEP_A0(Mat A, Vec x, Vec y)
{
  // Apply the operator: K (no transform) or M .
  palace::slepc::SlepcPEPSolver *slepc;
  palace::petsc::PetscParVector xx(x, true), yy(y, true);
  PetscFunctionBeginUser;

  PetscCall(MatShellGetContext(A, (void **)&slepc));
  MFEM_VERIFY(slepc, "Invalid PETSc shell matrix context for SLEPc!");
  {
    slepc->GetOpK()->Mult(xx, yy);
    yy.Scale(slepc->GetScalingDelta());
  }
  PetscFunctionReturn(0);
}

PetscErrorCode __mat_apply_PEP_A1(Mat A, Vec x, Vec y)
{
  // Apply the operator: C (no transform) or (C + 2σ M) .
  palace::slepc::SlepcPEPSolver *slepc;
  palace::petsc::PetscParVector xx(x, true), yy(y, true);
  PetscFunctionBeginUser;

  PetscCall(MatShellGetContext(A, (void **)&slepc));
  MFEM_VERIFY(slepc, "Invalid PETSc shell matrix context for SLEPc!");
  {
    slepc->GetOpC()->Mult(xx, yy);
    yy.Scale(slepc->GetScalingDelta() * slepc->GetScalingGamma());
  }
  PetscFunctionReturn(0);
}

PetscErrorCode __mat_apply_PEP_A2(Mat A, Vec x, Vec y)
{
  // Apply the operator: M (no transform) or (K + σ C + σ² M) .
  palace::slepc::SlepcPEPSolver *slepc;
  palace::petsc::PetscParVector xx(x, true), yy(y, true);
  PetscFunctionBeginUser;

  PetscCall(MatShellGetContext(A, (void **)&slepc));
  MFEM_VERIFY(slepc, "Invalid PETSc shell matrix context for SLEPc!");
  {
    slepc->GetOpM()->Mult(xx, yy);
    yy.Scale(slepc->GetScalingDelta() * slepc->GetScalingGamma() *
             slepc->GetScalingGamma());
  }
  PetscFunctionReturn(0);
}

PetscErrorCode __pc_apply_PEP(PC pc, Vec x, Vec y)
{
  // Solve the linear system associated with the generalized eigenvalue problem: y = M⁻¹ x ,
  // or shift-and-invert spectral transformation: y = P(σ)⁻¹ x . Enforces the
  // divergence-free constraint using the supplied projector.
  palace::slepc::SlepcPEPSolver *slepc;
  palace::petsc::PetscParVector xx(x, true), yy(y, true);
  PetscFunctionBeginUser;

  PetscCall(PCShellGetContext(pc, (void **)&slepc));
  MFEM_VERIFY(slepc, "Invalid PETSc shell PC context for SLEPc!");
  slepc->GetKspSolver()->Mult(xx, yy);
  if (!slepc->IsShiftInvert())
  {
    yy.Scale(1.0 / (slepc->GetScalingDelta() * slepc->GetScalingGamma() *
                    slepc->GetScalingGamma()));
  }
  else
  {
    yy.Scale(1.0 / slepc->GetScalingDelta());
  }

  // Debug
  // Mpi::Print(" Before projection: {:e}\n", yy);

  if (slepc->GetDivFreeSolver())
  {
    slepc->GetDivFreeSolver()->Mult(yy);
  }

  // Debug
  // Mpi::Print(" After projection: {:e}\n", yy);

  PetscFunctionReturn(0);
}

#endif
