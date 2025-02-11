// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "slepc.hpp"

#if defined(PALACE_WITH_SLEPC)

#include <algorithm>
#include <petsc.h>
#include <slepc.h>
#include <mfem.hpp>
#include "linalg/divfree.hpp"
#include "utils/communication.hpp"

static PetscErrorCode __mat_apply_EPS_A0(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_EPS_A1(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_EPS_B(Mat, Vec, Vec);
static PetscErrorCode __pc_apply_EPS(PC, Vec, Vec);
static PetscErrorCode __mat_apply_PEPLinear_L0(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_PEPLinear_L1(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_PEPLinear_B(Mat, Vec, Vec);
static PetscErrorCode __pc_apply_PEPLinear(PC, Vec, Vec);
static PetscErrorCode __mat_apply_PEP_A0(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_PEP_A1(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_PEP_A2(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_PEP_B(Mat, Vec, Vec);
static PetscErrorCode __pc_apply_PEP(PC, Vec, Vec);

namespace
{

inline PetscErrorCode FromPetscVec(Vec x, palace::ComplexVector &y, int block = 0,
                                   int nblocks = 1)
{
  PetscInt n;
  const PetscScalar *px;
  PetscMemType mtype;
  PetscCall(VecGetLocalSize(x, &n));
  MFEM_ASSERT(y.Size() * nblocks == n,
              "Invalid size mismatch for PETSc vector conversion!");
  PetscCall(VecGetArrayReadAndMemType(x, &px, &mtype));
  y.Set(px + block * n / nblocks, n / nblocks, PetscMemTypeDevice(mtype));
  PetscCall(VecRestoreArrayReadAndMemType(x, &px));
  return PETSC_SUCCESS;
}

inline PetscErrorCode FromPetscVec(Vec x, palace::ComplexVector &y1,
                                   palace::ComplexVector &y2)
{
  PetscInt n;
  const PetscScalar *px;
  PetscMemType mtype;
  PetscCall(VecGetLocalSize(x, &n));
  MFEM_ASSERT(y1.Size() == n / 2 && y2.Size() == n / 2,
              "Invalid size mismatch for PETSc vector conversion!");
  PetscCall(VecGetArrayReadAndMemType(x, &px, &mtype));
  y1.Set(px, n / 2, PetscMemTypeDevice(mtype));
  y2.Set(px + n / 2, n / 2, PetscMemTypeDevice(mtype));
  PetscCall(VecRestoreArrayReadAndMemType(x, &px));
  return PETSC_SUCCESS;
}

inline PetscErrorCode ToPetscVec(const palace::ComplexVector &x, Vec y, int block = 0,
                                 int nblocks = 1)
{
  PetscInt n;
  PetscScalar *py;
  PetscMemType mtype;
  PetscCall(VecGetLocalSize(y, &n));
  MFEM_ASSERT(x.Size() * nblocks == n,
              "Invalid size mismatch for PETSc vector conversion!");
  PetscCall(VecGetArrayWriteAndMemType(y, &py, &mtype));
  x.Get(py + block * n / nblocks, n / nblocks, PetscMemTypeDevice(mtype));
  PetscCall(VecRestoreArrayWriteAndMemType(y, &py));
  return PETSC_SUCCESS;
}

inline PetscErrorCode ToPetscVec(const palace::ComplexVector &x1,
                                 const palace::ComplexVector &x2, Vec y)
{
  PetscInt n;
  PetscScalar *py;
  PetscMemType mtype;
  PetscCall(VecGetLocalSize(y, &n));
  MFEM_ASSERT(x1.Size() == n / 2 && x2.Size() == n / 2,
              "Invalid size mismatch for PETSc vector conversion!");
  PetscCall(VecGetArrayWriteAndMemType(y, &py, &mtype));
  x1.Get(py, n / 2, PetscMemTypeDevice(mtype));
  x2.Get(py + n / 2, n / 2, PetscMemTypeDevice(mtype));
  PetscCall(VecRestoreArrayWriteAndMemType(y, &py));
  return PETSC_SUCCESS;
}

}  // namespace

namespace palace::slepc
{

namespace
{

inline PetscErrorCode ConfigurePetscDevice()
{
  // Tell PETSc to use the same CUDA or HIP device as MFEM.
  if (mfem::Device::Allows(mfem::Backend::CUDA_MASK))
  {
    PetscCall(PetscOptionsSetValue(NULL, "-use_gpu_aware_mpi",
                                   mfem::Device::GetGPUAwareMPI() ? "1" : "0"));
    PetscCall(PetscOptionsSetValue(NULL, "-device_select_cuda",
                                   std::to_string(mfem::Device::GetId()).c_str()));
    // PetscCall(PetscOptionsSetValue(NULL, "-bv_type", "svec"));
    // PetscCall(PetscOptionsSetValue(NULL, "-vec_type", "cuda"));
  }
  if (mfem::Device::Allows(mfem::Backend::HIP_MASK))
  {
    PetscCall(PetscOptionsSetValue(NULL, "-use_gpu_aware_mpi",
                                   mfem::Device::GetGPUAwareMPI() ? "1" : "0"));
    PetscCall(PetscOptionsSetValue(NULL, "-device_select_hip",
                                   std::to_string(mfem::Device::GetId()).c_str()));
    // PetscCall(PetscOptionsSetValue(NULL, "-bv_type", "svec"));
    // PetscCall(PetscOptionsSetValue(NULL, "-vec_type", "hip"));
  }
  return PETSC_SUCCESS;
}

inline VecType PetscVecType()
{
  if (mfem::Device::Allows(mfem::Backend::CUDA_MASK))
  {
    return VECCUDA;
  }
  if (mfem::Device::Allows(mfem::Backend::HIP_MASK))
  {
    return VECHIP;
  }
  return VECSTANDARD;
}

struct MatShellContext
{
  const ComplexOperator &A;
  ComplexVector &x, &y;
};

PetscErrorCode __mat_apply_shell(Mat A, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  MatShellContext *ctx;
  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context for SLEPc!");

  PetscCall(FromPetscVec(x, ctx->x));
  ctx->A.Mult(ctx->x, ctx->y);
  PetscCall(ToPetscVec(ctx->y, y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __mat_apply_transpose_shell(Mat A, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  MatShellContext *ctx;
  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context for SLEPc!");

  PetscCall(FromPetscVec(x, ctx->x));
  ctx->A.MultTranspose(ctx->x, ctx->y);
  PetscCall(ToPetscVec(ctx->y, y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __mat_apply_hermitian_transpose_shell(Mat A, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  MatShellContext *ctx;
  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context for SLEPc!");

  PetscCall(FromPetscVec(x, ctx->x));
  ctx->A.MultHermitianTranspose(ctx->x, ctx->y);
  PetscCall(ToPetscVec(ctx->y, y));

  PetscFunctionReturn(PETSC_SUCCESS);
};

inline void ConfigurePCShell(ST st, void *ctx, PetscErrorCode (*__pc_apply)(PC, Vec, Vec))
{
  KSP ksp;
  PC pc;
  PalacePetscCall(STGetKSP(st, &ksp));
  PalacePetscCall(KSPGetPC(ksp, &pc));
  PalacePetscCall(PCSetType(pc, PCSHELL));
  PalacePetscCall(PCShellSetContext(pc, ctx));
  PalacePetscCall(PCShellSetApply(pc, __pc_apply));
}

inline void ConfigureRG(RG rg, PetscReal lr, PetscReal ur, PetscReal li, PetscReal ui,
                        bool complement = false)
{
  PalacePetscCall(RGSetType(rg, RGINTERVAL));
  PalacePetscCall(RGIntervalSetEndpoints(rg, lr, ur, li, ui));
  if (complement)
  {
    PalacePetscCall(RGSetComplement(rg, PETSC_TRUE));
  }
}

}  // namespace

void Initialize(int &argc, char **&argv, const char rc_file[], const char help[])
{
  ConfigurePetscDevice();
  PalacePetscCall(SlepcInitialize(&argc, &argv, rc_file, help));

  // Remove default PETSc signal handling, since it can be confusing when the errors occur
  // not from within SLEPc/PETSc.
  PalacePetscCall(PetscPopSignalHandler());
}

void Initialize()
{
  ConfigurePetscDevice();
  PalacePetscCall(SlepcInitializeNoArguments());

  // Remove default PETSc signal handling, since it can be confusing when the errors occur
  // not from within SLEPc/PETSc.
  PalacePetscCall(PetscPopSignalHandler());
}

void Finalize()
{
  PalacePetscCall(SlepcFinalize());
}

PetscReal GetMaxSingularValue(MPI_Comm comm, const ComplexOperator &A, bool herm,
                              PetscReal tol, PetscInt max_it)
{
  // This method assumes the provided operator has the required operations for SLEPc's EPS
  // or SVD solvers, namely MATOP_MULT and MATOP_MULT_HERMITIAN_TRANSPOSE (if the matrix
  // is not Hermitian).
  MFEM_VERIFY(A.Height() == A.Width(), "Spectral norm requires a square matrix!");
  const PetscInt n = A.Height();
  ComplexVector x(n), y(n);
  x.UseDevice(true);
  y.UseDevice(true);
  MatShellContext ctx = {A, x, y};
  Mat A0;
  PalacePetscCall(
      MatCreateShell(comm, n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)&ctx, &A0));
  PalacePetscCall(MatShellSetOperation(A0, MATOP_MULT, (void (*)(void))__mat_apply_shell));
  PalacePetscCall(MatShellSetVecType(A0, PetscVecType()));
  if (herm)
  {
    EPS eps;
    PetscInt num_conv;
    PetscScalar eig;
    PalacePetscCall(EPSCreate(comm, &eps));
    PalacePetscCall(EPSSetOperators(eps, A0, nullptr));
    PalacePetscCall(EPSSetProblemType(eps, EPS_HEP));
    PalacePetscCall(EPSSetWhichEigenpairs(eps, EPS_LARGEST_MAGNITUDE));
    PalacePetscCall(EPSSetDimensions(eps, 1, PETSC_DEFAULT, PETSC_DEFAULT));
    PalacePetscCall(EPSSetTolerances(eps, tol, max_it));
    PalacePetscCall(EPSSolve(eps));
    PalacePetscCall(EPSGetConverged(eps, &num_conv));
    if (num_conv < 1)
    {
      Mpi::Warning(comm, "SLEPc EPS solve did not converge for maximum singular value!\n");
      eig = 0.0;
    }
    else
    {
      PalacePetscCall(EPSGetEigenvalue(eps, 0, &eig, nullptr));
      MFEM_VERIFY(PetscImaginaryPart(eig) == 0.0,
                  "Unexpected complex eigenvalue for Hermitian matrix (λ = " << eig
                                                                             << ")!");
    }
    PalacePetscCall(EPSDestroy(&eps));
    PalacePetscCall(MatDestroy(&A0));
    return PetscAbsScalar(eig);
  }
  else
  {
    PalacePetscCall(MatShellSetOperation(A0, MATOP_MULT_TRANSPOSE,
                                         (void (*)(void))__mat_apply_transpose_shell));
    PalacePetscCall(
        MatShellSetOperation(A0, MATOP_MULT_HERMITIAN_TRANSPOSE,
                             (void (*)(void))__mat_apply_hermitian_transpose_shell));

    SVD svd;
    PetscInt num_conv;
    PetscReal sigma;
    PalacePetscCall(SVDCreate(comm, &svd));
    PalacePetscCall(SVDSetOperators(svd, A0, nullptr));
    PalacePetscCall(SVDSetProblemType(svd, SVD_STANDARD));
    PalacePetscCall(SVDSetWhichSingularTriplets(svd, SVD_LARGEST));
    PalacePetscCall(SVDSetDimensions(svd, 1, PETSC_DEFAULT, PETSC_DEFAULT));
    PalacePetscCall(SVDSetTolerances(svd, tol, max_it));
    PalacePetscCall(SVDSolve(svd));
    PalacePetscCall(SVDGetConverged(svd, &num_conv));
    if (num_conv < 1)
    {
      Mpi::Warning(comm, "SLEPc SVD solve did not converge for maximum singular value!\n");
      sigma = 0.0;
    }
    else
    {
      PalacePetscCall(SVDGetSingularTriplet(svd, 0, &sigma, nullptr, nullptr));
    }
    PalacePetscCall(SVDDestroy(&svd));
    PalacePetscCall(MatDestroy(&A0));
    return sigma;
  }
}

// Eigensolver base class methods

SlepcEigenvalueSolver::SlepcEigenvalueSolver(int print) : print(print)
{
  sinvert = false;
  region = true;
  sigma = 0.0;
  gamma = delta = 1.0;

  opInv = nullptr;
  opProj = nullptr;
  opB = nullptr;

  B0 = nullptr;
  v0 = nullptr;

  cl_custom = false;
}

SlepcEigenvalueSolver::~SlepcEigenvalueSolver()
{
  PalacePetscCall(MatDestroy(&B0));
  PalacePetscCall(VecDestroy(&v0));
}

void SlepcEigenvalueSolver::SetOperators(const ComplexOperator &K, const ComplexOperator &M,
                                         EigenvalueSolver::ScaleType type)
{
  MFEM_ABORT("SetOperators not defined for base class SlepcEigenvalueSolver!");
}

void SlepcEigenvalueSolver::SetOperators(const ComplexOperator &K, const ComplexOperator &C,
                                         const ComplexOperator &M,
                                         EigenvalueSolver::ScaleType type)
{
  MFEM_ABORT("SetOperators not defined for base class SlepcEigenvalueSolver!");
}

void SlepcEigenvalueSolver::SetLinearSolver(const ComplexKspSolver &ksp)
{
  opInv = &ksp;
}

void SlepcEigenvalueSolver::SetDivFreeProjector(const DivFreeSolver<ComplexVector> &divfree)
{
  opProj = &divfree;
}

void SlepcEigenvalueSolver::SetBMat(const Operator &B)
{
  opB = &B;
}

void SlepcEigenvalueSolver::SetShiftInvert(std::complex<double> s, bool precond)
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
  sigma = s;  // Wait until solve time to call EPS/PEPSetTarget
  sinvert = true;
}

void SlepcEigenvalueSolver::SetOrthogonalization(bool mgs, bool cgs2)
{
  // The SLEPc default is CGS with refinement if needed.
  if (mgs || cgs2)
  {
    BV bv = GetBV();
    BVOrthogType type;
    BVOrthogRefineType refine;
    if (mgs)
    {
      type = BV_ORTHOG_MGS;
      refine = BV_ORTHOG_REFINE_NEVER;
    }
    else  // cgs2
    {
      type = BV_ORTHOG_CGS;
      refine = BV_ORTHOG_REFINE_ALWAYS;
    }
    PalacePetscCall(BVSetOrthogonalization(bv, type, refine, 1.0, BV_ORTHOG_BLOCK_GS));
  }
}

void SlepcEigenvalueSolver::Customize()
{
  // Configure the KSP object for non-preconditioned spectral transformations.
  PetscBool precond;
  ST st = GetST();
  PalacePetscCall(
      PetscObjectTypeCompare(reinterpret_cast<PetscObject>(st), STPRECOND, &precond));
  if (!precond)
  {
    KSP ksp;
    PalacePetscCall(STGetKSP(st, &ksp));
    PalacePetscCall(KSPSetType(ksp, KSPPREONLY));
  }

  // Configure the region based on the given target if necessary.
  if (sinvert && region)
  {
    if (PetscImaginaryPart(sigma) == 0.0)
    {
      PetscReal sr = PetscRealPart(sigma);
      if (sr > 0.0)
      {
        ConfigureRG(GetRG(), sr / gamma, mfem::infinity(), -mfem::infinity(),
                    mfem::infinity());
      }
      else if (sr < 0.0)
      {
        ConfigureRG(GetRG(), -mfem::infinity(), sr / gamma, -mfem::infinity(),
                    mfem::infinity());
      }
    }
    else if (PetscRealPart(sigma) == 0.0)
    {
      PetscReal si = PetscImaginaryPart(sigma);
      if (si > 0.0)
      {
        ConfigureRG(GetRG(), -mfem::infinity(), mfem::infinity(), si / gamma,
                    mfem::infinity());
      }
      else if (si < 0.0)
      {
        ConfigureRG(GetRG(), -mfem::infinity(), mfem::infinity(), -mfem::infinity(),
                    si / gamma);
      }
    }
    else
    {
      MFEM_ABORT("Shift-and-invert with general complex eigenvalue target is unsupported!");
    }
  }
}

PetscReal SlepcEigenvalueSolver::GetEigenvectorNorm(const ComplexVector &x,
                                                    ComplexVector &Bx) const
{
  if (opB)
  {
    return linalg::Norml2(GetComm(), x, *opB, Bx);
  }
  else
  {
    return linalg::Norml2(GetComm(), x);
  }
}

PetscReal SlepcEigenvalueSolver::GetError(int i, EigenvalueSolver::ErrorType type) const
{
  switch (type)
  {
    case ErrorType::ABSOLUTE:
      return res.get()[i];
    case ErrorType::RELATIVE:
      return res.get()[i] / PetscAbsScalar(GetEigenvalue(i));
    case ErrorType::BACKWARD:
      return res.get()[i] / GetBackwardScaling(GetEigenvalue(i));
  }
  return 0.0;
}

void SlepcEigenvalueSolver::RescaleEigenvectors(int num_eig)
{
  res = std::make_unique<PetscReal[]>(num_eig);
  xscale = std::make_unique<PetscReal[]>(num_eig);
  for (int i = 0; i < num_eig; i++)
  {
    xscale.get()[i] = 0.0;
    GetEigenvector(i, x1);
    xscale.get()[i] = 1.0 / GetEigenvectorNorm(x1, y1);
    res.get()[i] =
        GetResidualNorm(GetEigenvalue(i), x1, y1) / linalg::Norml2(GetComm(), x1);
  }
}

// EPS specific methods

SlepcEPSSolverBase::SlepcEPSSolverBase(MPI_Comm comm, int print, const std::string &prefix)
  : SlepcEigenvalueSolver(print)
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
  A0 = A1 = nullptr;
}

SlepcEPSSolverBase::~SlepcEPSSolverBase()
{
  PalacePetscCall(EPSDestroy(&eps));
  PalacePetscCall(MatDestroy(&A0));
  PalacePetscCall(MatDestroy(&A1));
}

void SlepcEPSSolverBase::SetNumModes(int num_eig, int num_vec)
{
  PalacePetscCall(EPSSetDimensions(eps, num_eig, (num_vec > 0) ? num_vec : PETSC_DEFAULT,
                                   PETSC_DEFAULT));
}

void SlepcEPSSolverBase::SetTol(PetscReal tol)
{
  PalacePetscCall(EPSSetTolerances(eps, tol, PETSC_DEFAULT));
  PalacePetscCall(EPSSetConvergenceTest(eps, EPS_CONV_REL));
  // PalacePetscCall(EPSSetTrackAll(eps, PETSC_TRUE));
  // PalacePetscCall(EPSSetTrueResidual(eps, PETSC_TRUE));
}

void SlepcEPSSolverBase::SetMaxIter(int max_it)
{
  PalacePetscCall(
      EPSSetTolerances(eps, PETSC_DEFAULT, (max_it > 0) ? max_it : PETSC_DEFAULT));
}

void SlepcEPSSolverBase::SetWhichEigenpairs(EigenvalueSolver::WhichType type)
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
  }
}

void SlepcEPSSolverBase::SetProblemType(SlepcEigenvalueSolver::ProblemType type)
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
    case ProblemType::HYPERBOLIC:
    case ProblemType::GYROSCOPIC:
      MFEM_ABORT("Problem type not implemented!");
      break;
  }
}

void SlepcEPSSolverBase::SetType(SlepcEigenvalueSolver::Type type)
{
  switch (type)
  {
    case Type::KRYLOVSCHUR:
      PalacePetscCall(EPSSetType(eps, EPSKRYLOVSCHUR));
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
    case Type::TOAR:
    case Type::STOAR:
    case Type::QARNOLDI:
      MFEM_ABORT("Eigenvalue solver type not implemented!");
      break;
  }
}

void SlepcEPSSolverBase::SetInitialSpace(const ComplexVector &v)
{
  MFEM_VERIFY(
      A0 && A1,
      "Must call SetOperators before using SetInitialSpace for SLEPc eigenvalue solver!");
  if (!v0)
  {
    PalacePetscCall(MatCreateVecs(A0, nullptr, &v0));
  }
  PalacePetscCall(ToPetscVec(v, v0));
  Vec is[1] = {v0};
  PalacePetscCall(EPSSetInitialSpace(eps, 1, is));
}

void SlepcEPSSolverBase::Customize()
{
  SlepcEigenvalueSolver::Customize();
  PalacePetscCall(EPSSetTarget(eps, sigma / gamma));
  if (!cl_custom)
  {
    PalacePetscCall(EPSSetFromOptions(eps));
    if (print > 0)
    {
      PetscOptionsView(nullptr, PETSC_VIEWER_STDOUT_(GetComm()));
      Mpi::Print(GetComm(), "\n");
    }
    cl_custom = true;
  }
}

int SlepcEPSSolverBase::Solve()
{
  MFEM_VERIFY(A0 && A1 && opInv, "Operators are not set for SlepcEPSSolverBase!");

  // Solve the eigenvalue problem.
  PetscInt num_conv;
  Customize();
  PalacePetscCall(EPSSolve(eps));
  PalacePetscCall(EPSGetConverged(eps, &num_conv));
  if (print > 0)
  {
    Mpi::Print(GetComm(), "\n");
    PalacePetscCall(EPSConvergedReasonView(eps, PETSC_VIEWER_STDOUT_(GetComm())));
    Mpi::Print(GetComm(),
               " Total number of linear systems solved: {:d}\n"
               " Total number of linear solver iterations: {:d}\n",
               opInv->NumTotalMult(), opInv->NumTotalMultIterations());
  }

  // Compute and store the eigenpair residuals.
  RescaleEigenvectors(num_conv);
  return (int)num_conv;
}

std::complex<double> SlepcEPSSolverBase::GetEigenvalue(int i) const
{
  PetscScalar l;
  PalacePetscCall(EPSGetEigenvalue(eps, i, &l, nullptr));
  return l * gamma;
}

void SlepcEPSSolverBase::GetEigenvector(int i, ComplexVector &x) const
{
  MFEM_VERIFY(
      v0,
      "Must call SetOperators before using GetEigenvector for SLEPc eigenvalue solver!");
  PalacePetscCall(EPSGetEigenvector(eps, i, v0, nullptr));
  PalacePetscCall(FromPetscVec(v0, x));
  if (xscale.get()[i] > 0.0)
  {
    x *= xscale.get()[i];
  }
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

SlepcEPSSolver::SlepcEPSSolver(MPI_Comm comm, int print, const std::string &prefix)
  : SlepcEPSSolverBase(comm, print, prefix)
{
  opK = opM = nullptr;
  normK = normM = 0.0;
}

void SlepcEPSSolver::SetOperators(const ComplexOperator &K, const ComplexOperator &M,
                                  EigenvalueSolver::ScaleType type)
{
  // Construct shell matrices for the scaled operators which define the generalized
  // eigenvalue problem.
  const bool first = (opK == nullptr);
  opK = &K;
  opM = &M;

  if (first)
  {
    const PetscInt n = opK->Height();
    PalacePetscCall(
        MatCreateShell(GetComm(), n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &A0));
    PalacePetscCall(
        MatCreateShell(GetComm(), n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &A1));
    PalacePetscCall(
        MatShellSetOperation(A0, MATOP_MULT, (void (*)(void))__mat_apply_EPS_A0));
    PalacePetscCall(
        MatShellSetOperation(A1, MATOP_MULT, (void (*)(void))__mat_apply_EPS_A1));
    PalacePetscCall(MatShellSetVecType(A0, PetscVecType()));
    PalacePetscCall(MatShellSetVecType(A1, PetscVecType()));
    PalacePetscCall(EPSSetOperators(eps, A0, A1));
  }

  if (first && type != ScaleType::NONE)
  {
    normK = linalg::SpectralNorm(GetComm(), *opK, opK->IsReal());
    normM = linalg::SpectralNorm(GetComm(), *opM, opM->IsReal());
    MFEM_VERIFY(normK >= 0.0 && normM >= 0.0, "Invalid matrix norms for EPS scaling!");
    if (normK > 0 && normM > 0.0)
    {
      gamma = normK / normM;  // Store γ² for linear problem
      delta = 2.0 / normK;
    }
  }

  // Set up workspace.
  if (!v0)
  {
    PalacePetscCall(MatCreateVecs(A0, nullptr, &v0));
  }
  x1.SetSize(opK->Height());
  y1.SetSize(opK->Height());
  x1.UseDevice(true);
  y1.UseDevice(true);

  // Configure linear solver for generalized problem or spectral transformation. This also
  // allows use of the divergence-free projector as a linear solve side-effect.
  if (first)
  {
    ConfigurePCShell(GetST(), (void *)this, __pc_apply_EPS);
  }
}

void SlepcEPSSolver::SetBMat(const Operator &B)
{
  SlepcEigenvalueSolver::SetBMat(B);

  const PetscInt n = B.Height();
  PalacePetscCall(
      MatCreateShell(GetComm(), n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &B0));
  PalacePetscCall(MatShellSetOperation(B0, MATOP_MULT, (void (*)(void))__mat_apply_EPS_B));
  PalacePetscCall(MatShellSetVecType(B0, PetscVecType()));

  BV bv = GetBV();
  PalacePetscCall(BVSetMatrix(bv, B0, PETSC_FALSE));
}

PetscReal SlepcEPSSolver::GetResidualNorm(PetscScalar l, const ComplexVector &x,
                                          ComplexVector &r) const
{
  // Compute the i-th eigenpair residual: || (K - λ M) x ||₂ for eigenvalue λ.
  opK->Mult(x, r);
  opM->AddMult(x, r, -l);
  return linalg::Norml2(GetComm(), r);
}

PetscReal SlepcEPSSolver::GetBackwardScaling(PetscScalar l) const
{
  // Make sure not to use norms from scaling as this can be confusing if they are different.
  // Note that SLEPc typically uses ||.||∞, not the 2-norm.
  if (normK <= 0.0)
  {
    normK = linalg::SpectralNorm(GetComm(), *opK, opK->IsReal());
  }
  if (normM <= 0.0)
  {
    normM = linalg::SpectralNorm(GetComm(), *opM, opM->IsReal());
  }
  return normK + PetscAbsScalar(l) * normM;
}

SlepcPEPLinearSolver::SlepcPEPLinearSolver(MPI_Comm comm, int print,
                                           const std::string &prefix)
  : SlepcEPSSolverBase(comm, print, prefix)
{
  opK = opC = opM = nullptr;
  normK = normC = normM = 0.0;
}

void SlepcPEPLinearSolver::SetOperators(const ComplexOperator &K, const ComplexOperator &C,
                                        const ComplexOperator &M,
                                        EigenvalueSolver::ScaleType type)
{
  // Construct shell matrices for the scaled linearized operators which define the block 2x2
  // eigenvalue problem.
  const bool first = (opK == nullptr);
  opK = &K;
  opC = &C;
  opM = &M;

  if (first)
  {
    const PetscInt n = opK->Height();
    PalacePetscCall(MatCreateShell(GetComm(), 2 * n, 2 * n, PETSC_DECIDE, PETSC_DECIDE,
                                   (void *)this, &A0));
    PalacePetscCall(MatCreateShell(GetComm(), 2 * n, 2 * n, PETSC_DECIDE, PETSC_DECIDE,
                                   (void *)this, &A1));
    PalacePetscCall(
        MatShellSetOperation(A0, MATOP_MULT, (void (*)(void))__mat_apply_PEPLinear_L0));
    PalacePetscCall(
        MatShellSetOperation(A1, MATOP_MULT, (void (*)(void))__mat_apply_PEPLinear_L1));
    PalacePetscCall(MatShellSetVecType(A0, PetscVecType()));
    PalacePetscCall(MatShellSetVecType(A1, PetscVecType()));
    PalacePetscCall(EPSSetOperators(eps, A0, A1));
  }

  if (first && type != ScaleType::NONE)
  {
    normK = linalg::SpectralNorm(GetComm(), *opK, opK->IsReal());
    normC = linalg::SpectralNorm(GetComm(), *opC, opC->IsReal());
    normM = linalg::SpectralNorm(GetComm(), *opM, opM->IsReal());
    MFEM_VERIFY(normK >= 0.0 && normC >= 0.0 && normM >= 0.0,
                "Invalid matrix norms for PEP scaling!");
    if (normK > 0 && normC > 0.0 && normM > 0.0)
    {
      gamma = std::sqrt(normK / normM);
      delta = 2.0 / (normK + gamma * normC);
    }
  }

  // Set up workspace.
  if (!v0)
  {
    PalacePetscCall(MatCreateVecs(A0, nullptr, &v0));
  }
  x1.SetSize(opK->Height());
  x2.SetSize(opK->Height());
  y1.SetSize(opK->Height());
  y2.SetSize(opK->Height());
  x1.UseDevice(true);
  x2.UseDevice(true);
  y1.UseDevice(true);
  y2.UseDevice(true);

  // Configure linear solver.
  if (first)
  {
    ConfigurePCShell(GetST(), (void *)this, __pc_apply_PEPLinear);
  }
}

void SlepcPEPLinearSolver::SetBMat(const Operator &B)
{
  SlepcEigenvalueSolver::SetBMat(B);

  const PetscInt n = B.Height();
  PalacePetscCall(MatCreateShell(GetComm(), 2 * n, 2 * n, PETSC_DECIDE, PETSC_DECIDE,
                                 (void *)this, &B0));
  PalacePetscCall(
      MatShellSetOperation(B0, MATOP_MULT, (void (*)(void))__mat_apply_PEPLinear_B));
  PalacePetscCall(MatShellSetVecType(B0, PetscVecType()));

  BV bv = GetBV();
  PalacePetscCall(BVSetMatrix(bv, B0, PETSC_FALSE));
}

void SlepcPEPLinearSolver::SetInitialSpace(const ComplexVector &v)
{
  MFEM_VERIFY(
      A0 && A1,
      "Must call SetOperators before using SetInitialSpace for SLEPc eigenvalue solver!");
  if (!v0)
  {
    PalacePetscCall(MatCreateVecs(A0, nullptr, &v0));
  }
  PalacePetscCall(ToPetscVec(v, v0, 0, 2));
  Vec is[1] = {v0};
  PalacePetscCall(EPSSetInitialSpace(eps, 1, is));
}

void SlepcPEPLinearSolver::GetEigenvector(int i, ComplexVector &x) const
{
  // Select the most accurate x for y = [x₁; x₂] from the linearized eigenvalue problem. Or,
  // just take x = x₁.
  MFEM_VERIFY(
      v0,
      "Must call SetOperators before using GetEigenvector for SLEPc eigenvalue solver!");
  PalacePetscCall(EPSGetEigenvector(eps, i, v0, nullptr));
  PalacePetscCall(FromPetscVec(v0, x, 0, 2));
  if (xscale.get()[i] > 0.0)
  {
    x *= xscale.get()[i];
  }
}

PetscReal SlepcPEPLinearSolver::GetResidualNorm(PetscScalar l, const ComplexVector &x,
                                                ComplexVector &r) const
{
  // Compute the i-th eigenpair residual: || P(λ) x ||₂ = || (K + λ C + λ² M) x ||₂ for
  // eigenvalue λ.
  opK->Mult(x, r);
  opC->AddMult(x, r, l);
  opM->AddMult(x, r, l * l);
  return linalg::Norml2(GetComm(), r);
}

PetscReal SlepcPEPLinearSolver::GetBackwardScaling(PetscScalar l) const
{
  // Make sure not to use norms from scaling as this can be confusing if they are different.
  // Note that SLEPc typically uses ||.||∞, not the 2-norm.
  if (normK <= 0.0)
  {
    normK = linalg::SpectralNorm(GetComm(), *opK, opK->IsReal());
  }
  if (normC <= 0.0)
  {
    normC = linalg::SpectralNorm(GetComm(), *opC, opC->IsReal());
  }
  if (normM <= 0.0)
  {
    normM = linalg::SpectralNorm(GetComm(), *opM, opM->IsReal());
  }
  PetscReal t = PetscAbsScalar(l);
  return normK + t * normC + t * t * normM;
}

// PEP specific methods

SlepcPEPSolverBase::SlepcPEPSolverBase(MPI_Comm comm, int print, const std::string &prefix)
  : SlepcEigenvalueSolver(print)
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
  PalacePetscCall(PEPDestroy(&pep));
  PalacePetscCall(MatDestroy(&A0));
  PalacePetscCall(MatDestroy(&A1));
  PalacePetscCall(MatDestroy(&A2));
}

void SlepcPEPSolverBase::SetNumModes(int num_eig, int num_vec)
{
  PalacePetscCall(PEPSetDimensions(pep, num_eig, (num_vec > 0) ? num_vec : PETSC_DEFAULT,
                                   PETSC_DEFAULT));
}

void SlepcPEPSolverBase::SetTol(PetscReal tol)
{
  PalacePetscCall(PEPSetTolerances(pep, tol, PETSC_DEFAULT));
  PalacePetscCall(PEPSetConvergenceTest(pep, PEP_CONV_REL));
  // PalacePetscCall(PEPSetTrackAll(pep, PETSC_TRUE));
}

void SlepcPEPSolverBase::SetMaxIter(int max_it)
{
  PalacePetscCall(
      PEPSetTolerances(pep, PETSC_DEFAULT, (max_it > 0) ? max_it : PETSC_DEFAULT));
}

void SlepcPEPSolverBase::SetWhichEigenpairs(EigenvalueSolver::WhichType type)
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
  }
}

void SlepcPEPSolverBase::SetProblemType(SlepcEigenvalueSolver::ProblemType type)
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
  }
}

void SlepcPEPSolverBase::SetType(SlepcEigenvalueSolver::Type type)
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
    case Type::KRYLOVSCHUR:
    case Type::POWER:
    case Type::SUBSPACE:
      MFEM_ABORT("Eigenvalue solver type not implemented!");
      break;
  }
}

void SlepcPEPSolverBase::SetInitialSpace(const ComplexVector &v)
{
  MFEM_VERIFY(
      A0 && A1 && A2,
      "Must call SetOperators before using SetInitialSpace for SLEPc eigenvalue solver!");
  if (!v0)
  {
    PalacePetscCall(MatCreateVecs(A0, nullptr, &v0));
  }
  PalacePetscCall(ToPetscVec(v, v0));
  Vec is[1] = {v0};
  PalacePetscCall(PEPSetInitialSpace(pep, 1, is));
}

void SlepcPEPSolverBase::Customize()
{
  SlepcEigenvalueSolver::Customize();
  PalacePetscCall(PEPSetTarget(pep, sigma / gamma));
  if (!cl_custom)
  {
    PalacePetscCall(PEPSetFromOptions(pep));
    if (print > 0)
    {
      PetscOptionsView(nullptr, PETSC_VIEWER_STDOUT_(GetComm()));
      Mpi::Print(GetComm(), "\n");
    }
    cl_custom = true;
  }
}

int SlepcPEPSolverBase::Solve()
{
  MFEM_VERIFY(A0 && A1 && A2 && opInv, "Operators are not set for SlepcPEPSolverBase!");

  // Solve the eigenvalue problem.
  PetscInt num_conv;
  Customize();
  PalacePetscCall(PEPSolve(pep));
  PalacePetscCall(PEPGetConverged(pep, &num_conv));
  if (print > 0)
  {
    Mpi::Print(GetComm(), "\n");
    PalacePetscCall(PEPConvergedReasonView(pep, PETSC_VIEWER_STDOUT_(GetComm())));
    Mpi::Print(GetComm(),
               " Total number of linear systems solved: {:d}\n"
               " Total number of linear solver iterations: {:d}\n",
               opInv->NumTotalMult(), opInv->NumTotalMultIterations());
  }

  // Compute and store the eigenpair residuals.
  RescaleEigenvectors(num_conv);
  return (int)num_conv;
}

std::complex<double> SlepcPEPSolverBase::GetEigenvalue(int i) const
{
  PetscScalar l;
  PalacePetscCall(PEPGetEigenpair(pep, i, &l, nullptr, nullptr, nullptr));
  return l * gamma;
}

void SlepcPEPSolverBase::GetEigenvector(int i, ComplexVector &x) const
{
  MFEM_VERIFY(
      v0,
      "Must call SetOperators before using GetEigenvector for SLEPc eigenvalue solver!");
  PalacePetscCall(PEPGetEigenpair(pep, i, nullptr, nullptr, v0, nullptr));
  PalacePetscCall(FromPetscVec(v0, x));
  if (xscale.get()[i] > 0.0)
  {
    x *= xscale.get()[i];
  }
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

SlepcPEPSolver::SlepcPEPSolver(MPI_Comm comm, int print, const std::string &prefix)
  : SlepcPEPSolverBase(comm, print, prefix)
{
  opK = opC = opM = nullptr;
  normK = normC = normM = 0.0;
}

void SlepcPEPSolver::SetOperators(const ComplexOperator &K, const ComplexOperator &C,
                                  const ComplexOperator &M,
                                  EigenvalueSolver::ScaleType type)
{
  // Construct shell matrices for the scaled operators which define the quadratic polynomial
  // eigenvalue problem.
  const bool first = (opK == nullptr);
  opK = &K;
  opC = &C;
  opM = &M;

  if (first)
  {
    const PetscInt n = opK->Height();
    PalacePetscCall(
        MatCreateShell(GetComm(), n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &A0));
    PalacePetscCall(
        MatCreateShell(GetComm(), n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &A1));
    PalacePetscCall(
        MatCreateShell(GetComm(), n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &A2));
    PalacePetscCall(
        MatShellSetOperation(A0, MATOP_MULT, (void (*)(void))__mat_apply_PEP_A0));
    PalacePetscCall(
        MatShellSetOperation(A1, MATOP_MULT, (void (*)(void))__mat_apply_PEP_A1));
    PalacePetscCall(
        MatShellSetOperation(A2, MATOP_MULT, (void (*)(void))__mat_apply_PEP_A2));
    PalacePetscCall(MatShellSetVecType(A0, PetscVecType()));
    PalacePetscCall(MatShellSetVecType(A1, PetscVecType()));
    PalacePetscCall(MatShellSetVecType(A2, PetscVecType()));
    Mat A[3] = {A0, A1, A2};
    PalacePetscCall(PEPSetOperators(pep, 3, A));
  }

  if (first && type != ScaleType::NONE)
  {
    normK = linalg::SpectralNorm(GetComm(), *opK, opK->IsReal());
    normC = linalg::SpectralNorm(GetComm(), *opC, opC->IsReal());
    normM = linalg::SpectralNorm(GetComm(), *opM, opM->IsReal());
    MFEM_VERIFY(normK >= 0.0 && normC >= 0.0 && normM >= 0.0,
                "Invalid matrix norms for PEP scaling!");
    if (normK > 0 && normC > 0.0 && normM > 0.0)
    {
      gamma = std::sqrt(normK / normM);
      delta = 2.0 / (normK + gamma * normC);
    }
  }

  // Set up workspace.
  if (!v0)
  {
    PalacePetscCall(MatCreateVecs(A0, nullptr, &v0));
  }
  x1.SetSize(opK->Height());
  y1.SetSize(opK->Height());

  // Configure linear solver.
  if (first)
  {
    ConfigurePCShell(GetST(), (void *)this, __pc_apply_PEP);
  }
}

void SlepcPEPSolver::SetBMat(const Operator &B)
{
  SlepcEigenvalueSolver::SetBMat(B);

  const PetscInt n = B.Height();
  PalacePetscCall(
      MatCreateShell(GetComm(), n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &B0));
  PalacePetscCall(MatShellSetOperation(B0, MATOP_MULT, (void (*)(void))__mat_apply_PEP_B));
  PalacePetscCall(MatShellSetVecType(B0, PetscVecType()));

  BV bv = GetBV();
  PalacePetscCall(BVSetMatrix(bv, B0, PETSC_FALSE));
}

PetscReal SlepcPEPSolver::GetResidualNorm(PetscScalar l, const ComplexVector &x,
                                          ComplexVector &r) const
{
  // Compute the i-th eigenpair residual: || P(λ) x ||₂ = || (K + λ C + λ² M) x ||₂ for
  // eigenvalue λ.
  opK->Mult(x, r);
  opC->AddMult(x, r, l);
  opM->AddMult(x, r, l * l);
  return linalg::Norml2(GetComm(), r);
}

PetscReal SlepcPEPSolver::GetBackwardScaling(PetscScalar l) const
{
  // Make sure not to use norms from scaling as this can be confusing if they are different.
  // Note that SLEPc typically uses ||.||∞, not Frobenius.
  if (normK <= 0.0)
  {
    normK = linalg::SpectralNorm(GetComm(), *opK, opK->IsReal());
  }
  if (normC <= 0.0)
  {
    normC = linalg::SpectralNorm(GetComm(), *opC, opC->IsReal());
  }
  if (normM <= 0.0)
  {
    normM = linalg::SpectralNorm(GetComm(), *opM, opM->IsReal());
  }
  PetscReal t = PetscAbsScalar(l);
  return normK + t * normC + t * t * normM;
}

}  // namespace palace::slepc

PetscErrorCode __mat_apply_EPS_A0(Mat A, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  palace::slepc::SlepcEPSSolver *ctx;
  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context for SLEPc!");

  PetscCall(FromPetscVec(x, ctx->x1));
  ctx->opK->Mult(ctx->x1, ctx->y1);
  ctx->y1 *= ctx->delta;
  PetscCall(ToPetscVec(ctx->y1, y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __mat_apply_EPS_A1(Mat A, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  palace::slepc::SlepcEPSSolver *ctx;
  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context for SLEPc!");

  PetscCall(FromPetscVec(x, ctx->x1));
  ctx->opM->Mult(ctx->x1, ctx->y1);
  ctx->y1 *= ctx->delta * ctx->gamma;
  PetscCall(ToPetscVec(ctx->y1, y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __mat_apply_EPS_B(Mat A, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  palace::slepc::SlepcEPSSolver *ctx;
  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context for SLEPc!");

  PetscCall(FromPetscVec(x, ctx->x1));
  ctx->opB->Mult(ctx->x1.Real(), ctx->y1.Real());
  ctx->opB->Mult(ctx->x1.Imag(), ctx->y1.Imag());
  ctx->y1 *= ctx->delta * ctx->gamma;
  PetscCall(ToPetscVec(ctx->y1, y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __pc_apply_EPS(PC pc, Vec x, Vec y)
{
  // Solve the linear system associated with the generalized eigenvalue problem: y =
  // M⁻¹ x, or shift-and-invert spectral transformation: y = (K - σ M)⁻¹ x . Enforces the
  // divergence-free constraint using the supplied projector.
  PetscFunctionBeginUser;
  palace::slepc::SlepcEPSSolver *ctx;
  PetscCall(PCShellGetContext(pc, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell PC context for SLEPc!");

  PetscCall(FromPetscVec(x, ctx->x1));
  ctx->opInv->Mult(ctx->x1, ctx->y1);
  if (!ctx->sinvert)
  {
    ctx->y1 *= 1.0 / (ctx->delta * ctx->gamma);
  }
  else
  {
    ctx->y1 *= 1.0 / ctx->delta;
  }
  if (ctx->opProj)
  {
    // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(ctx->GetComm(), ctx->y1));
    ctx->opProj->Mult(ctx->y1);
    // Mpi::Print(" After projection: {:e}\n", linalg::Norml2(ctx->GetComm(), ctx->y1));
  }
  PetscCall(ToPetscVec(ctx->y1, y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __mat_apply_PEPLinear_L0(Mat A, Vec x, Vec y)
{
  // Apply the linearized operator L₀ = [  0  I ]
  //                                    [ -K -C ] .
  PetscFunctionBeginUser;
  palace::slepc::SlepcPEPLinearSolver *ctx;
  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context for SLEPc!");
  PetscCall(FromPetscVec(x, ctx->x1, ctx->x2));
  ctx->y1 = ctx->x2;
  ctx->opC->Mult(ctx->x2, ctx->y2);
  ctx->y2 *= ctx->gamma;
  ctx->opK->AddMult(ctx->x1, ctx->y2, std::complex<double>(1.0, 0.0));
  ctx->y2 *= -ctx->delta;
  PetscCall(ToPetscVec(ctx->y1, ctx->y2, y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __mat_apply_PEPLinear_L1(Mat A, Vec x, Vec y)
{
  // Apply the linearized operator L₁ = [ I  0 ]
  //                                    [ 0  M ] .
  PetscFunctionBeginUser;
  palace::slepc::SlepcPEPLinearSolver *ctx;
  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context for SLEPc!");

  PetscCall(FromPetscVec(x, ctx->x1, ctx->x2));
  ctx->y1 = ctx->x1;
  ctx->opM->Mult(ctx->x2, ctx->y2);
  ctx->y2 *= ctx->delta * ctx->gamma * ctx->gamma;
  PetscCall(ToPetscVec(ctx->y1, ctx->y2, y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __mat_apply_PEPLinear_B(Mat A, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  palace::slepc::SlepcPEPLinearSolver *ctx;
  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context for SLEPc!");

  PetscCall(FromPetscVec(x, ctx->x1, ctx->x2));
  ctx->opB->Mult(ctx->x1.Real(), ctx->y1.Real());
  ctx->opB->Mult(ctx->x1.Imag(), ctx->y1.Imag());
  ctx->opB->Mult(ctx->x2.Real(), ctx->y2.Real());
  ctx->opB->Mult(ctx->x2.Imag(), ctx->y2.Imag());
  ctx->y1 *= ctx->delta * ctx->gamma * ctx->gamma;
  ctx->y2 *= ctx->delta * ctx->gamma * ctx->gamma;
  PetscCall(ToPetscVec(ctx->y1, ctx->y2, y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __pc_apply_PEPLinear(PC pc, Vec x, Vec y)
{
  // Solve the linear system associated with the generalized eigenvalue problem after
  // linearization: y = L₁⁻¹ x, or with the shift-and-invert spectral transformation:
  // y = (L₀ - σ L₁)⁻¹ x, with:
  //               L₀ = [  0  I ]    L₁ = [ I  0 ]
  //                    [ -K -C ] ,       [ 0  M ] .
  // Enforces the divergence-free constraint using the supplied projector.
  PetscFunctionBeginUser;
  palace::slepc::SlepcPEPLinearSolver *ctx;
  PetscCall(PCShellGetContext(pc, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell PC context for SLEPc!");

  PetscCall(FromPetscVec(x, ctx->x1, ctx->x2));
  if (!ctx->sinvert)
  {
    ctx->y1 = ctx->x1;
    if (ctx->opProj)
    {
      // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(ctx->GetComm(), ctx->y1));
      ctx->opProj->Mult(ctx->y1);
      // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(ctx->GetComm(), ctx->y1));
    }

    ctx->opInv->Mult(ctx->x2, ctx->y2);
    ctx->y2 *= 1.0 / (ctx->delta * ctx->gamma * ctx->gamma);
    if (ctx->opProj)
    {
      // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(ctx->GetComm(), ctx->y2));
      ctx->opProj->Mult(ctx->y2);
      // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(ctx->GetComm(), ctx->y2));
    }
  }
  else
  {
    ctx->y1.AXPBY(-ctx->sigma / (ctx->delta * ctx->gamma), ctx->x2, 0.0);  // Temporarily
    ctx->opK->AddMult(ctx->x1, ctx->y1, std::complex<double>(1.0, 0.0));
    ctx->opInv->Mult(ctx->y1, ctx->y2);
    if (ctx->opProj)
    {
      // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(ctx->GetComm(), ctx->y2));
      ctx->opProj->Mult(ctx->y2);
      // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(ctx->GetComm(), ctx->y2));
    }

    ctx->y1.AXPBYPCZ(ctx->gamma / ctx->sigma, ctx->y2, -ctx->gamma / ctx->sigma, ctx->x1,
                     0.0);
    if (ctx->opProj)
    {
      // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(ctx->GetComm(), ctx->y1));
      ctx->opProj->Mult(ctx->y1);
      // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(ctx->GetComm(), ctx->y1));
    }
  }
  PetscCall(ToPetscVec(ctx->y1, ctx->y2, y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __mat_apply_PEP_A0(Mat A, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  palace::slepc::SlepcPEPSolver *ctx;
  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context for SLEPc!");

  PetscCall(FromPetscVec(x, ctx->x1));
  ctx->opK->Mult(ctx->x1, ctx->y1);
  PetscCall(ToPetscVec(ctx->y1, y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __mat_apply_PEP_A1(Mat A, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  palace::slepc::SlepcPEPSolver *ctx;
  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context for SLEPc!");

  PetscCall(FromPetscVec(x, ctx->x1));
  ctx->opC->Mult(ctx->x1, ctx->y1);
  PetscCall(ToPetscVec(ctx->y1, y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __mat_apply_PEP_A2(Mat A, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  palace::slepc::SlepcPEPSolver *ctx;
  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context for SLEPc!");

  PetscCall(FromPetscVec(x, ctx->x1));
  ctx->opM->Mult(ctx->x1, ctx->y1);
  PetscCall(ToPetscVec(ctx->y1, y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __mat_apply_PEP_B(Mat A, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  palace::slepc::SlepcPEPSolver *ctx;
  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context for SLEPc!");

  PetscCall(FromPetscVec(x, ctx->x1));
  ctx->opB->Mult(ctx->x1.Real(), ctx->y1.Real());
  ctx->opB->Mult(ctx->x1.Imag(), ctx->y1.Imag());
  ctx->y1 *= ctx->delta * ctx->gamma;
  PetscCall(ToPetscVec(ctx->y1, y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __pc_apply_PEP(PC pc, Vec x, Vec y)
{
  // Solve the linear system associated with the generalized eigenvalue problem: y = M⁻¹ x,
  // or shift-and-invert spectral transformation: y = P(σ)⁻¹ x . Enforces the divergence-
  // free constraint using the supplied projector.
  PetscFunctionBeginUser;
  palace::slepc::SlepcPEPSolver *ctx;
  PetscCall(PCShellGetContext(pc, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell PC context for SLEPc!");

  PetscCall(FromPetscVec(x, ctx->x1));
  ctx->opInv->Mult(ctx->x1, ctx->y1);
  if (!ctx->sinvert)
  {
    ctx->y1 *= 1.0 / (ctx->delta * ctx->gamma * ctx->gamma);
  }
  else
  {
    ctx->y1 *= 1.0 / ctx->delta;
  }
  if (ctx->opProj)
  {
    // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(ctx->GetComm(), ctx->y1));
    ctx->opProj->Mult(ctx->y1);
    // Mpi::Print(" After projection: {:e}\n", linalg::Norml2(ctx->GetComm(), ctx->y1));
  }
  PetscCall(ToPetscVec(ctx->y1, y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

#endif
