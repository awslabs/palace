// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "slepc.hpp"

#if defined(PALACE_WITH_SLEPC)

#include <algorithm>
#include <petsc.h>
#include <slepc.h>
#include <mfem.hpp>
#include "linalg/divfree.hpp"
#include "linalg/nleps.hpp"
#include "linalg/rap.hpp"
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
// for NEP
static PetscErrorCode __mat_apply_NEP_A(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_NEP_J(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_NEP_B(Mat, Vec, Vec);
static PetscErrorCode __pc_apply_NEP(PC, Vec, Vec);
static PetscErrorCode __mat_apply_NLEIGS_term(Mat, Vec, Vec);
static PetscErrorCode __mat_duplicate_NLEIGS_term(Mat, MatDuplicateOption, Mat *);
static PetscErrorCode __mat_destroy_NLEIGS_term(Mat);
static PetscErrorCode __nleigs_singularities_callback(NEP, PetscInt *, PetscScalar *,
                                                      void *);
static PetscErrorCode __pc_apply_NLEIGS(PC, Vec, Vec);
static PetscErrorCode __form_NEP_function(NEP, PetscScalar, Mat, Mat, void *);
static PetscErrorCode __form_NEP_jacobian(NEP, PetscScalar, Mat, void *);

using namespace std::complex_literals;

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

// Eigensolver base class methods.

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

void SlepcEigenvalueSolver::SetLinearSolver(ComplexKspSolver &ksp)
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

// EPS specific methods.

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
  PalacePetscCall(EPSSetDimensions(eps, num_eig, num_vec, PETSC_DEFAULT));
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
  PalacePetscCall(EPSSetTolerances(eps, PETSC_DEFAULT, max_it));
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
    case ProblemType::GENERAL:
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
    case Type::SLP:
    case Type::NLEIGS:
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
    if (normK > 0 && normC >= 0.0 && normM > 0.0)
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
  if (opC)
  {
    opC->AddMult(x, r, l);
  }
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
  if (normC <= 0.0 && opC)
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

// PEP specific methods.

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
  PalacePetscCall(PEPSetDimensions(pep, num_eig, num_vec, PETSC_DEFAULT));
}

void SlepcPEPSolverBase::SetTol(PetscReal tol)
{
  PalacePetscCall(PEPSetTolerances(pep, tol, PETSC_DEFAULT));
  PalacePetscCall(PEPSetConvergenceTest(pep, PEP_CONV_REL));
  // PalacePetscCall(PEPSetTrackAll(pep, PETSC_TRUE));
}

void SlepcPEPSolverBase::SetMaxIter(int max_it)
{
  PalacePetscCall(PEPSetTolerances(pep, PETSC_DEFAULT, max_it));
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
    case ProblemType::GENERAL:
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
    case Type::SLP:
    case Type::NLEIGS:
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
    if (normK > 0 && normC >= 0.0 && normM > 0.0)
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
  if (opC)
  {
    opC->AddMult(x, r, l);
  }
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
  if (normC <= 0.0 && opC)
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

// NEP specific methods.

SlepcNEPSolverBase::SlepcNEPSolverBase(MPI_Comm comm, int print, const std::string &prefix)
  : SlepcEigenvalueSolver(print)
{
  PalacePetscCall(NEPCreate(comm, &nep));
  PalacePetscCall(NEPSetOptionsPrefix(nep, prefix.c_str()));
  if (print > 0)
  {
    std::string opts = "-nep_monitor";
    if (print > 2)
    {
      opts.append(" -nep_view");
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
  A = J = nullptr;
}

SlepcNEPSolverBase::~SlepcNEPSolverBase()
{
  PalacePetscCall(NEPDestroy(&nep));
  PalacePetscCall(MatDestroy(&A));
  PalacePetscCall(MatDestroy(&J));
}

void SlepcNEPSolverBase::SetNumModes(int num_eig, int num_vec)
{
  PalacePetscCall(NEPSetDimensions(nep, num_eig, num_vec, PETSC_DEFAULT));
}

void SlepcNEPSolverBase::SetTol(PetscReal tol)
{
  PalacePetscCall(NEPSetTolerances(nep, tol, PETSC_DEFAULT));
  PalacePetscCall(NEPSetConvergenceTest(nep, NEP_CONV_REL));
}

void SlepcNEPSolverBase::SetMaxIter(int max_it)
{
  PalacePetscCall(NEPSetTolerances(nep, PETSC_DEFAULT, max_it));
}

void SlepcNEPSolverBase::SetShiftInvert(std::complex<double> s, bool precond)
{
  sigma = s;  // Wait until solve time to call NEPSetTarget
  sinvert = false;
}

void SlepcNEPSolverBase::SetWhichEigenpairs(EigenvalueSolver::WhichType type)
{
  switch (type)
  {
    case WhichType::LARGEST_MAGNITUDE:
      PalacePetscCall(NEPSetWhichEigenpairs(nep, NEP_LARGEST_MAGNITUDE));
      region = false;
      break;
    case WhichType::SMALLEST_MAGNITUDE:
      PalacePetscCall(NEPSetWhichEigenpairs(nep, NEP_SMALLEST_MAGNITUDE));
      region = false;
      break;
    case WhichType::LARGEST_REAL:
      PalacePetscCall(NEPSetWhichEigenpairs(nep, NEP_LARGEST_REAL));
      break;
    case WhichType::SMALLEST_REAL:
      PalacePetscCall(NEPSetWhichEigenpairs(nep, NEP_SMALLEST_REAL));
      break;
    case WhichType::LARGEST_IMAGINARY:
      PalacePetscCall(NEPSetWhichEigenpairs(nep, NEP_LARGEST_IMAGINARY));
      break;
    case WhichType::SMALLEST_IMAGINARY:
      PalacePetscCall(NEPSetWhichEigenpairs(nep, NEP_SMALLEST_IMAGINARY));
      break;
    case WhichType::TARGET_MAGNITUDE:
      PalacePetscCall(NEPSetWhichEigenpairs(nep, NEP_TARGET_MAGNITUDE));
      region = false;
      break;
    case WhichType::TARGET_REAL:
      PalacePetscCall(NEPSetWhichEigenpairs(nep, NEP_TARGET_REAL));
      break;
    case WhichType::TARGET_IMAGINARY:
      PalacePetscCall(NEPSetWhichEigenpairs(nep, NEP_TARGET_IMAGINARY));
      break;
  }
}

void SlepcNEPSolverBase::SetProblemType(SlepcEigenvalueSolver::ProblemType type)
{
  switch (type)
  {
    case ProblemType::GENERAL:
      PalacePetscCall(NEPSetProblemType(nep, NEP_GENERAL));
      break;
    case ProblemType::HERMITIAN:
    case ProblemType::GEN_HERMITIAN:
    case ProblemType::NON_HERMITIAN:
    case ProblemType::GEN_INDEFINITE:
    case ProblemType::GEN_NON_HERMITIAN:
    case ProblemType::HYPERBOLIC:
    case ProblemType::GYROSCOPIC:
      MFEM_ABORT("Problem type not implemented!");
      break;
  }
}

void SlepcNEPSolverBase::SetType(SlepcEigenvalueSolver::Type type)
{
  switch (type)
  {
    case Type::SLP:
      PalacePetscCall(NEPSetType(nep, NEPSLP));
      break;
    case Type::NLEIGS:
      PalacePetscCall(NEPSetType(nep, NEPNLEIGS));
      break;
    case Type::KRYLOVSCHUR:
    case Type::POWER:
    case Type::SUBSPACE:
    case Type::JD:
    case Type::TOAR:
    case Type::STOAR:
    case Type::QARNOLDI:
      MFEM_ABORT("Eigenvalue solver type not implemented!");
      break;
  }
}

void SlepcNEPSolverBase::SetInitialSpace(const ComplexVector &v)
{
  MFEM_VERIFY(
      A && J,
      "Must call SetOperators before using SetInitialSpace for SLEPc eigenvalue solver!");
  if (!v0)
  {
    PalacePetscCall(MatCreateVecs(A, nullptr, &v0));
  }
  PalacePetscCall(ToPetscVec(v, v0));
  Vec is[1] = {v0};
  PalacePetscCall(NEPSetInitialSpace(nep, 1, is));
}

void SlepcNEPSolverBase::Customize()
{
  // Configure the region based on the given target if necessary.
  PalacePetscCall(NEPSetTarget(nep, sigma));
  if (!cl_custom)
  {
    PalacePetscCall(NEPSetFromOptions(nep));
    if (print > 0)
    {
      PetscOptionsView(nullptr, PETSC_VIEWER_STDOUT_(GetComm()));
      Mpi::Print(GetComm(), "\n");
    }
    cl_custom = true;
  }
}

int SlepcNEPSolverBase::Solve()
{
  MFEM_VERIFY(A && J && opInv, "Operators are not set for SlepcNEPSolverBase!");

  // Solve the eigenvalue problem.
  perm.reset();
  PetscInt num_conv;
  Customize();
  PalacePetscCall(NEPSolve(nep));
  PalacePetscCall(NEPGetConverged(nep, &num_conv));
  if (print > 0)
  {
    Mpi::Print(GetComm(), "\n");
    PalacePetscCall(NEPConvergedReasonView(nep, PETSC_VIEWER_STDOUT_(GetComm())));
    Mpi::Print(GetComm(),
               " Total number of linear systems solved: {:d}\n"
               " Total number of linear solver iterations: {:d}\n",
               opInv->NumTotalMult(), opInv->NumTotalMultIterations());
  }

  // Compute and store the ordered eigenpair residuals.
  const int nev = (int)num_conv;
  perm = std::make_unique<int[]>(nev);
  std::vector<std::complex<double>> eig(nev);
  for (int i = 0; i < nev; i++)
  {
    PetscScalar l;
    PalacePetscCall(NEPGetEigenpair(nep, i, &l, nullptr, nullptr, nullptr));
    eig[i] = l;
    perm[i] = i;
  }
  // Sort by ascending imaginary component.
  std::sort(perm.get(), perm.get() + nev,
            [&eig](auto l, auto r) { return eig[l].imag() < eig[r].imag(); });
  RescaleEigenvectors(nev);
  return nev;
}

std::complex<double> SlepcNEPSolverBase::GetEigenvalue(int i) const
{
  PetscScalar l;
  const int &j = perm.get()[i];
  PalacePetscCall(NEPGetEigenpair(nep, j, &l, nullptr, nullptr, nullptr));
  return l;
}

void SlepcNEPSolverBase::GetEigenvector(int i, ComplexVector &x) const
{
  MFEM_VERIFY(
      v0,
      "Must call SetOperators before using GetEigenvector for SLEPc eigenvalue solver!");
  const int &j = perm.get()[i];
  PalacePetscCall(NEPGetEigenpair(nep, j, nullptr, nullptr, v0, nullptr));
  PalacePetscCall(FromPetscVec(v0, x));
  if (xscale.get()[i] > 0.0)
  {
    x *= xscale.get()[i];
  }
}

BV SlepcNEPSolverBase::GetBV() const
{
  BV bv;
  PalacePetscCall(NEPGetBV(nep, &bv));
  return bv;
}

ST SlepcNEPSolverBase::GetST() const
{
  ST st = nullptr;
  // NEPGetST does not exist.
  return st;
}

RG SlepcNEPSolverBase::GetRG() const
{
  RG rg;
  PalacePetscCall(NEPGetRG(nep, &rg));
  return rg;
}

SlepcNEPSolver::SlepcNEPSolver(MPI_Comm comm, int print, const std::string &prefix)
  : SlepcNEPSolverBase(comm, print, prefix)
{
  opK = opC = opM = nullptr;
  normK = normC = normM = 0.0;
}

void SlepcNEPSolver::SetExtraSystemMatrix(
    std::function<std::unique_ptr<ComplexOperator>(double)> A2)
{
  funcA2 = A2;
}

std::unique_ptr<ComplexOperator> SlepcNEPSolver::BuildA2(std::complex<double> lambda) const
{
  if (use_complex_a2 && funcA2Complex)
  {
    return (*funcA2Complex)(lambda);
  }
  return (*funcA2)(std::abs(lambda.imag()));
}

void SlepcNEPSolver::SetExtraSystemMatrixComplex(
    std::function<std::unique_ptr<ComplexOperator>(std::complex<double>)> A2c)
{
  funcA2Complex = A2c;
}

void SlepcNEPSolver::SetPreconditionerUpdate(
    std::function<
        std::unique_ptr<ComplexOperator>(std::complex<double>, std::complex<double>,
                                         std::complex<double>, std::complex<double>)>
        P)
{
  funcP = P;
}

void SlepcNEPSolver::SetOperators(const ComplexOperator &K, const ComplexOperator &M,
                                  EigenvalueSolver::ScaleType type)
{
  // Construct shell matrices for the scaled operators which define the quadratic polynomial
  // eigenvalue problem.
  const bool first = (opK == nullptr);
  opK = &K;
  opM = &M;

  if (first)
  {
    const PetscInt n = opK->Height();
    PalacePetscCall(
        MatCreateShell(GetComm(), n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &A));
    PalacePetscCall(
        MatCreateShell(GetComm(), n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &J));
    PalacePetscCall(MatShellSetOperation(A, MATOP_MULT, (void (*)(void))__mat_apply_NEP_A));
    PalacePetscCall(MatShellSetOperation(J, MATOP_MULT, (void (*)(void))__mat_apply_NEP_J));
    PalacePetscCall(MatShellSetVecType(A, PetscVecType()));
    PalacePetscCall(MatShellSetVecType(J, PetscVecType()));
    PalacePetscCall(NEPSetFunction(nep, A, A, __form_NEP_function, NULL));
    PalacePetscCall(NEPSetJacobian(nep, J, __form_NEP_jacobian, NULL));
  }

  if (first && type != ScaleType::NONE)
  {
    normK = linalg::SpectralNorm(GetComm(), *opK, opK->IsReal());
    normM = linalg::SpectralNorm(GetComm(), *opM, opM->IsReal());
    MFEM_VERIFY(normK >= 0.0 && normM >= 0.0, "Invalid matrix norms for NEP scaling!");
    if (normK > 0 && normM > 0.0)
    {
      gamma = std::sqrt(normK / normM);
      delta = 2.0 / normK;
    }
  }

  // Set up workspace.
  if (!v0)
  {
    PalacePetscCall(MatCreateVecs(A, nullptr, &v0));
  }
  x1.SetSize(opK->Height());
  y1.SetSize(opK->Height());

  // Configure linear solver.
  if (first)
  {
    // SLP.
    PC pc;
    KSP ksp;
    EPS eps;
    PalacePetscCall(NEPSLPGetKSP(nep, &ksp));
    PalacePetscCall(KSPSetType(ksp, KSPPREONLY));
    PalacePetscCall(NEPSLPGetEPS(nep, &eps));
    PalacePetscCall(EPSSetType(eps, EPSKRYLOVSCHUR));
    PalacePetscCall(KSPGetPC(ksp, &pc));
    PalacePetscCall(PCSetType(pc, PCSHELL));
    PalacePetscCall(PCShellSetContext(pc, (void *)this));
    PalacePetscCall(PCShellSetApply(pc, __pc_apply_NEP));
  }
}

void SlepcNEPSolver::SetOperators(const ComplexOperator &K, const ComplexOperator &C,
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
        MatCreateShell(GetComm(), n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &A));
    PalacePetscCall(
        MatCreateShell(GetComm(), n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &J));
    PalacePetscCall(MatShellSetOperation(A, MATOP_MULT, (void (*)(void))__mat_apply_NEP_A));
    PalacePetscCall(MatShellSetOperation(J, MATOP_MULT, (void (*)(void))__mat_apply_NEP_J));
    PalacePetscCall(MatShellSetVecType(A, PetscVecType()));
    PalacePetscCall(MatShellSetVecType(J, PetscVecType()));
    PalacePetscCall(NEPSetFunction(nep, A, A, __form_NEP_function, NULL));
    PalacePetscCall(NEPSetJacobian(nep, J, __form_NEP_jacobian, NULL));
  }

  if (first && type != ScaleType::NONE)
  {
    normK = linalg::SpectralNorm(GetComm(), *opK, opK->IsReal());
    normC = linalg::SpectralNorm(GetComm(), *opC, opC->IsReal());
    normM = linalg::SpectralNorm(GetComm(), *opM, opM->IsReal());
    MFEM_VERIFY(normK >= 0.0 && normC >= 0.0 && normM >= 0.0,
                "Invalid matrix norms for NEP scaling!");
    if (normK > 0 && normC >= 0.0 && normM > 0.0)
    {
      gamma = std::sqrt(normK / normM);
      delta = 2.0 / (normK + gamma * normC);
    }
  }

  // Set up workspace.
  if (!v0)
  {
    PalacePetscCall(MatCreateVecs(A, nullptr, &v0));
  }
  x1.SetSize(opK->Height());
  y1.SetSize(opK->Height());

  // Configure linear solver.
  if (first)
  {
    // SLP.
    PC pc;
    KSP ksp;
    EPS eps;
    PalacePetscCall(NEPSLPGetKSP(nep, &ksp));
    PalacePetscCall(KSPSetType(ksp, KSPPREONLY));
    PalacePetscCall(NEPSLPGetEPS(nep, &eps));
    PalacePetscCall(EPSSetType(eps, EPSKRYLOVSCHUR));
    PalacePetscCall(KSPGetPC(ksp, &pc));
    PalacePetscCall(PCSetType(pc, PCSHELL));
    PalacePetscCall(PCShellSetContext(pc, (void *)this));
    PalacePetscCall(PCShellSetApply(pc, __pc_apply_NEP));
  }
}

void SlepcNEPSolver::SetBMat(const Operator &B)
{
  SlepcEigenvalueSolver::SetBMat(B);

  const PetscInt n = B.Height();
  PalacePetscCall(
      MatCreateShell(GetComm(), n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &B0));
  PalacePetscCall(MatShellSetOperation(B0, MATOP_MULT, (void (*)(void))__mat_apply_NEP_B));
  PalacePetscCall(MatShellSetVecType(B0, PetscVecType()));

  BV bv = GetBV();
  PalacePetscCall(BVSetMatrix(bv, B0, PETSC_FALSE));
}

PetscReal SlepcNEPSolver::GetResidualNorm(PetscScalar l, const ComplexVector &x,
                                          ComplexVector &r) const
{
  // Compute the i-th eigenpair residual: || P(λ) x ||₂ = || (K + λ C + λ² M) x ||₂ for
  // eigenvalue λ.
  opK->Mult(x, r);
  if (opC)
  {
    opC->AddMult(x, r, l);
  }
  opM->AddMult(x, r, l * l);
  if (funcA2 || funcA2Complex)
  {
    auto A2 = BuildA2(l);
    A2->AddMult(x, r, 1.0 + 0.0i);
  }
  return linalg::Norml2(GetComm(), r);
}

PetscReal SlepcNEPSolver::GetBackwardScaling(PetscScalar l) const
{
  // Make sure not to use norms from scaling as this can be confusing if they are different.
  // Note that SLEPc typically uses ||.||∞, not Frobenius.
  if (normK <= 0.0)
  {
    normK = linalg::SpectralNorm(GetComm(), *opK, opK->IsReal());
  }
  if (normC <= 0.0 && opC)
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

// SLEPc NEP NLEIGS split-form solver implementation. See header for design notes.

struct SingularityCtx
{
  std::vector<PetscScalar> values;
};

SlepcNEPNLEIGSSolver::SlepcNEPNLEIGSSolver(MPI_Comm comm, int print,
                                           const std::string &prefix)
  : SlepcNEPSolverBase(comm, print, prefix)
{
}

SlepcNEPNLEIGSSolver::~SlepcNEPNLEIGSSolver()
{
  // Order matters: tear down NEP first so it releases its references to the Mats and
  // FNs, then destroy the local Mats / FNs / contexts. The base destructor will
  // double-call NEPDestroy harmlessly (sets to nullptr).
  if (nep)
  {
    PalacePetscCall(NEPDestroy(&nep));
  }
  for (auto &m : term_mats_)
  {
    PalacePetscCall(MatDestroy(&m));
  }
  term_mats_.clear();
  for (auto &f : term_fns_)
  {
    PalacePetscCall(FNDestroy(&f));
  }
  term_fns_.clear();
  // term_ctxs_ destructed automatically.
}

void SlepcNEPNLEIGSSolver::SetOperators(const ComplexOperator &K, const ComplexOperator &C,
                                        const ComplexOperator &M, ScaleType type)
{
  // Reserve the first three split-form slots for the polynomial pencil:
  //   K · 1, C · λ, M · λ². Extra A2(λ) terms (wave-port, farfield, surface
  // conductivity) are appended via AddSplitFormTerm before Solve.
  opK_ref = &K;
  opC_ref = &C;
  opM_ref = &M;
  problem_size_ = K.Height();
  x1.SetSize(problem_size_);
  y1.SetSize(problem_size_);
  x1.UseDevice(true);
  y1.UseDevice(true);

  if (type != ScaleType::NONE)
  {
    normK = linalg::SpectralNorm(GetComm(), K, K.IsReal());
    normC = linalg::SpectralNorm(GetComm(), C, C.IsReal());
    normM = linalg::SpectralNorm(GetComm(), M, M.IsReal());
    if (normK > 0.0 && normM > 0.0)
    {
      gamma = std::sqrt(normK / normM);
      delta = 2.0 / (normK + gamma * normC);
    }
  }

  // Build FNs for the polynomial pencil: f₀(λ) = 1, f₁(λ) = λ, f₂(λ) = λ². FN_RATIONAL
  // stores polynomial coefficients in *highest-degree-first* order: [c_n,...,c_1,c_0]
  // (verified from SLEPc fnrational.c Horner loop).
  AddSplitFormTerm(opK_ref,
                   [&]()
                   {
                     FN fn;
                     PalacePetscCall(FNCreate(GetComm(), &fn));
                     PalacePetscCall(FNSetType(fn, FNRATIONAL));
                     PetscScalar p[1] = {1.0 + 0.0i};
                     PalacePetscCall(FNRationalSetNumerator(fn, 1, p));
                     return fn;
                   }());
  AddSplitFormTerm(opC_ref,
                   [&]()
                   {
                     FN fn;
                     PalacePetscCall(FNCreate(GetComm(), &fn));
                     PalacePetscCall(FNSetType(fn, FNRATIONAL));
                     PetscScalar p[2] = {1.0 + 0.0i, 0.0 + 0.0i};  // 1·λ + 0
                     PalacePetscCall(FNRationalSetNumerator(fn, 2, p));
                     return fn;
                   }());
  AddSplitFormTerm(opM_ref,
                   [&]()
                   {
                     FN fn;
                     PalacePetscCall(FNCreate(GetComm(), &fn));
                     PalacePetscCall(FNSetType(fn, FNRATIONAL));
                     PetscScalar p[3] = {1.0 + 0.0i, 0.0 + 0.0i,
                                         0.0 + 0.0i};  // 1·λ² + 0·λ + 0
                     PalacePetscCall(FNRationalSetNumerator(fn, 3, p));
                     return fn;
                   }());
}

void SlepcNEPNLEIGSSolver::SetBMat(const Operator &B)
{
  opB = &B;
}

void SlepcNEPNLEIGSSolver::AddSplitFormTerm(const ComplexOperator *op, FN fn)
{
  MFEM_VERIFY(op, "SlepcNEPNLEIGSSolver::AddSplitFormTerm: null operator!");
  MFEM_VERIFY(problem_size_ > 0, "Must call SetOperators before AddSplitFormTerm.");

  // Heap-allocate the per-term context so MatShellGetContext returns a stable pointer
  // throughout the NEP solve.
  auto ctx = std::make_unique<TermCtx>();
  ctx->op = op;
  ctx->x.SetSize(problem_size_);
  ctx->y.SetSize(problem_size_);
  ctx->x.UseDevice(true);
  ctx->y.UseDevice(true);

  Mat shell;
  PalacePetscCall(MatCreateShell(GetComm(), problem_size_, problem_size_, PETSC_DECIDE,
                                 PETSC_DECIDE, ctx.get(), &shell));
  PalacePetscCall(
      MatShellSetOperation(shell, MATOP_MULT, (void (*)(void))__mat_apply_NLEIGS_term));
  PalacePetscCall(MatShellSetOperation(shell, MATOP_DUPLICATE,
                                       (void (*)(void))__mat_duplicate_NLEIGS_term));
  PalacePetscCall(MatShellSetOperation(shell, MATOP_DESTROY,
                                       (void (*)(void))__mat_destroy_NLEIGS_term));
  PalacePetscCall(MatShellSetVecType(shell, PetscVecType()));

  term_ctxs_.push_back(std::move(ctx));
  term_mats_.push_back(shell);
  term_fns_.push_back(fn);
}

void SlepcNEPNLEIGSSolver::Customize()
{
  MFEM_VERIFY(!term_mats_.empty(), "SlepcNEPNLEIGSSolver: no split-form terms registered!");

  // Register the split-form operator with NLEIGS. The MatStructure flag is
  // DIFFERENT_NONZERO_PATTERN because shell matrices have no nonzero pattern; SLEPc
  // would conservatively assume this anyway.
  PalacePetscCall(NEPSetSplitOperator(nep, static_cast<PetscInt>(term_mats_.size()),
                                      term_mats_.data(), term_fns_.data(),
                                      DIFFERENT_NONZERO_PATTERN));

  // NLEIGS-specific tuning. SLEPc forbids full basis together with rational-Krylov
  // shifts (NEPNLEIGSSetRKShifts), so multiple shifts force full basis off. Note: on
  // Palace's split-form pencil full_basis = false is known to give poor / non-
  // converging results, so rk_shifts_ > 1 is effectively unusable — kept only for
  // experimentation.
  const bool full_basis = full_basis_ && (rk_shifts_ <= 1);
  if (full_basis_ && rk_shifts_ > 1)
  {
    Mpi::Warning(GetComm(), "NLEIGS: NLEIGSRKShifts > 1 is incompatible with full basis; "
                            "disabling full basis (expect degraded convergence).\n");
  }
  PalacePetscCall(NEPNLEIGSSetFullBasis(nep, full_basis ? PETSC_TRUE : PETSC_FALSE));
  PalacePetscCall(NEPNLEIGSSetInterpolation(nep, interp_tol_, interp_deg_));

  // Singularity set Ξ. SLEPc has two paths: (a) if NO singularities callback is
  // installed and the problem type is NEP_GENERAL, SLEPc runs its AAA algorithm to
  // discover the singularities automatically by sampling T(λ) on the target contour;
  // (b) explicit list via NEPNLEIGSSetSingularitiesFunction. When the caller provides
  // a singularity set (singularities_ non-empty — e.g. the wave-port branch points),
  // we install the callback to pass it; otherwise we install NOTHING and let SLEPc's
  // AAA auto-discovery run. (Leentvaar/Leja–Bagby rejects singularities at the origin,
  // so the farfield 1/(2λ) pole and TEM kₙ branch at λ=0 must not be passed verbatim.)
  static SingularityCtx ctx;  // file-scope singleton — only one NLEIGS solver runs
                              // at a time per Palace invocation.
  if (!singularities_.empty())
  {
    ctx.values.assign(singularities_.begin(), singularities_.end());
    PalacePetscCall(
        NEPNLEIGSSetSingularitiesFunction(nep, __nleigs_singularities_callback, &ctx));
    if (print > 0)
    {
      Mpi::Print(GetComm(), " NLEIGS: {} explicit singularit{} installed\n",
                 ctx.values.size(), ctx.values.size() == 1 ? "y" : "ies");
    }
  }
  else if (print > 0)
  {
    Mpi::Print(GetComm(), " NLEIGS: no singularities set — SLEPc AAA auto-discovery\n");
  }

  // Target region Σ. If user supplied a 4-tuple, build RGINTERVAL with those bounds;
  // otherwise auto-derive from sigma at solve time. We require sigma to have been
  // set via SetShiftInvert before Customize is called.
  RG rg;
  PalacePetscCall(NEPGetRG(nep, &rg));
  PalacePetscCall(RGSetType(rg, RGINTERVAL));
  if (region_.size() == 4)
  {
    PalacePetscCall(
        RGIntervalSetEndpoints(rg, region_[0], region_[1], region_[2], region_[3]));
  }
  else
  {
    // NLEIGS places Leja-Bagby interpolation nodes ON the RG boundary, and the
    // rational interpolant's eigenvalues at the boundary appear as "converged"
    // spurious modes. We bracket the band [target, target_upper] tightly along
    // the imaginary axis: the lower bound is the target itself (NO padding below).
    // For non-TEM wave ports, the cutoff frequency lies just below the band of
    // interest; extending the RG below the user's target risks placing
    // interpolation nodes near or below cutoff, where kₙ(λ) becomes purely real
    // (evanescent) and the rational fit struggles. The user is responsible for
    // selecting Target above the relevant cutoff frequency.
    //
    // Real-axis bandwidth ±target_upper allows substantial damping for low-Q
    // modes. Users override via NLEIGSRegion.
    const PetscReal target_im = sigma.imag();
    const PetscReal upper_im = (target_upper_ > 0.0) ? target_upper_ : 3.0 * target_im;
    MFEM_VERIFY(upper_im > target_im,
                "SlepcNEPNLEIGSSolver: NLEIGS region requires target_upper > target. "
                "Set TargetUpper in config or NLEIGSRegion explicitly.");
    const PetscReal lower_im = target_im;
    const PetscReal damp = upper_im;
    PalacePetscCall(RGIntervalSetEndpoints(rg, -damp, damp, lower_im, upper_im));
  }

  // Set target.
  PalacePetscCall(NEPSetTarget(nep, sigma));

  // Rational-Krylov shifts. With a single shift (default) NLEIGS reduces to
  // shift-and-invert Krylov-Schur at the target, which tends to pile converged Ritz
  // values onto one boundary point of the RG (the source of the spurious near-target
  // cluster for near-cutoff non-TEM ports). Spreading rk_shifts_ shifts across the
  // imaginary-axis band [target, target_upper] runs the rational-Krylov variant, which
  // resolves the interior modes instead of collapsing onto the target. Shifts are pure
  // imaginary (λ = i·ω), matching the physical mode axis; all inner KSPs remain routed
  // through the single target preconditioner (opInv) — the RK iteration tolerates the
  // inexact per-shift PC.
  if (rk_shifts_ > 1)
  {
    const PetscReal target_im = sigma.imag();
    const PetscReal upper_im = (target_upper_ > 0.0) ? target_upper_ : 3.0 * target_im;
    std::vector<PetscScalar> shifts(rk_shifts_);
    for (int i = 0; i < rk_shifts_; i++)
    {
      const double frac =
          (rk_shifts_ == 1) ? 0.0 : static_cast<double>(i) / (rk_shifts_ - 1);
      shifts[i] = std::complex<double>(0.0, target_im + frac * (upper_im - target_im));
    }
    PalacePetscCall(
        NEPNLEIGSSetRKShifts(nep, static_cast<PetscInt>(shifts.size()), shifts.data()));
    if (print > 0)
    {
      Mpi::Print(GetComm(), " NLEIGS rational-Krylov shifts: {} on Im[{:.3e}, {:.3e}]\n",
                 rk_shifts_, target_im, upper_im);
    }
  }

  // NLEIGS' inner KSPs default to KSPGMRES + PCLU which requires a sparse direct
  // factor (MatGetFactor) — Palace builds PETSc without MUMPS/SuperLU. Override
  // each inner KSP to use KSPPREONLY + PCSHELL routed through Palace's own
  // ComplexKspSolver (opInv). NLEIGSGetKSPs LAZILY ALLOCATES the inner KSPs on
  // first call; doing so here BEFORE NEPSetUp lets NLEIGS see our configured KSPs
  // when its divided-differences routine subsequently calls KSPSetUp.
  PetscInt nshifts = 0;
  KSP *ksps = nullptr;
  PalacePetscCall(NEPNLEIGSGetKSPs(nep, &nshifts, &ksps));
  ksp_ctxs_.clear();
  ksp_ctxs_.reserve(nshifts);
  for (PetscInt i = 0; i < nshifts; i++)
  {
    auto ctx = std::make_unique<KSPCtx>();
    ctx->solver = this;
    // The actual rational interpolation shift for this inner KSP is set by NLEIGS
    // after divided differences completes. Defer the shift extraction to the first
    // PC apply by reading KSPGetOperators(ksp, &T, &P) — T is the assembled
    // T(σ_i) at that point. For the simpler approach taken here, we use the NEP
    // target as a single-shift approximation; NLEIGS' Krylov method will then
    // iterate to convergence regardless.
    ctx->sigma_shift = sigma;
    PalacePetscCall(KSPSetType(ksps[i], KSPPREONLY));
    PC pc;
    PalacePetscCall(KSPGetPC(ksps[i], &pc));
    PalacePetscCall(PCSetType(pc, PCSHELL));
    PalacePetscCall(PCShellSetContext(pc, (void *)ctx.get()));
    PalacePetscCall(PCShellSetApply(pc, __pc_apply_NLEIGS));
    ksp_ctxs_.push_back(std::move(ctx));
  }
  if (print > 0)
  {
    Mpi::Print(GetComm(), " NLEIGS inner KSPs: {} shifts, PCSHELL → opInv\n", nshifts);
  }
  // Convergence test: NEP_CONV_REL (= ||r||/|λ|). NEP_CONV_NORM requires a matrix-norm
  // operation our split-form shell matrices do not provide. Spurious RG-boundary modes
  // (Leja–Bagby interpolation-node artifacts) can pass this test with large true
  // ||T(λ)x||; the reported Error (Abs.) column exposes them so they can be filtered.
  PalacePetscCall(NEPSetConvergenceTest(nep, NEP_CONV_REL));
  // Skip NEPSetFromOptions for NLEIGS — it would touch the inner-KSP options and
  // interfere with our explicit PCSHELL configuration above. The user can still set
  // -nep_* options via the command-line interface; those that don't conflict will
  // apply via NEPSolve's internal setup paths.
  cl_custom = true;
}

std::unique_ptr<ComplexOperator>
SlepcNEPNLEIGSSolver::BuildTAtShift(std::complex<double> sigma_i) const
{
  // T(σ_i) = K + σ·C + σ²·M + Σ_{j≥3} f_j(σ_i)·A_j. Assembled via BuildParSumOperator
  // over the same constant matrices that NLEIGS uses internally.
  std::vector<std::complex<double>> coeffs;
  std::vector<const ComplexParOperator *> ops;
  for (std::size_t k = 0; k < term_ctxs_.size(); k++)
  {
    PetscScalar f_val;
    PalacePetscCall(FNEvaluateFunction(term_fns_[k], sigma_i, &f_val));
    coeffs.push_back(f_val);
    ops.push_back(static_cast<const ComplexParOperator *>(term_ctxs_[k]->op));
  }
  return BuildParSumOperator(coeffs, ops);
}

std::unique_ptr<ComplexOperator>
SlepcNEPNLEIGSSolver::BuildPAtShift(std::complex<double> sigma_i) const
{
  // The funcP signature is funcP(a0, a1, a2, ω) building Palace's preconditioner for
  // the operator a0·K + a1·C + a2·M + (extra at frequency ω). Match the polynomial
  // coefficients to T(σ_i) = K + σ·C + σ²·M + ...; pass ω = |Im σ_i| so the bulk
  // contribution is evaluated at the right frequency.
  MFEM_VERIFY(funcP_, "NLEIGS BuildPAtShift: no preconditioner builder set!");
  // NLEIGS uses the closed-form (fit) wave-port BC, evaluated on the real-ω axis, so the
  // preconditioner BC frequency is the real |Im σ_i| (passed as a real-valued complex).
  return (*funcP_)(std::complex<double>(1.0, 0.0), sigma_i, sigma_i * sigma_i,
                   std::complex<double>(std::abs(sigma_i.imag()), 0.0));
}

PetscReal SlepcNEPNLEIGSSolver::GetResidualNorm(PetscScalar l, const ComplexVector &x,
                                                ComplexVector &r) const
{
  // T(l) x = (K + l·C + l²·M + Σ extra A2 terms) x. We evaluate the full split-form
  // residual by re-applying the per-term shells; SLEPc could provide this via
  // NEPComputeFunction but going through term_ctxs_ keeps the path explicit.
  MFEM_VERIFY(opK_ref && opM_ref, "Operators not set!");
  opK_ref->Mult(x, r);
  if (opC_ref)
  {
    opC_ref->AddMult(x, r, l);
  }
  opM_ref->AddMult(x, r, l * l);
  // Extra A2(λ) terms: indices 3..end of term_mats_ (slots 0,1,2 are K,C,M handled
  // above; the FN values for those are 1, λ, λ²).
  for (std::size_t k = 3; k < term_ctxs_.size(); k++)
  {
    PetscScalar f_val;
    PalacePetscCall(FNEvaluateFunction(term_fns_[k], l, &f_val));
    term_ctxs_[k]->op->AddMult(x, r, f_val);
  }
  return linalg::Norml2(GetComm(), r);
}

PetscReal SlepcNEPNLEIGSSolver::GetBackwardScaling(PetscScalar l) const
{
  // Backward scaling for split-form NLEIGS:
  //   denom = ||K|| + |λ|·||C|| + |λ|²·||M|| + Σ_{j>2} |f_j(λ)| · ||A_j||
  // The Σ over BC terms (wave-port, farfield) is critical because those operators
  // can dominate near the eigenvalue band — without them backward error is
  // artificially small (or undefined when polynomial part is small).
  if (normK <= 0.0 && opK_ref)
  {
    normK = linalg::SpectralNorm(GetComm(), *opK_ref, opK_ref->IsReal());
  }
  if (normC <= 0.0 && opC_ref)
  {
    normC = linalg::SpectralNorm(GetComm(), *opC_ref, opC_ref->IsReal());
  }
  if (normM <= 0.0 && opM_ref)
  {
    normM = linalg::SpectralNorm(GetComm(), *opM_ref, opM_ref->IsReal());
  }
  if (term_norms_.size() != term_ctxs_.size())
  {
    term_norms_.assign(term_ctxs_.size(), -1.0);
  }
  const PetscReal t = PetscAbsScalar(l);
  PetscReal denom = normK + t * normC + t * t * normM;
  // Add per-term A2 contributions (slots 3..end of term_ctxs_; slots 0,1,2 are K,C,M
  // already counted in normK/C/M above).
  for (std::size_t k = 3; k < term_ctxs_.size(); k++)
  {
    if (term_norms_[k] < 0.0)
    {
      term_norms_[k] =
          linalg::SpectralNorm(GetComm(), *term_ctxs_[k]->op, term_ctxs_[k]->op->IsReal());
    }
    PetscScalar f_val;
    PalacePetscCall(FNEvaluateFunction(term_fns_[k], l, &f_val));
    denom += PetscAbsScalar(f_val) * term_norms_[k];
  }
  return denom;
}

void SlepcNEPNLEIGSSolver::SetInitialSpace(const ComplexVector &v)
{
  MFEM_VERIFY(!term_mats_.empty(),
              "SlepcNEPNLEIGSSolver::SetInitialSpace: must register split-form terms "
              "first via SetOperators!");
  if (!v0)
  {
    // Use MatCreateVecs on one of our shell mats so the resulting Vec layout
    // matches what NEP/NLEIGS expects internally (size, distribution, type).
    PalacePetscCall(MatCreateVecs(term_mats_.front(), nullptr, &v0));
  }
  PalacePetscCall(ToPetscVec(v, v0));
  Vec is[1] = {v0};
  PalacePetscCall(NEPSetInitialSpace(nep, 1, is));
}

int SlepcNEPNLEIGSSolver::Solve()
{
  MFEM_VERIFY(opInv, "Linear solver not set for SlepcNEPNLEIGSSolver!");
  MFEM_VERIFY(!term_mats_.empty(),
              "SlepcNEPNLEIGSSolver::Solve: no split-form terms registered!");

  // SLEPc 3.21+ requires NEPSetSplitOperator to be called before NEPSetFromOptions
  // touches NLEIGS-specific options. Customize handles both in the right order.
  perm.reset();
  PetscInt num_conv;
  Customize();
  PalacePetscCall(NEPSolve(nep));
  PalacePetscCall(NEPGetConverged(nep, &num_conv));
  if (print > 0)
  {
    Mpi::Print(GetComm(), "\n");
    PalacePetscCall(NEPConvergedReasonView(nep, PETSC_VIEWER_STDOUT_(GetComm())));
    Mpi::Print(GetComm(),
               " Total number of linear systems solved: {:d}\n"
               " Total number of linear solver iterations: {:d}\n",
               opInv->NumTotalMult(), opInv->NumTotalMultIterations());
  }

  // Sort by ascending imaginary component, mirroring SlepcNEPSolverBase::Solve.
  const int nev = static_cast<int>(num_conv);
  perm = std::make_unique<int[]>(nev);
  std::vector<std::complex<double>> eig(nev);
  for (int i = 0; i < nev; i++)
  {
    PetscScalar l;
    PalacePetscCall(NEPGetEigenpair(nep, i, &l, nullptr, nullptr, nullptr));
    eig[i] = l;
    perm[i] = i;
  }
  std::sort(perm.get(), perm.get() + nev,
            [&eig](auto l, auto r) { return eig[l].imag() < eig[r].imag(); });
  RescaleEigenvectors(nev);
  return nev;
}

}  // namespace palace::slepc

PetscErrorCode __mat_apply_NLEIGS_term(Mat A, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  palace::slepc::SlepcNEPNLEIGSSolver::TermCtx *ctx;
  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx && ctx->op, "Invalid PETSc shell matrix context for NLEIGS term!");
  PetscCall(FromPetscVec(x, ctx->x));
  ctx->op->Mult(ctx->x, ctx->y);
  PetscCall(ToPetscVec(ctx->y, y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __mat_duplicate_NLEIGS_term(Mat A, MatDuplicateOption /*op*/, Mat *B)
{
  // Shallow duplicate for SLEPc NLEIGS' internal NEP_NLEIGS_MATSHELL wrapper. The
  // underlying ComplexOperator is immutable during the NLEIGS solve, so the
  // duplicate can point at the same TermCtx. The caller never modifies values; the
  // duplicate is used as one entry in a (matrix, scalar) sum tracked by SLEPc's
  // own shell wrapper. The Destroy operation is a no-op (the original TermCtx is
  // owned by SlepcNEPNLEIGSSolver and outlives any duplicate).
  PetscFunctionBeginUser;
  palace::slepc::SlepcNEPNLEIGSSolver::TermCtx *ctx;
  PetscCall(MatShellGetContext(A, (void **)&ctx));
  PetscInt m, n, M, N;
  PetscCall(MatGetLocalSize(A, &m, &n));
  PetscCall(MatGetSize(A, &M, &N));
  PetscCall(MatCreateShell(PetscObjectComm((PetscObject)A), m, n, M, N, (void *)ctx, B));
  PetscCall(MatShellSetOperation(*B, MATOP_MULT, (void (*)(void))__mat_apply_NLEIGS_term));
  PetscCall(MatShellSetOperation(*B, MATOP_DUPLICATE,
                                 (void (*)(void))__mat_duplicate_NLEIGS_term));
  PetscCall(
      MatShellSetOperation(*B, MATOP_DESTROY, (void (*)(void))__mat_destroy_NLEIGS_term));
  PetscCall(MatShellSetVecType(*B, palace::slepc::PetscVecType()));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __mat_destroy_NLEIGS_term(Mat /*A*/)
{
  // No-op destroy: the TermCtx pointed to by MatShellGetContext is owned by the
  // SlepcNEPNLEIGSSolver (heap-allocated unique_ptr in term_ctxs_). Both the
  // original split-form mats AND any duplicates SLEPc creates inside its NEP_NLEIGS
  // wrapper share the same TermCtx via shallow duplication, so we must NOT free it
  // here. The shell Mat itself is freed by the caller (PETSc) via the standard
  // PetscObject destruction path.
  PetscFunctionBeginUser;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __pc_apply_NLEIGS(PC pc, Vec x, Vec y)
{
  // NLEIGS inner KSPs apply T(σ_i)⁻¹ approximately at each rational interpolation
  // shift σ_i. Palace builds opInv ONCE upstream (eigensolver driver) with the
  // assembled A = K + iσ·C − σ²·M + A2(σ) at the eigenmode target σ; reusing that
  // single opInv as the preconditioner for ALL shifts is a tractable approximation
  // — NLEIGS' Krylov iteration tolerates an inexact PC. Reconfiguring opInv per
  // shift would require rebuilding the multigrid preconditioner with shift-σ_i
  // coarse-grid factorizations, which Palace's MfemWrapperSolver→SuperLU coarse
  // path doesn't support for shell-based operators.
  //
  // Net effect: the inner KSPs are PREONLY → PCSHELL → opInv->Mult(x, y). Krylov
  // dimension grows with NLEIGS' rational basis size, so even a fixed-σ PC
  // converges the outer problem.
  PetscFunctionBeginUser;
  palace::slepc::SlepcNEPNLEIGSSolver::KSPCtx *ctx;
  PetscCall(PCShellGetContext(pc, (void **)&ctx));
  MFEM_VERIFY(ctx && ctx->solver, "Invalid PETSc shell PC context for NLEIGS!");
  auto *solver = ctx->solver;
  auto *opInv = solver->GetOpInv();
  MFEM_VERIFY(opInv, "NLEIGS PC apply: solver-level KSP not set!");
  PetscInt vlen;
  PetscCall(VecGetLocalSize(x, &vlen));
  palace::ComplexVector xv, yv;
  xv.SetSize(vlen);
  yv.SetSize(vlen);
  xv.UseDevice(true);
  yv.UseDevice(true);
  PetscCall(FromPetscVec(x, xv));
  opInv->Mult(xv, yv);
  if (auto *opProj = solver->GetOpProj())
  {
    opProj->Mult(yv);
  }
  PetscCall(ToPetscVec(yv, y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __nleigs_singularities_callback(NEP, PetscInt *maxnp, PetscScalar *xi,
                                               void *vctx)
{
  PetscFunctionBeginUser;
  auto *ctx = static_cast<palace::slepc::SingularityCtx *>(vctx);
  // SLEPc passes *maxnp as the buffer capacity on input; write min(values.size(),
  // capacity) and update *maxnp. Truncation emits a warning; Leja-Bagby still picks
  // interpolation nodes in that case.
  PetscInt cap = *maxnp;
  PetscInt n = std::min<PetscInt>(static_cast<PetscInt>(ctx->values.size()), cap);
  for (PetscInt i = 0; i < n; i++)
  {
    xi[i] = ctx->values[i];
  }
  if (static_cast<PetscInt>(ctx->values.size()) > cap)
  {
    palace::Mpi::Warning(
        " NLEIGS singularity callback: {} singularities provided but SLEPc buffer "
        "capacity is {} — truncating\n",
        ctx->values.size(), cap);
  }
  *maxnp = n;
  PetscFunctionReturn(PETSC_SUCCESS);
}

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
  if (ctx->opC)
  {
    ctx->opC->Mult(ctx->x2, ctx->y2);
  }
  else
  {
    ctx->y2 = 0.0;
  }
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
  if (ctx->opC)
  {
    ctx->opC->Mult(ctx->x1, ctx->y1);
  }
  else
  {
    ctx->y1 = 0.0;
  }
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

PetscErrorCode __mat_apply_NEP_A(Mat A, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  palace::slepc::SlepcNEPSolver *ctx;
  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context for SLEPc!");
  PetscCall(FromPetscVec(x, ctx->x1));
  ctx->opA->Mult(ctx->x1, ctx->y1);
  PetscCall(ToPetscVec(ctx->y1, y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __mat_apply_NEP_J(Mat J, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  palace::slepc::SlepcNEPSolver *ctx;
  PetscCall(MatShellGetContext(J, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context for SLEPc!");
  PetscCall(FromPetscVec(x, ctx->x1));
  ctx->opJ->Mult(ctx->x1, ctx->y1);
  PetscCall(ToPetscVec(ctx->y1, y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __mat_apply_NEP_B(Mat A, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  palace::slepc::SlepcNEPSolver *ctx;
  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context for SLEPc!");
  PetscCall(FromPetscVec(x, ctx->x1));
  ctx->opB->Mult(ctx->x1.Real(), ctx->y1.Real());
  ctx->opB->Mult(ctx->x1.Imag(), ctx->y1.Imag());
  PetscCall(ToPetscVec(ctx->y1, y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __pc_apply_NEP(PC pc, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  palace::slepc::SlepcNEPSolver *ctx;
  PetscCall(PCShellGetContext(pc, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell PC context for SLEPc!");
  PetscCall(FromPetscVec(x, ctx->x1));
  // Updating PC for new λ is needed for SLP, but should not be done for NLEIGS.
  if (ctx->new_lambda && !ctx->first_pc)
  {
    if (ctx->lambda.imag() == 0.0)
      ctx->lambda = ctx->sigma;
    ctx->opA2_pc = ctx->BuildA2(ctx->lambda);
    ctx->opA_pc = palace::BuildParSumOperator(
        {1.0 + 0.0i, ctx->lambda, ctx->lambda * ctx->lambda, 1.0 + 0.0i},
        {ctx->opK, ctx->opC, ctx->opM, ctx->opA2_pc.get()}, true);
    ctx->opP_pc =
        (*ctx->funcP)(std::complex<double>(1.0, 0.0), ctx->lambda,
                      ctx->lambda * ctx->lambda, ctx->PreconditionerBCFreq(ctx->lambda));
    ctx->opInv->SetOperators(*ctx->opA_pc, *ctx->opP_pc);
    ctx->new_lambda = false;
  }
  else if (ctx->first_pc)
  {
    ctx->first_pc = false;
    ctx->new_lambda = false;
  }
  ctx->opInv->Mult(ctx->x1, ctx->y1);
  if (ctx->opProj)
  {
    // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(ctx->GetComm(), ctx->y1));
    ctx->opProj->Mult(ctx->y1);
    // Mpi::Print(" After projection: {:e}\n", linalg::Norml2(ctx->GetComm(), ctx->y1));
  }
  PetscCall(ToPetscVec(ctx->y1, y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __form_NEP_function(NEP nep, PetscScalar lambda, Mat fun, Mat B, void *ctx)
{
  PetscFunctionBeginUser;
  palace::slepc::SlepcNEPSolver *ctxF;
  PetscCall(MatShellGetContext(fun, (void **)&ctxF));
  // A(λ) = K + λ C + λ² M + A2(λ). Routes through BuildA2 so the holomorphic A2(λ)
  // overload is preferred when funcA2Complex is set; otherwise falls back to the
  // legacy A2(|Im λ|) form.
  ctxF->opA2 = ctxF->BuildA2(lambda);
  ctxF->opA = palace::BuildParSumOperator(
      {1.0 + 0.0i, lambda, lambda * lambda, 1.0 + 0.0i},
      {ctxF->opK, ctxF->opC, ctxF->opM, ctxF->opA2.get()}, true);
  ctxF->lambda = lambda;
  ctxF->new_lambda = true;  // flag to update the preconditioner in SLP
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __form_NEP_jacobian(NEP nep, PetscScalar lambda, Mat fun, void *ctx)
{
  PetscFunctionBeginUser;
  palace::slepc::SlepcNEPSolver *ctxF;
  PetscCall(MatShellGetContext(fun, (void **)&ctxF));
  // A(λ) = K + λ C + λ² M + A2(λ).
  // J(λ) = C + 2 λ M + A2'(λ).
  // FD perturbation direction matches the A2 evaluation path: when funcA2Complex is
  // set the perturbation is a complex λ-step (imaginary axis); otherwise it's a
  // real-ω step on |Im λ| (legacy form, equivalent on the imag axis).
  ctxF->opA2 = ctxF->BuildA2(lambda);
  const auto eps = std::sqrt(std::numeric_limits<double>::epsilon());
  if (ctxF->use_complex_a2 && ctxF->funcA2Complex)
  {
    const std::complex<double> step{0.0, eps * std::abs(lambda.imag())};
    ctxF->opA2p = (*ctxF->funcA2Complex)(lambda + step);
    ctxF->opAJ = palace::BuildParSumOperator({1.0 / step, -1.0 / step},
                                             {ctxF->opA2p.get(), ctxF->opA2.get()}, true);
  }
  else
  {
    ctxF->opA2p = (*ctxF->funcA2)(std::abs(lambda.imag()) * (1.0 + eps));
    const std::complex<double> denom{0.0, eps * std::abs(lambda.imag())};
    ctxF->opAJ = palace::BuildParSumOperator({1.0 / denom, -1.0 / denom},
                                             {ctxF->opA2p.get(), ctxF->opA2.get()}, true);
  }
  ctxF->opJ = palace::BuildParSumOperator(
      {0.0 + 0.0i, 1.0 + 0.0i, 2.0 * lambda, 1.0 + 0.0i},
      {ctxF->opK, ctxF->opC, ctxF->opM, ctxF->opAJ.get()}, true);
  PetscFunctionReturn(PETSC_SUCCESS);
}

#endif
