// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "slepc.hpp"

#if defined(PALACE_WITH_SLEPC)

#include <algorithm>
#include <petsc.h>
#include <slepc.h>
#include <mfem.hpp>
#include "linalg/divfree.hpp"
#include "linalg/rap.hpp" // test for ComplexParOperator
#include "linalg/orthog.hpp" // test for orthog
#include <Eigen/Dense> // test for deflation
#include "linalg/nleps.hpp" // test for newton interp?
#include "utils/communication.hpp"

static PetscErrorCode __mat_apply_EPS_A0(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_EPS_A1(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_EPS_B(Mat, Vec, Vec);
static PetscErrorCode __pc_apply_EPS(PC, Vec, Vec);
static PetscErrorCode __mat_apply_PEPLinear_L0(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_PEPLinear_L1(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_PEPLinear_B(Mat, Vec, Vec);
static PetscErrorCode __pc_apply_PEPLinear(PC, Vec, Vec);
static PetscErrorCode __mat_apply_PEPLinear2_L0(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_PEPLinear2_L1(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_PEPLinear2_B(Mat, Vec, Vec);
static PetscErrorCode __pc_apply_PEPLinear2(PC, Vec, Vec);
static PetscErrorCode __mat_apply_PEP_A0(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_PEP_A1(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_PEP_A2(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_PEP_B(Mat, Vec, Vec);
static PetscErrorCode __pc_apply_PEP(PC, Vec, Vec);
static PetscErrorCode __mat_apply_PEPCheb_A0(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_PEPCheb_A1(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_PEPCheb_A2(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_PEPCheb_AN(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_PEPCheb_B(Mat, Vec, Vec);
static PetscErrorCode __pc_apply_PEPCheb(PC, Vec, Vec);
static PetscErrorCode __mat_apply_PEPRat_A0(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_PEPRat_A1(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_PEPRat_A2(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_PEPRat_AN(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_PEPRat_B(Mat, Vec, Vec);
static PetscErrorCode __pc_apply_PEPRat(PC, Vec, Vec);
static PetscErrorCode __mat_apply_NEP_A(Mat, Vec, Vec);
static PetscErrorCode __mat_destroy_NEP_A(Mat);
static PetscErrorCode __mat_duplicate_NEP_A(Mat, MatDuplicateOption, Mat *);
static PetscErrorCode __mat_apply_NEP_J(Mat, Vec, Vec);  // remove?
// static PetscErrorCode __mat_destroy_NEP_J(Mat);
static PetscErrorCode __mat_duplicate_NEP_J(Mat, MatDuplicateOption, Mat *);
static PetscErrorCode __pc_apply_NEP(PC, Vec, Vec);
static PetscErrorCode __form_NEP_function(NEP, PetscScalar, Mat, Mat, void *);
static PetscErrorCode __form_NEP_jacobian(NEP, PetscScalar, Mat, void *);  // remove?
static PetscErrorCode __compute_singularities(NEP, PetscInt*, PetscScalar*, void*);

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
  sigma_max = 0.0;
  l0 = 0.0;
  has_A2 = false;
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

void SlepcEigenvalueSolver::SetOperators(SpaceOperator &space_op,
  const ComplexOperator &K, const ComplexOperator &M, EigenvalueSolver::ScaleType type)
{
MFEM_ABORT("SetOperators not defined for base class SlepcEigenvalueSolver!");
}

void SlepcEigenvalueSolver::SetOperators(const ComplexOperator &K, const ComplexOperator &C,
                                         const ComplexOperator &M,
                                         EigenvalueSolver::ScaleType type)
{
  MFEM_ABORT("SetOperators not defined for base class SlepcEigenvalueSolver!");
}

void SlepcEigenvalueSolver::SetOperators(SpaceOperator &space_op,
  const ComplexOperator &K, const ComplexOperator &C,
  const ComplexOperator &M,
  EigenvalueSolver::ScaleType type)
{
MFEM_ABORT("SetOperators not defined for base class SlepcEigenvalueSolver!");
}

void SlepcEigenvalueSolver::SetLinearA2Operators(
  const ComplexOperator &A2_0, const ComplexOperator &A2_1,
  const ComplexOperator &A2_2)
{
MFEM_ABORT("SetLinearA2Operators not defined for base class SlepcEigenvalueSolver!");
}

void SlepcEigenvalueSolver::SetA2Interpolation(const Interpolation &A2interp)
{
  MFEM_ABORT("SetA2Interpolation not defined for base class SlepcEigenvalueSolver!");
};

void SlepcEigenvalueSolver::SetLinearSolver(/*const*/ ComplexKspSolver &ksp)
{
  opInv = &ksp;
}

void SlepcEigenvalueSolver::SetIoData(const IoData &iodata)
{
  opIodata = &iodata;
}

void SlepcEigenvalueSolver::SetDivFreeProjector(const DivFreeSolver<ComplexVector> &divfree)
{
  opProj = &divfree;
}

void SlepcEigenvalueSolver::SetBMat(const Operator &B)
{
  opB = &B;
}

//void SlepcEigenvalueSolver::SetShiftInvert(std::complex<double> s, std::complex<double> l, std::complex<double> s_max, bool precond)
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
  //PalacePetscCall(STSetTransform(st, PETSC_TRUE)); //need false for chebyshev?
  PalacePetscCall(STSetMatMode(st, ST_MATMODE_SHELL));
  sigma = s;  // Wait until solve time to call EPS/PEPSetTarget
  sinvert = true;
  //l0 = l;
  //sigma_max = s_max;
  l0 = sigma;
  sigma_max = 2.0 * sigma; // don't love this..
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
    //Mpi::Print("\nCAREFUL, TESTING REGION WIDTH IN SlepcEigenvalueSolver ConfigureRG!!\n");
    //if (sigma_max == 0.0) sigma_max = 10.0 * sigma; // hack
    if (PetscImaginaryPart(sigma) == 0.0)
    {
      PetscReal sr = PetscRealPart(sigma);
      PetscReal sr_max = PetscRealPart(sigma_max);
      if (sr > 0.0)
      {
        ConfigureRG(GetRG(), sr / gamma, mfem::infinity(), -mfem::infinity(),
        //ConfigureRG(GetRG(), sr / gamma, sr_max / gamma, -mfem::infinity(),
                    mfem::infinity());
      }
      else if (sr < 0.0)
      {
        ConfigureRG(GetRG(), -mfem::infinity(), sr / gamma, -mfem::infinity(),
        //ConfigureRG(GetRG(), - sr_max / gamma, sr / gamma, -mfem::infinity(),
                    mfem::infinity());
      }
    }
    else if (PetscRealPart(sigma) == 0.0)
    {
      PetscReal si = PetscImaginaryPart(sigma);
      PetscReal si_max = PetscImaginaryPart(sigma_max);
      if (si > 0.0)
      {
        ConfigureRG(GetRG(), -mfem::infinity(), mfem::infinity(), si / gamma,
                    mfem::infinity());
        //            si_max / gamma);
      }
      else if (si < 0.0)
      {
        ConfigureRG(GetRG(), -mfem::infinity(), mfem::infinity(), -mfem::infinity(),
        //ConfigureRG(GetRG(), -mfem::infinity(), mfem::infinity(), - si_max / gamma,
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

// NEP specific methods

SlepcNEPSolver::SlepcNEPSolver(MPI_Comm comm, int print, const std::string &prefix)
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
  T = TJ = nullptr;  // remove TJ?
  opK = opC = opM = nullptr;
  normA = normJ = 0.0;  // remove normJ?
}

SlepcNEPSolver::~SlepcNEPSolver()
{
  PalacePetscCall(NEPDestroy(&nep));
  PalacePetscCall(MatDestroy(&T));
  PalacePetscCall(MatDestroy(&TJ));
}

void SlepcNEPSolver::SetNumModes(int num_eig, int num_vec)
{
  PalacePetscCall(NEPSetDimensions(nep, num_eig, (num_vec > 0) ? num_vec : PETSC_DEFAULT,
                                   PETSC_DEFAULT));
}

void SlepcNEPSolver::SetTol(PetscReal tol)
{
  PalacePetscCall(NEPSetTolerances(nep, tol, PETSC_DEFAULT));
  PalacePetscCall(NEPSetConvergenceTest(nep, NEP_CONV_REL));
}

void SlepcNEPSolver::SetMaxIter(int max_it)
{
  PalacePetscCall(
      NEPSetTolerances(nep, PETSC_DEFAULT, (max_it > 0) ? max_it : PETSC_DEFAULT));
}

void SlepcNEPSolver::SetWhichEigenpairs(EigenvalueSolver::WhichType type)
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

void SlepcNEPSolver::SetProblemType(SlepcEigenvalueSolver::ProblemType type)
{
  switch (type)
  {
    case ProblemType::GENERAL:
      PalacePetscCall(NEPSetProblemType(nep, NEP_GENERAL));
      break;
    case ProblemType::RATIONAL:
      PalacePetscCall(NEPSetProblemType(nep, NEP_RATIONAL));
      break;
    case ProblemType::HERMITIAN:
    case ProblemType::NON_HERMITIAN:
    case ProblemType::GEN_HERMITIAN:
    case ProblemType::GEN_INDEFINITE:
    case ProblemType::GEN_NON_HERMITIAN:
    case ProblemType::HYPERBOLIC:
    case ProblemType::GYROSCOPIC:
      MFEM_ABORT("Problem type not implemented!");
      break;
  }
}

void SlepcNEPSolver::SetType(SlepcEigenvalueSolver::Type type)
{
  switch (type)
  {
    case Type::RII:
      PalacePetscCall(NEPSetType(nep, NEPRII));
      break;
    case Type::SLP:
      PalacePetscCall(NEPSetType(nep, NEPSLP));
      break;
    case Type::NARNOLDI:
      PalacePetscCall(NEPSetType(nep, NEPNARNOLDI));
      break;
    case Type::CISS:
      PalacePetscCall(NEPSetType(nep, NEPCISS));
      break;
    case Type::INTERPOL:
      PalacePetscCall(NEPSetType(nep, NEPINTERPOL));
      break;
    case Type::NLEIGS:
      PalacePetscCall(NEPSetType(nep, NEPNLEIGS));
      PalacePetscCall(NEPNLEIGSSetFullBasis(nep, PETSC_TRUE));
      // PalacePetscCall(NEPNLEIGSSetInterpolation(nep, PETSC_DETERMINE,
      // PETSC_DETERMINE));//what makes sense? PetscReal tol; PetscInt degree;
      // PalacePetscCall(NEPNLEIGSGetInterpolation(nep, &tol, &degree));
      // std::cerr << "Default NLEIGS Interpolation settings tol: " << tol << " degree: " <<
      // degree << "\n"; PalacePetscCall(NEPNLEIGSSetInterpolation(nep, 1e-6, 10));//what
      // makes sense?
      //PalacePetscCall(NEPNLEIGSSetInterpolation(nep, 1e-6, 50));  // what makes sense?
      PalacePetscCall(NEPNLEIGSSetInterpolation(nep, 1e-6, 30));  // 5 is fast but results are bad
      //PalacePetscCall(NEPNLEIGSSetInterpolation(nep, 1e-8, 300));  // test tighter, takes much longer
      //PalacePetscCall(NEPNLEIGSSetSingularitiesFunction(nep, __compute_singularities, NULL));

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

void SlepcNEPSolver::SetInitialSpace(const ComplexVector &v)
{
  MFEM_VERIFY(
      T,  // T && TJ,
      "Must call SetOperators before using SetInitialSpace for SLEPc eigenvalue solver!");
  if (!v0)
  {
    PalacePetscCall(MatCreateVecs(T, nullptr, &v0));
  }
  PalacePetscCall(ToPetscVec(v, v0));
  Vec is[1] = {v0};
  PalacePetscCall(NEPSetInitialSpace(nep, 1, is));
}

//void SlepcNEPSolver::SetShiftInvert(std::complex<double> s, std::complex<double> l, std::complex<double> s_max, bool precond)
void SlepcNEPSolver::SetShiftInvert(std::complex<double> s, bool precond)
{
  sigma = s;  // Wait until solve time to call NEPSetTarget
  sinvert = true;
  //sigma_max = s_max;
  //l0 = l;
  sigma_max = s * 2.0;
  l0 = s;
}

void SlepcNEPSolver::Customize()
{
  // Configure the region based on the given target if necessary.
  // NLEIGS requires a bounded region (can't use infinity bounds).
  if (sinvert && region)
  {
    std::cerr << "NEP region sigma: " << sigma << "\n";
    const double upper_bound_fac = 10.0; //80.0;
    const double width = 1.0; //1.0 // On CPW ex, width 0.2 takes super long compared to 1
    if (PetscImaginaryPart(sigma) == 0.0)
    {
      PetscReal sr = PetscRealPart(sigma);
      double bound = upper_bound_fac * std::abs(sr) / gamma;  // what value to use? same bounds for real and imag?
      if (sr > 0.0)
      {
        ConfigureRG(GetRG(), sr / gamma, bound, -bound*width, bound*width);
      }
      else if (sr < 0.0)
      {
        ConfigureRG(GetRG(), -bound, sr / gamma, -bound*width, bound*width);
      }
    }
    else if (PetscRealPart(sigma) == 0.0)
    {
      PetscReal si = PetscImaginaryPart(sigma);
      double bound = upper_bound_fac * std::abs(si) / gamma;  // what value to use? same bounds for real and imag?
      if (si > 0.0)
      {
        ConfigureRG(GetRG(), -bound*width, bound*width, si / gamma, bound);
      }
      else if (si < 0.0)
      {
        ConfigureRG(GetRG(), -bound*width, bound*width, -bound, si / gamma);
      }
    }
    else
    {
      MFEM_ABORT("Shift-and-invert with general complex eigenvalue target is unsupported!");
    }
  }
  std::cerr << "NEPSetTarget with sigma: " << sigma << "!!\n";
  PalacePetscCall(NEPSetTarget(nep, sigma / gamma));
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

int SlepcNEPSolver::Solve()
{
  MFEM_VERIFY(T /*&& TJ*/ && opInv, "Operators are not set for SlepcNEPSolver!");
  Mpi::Print("!!!!! NEPSolver::Solve etc\n");
  // Solve the eigenvalue problem.
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

  // Compute and store the eigenpair residuals.
  RescaleEigenvectors(num_conv);
  return (int)num_conv;
}

std::complex<double> SlepcNEPSolver::GetEigenvalue(int i) const
{
  PetscScalar l;
  PalacePetscCall(NEPGetEigenpair(nep, i, &l, nullptr, nullptr, nullptr));
  return l * gamma;
}

void SlepcNEPSolver::GetEigenvector(int i, ComplexVector &x) const
{
  MFEM_VERIFY(
      v0,
      "Must call SetOperators before using GetEigenvector for SLEPc eigenvalue solver!");
  PalacePetscCall(NEPGetEigenpair(nep, i, nullptr, nullptr, v0, nullptr));
  PalacePetscCall(FromPetscVec(v0, x));
  if (xscale.get()[i] > 0.0)
  {
    x *= xscale.get()[i];
  }
}

BV SlepcNEPSolver::GetBV() const
{
  BV bv;
  PalacePetscCall(NEPGetBV(nep, &bv));
  return bv;
}

ST SlepcNEPSolver::GetST() const
{
  ST st;
  // PalacePetscCall(NEPGetST(nep, &st)); // does not exist... what to do about this?
  return st;
}

RG SlepcNEPSolver::GetRG() const
{
  RG rg;
  PalacePetscCall(NEPGetRG(nep, &rg));
  return rg;
}

// NEPSetOperators does not exist, we need to call NEPSetFunction and NEPSetJacobian
// with corresponding FormFunctions
void SlepcNEPSolver::SetOperators(SpaceOperator &space_op_ref, const ComplexOperator &K,
                                  const ComplexOperator &C, const ComplexOperator &M,
                                  EigenvalueSolver::ScaleType type)
{
  const bool first = (opK == nullptr);
  has_A2 = true;
  opK = &K;
  opC = &C;
  opM = &M;
  space_op = &space_op_ref;

  if (first)
  {
    const PetscInt n = opK->Height();
    PalacePetscCall(
        MatCreateShell(GetComm(), n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &T));
    PalacePetscCall(MatShellSetOperation(T, MATOP_MULT, (void (*)(void))__mat_apply_NEP_A));
    // PalacePetscCall(
    //     MatShellSetOperation(T, MATOP_GET_DIAGONAL, (void
    //     (*)(void))__mat_diagonal_NEP_A));
    // PalacePetscCall(
    //     MatShellSetOperation(T, MATOP_DESTROY, (void (*)(void))__mat_destroy_NEP_A));
    PalacePetscCall(
        MatShellSetOperation(T, MATOP_DUPLICATE, (void (*)(void))__mat_duplicate_NEP_A));
    PalacePetscCall(MatShellSetVecType(T, PetscVecType()));
    PalacePetscCall(NEPSetFunction(nep, T, T, __form_NEP_function, NULL));

    // Jacobian, remove later it not needed
    PalacePetscCall(
        MatCreateShell(GetComm(), n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &TJ));
    PalacePetscCall(
        MatShellSetOperation(TJ, MATOP_MULT, (void (*)(void))__mat_apply_NEP_J));
    // PalacePetscCall(
    //     MatShellSetOperation(T, MATOP_GET_DIAGONAL, (void
    //     (*)(void))__mat_diagonal_NEP_J));
    // PalacePetscCall(
    //     MatShellSetOperation(T, MATOP_DESTROY, (void (*)(void))__mat_destroy_NEP_A));
    PalacePetscCall(
        MatShellSetOperation(TJ, MATOP_DUPLICATE, (void (*)(void))__mat_duplicate_NEP_J));
    PalacePetscCall(MatShellSetVecType(TJ, PetscVecType()));
    PalacePetscCall(NEPSetJacobian(nep, TJ, __form_NEP_jacobian, NULL));
  }

  // Remove this?
  /*
  if (first && type != ScaleType::NONE)
  {
    normA = linalg::SpectralNorm(GetComm(), *opA, opA->IsReal());
    normJ = linalg::SpectralNorm(GetComm(), *opJ, opJ->IsReal());
    MFEM_VERIFY(normA >= 0.0 && normJ >= 0.0, "Invalid matrix norms for EPS scaling!");
    if (normA > 0 && normJ > 0.0)
    {
      gamma = normA / normJ;  // Store γ² for linear problem
      delta = 2.0 / normA;
    }
  }
  */

  // Set up workspace.
  if (!v0)
  {
    PalacePetscCall(MatCreateVecs(T, nullptr, &v0));
  }
  x1.SetSize(opK->Height());  //?
  y1.SetSize(opK->Height());  //?
  x1.UseDevice(true);
  y1.UseDevice(true);

  // Configure linear solver. Here or elsewhere?
  if (first)
  {
    PC pc;
    PetscInt nsolve;
    // NLEIGS
    /*
    KSP *ksp;
    PalacePetscCall(NEPNLEIGSGetKSPs(nep, &nsolve, &ksp));
    for (int i = 0; i < nsolve; i++)
    {
      PalacePetscCall(KSPSetType(ksp[i], KSPPREONLY));
      PalacePetscCall(KSPGetPC(ksp[i], &pc));
      PalacePetscCall(PCSetType(pc, PCSHELL));
      PalacePetscCall(PCShellSetContext(pc, (void *)this));
      PalacePetscCall(PCShellSetApply(pc, __pc_apply_NEP));
    }
    */
    // RII (needs Jacobian and target_magnitude)
    /*
    KSP ksp;
    PalacePetscCall(NEPRIISetLagPreconditioner(nep, 1));
    //PalacePetscCall(NEPRIISetConstCorrectionTol(nep, PETSC_TRUE));
    //PalacePetscCall(NEPRIISetDeflationThreshold(nep, 1.0e-6));//?
    PalacePetscCall(NEPRIIGetKSP(nep, &ksp));
    PalacePetscCall(KSPSetType(ksp, KSPPREONLY));
    PalacePetscCall(KSPGetPC(ksp, &pc));
    PalacePetscCall(PCSetType(pc, PCSHELL));
    PalacePetscCall(PCShellSetContext(pc, (void *)this));
    PalacePetscCall(PCShellSetApply(pc, __pc_apply_NEP));
    */
    // SLP
    /**/
    KSP ksp;
    EPS eps;//?
    PalacePetscCall(NEPSLPGetKSP(nep, &ksp));
    PalacePetscCall(KSPSetType(ksp, KSPPREONLY));
    PalacePetscCall(NEPSLPGetEPS(nep, &eps));
    PalacePetscCall(EPSSetType(eps, EPSKRYLOVSCHUR)); //??
    PalacePetscCall(KSPGetPC(ksp, &pc));
    PalacePetscCall(PCSetType(pc, PCSHELL));
    PalacePetscCall(PCShellSetContext(pc, (void *)this));
    PalacePetscCall(PCShellSetApply(pc, __pc_apply_NEP));
   /**/
  }
}

PetscReal SlepcNEPSolver::GetResidualNorm(PetscScalar l, const ComplexVector &x,
                                          ComplexVector &r) const
{
  Mpi::Print("SlepcNEPSolver GetResidualNorm with l: {:e}{:+e}i\n", l.real(), l.imag());
  // Compute the i-th eigenpair residual: || A(λ) x ||₂ for eigenvalue λ.
  auto A2 = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(l.imag()), Operator::DIAG_ZERO); //std:abs???
  auto A = space_op->GetSystemMatrix(std::complex<double>(1.0, 0.0), l, l * l, opK, opC,
                                     opM, A2.get());
  A->Mult(x, r);
  return linalg::Norml2(GetComm(), r);
}

PetscReal SlepcNEPSolver::GetBackwardScaling(PetscScalar l) const
{
  // Make sure not to use norms from scaling as this can be confusing if they are different.
  // Note that SLEPc typically uses ||.||∞, not the 2-norm.
  auto A2 = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(l.imag()), Operator::DIAG_ZERO); //std:abs???
  auto A = space_op->GetSystemMatrix(std::complex<double>(1.0, 0.0), l, l * l, opK, opC,
                                     opM, A2.get());
  normA = linalg::SpectralNorm(GetComm(), *A, A->IsReal());
  return normA;
}

///// END  /////

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
  Mpi::Print("!!!!! EPSSolver::Solve etc with has_A2: {}\n", has_A2);
  if (has_A2)
  {
    //Mpi::Print("Creating opA2 in EPSSolverBase::Solve\n");
    // Test
    //const auto dl = std::sqrt(std::numeric_limits<double>::epsilon());
    //opA2 = space_op->GetExtraSystemMatrix<palace::ComplexOperator>(std::abs(l0.real()), palace::Operator::DIAG_ZERO);
    //opA2p = space_op->GetExtraSystemMatrix<palace::ComplexOperator>(std::abs(l0.real()) * (1.0 + dl), palace::Operator::DIAG_ZERO);
    //opAJ = space_op->GetExtraSystemMatrixJacobian<palace::ComplexOperator>(dl * std::abs(l0.real()), 1, opA2p.get(), opA2.get());//, opA2.get());
    Mpi::Print("Creating rational interpolation ops in EPSSolverBase::Solve with sigma: {}, sigma_max: {}\n", sigma.imag(), sigma_max.imag());
    NewtonA2_deg2(sigma, sigma_max);
  }

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
  if (i < eigen_values.size())
  {
    return eigen_values[i]; // * gamma; ??
  }
  PetscScalar l;
  PalacePetscCall(EPSGetEigenvalue(eps, i, &l, nullptr));
  return l * gamma;
}

void SlepcEPSSolverBase::GetEigenvector(int i, ComplexVector &x) const
{
  if (i < eigen_vectors.size())
  {
    x = eigen_vectors[i];
    return;
  }
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

void SlepcEPSSolver::SetOperators(SpaceOperator &space_op_ref, const ComplexOperator &K, const ComplexOperator &M,
                                  EigenvalueSolver::ScaleType type)
{
  // Construct shell matrices for the scaled operators which define the generalized
  // eigenvalue problem.
  const bool first = (opK == nullptr);
  has_A2 = true;//false;//true;
  opK = &K;
  opM = &M;
  space_op = &space_op_ref;

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
  if (has_A2) {
  auto A2 = space_op->GetExtraSystemMatrix<ComplexOperator>(std::sqrt(l).real(), Operator::DIAG_ZERO); //std:abs???
    A2->AddMult(x, r, std::complex<double>(1.0, 0.0));
  }
  //auto A = space_op->GetSystemMatrix(std::complex<double>(1.0, 0.0), l, l * l, opK, nullptr,  opM, A2.get()); // check if equiv to opK, opM, opA2 above??
  //A->Mult(x, r);
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
  opA2_0 = opA2_1 = opA2_2 = nullptr;//test
  opA2interp = nullptr;
  normK = normC = normM = 0.0;
}

//void SlepcPEPLinearSolver::SetOperators(const ComplexOperator &K, const ComplexOperator &C,
//                                        const ComplexOperator &M,
//                                        EigenvalueSolver::ScaleType type)
void SlepcPEPLinearSolver::SetOperators(SpaceOperator &space_op_ref, const ComplexOperator &K,
                                  const ComplexOperator &C, const ComplexOperator &M,
                                  EigenvalueSolver::ScaleType type)
{
  // Construct shell matrices for the scaled linearized operators which define the block 2x2
  // eigenvalue problem.
  const bool first = (opK == nullptr);
  has_A2 = true;//false;//true;
  opK = &K;
  opC = &C;
  opM = &M;
  space_op = &space_op_ref;

  if (first)
  {
    const PetscInt n = opK->Height();
    PalacePetscCall(MatCreateShell(GetComm(), 2 * n, 2 * n, PETSC_DECIDE, PETSC_DECIDE,
                                   (void *)this, &A0));
    PalacePetscCall(MatCreateShell(GetComm(), 2 * n, 2 * n, PETSC_DECIDE, PETSC_DECIDE,
                                   (void *)this, &A1));
    PalacePetscCall(
        MatShellSetOperation(A0, MATOP_MULT, (void (*)(void))__mat_apply_PEPLinear_L0)); //test
    PalacePetscCall(
        MatShellSetOperation(A1, MATOP_MULT, (void (*)(void))__mat_apply_PEPLinear_L1)); //test
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
    ConfigurePCShell(GetST(), (void *)this, __pc_apply_PEPLinear); //test
  }
}

void SlepcPEPLinearSolver::SetLinearA2Operators(const ComplexOperator &A2_0, const ComplexOperator &A2_1, const ComplexOperator &A2_2)
{
  Mpi::Print("PEPLinearSolver begin\n");
  has_A2 = true; // should be false otherwise!
  opA2_0 = &A2_0;
  opA2_1 = &A2_1;
  opA2_2 = &A2_2;
  Mpi::Print("PEPLinearSolver end\n");
}

void SlepcPEPSolver::SetLinearA2Operators(const ComplexOperator &A2_0, const ComplexOperator &A2_1, const ComplexOperator &A2_2)
{
  Mpi::Print("PEPSolver begin\n");
  has_A2 = true;
  opA2_0 = &A2_0;
  opA2_1 = &A2_1;
  opA2_2 = &A2_2;
  Mpi::Print("PEPSolver end\n");
}

void SlepcPEPLinearSolver::SetA2Interpolation(const Interpolation &A2interp)
{
  Mpi::Print("PEPLinearSolver begin\n");
  has_A2 = true; // should be false otherwise!
  opA2interp =  &A2interp;
  Mpi::Print("PEPLinearSolver end\n");
}

void SlepcPEPSolver::SetA2Interpolation(const Interpolation &A2interp)
{
  Mpi::Print("PEPSolver begin\n");
  has_A2 = true; // should be false otherwise!
  opA2interp =  &A2interp;
  Mpi::Print("PEPLinearSolver end\n");
}

int SlepcPEPLinearSolver::Solve()
{
  MFEM_VERIFY(A0 && A1 && opInv, "Operators are not set for SlepcPEPSolverBase!"); // CHECK IF IT SHOULD BE A0 and A1 or something else??
  Mpi::Print("!!!!! PEPLinearSolver::Solve etc\n");
  if (has_A2)
  {
    Mpi::Print("Creating opA2 in PEPLinearSolverBase::Solve\n");
    // Test
    //const auto dl = std::sqrt(std::numeric_limits<double>::epsilon());
    //opA2 = space_op->GetExtraSystemMatrix<palace::ComplexOperator>(std::abs(l0.imag()), palace::Operator::DIAG_ZERO);
    //opA2p = space_op->GetExtraSystemMatrix<palace::ComplexOperator>(std::abs(l0.imag()) * (1.0 + dl), palace::Operator::DIAG_ZERO);
    //opAJ = space_op->GetExtraSystemMatrixJacobian<palace::ComplexOperator>(dl * std::abs(l0.imag()), 1, opA2p.get(), opA2.get());//, opA2.get());
    NewtonA2_deg2(sigma, sigma_max);
    //ChebyshevA2(2, sigma.imag(), sigma_max.imag());
  }

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

  // TEST TO REFINE EIGENVALUES WITH QUASI-NEWTON METHOD?
  int max_outer_it = 100;
  int max_inner_it = 5;
  double tol = 1e-6;
  ComplexVector v, u, w, c, w0, z;
  v.SetSize(opK->Height());
  u.SetSize(opK->Height());
  w.SetSize(opK->Height());
  c.SetSize(opK->Height());
  w0.SetSize(opK->Height());
  z.SetSize(opK->Height());
  v.UseDevice(true);
  u.UseDevice(true);
  w.UseDevice(true);
  c.UseDevice(true);
  w0.UseDevice(true);
  z.UseDevice(true);
  // Set arbitrary w0? or set arbitrary c and solve for w0?
  linalg::SetRandom(GetComm(), c);
  opInv->Mult(c, w0);
  const auto dl = std::sqrt(std::numeric_limits<double>::epsilon()); // TODO: define only once above
  space_op->GetWavePortOp().SetSuppressOutput(true); //suppressoutput!
  xscale = std::make_unique<PetscReal[]>(num_conv); // ??
  //eigen_values.resize(num_conv);
  //eigen_vectors.resize(num_conv);
  //Mpi::Print("Found {:d} eigenvalues in linear problem\n", num_conv);
  for (int i = 0; i < num_conv; i++)
  {
    std::complex<double> eig = GetEigenvalue(i);
    std::complex<double> init_eig = eig;
    xscale.get()[i] = 0.0; // annoying, maybe use rescale eigenvectors??
    GetEigenvector(i, v);
    //Mpi::Print("Refining eigenvalue #{:d}: {:e}{:+e}i\n", i, eig.real(), eig.imag());
    //
    // IF WE STICK WITH THIS, SHOULD PROBABLY MOVE REFINEMENT OUT OF SLEPC.CPP AND IN EIGENSOLVER.CPP OR ELSEWHERE?
    //
    /*
    // RII (similar to Algorithm 2 https://arxiv.org/pdf/1910.11712)
    int update_freq = 3;//test
    // update_freq=1  5, 14 its
    //1,      +7.219130e+00,      +3.798150e-01,      +9.516626e+00,      +1.261412e-10,      +2.701629e-07
    //2,      +7.761920e+00,      +3.774097e-01,      +1.029530e+01,      +4.620842e-10,      +9.906911e-07
    // update_freq=2  5, 9 its
    //1,      +7.219038e+00,      +3.797892e-01,      +9.517150e+00,      +1.457863e-11,      +3.122377e-08
    //2,      +7.760907e+00,      +3.661238e-01,      +1.061053e+01,      +2.696456e-11,      +5.781083e-08
    // update_freq=5  4, 8 its
    //1,      +7.219082e+00,      +3.798019e-01,      +9.516890e+00,      +2.215002e-10,      +4.743980e-07
    //2,      +7.759773e+00,      +3.710332e-01,      +1.046893e+01,      +3.102829e-10,      +6.652317e-07

    double res = 1.0;
    opA2 = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(eig.imag()), Operator::DIAG_ZERO);
    opA = space_op->GetSystemMatrix(std::complex<double>(1.0, 0.0), eig, eig * eig, opK, opC, opM, opA2.get());
        // WHEN eig =! eig.imag() the preconditioner matrix is no longer exact!!!!!!!!
    //opP = space_op->GetPreconditionerMatrix<ComplexOperator>(1.0, eig.imag(), -eig.imag() * eig.imag(), eig.imag());
    opP = space_op->GetPreconditionerMatrix<ComplexOperator>(std::complex<double>(1.0, 0.0), eig, eig * eig, eig.imag());
    opInv->SetOperators(*opA, *opP);

    opA->Mult(v, u);
    double init_res = linalg::Norml2(GetComm(), u);
    // Normalize eigenvector???
    v *= 1.0 / linalg::Norml2(GetComm(), v);

    int it = 0;
    while (it < max_outer_it)
    {
      int inner_it = 0;
      while (inner_it < max_inner_it)
      {
        auto A2 = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(eig.imag()), Operator::DIAG_ZERO); //std:abs???
        auto A = space_op->GetSystemMatrix(std::complex<double>(1.0, 0.0), eig, eig * eig, opK, opC, opM, A2.get());
        auto A2p = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(eig.imag()) * (1.0 + dl), Operator::DIAG_ZERO); // Maybe use this only for it=0 and after use prev A2?
        auto AJ = space_op->GetExtraSystemMatrixJacobian<ComplexOperator>(dl * std::abs(eig.imag()), 1, A2p.get(), A2.get(), A2.get());
        auto J = space_op->GetSystemMatrix(std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0), 2 * eig, opK, opC, opM, AJ.get());

        A->Mult(v, u);
        opInv->Mult(u, w);
        const std::complex<double> num = linalg::TransposeDot(GetComm(), w, v);//w.TransposeDot(v);
        J->Mult(v, u);
        opInv->Mult(u, w);
        const std::complex<double> den = linalg::TransposeDot(GetComm(), w, v);//w.TransposeDot(v);
        const auto mu = num / den;
        //Mpi::Print("\n RII inner iteration {:d}.{:d} eig: {:e}{:+e}i, |mu|/|eig|: {}\n\n", it, inner_it, eig.real(), eig.imag(), std::abs(mu)/std::abs(eig));

        if (std::abs(mu) < dl * std::abs(eig))
        {
          Mpi::Print("\nRII inner iteration {:d}.{:d}.{:d} reached tolerance with eig: {:e}{:+e}i\n\n", i, it, inner_it, eig.real(), eig.imag());
          break;
        }
        eig -= mu;
        inner_it++;
      }

      auto A2 = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(eig.imag()), Operator::DIAG_ZERO); //std:abs???
      auto A = space_op->GetSystemMatrix(std::complex<double>(1.0, 0.0), eig, eig * eig, opK, opC, opM, A2.get());
      A->Mult(v, u);
      res = linalg::Norml2(GetComm(), u);
      Mpi::Print(GetComm(), "RII outer iteration: {:d}, {:d}, eig: {:e}{:+e}i, res: {}\n", i, it, eig.real(), eig.imag(), res);
      if (res < tol)
      {
        eigen_values.push_back(eig);
        eigen_vectors.push_back(v);
        //eigen_values[i] = eig;
        //eigen_vectors[i] = v;
        Mpi::Print("\nRII for eig #{:d} reached tolerance in {:d} iterations with init res: {}, final res: {}, init_eig: {:e}{:+e}i, final eig: {:e}{:+e}i\n\n\n", i, it+1, init_res, res, init_eig.real(), init_eig.imag(), eig.real(), eig.imag());
        break;
      }

      // Update opInv or not??
      if (it % update_freq == 0)
      {
        Mpi::Print("Update opInv at outer it #{:d}\n", it);
        opA2 = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(eig.imag()), Operator::DIAG_ZERO);
        opA = space_op->GetSystemMatrix(std::complex<double>(1.0, 0.0), eig, eig * eig, opK, opC, opM, opA2.get());
            // WHEN eig =! eig.imag() the preconditioner matrix is no longer exact!!!!!!!!
        //opP = space_op->GetPreconditionerMatrix<ComplexOperator>(1.0, eig.imag(), -eig.imag() * eig.imag(), eig.imag());
        opP = space_op->GetPreconditionerMatrix<ComplexOperator>(std::complex<double>(1.0, 0.0), eig, eig * eig, eig.imag());
        opInv->SetOperators(*opA, *opP);
      }
      opInv->Mult(u, w);
      linalg::AXPBY(std::complex<double>(-1.0, 0.0), w, std::complex<double>(1.0, 0.0), v);
      v *= 1.0 / linalg::Norml2(GetComm(), v);
      it++;
      if (it == max_outer_it)
      {
        // We should instead set it to the lowest residual iteration, which MAY NOT be the last one
        eigen_values.push_back(eig);
        eigen_vectors.push_back(v);
        //eigen_values[i] = eig;
        //eigen_vectors[i] = v;
        Mpi::Print("\nRII for eig #{:d} reached max it with init res: {}, final res: {}, init_eig: {:e}{:+e}i, final eig: {:e}{:+e}i\n\n\n", i, init_res, res, init_eig.real(), init_eig.imag(), eig.real(), eig.imag());
      }
    }
    */

    // Waveguide adapter:
    // eig0: -1.114373e+00+2.118196e+01i
    // eig1: -1.110740e+00+2.277491e+01i
    /*
    // Quasi-Newton 2 from https://arxiv.org/pdf/1702.08492
    // Normalize eigenvector???
    v *= 1.0 / linalg::Norml2(GetComm(), v);

    opA2 = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(eig.imag()), Operator::DIAG_ZERO);
    opA = space_op->GetSystemMatrix(std::complex<double>(1.0, 0.0), eig, eig * eig, opK, opC, opM, opA2.get());
    // WHEN eig =! eig.imag() the preconditioner matrix is no longer exact!!!!!!!!
    //opP = space_op->GetPreconditionerMatrix<ComplexOperator>(1.0, eig.imag(), -eig.imag() * eig.imag(), eig.imag());
    opP = space_op->GetPreconditionerMatrix<ComplexOperator>(std::complex<double>(1.0, 0.0), eig, eig * eig, eig.imag());
    opInv->SetOperators(*opA, *opP);
    int update_freq = 3;//test
    // update_freq=1  43, 30 its
    //1,      +7.120867e+00,      +3.683617e-01,      +9.678516e+00,      +4.250972e-10,      +9.102878e-07
    //2,      +7.577326e+00,      +4.315980e-01,      +8.792450e+00,      +2.409751e-10,      +5.164587e-07
    // update_freq=2  43, 32 its
    //1,      +7.112863e+00,      +3.641351e-01,      +9.779582e+00,      +3.666994e-10,      +7.852252e-07
    //2,      +7.563062e+00,      +4.182341e-01,      +9.055476e+00,      +3.975758e-10,      +8.520614e-07
    // update_freq=5  51, 49 its
    //1,      +7.071786e+00,      +3.527732e-01,      +1.003560e+01,      +4.322906e-10,      +9.256083e-07
    //2,      +7.485256e+00,      +4.122635e-01,      +9.092000e+00,      +2.992645e-10,      +6.412705e-07

    double alpha = 1.0;
    double min_res = 1e6;
    int min_it = 0;
    double init_res = 1e6;
    double res = 1e6;
    double prev_res = 1e6;
    int it = 0;
    while (it < max_outer_it)
    {

      // Compute u = A * v and check residual
      auto A2n = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(eig.imag()), Operator::DIAG_ZERO); //std:abs???
      auto A = space_op->GetSystemMatrix(std::complex<double>(1.0, 0.0), eig, eig * eig, opK, opC, opM, A2n.get());

      A->Mult(v, u);
      //prev_res = res;
      res = linalg::Norml2(GetComm(), u) / linalg::Norml2(GetComm(), v);
      if (it == 0) init_res = res;
      if (res < min_res){min_res = res; min_it = it;}
      //min_res = std::min(res, min_res);
      Mpi::Print(GetComm(), "i: {}, it: {}, eig: {}, {}, res: {}, alpha: {}\n", i, it, eig.real(), eig.imag(), res, alpha);

      if (res < tol)
      {
        eigen_values.push_back(eig);
        eigen_vectors.push_back(v);
        break;
      }

      // Update safety factor?
      //if (res > prev_res) alpha = std::max(0.4, 0.9*alpha);
      //if (res < prev_res) alpha = std::min(1.0, 1.1*alpha);

      // Compute w = J * v
      opA2p = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(eig.imag()) * (1.0 + dl), Operator::DIAG_ZERO); // Maybe use this only for it=0 and after use prev A2?
      opAJ = space_op->GetExtraSystemMatrixJacobian<ComplexOperator>(dl * std::abs(eig.imag()), 1, opA2p.get(), A2n.get());
      opJ = space_op->GetSystemMatrix(std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0), 2 * eig, opK, opC, opM, opAJ.get());
      opJ->Mult(v, w);

      // Compute delta_eig = - dot(w0, u) / dot(w0, w)
      //std::complex<double> delta = - u.TransposeDot(w0) / w.TransposeDot(w0) * alpha; // safety factor? use some kind of momentum?
      std::complex<double> delta = - linalg::TransposeDot(GetComm(), u, w0) / linalg::TransposeDot(GetComm(), w, w0) * alpha; // safety factor? use some kind of momentum?

      // Update eigenvalue eig += delta_eig
      eig += delta;

      // Compute z = -(delta_eig * w + u)  TODO: reuse w instead of different variable z?
      linalg::AXPBYPCZ(-delta, w, std::complex<double>(-1.0, 0.0), u, std::complex<double>(0.0, 0.0), z);

      //  M (x_k+1 - x_k) = z
      // TEST UPDATE P and ksp every iteration?!
      if (it > 0 && it % update_freq == 0)
      {
        Mpi::Print("Update opInv\n");
        opA2 = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(eig.imag()), Operator::DIAG_ZERO);
        opA = space_op->GetSystemMatrix(std::complex<double>(1.0, 0.0), eig, eig * eig, opK, opC, opM, opA2.get());
            // WHEN eig =! eig.imag() the preconditioner matrix is no longer exact!!!!!!!!
        //opP = space_op->GetPreconditionerMatrix<ComplexOperator>(1.0, eig.imag(), -eig.imag() * eig.imag(), eig.imag());
        opP = space_op->GetPreconditionerMatrix<ComplexOperator>(std::complex<double>(1.0, 0.0), eig, eig * eig, eig.imag());
        opInv->SetOperators(*opA, *opP);
      }
      opInv->Mult(z, u);
      v += u;
      // Normalize eigenvector???
      v *= 1.0 / linalg::Norml2(GetComm(), v);

      it++;
      if (it == max_outer_it)
      {
        eigen_values.push_back(eig);
        eigen_vectors.push_back(v);
        break;
      }
    }
    Mpi::Print(GetComm(), "\n\n i: {}, init_res: {}, min_res: {}, min_it: {}\n\n", i, init_res, min_res, min_it);
    */
  }

  space_op->GetWavePortOp().SetSuppressOutput(false); //suppressoutput!

  // Compute and store the eigenpair residuals.
  RescaleEigenvectors(num_conv);
  return (int)num_conv;
}

void SlepcPEPLinearSolver::SetBMat(const Operator &B)
{
  SlepcEigenvalueSolver::SetBMat(B);

  const PetscInt n = B.Height();
  PalacePetscCall(MatCreateShell(GetComm(), 2 * n, 2 * n, PETSC_DECIDE, PETSC_DECIDE,
                                 (void *)this, &B0));
  PalacePetscCall(
      MatShellSetOperation(B0, MATOP_MULT, (void (*)(void))__mat_apply_PEPLinear_B));//test
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
  if (i < eigen_vectors.size())
  {
    x = eigen_vectors[i];
    return;
  }
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

  // SHOULDN'T BE NECESSARY TO DO THIS HERE...
  //const auto eps = std::sqrt(std::numeric_limits<double>::epsilon());
  //auto t_opA2 = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(l0.imag()), Operator::DIAG_ZERO);
  //auto t_opA2p = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(l0.imag()) * (1.0 + eps), Operator::DIAG_ZERO);
  //auto t_opAJ = space_op->GetExtraSystemMatrixJacobian<ComplexOperator>(eps * std::abs(l0.imag()), 1, t_opA2p.get(), t_opA2.get(), t_opA2.get());

  //Mpi::Print("GetResNorm sigma: {}, {}\n", sigma.real(), sigma.imag());
  //opK->Mult(x, r);
  //Mpi::Print("A2\n");
  //t_opA2->AddMult(x, r, std::complex<double>(1.0, 0.0));
  //Mpi::Print("J sigma\n");
  //t_opAJ->AddMult(x, r, -sigma);
  //Mpi::Print("J l.imag()\n");
  //t_opAJ->AddMult(x, r, std::complex<double>(0.0, l.imag()));
  //opC->AddMult(x, r, l);
  //opM->AddMult(x, r, l * l);
  //auto A2 = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(l.imag()), Operator::DIAG_ZERO); //std:abs???
  //A2->AddMult(x, r, 1.0); //test!?
  //Mpi::Print("GetResNorm done\n");

  Mpi::Print("PEPLinear GetResidualNorm with A2 for l: {:e}{:+e}i!\n", l.real(), l.imag());
  auto A2 = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(l.imag()), Operator::DIAG_ZERO); //std:abs???
  auto A = space_op->GetSystemMatrix(std::complex<double>(1.0, 0.0), l, l * l, opK, opC, opM, A2.get());
  A->Mult(x, r);
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
    case Type::LINEAR:
      PalacePetscCall(PEPSetType(pep, PEPLINEAR));
      //region = false;
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

int SlepcPEPSolverBase::Solve() // test
{
  MFEM_VERIFY(A0 && A1 && A2 && opInv, "Operators are not set for SlepcPEPSolverBase!");
  Mpi::Print("!!!!! PEPSolverBase::Solve etc\n");
  // Test
  if (has_A2)
  {
    const auto dl = std::sqrt(std::numeric_limits<double>::epsilon());
    opA2 = space_op->GetExtraSystemMatrix<palace::ComplexOperator>(std::abs(l0.imag()), palace::Operator::DIAG_ZERO);
    opA2p = space_op->GetExtraSystemMatrix<palace::ComplexOperator>(std::abs(l0.imag()) * (1.0 + dl), palace::Operator::DIAG_ZERO);
    opAJ = space_op->GetExtraSystemMatrixJacobian<palace::ComplexOperator>(dl * std::abs(l0.imag()), 1, opA2p.get(), opA2.get());
    //RationalA2(100, 1e-6);
    NewtonA2_deg2(sigma, sigma_max);
    //ChebyshevA2(2, sigma.imag(), sigma_max.imag());
  }

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

template <typename VecType, typename ScalarType>
inline void OrthogonalizeColumn(GmresSolverBase::OrthogType type, MPI_Comm comm,
                                const std::vector<VecType> &V, VecType &w, ScalarType *Rj,
                                int j)
{
  // Orthogonalize w against the leading j columns of V.
  switch (type)
  {
    case GmresSolverBase::OrthogType::MGS:
      linalg::OrthogonalizeColumnMGS(comm, V, w, Rj, j);
      break;
    case GmresSolverBase::OrthogType::CGS:
      linalg::OrthogonalizeColumnCGS(comm, V, w, Rj, j);
      break;
    case GmresSolverBase::OrthogType::CGS2:
      linalg::OrthogonalizeColumnCGS(comm, V, w, Rj, j, true);
      break;
  }
}

int SlepcPEPSolver::Solve() // test
{
  MFEM_VERIFY(A0 && A1 && A2 && opInv, "Operators are not set for SlepcPEPSolver!");
  Mpi::Print("!!!!! PEPSolver::Solve etc\n");
  // Test
  if (has_A2)
  {
    const auto dl = std::sqrt(std::numeric_limits<double>::epsilon());
    opA2 = space_op->GetExtraSystemMatrix<palace::ComplexOperator>(std::abs(l0.imag()), palace::Operator::DIAG_ZERO);
    opA2p = space_op->GetExtraSystemMatrix<palace::ComplexOperator>(std::abs(l0.imag()) * (1.0 + dl), palace::Operator::DIAG_ZERO);
    opAJ = space_op->GetExtraSystemMatrixJacobian<palace::ComplexOperator>(dl * std::abs(l0.imag()), 1, opA2p.get(), opA2.get());
    //RationalA2(100, 1e-6);
    NewtonA2_deg2(sigma, sigma_max);
    //ChebyshevA2(2, sigma.imag(), sigma_max.imag());
  }

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

  // Copied from PEPLinearSolve
 // TEST TO REFINE EIGENVALUES WITH QUASI-NEWTON METHOD?
 /*
  int max_outer_it = 100;//100;
  int max_inner_it = 5;//20;
  double tol = 1e-6;
  //double deflation_tol = 1e-8;//test, no effect if < tol
  ComplexVector v, vold, v1, u, w, c, w0, z;
  const int size = opK->Height();
  v.SetSize(size);
  vold.SetSize(size);
  v1.SetSize(size);
  u.SetSize(size);
  w.SetSize(size);
  c.SetSize(size);
  w0.SetSize(size);
  z.SetSize(size);
  v.UseDevice(true);
  vold.UseDevice(true);
  v1.UseDevice(true);
  u.UseDevice(true);
  w.UseDevice(true);
  c.UseDevice(true);
  w0.UseDevice(true);
  z.UseDevice(true);

  const auto dl = std::sqrt(std::numeric_limits<double>::epsilon()); // TODO: define only once above
  space_op->GetWavePortOp().SetSuppressOutput(true); //suppressoutput!
  xscale = std::make_unique<PetscReal[]>(num_conv); // ??

  Eigen::MatrixXcd Xeigen;//(size, num_conv); // for validity tests
  Eigen::VectorXcd v_eigen(size), c_eigen(size), w0_eigen(size), u_eigen(size), w_eigen(size); // for validity tests
  auto *cr = c.Real().Read();
  auto *ci = c.Imag().Read();
  for (int i = 0; i < size; i++) c_eigen(i) = std::complex<double>(cr[i], ci[i]);

  Eigen::MatrixXcd H;
  std::vector<ComplexVector> X;
  Eigen::VectorXcd u2, z2;


  int k = 0;
  while (k < num_conv)
  {
    // do we want to check c (defined above) inside the loop so it's random and thus avoids getting stuck on a "BAD" guess?
    // Set arbitrary c and solve for w0?
    linalg::SetRandom(space_op->GetComm(), c);
    std::srand((unsigned int) time(0)); // ?
    Eigen::VectorXcd c_eig = Eigen::VectorXcd::Random(k);
    Eigen::VectorXcd w0_eig;

    bool deflation = true;//false;
    std::complex<double> eig = GetEigenvalue(k);
    std::complex<double> eig_update = eig;
    std::complex<double> init_eig = eig;
    xscale.get()[k] = 0.0; // annoying, maybe use rescale eigenvectors??
    GetEigenvector(k, v);
    Eigen::VectorXcd v2 = Eigen::VectorXcd::Constant(k, 0.0);// how to initialize v2? zero or random?
    if (k > 0 && deflation)
    {
      std::srand((unsigned int) time(0));
      v2.setRandom();
      //std::cout << "random v2: " << v2 << "\n";
    }

    // Removed RII but we should test what works better when deflation is implemented
    // Quasi-Newton 2 from https://arxiv.org/pdf/1702.08492
    double norm_v = std::sqrt(linalg::Norml2(GetComm(), v, true) + v2.squaredNorm());
    v *= 1.0 / norm_v;
    v2 *= 1.0 / norm_v;
    norm_v = std::sqrt(linalg::Norml2(GetComm(), v, true) + v2.squaredNorm()); // should be 1...
    //std::cout << "v2 normalized: " << v2 << "\n";
    //linalg::Normalize(space_op->GetComm(), v);

    // validity tests
    auto *vr = v.Real().Read();
    auto *vi = v.Imag().Read();
    for (int i = 0; i < size; i++) v_eigen(i) = std::complex<double>(vr[i], vi[i]);


    // Deflation
    //https://arxiv.org/pdf/1910.11712
    // T*(l) = |T(l) U(l)|
    //         |A(l) B(l)|
    // |y1|= T*(l) |z1| ->
    // |y2|        |z2|
    // y1 = T(l) z1 + U(l) z2 = T(l) z1 + T(l)X(lI - H)^-1 z2
    // y2 = A(l) z1 + B(l) z2 = sum i=0..p l_i (XH^i)* z1 + sum i=0..p (XH^i)* X qi(l)z2
    // qi = sum j=0..i-1 l^j H^(i-j-1)

    // Linear solve |T(sigma) U(sigma)| |x1| = |b1|
    //              |A(sigma) B(sigma)| |x2|   |b2|
    // 1) S(sigma) = B(sigma) - A(sigma)*X*(sigma*I-H)^-1
    // 2) v = T(sigma)^-1 b1
    // 3) x2 = S(sigma)^-1 (b2 - A(sigma)v))
    // 4) x1 = v - X*(sigma*I-H)^-1 x2
    int p = 1; // minimality index

    opA2 = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(eig.imag()), Operator::DIAG_ZERO);
    opA = space_op->GetSystemMatrix(std::complex<double>(1.0, 0.0), eig, eig * eig, opK, opC, opM, opA2.get());
    opP = space_op->GetPreconditionerMatrix<ComplexOperator>(std::complex<double>(1.0, 0.0), eig, eig * eig, eig.imag());
    opInv->SetOperators(*opA, *opP);
    opInv->Mult(c, w0); // should determine w0_eig too!
    // update w0_eig
    if (k > 0 && deflation) //Effenberger
    {
    // w0 = T^1 c1 (done above with opInv)
    // w0_2 = SS^-1 (c2 - A*w0) where SS = B(sigma) - A(sigma) X S^-1. if p==1, B = 0 and A = X.adjoint() so SS = X^* X S^-1
    // w0_1 = w0 - X*S*w0_2
    w0_eig.conservativeResize(k);
    for (int j = 0; j < k; j++)
    {
      w0_eig(j) = c_eig(j) - linalg::Dot(GetComm(), w0, X[j]);
    }
    Eigen::MatrixXcd SS(k, k);
    for (int i = 0; i < k; i++)
    {
      for (int j = 0; j < k; j++)
      {
        SS(i,j) = linalg::Dot(GetComm(), X[i], X[j]);
      }
    }
    const auto S = (eig_update * Eigen::MatrixXcd::Identity(k, k) - H).inverse(); // use eig or something else?
    SS = - SS * S;
    w0_eig = SS.inverse() * w0_eig;
    // w0_1 = w0 - X*S*w0_eig
    const auto phi1 = S * w0_eig;
    ComplexVector phi2; phi2.SetSize(size); phi2.UseDevice(false); // keep on CPU?! Might want to implement methods in vector.cpp to do this
    phi2 = 0.0;
    auto *phi2r = phi2.Real().HostWrite();
    auto *phi2i = phi2.Imag().HostWrite();
    for (int j = 0; j < k; j++)
    {
      auto *XR = X[j].Real().Read();
      auto *XI = X[j].Imag().Read();
      for (int i = 0; i < size; i++)
      {
        phi2r[i] += phi1(j).real() * XR[i] - phi1(j).imag() * XI[i]; // transposedot
        phi2i[i] += phi1(j).imag() * XR[i] + phi1(j).real() * XI[i]; // transposedot
      }
    }
    linalg::AXPY(-1.0, phi2, w0);
    }
    auto *w0r = w0.Real().Read();
    auto *w0i = w0.Imag().Read();
    for (int i = 0; i < size; i++) w0_eigen(i) = std::complex<double>(w0r[i], w0i[i]);

    int update_freq = 4;// SHOULD REVISIT THIS AND FIGURE OUT BEST FREQUENCY around 4-5 seems good?

    double min_res = 1e6;
    int min_it = 0;
    double init_res = 1e6;
    double res = 1e6;
    int large_res = 0;
    int it = 0;
    while (it < max_outer_it)
    {
      // Compute u = A * v and check residual
      auto A2n = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(eig.imag()), Operator::DIAG_ZERO);
      auto A = space_op->GetSystemMatrix(std::complex<double>(1.0, 0.0), eig, eig * eig, opK, opC, opM, A2n.get());
      A->Mult(v, u);

      deflation = (k > 0); // could add a residual check?!

      // Effenberger deflation
      if (deflation)
      {
        // u1 = T(l) v1 + U(l) v2 = T(l) v1 + T(l)X(lI - H)^-1 v2
        const auto S = (eig * Eigen::MatrixXcd::Identity(k, k) - H).inverse();
        const auto phi1 = S * v2;
        ComplexVector phi2; phi2.SetSize(size); phi2.UseDevice(false); // keep on CPU?! Might want to implement methods in vector.cpp to do this
        phi2 = 0.0;
        auto *phi2r = phi2.Real().HostWrite();
        auto *phi2i = phi2.Imag().HostWrite();
        for (int j = 0; j < k; j++)
        {
          auto *XR = X[j].Real().Read();
          auto *XI = X[j].Imag().Read();
          for (int i = 0; i < size; i++)
          {
            phi2r[i] += phi1(j).real() * XR[i] - phi1(j).imag() * XI[i];
            phi2i[i] += phi1(j).imag() * XR[i] + phi1(j).real() * XI[i];
          }
        }
        // Test using Eigen
        Eigen::VectorXcd phi2_test(size), phi2_test1(size);
        phi2_test = Xeigen * phi1;
        auto *phi2r2 = phi2.Real().Read();
        auto *phi2i2 = phi2.Imag().Read();
        for (int i = 0; i < size; i++) phi2_test1(i) = std::complex<double>(phi2r2[i], phi2i2[i]);
        //std::cout << " diff phi2 eigen/palace: " << (phi2_test - phi2_test1).norm() << "\n";
        //if ((phi2_test - phi2_test1).norm() > 1e-8) std::cout << "\n\n\n\n DIFFERENT EIGEN AND PALACE phi2!!! \n\n\n\n";
        //exit(0);

        A->AddMult(phi2, u, 1.0);
        // u2 = A(l) v1 + B(l) v2 = sum i=0..p l_i (XH^i)* v1 + sum i=0..p (XH^i)* X qi(l) v2
        // if p = 1, A(l) = X.adjoint and B(l) = 0!
        u2.conservativeResize(k);
        Eigen::VectorXcd u2_0(k), u2_1(k), u2_2(k), u2_3(k);
        for (int j = 0; j < k; j++)
        {
          u2_0(j) = linalg::TransposeDot(GetComm(), v, X[j]); //??
          u2_1(j) = linalg::TransposeDot(GetComm(), X[j], v); // not sure, try switching X[j] and v and try using Dot instead?
          u2_2(j) = linalg::Dot(GetComm(), X[j], v);
          u2_3(j) = linalg::Dot(GetComm(), v, X[j]);
          u2(j) = linalg::Dot(GetComm(), v, X[j]);
        }
        // validity test
        Eigen::VectorXcd u2_test(k), u2_test2(k);
        for (int j = 0; j < k; j++) u2_test(j) = std::complex<double>(Xeigen.col(j).adjoint() * v_eigen);
        u2_test2 = Xeigen.adjoint() * v_eigen;
        //std::cout << " diff u2 test eigen : " << (u2_test2 - u2_test).norm() << "\n";
        //std::cout << " diff u2_0 eigen/palace: " << (u2_test - u2_0).norm() << "\n";
        //std::cout << " diff u2_1 eigen/palace: " << (u2_test - u2_1).norm() << "\n";
        //std::cout << " diff u2_2 eigen/palace: " << (u2_test - u2_2).norm() << "\n";
        //std::cout << " diff u2_3 eigen/palace: " << (u2_test - u2_3).norm() << "\n";
        // std::cout << " diff u2 eigen/palace: " << (u2_test - u2).norm() << "\n";
        //std::cout << "u2_test: " << u2_test << " u2: " << u2 << "\n";
        //if ((u2_test - u2).norm() > 1e-8) std::cout << "\n\n\n\n DIFFERENT EIGEN AND PALACE u2!!! \n\n\n\n";
        //exit(0);
      }
      auto *ur = u.Real().Read();
      auto *ui = u.Imag().Read();
      for (int i = 0; i < size; i++) u_eigen(i) = std::complex<double>(ur[i], ui[i]);

      if (deflation)
      {
        //Mpi::Print("L2289 squared norm(u): {}, norm(u2): {}, norm(v): {}, norm(v2): {}\n", linalg::Norml2(GetComm(), u, true), u2.squaredNorm(), linalg::Norml2(GetComm(), v, true), v2.squaredNorm());
        res = std::sqrt(linalg::Norml2(GetComm(), u, true) + u2.squaredNorm()) / norm_v;
        // val test
        double res_eigen = std::sqrt(u_eigen.squaredNorm() + u2.squaredNorm()) / norm_v;
        //std::cout << " diff res eigen/palace: " << (res_eigen - res) << "\n";
        //if (std::abs(res_eigen - res) > 1e-8) std::cout << "\n\n\n\n DIFFERENT EIGEN AND PALACE res!!! \n\n\n\n";
      }
      else
      {
        res = linalg::Norml2(space_op->GetComm(), u) / linalg::Norml2(space_op->GetComm(), v);
      }
      if (it == 0) init_res = res;
      if (res < min_res){min_res = res; min_it = it;}
      //min_res = std::min(res, min_res);
      Mpi::Print(space_op->GetComm(), "k: {}, it: {}, eig: {}, {}, res: {}\n", k, it, eig.real(), eig.imag(), res);

      if (res < tol) // In addition to res < tol, do we also want to check if it's within the range (i.e. greater than the target)???
      {
        // Update the invariant pair with normalization.
        const auto scale = linalg::Norml2(GetComm(), v);
        v *= 1.0 / scale;
        X.push_back(v);
        Xeigen.conservativeResize(size, k + 1);
        auto *vr = v.Real().Read();
        auto *vi = v.Imag().Read();
        for (int i = 0; i < size; i++) Xeigen(i, k) = std::complex<double>(vr[i], vi[i]);
        H.conservativeResizeLike(Eigen::MatrixXd::Zero(k + 1, k + 1));
        H.col(k).head(k) = v2 / scale;
        H(k, k) = eig;
        // Also store here?
        eigen_values.push_back(eig);
        eigen_vectors.push_back(v);
        k++; // only increment eigenpairs when a converged pair is found!
        break;
      }
      // Stop if large residual for 10 consecutive iterations
      large_res = (res > 0.9) ? large_res + 1 : 0;
      if (large_res > 10)
      {
        Mpi::Print(GetComm(), "k: {}, Newton not converging after {} iterations, restarting\n", k, it);
        break;
      }

      // Compute w = J * v
      opA2p = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(eig.imag()) * (1.0 + dl), Operator::DIAG_ZERO); // Maybe use this only for it=0 and after use prev A2?
      //opAJ = space_op->GetExtraSystemMatrixJacobian<ComplexOperator>(dl * std::abs(eig.imag()), 1, opA2p.get(), A2n.get());
      std::complex<double> denom = dl * std::abs(eig.imag());
      opAJ = space_op->GetDividedDifferenceMatrix<ComplexOperator>(denom, opA2p.get(), A2n.get(), Operator::DIAG_ZERO);
      opJ = space_op->GetSystemMatrix(std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0), 2 * eig, opK, opC, opM, opAJ.get());
      opJ->Mult(v, w);

      if (deflation)
      {
        // w1 = T'(l) v1 + U'(l) v2 = T'(l) v1 + T'(l)XS v2 - T(l)XS^2 v2
        const auto S = (eig * Eigen::MatrixXcd::Identity(k, k) - H).inverse(); // should re-use the one defined above!!
        const auto phi1 = S * v2;
        const auto phi2 = S * phi1;
        //std::cout << "phi2.size(): " << phi2.size() << " phi2: " << phi2 << "\n";
        ComplexVector phi3; phi3.SetSize(size); phi3.UseDevice(false); // keep on CPU?! Might want to implement methods in vector.cpp to do this
        ComplexVector phi4; phi4.SetSize(size); phi4.UseDevice(false); // keep on CPU?! Might want to implement methods in vector.cpp to do this
        phi3 = 0.0;
        phi4 = 0.0;
        auto *phi3r = phi3.Real().HostWrite();
        auto *phi3i = phi3.Imag().HostWrite();
        auto *phi4r = phi4.Real().HostWrite();
        auto *phi4i = phi4.Imag().HostWrite();
        for (int j = 0; j < k; j++)
        {
          auto *XR = X[j].Real().Read();
          auto *XI = X[j].Imag().Read();
          for (int i = 0; i < size; i++)
          {
            phi3r[i] += phi1(j).real() * XR[i] - phi1(j).imag() * XI[i];
            phi3i[i] += phi1(j).imag() * XR[i] + phi1(j).real() * XI[i];
            phi4r[i] += phi2(j).real() * XR[i] - phi2(j).imag() * XI[i];
            phi4i[i] += phi2(j).imag() * XR[i] + phi2(j).real() * XI[i];
          }
        }
        // Tests for validity
        Eigen::VectorXcd phi3_test = Xeigen * S * v2;
        Eigen::VectorXcd phi4_test = Xeigen * S * S * v2;
        Eigen::VectorXcd phi3_test1(size), phi4_test1(size);
        auto *phi3r1 = phi3.Real().Read();
        auto *phi3i1 = phi3.Imag().Read();
        auto *phi4r1 = phi4.Real().Read();
        auto *phi4i1 = phi4.Imag().Read();
        for (int i = 0; i < size; i++)
        {
          phi3_test1(i) = std::complex<double>(phi3r1[i], phi3i1[i]);
          phi4_test1(i) = std::complex<double>(phi4r1[i], phi4i1[i]);
        }
        //std::cout << " diff phi3 eigen/palace: " << (phi3_test - phi3_test1).norm() << "\n";
        //std::cout << " diff phi4 eigen/palace: " << (phi4_test - phi4_test1).norm() << "\n";
        //if ((phi3_test - phi3_test1).norm() > 1e-8) std::cout << "\n\n\n\n DIFFERENT EIGEN AND PALACE phi3!!! \n\n\n\n";
        //if ((phi4_test - phi4_test1).norm() > 1e-8) std::cout << "\n\n\n\n DIFFERENT EIGEN AND PALACE phi4!!! \n\n\n\n";
        opJ->AddMult(phi3, w, 1.0);
        A->AddMult(phi4, w, -1.0);
        // with p = 1, A'(l) and B'(l) are zero so w2 is zero?
      }
      // validity test
      auto *wr = w.Real().Read();
      auto *wi = w.Imag().Read();
      for (int i = 0; i < size; i++) w_eigen(i) = std::complex<double>(wr[i], wi[i]);

      // Compute delta = - dot(w0, u) / dot(w0, w)
      std::complex<double> delta;
      if (deflation)
      {
        // test, not sure
        //std::complex<double> u2_w0 = u2.adjoint() * w0_eig.head(k);
        std::complex<double> u2_w0 = w0_eig.adjoint() * u2;
        //Mpi::Print("L2364 u2_w0: {}+{}i\n", u2_w0.real(), u2_w0.imag());
        //delta = - (linalg::TransposeDot(space_op->GetComm(), u, w0) + u2_w0) / linalg::TransposeDot(space_op->GetComm(), w, w0); // maybe wrong?!
        std::complex<double> delta_0 = - (linalg::TransposeDot(space_op->GetComm(), u, w0) + u2_w0) / linalg::TransposeDot(space_op->GetComm(), w, w0);
        std::complex<double> delta_1 = - (linalg::TransposeDot(space_op->GetComm(), w0, u) + u2_w0) / linalg::TransposeDot(space_op->GetComm(), w0, w);
        std::complex<double> delta_2 = - (linalg::Dot(space_op->GetComm(), u, w0) + u2_w0) / linalg::Dot(space_op->GetComm(), w, w0); // this is right?!?!!?
        std::complex<double> delta_3 = - (linalg::Dot(space_op->GetComm(), w0, u) + u2_w0) / linalg::Dot(space_op->GetComm(), w0, w);
        delta = - (linalg::Dot(space_op->GetComm(), u, w0) + u2_w0) / linalg::Dot(space_op->GetComm(), w, w0);
        // validity test
        std::complex<double> num = w0_eigen.adjoint() * u_eigen;
        std::complex<double> den = w0_eigen.adjoint() * w_eigen;
        num += u2_w0;
        std::complex<double> delta_eigen = -num/den;
        //std::cout << " diff delta0 eigen/palace: " << (delta_eigen - delta_0) << "\n";
        //std::cout << " diff delta1 eigen/palace: " << (delta_eigen - delta_1) << "\n";
        //std::cout << " diff delta2 eigen/palace: " << (delta_eigen - delta_2) << "\n";
        //std::cout << " diff delta3 eigen/palace: " << (delta_eigen - delta_3) << "\n";
        //std::cout << " diff delta eigen/palace: " << (delta_eigen - delta) << "\n";
        //if (std::abs(delta_eigen - delta) > 1e-8) std::cout << "\n\n\n\n DIFFERENT EIGEN AND PALACE delta!!! \n\n\n\n";
      }
      else
      {
        //delta = - linalg::TransposeDot(space_op->GetComm(), u, w0) / linalg::TransposeDot(space_op->GetComm(), w, w0);
        delta = - linalg::Dot(space_op->GetComm(), u, w0) / linalg::Dot(space_op->GetComm(), w, w0);
      }
      // Update eigenvalue eig += delta
      eig += delta;

      // Compute z = -(delta * w + u)  TODO: reuse w instead of different variable z?
      z.AXPBYPCZ(-delta, w, std::complex<double>(-1.0, 0.0), u, std::complex<double>(0.0, 0.0));
      // z2 = -delta w2 - u2 = -u2 since w2 is zero?

      //  M (x_k+1 - x_k) = z
      // TEST UPDATE P and ksp every x iteration?!
      if (it > 0 && it % update_freq == 0)
      {
        eig_update = eig;
        //Mpi::Print("Update opInv\n");
        opA2 = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(eig.imag()), Operator::DIAG_ZERO);
        opA = space_op->GetSystemMatrix(std::complex<double>(1.0, 0.0), eig, eig * eig, opK, opC, opM, opA2.get());
        opP = space_op->GetPreconditionerMatrix<ComplexOperator>(std::complex<double>(1.0, 0.0), eig, eig * eig, eig.imag());
        opInv->SetOperators(*opA, *opP);
        opInv->Mult(c, w0);
        // update w0_eig
        if (k > 0 && deflation) //Effenberger
        {
        // w0 = T^1 c1 (done above with opInv)
        // w0_2 = SS^-1 (c2 - A*w0) where SS = B(sigma) - A(sigma) X S^-1. if p==1, B = 0 and A = X.adjoint() so SS = X^* X S^-1
        // w0_1 = w0 - X*S*w0_2

        w0_eig.conservativeResize(k);
        for (int j = 0; j < k; j++)
        {
          w0_eig(j) = c_eig(j) - linalg::Dot(GetComm(), w0, X[j]);
        }
        Eigen::MatrixXcd SS(k, k);
        for (int i = 0; i < k; i++)
        {
          for (int j = 0; j < k; j++)
          {
            SS(i,j) = linalg::Dot(GetComm(), X[i], X[j]);
          }
        }
        const auto S = (eig * Eigen::MatrixXcd::Identity(k, k) - H).inverse(); // use eig or something else?
        SS = - SS * S;
        w0_eig = SS.inverse() * w0_eig;
        // w0_1 = w0 - X*S*w0_eig
        const auto phi1 = S * w0_eig;
        ComplexVector phi2; phi2.SetSize(size); phi2.UseDevice(false); // keep on CPU?! Might want to implement methods in vector.cpp to do this
        phi2 = 0.0;
        auto *phi2r = phi2.Real().HostWrite();
        auto *phi2i = phi2.Imag().HostWrite();
        for (int j = 0; j < k; j++)
        {
          auto *XR = X[j].Real().Read();
          auto *XI = X[j].Imag().Read();
          for (int i = 0; i < size; i++)
          {
            phi2r[i] += phi1(j).real() * XR[i] - phi1(j).imag() * XI[i];
            phi2i[i] += phi1(j).imag() * XR[i] + phi1(j).real() * XI[i];
          }
        }
        linalg::AXPY(-1.0, phi2, w0);
        }
        auto *w0r = w0.Real().Read();
        auto *w0i = w0.Imag().Read();
        for (int i = 0; i < size; i++) w0_eigen(i) = std::complex<double>(w0r[i], w0i[i]);
      }
      opInv->Mult(z, u);
      // validity test
      auto *ur2 = u.Real().Read();
      auto *ui2 = u.Imag().Read();
      for (int i = 0; i < size; i++) u_eigen(i) = std::complex<double>(ur2[i], ui2[i]);
      if (k > 0 && deflation) //Effenberger
      {
        // u = T^1 z1 (done above with opInv)
        // u2 = SS^-1 (z2 - A*u) where SS = B(sigma) - A(sigma) X S^-1. if p==1, B = 0 and A = X.adjoint() so SS = X^* X S^-1
        // u1 = u - X*S*u2
        Eigen::VectorXcd z2 = -u2;
        u2.conservativeResize(k);
        Eigen::VectorXcd u2_0(k), u2_1(k), u2_2(k), u2_3(k);
        for (int j = 0; j < k; j++)
        {
          u2_0(j) = linalg::TransposeDot(GetComm(), u, X[j]); //??
          u2_1(j) = linalg::TransposeDot(GetComm(), X[j], u); // not sure, try switching X[j] and v and try using Dot instead?
          u2_2(j) = linalg::Dot(GetComm(), X[j], u);
          u2_3(j) = linalg::Dot(GetComm(), u, X[j]);
          u2(j) = linalg::Dot(GetComm(), u, X[j]);
        }
        // validity test
        Eigen::VectorXcd u2_test(k), u2_test2(k);
        for (int j = 0; j < k; j++) u2_test(j) = std::complex<double>(Xeigen.col(j).adjoint() * u_eigen);
        u2_test2 = Xeigen.adjoint() * u_eigen;
        //std::cout << " diff u2 test eigen : " << (u2_test2 - u2_test).norm() << "\n";
        //std::cout << " diff u2_0 eigen/palace: " << (u2_test - u2_0).norm() << "\n";
        //std::cout << " diff u2_1 eigen/palace: " << (u2_test - u2_1).norm() << "\n";
        //std::cout << " diff u2_2 eigen/palace: " << (u2_test - u2_2).norm() << "\n";
        //std::cout << " diff u2_3 eigen/palace: " << (u2_test - u2_3).norm() << "\n";
        //std::cout << " diff u2 eigen/palace: " << (u2_test - u2).norm() << "\n";
        //if ((u2_test - u2).norm() > 1e-8) std::cout << "\n\n\n\n DIFFERENT EIGEN AND PALACE u2 solve!!! \n\n\n\n";
        //std::cout << "u2_test: " << u2_test << " u2: " << u2 << "\n";
        z2 = z2 - u2;
        Eigen::MatrixXcd SS(k, k);
        Eigen::MatrixXcd SS_0(k, k), SS_1(k, k), SS_2(k, k), SS_3(k, k);
        Eigen::MatrixXcd SS_test(k, k), SS_test2(k, k);
        SS_test2 = Xeigen.adjoint() * Xeigen;
        for (int i = 0; i < k; i++)
        {
          for (int j = 0; j < k; j++)
          {
            SS_test(i,j) = std::complex<double>(Xeigen.col(i).adjoint() * Xeigen.col(j));
            SS_0(i,j) = linalg::TransposeDot(GetComm(), X[j], X[i]);
            SS_1(i,j) = linalg::TransposeDot(GetComm(), X[i], X[j]);
            SS_2(i,j) = linalg::Dot(GetComm(), X[j], X[i]);
            SS_3(i,j) = linalg::Dot(GetComm(), X[i], X[j]);
            SS(i,j) = linalg::Dot(GetComm(), X[i], X[j]); // 2 or 3??
          }
        }
        //std::cout << " diff SS test eigen : " << (SS_test2 - SS_test).norm() << "\n";
        //std::cout << " diff SS_0 eigen/palace: " << (SS_test - SS_0).norm() << "\n";
        //std::cout << " diff SS_1 eigen/palace: " << (SS_test - SS_1).norm() << "\n";
        //std::cout << " diff SS_2 eigen/palace: " << (SS_test - SS_2).norm() << "\n";
        //std::cout << " diff SS_3 eigen/palace: " << (SS_test - SS_3).norm() << "\n";
        //std::cout << " diff SS eigen/palace: " << (SS_test - SS).norm() << "\n";
        //if ((SS_test - SS).norm() > 1e-8) std::cout << "\n\n\n\n DIFFERENT EIGEN AND PALACE SS!!! \n\n\n\n";

        const auto S = (eig_update * Eigen::MatrixXcd::Identity(k, k) - H).inverse(); // use eig or something else?
        SS = - SS * S;
        u2 = SS.inverse() * z2;
        // u1 = u - X*S*u2
        const auto phi1 = S * u2;
        ComplexVector phi2; phi2.SetSize(size); phi2.UseDevice(false); // keep on CPU?! Might want to implement methods in vector.cpp to do this
        phi2 = 0.0;
        auto *phi2r = phi2.Real().HostWrite();
        auto *phi2i = phi2.Imag().HostWrite();
        for (int j = 0; j < k; j++)
        {
          auto *XR = X[j].Real().Read();
          auto *XI = X[j].Imag().Read();
          for (int i = 0; i < size; i++)
          {
            phi2r[i] += phi1(j).real() * XR[i] - phi1(j).imag() * XI[i]; // transposedot
            phi2i[i] += phi1(j).imag() * XR[i] + phi1(j).real() * XI[i]; // transposedot
          }
        }
        // Test using Eigen
        Eigen::VectorXcd phi2_test(size), phi2_test1(size);
        phi2_test = Xeigen * phi1;
        auto *phi2r2 = phi2.Real().Read();
        auto *phi2i2 = phi2.Imag().Read();
        for (int i = 0; i < size; i++) phi2_test1(i) = std::complex<double>(phi2r2[i], phi2i2[i]);
        //std::cout << " diff phi2 solve eigen/palace: " << (phi2_test - phi2_test1).norm() << "\n";
        //if ((phi2_test - phi2_test1).norm() > 1e-8) std::cout << "\n\n\n\n DIFFERENT EIGEN AND PALACE phi2 solve!!! \n\n\n\n";
        linalg::AXPY(-1.0, phi2, u);
        v2 += u2; // update v2?
        //std::cout << "updated v2: " << v2 << "\n";
      }
      v += u;
      norm_v = std::sqrt(linalg::Norml2(GetComm(), v, true) + v2.squaredNorm());
      v *= 1.0 / norm_v;
      v2 *= 1.0 / norm_v;
      norm_v = std::sqrt(linalg::Norml2(GetComm(), v, true) + v2.squaredNorm()); // should be 1...
      //std::cout << "updated v2 normalized: " << v2 << "\n";
      auto *vr = v.Real().Read();
      auto *vi = v.Imag().Read();
      for (int i = 0; i < size; i++) v_eigen(i) = std::complex<double>(vr[i], vi[i]);

      it++;
      if (it == max_outer_it) // ACTUALLY, WE DO NOT WANT TO SAVE UNCONVERGED EIGENPAIRS, CAN MESS UP FUTURE DEFLATION ITERATIONS
      {
        Mpi::Print("Quasi-Newton did not converge in {} iterations\n", max_outer_it);
      }
    }
    Mpi::Print(space_op->GetComm(), "\n\n i: {}, init_res: {}, min_res: {}, min_it: {}\n\n", k, init_res, min_res, min_it);
  }

  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eps;
  eps.compute(H);
  //std::cout << "eps.eigenvalues: " << eps.eigenvalues() << "\n";
  // check that eps.eigenvalues().size() == eigen_values.size()
  std::vector<int> order(eigen_values.size()), order_eigen(eps.eigenvalues().size()), order2(eigen_values.size());
  std::iota(order.begin(), order.end(), 0);
  std::iota(order_eigen.begin(), order_eigen.end(), 0);
  std::iota(order2.begin(), order2.end(), 0);
  std::sort(order.begin(), order.end(),
            [&eigen_values = this->eigen_values](auto l, auto r) { return eigen_values[l].imag() < eigen_values[r].imag(); });
  std::sort(order_eigen.begin(), order_eigen.end(),
            [&epseig = eps.eigenvalues()](auto l, auto r) { return epseig(l).imag() < epseig(r).imag(); });
  std::sort(order2.begin(), order2.end(), [&order](auto l, auto r) { return order[l] < order[r]; });
  for(int k = 0; k < num_conv; k++) Mpi::Print("k: {}, order: {}, order2: {}, eig: {}+{}i\n", k, order[k], order2[k], eigen_values[k].real(), eigen_values[k].imag());//std::cout << "k: " << k << " order: " << order[k] << " eig: " << eigen_values[k] << "\n";
  for(int k = 0; k < num_conv; k++) Mpi::Print("k: {}, ordereig: {}, eig: {}+{}i\n", k, order_eigen[k], eps.eigenvalues()(k).real(), eps.eigenvalues()(k).imag());//std::cout << "k: " << k << " ordereig: " << order_eigen[k] << " eig: " << eps.eigenvalues()(k) << "\n";
  Eigen::MatrixXcd Xeig = eps.eigenvectors();
  std::vector<Eigen::VectorXcd> sorted_Xeig;
  std::vector<ComplexVector> sorted_eigen_vectors;
  std::vector<std::complex<double>> sorted_eigen_values;
  for (int k = 0; k < num_conv; k++)
  {
    //Mpi::Print("k: {}, eig: {}+{}i, epseig: {}+{}i\n", k, eigen_values[order[k]].real(), eigen_values[order[k]].imag(), eps.eigenvalues()(order_eigen[k]).real(), eps.eigenvalues()(order_eigen[k]).imag());
    //std::cout << "k: " << k << " eig: " << eigen_values[order[k]] << " epseig: " << eps.eigenvalues()(order_eigen[k]) << "\n";
    sorted_Xeig.push_back(Xeig.col(order_eigen[k]));//should be col?
    //sorted_eigen_vectors.push_back(eigen_vectors[order[k]]);
    //sorted_eigen_values.push_back(eigen_values[order[k]]);

    // original order
    //sorted_Xeig.push_back(Xeig.col(k));
    sorted_eigen_vectors.push_back(eigen_vectors[k]);
    sorted_eigen_values.push_back(eigen_values[k]);
  }

  //Xeigen *= eps.eigenvectors(); // [n x k] x [k x k] = [n x k]
  //Xeigen *= Xeigen.colwise().norm().cwiseInverse().asDiagonal(); //
  for(int k = 0; k < num_conv; k++)
  {
    double min_res = 1e6; int min_k = -1;
    ComplexVector eigv_best; eigv_best.SetSize(size); eigv_best.UseDevice(false);
    //for(int k2 = 0; k2 < num_conv; k2++) // find which Xeig columns leads to lowest residual
    //{
    int k2 = order2[k];//order_eigen[order[k]];
    //const size_t k_eigen = order_eigen[k]; // double check
    ComplexVector eigv; eigv.SetSize(size); eigv.UseDevice(false); // keep on CPU?! Might want to implement methods in vector.cpp to do this
    eigv = 0.0;
    auto *eigvr = eigv.Real().HostWrite();
    auto *eigvi = eigv.Imag().HostWrite();
    for (int j = 0; j < num_conv; j++)
    {
      //const size_t idx = order[j]; // double check;
      //const size_t j_eigen = order_eigen[j]; // double check
      //auto *XR = X[j].Real().Read();
      //auto *XI = X[j].Imag().Read();
      auto *XR = sorted_eigen_vectors[j].Real().Read();
      auto *XI = sorted_eigen_vectors[j].Imag().Read();
      for (int i = 0; i < size; i++)
      {
        //eigvr[i] += Xeig(j,k).real() * XR[i] - Xeig(j,k).imag() * XI[i];
        //eigvi[i] += Xeig(j,k).imag() * XR[i] + Xeig(j,k).real() * XI[i];
        eigvr[i] += sorted_Xeig[k2](j).real() * XR[i] - sorted_Xeig[k2](j).imag() * XI[i];
        eigvi[i] += sorted_Xeig[k2](j).imag() * XR[i] + sorted_Xeig[k2](j).real() * XI[i];
      }
    }
    const auto scale = linalg::Norml2(GetComm(), eigv);
    //std::cout << "k: " << k << " scale: " << scale << "\n";

    // evaluate residual
    std::complex<double> eig = sorted_eigen_values[k];
    auto A2n = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(eig.imag()), Operator::DIAG_ZERO); //std:abs???
    auto A = space_op->GetSystemMatrix(std::complex<double>(1.0, 0.0), eig, eig * eig, opK, opC, opM, A2n.get());
    A->Mult(eigv, u);
    double res = linalg::Norml2(GetComm(), u) / scale;
    if (res < min_res) {min_res = res; min_k = k2; eigv_best = eigv;}
    //} // test loop k2
    Mpi::Print("k: {}, min_res: {} at k: {}\n", k, min_res, min_k);
    eigen_vectors[k] = eigv_best; // double check k. Should reorder?!
    eigen_values[k] = sorted_eigen_values[k];
  }
    space_op->GetWavePortOp().SetSuppressOutput(false); //suppressoutput!
  */

  // Compute and store the eigenpair residuals.
  RescaleEigenvectors(num_conv);
  return (int)num_conv;
}

std::complex<double> SlepcPEPSolverBase::GetEigenvalue(int i) const
{
  if (i < eigen_values.size())
  {
    return eigen_values[i]; // * gamma; ??
  }
  PetscScalar l;
  PalacePetscCall(PEPGetEigenpair(pep, i, &l, nullptr, nullptr, nullptr));
  return l * gamma;
}

void SlepcPEPSolverBase::GetEigenvector(int i, ComplexVector &x) const
{
  if (i < eigen_vectors.size())
  {
    x = eigen_vectors[i];
    return;
  }
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
  opA2_0 = opA2_1 = opA2_2 = nullptr;//test
  opA2interp = nullptr;
  normK = normC = normM = 0.0;
}


// Write a function that builds the b vector instead?
// Do we need beta or should it just be 1 when there are no singularities?
double b_coeff(double lambda, int j, std::vector<double> &sigma_vec, std::vector<double> &beta)
{
  if (j == 0)
  {
    return 1.0; // / beta[0];
  }
  else
  {
    return (lambda - sigma_vec[j-1]) * b_coeff(lambda, j-1, sigma_vec, beta);// / beta[j];
  }
}

void SlepcPEPSolverBase::NewtonA2_deg2(PetscScalar target_min, PetscScalar target_max)
{
  int npoints = 3;//10;
  // interpolation points
  //std::vector<double> xs(npoints);
  xs.resize(npoints);
  for (int j = 0; j < npoints; j++) xs[j] = target_min + j * (target_max - target_min) / (npoints - 1);
  // Divided difference matrices (order 0 -> A2, order > 0 divided differences)
  //std::vector<std::vector<std::unique_ptr<ComplexOperator>>> D_j(npoints);
  D_j.resize(npoints);
  for (int k = 0; k < npoints; k++) // Order
  {
    for (int j = 0; j < npoints - k; j++)
    {
      if (k == 0)
      {
        auto A2i = space_op->GetExtraSystemMatrix<ComplexOperator>(xs[j].imag(), Operator::DIAG_ZERO);
        D_j[k].push_back(std::move(A2i));
      }
      else
      {
        //std::complex<double> denom = (xs[j+k]-xs[j]).imag();
        std::complex<double> denom = (xs[j+k]-xs[j]);
        auto A2dd = space_op->GetDividedDifferenceMatrix<ComplexOperator>(denom, D_j[k-1][j+1].get(), D_j[k-1][j].get(), Operator::DIAG_ZERO);
        D_j[k].push_back(std::move(A2dd));
      }
    }
  }
}

void SlepcEPSSolverBase::NewtonA2_deg2(PetscScalar target_min, PetscScalar target_max)
{
  int npoints = 3;//10;
  // interpolation points
  //std::vector<double> xs(npoints);
  xs.resize(npoints);
  for (int j = 0; j < npoints; j++) xs[j] = target_min + j * (target_max - target_min) / (npoints - 1);
  // Divided difference matrices (order 0 -> A2, order > 0 divided differences)
  //std::vector<std::vector<std::unique_ptr<ComplexOperator>>> D_j(npoints);
  D_j.resize(npoints);
  for (int k = 0; k < npoints; k++) // Order
  {
    for (int j = 0; j < npoints - k; j++)
    {
      if (k == 0)
      {
        auto A2i = space_op->GetExtraSystemMatrix<ComplexOperator>(xs[j].imag(), Operator::DIAG_ZERO);
        D_j[k].push_back(std::move(A2i));
      }
      else
      {
        //std::complex<double> denom = (xs[j+k]-xs[j]).imag();
        std::complex<double> denom = (xs[j+k]-xs[j]);
        auto A2dd = space_op->GetDividedDifferenceMatrix<ComplexOperator>(denom, D_j[k-1][j+1].get(), D_j[k-1][j].get(), Operator::DIAG_ZERO);
        D_j[k].push_back(std::move(A2dd));
      }
    }
  }
}

std::complex<double> b_coeff(
  std::complex<double> lambda, int j,
  std::vector<std::complex<double>> &sigma,
  std::vector<std::complex<double>> &xi,
  std::vector<std::complex<double>> &beta)
{
  if (j == 0)
  {
    return std::complex<double>(1.0, 0.0);
  }
  else
  {
    return (lambda - sigma[j-1]) / (beta[j] * (1.0 - lambda / xi[j])) * b_coeff(lambda, j-1, sigma, xi, beta);
  }
}
/*
void SlepcPEPSolverBase::RationalA2_deg2(
  const std::vector<std::complex<double>> &s,
  const std::vector<std::complex<double>> &xi,
  const std::vector<std::complex<double>> &beta)
{
  int npoints = s.size();

  D_j.resize(npoints);
  for (int k = 0; k < npoints; k++) // Order
  {
    for (int j = 0; j < npoints - k; j++)
    {
      if (k == 0)
      {
        auto A2i = space_op->GetExtraSystemMatrix<ComplexOperator>(s[j].imag(), Operator::DIAG_ZERO);
        D_j[k].push_back(std::move(A2i));
      }
      else
      {
        std::complex<double> denom = b_coeff()
        auto A2dd = space_op->GetDividedDifferenceMatrix<ComplexOperator>(denom, D_j[k-1][j+1].get(), D_j[k-1][j].get(), Operator::DIAG_ZERO);
        D_j[k].push_back(std::move(A2dd));
      }
    }
  }
}
  */

void SlepcPEPSolverBase::RationalA2(int degree, double tol)
{
  // R_d(λ) = sum_j b_j(λ) D_j
  // b_0(λ) = 1, b_j(λ) = (λ - σ_j-1) b_j-1 / ß_j
  // D_0 = ß_0 T(σ_0), D_j = (T_(σ_j) - R_j-1(σ_j)) / b_j(σ_j)
  // ß_j chosen so that max |b_j(λ)| = 1
  // d chosen so ||D_d|| / ||D_0|| < tol


  // For now let's just use 100 points spaced linearly in [sigma, 10*sigma]
  const int max_points = 100;
  const int max_sigma = 10;
  std::vector<double> sigma_vec(max_points), beta(max_points), nrs(max_points), s(max_points);
  double sigma_0 = sigma.imag();
  double dsigma = (max_sigma * sigma_0 - sigma_0) / double(max_points - 1);
  for (int j = 0; j < max_points; j++)
  {
    sigma_vec[j] = sigma_0 + double(j) * dsigma;
    nrs[j] = 1.0;
  }

  /*
  int ddmaxit = max_points;
  s[0] = sigma_vec[0];
  beta[0] = 1.0;
  for (int k=1; k < ddmaxit; k++)
  {
    double maxnrs = 0.0;
    for (int i = 0; i < max_points; i++)
    {
      nrs[i] *= (sigma_vec[i] - s[k-1])/beta[k-1];
      if (std::abs(nrs[i]) > maxnrs)
      {
        maxnrs = std::abs(nrs[i]);
        s[k] = sigma_vec[i];
      }
    }
    beta[k] = maxnrs;
    Mpi::Print("beta[{}]: {}\n", k, beta[k]);
  }
  */

  Mpi::Print("Computing A2 with sigma_0: {}\n", sigma_0);
  auto A2_0 = space_op->GetExtraSystemMatrix<ComplexOperator>(sigma_0, palace::Operator::DIAG_ZERO);
  int h = A2_0->Height(), w = A2_0->Width();
  double norm_0 = linalg::SpectralNorm(GetComm(), *A2_0);
  Mpi::Print("norm_D0: {}\n", norm_0);
  std::vector<std::unique_ptr<ComplexOperator>> D_j;
  std::vector<std::unique_ptr<ComplexOperator>> R_j;
  std::vector<std::unique_ptr<ComplexOperator>> A2_j; //needed?
  D_j.reserve(max_points);//not sure needed
  R_j.reserve(max_points);//not sure needed
  A2_j.reserve(max_points);//not sure needed
  D_j.push_back(std::move(A2_0));
  //D_j.emplace_back(std::move(A2_0));
  Mpi::Print("Computing norm_D0 using D_j\n");
  norm_0 = linalg::SpectralNorm(GetComm(), *D_j[0]);
  Mpi::Print("norm_D0 using D_j: {}\n", norm_0);
  for (int j = 1; j < max_points; j++)
  {
    Mpi::Print("j: {}\n", j);
    double b_j = b_coeff(sigma_vec[j], j, sigma_vec, beta);
    Mpi::Print("Computed b[{}]: {}\n", j, b_j);
    auto A2 = space_op->GetExtraSystemMatrix<ComplexOperator>(sigma_vec[j], palace::Operator::DIAG_ZERO);
    Mpi::Print("Computed A2\n");
    A2_j.push_back(std::move(A2));
    //A2_j.emplace_back(std::move(A2));
    Mpi::Print("Moved A2 into A2_j\n");
    auto sumr = std::make_unique<SumOperator>(h, w);
    auto sumi = std::make_unique<SumOperator>(h, w);
    for (int i = 0; i < j; i++) // i < j or j-1?
    {
      double b_i = b_coeff(sigma_vec[i], i, sigma_vec, beta);
      Mpi::Print("b[{}]: {}\n", i, b_i);
      const auto *PtAP_D = (D_j[i].get()) ? dynamic_cast<const ComplexParOperator *>(D_j[i].get()) : nullptr;
      Mpi::Print("Use PtAP_D with h x w: {} x {}\n", PtAP_D->LocalOperator().Height(), PtAP_D->LocalOperator().Width());
      if (PtAP_D && b_i != 0.0)
      {
        Mpi::Print("if PtAP_D and b_i != 0.0\n");
        //sumr->AddOperator(*D_j[i]->Real(), b_i);
        sumr->AddOperator(*PtAP_D->LocalOperator().Real(), b_i); //this isn't working, figure out why...
        Mpi::Print("Added D_j[{}]->Real to sumr\n", i);
        //sumi->AddOperator(*D_j[i]->Imag(), b_i);
        sumi->AddOperator(*PtAP_D->LocalOperator().Imag(), b_i);
        Mpi::Print("Added D_j[{}]->Imag to sumi\n", i);
      }
    }
    Mpi::Print("Computed sumr, sumi\n");
    auto R = std::make_unique<ComplexParOperator>(std::move(sumr), std::move(sumi), space_op->GetNDSpace());
    Mpi::Print("Built R\n");
    R_j.push_back(std::move(R));
    //R_j.emplace_back(std::move(R));
    Mpi::Print("Moved R into R_j\n");
    auto D = space_op->GetExtraSystemMatrixJacobian<ComplexOperator>(b_j, 1, A2_j[j-1].get(), R_j[j-1].get());
    Mpi::Print("Computed D\n");
    double norm = linalg::SpectralNorm(GetComm(), *D);
    Mpi::Print("Norm D: {}\n", norm);
    D_j.push_back(std::move(D));
    //D_j.emplace_back(std::move(D));
    Mpi::Print("Moved D into D_j\n");
  }
  std::exit(0);
}

//
// Convert T_k((2*lambda - (a+b))/(b-a)) to monomial form in lambda
std::vector<double> chebyshevPolynomialToMonomial(int k, double scale, double shift) {
    std::vector<double> coeffs(k + 1, 0.0);

    if (k == 0) {
        coeffs[0] = 1.0;  // T_0 = 1
        return coeffs;
    }

    if (k == 1) {
        // T_1(x) = x = (2*lambda - (a+b))/(b-a) = scale*lambda + shift * scale
        coeffs[0] = shift * scale;  // Constant term
        coeffs[1] = scale;       // Linear term
        return coeffs;
    }

    // Use recurrence: T_{k+1}(x) = 2*x*T_k(x) - T_{k-1}(x)
    std::vector<double> T_prev2 = chebyshevPolynomialToMonomial(k-2, scale, shift);
    std::vector<double> T_prev1 = chebyshevPolynomialToMonomial(k-1, scale, shift);

    // T_k = 2*x*T_{k-1} - T_{k-2}
    // where x = (2*lambda - (a+b))/(b-a) = scale * lambda + shift * scale

    // First compute 2*x*T_{k-1}
    std::vector<double> temp(k + 1, 0.0);
    for (int i = 0; i < T_prev1.size(); ++i) {
        // Multiply T_{k-1} by 2*x = 2*(2*lambda - (a+b))/(b-a) = 2 * scale * lambda + 2 * scale * shift
        if (i + 1 < temp.size()) {
            temp[i + 1] += T_prev1[i] * 2.0 * scale;  // lambda term
        }
        temp[i] += T_prev1[i] * 2 * scale * shift;  // constant term
    }

    // Subtract T_{k-2}
    for (int i = 0; i < T_prev2.size() && i < temp.size(); ++i) {
        temp[i] -= T_prev2[i];
    }

    return temp;
}

// Convert Chebyshev coefficients to monomial coefficients on original interval [a, b]
std::vector<double> chebyshevToMonomialOriginalRange(const std::vector<double>& cheb_coeffs, double a, double b) {

    int n = cheb_coeffs.size();
    std::vector<double> monomial_coeffs(n, 0.0);

    // Transformation parameters: x = (2*lambda - (a+b))/(b-a)
    // So lambda = ((b-a)*x + (a+b))/2
    double scale = 2.0 / (b - a);
    double shift = -(a + b) / 2.0;

    // Convert each Chebyshev polynomial T_k(x) to monomial form in lambda
    for (int k = 0; k < n; ++k) {
        if (std::abs(cheb_coeffs[k]) < 1e-15) continue;

        // Get monomial coefficients for T_k(x) where x = (2*lambda - (a+b))/(b-a)
        std::vector<double> tk_coeffs = chebyshevPolynomialToMonomial(k, scale, shift);

        // Add contribution to final polynomial
        for (int j = 0; j < tk_coeffs.size() && j < n; ++j) {
            monomial_coeffs[j] += cheb_coeffs[k] * tk_coeffs[j];
        }
    }

    return monomial_coeffs;
}

// l in [a, b], x in [-1, 1]
// x = (2 * (l - a) / (b - a)) - 1
// x = (2 * l - 2 * a - b + a) / (b - a) = (2 * l - (a + b)) / (b - a) = 2 * l / (b - a) - (a + b)/(b - a) = l * scale + shift
// scale = 2 / (b - a), shift = - (a + b)/(b - a)
// or x = scale * (l + shift) with scale = 2 / (b - a) and shift = - (a + b) / 2
void SlepcPEPSolverBase::ChebyshevA2(int degree, double target_min, double target_max)
{
  Mpi::Print("\n\n Now testing Chebyshev interpolation!\n\n");
  double b = target_max, a = target_min;
  std::vector<double> nodes(degree+1);
  std::vector<std::vector<double>> P_coeffs(degree+1);
  for (int i = 0; i <= degree; i++)
  {
    nodes[i] = 0.5 * (a + b) + 0.5 * (b - a) * std::cos((i + 0.5) * M_PI / (degree + 1.0));
    Mpi::Print("Chebyshev node[{}]: {}\n", i, nodes[i]);
    auto A2i = space_op->GetExtraSystemMatrix<ComplexOperator>(nodes[i], Operator::DIAG_ZERO);
    opA2_sample.push_back(std::move(A2i));
  }
  for (std::size_t i = 0; i <= degree; i++)
  {
    const double p = (i > 0) ? 2.0 : 1.0;
    for (std::size_t j = 0; j <= degree; j++)
    {
      const double xj = (2.0 * (nodes[j] - a) / (b - a)) - 1.0;
      const double coeff = p / (degree + 1.0) * std::cos(i * std::acos(xj));
      P_coeffs[i].push_back(coeff);
      //Mpi::Print("Polynomial coeff {}, {}: {}\n", i, j, coeff);
    }
    auto cheb = space_op->GetExtraSystemMatrixSum2(P_coeffs[i], opA2_sample, Operator::DIAG_ZERO);
    opA2_cheb.push_back(std::move(cheb));
  }

  std::vector<double> test_p = {nodes[0], nodes[degree], nodes[0]*0.9999, nodes[0]*0.99, nodes[0]*0.9, nodes[degree]*1.01, nodes[degree]*1.1, target_min, target_min*1.1, target_min*1.5, target_min*2};//, target_min*3.01, target_min*5, target_min*6.01, target_min*9};
  //std::vector<double> total_res(degree+1, 0.0);
  double total_res = 0.0;
  for (auto l : test_p)
  {
    auto A2_l = space_op->GetExtraSystemMatrix<ComplexOperator>(l, Operator::DIAG_ZERO);
    //auto normA2_l = linalg::SpectralNorm(GetComm(), *A2_l, A2_l->IsReal());
    ComplexVector tt(A2_l->Height());
    tt = std::complex<double>(1.23, 1.23);
    ComplexVector x1(A2_l->Height()), x2(A2_l->Height()), diff(A2_l->Height());
    x1 = 0.0; x2 = 0.0;
    A2_l->Mult(tt, x1);
    double normx1 = linalg::Norml2(space_op->GetComm(), x1);
    const double xj = (2.0 * (l - a) / (b - a)) - 1.0;
    Mpi::Print("test lambda: {} ({} in [-1,1]), with {} interpolation points\n", l, xj, degree+1);

    int k = degree;
    //for (int k = 0; k <= degree; k++)
    {
      x2 = 0.0;
      for(int j = 0; j <= k; j++)
      {
        opA2_cheb[j]->AddMult(tt, x2, std::cos(j * std::acos(xj)));
      }
      diff = 0.0; linalg::AXPBYPCZ(1.0, x1, -1.0, x2, 0.0, diff);
      double res = linalg::Norml2(space_op->GetComm(), diff);
      //Mpi::Print("Order {} res: {}, res/normx1: {}, \n\n", k, res, res/normx1);
      Mpi::Print("Order {} res: {}, res/normx1: {}, \n\n", degree, res, res/normx1);
      total_res += res/normx1;
    }

    // Monomial basis test - only works for order 2 for now
    // x = (2 * l - (a + b))/(b - a) = 2 * l/(b-a) - (a+b)/(b-a) = alpha * l + beta
    // alpha = 2 / (b - a), beta = -(a+b)/(b-a)
    if (k == 2)
    {
      double alpha = 2.0 / (b-a);
      double beta = -(a+b) / (b-a);
      x2 = 0.0;
      // a0 = c0 + c1 beta + c2 (2 beta^2 - 1)
      opA2_cheb[0]->AddMult(tt, x2, 1.0);
      opA2_cheb[1]->AddMult(tt, x2, beta);
      opA2_cheb[2]->AddMult(tt, x2, 2.0 * beta * beta - 1.0);
      // a1 = c1 alpha + c2 4 alpha beta
      opA2_cheb[1]->AddMult(tt, x2, alpha * l);
      opA2_cheb[2]->AddMult(tt, x2, 4.0 * alpha * beta * l);
      // a2 = 2 c2 alpha^2
      opA2_cheb[2]->AddMult(tt, x2, 2.0 * alpha * alpha * l * l);
      diff = 0.0; linalg::AXPBYPCZ(1.0, x1, -1.0, x2, 0.0, diff);
      double res = linalg::Norml2(space_op->GetComm(), diff);
      //Mpi::Print("Order {} res: {}, res/normx1: {}, \n\n", k, res, res/normx1);
      Mpi::Print("Order {} monomial res: {}, res/normx1: {}, \n\n", degree, res, res/normx1);

      cheb_to_mono.resize(3);
      cheb_to_mono[0].push_back(1.0);
      cheb_to_mono[0].push_back(beta);
      cheb_to_mono[0].push_back(2.0 * beta * beta - 1.0);
      cheb_to_mono[1].push_back(0.0);
      cheb_to_mono[1].push_back(alpha);
      cheb_to_mono[1].push_back(4.0 * alpha * beta);
      cheb_to_mono[2].push_back(0.0);
      cheb_to_mono[2].push_back(0.0);
      cheb_to_mono[2].push_back(2.0 * alpha * alpha);
    }



  }
  //for (int k = 0; k <= degree; k++)
  Mpi::Print("Cheb Order {}, total_res: {}\n", degree, total_res / test_p.size());
  //std::exit(0);
}

//void SlepcPEPSolver::SetOperators(const ComplexOperator &K, const ComplexOperator &C,
//                                  const ComplexOperator &M,
//                                  EigenvalueSolver::ScaleType type)
void SlepcPEPSolver::SetOperators(SpaceOperator &space_op_ref, const ComplexOperator &K,
                                  const ComplexOperator &C, const ComplexOperator &M,
                                  EigenvalueSolver::ScaleType type)
{
  // Construct shell matrices for the scaled operators which define the quadratic polynomial
  // eigenvalue problem.
  const bool first = (opK == nullptr);
  has_A2 = true;//false;//true;
  opK = &K;
  opC = &C;
  opM = &M;
  space_op = &space_op_ref;

  // Get Chebyshev polynomial degree?

  if (first)
  {
    //PalacePetscCall(PEPSetBasis(pep, PEP_BASIS_CHEBYSHEV1));
    const PetscInt n = opK->Height();
    PalacePetscCall(
        MatCreateShell(GetComm(), n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &A0));
    PalacePetscCall(
        MatCreateShell(GetComm(), n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &A1));
    PalacePetscCall(
        MatCreateShell(GetComm(), n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &A2));
    PalacePetscCall(
        //MatShellSetOperation(A0, MATOP_MULT, (void (*)(void))__mat_apply_PEP_A0));
        //MatShellSetOperation(A0, MATOP_MULT, (void (*)(void))__mat_apply_PEPCheb_A0));
        MatShellSetOperation(A0, MATOP_MULT, (void (*)(void))__mat_apply_PEPRat_A0));
    PalacePetscCall(
        //MatShellSetOperation(A1, MATOP_MULT, (void (*)(void))__mat_apply_PEP_A1));
        //MatShellSetOperation(A1, MATOP_MULT, (void (*)(void))__mat_apply_PEPCheb_A1));
        MatShellSetOperation(A1, MATOP_MULT, (void (*)(void))__mat_apply_PEPRat_A1));
    PalacePetscCall(
        //MatShellSetOperation(A2, MATOP_MULT, (void (*)(void))__mat_apply_PEP_A2));
        //MatShellSetOperation(A2, MATOP_MULT, (void (*)(void))__mat_apply_PEPCheb_A2));
        MatShellSetOperation(A2, MATOP_MULT, (void (*)(void))__mat_apply_PEPRat_A2));
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
    //ConfigurePCShell(GetST(), (void *)this, __pc_apply_PEP);
    //ConfigurePCShell(GetST(), (void *)this, __pc_apply_PEPCheb); //different or is it the same?
    ConfigurePCShell(GetST(), (void *)this, __pc_apply_PEPRat); //different or is it the same?
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
  Mpi::Print("PEPSolver GetResidualNorm with A2 for l: {:e}{:+e}i!\n", l.real(), l.imag());
  // Compute the i-th eigenpair residual: || P(λ) x ||₂ = || (K + λ C + λ² M) x ||₂ for
  // eigenvalue λ.
  auto A2 = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(l.imag()), Operator::DIAG_ZERO); //std:abs???
  opK->Mult(x, r);
  opC->AddMult(x, r, l);
  opM->AddMult(x, r, l * l);
  A2->AddMult(x, r, 1.0);//std::complex<double>(1.0, 0.0));
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
  if (ctx->has_A2) ctx->opA2->AddMult(ctx->x1, ctx->y1, std::complex<double>(1.0, 0.0));
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
  // Not sure why but this function appears to never be called?!
  std::exit(0);
  PetscCall(FromPetscVec(x, ctx->x1, ctx->x2));
  ctx->y1 = ctx->x2;
  ctx->opC->Mult(ctx->x2, ctx->y2);
  ctx->D_j[1][0]->AddMult(ctx->x2, ctx->y2, std::complex<double>(1.0, 0.0));
  ctx->D_j[2][0]->AddMult(ctx->x2, ctx->y2, -(ctx->xs[0] + ctx->xs[1]));
  //ctx->opAJ->AddMult(ctx->x2, ctx->y2, std::complex<double>(1.0, 0.0)); // TEST??? (OR SHOULD ONLY BE IMAG PART OF LAMBDA?!)
  ctx->y2 *= ctx->gamma;
  ctx->opK->AddMult(ctx->x1, ctx->y2, std::complex<double>(1.0, 0.0));
  ctx->D_j[0][0]->AddMult(ctx->x1, ctx->y2, std::complex<double>(1.0, 0.0));
  ctx->D_j[1][0]->AddMult(ctx->x1, ctx->y2, -ctx->xs[0]);
  ctx->D_j[2][0]->AddMult(ctx->x1, ctx->y2, ctx->xs[0]*ctx->xs[1]);
  //if (ctx->has_A2) ctx->opA2->AddMult(ctx->x1, ctx->y2, std::complex<double>(1.0, 0.0)); // TEST???
  //ctx->opAJ->AddMult(ctx->x1, ctx->y2, -ctx->l0); //TEST??? l0 like this?
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
  ctx->D_j[2][0]->AddMult(ctx->x2, ctx->y2, std::complex<double>(1.0, 0.0));
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
    //std::cerr << "__pc_apply_PEPLinear !sinvert\n";
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

  ctx->D_j[0][0]->AddMult(ctx->x1, ctx->y1, std::complex<double>(1.0, 0.0));
  ctx->D_j[1][0]->AddMult(ctx->x1, ctx->y1, -ctx->xs[0]);
  ctx->D_j[2][0]->AddMult(ctx->x1, ctx->y1, ctx->xs[0]*ctx->xs[1]);

    //if (ctx->has_A2) ctx->opA2->AddMult(ctx->x1, ctx->y1, std::complex<double>(1.0, 0.0)); // TEST???
    //ctx->opAJ->AddMult(ctx->x1, ctx->y1, -ctx->l0); // TEST??? for this to make sense maybe J needs to be in opInv!?
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


PetscErrorCode __mat_apply_PEPLinear2_L0(Mat A, Vec x, Vec y)
{
  // Apply the linearized operator L₀ = [ -K  0 ]
  //                                    [  0  I ] .
  PetscFunctionBeginUser;
  palace::slepc::SlepcPEPLinearSolver *ctx;
  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context for SLEPc!");

  PetscCall(FromPetscVec(x, ctx->x1, ctx->x2));
  ctx->opK->Mult(ctx->x1, ctx->y1);
  ctx->opA2->AddMult(ctx->x1, ctx->y1, std::complex<double>(1.0, 0.0));
  ctx->opJ->AddMult(ctx->x1, ctx->y2, -ctx->l0);
  ctx->y1 *= -1.0;
  ctx->y2 = ctx->x2;
  // TODO: FIGURE OUT SCALING!!
  PetscCall(ToPetscVec(ctx->y1, ctx->y2, y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __mat_apply_PEPLinear2_L1(Mat A, Vec x, Vec y)
{
  // Apply the linearized operator L₁ = [ C  M ]
  //                                    [ I  0 ] .
  PetscFunctionBeginUser;
  palace::slepc::SlepcPEPLinearSolver *ctx;
  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context for SLEPc!");

  PetscCall(FromPetscVec(x, ctx->x1, ctx->x2));
  ctx->opC->Mult(ctx->x1, ctx->y1);
  //ctx->opJ->AddMult(ctx->x1, ctx->y1, std::complex<double>(1.0, 0.0)); // TEST
  ctx->opM->AddMult(ctx->x2, ctx->y1, std::complex<double>(1.0, 0.0));
  ctx->y2 = ctx->x1;
  // TODO: FIGURE OUT SCALING!!
  //ctx->y1 = ctx->x1;
  //ctx->opM->Mult(ctx->x2, ctx->y2);
  //ctx->y2 *= ctx->delta * ctx->gamma * ctx->gamma;
  PetscCall(ToPetscVec(ctx->y1, ctx->y2, y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __mat_apply_PEPLinear2_B(Mat A, Vec x, Vec y)
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

PetscErrorCode __pc_apply_PEPLinear2(PC pc, Vec x, Vec y)
{
  // Solve the linear system associated with the generalized eigenvalue problem after
  // linearization: y = L₁⁻¹ x, or with the shift-and-invert spectral transformation:
  // y = (L₀ - σ L₁)⁻¹ x, with:
  //               L₀ = [ -K  0 ]    L₁ = [ C  M ]
  //                    [  0  I ] ,       [ I  0 ] .
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
    ctx->y1.AXPBY(-ctx->sigma / (ctx->delta * ctx->gamma), ctx->x1, 0.0);
    ctx->opK->AddMult(ctx->x2, ctx->y1, std::complex<double>(1.0, 0.0));
    ctx->opA2->AddMult(ctx->x2, ctx->y1, std::complex<double>(1.0, 0.0));
    ctx->opAJ->AddMult(ctx->x2, ctx->y1, -ctx->l0);
    ctx->opC->AddMult(ctx->x2, ctx->y1, ctx->sigma);
    ctx->opAJ->AddMult(ctx->x2, ctx->y1, ctx->sigma);
    ctx->opInv->Mult(ctx->y1, ctx->y2);
    if (ctx->opProj)
    {
      // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(ctx->GetComm(), ctx->y2));
      ctx->opProj->Mult(ctx->y2);
      // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(ctx->GetComm(), ctx->y2));
    }
    ctx->y1.AXPBYPCZ(ctx->gamma / ctx->sigma, ctx->y2, -ctx->gamma / ctx->sigma, ctx->x2,
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
  std::cerr << "__mat_apply_PEP_A0\n";
  std::exit(0); // THIS FUNCTION IS NEVER CALLED?!? WHY?!?!?
  PetscCall(FromPetscVec(x, ctx->x1));
  ctx->opK->Mult(ctx->x1, ctx->y1);
  //ctx->opA2->AddMult(ctx->x1, ctx->y1, std::complex<double>(1.0, 0.0)); //TEST???
  //ctx->opAJ->AddMult(ctx->x1, ctx->y1, -ctx->l0); //TEST???
  PetscCall(ToPetscVec(ctx->y1, y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __mat_apply_PEP_A1(Mat A, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  palace::slepc::SlepcPEPSolver *ctx;
  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context for SLEPc!");
  std::cerr << "__mat_apply_PEP_A1\n";
  PetscCall(FromPetscVec(x, ctx->x1));
  ctx->opC->Mult(ctx->x1, ctx->y1);
  //ctx->opAJ->AddMult(ctx->x1, ctx->y1, std::complex<double>(1.0, 0.0)); //TEST???
  PetscCall(ToPetscVec(ctx->y1, y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __mat_apply_PEP_A2(Mat A, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  palace::slepc::SlepcPEPSolver *ctx;
  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context for SLEPc!");
  std::cerr << "__mat_apply_PEP_A2\n";
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
  std::cerr << "__mat_apply_PEP_B\n";
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
  std::cerr << "__pc_apply_PEP\n";
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

PetscErrorCode __mat_apply_PEPCheb_A0(Mat A, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  palace::slepc::SlepcPEPSolver *ctx;
  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context for SLEPc!");
  std::cerr << "__mat_apply_PEPCheb_A0\n";
  std::exit(0); // THIS FUNCTION IS NEVER CALLED?!? WHY?!?!?
  PetscCall(FromPetscVec(x, ctx->x1));
  ctx->opK->Mult(ctx->x1, ctx->y1);
  ctx->opM->AddMult(ctx->x1, ctx->y1, std::complex<double>(0.5, 0.0));
  if (ctx->opA2_cheb.size() > 0) ctx->opA2_cheb[0]->AddMult(ctx->x1, ctx->y1, std::complex<double>(1.0, 0.0)); //coeff??
  //ctx->opA2->AddMult(ctx->x1, ctx->y1, std::complex<double>(1.0, 0.0)); //TEST???
  //ctx->opAJ->AddMult(ctx->x1, ctx->y1, -ctx->l0); //TEST???
  PetscCall(ToPetscVec(ctx->y1, y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __mat_apply_PEPCheb_A1(Mat A, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  palace::slepc::SlepcPEPSolver *ctx;
  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context for SLEPc!");
  std::cerr << "__mat_apply_PEPCheb_A1\n";
  PetscCall(FromPetscVec(x, ctx->x1));
  ctx->opC->Mult(ctx->x1, ctx->y1);
  if (ctx->opA2_cheb.size() > 1) ctx->opA2_cheb[1]->AddMult(ctx->x1, ctx->y1, std::complex<double>(1.0, 0.0)); // coeff??
  //ctx->opAJ->AddMult(ctx->x1, ctx->y1, std::complex<double>(1.0, 0.0)); //TEST???
  PetscCall(ToPetscVec(ctx->y1, y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __mat_apply_PEPCheb_A2(Mat A, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  palace::slepc::SlepcPEPSolver *ctx;
  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context for SLEPc!");
  std::cerr << "__mat_apply_PEPCheb_A2\n";
  PetscCall(FromPetscVec(x, ctx->x1));
  ctx->opM->Mult(ctx->x1, ctx->y1);
  ctx->y1 *= std::complex<double>(0.5, 0.0); //
  if (ctx->opA2_cheb.size() > 2) ctx->opA2_cheb[2]->AddMult(ctx->x1, ctx->y1, std::complex<double>(1.0, 0.0)); // coeff??
  PetscCall(ToPetscVec(ctx->y1, y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __mat_apply_PEPCheb_B(Mat A, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  palace::slepc::SlepcPEPSolver *ctx;
  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context for SLEPc!");
  std::cerr << "__mat_apply_PEPCheb_B\n";
  PetscCall(FromPetscVec(x, ctx->x1));
  ctx->opB->Mult(ctx->x1.Real(), ctx->y1.Real());
  ctx->opB->Mult(ctx->x1.Imag(), ctx->y1.Imag());
  ctx->y1 *= ctx->delta * ctx->gamma;
  PetscCall(ToPetscVec(ctx->y1, y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __pc_apply_PEPCheb(PC pc, Vec x, Vec y)
{
  // Solve the linear system associated with the generalized eigenvalue problem: y = M⁻¹ x,
  // or shift-and-invert spectral transformation: y = P(σ)⁻¹ x . Enforces the divergence-
  // free constraint using the supplied projector.
  PetscFunctionBeginUser;
  palace::slepc::SlepcPEPSolver *ctx;
  PetscCall(PCShellGetContext(pc, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell PC context for SLEPc!");
  std::cerr << "__pc_apply_PEPCheb\n";
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

PetscErrorCode __mat_apply_PEPRat_A0(Mat A, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  palace::slepc::SlepcPEPSolver *ctx;
  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context for SLEPc!");
  std::cerr << "__mat_apply_PEPRat_A0\n";
  std::exit(0); // THIS FUNCTION IS NEVER CALLED?!? WHY?!?!?
  PetscCall(FromPetscVec(x, ctx->x1));
  ctx->opK->Mult(ctx->x1, ctx->y1);
  if (ctx->has_A2)
  {
    ctx->opA2interp->AddMult(0, ctx->x1, ctx->y1, std::complex<double>(1.0, 0.0));
    //ctx->opA2_0->AddMult(ctx->x1, ctx->y1, std::complex<double>(1.0, 0.0));
    //ctx->opA2_1->AddMult(ctx->x1, ctx->y1, -ctx->xs[0]);
    //ctx->opA2_2->AddMult(ctx->x1, ctx->y1, ctx->xs[0]*ctx->xs[1]);
    //for (int i = 0; i < 3; i++) {if (ctx->cheb_to_mono[0][i] != 0.0) ctx->opA2_cheb[i]->AddMult(ctx->x1, ctx->y1, ctx->cheb_to_mono[0][i]);}
  }
  //ctx->D_j[0][0]->AddMult(ctx->x1, ctx->y1, std::complex<double>(1.0, 0.0));
  //ctx->D_j[1][0]->AddMult(ctx->x1, ctx->y1, -ctx->xs[0]);
  //ctx->D_j[2][0]->AddMult(ctx->x1, ctx->y1, ctx->xs[0]*ctx->xs[1]);
  PetscCall(ToPetscVec(ctx->y1, y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __mat_apply_PEPRat_A1(Mat A, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  palace::slepc::SlepcPEPSolver *ctx;
  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context for SLEPc!");
  std::cerr << "__mat_apply_PEPRat_A1\n";
  PetscCall(FromPetscVec(x, ctx->x1));
  ctx->opC->Mult(ctx->x1, ctx->y1);
  if (ctx->has_A2)
  {
    std::cerr << "__mat_apply_PEPRat_A1 with opA2\n";
    ctx->opA2interp->AddMult(1, ctx->x1, ctx->y1, std::complex<double>(1.0, 0.0));
    //ctx->opA2_1->AddMult(ctx->x1, ctx->y1, std::complex<double>(1.0, 0.0));
    //ctx->opA2_2->AddMult(ctx->x1, ctx->y1, -(ctx->xs[0] + ctx->xs[1]));
    //for (int i = 0; i < 3; i++) {if (ctx->cheb_to_mono[1][i] != 0.0) ctx->opA2_cheb[i]->AddMult(ctx->x1, ctx->y1, ctx->cheb_to_mono[1][i]);}
  }
  //ctx->D_j[1][0]->AddMult(ctx->x1, ctx->y1, std::complex<double>(1.0, 0.0));
  //ctx->D_j[2][0]->AddMult(ctx->x1, ctx->y1, -(ctx->xs[0] + ctx->xs[1]));
  PetscCall(ToPetscVec(ctx->y1, y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __mat_apply_PEPRat_A2(Mat A, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  palace::slepc::SlepcPEPSolver *ctx;
  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context for SLEPc!");
  std::cerr << "__mat_apply_PEPRat_A2\n";
  PetscCall(FromPetscVec(x, ctx->x1));
  ctx->opM->Mult(ctx->x1, ctx->y1);
  if (ctx->has_A2)
  {
    std::cerr << "__mat_apply_PEPRat_A2 with opA2\n";
    ctx->opA2interp->AddMult(2, ctx->x1, ctx->y1, std::complex<double>(1.0, 0.0));
    //ctx->opA2_2->AddMult(ctx->x1, ctx->y1, std::complex<double>(1.0, 0.0));
    //for (int i = 0; i < 3; i++) {if (ctx->cheb_to_mono[2][i] != 0.0) ctx->opA2_cheb[i]->AddMult(ctx->x1, ctx->y1, ctx->cheb_to_mono[2][i]);}
  }
  //ctx->D_j[2][0]->AddMult(ctx->x1, ctx->y1, std::complex<double>(1.0, 0.0));
  PetscCall(ToPetscVec(ctx->y1, y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __mat_apply_PEPRat_B(Mat A, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  palace::slepc::SlepcPEPSolver *ctx;
  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context for SLEPc!");
  std::cerr << "__mat_apply_PEPRat_B\n";
  PetscCall(FromPetscVec(x, ctx->x1));
  ctx->opB->Mult(ctx->x1.Real(), ctx->y1.Real());
  ctx->opB->Mult(ctx->x1.Imag(), ctx->y1.Imag());
  ctx->y1 *= ctx->delta * ctx->gamma;
  PetscCall(ToPetscVec(ctx->y1, y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __pc_apply_PEPRat(PC pc, Vec x, Vec y)
{
  // Solve the linear system associated with the generalized eigenvalue problem: y = M⁻¹ x,
  // or shift-and-invert spectral transformation: y = P(σ)⁻¹ x . Enforces the divergence-
  // free constraint using the supplied projector.
  PetscFunctionBeginUser;
  palace::slepc::SlepcPEPSolver *ctx;
  PetscCall(PCShellGetContext(pc, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell PC context for SLEPc!");
  std::cerr << "__pc_apply_PEPRat\n";
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
  //std::cerr << "In  __mat_apply_NEP_A with lambda: " << ctx->lambda_test << "\n";
  // PetscCall(FromPetscVec(x, ctx->x1));
  // std::cerr << "In  __mat_apply_NEP_A after frompetscvec\n";
  // std::cerr << "In  __mat_apply_NEP_A ctx->x1.Size(): " << ctx->x1.Size() << "
  // ctx->y1.Size(): " << ctx->y1.Size() << "\n";
  palace::ComplexVector t1(ctx->opA->Height()), t2(ctx->opA->Height());
  t1.UseDevice(true);
  t1 = 0.0;
  t2.UseDevice(true);
  t2 = 0.0;
  PetscCall(FromPetscVec(x, t1));
  ctx->opA->Mult(t1, t2);  // NEED TO FIGURE OUT WHY IT FAILS WITH x1, y1?!?!!?!?!?!?
  PetscCall(ToPetscVec(t2, y));
  // ctx->opA->Mult(ctx->x1, ctx->y1);
  // std::cerr << "In  __mat_apply_NEP_A after mult\n";
  //  ctx->y1 *= ctx->delta; //needed?
  // PetscCall(ToPetscVec(ctx->y1, y));
  // std::cerr << "In  __mat_apply_NEP_A after topetscvec\n";
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __mat_destroy_NEP_A(Mat A)
{
  PetscFunctionBeginUser;
  palace::slepc::SlepcNEPSolver *ctx;
  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context for SLEPc!");
  PetscCall(PetscFree(ctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __mat_duplicate_NEP_A(Mat A, MatDuplicateOption op, Mat *B)
{
  PetscFunctionBeginUser;
  palace::slepc::SlepcNEPSolver *ctxA, *ctxB;
  PetscInt m, n, M, N;
  MPI_Comm comm;
  PetscCall(MatShellGetContext(A, (void **)&ctxA));
  MFEM_VERIFY(ctxA, "Invalid PETSc shell matrix context for SLEPc!");
  PetscCall(MatGetSize(A, &M, &N));
  PetscCall(MatGetLocalSize(A, &m, &n));
  std::cerr << "In  __mat_duplicate_NEP_A \n";
  PetscCall(PetscNew(&ctxB));
  ctxB->space_op = ctxA->space_op;
  ctxB->opK = ctxA->opK;
  ctxB->opC = ctxA->opC;
  ctxB->opM = ctxA->opM;
  ctxB->lambda_test = ctxA->lambda_test;  // needed?
  ctxB->sigma = ctxA->sigma;
  ctxB->opInv = ctxA->opInv;

  std::complex<double> lambda = ctxB->lambda_test;// + ctxB->sigma; // TEST???
  if (lambda.imag() == 0.0) lambda += ctxB->sigma; // WHY
  std::cerr << "In __mat_duplicate_NEP_A ctxA->sigma: " << ctxA->sigma.real() << "+" << ctxA->sigma.imag() <<"i\n";
  std::cerr << "In __mat_duplicate_NEP_A lambda: " << lambda.real() << "+" << lambda.imag() <<"i\n";
  // why recalculate???? can't we just do ctxB->opA = ctxA->opA??
  ctxB->opA2 = ctxB->space_op->GetExtraSystemMatrix<palace::ComplexOperator>(
      std::abs(lambda.imag()), palace::Operator::DIAG_ZERO); //std:abs???
  ctxB->opA = ctxB->space_op->GetSystemMatrix(std::complex<double>(1.0, 0.0), lambda,
                                              lambda * lambda, ctxB->opK, ctxB->opC,
                                              ctxB->opM, ctxB->opA2.get());
  std::cerr << "In __mat_duplicate_NEP_A done calculating opA2 and opA\n";
  ctxB->x1.SetSize(ctxB->opK->Height());  //?
  ctxB->y1.SetSize(ctxB->opK->Height());  //?
  ctxB->x1.UseDevice(true);
  ctxB->y1.UseDevice(true);
  // need anything else?

  PetscCall(PetscObjectGetComm((PetscObject)A, &comm));
  PetscCall(MatCreateShell(comm, m, n, M, N, (void *)ctxB, B));
  PetscCall(MatShellSetOperation(*B, MATOP_MULT, (void (*)(void))__mat_apply_NEP_A));
  PetscCall(MatShellSetOperation(*B, MATOP_DESTROY, (void (*)(void))__mat_destroy_NEP_A));
  PetscCall(
      MatShellSetOperation(*B, MATOP_DUPLICATE, (void (*)(void))__mat_duplicate_NEP_A));
  PetscCall(MatShellSetVecType(*B, palace::slepc::PetscVecType()));
  PetscFunctionReturn(PETSC_SUCCESS);
}

// Remove if we don't need J
PetscErrorCode __mat_duplicate_NEP_J(Mat A, MatDuplicateOption op, Mat *B)
{
  PetscFunctionBeginUser;
  palace::slepc::SlepcNEPSolver *ctxA, *ctxB;
  PetscInt m, n, M, N;
  MPI_Comm comm;
  PetscCall(MatShellGetContext(A, (void **)&ctxA));
  MFEM_VERIFY(ctxA, "Invalid PETSc shell matrix context for SLEPc!");
  PetscCall(MatGetSize(A, &M, &N));
  PetscCall(MatGetLocalSize(A, &m, &n));
std::cerr << "In  __mat_duplicate_NEP_J \n";
  PetscCall(PetscNew(&ctxB));
  ctxB->space_op = ctxA->space_op;
  ctxB->opK = ctxA->opK;
  ctxB->opC = ctxA->opC;
  ctxB->opM = ctxA->opM;
  ctxB->lambda_J = ctxA->lambda_J;  // needed?
  ctxB->opInv = ctxA->opInv;
  ctxB->sigma = ctxA->sigma; //??
  // why recalculate???? can't we just do ctxB->opJ = ctxA->opJ??
  std::complex<double> lambda = ctxB->lambda_test;
  std::cerr << "In __mat_duplicate_NEP_J ctxA->sigma: " << ctxA->sigma.real() << "+" << ctxA->sigma.imag() <<"i\n";
  std::cerr << "In __mat_duplicate_NEP_J lambda: " << lambda.real() << "+" << lambda.imag() <<"i\n";
  ctxB->opA2 = ctxB->space_op->GetExtraSystemMatrix<palace::ComplexOperator>(
      std::abs(lambda.imag()), palace::Operator::DIAG_ZERO); //std:abs???
  const auto eps = std::sqrt(std::numeric_limits<double>::epsilon());
  ctxB->opA2p = ctxB->space_op->GetExtraSystemMatrix<palace::ComplexOperator>(
      std::abs(lambda.imag()) * (1.0 + eps), palace::Operator::DIAG_ZERO);//std:abs???
  ctxB->opAJ = ctxB->space_op->GetExtraSystemMatrixJacobian<palace::ComplexOperator>(eps * std::abs(lambda.imag()), 1, ctxB->opA2p.get(), ctxB->opA2.get());//, ctxB->opA2.get()); // third operator not used?
  ctxB->opJ = ctxB->space_op->GetSystemMatrix(
      std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0), 2.0 * lambda,
      ctxB->opK, ctxB->opC, ctxB->opM, ctxB->opAJ.get());

  ctxB->x1.SetSize(ctxB->opK->Height());  //?
  ctxB->y1.SetSize(ctxB->opK->Height());  //?
  ctxB->x1.UseDevice(true);
  ctxB->y1.UseDevice(true);
  // need anything else?

  PetscCall(PetscObjectGetComm((PetscObject)A, &comm));
  PetscCall(MatCreateShell(comm, m, n, M, N, (void *)ctxB, B));
  PetscCall(MatShellSetOperation(*B, MATOP_MULT, (void (*)(void))__mat_apply_NEP_J));
  PetscCall(MatShellSetOperation(*B, MATOP_DESTROY, (void (*)(void))__mat_destroy_NEP_A));
  PetscCall(
      MatShellSetOperation(*B, MATOP_DUPLICATE, (void (*)(void))__mat_duplicate_NEP_J));
  PetscCall(MatShellSetVecType(*B, palace::slepc::PetscVecType()));
  PetscFunctionReturn(PETSC_SUCCESS);
}

// Remove if we don't need J
PetscErrorCode __mat_apply_NEP_J(Mat J, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  palace::slepc::SlepcNEPSolver *ctx;
  PetscCall(MatShellGetContext(J, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context for SLEPc!");
std::cerr << "In  __mat_apply__NEP_J \n";
  PetscCall(FromPetscVec(x, ctx->x1));
  ctx->opJ->Mult(ctx->x1, ctx->y1);
  // ctx->y1 *= ctx->delta * ctx->gamma; //needed?
  PetscCall(ToPetscVec(ctx->y1, y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __pc_apply_NEP(PC pc, Vec x, Vec y)
{
  //
  PetscFunctionBeginUser;
  palace::slepc::SlepcNEPSolver *ctx;
  PetscCall(PCShellGetContext(pc, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell PC context for SLEPc!");

  PetscCall(FromPetscVec(x, ctx->x1));
  /**/
  if (ctx->new_lambda)
  {
    if (ctx->lambda_test.imag() == 0.0) ctx->lambda_test = ctx->sigma;
    std::cerr << "__pc_apply_NEP with new ctx->lambda_test: " << ctx->lambda_test.real() << "+" << ctx->lambda_test.imag() <<"i\n";
    std::cerr << "__pc_apply_NEP with new ctx->lambda_J: " << ctx->lambda_J.real() << "+" << ctx->lambda_J.imag() <<"i\n";
    //auto opA2_pc = ctx->space_op->GetExtraSystemMatrix<palace::ComplexOperator>(std::abs(ctx->lambda_test.imag()), palace::Operator::DIAG_ZERO);
    //auto opA_pc = ctx->space_op->GetSystemMatrix(std::complex<double>(1.0, 0.0), ctx->lambda_test, ctx->lambda_test * ctx->lambda_test, ctx->opK, ctx->opC, ctx->opM, opA2_pc.get());
    //auto opP_pc = ctx->space_op->GetPreconditionerMatrix<palace::ComplexOperator>(std::complex<double>(1.0, 0.0), ctx->lambda_test, ctx->lambda_test * ctx->lambda_test, ctx->lambda_test.imag());
    ctx->opA2_pc = ctx->space_op->GetExtraSystemMatrix<palace::ComplexOperator>(std::abs(ctx->lambda_test.imag()), palace::Operator::DIAG_ZERO);
    ctx->opA_pc = ctx->space_op->GetSystemMatrix(std::complex<double>(1.0, 0.0), ctx->lambda_test, ctx->lambda_test * ctx->lambda_test, ctx->opK, ctx->opC, ctx->opM, ctx->opA2_pc.get());
    ctx->opP_pc = ctx->space_op->GetPreconditionerMatrix<palace::ComplexOperator>(std::complex<double>(1.0, 0.0), ctx->lambda_test, ctx->lambda_test * ctx->lambda_test, ctx->lambda_test.imag());
    ctx->new_lambda = false;
  }

  // If I call ctx->opInv->SetOperators I get a segfault, so I'm creating a new ksp object here for testing
  // This is inefficient and expensive though, so we'd want to do something else instead...
  //std::cerr << "__pc_apply_NEP make unique ksp\n";
  // This is needed for Newton-based solvers:
    auto ksp = std::make_unique<palace::ComplexKspSolver>(*ctx->opIodata, ctx->space_op->GetNDSpaces(), &ctx->space_op->GetH1Spaces());
  ksp->SetOperators(*ctx->opA_pc, *ctx->opP_pc);
  ksp->Mult(ctx->x1, ctx->y1);
 /**/
  // For NLEIGS we don't want to change the operators
  //ctx->opInv->Mult(ctx->x1, ctx->y1);

  //std::cerr << "__pc_apply_NEP after Mult\n";
  /**/
  if (!ctx->sinvert)
  {
    ctx->y1 *= 1.0 / (ctx->delta * ctx->gamma);
  }
  else
  {
    ctx->y1 *= 1.0 / ctx->delta;
  }
  /**/
  if (ctx->opProj)
  {
    // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(ctx->GetComm(), ctx->y1));
    ctx->opProj->Mult(ctx->y1);
    // Mpi::Print(" After projection: {:e}\n", linalg::Norml2(ctx->GetComm(), ctx->y1));
  }

  PetscCall(ToPetscVec(ctx->y1, y));
std::cerr << "__pc_apply_NEP done \n";
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __form_NEP_function(NEP nep, PetscScalar lambda, Mat fun, Mat B, void *ctx)
{
  PetscFunctionBeginUser;
  std::cerr << "In SlepcNepSolver FormFunction with lambda: " << lambda << "\n";

  palace::slepc::SlepcNEPSolver *ctxF;
  PetscCall(MatShellGetContext(fun, (void **)&ctxF));

  // A(λ) = K + λ C + λ² M + A2(Im{λ}).
  ctxF->opA2 = ctxF->space_op->GetExtraSystemMatrix<palace::ComplexOperator>(
      std::abs(lambda.imag()), palace::Operator::DIAG_ZERO); //std:abs???
  ctxF->opA = ctxF->space_op->GetSystemMatrix(std::complex<double>(1.0, 0.0), lambda,
                                              lambda * lambda, ctxF->opK, ctxF->opC,
                                              ctxF->opM, ctxF->opA2.get());
  //ctxF->opP = ctxF->space_op->GetPreconditionerMatrix<palace::ComplexOperator>(std::complex<double>(1.0, 0.0), lambda, lambda * lambda, lambda.imag());

  ctxF->lambda_test = lambda;  // needed for duplication?
  ctxF->new_lambda = true;

  PetscFunctionReturn(PETSC_SUCCESS);
}

// May not need FormJacobian depending on choice on solver type
PetscErrorCode __form_NEP_jacobian(NEP nep, PetscScalar lambda, Mat fun, void *ctx)
{
  // A(λ) = K + λ C + λ² M + A2(Im{λ}).
  // A'(λ) = C + 2 λ M + A2'(Im{λ}.
  PetscFunctionBeginUser;
  std::cerr << "In SlepcNepSolver FormJacobian with lambda: " << lambda << "\n";
  palace::slepc::SlepcNEPSolver *ctxF;
  PetscCall(MatShellGetContext(fun, (void **)&ctxF));
  ctxF->opA2 = ctxF->space_op->GetExtraSystemMatrix<palace::ComplexOperator>(std::abs(lambda.imag()), palace::Operator::DIAG_ZERO);//std:abs???
  const auto eps = std::sqrt(std::numeric_limits<double>::epsilon());
  ctxF->opA2p = ctxF->space_op->GetExtraSystemMatrix<palace::ComplexOperator>(std::abs(lambda.imag()) * (1.0 + eps), palace::Operator::DIAG_ZERO);//std:abs???
  ctxF->opAJ = ctxF->space_op->GetExtraSystemMatrixJacobian<palace::ComplexOperator>(eps * std::abs(lambda.imag()), 1, ctxF->opA2p.get(), ctxF->opA2.get());//, ctxF->opA2.get()); // third operator not used?
  ctxF->opJ = ctxF->space_op->GetSystemMatrix(std::complex<double>(0.0, 0.0),
                                              std::complex<double>(1.0, 0.0), 2.0 * lambda,
                                              ctxF->opK, ctxF->opC, ctxF->opM, ctxF->opAJ.get());
  ctxF->lambda_J = lambda;      // needed for duplication?
  std::cerr << "Leaving SlepcNepSolver FormJacobian\n";
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __compute_singularities(NEP nep, PetscInt *maxnp, PetscScalar *xi, void *pt)
{
  PetscReal h;
  PetscInt i; // JUST A TEST
  PetscFunctionBeginUser;
  std::cerr << "In SlepcNepSolver compute_singularities maxnp: " << *maxnp << "\n";
  h = 11.0/(*maxnp-1);
   xi[0] = -1e-5; xi[*maxnp-1] = -1e+6;
  for (i=1;i<*maxnp-1;i++) xi[i] = -PetscPowReal(10,-5+h*i);
  PetscFunctionReturn(PETSC_SUCCESS);
}

#endif
