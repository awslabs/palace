// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "slepc.hpp"

#if defined(PALACE_WITH_SLEPC)

#include <algorithm>
#include <cmath>
#include <petsc.h>
#include <slepc.h>
#include <mfem.hpp>
#include "linalg/divfree.hpp"
#include "linalg/nleps.hpp"
#include "linalg/rap.hpp"
#include "linalg/solverstats.hpp"
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

PetscErrorCode __eps_region_monitor(EPS eps, PetscInt its, PetscInt nconv,
                                    PetscScalar *eigr, PetscScalar *eigi, PetscReal *errest,
                                    PetscInt nest, void *ctx)
{
  PetscFunctionBeginUser;
  PetscReal tol;
  PetscCall(EPSGetTolerances(eps, &tol, nullptr));
  static_cast<SlepcEigenvalueSolver *>(ctx)->CountRegionExcluded(tol, eigr, eigi, errest,
                                                                 static_cast<int>(nest));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __pep_region_monitor(PEP pep, PetscInt its, PetscInt nconv,
                                    PetscScalar *eigr, PetscScalar *eigi, PetscReal *errest,
                                    PetscInt nest, void *ctx)
{
  PetscFunctionBeginUser;
  PetscReal tol;
  PetscCall(PEPGetTolerances(pep, &tol, nullptr));
  static_cast<SlepcEigenvalueSolver *>(ctx)->CountRegionExcluded(tol, eigr, eigi, errest,
                                                                 static_cast<int>(nest));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __eps_relaxation_monitor(EPS eps, PetscInt its, PetscInt nconv,
                                        PetscScalar *eigr, PetscScalar *eigi,
                                        PetscReal *errest, PetscInt nest, void *ctx)
{
  PetscFunctionBeginUser;
  PetscInt nev;
  PetscCall(EPSGetDimensions(eps, &nev, nullptr, nullptr));
  static_cast<SlepcEPSSolverBase *>(ctx)->UpdateLinearSolverTol(
      static_cast<int>(nconv), static_cast<int>(nev), eigr, eigi, errest,
      static_cast<int>(nest));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __eps_divfree_monitor(EPS eps, PetscInt its, PetscInt nconv,
                                     PetscScalar *eigr, PetscScalar *eigi,
                                     PetscReal *errest, PetscInt nest, void *ctx)
{
  PetscFunctionBeginUser;
  PetscInt nev;
  PetscCall(EPSGetDimensions(eps, &nev, nullptr, nullptr));
  static_cast<SlepcEPSSolverBase *>(ctx)->CheckDivFreeTol(static_cast<int>(nev), eigr, eigi,
                                                          static_cast<int>(nest));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __eps_convergence_test(EPS eps, PetscScalar eigr, PetscScalar eigi,
                                      PetscReal res, PetscReal *errest, void *ctx)
{
  // The error estimate relative to the eigenvalue (as EPS_CONV_REL; converted to the
  // eigenvalue with a target apart from the shift), for a Ritz value which also settled
  // since the previous restart.
  PetscFunctionBeginUser;
  PetscCall(EPSConvergedRelative(eps, eigr, eigi, res, errest, nullptr));
  auto *solver = static_cast<SlepcEPSSolverBase *>(ctx);
  *errest *= solver->GetTargetErrorScale(eigr);
  solver->CheckRitzValue(eigr, errest);
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode __eps_stopping_test(EPS eps, PetscInt its, PetscInt max_it, PetscInt nconv,
                                   PetscInt nev, EPSConvergedReason *reason, void *ctx)
{
  // SLEPc's default stopping test, which also ends a convergence check of the Ritz values,
  // and stops a solve to be repeated with a smaller divergence-free projection tolerance.
  PetscFunctionBeginUser;
  PetscCall(EPSStoppingBasic(eps, its, max_it, nconv, nev, reason, nullptr));
  auto *solver = static_cast<SlepcEPSSolverBase *>(ctx);
  solver->EndRitzValueCheck();
  if (solver->StopForDivFreeTol())
  {
    *reason = EPS_CONVERGED_USER;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

// Whether the relative tolerance tol of the divergence-free projection of the
// shift-and-invert linear solve results is too large for the largest ratio
// r = |λ - σ| / |σ| of the wanted eigenvalues λ for the shift σ, and the tolerance which
// then replaces it (GetDivFreeSolveTol). Gradients are eigenvectors for λ = 0, so the
// shift-and-invert operator amplifies them r times more than a wanted eigenvector. A
// projection to the relative tolerance tol leaves about tol times the gradient components
// of each result, which the next operator application amplifies r times more than the
// wanted components: the gradient components of the Krylov basis stay negligible only for
// tol r << 1, and from tol r ~ 1 on they pollute the Ritz values, which then converge to
// spurious eigenpairs with small backward errors. A tolerance is kept while tol r <= 0.1;
// otherwise it becomes 0.01 / r, with a margin for Ritz values which move farther from σ,
// also below the tolerance of the projector (Solver.Linear.DivFreeTol), which only sets
// the starting one if larger than 1e-3.
bool DivFreeTolTooLarge(double tol, double r)
{
  constexpr double max_product = 0.1;
  return tol * r > max_product;
}

double GetDivFreeSolveTol(double tol, double r)
{
  constexpr double product = 0.01;
  return DivFreeTolTooLarge(tol, r) ? product / r : tol;
}

}  // namespace

void ConfigureQualityFactorRegion(RG rg, PetscReal s, double q_min)
{
  // For λ = iω, the quality factor of a damped mode Q = |ω| / (2 Im{ω}) >= q_min if
  // Re{ω} / (2 Im{ω}) = Im{λ} / (-2 Re{λ}) >= (q_min² - 1/4)^½. Growing modes (Re{λ} > 0)
  // are kept for |Q| >= q_grow, which includes lossless modes computed with a small error
  // of either sign. The region is a trapezoid from Im{λ} = s up to far above any eigenvalue
  // near the target.
  MFEM_VERIFY(s > 0.0 && q_min > 0.5,
              "Invalid target or minimum quality factor for the eigenvalue region!");
  constexpr double q_grow = 10.0;
  const PetscReal a = 0.5 / std::sqrt(q_min * q_min - 0.25);
  const PetscReal b = 0.5 / std::sqrt(q_grow * q_grow - 0.25);
  const PetscReal h = 1.0e3 * std::max(PetscReal(1.0), s);
  PetscScalar v[4] = {PetscScalar(-a * s, s), PetscScalar(b * s, s), PetscScalar(b * h, h),
                      PetscScalar(-a * h, h)};
  PalacePetscCall(RGSetType(rg, RGPOLYGON));
  PalacePetscCall(RGPolygonSetVertices(rg, 4, v, nullptr));
}

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

// Eigensolver base class methods.

SlepcEigenvalueSolver::SlepcEigenvalueSolver(int print) : print(print)
{
  sinvert = false;
  region = true;
  sigma = tau = 0.0;
  has_tau = st_shift_set = false;
  gamma = delta = 1.0;

  opInv = nullptr;
  opProj = nullptr;
  opB = nullptr;

  B0 = nullptr;
  v0 = nullptr;

  cl_custom = false;
  q_min = 0.0;
  q_region = q_region_monitor = false;
  q_region_count[0] = q_region_count[1] = q_region_count[2] = 0;
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
  sigma = tau = s;  // Wait until solve time to call EPS/PEPSetTarget
  has_tau = false;
  sinvert = true;
}

void SlepcEigenvalueSolver::SetTarget(std::complex<double> t)
{
  MFEM_VERIFY(sinvert && SupportsTarget(),
              "A target apart from the shift requires the shift-and-invert spectral "
              "transformation of a solver which supports it!");
  tau = t;
  has_tau = (tau != sigma);
}

PetscReal SlepcEigenvalueSolver::GetTargetErrorScale(PetscScalar theta) const
{
  if (!has_tau || theta == 0.0)
  {
    return 1.0;
  }
  // λ - σ = 1 / θ in the scaled variables of the solver.
  const PetscScalar l = sigma / gamma + 1.0 / theta;
  const PetscReal a = PetscAbsScalar(l);
  return (a > 0.0) ? std::max(PetscReal(1.0), 1.0 / (PetscAbsScalar(theta) * a)) : 1.0;
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
    // With refinement, a breakdown (linear dependence) is flagged when the second pass
    // reduces the norm below eta times that after the first one: eta = 1 would flag
    // almost every vector of a nearly invariant Krylov space (for example after a warm
    // start), so use the default eta = 1/√2.
    PalacePetscCall(
        BVSetOrthogonalization(bv, type, refine, PETSC_DETERMINE, BV_ORTHOG_BLOCK_GS));
  }
}

void SlepcEigenvalueSolver::CountRegionExcluded(PetscReal tol, const PetscScalar *eigr,
                                                const PetscScalar *eigi,
                                                const PetscReal *errest, int nest)
{
  // The region applies to the back-transformed eigenvalues (as in the sorting criterion of
  // the solver).
  ST st = GetST();
  RG rg = GetRG();
  q_region_count[0] = nest;
  q_region_count[1] = q_region_count[2] = 0;
  for (int i = 0; i < nest; i++)
  {
    PetscScalar er = eigr[i], ei = eigi[i];
    PetscInt inside;
    PalacePetscCall(STBackTransform(st, 1, &er, &ei));
    PalacePetscCall(RGCheckInside(rg, 1, &er, &ei, &inside));
    if (inside < 0)
    {
      q_region_count[1]++;
      q_region_count[2] += (errest[i] < tol);
    }
  }
}

void SlepcEigenvalueSolver::PrintRegionExcluded() const
{
  Mpi::Print(GetComm(),
             " Search region Q >= {:.3g}: {:d} of {:d} Ritz values outside in the last "
             "iteration ({:d} with error estimates below the tolerance)\n",
             q_min, q_region_count[1], q_region_count[0], q_region_count[2]);
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
  q_region = false;
  if (sinvert && region)
  {
    if (PetscImaginaryPart(tau) == 0.0)
    {
      PetscReal sr = PetscRealPart(tau);
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
    else if (PetscRealPart(tau) == 0.0)
    {
      PetscReal si = PetscImaginaryPart(tau);
      if (si > 0.0 && q_min > 0.5)
      {
        ConfigureQualityFactorRegion(GetRG(), si / gamma, q_min);
        q_region = true;
      }
      else if (si > 0.0)
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
  relax_c = relax_lin_tol = relax_theta_max = 0.0;
  relax_monitor = false;
}

SlepcEPSSolverBase::~SlepcEPSSolverBase()
{
  PalacePetscCall(EPSDestroy(&eps));
  PalacePetscCall(MatDestroy(&A0));
  PalacePetscCall(MatDestroy(&A1));
  PalacePetscCall(VecDestroy(&v_init));
}

void SlepcEPSSolverBase::SaveInitialSpace(Vec v)
{
  if (!v_init)
  {
    PalacePetscCall(VecDuplicate(v, &v_init));
  }
  PalacePetscCall(VecCopy(v, v_init));
  has_init = true;
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
  SaveInitialSpace(is[0]);
}

void SlepcEPSSolverBase::ResetInitialSpace()
{
  if (has_init)
  {
    Vec is[1] = {v_init};
    PalacePetscCall(EPSSetInitialSpace(eps, 1, is));
  }
}

void SlepcEPSSolverBase::Customize()
{
  SlepcEigenvalueSolver::Customize();
  PalacePetscCall(EPSSetTarget(eps, tau / gamma));
  if (has_tau || st_shift_set)
  {
    // The target only sets the default shift of the spectral transformation.
    PalacePetscCall(STSetShift(GetST(), sigma / gamma));
    st_shift_set = true;
  }
  if (!cl_custom)
  {
#if defined(SLEPC_HAVE_EPS_KRYLOVSCHUR_CHECK_INTERVAL)
    // Convergence checks within the Krylov-Schur cycles (SLEPc patch, see
    // SetConvergenceCheckInterval); command-line options
    // (-eps_krylovschur_check_interval) override it.
    PalacePetscCall(EPSKrylovSchurSetCheckInterval(eps, check_interval));
#endif
    PalacePetscCall(EPSSetFromOptions(eps));
    if (print > 0)
    {
      PetscOptionsView(nullptr, PETSC_VIEWER_STDOUT_(GetComm()));
      Mpi::Print(GetComm(), "\n");
    }
    cl_custom = true;
  }
}

void SlepcEPSSolverBase::CheckRitzValue(PetscScalar theta, PetscReal *errest)
{
  // SLEPc's error estimate for a Ritz pair (θ, z) of the transformed problem assumes exact
  // applications of the spectral transformation, so it does not see the error of iterative
  // linear solves in the Ritz value, which can exceed the tolerance by orders of magnitude
  // for an eigenpair converging as its Ritz value appears: after a restart which purged the
  // Ritz values outside of the search region, or in the first Krylov cycle of a
  // warm-started solve. A Ritz value with an error estimate below the tolerance converges
  // only if the nearest Ritz value of the previous restart is the same eigenpair: within a
  // relative distance max(Tol, min(1e-3, 0.1 d₂)) of it, for the relative distance d₂ of
  // the second nearest one. Otherwise (or without a previous restart, except in the first
  // Krylov cycle of a solve which is not warm-started) its error estimate becomes that
  // relative distance, unless it was kept unconverged for that reason at the two previous
  // restarts already. The convergence test is called at each restart for every unconverged
  // Ritz value, in order (EPSSetTrackAll), and only those before the first unconverged one
  // converge. For well-separated Ritz values the distance allowed is 1e-3: the test holds
  // back the Ritz values which just appeared, not ones which move slowly between restarts.
  constexpr double max_distance = 1.0e-3, max_fraction = 0.1;
  constexpr int max_unsettled = 2;
  PetscInt its;
  PalacePetscCall(EPSGetIterationNumber(eps, &its));
  if (its != ritz_its)
  {
    ritz_its = its;
    ritz_prev.swap(ritz_last);
    ritz_last.clear();
    ritz_curr.clear();
    ritz_blocked = false;
  }
  PetscReal tol;
  PalacePetscCall(EPSGetTolerances(eps, &tol, nullptr));
  int unsettled = 0;
  if (*errest < tol && (its > 1 || warm_start))
  {
    PetscReal d1 = mfem::infinity(), d2 = mfem::infinity();
    int count = 0;
    for (const auto &[theta_prev, count_prev] : ritz_prev)
    {
      const PetscReal d = PetscAbsScalar(theta - theta_prev);
      if (d < d1)
      {
        d2 = d1;
        d1 = d;
        count = count_prev;
      }
      else if (d < d2)
      {
        d2 = d;
      }
    }
    const PetscReal a = PetscAbsScalar(theta);
    d1 /= a;
    d2 /= a;
    if (d1 > std::max(static_cast<double>(tol), std::min(max_distance, max_fraction * d2)))
    {
      if (count < max_unsettled)
      {
        *errest = std::min(d1, PetscReal(1.0));
        unsettled = count + 1;
        ritz_num_delayed += !ritz_blocked;
      }
      else
      {
        if (!ritz_blocked)
        {
          ritz_num_unsettled++;
          ritz_unsettled.push_back(theta);
        }
      }
    }
  }
  ritz_blocked = ritz_blocked || !(*errest < tol);
  ritz_curr.emplace_back(theta, unsettled);
}

void SlepcEPSSolverBase::EndRitzValueCheck()
{
  ritz_last.swap(ritz_curr);
  ritz_curr.clear();
  ritz_blocked = false;
}

void SlepcEPSSolverBase::UpdateLinearSolverTol(int nconv, int nev, const PetscScalar *eigr,
                                               const PetscScalar *eigi,
                                               const PetscReal *errest, int nest)
{
  // Smallest relaxed tolerance over the wanted unconverged eigenpairs, from their error
  // estimates, the ratio of the largest transformed eigenvalue so far to theirs, and their
  // relative separation from the other Ritz values.
  PetscReal eig_tol;
  PalacePetscCall(EPSGetTolerances(eps, &eig_tol, nullptr));
  for (int i = 0; i < nest; i++)
  {
    relax_theta_max = std::max(relax_theta_max, SlepcAbsEigenvalue(eigr[i], eigi[i]));
  }
  double tol = mfem::infinity(), r_max = 0.0, ratio_max = 1.0;
  for (int k = nconv; k < std::min(nev, nest); k++)
  {
    const PetscReal theta = SlepcAbsEigenvalue(eigr[k], eigi[k]);
    PetscReal gap = mfem::infinity();
    for (int j = 0; j < nest; j++)
    {
      if (j != k)
      {
        gap = std::min(gap, SlepcAbsEigenvalue(eigr[j] - eigr[k], eigi[j] - eigi[k]));
      }
    }
    const double rel_gap = (theta > 0.0) ? gap / theta : 0.0;
    const double ratio = (theta > 0.0) ? relax_theta_max / theta : mfem::infinity();
    tol = std::min(tol, GetRelaxedLinearSolverTol(relax_c, eig_tol, relax_lin_tol,
                                                  errest[k], ratio, rel_gap));
    r_max = std::max(r_max, static_cast<double>(errest[k]));
    ratio_max = std::max(ratio_max, ratio);
  }
  if (!(tol < mfem::infinity()) || relax_suspend)
  {
    tol = relax_lin_tol;
  }
  relax_max_used = std::max(relax_max_used, tol);
  if (tol != opInv->GetRelTol())
  {
    opInv->SetRelTol(tol);
    if (print > 0)
    {
      Mpi::Print(GetComm(),
                 " Linear solver tolerance {:.3e} (wanted error estimates <= {:.3e}, "
                 "eigenvalue ratio <= {:.1f})\n",
                 tol, r_max, ratio_max);
    }
  }
}

void SlepcEPSSolverBase::CheckDivFreeTol(int nev, const PetscScalar *eigr,
                                         const PetscScalar *eigi, int nest)
{
  // The gradients have the transformed eigenvalue -1 / σ, so the ratio |λ - σ| / |σ| of a
  // Ritz value θ of the transformed problem is 1 / (|θ| |σ|), in the scaled variables. Once
  // the ratio of a wanted one requires a smaller tolerance (see GetDivFreeSolveTol), the
  // solve stops at the end of the next Krylov cycle, which skips its linear solves, unless
  // it is the last one allowed. The maximum over the processes makes the decision the
  // same on all of them.
  if (divfree_stop > 0.0)
  {
    return;
  }
  const PetscReal s = PetscAbsScalar(sigma) / gamma;
  double r = 0.0;
  for (int k = 0; k < std::min(nev, nest); k++)
  {
    const PetscReal theta = SlepcAbsEigenvalue(eigr[k], eigi[k]);
    r = std::max(r, (theta * s > 0.0) ? 1.0 / (theta * s) : mfem::infinity());
  }
  Mpi::GlobalMax(1, &r, GetComm());
  divfree_ratio = std::max(divfree_ratio, r);
  divfree_ratio_last = r;
  if (!divfree_final && DivFreeTolTooLarge(divfree_tol, r))
  {
    divfree_stop = r;
  }
}

int SlepcEPSSolverBase::Solve()
{
  MFEM_VERIFY(A0 && A1 && opInv, "Operators are not set for SlepcEPSSolverBase!");

  // Solve the eigenvalue problem. With the linear solver relaxation, a monitor updates the
  // linear solver tolerance after each restart from the error estimates of all Ritz values
  // (EPSSetTrackAll, which does not change the convergence test), and the tolerance is
  // restored after the solve.
  PetscInt num_conv;
  Customize();
  relax_lin_tol = opInv->GetRelTol();
  relax_theta_max = 0.0;
  // With a tighter tolerance for the first Krylov cycle, the monitor restores the linear
  // solver tolerance (or the relaxed one) after the first restart.
  const bool first_cycle = (first_cycle_tol > 0.0 && first_cycle_tol < relax_lin_tol);
  if ((relax_c > 0.0 || first_cycle) && !relax_monitor)
  {
    PalacePetscCall(EPSMonitorSet(eps, __eps_relaxation_monitor, (void *)this, nullptr));
    PalacePetscCall(EPSSetTrackAll(eps, PETSC_TRUE));
    relax_monitor = true;
  }
  if (q_region && print > 0 && !q_region_monitor)
  {
    PalacePetscCall(EPSMonitorSet(eps, __eps_region_monitor, (void *)this, nullptr));
    PalacePetscCall(EPSSetTrackAll(eps, PETSC_TRUE));
    q_region_monitor = true;
  }
  q_region_count[0] = q_region_count[1] = q_region_count[2] = 0;

  // Tolerance of the divergence-free projection of the linear solve results: 1e-3 (or the
  // one of the projector, if larger), or a smaller one if the eigenvalue estimates require
  // it, and a monitor checks it against the wanted Ritz values after each restart (see
  // CheckDivFreeTol and GetDivFreeSolveTol).
  const bool divfree_check = (opProj && sinvert);
  if (divfree_check && !divfree_monitor)
  {
    PalacePetscCall(EPSMonitorSet(eps, __eps_divfree_monitor, (void *)this, nullptr));
    divfree_monitor = true;
  }
  PetscInt divfree_nev;
  PalacePetscCall(EPSGetDimensions(eps, &divfree_nev, nullptr, nullptr));
  auto GetDivFreeRatio = [this](std::complex<double> l)
  {
    const std::complex<double> s(sigma);
    return (std::abs(s) > 0.0) ? std::abs(l - s) / std::abs(s) : mfem::infinity();
  };
  divfree_tol = opProj ? std::max(opProj->GetTol(), 1.0e-3) : 1.0e-3;
  divfree_stop = divfree_ratio = divfree_ratio_last = 0.0;
  divfree_repeats = 0;
  divfree_final = false;
  if (divfree_check && !divfree_estimates.empty())
  {
    const std::size_t n = std::min(divfree_estimates.size(), std::size_t(divfree_nev));
    double r = 0.0;
    for (std::size_t i = 0; i < n; i++)
    {
      r = std::max(r, GetDivFreeRatio(divfree_estimates[i]));
    }
    const double tol0 = divfree_tol;
    divfree_tol = GetDivFreeSolveTol(divfree_tol, r);
    if (divfree_tol < tol0)
    {
      Mpi::Print(GetComm(),
                 " Divergence-free projection tolerance {:.1e} for the linear solve "
                 "results (eigenvalue estimates at |λ - σ| / |σ| <= {:.3g})\n",
                 divfree_tol, r);
    }
  }
  PalacePetscCall(
      EPSSetConvergenceTestFunction(eps, __eps_convergence_test, (void *)this, nullptr));
  PalacePetscCall(
      EPSSetStoppingTestFunction(eps, __eps_stopping_test, (void *)this, nullptr));
  PalacePetscCall(EPSSetTrackAll(eps, PETSC_TRUE));
  ritz_prev.clear();
  ritz_last.clear();
  ritz_curr.clear();
  ritz_unsettled.clear();
  ritz_its = -1;
  ritz_num_delayed = ritz_num_unsettled = 0;
  relax_max_used = relax_lin_tol;
  relax_suspend = false;
  // A relaxed linear solve which converges falsely is repeated at the strict tolerance.
  opInv->SetFallbackRelTol((relax_c > 0.0) ? relax_lin_tol : 0.0);
  if (first_cycle)
  {
    opInv->SetRelTol(first_cycle_tol);
    if (print > 0)
    {
      Mpi::Print(GetComm(), " Linear solver tolerance {:.3e} for the first Krylov cycle\n",
                 first_cycle_tol);
    }
  }

  // A solve stopped by the divergence-free projection check is repeated from the same
  // initial space with the tolerance for the ratio |λ - σ| / |σ| which stopped it, at most
  // twice: the second repeat also takes the tolerance of the projector
  // (Solver.Linear.DivFreeTol) if smaller, and is not stopped (see the warning at the end
  // of the solve). Returns the number of iterations of the stopped solves.
  auto RepeatForDivFreeTol = [&]()
  {
    constexpr int max_repeats = 2;
    PetscInt stopped_its = 0;
    while (divfree_stop > 0.0)
    {
      PetscInt its;
      PalacePetscCall(EPSGetIterationNumber(eps, &its));
      stopped_its += its;
      const double r = divfree_stop;
      divfree_repeats++;
      divfree_final = (divfree_repeats >= max_repeats);
      divfree_tol = GetDivFreeSolveTol(divfree_tol, r);
      if (divfree_final)
      {
        divfree_tol = std::min(divfree_tol, opProj->GetTol());
      }
      Mpi::Print(GetComm(),
                 "\n Repeating the eigenvalue solve with the divergence-free projection "
                 "tolerance {:.1e} (wanted Ritz values at |λ - σ| / |σ| = {:.3g} after "
                 "{:d} restart{})\n",
                 divfree_tol, r, its, (its == 1) ? "" : "s");
      if (has_init)
      {
        Vec is[1] = {v_init};
        PalacePetscCall(EPSSetInitialSpace(eps, 1, is));
      }
      ritz_prev.clear();
      ritz_last.clear();
      ritz_curr.clear();
      ritz_unsettled.clear();
      ritz_its = -1;
      ritz_num_delayed = ritz_num_unsettled = 0;
      relax_theta_max = 0.0;
      relax_max_used = relax_lin_tol;
      divfree_stop = divfree_ratio = divfree_ratio_last = 0.0;
      opInv->SetRelTol(first_cycle ? first_cycle_tol : relax_lin_tol);
      PalacePetscCall(EPSSolve(eps));
    }
    return stopped_its;
  };
  PalacePetscCall(EPSSolve(eps));
  const PetscInt divfree_its = RepeatForDivFreeTol();
  opInv->SetRelTol(relax_lin_tol);
  PalacePetscCall(EPSGetConverged(eps, &num_conv));
  {
    PetscInt num_it, num_req;
    PalacePetscCall(EPSGetIterationNumber(eps, &num_it));
    PalacePetscCall(EPSGetDimensions(eps, &num_req, nullptr, nullptr));
    perf::RecordEigenSolve(opInv->GetRole(), static_cast<int>(num_it + divfree_its),
                           static_cast<int>(num_conv), static_cast<int>(num_req));
  }
  if (print > 0)
  {
    Mpi::Print(GetComm(), "\n");
    PalacePetscCall(EPSConvergedReasonView(eps, PETSC_VIEWER_STDOUT_(GetComm())));
    Mpi::Print(GetComm(),
               " Total number of linear systems solved: {:d}\n"
               " Total number of linear solver iterations: {:d}\n",
               opInv->NumTotalMult(), opInv->NumTotalMultIterations());
    if (q_region && q_region_monitor)
    {
      PrintRegionExcluded();
    }
    PrintRitzValueChecks();
  }

  // Compute and store the eigenpair residuals.
  RescaleEigenvectors(num_conv);

  // A solve with relaxed linear solves is repeated from the same initial space with the
  // strict linear solver tolerance if fewer eigenpairs than requested converged, if a
  // requested eigenpair was accepted with a Ritz value which did not settle between
  // restarts (see CheckRitzValue), or if one has a backward error above twice the tolerance
  // scaled by |λ - σ| / |σ|. The gap-aware bound of the relaxation keeps the error of the
  // others at the level of the strict solves, whose backward errors themselves can exceed
  // the tolerance several times.
  if (relax_c > 0.0 && relax_max_used > relax_lin_tol)
  {
    PetscInt nev;
    PetscReal eig_tol;
    PalacePetscCall(EPSGetDimensions(eps, &nev, nullptr, nullptr));
    PalacePetscCall(EPSGetTolerances(eps, &eig_tol, nullptr));
    PetscInt num_unsettled = 0;
    double worst = 0.0;
    for (PetscInt i = 0; i < std::min(nev, num_conv); i++)
    {
      // The backward error of strict solves grows with the distance of the eigenvalue from
      // the shift relative to the shift (the residuals of the linear solves enter the
      // residual of the eigenpair amplified by about |λ - σ| / |σ|): only a backward error
      // beyond twice that times the tolerance is a sign of relaxed solves gone wrong.
      const std::complex<double> l = GetEigenvalue(static_cast<int>(i));
      const double amplification =
          (std::abs(sigma) > 0.0)
              ? std::max(1.0, std::abs(l - std::complex<double>(sigma)) / std::abs(sigma))
              : 1.0;
      worst = std::max(worst, GetError(static_cast<int>(i), ErrorType::BACKWARD) /
                                  (amplification * eig_tol));
    }
    constexpr double max_error = 2.0;
    if (num_conv >= nev && !ritz_unsettled.empty())
    {
      PetscScalar target;
      PalacePetscCall(EPSGetTarget(eps, &target));
      for (PetscInt i = 0; i < nev; i++)
      {
        PetscScalar l;
        PalacePetscCall(EPSGetEigenvalue(eps, i, &l, nullptr));
        const PetscScalar theta = 1.0 / (l - target);
        for (const auto &theta_k : ritz_unsettled)
        {
          if (PetscAbsScalar(theta - theta_k) <= 1.0e-8 * PetscAbsScalar(theta_k))
          {
            num_unsettled++;
            break;
          }
        }
      }
    }
    if (num_conv < nev || num_unsettled > 0 || worst > max_error)
    {
      Mpi::Print(GetComm(),
                 " Repeating the eigenvalue solve with the linear solver tolerance {:.1e} "
                 "(after relaxed linear solves, {})\n",
                 relax_lin_tol,
                 (num_conv < nev)
                     ? fmt::format("{:d} of {:d} eigenpairs converged", num_conv, nev)
                 : (num_unsettled > 0)
                     ? fmt::format("{:d} requested eigenpair{} not settled", num_unsettled,
                                   (num_unsettled == 1) ? "" : "s")
                     : fmt::format("backward error {:.1f} times the tolerance scaled by "
                                   "|λ - σ| / |σ|",
                                   worst));
      if (has_init)
      {
        Vec is[1] = {v_init};
        PalacePetscCall(EPSSetInitialSpace(eps, 1, is));
      }
      const int mult0 = opInv->NumTotalMult(), it0 = opInv->NumTotalMultIterations();
      relax_suspend = true;
      ritz_prev.clear();
      ritz_last.clear();
      ritz_curr.clear();
      ritz_unsettled.clear();
      ritz_its = -1;
      ritz_num_delayed = ritz_num_unsettled = 0;
      if (first_cycle)
      {
        opInv->SetRelTol(first_cycle_tol);
      }
      PalacePetscCall(EPSSolve(eps));
      RepeatForDivFreeTol();
      relax_suspend = false;
      opInv->SetRelTol(relax_lin_tol);
      PalacePetscCall(EPSGetConverged(eps, &num_conv));
      if (print > 0)
      {
        Mpi::Print(GetComm(), "\n");
        PalacePetscCall(EPSConvergedReasonView(eps, PETSC_VIEWER_STDOUT_(GetComm())));
        Mpi::Print(GetComm(),
                   " Total number of linear systems solved: {:d}\n"
                   " Total number of linear solver iterations: {:d}\n",
                   opInv->NumTotalMult() - mult0, opInv->NumTotalMultIterations() - it0);
        PrintRitzValueChecks();
      }
      RescaleEigenvectors(num_conv);
    }
  }
  opInv->SetFallbackRelTol(0.0);
  if (divfree_check)
  {
    if (print > 0)
    {
      Mpi::Print(GetComm(),
                 " Divergence-free projection tolerance {:.1e} for the linear solve "
                 "results (wanted Ritz values at |λ - σ| / |σ| <= {:.3g})\n",
                 divfree_tol, divfree_ratio);
    }
    if (DivFreeTolTooLarge(divfree_tol, divfree_ratio_last))
    {
      // Only the last solve allowed can end with a tolerance too large for its eigenvalues.
      Mpi::Warning(GetComm(),
                   "The divergence-free projection tolerance {:.1e} of the eigenvalue "
                   "solve is too large for its wanted eigenvalues at |λ - σ| / |σ| = "
                   "{:.3g} (after {:d} repeated solves): gradient components may pollute "
                   "the eigenpairs!\nDecrease "
                   "config[\"Solver\"][\"Linear\"][\"DivFreeTol\"], or move the target "
                   "closer to the eigenvalues.\n",
                   divfree_tol, divfree_ratio_last, divfree_repeats);
    }
    if (num_conv > 0)
    {
      divfree_estimates.clear();
      for (PetscInt i = 0; i < std::min(divfree_nev, num_conv); i++)
      {
        divfree_estimates.push_back(GetEigenvalue(static_cast<int>(i)));
      }
    }
  }
  return static_cast<int>(num_conv);
}

void SlepcEPSSolverBase::PrintRitzValueChecks() const
{
  if (ritz_num_delayed > 0 || ritz_num_unsettled > 0)
  {
    Mpi::Print(GetComm(),
               " Convergence delayed {:d} time{} for Ritz values not settled since the "
               "previous restart{}\n",
               ritz_num_delayed, (ritz_num_delayed == 1) ? "" : "s",
               (ritz_num_unsettled > 0)
                   ? fmt::format(" ({:d} accepted unsettled)", ritz_num_unsettled)
                   : "");
  }
}

int SlepcEPSSolverBase::GetNumIterations() const
{
  PetscInt num_it;
  PalacePetscCall(EPSGetIterationNumber(eps, &num_it));
  return static_cast<int>(num_it);
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
  SaveInitialSpace(is[0]);
}

void SlepcPEPLinearSolver::SetInitialSpaceBlocks(const ComplexVector &v,
                                                 const ComplexVector &v_lambda)
{
  MFEM_VERIFY(
      A0 && A1,
      "Must call SetOperators before using SetInitialSpace for SLEPc eigenvalue solver!");
  if (!v0)
  {
    PalacePetscCall(MatCreateVecs(A0, nullptr, &v0));
  }
  // The second block of the eigenvectors of the linearization is (λ / γ) x.
  ComplexVector v2(v_lambda);
  v2 *= 1.0 / gamma;
  PalacePetscCall(ToPetscVec(v, v2, v0));
  Vec is[1] = {v0};
  PalacePetscCall(EPSSetInitialSpace(eps, 1, is));
  SaveInitialSpace(is[0]);
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
  MFEM_VERIFY(!has_tau, "PEP solvers do not support a target apart from the shift!");
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
  if (q_region && print > 0 && !q_region_monitor)
  {
    PalacePetscCall(PEPMonitorSet(pep, __pep_region_monitor, (void *)this, nullptr));
    PalacePetscCall(PEPSetTrackAll(pep, PETSC_TRUE));
    q_region_monitor = true;
  }
  q_region_count[0] = q_region_count[1] = q_region_count[2] = 0;
  PalacePetscCall(PEPSolve(pep));
  PalacePetscCall(PEPGetConverged(pep, &num_conv));
  {
    PetscInt num_it, num_req;
    PalacePetscCall(PEPGetIterationNumber(pep, &num_it));
    PalacePetscCall(PEPGetDimensions(pep, &num_req, nullptr, nullptr));
    perf::RecordEigenSolve(opInv->GetRole(), static_cast<int>(num_it),
                           static_cast<int>(num_conv), static_cast<int>(num_req));
  }
  if (print > 0)
  {
    Mpi::Print(GetComm(), "\n");
    PalacePetscCall(PEPConvergedReasonView(pep, PETSC_VIEWER_STDOUT_(GetComm())));
    Mpi::Print(GetComm(),
               " Total number of linear systems solved: {:d}\n"
               " Total number of linear solver iterations: {:d}\n",
               opInv->NumTotalMult(), opInv->NumTotalMultIterations());
    if (q_region && q_region_monitor)
    {
      PrintRegionExcluded();
    }
  }

  // Compute and store the eigenpair residuals.
  RescaleEigenvectors(num_conv);
  return static_cast<int>(num_conv);
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
  sigma = tau = s;  // Wait until solve time to call NEPSetTarget
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
  {
    PetscInt num_it, num_req;
    PalacePetscCall(NEPGetIterationNumber(nep, &num_it));
    PalacePetscCall(NEPGetDimensions(nep, &num_req, nullptr, nullptr));
    perf::RecordEigenSolve(opInv->GetRole(), static_cast<int>(num_it),
                           static_cast<int>(num_conv), static_cast<int>(num_req));
  }
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
    std::function<std::unique_ptr<ComplexOperator>(std::complex<double>)> A2)
{
  funcA2 = A2;
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
  if (funcA2)
  {
    auto A2 = (*funcA2)(l);
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

  if (ctx->divfree_stop > 0.0)
  {
    // The solve stops at the end of the Krylov cycle, to be repeated with a smaller
    // divergence-free projection tolerance: skip its remaining linear solves.
    PetscCall(VecCopy(x, y));
    PetscFunctionReturn(PETSC_SUCCESS);
  }
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
    ctx->opProj->Mult(ctx->y1, ctx->divfree_tol);
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

  if (ctx->divfree_stop > 0.0)
  {
    // The solve stops at the end of the Krylov cycle (see __pc_apply_EPS).
    PetscCall(VecCopy(x, y));
    PetscFunctionReturn(PETSC_SUCCESS);
  }
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
      ctx->opProj->Mult(ctx->y2, ctx->divfree_tol);
    }

    ctx->y1.AXPBYPCZ(ctx->gamma / ctx->sigma, ctx->y2, -ctx->gamma / ctx->sigma, ctx->x1,
                     0.0);
    if (ctx->opProj)
    {
      // Without this projection, gradients are eigenvectors with λ = 0, which dominate the
      // targeted eigenvalues when these are farther from σ than σ from 0. The divergence of
      // y1 only comes from those of y2 and of the projected basis, and the tolerance is the
      // same as for y2 (see GetDivFreeSolveTol).
      ctx->opProj->Mult(ctx->y1, ctx->divfree_tol);
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
    ctx->opA2_pc = (*ctx->funcA2)(ctx->lambda);
    ctx->opA_pc =
        palace::BuildOperatorWithA2({1.0 + 0.0i, ctx->lambda, ctx->lambda * ctx->lambda},
                                    {ctx->opK, ctx->opC, ctx->opM}, ctx->opA2_pc.get());
    ctx->opP_pc = (*ctx->funcP)(std::complex<double>(1.0, 0.0), ctx->lambda,
                                ctx->lambda * ctx->lambda,
                                ctx->lambda / std::complex<double>(0.0, 1.0));
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
  // A(λ) = K + λ C + λ² M + A2(λ).
  ctxF->opA2 = (*ctxF->funcA2)(lambda);
  ctxF->opA =
      palace::BuildOperatorWithA2({1.0 + 0.0i, lambda, lambda * lambda},
                                  {ctxF->opK, ctxF->opC, ctxF->opM}, ctxF->opA2.get());
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
  // The A2 finite difference is a matrix-free sum (each A2 may itself be a matrix-free
  // SumComplexOperator carrying the wave-port correction).
  ctxF->opA2_jac = (*ctxF->funcA2)(lambda);
  const auto eps = std::sqrt(std::numeric_limits<double>::epsilon());
  ctxF->opA2p = (*ctxF->funcA2)(lambda * (1.0 + eps));
  std::complex<double> denom = eps * lambda;
  auto opAJ = std::make_unique<palace::SumComplexOperator>(ctxF->opA2p->Height(),
                                                           ctxF->opA2p->Width());
  opAJ->AddOperator(*ctxF->opA2p, 1.0 / denom);
  opAJ->AddOperator(*ctxF->opA2_jac, -1.0 / denom);
  ctxF->opAJ = std::move(opAJ);
  ctxF->opJ =
      palace::BuildOperatorWithA2({0.0 + 0.0i, 1.0 + 0.0i, 2.0 * lambda},
                                  {ctxF->opK, ctxF->opC, ctxF->opM}, ctxF->opAJ.get());
  PetscFunctionReturn(PETSC_SUCCESS);
}

#endif
