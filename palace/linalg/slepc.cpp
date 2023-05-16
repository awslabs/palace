// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "slepc.hpp"

#if defined(PALACE_WITH_SLEPC)

#include <algorithm>
#include <petsc.h>
#include <slepc.h>
#include <mfem.hpp>
#include "linalg/divfree.hpp"
#include "linalg/ksp.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

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

PetscReal GetMaxSingularValue(MPI_Comm comm, const ComplexOperator &A, bool herm,
                              PetscReal tol, PetscInt max_it)
{
  // This method assumes the provided operator has the required operations for SLEPc's EPS
  // or SVD solvers, namely MATOP_MULT and MATOP_MULT_HERMITIAN_TRANSPOSE (if the matrix
  // is not Hermitian).
  Mat A0;
  PetscInt n = A.Height() / 2;
  PalacePetscCall(MatCreateShell(comm, n, n, PETSC_DECIDE, PETSC_DECIDE, nullptr, &A0));
  ComplexVector x(A.Height()), y(A.Height());
  auto __mat_apply_shell = [&A, &x, &y](Mat, Vec x0, Vec y0) -> PetscErrorCode
  {
    PetscFunctionBeginUser;
    PetscInt n;
    PetscCall(VecGetLocalSize(x0, &n));

    const PetscScalar *px0;
    PetscCall(VecGetArrayRead(x0, &px0));
    x.Set(px0, n);
    PetscCall(VecRestoreArrayRead(x0, &px0));

    A.Mult(x, y);

    PetscScalar *py0;
    PetscCall(VecGetArrayWrite(y0, &py0));
    y.Get(py0, n);
    PetscCall(VecRestoreArrayWrite(y0, &py0));

    PetscFunctionReturn(0);
  };
  auto __mat_apply_transpose_shell = [&A, &x, &y](Mat, Vec x0, Vec y0) -> PetscErrorCode
  {
    PetscFunctionBeginUser;
    PetscInt n;
    PetscCall(VecGetLocalSize(x0, &n));

    const PetscScalar *px0;
    PetscCall(VecGetArrayRead(x0, &px0));
    x.Set(px0, n);
    PetscCall(VecRestoreArrayRead(x0, &px0));

    A.MultTranspose(x, y);

    PetscScalar *py0;
    PetscCall(VecGetArrayWrite(y0, &py0));
    y.Get(py0, n);
    PetscCall(VecRestoreArrayWrite(y0, &py0));

    PetscFunctionReturn(0);
  };
  auto __mat_apply_hermitian_transpose_shell = [&A, &x, &y](Mat, Vec x0,
                                                            Vec y0) -> PetscErrorCode
  {
    PetscFunctionBeginUser;
    PetscInt n;
    PetscCall(VecGetLocalSize(x0, &n));

    const PetscScalar *px0;
    PetscCall(VecGetArrayRead(x0, &px0));
    x.Set(px0, n);
    PetscCall(VecRestoreArrayRead(x0, &px0));

    A.MultHermitianTranspose(x, y);

    PetscScalar *py0;
    PetscCall(VecGetArrayWrite(y0, &py0));
    y.Get(py0, n);
    PetscCall(VecRestoreArrayWrite(y0, &py0));

    PetscFunctionReturn(0);
  };
  PalacePetscCall(MatShellSetOperation(A0, MATOP_MULT, (void (*)()) & __mat_apply_shell));

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
                                         (void (*)()) & __mat_apply_transpose_shell));
    PalacePetscCall(
        MatShellSetOperation(A0, MATOP_MULT_HERMITIAN_TRANSPOSE,
                             (void (*)()) & __mat_apply_hermitian_transpose_shell));

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

SlepcEigenSolver::SlepcEigenSolver(int print) : print(print)
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

SlepcEigenSolver::~SlepcEigenSolver()
{
  PalacePetscCall(MatDestroy(&B0));
  PalacePetscCall(VecDestroy(&v0));
}

void SlepcEigenSolver::SetOperators(const ComplexOperator &K, const ComplexOperator &M,
                                    EigenvalueSolver::ScaleType type)
{
  MFEM_ABORT("SetOperators not defined for base class SlepcEigenSolver!");
}

void SlepcEigenSolver::SetOperators(const ComplexOperator &K, const ComplexOperator &C,
                                    const ComplexOperator &M,
                                    EigenvalueSolver::ScaleType type)
{
  MFEM_ABORT("SetOperators not defined for base class SlepcEigenSolver!");
}

void SlepcEigenSolver::SetLinearSolver(const ComplexKspSolver &ksp)
{
  opInv = &ksp;
}

void SlepcEigenSolver::SetDivFreeProjector(const DivFreeSolver &divfree)
{
  opProj = &divfree;
}

void SlepcEigenSolver::SetBMat(const Operator &B)
{
  opB = &B;
}

void SlepcEigenSolver::SetShiftInvert(PetscScalar s, bool precond)
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
        SetRegion(PetscRealPart(sigma) / gamma - PETSC_i * mfem::infinity(),
                  mfem::infinity() + PETSC_i * mfem::infinity());
      }
      else if (PetscRealPart(sigma) < 0.0)
      {
        SetRegion(-mfem::infinity() - PETSC_i * mfem::infinity(),
                  PetscRealPart(sigma) / gamma + PETSC_i * mfem::infinity());
      }
    }
    else if (PetscRealPart(sigma) == 0.0)
    {
      if (PetscImaginaryPart(sigma) > 0.0)
      {
        SetRegion(-mfem::infinity() + PETSC_i * PetscImaginaryPart(sigma) / gamma,
                  mfem::infinity() + PETSC_i * mfem::infinity());
      }
      else if (PetscImaginaryPart(sigma) < 0.0)
      {
        SetRegion(-mfem::infinity() - PETSC_i * mfem::infinity(),
                  PetscImaginaryPart(sigma) / gamma + PETSC_i * mfem::infinity());
      }
    }
    else
    {
      MFEM_ABORT("Shift-and-invert with general complex eigenvalue target "
                 "is unsupported!");
    }
  }
}

void SlepcEigenSolver::SetRegion(PetscScalar lower_left, PetscScalar upper_right,
                                 bool complement)
{
  RG rg = GetRG();
  PalacePetscCall(RGSetType(rg, RGINTERVAL));
  PalacePetscCall(RGIntervalSetEndpoints(
      rg, PetscRealPart(lower_left), PetscRealPart(upper_right),
      PetscImaginaryPart(lower_left), PetscImaginaryPart(upper_right)));
  if (complement)
  {
    PalacePetscCall(RGSetComplement(rg, PETSC_TRUE));
  }
}

PetscScalar SlepcEigenSolver::GetBackTransform(PetscScalar l) const
{
  return gamma * l;
}

PetscReal SlepcEigenSolver::GetError(int i, EigenvalueSolver::ErrorType type) const
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

// EPS specific methods

SlepcEPSSolverBase::SlepcEPSSolverBase(MPI_Comm comm, int print, const std::string &prefix)
  : SlepcEigenSolver(print)
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
    case ProblemType::HYPERBOLIC:
    case ProblemType::GYROSCOPIC:
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

  PetscInt n;
  PalacePetscCall(VecGetLocalSize(v0, &n));
  MFEM_VERIFY(v.Size() == 2 * n,
              "Invalid size mismatch for provided initial space vector!");

  PetscScalar *pv0;
  PalacePetscCall(VecGetArrayWrite(v0, &pv0));
  v.Get(pv0, n);
  PalacePetscCall(VecRestoreArrayWrite(v0, &pv0));

  Vec is[1] = {v0};
  PalacePetscCall(EPSSetInitialSpace(eps, 1, is));
}

void SlepcEPSSolverBase::Customize()
{
  SlepcEigenSolver::Customize();
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
               opInv->NumTotalMult(), opInv->NumTotalMultIter());
  }

  // Compute and store the eigenpair residuals.
  res = std::make_unique<PetscReal>(num_conv);
  for (int i = 0; i < num_conv; i++)
  {
    res.get()[i] = GetResidualNorm(i);
  }
  return (int)num_conv;
}

PetscScalar SlepcEPSSolverBase::GetEigenvalue(int i) const
{
  PetscScalar l;
  PalacePetscCall(EPSGetEigenvalue(eps, i, &l, nullptr));
  return GetBackTransform(l);
}

void SlepcEPSSolverBase::GetEigenvector(int i, ComplexVector &x) const
{
  MFEM_VERIFY(
      v0,
      "Must call SetOperators before using GetEigenvector for SLEPc eigenvalue solver!");
  PalacePetscCall(EPSGetEigenvector(eps, i, v0, nullptr));

  PetscInt n;
  PalacePetscCall(VecGetLocalSize(v0, &n));
  MFEM_VERIFY(x.Size() == 2 * n, "Invalid size mismatch for provided eigenvector!");

  const PetscScalar *pv0;
  PalacePetscCall(VecGetArrayRead(v0, &pv0));
  x.Set(pv0, n);
  PalacePetscCall(VecRestoreArrayRead(v0, &pv0));
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
  bool first = (opK == nullptr);
  opK = &K;
  opM = &M;

  PetscInt n = opK->Height() / 2;
  PalacePetscCall(
      MatCreateShell(GetComm(), n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &A0));
  PalacePetscCall(
      MatCreateShell(GetComm(), n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &A1));
  auto __mat_apply_shell_A0 = [this](Mat, Vec x_, Vec y_) -> PetscErrorCode
  {
    PetscFunctionBeginUser;
    PetscInt n;
    PetscCall(VecGetLocalSize(x_, &n));

    const PetscScalar *px_;
    PetscCall(VecGetArrayRead(x_, &px_));
    x.Set(px_, n);
    PetscCall(VecRestoreArrayRead(x_, &px_));

    opK->Mult(x, y);

    PetscScalar *py_;
    PetscCall(VecGetArrayWrite(y_, &py_));
    y.Get(py_, n);
    PetscCall(VecRestoreArrayWrite(y_, &py_));

    PetscFunctionReturn(0);
  };
  auto __mat_apply_shell_A1 = [this](Mat, Vec x_, Vec y_) -> PetscErrorCode
  {
    PetscFunctionBeginUser;
    PetscInt n;
    PetscCall(VecGetLocalSize(x_, &n));

    const PetscScalar *px_;
    PetscCall(VecGetArrayRead(x_, &px_));
    x.Set(px_, n);
    PetscCall(VecRestoreArrayRead(x_, &px_));

    opM->Mult(x, y);

    PetscScalar *py_;
    PetscCall(VecGetArrayWrite(y_, &py_));
    y.Get(py_, n);
    PetscCall(VecRestoreArrayWrite(y_, &py_));

    PetscFunctionReturn(0);
  };
  PalacePetscCall(
      MatShellSetOperation(A0, MATOP_MULT, (void (*)()) & __mat_apply_shell_A0));
  PalacePetscCall(
      MatShellSetOperation(A1, MATOP_MULT, (void (*)()) & __mat_apply_shell_A1));
  PalacePetscCall(EPSSetOperators(eps, A0, A1));

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
  x.SetSize(opK->Height());
  y.SetSize(opK->Height());

  // Configure linear solver for generalized problem or spectral transformation. This also
  // allows use of the divergence-free projector as a linear solve side-effect.
  if (first)
  {
    ConfigurePCShell();
  }
}

void SlepcEPSSolver::SetBMat(const Operator &B)
{
  SlepcEigenSolver::SetBMat(B);

  PetscInt n = B.Height();
  PalacePetscCall(
      MatCreateShell(GetComm(), n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &B0));
  auto __mat_apply_shell = [this](Mat, Vec x_, Vec y_) -> PetscErrorCode
  {
    PetscFunctionBeginUser;
    PetscInt n;
    PetscCall(VecGetLocalSize(x_, &n));

    const PetscScalar *px_;
    PetscCall(VecGetArrayRead(x_, &px_));
    x.Set(px_, n);
    PetscCall(VecRestoreArrayRead(x_, &px_));

    opB->Mult(x.Real(), y.Real());
    opB->Mult(x.Imag(), y.Imag());
    y *= delta * gamma;

    PetscScalar *py_;
    PetscCall(VecGetArrayWrite(y_, &py_));
    y.Get(py_, n);
    PetscCall(VecRestoreArrayWrite(y_, &py_));

    PetscFunctionReturn(0);
  };
  PalacePetscCall(MatShellSetOperation(B0, MATOP_MULT, (void (*)()) & __mat_apply_shell));

  BV bv = GetBV();
  PalacePetscCall(BVSetMatrix(bv, B0, PETSC_FALSE));
}

void SlepcEPSSolver::ConfigurePCShell()
{
  auto __pc_apply = [this](PC, Vec x_, Vec y_) -> PetscErrorCode
  {
    // Solve the linear system associated with the generalized eigenvalue problem: y =
    // M⁻¹ x, or shift-and-invert spectral transformation: y = (K - σ M)⁻¹ x . Enforces the
    // divergence-free constraint using the supplied projector.
    PetscFunctionBeginUser;
    PetscInt n;
    PetscCall(VecGetLocalSize(x_, &n));

    const PetscScalar *px_;
    PetscCall(VecGetArrayRead(x_, &px_));
    x.Set(px_, n);
    PetscCall(VecRestoreArrayRead(x_, &px_));

    opInv->Mult(x, y);
    if (!sinvert)
    {
      y *= 1.0 / (delta * gamma);
    }
    else
    {
      y *= 1.0 / delta;
    }
    if (opProj)
    {
      // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(GetComm(), y));
      opProj->Mult(y);
      // Mpi::Print(" After projection: {:e}\n", linalg::Norml2(GetComm(), y));
    }

    PetscScalar *py_;
    PetscCall(VecGetArrayWrite(y_, &py_));
    y.Get(py_, n);
    PetscCall(VecRestoreArrayWrite(y_, &py_));

    PetscFunctionReturn(0);
  };

  KSP ksp;
  PC pc;
  ST st = GetST();
  PalacePetscCall(STGetKSP(st, &ksp));
  PalacePetscCall(KSPGetPC(ksp, &pc));
  PalacePetscCall(PCSetType(pc, PCSHELL));
  PalacePetscCall(PCShellSetApply(pc, (PetscErrorCode(*)(PC, Vec, Vec)) & __pc_apply));
}

PetscReal SlepcEPSSolver::GetResidualNorm(int i) const
{
  // Compute the i-th eigenpair residual: || (K - λ M) x ||₂ for eigenvalue λ.
  PetscScalar l = GetEigenvalue(i);
  GetEigenvector(i, x);
  opK->Mult(x, y);
  opM->AddMult(x, y, -l);
  return linalg::Norml2(GetComm(), y);
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
  bool first = (opK == nullptr);
  opK = &K;
  opC = &C;
  opM = &M;

  PetscInt n = opK->Height() / 2;
  PalacePetscCall(MatCreateShell(GetComm(), 2 * n, 2 * n, PETSC_DECIDE, PETSC_DECIDE,
                                 (void *)this, &A0));
  PalacePetscCall(MatCreateShell(GetComm(), 2 * n, 2 * n, PETSC_DECIDE, PETSC_DECIDE,
                                 (void *)this, &A1));
  auto __mat_apply_shell_A0 = [this](Mat, Vec x_, Vec y_) -> PetscErrorCode
  {
    // Apply the linearized operator L₀ = [  0  I ]
    //                                    [ -K -C ] .
    PetscFunctionBeginUser;
    PetscInt n;
    PetscCall(VecGetLocalSize(x_, &n));

    const PetscScalar *px_;
    PetscCall(VecGetArrayRead(x_, &px_));
    x1.Set(px_, n / 2);
    x2.Set(px_ + n / 2, n / 2);
    PetscCall(VecRestoreArrayRead(x_, &px_));

    y1 = x2;
    opC->Mult(x2, y2);
    y2 *= gamma;
    opK->AddMult(x1, y2, std::complex<double>(1.0, 0.0));
    y2 *= -delta;

    PetscScalar *py_;
    PetscCall(VecGetArrayWrite(y_, &py_));
    y1.Get(py_, n / 2);
    y2.Get(py_ + n / 2, n / 2);
    PetscCall(VecRestoreArrayWrite(y_, &py_));

    PetscFunctionReturn(0);
  };
  auto __mat_apply_shell_A1 = [this](Mat, Vec x_, Vec y_) -> PetscErrorCode
  {
    // Apply the linearized operator L₁ = [ I  0 ]
    //                                    [ 0  M ] .
    PetscFunctionBeginUser;
    PetscInt n;
    PetscCall(VecGetLocalSize(x_, &n));

    const PetscScalar *px_;
    PetscCall(VecGetArrayRead(x_, &px_));
    x1.Set(px_, n / 2);
    x2.Set(px_ + n / 2, n / 2);
    PetscCall(VecRestoreArrayRead(x_, &px_));

    y1 = x1;
    opM->Mult(x2, y2);
    y2 *= delta * gamma * gamma;

    PetscScalar *py_;
    PetscCall(VecGetArrayWrite(y_, &py_));
    y1.Get(py_, n / 2);
    y2.Get(py_ + n / 2, n / 2);
    PetscCall(VecRestoreArrayWrite(y_, &py_));

    PetscFunctionReturn(0);
  };
  PalacePetscCall(
      MatShellSetOperation(A0, MATOP_MULT, (void (*)()) & __mat_apply_shell_A0));
  PalacePetscCall(
      MatShellSetOperation(A1, MATOP_MULT, (void (*)()) & __mat_apply_shell_A1));
  PalacePetscCall(EPSSetOperators(eps, A0, A1));

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

  // Configure linear solver.
  if (first)
  {
    ConfigurePCShell();
  }
}

void SlepcPEPLinearSolver::SetBMat(const Operator &B)
{
  SlepcEigenSolver::SetBMat(B);

  PetscInt n = B.Height();
  PalacePetscCall(MatCreateShell(GetComm(), 2 * n, 2 * n, PETSC_DECIDE, PETSC_DECIDE,
                                 (void *)this, &B0));
  auto __mat_apply_shell = [this](Mat, Vec x_, Vec y_) -> PetscErrorCode
  {
    PetscFunctionBeginUser;
    PetscInt n;
    PetscCall(VecGetLocalSize(x_, &n));

    const PetscScalar *px_;
    PetscCall(VecGetArrayRead(x_, &px_));
    x1.Set(px_, n / 2);
    x2.Set(px_ + n / 2, n / 2);
    PetscCall(VecRestoreArrayRead(x_, &px_));

    opB->Mult(x1.Real(), y1.Real());
    opB->Mult(x1.Imag(), y1.Imag());
    opB->Mult(x2.Real(), y2.Real());
    opB->Mult(x2.Imag(), y2.Imag());
    y1 *= delta * gamma * gamma;
    y2 *= delta * gamma * gamma;

    PetscScalar *py_;
    PetscCall(VecGetArrayWrite(y_, &py_));
    y1.Get(py_, n / 2);
    y2.Get(py_ + n / 2, n / 2);
    PetscCall(VecRestoreArrayWrite(y_, &py_));

    PetscFunctionReturn(0);
  };
  PalacePetscCall(MatShellSetOperation(B0, MATOP_MULT, (void (*)()) & __mat_apply_shell));

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

  PetscInt n;
  PalacePetscCall(VecGetLocalSize(v0, &n));
  MFEM_VERIFY(2 * v.Size() == 2 * n,
              "Invalid size mismatch for provided initial space vector!");

  PetscScalar *pv0;
  PalacePetscCall(VecGetArrayWrite(v0, &pv0));
  v.Get(pv0, n / 2);
  std::fill(pv0 + n / 2, pv0 + n, 0.0);
  PalacePetscCall(VecRestoreArrayWrite(v0, &pv0));

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

  PetscInt n;
  PalacePetscCall(VecGetLocalSize(v0, &n));
  MFEM_VERIFY(2 * x1.Size() == 2 * n, "Invalid size mismatch for provided eigenvector!");

  const PetscScalar *pv0;
  PalacePetscCall(VecGetArrayRead(v0, &pv0));
  x1.Set(pv0, n / 2);
  PalacePetscCall(VecRestoreArrayRead(v0, &pv0));
  if (opB)
  {
    linalg::Normalize(GetComm(), x1, *opB, y1);
  }
  else
  {
    linalg::Normalize(GetComm(), x1);
  }
}

void SlepcPEPLinearSolver::ConfigurePCShell()
{
  auto __pc_apply = [this](PC, Vec x_, Vec y_) -> PetscErrorCode
  {
    // Solve the linear system associated with the generalized eigenvalue problem after
    // linearization: y = L₁⁻¹ x, or with the shift-and-invert spectral transformation:
    // y = (L₀ - σ L₁)⁻¹ x, with:
    //               L₀ = [  0  I ]    L₁ = [ I  0 ]
    //                    [ -K -C ] ,       [ 0  M ] .
    // Enforces the divergence-free constraint using the supplied projector.
    PetscFunctionBeginUser;
    PetscInt n;
    PetscCall(VecGetLocalSize(x_, &n));

    const PetscScalar *px_;
    PetscCall(VecGetArrayRead(x_, &px_));
    x1.Set(px_, n / 2);
    x2.Set(px_ + n / 2, n / 2);
    PetscCall(VecRestoreArrayRead(x_, &px_));

    if (!sinvert)
    {
      y1 = x1;
      if (opProj)
      {
        // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(GetComm(), y1));
        opProj->Mult(y1);
        // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(GetComm(), y1));
      }

      opInv->Mult(x2, y2);
      y2 *= 1.0 / (delta * gamma * gamma);
      if (opProj)
      {
        // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(GetComm(), y2));
        opProj->Mult(y2);
        // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(GetComm(), y2));
      }
    }
    else
    {
      y1.AXPBY(-sigma / (delta * gamma), x2, 0.0);  // Temporarily
      opK->AddMult(x1, y1, std::complex<double>(1.0, 0.0));
      opInv->Mult(y1, y2);
      if (opProj)
      {
        // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(GetComm(), y2));
        opProj->Mult(y2);
        // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(GetComm(), y2));
      }

      y1.AXPBYPCZ(gamma / sigma, y2, -gamma / sigma, x1, 0.0);
      if (opProj)
      {
        // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(GetComm(), y1));
        opProj->Mult(y1);
        // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(GetComm(), y1));
      }
    }

    PetscScalar *py_;
    PetscCall(VecGetArrayWrite(y_, &py_));
    y1.Get(py_, n / 2);
    y2.Get(py_ + n / 2, n / 2);
    PetscCall(VecRestoreArrayWrite(y_, &py_));

    PetscFunctionReturn(0);
  };

  KSP ksp;
  PC pc;
  ST st = GetST();
  PalacePetscCall(STGetKSP(st, &ksp));
  PalacePetscCall(KSPGetPC(ksp, &pc));
  PalacePetscCall(PCSetType(pc, PCSHELL));
  PalacePetscCall(PCShellSetApply(pc, (PetscErrorCode(*)(PC, Vec, Vec)) & __pc_apply));
}

PetscReal SlepcPEPLinearSolver::GetResidualNorm(int i) const
{
  // Compute the i-th eigenpair residual: || P(λ) x ||₂ = || (K + λ C + λ² M) x ||₂ for
  // eigenvalue λ.
  PetscScalar l = GetEigenvalue(i);
  GetEigenvector(i, x1);
  opK->Mult(x1, y1);
  opC->AddMult(x1, y1, l);
  opM->AddMult(x1, y1, l * l);
  return linalg::Norml2(GetComm(), y1);
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
  : SlepcEigenSolver(print)
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

  PetscInt n;
  PalacePetscCall(VecGetLocalSize(v0, &n));
  MFEM_VERIFY(v.Size() == 2 * n,
              "Invalid size mismatch for provided initial space vector!");

  PetscScalar *pv0;
  PalacePetscCall(VecGetArrayWrite(v0, &pv0));
  v.Get(pv0, n);
  PalacePetscCall(VecRestoreArrayWrite(v0, &pv0));

  Vec is[1] = {v0};
  PalacePetscCall(PEPSetInitialSpace(pep, 1, is));
}

void SlepcPEPSolverBase::Customize()
{
  SlepcEigenSolver::Customize();
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
               opInv->NumTotalMult(), opInv->NumTotalMultIter());
  }

  // Compute and store the eigenpair residuals.
  res = std::make_unique<PetscReal>(num_conv);
  for (int i = 0; i < num_conv; i++)
  {
    res.get()[i] = GetResidualNorm(i);
  }
  return (int)num_conv;
}

PetscScalar SlepcPEPSolverBase::GetEigenvalue(int i) const
{
  PetscScalar l;
  PalacePetscCall(PEPGetEigenpair(pep, i, &l, nullptr, nullptr, nullptr));
  return GetBackTransform(l);
}

void SlepcPEPSolverBase::GetEigenvector(int i, ComplexVector &x) const
{
  MFEM_VERIFY(
      v0,
      "Must call SetOperators before using GetEigenvector for SLEPc eigenvalue solver!");
  PalacePetscCall(PEPGetEigenpair(pep, i, nullptr, nullptr, v0, nullptr));

  PetscInt n;
  PalacePetscCall(VecGetLocalSize(v0, &n));
  MFEM_VERIFY(x.Size() == 2 * n, "Invalid size mismatch for provided eigenvector!");

  const PetscScalar *pv0;
  PalacePetscCall(VecGetArrayRead(v0, &pv0));
  x.Set(pv0, n);
  PalacePetscCall(VecRestoreArrayRead(v0, &pv0));
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
  bool first = (opK == nullptr);
  opK = &K;
  opC = &C;
  opM = &M;

  PetscInt n = opK->Height() / 2;
  PalacePetscCall(
      MatCreateShell(GetComm(), n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &A0));
  PalacePetscCall(
      MatCreateShell(GetComm(), n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &A1));
  PalacePetscCall(
      MatCreateShell(GetComm(), n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &A2));
  auto __mat_apply_shell_A0 = [this](Mat, Vec x_, Vec y_) -> PetscErrorCode
  {
    PetscFunctionBeginUser;
    PetscInt n;
    PetscCall(VecGetLocalSize(x_, &n));

    const PetscScalar *px_;
    PetscCall(VecGetArrayRead(x_, &px_));
    x.Set(px_, n);
    PetscCall(VecRestoreArrayRead(x_, &px_));

    opK->Mult(x, y);

    PetscScalar *py_;
    PetscCall(VecGetArrayWrite(y_, &py_));
    y.Get(py_, n);
    PetscCall(VecRestoreArrayWrite(y_, &py_));

    PetscFunctionReturn(0);
  };
  auto __mat_apply_shell_A1 = [this](Mat, Vec x_, Vec y_) -> PetscErrorCode
  {
    PetscFunctionBeginUser;
    PetscInt n;
    PetscCall(VecGetLocalSize(x_, &n));

    const PetscScalar *px_;
    PetscCall(VecGetArrayRead(x_, &px_));
    x.Set(px_, n);
    PetscCall(VecRestoreArrayRead(x_, &px_));

    opC->Mult(x, y);

    PetscScalar *py_;
    PetscCall(VecGetArrayWrite(y_, &py_));
    y.Get(py_, n);
    PetscCall(VecRestoreArrayWrite(y_, &py_));

    PetscFunctionReturn(0);
  };
  auto __mat_apply_shell_A2 = [this](Mat, Vec x_, Vec y_) -> PetscErrorCode
  {
    PetscFunctionBeginUser;
    PetscInt n;
    PetscCall(VecGetLocalSize(x_, &n));

    const PetscScalar *px_;
    PetscCall(VecGetArrayRead(x_, &px_));
    x.Set(px_, n);
    PetscCall(VecRestoreArrayRead(x_, &px_));

    opM->Mult(x, y);

    PetscScalar *py_;
    PetscCall(VecGetArrayWrite(y_, &py_));
    y.Get(py_, n);
    PetscCall(VecRestoreArrayWrite(y_, &py_));

    PetscFunctionReturn(0);
  };
  PalacePetscCall(
      MatShellSetOperation(A0, MATOP_MULT, (void (*)()) & __mat_apply_shell_A0));
  PalacePetscCall(
      MatShellSetOperation(A1, MATOP_MULT, (void (*)()) & __mat_apply_shell_A1));
  PalacePetscCall(
      MatShellSetOperation(A2, MATOP_MULT, (void (*)()) & __mat_apply_shell_A2));
  Mat A[3] = {A0, A1, A2};
  PalacePetscCall(PEPSetOperators(pep, 3, A));

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
  x.SetSize(opK->Height());
  y.SetSize(opK->Height());

  // Configure linear solver.
  if (first)
  {
    ConfigurePCShell();
  }
}

void SlepcPEPSolver::SetBMat(const Operator &B)
{
  SlepcEigenSolver::SetBMat(B);

  PetscInt n = B.Height();
  PalacePetscCall(
      MatCreateShell(GetComm(), n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &B0));
  auto __mat_apply_shell = [this](Mat, Vec x_, Vec y_) -> PetscErrorCode
  {
    PetscFunctionBeginUser;
    PetscInt n;
    PetscCall(VecGetLocalSize(x_, &n));

    const PetscScalar *px_;
    PetscCall(VecGetArrayRead(x_, &px_));
    x.Set(px_, n);
    PetscCall(VecRestoreArrayRead(x_, &px_));

    opB->Mult(x.Real(), y.Real());
    opB->Mult(x.Imag(), y.Imag());
    y *= delta * gamma;

    PetscScalar *py_;
    PetscCall(VecGetArrayWrite(y_, &py_));
    y.Get(py_, n);
    PetscCall(VecRestoreArrayWrite(y_, &py_));

    PetscFunctionReturn(0);
  };
  PalacePetscCall(MatShellSetOperation(B0, MATOP_MULT, (void (*)()) & __mat_apply_shell));

  BV bv = GetBV();
  PalacePetscCall(BVSetMatrix(bv, B0, PETSC_FALSE));
}

void SlepcPEPSolver::ConfigurePCShell()
{
  auto __pc_apply = [this](PC, Vec x_, Vec y_) -> PetscErrorCode
  {
    // Solve the linear system associated with the generalized eigenvalue problem: y =
    // M⁻¹ x, or shift-and-invert spectral transformation: y = P(σ)⁻¹ x . Enforces the
    // divergence-free constraint using the supplied projector.
    PetscFunctionBeginUser;
    PetscInt n;
    PetscCall(VecGetLocalSize(x_, &n));

    const PetscScalar *px_;
    PetscCall(VecGetArrayRead(x_, &px_));
    x.Set(px_, n);
    PetscCall(VecRestoreArrayRead(x_, &px_));

    opInv->Mult(x, y);
    if (!sinvert)
    {
      y *= 1.0 / (delta * gamma * gamma);
    }
    else
    {
      y *= 1.0 / delta;
    }
    if (opProj)
    {
      // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(GetComm(), y));
      opProj->Mult(y);
      // Mpi::Print(" After projection: {:e}\n", linalg::Norml2(GetComm(), y));
    }

    PetscScalar *py_;
    PetscCall(VecGetArrayWrite(y_, &py_));
    y.Get(py_, n);
    PetscCall(VecRestoreArrayWrite(y_, &py_));

    PetscFunctionReturn(0);
  };

  KSP ksp;
  PC pc;
  ST st = GetST();
  PalacePetscCall(STGetKSP(st, &ksp));
  PalacePetscCall(KSPGetPC(ksp, &pc));
  PalacePetscCall(PCSetType(pc, PCSHELL));
  PalacePetscCall(PCShellSetApply(pc, (PetscErrorCode(*)(PC, Vec, Vec)) & __pc_apply));
}

PetscReal SlepcPEPSolver::GetResidualNorm(int i) const
{
  // Compute the i-th eigenpair residual: || P(λ) x ||₂ = || (K + λ C + λ² M) x ||₂ for
  // eigenvalue λ.
  PetscScalar l = GetEigenvalue(i);
  GetEigenvector(i, x);
  opK->Mult(x, y);
  opC->AddMult(x, y, l);
  opM->AddMult(x, y, l * l);
  return linalg::Norml2(GetComm(), y);
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

#endif
