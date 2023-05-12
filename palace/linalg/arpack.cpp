// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "arpack.hpp"

#if 0  // XX TODO DISABLE ARPACK FOR NOW

#if defined(PALACE_WITH_ARPACK)

#if defined(__GNUC__) && defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wc99-extensions"
#endif

#include <algorithm>
#include <string>
#include <mfem.hpp>
// clang-format off
#include <parpack.hpp>  // ARPACK headers
#include <debug_c.hpp>
#include <stat_c.hpp>
// clang-format on
#include "linalg/divfree.hpp"
#include "linalg/ksp.hpp"
#include "linalg/petsc.hpp"
#include "utils/communication.hpp"

namespace palace::arpack
{

// Base class methods

ArpackEigenSolver::ArpackEigenSolver(int print_lvl)
{
  // Initialization.
  print = print_lvl;
  info = 0;
  nev = ncv = 0;
  rtol = 0.0;
  max_it = 0;
  which_option = ::arpack::which::largest_magnitude;
  sinvert = false;
  sigma = 0.0;
  gamma = delta = 1.0;

  eig = nullptr;
  perm = nullptr;
  V = nullptr;
  res = nullptr;
  r = nullptr;
  opInv = nullptr;
  opProj = nullptr;
  opB = nullptr;

  // Configure debugging output.
  a_int logfill = 6, ndigit = -6, mgetv0 = 0;
  a_int _aupd = (print_lvl > 2) ? 1 : 0,
        _aup2 = (print_lvl > 2) ? 2 : ((print_lvl > 0) ? 1 : 0), _aitr = 0, _eigh = 0,
        _gets = 0, _apps = 0, _eupd = 0;
  debug_c(logfill, ndigit, mgetv0, _aupd, _aup2, _aitr, _eigh, _gets, _apps, _eupd, _aupd,
          _aup2, _aitr, _eigh, _gets, _apps, _eupd, _aupd, _aup2, _aitr, _eigh, _gets,
          _apps, _eupd);
  cstatn_c();
}

ArpackEigenSolver::~ArpackEigenSolver()
{
  delete[] eig;
  delete[] perm;
  delete[] res;
  delete V;
  delete r;
}

void ArpackEigenSolver::SetOperators(const petsc::PetscParMatrix &K,
                                     const petsc::PetscParMatrix &M,
                                     EigenSolverBase::ScaleType type)
{
  MFEM_ABORT("SetOperators not defined for base class ArpackEigenSolver!");
}

void ArpackEigenSolver::SetOperators(const petsc::PetscParMatrix &K,
                                     const petsc::PetscParMatrix &C,
                                     const petsc::PetscParMatrix &M,
                                     EigenSolverBase::ScaleType type)
{
  MFEM_ABORT("SetOperators not defined for base class ArpackEigenSolver!");
}

void ArpackEigenSolver::SetLinearSolver(const KspSolver &ksp)
{
  opInv = &ksp;
}

void ArpackEigenSolver::SetProjector(const DivFreeSolver &divfree)
{
  opProj = &divfree;
}

void ArpackEigenSolver::SetBMat(const petsc::PetscParMatrix &B)
{
  MFEM_VERIFY(!opB || opB->GetNumRows() == B.GetNumRows(),
              "Invalid modification of eigenvalue problem size!");
  opB = &B;
}

void ArpackEigenSolver::SetNumModes(int numeig, int numvec)
{
  if (nev > 0 && numeig != nev)
  {
    delete[] eig;
    delete[] perm;
    delete[] res;
    eig = nullptr;
    perm = nullptr;
    res = nullptr;
  }
  if (ncv > 0 && numvec != ncv)
  {
    delete V;
    V = nullptr;
  }
  nev = numeig;
  ncv = (numvec > 0) ? numvec : std::max(20, 2 * nev + 1);  // Default from SLEPc
}

void ArpackEigenSolver::SetTol(double tol)
{
  rtol = tol;
}

void ArpackEigenSolver::SetMaxIter(int maxits)
{
  max_it = maxits;
}

void ArpackEigenSolver::SetWhichEigenpairs(EigenSolverBase::WhichType type)
{
  switch (type)
  {
    case WhichType::LARGEST_MAGNITUDE:
    case WhichType::TARGET_MAGNITUDE:
      which_option = ::arpack::which::largest_magnitude;
      break;
    case WhichType::SMALLEST_MAGNITUDE:
      which_option = ::arpack::which::smallest_magnitude;
      break;
    case WhichType::LARGEST_REAL:
      which_option = ::arpack::which::largest_real;
      break;
    case WhichType::SMALLEST_REAL:
      which_option = ::arpack::which::smallest_real;
      break;
    case WhichType::LARGEST_IMAGINARY:
      which_option = ::arpack::which::largest_imaginary;
      break;
    case WhichType::SMALLEST_IMAGINARY:
      which_option = ::arpack::which::smallest_imaginary;
      break;
    case WhichType::TARGET_REAL:
    case WhichType::TARGET_IMAGINARY:
      MFEM_ABORT("ARPACK eigenvalue solver does not implement TARGET_REAL or "
                 "TARGET_IMAGINARY for SetWhichEigenpairs!");
      break;
  }
}

void ArpackEigenSolver::SetShiftInvert(double tr, double ti, bool precond)
{
  MFEM_VERIFY(!precond, "ARPACK eigenvalue solver does not support preconditioned "
                        "spectral transformation option!");
  sigma = tr + PETSC_i * ti;
  sinvert = true;
}

void ArpackEigenSolver::SetInitialSpace(const petsc::PetscParVector &v)
{
  if (!r)
  {
    r = new petsc::PetscParVector(v);
  }
  else
  {
    MFEM_VERIFY(v.GetSize() == r->GetSize(),
                "Invalid modification of eigenvalue problem size!");
    r->Copy(v);
  }
  info = 1;
}

int ArpackEigenSolver::SolveInternal(petsc::PetscParVector &r_, petsc::PetscDenseMatrix &V_,
                                     PetscScalar *eig_, int *perm_)
{
  MPI_Comm comm;
  MPI_Fint fcomm;
  a_int ido, info_ = (a_int)info;
  a_int iparam[11] = {0}, ipntr[14] = {0};
  a_int n, nev_, ncv_;
  ::arpack::bmat bmat_option =
      (opB) ? ::arpack::bmat::generalized : ::arpack::bmat::identity;
  PetscScalar *workd, *workl;
  double *rwork;
  a_int lworkl;

  comm = r_.GetComm();
  fcomm = MPI_Comm_c2f(comm);
  iparam[0] = 1;                // Exact shifts
  iparam[2] = (a_int)max_it;    // Maximum number of Arnoldi iterations
  iparam[3] = 1;                // Block size
  iparam[4] = 0;                // Number of converged Ritz values
  iparam[6] = sinvert ? 3 : 1;  // Problem mode

  // Set problem sizes. The cast to int should always be safe because this is a local size.
  n = (a_int)r_.GetSize();
  nev_ = (a_int)nev;
  ncv_ = (a_int)ncv;

  // Allocate work arrays.
  lworkl = 3 * ncv_ * ncv_ + 5 * ncv_;
  workd = new PetscScalar[3 * n];
  workl = new PetscScalar[lworkl];
  rwork = new double[ncv_];

  PetscScalar *pr_ = r_.GetArray();
  PetscScalar *pV_ = V_.GetArray();
  petsc::PetscParVector x(comm, n, PETSC_DECIDE, nullptr);
  petsc::PetscParVector y(comm, n, PETSC_DECIDE, nullptr);

  // Begin RCI loop.
  ido = 0;
  while (true)
  {
    // Call complex problem driver.
    naupd(fcomm, ido, bmat_option, n, which_option, nev_, rtol, pr_, ncv_, pV_, n, iparam,
          ipntr, workd, workl, lworkl, rwork, info_);
    CheckInfoAUPD(info_);

    // We never use pre-computed B * x in workd[ipntr[2]-1].
    x.PlaceArray(&workd[ipntr[0] - 1]);
    y.PlaceArray(&workd[ipntr[1] - 1]);
    if (ido == 1 || ido == -1)
    {
      ApplyOp(x, y);
    }
    else if (ido == 2)
    {
      ApplyOpB(x, y);
    }
    else if (ido == 99)
    {
      break;
    }
    else
    {
      MFEM_ABORT("Internal error in ARPACK RCI interface!");
    }
    x.ResetArray();
    y.ResetArray();
  }

  // Print some log information.
  int niter = (int)iparam[2];
  int nconv = (int)iparam[4];
  if (print > 0)
  {
    Mpi::Print(comm,
               "\n ARPACK {} eigensolve {} ({:d} eigenpairs); iterations {:d}\n"
               " Total number of linear systems solved: {:d}\n"
               " Total number of linear solver iterations: {:d}\n",
               GetName(), (nconv >= nev_) ? "converged" : "finished", nconv, niter,
               opInv->GetTotalNumMult(), opInv->GetTotalNumIter());
  }
  if (nconv < nev_)
  {
    Mpi::Warning(
        comm, "ARPACK eigenvalue solver found only {:d} of requested {:d} eigenvalues!\n",
        nconv, nev_);
  }

  // Postprocess eigenvalues and eigenvectors.
  a_int rvec = 1;
  ::arpack::howmny howmny_option = ::arpack::howmny::ritz_vectors;
  a_int *select;
  PetscScalar *workev;

  // Allocate eigenvalue storage and work arrays.
  select = new a_int[ncv_];
  workev = new PetscScalar[2 * ncv_];

  // Call complex problem driver.
  PetscScalar sigma_ = sigma / gamma;
  neupd(fcomm, rvec, howmny_option, select, eig_, pV_, n, sigma_, workev, bmat_option, n,
        which_option, nev_, rtol, pr_, ncv_, pV_, n, iparam, ipntr, workd, workl, lworkl,
        rwork, info_);
  CheckInfoEUPD(info_);

  // Unscale and properly sort the eigenvalues.
  auto CompareReal = [&eig_](const int &l, const int &r)
  { return PetscRealPart(eig_[l]) < PetscRealPart(eig_[r]); };
  auto CompareImag = [&eig_](const int &l, const int &r)
  { return PetscImaginaryPart(eig_[l]) < PetscImaginaryPart(eig_[r]); };
  auto CompareAbs = [&eig_](const int &l, const int &r)
  { return PetscAbsScalar(eig_[l]) < PetscAbsScalar(eig_[r]); };
  for (int i = 0; i < nev_; i++)
  {
    eig_[i] = eig_[i] * gamma;
    perm_[i] = i;
  }
  if (which_option == ::arpack::which::largest_real ||
      which_option == ::arpack::which::smallest_real)
  {
    std::sort(perm_, perm_ + nev_, CompareReal);
  }
  else if (which_option == ::arpack::which::largest_imaginary ||
           which_option == ::arpack::which::smallest_imaginary)
  {
    std::sort(perm_, perm_ + nev_, CompareImag);
  }
  else
  {
    std::sort(perm_, perm_ + nev_, CompareAbs);
  }

  // Cleanup.
  r_.RestoreArray(pr_);
  V_.RestoreArray(pV_);
  delete[] select;
  delete[] workev;
  delete[] workd;
  delete[] workl;
  delete[] rwork;

  return nconv;
}

void ArpackEigenSolver::CheckParameters() const
{
  MFEM_VERIFY(nev > 0, "Number of requested modes is not positive!");
  MFEM_VERIFY(rtol > 0.0, "Eigensolver tolerance is not positive!");
  MFEM_VERIFY(opInv, "No linear solver provided for operator!");
}

void ArpackEigenSolver::CheckInfoAUPD(int info) const
{
  if (info != 0)
  {
    std::string msg = "ARPACK pznaupd error: ";
    switch (info)
    {
      case 1:
        msg += "Maximum number of iterations taken, all possible eigenvalues "
               "have been found!";
        break;
      case 2:
        msg += "No longer an informational error (deprecated starting with "
               "release 2 of ARPACK)!";
        break;
      case 3:
        msg += "No shifts could be applied during a cycle of the Implicitly "
               "restarted Arnoldi iteration!";
        break;
      case -1:
        msg += "N must be positive!";
        break;
      case -2:
        msg += "NEV must be positive!";
        break;
      case -3:
        msg += "NCV-NEV >= 2 and less than or equal to N!";
        break;
      case -4:
        msg += "The maximum number of Arnoldi update iterations allowed must "
               "be greater than zero!";
        break;
      case -5:
        msg += "WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'";
        break;
      case -6:
        msg += "BMAT must be one of 'I' or 'G'!";
        break;
      case -7:
        msg += "Length of private work array WORKL is not sufficient!";
        break;
      case -8:
        msg += "Error return from LAPACK eigenvalue calculation!";
        break;
      case -9:
        msg += "Starting vector is zero!";
        break;
      case -10:
        msg += "IPARAM(7) must be 1, 2, or 3!";
        break;
      case -11:
        msg += "IPARAM(7) = 1 and BMAT = 'G' are incompatible!";
        break;
      case -12:
        msg += "IPARAM(1) must be equal to 0 or 1!";
        break;
      case -9999:
        msg += "Could not build an Arnoldi factorization!";
        break;
      default:
        msg += "Unknown ARPACK error message!";
        break;
    }
    MFEM_ABORT(msg.c_str());
  }
}

void ArpackEigenSolver::CheckInfoEUPD(int info) const
{
  if (info != 0)
  {
    std::string msg = "ARPACK pzneupd error: ";
    switch (info)
    {
      case 1:
        msg += "The Schur form computed by LAPACK routine csheqr could not "
               "be reordered by LAPACK routine ztrsen!";
        break;
      case -1:
        msg += "N must be positive!";
        break;
      case -2:
        msg += "NEV must be positive!";
        break;
      case -3:
        msg += "NCV-NEV >= 2 and less than or equal to N!";
        break;
      case -4:
        msg += "The maximum number of Arnoldi update iterations allowed must "
               "be greater than zero!";
        break;
      case -5:
        msg += "WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'";
        break;
      case -6:
        msg += "BMAT must be one of 'I' or 'G'!";
        break;
      case -7:
        msg += "Length of private work array WORKL is not sufficient!";
        break;
      case -8:
        msg += "Error return from LAPACK eigenvalue calculation!";
        break;
      case -9:
        msg += "Error return from calculation of eigenvectors!";
        break;
      case -10:
        msg += "IPARAM(7) must be 1, 2, or 3!";
        break;
      case -11:
        msg += "IPARAM(7) = 1 and BMAT = 'G' are incompatible!";
        break;
      case -12:
        msg += "HOWMNY = 'S' not yet implemented!";
        break;
      case -13:
        msg += "HOWMNY must be one of 'A' or 'P' if RVEC = true!";
        break;
      case -14:
        msg += "PZNAUPD did not find any eigenvalues to sufficient accuracy!";
        break;
      case -15:
        msg += "ZNEUPD got a different count of the number of converged Ritz "
               "values than ZNAUPD got!";
        break;
      default:
        msg += "Unknown ARPACK error message!";
        break;
    }
    MFEM_ABORT(msg.c_str());
  }
}

void ArpackEigenSolver::GetEigenvalue(int i, double &eigr, double &eigi) const
{
  MFEM_VERIFY(eig && i >= 0 && i < nev,
              "Out of range eigenpair requested (i = " << i << ", nev = " << nev << ")!");
  const int &j = perm[i];
  eigr = PetscRealPart(eig[j]);
  eigi = PetscImaginaryPart(eig[j]);
}

void ArpackEigenSolver::GetEigenvector(int i, petsc::PetscParVector &x) const
{
  MFEM_VERIFY(eig && i >= 0 && i < nev,
              "Out of range eigenpair requested (i = " << i << ", nev = " << nev << ")!");
  const int &j = perm[i];
  const petsc::PetscParVector v = V->GetColumnRead(j);
  x.Copy(v);
  V->RestoreColumnRead(j, v);
}

void ArpackEigenSolver::GetError(int i, EigenSolverBase::ErrorType type, double &err) const
{
  MFEM_VERIFY(eig && i >= 0 && i < nev,
              "Out of range eigenpair requested (i = " << i << ", nev = " << nev << ")!");
  const int &j = perm[i];
  if (res[j] <= 0.0)
  {
    const petsc::PetscParVector v = V->GetColumnRead(j);
    GetResidual(eig[j], v, *r);
    res[j] = r->Norml2() / v.Norml2();
    V->RestoreColumnRead(j, v);
  }
  switch (type)
  {
    case ErrorType::ABSOLUTE:
      err = res[j];
      break;
    case ErrorType::RELATIVE:
      err = res[j] / PetscAbsScalar(eig[j]);
      break;
    case ErrorType::BACKWARD:
      err = res[j] / GetBackwardScaling(eig[j]);
      break;
  }
}

// EPS specific methods

ArpackEPSSolver::ArpackEPSSolver(int print_lvl) : ArpackEigenSolver(print_lvl)
{
  opK = opM = nullptr;
  normK = normM = 0.0;
  z = nullptr;
}

ArpackEPSSolver::~ArpackEPSSolver()
{
  delete z;
}

void ArpackEPSSolver::SetOperators(const petsc::PetscParMatrix &K,
                                   const petsc::PetscParMatrix &M,
                                   EigenSolverBase::ScaleType type)
{
  MFEM_VERIFY(!opK || opK->GetNumRows() == K.GetNumRows(),
              "Invalid modification of eigenvalue problem size!");
  bool first = (opK == nullptr);
  opK = &K;
  opM = &M;
  if (first && type != ScaleType::NONE)
  {
    normK = opK->Norm2();
    normM = opM->Norm2();
    MFEM_VERIFY(normK >= 0.0 && normM >= 0.0, "Invalid matrix norms for EPS scaling!");
    if (normK > 0 && normM > 0.0)
    {
      gamma = normK / normM;  // Store γ² for linear problem
      delta = 2.0 / normK;
    }
  }

  // Set up workspace.
  if (!z)
  {
    z = new petsc::PetscParVector(K);
  }
}

int ArpackEPSSolver::Solve()
{
  // Check input parameters.
  CheckParameters();
  MFEM_VERIFY(opK && opM, "Operators are not set for ArpackEPSSolver!");

  // Set some defaults (default maximum iterations from SLEPc).
  PetscInt n = opK->GetNumRows(), N = opK->GetGlobalNumRows();
  if (ncv > N)
  {
    ncv = (int)N;
  }
  if (max_it <= 0)
  {
    max_it = std::max(300, (int)(2 * N / ncv));
  }

  // Initialize if user did not provide an initial space.
  if (!r)
  {
    info = 0;
    r = new petsc::PetscParVector(*opK);
  }
  if (!info)
  {
    r->SetZero();
  }

  // Allocate Arnoldi basis for the problem.
  if (!V)
  {
    V = new petsc::PetscDenseMatrix(opK->GetComm(), n, PETSC_DECIDE, PETSC_DECIDE, ncv,
                                    nullptr);
  }

  // Cache residual norms when calculated later on.
  if (!eig)
  {
    eig = new PetscScalar[nev + 1];
    perm = new int[nev + 1];
    res = new double[nev + 1];
  }
  for (int i = 0; i < nev + 1; i++)
  {
    res[i] = -1.0;
  }

  // Solve the generalized eigenvalue problem.
  int nconv = SolveInternal(*r, *V, eig, perm);

  // Reset for next solve.
  info = 0;
  return nconv;
}

void ArpackEPSSolver::ApplyOp(const petsc::PetscParVector &x,
                              petsc::PetscParVector &y) const
{
  // Case 1: No spectral transformation (opInv = M⁻¹)
  //               y = M⁻¹ K x .
  // Case 2: Shift-and-invert spectral transformation (opInv = (K - σ M)⁻¹)
  //               y = (K - σ M)⁻¹ M x .
  if (!sinvert)
  {
    opK->Mult(x, *z);
    opInv->Mult(*z, y);
    y.Scale(1.0 / gamma);
  }
  else
  {
    opM->Mult(x, *z);
    opInv->Mult(*z, y);
    y.Scale(gamma);
  }

  // Debug
  // Mpi::Print(" Before projection: {:e}\n", y.Norml2());

  if (opProj)
  {
    opProj->Mult(y);
  }

  // Debug
  // Mpi::Print(" After projection: {:e}\n", y.Norml2());
}

void ArpackEPSSolver::ApplyOpB(const petsc::PetscParVector &x,
                               petsc::PetscParVector &y) const
{
  MFEM_VERIFY(opB, "No B operator for weighted inner product in ARPACK solve!");
  opB->Mult(x, y);
  y.Scale(delta * gamma);
}

void ArpackEPSSolver::GetResidual(PetscScalar l, const petsc::PetscParVector &x,
                                  petsc::PetscParVector &r) const
{
  // r = (K - λ M) x for eigenvalue λ.
  opM->Mult(x, r);
  r.Scale(-l);
  opK->MultAdd(x, r);
}

double ArpackEPSSolver::GetBackwardScaling(PetscScalar l) const
{
  // Make sure not to use norms from scaling as this can be confusing if they are different.
  // Note that SLEPc uses ||.||∞, not the 2-norm.
  if (normK <= 0.0)
  {
    normK = opK->Norm2();
  }
  if (normM <= 0.0)
  {
    normM = opM->Norm2();
  }
  return normK + PetscAbsScalar(l) * normM;
}

// PEP specific methods

ArpackPEPSolver::ArpackPEPSolver(int print_lvl) : ArpackEigenSolver(print_lvl)
{
  opK = opC = opM = nullptr;
  normK = normC = normM = 0.0;
  x1 = x2 = y1 = y2 = z = nullptr;
}

ArpackPEPSolver::~ArpackPEPSolver()
{
  delete x1;
  delete x2;
  delete y1;
  delete y2;
  delete z;
}

void ArpackPEPSolver::SetOperators(const petsc::PetscParMatrix &K,
                                   const petsc::PetscParMatrix &C,
                                   const petsc::PetscParMatrix &M,
                                   EigenSolverBase::ScaleType type)
{
  MFEM_VERIFY(!opK || opK->GetNumRows() == K.GetNumRows(),
              "Invalid modification of eigenvalue problem size!");
  bool first = (opK == nullptr);
  opK = &K;
  opC = &C;
  opM = &M;
  if (first && type != ScaleType::NONE)
  {
    normK = opK->Norm2();
    normC = opC->Norm2();
    normM = opM->Norm2();
    MFEM_VERIFY(normK >= 0.0 && normC >= 0.0 && normM >= 0.0,
                "Invalid matrix norms for PEP scaling!");
    if (normK > 0 && normC > 0.0 && normM > 0.0)
    {
      gamma = std::sqrt(normK / normM);
      delta = 2.0 / (normK + gamma * normC);
    }
  }

  // Set up workspace.
  if (!z)
  {
    MPI_Comm comm = K.GetComm();
    PetscInt n = K.GetNumRows();
    delete x1;
    delete x2;
    delete y1;
    delete y2;
    delete z;
    x1 = new petsc::PetscParVector(comm, n, PETSC_DECIDE, nullptr);
    x2 = new petsc::PetscParVector(comm, n, PETSC_DECIDE, nullptr);
    y1 = new petsc::PetscParVector(comm, n, PETSC_DECIDE, nullptr);
    y2 = new petsc::PetscParVector(comm, n, PETSC_DECIDE, nullptr);
    z = new petsc::PetscParVector(K);
  }
}

int ArpackPEPSolver::Solve()
{
  // Check input parameters.
  CheckParameters();
  MFEM_VERIFY(opK && opC && opM, "Operators are not set for ArpackPEPSolver!");

  // Set some defaults (from SLEPc ARPACK interface). The problem size is the size of the
  // 2x2 block linearized problem.
  PetscInt n = opK->GetNumRows(), N = opK->GetGlobalNumRows();
  if (ncv > 2 * N)
  {
    ncv = 2 * (int)N;
  }
  if (max_it <= 0)
  {
    max_it = std::max(300, 4 * (int)(N / ncv));
  }

  // Initialize if user did not provide an initial space.
  if (!r)
  {
    info = 0;
    r = new petsc::PetscParVector(*opK);
  }
  if (!info)
  {
    r->SetZero();
  }
  petsc::PetscParVector *s = new petsc::PetscParVector(opK->GetComm(), 2 * n, PETSC_DECIDE);
  PetscScalar *ps = GetBlocks(*s, *x1, *x2);
  x1->Copy(*r);
  x2->SetZero();  // Second block initialized to zero even with initial guess
  RestoreBlocks(ps, *s, *x1, *x2);

  // Allocate Arnoldi basis for original and linearized problem.
  if (!V)
  {
    V = new petsc::PetscDenseMatrix(opK->GetComm(), n, PETSC_DECIDE, PETSC_DECIDE, ncv,
                                    nullptr);
  }
  petsc::PetscDenseMatrix *W = new petsc::PetscDenseMatrix(
      opK->GetComm(), 2 * n, PETSC_DECIDE, PETSC_DECIDE, ncv, nullptr);

  // Cache residual norms when calculated later on.
  if (!eig)
  {
    eig = new PetscScalar[nev + 1];
    perm = new int[nev + 1];
    res = new double[nev + 1];
  }
  for (int i = 0; i < nev + 1; i++)
  {
    res[i] = -1.0;
  }

  // Solve the linearized eigenvalue problem.
  int nconv = SolveInternal(*s, *W, eig, perm);

  // Eigenvector extraction from the linearized eigenproblem.
  for (int i = 0; i < nev; i++)
  {
    petsc::PetscParVector w = W->GetColumn(i);
    petsc::PetscParVector v = V->GetColumn(i);
    ExtractEigenvector(eig[i], w, v);
    W->RestoreColumn(i, w);
    V->RestoreColumn(i, v);
  }

  // Cleanup auxiliary basis and residual vector.
  delete W;
  delete s;

  // Reset for next solve.
  info = 0;
  return nconv;
}

void ArpackPEPSolver::ApplyOp(const petsc::PetscParVector &x,
                              petsc::PetscParVector &y) const
{
  // Case 1: No spectral transformation (opInv = M⁻¹)
  //               y = L₁⁻¹ L₀ x .
  // Case 2: Shift-and-invert spectral transformation (opInv = P(σ)⁻¹)
  //               y = (L₀ - σ L₁)⁻¹ L₁ x .
  // With:
  //               L₀ = [ -K  0 ]    L₁ = [ C  M ]
  //                    [  0  M ] ,       [ M  0 ] .
  PetscScalar *px = GetBlocks(const_cast<petsc::PetscParVector &>(x), *x1, *x2);
  PetscScalar *py = GetBlocks(y, *y1, *y2);
  if (!sinvert)
  {
    opC->Mult(*x2, *z);
    z->Scale(gamma);
    opK->MultAdd(*x1, *z);
    opInv->Mult(*z, *y2);
    y2->Scale(-1.0 / (gamma * gamma));
    if (opProj)
    {
      opProj->Mult(*y2);
    }
    y1->Copy(*x2);
    if (opProj)
    {
      opProj->Mult(*y1);
    }
  }
  else
  {
    y1->AXPBYPCZ(sigma, *x1, gamma, *x2, 0.0);  // Just temporarily
    opM->Mult(*y1, *z);
    opC->MultAdd(*x1, *z);
    z->Scale(-gamma);
    opInv->Mult(*z, *y1);

    // Debug
    // Mpi::Print(" Before projection: {:e}\n", y1->Norml2());

    if (opProj)
    {
      opProj->Mult(*y1);
    }

    // Debug
    // Mpi::Print(" After projection: {:e}\n", y1->Norml2());

    y2->AXPBYPCZ(sigma / gamma, *y1, 1.0, *x1, 0.0);

    // Debug
    // Mpi::Print(" Before projection: {:e}\n", y2->Norml2());

    if (opProj)
    {
      opProj->Mult(*y2);
    }

    // Debug
    // Mpi::Print(" After projection: {:e}\n", y2->Norml2());
  }
  RestoreBlocks(px, const_cast<petsc::PetscParVector &>(x), *x1, *x2);
  RestoreBlocks(py, y, *y1, *y2);
}

void ArpackPEPSolver::ApplyOpB(const petsc::PetscParVector &x,
                               petsc::PetscParVector &y) const
{
  MFEM_VERIFY(opB, "No B operator for weighted inner product in ARPACK solve!");
  PetscScalar *px = GetBlocks(const_cast<petsc::PetscParVector &>(x), *x1, *x2);
  PetscScalar *py = GetBlocks(y, *y1, *y2);
  opB->Mult(*x1, *y1);
  opB->Mult(*x2, *y2);
  y1->Scale(delta * gamma * gamma);
  y2->Scale(delta * gamma * gamma);
  RestoreBlocks(px, const_cast<petsc::PetscParVector &>(x), *x1, *x2);
  RestoreBlocks(py, y, *y1, *y2);
}

void ArpackPEPSolver::GetResidual(PetscScalar l, const petsc::PetscParVector &x,
                                  petsc::PetscParVector &r) const
{
  // r = P(λ) x = (K + λ C + λ² M) x for eigenvalue λ.
  opM->Mult(x, r);
  r.Scale(l);
  opC->MultAdd(x, r);
  r.Scale(l);
  opK->MultAdd(x, r);
}

double ArpackPEPSolver::GetBackwardScaling(PetscScalar l) const
{
  // Make sure not to use norms from scaling as this can be confusing if they are different.
  // Note that SLEPc uses ||.||∞, not the 2-norm.
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
  double t = PetscAbsScalar(l);
  return normK + t * normC + t * t * normM;
}

void ArpackPEPSolver::ExtractEigenvector(PetscScalar l, petsc::PetscParVector &y,
                                         petsc::PetscParVector &x)
{
  // Select the most accurate x for y = [x₁; x₂] from the linearized eigenvalue problem.
  PetscScalar *py = GetBlocks(y, *y1, *y2);
  {
    if (opB)
    {
      y1->Normalize(*opB, *r);
    }
    else
    {
      y1->Normalize();
    }
    x.Copy(*y1);
  }
  RestoreBlocks(py, y, *y1, *y2);
}

PetscScalar *ArpackPEPSolver::GetBlocks(petsc::PetscParVector &v, petsc::PetscParVector &v1,
                                        petsc::PetscParVector &v2) const
{
  PetscInt n1 = v1.GetSize(), n2 = v2.GetSize();
  MFEM_VERIFY(n1 + n2 == v.GetSize(), "Unexpected size in PEP linearization!");
  PetscScalar *pv = v.GetArray();
  v1.PlaceArray(pv);
  v2.PlaceArray(pv + n1);
  return pv;
}

void ArpackPEPSolver::RestoreBlocks(PetscScalar *pv, petsc::PetscParVector &v,
                                    petsc::PetscParVector &v1,
                                    petsc::PetscParVector &v2) const
{
  v1.ResetArray();
  v2.ResetArray();
  v.RestoreArray(pv);
}

}  // namespace palace::arpack

#if defined(__GNUC__) && defined(__clang__)
#pragma GCC diagnostic pop
#endif

#endif

#endif
