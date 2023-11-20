// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "arpack.hpp"

#if defined(PALACE_WITH_ARPACK)

#include <algorithm>
#include <string>
#include <mfem.hpp>
// clang-format off
#include <parpack.hpp>  // ARPACK headers
#include <debug_c.hpp>
#include <stat_c.hpp>
// clang-format on
#include "linalg/divfree.hpp"
#include "utils/communication.hpp"

namespace
{

void CheckInfoAUPD(a_int info)
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

void CheckInfoEUPD(a_int info)
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

}  // namespace

namespace palace::arpack
{

// Base class methods

ArpackEigenvalueSolver::ArpackEigenvalueSolver(MPI_Comm comm, int print)
  : comm(comm), print(print)
{
  // Initialization.
  info = 0;
  nev = ncv = n = 0;
  rtol = 0.0;
  arpack_it = 0;
  which_type = WhichType::LARGEST_MAGNITUDE;
  gamma = delta = 1.0;
  sinvert = false;
  sigma = 0.0;

  opInv = nullptr;
  opProj = nullptr;
  opB = nullptr;

  // Configure debugging output.
  a_int logfill = 6, ndigit = -6, mgetv0 = 0;
  a_int _aupd = (print > 2) ? 1 : 0, _aup2 = (print > 2) ? 2 : ((print > 0) ? 1 : 0),
        _aitr = 0, _eigh = 0, _gets = 0, _apps = 0, _eupd = 0;
  debug_c(logfill, ndigit, mgetv0, _aupd, _aup2, _aitr, _eigh, _gets, _apps, _eupd, _aupd,
          _aup2, _aitr, _eigh, _gets, _apps, _eupd, _aupd, _aup2, _aitr, _eigh, _gets,
          _apps, _eupd);
  cstatn_c();
}

void ArpackEigenvalueSolver::SetOperators(const ComplexOperator &K,
                                          const ComplexOperator &M,
                                          EigenvalueSolver::ScaleType type)
{
  MFEM_ABORT("SetOperators not defined for base class ArpackEigenvalueSolver!");
}

void ArpackEigenvalueSolver::SetOperators(const ComplexOperator &K,
                                          const ComplexOperator &C,
                                          const ComplexOperator &M,
                                          EigenvalueSolver::ScaleType type)
{
  MFEM_ABORT("SetOperators not defined for base class ArpackEigenvalueSolver!");
}

void ArpackEigenvalueSolver::SetLinearSolver(const ComplexKspSolver &ksp)
{
  opInv = &ksp;
}

void ArpackEigenvalueSolver::SetDivFreeProjector(const DivFreeSolver &divfree)
{
  opProj = &divfree;
}

void ArpackEigenvalueSolver::SetBMat(const Operator &B)
{
  MFEM_VERIFY(!opB || opB->Height() == B.Height(),
              "Invalid modification of eigenvalue problem size!");
  opB = &B;
}

void ArpackEigenvalueSolver::SetNumModes(int num_eig, int num_vec)
{
  if (nev > 0 && num_eig != nev)
  {
    eig.reset();
    perm.reset();
    res.reset();
  }
  if (ncv > 0 && num_vec != ncv)
  {
    V.reset();
  }
  nev = num_eig;
  ncv = (num_vec > 0) ? num_vec : std::max(20, 2 * nev + 1);  // Default from SLEPc
}

void ArpackEigenvalueSolver::SetTol(double tol)
{
  rtol = tol;
}

void ArpackEigenvalueSolver::SetMaxIter(int max_it)
{
  arpack_it = max_it;
}

void ArpackEigenvalueSolver::SetWhichEigenpairs(EigenvalueSolver::WhichType type)
{
  which_type = type;
}

void ArpackEigenvalueSolver::SetShiftInvert(std::complex<double> s, bool precond)
{
  MFEM_VERIFY(!precond, "ARPACK eigenvalue solver does not support preconditioned "
                        "spectral transformation option!");
  sigma = s;
  sinvert = true;
}

void ArpackEigenvalueSolver::SetInitialSpace(const ComplexVector &v)
{
  MFEM_VERIFY(
      n > 0,
      "Must call SetOperators before using SetInitialSpace for ARPACK eigenvalue solver!");
  if (!r)
  {
    r = std::make_unique<std::complex<double>[]>(n);
  }
  MFEM_VERIFY(v.Size() == n, "Invalid size mismatch for provided initial space vector!");
  v.Get(r.get(), n);
  info = 1;
}

int ArpackEigenvalueSolver::SolveInternal(int n, std::complex<double> *r,
                                          std::complex<double> *V,
                                          std::complex<double> *eig, int *perm)
{
  MPI_Fint fcomm = MPI_Comm_c2f(comm);
  a_int iparam[11] = {0};
  iparam[0] = 1;                 // Exact shifts
  iparam[2] = (a_int)arpack_it;  // Maximum number of Arnoldi iterations
  iparam[3] = 1;                 // Block size
  iparam[4] = 0;                 // Number of converged Ritz values
  iparam[6] = sinvert ? 3 : 1;   // Problem mode

  ::arpack::bmat bmat_option =
      (opB) ? ::arpack::bmat::generalized : ::arpack::bmat::identity;

  ::arpack::which which_option = ::arpack::which::largest_magnitude;
  switch (which_type)
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

  // Allocate work arrays.
  a_int lworkl = 3 * ncv * ncv + 5 * ncv;
  auto workd = std::make_unique<std::complex<double>[]>(3 * n);
  auto workl = std::make_unique<std::complex<double>[]>(lworkl);
  auto rwork = std::make_unique<double[]>(ncv);

  // Begin RCI loop.
  a_int ido = 0, ainfo = (a_int)info, ipntr[14] = {0};
  while (true)
  {
    // Call complex problem driver.
    naupd(fcomm, ido, bmat_option, (a_int)n, which_option, (a_int)nev, rtol, r, (a_int)ncv,
          V, (a_int)n, iparam, ipntr, workd.get(), workl.get(), lworkl, rwork.get(), ainfo);
    CheckInfoAUPD(ainfo);

    // We never use pre-computed B * x in workd[ipntr[2] - 1].
    if (ido == 1 || ido == -1)
    {
      ApplyOp(&workd.get()[ipntr[0] - 1], &workd.get()[ipntr[1] - 1]);
    }
    else if (ido == 2)
    {
      ApplyOpB(&workd.get()[ipntr[0] - 1], &workd.get()[ipntr[1] - 1]);
    }
    else if (ido == 99)
    {
      break;
    }
    else
    {
      MFEM_ABORT("Internal error in ARPACK RCI interface!");
    }
  }

  // Print some log information.
  int num_it = (int)iparam[2];
  int num_conv = (int)iparam[4];
  if (print > 0)
  {
    Mpi::Print(comm,
               "\n ARPACK {} eigensolve {} ({:d} eigenpairs); iterations {:d}\n"
               " Total number of linear systems solved: {:d}\n"
               " Total number of linear solver iterations: {:d}\n",
               GetName(), (num_conv >= nev) ? "converged" : "finished", num_conv, num_it,
               opInv->NumTotalMult(), opInv->NumTotalMultIterations());
  }
  if (num_conv < nev)
  {
    Mpi::Warning(
        comm, "ARPACK eigenvalue solver found only {:d} of requested {:d} eigenvalues!\n",
        num_conv, nev);
  }

  // Postprocess eigenvalues and eigenvectors.
  a_int rvec = 1;
  ::arpack::howmny howmny_option = ::arpack::howmny::ritz_vectors;

  // Allocate eigenvalue storage and work arrays.
  auto select = std::make_unique<a_int[]>(ncv);
  auto workev = std::make_unique<std::complex<double>[]>(2 * ncv);

  // Call complex problem driver.
  neupd(fcomm, rvec, howmny_option, select.get(), eig, V, (a_int)n, sigma / gamma,
        workev.get(), bmat_option, (a_int)n, which_option, (a_int)nev, rtol, r, (a_int)ncv,
        V, (a_int)n, iparam, ipntr, workd.get(), workl.get(), lworkl, rwork.get(), ainfo);
  CheckInfoEUPD(ainfo);

  // Unscale and properly sort the eigenvalues.
  auto CompareReal = [&eig](const int &l, const int &r)
  { return eig[l].real() < eig[r].real(); };
  auto CompareImag = [&eig](const int &l, const int &r)
  { return eig[l].imag() < eig[r].imag(); };
  auto CompareAbs = [&eig](const int &l, const int &r)
  { return std::abs(eig[l]) < std::abs(eig[r]); };
  for (int i = 0; i < nev; i++)
  {
    eig[i] = eig[i] * gamma;
    perm[i] = i;
  }
  if (which_option == ::arpack::which::largest_real ||
      which_option == ::arpack::which::smallest_real)
  {
    std::sort(perm, perm + nev, CompareReal);
  }
  else if (which_option == ::arpack::which::largest_imaginary ||
           which_option == ::arpack::which::smallest_imaginary)
  {
    std::sort(perm, perm + nev, CompareImag);
  }
  else
  {
    std::sort(perm, perm + nev, CompareAbs);
  }

  return num_conv;
}

void ArpackEigenvalueSolver::CheckParameters() const
{
  MFEM_VERIFY(n > 0, "Operators are not set for ARPACK eigenvalue solver!");
  MFEM_VERIFY(nev > 0, "Number of requested modes is not positive!");
  MFEM_VERIFY(rtol > 0.0, "Eigensolver tolerance is not positive!");
  MFEM_VERIFY(opInv, "No linear solver provided for operator!");
}

std::complex<double> ArpackEigenvalueSolver::GetEigenvalue(int i) const
{
  MFEM_VERIFY(eig && i >= 0 && i < nev,
              "Out of range eigenpair requested (i = " << i << ", nev = " << nev << ")!");
  const int &j = perm.get()[i];
  return eig.get()[j];
}

void ArpackEigenvalueSolver::GetEigenvector(int i, ComplexVector &x) const
{
  MFEM_VERIFY(eig && i >= 0 && i < nev,
              "Out of range eigenpair requested (i = " << i << ", nev = " << nev << ")!");
  MFEM_VERIFY(x.Size() == n, "Invalid size mismatch for provided eigenvector!");
  const int &j = perm.get()[i];
  x.Set(V.get() + j * n, n);
}

double ArpackEigenvalueSolver::GetError(int i, EigenvalueSolver::ErrorType type) const
{
  MFEM_VERIFY(eig && i >= 0 && i < nev,
              "Out of range eigenpair requested (i = " << i << ", nev = " << nev << ")!");
  const int &j = perm.get()[i];
  switch (type)
  {
    case ErrorType::ABSOLUTE:
      return res.get()[j];
    case ErrorType::RELATIVE:
      return res.get()[j] / std::abs(eig.get()[j]);
    case ErrorType::BACKWARD:
      return res.get()[j] / GetBackwardScaling(eig.get()[j]);
  }
  return 0.0;
}

// EPS specific methods

ArpackEPSSolver::ArpackEPSSolver(MPI_Comm comm, int print)
  : ArpackEigenvalueSolver(comm, print)
{
  opK = opM = nullptr;
  normK = normM = 0.0;
}

void ArpackEPSSolver::SetOperators(const ComplexOperator &K, const ComplexOperator &M,
                                   EigenvalueSolver::ScaleType type)
{
  MFEM_VERIFY(!opK || opK->Height() == K.Height(),
              "Invalid modification of eigenvalue problem size!");
  bool first = (opK == nullptr);
  opK = &K;
  opM = &M;
  if (first && type != ScaleType::NONE)
  {
    normK = linalg::SpectralNorm(comm, *opK, opK->IsReal());
    normM = linalg::SpectralNorm(comm, *opM, opM->IsReal());
    MFEM_VERIFY(normK >= 0.0 && normM >= 0.0, "Invalid matrix norms for EPS scaling!");
    if (normK > 0 && normM > 0.0)
    {
      gamma = normK / normM;  // Store γ² for linear problem
      delta = 2.0 / normK;
    }
  }

  // Set up workspace.
  x.SetSize(opK->Height());
  y.SetSize(opK->Height());
  z.SetSize(opK->Height());
  n = opK->Height();
}

int ArpackEPSSolver::Solve()
{
  // Set some defaults (default maximum iterations from SLEPc).
  CheckParameters();
  HYPRE_BigInt N = linalg::GlobalSize(comm, z);
  if (ncv > N)
  {
    ncv = mfem::internal::to_int(N);
  }
  if (arpack_it <= 0)
  {
    arpack_it = std::max(300, mfem::internal::to_int(2 * N / ncv));
  }

  // Initialize if user did not provide an initial space.
  if (!r)
  {
    r = std::make_unique<std::complex<double>[]>(n);
    info = 0;
  }
  if (!info)
  {
    std::fill(r.get(), r.get() + n, 0.0);
  }

  // Allocate Arnoldi basis for the problem.
  if (!V)
  {
    V = std::make_unique<std::complex<double>[]>(n * ncv);
  }

  // Allocate storage for eigenvalues and residual norms.
  if (!eig)
  {
    eig = std::make_unique<std::complex<double>[]>(nev + 1);
    perm = std::make_unique<int[]>(nev);
    res = std::make_unique<double[]>(nev);
  }

  // Solve the generalized eigenvalue problem.
  int num_conv = SolveInternal(n, r.get(), V.get(), eig.get(), perm.get());

  // Compute the eigenpair residuals: || (K - λ M) x ||₂ for eigenvalue λ.
  for (int i = 0; i < nev; i++)
  {
    const std::complex<double> l = eig.get()[i];
    x.Set(V.get() + i * n, n);
    opK->Mult(x, y);
    opM->AddMult(x, y, -l);
    res.get()[i] = linalg::Norml2(comm, y);
  }

  // Reset for next solve.
  info = 0;
  return num_conv;
}

void ArpackEPSSolver::ApplyOp(const std::complex<double> *px,
                              std::complex<double> *py) const
{
  // Case 1: No spectral transformation (opInv = M⁻¹)
  //               y = M⁻¹ K x .
  // Case 2: Shift-and-invert spectral transformation (opInv = (K - σ M)⁻¹)
  //               y = (K - σ M)⁻¹ M x .
  x.Set(px, n);
  if (!sinvert)
  {
    opK->Mult(x, z);
    opInv->Mult(z, y);
    y *= 1.0 / gamma;
  }
  else
  {
    opM->Mult(x, z);
    opInv->Mult(z, y);
    y *= gamma;
  }
  if (opProj)
  {
    // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(comm, y));
    opProj->Mult(y);
    // Mpi::Print(" After projection: {:e}\n", linalg::Norml2(comm, y));
  }
  y.Get(py, n);
}

void ArpackEPSSolver::ApplyOpB(const std::complex<double> *px,
                               std::complex<double> *py) const
{
  MFEM_VERIFY(opB, "No B operator for weighted inner product in ARPACK solve!");
  x.Set(px, n);
  opB->Mult(x.Real(), y.Real());
  opB->Mult(x.Imag(), y.Imag());
  y *= delta * gamma;
  y.Get(py, n);
}

double ArpackEPSSolver::GetBackwardScaling(std::complex<double> l) const
{
  // Make sure not to use norms from scaling as this can be confusing if they are different.
  // Note that SLEPc uses ||.||∞, not the 2-norm.
  if (normK <= 0.0)
  {
    normK = linalg::SpectralNorm(comm, *opK, opK->IsReal());
  }
  if (normM <= 0.0)
  {
    normM = linalg::SpectralNorm(comm, *opM, opM->IsReal());
  }
  return normK + std::abs(l) * normM;
}

// PEP specific methods

ArpackPEPSolver::ArpackPEPSolver(MPI_Comm comm, int print)
  : ArpackEigenvalueSolver(comm, print)
{
  opK = opC = opM = nullptr;
  normK = normC = normM = 0.0;
}

void ArpackPEPSolver::SetOperators(const ComplexOperator &K, const ComplexOperator &C,
                                   const ComplexOperator &M,
                                   EigenvalueSolver::ScaleType type)
{
  MFEM_VERIFY(!opK || opK->Height() == K.Height(),
              "Invalid modification of eigenvalue problem size!");
  bool first = (opK == nullptr);
  opK = &K;
  opC = &C;
  opM = &M;
  if (first && type != ScaleType::NONE)
  {
    normK = linalg::SpectralNorm(comm, *opK, opK->IsReal());
    normC = linalg::SpectralNorm(comm, *opC, opC->IsReal());
    normM = linalg::SpectralNorm(comm, *opM, opM->IsReal());
    MFEM_VERIFY(normK >= 0.0 && normC >= 0.0 && normM >= 0.0,
                "Invalid matrix norms for PEP scaling!");
    if (normK > 0 && normC > 0.0 && normM > 0.0)
    {
      gamma = std::sqrt(normK / normM);
      delta = 2.0 / (normK + gamma * normC);
    }
  }

  // Set up workspace.
  x1.SetSize(opK->Height());
  x2.SetSize(opK->Height());
  y1.SetSize(opK->Height());
  y2.SetSize(opK->Height());
  z.SetSize(opK->Height());
  n = opK->Height();
}

int ArpackPEPSolver::Solve()
{
  // Set some defaults (from SLEPc ARPACK interface). The problem size is the size of the
  // 2x2 block linearized problem.
  CheckParameters();
  HYPRE_BigInt N = linalg::GlobalSize(comm, z);
  if (ncv > 2 * N)
  {
    ncv = mfem::internal::to_int(2 * N);
  }
  if (arpack_it <= 0)
  {
    arpack_it = std::max(300, mfem::internal::to_int(4 * N / ncv));
  }

  // Initialize if user did not provide an initial space.
  if (!r)
  {
    r = std::make_unique<std::complex<double>[]>(n);
    info = 0;
  }
  if (!info)
  {
    std::fill(r.get(), r.get() + n, 0.0);
  }
  auto s = std::make_unique<std::complex<double>[]>(2 * n);
  std::copy(r.get(), r.get() + n, s.get());
  std::fill(s.get() + n, s.get() + 2 * n, 0.0);

  // Allocate Arnoldi basis for original and linearized problem.
  if (!V)
  {
    V = std::make_unique<std::complex<double>[]>(n * ncv);
  }
  auto W = std::make_unique<std::complex<double>[]>(2 * n * ncv);

  // Allocate storage for eigenvalues and residual norms.
  if (!eig)
  {
    eig = std::make_unique<std::complex<double>[]>(nev + 1);
    perm = std::make_unique<int[]>(nev + 1);
    res = std::make_unique<double[]>(nev + 1);
  }

  // Solve the linearized eigenvalue problem.
  int num_conv = SolveInternal(2 * n, s.get(), W.get(), eig.get(), perm.get());

  // Extract the eigenvector from the linearized problem and compute the eigenpair
  // residuals: || P(λ) x ||₂ = || (K + λ C + λ² M) x ||₂ for eigenvalue λ.
  for (int i = 0; i < nev; i++)
  {
    const std::complex<double> &l = eig.get()[i];
    ExtractEigenvector(l, W.get() + i * 2 * n, V.get() + i * n);
    x1.Set(V.get() + i * n, n);
    opK->Mult(x1, y1);
    opC->AddMult(x1, y1, l);
    opM->AddMult(x1, y1, l * l);
    res.get()[i] = linalg::Norml2(comm, y1);
  }

  // Reset for next solve.
  info = 0;
  return num_conv;
}

void ArpackPEPSolver::ApplyOp(const std::complex<double> *px,
                              std::complex<double> *py) const
{
  // Case 1: No spectral transformation (opInv = M⁻¹)
  //               y = L₁⁻¹ L₀ x .
  // Case 2: Shift-and-invert spectral transformation (opInv = P(σ)⁻¹)
  //               y = (L₀ - σ L₁)⁻¹ L₁ x .
  // With:
  //               L₀ = [ -K  0 ]    L₁ = [ C  M ]
  //                    [  0  M ] ,       [ M  0 ] .
  x1.Set(px, n);
  x2.Set(px + n, n);
  if (!sinvert)
  {
    y1 = x2;
    if (opProj)
    {
      // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(comm, y1));
      opProj->Mult(y1);
      // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(comm, y1));
    }

    opK->Mult(x1, z);
    opC->AddMult(x2, z, std::complex<double>(gamma, 0.0));
    opInv->Mult(z, y2);
    y2 *= -1.0 / (gamma * gamma);
    if (opProj)
    {
      // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(comm, y2));
      opProj->Mult(y2);
      // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(comm, y2));
    }
  }
  else
  {
    y2.AXPBYPCZ(sigma, x1, gamma, x2, 0.0);  // Just temporarily
    opM->Mult(y2, z);
    opC->AddMult(x1, z, std::complex<double>(1.0, 0.0));
    opInv->Mult(z, y1);
    y1 *= -gamma;
    if (opProj)
    {
      // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(comm, y1));
      opProj->Mult(y1);
      // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(comm, y1));
    }

    y2.AXPBYPCZ(sigma / gamma, y1, 1.0, x1, 0.0);
    if (opProj)
    {
      // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(comm, y2));
      opProj->Mult(y2);
      // Mpi::Print(" Before projection: {:e}\n", linalg::Norml2(comm, y2));
    }
  }
  y1.Get(py, n);
  y2.Get(py + n, n);
}

void ArpackPEPSolver::ApplyOpB(const std::complex<double> *px,
                               std::complex<double> *py) const
{
  MFEM_VERIFY(opB, "No B operator for weighted inner product in ARPACK solve!");
  x1.Set(px, n);
  x2.Set(px + n, n);
  opB->Mult(x1.Real(), y1.Real());
  opB->Mult(x1.Imag(), y1.Imag());
  opB->Mult(x2.Real(), y2.Real());
  opB->Mult(x2.Imag(), y2.Imag());
  y1 *= delta * gamma * gamma;
  y2 *= delta * gamma * gamma;
  y1.Get(py, n);
  y2.Get(py + n, n);
}

double ArpackPEPSolver::GetBackwardScaling(std::complex<double> l) const
{
  // Make sure not to use norms from scaling as this can be confusing if they are different.
  // Note that SLEPc uses ||.||∞, not the 2-norm.
  if (normK <= 0.0)
  {
    normK = linalg::SpectralNorm(comm, *opK, opK->IsReal());
  }
  if (normC <= 0.0)
  {
    normC = linalg::SpectralNorm(comm, *opC, opC->IsReal());
  }
  if (normM <= 0.0)
  {
    normM = linalg::SpectralNorm(comm, *opM, opM->IsReal());
  }
  double t = std::abs(l);
  return normK + t * normC + t * t * normM;
}

void ArpackPEPSolver::ExtractEigenvector(std::complex<double> l,
                                         const std::complex<double> *py,
                                         std::complex<double> *px) const
{
  // Select the most accurate x for y = [x₁; x₂] from the linearized eigenvalue problem. Or,
  // just take x = x₁.
  x1.Set(py, n);
  if (opB)
  {
    linalg::Normalize(comm, x1, *opB, y1);
  }
  else
  {
    linalg::Normalize(comm, x1);
  }
  x1.Get(px, n);
}

}  // namespace palace::arpack

#endif
