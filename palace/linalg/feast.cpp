// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "feast.hpp"

#if defined(PALACE_WITH_SLEPC)

#include <memory>
#include <string>
#include <vector>
#include <petsc.h>
#include <slepc.h>
#include <slepcblaslapack.h>
#include "linalg/divfree.hpp"
#include "linalg/ksp.hpp"
#include "linalg/pc.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

static PetscErrorCode __mat_apply_FEAST_EPS(Mat, Vec, Vec);
static PetscErrorCode __mat_apply_FEAST_PEP(Mat, Vec, Vec);

namespace palace::feast
{

namespace internal
{

// Linear solver helper class

class FeastLinearSolver
{
public:
  PetscScalar zk, wk;
  KspSolver ksp;
  KspPreconditioner pc;
  const petsc::PetscParMatrix *opK, *opC, *opM;  // Reference to EVP operators (not owned)

private:
  SpaceOperator &spaceop;  // Reference to spatial discretization (not owned)
  std::unique_ptr<petsc::PetscParMatrix> A;
  std::vector<std::unique_ptr<mfem::Operator>> P, AuxP;

public:
  FeastLinearSolver(int k, MPI_Comm comm, const IoData &iodata, SpaceOperator &sp)
    : zk(0.0), wk(0.0), ksp(comm, iodata, "ksp" + std::to_string(k + 1) + "_"),
      pc(iodata, sp.GetDbcMarker(), sp.GetNDSpaces(), &sp.GetH1Spaces()), spaceop(sp)
  {
    ksp.SetTabLevel(1);
    ksp.SetPrintOptions(false);
    ksp.SetPreconditioner(pc);
    opK = opC = opM = nullptr;
  }

  void SetOperators(PetscScalar z, PetscScalar w, const petsc::PetscParMatrix &K,
                    const petsc::PetscParMatrix &M)
  {
    zk = z;
    wk = w;
    opK = &K;
    opM = &M;
    {
      Mat A_;
      MPI_Comm comm = K.GetComm();
      PetscInt n = K.GetNumRows();
      PalacePetscCall(
          MatCreateShell(comm, n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &A_));
      PalacePetscCall(
          MatShellSetOperation(A_, MATOP_MULT,
                               (void (*)()) static_cast<PetscErrorCode (*)(Mat, Vec, Vec)>(
                                   &__mat_apply_FEAST_EPS)));
      A = std::make_unique<petsc::PetscParMatrix>(A_, false);  // Inherits the PETSc Mat
      ksp.SetOperator(*A);
    }
    const double sigma = PetscSqrtReal(PetscAbsScalar(zk));
    constexpr bool print = false;
    spaceop.GetPreconditionerMatrix(sigma, P, AuxP, print);
    pc.SetOperator(P, &AuxP);
  }

  void SetOperators(PetscScalar z, PetscScalar w, const petsc::PetscParMatrix &K,
                    const petsc::PetscParMatrix &C, const petsc::PetscParMatrix &M,
                    KspPreconditioner *op = nullptr)
  {
    zk = z;
    wk = w;
    opK = &K;
    opC = &C;
    opM = &M;
    {
      Mat A_;
      MPI_Comm comm = K.GetComm();
      PetscInt n = K.GetNumRows();
      PalacePetscCall(
          MatCreateShell(comm, n, n, PETSC_DECIDE, PETSC_DECIDE, (void *)this, &A_));
      PalacePetscCall(
          MatShellSetOperation(A_, MATOP_MULT,
                               (void (*)()) static_cast<PetscErrorCode (*)(Mat, Vec, Vec)>(
                                   &__mat_apply_FEAST_PEP)));
      A = std::make_unique<petsc::PetscParMatrix>(A_, false);  // Inherits the PETSc Mat
      ksp.SetOperator(*A);
    }
    const double sigma = PetscAbsScalar(zk);
    constexpr bool print = false;
    spaceop.GetPreconditionerMatrix(sigma, P, AuxP, print);
    pc.SetOperator(P, &AuxP);
  }

  void Mult(const PetscScalar *eig, const petsc::PetscDenseMatrix &X,
            const petsc::PetscDenseMatrix &R, petsc::PetscDenseMatrix &Q,
            petsc::PetscParVector &v, bool *converged, PetscReal gamma) const
  {
    // Solve P(zₖ) Qₖ = R, Q += wₖ (X - Qₖ) (zₖ I - Λ)⁻¹ (residual-inverse iteration). Note:
    // Q may have g.t. m0 columns, but we just use the first m0 for the result (X should
    // have exactly m0 columns).
    PetscInt m0 = X.GetGlobalNumCols();
    PetscInt M = Q.GetGlobalNumCols() / (2 * m0);
    MFEM_VERIFY(M == 1 || M == 2,
                "FEAST eigensolver only supports up to 2 subspace moments!");
    for (PetscInt j = 0; j < m0; j++)
    {
      const petsc::PetscParVector x = X.GetColumnRead(j);
      if (converged && converged[j])
      {
        // When R[j] is converged, Q[j] += wₖ/(zₖ - λₖ) X[j] (with Qₖ[j] = 0) .
        v.AXPBY(wk / (zk / gamma - eig[j]), x, 0.0);
      }
      else
      {
        const petsc::PetscParVector r = R.GetColumnRead(j);
        ksp.Mult(r, v);
        v.AXPBY(wk / (zk / gamma - eig[j]), x, -wk / (zk / gamma - eig[j]));
        R.RestoreColumnRead(j, r);
      }
      X.RestoreColumnRead(j, x);

      petsc::PetscParVector q = Q.GetColumn(j);
      q.AXPY(1.0, v);
      Q.RestoreColumn(j, q);
      if (M > 1)
      {
        petsc::PetscParVector q = Q.GetColumn(j + m0);
        q.AXPY(zk / gamma, v);
        Q.RestoreColumn(j + m0, q);
      }
    }
  }

  PetscScalar Mult(const petsc::PetscDenseMatrix &X, petsc::PetscParVector &r,
                   petsc::PetscParVector &v) const
  {
    // Solve P(zₖ) Qₖ = P'(zₖ) X, sum += wₖ tr(Xᵀ Qₖ) for estimating the eigenvalue count
    // inside of the contour.
    PetscInt m0 = X.GetGlobalNumCols();
    PetscScalar sum = 0.0;
    for (PetscInt j = 0; j < m0; j++)
    {
      const petsc::PetscParVector x = X.GetColumnRead(j);
      opM->Mult(x, r);
      if (opC)
      {
        r.Scale(zk);
        opC->MultAdd(x, r);
      }
      ksp.Mult(r, v);
      sum += x.TransposeDot(v);
      X.RestoreColumnRead(j, x);
    }
    return wk * sum;
  }
};

}  // namespace internal

// Base class methods

FeastEigenSolver::FeastEigenSolver(MPI_Comm comm, const IoData &iodata,
                                   SpaceOperator &spaceop, int np, int print_lvl)
{
  // Initialization.
  print = print_lvl;
  info = 0;
  nev = m0 = mQ = 0;
  M = iodata.solver.eigenmode.feast_moments;
  MFEM_VERIFY(M == 1 || M == 2,
              "FEAST eigensolver only supports up to 2 subspace moments!");
  rtol = 0.0;
  max_it = 0;
  gamma = delta = 1.0;
  bl = tr = 0.0;
  real_threshold = imag_threshold = false;

  eig = nullptr;
  perm = nullptr;
  X = nullptr;
  res = nullptr;
  r0 = nullptr;
  opProj = nullptr;
  opB = nullptr;

  // Construct the linear solvers for each quadrature point.
  opInv.reserve(np);
  for (int k = 0; k < np; k++)
  {
    opInv.emplace_back(k, comm, iodata, spaceop);
  }
}

FeastEigenSolver::~FeastEigenSolver()
{
  delete[] eig;
  delete[] perm;
  delete[] res;
  delete X;
  delete r0;
}

void FeastEigenSolver::SetOperators(const petsc::PetscParMatrix &K,
                                    const petsc::PetscParMatrix &M,
                                    EigenSolverBase::ScaleType type)
{
  MFEM_ABORT("SetOperators not defined for base class FeastEigenSolver!");
}

void FeastEigenSolver::SetOperators(const petsc::PetscParMatrix &K,
                                    const petsc::PetscParMatrix &C,
                                    const petsc::PetscParMatrix &M,
                                    EigenSolverBase::ScaleType type)
{
  MFEM_ABORT("SetOperators not defined for base class FeastEigenSolver!");
}

void FeastEigenSolver::SetProjector(const DivFreeSolver &divfree)
{
  opProj = &divfree;
}

void FeastEigenSolver::SetNumModes(int numeig, int numvec)
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
  if (m0 > 0 && numvec != m0)
  {
    delete X;
    X = nullptr;
  }
  nev = numeig;
  if (numvec > 0)
  {
    m0 = numvec;
  }
  else
  {
    if (nev <= 3)
    {
      m0 = std::max(nev + 2, 2 * nev);  // Just a guess for subspace dimension
    }
    else
    {
      m0 = std::max(nev + 3, nev + (nev + 1) / 2);
    }
  }
  mQ = 2 * M * m0;  // Real-valued basis splitting leads to factor of 2
}

void FeastEigenSolver::SetTol(double tol)
{
  rtol = tol;
}

void FeastEigenSolver::SetMaxIter(int maxits)
{
  max_it = maxits;
}

void FeastEigenSolver::SetContour(double blr, double bli, double trr, double tri,
                                  bool filter_small_real, bool filter_small_imag)
{
  MFEM_VERIFY(blr <= trr && bli <= tri,
              "Integration contour must be defined by bottom-left and top-right "
              "points in the complex plane!");
  bl = blr + PETSC_i * bli;
  tr = trr + PETSC_i * tri;
  real_threshold = filter_small_real;
  imag_threshold = filter_small_imag;
}

void FeastEigenSolver::SetBMat(const petsc::PetscParMatrix &B)
{
  opB = &B;
}

void FeastEigenSolver::SetInitialSpace(const petsc::PetscParVector &v)
{
  if (!r0)
  {
    r0 = new petsc::PetscParVector(v);
  }
  else
  {
    MFEM_VERIFY(v.GetSize() == r0->GetSize(),
                "Invalid modification of eigenvalue problem size!");
    r0->Copy(v);
  }
  info = 1;
}

int FeastEigenSolver::SolveInternal(RG rg)
{
  // Allocate space for subspace and residuals. R is constructed with mQ columns for
  // computing products of form R = A Q during projection.
  MFEM_VERIFY(X && X->GetGlobalNumCols() == m0,
              "Unexpected number of eigenvector columns in FEAST solver!");
  MPI_Comm comm = X->GetComm();
  PetscInt n = X->GetNumRows();
  petsc::PetscDenseMatrix R(comm, n, PETSC_DECIDE, PETSC_DECIDE, mQ, nullptr);
  petsc::PetscDenseMatrix Q(comm, n, PETSC_DECIDE, PETSC_DECIDE, mQ, nullptr);

  // Allocate other workspace variables.
  PetscInt *inside = new PetscInt[m0];
  bool *converged = new bool[m0];
  if (!eig)
  {
    eig = new PetscScalar[m0];
    perm = new PetscInt[m0];
    res = new PetscReal[m0];
  }
  for (PetscInt j = 0; j < m0; j++)
  {
    res[j] = -1.0;
  }
  mfem::Vector qr(n), qi(n);

#if 0
  // XX TODO: Stochastic estimates
  bool est_stochastic = true;
  if (est_stochastic)
  {
    X->SetRandomReal(0, m0);
    if (info)
    {
      for (PetscInt j = 0; j < m0; j++)
      {
        // Ensure homogeneous Dirichlet BC are satisfied by the subspace.
        petsc::PetscParVector x = X->GetColumn(j);
        x.PointwiseMult(*r0, false);
        X->RestoreColumn(j, x);
      }
    }
    X->SetRandomSign(0, m0, true);

    PetscScalar sum = 0;
    petsc::PetscParVector r = R.GetColumn(0);  // Just for workspace
    for (const auto &op : opInv)
    {
      sum += op.Mult(*X, r, *r0);
    }
    R.RestoreColumn(0, r);
    PetscInt m = (PetscInt)PetscCeilReal(PetscAbsScalar(sum)/(PetscReal)m0);

    // Debug
    Mpi::Print("Eigenvalue estimate: {:d}\n", m);
  }
#endif

  // Initialize the subspace.
  Q.SetRandom(0, mQ / 2);
  if (info)
  {
    petsc::PetscParVector q = Q.GetColumn(0);
    q.Copy(*r0);
    Q.RestoreColumn(0, q);
    for (PetscInt j = 1; j < mQ / 2; j++)
    {
      // Ensure homogeneous Dirichlet BC are satisfied by the starting subspace.
      petsc::PetscParVector q = Q.GetColumn(j);
      q.PointwiseMult(*r0, false);
      Q.RestoreColumn(j, q);
    }
  }

  // Begin main FEAST loop.
  int it = 0, nconv, ninside;
  while (true)
  {
    // Orthonormalize the (real-valued) basis Q.
    {
      bool mgs = false, cgs2 = true;
      for (PetscInt j = 0; j < mQ / 2; j++)
      {
        petsc::PetscParVector q1 = Q.GetColumn(j);
        q1.GetToVectors(qr, qi);
        if (opProj)
        {
          opProj->Mult(qr);
          opProj->Mult(qi);
        }
        q1.SetFromVector(qr);
        Q.RestoreColumn(j, q1);

        petsc::PetscParVector q2 = Q.GetColumn(j + mQ / 2);
        q2.SetFromVector(qi);
        Q.RestoreColumn(j + mQ / 2, q2);
      }
      for (PetscInt j = 0; j < mQ; j++)
      {
        if (opB)
        {
          Q.OrthonormalizeColumn(j, mgs, cgs2, *opB, *r0);
        }
        else
        {
          Q.OrthonormalizeColumn(j, mgs, cgs2);
        }
      }
    }

    // Form and solve the projected EVP. Select the m0 best eigenpair candidates and
    // reconstruct the full-dimensional eigenvectors.
    SolveProjectedProblem(Q, R, *X, eig);

    // Update the eigenpair residuals and check convergence. Residual calculation and
    // convergence tests occur in the unscaled space.
    nconv = ninside = 0;
    bool check = true;
    PetscReal rmin = mfem::infinity(), rmax = 0.0;
    PetscInt jmin = -1, jmax = -1;
    if (rg)
    {
      PalacePetscCall(RGCheckInside(rg, m0, eig, nullptr, inside));
    }
    else
    {
      for (PetscInt j = 0; j < m0; j++)
      {
        inside[j] = true;
      }
    }
    for (PetscInt j = 0; j < m0; j++)
    {
      PetscScalar sigma = eig[j] * gamma;
      petsc::PetscParVector x = X->GetColumn(j);
      petsc::PetscParVector r = R.GetColumn(j);
      if (opB)
      {
        x.Normalize(*opB, *r0);
      }
      else
      {
        x.Normalize();
      }
      GetResidual(sigma, x, r);
      PetscReal res = r.Norml2() / (x.Norml2() * PetscAbsScalar(sigma));
      // PetscReal res = r.Norml2()/x.Norml2();
      X->RestoreColumn(j, x);
      R.RestoreColumn(j, r);
      if (res < rtol)
      {
        // Mark converged even for eigenvalues outside the contour.
        converged[j] = true;
        nconv++;
        if (res > rmax)
        {
          rmax = res;
          jmax = j;
        }
      }
      else
      {
        converged[j] = false;
        if (res < rmin)
        {
          rmin = res;
          jmin = j;
        }
      }
      if (inside[j] >= 0)
      {
        ninside++;
        if (!converged[j])
        {
          check = false;  // Only finish when inside eigenvalues are converged
        }
      }

      // Debug
      // Mpi::Print(comm, " res[{:d}] = {:e} (eig = {:+e}{:+e}i, inside = {:d})\n",
      //            j, res, PetscRealPart(sigma),
      //            PetscImaginaryPart(sigma), inside[j]);
    }
    if (print > 0)
    {
      if (ninside > 0 || nconv > 0)
      {
        if (jmin >= 0)
        {
          Mpi::Print(comm,
                     "  {:d} FEAST inside={:d} converged={:d} first "
                     "unconverged value (error) {:+.3e}{:+.3e}i ({:.6e})\n",
                     it, ninside, nconv, PetscRealPart(eig[jmin] * gamma),
                     PetscImaginaryPart(eig[jmin] * gamma), rmin);
        }
        else
        {
          Mpi::Print(comm,
                     "  {:d} FEAST inside={:d} converged={:d} last "
                     "converged value (error) {:+.3e}{:+.3e}i ({:.6e})\n",
                     it, ninside, nconv, PetscRealPart(eig[jmax] * gamma),
                     PetscImaginaryPart(eig[jmax] * gamma), rmax);
        }
      }
      else
      {
        Mpi::Print(comm, "  {:d} FEAST inside=0\n", it);
      }
    }
    // Check convergence: All inside must be converged + any outside if user specified nev
    // too large.
    if ((check && nconv >= nev) || it == max_it)
    {
      break;
    }

    // Update subspace with contour integral (accumulates to first M*m0 columns of Q).
    Q.Scale(0.0);
    for (const auto &op : opInv)
    {
      op.Mult(eig, *X, R, Q, *r0, converged, gamma);
    }
    it++;
  }

  // Print some log information.
  if (print > 0)
  {
    Mpi::Print(comm,
               "\n FEAST {} eigensolve {} ({:d} eigenpairs); iterations {:d}\n"
               " Total number of linear systems solved: {:d}\n"
               " Total number of linear solver iterations: {:d}\n",
               GetName(), (it == max_it) ? "finished" : "converged", nconv, it,
               GetTotalKspMult(), GetTotalKspIter());
  }
  if (it == max_it)
  {
    Mpi::Warning(comm,
                 "FEAST eigenvalue solver reached maximum {:d} "
                 "iterations!\nFound {:d} converged eigenvales of requested {:d}!\n",
                 it, nconv, nev);
  }

  // Unscale and sort the eigenvalues in ascending order.
  auto CompareAbs = [converged, this](const PetscInt &l, const PetscInt &r)
  {
    if (!converged[l] && converged[r])
    {
      return false;
    }
    else if (converged[l] && !converged[r])
    {
      return true;
    }
    return (PetscAbsScalar(eig[l]) < PetscAbsScalar(eig[r]));
  };
  for (PetscInt j = 0; j < m0; j++)
  {
    eig[j] = eig[j] * gamma;
    perm[j] = j;
  }
  std::sort(perm, perm + m0, CompareAbs);

  // Cleanup.
  delete[] inside;
  delete[] converged;

  // Reset for next solve.
  info = 0;
  return nconv;
}

void FeastEigenSolver::CheckParameters()
{
  MFEM_VERIFY(nev > 0, "Number of requested modes is not positive!");
  MFEM_VERIFY(rtol > 0.0, "Eigensolver tolerance is not positive!");
  MFEM_VERIFY(!(bl == 0.0 && tr == 0.0), "Integration contour has not been defined!");
  if (max_it <= 0)
  {
    max_it = 15;
  }
}

RG FeastEigenSolver::ConfigureRG(PetscScalar *&z, PetscScalar *&w)
{
  int np = static_cast<int>(opInv.size());
  if (np == 1)
  {
    z = new PetscScalar[np];
    w = new PetscScalar[np];
    z[0] = 0.5 * (bl + tr) / gamma;  // User should pass in bl = tr = target
    w[0] = 1.0;
    return nullptr;
  }
  else
  {
    RG rg;
    PalacePetscCall(RGCreate(PETSC_COMM_SELF, &rg));
    MFEM_VERIFY(PetscRealPart(tr - bl) > 0.0 && PetscImaginaryPart(tr - bl) > 0.0,
                "Contour must have nonzero and finite aspect ratio!");
    PetscScalar c = 0.5 * (bl + tr) / gamma;
    PetscReal r = 0.5 * PetscRealPart(tr - bl) / gamma;
    PetscReal vscale = 0.5 * PetscImaginaryPart(tr - bl) / (r * gamma);
    PalacePetscCall(RGSetType(rg, RGELLIPSE));
    PalacePetscCall(RGEllipseSetParameters(rg, c, r, vscale));
    // MFEM_VERIFY(opInv.size() % 4 == 0,
    //             "Number of contour quadrature points for rectangular region
    //             must be evenly divisible by 4!");
    // PalacePetscCall(RGSetType(rg, RGINTERVAL));
    // PalacePetscCall(RGIntervalSetEndpoints(rg, PetscRealPart(bl)/gamma,
    // PetscRealPart(tr)/gamma, // PetscImaginaryPart(bl)/gamma,
    // PetscImaginaryPart(tr)/gamma));

    z = new PetscScalar[np];
    w = new PetscScalar[np];
    if (PetscImaginaryPart(c) == 0.0 || PetscRealPart(c) == 0.0)
    {
      // Contour is symmetric about an axis and we place the first quadrature point at θ
      // = -π/2 (imaginary-axis symmetry) or θ = π (real-axis symmetry).
      PetscReal shift = (PetscRealPart(c) == 0.0) ? -0.5 * PETSC_PI : PETSC_PI;
      for (int k = 0; k < np; k++)
      {
        PetscReal theta = 2.0 * PETSC_PI * k / (PetscReal)np + shift;
        z[k] = c + r * (PetscCosReal(theta) + PETSC_i * vscale * PetscSinReal(theta));
        w[k] = r * (vscale * PetscCosReal(theta) + PETSC_i * PetscSinReal(theta)) /
               (PetscReal)np;
      }
    }
    else
    {
      PetscScalar *zn = new PetscScalar[np];
      PalacePetscCall(RGComputeQuadrature(rg, RG_QUADRULE_TRAPEZOIDAL, np, z, zn, w));
      delete[] zn;
    }
    return rg;
  }
}

PetscInt *FeastEigenSolver::SortEigenvalues(const PetscScalar *eig_, PetscInt m) const
{
  PetscReal rthresh = (real_threshold) ? 0.01 * PetscRealPart(bl) / gamma : 0.0;
  PetscReal ithresh = (imag_threshold) ? 0.01 * PetscImaginaryPart(bl) / gamma : 0.0;
  PetscScalar target = 0.5 * (bl + tr) / gamma;
  PetscReal vscale =
      (bl == tr) ? 1.0 : PetscImaginaryPart(tr - bl) / PetscRealPart(tr - bl);
  auto CompareTargetAbs =
      [eig_, rthresh, ithresh, target, vscale](const PetscInt &l, const PetscInt &r)
  {
    PetscReal lr = PetscAbsReal(PetscRealPart(eig_[l]));
    PetscReal li = PetscAbsReal(PetscImaginaryPart(eig_[l]));
    PetscReal rr = PetscAbsReal(PetscRealPart(eig_[r]));
    PetscReal ri = PetscAbsReal(PetscImaginaryPart(eig_[r]));
    if ((li < ithresh && ri >= ithresh) || (lr < rthresh && rr >= rthresh))
    {
      return false;
    }
    else if ((li >= ithresh && ri < ithresh) || (lr >= rthresh && rr < rthresh))
    {
      return true;
    }
    PetscScalar dl = eig_[l] - target;
    PetscScalar dr = eig_[r] - target;
    PetscReal vl = PetscRealPart(dl) * PetscRealPart(dl) +
                   PetscImaginaryPart(dl) * PetscImaginaryPart(dl) / (vscale * vscale);
    PetscReal vr = PetscRealPart(dr) * PetscRealPart(dr) +
                   PetscImaginaryPart(dr) * PetscImaginaryPart(dr) / (vscale * vscale);
    return (vl < vr);
  };
  PetscInt *perm_ = new PetscInt[m];
  for (PetscInt i = 0; i < m; i++)
  {
    perm_[i] = i;
  }
  std::sort(perm_, perm_ + m, CompareTargetAbs);
  return perm_;
}

void FeastEigenSolver::BVMatProjectInternal(const petsc::PetscDenseMatrix &Q,
                                            const petsc::PetscParMatrix &A,
                                            petsc::PetscDenseMatrix &Ar,
                                            petsc::PetscDenseMatrix &R,
                                            PetscReal scale) const
{
  // Compute Ar = Qᴴ A Q. We assume Q is real and thus the result is complex symmetric if A
  // is symmetric. Ar is replicated across all processes(sequential mQ x mQ matrix).
  MFEM_VERIFY(A.GetSymmetric() && Ar.GetSymmetric(),
              "BVMatProjectInternal is specialized for symmetric matrices!");
  MFEM_VERIFY(Q.GetGlobalNumCols() == mQ && R.GetGlobalNumCols() == mQ &&
                  Ar.GetNumRows() == mQ && Ar.GetNumCols() == mQ,
              "Unexpected number of basis columns in FEAST solver!");
  mfem::Vector qr(Q.GetNumRows());
  for (PetscInt j = 0; j < mQ; j++)
  {
    const petsc::PetscParVector q = Q.GetColumnRead(j);
    petsc::PetscParVector r = R.GetColumn(j);
    q.GetToVector(qr);
    A.Mult(qr, r);
    Q.RestoreColumnRead(j, q);
    R.RestoreColumn(j, r);
  }
  PetscInt n = A.GetNumRows();
  const PetscScalar *pQ = Q.GetArrayRead(), *pR = R.GetArrayRead();
  petsc::PetscDenseMatrix locQ(n, mQ, const_cast<PetscScalar *>(pQ));
  petsc::PetscDenseMatrix locR(n, mQ, const_cast<PetscScalar *>(pR));
  locQ.MatTransposeMult(locR, Ar);  // Qᴴ = Qᵀ
  Q.RestoreArrayRead(pQ);
  R.RestoreArrayRead(pR);

  // Global reduction over all processes.
  PetscScalar *pAr = Ar.GetArray();
  Mpi::GlobalSum(mQ * mQ, pAr, Q.GetComm());
  Ar.RestoreArray(pAr);
  Ar.Scale(scale);
}

int FeastEigenSolver::GetTotalKspMult() const
{
  int ksp_mult = 0;
  for (const auto &op : opInv)
  {
    ksp_mult += op.ksp.GetTotalNumMult();
  }
  return ksp_mult;
}

int FeastEigenSolver::GetTotalKspIter() const
{
  int ksp_it = 0;
  for (const auto &op : opInv)
  {
    ksp_it += op.ksp.GetTotalNumIter();
  }
  return ksp_it;
}

void FeastEigenSolver::GetEigenvalue(int i, double &eigr, double &eigi) const
{
  MFEM_VERIFY(eig && i >= 0 && i < m0,
              "Out of range eigenpair requested (i = " << i << ", m0 = " << m0 << ")!");
  const int &j = perm[i];
  eigr = PetscRealPart(eig[j]);
  eigi = PetscImaginaryPart(eig[j]);
}

void FeastEigenSolver::GetEigenvector(int i, petsc::PetscParVector &v) const
{
  MFEM_VERIFY(eig && i >= 0 && i < m0,
              "Out of range eigenpair requested (i = " << i << ", m0 = " << m0 << ")!");
  const int &j = perm[i];
  const petsc::PetscParVector x = X->GetColumnRead(j);
  v.Copy(x);
  X->RestoreColumnRead(j, x);
}

void FeastEigenSolver::GetError(int i, EigenSolverBase::ErrorType type, double &err) const
{
  MFEM_VERIFY(eig && i >= 0 && i < m0,
              "Out of range eigenpair requested (i = " << i << ", m0 = " << m0 << ")!");
  const int &j = perm[i];
  if (res[j] <= 0.0)
  {
    const petsc::PetscParVector x = X->GetColumnRead(j);
    GetResidual(eig[j], x, *r0);
    res[j] = r0->Norml2() / x.Norml2();
    X->RestoreColumnRead(j, x);
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
    default:
      MFEM_ABORT("Eigenpair error type not implemented!");
      break;
  }
}

// EPS specific methods

FeastEPSSolver::FeastEPSSolver(MPI_Comm comm, const IoData &iodata, SpaceOperator &spaceop,
                               int np, int print_lvl)
  : FeastEigenSolver(comm, iodata, spaceop, np, print_lvl)
{
  opK = opM = nullptr;
  normK = normM = 0.0;
  AQ = BQ = XQ = XQ0 = nullptr;
}

void FeastEPSSolver::SetOperators(const petsc::PetscParMatrix &K,
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
    MFEM_VERIFY(normK > 0.0 && normM > 0.0, "Invalid matrix norms for EPS scaling!");
    gamma = normK / normM;  // Store γ² for linear problem
    delta = 2.0 / normK;
  }
}

int FeastEPSSolver::Solve()
{
  // Check inputs.
  CheckParameters();
  MFEM_VERIFY(opK && opM, "Operators are not set for FeastEPSSolver!");

  // Allocate storage for eigenvectors.
  MPI_Comm comm = opK->GetComm();
  if (!X)
  {
    X = new petsc::PetscDenseMatrix(comm, opK->GetNumRows(), PETSC_DECIDE, PETSC_DECIDE, m0,
                                    nullptr);
  }
  if (!r0)
  {
    r0 = new petsc::PetscParVector(*opK);
  }

  // Allocate sequential matrices for the projected generalized eigenvalue problems at each
  // iteration, and associated eigenvectors.
  AQ = new petsc::PetscDenseMatrix(mQ, mQ, nullptr);
  AQ->CopySymmetry(*opK);
  BQ = new petsc::PetscDenseMatrix(mQ, mQ, nullptr);
  BQ->CopySymmetry(*opM);
  XQ = new petsc::PetscDenseMatrix(mQ, mQ, nullptr);
  XQ0 = new petsc::PetscDenseMatrix(mQ, m0, nullptr);

  // Create region object for integration contour and configure the linear solvers at each
  // integration point. The linear solves use the unscaled space.
  PetscScalar *z, *w;
  RG rg = ConfigureRG(z, w);
  Mpi::Print(comm, "Quadrature points for FEAST contour\n");
  for (int k = 0; k < static_cast<int>(opInv.size()); k++)
  {
    Mpi::Print(comm, " {:d}: zₖ = {:+.3e}{:+3e}i\n", k + 1, PetscRealPart(z[k]) * gamma,
               PetscImaginaryPart(z[k]) * gamma);
    opInv[k].SetOperators(z[k] * gamma, w[k], *opK, *opM);
#if 0
    int l = 0;
    for (; l < k; l++)
    {
      constexpr double atol = 1.0e-9;
      if (PetscAbsReal(PetscAbsScalar(z[k]) - PetscAbsScalar(z[l])) < atol)
      {
        // Reuse preconditioner assembled for contour point with same real magnitude.
        opInv[k].SetOperators(z[k] * gamma, w[k], *opK, *opM, opInv[l].pc);
        break;
      }
    }
    if (l == k)
    {
      opInv[k].SetOperators(z[k] * gamma, w[k], *opK, *opM);
    }
#endif
  }
  Mpi::Print(comm, "\n");
  delete[] z;
  delete[] w;

  // Solve the quadratic eigenvalue problem.
  int nconv = SolveInternal(rg);

  // Cleanup.
  PalacePetscCall(RGDestroy(&rg));
  delete AQ;
  delete BQ;
  delete XQ;
  delete XQ0;

  return nconv;
}

void FeastEPSSolver::SolveProjectedProblem(const petsc::PetscDenseMatrix &Q_,
                                           petsc::PetscDenseMatrix &R_,
                                           petsc::PetscDenseMatrix &X_, PetscScalar *eig_)
{
  // Form mQ x mQ projected matrices.
  // AQ->Scale(0.0);
  // BQ->Scale(0.0);
  BVMatProjectInternal(Q_, *opK, *AQ, R_, delta);
  BVMatProjectInternal(Q_, *opM, *BQ, R_, delta * gamma);

  // Solve projected EVP using LAPACK wrapper.
  PetscBLASInt info, n, lwork, lrwork;
  PetscScalar *work, *alpha, *beta;
  PetscReal *rwork;
  PetscBLASIntCast(mQ, &n);
  lwork = 2 * n;
  lrwork = 8 * n;
  work = new PetscScalar[lwork];
  rwork = new PetscReal[lrwork];
  alpha = new PetscScalar[n];
  beta = new PetscScalar[n];

  PetscScalar *pAQ = AQ->GetArray();
  PetscScalar *pBQ = BQ->GetArray();
  PetscScalar *pXQ = XQ->GetArray();
  LAPACKggev_("N", "V", &n, pAQ, &n, pBQ, &n, alpha, beta, nullptr, &n, pXQ, &n, work,
              &lwork, rwork, &info);
  AQ->RestoreArray(pAQ);
  BQ->RestoreArray(pBQ);
  XQ->RestoreArray(pXQ);

  // Sort eigenpairs by distance to center.
  for (PetscBLASInt i = 0; i < n; i++)
  {
    alpha[i] /= beta[i];
  }

  // Debug
  // Mpi::Print(Q_.GetComm(), "Before sort, eigenvalues:\n");
  // for (PetscBLASInt i = 0; i < n; i++)
  // {
  //   Mpi::Print(Q_.GetComm(), " {:+e}{:+e}i\n",
  //               PetscRealPart(alpha[i]*gamma),
  //               PetscImaginaryPart(alpha[i]*gamma));
  // }

  PetscInt *sort = SortEigenvalues(alpha, n);
  for (PetscInt i = 0; i < m0; i++)
  {
    eig_[i] = alpha[sort[i]];
    const petsc::PetscParVector xq = XQ->GetColumnRead(sort[i]);
    petsc::PetscParVector xq0 = XQ0->GetColumn(i);
    xq0.Copy(xq);
    XQ->RestoreColumnRead(sort[i], xq);
    XQ0->RestoreColumn(i, xq0);
  }

  // Cleanup.
  delete[] sort;
  delete[] work;
  delete[] rwork;
  delete[] alpha;
  delete[] beta;

  // Reconstruct the first m0 high-dimensional eigenvectors.
  const PetscScalar *pQ = Q_.GetArrayRead();
  PetscScalar *pX = X_.GetArray();
  petsc::PetscDenseMatrix locQ(X_.GetNumRows(), mQ, const_cast<PetscScalar *>(pQ));
  petsc::PetscDenseMatrix locX(X_.GetNumRows(), m0, pX);
  locQ.MatMult(*XQ0, locX);
  Q_.RestoreArrayRead(pQ);
  X_.RestoreArray(pX);
}

void FeastEPSSolver::GetResidual(PetscScalar eig_, const petsc::PetscParVector &x_,
                                 petsc::PetscParVector &r_) const
{
  // r = (K - λ M) x for eigenvalue λ.
  opM->Mult(x_, r_);
  r_.Scale(-eig_);
  opK->MultAdd(x_, r_);
}

PetscReal FeastEPSSolver::GetBackwardScaling(PetscScalar eig_) const
{
  // Make sure not to use norms from scaling as this can be confusing if they are different.
  // Note that SLEPc uses ||.||∞, not Frobenius.
  if (normK <= 0.0)
  {
    normK = opK->NormInf();
  }
  if (normM <= 0.0)
  {
    normM = opM->NormInf();
  }
  return normK + PetscAbsScalar(eig_) * normM;
}

// PEP specific methods

FeastPEPSolver::FeastPEPSolver(MPI_Comm comm, const IoData &iodata, SpaceOperator &spaceop,
                               int np, int print_lvl)
  : FeastEigenSolver(comm, iodata, spaceop, np, print_lvl)
{
  opK = opC = opM = nullptr;
  normK = normC = normM = 0.0;
  AQ = BQ = AQ0 = XQ = XQ0 = nullptr;
}

void FeastPEPSolver::SetOperators(const petsc::PetscParMatrix &K,
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
    MFEM_VERIFY(normK > 0.0 && normC > 0.0 && normM > 0.0,
                "Invalid matrix norms for PEP scaling!");
    gamma = std::sqrt(normK / normM);
    delta = 2.0 / (normK + gamma * normC);
  }
}

int FeastPEPSolver::Solve()
{
  // Check inputs.
  CheckParameters();
  MFEM_VERIFY(opK && opC && opM, "Operators are not set for FeastPEPSolver!");

  // Allocate storage for eigenvectors.
  MPI_Comm comm = opK->GetComm();
  if (!X)
  {
    X = new petsc::PetscDenseMatrix(comm, opK->GetNumRows(), PETSC_DECIDE, PETSC_DECIDE, m0,
                                    nullptr);
  }
  if (!r0)
  {
    r0 = new petsc::PetscParVector(*opK);
  }

  // Allocate sequential matrices for the projected linearized generalized eigenvalue
  // problems at each iteration, and associated eigenvectors.
  AQ = new petsc::PetscDenseMatrix(2 * mQ, 2 * mQ, nullptr);
  BQ = new petsc::PetscDenseMatrix(2 * mQ, 2 * mQ, nullptr);
  AQ0 = new petsc::PetscDenseMatrix(mQ, mQ, nullptr);
  AQ0->SetSymmetric(opK->GetSymmetric() && opC->GetSymmetric() && opM->GetSymmetric());
  XQ = new petsc::PetscDenseMatrix(2 * mQ, 2 * mQ, nullptr);
  XQ0 = new petsc::PetscDenseMatrix(mQ, m0, nullptr);

  // Create region object for integration contour and configure the linear solvers at each
  // integration point. The linear solves use the unscaled space.
  PetscScalar *z, *w;
  RG rg = ConfigureRG(z, w);
  Mpi::Print(comm, "Quadrature points for FEAST contour\n");
  for (int k = 0; k < static_cast<int>(opInv.size()); k++)
  {
    Mpi::Print(comm, " {:d}: zₖ = {:+.3e}{:+.3e}i\n", k + 1, PetscRealPart(z[k]) * gamma,
               PetscImaginaryPart(z[k]) * gamma);
    opInv[k].SetOperators(z[k] * gamma, w[k], *opK, *opC, *opM);
#if 0
    int l = 0;
    for (; l < k; l++)
    {
      constexpr double atol = 1.0e-9;
      if (PetscAbsReal(PetscAbsScalar(z[k]) - PetscAbsScalar(z[l])) < atol)
      {
        // Reuse preconditioner assembled for contour point with same real magnitude.
        opInv[k].SetOperators(z[k] * gamma, w[k], *opK, *opC, *opM, opInv[l].pc);
        break;
      }
    }
    if (l == k)
    {
      opInv[k].SetOperators(z[k] * gamma, w[k], *opK, *opC, *opM);
    }
#endif
  }
  Mpi::Print(comm, "\n");
  delete[] z;
  delete[] w;

  // Solve the quadratic eigenvalue problem.
  int nconv = SolveInternal(rg);

  // Cleanup.
  PalacePetscCall(RGDestroy(&rg));
  delete AQ;
  delete BQ;
  delete AQ0;
  delete XQ;
  delete XQ0;

  return nconv;
}

void FeastPEPSolver::SolveProjectedProblem(const petsc::PetscDenseMatrix &Q_,
                                           petsc::PetscDenseMatrix &R_,
                                           petsc::PetscDenseMatrix &X_, PetscScalar *eig_)
{
  // Form mQ x mQ projected matrices and construct the canonincal linearization:
  //               L₀ = [  0   I ]    L₁ = [ I  0 ]
  //                    [ -K  -C ] ,       [ 0  M ] .
  AQ->Scale(0.0);
  BQ->Scale(0.0);
  PetscScalar *pAQ = AQ->GetArray();
  PetscScalar *pBQ = BQ->GetArray();
  for (PetscInt i = 0; i < mQ; i++)
  {
    pAQ[i + 2 * mQ * (i + mQ)] = 1.0;
  }
  {
    // AQ0->Scale(0.0);
    BVMatProjectInternal(Q_, *opK, *AQ0, R_, delta);

    const PetscScalar *pAQ0 = AQ0->GetArrayRead();
    for (PetscInt j = 0; j < mQ; j++)
    {
      for (PetscInt i = 0; i < mQ; i++)
      {
        pAQ[i + mQ + 2 * mQ * j] = -pAQ0[i + mQ * j];
      }
    }
    AQ0->RestoreArrayRead(pAQ0);
  }
  {
    // AQ0->Scale(0.0);
    BVMatProjectInternal(Q_, *opC, *AQ0, R_, delta * gamma);

    const PetscScalar *pAQ0 = AQ0->GetArrayRead();
    for (PetscInt j = 0; j < mQ; j++)
    {
      for (PetscInt i = 0; i < mQ; i++)
      {
        pAQ[i + mQ + 2 * mQ * (j + mQ)] = -pAQ0[i + mQ * j];
      }
    }
    AQ0->RestoreArrayRead(pAQ0);
  }
  for (PetscInt i = 0; i < mQ; i++)
  {
    pBQ[i + 2 * mQ * i] = 1.0;
  }
  {
    // AQ0->Scale(0.0);
    BVMatProjectInternal(Q_, *opM, *AQ0, R_, delta * gamma * gamma);

    const PetscScalar *pAQ0 = AQ0->GetArrayRead();
    for (PetscInt j = 0; j < mQ; j++)
    {
      PalacePetscCall(PetscArraycpy(pBQ + mQ + 2 * mQ * (j + mQ), pAQ0 + mQ * j, mQ));
    }
    AQ0->RestoreArrayRead(pAQ0);
  }

  // Solve projected EVP using LAPACK wrapper.
  PetscBLASInt info, n, lwork, lrwork;
  PetscScalar *work, *alpha, *beta;
  PetscReal *rwork;
  PetscBLASIntCast(2 * mQ, &n);
  lwork = 2 * n;
  lrwork = 8 * n;
  work = new PetscScalar[lwork];
  rwork = new PetscReal[lrwork];
  alpha = new PetscScalar[n];
  beta = new PetscScalar[n];

  PetscScalar *pXQ = XQ->GetArray();
  LAPACKggev_("N", "V", &n, pAQ, &n, pBQ, &n, alpha, beta, nullptr, &n, pXQ, &n, work,
              &lwork, rwork, &info);
  AQ->RestoreArray(pAQ);
  BQ->RestoreArray(pBQ);
  XQ->RestoreArray(pXQ);

  // Sort eigenpairs by distance to center. From the linearization, we extract the
  // eigenvectors from the top block and normalize later on.
  for (PetscBLASInt i = 0; i < n; i++)
  {
    alpha[i] /= beta[i];
  }

  // Debug
  // Mpi::Print(Q_.GetComm(), "Before sort, eigenvalues:\n");
  // for (PetscBLASInt i = 0; i < n; i++)
  // {
  //   Mpi::Print(Q_.GetComm(), " {:+e}{:+e}i\n",
  //               PetscRealPart(alpha[i]*gamma),
  //               PetscImaginaryPart(alpha[i]*gamma));
  // }

  PetscInt *sort = SortEigenvalues(alpha, n);
  for (PetscInt i = 0; i < m0; i++)
  {
    eig_[i] = alpha[sort[i]];
    const PetscScalar *pXQ = XQ->GetArrayRead();
    PetscScalar *pXQ0 = XQ0->GetArray();
    PalacePetscCall(PetscArraycpy(pXQ0 + mQ * i, pXQ + 2 * mQ * sort[i], mQ));
    XQ->RestoreArrayRead(pXQ);
    XQ0->RestoreArray(pXQ0);
  }

  // Cleanup.
  delete[] sort;
  delete[] work;
  delete[] rwork;
  delete[] alpha;
  delete[] beta;

  // Reconstruct the first m0 high-dimensional eigenvectors.
  const PetscScalar *pQ = Q_.GetArrayRead();
  PetscScalar *pX = X_.GetArray();
  petsc::PetscDenseMatrix locQ(X_.GetNumRows(), mQ, const_cast<PetscScalar *>(pQ));
  petsc::PetscDenseMatrix locX(X_.GetNumRows(), m0, pX);
  locQ.MatMult(*XQ0, locX);
  Q_.RestoreArrayRead(pQ);
  X_.RestoreArray(pX);
}

void FeastPEPSolver::GetResidual(PetscScalar eig_, const petsc::PetscParVector &x_,
                                 petsc::PetscParVector &r_) const
{
  // r = P(λ) x = (K + λ C + λ² M) x for eigenvalue λ.
  opM->Mult(x_, r_);
  r_.Scale(eig_);
  opC->MultAdd(x_, r_);
  r_.Scale(eig_);
  opK->MultAdd(x_, r_);
}

PetscReal FeastPEPSolver::GetBackwardScaling(PetscScalar eig_) const
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
  PetscReal t = PetscAbsScalar(eig_);
  return normK + t * normC + t * t * normM;
}

}  // namespace palace::feast

PetscErrorCode __mat_apply_FEAST_EPS(Mat A, Vec x, Vec y)
{
  // Apply the operator: K - zₖ M .
  palace::feast::internal::FeastLinearSolver *feast;
  palace::petsc::PetscParVector xx(x, true), yy(y, true);
  PetscFunctionBeginUser;

  PetscCall(MatShellGetContext(A, (void **)&feast));
  MFEM_VERIFY(feast, "Invalid PETSc shell matrix context for FEAST!");
  {
    feast->opM->Mult(xx, yy);
    yy.Scale(-feast->zk);
    feast->opK->MultAdd(xx, yy);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode __mat_apply_FEAST_PEP(Mat A, Vec x, Vec y)
{
  // Apply the operator: K + zₖ C + zₖ² M .
  palace::feast::internal::FeastLinearSolver *feast;
  palace::petsc::PetscParVector xx(x, true), yy(y, true);
  PetscFunctionBeginUser;

  PetscCall(MatShellGetContext(A, (void **)&feast));
  MFEM_VERIFY(feast, "Invalid PETSc shell matrix context for FEAST!");
  {
    feast->opM->Mult(xx, yy);
    yy.Scale(feast->zk);
    feast->opC->MultAdd(xx, yy);
    yy.Scale(feast->zk);
    feast->opK->MultAdd(xx, yy);
  }
  PetscFunctionReturn(0);
}

#endif
