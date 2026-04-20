// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "nleps.hpp"

#include <algorithm>
#include <Eigen/Dense>
#include <mfem.hpp>
#include "linalg/divfree.hpp"
#include "linalg/rap.hpp"
#include "utils/communication.hpp"

namespace palace
{

using namespace std::complex_literals;
// Base class methods.

NonLinearEigenvalueSolver::NonLinearEigenvalueSolver(MPI_Comm comm, int print)
  : comm(comm), print(print)
{
  // Initialization.
  nleps_it = 0;

  gamma = delta = 1.0;
  sigma = 0.0;

  opInv = nullptr;
  opProj = nullptr;
  opB = nullptr;
}

void NonLinearEigenvalueSolver::SetLinearSolver(ComplexKspSolver &ksp)
{
  opInv = &ksp;
}

void NonLinearEigenvalueSolver::SetDivFreeProjector(
    const DivFreeSolver<ComplexVector> &divfree)
{
  opProj = &divfree;
}

void NonLinearEigenvalueSolver::SetBMat(const Operator &B)
{
  opB = &B;
}

void NonLinearEigenvalueSolver::SetNumModes(int num_eig, int num_vec)
{
  if (nev > 0 && num_eig != nev)
  {
    res.reset();
    xscale.reset();
    perm.reset();
  }
  nev = num_eig;
}

void NonLinearEigenvalueSolver::SetTol(double tol)
{
  rtol = tol;
}

void NonLinearEigenvalueSolver::SetMaxIter(int max_it)
{
  nleps_it = max_it;
}

void NonLinearEigenvalueSolver::SetWhichEigenpairs(EigenvalueSolver::WhichType type)
{
  which_type = type;
}

void NonLinearEigenvalueSolver::SetShiftInvert(std::complex<double> s, bool precond)
{
  MFEM_VERIFY(!precond, "Nonlinear eigenvalue solver does not support preconditioned "
                        "spectral transformation option!");
  sigma = s;
  sinvert = true;
}

void NonLinearEigenvalueSolver::SetInitialSpace(const ComplexVector &v)
{
  MFEM_ABORT("SetInitialSpace not defined for base class NonLinearEigenvalueSolver!");
}

std::complex<double> NonLinearEigenvalueSolver::GetEigenvalue(int i) const
{
  MFEM_VERIFY(i >= 0 && i < nev,
              "Out of range eigenpair requested (i = " << i << ", nev = " << nev << ")!");
  const int &j = perm.get()[i];
  return eigenvalues[j];
}

void NonLinearEigenvalueSolver::GetEigenvector(int i, ComplexVector &x) const
{
  MFEM_VERIFY(i >= 0 && i < nev,
              "Out of range eigenpair requested (i = " << i << ", nev = " << nev << ")!");
  MFEM_VERIFY(x.Size() == n, "Invalid size mismatch for provided eigenvector!");
  const int &j = perm.get()[i];
  x = eigenvectors[j];
  if (xscale.get()[j] > 0.0)
  {
    x *= xscale.get()[j];
  }
}

double NonLinearEigenvalueSolver::GetEigenvectorNorm(const ComplexVector &x,
                                                     ComplexVector &Bx) const
{
  if (opB)
  {
    return linalg::Norml2(comm, x, *opB, Bx);
  }
  else
  {
    return linalg::Norml2(comm, x);
  }
}

double NonLinearEigenvalueSolver::GetError(int i, EigenvalueSolver::ErrorType type) const
{
  MFEM_VERIFY(i >= 0 && i < nev,
              "Out of range eigenpair requested (i = " << i << ", nev = " << nev << ")!");
  const int &j = perm.get()[i];
  switch (type)
  {
    case ErrorType::ABSOLUTE:
      return res.get()[j];
    case ErrorType::RELATIVE:
      return res.get()[j] / std::abs(eigenvalues[j]);
    case ErrorType::BACKWARD:
      return res.get()[j] / GetBackwardScaling(eigenvalues[j]);
  }
  return 0.0;
}

void NonLinearEigenvalueSolver::RescaleEigenvectors(int num_eig)
{
  res = std::make_unique<double[]>(num_eig);
  xscale = std::make_unique<double[]>(num_eig);
  for (int i = 0; i < num_eig; i++)
  {
    x1 = eigenvectors[i];
    xscale.get()[i] = 1.0 / GetEigenvectorNorm(x1, y1);
    res.get()[i] = GetResidualNorm(eigenvalues[i], x1, y1) / linalg::Norml2(comm, x1);
  }
}

// Quasi-Newton specific methods.
QuasiNewtonSolver::QuasiNewtonSolver(MPI_Comm comm,
                                     std::unique_ptr<EigenvalueSolver> linear_eigensolver,
                                     int num_conv, int print, bool refine)
  : NonLinearEigenvalueSolver(comm, print),
    linear_eigensolver_(std::move(linear_eigensolver)), nev_linear(num_conv),
    refine_nonlinear(refine)
{
  opK = opC = opM = nullptr;
  normK = normC = normM = 0.0;
}

// Set the update frequency of the preconditioner.
void QuasiNewtonSolver::SetPreconditionerLag(int preconditioner_update_freq,
                                             double preconditioner_update_tol)
{
  preconditioner_lag = preconditioner_update_freq;
  preconditioner_tol = preconditioner_update_tol;
}

// Set the maximum number of restarts with the same initial guess.
void QuasiNewtonSolver::SetMaxRestart(int max_num_restart)
{
  max_restart = max_num_restart;
}

void QuasiNewtonSolver::SetExtraSystemMatrix(
    std::function<std::unique_ptr<ComplexOperator>(double)> A2)
{
  funcA2 = A2;
}

void QuasiNewtonSolver::SetPreconditionerUpdate(
    std::function<std::unique_ptr<ComplexOperator>(
        std::complex<double>, std::complex<double>, std::complex<double>, double)>
        P)
{
  funcP = P;
}

void QuasiNewtonSolver::SetOperators(const ComplexOperator &K, const ComplexOperator &M,
                                     EigenvalueSolver::ScaleType type)
{
  MFEM_VERIFY(!opK || K.Height() == n, "Invalid modification of eigenvalue problem size!");
  bool first = (opK == nullptr);
  opK = &K;
  opM = &M;

  if (first && type != ScaleType::NONE)
  {
    normK = linalg::SpectralNorm(comm, *opK, opK->IsReal());
    normM = linalg::SpectralNorm(comm, *opM, opM->IsReal());
    MFEM_VERIFY(normK >= 0.0 && normM >= 0.0,
                "Invalid matrix norms for Quasi-Newton scaling!");
    if (normK > 0 && normM > 0.0)
    {
      gamma = std::sqrt(normK / normM);
      delta = 2.0 / (normK);
    }
  }

  n = opK->Height();

  // Set up workspace.
  x1.SetSize(opK->Height());
  y1.SetSize(opK->Height());
  x1.UseDevice(true);
  y1.UseDevice(true);
}

void QuasiNewtonSolver::SetOperators(const ComplexOperator &K, const ComplexOperator &C,
                                     const ComplexOperator &M,
                                     EigenvalueSolver::ScaleType type)
{
  MFEM_VERIFY(!opK || K.Height() == n, "Invalid modification of eigenvalue problem size!");
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
                "Invalid matrix norms for Quasi-Newton scaling!");
    if (normK > 0 && normC >= 0.0 && normM > 0.0)
    {
      gamma = std::sqrt(normK / normM);
      delta = 2.0 / (normK + gamma * normC);
    }
  }

  n = opK->Height();

  // Set up workspace.
  x1.SetSize(opK->Height());
  y1.SetSize(opK->Height());
  x1.UseDevice(true);
  y1.UseDevice(true);
}

void QuasiNewtonSolver::SetInitialGuess()
{

  MFEM_VERIFY(n > 0, "Must call SetOperators before using SetInitialguess for nonlinear "
                     "eigenvalue solver!");
  MFEM_VERIFY(nev > 0, "Must call SetNumModes before using SetInitialguess for nonlinear "
                       "eigenvalue solver!");

  // Get eigenmodes initial guesses from linear eigensolver
  eigenvalues.resize(nev_linear);
  eigenvectors.resize(nev_linear);
  for (int i = 0; i < nev_linear; i++)
  {
    eigenvalues[i] = linear_eigensolver_->GetEigenvalue(i);
    linear_eigensolver_->GetEigenvector(i, x1);
    eigenvectors[i] = x1;
  }

  // Compute errors.
  RescaleEigenvectors(nev_linear);

  // Initialize eigenpairs ordering.
  perm = std::make_unique<int[]>(nev_linear);
  std::iota(perm.get(), perm.get() + nev_linear, 0);

  // Early return if nonlinear Newton won't be used.
  if (!refine_nonlinear)
    return;

  // If the number of initial guesses is greater than the number of requested modes
  // de-prioritize the initial guesses that have larger errors.
  std::vector<std::size_t> indices(nev_linear);
  std::iota(indices.begin(), indices.end(), 0);
  if (nev_linear > nev)
  {
    double min_error = res.get()[0];
    for (int i = 0; i < nev_linear; i++)
    {
      min_error = std::min(min_error, res.get()[i]);
    }
    const double threshold = 100.0 * min_error;
    std::sort(indices.begin(), indices.end(),
              [&](const auto i, const auto j)
              {
                if (res.get()[i] < threshold && res.get()[j] > threshold)
                {
                  return true;
                }
                else if (res.get()[i] > threshold && res.get()[j] < threshold)
                {
                  return false;
                }
                else
                {
                  return eigenvalues[i].imag() < eigenvalues[j].imag();
                }
              });
  }
  for (int i = 0; i < nev_linear; i++)
  {
    eigenvalues[i] = linear_eigensolver_->GetEigenvalue(indices[i]);
    linear_eigensolver_->GetEigenvector(indices[i], x1);
    linalg::NormalizePhase(comm, x1);
    eigenvectors[i] = x1;
  }

  // Get ordering of the eigenpairs.
  std::sort(perm.get(), perm.get() + nev_linear, [&eig = this->eigenvalues](auto l, auto r)
            { return eig[l].imag() < eig[r].imag(); });
}

namespace
{
// Multiply an (n x k) matrix (vector of size k of ComplexVectors of size n) by a vector of
// size k, returning a ComplexVector of size n.
ComplexVector MatVecMult(const std::vector<ComplexVector> &X, const Eigen::VectorXcd &y)
{
  // Cast to avoid compiler warnings about types.
  MFEM_ASSERT(static_cast<Eigen::Index>(X.size()) == y.size(),
              "Mismatch in dimension of input vectors!");
  const int k = X.size();
  const int n = X[0].Size();
  const bool use_dev = X[0].UseDevice();
  ComplexVector z;
  z.SetSize(n);
  z.UseDevice(use_dev);
  z = 0.0;
  for (int j = 0; j < k; j++)
  {
    linalg::AXPBYPCZ(y(j).real(), X[j].Real(), -y(j).imag(), X[j].Imag(), 1.0, z.Real());
    linalg::AXPBYPCZ(y(j).imag(), X[j].Real(), y(j).real(), X[j].Imag(), 1.0, z.Imag());
  }
  return z;
}

}  // namespace

int QuasiNewtonSolver::Solve()
{
  // Quasi-Newton method for nonlinear eigenvalue problems.
  // Reference: Jarlebring, Koskela, Mele, Disguised and new quasi-Newton methods for
  //            nonlinear eigenvalue problems, Numerical Algorithms (2018).
  // Using the deflation scheme used by SLEPc's NEP solver with minimality index set to 1.
  // Reference: Effenberger, Robust successive computation of eigenpairs for nonlinear
  //            eigenvalue problems, SIAM J. Matrix Anal. Appl. (2013).
  // The deflation scheme solves an extended problem of size n + k, where n is the original
  // problem size and k is the number of converged eigenpairs. The extended operators are
  // never explicitly constructed and two separate vectors of length n and k are used to
  // store the extended solution: [v, v2] where v is a ComplexVector distributed across all
  // processes and v2 is an Eigen::VectorXcd stored redundantly on all processes.

  // Set initial guess from linear eigensolver.
  SetInitialGuess();

  // Return early if not refining the eigenmodes with Newton.
  if (!refine_nonlinear)
  {
    nev = eigenvalues.size();
    return nev;
  }

  // u/u2 hold the residual at the committed iterate; du/du2 the Newton correction.
  // v_trial/v2_trial are scratch for line-search trial normalization (v itself cannot
  // be rolled back once normalized).
  ComplexVector v, u, w, c, w0, z, du, v_trial;
  v.SetSize(n);
  u.SetSize(n);
  w.SetSize(n);
  c.SetSize(n);
  w0.SetSize(n);
  z.SetSize(n);
  du.SetSize(n);
  v_trial.SetSize(n);
  v.UseDevice(true);
  u.UseDevice(true);
  w.UseDevice(true);
  c.UseDevice(true);
  w0.UseDevice(true);
  z.UseDevice(true);
  du.UseDevice(true);
  v_trial.UseDevice(true);

  // Eigen Matrix/Vectors for extended operator of size k.
  Eigen::MatrixXcd H;
  Eigen::VectorXcd u2, z2, c2, w2, v2, du2, v2_trial;

  // Storage for eigenpairs.
  std::vector<ComplexVector> X;
  std::vector<std::complex<double>> eigs;
  X.reserve(nev);

  // Set defaults.
  if (nleps_it <= 0)
  {
    nleps_it = 100;
  }

  // Delta used in to compute divided difference Jacobian.
  const auto delta = std::sqrt(std::numeric_limits<double>::epsilon());

  // Set a seed and distribution for random Eigen vectors to ensure the same values on all
  // ranks.
  unsigned int seed = 111;
  std::mt19937 gen(seed);
  std::uniform_real_distribution<> dist_real(-1.0, 1.0);

  const int num_init_guess = eigenvalues.size();
  std::uniform_int_distribution<int> dist_int(0, num_init_guess - 1);

  // Save the user-specified solver tolerances.
  const double ksp_rel_tol = opInv->GetRelTol();
  const double inexact_tol = std::sqrt(rtol);

  int k = 0, restart = 0, guess_idx = 0;
  while (k < nev)
  {
    // If > max_restart with the same initial guess, skip to next initial guess.
    // If we tried all initial guesses and the random guess, end search even if k < nev.
    if (restart > max_restart)
    {
      if (guess_idx < num_init_guess)
      {
        guess_idx++;
        restart = 0;
      }
      else
      {
        break;
      }
    }

    // Set the eigenpair estimate to the initial guess.
    std::complex<double> eig, eig_opInv;
    if (guess_idx < num_init_guess)
    {
      eig = eigenvalues[guess_idx];
      v = eigenvectors[guess_idx];
    }
    else
    {
      eig = sigma;
      if (num_init_guess < 3)
      {
        // Set purely random vector.
        linalg::SetRandom(GetComm(), v);
        linalg::SetSubVector(
            v, *dynamic_cast<const ComplexParOperator *>(opK)->GetEssentialTrueDofs(), 0.0);
      }
      else
      {
        // Set random vector as the average of two distinct randomly-chosen initial guesses.
        int i1 = dist_int(gen);
        int i2;
        do
        {
          i2 = dist_int(gen);
        } while (i2 == i1);
        v.AXPBYPCZ(0.5, eigenvectors[i1], 0.5, eigenvectors[i2], 0.0);
      }
    }
    eig_opInv = eig;  // eigenvalue estimate used in the (lagged) preconditioner

    // Set the "random" c vector and the deflation component of the eigenpair initial guess.
    linalg::SetRandom(GetComm(), c, seed);  // set seed for deterministic behavior
    c2.conservativeResize(k);
    v2.conservativeResize(k);
    for (int i = 0; i < k; i++)
    {
      c2(i) = std::complex<double>(dist_real(gen), dist_real(gen));
      v2(i) = std::complex<double>(dist_real(gen), dist_real(gen));
    }

    // Normalize random c vector.
    double norm_c = std::sqrt(std::abs(linalg::Dot(GetComm(), c, c)) + c2.squaredNorm());
    c *= 1.0 / norm_c;
    c2 *= 1.0 / norm_c;

    // Normalize eigenvector estimate.
    double norm_v = std::sqrt(std::abs(linalg::Dot(GetComm(), v, v)) + v2.squaredNorm());
    v *= 1.0 / norm_v;
    v2 *= 1.0 / norm_v;

    // Set the linear solver operators.
    opA2 = (*funcA2)(std::abs(eig.imag()));
    opA = BuildParSumOperator({1.0 + 0.0i, eig, eig * eig, 1.0 + 0.0i},
                              {opK, opC, opM, opA2.get()}, true);
    opP = (*funcP)(1.0 + 0.0i, eig, eig * eig, eig.imag());
    opInv->SetOperators(*opA, *opP);
    opInv->SetAbsTol(1.0e-12);

    // Linear solve with the extended operator of the deflated problem.
    auto deflated_solve = [&](const ComplexVector &b1, const Eigen::VectorXcd &b2,
                              ComplexVector &x1, Eigen::VectorXcd &x2)
    {
      // Solve the block linear system
      // |T(σ) U(σ)| |x1| = |b1|
      // |A(σ) B(σ)| |x2|   |b2|
      // x1 = T^-1 b1
      // x2 = SS^-1 (b2 - A x1) where SS = (B - A T^-1 U) = - X^* X S^-1
      // x1 = x1 - X S x2
      opInv->Mult(b1, x1);
      if (k == 0)  // no deflation
      {
        return;
      }
      x2.conservativeResize(k);
      for (int j = 0; j < k; j++)
      {
        x2(j) = b2(j) - linalg::Dot(GetComm(), x1, X[j]);
      }
      Eigen::MatrixXcd SS(k, k);
      for (int i = 0; i < k; i++)
      {
        for (int j = 0; j < k; j++)
        {
          SS(i, j) = linalg::Dot(GetComm(), X[i], X[j]);
        }
      }
      const Eigen::MatrixXcd S = eig_opInv * Eigen::MatrixXcd::Identity(k, k) - H;
      SS = -S.fullPivLu().solve(SS);
      x2 = SS.fullPivLu().solve(x2);
      const ComplexVector XSx2 = MatVecMult(X, S.fullPivLu().solve(x2));
      linalg::AXPY(-1.0, XSx2, x1);
    };

    // Compute w0 = T^-1 c and normalize it. The w0 vector is only used as a
    // projection direction for the eigenvalue correction, so moderate accuracy suffices.
    opInv->SetRelTol(std::max(ksp_rel_tol, inexact_tol));
    deflated_solve(c, c2, w0, w2);
    double norm_w0 = std::sqrt(std::abs(linalg::Dot(GetComm(), w0, w0)) + w2.squaredNorm());
    w0 *= 1.0 / norm_w0;
    w2 *= 1.0 / norm_w0;

    // Evaluate the deflated residual r = T(lam) vv + T(lam) X (lam I - H)^-1 vv2, with
    // rr2 = X^* vv. A2_out returns the built A2 operator so the caller can hold onto it
    // and skip re-assembling at the same lam.
    auto compute_residual = [this, &k, &H,
                             &X](std::complex<double> lam, const ComplexVector &vv,
                                 const Eigen::VectorXcd &vv2, ComplexVector &rr,
                                 Eigen::VectorXcd &rr2,
                                 std::unique_ptr<ComplexOperator> &A2_out) -> double
    {
      A2_out = (*funcA2)(std::abs(lam.imag()));
      auto A = BuildParSumOperator({1.0 + 0.0i, lam, lam * lam, 1.0 + 0.0i},
                                   {opK, opC, opM, A2_out.get()}, true);
      A->Mult(vv, rr);
      if (k > 0)
      {
        const Eigen::MatrixXcd S = lam * Eigen::MatrixXcd::Identity(k, k) - H;
        const ComplexVector XSvv2 = MatVecMult(X, S.fullPivLu().solve(vv2));
        A->AddMult(XSvv2, rr, 1.0);
        rr2.conservativeResize(k);
        for (int j = 0; j < k; j++)
        {
          rr2(j) = linalg::Dot(GetComm(), vv, X[j]);
        }
      }
      else
      {
        rr2.resize(0);
      }
      return std::sqrt(std::abs(linalg::Dot(GetComm(), rr, rr)) + rr2.squaredNorm());
    };

    // Standard Armijo constants. At the last backtrack we commit the smallest alpha
    // anyway; the outer divergence check restarts the guess if the step was truly bad.
    constexpr double armijo_c = 1.0e-4;
    constexpr double backtrack_factor = 0.5;
    constexpr int max_backtrack = 10;

    // A2 operator at the committed eig, carried across iterations so the Jacobian can
    // reuse the assembly from the accepted trial.
    std::unique_ptr<ComplexOperator> A2n;
    double res = compute_residual(eig, v, v2, u, u2, A2n);

    int it = 0, diverged_it = 0;
    while (it < nleps_it)
    {
      if (print > 0)
      {
        Mpi::Print(GetComm(),
                   "{:d} NLEPS (nconv={:d}, restart={:d}) residual norm {:.6e}\n", it, k,
                   restart, res);
      }

      // End if residual below tolerance and eigenvalue above the target.
      if (res < rtol)
      {
        if (print > 0)
        {
          Mpi::Print(GetComm(),
                     "Eigenvalue {:d}, Quasi-Newton converged in {:d} iterations "
                     "({:.3e}{:+.3e}i).\n",
                     k, it, eig.real(), eig.imag());
        }
        // Update the invariant pair with normalization.
        const auto scale = linalg::Norml2(GetComm(), v);
        v *= 1.0 / scale;
        eigs.resize(k + 1);
        eigs[k] = eig;
        X.resize(k + 1);
        X[k] = v;
        H.conservativeResizeLike(Eigen::MatrixXd::Zero(k + 1, k + 1));
        H.col(k).head(k) = v2 / scale;
        H(k, k) = eig;
        k++;
        // If the eigenvalue is inside the desired range, increment initial guess index
        // Otherwise, use the same initial guess again and increment number of desired
        // eigenvalues.
        if (eig.imag() > sigma.imag())
        {
          guess_idx++;
        }
        else
        {
          nev++;
        }
        restart = 0;  // reset restart counter
        break;
      }
      // Stop if large residual for 10 consecutive iterations.
      diverged_it = (res > 0.9) ? diverged_it + 1 : 0;
      if (diverged_it > 10)
      {
        if (print > 0)
        {
          Mpi::Print(GetComm(),
                     "Eigenvalue {:d}, Quasi-Newton not converging after {:d} iterations, "
                     "restarting.\n",
                     k, it);
        }
        restart++;
        break;
      }

      // Compute w = J * v.
      auto opA2p = (*funcA2)(std::abs(eig.imag()) * (1.0 + delta));
      const std::complex<double> denom =
          std::complex<double>(0.0, delta * std::abs(eig.imag()));
      std::unique_ptr<ComplexOperator> opAJ =
          BuildParSumOperator({1.0 / denom, -1.0 / denom}, {opA2p.get(), A2n.get()}, true);
      auto opJ = BuildParSumOperator({0.0 + 0.0i, 1.0 + 0.0i, 2.0 * eig, 1.0 + 0.0i},
                                     {opK, opC, opM, opAJ.get()}, true);
      opJ->Mult(v, w);
      if (k > 0)  // Deflation
      {
        // w1 = T'(l) v1 + U'(l) v2 = T'(l) v1 + T'(l)XS v2 - T(l)XS^2 v2. Scoping T(l)
        // here lets the line search overwrite A2n freely; with no deflation we skip it.
        auto A = BuildParSumOperator({1.0 + 0.0i, eig, eig * eig, 1.0 + 0.0i},
                                     {opK, opC, opM, A2n.get()}, true);
        const Eigen::MatrixXcd S = eig * Eigen::MatrixXcd::Identity(k, k) - H;
        const Eigen::VectorXcd Sv2 = S.fullPivLu().solve(v2);
        const ComplexVector XSv2 = MatVecMult(X, Sv2);
        const ComplexVector XSSv2 = MatVecMult(X, S.fullPivLu().solve(Sv2));
        opJ->AddMult(XSv2, w, 1.0);
        A->AddMult(XSSv2, w, -1.0);
      }

      // Undamped Newton step for the eigenvalue; the line search damps it.
      const std::complex<double> u2_w0 = std::complex<double>(w2.adjoint() * u2);
      const std::complex<double> delta_eig =
          -(linalg::Dot(GetComm(), u, w0) + u2_w0) / linalg::Dot(GetComm(), w, w0);
      z.AXPBYPCZ(-delta_eig, w, -1.0, u, 0.0);
      z2 = -u2;

      // Inexact Newton: loosen the linear solve tolerance when the outer residual is
      // large, to avoid over-solving when T(σ) is nearly singular.
      opInv->SetRelTol(std::max(ksp_rel_tol, std::min(inexact_tol, res)));
      deflated_solve(z, z2, du, du2);

      // Armijo backtracking on the coupled (eig, v) step. Without damping, Newton
      // overshoots when the linear-eigensolver seed is outside the basin or <w0, w> is
      // near-singular — observed on adapter/hybrid mode 3 at NP >= 32.
      double alpha = 1.0;
      int bt = 0;
      for (; bt < max_backtrack; bt++)
      {
        const std::complex<double> eig_trial = eig + alpha * delta_eig;

        v_trial.AXPBYPCZ(1.0, v, alpha, du, 0.0);
        v2_trial = v2 + alpha * du2;
        const double norm_v_trial = std::sqrt(
            std::abs(linalg::Dot(GetComm(), v_trial, v_trial)) + v2_trial.squaredNorm());
        v_trial *= 1.0 / norm_v_trial;
        v2_trial *= 1.0 / norm_v_trial;

        // In-place writes into u, u2, A2n are safe: u/u2 were consumed into z above,
        // and no outer reference to A2n outlives this loop.
        const double res_trial = compute_residual(eig_trial, v_trial, v2_trial, u, u2, A2n);

        if (res_trial <= (1.0 - armijo_c * alpha) * res || bt == max_backtrack - 1)
        {
          std::swap(v, v_trial);
          std::swap(v2, v2_trial);
          eig = eig_trial;
          res = res_trial;
          break;
        }
        alpha *= backtrack_factor;
      }
      if (print > 0 && bt > 0)
      {
        Mpi::Print(GetComm(),
                   "   NLEPS Armijo backtracks={:d}, alpha={:.3e}, res_new={:.6e}\n", bt,
                   alpha, res);
      }

      // Lagged preconditioner update using the committed (post-line-search) eig.
      // Updating too close to an eigenvalue can cause numerical instability, so lag
      // by preconditioner_lag iterations.
      if (it > 0 && it % preconditioner_lag == 0 && res > preconditioner_tol)
      {
        eig_opInv = eig;
        opA2 = (*funcA2)(std::abs(eig_opInv.imag()));
        opA =
            BuildParSumOperator({1.0 + 0.0i, eig_opInv, eig_opInv * eig_opInv, 1.0 + 0.0i},
                                {opK, opC, opM, opA2.get()}, true);
        opP = (*funcP)(1.0 + 0.0i, eig_opInv, eig_opInv * eig_opInv, eig_opInv.imag());
        opInv->SetOperators(*opA, *opP);
        // Recompute w0 and normalize.
        opInv->SetRelTol(std::max(ksp_rel_tol, inexact_tol));
        deflated_solve(c, c2, w0, w2);
        double norm_w0 =
            std::sqrt(std::abs(linalg::Dot(GetComm(), w0, w0)) + w2.squaredNorm());
        w0 *= 1.0 / norm_w0;
        w2 *= 1.0 / norm_w0;
      }

      it++;
      if (it == nleps_it)
      {
        if (print > 0)
        {
          Mpi::Print(GetComm(),
                     "Eigenvalue {:d}, Quasi-Newton did not converge in {:d} iterations, "
                     "restarting.\n",
                     k, nleps_it);
        }
        restart++;
      }
    }
  }
  nev = k;  // in case some guesses did not converge

  // Eigenpair extraction from the invariant pair (X, H).
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eps;
  eps.compute(H);
  // H eigenvectors are ordered arbitrarily, need to match them to order of X.
  std::vector<int> order(nev), order_eigen(nev), order2(nev);
  std::iota(order.begin(), order.end(), 0);
  std::iota(order_eigen.begin(), order_eigen.end(), 0);
  std::iota(order2.begin(), order2.end(), 0);
  std::sort(order.begin(), order.end(),
            [&](auto l, auto r) { return eigs[l].imag() < eigs[r].imag(); });
  std::sort(order_eigen.begin(), order_eigen.end(),
            [&epseig = eps.eigenvalues()](auto l, auto r)
            { return epseig(l).imag() < epseig(r).imag(); });
  std::sort(order2.begin(), order2.end(),
            [&](auto l, auto r) { return order[l] < order[r]; });

  // Sort Eigen eigenvectors.
  std::vector<Eigen::VectorXcd> Xeig;
  for (int i = 0; i < nev; i++)
  {
    Xeig.push_back(eps.eigenvectors().col(order_eigen[i]));
  }

  // Recover the eigenvectors in the target range.
  eigenvalues.clear();
  eigenvectors.clear();
  for (int i = 0; i < nev; i++)
  {
    if (eigs[i].imag() > sigma.imag())
    {
      ComplexVector eigv = MatVecMult(X, Xeig[order2[i]]);
      eigenvalues.push_back(eigs[i]);
      eigenvectors.push_back(eigv);
    }
  }
  nev = eigenvalues.size();

  // Get ordering of the reported eigenpairs.
  perm = std::make_unique<int[]>(nev);
  std::iota(perm.get(), perm.get() + nev, 0);
  std::sort(perm.get(), perm.get() + nev, [&eig = this->eigenvalues](auto l, auto r)
            { return eig[l].imag() < eig[r].imag(); });

  // Compute the eigenpair residuals for eigenvalue λ.
  RescaleEigenvectors(nev);

  return nev;
}

double QuasiNewtonSolver::GetResidualNorm(std::complex<double> l, const ComplexVector &x,
                                          ComplexVector &r) const
{
  // Compute the i-th eigenpair residual: || P(λ) x ||₂ = || (K + λ C + λ² M + A2(λ)) x ||₂
  // for eigenvalue λ.
  opK->Mult(x, r);
  if (opC)
  {
    opC->AddMult(x, r, l);
  }
  opM->AddMult(x, r, l * l);
  auto A2 = (*funcA2)(std::abs(l.imag()));
  A2->AddMult(x, r, 1.0);
  return linalg::Norml2(comm, r);
}

double QuasiNewtonSolver::GetBackwardScaling(std::complex<double> l) const
{
  // Make sure not to use norms from scaling as this can be confusing if they are different.
  // Note that SLEPc uses ||.||∞, not the 2-norm.
  if (normK <= 0.0)
  {
    normK = linalg::SpectralNorm(comm, *opK, opK->IsReal());
  }
  if (normC <= 0.0 && opC)
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

NewtonInterpolationOperator::NewtonInterpolationOperator(
    std::function<std::unique_ptr<ComplexOperator>(double)> funcA2, int size)
  : funcA2(funcA2)
{
  rhs.SetSize(size);
  rhs.UseDevice(true);
}

// Compute the elementary symmetric polynomial. Used to convert from Newton to monomial
// basis.
template <typename ScalarType>
ScalarType elementarySymmetric(const std::vector<ScalarType> &points, int k, int n)
{
  if (k == 0)
  {
    return 1.0;
  }
  if (k > n || k < 0 || n == 0)
  {
    return 0.0;
  }
  return elementarySymmetric(points, k, n - 1) +
         points[n - 1] * elementarySymmetric(points, k - 1, n - 1);
}

void NewtonInterpolationOperator::Interpolate(const std::complex<double> sigma_min,
                                              const std::complex<double> sigma_max)
{
  // Reset operators and sample points each time Interpolate is called.
  ops.clear();
  ops.resize(num_points);
  points.clear();
  points.resize(num_points);

  // Linearly spaced sample points.
  for (int j = 0; j < num_points; j++)
  {
    points[j] = sigma_min + (double)j * (sigma_max - sigma_min) / (double)(num_points - 1);
  }

  // Build divided difference matrices.
  for (int k = 0; k < num_points; k++)
  {
    for (int j = 0; j < num_points - k; j++)
    {
      if (k == 0)
      {
        auto A2j = (funcA2)(points[j].imag());
        ops[k].push_back(std::move(A2j));
      }
      else
      {
        std::complex<double> denom = points[j + k] - points[j];
        auto A2dd =
            BuildParSumOperator({1.0 / denom, -1.0 / denom},
                                {ops[k - 1][j + 1].get(), ops[k - 1][j].get()}, true);
        ops[k].push_back(std::move(A2dd));
      }
    }
  }

  // Compute monomial coefficients as a function of the Newton polynomial coefficients.
  coeffs.clear();
  coeffs.assign(num_points, std::vector<std::complex<double>>(num_points, 0.0));
  for (int k = 0; k < num_points; k++)
  {
    for (int j = k; j < num_points; j++)
    {
      double sign = ((j - k) % 2 == 0) ? 1 : -1;
      coeffs[k][j] = sign * elementarySymmetric(points, j - k, j);
    }
  }
}

std::unique_ptr<ComplexOperator>
NewtonInterpolationOperator::GetInterpolationOperator(int order) const
{
  MFEM_VERIFY(order >= 0 && order < num_points,
              "Order must be greater than or equal to 0 and smaller than the number of "
              "interpolation points!");
  return BuildParSumOperator({coeffs[order][0], coeffs[order][1], coeffs[order][2]},
                             {ops[0][0].get(), ops[1][0].get(), ops[2][0].get()}, true);
}

void NewtonInterpolationOperator::Mult(int order, const ComplexVector &x,
                                       ComplexVector &y) const
{
  MFEM_VERIFY(order >= 0 && order < num_points,
              "Order must be greater than or equal to 0 and smaller than the number of "
              "interpolation points!");

  y = 0.0;
  for (int j = 0; j < num_points; j++)
  {
    if (coeffs[order][j] != 0.0)
    {
      ops[j][0]->AddMult(x, y, coeffs[order][j]);
    }
  }
}

void NewtonInterpolationOperator::AddMult(int order, const ComplexVector &x,
                                          ComplexVector &y, std::complex<double> a) const
{
  this->Mult(order, x, rhs);
  rhs *= a;
  y += rhs;
}

}  // namespace palace
