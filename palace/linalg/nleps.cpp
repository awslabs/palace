// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "nleps.hpp"

#include <algorithm>
#include <string>
#include <Eigen/Dense>
#include <mfem.hpp>
#include "linalg/divfree.hpp"
#include "utils/communication.hpp"

namespace palace
{

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

void NonLinearEigenvalueSolver::SetOperators(const ComplexOperator &K,
                                             const ComplexOperator &M,
                                             EigenvalueSolver::ScaleType type)
{
  MFEM_ABORT("SetOperators not defined for base class NonLinearEigenvalueSolver!");
}

void NonLinearEigenvalueSolver::SetOperators(const ComplexOperator &K,
                                             const ComplexOperator &C,
                                             const ComplexOperator &M,
                                             EigenvalueSolver::ScaleType type)
{
  MFEM_ABORT("SetOperators not defined for base class NonLinearEigenvalueSolver!");
}

void NonLinearEigenvalueSolver::SetOperators(SpaceOperator &space_op,
                                             const ComplexOperator &K,
                                             const ComplexOperator &C,
                                             const ComplexOperator &M,
                                             EigenvalueSolver::ScaleType type)
{
  MFEM_ABORT("SetOperators not defined for base class NonLinearEigenvalueSolver!");
}

void NonLinearEigenvalueSolver::SetNLInterpolation(const Interpolation &interp)
{
  MFEM_ABORT("SetNLInterpolation not defined for base class NonLinearEigenvalueSolver!");
}

void NonLinearEigenvalueSolver::SetPreconditionerLag(int preconditioner_update_freq,
                                                     double preconditioner_update_tol)
{
  MFEM_ABORT("SetPreconditionerLag not defined for base class NonLinearEigenvalueSolver!");
}

void NonLinearEigenvalueSolver::SetMaxRestart(int max_num_restart)
{
  MFEM_ABORT("SetMaxRestart not defined for base class NonLinearEigenvalueSolver!");
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

void NonLinearEigenvalueSolver::SetInitialGuess(
    const std::vector<std::complex<double>> &init_eig,
    const std::vector<ComplexVector> &init_V)
{
  MFEM_VERIFY(n > 0, "Must call SetOperators before using SetInitialguess for nonlinear "
                     "eigenvalue solver!");
  MFEM_VERIFY(
      init_eig.size() == init_V.size(),
      "SetInitialGuess requires the same number of eigenvalues and eigenvectors guesses!");

  init_eigenvalues.resize(init_eig.size());
  init_eigenvectors.resize(init_eig.size());
  for (int i = 0; i < init_eig.size(); i++)
  {
    init_eigenvalues[i] = init_eig[i];
    init_eigenvectors[i] = init_V[i];
  }
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
QuasiNewtonSolver::QuasiNewtonSolver(MPI_Comm comm, int print)
  : NonLinearEigenvalueSolver(comm, print)
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

void QuasiNewtonSolver::SetOperators(SpaceOperator &space_op_ref, const ComplexOperator &K,
                                     const ComplexOperator &C, const ComplexOperator &M,
                                     EigenvalueSolver::ScaleType type)
{
  MFEM_VERIFY(!opK || K.Height() == n, "Invalid modification of eigenvalue problem size!");
  bool first = (opK == nullptr);
  opK = &K;
  opC = &C;
  opM = &M;
  space_op = &space_op_ref;

  if (first && type != ScaleType::NONE)
  {
    normK = linalg::SpectralNorm(comm, *opK, opK->IsReal());
    if (opC)
    {
      normC = linalg::SpectralNorm(comm, *opC, opC->IsReal());
    }
    normM = linalg::SpectralNorm(comm, *opM, opM->IsReal());
    MFEM_VERIFY(normK >= 0.0 && normC >= 0.0 && normM >= 0.0,
                "Invalid matrix norms for Quasi-Newton scaling!");
    if (normK > 0 && normC > 0.0 && normM > 0.0)
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

namespace
{
// Multiply an (n x k) matrix (vector of size k of ComplexVectors of size n) by a vector of
// size k, returning a ComplexVector of size n.
ComplexVector MatVecMult(const std::vector<ComplexVector> &X, const Eigen::VectorXcd &y)
{
  MFEM_ASSERT(X.size() == y.size(), "Mismatch in dimension of input vectors!");
  const size_t k = X.size();
  const size_t n = X[0].Size();
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

  // Palace ComplexVectors of size n.
  ComplexVector v, u, w, c, w0, z;
  v.SetSize(n);
  u.SetSize(n);
  w.SetSize(n);
  c.SetSize(n);
  w0.SetSize(n);
  z.SetSize(n);
  v.UseDevice(true);
  u.UseDevice(true);
  w.UseDevice(true);
  c.UseDevice(true);
  w0.UseDevice(true);
  z.UseDevice(true);

  // Eigen Matrix/Vectors for extended operator of size k.
  Eigen::MatrixXcd H;
  Eigen::VectorXcd u2, z2, c2, w2, v2;

  // Storage for eigenpairs.
  std::vector<ComplexVector> X;
  std::vector<std::complex<double>> eigs;

  // Reset any previously computed eigenpairs.
  eigenvalues.clear();
  eigenvectors.clear();
  eigenvectors.reserve(nev);
  X.reserve(nev);

  // Set defaults.
  if (nleps_it <= 0)
  {
    nleps_it = 100;
  }

  // Delta used in to compute divided difference Jacobian.
  const auto delta = std::sqrt(std::numeric_limits<double>::epsilon());

  // Suppress wave port output during Newton iterations.
  space_op->GetWavePortOp().SetSuppressOutput(true);

  // Set a seed and distribution for random Eigen vectors to ensure the same values on all
  // ranks.
  unsigned int seed = nev;
  std::mt19937 gen(seed);
  std::uniform_real_distribution<> dis(-1.0, 1.0);

  const int num_init_guess = init_eigenvalues.size();
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
      eig = init_eigenvalues[guess_idx];
      v = init_eigenvectors[guess_idx];
    }
    else
    {
      eig = sigma;
      space_op->GetRandomInitialVector(v);
    }
    eig_opInv = eig;  // eigenvalue estimate used in the (lagged) preconditioner

    // Set the "random" c vector and the deflation component of the eigenpair initial guess.
    linalg::SetRandom(GetComm(), c, seed);  // set seed for deterministic behavior
    c2.conservativeResize(k);
    v2.conservativeResize(k);
    for (int i = 0; i < k; i++)
    {
      c2(i) = std::complex<double>(dis(gen), dis(gen));
      v2(i) = std::complex<double>(dis(gen), dis(gen));
    }

    // Normalize random c vector.
    double norm_c = std::sqrt(linalg::Norml2(GetComm(), c, true) + c2.squaredNorm());
    c *= 1.0 / norm_c;
    c2 *= 1.0 / norm_c;

    // Normalize eigenvector estimate.
    double norm_v = std::sqrt(linalg::Norml2(GetComm(), v, true) + v2.squaredNorm());
    v *= 1.0 / norm_v;
    v2 *= 1.0 / norm_v;

    // Set the linear solver operators.
    opA2 = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(eig.imag()),
                                                           Operator::DIAG_ZERO);
    opA = space_op->GetSystemMatrix(std::complex<double>(1.0, 0.0), eig, eig * eig, opK,
                                    opC, opM, opA2.get());
    opP = space_op->GetPreconditionerMatrix<ComplexOperator>(std::complex<double>(1.0, 0.0),
                                                             eig, eig * eig, eig.imag());
    opInv->SetOperators(*opA, *opP);

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

    // Compute w0 = T^-1 c and normalize it.
    deflated_solve(c, c2, w0, w2);
    double norm_w0 = std::sqrt(linalg::Norml2(GetComm(), w0, true) + w2.squaredNorm());
    w0 *= 1.0 / norm_w0;
    w2 *= 1.0 / norm_w0;

    // Newton iterations.
    double res = mfem::infinity();
    int it = 0, diverged_it = 0;
    while (it < nleps_it)
    {
      // Compute u = A * v.
      auto A2n = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(eig.imag()),
                                                                 Operator::DIAG_ZERO);
      auto A = space_op->GetSystemMatrix(std::complex<double>(1.0, 0.0), eig, eig * eig,
                                         opK, opC, opM, A2n.get());
      A->Mult(v, u);
      if (k > 0)  // Deflation
      {
        // u1 = T(l) v1 + U(l) v2 = T(l) v1 + T(l)X(lI - H)^-1 v2.
        const Eigen::MatrixXcd S = eig * Eigen::MatrixXcd::Identity(k, k) - H;
        const ComplexVector XSv2 = MatVecMult(X, S.fullPivLu().solve(v2));
        A->AddMult(XSv2, u, 1.0);
        // u2 = X^* v1.
        u2.conservativeResize(k);
        for (int j = 0; j < k; j++)
        {
          u2(j) = linalg::Dot(GetComm(), v, X[j]);
        }
      }

      // Compute residual.
      res = std::sqrt(linalg::Norml2(GetComm(), u, true) + u2.squaredNorm());
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
      auto opA2p = space_op->GetExtraSystemMatrix<ComplexOperator>(
          std::abs(eig.imag()) * (1.0 + delta), Operator::DIAG_ZERO);
      const std::complex<double> denom =
          std::complex<double>(0.0, delta * std::abs(eig.imag()));
      auto opAJ = space_op->GetDividedDifferenceMatrix<ComplexOperator>(
          denom, opA2p.get(), A2n.get(), Operator::DIAG_ZERO);
      auto opJ = space_op->GetSystemMatrix(
          std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0),
          std::complex<double>(2.0, 0.0) * eig, opK, opC, opM, opAJ.get());
      opJ->Mult(v, w);
      if (k > 0)  // Deflation
      {
        // w1 = T'(l) v1 + U'(l) v2 = T'(l) v1 + T'(l)XS v2 - T(l)XS^2 v2.
        const Eigen::MatrixXcd S = eig * Eigen::MatrixXcd::Identity(k, k) - H;
        const Eigen::VectorXcd Sv2 = S.fullPivLu().solve(v2);
        const ComplexVector XSv2 = MatVecMult(X, Sv2);
        const ComplexVector XSSv2 = MatVecMult(X, S.fullPivLu().solve(Sv2));
        opJ->AddMult(XSv2, w, 1.0);
        A->AddMult(XSSv2, w, -1.0);
      }

      // Compute delta = - dot(w0, u) / dot(w0, w).
      const std::complex<double> u2_w0 = std::complex<double>(w2.adjoint() * u2);
      const std::complex<double> delta =
          -(linalg::Dot(GetComm(), u, w0) + u2_w0) / linalg::Dot(GetComm(), w, w0);

      // Update eigenvalue.
      eig += delta;

      // Compute z = -(delta * w + u).
      z.AXPBYPCZ(-delta, w, -1.0, u, 0.0);
      z2 = -u2;

      //  Update preconditioner if needed. Updating the preconditioner as infrequently as
      //  possible gives the best performance and robustness. Updating the preconditioner
      //  close to an eigenvalue can lead to numerical instability.
      if (it > 0 && it % preconditioner_lag == 0 && res > preconditioner_tol)
      {
        eig_opInv = eig;
        opA2 = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(eig.imag()),
                                                               Operator::DIAG_ZERO);
        opA = space_op->GetSystemMatrix(std::complex<double>(1.0, 0.0), eig, eig * eig, opK,
                                        opC, opM, opA2.get());
        opP = space_op->GetPreconditionerMatrix<ComplexOperator>(
            std::complex<double>(1.0, 0.0), eig, eig * eig, eig.imag());
        opInv->SetOperators(*opA, *opP);
        // Recompute w0 and normalize.
        deflated_solve(c, c2, w0, w2);
        double norm_w0 = std::sqrt(linalg::Norml2(GetComm(), w0, true) + w2.squaredNorm());
        w0 *= 1.0 / norm_w0;
        w2 *= 1.0 / norm_w0;
      }

      // Solve M (v_k+1 - v_k) = z.
      deflated_solve(z, z2, u, u2);

      // Update and normalize eigenvector estimate.
      v += u;
      v2 += u2;
      norm_v = std::sqrt(linalg::Norml2(GetComm(), v, true) + v2.squaredNorm());
      v *= 1.0 / norm_v;
      v2 *= 1.0 / norm_v;

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

  space_op->GetWavePortOp().SetSuppressOutput(false);

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
  auto A2 = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(l.imag()),
                                                            Operator::DIAG_ZERO);
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

NewtonInterpolationOperator::NewtonInterpolationOperator(SpaceOperator &space_op)
  : space_op(&space_op)
{
  rhs.SetSize(space_op.GetNDSpace().GetTrueVSize());
  rhs.UseDevice(true);
}

// Compute the elementary symmetric polynomial. Used to convert from Newton to monomial
// basis.
template <typename ScalarType>
ScalarType elementarySymmetric(const std::vector<ScalarType> &points, const int k,
                               const int n)
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

void NewtonInterpolationOperator::Interpolate(const int order,
                                              const std::complex<double> sigma_min,
                                              const std::complex<double> sigma_max)
{
  MFEM_VERIFY(order >= 0, "Interpolation order must be greater than or equal to 0!");

  // Reset operators and sample points each time Interpolate is called.
  num_points = order + 1;
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
        auto A2j = space_op->GetExtraSystemMatrix<ComplexOperator>(points[j].imag(),
                                                                   Operator::DIAG_ZERO);
        ops[k].push_back(std::move(A2j));
      }
      else
      {
        std::complex<double> denom = points[j + k] - points[j];
        auto A2dd = space_op->GetDividedDifferenceMatrix<ComplexOperator>(
            denom, ops[k - 1][j + 1].get(), ops[k - 1][j].get(), Operator::DIAG_ZERO);
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

void NewtonInterpolationOperator::Mult(const int order, const ComplexVector &x,
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

void NewtonInterpolationOperator::AddMult(const int order, const ComplexVector &x,
                                          ComplexVector &y, std::complex<double> a) const
{
  this->Mult(order, x, rhs);
  rhs *= a;
  y += rhs;
}

}  // namespace palace
