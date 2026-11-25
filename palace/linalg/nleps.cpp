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
  preconditioner_lag = std::max(preconditioner_update_freq, 1);  // 0 as 1, every iteration
  preconditioner_tol = preconditioner_update_tol;
}

// Set the maximum number of restarts with the same initial guess.
void QuasiNewtonSolver::SetMaxRestart(int max_num_restart)
{
  max_restart = max_num_restart;
}

void QuasiNewtonSolver::SetExtraSystemMatrix(
    std::function<std::unique_ptr<ComplexOperator>(std::complex<double>)> A2)
{
  funcA2 = A2;
}

void QuasiNewtonSolver::SetPreconditionerUpdate(
    std::function<
        std::unique_ptr<ComplexOperator>(std::complex<double>, std::complex<double>,
                                         std::complex<double>, std::complex<double>)>
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

  // Get eigenmodes initial guesses from linear eigensolver, with λ = iω = i √μ for a
  // linear eigenvalue problem in μ = -λ² = ω².
  auto GetLinearEigenvalue = [this](int i)
  {
    const std::complex<double> l = linear_eigensolver_->GetEigenvalue(i);
    return linear_squared ? std::complex<double>(0.0, 1.0) * std::sqrt(l) : l;
  };
  eigenvalues.resize(nev_linear);
  eigenvectors.resize(nev_linear);
  for (int i = 0; i < nev_linear; i++)
  {
    eigenvalues[i] = GetLinearEigenvalue(i);
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
    const double min_error = *std::min_element(res.get(), res.get() + nev_linear);
    const double threshold = 100.0 * min_error;
    // One partition key keeps a strict weak ordering when residuals equal the threshold.
    std::ranges::sort(indices,
                      [&](const auto i, const auto j)
                      {
                        const bool i_small = res.get()[i] <= threshold;
                        const bool j_small = res.get()[j] <= threshold;
                        if (i_small != j_small)
                        {
                          return i_small;
                        }
                        return eigenvalues[i].imag() < eigenvalues[j].imag();
                      });
  }
  for (int i = 0; i < nev_linear; i++)
  {
    eigenvalues[i] = GetLinearEigenvalue(indices[i]);
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
  // problem size and k is the number of deflated eigenpairs, those of the minimal invariant
  // pair (X, H). The extended operators are never explicitly constructed and two separate
  // vectors of length n and k are used to store the extended solution: [v, v2] where v is a
  // ComplexVector distributed across all processes and v2 is an Eigen::VectorXcd stored
  // redundantly on all processes. An extended eigenpair (λ, [v, v2]) gives the eigenpair
  // (λ, x) of T, x = v + X (λ I - H)⁻¹ v2, and Xᴴ v = 0: v is the part of x orthogonal to
  // the deflated eigenvectors.
  //
  // The initial guesses are refined in the order of SetInitialGuess, with restarts, until
  // the requested number of eigenvalues is found above the target. The iteration from a
  // guess converges to an eigenpair which need not be the one the guess approximates, and
  // the guesses refined first are the most accurate ones, not the nearest: the eigenvalues
  // found can be farther from the target than those of the guesses left. Each guess left is
  // then refined once, nearest the target first, while it is nearer than the requested
  // number of eigenvalues found (see Nearer below), and the nearest of those found are
  // returned.

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
  // be rolled back once normalized). x holds the eigenvector of a converged iterate.
  ComplexVector v, u, w, c, w0, z, du, v_trial, x;
  v.SetSize(n);
  u.SetSize(n);
  w.SetSize(n);
  c.SetSize(n);
  w0.SetSize(n);
  z.SetSize(n);
  du.SetSize(n);
  v_trial.SetSize(n);
  x.SetSize(n);
  v.UseDevice(true);
  u.UseDevice(true);
  w.UseDevice(true);
  c.UseDevice(true);
  w0.UseDevice(true);
  z.UseDevice(true);
  du.UseDevice(true);
  v_trial.UseDevice(true);
  x.UseDevice(true);

  // Eigen Matrix/Vectors for extended operator of size k, and the Gram matrix G = Xᴴ X of
  // the deflated eigenvectors with its factorization.
  Eigen::MatrixXcd H, G;
  Eigen::FullPivLU<Eigen::MatrixXcd> G_lu;
  Eigen::VectorXcd u2, z2, c2, w2, v2, du2, v2_trial;

  // Storage for the deflated eigenvectors, and for the eigenpairs found above the target:
  // the unit eigenvector of a deflated eigenpair is X z, with its coefficients z in the
  // first columns of X, and that of an eigenpair which is not deflated is found_vec[vec].
  struct FoundEigenpair
  {
    std::complex<double> eig;
    Eigen::VectorXcd z;
    int vec;
  };
  std::vector<ComplexVector> X;
  X.reserve(nev);
  std::vector<FoundEigenpair> found;
  std::vector<ComplexVector> found_vec;

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

  // An eigenvector computed to the tolerance which lies in the span of the deflated
  // eigenvectors has a part orthogonal to them at the level of the tolerance (relative to
  // its norm), and one of order one otherwise: the two are told apart at √tol. Eigenvalues
  // within √tol of each other are taken as the same eigenvalue.
  const double span_tol = std::sqrt(rtol);

  // Whether initial guess i is nearer the target than λₙ, the num_requested-th nearest
  // eigenvalue found above it: above the target and below λₙ, by more than the tolerance of
  // equal eigenvalues (the guess of the partner of a degenerate λₙ is not nearer), and
  // nearer the target than λₙ in the complex plane too. The initial guesses solve an
  // approximation of T(λ) near the target frequencies, less accurate away from them: a
  // guess farther from the target than λₙ in the complex plane, that of a more strongly
  // damped mode, is the least likely to converge, and each iteration costs a linear solve.
  const int num_requested = nev;
  auto Nearer = [&](int i)
  {
    std::vector<std::complex<double>> eig_found;
    for (const auto &f : found)
    {
      eig_found.push_back(f.eig);
    }
    std::ranges::nth_element(eig_found, eig_found.begin() + (num_requested - 1),
                             [](auto l, auto r) { return l.imag() < r.imag(); });
    const auto eig_n = eig_found[num_requested - 1];
    return eigenvalues[i].imag() > sigma.imag() &&
           eigenvalues[i].imag() < eig_n.imag() - span_tol * std::abs(eig_n) &&
           std::abs(eigenvalues[i] - sigma) < std::abs(eig_n - sigma);
  };

  // The initial guesses left when num_requested eigenvalues are found above the target,
  // nearest the target first: each is refined once if it is nearer than those found.
  std::vector<int> nearer;
  std::size_t nearer_idx = 0;
  bool refine_nearer = false;

  int k = 0, restart = 0, guess_idx = 0;
  while (true)
  {
    if (!refine_nearer)
    {
      // If > max_restart with the same initial guess, skip to next initial guess.
      // If we tried all initial guesses and the random guess, end search even if fewer
      // eigenpairs than requested were found.
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
      if (static_cast<int>(found.size()) >= num_requested)
      {
        refine_nearer = true;
        nearer.resize(std::max(num_init_guess - guess_idx, 0));
        std::iota(nearer.begin(), nearer.end(), guess_idx);
        std::ranges::stable_sort(nearer, [&](const auto i, const auto j)
                                 { return eigenvalues[i].imag() < eigenvalues[j].imag(); });
      }
    }
    if (refine_nearer)
    {
      while (nearer_idx < nearer.size() && !Nearer(nearer[nearer_idx]))
      {
        nearer_idx++;
      }
      if (nearer_idx == nearer.size())
      {
        break;
      }
      restart = 0;
    }

    // The initial guess of this attempt (num_init_guess for the random guess).
    const int guess = refine_nearer ? nearer[nearer_idx++] : guess_idx;

    // Set the eigenpair estimate to the initial guess.
    std::complex<double> eig, eig_opInv;
    if (guess < num_init_guess)
    {
      eig = eigenvalues[guess];
      v = eigenvectors[guess];
    }
    else
    {
      eig = sigma;
      if (num_init_guess < 3)
      {
        // Set purely random vector, zeroing the essential-BC true dofs. On a
        // partitioned mesh GetEssentialTrueDofs() returns nullptr on ranks whose
        // local partition owns no essential (Dirichlet/PEC) dofs, in which case
        // there is nothing to zero on this rank.
        linalg::SetRandom(GetComm(), v);
        if (const auto *ess_tdof_list =
                dynamic_cast<const ComplexParOperator *>(opK)->GetEssentialTrueDofs())
        {
          linalg::SetSubVector(v, *ess_tdof_list, 0.0);
        }
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

    // With no deflated eigenpair, a new attempt from the same initial guess repeats this
    // one (the same guess and c vector, and no deflation component), except for the random
    // guess made of two initial guesses: the search moves on to the next guess instead.
    const bool same_attempt = (k == 0) && (guess < num_init_guess || num_init_guess < 3);
    const char *next_attempt = refine_nearer ? "next initial guess" : "restarting";

    // Set the "random" c vector and the deflation component of the eigenpair initial guess.
    // The c vector has a seed on each process, from its rank, for deterministic behavior:
    // with one seed on every process, each would hold the same entries (the first of them
    // 8.8e-4 after srand(111) on macOS).
    linalg::SetRandom(GetComm(), c);
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
    opA2 = (*funcA2)(eig);
    opA = BuildOperatorWithA2({1.0 + 0.0i, eig, eig * eig}, {opK, opC, opM}, opA2.get());
    opP = (*funcP)(1.0 + 0.0i, eig, eig * eig,
                   eig / std::complex<double>(0.0, 1.0));  // ω = λ/i
    opInv->SetOperators(*opA, *opP);
    opInv->SetAbsTol(1.0e-12);

    // Linear solve with the extended operator of the deflated problem.
    auto deflated_solve = [&](const ComplexVector &b1, const Eigen::VectorXcd &b2,
                              ComplexVector &x1, Eigen::VectorXcd &x2)
    {
      // Solve the block linear system
      // |T(σ)  T(σ) X S⁻¹| |x1| = |b1|
      // |  Xᴴ         0  | |x2|   |b2|
      // with S = σ I - H and G = Xᴴ X:
      // x1 = T⁻¹ b1 + X G⁻¹ (b2 - Xᴴ T⁻¹ b1)
      // x2 = -S G⁻¹ (b2 - Xᴴ T⁻¹ b1)
      opInv->Mult(b1, x1);
      if (k == 0)  // no deflation
      {
        return;
      }
      Eigen::VectorXcd y(k);
      for (int j = 0; j < k; j++)
      {
        y(j) = b2(j) - linalg::Dot(GetComm(), x1, X[j]);
      }
      y = G_lu.solve(y);
      x2 = -(eig_opInv * Eigen::MatrixXcd::Identity(k, k) - H) * y;
      linalg::AXPY(1.0, MatVecMult(X, y), x1);
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
      A2_out = (*funcA2)(lam);
      auto A =
          BuildOperatorWithA2({1.0 + 0.0i, lam, lam * lam}, {opK, opC, opM}, A2_out.get());
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
    double res_check = res;  // Residual at the previous preconditioner update check

    int it = 0, diverged_it = 0, stalled_windows = 0;
    while (it < nleps_it)
    {
      if (print > 0)
      {
        Mpi::Print(GetComm(),
                   "{:d} NLEPS (nconv={:d}, restart={:d}) residual norm {:.6e}\n", it, k,
                   restart, res);
      }

      // End if residual below tolerance, for the extended eigenpair and for the eigenpair
      // (λ, x) of T itself, T(λ) x = u: with deflation, x is not normalized, and ‖x‖ < 1
      // when x lies near the span of the deflated eigenvectors.
      double norm_x = 0.0;
      Eigen::VectorXcd Sv2;  // (λ I - H)⁻¹ v2
      if (res < rtol)
      {
        x = v;
        if (k > 0)
        {
          const Eigen::MatrixXcd S = eig * Eigen::MatrixXcd::Identity(k, k) - H;
          Sv2 = S.fullPivLu().solve(v2);
          linalg::AXPY(1.0, MatVecMult(X, Sv2), x);
        }
        norm_x = linalg::Norml2(GetComm(), x);
      }
      if (res < rtol && linalg::Norml2(GetComm(), u) < rtol * norm_x)
      {
        const bool in_range = (eig.imag() > sigma.imag());
        // The eigenpair is deflated if its eigenvector is not in the span of the deflated
        // eigenvectors (eigenvectors of distinct eigenvalues of a nonlinear eigenvalue
        // problem can be linearly dependent). In that span, the invariant pair extended
        // with it would not be minimal: its new column v / ‖v‖ would be rounding error,
        // with the entries of H in 1 / ‖v‖, and the deflated problem of the next
        // eigenpairs would be wrong.
        const double scale = linalg::Norml2(GetComm(), v);
        const bool deflate = (scale >= span_tol * norm_x);
        // An eigenpair which is not deflated can be found again: it is one found before
        // if its eigenvalue and eigenvector are those of one.
        bool found_before = false;
        if (in_range && !deflate)
        {
          Eigen::VectorXcd Xx;  // Xᴴ x, for the eigenvectors of the deflated eigenpairs
          for (const auto &f : found)
          {
            if (std::abs(f.eig - eig) > span_tol * std::abs(f.eig))
            {
              continue;
            }
            std::complex<double> overlap;
            if (f.vec >= 0)
            {
              overlap = linalg::Dot(GetComm(), x, found_vec[f.vec]);
            }
            else
            {
              if (Xx.size() == 0)
              {
                Xx.resize(k);
                for (int j = 0; j < k; j++)
                {
                  Xx(j) = linalg::Dot(GetComm(), x, X[j]);
                }
              }
              overlap = f.z.dot(Xx.head(f.z.size()));
            }
            if (std::abs(overlap) >= 0.99 * norm_x)
            {
              found_before = true;
              break;
            }
          }
        }
        if (print > 0)
        {
          Mpi::Print(GetComm(),
                     "Eigenvalue {:d}, Quasi-Newton converged in {:d} iterations "
                     "({:.3e}{:+.3e}i){}.\n",
                     k, it, eig.real(), eig.imag(),
                     found_before ? ", found before"
                     : !deflate   ? ", in the span of the deflated eigenvectors"
                                  : "");
        }
        if (in_range && !found_before)
        {
          auto &f = found.emplace_back();
          f.eig = eig;
          if (deflate)
          {
            // The coefficients of x = v + X (λ I - H)⁻¹ v2 in X, with v the new column.
            f.z.resize(k + 1);
            f.z.head(k) = Sv2;
            f.z(k) = scale;
            f.z *= 1.0 / norm_x;
            f.vec = -1;
          }
          else
          {
            found_vec.push_back(x);
            found_vec.back() *= 1.0 / norm_x;
            f.vec = static_cast<int>(found_vec.size()) - 1;
          }
        }
        if (deflate)
        {
          // Update the invariant pair with normalization.
          v *= 1.0 / scale;
          X.resize(k + 1);
          X[k] = v;
          H.conservativeResizeLike(Eigen::MatrixXd::Zero(k + 1, k + 1));
          H.col(k).head(k) = v2 / scale;
          H(k, k) = eig;
          G.conservativeResize(k + 1, k + 1);
          for (int j = 0; j <= k; j++)
          {
            G(j, k) = linalg::Dot(GetComm(), X[k], X[j]);  // X[j]ᴴ X[k]
            G(k, j) = std::conj(G(j, k));
          }
          G_lu.compute(G);
          k++;
        }
        if (!refine_nearer)
        {
          // If the eigenvalue is new and inside the desired range, increment initial guess
          // index. Otherwise, use the same initial guess again: in a new attempt once the
          // eigenpair is deflated, and else in a restart.
          if (in_range && !found_before)
          {
            guess_idx++;
            restart = 0;
          }
          else if (deflate)
          {
            restart = 0;
          }
          else
          {
            restart++;
          }
        }
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
                     "{}.\n",
                     k, it, next_attempt);
        }
        restart = same_attempt ? max_restart + 1 : restart + 1;
        break;
      }

      // Compute w = J * v. The A2 finite difference is a matrix-free sum (each A2 may
      // itself be a matrix-free SumComplexOperator carrying the wave-port correction).
      auto opA2p = (*funcA2)(eig * (1.0 + delta));
      const std::complex<double> denom = delta * eig;
      auto opAJ = std::make_unique<SumComplexOperator>(opA2p->Height(), opA2p->Width());
      opAJ->AddOperator(*opA2p, 1.0 / denom);
      opAJ->AddOperator(*A2n, -1.0 / denom);
      auto opJ = BuildOperatorWithA2({0.0 + 0.0i, 1.0 + 0.0i, 2.0 * eig}, {opK, opC, opM},
                                     opAJ.get());
      opJ->Mult(v, w);
      if (k > 0)  // Deflation
      {
        // w1 = T'(l) v1 + U'(l) v2 = T'(l) v1 + T'(l)XS v2 - T(l)XS^2 v2. Scoping T(l)
        // here lets the line search overwrite A2n freely; with no deflation we skip it.
        auto A =
            BuildOperatorWithA2({1.0 + 0.0i, eig, eig * eig}, {opK, opC, opM}, A2n.get());
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
      // by preconditioner_lag iterations. Below preconditioner_tol, the operator of the
      // last update is kept while the iteration converges fast enough with it: it then
      // converges linearly, at a rate set by the distance of eig_opInv from the eigenvalue.
      // It is updated if, at the rate since the previous check, the residual would not
      // reach the tolerance within the next preconditioner_lag iterations.
      bool update = false;
      if (it > 0 && it % preconditioner_lag == 0)
      {
        update = (res > preconditioner_tol) || (res * (res / res_check) > rtol);
        // A window of preconditioner_lag iterations stalls if it reduces the residual by
        // less than 10% (below preconditioner_tol, the preconditioner is updated after such
        // a window). Only the attempts from the nearer initial guesses, once num_requested
        // eigenvalues are found, are ended when they stall: until then, an attempt which
        // stalls for a while can still converge to an eigenvalue the next ones would miss.
        stalled_windows =
            (refine_nearer && res > 0.9 * res_check) ? stalled_windows + 1 : 0;
        res_check = res;
      }
      // An attempt from a nearer initial guess which stalls in two windows in a row is
      // ended.
      if (stalled_windows >= 2)
      {
        if (print > 0)
        {
          Mpi::Print(GetComm(),
                     "Eigenvalue {:d}, Quasi-Newton stalled after {:d} iterations, {}.\n",
                     k, it + 1, next_attempt);
        }
        break;
      }
      if (update)
      {
        eig_opInv = eig;
        opA2 = (*funcA2)(eig_opInv);
        opA = BuildOperatorWithA2({1.0 + 0.0i, eig_opInv, eig_opInv * eig_opInv},
                                  {opK, opC, opM}, opA2.get());
        opP = (*funcP)(1.0 + 0.0i, eig_opInv, eig_opInv * eig_opInv,
                       eig_opInv / std::complex<double>(0.0, 1.0));
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
                     "{}.\n",
                     k, nleps_it, next_attempt);
        }
        restart = same_attempt ? max_restart + 1 : restart + 1;
      }
    }
  }

  // The requested eigenpairs are the num_requested nearest the target of those found above
  // it, and all those found with the eigenvalue of the last one, so that a degenerate
  // eigenvalue is not split.
  std::vector<std::size_t> order(found.size());
  std::iota(order.begin(), order.end(), 0);
  std::ranges::stable_sort(order, [&](auto l, auto r)
                           { return found[l].eig.imag() < found[r].eig.imag(); });
  std::size_t num_return = std::min(order.size(), static_cast<std::size_t>(num_requested));
  while (num_return < order.size() &&
         std::abs(found[order[num_return]].eig - found[order[num_return - 1]].eig) <=
             span_tol * std::abs(found[order[num_return - 1]].eig))
  {
    num_return++;
  }
  order.resize(num_return);
  eigenvalues.clear();
  eigenvectors.clear();
  for (const auto i : order)
  {
    eigenvalues.push_back(found[i].eig);
    if (found[i].vec >= 0)
    {
      eigenvectors.push_back(std::move(found_vec[found[i].vec]));
    }
    else
    {
      Eigen::VectorXcd z_k = Eigen::VectorXcd::Zero(k);
      z_k.head(found[i].z.size()) = found[i].z;
      eigenvectors.push_back(MatVecMult(X, z_k));
    }
  }
  nev = eigenvalues.size();

  // The reported eigenpairs are ordered by their imaginary part.
  perm = std::make_unique<int[]>(nev);
  std::iota(perm.get(), perm.get() + nev, 0);

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
  auto A2 = (*funcA2)(l);
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
    std::function<std::unique_ptr<ComplexOperator>(std::complex<double>)> funcA2, int size)
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
  frozen_ops.clear();

  // Linearly spaced sample points.
  const auto n_intervals = static_cast<double>(num_points - 1);
  for (int j = 0; j < num_points; j++)
  {
    points[j] = sigma_min + static_cast<double>(j) * (sigma_max - sigma_min) / n_intervals;
  }

  // Build divided difference matrices.
  for (int k = 0; k < num_points; k++)
  {
    for (int j = 0; j < num_points - k; j++)
    {
      if (k == 0)
      {
        auto A2j = (funcA2)(points[j]);
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

  // Assemble the single per-order (coefficient, operator) representation read by
  // GetInterpolationOperator, Mult, and AddMult: the order-th pencil coefficient operator
  // is Σₜ pencil_coeffs[order][t]·pencil_ops[order][t]. The base terms are the
  // Newton→monomial combination of the leading divided-difference operators;
  // AddFrozenPole appends one correction term per frozen pole.
  pencil_coeffs.assign(num_points, {});
  pencil_ops.assign(num_points, {});
  for (int order = 0; order < num_points; order++)
  {
    for (int j = 0; j < num_points; j++)
    {
      pencil_coeffs[order].push_back(coeffs[order][j]);
      pencil_ops[order].push_back(ops[j][0].get());
    }
  }
}

std::vector<std::complex<double>>
NewtonInterpolationOperator::ComputeLeadingDividedDifferences(
    const std::function<std::complex<double>(std::complex<double>)> &f) const
{
  // In-place leading divided differences: after pass k, dd[j] = f[x_{j-k}, ..., x_j] for
  // j >= k, so on completion dd[j] = f[x_0, ..., x_j] — the scalar analogue of the
  // operator divided differences ops[j][0], with the same denominators
  // points[j] - points[j - k] as the operator DAG in Interpolate().
  std::vector<std::complex<double>> dd(num_points);
  for (int j = 0; j < num_points; j++)
  {
    dd[j] = f(points[j]);
  }
  for (int k = 1; k < num_points; k++)
  {
    for (int j = num_points - 1; j >= k; j--)
    {
      dd[j] = (dd[j] - dd[j - 1]) / (points[j] - points[j - k]);
    }
  }
  return dd;
}

std::vector<std::complex<double>> NewtonInterpolationOperator::ComputeMonomialCoefficients(
    const std::vector<std::complex<double>> &dd) const
{
  // q[order] = Σⱼ coeffs[order][j]·dd[j], the same Newton→monomial combination the
  // operator path applies to ops[j][0], so the scalar Newton interpolant through the
  // nodes is p(λ) = Σ q[order]·λ^order.
  std::vector<std::complex<double>> q(num_points, 0.0);
  for (int order = 0; order < num_points; order++)
  {
    for (int j = 0; j < num_points; j++)
    {
      q[order] += coeffs[order][j] * dd[j];
    }
  }
  return q;
}

void NewtonInterpolationOperator::AddFrozenPole(
    std::unique_ptr<ComplexOperator> M,
    const std::function<std::complex<double>(std::complex<double>)> &f,
    std::complex<double> lambda_target)
{
  if (!M)
  {
    return;
  }
  MFEM_VERIFY(!points.empty() && !coeffs.empty(),
              "AddFrozenPole requires a prior Interpolate() call!");

  // Per-order multiplier of M: remove the interpolated pole at every order (the scalar
  // Newton interpolant of f through the same nodes as the operator path), then re-add the
  // pole frozen at the target into order 0. By linearity of the interpolation, this
  // yields the smooth remainder pencil with the pole frozen in K. The correction terms
  // are appended to the shared pencil representation, so GetInterpolationOperator, Mult,
  // and AddMult all see them.
  const auto q = ComputeMonomialCoefficients(ComputeLeadingDividedDifferences(f));
  const ComplexOperator *M_raw = M.get();
  frozen_ops.push_back(std::move(M));
  for (int order = 0; order < num_points; order++)
  {
    std::complex<double> corr = -q[order];
    if (order == 0)
    {
      corr += f(lambda_target);
    }
    pencil_coeffs[order].push_back(corr);
    pencil_ops[order].push_back(M_raw);
  }
}

bool NewtonInterpolationOperator::DetermineFrozen(
    const std::function<std::complex<double>(std::complex<double>)> &f_full,
    const std::function<std::complex<double>(std::complex<double>)> &f_frozen,
    std::complex<double> lambda_target, double &fit_err, double &freeze_err) const
{
  MFEM_VERIFY(!points.empty() && !coeffs.empty(),
              "DetermineFrozen requires a prior Interpolate() call!");

  // Monomial coefficients of the scalar Newton interpolants through the same nodes as the
  // operator path: the full coefficient (fit branch) and its smooth part f_full - f_frozen
  // (what AddFrozenPole(M, f_frozen, lambda_target) leaves interpolated, exact when that
  // part is a polynomial of degree < num_points).
  const auto q_full = ComputeMonomialCoefficients(ComputeLeadingDividedDifferences(f_full));
  const auto q_smooth = ComputeMonomialCoefficients(ComputeLeadingDividedDifferences(
      [&](std::complex<double> lambda) { return f_full(lambda) - f_frozen(lambda); }));
  auto eval = [&](const std::vector<std::complex<double>> &q, std::complex<double> lambda)
  {
    // Horner evaluation of Σ q[order]·λ^order.
    std::complex<double> val = 0.0;
    for (int order = num_points - 1; order >= 0; order--)
    {
      val = val * lambda + q[order];
    }
    return val;
  };

  // Max relative error over a dense sample of the window. The freeze branch is what
  // AddFrozenPole(M, f_frozen, lambda_target) produces: the interpolant of the smooth part
  // plus the frozen value f_frozen(lambda_target).
  const std::complex<double> frozen_target = f_frozen(lambda_target);
  constexpr int n_samples = 101;
  const auto n_intervals = static_cast<double>(n_samples - 1);
  fit_err = freeze_err = 0.0;
  for (int i = 0; i < n_samples; i++)
  {
    const std::complex<double> lambda =
        points.front() +
        static_cast<double>(i) * (points.back() - points.front()) / n_intervals;
    const std::complex<double> f_exact = f_full(lambda);
    const double scale = std::max(std::abs(f_exact), 1.0e-300);
    fit_err = std::max(fit_err, std::abs(eval(q_full, lambda) - f_exact) / scale);
    freeze_err = std::max(
        freeze_err, std::abs(eval(q_smooth, lambda) + frozen_target - f_exact) / scale);
  }
  return freeze_err < fit_err;
}

std::unique_ptr<ComplexOperator>
NewtonInterpolationOperator::GetInterpolationOperator(int order) const
{
  MFEM_VERIFY(order >= 0 && order < num_points,
              "Order must be greater than or equal to 0 and smaller than the number of "
              "interpolation points!");
  return BuildParSumOperator(pencil_coeffs[order], pencil_ops[order], true);
}

void NewtonInterpolationOperator::Mult(int order, const ComplexVector &x,
                                       ComplexVector &y) const
{
  MFEM_VERIFY(order >= 0 && order < num_points,
              "Order must be greater than or equal to 0 and smaller than the number of "
              "interpolation points!");

  y = 0.0;
  for (std::size_t t = 0; t < pencil_ops[order].size(); t++)
  {
    if (pencil_coeffs[order][t] != 0.0)
    {
      pencil_ops[order][t]->AddMult(x, y, pencil_coeffs[order][t]);
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
