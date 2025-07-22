// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "nleps.hpp"

#include <algorithm>
#include <Eigen/Dense> // test for deflation
#include <string>
#include <mfem.hpp>
#include <mfem/general/forall.hpp>
#include "linalg/divfree.hpp"
#include "utils/communication.hpp"


namespace palace::nleps
{

// Base class methods.

NonLinearEigenvalueSolver::NonLinearEigenvalueSolver(MPI_Comm comm, int print)
  : comm(comm), print(print)
{
  // Initialization.
  nleps_it = 0;

  gamma = delta = 1.0; // get that to work!?
  sigma = 0.0;

  opInv = nullptr;
  opProj = nullptr;
  opB = nullptr; // figure out whether we should keep or not
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

void NonLinearEigenvalueSolver::SetLinearSolver(/*const*/ ComplexKspSolver &ksp)
{
  opInv = &ksp;
}

void NonLinearEigenvalueSolver::SetIoData(const IoData &iodata)
{
  opIodata = &iodata;
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
  }
  if (num_eig > 0)
  {
    eigenvalues.resize(num_eig);
    eigenvectors.resize(num_eig);
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
  MFEM_VERIFY(
      n > 0,
      "Must call SetOperators before using SetInitialSpace for nonlinear eigenvalue solver!");
  if (!r)
  {
    r = std::make_unique<std::complex<double>[]>(n);
  }
  MFEM_VERIFY(v.Size() == n, "Invalid size mismatch for provided initial space vector!");
  v.Get(r.get(), n, false);
}

void NonLinearEigenvalueSolver::SetInitialGuess(const std::vector<std::complex<double>> init_eig, const std::vector<ComplexVector> &init_V)
{
  MFEM_VERIFY(
      n > 0,
      "Must call SetOperators before using SetInitialguess for nonlinear eigenvalue solver!");
  MFEM_VERIFY(init_eig.size() == init_V.size(), "SetInitialGuess requires the same number of eigenvalues and eigenvectors\n");

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
  //MFEM_VERIFY(eig && i >= 0 && i < nev,
  MFEM_VERIFY(i >= 0 && i < nev,
              "Out of range eigenpair requested (i = " << i << ", nev = " << nev << ")!");
  const int &j = perm.get()[i];
  return eigenvalues[j];//eig.get()[j];
}

void NonLinearEigenvalueSolver::GetEigenvector(int i, ComplexVector &x) const
{
  //MFEM_VERIFY(eig && i >= 0 && i < nev,
  MFEM_VERIFY(i >= 0 && i < nev,
              "Out of range eigenpair requested (i = " << i << ", nev = " << nev << ")!");
  MFEM_VERIFY(x.Size() == n, "Invalid size mismatch for provided eigenvector!");
  const int &j = perm.get()[i];
  //x.Set(V.get() + j * n, n, false);
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
 // MFEM_VERIFY(eig && i >= 0 && i < nev,
  MFEM_VERIFY(i >= 0 && i < nev,
              "Out of range eigenpair requested (i = " << i << ", nev = " << nev << ")!");
  const int &j = perm.get()[i];
  switch (type)
  {
    case ErrorType::ABSOLUTE:
      return res.get()[j];
    case ErrorType::RELATIVE:
      //return res.get()[j] / std::abs(eig.get()[j]);
      return res.get()[j] / std::abs(eigenvalues[j]);
    case ErrorType::BACKWARD:
      //return res.get()[j] / GetBackwardScaling(eig.get()[j]);
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
    //x1.Set(V.get() + i * n, n, false);
    x1 = eigenvectors[i];
    xscale.get()[i] = 1.0 / GetEigenvectorNorm(x1, y1);
    //res.get()[i] = GetResidualNorm(eig.get()[i], x1, y1) / linalg::Norml2(comm, x1);
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

void QuasiNewtonSolver::SetOperators(SpaceOperator &space_op_ref,
                                    const ComplexOperator &K, const ComplexOperator &C,
                                    const ComplexOperator &M, EigenvalueSolver::ScaleType type)
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
    normC = linalg::SpectralNorm(comm, *opC, opC->IsReal());
    normM = linalg::SpectralNorm(comm, *opM, opM->IsReal());
    MFEM_VERIFY(normK >= 0.0 && normC >= 0.0 && normM >= 0.0,
                "Invalid matrix norms for Quasi-Newton scaling!");
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
  //z1.SetSize(opK->Height());
  x1.UseDevice(true);
  x2.UseDevice(true);
  y1.UseDevice(true);
  y2.UseDevice(true);
  //z1.UseDevice(true);
  n = opK->Height();
}

namespace
{
  // z = X * y where X is a vector of size k or ComplexVector's of size n, y is a vector of size k
  // z will be a ComplexVector of size n
  //ComplexVector MatVecMult(const std::vector<ComplexVector> &X, const std::vector<std::complex<double>> y, bool on_dev)
  ComplexVector MatVecMult(const std::vector<ComplexVector> &X, const Eigen::VectorXcd &y, bool on_dev)
  {
    MFEM_ASSERT(X.size() == y.size(), "Mismatch in dimension of input vectors!");
    const size_t k = X.size();
    const size_t n = X[0].Size();
    ComplexVector z; z.SetSize(n); z.UseDevice(on_dev); // not sure about on_dev?
    z = 0.0;
    auto *zr = z.Real().Write(on_dev);//z.Real().HostWrite(); // not sure about on_dev
    auto *zi = z.Imag().Write(on_dev);//HostWrite();
    for (int j = 0; j < k; j++)
    {
      auto *XR = X[j].Real().Read(on_dev); //not sure about on_dev. Need to figure out if X is on device or not?
      auto *XI = X[j].Imag().Read(on_dev);
      mfem::forall_switch(on_dev, n,
                          [=] MFEM_HOST_DEVICE(int i)
                          {
                            zr[i] += y(j).real() * XR[i] - y(j).imag() * XI[i];
                            zi[i] += y(j).imag() * XR[i] + y(j).real() * XI[i];
                          });
    }
    return z;
  }

} // test


int QuasiNewtonSolver::Solve()
{
  // Quasi-Newton method for nonlinear eigenvalue problems.
  // Reference: Jarlebring, Koskela, Mele, Disguised and new quasi-Newton methods for nonlinear
  //            eigenvalue problems, Numerical Algorithms (2018).
  // Using the deflation scheme used by SLEPc's NEP solver with minimality index set to 1.
  // Reference: Effenberger, Robust successive computation of eigenpairs for nonlinear
  //            eigenvalue problems, SIAM J. Matrix Anal. Appl. (2013).
  // The deflation scheme solves an extended problem of size n + k, where n is the original
  // problem size and k is the number of converged eigenpairs. The extended operators are never
  // explicitly constructed and two separate vectors of length n and k are used to store the
  // extended solution: [v, v2] where v is a ComplexVector distributed across all processes
  // and v2 is an Eigen::VectorXcd stored redundantly on all processes.

  // Use x1, x2, y1, y2 instead of some of these?!
  ComplexVector v, u, w, c, w0, z;
  std::vector<ComplexVector> X;
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

  const int update_freq = 4;// SHOULD REVISIT THIS AND FIGURE OUT BEST FREQUENCY around 4-5 seems good?

  perm = std::make_unique<int[]>(nev); // use this or order below??

  const auto dl = std::sqrt(std::numeric_limits<double>::epsilon()); // TODO: define only once above
  space_op->GetWavePortOp().SetSuppressOutput(true); //suppressoutput?

  Eigen::MatrixXcd H;
  Eigen::VectorXcd u2, z2;

  int k = 0;
  while (k < nev)
  {
    // do we want set c inside the loop so it's random and thus avoids getting stuck on a "BAD" guess?
    linalg::SetRandom(GetComm(), c); // can it be made deterministc? seed = k?
    std::srand((unsigned int) time(0)); // ?don't like this?! ensures all processes have the same but it's ugly
    Eigen::VectorXcd c2 = Eigen::VectorXcd::Random(k); // is this one already deterministic?
    Eigen::VectorXcd w2;

    // Set the initial guess.
    std::complex<double> eig, eig_opInv;
    if (k < init_eigenvalues.size())
    {
      eig = init_eigenvalues[k];
      v = init_eigenvectors[k];
    }
    else if (init_eigenvalues.size() > 0)
    {
      const int last_idx = init_eigenvalues.size() - 1;
      eig = init_eigenvalues[last_idx];
      v = init_eigenvectors[last_idx];
    }
    else
    {
      eig = sigma;
      linalg::SetRandom(GetComm(), v); // pass a seed or not?
    }
    eig_opInv = eig;
    // The deflation component is always random.
    std::srand((unsigned int) time(0)); // don't like this!?
    Eigen::VectorXcd v2 = Eigen::VectorXcd::Random(k);

    // Normalize eigenvector estimate.
    double norm_v = std::sqrt(linalg::Norml2(GetComm(), v, true) + v2.squaredNorm());
    v *= 1.0 / norm_v;
    v2 *= 1.0 / norm_v;

    opA2 = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(eig.imag()), Operator::DIAG_ZERO);
    opA = space_op->GetSystemMatrix(std::complex<double>(1.0, 0.0), eig, eig * eig, opK, opC, opM, opA2.get());
    opP = space_op->GetPreconditionerMatrix<ComplexOperator>(std::complex<double>(1.0, 0.0), eig, eig * eig, eig.imag());
    opInv->SetOperators(*opA, *opP);

    auto deflated_solve = [&](const std::complex<double> lambda, const ComplexVector &b1, const Eigen::VectorXcd &b2, ComplexVector &x1, Eigen::VectorXcd &x2)
    {
      // Solve the block linear system
      // |T(σ) U(σ)| |x1| = |b1|
      // |A(σ) B(σ)| |x2|   |b2|
      // x1 = T^-1 b1
      // x2 = SS^-1 (b2 - A x1) where SS = (B - A T^-1 U) = - X^* X S^-1
      // x1 = x1 - X S x2
      opInv->Mult(b1, x1);
      if (k == 0)
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
          SS(i,j) = linalg::Dot(GetComm(), X[i], X[j]);
        }
      }
      const auto S = (lambda * Eigen::MatrixXcd::Identity(k, k) - H).inverse();
      SS = - SS * S;
      x2 = SS.inverse() * x2;
      ComplexVector XSx2 = MatVecMult(X, S * x2, true);
      linalg::AXPY(-1.0, XSx2, x1);
    };
    deflated_solve(eig_opInv, c, c2, w0, w2);

    // Newton iterations.
    double res = mfem::infinity();
    int it = 0, diverged_it = 0;
    while (it < nleps_it)
    {
      // Compute u = A * v.
      auto A2n = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(eig.imag()), Operator::DIAG_ZERO);
      auto A = space_op->GetSystemMatrix(std::complex<double>(1.0, 0.0), eig, eig * eig, opK, opC, opM, A2n.get());
      A->Mult(v, u);
      if (k > 0) // Deflation
      {
        // u1 = T(l) v1 + U(l) v2 = T(l) v1 + T(l)X(lI - H)^-1 v2
        const auto S = (eig * Eigen::MatrixXcd::Identity(k, k) - H).inverse();
        ComplexVector XSv2 = MatVecMult(X, S * v2, true);
        A->AddMult(XSv2, u, 1.0);
        // u2 = X^* v1
        u2.conservativeResize(k);
        for (int j = 0; j < k; j++)
        {
          u2(j) = linalg::Dot(GetComm(), v, X[j]);
        }
      }

      // Compute residual.
      res = std::sqrt(linalg::Norml2(GetComm(), u, true) + u2.squaredNorm());
      Mpi::Print(GetComm(), "k: {}, it: {}, eig: {}, {}, res: {}\n", k, it, eig.real(), eig.imag(), res); // only print if print > 0? >1 ????

      // End if residual below tolerance and eigenvalue above the target.
      if (res < rtol)
      {
        if (eig.imag() > sigma.imag())
        {
          // Update the invariant pair with normalization.
          const auto scale = linalg::Norml2(GetComm(), v);
          v *= 1.0 / scale;
          X.push_back(v);
          H.conservativeResizeLike(Eigen::MatrixXd::Zero(k + 1, k + 1));
          H.col(k).head(k) = v2 / scale;
          H(k, k) = eig;
          eigenvalues[k] = eig;
          k++;
        }
        break;
      }
      // Stop if large residual for 10 consecutive iterations.
      diverged_it = (res > 0.9) ? diverged_it + 1 : 0;
      if (diverged_it > 10)
      {
        Mpi::Print(GetComm(), "Eigenvalue {}, Quasi-Newton not converging after {} iterations, restarting.\n", k, it);
        break;
      }

      // Compute w = J * v.
      auto opA2p = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(eig.imag()) * (1.0 + dl), Operator::DIAG_ZERO);
      std::complex<double> denom = dl * std::abs(eig.imag());
      auto opAJ = space_op->GetDividedDifferenceMatrix<ComplexOperator>(denom, opA2p.get(), A2n.get(), Operator::DIAG_ZERO);
      auto opJ = space_op->GetSystemMatrix(std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0), std::complex<double>(2.0, 0.0) * eig, opK, opC, opM, opAJ.get());
      opJ->Mult(v, w);
      if (k > 0) // Deflation
      {
        // w1 = T'(l) v1 + U'(l) v2 = T'(l) v1 + T'(l)XS v2 - T(l)XS^2 v2
        const auto S = (eig * Eigen::MatrixXcd::Identity(k, k) - H).inverse(); // should re-use the one defined above!!
        ComplexVector XSv2 = MatVecMult(X, S * v2, true);
        ComplexVector XSSv2 = MatVecMult(X, S * S * v2, true);
        opJ->AddMult(XSv2, w, 1.0);
        A->AddMult(XSSv2, w, -1.0);
      }

      // Compute delta = - dot(w0, u) / dot(w0, w).
      std::complex<double> delta;
      std::complex<double> u2_w0(0.0, 0.0);
      if (k > 0)
      {
        u2_w0 = std::complex<double>(w2.adjoint() * u2);
      }
      delta = - (linalg::Dot(GetComm(), u, w0) + u2_w0) / linalg::Dot(GetComm(), w, w0);

      // Update eigenvalue.
      eig += delta;

      // Compute z = -(delta * w + u).
      z.AXPBYPCZ(-delta, w, std::complex<double>(-1.0, 0.0), u, std::complex<double>(0.0, 0.0));
      z2 = -u2;

      //  M (x_k+1 - x_k) = z.
      if (it > 0 && it % update_freq == 0)
      {
        eig_opInv = eig;
        opA2 = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(eig.imag()), Operator::DIAG_ZERO);
        opA = space_op->GetSystemMatrix(std::complex<double>(1.0, 0.0), eig, eig * eig, opK, opC, opM, opA2.get());
        opP = space_op->GetPreconditionerMatrix<ComplexOperator>(std::complex<double>(1.0, 0.0), eig, eig * eig, eig.imag());
        opInv->SetOperators(*opA, *opP);
        deflated_solve(eig_opInv, c, c2, w0, w2);
      }
      deflated_solve(eig_opInv, z, z2, u, u2);

      // Update and normalize eigenvector estimate.
      v2 += u2;
      v += u;
      norm_v = std::sqrt(linalg::Norml2(GetComm(), v, true) + v2.squaredNorm());
      v *= 1.0 / norm_v;
      v2 *= 1.0 / norm_v;

      it++;
      if (it == nleps_it)
      {
        Mpi::Print(GetComm(), "Eigenvalue {}, Quasi-Newton did not converge in {} iterations, restarting.\n", k, nleps_it); // only if print > 0? > 1??
      }
    }
    Mpi::Print(GetComm(), "Eigenvalue {}, Quasi-Newton converged in {} iterations.\n", k, it); // only if print > 0?
  }

  // Eigenpair extraction from the invariant pair (X, H).
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eps;
  eps.compute(H);

  std::vector<int> order(eigenvalues.size()), order_eigen(eps.eigenvalues().size()), order2(eigenvalues.size());
  std::iota(order.begin(), order.end(), 0);
  std::iota(order_eigen.begin(), order_eigen.end(), 0);
  std::iota(order2.begin(), order2.end(), 0);
  std::sort(order.begin(), order.end(),
            [&eigen_values = this->eigenvalues](auto l, auto r) { return eigen_values[l].imag() < eigen_values[r].imag(); });
  std::sort(order_eigen.begin(), order_eigen.end(),
            [&epseig = eps.eigenvalues()](auto l, auto r) { return epseig(l).imag() < epseig(r).imag(); });
  std::sort(order2.begin(), order2.end(), [&order](auto l, auto r) { return order[l] < order[r]; });

  // Sort Eigen eigenvectors
  std::vector<Eigen::VectorXcd> Xeig;
  for (int k = 0; k < nev; k++)
  {
    perm[k] = order[k]; // stupid, only use one!
    Xeig.push_back(eps.eigenvectors().col(order_eigen[k]));
  }

  for(int k = 0; k < nev; k++)
  {
    ComplexVector eigv = MatVecMult(X, Xeig[order2[k]], true);
    eigenvectors[k] = eigv;
  }

  // Compute the eigenpair residuals for eigenvalue λ.
  RescaleEigenvectors(nev);

  space_op->GetWavePortOp().SetSuppressOutput(false); //reset?
  //
  return nev;
}

//void QuasiNewtonSolver::ApplyOp(const std::complex<double> *px,
//                              std::complex<double> *py) const
//{
//
//}

//void QuasiNewtonSolver::ApplyOpB(const std::complex<double> *px,
//                               std::complex<double> *py) const
//{
//  MFEM_VERIFY(opB, "No B operator for weighted inner product in QuasiNewton solve!");
//  x1.Set(px, n, false);
//  x2.Set(px + n, n, false);
//  opB->Mult(x1.Real(), y1.Real());
//  opB->Mult(x1.Imag(), y1.Imag());
//  opB->Mult(x2.Real(), y2.Real());
//  opB->Mult(x2.Imag(), y2.Imag());
//  y1 *= delta * gamma * gamma;
//  y2 *= delta * gamma * gamma;
//  y1.Get(py, n, false);
//  y2.Get(py + n, n, false);
//}

double QuasiNewtonSolver::GetResidualNorm(std::complex<double> l, const ComplexVector &x,
                                        ComplexVector &r) const
{
  // Compute the i-th eigenpair residual: || P(λ) x ||₂ = || (K + λ C + λ² M + A2(λ)) x ||₂ for
  // eigenvalue λ.
  opK->Mult(x, r);
  opC->AddMult(x, r, l);
  opM->AddMult(x, r, l * l);
  // if (has_A2) ???
  auto A2 = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(l.imag()), Operator::DIAG_ZERO);
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

}  // namespace palace::nleps
