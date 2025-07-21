// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "nleps.hpp"

#include <algorithm>
#include <Eigen/Dense> // test for deflation
#include <string>
#include <mfem.hpp>
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
      return res.get()[j] / std::abs(eig.get()[j]);
    case ErrorType::BACKWARD:
      return res.get()[j] / GetBackwardScaling(eig.get()[j]);
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
    res.get()[i] = GetResidualNorm(eig.get()[i], x1, y1) / linalg::Norml2(comm, x1);
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



int QuasiNewtonSolver::Solve()
{
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

  perm = std::make_unique<int[]>(nev); // use this or order below??

  const auto dl = std::sqrt(std::numeric_limits<double>::epsilon()); // TODO: define only once above
  space_op->GetWavePortOp().SetSuppressOutput(true); //suppressoutput?

  Eigen::MatrixXcd H;
  std::vector<ComplexVector> X; // if this is always the same as eigenvectors, use only one of the two!
  Eigen::VectorXcd u2, z2;
  Mpi::Print("nev: {}\n", nev);
  int k = 0;
  while (k < nev)
  {
    Mpi::Print("inside while k < nev loop, k: {}\n", k);
    // do we want set c inside the loop so it's random and thus avoids getting stuck on a "BAD" guess?
    // Set arbitrary c and solve for w0?
    linalg::SetRandom(GetComm(), c);
    std::srand((unsigned int) time(0)); // ?
    Eigen::VectorXcd c_eig = Eigen::VectorXcd::Random(k);
    Eigen::VectorXcd w0_eig;

    bool deflation = true;//false;
    std::complex<double> eig;
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
      linalg::SetRandom(GetComm(), v);
    }
    std::complex<double> eig_update = eig;

    Eigen::VectorXcd v2 = Eigen::VectorXcd::Constant(k, 0.0);// how to initialize v2? zero or random?
    if (k > 0 && deflation)
    {
      std::srand((unsigned int) time(0));
      v2.setRandom();
    }

    // Removed RII but we should test what works better when deflation is implemented
    // Quasi-Newton 2 from https://arxiv.org/pdf/1702.08492
    double norm_v = std::sqrt(linalg::Norml2(GetComm(), v, true) + v2.squaredNorm());
    v *= 1.0 / norm_v;
    v2 *= 1.0 / norm_v;
    norm_v = std::sqrt(linalg::Norml2(GetComm(), v, true) + v2.squaredNorm()); // should be 1...

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
    opInv->Mult(c, w0);
    // update w0_eig // this whole deflated solve should be a separate function since it's being called a few times!
    // deflated_solve(c_eig, w0, w0_eig, k, lambda) // uses X and H but those should be class data members
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
      ComplexVector phi2; phi2.SetSize(n); phi2.UseDevice(false); // keep on CPU?! Might want to implement methods in vector.cpp to do this
      phi2 = 0.0;
      auto *phi2r = phi2.Real().HostWrite();
      auto *phi2i = phi2.Imag().HostWrite();
      for (int j = 0; j < k; j++)
      {
        auto *XR = X[j].Real().Read();
        auto *XI = X[j].Imag().Read();
        for (int i = 0; i < n; i++)
        {
          phi2r[i] += phi1(j).real() * XR[i] - phi1(j).imag() * XI[i]; // transposedot
          phi2i[i] += phi1(j).imag() * XR[i] + phi1(j).real() * XI[i]; // transposedot
        }
      }
      linalg::AXPY(-1.0, phi2, w0);
    }

    int update_freq = 4;// SHOULD REVISIT THIS AND FIGURE OUT BEST FREQUENCY around 4-5 seems good?

    double min_res = 1e6;
    int min_it = 0;
    double init_res = 1e6;
    double res = 1e6;
    int large_res = 0;
    int it = 0;
    while (it < nleps_it)
    {
      Mpi::Print("inside while it < nleps_it lool it: {}\n", it);
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
        ComplexVector phi2; phi2.SetSize(n); phi2.UseDevice(false); // keep on CPU?! Might want to implement methods in vector.cpp to do this
        phi2 = 0.0;
        auto *phi2r = phi2.Real().HostWrite();
        auto *phi2i = phi2.Imag().HostWrite();
        for (int j = 0; j < k; j++)
        {
          auto *XR = X[j].Real().Read();
          auto *XI = X[j].Imag().Read();
          for (int i = 0; i < n; i++)
          {
            phi2r[i] += phi1(j).real() * XR[i] - phi1(j).imag() * XI[i];
            phi2i[i] += phi1(j).imag() * XR[i] + phi1(j).real() * XI[i];
          }
        }

        A->AddMult(phi2, u, 1.0);
        // u2 = A(l) v1 + B(l) v2 = sum i=0..p l_i (XH^i)* v1 + sum i=0..p (XH^i)* X qi(l) v2
        // if p = 1, A(l) = X.adjoint and B(l) = 0!
        u2.conservativeResize(k);
        for (int j = 0; j < k; j++)
        {
          u2(j) = linalg::Dot(GetComm(), v, X[j]);
        }
      }

      if (deflation)
      {
        res = std::sqrt(linalg::Norml2(GetComm(), u, true) + u2.squaredNorm()) / norm_v;
      }
      else
      {
        res = linalg::Norml2(GetComm(), u) / linalg::Norml2(GetComm(), v);
      }
      if (it == 0) init_res = res;
      if (res < min_res){min_res = res; min_it = it;}
      //min_res = std::min(res, min_res);
      Mpi::Print(GetComm(), "k: {}, it: {}, eig: {}, {}, res: {}\n", k, it, eig.real(), eig.imag(), res); // only print if print > 0? >1 ????

      if (res < rtol) // In addition to res < tol, do we also want to check if it's within the range (i.e. greater than the target)???
      {
        // Update the invariant pair with normalization.
        Mpi::Print("L463\n");
        const auto scale = linalg::Norml2(GetComm(), v);
        v *= 1.0 / scale;
        X.push_back(v);
        H.conservativeResizeLike(Eigen::MatrixXd::Zero(k + 1, k + 1));
        H.col(k).head(k) = v2 / scale;
        H(k, k) = eig;
        // Also store here?
        Mpi::Print("L471 eigenvalues.size(): {}\n", eigenvalues.size());
        eigenvalues[k] = eig;
        Mpi::Print("L473 eigenvectors.size(): {}\n", eigenvectors.size());
        eigenvectors[k] = v;
        k++; // only increment eigenpairs when a converged pair is found!
        break;
      }
      // Stop if large residual for 10 consecutive iterations
      large_res = (res > 0.9) ? large_res + 1 : 0;
      if (large_res > 10)
      {
        Mpi::Print(GetComm(), "k: {}, Quasi-Newton not converging after {} iterations, restarting\n", k, it);
        break;
      }

      // Compute w = J * v
      auto opA2p = space_op->GetExtraSystemMatrix<ComplexOperator>(std::abs(eig.imag()) * (1.0 + dl), Operator::DIAG_ZERO);
      std::complex<double> denom = dl * std::abs(eig.imag());
      auto opAJ = space_op->GetDividedDifferenceMatrix<ComplexOperator>(denom, opA2p.get(), A2n.get(), Operator::DIAG_ZERO);
      auto opJ = space_op->GetSystemMatrix(std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0), std::complex<double>(2.0, 0.0) * eig, opK, opC, opM, opAJ.get());
      opJ->Mult(v, w);
      // Create deflated Jac multiply function?
      // jacmult(w, v2, k, H) // assume X is a class data member? maybe H too?
      if (deflation)
      {
        // w1 = T'(l) v1 + U'(l) v2 = T'(l) v1 + T'(l)XS v2 - T(l)XS^2 v2
        const auto S = (eig * Eigen::MatrixXcd::Identity(k, k) - H).inverse(); // should re-use the one defined above!!
        const auto phi1 = S * v2;
        const auto phi2 = S * phi1;
        ComplexVector phi3; phi3.SetSize(n); phi3.UseDevice(false); // keep on CPU?! Might want to implement methods in vector.cpp to do this
        ComplexVector phi4; phi4.SetSize(n); phi4.UseDevice(false); // keep on CPU?! Might want to implement methods in vector.cpp to do this
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
          for (int i = 0; i < n; i++)
          {
            phi3r[i] += phi1(j).real() * XR[i] - phi1(j).imag() * XI[i];
            phi3i[i] += phi1(j).imag() * XR[i] + phi1(j).real() * XI[i];
            phi4r[i] += phi2(j).real() * XR[i] - phi2(j).imag() * XI[i];
            phi4i[i] += phi2(j).imag() * XR[i] + phi2(j).real() * XI[i];
          }
        }
        opJ->AddMult(phi3, w, 1.0);
        A->AddMult(phi4, w, -1.0);
        // with p = 1, A'(l) and B'(l) are zero so w2 is zero?
      }

      // Compute delta = - dot(w0, u) / dot(w0, w)
      std::complex<double> delta;
      std::complex<double> u2_w0(0.0, 0.0);
      if (deflation)
      {
        u2_w0 = std::complex<double>(w0_eig.adjoint() * u2);
      }
      delta = - (linalg::Dot(GetComm(), u, w0) + u2_w0) / linalg::Dot(GetComm(), w, w0);

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
        // use deflated solve function
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
          ComplexVector phi2; phi2.SetSize(n); phi2.UseDevice(false); // keep on CPU?! Might want to implement methods in vector.cpp to do this
          phi2 = 0.0;
          auto *phi2r = phi2.Real().HostWrite();
          auto *phi2i = phi2.Imag().HostWrite();
          for (int j = 0; j < k; j++)
          {
            auto *XR = X[j].Real().Read();
            auto *XI = X[j].Imag().Read();
            for (int i = 0; i < n; i++)
            {
              phi2r[i] += phi1(j).real() * XR[i] - phi1(j).imag() * XI[i];
              phi2i[i] += phi1(j).imag() * XR[i] + phi1(j).real() * XI[i];
            }
          }
          linalg::AXPY(-1.0, phi2, w0);
        }
      }
      opInv->Mult(z, u);
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
          u2(j) = linalg::Dot(GetComm(), u, X[j]);
        }
        z2 = z2 - u2;
        Eigen::MatrixXcd SS(k, k);
        for (int i = 0; i < k; i++)
        {
          for (int j = 0; j < k; j++)
          {
            SS(i,j) = linalg::Dot(GetComm(), X[i], X[j]); // 2 or 3??
          }
        }
        const auto S = (eig_update * Eigen::MatrixXcd::Identity(k, k) - H).inverse(); // use eig or something else?
        SS = - SS * S;
        u2 = SS.inverse() * z2;
        // u1 = u - X*S*u2
        const auto phi1 = S * u2;
        ComplexVector phi2; phi2.SetSize(n); phi2.UseDevice(false); // keep on CPU?! Might want to implement methods in vector.cpp to do this
        phi2 = 0.0;
        auto *phi2r = phi2.Real().HostWrite();
        auto *phi2i = phi2.Imag().HostWrite();
        for (int j = 0; j < k; j++)
        {
          auto *XR = X[j].Real().Read();
          auto *XI = X[j].Imag().Read();
          for (int i = 0; i < n; i++)
          {
            phi2r[i] += phi1(j).real() * XR[i] - phi1(j).imag() * XI[i]; // transposedot
            phi2i[i] += phi1(j).imag() * XR[i] + phi1(j).real() * XI[i]; // transposedot
          }
        }
        linalg::AXPY(-1.0, phi2, u);
        v2 += u2;
      }
      v += u;
      norm_v = std::sqrt(linalg::Norml2(GetComm(), v, true) + v2.squaredNorm());
      v *= 1.0 / norm_v;
      v2 *= 1.0 / norm_v;
      norm_v = std::sqrt(linalg::Norml2(GetComm(), v, true) + v2.squaredNorm()); // should be 1...

      it++;
      if (it == nleps_it) // ACTUALLY, WE DO NOT WANT TO SAVE UNCONVERGED EIGENPAIRS, CAN MESS UP FUTURE DEFLATION ITERATIONS
      {
        Mpi::Print("Quasi-Newton did not converge in {} iterations\n", nleps_it);
      }
    }
    Mpi::Print(GetComm(), "\n\n i: {}, init_res: {}, min_res: {}, min_it: {}\n\n", k, init_res, min_res, min_it);
  }
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

  Eigen::MatrixXcd Xeig = eps.eigenvectors();
  std::vector<Eigen::VectorXcd> sorted_Xeig;

  for (int k = 0; k < nev; k++)
  {
    perm[k] = order[k]; // stupid, only use one!
    sorted_Xeig.push_back(Xeig.col(order_eigen[k]));
  }

  for(int k = 0; k < nev; k++)
  {
    int k2 = order2[k];
    ComplexVector eigv; eigv.SetSize(n); eigv.UseDevice(false); // keep on CPU?! Might want to implement methods in vector.cpp to do this
    eigv = 0.0;
    auto *eigvr = eigv.Real().HostWrite();
    auto *eigvi = eigv.Imag().HostWrite();
    for (int j = 0; j < nev; j++)
    {
      auto *XR = X[j].Real().Read();
      auto *XI = X[j].Imag().Read();
      for (int i = 0; i < n; i++)
      {
        //eigvr[i] += Xeig(j,k).real() * XR[i] - Xeig(j,k).imag() * XI[i];
        //eigvi[i] += Xeig(j,k).imag() * XR[i] + Xeig(j,k).real() * XI[i];
        eigvr[i] += sorted_Xeig[k2](j).real() * XR[i] - sorted_Xeig[k2](j).imag() * XI[i];
        eigvi[i] += sorted_Xeig[k2](j).imag() * XR[i] + sorted_Xeig[k2](j).real() * XI[i];
      }
    }
    eigenvectors[k] = eigv;
  }


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
