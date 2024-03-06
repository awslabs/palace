// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "romoperator.hpp"

#include <algorithm>
#include <limits>
#include <numeric>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
#include <mfem.hpp>
#include "linalg/orthog.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

// Eigen does not provide a complex-valued genearlized eigenvalue solver, so we use LAPACK
// for this.
extern "C"
{
  void zggev_(char *, char *, int *, std::complex<double> *, int *, std::complex<double> *,
              int *, std::complex<double> *, std::complex<double> *, std::complex<double> *,
              int *, std::complex<double> *, int *, std::complex<double> *, int *, double *,
              int *);
}

namespace palace
{

using namespace std::complex_literals;

namespace
{

constexpr auto ORTHOG_TOL = 1.0e-12;

template <typename VecType, typename ScalarType>
inline void OrthogonalizeColumn(GmresSolverBase::OrthogType type, MPI_Comm comm,
                                const std::vector<VecType> &V, VecType &w, ScalarType *Rj,
                                int j)
{
  // Orthogonalize w against the leading j columns of V.
  switch (type)
  {
    case GmresSolverBase::OrthogType::MGS:
      linalg::OrthogonalizeColumnMGS(comm, V, w, Rj, j);
      break;
    case GmresSolverBase::OrthogType::CGS:
      linalg::OrthogonalizeColumnCGS(comm, V, w, Rj, j);
      break;
    case GmresSolverBase::OrthogType::CGS2:
      linalg::OrthogonalizeColumnCGS(comm, V, w, Rj, j, true);
      break;
  }
}

template <typename MatType>
inline void ProjectMatInternal(MPI_Comm comm, const std::vector<Vector> &V,
                               const ComplexOperator &A, MatType &Ar, ComplexVector &r,
                               int n0)
{
  // Update Ar = Vᴴ A V for the new basis dimension n0 -> n. V is real and thus the result
  // is complex symmetric if A is symmetric (which we assume is the case). Ar is replicated
  // across all processes as a sequential n x n matrix.
  const auto n = Ar.rows();
  MFEM_VERIFY(n0 < n, "Invalid dimensions in PROM matrix projection!");
  for (int j = n0; j < n; j++)
  {
    // Fill block of Vᴴ A V = [  | Vᴴ A vj ] . We can optimize the matrix-vector product
    // since the columns of V are real.
    MFEM_VERIFY(A.HasReal() || A.HasImag(),
                "Invalid zero ComplexOperator for PROM matrix projection!");
    if (A.HasReal())
    {
      A.Real()->Mult(V[j], r.Real());
    }
    if (A.HasImag())
    {
      A.Imag()->Mult(V[j], r.Imag());
    }
    for (int i = 0; i < n; i++)
    {
      Ar(i, j).real(A.HasReal() ? V[i] * r.Real() : 0.0);  // Local inner product
      Ar(i, j).imag(A.HasImag() ? V[i] * r.Imag() : 0.0);
    }
  }
  Mpi::GlobalSum((n - n0) * n, Ar.data() + n0 * n, comm);

  // Fill lower block of Vᴴ A V = [ ____________  |  ]
  //                              [ vjᴴ A V[1:n0] |  ] .
  for (int j = 0; j < n0; j++)
  {
    for (int i = n0; i < n; i++)
    {
      Ar(i, j) = Ar(j, i);
    }
  }
}

template <typename VecType>
inline void ProjectVecInternal(MPI_Comm comm, const std::vector<Vector> &V,
                               const ComplexVector &b, VecType &br, int n0)
{
  // Update br = Vᴴ b for the new basis dimension n0 -> n. br is replicated across all
  // processes as a sequential n-dimensional vector.
  const auto n = br.size();
  MFEM_VERIFY(n0 < n, "Invalid dimensions in PROM vector projection!");
  for (int i = n0; i < n; i++)
  {
    br(i).real(V[i] * b.Real());  // Local inner product
    br(i).imag(V[i] * b.Imag());
  }
  Mpi::GlobalSum(n - n0, br.data() + n0, comm);
}

template <typename MatType, typename VecType>
inline void ComputeMRI(const MatType &R, VecType &q)
{
  // Compute the coefficients of the minimal rational interpolation (MRI):
  // u = [sum_i u_i q_i / (z - z_i)] / [sum_i q_i / (z - z_i)]. The coefficients are given
  // by the right singular vector of R corresponding to the minimum singular value.
  const auto S = R.rows();
  MFEM_ASSERT(S > 0 && R.cols() == S, "Invalid dimension mismatch when computing MRI!");
  // For Eigen = v3.4.0 (latest tagged release as of 10/2023)
  Eigen::JacobiSVD<Eigen::MatrixXcd> svd;
  svd.compute(R, Eigen::ComputeFullV);
  // For Eigen > v3.4.0 (GitLab repo is at v3.4.90 as of 10/2023)
  // Eigen::JacobiSVD<Eigen::MatrixXcd, Eigen::ComputeFullV> svd;
  // svd.compute(R);
  const auto &sigma = svd.singularValues();
  auto m = S - 1;
  while (m > 0 && sigma[m] < ORTHOG_TOL * sigma[0])
  {
    Mpi::Warning("Minimal rational interpolation encountered rank-deficient matrix: "
                 "σ[{:d}] = {:.3e} (σ[0] = {:.3e})!\n",
                 m, sigma[m], sigma[0]);
    m--;
  }
  q = svd.matrixV().col(m);
}

template <typename MatType, typename VecType>
inline void ZGGEV(MatType &A, MatType &B, VecType &D, MatType &VR)
{
  // Wrapper for LAPACK's (z)ggev. A and B are overwritten by their Schur decompositions.
  MFEM_VERIFY(A.rows() == A.cols() && B.rows() == B.cols() && A.rows() == B.rows(),
              "Generalized eigenvalue problem expects A, B matrices to be square and have "
              "same dimensions!");
  char jobvl = 'N', jobvr = 'V';
  int n = static_cast<int>(A.rows()), lwork = 2 * n;
  std::vector<std::complex<double>> alpha(n), beta(n), work(lwork);
  std::vector<double> rwork(8 * n);
  MatType VL(0, 0);
  VR.resize(n, n);
  int info = 0;

  zggev_(&jobvl, &jobvr, &n, A.data(), &n, B.data(), &n, alpha.data(), beta.data(),
         VL.data(), &n, VR.data(), &n, work.data(), &lwork, rwork.data(), &info);
  MFEM_VERIFY(info == 0, "ZGGEV failed with info = " << info << "!");

  // Postprocess the eigenvalues and eigenvectors (return unit 2-norm eigenvectors).
  D.resize(n);
  for (int i = 0; i < n; i++)
  {
    D(i) = (beta[i] == 0.0)
               ? ((alpha[i] == 0.0) ? std::numeric_limits<std::complex<double>>::quiet_NaN()
                                    : mfem::infinity())
               : alpha[i] / beta[i];
    VR.col(i) /= VR.col(i).norm();
  }
}

template <typename VecType>
void ProlongatePROMSolution(std::size_t n, const std::vector<Vector> &V, const VecType &y,
                            ComplexVector &u)
{
  u = 0.0;
  for (std::size_t j = 0; j < n; j += 2)
  {
    if (j + 1 < n)
    {
      linalg::AXPBYPCZ(y(j).real(), V[j], y(j + 1).real(), V[j + 1], 1.0, u.Real());
      linalg::AXPBYPCZ(y(j).imag(), V[j], y(j + 1).imag(), V[j + 1], 1.0, u.Imag());
    }
    else
    {
      linalg::AXPY(y(j).real(), V[j], u.Real());
      linalg::AXPY(y(j).imag(), V[j], u.Imag());
    }
  }
}

}  // namespace

RomOperator::RomOperator(const IoData &iodata, SpaceOperator &spaceop, int max_size)
  : spaceop(spaceop)
{
  // Construct the system matrices defining the linear operator. PEC boundaries are handled
  // simply by setting diagonal entries of the system matrix for the corresponding dofs.
  // Because the Dirichlet BC is always homogenous, no special elimination is required on
  // the RHS. The damping matrix may be nullptr.
  K = spaceop.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
  C = spaceop.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  M = spaceop.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  MFEM_VERIFY(K && M, "Invalid empty HDM matrices when constructing PROM!");

  // Set up RHS vector (linear in frequency part) for the incident field at port boundaries,
  // and the vector for the solution, which satisfies the Dirichlet (PEC) BC.
  if (!spaceop.GetExcitationVector1(RHS1))
  {
    RHS1.SetSize(0);
  }
  has_A2 = has_RHS2 = true;

  // Initialize working vector storage.
  r.SetSize(K->Height());

  // Set up the linear solver and set operators but don't set the operators yet (this will
  // be done during an HDM solve at a given parameter point). The preconditioner for the
  // complex linear system is constructed from a real approximation to the complex system
  // matrix.
  ksp = std::make_unique<ComplexKspSolver>(iodata, spaceop.GetNDSpaces(),
                                           &spaceop.GetH1Spaces());

  // The initial PROM basis is empty. The provided maximum dimension is the number of sample
  // points (2 basis vectors per point). Basis orthogonalization method is configured using
  // GMRES/FGMRES settings.
  MFEM_VERIFY(max_size > 0, "Reduced order basis storage must have > 0 columns!");
  V.resize(2 * max_size, Vector());
  Q.resize(max_size, ComplexVector());
  dim_V = dim_Q = 0;
  switch (iodata.solver.linear.gs_orthog_type)
  {
    case config::LinearSolverData::OrthogType::MGS:
      orthog_type = GmresSolverBase::OrthogType::MGS;
      break;
    case config::LinearSolverData::OrthogType::CGS:
      orthog_type = GmresSolverBase::OrthogType::CGS;
      break;
    case config::LinearSolverData::OrthogType::CGS2:
      orthog_type = GmresSolverBase::OrthogType::CGS2;
      break;
  }
}

void RomOperator::SolveHDM(double omega, ComplexVector &u)
{
  // Compute HDM solution at the given frequency. The system matrix, A = K + iω C - ω² M +
  // A2(ω) is built by summing the underlying operator contributions.
  BlockTimer bt0(Timer::CONSTRUCT);
  auto A2 = spaceop.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO);
  has_A2 = (A2 != nullptr);
  auto A = spaceop.GetSystemMatrix(std::complex<double>(1.0, 0.0), 1i * omega,
                                   std::complex<double>(-omega * omega, 0.0), K.get(),
                                   C.get(), M.get(), A2.get());
  auto P =
      spaceop.GetPreconditionerMatrix<ComplexOperator>(1.0, omega, -omega * omega, omega);
  ksp->SetOperators(*A, *P);

  // The HDM excitation vector is computed as RHS = iω RHS1 + RHS2(ω).
  Mpi::Print("\n");
  if (has_RHS2)
  {
    has_RHS2 = spaceop.GetExcitationVector2(omega, r);
  }
  else
  {
    r = 0.0;
  }
  if (RHS1.Size())
  {
    r.Add(1i * omega, RHS1);
  }

  // Solve the linear system.
  BlockTimer bt1(Timer::SOLVE);
  ksp->Mult(r, u);
}

void RomOperator::UpdatePROM(double omega, ComplexVector &u)
{
  // Update V. The basis is always real (each complex solution adds two basis vectors if it
  // has a nonzero real and imaginary parts).
  MPI_Comm comm = spaceop.GetComm();
  const double normr = linalg::Norml2(comm, u.Real());
  const double normi = linalg::Norml2(comm, u.Imag());
  const bool has_real = (normr > ORTHOG_TOL * std::sqrt(normr * normr + normi * normi));
  const bool has_imag = (normi > ORTHOG_TOL * std::sqrt(normr * normr + normi * normi));
  MFEM_VERIFY(dim_V + has_real + has_imag <= V.size(),
              "Unable to increase basis storage size, increase maximum number of vectors!");
  const std::size_t dim_V0 = dim_V;
  std::vector<double> H(dim_V + has_real + has_imag);
  if (has_real)
  {
    V[dim_V] = u.Real();
    OrthogonalizeColumn(orthog_type, comm, V, V[dim_V], H.data(), dim_V);
    H[dim_V] = linalg::Norml2(comm, V[dim_V]);
    V[dim_V] *= 1.0 / H[dim_V];
    dim_V++;
  }
  if (has_imag)
  {
    V[dim_V] = u.Imag();
    OrthogonalizeColumn(orthog_type, comm, V, V[dim_V], H.data(), dim_V);
    H[dim_V] = linalg::Norml2(comm, V[dim_V]);
    V[dim_V] *= 1.0 / H[dim_V];
    dim_V++;
  }

  // Update reduced-order operators. Resize preserves the upper dim0 x dim0 block of each
  // matrix and first dim0 entries of each vector and the projection uses the values
  // computed for the unchanged basis vectors.
  Kr.conservativeResize(dim_V, dim_V);
  ProjectMatInternal(comm, V, *K, Kr, r, dim_V0);
  if (C)
  {
    Cr.conservativeResize(dim_V, dim_V);
    ProjectMatInternal(comm, V, *C, Cr, r, dim_V0);
  }
  Mr.conservativeResize(dim_V, dim_V);
  ProjectMatInternal(comm, V, *M, Mr, r, dim_V0);
  Ar.resize(dim_V, dim_V);
  if (RHS1.Size())
  {
    RHS1r.conservativeResize(dim_V);
    ProjectVecInternal(comm, V, RHS1, RHS1r, dim_V0);
  }
  RHSr.resize(dim_V);

  // Compute the coefficients for the minimal rational interpolation of the state u used
  // as an error indicator. The complex-valued snapshot matrix U = [{u_i}] is stored by its
  // QR decomposition.
  MFEM_VERIFY(dim_Q + 1 <= Q.size(),
              "Unable to increase basis storage size, increase maximum number of vectors!");
  R.conservativeResizeLike(Eigen::MatrixXd::Zero(dim_Q + 1, dim_Q + 1));
  Q[dim_Q] = u;
  OrthogonalizeColumn(orthog_type, comm, Q, Q[dim_Q], R.col(dim_Q).data(), dim_Q);
  R(dim_Q, dim_Q) = linalg::Norml2(comm, Q[dim_Q]);
  Q[dim_Q] *= 1.0 / R(dim_Q, dim_Q);
  dim_Q++;
  ComputeMRI(R, q);
  // if (Mpi::Root(comm))
  // {
  //   std::cout << "MRI (S = " << dim_Q << "):\n"
  //   std::cout << "R =\n" << R << "\n";
  //   std::cout << "q =\n" << q << "\n";
  // }
  z.push_back(omega);
}

void RomOperator::SolvePROM(double omega, ComplexVector &u)
{
  // Assemble the PROM linear system at the given frequency. The PROM system is defined by
  // the matrix Aᵣ(ω) = Kᵣ + iω Cᵣ - ω² Mᵣ + Vᴴ A2 V(ω) and source vector RHSᵣ(ω) = iω RHS1ᵣ
  // + Vᴴ RHS2(ω). A2(ω) and RHS2(ω) are constructed only if required and are only nonzero
  // on boundaries, will be empty if not needed.
  if (has_A2)
  {
    auto A2 = spaceop.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO);
    ProjectMatInternal(spaceop.GetComm(), V, *A2, Ar, r, 0);
  }
  else
  {
    Ar.setZero();
  }
  Ar += Kr;
  if (Cr.rows() > 0)
  {
    Ar += (1i * omega) * Cr;
  }
  Ar += (-omega * omega) * Mr;

  if (has_RHS2)
  {
    spaceop.GetExcitationVector2(omega, r);
    ProjectVecInternal(spaceop.GetComm(), V, r, RHSr, 0);
  }
  else
  {
    RHSr.setZero();
  }
  if (RHS1.Size())
  {
    RHSr += (1i * omega) * RHS1r;
  }

  // Compute PROM solution at the given frequency and expand into high-dimensional space.
  // The PROM is solved on every process so the matrix-vector product for vector expansion
  // does not require communication.
  if constexpr (false)
  {
    // LDLT solve
    RHSr = Ar.ldlt().solve(RHSr);
    // RHSr = Ar.selfadjointView<Eigen::Lower>().ldlt().solve(RHSr);
  }
  else
  {
    // LU solve
    RHSr = Ar.partialPivLu().solve(RHSr);
  }
  ProlongatePROMSolution(dim_V, V, RHSr, u);
}

double RomOperator::FindMaxError(double delta) const
{
  // Return an estimate for argmax_z ||u(z) - V y(z)|| as argmin_z |Q(z)| with Q(z) =
  // sum_i q_z / (z - z_i) (denominator of the barycentric interpolation of u). The roots of
  // Q are given analytically as the solution to an S + 1 dimensional eigenvalue problem.
  const auto S = dim_Q;
  MFEM_VERIFY(S >= 2, "Maximum error can only be found once two sample points have been "
                      "added to the PROM to define the parameter domain!");
  double start = *std::min_element(z.begin(), z.end());
  double end = *std::max_element(z.begin(), z.end());
  Eigen::Map<const Eigen::VectorXd> z_map(z.data(), S);

  Eigen::MatrixXcd A = Eigen::MatrixXcd::Zero(S + 1, S + 1);
  A.col(0).tail(S) = Eigen::VectorXcd::Ones(S);
  A.row(0).tail(S) = q;
  A.diagonal().tail(S) = z_map.array();

  Eigen::MatrixXcd B = Eigen::MatrixXcd::Identity(S + 1, S + 1);
  B(0, 0) = 0.0;

  Eigen::VectorXcd D;
  Eigen::MatrixXcd X;
  ZGGEV(A, B, D, X);

  // If there are multiple roots in [start, end], pick the one furthest from the existing
  // set of samples.
  double dist_star = 0.0, z_star = 0.0;
  for (auto d : D)
  {
    if (std::real(d) >= start && std::real(d) <= end)
    {
      const double dist = (z_map.array() - std::real(d)).abs().maxCoeff();
      if (dist > dist_star)
      {
        z_star = std::real(d);
        dist_star = dist;
      }
    }
  }
  if (z_star > 0.0)
  {
    return z_star;
  }

  // XX TODO DEBUG: FALLBACK TO SAMPLING ON GRID WITH NO POLES
  Mpi::Print("\nFalling back to sampling z* on grid [{:.3e}, {:.3e}]\n", start, end);

  // Fall back to sampling Q on discrete points if no roots exist in [start, end].
  double Q_star = mfem::infinity();
  if (delta < 0.0)
  {
    delta = -delta;
  }
  while (start <= end)
  {
    const double Q = std::abs((q.array() / (z_map.array() - start)).sum());
    if (Q < Q_star)
    {
      z_star = start;
      Q_star = Q;
    }
    start += delta;
  }
  MFEM_VERIFY(z_star > 0.0, "Could not locate a maximum error in the range ["
                                << start << ", " << end << "]!");
  return z_star;
}

namespace
{

template <typename F, typename VecType>
void MSLP(int n, F EvalFunction, std::complex<double> &lambda, VecType &x,
          double tol = 1.0e-9, int max_it = 100)
{
  MFEM_VERIFY(x.size() == n,
              "Must provide an initial guess for the eigenvector x in MSLP solver!");

  using MatType = Eigen::MatrixXcd;
  MatType T(n, n), dT(n, n), X(n, n);
  VecType r(n), mu(n);
  Eigen::ComplexEigenSolver<MatType> eps;

  int it = 0;
  while (it < max_it)
  {
    // Check convergence.
    EvalFunction(lambda, T, dT, true, (it == 0));
    r = T * x;

    // // XX TODO DEBUG WIP
    // Mpi::Print("MSLP iteration {:d}, l = {:e}{:+e}i, ||r|| = {:e}, ||T|| = {:e}\n", it,
    //            lambda.real(), lambda.imag(), r.norm(), T.norm());

    double res = r.norm() / (T.norm() * x.norm());
    if (res < tol)
    {
      break;
    }

    // Set up and solve the linear EVP.
    ZGGEV(T, dT, mu, X);

    // Update eigenpair estimates.
    const auto i =
        std::distance(mu.begin(), std::min_element(mu.begin(), mu.end(),
                                                   [](auto l, auto r)
                                                   { return std::abs(l) < std::abs(r); }));
    lambda -= mu(i);
    x = X.col(i);
    it++;
  }
  // if (it == max_it)
  // {
  //   EvalFunction(lambda, T, dT, false, false);
  //   r = T * x;
  //   Mpi::Warning(
  //       "MSLP solver did not converge, ||Tx|| / ||T|| ||x|| = {:.3e} (tol = {:.3e})!\n",
  //       r.norm() / (T.norm() * x.norm()), tol);
  // }
}

template <typename F, typename MatType, typename VecType>
void SolveNEP(int n, int num_eig, std::complex<double> sigma, F EvalFunction, VecType &D,
              MatType &X)
{
  // This is the deflation scheme used by SLEPc's NEP solver.
  // Reference: Effenberger, Robust successive computation of eigenpairs for nonlinear
  //                eigenvalue problems, SIAM J. Matrix Anal. Appl. (2013).
  MatType T(n, n), dT(n, n), H;
  for (int k = 0; k < num_eig; k++)
  {
    // Precompute all required powers of H. The minimality index of (X, H) is 1 if
    // the eigenvectors are linearly independent.
    int p = 1;
    if constexpr (false)
    {
      if (k > 0)
      {
        MatType Vl = X;
        Eigen::ColPivHouseholderQR<MatType> qr;
        qr.compute(Vl);
        while (qr.rank() < k)
        {
          Vl.conservativeResize((p + 1) * n, k);
          Vl.bottomRows(n) = Vl.block((p - 1) * n, 0, n, k) * H;
          qr.compute(Vl);
          p++;
        }
      }
    }
    std::vector<MatType> HH(p);
    HH[0] = MatType::Identity(k, k);
    for (int i = 1; i < p; i++)
    {
      HH[i] = HH[i - 1] * H;
    }

    auto EvalDeflated =
        [&](std::complex<double> l, MatType &Tp, MatType &dTp, bool jacobian, bool first)
    {
      if (k == 0)
      {
        EvalFunction(l, Tp, dTp, jacobian, first);
        return;
      }

      // Compute the extended operators of the deflated problem, with explicit computation
      // of S = (λ I - H)⁻¹ for U(λ) and U'(λ). When constructing S, for some reason the
      // matrix inverse works more nicely than triangularView<Eigen::Upper>().solve().
      EvalFunction(l, T, dT, jacobian, first);
      const auto S = (l * MatType::Identity(k, k) - H).inverse();
      Tp.topLeftCorner(n, n) = T;
      Tp.topRightCorner(n, k) = T * (X * S);

      // Second row: A(λ) x + B(λ) t = 0.
      Tp.bottomRows(k) = MatType::Zero(k, n + k);
      MatType XH(n, k);
      for (int i = 0; i < p; i++)
      {
        XH = (X * HH[i]).adjoint();
        Tp.bottomLeftCorner(k, n) += std::pow(l, i) * XH;
        if (i > 0)
        {
          MatType qi = MatType::Zero(k, k);
          for (int j = 0; j <= i - 1; j++)
          {
            qi += std::pow(l, j) * HH[i - j - 1];
          }
          Tp.bottomRightCorner(k, k) += XH * (X * qi);
        }
      }

      if (jacobian)
      {
        dTp.topLeftCorner(n, n) = dT;
        dTp.topRightCorner(n, k) = dT * X * S;
        dTp.topRightCorner(n, k) -= Tp.topRightCorner(n, k) * S;
        dTp.bottomRows(k) = MatType::Zero(k, n + k);
        for (int i = 1; i < p; i++)
        {
          XH = (X * HH[i]).adjoint();
          dTp.bottomLeftCorner(k, n) += (std::pow(l, i - 1) * (double)i) * XH;
          if (i > 1)
          {
            MatType qi = MatType::Zero(k, k);
            for (int j = 1; j <= i - 1; j++)
            {
              qi += (std::pow(l, j - 1) * (double)j) * HH[i - j - 1];
            }
            dTp.bottomRightCorner(k, k) += XH * (X * qi);
          }
        }
      }
    };

    // Solve the deflated NEP with initial guess σ.
    auto lambda = sigma;
    VecType x = VecType::Random(n + k);
    MSLP(n + k, EvalDeflated, lambda, x);

    // XX TODO DEBUG WIP
    Mpi::Print("Eigenvalue {:d}/{:d}, l = {:e}{:+e}i\n", k + 1, num_eig, lambda.real(),
               lambda.imag());

    // Update the invariant pair with normalization. This seems to work better than taking
    // the pair (X, H) and normalizing it via a QR decomposition of X.
    X.conservativeResize(n, k + 1);
    H.conservativeResizeLike(Eigen::MatrixXd::Zero(k + 1, k + 1));
    const auto scale = x.head(n).norm();
    X.col(k) = x.head(n) / scale;
    H.col(k).head(k) = x.tail(k) / scale;
    H(k, k) = lambda;
  }

  // Eigenpair extraction from the invariant pair (X, H).
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eps;
  eps.compute(H);
  D = eps.eigenvalues();
  X *= eps.eigenvectors();
  X *= X.colwise().norm().cwiseInverse().asDiagonal();
}

template <typename F, typename MatType, typename VecType>
void SolveNEP2(int n, int num_eig, std::complex<double> sigma, F EvalFunction, VecType &D,
               MatType &X)
{
  // Deflate the converged eigenvectors explicitly (works as long as eigenvectors are
  // independent).
  MatType QQ;
  D.resize(num_eig);
  X.resize(n, num_eig);
  for (int k = 0; k < num_eig; k++)
  {
    auto EvalDeflated =
        [&](std::complex<double> l, MatType &T, MatType &dT, bool jacobian, bool first)
    {
      EvalFunction(l, T, dT, jacobian, first);
      if (k == 0)
      {
        return;
      }
      const double penalty = 1.0e3 * T.norm();
      T += penalty * QQ;
    };

    // Solve the deflated NEP with initial guess σ.
    auto lambda = sigma;
    VecType x = VecType::Random(n);
    MSLP(n, EvalDeflated, lambda, x);
    D(k) = lambda;
    X.col(k) = x;

    // XX TODO DEBUG WIP
    Mpi::Print("Eigenvalue {:d}/{:d}, l = {:e}{:+e}i\n", k + 1, num_eig, D(k).real(),
               D(k).imag());

    // Update basis for deflation (Eigen QR returns the full Q).
    Eigen::ColPivHouseholderQR<MatType> qr(X.leftCols(k + 1));
    const MatType Q = qr.householderQ() * MatType::Identity(n, k + 1);
    QQ = Q * Q.adjoint();
  }
}

}  // namespace

std::vector<std::complex<double>> RomOperator::ComputeEigenvalueEstimates()
{

  // XX TODO ADD EIGENVECTORS TO COMPUTE HDM RESIDUAL ERROR?

  MFEM_VERIFY(dim_V > 0,
              "Eigenvalue estimates are only available for a PROM with nonzero dimension!");
  Eigen::VectorXcd omega;
  Eigen::MatrixXcd X;
  if (!has_A2)
  {
    if (Cr.rows() == 0)
    {
      // Linear generalized EVP: M⁻¹ K x = μ x. Linear EVP has eigenvalue μ = -λ² = ω².
      Eigen::MatrixXcd cKr(Kr), cMr(Mr);
      ZGGEV(cKr, cMr, omega, X);
      for (std::size_t i = 0; i < dim_V; i++)
      {
        omega(i) = std::sqrt(omega(i));
      }
    }
    else
    {
      // Quadratic EVP: P(λ) x = (K + λ C + λ² M) x = 0, solved via linearization.
      Eigen::MatrixXcd L0 = Eigen::MatrixXcd::Zero(2 * dim_V, 2 * dim_V);
      L0.topRightCorner(dim_V, dim_V) = Eigen::MatrixXcd::Identity(dim_V, dim_V);
      L0.bottomLeftCorner(dim_V, dim_V) = -Kr;
      L0.bottomRightCorner(dim_V, dim_V) = -Cr;
      Eigen::MatrixXcd L1 = Eigen::MatrixXcd::Zero(2 * dim_V, 2 * dim_V);
      L1.topLeftCorner(dim_V, dim_V) = Eigen::MatrixXcd::Identity(dim_V, dim_V);
      L1.bottomRightCorner(dim_V, dim_V) = Mr;
      ZGGEV(L0, L1, omega, X);
      for (std::size_t i = 0; i < 2 * dim_V; i++)
      {
        // if (omega(i).imag() >= -0.1 * std::abs(omega(i).real()))
        {
          omega(i) /= 1i;
        }
        // else
        // {

        //   // XX TODO: For dynamics, do we really want to do this?

        //   // Ignore eigenpairs outside the desired range.
        //   omega(i) = mfem::infinity();
        // }
      }
    }
  }
  else
  {
    // General nonlinear EVP: T(λ) x = (K + λ C + λ² M + A2(Im{λ})) x = 0. If C != 0, the
    // problem is at least quadratic. The all processes solve the eigenvalue problem
    // together.
    std::complex<double> l_prev;
    Eigen::MatrixXcd Ar_prev;
    auto EvalFunction = [&l_prev, &Ar_prev, this](std::complex<double> l,
                                                  Eigen::MatrixXcd &T, Eigen::MatrixXcd &dT,
                                                  bool jacobian, bool first)
    {
      // Assemble T(λ) = K + λ C + λ² M + A2(Im{λ}) and T'(λ) = C + 2 λ M + A2'(Im{λ}) .
      if (has_A2)
      {
        auto A2 = spaceop.GetExtraSystemMatrix<ComplexOperator>(std::abs(l.imag()),
                                                                Operator::DIAG_ZERO);
        ProjectMatInternal(spaceop.GetComm(), V, *A2, Ar, r, 0);
        T = Ar;
      }
      else
      {
        T.setZero();
      }
      T += Kr;
      if (Cr.rows() > 0)
      {
        T += l * Cr;
      }
      T += (l * l) * Mr;
      if (jacobian)
      {
        if (has_A2)
        {
          // Evaluate A2' by finite differencing, and reuse the existing evaluation from the
          // residual evaluation.
          if (!first && Ar_prev.rows() > 0)
          {
            dT = Ar_prev;
            dT -= Ar;
            dT * -1.0 / (l_prev - l);
          }
          else
          {
            const auto eps = std::sqrt(std::numeric_limits<double>::epsilon());
            auto A2 = spaceop.GetExtraSystemMatrix<ComplexOperator>(
                std::abs(l.imag()) * (1.0 + eps), Operator::DIAG_ZERO);
            ProjectMatInternal(spaceop.GetComm(), V, *A2, dT, r, 0);
            dT -= Ar;
            dT *= 1.0 / (eps * std::abs(l.imag()));
          }
        }
        else
        {
          dT.setZero();
        }
        if (Cr.rows() > 0)
        {
          dT += Cr;
        }
        dT += (2.0 * l) * Mr;
      }
      if (has_A2)
      {
        l_prev = l;
        Ar_prev = Ar;
      }
    };

    // For the initial guess, we choose the second smallest entry of the sample points. We
    // could alternatively use the center of the sample range.
    double sigma = *std::min_element(z.begin(), z.end());
    // sigma = *std::min_element(z.begin(), z.end(),
    //                           [sigma](auto l, auto r)
    //                           {
    //                             if (l == sigma)
    //                               return false;
    //                             if (r == sigma)
    //                               return true;
    //                             return l < r;
    //                           });
    sigma = 0.5 * (sigma + *std::max_element(z.begin(), z.end()));
    const std::size_t num_eig = dim_V;  // XX TODO SPECIFY AT RUNTIME WITH TOLERANCES
    if constexpr (false)
    {
      // Variant with explicit deflation of eigenbasis.
      SolveNEP2(dim_V, num_eig, 1i * sigma, EvalFunction, omega, X);
    }
    else
    {
      // Variant with extended problem deflation from Effenberger.
      SolveNEP(dim_V, num_eig, 1i * sigma, EvalFunction, omega, X);
    }
    for (std::size_t i = 0; i < num_eig; i++)
    {
      // if (omega(i).imag() >= -0.1 * std::abs(omega(i).real()))
      {
        omega(i) /= 1i;
      }
      // else
      // {

      //   // XX TODO: For dynamics, do we really want to do this?

      //   // Ignore eigenpairs outside the desired range.
      //   omega(i) = mfem::infinity();
      // }
    }
  }

  // XX TODO: Evaluate A2(ω) in the complex plane (upgrade models)

  // Compute HDM eigenvectors and residual norms.
  std::vector<std::size_t> perm(omega.size());
  std::iota(perm.begin(), perm.end(), 0);
  std::sort(perm.begin(), perm.end(),
            [&omega](auto l, auto r) { return std::abs(omega(l)) < std::abs(omega(r)); });


  //XX TODO DEBUG
  Mpi::Print("\n\nEigenpair residuals:\n");


  ComplexVector u(r.Size());
  for (int i = 0; i < omega.size(); i++)
  {
    // Compute HDM eigenmode.
    ProlongatePROMSolution(dim_V, V, X.col(i), u);
    linalg::Normalize(spaceop.GetComm(), u);

    // Evaluate the HDM eigenpair residual.
    std::unique_ptr<ComplexOperator> A2;
    if (has_A2)
    {
      A2 = spaceop.GetExtraSystemMatrix<ComplexOperator>(omega(i).real(),
                                                         Operator::DIAG_ZERO);
    }
    auto A =
        spaceop.GetSystemMatrix(std::complex<double>(1.0, 0.0), 1i * omega(i),
                                -omega(i) * omega(i), K.get(), C.get(), M.get(), A2.get());
    A->Mult(u, r);
    double res = linalg::Norml2(spaceop.GetComm(), r);


    //XX TODO DEBUG
    Mpi::Print("omega = {:e}{:+e}, ||r|| = {:e}\n",
               omega(i).real(), omega(i).imag(), res);


    // XX TODO WIP.... (scaling, storage, etc.)



  }

  // XX TODO FILTERING...

  return {omega.begin(), omega.end()};
}

}  // namespace palace
