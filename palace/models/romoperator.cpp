// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "romoperator.hpp"

#include <Eigen/SVD>
#include <mfem.hpp>
#include "linalg/orthog.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

// Eigen does not provide a complex-valued generalized eigenvalue solver, so we use LAPACK
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

inline void ProjectMatInternal(MPI_Comm comm, const std::vector<Vector> &V,
                               const ComplexOperator &A, Eigen::MatrixXcd &Ar,
                               ComplexVector &r, int n0)
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
    MFEM_VERIFY(A.Real() || A.Imag(),
                "Invalid zero ComplexOperator for PROM matrix projection!");
    if (A.Real())
    {
      A.Real()->Mult(V[j], r.Real());
    }
    if (A.Imag())
    {
      A.Imag()->Mult(V[j], r.Imag());
    }
    for (int i = 0; i < n; i++)
    {
      Ar(i, j).real(A.Real() ? V[i] * r.Real() : 0.0);  // Local inner product
      Ar(i, j).imag(A.Imag() ? V[i] * r.Imag() : 0.0);
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

inline void ProjectVecInternal(MPI_Comm comm, const std::vector<Vector> &V,
                               const ComplexVector &b, Eigen::VectorXcd &br, int n0)
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

inline void ComputeMRI(const Eigen::MatrixXcd &R, Eigen::VectorXcd &q)
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
inline void ProlongatePROMSolution(std::size_t n, const std::vector<Vector> &V,
                                   const VecType &y, ComplexVector &u)
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

MinimalRationalInterpolation::MinimalRationalInterpolation(int max_size)
{
  Q.resize(max_size, ComplexVector());
}

void MinimalRationalInterpolation::AddSolutionSample(
    double omega, const ComplexVector &u, const SpaceOperator &space_op,
    GmresSolverBase::OrthogType orthog_type)
{
  MPI_Comm comm = space_op.GetComm();

  // Compute the coefficients for the minimal rational interpolation of the state u used
  // as an error indicator. The complex-valued snapshot matrix U = [{u_i, (iω) u_i}] is
  // stored by its QR decomposition.
  MFEM_VERIFY(dim_Q + 1 <= Q.size(),
              "Unable to increase basis storage size, increase maximum number of vectors!");
  R.conservativeResizeLike(Eigen::MatrixXd::Zero(dim_Q + 1, dim_Q + 1));
  {
    std::vector<const ComplexVector *> blocks = {&u, &u};
    std::vector<std::complex<double>> s = {1.0, 1i * omega};
    Q[dim_Q].SetSize(2 * u.Size());
    Q[dim_Q].UseDevice(true);
    Q[dim_Q].SetBlocks(blocks, s);
  }
  OrthogonalizeColumn(orthog_type, comm, Q, Q[dim_Q], R.col(dim_Q).data(), dim_Q);
  R(dim_Q, dim_Q) = linalg::Norml2(comm, Q[dim_Q]);
  Q[dim_Q] *= 1.0 / R(dim_Q, dim_Q);
  dim_Q++;
  ComputeMRI(R, q);
  if constexpr (false)
  {
    Mpi::Print("MRI (S = {}):\nR = {}\nq = {}", dim_Q, R, q);
  }
  z.push_back(omega);
}

std::vector<double> MinimalRationalInterpolation::FindMaxError(int N) const
{
  // Return an estimate for argmax_z ||u(z) - V y(z)|| as argmin_z |Q(z)| with Q(z) =
  // sum_i q_z / (z - z_i) (denominator of the barycentric interpolation of u). The roots of
  // Q are given analytically as the solution to an S + 1 dimensional eigenvalue problem.
  BlockTimer bt(Timer::CONSTRUCT_PROM);
  const auto S = dim_Q;
  MFEM_VERIFY(S >= 2, "Maximum error can only be found once two sample points have been "
                      "added to the PROM to define the parameter domain!");
  double start = *std::min_element(z.begin(), z.end());
  double end = *std::max_element(z.begin(), z.end());
  Eigen::Map<const Eigen::VectorXd> z_map(z.data(), S);
  std::vector<std::complex<double>> z_star(N, 0.0);

  // XX TODO: For now, we explicitly minimize Q on the real line since we don't allow
  //          samples at complex-valued points (yet).

  // Eigen::MatrixXcd A = Eigen::MatrixXcd::Zero(S + 1, S + 1);
  // A.diagonal().head(S) = z_map.array();
  // A.row(S).head(S) = q;
  // A.col(S).head(S) = Eigen::VectorXcd::Ones(S);

  // Eigen::MatrixXcd B = Eigen::MatrixXcd::Identity(S + 1, S + 1);
  // B(S, S) = 0.0;

  // Eigen::VectorXcd D;
  // Eigen::MatrixXcd X;
  // ZGGEV(A, B, D, X);

  // // If there are multiple roots in [start, end], pick the ones furthest from the
  // // existing set of samples.
  // {
  //   std::vector<double> dist_star(N, 0.0);
  //   for (auto d : D)
  //   {
  //     if (std::real(d) < start || std::real(d) > end)
  //     {
  //       continue;
  //     }
  //     const double dist = (z_map.array() - std::real(d)).abs().maxCoeff();
  //     for (int i = 0; i < N; i++)
  //     {
  //       if (dist > dist_star[i])
  //       {
  //         for (int j = i + 1; j < N; j++)
  //         {
  //           z_star[j] = z_star[j - 1];
  //           dist_star[j] = dist_star[j - 1];
  //         }
  //         z_star[i] = start;
  //         dist_star[i] = dist;
  //       }
  //     }
  //   }
  // }

  // Fall back to sampling Q on discrete points if no roots exist in [start, end].
  if (std::abs(z_star[0]) == 0.0)
  {
    const auto delta = (end - start) / 1.0e6;
    std::vector<double> Q_star(N, mfem::infinity());
    while (start <= end)
    {
      const double Q = std::abs((q.array() / (z_map.array() - start)).sum());
      for (int i = 0; i < N; i++)
      {
        if (Q < Q_star[i])
        {
          for (int j = i + 1; j < N; j++)
          {
            z_star[j] = z_star[j - 1];
            Q_star[j] = Q_star[j - 1];
          }
          z_star[i] = start;
          Q_star[i] = Q;
        }
      }
      start += delta;
    }
    MFEM_VERIFY(
        N == 0 || std::abs(z_star[0]) > 0.0,
        fmt::format("Could not locate a maximum error in the range [{}, {}]!", start, end));
  }
  std::vector<double> vals(z_star.size());
  std::transform(z_star.begin(), z_star.end(), vals.begin(),
                 [](std::complex<double> z) { return std::real(z); });
  return vals;
}

RomOperator::RomOperator(const IoData &iodata, SpaceOperator &space_op,
                         int max_size_per_excitation)
  : space_op(space_op)
{
  // Construct the system matrices defining the linear operator. PEC boundaries are handled
  // simply by setting diagonal entries of the system matrix for the corresponding dofs.
  // Because the Dirichlet BC is always homogeneous, no special elimination is required on
  // the RHS. The damping matrix may be nullptr.
  K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
  C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  MFEM_VERIFY(K && M, "Invalid empty HDM matrices when constructing PROM!");

  // Initialize working vector storage.
  r.SetSize(K->Height());
  r.UseDevice(true);

  // Set up the linear solver and set operators but don't set the operators yet (this will
  // be done during an HDM solve at a given parameter point). The preconditioner for the
  // complex linear system is constructed from a real approximation to the complex system
  // matrix.
  ksp = std::make_unique<ComplexKspSolver>(iodata, space_op.GetNDSpaces(),
                                           &space_op.GetH1Spaces());

  auto excitation_helper = space_op.GetPortExcitations();

  // The initial PROM basis is empty. The provided maximum dimension is the number of sample
  // points (2 basis vectors per point). Basis orthogonalization method is configured using
  // GMRES/FGMRES settings.
  MFEM_VERIFY(max_size_per_excitation * excitation_helper.Size() > 0,
              "Reduced order basis storage must have > 0 columns!");
  V.resize(2 * max_size_per_excitation * excitation_helper.Size(), Vector());
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
  // Set up MRI.
  for (const auto &[excitation_idx, data] : excitation_helper)
  {
    mri.emplace(excitation_idx, MinimalRationalInterpolation(max_size_per_excitation));
  }
}

void RomOperator::SetExcitationIndex(int excitation_idx)
{
  // Set up RHS vector (linear in frequency part) for the incident field at port boundaries,
  // and the vector for the solution, which satisfies the Dirichlet (PEC) BC.
  excitation_idx_cache = excitation_idx;
  has_RHS1 = space_op.GetExcitationVector1(excitation_idx_cache, RHS1);
  if (!has_RHS1)
  {
    RHS1.SetSize(0);
  }
  else
  {
    // Project RHS1 to RHS1r with current PROM.
    if (dim_V > 0)
    {
      auto comm = space_op.GetComm();
      RHS1r.conservativeResize(dim_V);
      ProjectVecInternal(comm, V, RHS1, RHS1r, 0);
    }
  }
}

void RomOperator::SolveHDM(int excitation_idx, double omega, ComplexVector &u)
{
  if (excitation_idx_cache != excitation_idx)
  {
    SetExcitationIndex(excitation_idx);
  }
  // Compute HDM solution at the given frequency. The system matrix, A = K + iω C - ω² M +
  // A2(ω) is built by summing the underlying operator contributions.
  A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO);
  has_A2 = (A2 != nullptr);
  auto A = space_op.GetSystemMatrix(std::complex<double>(1.0, 0.0), 1i * omega,
                                    std::complex<double>(-omega * omega, 0.0), K.get(),
                                    C.get(), M.get(), A2.get());
  auto P =
      space_op.GetPreconditionerMatrix<ComplexOperator>(1.0, omega, -omega * omega, omega);
  ksp->SetOperators(*A, *P);

  // The HDM excitation vector is computed as RHS = iω RHS1 + RHS2(ω).
  Mpi::Print("\n");
  if (has_RHS2)
  {
    has_RHS2 = space_op.GetExcitationVector2(excitation_idx, omega, r);
  }
  else
  {
    r = 0.0;
  }
  if (has_RHS1)
  {
    r.Add(1i * omega, RHS1);
  }

  // Solve the linear system.
  ksp->Mult(r, u);
}

void RomOperator::UpdatePROM(const ComplexVector &u)
{

  // Update V. The basis is always real (each complex solution adds two basis vectors if it
  // has a nonzero real and imaginary parts).
  BlockTimer bt(Timer::CONSTRUCT_PROM);
  MPI_Comm comm = space_op.GetComm();
  const double normr = linalg::Norml2(comm, u.Real());
  const double normi = linalg::Norml2(comm, u.Imag());
  const bool has_real = (normr > ORTHOG_TOL * std::sqrt(normr * normr + normi * normi));
  const bool has_imag = (normi > ORTHOG_TOL * std::sqrt(normr * normr + normi * normi));
  MFEM_VERIFY(dim_V + has_real + has_imag <= V.size(),
              "Unable to increase basis storage size, increase maximum number of vectors!");
  const std::size_t dim_V0 = dim_V;
  std::vector<double> H(dim_V + static_cast<std::size_t>(has_real) +
                        static_cast<std::size_t>(has_imag));
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
}

void RomOperator::UpdateMRI(int excitation_idx, double omega, const ComplexVector &u)
{
  BlockTimer bt(Timer::CONSTRUCT_PROM);
  mri.at(excitation_idx).AddSolutionSample(omega, u, space_op, orthog_type);
  z.push_back(omega);// test for eigs
}

void RomOperator::SolvePROM(int excitation_idx, double omega, ComplexVector &u)
{
  if (excitation_idx_cache != excitation_idx)
  {
    SetExcitationIndex(excitation_idx);
  }

  // Assemble the PROM linear system at the given frequency. The PROM system is defined by
  // the matrix Aᵣ(ω) = Kᵣ + iω Cᵣ - ω² Mᵣ + Vᴴ A2 V(ω) and source vector RHSᵣ(ω) =
  // iω RHS1ᵣ + Vᴴ RHS2(ω). A2(ω) and RHS2(ω) are constructed only if required and are
  // only nonzero on boundaries, will be empty if not needed.
  if (has_A2 && Ar.rows() > 0)
  {
    A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO);
    ProjectMatInternal(space_op.GetComm(), V, *A2, Ar, r, 0);
  }
  else
  {
    Ar.setZero();
  }
  Ar += Kr;
  if (C)
  {
    Ar += (1i * omega) * Cr;
  }
  Ar += (-omega * omega) * Mr;

  if (has_RHS2 && RHSr.size() > 0)
  {
    space_op.GetExcitationVector2(excitation_idx, omega, RHS2);
    ProjectVecInternal(space_op.GetComm(), V, RHS2, RHSr, 0);
  }
  else
  {
    RHSr.setZero();
  }
  if (has_RHS1)
  {
    RHSr += (1i * omega) * RHS1r;
  }

  // Compute PROM solution at the given frequency and expand into high-dimensional space.
  // The PROM is solved on every process so the matrix-vector product for vector expansion
  // does not require communication.
  BlockTimer bt(Timer::SOLVE_PROM);
  if constexpr (false)
  {
    // LDLT solve.
    RHSr = Ar.ldlt().solve(RHSr);
    RHSr = Ar.selfadjointView<Eigen::Lower>().ldlt().solve(RHSr);
  }
  else
  {
    // LU solve.
    RHSr = Ar.partialPivLu().solve(RHSr);
  }
  ProlongatePROMSolution(dim_V, V, RHSr, u);
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
  Mpi::Print("In Solve NEP num_eig: {}\n", num_eig);
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
    Mpi::Print("Call MSLP\n");
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


namespace
{

enum class BasisType
{
  MONOMIAL,
  CHEBYSHEV_1
};

struct Eigenpairs
{
  Eigen::VectorXcd D;
  Eigen::MatrixXcd X;
};

Eigenpairs SolvePEP(const std::vector<const Eigen::MatrixXcd *> &P, BasisType type)
{
  // Solve the polynomial EVP: P(λ) x = (∑ᵢ λⁱ Pᵢ) x = 0 by linearization. If provided in
  // the Chebyshev basis instead of the monomial basis, the problem is: P(λ) x =
  // (∑ᵢ Tᵢ(λ) Pᵢ) x = 0, where Tᵢ are the Chebyshev polynomials of the first or second
  // kind.
  const std::size_t d = P.size() - 1;
  const std::size_t n = [&P]()
  {
    for (const auto Pi : P)
    {
      if (Pi)
      {
        return Pi->rows();
      }
    }
    return 0l;
  }();
  MFEM_VERIFY(d > 0 && n > 0 && P[d], "Invalid inputs for polynomial eigenvalue problem!");

  Eigen::MatrixXcd L0 = Eigen::MatrixXcd::Zero(d * n, d * n),
                   L1 = Eigen::MatrixXcd::Zero(d * n, d * n);
  switch (type)
  {
    case BasisType::MONOMIAL:
      for (std::size_t i = 0; i < d - 1; i++)
      {
        L0.block(i * n, (i + 1) * n, n, n) = Eigen::MatrixXcd::Identity(n, n);
      }
      for (std::size_t i = 0; i < d; i++)
      {
        if (P[i])
        {
          L0.block((d - 1) * n, i * n, n, n) = -(*P[i]);
        }
      }
      for (std::size_t i = 0; i < d - 1; i++)
      {
        L1.block(i * n, i * n, n, n) = Eigen::MatrixXcd::Identity(n, n);
      }
      L1.block((d - 1) * n, (d - 1) * n, n, n) = *P[d];
      break;
    case BasisType::CHEBYSHEV_1:
      for (std::size_t i = 0; i < d - 1; i++)
      {
        L0.block(i * n, (i + 1) * n, n, n) = Eigen::MatrixXcd::Identity(n, n);
        if (i > 0)
        {
          L0.block(i * n, (i - 1) * n, n, n) = Eigen::MatrixXcd::Identity(n, n);
        }
      }
      for (std::size_t i = 0; i < d; i++)
      {
        if (P[i])
        {
          L0.block((d - 1) * n, i * n, n, n) = -(*P[i]);
        }
      }
      L0.block((d - 1) * n, (d - 2) * n, n, n) += *P[d];
      for (std::size_t i = 0; i < d - 1; i++)
      {
        L1.block(i * n, i * n, n, n) =
            (i > 0 ? 2.0 : 1.0) * Eigen::MatrixXcd::Identity(n, n);
      }
      L1.block((d - 1) * n, (d - 1) * n, n, n) = 2.0 * (*P[d]);
      break;
  }

  // Solve the eigenvalue problem.
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eps;
  eps.compute(L1.partialPivLu().solve(L0));
  return {eps.eigenvalues(), eps.eigenvectors()};
}

std::vector<double> Chebyshev1Nodes(int d, double a, double b)
{
  // Chebyshev nodes of the first kind.
  std::vector<double> nodes(d + 1);
  for (int i = 0; i <= d; i++)
  {
    nodes[i] = 0.5 * (a + b) + 0.5 * (b - a) * std::cos((i + 0.5) * M_PI / (d + 1.0));
  }
  return nodes;
}

template <typename MatType>
void DCTIII(double a, double b, const std::vector<double> &nodes,
            const std::vector<MatType> &T, std::vector<MatType> &P)
{
  // Inverse transform for Chebyshev nodes of the first kind.
  const std::size_t d = T.size() - 1;
  for (std::size_t i = 0; i <= d; i++)
  {
    P[i] = MatType::Zero(T[0].rows(), T[0].cols());
    const double p = (i > 0) ? 2.0 : 1.0;
    for (std::size_t j = 0; j <= d; j++)
    {
      const double xj = (2.0 * (nodes[j] - a) / (b - a)) - 1.0;
      P[i] += (p / (d + 1.0) * std::cos(i * std::acos(xj))) * T[j];
    }

    // // XX TODO UNCLEAR IF NEED DOMAIN [-1, 1] -> [a, b] SCALING HERE TOO?
    // P[i] = (1.0 / (d + 1.0)) * T[0];
    // for (std::size_t j = 1; j <= d; j++)
    // {
    //   P[i] += (2.0 / (d + 1.0) * std::cos((i + 0.5) * j * M_PI / (d + 1.0))) * T[j];
    // }
  }
}

}  // namespace


std::vector<std::complex<double>> RomOperator::ComputeEigenvalueEstimates() /*const*/
{
  // XX TODO: Not yet implemented
  //MFEM_ABORT("Eigenvalue estimates for PROM operators are not yet implemented!");
  //return {};
  Mpi::Print("PROM ComputeEigenvalueEstimates!\n");
  MFEM_VERIFY(dim_V > 0,
              "Eigenvalue estimates are only available for a PROM with nonzero dimension!");
  Eigen::VectorXcd omega;
  Eigen::MatrixXcd X;
  if (!has_A2)
  {
    if (Cr.rows() == 0)
    {
      Mpi::Print("Linear generalized EVP\n");
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
      Mpi::Print("Quadratic EVP\n");
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
    Mpi::Print("General nonlinear EVP\n");
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
        auto A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(std::abs(l.imag()),
                                                                Operator::DIAG_ZERO);
        ProjectMatInternal(space_op.GetComm(), V, *A2, Ar, r, 0);
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
            auto A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(
                std::abs(l.imag()) * (1.0 + eps), Operator::DIAG_ZERO);
            ProjectMatInternal(space_op.GetComm(), V, *A2, dT, r, 0);
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
      Mpi::Print("Call SolveNEP with sigma: {}\n", sigma);
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
    linalg::Normalize(space_op.GetComm(), u);

    // Evaluate the HDM eigenpair residual.
    std::unique_ptr<ComplexOperator> A2;
    if (has_A2)
    {
      A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(omega(i).real(),
                                                         Operator::DIAG_ZERO);
    }
    auto A =
        space_op.GetSystemMatrix(std::complex<double>(1.0, 0.0), 1i * omega(i),
                                -omega(i) * omega(i), K.get(), C.get(), M.get(), A2.get());
    A->Mult(u, r);
    double res = linalg::Norml2(space_op.GetComm(), r);


    //XX TODO DEBUG
    Mpi::Print("omega = {:e}{:+e}, ||r|| = {:e}\n",
               omega(i).real(), omega(i).imag(), res);


    // XX TODO WIP.... (scaling, storage, etc.)
  }

  // XX TODO FILTERING...

  return {omega.begin(), omega.end()};
}



std::vector<std::complex<double>> RomOperator::ComputeEigenvalueEstimates2(double start,
                                                                          double end)
{

  // XX TODO ADD EIGENVECTORS TO COMPUTE HDM RESIDUAL ERROR?

  MFEM_VERIFY(dim_V > 0,
              "Eigenvalue estimates are only available for a PROM with nonzero dimension!");
  if (!has_A2)
  {
    if (Cr.rows() == 0)
    {
      // Linear generalized EVP: M⁻¹ K x = ω² x (Eigen does not support complex-valued
      // generalized EVPs).
      Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eps;
      eps.compute(Mr.partialPivLu().solve(Kr));
      auto D = eps.eigenvalues();
      std::vector<std::complex<double>> omega(D.size());
      for (auto eig : D)
      {
        omega.push_back(std::sqrt(eig));
      }
      std::sort(omega.begin(), omega.end(),
                [](auto l, auto r) { return std::abs(l) < std::abs(r); });
      return omega;
    }
    else
    {
      // Quadratic EVP: P(ω) x = (K + iω C - ω² M) x = 0 , solved via linearization.
      const Eigen::MatrixXcd iCr = 1i * Cr;
      const Eigen::MatrixXcd nMr = -Mr;
      std::vector<const Eigen::MatrixXcd *> P = {&Kr, &iCr, &nMr};
      auto [D, X] = SolvePEP(P, BasisType::MONOMIAL);
      std::vector<std::complex<double>> omega;
      omega.reserve(D.size());
      for (auto eig : D)
      {
        if (eig.real() >= -0.1 * std::abs(eig.imag()))
        {
          omega.push_back(eig);
        }
      }
      std::sort(omega.begin(), omega.end(),
                [](auto l, auto r) { return std::abs(l) < std::abs(r); });
      return omega;
    }
  }
  else
  {
    // General nonlinear EVP: T(ω) x = (K + iω C - ω² M + A2(Re{ω})) x = 0 . We first
    // interpolate A2 as a polynomial function of ω, using the first d + 1 Chebyshev
    // polynomials, and then solve the polynomial EVP via linearization.
    constexpr std::size_t d = 4;
    const auto nodes = Chebyshev1Nodes(d, start, end);
    MFEM_VERIFY(nodes.size() == d + 1,
                "Unexpected number of Chebyshev interpolation nodes for order " << d
                                                                                << "!");

    // Get Chebyshev polynomial coefficients by inverse DCT.
    std::vector<Eigen::MatrixXcd> T2(d + 1), P2(d + 1);
    for (std::size_t i = 0; i <= d; i++)
    {
      auto A2 =
          space_op.GetExtraSystemMatrix<ComplexOperator>(nodes[i], Operator::DIAG_ZERO);
      T2[i].resize(dim_V, dim_V);
      ProjectMatInternal(space_op.GetComm(), V, *A2, T2[i], r, 0);

      // //XX TODO TEST INTERPOLATION EXACT
      // T2[i] = Kr;
      // // if (Cr.rows() > 0)
      // // {
      // //   T2[i] += 1i * nodes[i] * Cr;
      // // }
      // T2[i] -= nodes[i] * nodes[i] * Mr;
    }
    DCTIII(start, end, nodes, T2, P2);

    // Add contributions of the quadratic EVP in the Chebyshev basis.
    const std::size_t d_min = 2;
    if (d < d_min)
    {
      const auto P2b(P2);
      P2.resize(d_min + 1);
      for (std::size_t i = 0; i <= d; i++)
      {
        P2[i] = P2b[i];
      }
      for (std::size_t i = d + 1; i <= d_min; i++)
      {
        P2[i] = Eigen::MatrixXcd::Zero(dim_V, dim_V);
      }
    }
    P2[0] += Kr;
    P2[0] += -0.5 * Mr;
    if (Cr.rows() > 0)
    {
      P2[1] += 1i * Cr;
    }
    P2[2] += -0.5 * Mr;

    // P2[0] = Kr;  //XX TODO FOR TESTING...
    // P2[0] += -0.5 * Mr;
    // P2[1].setZero();
    // P2[2] = -0.5 * Mr;

    // // XX TODO WIP DEBUG
    // if (Mpi::Root(Mpi::World()))
    // {
    //   for (std::size_t i = 0; i <= std::max(d, d_min); i++)
    //   {
    //     std::cout << "P[" << i << "]:\n" << P2[i] << "\n";
    //   }
    // }

    // Solve the polynomial EVP.
    std::vector<const Eigen::MatrixXcd *> P(std::max(d, d_min) + 1);
    for (std::size_t i = 0; i <= std::max(d, d_min); i++)
    {
      P[i] = &P2[i];
    }
    auto [D, X] = SolvePEP(P, BasisType::CHEBYSHEV_1);
    std::vector<std::complex<double>> omega;
    omega.reserve(D.size());
    for (auto eig : D)
    {
      if (eig.real() >= -0.1 * std::abs(eig.imag()))
      {
        omega.push_back(eig);
      }
    }
    std::sort(omega.begin(), omega.end(),
              [](auto l, auto r) { return std::abs(l) < std::abs(r); });

    //XX TODO DEBUG
    Mpi::Print("\n\nEigenpair residuals:\n");

    ComplexVector u(r.Size());
    for (int i = 0; i < omega.size(); i++)
    {
      // Compute HDM eigenmode.
      ProlongatePROMSolution(dim_V, V, X.col(i), u);
      linalg::Normalize(space_op.GetComm(), u);

      // Evaluate the HDM eigenpair residual.
      std::unique_ptr<ComplexOperator> A2;
      if (has_A2)
      {
        A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(omega[i].real(),
                                                         Operator::DIAG_ZERO);
      }
      auto A =
          space_op.GetSystemMatrix(std::complex<double>(1.0, 0.0), 1i * omega[i],
                                -omega[i] * omega[i], K.get(), C.get(), M.get(), A2.get());
      A->Mult(u, r);
      double res = linalg::Norml2(space_op.GetComm(), r);


      //XX TODO DEBUG
      Mpi::Print("omega = {:e}{:+e}, ||r|| = {:e}\n",
               omega[i].real(), omega[i].imag(), res);


    // XX TODO WIP.... (scaling, storage, etc.)
    }

    return omega;
  }
}

}  // namespace palace
