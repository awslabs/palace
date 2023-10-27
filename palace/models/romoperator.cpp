// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "romoperator.hpp"

#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
#include <mfem.hpp>
#include "linalg/orthog.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

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
  // u = [sum_s u_s q_s / (z - z_s)] / [sum_s q_s / (z - z_s)]. The coefficients are given
  // by the right singular vector of R corresponding to the minimum singular value.
  // solution to the
  const auto S = R.rows();
  MFEM_ASSERT(S > 0 && R.cols() == S, "Invalid dimension mismatch when computing MRI!");
  // For Eigen = v3.4.0 (latest tagged release as of 10/2023)
  Eigen::JacobiSVD<Eigen::MatrixXcd> svd;
  svd.compute(R, Eigen::ComputeFullV | Eigen::NoQRPreconditioner);
  // For Eigen > v3.4.0 (GitLab repo is at v3.4.90 as of 10/2023)
  // Eigen::JacobiSVD<Eigen::MatrixXcd, Eigen::ComputeFullV | Eigen::NoQRPreconditioner>
  // svd; svd.compute(R);
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
  A2 = spaceop.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO);
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
  // as an error indicator. The complex-valued snapshot matrix U = [{u_s}] is stored by its
  // QR decomposition.
  MFEM_VERIFY(dim_Q + 1 <= Q.size(),
              "Unable to increase basis storage size, increase maximum number of vectors!");
  R.conservativeResizeLike(Eigen::MatrixXd::Zero(dim_Q + 1, dim_Q + 1));
  z.conservativeResize(dim_Q + 1);
  z(dim_Q) = omega;
  {
    Q[dim_Q] = u;
    OrthogonalizeColumn(orthog_type, comm, Q, Q[dim_Q], R.col(dim_Q).data(), dim_Q);
    R(dim_Q, dim_Q) = linalg::Norml2(comm, Q[dim_Q]);
    Q[dim_Q] *= 1.0 / R(dim_Q, dim_Q);
    dim_Q++;
  }
  ComputeMRI(R, q);
  // if (Mpi::Root(comm))
  // {
  //   std::cout << "MRI (S = " << dim_Q << "):\n"
  //   std::cout << "R =\n" << R << "\n";
  //   std::cout << "q =\n" << q << "\n";
  // }

  // Update the set of sampled parameter points.
  S.push_back(omega);
}

void RomOperator::SolvePROM(double omega, ComplexVector &u)
{
  // Assemble the PROM linear system at the given frequency. The PROM system is defined by
  // the matrix Aᵣ(ω) = Kᵣ + iω Cᵣ - ω² Mᵣ + Vᴴ A2 V(ω) and source vector RHSᵣ(ω) = iω RHS1ᵣ
  // + Vᴴ RHS2(ω). A2(ω) and RHS2(ω) are constructed only if required and are only nonzero
  // on boundaries, will be empty if not needed.
  if (has_A2)
  {
    A2 = spaceop.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO);
    ProjectMatInternal(spaceop.GetComm(), V, *A2, Ar, r, 0);
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

  if (has_RHS2)
  {
    spaceop.GetExcitationVector2(omega, RHS2);
    ProjectVecInternal(spaceop.GetComm(), V, RHS2, RHSr, 0);
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
  u = 0.0;
  for (std::size_t j = 0; j < dim_V; j += 2)
  {
    if (j + 1 < dim_V)
    {
      linalg::AXPBYPCZ(RHSr(j).real(), V[j], RHSr(j + 1).real(), V[j + 1], 1.0, u.Real());
      linalg::AXPBYPCZ(RHSr(j).imag(), V[j], RHSr(j + 1).imag(), V[j + 1], 1.0, u.Imag());
    }
    else
    {
      linalg::AXPY(RHSr(j).real(), V[j], u.Real());
      linalg::AXPY(RHSr(j).imag(), V[j], u.Imag());
    }
  }
}

double RomOperator::FindMaxError(double start, double delta, int num_steps) const
{
  // Return an estimate for argmax_z ||u(z) - V y(z)|| as argmin_z |Q(z)| with
  // Q(z) = sum_s q_z / (z - z_s) (denominator of the barycentric interpolation of u).
  if (delta < 0.0)
  {
    start = start + (num_steps - 1) * delta;
    delta = -delta;
  }
  double omega_star = 0.0, Q_star = mfem::infinity();
  for (int step = 0; step < num_steps; step++)
  {
    const double omega = start + step * delta;
    const double Q = std::abs((q.array() / (z.array() - omega)).sum());
    if (Q < Q_star)
    {
      omega_star = omega;
      Q_star = Q;
    }
  }
  MFEM_VERIFY(omega_star > 0.0, "Unable to find location for maximum error!");
  return omega_star;
}

std::vector<std::complex<double>> RomOperator::ComputeEigenvalueEstimates() const
{
  if (!has_A2)
  {
    if (Cr.rows() == 0)
    {
      // Linear generalized EVP: M⁻¹ K x = λ x (Eigen does not support complex-valued
      // generalized EVPs).
      Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eps;
      eps.compute(Mr.partialPivLu().solve(Kr), false);
      const auto lambda = eps.eigenvalues();
      return {lambda.begin(), lambda.end()};
    }
    else
    {
      // Quadratic EVP: P(λ) x = (K + λ C + λ² M) x = 0 , solved via linearization.
      Eigen::MatrixXcd L0 = Eigen::MatrixXcd::Zero(2 * dim_V, 2 * dim_V);
      L0.topRightCorner(dim_V, dim_V) = Eigen::MatrixXcd::Identity(dim_V, dim_V);
      L0.bottomLeftCorner(dim_V, dim_V) = -Kr;
      L0.bottomRightCorner(dim_V, dim_V) = -Cr;

      Eigen::MatrixXcd L1 = Eigen::MatrixXcd::Zero(2 * dim_V, 2 * dim_V);
      L1.topLeftCorner(dim_V, dim_V) = Eigen::MatrixXcd::Identity(dim_V, dim_V);
      L1.bottomRightCorner(dim_V, dim_V) = Mr;

      Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eps;
      eps.compute(L1.partialPivLu().solve(L0), false);
      const auto lambda = eps.eigenvalues();
      return {lambda.begin(), lambda.end()};
    }
  }
  else
  {
    // General nonlinear EVP: T(λ) x = (K + λ C + λ² M + A2(λ)) x = 0 .

    // XX TODO: Not yet implemented
    MFEM_ABORT("NEVP estimates for PROM operators are not yet implemented!");

    return {};
  }
}

}  // namespace palace
