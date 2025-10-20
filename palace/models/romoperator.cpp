// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "romoperator.hpp"

#include <iostream>
#include <Eigen/SVD>
#include <fmt/format.h>
#include <mfem.hpp>
#include <nlohmann/json.hpp>
#include "drivers/drivensolver.hpp"
#include "fem/bilinearform.hpp"
#include "fem/coefficient.hpp"
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "linalg/operator.hpp"
#include "linalg/orthog.hpp"
#include "linalg/rap.hpp"
#include "models/materialoperator.hpp"
#include "models/postoperator.hpp"
#include "models/postoperatorcsv.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/filesystem.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/tablecsv.hpp"
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

// void AddIntegrators(BilinearForm &a, const MaterialPropertyCoefficient *df,
//                     const MaterialPropertyCoefficient *f,
//                     const MaterialPropertyCoefficient *dfb,
//                     const MaterialPropertyCoefficient *fb,
//                     const MaterialPropertyCoefficient *fp, bool assemble_q_data = false)
// {
//   if (df && !df->empty() && f && !f->empty())
//   {
//     a.AddDomainIntegrator<CurlCurlMassIntegrator>(*df, *f);
//   }
//   else
//   {
//     if (df && !df->empty())
//     {
//       a.AddDomainIntegrator<CurlCurlIntegrator>(*df);
//     }
//     if (f && !f->empty())
//     {
//       a.AddDomainIntegrator<VectorFEMassIntegrator>(*f);
//     }
//   }
//   if (dfb && !dfb->empty() && fb && !fb->empty())
//   {
//     a.AddBoundaryIntegrator<CurlCurlMassIntegrator>(*dfb, *fb);
//   }
//   else
//   {
//     if (dfb && !dfb->empty())
//     {
//       a.AddBoundaryIntegrator<CurlCurlIntegrator>(*dfb);
//     }
//     if (fb && !fb->empty())
//     {
//       a.AddBoundaryIntegrator<VectorFEMassIntegrator>(*fb);
//     }
//   }
//   if (fp && !fp->empty())
//   {
//     a.AddDomainIntegrator<MixedVectorWeakCurlIntegrator>(*fp, true);
//     a.AddDomainIntegrator<MixedVectorCurlIntegrator>(*fp);
//   }
//   if (assemble_q_data)
//   {
//     a.AssembleQuadratureData();
//   }
// }

// auto AssembleOperator(const FiniteElementSpace &fespace,
//                       const MaterialPropertyCoefficient *df,
//                       const MaterialPropertyCoefficient *f,
//                       const MaterialPropertyCoefficient *dfb,
//                       const MaterialPropertyCoefficient *fb,
//                       const MaterialPropertyCoefficient *fp, bool skip_zeros = false,
//                       bool assemble_q_data = false)
// {
//   BilinearForm a(fespace);
//   AddIntegrators(a, df, f, dfb, fb, fp, assemble_q_data);
//   return a.Assemble(skip_zeros);
// }

constexpr auto ORTHOG_TOL = 1.0e-12;

template <typename VecType, typename ScalarType>
inline void OrthogonalizeColumn(Orthogonalization type, MPI_Comm comm,
                                const std::vector<VecType> &V, VecType &w, ScalarType *Rj,
                                int j)
{
  // Orthogonalize w against the leading j columns of V.
  switch (type)
  {
    case Orthogonalization::MGS:
      linalg::OrthogonalizeColumnMGS(comm, V, w, Rj, j);
      break;
    case Orthogonalization::CGS:
      linalg::OrthogonalizeColumnCGS(comm, V, w, Rj, j);
      break;
    case Orthogonalization::CGS2:
      linalg::OrthogonalizeColumnCGS(comm, V, w, Rj, j, true);
      break;
  }
}

// template <typename VecType, typename ScalarType>
// inline void OrthogonalizeColumnMGS(MPI_Comm comm, const std::vector<VecType> &V, VecType
// &w,
//                                    ScalarType *H, int m)
// {
//   MFEM_ASSERT(static_cast<std::size_t>(m) <= V.size(),
//               "Out of bounds number of columns for MGS orthogonalization!");
//   for (int j = 0; j < m; j++)
//   {
//     H[j] = linalg::Dot(comm, w, V[j]);  // Global inner product
//     w.Add(-H[j], V[j]);
//   }
// }

template <typename VecType, typename ScalarType>
inline void OrthogonalizeColumnNewMGS(MPI_Comm comm, const Operator &M_matrix,
                                      const std::vector<VecType> &V, VecType &w,
                                      ScalarType *Rj, int m)
{

  for (int j = 0; j < m; j++)
  {
    VecType mass_V_j;
    mass_V_j.SetSize(V[j].Size());
    mass_V_j.UseDevice(true);

    // M_matrix.Mult(V[j], mass_V_j);
    Rj[j] = linalg::Dot(comm, w, V[j]);  // Global inner product
    Mpi::Print("Rj[j] in OrthogonalizeColumnNew with j={}: {}\n", j, Rj[j]);
    w.Add(-Rj[j], V[j]);
  }
}

template <typename VecType, typename ScalarType>
inline void OrthogonalizeColumnNewCGS(MPI_Comm comm, const Operator &M_matrix,
                                      const std::vector<VecType> &V, VecType &w,
                                      ScalarType *H, int m, bool refine = false)
{
  MFEM_ASSERT(static_cast<std::size_t>(m) <= V.size(),
              "Out of bounds number of columns for CGS orthogonalization!");
  if (m == 0)
  {
    return;
  }
  VecType mass_V_j;
  mass_V_j.SetSize(V[0].Size());
  mass_V_j.UseDevice(true);

  for (int j = 0; j < m; j++)
  {
    M_matrix.Mult(V[j], mass_V_j);
    H[j] = w * mass_V_j;  // Local inner product
  }
  Mpi::GlobalSum(m, H, comm);
  for (int j = 0; j < m; j++)
  {
    w.Add(-H[j], V[j]);
  }
  if (refine)
  {
    std::vector<ScalarType> dH(m);
    for (int j = 0; j < m; j++)
    {
      M_matrix.Mult(V[j], mass_V_j);
      dH[j] = w * mass_V_j;  // Local inner product
    }
    Mpi::GlobalSum(m, dH.data(), comm);
    for (int j = 0; j < m; j++)
    {
      H[j] += dH[j];
      w.Add(-dH[j], V[j]);
    }
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

MinimalRationalInterpolation::MinimalRationalInterpolation(std::size_t max_size)
{
  z.reserve(max_size);
  Q.resize(max_size, ComplexVector());
}

void MinimalRationalInterpolation::AddSolutionSample(double omega, const ComplexVector &u,
                                                     MPI_Comm comm,
                                                     Orthogonalization orthog_type)
{
  // Compute the coefficients for the minimal rational interpolation of the state u used
  // as an error indicator. The complex-valued snapshot matrix U = [{u_i, (iω) u_i}] is
  // stored by its QR decomposition.
  MFEM_VERIFY(dim_Q + 1 <= Q.size(),
              "Unable to increase basis storage size, increase maximum number of vectors!");
  R.conservativeResizeLike(Eigen::MatrixXd::Zero(dim_Q + 1, dim_Q + 1));
  {
    std::vector<const ComplexVector *> blocks = {&u, &u};
    std::vector<std::complex<double>> s = {1.0 + 0i, 1i * omega};
    Q[dim_Q].SetSize(2 * u.Size());
    Q[dim_Q].UseDevice(true);
    Q[dim_Q].SetBlocks(blocks, s);
  }
  OrthogonalizeColumn(Orthogonalization::CGS2, comm, Q, Q[dim_Q], R.col(dim_Q).data(),
                      dim_Q);
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

std::vector<double> MinimalRationalInterpolation::FindMaxError(std::size_t N) const
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
  //
  // TODO: Currently we always use this brute for sampling and we could optimize this more.
  // Also, the case of N>1 samples is not very useful below. It will typically give us
  // multiple sample points right next to each other in the same local maximum, rather than
  // N separate local maxima.

  // We could use priority queue here to keep the N lowest values. However, we don't use
  // std::priority_queue class since we want to have access to the vector and also binary
  // tree structure of heap class as rebalancing is excessive overhead for tiny size N.
  using q_t = std::pair<std::complex<double>, double>;
  std::vector<q_t> queue{};
  queue.reserve(N);

  const std::size_t nr_sample = 1.0e6;  // must be >= N
  MFEM_VERIFY(N < nr_sample,
              fmt::format("Number of location of error maximum N={} needs to be less than "
                          "the fine sampling grid nr_sample={}.",
                          N, nr_sample));
  const auto delta = (end - start) / nr_sample;
  for (double z_sample = start; z_sample <= end; z_sample += delta)
  {
    const double Q_sample = std::abs((q.array() / (z_map.array() - z_sample)).sum());

    bool partial_full = (queue.size() < N);
    if (partial_full || Q_sample < queue.back().second)
    {
      auto it_loc = std::upper_bound(queue.begin(), queue.end(), Q_sample,
                                     [](double q, const q_t &p2) { return q < p2.second; });
      queue.insert(it_loc, std::make_pair(z_sample, Q_sample));
      if (!partial_full)
      {
        queue.pop_back();
      }
    }
  }
  MFEM_VERIFY(queue.size() == N,
              fmt::format("Internal failure: queue should be size should be N={} (got {})",
                          N, queue.size()));

  std::vector<double> vals(z_star.size());
  std::transform(queue.begin(), queue.end(), vals.begin(),
                 [](const q_t &p) { return p.first.real(); });
  return vals;
}

RomOperator::RomOperator(const IoData &iodata, SpaceOperator &space_op,
                         std::size_t max_size_per_excitation)
  : space_op(space_op), orthog_type(Orthogonalization::MGS)
{
  // Construct the system matrices defining the linear operator. PEC boundaries are handled
  // simply by setting diagonal entries of the system matrix for the corresponding dofs.
  // Because the Dirichlet BC is always homogeneous, no special elimination is required on
  // the RHS. The damping matrix may be nullptr.
  K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
  C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  MFEM_VERIFY(K && M, "Invalid empty HDM matrices when constructing PROM!");

  Kinnerproduct = BuildParSumOperator({1.0 + 0i}, {K.get()}, true);
  Cinnerproduct = BuildParSumOperator({1.0 + 0i}, {C.get()}, true);
  Minnerproduct = BuildParSumOperator({1.0 + 0i}, {M.get()}, true);

  BilinearForm m_ij_bulk(space_op.GetNDSpace());
  m_ij_bulk.AddDomainIntegrator<VectorFEMassIntegrator>();
  M_inner_product_weight_base = m_ij_bulk.FullAssemble(false);

  M_inner_product_weight_ksp =
      std::make_unique<ParOperator>(*M_inner_product_weight_base, space_op.GetNDSpace());
  M_inner_product_weight_ksp->SetEssentialTrueDofs(space_op.GetNDDbcTDofLists().back(),
                                                   Operator::DIAG_ONE);

  for (auto &[port_idx, port] : space_op.GetLumpedPortOp())
  {
    const auto &mesh = space_op.GetNDSpace().GetParMesh();

    mfem::Array<int> attr_marker;

    mfem::Array<int> attr_list;
    for (const auto &elem : port.elems)
    {
      attr_list.Append(elem->GetAttrList());
    }
    int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
    mesh::AttrToMarker(bdr_attr_max, attr_list, attr_marker);
    SumVectorCoefficient fb(space_op.GetMesh().SpaceDimension());
    for (const auto &elem : port.elems)
    {
      const double E_inc = 1.0 / std::sqrt(elem->GetGeometryWidth() *
                                           elem->GetGeometryLength() * port.elems.size());
      fb.AddCoefficient(elem->GetModeCoefficient(E_inc));
    }
    auto s = std::make_unique<mfem::LinearForm>(&space_op.GetNDSpace().Get());
    s->AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fb), attr_marker);
    s->UseFastAssembly(false);
    s->UseDevice(false);
    s->Assemble();
    s->UseDevice(true);

    port_forms.emplace_back(std::move(s));
  }

  // Initialize working vector storage.
  r.SetSize(K->Height());
  r.UseDevice(true);

  // Set up the linear solver and set operators but don't set the operators yet (this will
  // be done during an HDM solve at a given parameter point). The preconditioner for the
  // complex linear system is constructed from a real approximation to the complex system
  // matrix.
  ksp = std::make_unique<ComplexKspSolver>(iodata, space_op.GetNDSpaces(),
                                           &space_op.GetH1Spaces());

  ksp_rom =
      std::make_unique<KspSolver>(iodata, space_op.GetNDSpaces(), &space_op.GetH1Spaces());

  MFEM_VERIFY(max_size_per_excitation > 0, "Reduced order basis must have > 0 size!");

  auto max_prom_size = 2 * max_size_per_excitation * space_op.GetPortExcitations().Size();
  if (iodata.solver.driven.adaptive_circuit_synthesis)
  {
    max_prom_size += space_op.GetLumpedPortOp().Size();  // Lumped ports are real fields
  }

  // Reserve empty vectors.
  V.reserve(max_prom_size);
  v_node_label.reserve(max_prom_size);

  // Set up MinimalRationalInterpolation.
  for (const auto &[excitation_idx, data] : space_op.GetPortExcitations())
  {
    mri.emplace(excitation_idx, MinimalRationalInterpolation(max_size_per_excitation));
  }
}

void RomOperator::SetExcitationIndex(int excitation_idx)
{
  // Return if cached. Ctor constructs with excitation_idx_cache = 0 which is not a valid
  // excitation index, so this is triggered the first time it is called in drivensolver.
  if (excitation_idx_cache == excitation_idx)
  {
    return;
  }

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
    auto dim_V = V.size();
    if (dim_V > 0)
    {
      MPI_Comm comm = space_op.GetComm();
      RHS1r.conservativeResize(dim_V);
      ProjectVecInternal(comm, V, RHS1, RHS1r, 0);
    }
  }
}

void RomOperator::SolveHDM(int excitation_idx, double omega, ComplexVector &u)
{
  SetExcitationIndex(excitation_idx);

  // Compute HDM solution at the given frequency. The system matrix, A = K + iω C - ω² M +
  // A2(ω) is built by summing the underlying operator contributions.
  A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO);
  has_A2 = (A2 != nullptr);
  auto A = space_op.GetSystemMatrix(std::complex<double>(1.0, 0.0), 1i * omega,
                                    std::complex<double>(-omega * omega, 0.0), K.get(),
                                    C.get(), M.get(), A2.get());
  auto P = space_op.GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, 1i * omega,
                                                             -omega * omega + 0.0i, omega);
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
void RomOperator::UpdatePROMPortHack(const ComplexVector &port_form_vector,
                                     std::string_view node_label)
{
  // ComplexVector out;
  // out.UseDevice(true);
  // out.SetSize(port_form_vector.Size());
  // out = 0.0;
  // // auto P = space_op.GetPreconditionerMatrix<Operator, double>(0.0, 0.0, 1.0, 0.0);
  // // ksp_rom->SetOperators(*M_inner_product_weight_ksp, *P);
  // // ksp_rom->Mult(port_form_vector.Real(), out.Real());

  // mfem::GridFunction E_port_vec_lifted(&space_op.GetNDSpace().Get());
  // E_port_vec_lifted = 0.0;

  // mfem::SparseMatrix A;
  // mfem::Vector B, X;

  // mfem::Array<int> ess_dof;

  // mfem::BilinearForm m_ij_bulk(&space_op.GetNDSpace().Get());
  // m_ij_bulk.AddDomainIntegrator(new mfem::VectorFEMassIntegrator());
  // m_ij_bulk.Assemble(true);

  // m_ij_bulk.FormLinearSystem(ess_dof, E_port_vec_lifted, *port_forms.at(0), A, X, B);

  // mfem::CGSolver cg(MPI_COMM_WORLD);
  // cg.SetRelTol(1e-12);
  // cg.SetMaxIter(2000);
  // cg.SetPrintLevel(1);
  // cg.SetOperator(A);
  // cg.Mult(B, out.Real());

  // // m_ij_bulk.RecoverFEMSolution(X, E_port_vec.Real(), E_port_vec_lifted);

  // auto dot = port_form_vector.Dot(out);
  // Mpi::Print("\nSurface Form inner with Domain Primary ({}): {},{}\n", node_label,
  //            dot.real(), dot.imag());

  // for (int i = 0; i < port_forms.size(); i++)
  // {
  //   auto dot3 = (*port_forms.at(i)) * out.Real();
  //   Mpi::GlobalSum(1, &dot3, space_op.GetComm());
  //   Mpi::Print("Suface g_ij with u_surf, G^-1 u_surf: {}\n", dot3);
  // }

  // ComplexVector out2;
  // out2.UseDevice(true);
  // out2.SetSize(port_form_vector.Size());
  // out2 = 0.0;

  // M_inner_product_weight_base->Mult(out.Real(), out2.Real());
  // auto dot2 = out2.Dot(out);

  // Mpi::Print("\nSurface Form inner with Domain Primary ({}): {},{}\n", node_label,
  //            dot2.real(), dot2.imag());

  // for (int i = 0; i < port_forms.size(); i++)
  // {
  //   auto dot3 = (*port_forms.at(i)) * out2.Real();
  //   Mpi::GlobalSum(1, &dot3, space_op.GetComm());
  //   Mpi::Print("Suface g_ij with u_surf, G^-1 u_surf: {}\n", dot3);
  // }

  // // out *= 1.0 / std::sqrt(dot2);

  UpdatePROM(port_form_vector, node_label);
}

void RomOperator::AddLumpedPortModesForSynthesis(const IoData &iodata)
{
  auto &lumped_port_op = space_op.GetLumpedPortOp();
  for (const auto &[port_idx, port_data] : lumped_port_op)
  {
    ComplexVector port_excitation_E;
    port_excitation_E.UseDevice(true);
    space_op.GetLumpedPortExcitationVector(port_idx, port_excitation_E, true);
    UpdatePROM(port_excitation_E, fmt::format("port_{:d}", port_idx));
  }

  // Debug Print ROM Matrices
  if constexpr (true)  // ignore
  {
    fs::path folder_tmp = fs::path(iodata.problem.output) / "prom_port_debug";
    fs::create_directories(folder_tmp);
    PrintPROMMatrices(iodata.units, folder_tmp);
  }

  const auto &orth_R = GetRomOrthogonalityMatrix();
  MFEM_VERIFY(orth_R.isDiagonal(), "Lumped port modes should have exactly zero overlap.");

  // WIP Tests

  ComplexVector out;
  out.UseDevice(true);
  out.SetSize(port_form_vector.Size());
  out = 0.0;
  // auto P = space_op.GetPreconditionerMatrix<Operator, double>(0.0, 0.0, 1.0, 0.0);
  // auto ksp_rom =
  //     std::make_unique<KspSolver>(iodata, space_op.GetNDSpaces(),
  //     &space_op.GetH1Spaces());
  // ksp_rom->SetOperators(*M_inner_product_weight_ksp, *P);
  // ksp_rom->Mult(port_form_vector.Real(), out.Real());

  mfem::GridFunction E_port_vec_lifted(&space_op.GetNDSpace().Get());
  E_port_vec_lifted = 0.0;

  mfem::SparseMatrix A;
  mfem::Vector B, X;

  mfem::Array<int> ess_dof;

  mfem::BilinearForm m_ij_bulk(&space_op.GetNDSpace().Get());
  m_ij_bulk.AddDomainIntegrator(new mfem::VectorFEMassIntegrator());
  m_ij_bulk.Assemble(true);

  m_ij_bulk.FormLinearSystem(ess_dof, E_port_vec_lifted, *port_forms.at(0), A, X, B);

  mfem::CGSolver cg(MPI_COMM_WORLD);
  cg.SetRelTol(1e-12);
  cg.SetMaxIter(2000);
  cg.SetPrintLevel(1);
  cg.SetOperator(A);
  cg.SetOperator(A);
  cg.Mult(B, out.Real());

  // m_ij_bulk.RecoverFEMSolution(X, E_port_vec.Real(), E_port_vec_lifted);

  auto dot = port_form_vector.Dot(out);
  Mpi::Print("\nSurface Form inner with Domain Primary ({}): {},{}\n", node_label,
             dot.real(), dot.imag());

  for (int i = 0; i < port_forms.size(); i++)
  {
    auto dot3 = (*port_forms.at(i)) * out.Real();
    Mpi::GlobalSum(1, &dot3, space_op.GetComm());
    Mpi::Print("Suface g_ij with u_surf, G^-1 u_surf: {}\n", dot3);
  }

  ComplexVector out2;
  out2.UseDevice(true);
  out2.SetSize(port_form_vector.Size());
  out2 = 0.0;

  M_inner_product_weight_base->Mult(out.Real(), out2.Real());
  auto dot2 = out2.Dot(out);

  Mpi::Print("\nSurface Form inner with Domain Primary ({}): {},{}\n", node_label,
             dot2.real(), dot2.imag());

  for (int i = 0; i < port_forms.size(); i++)
  {
    auto dot3 = (*port_forms.at(i)) * out2.Real();
    Mpi::GlobalSum(1, &dot3, space_op.GetComm());
    Mpi::Print("Suface g_ij with u_surf, G^-1 u_surf: {}\n", dot3);
  }

  // out *= 1.0 / std::sqrt(dot2);

  UpdatePROM(out, node_label);
}

void RomOperator::UpdatePROM(const ComplexVector &u, std::string_view node_label)
{
  // Update PROM basis V. The basis is always real (each complex solution adds two basis
  // vectors, if it has a nonzero real and imaginary parts).
  BlockTimer bt(Timer::CONSTRUCT_PROM);
  MPI_Comm comm = space_op.GetComm();

  const auto norm_re = linalg::Norml2(comm, u.Real());
  const auto norm_im = linalg::Norml2(comm, u.Imag());
  const auto norm_tol = ORTHOG_TOL * std::sqrt(norm_re * norm_re + norm_im * norm_im);
  const bool has_real = (norm_re > norm_tol);
  const bool has_imag = (norm_im > norm_tol);

  const std::size_t dim_V_old = V.size();
  const std::size_t dim_V_new =
      V.size() + static_cast<std::size_t>(has_real) + static_cast<std::size_t>(has_imag);

  orth_R.conservativeResizeLike(Eigen::MatrixXd::Zero(dim_V_new, dim_V_new));

  auto add_real_vector_to_basis = [this](const Vector &vector)
  {
    auto dim_V = V.size();
    auto &v = V.emplace_back(vector);
    OrthogonalizeColumnNewCGS(space_op.GetComm(), *M_inner_product_weight_base, V, v,
                              orth_R.col(dim_V).data(), dim_V, true);

    Vector mass_V_j;
    mass_V_j.SetSize(V.back().Size());
    mass_V_j.UseDevice(true);

    M_inner_product_weight_base->Mult(V.back(), mass_V_j);
    orth_R(dim_V, dim_V) = std::sqrt(std::abs(
        linalg::Dot(space_op.GetComm(), V.back(), mass_V_j)));  // Global inner product

    // orth_R(dim_V, dim_V) = linalg::Norml2(space_op.GetComm(), v);
    v *= 1.0 / orth_R(dim_V, dim_V);
  };

  port_form_overlap.conservativeResize(dim_V_new, port_forms.size());
  if (has_real)
  {
    add_real_vector_to_basis(u.Real());
    v_node_label.emplace_back(fmt::format("{}_re", node_label));

    auto dim_V = V.size() - 1;
    for (int i = 0; i < port_forms.size(); i++)
    {
      auto dot = (*port_forms[i]) * (V.back());
      Mpi::GlobalSum(1, &dot, space_op.GetComm());
      port_form_overlap(dim_V, i) = dot;
    }
  }
  if (has_imag)
  {
    add_real_vector_to_basis(u.Imag());
    v_node_label.emplace_back(fmt::format("{}_im", node_label));

    auto dim_V = V.size() - 1;
    for (int i = 0; i < port_forms.size(); i++)
    {
      auto dot = (*port_forms.at(i)) * (V.back());
      Mpi::GlobalSum(1, &dot, space_op.GetComm());
      port_form_overlap(dim_V, i) = dot;
    }
  }

  // Update reduced-order operators. Resize preserves the upper dim0 x dim0 block of each
  // matrix and first dim0 entries of each vector and the projection uses the values
  // computed for the unchanged basis vectors.

  // mfem::Array<int> ess_dof;
  // dynamic_cast<ComplexParOperator*>(C.get())->SetEssentialTrueDofs(ess_dof);

  // auto Ckeep = dynamic_cast<ComplexParOperator*>(C.get());
  // Ckeep->SetEssentialTrueDofs(ess_dof, Operator::DiagonalPolicy::DIAG_KEEP);

  // // Update reduced-order operators.
  // Kr.conservativeResize(dim_V_new, dim_V_new);
  // ProjectMatInternal(comm, V, *K, Kr, r, dim_V_old);
  // if (C)
  // {
  //   Cr.conservativeResize(dim_V_new, dim_V_new);
  //   ProjectMatInternal(comm, V, *C, Cr, r, dim_V_old);
  // }
  // Mr.conservativeResize(dim_V_new, dim_V_new);
  // ProjectMatInternal(comm, V, *M, Mr, r, dim_V_old);
  Kr.conservativeResize(dim_V_new, dim_V_new);
  ProjectMatInternal(comm, V, *Kinnerproduct, Kr, r, dim_V_old);
  if (C)
  {
    Cr.conservativeResize(dim_V_new, dim_V_new);
    ProjectMatInternal(comm, V, *Cinnerproduct, Cr, r, dim_V_old);
  }
  Mr.conservativeResize(dim_V_new, dim_V_new);
  ProjectMatInternal(comm, V, *Minnerproduct, Mr, r, dim_V_old);
  if (RHS1.Size())
  {
    RHS1r.conservativeResize(dim_V_new);
    ProjectVecInternal(comm, V, RHS1, RHS1r, dim_V_old);
  }
}

void RomOperator::UpdateMRI(int excitation_idx, double omega, const ComplexVector &u)
{
  BlockTimer bt(Timer::CONSTRUCT_PROM);
  mri.at(excitation_idx).AddSolutionSample(omega, u, space_op.GetComm(), orthog_type);
}

void RomOperator::SolvePROM(int excitation_idx, double omega, ComplexVector &u)
{
  SetExcitationIndex(excitation_idx);

  // Assemble the PROM linear system at the given frequency. The PROM system is defined by
  // the matrix Aᵣ(ω) = Kᵣ + iω Cᵣ - ω² Mᵣ + Vᴴ A2 V(ω) and source vector RHSᵣ(ω) =
  // iω RHS1ᵣ + Vᴴ RHS2(ω). A2(ω) and RHS2(ω) are constructed only if required and are
  // only nonzero on boundaries, will be empty if not needed.

  // No basis states ill-defined: return zero vector to match current behaviour.
  if (V.size() == 0)
  {
    u = 0.0;
    return;
  }

  Ar.resize(V.size(), V.size());
  RHSr.resize(V.size());

  if (has_A2)
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

  if (has_RHS2)
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
  ProlongatePROMSolution(V.size(), V, RHSr, u);
}

std::vector<std::complex<double>> RomOperator::ComputeEigenvalueEstimates() const
{
  // XX TODO: Not yet implemented
  MFEM_ABORT("Eigenvalue estimates for PROM operators are not yet implemented!");
  return {};
}

void RomOperator::PrintPROMMatrices(const Units &units, const fs::path &post_dir) const
{
  BlockTimer bt0(Timer::POSTPRO);
  Mpi::Print(" Printing PROM Matrices to disk.\n");

  if (!Mpi::Root(space_op.GetComm()))
  {
    return;
  }
  auto print_table = [post_dir, labels = this->v_node_label](const Eigen::MatrixXd &mat,
                                                             std::string_view filename)
  {
    // MFEM_VERIFY(labels.size() == mat.cols(), "Inconsistent PROM size!");
    auto out = TableWithCSVFile(post_dir / filename);
    out.table.col_options.float_precision = 17;
    for (long i = 0; i < mat.cols(); i++)
    {
      out.table.insert(labels[i], labels[i]);
      auto &col = out.table[labels[i]];
      for (long j = 0; j < mat.rows(); j++)
      {
        col << mat(j, i);
      }
    }
    out.WriteFullTableTrunc();
  };

  // De-normalize PROM matrices voltages (both port and sampled). Define so that 1.0 on port
  // i corresponds to full (un-normalized solution), so you can use Linv, Rinv, C directly.
  //
  // In more detail: there are two normalizations  associated with ports: alpha & beta in
  // the notation of Marks and Williams  "A General Waveguide Circuit Theory" [J. Res. Natl.
  // Inst. Stand. Technol. 97, 533 (1992)]. "alpha" is the normalization power through the
  // port mode. Palace defines the input power to be normalized to 1 (see
  // LumpedPortData::GetExcitationPower). "beta" is the normalization convention of the port
  // voltage \beta * v_0 and current i_0 / \beta^* at constant power.

  // auto v_d = orth_R.diagonal().cwiseInverse().asDiagonal();
  auto v_d = Eigen::MatrixXcd::Identity(orth_R.rows(), orth_R.cols());
  // Note: When checking for imaginary parts, it is better to do this for K,C,M as this is a
  // nullptr check. Kr, Mr, Cr would require a numerical check on imag elements.

  print_table(port_form_overlap.real(), "port_form_overlap_re.csv");
  print_table(port_form_overlap.imag(), "port_form_overlap_im.csv");

  auto unit_henry = units.GetScaleFactor<Units::ValueType::INDUCTANCE>();
  auto m_Linv = ((1.0 / unit_henry) * v_d) * Kr * v_d;
  print_table(m_Linv.real(), "rom-Linv-re.csv");
  if (K->Imag())
  {
    print_table(m_Linv.imag(), "rom-Linv-im.csv");
  }

  auto unit_farad = units.GetScaleFactor<Units::ValueType::CAPACITANCE>();
  auto m_C = (unit_farad * v_d) * Mr * v_d;
  print_table(m_C.real(), "rom-C-re.csv");
  if (M->Imag())
  {
    print_table(m_C.imag(), "rom-C-im.csv");
  }

  // C & Cr are optional in UpdatePROM so follow this here. In practice, Cr always exists
  // since we need dissipative ports for a driven response, but this may change.
  if (C)
  {
    auto unit_ohm = units.GetScaleFactor<Units::ValueType::IMPEDANCE>();
    auto m_Rinv = ((1.0 / unit_ohm) * v_d) * Cr * v_d;
    print_table(m_Rinv.real(), "rom-Rinv-re.csv");
    if (C->Imag())
    {
      print_table(m_Rinv.imag(), "rom-Rinv-im.csv");
    }
  }

  // Print orth-R. Don't divide by diagonal to keep state normalization info.
  print_table(orth_R, "rom-orthogonalization-matrix-R.csv");

  // WIP: Compute and print boundary overlap matrix

  auto &mesh = space_op.GetNDSpace().GetParMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;

  SumCoefficient fb_i{};
  SumVectorCoefficient fb_v(mesh.SpaceDimension());

  std::vector<ComplexVector> port_v_primary(space_op.GetLumpedPortOp().Size());

  int port_idx = 0;
  for (const auto &[port_i, port] : space_op.GetLumpedPortOp())
  {
    port_v_primary.at(port_idx).SetSize(space_op.GetNDSpace().GetTrueVSize());
    port_v_primary.at(port_idx).UseDevice(true);
    port_v_primary.at(port_idx) = 0.0;

    mfem::Array<int> attr_list;

    SumVectorCoefficient fb(mesh.SpaceDimension());
    for (const auto &elem : port.elems)
    {
      auto port_element_scale_factor =
          1.0 / std::sqrt(elem->GetGeometryWidth() * elem->GetGeometryLength() *
                          port.elems.size());

      attr_list.Append(elem->GetAttrList());
      fb_v.AddCoefficient(elem->GetModeCoefficient(port_element_scale_factor));
      fb_i.AddCoefficient(
          std::make_unique<RestrictedCoefficient<mfem::ConstantCoefficient>>(
              elem->GetAttrList(), port_element_scale_factor));
      fb.AddCoefficient(elem->GetModeCoefficient(port_element_scale_factor));
    }

    {
      mesh::AttrToMarker(bdr_attr_max, attr_list, attr_marker);

      GridFunction E_port_vec(space_op.GetNDSpace());
      E_port_vec = 0.0;
      E_port_vec.Real().ProjectBdrCoefficientTangent(fb, attr_marker);
      E_port_vec.Real().GetTrueDofs(port_v_primary.at(port_idx).Real());
    }
    port_idx++;
  }

  mfem::BilinearForm g_ij_ports_id(&space_op.GetNDSpace().Get());
  g_ij_ports_id.AddBoundaryIntegrator(new mfem::VectorFEMassIntegrator(fb_i));
  g_ij_ports_id.Assemble();

  mfem::BilinearForm g_ij_ports_vec(&space_op.GetNDSpace().Get());
  g_ij_ports_vec.AddBoundaryIntegrator(new mfem::VectorFEMassIntegrator(fb_v));
  g_ij_ports_vec.Assemble();

  Eigen::MatrixXd boundary_overlap_i(V.size(), V.size());
  Eigen::MatrixXd boundary_overlap_v(V.size(), V.size());
  Eigen::MatrixXd boundary_overlap_primary(V.size(), V.size());
  Eigen::MatrixXd direct_orthogonality_check(V.size(), V.size());

  boundary_overlap_primary.setZero();

  for (long i = 0; i < V.size(); i++)
  {
    for (long j = 0; j < space_op.GetLumpedPortOp().Size(); j++)
    {
      boundary_overlap_i(i, j) = g_ij_ports_id.InnerProduct(V[i], port_v_primary[j].Real());
      boundary_overlap_v(i, j) =
          g_ij_ports_vec.InnerProduct(V[i], port_v_primary[j].Real());
    }
    for (long j = space_op.GetLumpedPortOp().Size(); j < V.size(); j++)
    {
      boundary_overlap_i(i, j) = g_ij_ports_id.InnerProduct(V[i], V[j]);
      boundary_overlap_v(i, j) = g_ij_ports_vec.InnerProduct(V[i], V[j]);
    }
  }

  for (long i = 0; i < space_op.GetLumpedPortOp().Size(); i++)
  {
    for (long j = 0; j < space_op.GetLumpedPortOp().Size(); j++)
    {
      boundary_overlap_primary(i, j) =
          g_ij_ports_id.InnerProduct(port_v_primary[i].Real(), port_v_primary[j].Real());
    }
  }
  Vector mass_V_j;
  mass_V_j.SetSize(V[0].Size());
  mass_V_j.UseDevice(true);

  for (long i = 0; i < V.size(); i++)
  {
    for (long j = 0; j < V.size(); j++)
    {
      M_inner_product_weight_base->Mult(V[j], mass_V_j);
      direct_orthogonality_check(i, j) = V[i] * mass_V_j;
    }
  }

  print_table(boundary_overlap_i, "boundary_overlap_i.csv");
  print_table(boundary_overlap_v, "boundary_overlap_v.csv");
  print_table(boundary_overlap_primary, "boundary_overlap_primary.csv");
  print_table(direct_orthogonality_check, "direct_orthogonality_check.csv");
}

}  // namespace palace
