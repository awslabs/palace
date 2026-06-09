// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "romoperator.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <string_view>
#include <tuple>
#include <utility>
#include <mfem.hpp>
#include "fem/bilinearform.hpp"
#include "fem/integrator.hpp"
#include "linalg/operator.hpp"
#include "linalg/orthog.hpp"
#include "linalg/rap.hpp"
#include "models/materialoperator.hpp"
#include "models/spaceoperator.hpp"
#include "utils/aaa.hpp"
#include "utils/communication.hpp"
#include "utils/filesystem.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/tablecsv.hpp"
#include "utils/timer.hpp"
#include "utils/units.hpp"

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

template <typename VecType, typename ScalarType,
          typename InnerProductW = linalg::IdentityInnerProduct>
inline void OrthogonalizeColumn(Orthogonalization type, MPI_Comm comm,
                                const std::vector<VecType> &V, VecType &w, ScalarType *Rj,
                                std::size_t j, const InnerProductW &dot_op = {})
{
  // Orthogonalize w against the leading j columns of V.
  switch (type)
  {
    case Orthogonalization::MGS:
      linalg::OrthogonalizeColumnMGS(comm, V, w, Rj, j, dot_op);
      break;
    case Orthogonalization::CGS:
      linalg::OrthogonalizeColumnCGS(comm, V, w, Rj, j, false, dot_op);
      break;
    case Orthogonalization::CGS2:
      linalg::OrthogonalizeColumnCGS(comm, V, w, Rj, j, true, dot_op);
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

// Evaluate the barycentric Lagrange interpolant of a vector-valued function f(ω) at the
// given ω, given samples {(ωₖ, fₖ)}. The samples map is keyed by ω (sorted by std::map),
// and the values are fixed-dimension Eigen::VectorXcd (the projected RHS2 in our use).
//
// Returns the result by value. If `omega` matches a sample point exactly to round-off,
// returns that sample directly to avoid 0/0. Tiny sample sets (S=1) reduce to constant
// extrapolation; S=2 to linear. Weights are recomputed at every call (S² in the number
// of samples, negligible compared to the linear-in-n cost of the result).
inline Eigen::VectorXcd
BarycentricInterpolate(const std::map<double, Eigen::VectorXcd> &samples, double omega)
{
  MFEM_VERIFY(!samples.empty(), "BarycentricInterpolate requires at least one sample!");
  if (samples.size() == 1)
  {
    return samples.begin()->second;
  }
  // Detect exact-match sample (within a relative tolerance proportional to the spread
  // of the sample frequencies, since ω⋆ values are inserted from a finite list of
  // double-precision frequencies and may compare exactly).
  double omega_max = std::abs(samples.rbegin()->first);
  double match_tol = 1.0e-14 * std::max(omega_max, 1.0);
  for (const auto &[w, v] : samples)
  {
    if (std::abs(w - omega) <= match_tol)
    {
      return v;
    }
  }
  // Recompute barycentric weights wₖ = 1 / Π_{j≠k} (ωₖ − ωⱼ).
  const std::size_t S = samples.size();
  std::vector<double> w_pts;
  std::vector<const Eigen::VectorXcd *> v_pts;
  w_pts.reserve(S);
  v_pts.reserve(S);
  for (const auto &[w, v] : samples)
  {
    w_pts.push_back(w);
    v_pts.push_back(&v);
  }
  std::vector<double> bary(S);
  for (std::size_t k = 0; k < S; k++)
  {
    double prod = 1.0;
    for (std::size_t j = 0; j < S; j++)
    {
      if (j == k)
      {
        continue;
      }
      prod *= (w_pts[k] - w_pts[j]);
    }
    bary[k] = 1.0 / prod;
  }
  // f(ω) = Σ (wₖ/(ω-ωₖ)) fₖ / Σ (wₖ/(ω-ωₖ)).
  Eigen::VectorXcd num = Eigen::VectorXcd::Zero(v_pts[0]->size());
  std::complex<double> denom{0.0, 0.0};
  for (std::size_t k = 0; k < S; k++)
  {
    double coeff = bary[k] / (omega - w_pts[k]);
    num += coeff * (*v_pts[k]);
    denom += coeff;
  }
  return num / denom;
}

inline void ComputeMRI(const Eigen::MatrixXcd &R, Eigen::VectorXcd &q)
{
  // Compute the coefficients of the minimal rational interpolation (MRI):
  // u = [sum_i u_i q_i / (z - z_i)] / [sum_i q_i / (z - z_i)]. The coefficients are given
  // by the right singular vector of R corresponding to the minimum singular value.
  const auto S = R.rows();
  MFEM_ASSERT(S > 0 && R.cols() == S, "Invalid dimension mismatch when computing MRI!");
  Eigen::JacobiSVD<Eigen::MatrixXcd, Eigen::ComputeFullV> svd;
  svd.compute(R);
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
    std::vector<std::complex<double>> s = {{1.0, 0.0}, {0.0, omega}};
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

  // std::vector<std::complex<double>> z_star(N, 0.0);
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

  // Number of uniformly spaced sample points for brute-force minimization of |Q(z)|. TODO:
  // Consider making configurable or scaling with the parameter range, so that it is e.g.
  // always 1 kHz.
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

  std::vector<double> vals(N);
  std::transform(queue.begin(), queue.end(), vals.begin(),
                 [](const q_t &p) { return p.first.real(); });
  return vals;
}

// Hybrid inner-product matrix. This is made from MassIntegrator of the domain and port
// boundaries summed together. However:
// - We zero out the domain mass matrix on the dof of the boundary, leaving on the boundary
//   mass matrix.
// - We weight the mass matrix by 1 / \eta with reference impedance \vert Z_R \vert = 1, so
//   that power orthogonality of modes is enforced.
// - We don't weight by material coefficients so that is fully real and corresponds to to
//   full space overlap (except excised bulk part).
//
// Zero out port dofs of bulk in ctor.
HybridBulkBoundaryOperator::HybridBulkBoundaryOperator(
    const SpaceOperator &space_op, DomainOrthogonalizationWeight domain_orthog_type)
{
  const auto &mat_op = space_op.GetMaterialOp();

  // Port attrs: To zero out true dof corresponding to attrs in bulk
  mfem::Array<int> port_attr_list_local{};

  // Ports:
  BilinearForm w_port(space_op.GetNDSpace());
  MaterialPropertyCoefficient fb_port(mat_op.MaxCeedBdrAttribute());

  for (const auto &[idx, data] : space_op.GetLumpedPortOp())
  {
    for (const auto &elem : data.elems)
    {
      // Want to add eta corresponding to a nominal Z_r = 1
      double eta_norm = data.GetToSquare(*elem);
      fb_port.AddMaterialProperty(data.mat_op.GetCeedBdrAttributes(elem->GetAttrList()),
                                  1.0 / eta_norm);
      port_attr_list_local.Append(elem->GetAttrList());
    }
  }
  // Add wave-port boundaries. The waveport modal field is added to the basis as a
  // synthesis port mode (see RomOperator::AddWavePortModesForSynthesis); to enforce
  // power orthogonality consistently the waveport boundary contributes the same kind of
  // overlap integral as the lumped ports. Using a unit weight here corresponds to a
  // reference impedance Z_R = 1 (in internal units) for the waveport, matching the
  // lumped-port convention. The exact value only affects diagonal scaling of synthesis
  // matrices, not their off-diagonal structure or the recovered S-parameters.
  for (const auto &[idx, data] : space_op.GetWavePortOp())
  {
    fb_port.AddMaterialProperty(mat_op.GetCeedBdrAttributes(data.GetAttrList()), 1.0);
    port_attr_list_local.Append(data.GetAttrList());
  }
  // Need to check this as this per MPI rank. Ranks where the material property is empty
  // should not add this integrator.
  if (!fb_port.empty())
  {
    w_port.AddBoundaryIntegrator<VectorFEMassIntegrator>(fb_port);
  }
  auto w_port_assemble = w_port.Assemble(false);
  W_inner_product_weight_port =
      std::make_unique<ParOperator>(std::move(w_port_assemble), space_op.GetNDSpace());

  // Convert port_attr_list into essential tdof
  int bdr_attr_max = (space_op.GetMesh().Get().bdr_attributes.Size() != 0)
                         ? space_op.GetMesh().Get().bdr_attributes.Max()
                         : 0;
  auto port_attr_marker = mesh::AttrToMarker(bdr_attr_max, port_attr_list_local);
  space_op.GetNDSpace().Get().GetEssentialTrueDofs(port_attr_marker, port_tdof_list);

  // Bulk, based on configuration (skip for FE_BASIS_IDENTITY)
  if (domain_orthog_type == DomainOrthogonalizationWeight::ENERGY ||
      domain_orthog_type == DomainOrthogonalizationWeight::SPACE_OVERLAP)
  {
    MaterialPropertyCoefficient epsilon_func = [&mat_op, &domain_orthog_type]()
    {
      if (domain_orthog_type == DomainOrthogonalizationWeight::ENERGY)
      {
        return MaterialPropertyCoefficient{mat_op.GetAttributeToMaterial(),
                                           mat_op.GetPermittivityReal()};
      }
      // SPACE_OVERLAP: Integrate \int dx E(x) E(x)
      // Use Palace existing palace machinery, but make a trivial bulk material.
      MaterialPropertyCoefficient eps_func_local(mat_op.MaxCeedAttribute());
      const auto &eps_ref = mat_op.GetPermittivityReal();
      mfem::DenseTensor eps_id(eps_ref.SizeI(), eps_ref.SizeJ(), eps_ref.SizeK());
      eps_id = 0.0;
      for (int k = 0; k < eps_id.SizeK(); k++)
      {
        for (int i = 0; i < eps_id.SizeI(); i++)
        {
          eps_id(i, i, k) = 1.0;
        }
      }
      eps_func_local.AddCoefficient(mat_op.GetAttributeToMaterial(), eps_id, 1.0);
      return eps_func_local;
    }();

    BilinearForm w_bulk(space_op.GetNDSpace());
    // Need to check this as this per MPI rank. Ranks where the material property is empty
    // should not add this integrator.
    if (!epsilon_func.empty())
    {
      w_bulk.AddDomainIntegrator<VectorFEMassIntegrator>(epsilon_func);
    }
    auto w_bulk_assemble = w_bulk.Assemble(false);
    W_inner_product_weight_bulk =
        std::make_unique<ParOperator>(std::move(w_bulk_assemble), space_op.GetNDSpace());
  }
  // Zero out port dofs of bulk in validation.
  validate_operators_zero_bulk_tdof();
}

RomOperator::RomOperator(const IoData &iodata, SpaceOperator &space_op,
                         std::size_t max_size_per_excitation)
  : space_op(space_op), orthog_type(iodata.solver.driven.adaptive_solver_gs_orthog_type)
{
  // Construct the system matrices defining the linear operator. PEC boundaries are
  // handled simply by setting diagonal entries of the system matrix for the corresponding
  // dofs. Because the Dirichlet BC is always homogeneous, no special elimination is
  // required on the RHS. The damping matrix may be nullptr.
  K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
  C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  MFEM_VERIFY(K && M, "Invalid empty HDM matrices when constructing PROM!");

  // Per-port boundary masses for wave ports (ω-independent). The wave-port contribution
  // to A(ω) is then assembled at each ω as Σ_p kₙ,p(ω)·M^(p)_{μ⁻¹,r}, where M^(p)_r is
  // projected onto the basis only when the basis grows. This avoids HDM-scale work in
  // the online phase. See WavePortOperator::AddBoundaryMassBdrCoefficients and
  // SpaceOperator::GetWavePortBoundaryMassMatrix.
  for (const auto &[port_idx, port_data] : space_op.GetWavePortOp())
  {
    auto Mp = space_op.GetWavePortBoundaryMassMatrix<ComplexOperator>(port_idx,
                                                                      Operator::DIAG_ZERO);
    if (Mp)
    {
      Mwp_p.emplace(port_idx, std::move(Mp));
    }
  }
  // Detect whether GetExtraSystemMatrix has any non-wave-port contributions (e.g.
  // second-order farfield, surface conductivity). The probe frequency does not matter
  // since these BCs are nonzero whenever configured.
  {
    auto A2_other_probe = space_op.GetExtraSystemMatrix<ComplexOperator>(
        1.0, Operator::DIAG_ZERO, /*include_wave_ports=*/false);
    has_other_A2 = (A2_other_probe != nullptr);
  }

  // Capture sweep band and synthesis tolerance for later use in
  // CalculateNormalizedPROMMatrices (polynomial fit residual + bound). After
  // config::Nondimensionalize (utils/configfile.cpp), sample_f stores
  // 2π·f_nondim — i.e. ω in nondimensional units, not the linear frequency. Use it
  // directly. Default the synthesis tol to adaptive_tol when the user did not set one
  // explicitly.
  const auto &sample_f = iodata.solver.driven.sample_f;
  if (!sample_f.empty())
  {
    sweep_omega_min = *std::min_element(sample_f.begin(), sample_f.end());
    sweep_omega_max = *std::max_element(sample_f.begin(), sample_f.end());
  }
  waveport_synthesis_tol = iodata.solver.driven.waveport_synthesis_tol > 0.0
                               ? iodata.solver.driven.waveport_synthesis_tol
                               : iodata.solver.driven.adaptive_tol;
  waveport_synthesis_order_max = iodata.solver.driven.waveport_synthesis_order_max;
  waveport_synthesis_rank_tol = iodata.solver.driven.waveport_synthesis_rank_tol;
  waveport_synthesis_force = iodata.solver.driven.waveport_synthesis_force;

  // Initialize working vector storage.
  r.SetSize(K->Height());
  r.UseDevice(true);

  // Set up the linear solver and set operators but don't set the operators yet (this will
  // be done during an HDM solve at a given parameter point). The preconditioner for the
  // complex linear system is constructed from a real approximation to the complex system
  // matrix.
  ksp = std::make_unique<ComplexKspSolver>(
      iodata.solver.linear, GetPreconditionerMatrixSymmetry(iodata), iodata.problem.verbose,
      space_op.GetNDSpaces(), &space_op.GetH1Spaces());

  MFEM_VERIFY(max_size_per_excitation > 0, "Reduced order basis must have > 0 size!");

  auto max_prom_size = 2 * max_size_per_excitation * space_op.GetPortExcitations().Size();
  if (iodata.solver.driven.adaptive_circuit_synthesis)
  {
    // Each lumped port included in synthesis contributes one real basis vector; ports
    // flagged out via IncludeInSynthesis = false add nothing. Reserve against the
    // included count, not the total port count, to avoid over-reserving basis storage
    // (one full-FE-space vector per excluded port would otherwise be reserved).
    max_prom_size += NumSynthesisPortModes();
    // Wave-port modes are added once per INCLUDED port (ports with IncludeInSynthesis =
    // false add nothing). The seeded basis vector at the reference frequency is generally
    // complex (mode field has both real and imaginary parts), so reserve up to two slots
    // per included port.
    max_prom_size += 2 * NumSynthesisWavePortModes();

    // Build inner-product weight matrix.
    weight_op_W = HybridBulkBoundaryOperator{
        space_op, iodata.solver.driven.adaptive_circuit_synthesis_domain_orthog};
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

// TODO: Add config-only constructor for unit testing once the new PROM
// circuit synthesis parameters (adaptive_circuit_synthesis,
// adaptive_circuit_synthesis_domain_orthog) are factored out of IoData.

void RomOperator::SetExcitationIndex(int excitation_idx)
{
  // Return if cached. Ctor constructs with excitation_idx_cache = 0 which is not a valid
  // excitation index, so this is triggered the first time it is called in
  // drivensolver.cpp.
  if (excitation_idx_cache == excitation_idx)
  {
    return;
  }

  // Set up RHS vector (linear in frequency part) for the incident field at port
  // boundaries, and the vector for the solution, which satisfies the Dirichlet (PEC) BC.
  excitation_idx_cache = excitation_idx;
  // Reset has_RHS2 so SolveHDM re-checks, since it may differ per excited port.
  has_RHS2 = true;
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

  // The HDM excitation vector is computed as RHS = iω RHS1 + RHS2(ω). When RHS2 is
  // present (excited wave port), capture an HDM-side copy at this ω⋆ for later
  // interpolation in SolvePROM. The HDM solution u itself is added to the basis by
  // UpdatePROM, but the projected RHS2 cannot be reconstructed from u alone, so the
  // sample must be cached here.
  Mpi::Print("\n");
  if (has_RHS2)
  {
    has_RHS2 = space_op.GetExcitationVector2(excitation_idx, omega, r);
    if (has_RHS2)
    {
      auto &cache = RHS2_hdm_samples[excitation_idx][omega];
      cache.SetSize(r.Size());
      cache.UseDevice(true);
      cache = r;
      // Project into the current basis (lazy — may be empty if this is the first
      // sample). Extension for newly added basis vectors happens in UpdatePROM.
      auto dim_V = V.size();
      auto &cache_r = RHS2_r_samples[excitation_idx][omega];
      cache_r.resize(dim_V);
      if (dim_V > 0)
      {
        ProjectVecInternal(space_op.GetComm(), V, cache, cache_r, 0);
      }
    }
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

std::size_t RomOperator::NumSynthesisPortModes() const
{
  // Each lumped port included in synthesis contributes exactly one real basis vector
  // (lumped port fields are real). Ports flagged out via IncludeInSynthesis = false
  // contribute nothing. Keep this the single source of truth for the port-mode count so
  // the reservation and the port-port block scaling cannot drift apart.
  std::size_t n = 0;
  for (const auto &[port_idx, port_data] : space_op.GetLumpedPortOp())
  {
    if (port_data.include_in_synthesis)
    {
      n++;
    }
  }
  return n;
}

std::size_t RomOperator::NumSynthesisWavePortModes() const
{
  // Each wave port included in synthesis contributes one port mode (the modal field at the
  // reference frequency). Ports flagged out via IncludeInSynthesis = false contribute
  // nothing. The mode is generally complex, so the caller reserves up to two basis vectors
  // (real + imaginary) per included port; see AddWavePortModesForSynthesis.
  std::size_t n = 0;
  for (const auto &[port_idx, port_data] : space_op.GetWavePortOp())
  {
    if (port_data.include_in_synthesis)
    {
      n++;
    }
  }
  return n;
}

void RomOperator::AddLumpedPortModesForSynthesis()
{
  // Add modes for lumped port to use them a circuit matrices.
  //
  // The excitation vector that we expect to add to the PROM is just the GridFunction
  // (primary vector) E_t which is the tangential electric field associated with that
  // port. The field is normalized according with an effective reference impedance of
  // \vert Z_R
  // \vert = 1, see SpaceOp::GetLumpedPortExcitationVectorPrimaryEt &
  // LumpedPortData::GetExcitationFieldEtNormSqWithUnityZR().
  //
  // The lumped ports currently implemented (rectangular and coax) are purely real.
  //
  // The hybrid weight matrix used for normalization weight_op_W will guarantee that the
  // generic vectors added to the ROM will be orthogonal with respect to the boundary
  // bilinear v_rom * g_port_boundary * e_t = 0 (unless v_rom == e_t). If we used a
  // conventional mass matrix of the space, this orthogonality would not be true as DoFs
  // of the boundary also contribute the bulk integral ('leak into the bulk').
  //
  // To have a sensible scattering matrix, we ensure port modes are orthogonal. Ports that
  // neighbor each other in space could fail this since they share DoFs on the finite
  // element mesh, even if they would be orthogonal in continuous space.

  ComplexVector vec;  // Workspace vector:  UpdatePROM interface requires ComplexVectors
  vec.SetSize(space_op.GetNDSpace().GetTrueVSize());
  vec.UseDevice(true);
  vec = 0.0;

  for (const auto &[port_idx, port_data] : space_op.GetLumpedPortOp())
  {
    if (!port_data.include_in_synthesis)
    {
      // The boundary condition for this port is still applied (see
      // LumpedPortOperator), but no port-mode vector is added to the PROM basis.
      // Excited ports always have include_in_synthesis = true (enforced by the config
      // parser) so the excitation vector is never silently dropped here.
      continue;
    }
    space_op.GetLumpedPortExcitationVectorPrimaryEt(port_idx, vec);
    UpdatePROM(vec, fmt::format("port_{:d}", port_idx));
  }

  // Check that the ports don't have any overlap. These should be exactly zero, if the ports
  // are distinct. But add this here for future case where lumped ports could on same attrs.
  // Mix orthogonalization error and reduction error.
  auto diag_tol = std::max(ORTHOG_TOL, std::numeric_limits<double>::epsilon() *
                                           std::sqrt(Mpi::Size(space_op.GetComm())));
  MFEM_VERIFY(orth_R.isDiagonal(diag_tol),
              "Lumped port fields on the mesh should have exactly zero overlap. This may "
              "be non-zero if attributes share edges.");
}

void RomOperator::AddWavePortModesForSynthesis(double omega_ref)
{
  // Add the modal field of each wave port to the basis as a port mode for synthesis,
  // analogous to AddLumpedPortModesForSynthesis. The modal field is evaluated at a
  // single reference frequency (typically the band centre) and projected onto the
  // parent ND space, restricted to the port boundary attributes.
  //
  // The mode field is generally complex; UpdatePROM splits it into real and imaginary
  // basis vectors as needed (only those that pass the orthogonalisation tolerance are
  // retained). The hybrid weight matrix W (now extended to include wave-port boundary
  // mass, see HybridBulkBoundaryOperator::HybridBulkBoundaryOperator) enforces the
  // boundary-overlap orthogonality condition for synthesis.
  ComplexVector vec;
  vec.SetSize(space_op.GetNDSpace().GetTrueVSize());
  vec.UseDevice(true);
  vec = 0.0;

  for (const auto &[port_idx, port_data] : space_op.GetWavePortOp())
  {
    if (!port_data.include_in_synthesis)
    {
      // The boundary condition for this port is still applied (see WavePortOperator),
      // but no port-mode vector is added to the PROM basis. Excited ports always have
      // include_in_synthesis = true (enforced by the config parser) so the excitation
      // vector is never silently dropped here.
      continue;
    }
    space_op.GetWavePortFieldVectorPrimaryEt(port_idx, omega_ref, vec);
    UpdatePROM(vec, fmt::format("waveport_{:d}", port_idx));
  }
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
  std::size_t dim_V_new = V.size() + std::size_t{has_real} + std::size_t{has_imag};

  orth_R.conservativeResizeLike(Eigen::MatrixXd::Zero(dim_V_new, dim_V_new));

  // Small lambda to add vector to basis. Lambda returns a bool, which is false when the new
  // vector is below the linear dependence tolerance. The MFEM_VERIFY happens after the
  // lambda. This is done for MPI syncronization reasons. If the MFEM_VERIFY throws / aborts
  // on one rank why others are in an inconsistent state, this could lead to MPI
  // communication failures in the memory unwinding and global MPI hard crash (e.g. in the
  // unit tests).
  auto add_real_vector_to_basis = [this](const Vector &vector,
                                         std::string_view node_label) -> bool
  {
    auto dim_V = V.size();
    MFEM_VERIFY(dim_V < V.capacity(),
                "PROM basis storage exceeded. Please increase maximum number of prom "
                "vector per excitation using AdaptiveMaxSamples.");
    auto &v = V.emplace_back(vector);
    double pre_norm;
    if (weight_op_W.has_value())
    {
      auto pre_norm_sq = weight_op_W->InnerProduct(space_op.GetComm(), v, v, r.Real());
      pre_norm = std::sqrt(std::abs(pre_norm_sq));
      OrthogonalizeColumn(
          orthog_type, space_op.GetComm(), V, v, orth_R.col(dim_V).data(), dim_V,
          [&W = *(this->weight_op_W), &r = this->r](const Vector &x, const Vector &y)
          { return W.InnerProduct(x, y, r.Real()); });
      auto norm_sq = weight_op_W->InnerProduct(space_op.GetComm(), v, v, r.Real());
      orth_R(dim_V, dim_V) = std::sqrt(std::abs(norm_sq));
    }
    else
    {
      pre_norm = linalg::Norml2(space_op.GetComm(), v);
      OrthogonalizeColumn(orthog_type, space_op.GetComm(), V, v, orth_R.col(dim_V).data(),
                          dim_V);
      orth_R(dim_V, dim_V) = linalg::Norml2(space_op.GetComm(), v);
    }

    if (orth_R(dim_V, dim_V) <= ORTHOG_TOL * pre_norm)
    {
      return false;
    }

    v *= 1.0 / orth_R(dim_V, dim_V);
    v_node_label.emplace_back(node_label);
    return true;
  };

  if (has_real && !add_real_vector_to_basis(u.Real(), fmt::format("{}_re", node_label)))
  {
    MFEM_ABORT("Linearly dependent vector added to PROM basis. This indicates a "
               "convergence issue or a code error (the same vector was added multiple "
               "times accidentally).");
  }
  if (has_imag && !add_real_vector_to_basis(u.Imag(), fmt::format("{}_im", node_label)))
  {
    MFEM_ABORT("Linearly dependent vector added to PROM basis. This indicates a "
               "convergence issue or a code error (the same vector was added multiple "
               "times accidentally).");
  }

  if (dim_V_new == dim_V_old)
  {
    return;
  }

  // Update reduced-order operators. Resize preserves the upper dim0 x dim0 block of each
  // matrix and first dim0 entries of each vector and the projection uses the values
  // computed for the unchanged basis vectors.
  Kr.conservativeResize(dim_V_new, dim_V_new);
  ProjectMatInternal(comm, V, *K, Kr, r, dim_V_old);
  if (C)
  {
    Cr.conservativeResize(dim_V_new, dim_V_new);
    ProjectMatInternal(comm, V, *C, Cr, r, dim_V_old);
  }
  Mr.conservativeResize(dim_V_new, dim_V_new);
  ProjectMatInternal(comm, V, *M, Mr, r, dim_V_old);
  // Per-port wave-port masses. M^(p)_r is initialized lazily so the map only contains
  // entries when the per-port HDM operator was non-null on this rank.
  for (auto &[port_idx, Mp_hdm] : Mwp_p)
  {
    auto &Mp_r = Mwp_p_r[port_idx];
    Mp_r.conservativeResize(dim_V_new, dim_V_new);
    ProjectMatInternal(comm, V, *Mp_hdm, Mp_r, r, dim_V_old);
  }
  if (RHS1.Size())
  {
    RHS1r.conservativeResize(dim_V_new);
    ProjectVecInternal(comm, V, RHS1, RHS1r, dim_V_old);
  }
  // Extend the cached projected RHS2 samples (one per offline (excitation, ω⋆) pair) to
  // cover the new basis vectors. The underlying HDM-side RHS2 sample is preserved so we
  // can compute the new entries Vᴴ·RHS2_hdm only for the newly added basis vectors.
  for (auto &[exc_idx, samples_r] : RHS2_r_samples)
  {
    auto &samples_hdm = RHS2_hdm_samples.at(exc_idx);
    for (auto &[omega_star, sample_r] : samples_r)
    {
      auto &sample_hdm = samples_hdm.at(omega_star);
      sample_r.conservativeResize(dim_V_new);
      ProjectVecInternal(comm, V, sample_hdm, sample_r, dim_V_old);
    }
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
  if (V.empty())
  {
    u = 0.0;
    return;
  }

  Ar.resize(V.size(), V.size());
  RHSr.resize(V.size());

  // Slow path: any remaining ω-nonlinear A2 contributors that we have not factored out.
  // Currently this is second-order farfield boundaries and surface conductivity. The
  // wave-port contribution is excluded here and applied below via the per-port factored
  // form.
  if (has_other_A2)
  {
    A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO,
                                                        /*include_wave_ports=*/false);
    if (A2)
    {
      ProjectMatInternal(space_op.GetComm(), V, *A2, Ar, r, 0);
    }
    else
    {
      Ar.setZero();
    }
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
  // Wave-port contribution: A_wp(ω) = i·Σ_p kₙ,p(ω)·M^(p)_{μ⁻¹}. The Mwp_p HDM
  // operators are purely imaginary (constructed with the boundary mass in the imaginary
  // slot), so multiplying their projection by the real scalar kₙ recovers the i·kₙ·M
  // contribution to the system matrix. The cross-section EVP for kₙ,p(ω) is small
  // relative to a full HDM solve and is the only ω-dependent work here.
  for (const auto &[port_idx, Mp_r] : Mwp_p_r)
  {
    double kn = space_op.GetWavePortOp().GetWavePortKn(port_idx, omega);
    Ar += std::complex<double>(kn, 0.0) * Mp_r;
  }

  // Wave-port (or any other RHS2) excitation: interpolate the cached projected
  // RHS2_r(ω⋆) samples in ω rather than re-assembling and re-projecting RHS2 from the
  // full HDM at each online frequency. The samples were captured during the offline
  // greedy SolveHDM calls and rolled forward in UpdatePROM. Empty cache means RHS2 was
  // zero at every offline sample, so it's zero here too.
  auto rhs2_it = RHS2_r_samples.find(excitation_idx);
  if (rhs2_it != RHS2_r_samples.end() && !rhs2_it->second.empty())
  {
    auto rhs2_interp = BarycentricInterpolate(rhs2_it->second, omega);
    MFEM_VERIFY(rhs2_interp.size() == static_cast<long>(V.size()),
                "Cached RHS2 projection size does not match the current PROM basis "
                "size. UpdatePROM is responsible for extending all cached samples "
                "whenever the basis grows; this indicates a bookkeeping bug.");
    RHSr = rhs2_interp;
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
    // QR solve, for maximal stability. The small system is cheap to compute but can be
    // numerically poorly conditioned to due the splitting of HDM solutions into Re and Im
    // into separate columns.
    RHSr = Ar.fullPivHouseholderQr().solve(RHSr);
  }
  ProlongatePROMSolution(V.size(), V, RHSr, u);
}

std::vector<std::complex<double>> RomOperator::ComputeEigenvalueEstimates() const
{
  // Newton fixed-point eigenvalues of the projected pencil with **exact** wave-port
  // dispersion:
  //
  //   (Kr + iω Cr − ω² Mr + Σₚ kₙ,ₚ(ω) · Mp_r) v = 0,
  //
  // computed by seeding from the polynomial-fit pencil (companion-form QEP) and
  // refining each seed with the dispersion correction (kₙ_exact − kₙ_poly)·Mp_r
  // evaluated at the current iterate. Returned ω values are in Palace internal
  // nondimensional units (rad / tc).
  //
  // Used as a self-test for the regime-2 augmented L⁻¹/R⁻¹/C realisation:
  // eigenvalues of the augmented pencil (computed externally from the printed
  // matrices) should match these to machine precision when the AAA fit residual is
  // below tolerance. Mismatch → bug in the augmentation, not the underlying
  // physics. For augmentation-disabled (regime 1) configurations, the seed and
  // refined values can disagree by the polynomial fit residual — that's expected.
  if (Kr.rows() == 0)
  {
    return {};
  }
  const long n = Kr.rows();

  // Per-port polynomial fit, identical to CalculateNormalizedPROMMatrices.
  std::map<int, Eigen::Vector3d> poly_coeffs;
  std::map<int, std::pair<double, double>> kn2_coeffs;  // (γ, α_kn²)
  if (sweep_omega_max > sweep_omega_min && !Mwp_p_r.empty())
  {
    constexpr int n_fit = 30;
    Eigen::VectorXd ws(n_fit);
    Eigen::MatrixXd vand(n_fit, 3), vand_kn2(n_fit, 2);
    const double w_lo = sweep_omega_min;
    const double w_hi = sweep_omega_max;
    double w_max = 0.0;
    for (int i = 0; i < n_fit; i++)
    {
      ws(i) = w_lo + (w_hi - w_lo) * i / (n_fit - 1.0);
      vand(i, 0) = 1.0;
      vand(i, 1) = ws(i);
      vand(i, 2) = ws(i) * ws(i);
      w_max = std::max(w_max, std::abs(ws(i)));
    }
    for (int i = 0; i < n_fit; i++)
    {
      vand_kn2(i, 0) = 1.0;
      vand_kn2(i, 1) = (ws(i) / w_max) * (ws(i) / w_max);
    }
    auto qr = vand.colPivHouseholderQr();
    auto qr_kn2 = vand_kn2.colPivHouseholderQr();
    for (auto &[port_idx, Mp_r] : Mwp_p_r)
    {
      Eigen::VectorXd ys(n_fit), ys2(n_fit);
      for (int i = 0; i < n_fit; i++)
      {
        ys(i) = space_op.GetWavePortOp().GetWavePortKn(port_idx, ws(i));
        ys2(i) = ys(i) * ys(i);
      }
      poly_coeffs[port_idx] = qr.solve(ys);
      Eigen::Vector2d c2 = qr_kn2.solve(ys2);
      kn2_coeffs[port_idx] = {c2(1) / (w_max * w_max), c2(0)};
    }
  }

  // Polynomial-pencil matrices.
  Eigen::MatrixXcd Kr_eff = Kr;
  Eigen::MatrixXcd Cr_eff = (Cr.rows() == n) ? Cr : Eigen::MatrixXcd::Zero(n, n);
  Eigen::MatrixXcd Mr_eff = Mr;
  for (auto &[port_idx, Mp_r] : Mwp_p_r)
  {
    const auto &c = poly_coeffs.at(port_idx);
    Kr_eff += std::complex<double>(c(0), 0.0) * Mp_r;
    Cr_eff += std::complex<double>(0.0, -c(1)) * Mp_r;
    Mr_eff += std::complex<double>(-c(2), 0.0) * Mp_r;
  }

  auto solve_qep = [&](const Eigen::MatrixXcd &K_in)
  {
    Eigen::MatrixXcd A = Eigen::MatrixXcd::Zero(2 * n, 2 * n);
    Eigen::MatrixXcd B = Eigen::MatrixXcd::Zero(2 * n, 2 * n);
    A.topRightCorner(n, n) = Eigen::MatrixXcd::Identity(n, n);
    A.bottomLeftCorner(n, n) = K_in;
    A.bottomRightCorner(n, n) = std::complex<double>(0.0, 1.0) * Cr_eff;
    B.topLeftCorner(n, n) = Eigen::MatrixXcd::Identity(n, n);
    B.bottomRightCorner(n, n) = Mr_eff;
    Eigen::VectorXcd D;
    Eigen::MatrixXcd VR;
    ZGGEV(A, B, D, VR);
    return D;
  };

  Eigen::VectorXcd seeds = solve_qep(Kr_eff);

  std::vector<std::complex<double>> refined;
  if (Mwp_p_r.empty())
  {
    for (int i = 0; i < seeds.size(); i++)
    {
      if (std::isfinite(seeds(i).real()) && std::isfinite(seeds(i).imag()) &&
          seeds(i).real() > 0.0)
      {
        refined.push_back(seeds(i));
      }
    }
  }
  else
  {
    constexpr int max_iter = 30;
    constexpr double rel_tol = 1.0e-12;
    for (int s = 0; s < seeds.size(); s++)
    {
      std::complex<double> w = seeds(s);
      if (!(std::isfinite(w.real()) && std::isfinite(w.imag()) && w.real() > 0.0))
      {
        continue;
      }
      bool converged = false;
      for (int it = 0; it < max_iter; it++)
      {
        Eigen::MatrixXcd K_corr = Kr_eff;
        for (auto &[port_idx, Mp_r] : Mwp_p_r)
        {
          // kₙ from the canonical kₙ² fit, evaluated at complex ω so the iteration
          // works above and below cutoff (where kₙ becomes imaginary).
          const auto &[gamma, alpha_kn2] = kn2_coeffs.at(port_idx);
          std::complex<double> kn_sq = gamma * w * w + alpha_kn2;
          std::complex<double> kn_exact = std::sqrt(kn_sq);
          const auto &c = poly_coeffs.at(port_idx);
          std::complex<double> kn_poly = c(0) + c(1) * w + c(2) * w * w;
          K_corr += (kn_exact - kn_poly) * Mp_r;
        }
        Eigen::VectorXcd D = solve_qep(K_corr);
        std::complex<double> best = w;
        double best_dist = std::numeric_limits<double>::infinity();
        for (int i = 0; i < D.size(); i++)
        {
          if (!(std::isfinite(D(i).real()) && std::isfinite(D(i).imag()) &&
                D(i).real() > 0.0))
          {
            continue;
          }
          double d = std::abs(D(i) - w);
          if (d < best_dist)
          {
            best_dist = d;
            best = D(i);
          }
        }
        double rel = std::abs(best - w) / std::max(std::abs(w), 1.0);
        w = best;
        if (rel < rel_tol)
        {
          converged = true;
          break;
        }
      }
      if (converged)
      {
        bool duplicate = false;
        for (const auto &existing : refined)
        {
          if (std::abs(existing - w) / std::max(std::abs(w), 1.0) < 1.0e-9)
          {
            duplicate = true;
            break;
          }
        }
        if (!duplicate)
        {
          refined.push_back(w);
        }
      }
    }
  }

  std::sort(refined.begin(), refined.end(),
            [](const std::complex<double> &a, const std::complex<double> &b)
            { return a.real() < b.real(); });
  return refined;
}

RomOperator::WavePortRegime
RomOperator::SelectWavePortRegime(int port_idx, double rel_err, bool meets_tol) const
{
  // AUTO: polynomial if residual meets tolerance, else augmented. Force settings
  // override but emit a warning when the user-requested regime is incompatible with
  // the computed residual.
  switch (waveport_synthesis_force)
  {
    case WavePortSynthesisRegime::POLYNOMIAL:
      if (!meets_tol)
      {
        Mpi::Warning(
            "Wave port {:d}: WavePortSynthesisForce=Polynomial, but the order-2 fit "
            "residual {:.3e} exceeds WavePortSynthesisTol={:.3e}. Proceeding with "
            "the larger error per user request.\n",
            port_idx, rel_err, waveport_synthesis_tol);
      }
      return WavePortRegime::Polynomial;
    case WavePortSynthesisRegime::AUGMENTED:
      return WavePortRegime::Augmented;
    case WavePortSynthesisRegime::AUTO:
    default:
      return meets_tol ? WavePortRegime::Polynomial : WavePortRegime::Augmented;
  }
}

RomOperator::WavePortDispersionFit
RomOperator::FitWavePortDispersion(int port_idx, const Eigen::MatrixXcd &Mp_r) const
{
  // Sample kₙ,p on the sweep band, fit a quadratic, optionally augment with AAA on the
  // residual. Returns a packed result describing the chosen regime and the fit data.
  WavePortDispersionFit fit;
  fit.port_idx = port_idx;

  constexpr int n_fit = 30;
  constexpr int n_dense = 1000;
  const double w_lo = sweep_omega_min;
  const double w_hi = sweep_omega_max;
  auto sample_omega = [w_lo, w_hi](int n)
  {
    std::vector<double> ws(n);
    for (int i = 0; i < n; i++)
    {
      ws[i] = w_lo + (w_hi - w_lo) * i / std::max(n - 1, 1);
    }
    return ws;
  };
  auto fit_omegas = sample_omega(n_fit);
  auto dense_omegas = sample_omega(n_dense);

  // Sample kₙ,p on the fit grid. Each call triggers (or hits the cache of)
  // WavePortData::Initialize. Cheap because the cross-section EVP is small.
  Eigen::VectorXd y_fit(n_fit);
  for (int i = 0; i < n_fit; i++)
  {
    y_fit(i) = space_op.GetWavePortOp().GetWavePortKn(port_idx, fit_omegas[i]);
  }
  // LSQ polynomial fit at order 2 in ω. Higher orders cannot be absorbed into the
  // K + iωC − ω²M synthesis structure (cf. design notes).
  Eigen::MatrixXd vandermonde(n_fit, 3);
  for (int i = 0; i < n_fit; i++)
  {
    const double w = fit_omegas[i];
    vandermonde(i, 0) = 1.0;
    vandermonde(i, 1) = w;
    vandermonde(i, 2) = w * w;
  }
  Eigen::Vector3d coeffs = vandermonde.colPivHouseholderQr().solve(y_fit);
  fit.alpha0 = coeffs(0);
  fit.alpha1 = coeffs(1);
  fit.alpha2 = coeffs(2);

  // Residual on the dense grid: relative error on kₙ. The 1/(ωμ) factor in the
  // wave-port admittance Y_p(ω) = kₙ,p(ω)/(iωμ) cancels in the ratio, so kₙ is the
  // right proxy for the user-visible synthesis accuracy (no HDM data required).
  double max_rel = 0.0, max_abs_truth = 0.0;
  for (int i = 0; i < n_dense; i++)
  {
    const double w = dense_omegas[i];
    const double truth = space_op.GetWavePortOp().GetWavePortKn(port_idx, w);
    const double poly = fit.alpha0 + fit.alpha1 * w + fit.alpha2 * w * w;
    max_abs_truth = std::max(max_abs_truth, std::abs(truth));
    max_rel = std::max(max_rel, std::abs(poly - truth));
  }
  fit.rel_err_polynomial = (max_abs_truth > 0.0) ? max_rel / max_abs_truth : 0.0;
  fit.rel_err_augmented = fit.rel_err_polynomial;
  const bool meets_tol = (fit.rel_err_polynomial <= waveport_synthesis_tol);

  fit.regime = SelectWavePortRegime(port_idx, fit.rel_err_polynomial, meets_tol);

  if (fit.regime == WavePortRegime::Polynomial)
  {
    if (meets_tol)
    {
      Mpi::Print(
          " Wave port {:d}: polynomial synthesis residual {:.3e} (tol {:.3e}, α₀={:.3e}, "
          "α₁={:.3e}, α₂={:.3e})\n",
          port_idx, fit.rel_err_polynomial, waveport_synthesis_tol, fit.alpha0,
          fit.alpha1, fit.alpha2);
    }
    return fit;
  }

  // Augmented regime: AAA rational fit on the polynomial residual δkₙ(ω). Cap the
  // pole count at waveport_synthesis_order_max.
  Eigen::VectorXcd z_aaa(n_fit), F_aaa(n_fit);
  for (int i = 0; i < n_fit; i++)
  {
    const double w = fit_omegas[i];
    z_aaa(i) = w;
    F_aaa(i) = y_fit(i) - (fit.alpha0 + fit.alpha1 * w + fit.alpha2 * w * w);
  }
  const double aaa_tol_rel = (max_abs_truth > 0.0)
                                 ? waveport_synthesis_tol * max_abs_truth /
                                       std::max(F_aaa.cwiseAbs().maxCoeff(), 1.0e-300)
                                 : waveport_synthesis_tol;
  auto aaa = utils::RunAAA(z_aaa, F_aaa, aaa_tol_rel,
                           std::max<std::size_t>(waveport_synthesis_order_max, 1));
  auto pr = utils::AAAToPoleResidue(aaa);

  // Fold the asymptote d into α₀ — equivalent and avoids carrying a separate constant
  // offset through the augmented block.
  const double aaa_d = pr.d.real();
  fit.alpha0 += aaa_d;

  // Re-evaluate residual on the dense grid using the augmented model.
  max_rel = 0.0;
  for (int i = 0; i < n_dense; i++)
  {
    const double w = dense_omegas[i];
    const double truth = space_op.GetWavePortOp().GetWavePortKn(port_idx, w);
    double aug = fit.alpha0 + fit.alpha1 * w + fit.alpha2 * w * w;
    for (long k = 0; k < pr.poles.size(); k++)
    {
      std::complex<double> wc = w;
      aug += (pr.residues(k) / (wc - pr.poles(k))).real();
    }
    max_rel = std::max(max_rel, std::abs(aug - truth));
  }
  fit.rel_err_augmented = (max_abs_truth > 0.0) ? max_rel / max_abs_truth : 0.0;

  // SVD of M_proj — keep every singular vector with σⱼ / σ_max above the rank
  // tolerance as a coupling direction. Rank is typically 1 (real-only mode field) or
  // 2 (complex mode field split into real+imag basis vectors).
  WavePortAuxBlock blk;
  blk.port_idx = port_idx;
  Eigen::MatrixXd M_proj = Mp_r.imag().eval();
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(M_proj, Eigen::ComputeThinU);
  const double sigma_max = svd.singularValues()(0);
  for (long j = 0; j < svd.singularValues().size(); j++)
  {
    double s = svd.singularValues()(j);
    if (sigma_max > 0.0 && s / sigma_max > waveport_synthesis_rank_tol)
    {
      blk.sigmas.push_back(s);
      blk.u_dirs.push_back(svd.matrixU().col(j));
    }
  }
  for (long k = 0; k < pr.poles.size(); k++)
  {
    blk.poles.push_back(pr.poles(k));
    blk.residues.push_back(pr.residues(k));
  }
  std::size_t aux_per_port = pr.poles.size() * blk.sigmas.size();
  std::size_t rank_used = blk.sigmas.size();
  Mpi::Print(" Wave port {:d}: augmented synthesis residual {:.3e} → {:.3e} "
             "(tol {:.3e}, {:d} pole{} × rank-{:d} mass = +{:d} aux state{}, "
             "α₀={:.3e}, α₁={:.3e}, α₂={:.3e}, d={:.3e})\n",
             port_idx, fit.rel_err_polynomial, fit.rel_err_augmented,
             waveport_synthesis_tol, pr.poles.size(), pr.poles.size() == 1 ? "" : "s",
             rank_used, aux_per_port, aux_per_port == 1 ? "" : "s", fit.alpha0,
             fit.alpha1, fit.alpha2, aaa_d);
  fit.aux = std::move(blk);
  return fit;
}

void RomOperator::ApplyPolynomialFitCorrections(const WavePortDispersionFit &fit,
                                                const Eigen::MatrixXcd &Mp_r,
                                                Eigen::MatrixXcd &Kr_corr,
                                                Eigen::MatrixXcd &Cr_corr,
                                                Eigen::MatrixXcd &Mr_corr)
{
  // Mp_r is purely imaginary (= i·M_proj). Multiplying by α₀ folds α₀·M_proj into the
  // imaginary part of Kr (which becomes Im(L⁻¹)); the i·α₁·M_proj term lives in the
  // real part of Cr (Re(R⁻¹)) and -α₂·M_proj·v adds to the imaginary part of Mr
  // (Im(C)). Signs are set so the sum recovers i·(α₀+α₁ω+α₂ω²)·M_proj·v in Aᵣ(ω).
  Kr_corr += std::complex<double>(fit.alpha0, 0.0) * Mp_r;
  Cr_corr += std::complex<double>(0.0, -fit.alpha1) * Mp_r;
  Mr_corr += std::complex<double>(-fit.alpha2, 0.0) * Mp_r;
}

RomOperator::AugmentedPencil RomOperator::BuildAugmentedPencil(
    const Eigen::MatrixXcd &Kr_total, const Eigen::MatrixXcd &Cr_total,
    const Eigen::MatrixXcd &Mr_total, const std::vector<WavePortAuxBlock> &aux_blocks,
    std::vector<std::string> &aux_labels)
{
  // Per pole–residue (pₖ, rₖ) at port p, with projected mass M_proj = Σⱼ σⱼ uⱼuⱼᵀ:
  // For each kept singular vector j, allocate one aux state sₖⱼ. The augmented pencil
  //   ⎡ K  K_va ⎤   ⎡ 0  0  ⎤    ⎡ 0  0 ⎤
  //   ⎢       ⎥+iω⎢      ⎥−ω²⎢     ⎥
  //   ⎣ K_avᵀ K_aa⎦   ⎣ 0 C_aa⎦    ⎣ 0  0 ⎦
  // with Kr_aa = -pₖ, Cr_aa = -i, K_va = a·uⱼ, K_avᵀ = a·uⱼᵀ symmetric. Eliminating
  // sₖⱼ gives the pencil contribution -a²·uⱼuⱼᵀ·v/(ω - pₖ). Choosing a² = -i·rₖ·σⱼ
  // makes the sum over j equal +i·rₖ·M_proj·v/(ω - pₖ), the desired rₖ-residue
  // contribution to the unaugmented wave-port pencil i·kₙ(ω)·Mp_r·v. Both K_va and
  // K_avᵀ get the same coupling so the augmented pencil is complex-symmetric (not
  // Hermitian — downstream eigensolvers handle this fine).
  std::size_t n_aux_total = 0;
  for (const auto &blk : aux_blocks)
  {
    n_aux_total += blk.poles.size() * blk.sigmas.size();
  }
  const long n_v = Kr_total.rows();
  const long n_aug = n_v + static_cast<long>(n_aux_total);
  AugmentedPencil aug;
  aug.Kr = Eigen::MatrixXcd::Zero(n_aug, n_aug);
  aug.Cr = Eigen::MatrixXcd::Zero(n_aug, n_aug);
  aug.Mr = Eigen::MatrixXcd::Zero(n_aug, n_aug);
  aug.Kr.topLeftCorner(n_v, n_v) = Kr_total;
  aug.Cr.topLeftCorner(n_v, n_v) = Cr_total;
  aug.Mr.topLeftCorner(n_v, n_v) = Mr_total;
  long aux_row = n_v;
  for (const auto &blk : aux_blocks)
  {
    for (std::size_t k = 0; k < blk.poles.size(); k++)
    {
      auto pk = blk.poles[k];
      auto rk = blk.residues[k];
      for (std::size_t j = 0; j < blk.sigmas.size(); j++)
      {
        std::complex<double> coupling =
            std::sqrt(std::complex<double>(0.0, -1.0) * rk * blk.sigmas[j]);
        aug.Kr(aux_row, aux_row) = -pk;
        aug.Cr(aux_row, aux_row) = std::complex<double>(0.0, -1.0);
        for (long i = 0; i < n_v; i++)
        {
          aug.Kr(i, aux_row) = coupling * blk.u_dirs[j](i);
          aug.Kr(aux_row, i) = coupling * blk.u_dirs[j](i);
        }
        aux_labels.push_back(
            fmt::format("waveport_{:d}_p{:d}d{:d}", blk.port_idx, k, j));
        aux_row++;
      }
    }
  }
  return aug;
}

RomOperator::NormalizedMatrices
RomOperator::CalculateNormalizedPROMMatrices(const Units &units) const
{
  using mat_t = Eigen::MatrixXcd;
  NormalizedMatrices out;

  // Port-row/column scaling. The port block of the synthesised admittance matrices is
  // independent of HDM solutions; for HDM rows we leave v_conc = 1 (orth_R(j,j) can be
  // tiny for nearly degenerate vectors and would make the matrices ill-conditioned).
  // Lumped ports come first in the basis (real-valued, in order); wave-port basis
  // vectors follow (real and imag parts each count separately, see
  // AddWavePortModesForSynthesis → UpdatePROM).
  Eigen::VectorXd v_conc = Eigen::VectorXd::Ones(GetReducedDimension());

  // Only ports with include_in_synthesis contribute a basis row, so the number of leading
  // lumped-port rows is NumSynthesisPortModes(), NOT the total number of lumped ports.
  const long n_port_modes = static_cast<long>(NumSynthesisPortModes());
  MFEM_ASSERT(n_port_modes <= GetReducedDimension(),
              "More lumped port modes than PROM basis vectors; basis is inconsistent!");
  // Wave-port basis vectors come right after the lumped-port rows; count how many were
  // actually added (real + imaginary parts each count separately, see
  // AddWavePortModesForSynthesis → UpdatePROM).
  long n_waveport_rows = 0;
  if (!Mwp_p_r.empty())
  {
    for (long j = n_port_modes; j < static_cast<long>(v_node_label.size()); j++)
    {
      if (v_node_label[j].rfind("waveport_", 0) == 0)
      {
        n_waveport_rows++;
      }
      else
      {
        break;
      }
    }
  }
  for (long j = 0; j < n_port_modes + n_waveport_rows; j++)
  {
    v_conc[j] = orth_R(j, j);
  }

  // Wave-port dispersion handling. Online, SolvePROM assembles
  //   Ar(ω) = Kr + iω Cr − ω² Mr + Σₚ kₙ,ₚ(ω)·Mp_r,
  // where Mp_r is purely imaginary. For circuit synthesis we replace kₙ,ₚ(ω) with an
  // approximation that fits inside a quadratic-in-ω pencil:
  //
  //   • Polynomial regime: kₙ,ₚ(ω) ≈ α₀ + α₁ω + α₂ω², absorbed into Kr/Cr/Mr.
  //   • Augmented regime: residual δkₙ(ω) = kₙ(ω) − (α₀+α₁ω+α₂ω²) fit by an AAA
  //     rational expansion d + Σₖ rₖ/(ω − pₖ); d folds into α₀, each pole–residue pair
  //     becomes auxiliary scalar states appended as new rows/columns to L⁻¹/R⁻¹/C.
  //
  // FitWavePortDispersion does the per-port sampling, fitting and regime selection;
  // ApplyPolynomialFitCorrections accumulates the α-contributions; BuildAugmentedPencil
  // extends the pencil with the aux states from collected fit results.
  Eigen::MatrixXcd Kr_corr = Eigen::MatrixXcd::Zero(Kr.rows(), Kr.cols());
  Eigen::MatrixXcd Cr_corr = Eigen::MatrixXcd::Zero(Kr.rows(), Kr.cols());
  Eigen::MatrixXcd Mr_corr = Eigen::MatrixXcd::Zero(Kr.rows(), Kr.cols());
  std::vector<WavePortAuxBlock> aux_blocks;
  if (!Mwp_p_r.empty() && sweep_omega_max > sweep_omega_min)
  {
    for (auto &[port_idx, Mp_r] : Mwp_p_r)
    {
      auto fit = FitWavePortDispersion(port_idx, Mp_r);
      ApplyPolynomialFitCorrections(fit, Mp_r, Kr_corr, Cr_corr, Mr_corr);
      if (fit.aux)
      {
        aux_blocks.push_back(std::move(*fit.aux));
      }
    }
  }

  // Polynomial-only matrices (basis dim n × n).
  Eigen::MatrixXcd Kr_total = Kr + Kr_corr;
  Eigen::MatrixXcd Mr_total = Mr + Mr_corr;
  Eigen::MatrixXcd Cr_total = Cr_corr;
  if (C)
  {
    Cr_total += Cr;
  }

  auto aug = BuildAugmentedPencil(Kr_total, Cr_total, Mr_total, aux_blocks,
                                  out.aux_labels);

  // v_d port-row scaling: extend with 1's for aux rows (no port-impedance scaling on
  // aux states — they're internal circuit nodes).
  const long n_v = Kr_total.rows();
  const long n_aug = aug.Kr.rows();
  Eigen::VectorXcd v_conc_aug = Eigen::VectorXcd::Ones(n_aug);
  for (long j = 0; j < n_v; j++)
  {
    v_conc_aug(j) = v_conc(j);
  }
  auto v_d_aug = v_conc_aug.asDiagonal();

  auto unit_henry_inv = 1.0 / units.GetScaleFactor<Units::ValueType::INDUCTANCE>();
  out.L_inv =
      std::make_unique<mat_t>((unit_henry_inv * v_d_aug * aug.Kr * v_d_aug).eval());

  auto unit_farad = units.GetScaleFactor<Units::ValueType::CAPACITANCE>();
  out.C = std::make_unique<mat_t>((unit_farad * v_d_aug * aug.Mr * v_d_aug).eval());

  // Emit R⁻¹ whenever there's any dissipative contribution: lumped resistance, surface
  // conductivity, wave-port α₁, or the aux state diagonal (regime 2).
  const bool has_R_inv = (aug.Cr.cwiseAbs().maxCoeff() > 0.0);
  if (has_R_inv)
  {
    auto unit_ohm_inv = 1.0 / units.GetScaleFactor<Units::ValueType::IMPEDANCE>();
    out.R_inv =
        std::make_unique<mat_t>((unit_ohm_inv * v_d_aug * aug.Cr * v_d_aug).eval());
  }

  return out;
}

void RomOperator::PrintPROMMatrices(const Units &units, const fs::path &post_dir) const
{
  BlockTimer bt0(Timer::POSTPRO);
  Mpi::Print(" Printing PROM Matrices to disk.\n");

  // Build the synthesised matrices on every rank, since the polynomial-fit step inside
  // calls into WavePortOperator::GetWavePortKn which triggers the cross-section EVP —
  // a collective operation that would deadlock if only the root participated. The
  // resulting matrices are replicated; only the root will write them to disk below.
  auto matrices = CalculateNormalizedPROMMatrices(units);
  const auto &inductance_L_inv = matrices.L_inv;
  const auto &resistance_R_inv = matrices.R_inv;
  const auto &capacitance_C = matrices.C;

  if (!Mpi::Root(space_op.GetComm()))
  {
    return;
  }
  // Row/column labels for the augmented matrices: basis labels followed by any
  // regime-2 aux state labels. The synthesised matrices have one row per basis
  // vector + one row per aux state (waveport_<idx>_aux_<k>).
  std::vector<std::string> labels = v_node_label;
  for (const auto &lab : matrices.aux_labels)
  {
    labels.push_back(lab);
  }
  auto print_table =
      [post_dir, labels](const Eigen::MatrixXd &mat, std::string_view filename)
  {
    MFEM_VERIFY((labels.size() == mat.cols()) && (labels.size() == mat.rows()),
                "Inconsistent PROM size!");

    auto out = TableWithCSVFile(post_dir / filename);
    out.table.col_options.float_precision = 17;
    for (long i = 0; i < mat.cols(); i++)
    {
      out.table.insert(labels[i], labels[i]);
      auto &col = out.table[i];
      for (long j = 0; j < mat.rows(); j++)
      {
        col << mat(j, i);
      }
    }
    out.WriteFullTableTrunc();
  };

  // For each synthesised matrix, emit the real part if its content is nonzero, and
  // similarly for the imaginary part. The HDM nullptr (e.g. K->Imag()) is not a
  // reliable proxy: the wave-port polynomial-fit corrections inject α₀ into Im(L⁻¹)
  // and α₂ into Im(C) even when the underlying HDM K/M have no imaginary part, and
  // wave-port dissipation fills R⁻¹ even when no Palace damping matrix exists. Gate
  // exclusively on the synthesised content so nothing is silently dropped.
  auto print_if_nonzero = [&](const Eigen::MatrixXd &mat, std::string_view filename)
  {
    if (mat.cwiseAbs().maxCoeff() > 0.0)
    {
      print_table(mat, filename);
    }
  };
  print_if_nonzero(inductance_L_inv->real(), "rom-Linv-re.csv");
  print_if_nonzero(inductance_L_inv->imag(), "rom-Linv-im.csv");
  print_if_nonzero(capacitance_C->real(), "rom-C-re.csv");
  print_if_nonzero(capacitance_C->imag(), "rom-C-im.csv");
  if (resistance_R_inv)
  {
    print_if_nonzero(resistance_R_inv->real(), "rom-Rinv-re.csv");
    print_if_nonzero(resistance_R_inv->imag(), "rom-Rinv-im.csv");
  }

  // Print orth-R. Don't divide by diagonal to keep state normalization info.
  // Pad with identity for aux states (regime 2): aux rows are synthetic circuit
  // nodes with no basis-orthogonality content; identity preserves the diagonal
  // form expected by downstream consumers.
  Eigen::MatrixXd orth_R_padded = Eigen::MatrixXd::Identity(labels.size(), labels.size());
  orth_R_padded.topLeftCorner(orth_R.rows(), orth_R.cols()) = orth_R;
  print_table(orth_R_padded, "rom-orthogonalization-matrix-R.csv");
}

}  // namespace palace
