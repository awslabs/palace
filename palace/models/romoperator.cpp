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
    // Wave-port modes are added once per port (one mode per port today). The seeded
    // basis vector at the reference frequency is generally complex (mode field has both
    // real and imaginary parts), so reserve up to two slots per port.
    max_prom_size += 2 * space_op.GetWavePortOp().Size();

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
  auto P = space_op.GetPreconditionerMatrix<ComplexOperator>(
      1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i, omega + 0.0i);
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
    space_op.GetLumpedPortExcitationVectorPrimaryEt(port_idx, vec, true);
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
    space_op.GetWavePortFieldVectorPrimaryEt(port_idx, omega_ref, vec, true);
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
    if (rhs2_interp.size() == static_cast<long>(V.size()))
    {
      RHSr = rhs2_interp;
    }
    else
    {
      // Stale (basis grew since this excitation was last touched in UpdatePROM and the
      // cached projection wasn't extended). This shouldn't happen in normal flow but
      // fall back to a safe direct compute.
      space_op.GetExcitationVector2(excitation_idx, omega, RHS2);
      ProjectVecInternal(space_op.GetComm(), V, RHS2, RHSr, 0);
    }
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
  // XX TODO: Not yet implemented
  MFEM_ABORT("Eigenvalue estimates for PROM operators are not yet implemented!");
  return {};
}

std::tuple<std::unique_ptr<Eigen::MatrixXcd>, std::unique_ptr<Eigen::MatrixXcd>,
           std::unique_ptr<Eigen::MatrixXcd>>
RomOperator::CalculateNormalizedPROMMatrices(const Units &units) const
{
  using mat_t = Eigen::MatrixXcd;
  std::unique_ptr<mat_t> inductance_L_inv = {};
  std::unique_ptr<mat_t> resistance_R_inv = {};
  std::unique_ptr<mat_t> capacitance_C = {};

  Eigen::VectorXd v_conc = Eigen::VectorXd::Ones(GetReducedDimension());

  // We will only scale the rows and columns corresponding to ports by orth_R(j,j). This
  // will correctly give back the port-port block of the circuit matrices, which is
  // independent of the HDM solutions. The normalization of the rows and columns of the HDM
  // solutions (synthesized nodes) are not relevant for port quantities. We choose v_conc =
  // 1.0, since for nearly degenerate vectors orth_R(j,j) could be tiny, which may lead to
  // poor condition of the circuit matrices.
  //
  // Lumped ports are real, at the beginning and in order. Only ports with
  // include_in_synthesis contribute a basis row, so the number of leading lumped-port
  // rows is NumSynthesisPortModes(), NOT the total number of lumped ports. Wave-port
  // modes follow the lumped ports in the basis (added by AddWavePortModesForSynthesis);
  // they may have both real and imaginary parts. The synthesised matrices are scaled
  // symmetrically by orth_R for those rows as well so the port-block diagonal entries
  // take their physical port-impedance value.
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

  // Lazy diagonal representation
  auto v_d = v_conc.asDiagonal();

  // Polynomial-fit corrections for wave-port dispersion. The wave-port contribution
  // i·kₙ,p(ω)·Mp_r enters the system matrix; its polynomial expansion
  // kₙ,p(ω) ≈ α₀ + α₁ω + α₂ω² is absorbed into Kr/Cr/Mr as documented in the design
  // notes (cf. prom-waveport.md):
  //   i α₀ Mp_r           → Im(Kr)        → Im(L⁻¹)
  //   α₁ Mp_r             → Re(Cr)        → Re(R⁻¹)
  //   −i α₂ Mp_r          → Im(Mr)        → Im(C)
  // (The signs follow Ar = Kr + iωCr − ω²Mr.)
  // Higher-order terms cannot be absorbed and require the augmented state space (regime
  // 2, future). For now: fit at order ≤ 2 and warn if the residual exceeds tolerance.
  Eigen::MatrixXcd Kr_corr = Eigen::MatrixXcd::Zero(Kr.rows(), Kr.cols());
  Eigen::MatrixXcd Cr_corr = Eigen::MatrixXcd::Zero(Kr.rows(), Kr.cols());
  Eigen::MatrixXcd Mr_corr = Eigen::MatrixXcd::Zero(Kr.rows(), Kr.cols());
  if (!Mwp_p_r.empty() && sweep_omega_max > sweep_omega_min)
  {
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

    const bool force_polynomial =
        (waveport_synthesis_force == WavePortSynthesisRegime::POLYNOMIAL);
    const bool force_augmented =
        (waveport_synthesis_force == WavePortSynthesisRegime::AUGMENTED);

    for (auto &[port_idx, Mp_r] : Mwp_p_r)
    {
      // Sample kₙ,p at fit and dense grids. Each call triggers (or hits the cache of)
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
      const double alpha0 = coeffs(0);
      const double alpha1 = coeffs(1);
      const double alpha2 = coeffs(2);

      // Evaluate residual on the dense grid: relative error on the per-port wave
      // admittance Y_p(ω) = kₙ,p(ω)/(iωμ). The 1/(ωμ) factor cancels in the ratio, so
      // we use kₙ directly. This is the cheap proxy for the user-visible synthesis
      // accuracy (no HDM data required).
      double max_rel = 0.0, max_abs_truth = 0.0;
      for (int i = 0; i < n_dense; i++)
      {
        const double w = dense_omegas[i];
        const double truth = space_op.GetWavePortOp().GetWavePortKn(port_idx, w);
        const double fit = alpha0 + alpha1 * w + alpha2 * w * w;
        max_abs_truth = std::max(max_abs_truth, std::abs(truth));
        max_rel = std::max(max_rel, std::abs(fit - truth));
      }
      const double rel_err = (max_abs_truth > 0.0) ? max_rel / max_abs_truth : 0.0;
      const bool meets_tol = (rel_err <= waveport_synthesis_tol);

      // Regime auto-dispatch. Today only regime 1 (polynomial absorption into L/R/C) is
      // implemented; regime 2 (augmented LC ladder) is reserved for cutoff-spanning
      // ports where the polynomial fit cannot meet the tolerance. Until regime 2 is
      // implemented, emit a clear warning when the polynomial path is the wrong choice
      // so the user knows their synthesis is approximate. The dispatch variable is
      // currently informational; once regime 2 lands it will gate the augmented-state
      // matrix construction below.
      WavePortSynthesisRegime regime_used = WavePortSynthesisRegime::POLYNOMIAL;
      (void)regime_used;  // suppress unused-variable until regime 2 lands
      if (force_augmented)
      {
        Mpi::Warning(
            "Wave port {:d}: WavePortSynthesisForce=Augmented requested but the "
            "augmented state-space (LC-ladder) synthesis is not yet implemented in this "
            "build. Falling back to polynomial fit (regime 1). Polynomial residual = "
            "{:.3e}, tol = {:.3e}.\n",
            port_idx, rel_err, waveport_synthesis_tol);
      }
      else if (force_polynomial && !meets_tol)
      {
        Mpi::Warning(
            "Wave port {:d}: WavePortSynthesisForce=Polynomial, but the order-2 fit "
            "residual {:.3e} exceeds WavePortSynthesisTol={:.3e}. Proceeding with the "
            "larger error per user request.\n",
            port_idx, rel_err, waveport_synthesis_tol);
      }
      else if (!meets_tol)
      {
        // AUTO regime, polynomial insufficient. Most likely a closed-waveguide port
        // operating near a modal cutoff (cf. prom-waveport.md "Case B"). The synthesis
        // output remains usable but accuracy is bounded by the residual reported here.
        Mpi::Warning(
            "Wave port {:d}: order-2 polynomial fit residual {:.3e} exceeds "
            "WavePortSynthesisTol={:.3e}. This typically indicates a closed-waveguide "
            "port operating near a modal cutoff. Augmented state-space (LC-ladder) "
            "synthesis would be needed for tighter accuracy but is not yet implemented; "
            "proceeding with the order-2 polynomial fit.\n",
            port_idx, rel_err, waveport_synthesis_tol);
      }
      else
      {
        Mpi::Print(" Wave port {:d}: polynomial synthesis residual {:.3e} (tol {:.3e}, "
                   "α₀={:.3e}, α₁={:.3e}, α₂={:.3e})\n",
                   port_idx, rel_err, waveport_synthesis_tol, alpha0, alpha1, alpha2);
      }

      // Apply absorbed contributions. Online, SolvePROM assembles
      //   Ar = Kr + iωCr − ω²Mr + Σ_p kₙ,p(ω) · Mp_r,
      // where Mp_r is purely imaginary (Re=0, Im = projected boundary mass M_proj).
      // Substituting kₙ,p(ω) ≈ α₀ + α₁ω + α₂ω² and matching ω-powers gives the
      // effective folded matrices:
      //   Kr_eff = Kr + α₀·Mp_r              (Im(Kr) gains α₀·M_proj)
      //   Cr_eff = Cr + (−i α₁)·Mp_r         (Re(Cr) gains α₁·M_proj)
      //   Mr_eff = Mr + (−α₂)·Mp_r           (Im(Mr) gains −α₂·M_proj)
      // The synthesis output then maps Kr→L⁻¹, Cr→R⁻¹, Mr→C with appropriate units.
      Kr_corr += std::complex<double>(alpha0, 0.0) * Mp_r;
      Cr_corr += std::complex<double>(0.0, -alpha1) * Mp_r;
      Mr_corr += std::complex<double>(-alpha2, 0.0) * Mp_r;
    }
  }

  auto unit_henry_inv = 1.0 / units.GetScaleFactor<Units::ValueType::INDUCTANCE>();
  inductance_L_inv =
      std::make_unique<mat_t>(((unit_henry_inv * v_d) * (Kr + Kr_corr) * v_d).eval());

  auto unit_farad = units.GetScaleFactor<Units::ValueType::CAPACITANCE>();
  capacitance_C =
      std::make_unique<mat_t>(((unit_farad * v_d) * (Mr + Mr_corr) * v_d).eval());

  // C & Cr are optional in UpdatePROM. In practice Cr exists whenever there is a
  // dissipative contribution: a lumped resistive port, surface conductivity, or — newly
  // — a wave port (its polynomial-fit α₁ term contributes to Cr_corr even when no
  // lumped C matrix exists). Emit R⁻¹ whenever any of these contributes.
  const bool has_R_inv = C || (Cr_corr.cwiseAbs().maxCoeff() > 0.0);
  if (has_R_inv)
  {
    auto unit_ohm_inv = 1.0 / units.GetScaleFactor<Units::ValueType::IMPEDANCE>();
    Eigen::MatrixXcd Cr_total = Cr_corr;
    if (C)
    {
      Cr_total += Cr;
    }
    resistance_R_inv =
        std::make_unique<mat_t>(((unit_ohm_inv * v_d) * Cr_total * v_d).eval());
  }

  return std::make_tuple(std::move(inductance_L_inv), std::move(resistance_R_inv),
                         std::move(capacitance_C));
}

void RomOperator::PrintPROMMatrices(const Units &units, const fs::path &post_dir) const
{
  BlockTimer bt0(Timer::POSTPRO);
  Mpi::Print(" Printing PROM Matrices to disk.\n");

  // Build the synthesised matrices on every rank, since the polynomial-fit step inside
  // calls into WavePortOperator::GetWavePortKn which triggers the cross-section EVP —
  // a collective operation that would deadlock if only the root participated. The
  // resulting matrices are replicated; only the root will write them to disk below.
  const auto [inductance_L_inv, resistance_R_inv, capacitance_C] =
      CalculateNormalizedPROMMatrices(units);

  if (!Mpi::Root(space_op.GetComm()))
  {
    return;
  }
  auto print_table = [post_dir, labels = this->v_node_label](const Eigen::MatrixXd &mat,
                                                             std::string_view filename)
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

  // Note: When checking for imaginary parts, it is better to do this for K,C,M as this is
  // a nullptr check. Kr, Mr, Cr would require a numerical check on imag elements.

  if (K->Real())
  {
    print_table(inductance_L_inv->real(), "rom-Linv-re.csv");
  }
  if (K->Imag())
  {
    print_table(inductance_L_inv->imag(), "rom-Linv-im.csv");
  }

  if (M->Real())
  {
    print_table(capacitance_C->real(), "rom-C-re.csv");
  }
  if (M->Imag())
  {
    print_table(capacitance_C->imag(), "rom-C-im.csv");
  }

  // R⁻¹ is emitted whenever there is a dissipative contribution (lumped resistance,
  // surface conductivity, or wave-port polynomial-fit α₁). When `C` (the Palace damping
  // matrix) is null, we still emit if `resistance_R_inv` is non-null — this is the case
  // for wave-port-only systems where the dissipation comes through Cr_corr.
  if (resistance_R_inv)
  {
    const bool has_real = (resistance_R_inv->real().cwiseAbs().maxCoeff() > 0.0);
    const bool has_imag = (resistance_R_inv->imag().cwiseAbs().maxCoeff() > 0.0);
    if (has_real)
    {
      print_table(resistance_R_inv->real(), "rom-Rinv-re.csv");
    }
    if (has_imag)
    {
      print_table(resistance_R_inv->imag(), "rom-Rinv-im.csv");
    }
  }

  // Print orth-R. Don't divide by diagonal to keep state normalization info.
  print_table(orth_R, "rom-orthogonalization-matrix-R.csv");
}

}  // namespace palace
