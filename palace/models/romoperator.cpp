// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "romoperator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <string_view>
#include <tuple>
#include <utility>
#include <mfem.hpp>
#include "fem/bilinearform.hpp"
#include "fem/integrator.hpp"
#include "linalg/operator.hpp"
#include "linalg/orthog.hpp"
#include "linalg/rap.hpp"
#include "models/floquetportoperator.hpp"
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
constexpr std::size_t WAVEPORT_AAA_ORDER_MAX = 12;
constexpr double WAVEPORT_SYNTHESIS_RANK_TOL_MAX = 1.0e-6;
// Cross-section EVP tolerances for sampling kₙ(ω)/M(ω) during synthesis, decoupled from the
// user's EigenTol/KSPTol so the fit is not resolved below the port-mode accuracy floor
// (fitting eigensolver noise makes the synthesized pencil MPI-partition dependent). The fit
// and rank tolerances below are floored at WAVEPORT_SYNTHESIS_EIG_TOL to match.
constexpr double WAVEPORT_SYNTHESIS_EIG_TOL = 1.0e-11;
constexpr double WAVEPORT_SYNTHESIS_KSP_TOL = 1.0e-12;
// Wave-port modal-correction synthesis subspace: band samples used to build the per-port
// subspace Q, and a rank cap (r(r+1)/2 scalar fits follow; cap guards a bad tol choice).
constexpr int WAVEPORT_SYNTHESIS_SUBSPACE_SAMPLES = 12;
constexpr int WAVEPORT_SYNTHESIS_SUBSPACE_RANK_MAX = 8;
// Rank floor for signed auxiliary residue factorizations. The previous 3e-3 cutoff dropped
// physical directions in modal-correction and boundary-mass residues, changing the realized
// Schur complement by percent-level amounts. A 1e-6 floor retains those directions while
// avoiding rank decisions in the numerical tail when AdaptiveTol is much tighter.
constexpr double WAVEPORT_SYNTHESIS_AUX_RANK_TOL = 1.0e-6;
// Significance floor for a modal-correction coupling part relative to its pair norm ‖S_pq‖.
// The modal subspace is real up to an arbitrary per-partition global phase, so Im(Q_p Q_qᵀ)
// is machine-epsilon phase noise (~1e-10·‖S_pq‖). Fitting it spends aux states on noise
// whose SVD rank flips across partitions, so skip any part below this floor as numerically
// zero.
constexpr double WAVEPORT_SYNTHESIS_MODAL_PART_TOL = 1.0e-6;
// Drop synthesized eigenvalues whose eigenvector energy in the basis rows falls below this
// fraction; aux states inject roots at their pole frequencies that live in the aux rows.
constexpr double SYNTHESIS_EIG_BASIS_FRAC_MIN = 1.0e-2;
// Synthesized-eigenvalue filter: the augmented realization's aux states (zero-capacitance
// rows, cond(C) = ∞) produce spurious near-critically-damped roots at Q ≲ 0.5 that carry
// no physical content.
constexpr double SYNTHESIS_EIG_Q_MIN = 0.5;

// Index of `target` in `labels`, or -1 when absent. Used to address rows of the
// synthesized matrices by their node label.
inline long LabelIndex(const std::vector<std::string> &labels, const std::string &target)
{
  auto it = std::find(labels.begin(), labels.end(), target);
  return (it == labels.end()) ? -1 : static_cast<long>(std::distance(labels.begin(), it));
}

// Shared sampling grids for the synthesis dispersion fits. Fit grids use
// Chebyshev–Gauss–Lobatto nodes (endpoints included, well-conditioned LSQ); dense
// validation grids use Chebyshev–Gauss (interior) nodes, which interlace the CGL fit
// nodes so validation points never coincide with fit points.
inline std::vector<double> SampleChebyshevLobatto(double w_lo, double w_hi, int n)
{
  const double w_mid = 0.5 * (w_lo + w_hi), w_half = 0.5 * (w_hi - w_lo);
  std::vector<double> ws(n);
  if (n == 1)
  {
    ws[0] = w_mid;
    return ws;
  }
  for (int i = 0; i < n; i++)
  {
    ws[i] = w_mid - w_half * std::cos(M_PI * i / (n - 1));
  }
  return ws;
}

inline std::vector<double> SampleChebyshevGauss(double w_lo, double w_hi, int n)
{
  const double w_mid = 0.5 * (w_lo + w_hi), w_half = 0.5 * (w_hi - w_lo);
  std::vector<double> ws(n);
  for (int j = 0; j < n; j++)
  {
    ws[j] = w_mid - w_half * std::cos(M_PI * (2 * j + 1) / (2.0 * n));
  }
  return ws;
}

// Complex LSQ quadratic fit y(ω) ≈ c0 + c1·ω + c2·ω² on the given nodes (real inputs
// come through with zero imaginary parts and produce real coefficients).
inline Eigen::Vector3cd FitQuadratic(const std::vector<double> &omegas,
                                     const Eigen::VectorXcd &y)
{
  const auto n = static_cast<long>(omegas.size());
  Eigen::MatrixXcd vandermonde(n, 3);
  for (long i = 0; i < n; i++)
  {
    const double w = omegas[static_cast<std::size_t>(i)];
    vandermonde(i, 0) = 1.0;
    vandermonde(i, 1) = w;
    vandermonde(i, 2) = w * w;
  }
  return vandermonde.colPivHouseholderQr().solve(y);
}

// AAA rational fit of a residual sampled on the fit grid, with the synthesis tolerance
// converted from "relative to the full function" (scale = max |truth| on the grid) to
// AAA's internal "relative to max |residual|" convention.
inline utils::AAAPoleResidue FitResidualPoles(const std::vector<double> &omegas,
                                              const Eigen::VectorXcd &residual, double tol,
                                              double truth_scale, std::size_t order_max)
{
  const auto n = static_cast<long>(omegas.size());
  Eigen::VectorXcd z(n);
  for (long i = 0; i < n; i++)
  {
    z(i) = omegas[static_cast<std::size_t>(i)];
  }
  const double aaa_tol_rel =
      (truth_scale > 0.0)
          ? tol * truth_scale / std::max(residual.cwiseAbs().maxCoeff(), 1.0e-300)
          : tol;
  auto aaa = utils::RunAAA(z, residual, aaa_tol_rel, std::max<std::size_t>(order_max, 1));
  return utils::AAAToPoleResidue(aaa);
}

// Unpack an AAAPoleResidue into the std::vector storage used by WavePortAuxBlock.
inline std::pair<std::vector<std::complex<double>>, std::vector<std::complex<double>>>
ToPoleResidueVectors(const utils::AAAPoleResidue &pr)
{
  std::vector<std::complex<double>> poles(pr.poles.data(),
                                          pr.poles.data() + pr.poles.size());
  std::vector<std::complex<double>> residues(pr.residues.data(),
                                             pr.residues.data() + pr.residues.size());
  return {std::move(poles), std::move(residues)};
}

// Fold the centered polynomial quotient from an AAA conversion into the quadratic
// synthesis pencil. Higher-degree quotients cannot be represented by K + iωC − ω²M and
// must not be silently discarded.
inline void FoldAAAPolynomialIntoQuadratic(const utils::AAAPoleResidue &pr,
                                           std::complex<double> &a0,
                                           std::complex<double> &a1,
                                           std::complex<double> &a2)
{
  MFEM_VERIFY(pr.polynomial.size() <= 3,
              "AAA polynomial quotient has degree "
                  << pr.polynomial.size() - 1
                  << ", which cannot be represented by a quadratic synthesis pencil!");
  MFEM_VERIFY(pr.polynomial_scale > 0.0,
              "AAA polynomial quotient must have a positive affine scale!");

  const std::complex<double> q0 = pr.polynomial.size() > 0 ? pr.polynomial(0) : 0.0;
  const std::complex<double> q1 = pr.polynomial.size() > 1 ? pr.polynomial(1) : 0.0;
  const std::complex<double> q2 = pr.polynomial.size() > 2 ? pr.polynomial(2) : 0.0;
  const double inverse_scale = 1.0 / pr.polynomial_scale;

  // AAA stores q in x = (ω-ωc)/s. Equivalently, x = ω/s-ω̂c with the nondimensional
  // center ω̂c = ωc/s. Compose q0+q1*x+q2*x² analytically into powers of the physical
  // omega used by the existing synthesis pencil; do not first expand it into an
  // ill-conditioned temporary vector of physical-omega monomial coefficients.
  const std::complex<double> nondimensional_center = pr.polynomial_center * inverse_scale;
  a0 +=
      q0 - q1 * nondimensional_center + q2 * nondimensional_center * nondimensional_center;
  a1 += (q1 - 2.0 * q2 * nondimensional_center) * inverse_scale;
  a2 += q2 * inverse_scale * inverse_scale;
}

// Evaluate the complex synthesis model α₀ + α₁ω + α₂ω² + Σₖ rₖ/(ω − pₖ).
inline std::complex<double> EvaluatePolyPlusPoles(std::complex<double> a0,
                                                  std::complex<double> a1,
                                                  std::complex<double> a2,
                                                  const utils::AAAPoleResidue &pr, double w)
{
  std::complex<double> val = a0 + a1 * w + a2 * w * w;
  for (long k = 0; k < pr.poles.size(); k++)
  {
    val += pr.residues(k) / (std::complex<double>(w, 0.0) - pr.poles(k));
  }
  return val;
}

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
                               ComplexVector &r, int n0, bool symmetric)
{
  // Update Ar = Vᴴ A V for the new basis dimension n0 -> n. V is real. Ar is replicated
  // across all processes as a sequential n x n matrix.
  const auto n = Ar.rows();
  MFEM_VERIFY(n0 < n, "Invalid dimensions in PROM matrix projection!");

  // Compute the right block: columns [n0, n), all rows [0, n).
  for (int j = n0; j < n; j++)
  {
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

  if (symmetric)
  {
    for (int j = 0; j < n0; j++)
    {
      for (int i = n0; i < n; i++)
      {
        Ar(i, j) = Ar(j, i);
      }
    }
    return;
  }

  // Compute the lower-left block directly: rows [n0, n), columns [0, n0).
  // No symmetry assumption — works for Hermitian, anti-symmetric, or general operators.
  for (int j = 0; j < n0; j++)
  {
    if (A.Real())
    {
      A.Real()->Mult(V[j], r.Real());
    }
    if (A.Imag())
    {
      A.Imag()->Mult(V[j], r.Imag());
    }
    for (int i = n0; i < n; i++)
    {
      Ar(i, j).real(A.Real() ? V[i] * r.Real() : 0.0);
      Ar(i, j).imag(A.Imag() ? V[i] * r.Imag() : 0.0);
    }
    Mpi::GlobalSum(n - n0, Ar.data() + j * n + n0, comm);
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
  Eigen::JacobiSVD<Eigen::MatrixXcd, Eigen::ComputeFullV> svd;
  svd.compute(R);
  const auto &sigma = svd.singularValues();
  auto m = S - 1;
  while (m > 0 && sigma[m] < ORTHOG_TOL * sigma[0])
  {
    Mpi::Warning("Minimal rational interpolation encountered a rank-deficient matrix: "
                 "σ[{:d}] = {:.3e} (σ[0] = {:.3e}). This can indicate that the "
                 "adaptive interpolation is near the accuracy limit of the HDM solves; "
                 "if adaptive convergence is poor, try tightening "
                 "config[\"Solver\"][\"Linear\"][\"Tol\"] or using a looser "
                 "config[\"Solver\"][\"Driven\"][\"AdaptiveTol\"].\n",
                 m, sigma[m], sigma[0]);
    m--;
  }
  q = svd.matrixV().col(m);
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
  z.push_back(omega);
}

std::vector<double> MinimalRationalInterpolation::FindMaxError(std::size_t N) const
{
  // Return an estimate for argmax_z ||u(z) - V y(z)|| as argmin_z |Q(z)| with Q(z) =
  // sum_i q_z / (z - z_i) (denominator of the barycentric interpolation of u).
  BlockTimer bt(Timer::CONSTRUCT_PROM);
  const auto S = dim_Q;
  MFEM_VERIFY(S >= 2, "Maximum error can only be found once two sample points have been "
                      "added to the PROM to define the parameter domain!");
  double start = *std::min_element(z.begin(), z.end());
  double end = *std::max_element(z.begin(), z.end());
  Eigen::Map<const Eigen::VectorXd> z_map(z.data(), S);

  // Sample Q on discrete points. The case of N>1 samples is not very useful below. It will
  // typically give us multiple sample points right next to each other in the same local
  // maximum, rather than N separate local maxima.

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
  // to A(ω) is then assembled at each ω as Σ_p k_{n,p}(ω)·M_{μ⁻¹,r,p}, where M_{r,p} is
  // projected onto the basis only when the basis grows. Cache masses for every active port
  // and for inactive ports included as unloaded synthesis nodes: the latter do not
  // contribute to A(ω), but their mass is required for node scaling and reference data.
  // This avoids HDM-scale work in the online phase. See
  // WavePortOperator::AddBoundaryMassBdrCoefficients and
  // SpaceOperator::GetWavePortBoundaryMassMatrix.
  for (const auto &[port_idx, port_data] : space_op.GetWavePortOp())
  {
    if (!port_data.active && !port_data.include_in_synthesis)
    {
      continue;
    }
    auto Mp = space_op.GetWavePortBoundaryMassMatrix<ComplexOperator>(port_idx,
                                                                      Operator::DIAG_ZERO);
    if (Mp)
    {
      Mwp_p.emplace(port_idx, std::move(Mp));
    }
  }
  // Detect whether GetExtraSystemMatrix has any non-wave-port contributions (e.g.
  // second-order farfield, surface conductivity, or Floquet Robin terms).
  {
    auto A2_other_probe = space_op.GetExtraSystemMatrix<ComplexOperator>(
        1.0, Operator::DIAG_ZERO, /*include_wave_ports=*/false);
    // The probe stamps every non-wave-port frequency-dependent term, including the
    // Floquet Robin BC, so a non-null probe is the complete condition.
    has_other_A2 = (A2_other_probe != nullptr);
  }

  // Cache the ω-independent boundary masses for the other frequency-dependent BCs so they
  // can be folded into circuit synthesis (projected per basis growth in UpdatePROM, fit
  // in CalculateNormalizedPROMMatrices). On the imaginary slot, matching Mwp_p.
  //   - 2nd-order farfield ABC: M_ff with full term i·(0.5/ω)·M_ff.
  //   - Surface conductivity: one boundary mass per active attribute group, term
  //     (i·ω/Z_g(ω))·A_σ_g.
  M_ff_ = space_op.GetFarfieldBoundaryCurlCurlMatrix<ComplexOperator>(Operator::DIAG_ZERO,
                                                                      /*imag_slot=*/true);
  {
    const auto &surf_op = space_op.GetSurfaceConductivityOp();
    Asig_g_.resize(surf_op.Size());
    Asig_g_r.resize(surf_op.Size());
    for (std::size_t g = 0; g < surf_op.Size(); g++)
    {
      if (surf_op.IsActive(g))
      {
        Asig_g_[g] = space_op.GetSurfaceConductivityBoundaryMatrix<ComplexOperator>(
            static_cast<int>(g), Operator::DIAG_ZERO);
      }
    }
  }

  // Per-port Floquet Robin µ⁻¹ boundary mass (imaginary slot, matching the convention).
  // The online scalar is γ₀,p(ω) = sqrt(max(0, ω²µε − kF²)), evaluated in closed form.
  for (const auto &[port_idx, port] : space_op.GetFloquetPortOp())
  {
    auto Mp = space_op.GetFloquetRobinBoundaryMassMatrix<ComplexOperator>(
        port_idx, Operator::DIAG_ZERO);
    if (Mp)
    {
      M_floquet_p_.emplace(port_idx, std::move(Mp));
    }
  }

  // Per-boundary rational surface impedance mass (imaginary slot, matching the
  // convention). The online scalar is f(ω) = g(iω)/i with g(s) = s·D(s)/N(s) the Robin
  // coefficient, evaluated in closed form; the synthesis fit is poly+AAA.
  {
    const auto &surf_rz_op = space_op.GetRationalImpedanceOp();
    Arz_b_.resize(surf_rz_op.GetNumBoundaries());
    Arz_b_r.resize(surf_rz_op.GetNumBoundaries());
    for (int b = 0; b < surf_rz_op.GetNumBoundaries(); b++)
    {
      Arz_b_[b] = space_op.GetRationalImpedanceBoundaryMassMatrix<ComplexOperator>(
          b, Operator::DIAG_ZERO, /*imag_slot=*/true);
    }
  }

  // Capture sweep band and synthesis controls for later use in
  // CalculateNormalizedPROMMatrices. After config::Nondimensionalize, sample_f stores
  // 2π·f_nondim, i.e. ω in nondimensional units. Use AdaptiveTol for the scalar fit
  // tolerance. Keep the SVD rank cutoff conservative for loose sweeps, but tighten it
  // automatically when the requested adaptive tolerance is below 1e-6.
  const auto &sample_f = iodata.solver.driven.sample_f;
  if (!sample_f.empty())
  {
    sweep_omega_samples = sample_f;
    sweep_omega_min = *std::min_element(sample_f.begin(), sample_f.end());
    sweep_omega_max = *std::max_element(sample_f.begin(), sample_f.end());
  }
  // Floor the fit/rank tolerances at the synthesis EVP accuracy floor: resolving finer than
  // the port modes are solved just chases eigensolver noise (see
  // WAVEPORT_SYNTHESIS_EIG_TOL).
  waveport_synthesis_tol =
      std::max(iodata.solver.driven.adaptive_tol, WAVEPORT_SYNTHESIS_EIG_TOL);
  waveport_synthesis_order_max = WAVEPORT_AAA_ORDER_MAX;
  waveport_synthesis_rank_tol = std::clamp(
      waveport_synthesis_tol, WAVEPORT_SYNTHESIS_EIG_TOL, WAVEPORT_SYNTHESIS_RANK_TOL_MAX);

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

  // Compute HDM solution at the given frequency. The sparse A2 is stored as a member for
  // PROM projection. The full frequency-dependent operator (A2 + low-rank Floquet DtN) is
  // built locally for the HDM system matrix.
  A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO);
  has_A2 = (A2 != nullptr);
  auto A2_full = space_op.GetExtraSystemOperator(omega, Operator::DIAG_ZERO);
  auto A = space_op.GetSystemMatrix(std::complex<double>(1.0, 0.0), 1i * omega,
                                    std::complex<double>(-omega * omega, 0.0), K.get(),
                                    C.get(), M.get(), A2_full.get());
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
  // single reference frequency (typically the band center) and projected onto the
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
  ProjectMatInternal(comm, V, *K, Kr, r, dim_V_old, true);
  if (C)
  {
    Cr.conservativeResize(dim_V_new, dim_V_new);
    ProjectMatInternal(comm, V, *C, Cr, r, dim_V_old,
                       !space_op.GetMaterialOp().HasFloquetFrequencyScaling());
  }
  Mr.conservativeResize(dim_V_new, dim_V_new);
  ProjectMatInternal(comm, V, *M, Mr, r, dim_V_old, true);
  // Per-port wave-port masses. M_{r,p} is initialized lazily so the map only contains
  // entries when the per-port HDM operator was non-null on this rank.
  for (auto &[port_idx, Mp_hdm] : Mwp_p)
  {
    auto &Mp_r = Mwp_p_r[port_idx];
    Mp_r.conservativeResize(dim_V_new, dim_V_new);
    ProjectMatInternal(comm, V, *Mp_hdm, Mp_r, r, dim_V_old, true);
  }
  // Other frequency-dependent BC boundary masses folded into circuit synthesis: the
  // 2nd-order farfield ABC (M_ff_) and each surface-conductivity group (Asig_g_). Projected
  // like the wave-port masses so the synthesis path treats them uniformly.
  if (M_ff_)
  {
    M_ff_r.conservativeResize(dim_V_new, dim_V_new);
    ProjectMatInternal(comm, V, *M_ff_, M_ff_r, r, dim_V_old, true);
  }
  for (std::size_t g = 0; g < Asig_g_.size(); g++)
  {
    if (Asig_g_[g])
    {
      Asig_g_r[g].conservativeResize(dim_V_new, dim_V_new);
      ProjectMatInternal(comm, V, *Asig_g_[g], Asig_g_r[g], r, dim_V_old, true);
    }
  }
  for (std::size_t b = 0; b < Arz_b_.size(); b++)
  {
    if (Arz_b_[b])
    {
      Arz_b_r[b].conservativeResize(dim_V_new, dim_V_new);
      ProjectMatInternal(comm, V, *Arz_b_[b], Arz_b_r[b], r, dim_V_old, true);
    }
  }
  // Per-port Floquet Robin boundary mass projection (same pattern as wave-port masses).
  for (auto &[port_idx, Mp_r] : M_floquet_p_r)
  {
    Mp_r.conservativeResize(dim_V_new, dim_V_new);
    ProjectMatInternal(comm, V, *M_floquet_p_.at(port_idx), Mp_r, r, dim_V_old, true);
  }
  // Initialize map entries for new ports that haven't been projected yet (first call).
  for (const auto &[port_idx, Mp_hdm] : M_floquet_p_)
  {
    if (M_floquet_p_r.find(port_idx) == M_floquet_p_r.end())
    {
      auto &Mp_r = M_floquet_p_r[port_idx];
      Mp_r.resize(dim_V_new, dim_V_new);
      ProjectMatInternal(comm, V, *Mp_hdm, Mp_r, r, 0, true);
    }
  }
  if (RHS1.Size())
  {
    RHS1r.conservativeResize(dim_V_new);
    ProjectVecInternal(comm, V, RHS1, RHS1r, dim_V_old);
  }

  // Update reduced Floquet port projection vectors. For F = Σ g_k v_k v_k^H,
  // the PROM contribution is Σ g_k (V^T v_k) (v_k^H V) = Σ g_k (V^T v_k) conj(V^T v_k)^T.
  // Since V is real: V^T v_k is stored as vk_V, conj(V^T v_k) as Vh_cvk.
  if (dim_V_old == 0)
  {
    // First time: enumerate all Floquet modes with for_dtn flag.
    floquet_reduced.clear();
    for (const auto &[port_idx, port] : space_op.GetFloquetPortOp())
    {
      for (const auto &order : port.GetOrders())
      {
        if (!HasFlag(order.use, FloquetModeUse::Dtn))
        {
          continue;
        }
        for (bool is_te : {true, false})
        {
          ReducedFloquetMode rm;
          rm.port_idx = port_idx;
          rm.order = &order;
          rm.is_te = is_te;
          rm.vk_V.resize(dim_V_new);
          rm.Vh_cvk.resize(dim_V_new);
          floquet_reduced.push_back(std::move(rm));
        }
      }
    }
  }
  for (auto &rm : floquet_reduced)
  {
    // Compute V^T v_k for new basis vectors [dim_V_old, dim_V_new).
    rm.vk_V.conservativeResize(dim_V_new);
    rm.Vh_cvk.conservativeResize(dim_V_new);
    for (int i = dim_V_old; i < dim_V_new; i++)
    {
      // V[i] is real. v_k is complex. V[i]^T v_k = (V[i] · v_k_real) + j(V[i] · v_k_imag).
      double dr = V[i] * rm.order->v[rm.is_te ? 0 : 1].Real();
      double di = V[i] * rm.order->v[rm.is_te ? 0 : 1].Imag();
      Mpi::GlobalSum(1, &dr, comm);
      Mpi::GlobalSum(1, &di, comm);
      std::complex<double> vt_vi(dr, di);
      rm.vk_V(i) = vt_vi;               // v_k^T V[i]
      rm.Vh_cvk(i) = std::conj(vt_vi);  // V[i]^H conj(v_k) = conj(V[i]^T v_k)
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

  // Refresh Floquet port state for this frequency once, up front: the factored Robin
  // term below reads gamma0 and the low-rank DtN correction reads gamma_sq (and, when
  // k_F scales with frequency, the recomputed mode vectors).
  if (!space_op.GetFloquetPortOp().Empty())
  {
    space_op.GetFloquetPortOp().Initialize(omega);
  }

  // Other ω-nonlinear A2 contributors (second-order farfield ABC and surface conductivity).
  // These are applied in factored form: their ω-independent boundary masses (M_ff_r,
  // Asig_g_r) were projected onto the basis once in UpdatePROM, exactly like the wave-port
  // masses, so the online cost is a per-ω scalar times an n×n matrix add — no per-ω HDM-
  // scale assembly or reprojection. This is algebraically identical to projecting the full
  // A2(ω) here (the scalar is uniform per boundary group, so it commutes with the
  // projection) and matches the HDM stamping to round-off.
  //
  // Robustness: the structural check below requires every factored operator we hold to be
  // sized to the current basis, but it cannot know whether the factored set is COMPLETE
  // (an ω-dependent non-wave-port BC in GetExtraSystemMatrix without a factored operator —
  // e.g. Floquet Robin terms — would be silently dropped). So on the first factored online
  // solve we additionally verify the factored Aᵣ contribution against the full HDM
  // projection; on any mismatch we latch other_A2_factored_ok = false and use the slow
  // fallback for the rest of the sweep.
  auto apply_factored_other_A2 = [&]()
  {
    // Factored 2nd-order farfield ABC: A2_ff(ω) = i·(0.5/ω)·M_ff. M_ff_r carries the
    // boundary mass on the imaginary slot (the i), so scaling by the real scalar 0.5/ω
    // reproduces the full contribution.
    if (M_ff_ && M_ff_r.rows() == static_cast<long>(V.size()))
    {
      Ar += std::complex<double>(0.5 / omega, 0.0) * M_ff_r;
    }
    // Factored surface conductivity, per active group: A2_σ,g(ω) = (i·ω/Z_g(ω))·A_σ,g =
    // EvaluateScalar(g,ω)·A_σ,g. Asig_g_r[g] carries A_σ,g on the imaginary slot, so the
    // scalar here is EvaluateScalar/i to avoid double-counting the i (matching the
    // synthesis convention in CalculateNormalizedPROMMatrices). EvaluateScalar is closed
    // form (skin depth + optional finite-thickness correction) — a few transcendental ops
    // per group, negligible versus the reduced solve. No AAA needed online.
    const auto &surf_op = space_op.GetSurfaceConductivityOp();
    for (std::size_t g = 0; g < Asig_g_.size(); g++)
    {
      if (Asig_g_[g] && Asig_g_r[g].rows() == static_cast<long>(V.size()))
      {
        const std::complex<double> s =
            surf_op.EvaluateScalar(g, std::complex<double>(omega, 0.0)) /
            std::complex<double>(0.0, 1.0);
        Ar += s * Asig_g_r[g];
      }
    }
    // Factored rational surface impedance, per boundary: A2_rz,b(ω) = g(iω)·M_b =
    // i·(g(iω)/i)·M_b. Arz_b_r[b] carries M_b on the imaginary slot, so the scalar here
    // is EvalRobinCoefficient/i (matching the synthesis convention). Closed form (two
    // Horner evaluations per boundary), no AAA needed online.
    const auto &surf_rz_op = space_op.GetRationalImpedanceOp();
    for (std::size_t b = 0; b < Arz_b_.size(); b++)
    {
      if (Arz_b_[b] && Arz_b_r[b].rows() == static_cast<long>(V.size()))
      {
        const std::complex<double> s =
            surf_rz_op.EvalRobinCoefficient(static_cast<int>(b),
                                            std::complex<double>(0.0, omega)) /
            std::complex<double>(0.0, 1.0);
        Ar += s * Arz_b_r[b];
      }
    }
    // Factored Floquet port Robin BC: A2_floquet,p(ω) = i·γ₀,p(ω)·M_floquet_p.
    // M_floquet_p_r carries the µ⁻¹ boundary mass on the imaginary slot (the i), so
    // the scalar multiplier is γ₀ (real, refreshed by the Initialize(omega) at the top
    // of SolvePROM).
    for (const auto &[port_idx, Mp_r] : M_floquet_p_r)
    {
      if (Mp_r.rows() == static_cast<long>(V.size()))
      {
        const double gamma0 = space_op.GetFloquetPortOp().GetPort(port_idx).GetGamma0();
        Ar += std::complex<double>(gamma0, 0.0) * Mp_r;
      }
    }
  };

  // Structural precondition for the factored path: every factored operator we hold must be
  // sized to the current basis, and we must hold at least one (else has_other_A2 came from
  // a BC we don't factor).
  bool other_A2_factored = false;
  if (has_other_A2 && other_A2_factored_ok)
  {
    const long n = static_cast<long>(V.size());
    bool any_factored = false, all_present = true;
    if (M_ff_)
    {
      (M_ff_r.rows() == n) ? (any_factored = true) : (all_present = false);
    }
    for (std::size_t g = 0; g < Asig_g_.size(); g++)
    {
      if (Asig_g_[g])
      {
        (Asig_g_r[g].rows() == n) ? (any_factored = true) : (all_present = false);
      }
    }
    for (std::size_t b = 0; b < Arz_b_.size(); b++)
    {
      if (Arz_b_[b])
      {
        (Arz_b_r[b].rows() == n) ? (any_factored = true) : (all_present = false);
      }
    }
    for (const auto &[port_idx, Mp_r] : M_floquet_p_r)
    {
      (Mp_r.rows() == n) ? (any_factored = true) : (all_present = false);
    }
    other_A2_factored = any_factored && all_present;
  }

  Ar.setZero();
  if (has_other_A2 && other_A2_factored)
  {
    apply_factored_other_A2();
    if (!other_A2_self_checked)
    {
      // One-time correctness self-check: compare the factored contribution to the full HDM
      // projection of A2_other(ω) at this frequency. Cheap (runs once per RomOperator).
      Eigen::MatrixXcd Ar_factored = Ar;
      Eigen::MatrixXcd Ar_hdm = Eigen::MatrixXcd::Zero(V.size(), V.size());
      A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO,
                                                          /*include_wave_ports=*/false);
      if (A2)
      {
        ProjectMatInternal(space_op.GetComm(), V, *A2, Ar_hdm, r, 0, true);
      }
      const double err = (Ar_factored - Ar_hdm).cwiseAbs().maxCoeff();
      const double ref =
          std::max(Ar_hdm.cwiseAbs().maxCoeff(), Ar_factored.cwiseAbs().maxCoeff());
      if (ref == 0.0)
      {
        // Both contributions vanish at this frequency (degenerate compare, e.g. all
        // frequency-dependent terms zero at the first sweep point): nothing learned,
        // retry the check at the next online frequency.
      }
      else
      {
        other_A2_self_checked = true;
        if (err / ref > 1.0e-9)
        {
          other_A2_factored_ok = false;
          Ar = Ar_hdm;  // Use the trusted HDM projection for this solve.
          Mpi::Warning(
              "Factored online A2 (farfield ABC, surface conductivity, rational "
              "impedance, Floquet Robin) disagrees with the full operator "
              "(rel. err {:.3e})!\n"
              "Reverting to the per-frequency assembled A2 for the remaining sweep. "
              "This indicates an ω-dependent boundary condition not covered by the "
              "factored path.\n",
              err / ref);
        }
      }
    }
  }
  else if (has_other_A2)
  {
    // Slow fallback: reassemble and reproject the full non-wave-port A2(ω) per ω.
    A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO,
                                                        /*include_wave_ports=*/false);
    if (A2)
    {
      ProjectMatInternal(space_op.GetComm(), V, *A2, Ar, r, 0, true);
    }
  }
  Ar += Kr;
  if (C)
  {
    Ar += (1i * omega) * Cr;
  }
  Ar += (-omega * omega) * Mr;
  // Wave-port contribution: A_wp(ω) = i·Σ_p k_{n,p}(ω)·M_{μ⁻¹,p}. GetWavePortKn re-solves
  // the per-port cross-section EVP at this ω and refreshes the modal post-processing state
  // used by MeasureWavePorts for S-parameters and power.
  {
    BlockTimer bt(Timer::WAVE_PORT);
    for (const auto &[port_idx, Mp_r] : Mwp_p_r)
    {
      const auto &port_data = space_op.GetWavePortOp().GetPort(port_idx);
      if (!port_data.active)
      {
        continue;
      }
      const double kn = space_op.GetWavePortOp().GetWavePortKn(port_idx, omega);
      Ar += std::complex<double>(kn, 0.0) * Mp_r;
    }

    // Add the low-rank wave-port modal correction Wᵣ(ω) = Σ_k g_k(ω) (Vᵀs_k)(Vᵀs_k)ᵀ, the
    // Galerkin projection of the complex-symmetric W = Σ_k g_k s_k s_kᵀ that the uniform
    // path applies via GetExtraSystemOperator (V is real, so Vᴴ = Vᵀ). Reassembled per ω
    // since the modal fields and reactions change with frequency. Without it Aᵣ carries
    // only i·kₙ·M while the reduced RHS carries the full modal n×H, breaking unitarity for
    // reactive/TM modes; for TEM modes s_full = s_scalar and Wᵣ ≡ 0.
    auto wp_terms = space_op.GetModalCorrectionTerms(omega);
    Eigen::VectorXcd sV(V.size());
    for (auto &term : wp_terms)
    {
      ProjectVecInternal(space_op.GetComm(), V, *term.s, sV, 0);
      Ar.noalias() += term.g * sV * sV.transpose();
    }
  }

  // Add low-rank Floquet port DtN correction: Fᵣ = Σ g_k(ω) (V^T v_k) conj(V^T v_k)^T
  // (per-frequency state refreshed by the Initialize(omega) at the top of SolvePROM).

  // When k_F scales with frequency, the mode vectors v_k change (polarization rotation).
  // Reproject onto the PROM basis V for the current frequency.
  if (space_op.GetMaterialOp().HasFloquetFrequencyScaling())
  {
    MPI_Comm comm = space_op.GetComm();
    auto dim_V = static_cast<int>(V.size());
    for (auto &rm : floquet_reduced)
    {
      for (int i = 0; i < dim_V; i++)
      {
        double dr = V[i] * rm.order->v[rm.is_te ? 0 : 1].Real();
        double di = V[i] * rm.order->v[rm.is_te ? 0 : 1].Imag();
        Mpi::GlobalSum(1, &dr, comm);
        Mpi::GlobalSum(1, &di, comm);
        std::complex<double> vt_vi(dr, di);
        rm.vk_V(i) = vt_vi;
        rm.Vh_cvk(i) = std::conj(vt_vi);
      }
    }
  }

  for (const auto &rm : floquet_reduced)
  {
    const auto &port = space_op.GetFloquetPortOp().GetPort(rm.port_idx);
    auto g = port.ComputeDtNCorrectionCoeff(*rm.order, rm.is_te);
    if (g != 0.0)
    {
      Ar.noalias() += g * rm.vk_V * rm.Vh_cvk.transpose();
    }
  }

  if (has_RHS2)
  {
    // NOTE: this per-ω HDM-scale assembly + projection of RHS2(ω) is intentional, not an
    // oversight. RHS2 depends on the wave-port modal fields, whose per-ω refresh (via the
    // cross-section EVP triggered above) is required anyway for correct S-parameter
    // post-processing; a cached-sample interpolation of the projected RHS2 was tried and
    // reverted for exactly this reason.
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
  // QR solve, for maximal stability. The small system is cheap to compute but can be
  // numerically poorly conditioned to due the splitting of HDM solutions into Re and Im
  // into separate columns.
  RHSr = Ar.fullPivHouseholderQr().solve(RHSr);
  ProlongatePROMSolution(V.size(), V, RHSr, u);
}

RomOperator::WavePortDispersionFit
RomOperator::FitWavePortDispersion(int port_idx, const Eigen::MatrixXcd &Mp_r) const
{
  // Sample kₙ,p on the sweep band, fit a quadratic, optionally augment with AAA on the
  // residual. Returns a packed result describing the chosen regime and the fit data.
  WavePortDispersionFit fit;
  fit.port_idx = port_idx;

  // kₙ(ω) is sampled by re-solving the per-port cross-section EVP at each ω (the
  // dominant cost here), so the sample counts are deliberately small. The model is at most
  // an order-2 polynomial plus a few AAA poles, so a modest fit grid suffices, and a
  // smooth analytic kₙ(ω) needs only a handful of out-of-sample points to bound the fit
  // residual. n_fit is AAA's candidate pool (it greedily selects support points FROM these
  // samples — it does not sample new ω), so it must comfortably exceed the pole budget.
  //
  // Both grids are Chebyshev–Gauss–Lobatto nodes on [w_lo, w_hi].
  const int n_fit = std::max(12, 2 * static_cast<int>(waveport_synthesis_order_max) + 4);
  const int n_dense = 2 * n_fit;
  auto fit_omegas = SampleChebyshevLobatto(sweep_omega_min, sweep_omega_max, n_fit);
  auto dense_omegas = SampleChebyshevGauss(sweep_omega_min, sweep_omega_max, n_dense);

  // Sample kₙ,p on the fit grid. Each call triggers (or hits the cache of)
  // WavePortData::Initialize. Cheap because the cross-section EVP is small.
  Eigen::VectorXcd y_fit(n_fit);
  for (int i = 0; i < n_fit; i++)
  {
    y_fit(i) = space_op.GetWavePortOp().GetWavePortKn(port_idx, fit_omegas[i]);
  }
  // LSQ polynomial fit at order 2 in ω. Higher orders cannot be absorbed into the
  // K + iωC − ω²M synthesis structure (cf. design notes). kₙ is real, so the complex
  // LSQ returns real coefficients (imaginary parts exactly zero).
  Eigen::Vector3cd coeffs = FitQuadratic(fit_omegas, y_fit);
  fit.alpha0 = coeffs(0).real();
  fit.alpha1 = coeffs(1).real();
  fit.alpha2 = coeffs(2).real();

  // Residual on the dense grid: relative error on kₙ. The 1/(ωμ) factor in the
  // wave-port admittance Y_p(ω) = kₙ,p(ω)/(iωμ) cancels in the ratio, so kₙ is the
  // right proxy for the user-visible synthesis accuracy (no HDM data required).
  double max_rel = 0.0, max_abs_truth = 0.0;
  for (int i = 0; i < n_dense; i++)
  {
    const double w = dense_omegas[i];
    const double truth = space_op.GetWavePortOp().GetWavePortKn(port_idx, w);
    // fit.aux is still empty here, so EvaluateWavePortKnFit returns the bare polynomial.
    const double poly = EvaluateWavePortKnFit(fit, w);
    max_abs_truth = std::max(max_abs_truth, std::abs(truth));
    max_rel = std::max(max_rel, std::abs(poly - truth));
  }
  fit.rel_err_polynomial = (max_abs_truth > 0.0) ? max_rel / max_abs_truth : 0.0;
  fit.rel_err_augmented = fit.rel_err_polynomial;
  const bool meets_tol = (fit.rel_err_polynomial <= waveport_synthesis_tol);

  fit.regime = meets_tol ? WavePortRegime::Polynomial : WavePortRegime::Augmented;

  if (fit.regime == WavePortRegime::Polynomial)
  {
    if (meets_tol)
    {
      Mpi::Print(
          " Wave port {:d}: polynomial synthesis residual {:.3e} (tol {:.3e}, α₀={:.3e}, "
          "α₁={:.3e}, α₂={:.3e})\n",
          port_idx, fit.rel_err_polynomial, waveport_synthesis_tol, fit.alpha0, fit.alpha1,
          fit.alpha2);
    }
    return fit;
  }

  // Augmented regime: AAA rational fit on the polynomial residual δkₙ(ω). Cap the
  // pole count at the internal guard waveport_synthesis_order_max.
  Eigen::VectorXcd F_aaa(n_fit);
  for (int i = 0; i < n_fit; i++)
  {
    const double w = fit_omegas[i];
    F_aaa(i) = y_fit(i) - (fit.alpha0 + fit.alpha1 * w + fit.alpha2 * w * w);
  }
  auto pr = FitResidualPoles(fit_omegas, F_aaa, waveport_synthesis_tol, max_abs_truth,
                             waveport_synthesis_order_max);

  // Fold any polynomial quotient from the barycentric-to-pole conversion into the base
  // quadratic. Real wave-port samples produce a real quotient up to roundoff.
  std::complex<double> alpha0 = fit.alpha0, alpha1 = fit.alpha1, alpha2 = fit.alpha2;
  FoldAAAPolynomialIntoQuadratic(pr, alpha0, alpha1, alpha2);
  fit.alpha0 = alpha0.real();
  fit.alpha1 = alpha1.real();
  fit.alpha2 = alpha2.real();
  const double aaa_q0 = pr.polynomial.size() > 0 ? pr.polynomial(0).real() : 0.0;

  auto [poles, residues] = ToPoleResidueVectors(pr);

  WavePortAuxBlock blk;
  blk.port_idx = port_idx;
  AddAuxBlockDirections(blk, Mp_r, waveport_synthesis_rank_tol);
  blk.poles = std::move(poles);
  blk.residues = std::move(residues);
  // Attach the aux block now so EvaluateWavePortKnFit picks up the AAA poles for the
  // augmented-model residual check below.
  fit.aux = std::move(blk);

  // Re-evaluate residual on the dense grid using the augmented model, via the shared
  // evaluator (α-polynomial + Σ Re(rₖ/(ω−pₖ))).
  max_rel = 0.0;
  for (int i = 0; i < n_dense; i++)
  {
    const double w = dense_omegas[i];
    const double truth = space_op.GetWavePortOp().GetWavePortKn(port_idx, w);
    max_rel = std::max(max_rel, std::abs(EvaluateWavePortKnFit(fit, w) - truth));
  }
  fit.rel_err_augmented = (max_abs_truth > 0.0) ? max_rel / max_abs_truth : 0.0;

  const std::size_t n_poles = fit.aux->poles.size();
  const std::size_t rank_used = fit.aux->weights.size();
  const std::size_t aux_per_port = n_poles * rank_used;
  Mpi::Print(" Wave port {:d}: augmented synthesis residual {:.3e} → {:.3e} "
             "(tol {:.3e}, {:d} pole{} × rank-{:d} mass = +{:d} aux state{}, "
             "α₀={:.3e}, α₁={:.3e}, α₂={:.3e}, q₀={:.3e})\n",
             port_idx, fit.rel_err_polynomial, fit.rel_err_augmented,
             waveport_synthesis_tol, n_poles, n_poles == 1 ? "" : "s", rank_used,
             aux_per_port, aux_per_port == 1 ? "" : "s", fit.alpha0, fit.alpha1, fit.alpha2,
             aaa_q0);
  if (fit.rel_err_augmented > waveport_synthesis_tol)
  {
    Mpi::Warning("Wave port {:d}: augmented synthesis residual {:.3e} exceeds "
                 "AdaptiveTol={:.3e}; synthesized L/R/C matrices may be less accurate "
                 "for this port.\n",
                 port_idx, fit.rel_err_augmented, waveport_synthesis_tol);
  }
  return fit;
}

double RomOperator::EvaluateWavePortKnFit(const WavePortDispersionFit &fit, double omega)
{
  // Real wave-port dispersion model: kₙ(ω) ≈ α₀ + α₁ω + α₂ω² + Σₖ Re(rₖ/(ω−pₖ)). The pole
  // sum is present only in the Augmented regime (fit.aux holds the AAA poles/residues); in
  // the Polynomial regime fit.aux is empty and this reduces to the quadratic. This is the
  // analytic continuation evaluated on the real axis — the same expression the synthesis
  // pencil realizes via Kr/Cr/Mr (α-part) plus aux states (pole part).
  double kn = fit.alpha0 + fit.alpha1 * omega + fit.alpha2 * omega * omega;
  if (fit.aux)
  {
    const std::complex<double> wc(omega, 0.0);
    for (std::size_t k = 0; k < fit.aux->poles.size(); k++)
    {
      kn += (fit.aux->residues[k] / (wc - fit.aux->poles[k])).real();
    }
  }
  return kn;
}

void RomOperator::ApplyPolynomialFitCorrections(const WavePortDispersionFit &fit,
                                                const Eigen::MatrixXcd &Mp_r,
                                                Eigen::MatrixXcd &Kr_corr,
                                                Eigen::MatrixXcd &Cr_corr,
                                                Eigen::MatrixXcd &Mr_corr)
{
  // Mp_r is purely imaginary (= i·M_proj). Multiplying by α₀ folds α₀·M_proj into the
  // imaginary part of Kr (which becomes Im(L⁻¹)); +α₁·M_proj appears in the real part
  // of Cr (Re(R⁻¹), since iω·(α₁·M_proj) supplies the i); and -α₂·M_proj adds to the
  // imaginary part of Mr (Im(C)). Signs are set so the sum recovers
  // i·(α₀+α₁ω+α₂ω²)·M_proj·v in Aᵣ(ω).
  Kr_corr += std::complex<double>(fit.alpha0, 0.0) * Mp_r;
  Cr_corr += std::complex<double>(0.0, -fit.alpha1) * Mp_r;
  Mr_corr += std::complex<double>(-fit.alpha2, 0.0) * Mp_r;
}

void RomOperator::ApplyComplexPolynomialFitCorrections(
    std::complex<double> alpha0, std::complex<double> alpha1, std::complex<double> alpha2,
    const Eigen::MatrixXcd &Mp_r, Eigen::MatrixXcd &Kr_corr, Eigen::MatrixXcd &Cr_corr,
    Eigen::MatrixXcd &Mr_corr)
{
  // Complex-α generalization of ApplyPolynomialFitCorrections. Mp_r = i·M_proj. The target
  // contribution to Aᵣ(ω) is i·(α₀+α₁ω+α₂ω²)·M_proj·v = (α₀+α₁ω+α₂ω²)·Mp_r·v, matched into
  // K + iωC − ω²M by: α₀ → Kr (constant), α₁ → Cr via iω·(α₁/i)·Mp_r so Cr_corr gets
  // -i·α₁·Mp_r, and α₂ → Mr via −ω²·(−α₂)·Mp_r so Mr_corr gets −α₂·Mp_r. (Identical sign
  // structure to the real-α version, now carrying the full complex coefficients.)
  Kr_corr += alpha0 * Mp_r;
  Cr_corr += std::complex<double>(0.0, -1.0) * alpha1 * Mp_r;
  Mr_corr += (-alpha2) * Mp_r;
}

bool RomOperator::AddAuxBlockDirections(WavePortAuxBlock &blk, const Eigen::MatrixXcd &Mp_r,
                                        double rank_tol)
{
  // Mp_r = i M_proj by convention. M_proj is real symmetric, but modal-correction
  // matrices are generally indefinite (for example q_p q_qᵀ + q_q q_pᵀ). An SVD with
  // unsigned singular values would realize |M_proj| rather than M_proj in the auxiliary
  // Schur complement. Use a signed self-adjoint eigendecomposition instead:
  //                         M_proj = Σ_j λ_j u_j u_jᵀ.
  const double matrix_scale = std::max(Mp_r.cwiseAbs().maxCoeff(), 1.0e-300);
  const double real_rel = Mp_r.real().cwiseAbs().maxCoeff() / matrix_scale;
  MFEM_VERIFY(real_rel <= 1.0e-10,
              "Synthesis residue direction matrix must use the purely imaginary slot "
                  << "(relative real part " << real_rel << ")!");
  Eigen::MatrixXd M_proj = Mp_r.imag().eval();
  const double sym_scale = std::max(M_proj.cwiseAbs().maxCoeff(), 1.0e-300);
  const double skew_rel = (M_proj - M_proj.transpose()).cwiseAbs().maxCoeff() / sym_scale;
  MFEM_VERIFY(skew_rel <= 1.0e-10, "Synthesis residue direction matrix must be symmetric "
                                       << "(relative skew part " << skew_rel << ")!");
  M_proj = 0.5 * (M_proj + M_proj.transpose()).eval();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(M_proj);
  MFEM_VERIFY(eig.info() == Eigen::Success,
              "Failed symmetric eigendecomposition of synthesis residue matrix!");
  if (eig.eigenvalues().size() == 0)
  {
    return false;
  }
  const double weight_max = eig.eigenvalues().cwiseAbs().maxCoeff();
  // Do not let an extremely tight fit tolerance resolve the partition-dependent numerical
  // tail, but retain physical modal-correction directions below the old 3e-3 floor.
  const double eff_tol = std::max(rank_tol, WAVEPORT_SYNTHESIS_AUX_RANK_TOL);
  // SelfAdjointEigenSolver orders eigenvalues increasingly. Iterate from largest absolute
  // weight only for deterministic labels; the factorization itself is order independent.
  std::vector<long> order(static_cast<std::size_t>(eig.eigenvalues().size()));
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(
      order.begin(), order.end(), [&](long a, long b)
      { return std::abs(eig.eigenvalues()(a)) > std::abs(eig.eigenvalues()(b)); });
  for (long j : order)
  {
    const double weight = eig.eigenvalues()(j);
    if (weight_max > 0.0 && std::abs(weight) / weight_max > eff_tol)
    {
      blk.weights.push_back(weight);
      blk.u_dirs.push_back(eig.eigenvectors().col(j));
    }
  }
  return !blk.weights.empty();
}

std::optional<RomOperator::WavePortAuxBlock>
RomOperator::MakeAuxBlock(std::string label, const Eigen::MatrixXcd &Mp_r,
                          const std::vector<std::complex<double>> &poles,
                          const std::vector<std::complex<double>> &residues,
                          double rank_tol)
{
  // Build an aux block from a projected boundary mass Mp_r (purely imaginary = i·M_proj)
  // and a pole-residue list: eigendecompose the real symmetric M_proj into signed coupling
  // directions, then
  // attach the poles. The label prefixes the aux-state row names in the synthesized
  // matrices; requiring it here (rather than patching blk.label afterwards) prevents an
  // unlabeled non-wave-port block from colliding with the waveport_<idx> label space.
  MFEM_VERIFY(!label.empty(), "MakeAuxBlock requires a nonempty aux-state label!");
  if (poles.empty())
  {
    return std::nullopt;
  }
  WavePortAuxBlock blk;
  blk.port_idx = -1;
  blk.label = std::move(label);
  AddAuxBlockDirections(blk, Mp_r, rank_tol);
  if (blk.weights.empty())
  {
    return std::nullopt;
  }
  blk.poles = poles;
  blk.residues = residues;
  return blk;
}

RomOperator::WavePortDispersionFit RomOperator::FitScalarDispersion(
    const std::string &label, const Eigen::MatrixXcd &Mp_r,
    const std::function<std::complex<double>(std::complex<double>)> &f,
    bool allow_augment) const
{
  // Generalized scalar-dispersion fit for a non-wave-port frequency-dependent BC. Mirrors
  // FitWavePortDispersion but samples an arbitrary (generally complex) scalar f(ω) instead
  // of the real wave-port kₙ(ω), and carries complex polynomial coefficients. The aux block
  // (when augmenting) is built from Mp_r via MakeAuxBlock. NOTE: this routine only sets the
  // fit metadata (regime + residuals) and the aux block; the polynomial part is applied by
  // the caller via ApplyComplexPolynomialFitCorrections using the complex coefficients
  // stored in alpha0c/alpha1c/alpha2c below (returned through the fit's real α fields is
  // not possible since those are real, so the caller reads the complex coeffs from the
  // returned struct's *_c members — see WavePortDispersionFit).
  WavePortDispersionFit fit;
  fit.port_idx = -1;  // not a wave port

  constexpr int n_fit = 30;
  constexpr int n_dense = 200;
  // Chebyshev nodes match FitWavePortDispersion: better-conditioned LSQ than uniform
  // sampling, and validation nodes interlace the fit nodes.
  auto fit_omegas = SampleChebyshevLobatto(sweep_omega_min, sweep_omega_max, n_fit);
  auto dense_omegas = SampleChebyshevGauss(sweep_omega_min, sweep_omega_max, n_dense);

  // Complex LSQ polynomial fit of f(ω) ≈ c0 + c1 ω + c2 ω² at order 2.
  Eigen::VectorXcd y(n_fit);
  for (int i = 0; i < n_fit; i++)
  {
    y(i) = f(std::complex<double>(fit_omegas[i], 0.0));
  }
  Eigen::Vector3cd c = FitQuadratic(fit_omegas, y);
  fit.alpha0c = c(0);
  fit.alpha1c = c(1);
  fit.alpha2c = c(2);

  // Dense-grid residual (absolute, relative to max |f|).
  double max_rel = 0.0, max_abs_truth = 0.0;
  for (int i = 0; i < n_dense; i++)
  {
    const double w = dense_omegas[i];
    const std::complex<double> truth = f(std::complex<double>(w, 0.0));
    const std::complex<double> poly = c(0) + c(1) * w + c(2) * w * w;
    max_abs_truth = std::max(max_abs_truth, std::abs(truth));
    max_rel = std::max(max_rel, std::abs(poly - truth));
  }
  fit.rel_err_polynomial = (max_abs_truth > 0.0) ? max_rel / max_abs_truth : 0.0;
  fit.rel_err_augmented = fit.rel_err_polynomial;
  const bool meets_tol = (fit.rel_err_polynomial <= waveport_synthesis_tol);

  if (!allow_augment || meets_tol)
  {
    fit.regime = WavePortRegime::Polynomial;
    Mpi::Print(" {}: polynomial synthesis residual {:.3e} (tol {:.3e})\n", label,
               fit.rel_err_polynomial, waveport_synthesis_tol);
    return fit;
  }

  // Augmented regime: AAA on the complex polynomial residual.
  Eigen::VectorXcd F_aaa(n_fit);
  for (int i = 0; i < n_fit; i++)
  {
    const double w = fit_omegas[i];
    F_aaa(i) = y(i) - (c(0) + c(1) * w + c(2) * w * w);
  }
  auto pr = FitResidualPoles(fit_omegas, F_aaa, waveport_synthesis_tol, max_abs_truth,
                             waveport_synthesis_order_max);
  FoldAAAPolynomialIntoQuadratic(pr, fit.alpha0c, fit.alpha1c, fit.alpha2c);

  max_rel = 0.0;
  for (int i = 0; i < n_dense; i++)
  {
    const double w = dense_omegas[i];
    const std::complex<double> truth = f(std::complex<double>(w, 0.0));
    max_rel = std::max(max_rel, std::abs(EvaluatePolyPlusPoles(fit.alpha0c, fit.alpha1c,
                                                               fit.alpha2c, pr, w) -
                                         truth));
  }
  fit.rel_err_augmented = (max_abs_truth > 0.0) ? max_rel / max_abs_truth : 0.0;

  auto [poles, residues] = ToPoleResidueVectors(pr);
  fit.aux = MakeAuxBlock(label, Mp_r, poles, residues, waveport_synthesis_rank_tol);
  fit.regime = WavePortRegime::Augmented;
  Mpi::Print(
      " {}: augmented synthesis residual {:.3e} → {:.3e} (tol {:.3e}, {:d} pole{})\n",
      label, fit.rel_err_polynomial, fit.rel_err_augmented, waveport_synthesis_tol,
      pr.poles.size(), pr.poles.size() == 1 ? "" : "s");
  return fit;
}

RomOperator::WavePortDispersionFit
RomOperator::FitRationalImpedanceDispersion(const std::string &label,
                                            const Eigen::MatrixXcd &Mp_r, int idx) const
{
  // Exact split of the rational Robin coefficient g(s) = s·D(s)/N(s) = P(s) + R(s)/N(s)
  // (long division, precomputed by SurfaceRationalImpedanceOperator). The full term is
  // g(iω)·M_b = i·f(ω)·M_b with f(ω) = g(iω)/i, so the polynomial part contributes
  // f_P(ω) = P(iω)/i = P₀/i + P₁ω + i·P₂ω² + ... — exactly representable in the pencil
  // for deg(P) <= 2 (guaranteed for any passive N/D). The strictly proper remainder
  // f_R(ω) = (R/N)(iω)/i is itself rational in ω with the poles of g (the zeros of Zs
  // mapped to the ω plane), so AAA on f_R directly reproduces it to rounding.
  const auto &surf_rz_op = space_op.GetRationalImpedanceOp();
  const int deg_P = surf_rz_op.GetRobinQuotientDegree(idx);
  if (deg_P > 2)
  {
    // Non-passive input (|deg N - deg D| > 1); fall back to the generic recipe.
    auto f = [&surf_rz_op, idx](std::complex<double> omega) -> std::complex<double>
    {
      return surf_rz_op.EvalRobinCoefficient(idx, std::complex<double>(0.0, 1.0) * omega) /
             std::complex<double>(0.0, 1.0);
    };
    return FitScalarDispersion(label, Mp_r, f, /*allow_augment=*/true);
  }

  WavePortDispersionFit fit;
  fit.port_idx = -1;

  // Polynomial part, exact: with P(s) = Σ_k P_k s^k, f_P(ω) = Σ_k P_k i^{k-1} ω^k.
  const auto &P = surf_rz_op.GetRobinQuotient(idx);
  const auto n_P = static_cast<int>(P.size());
  auto P_coeff = [&](int k) -> double  // P_k, coefficient of s^k (stored high-first)
  { return (k < n_P) ? P[n_P - 1 - k] : 0.0; };
  fit.alpha0c = P_coeff(0) / std::complex<double>(0.0, 1.0);
  fit.alpha1c = P_coeff(1);
  fit.alpha2c = P_coeff(2) * std::complex<double>(0.0, 1.0);

  constexpr int n_fit = 30;
  constexpr int n_dense = 200;
  auto fit_omegas = SampleChebyshevLobatto(sweep_omega_min, sweep_omega_max, n_fit);
  auto dense_omegas = SampleChebyshevGauss(sweep_omega_min, sweep_omega_max, n_dense);
  auto f_rem = [&surf_rz_op, idx](double w) -> std::complex<double>
  {
    return surf_rz_op.EvalRobinRemainder(idx, std::complex<double>(0.0, w)) /
           std::complex<double>(0.0, 1.0);
  };
  auto f_full = [&surf_rz_op, idx](double w) -> std::complex<double>
  {
    return surf_rz_op.EvalRobinCoefficient(idx, std::complex<double>(0.0, w)) /
           std::complex<double>(0.0, 1.0);
  };

  // AAA directly on the strictly proper remainder (no LSQ residual step — the
  // polynomial part is already exact from the long division).
  Eigen::VectorXcd F_aaa(n_fit);
  double max_abs_truth = 0.0;
  for (int i = 0; i < n_fit; i++)
  {
    F_aaa(i) = f_rem(fit_omegas[i]);
    max_abs_truth = std::max(max_abs_truth, std::abs(f_full(fit_omegas[i])));
  }
  auto pr = FitResidualPoles(fit_omegas, F_aaa, waveport_synthesis_tol, max_abs_truth,
                             waveport_synthesis_order_max);
  // The exact remainder is strictly proper, but preserve any small polynomial quotient
  // introduced by the numerical AAA conversion instead of silently dropping it.
  FoldAAAPolynomialIntoQuadratic(pr, fit.alpha0c, fit.alpha1c, fit.alpha2c);

  // Dense-grid residual of the full reconstruction (polynomial part + pole-residue).
  double max_rel = 0.0;
  for (int i = 0; i < n_dense; i++)
  {
    const double w = dense_omegas[i];
    max_rel = std::max(max_rel, std::abs(EvaluatePolyPlusPoles(fit.alpha0c, fit.alpha1c,
                                                               fit.alpha2c, pr, w) -
                                         f_full(w)));
  }
  fit.rel_err_polynomial = fit.rel_err_augmented =
      (max_abs_truth > 0.0) ? max_rel / max_abs_truth : 0.0;

  auto [poles, residues] = ToPoleResidueVectors(pr);
  fit.aux = MakeAuxBlock(label, Mp_r, poles, residues, waveport_synthesis_rank_tol);
  fit.regime = fit.aux ? WavePortRegime::Augmented : WavePortRegime::Polynomial;
  Mpi::Print(" {}: exact-split synthesis residual {:.3e} (tol {:.3e}, {:d} pole{})\n",
             label, fit.rel_err_augmented, waveport_synthesis_tol, pr.poles.size(),
             pr.poles.size() == 1 ? "" : "s");
  if (fit.rel_err_augmented > waveport_synthesis_tol)
  {
    Mpi::Warning("{}: rational impedance synthesis residual {:.3e} exceeds tolerance "
                 "{:.3e}!\n",
                 label, fit.rel_err_augmented, waveport_synthesis_tol);
  }
  return fit;
}

RomOperator::AugmentedPencil RomOperator::BuildAugmentedPencil(
    const Eigen::MatrixXcd &Kr_total, const Eigen::MatrixXcd &Cr_total,
    const Eigen::MatrixXcd &Mr_total, const std::vector<WavePortAuxBlock> &aux_blocks,
    std::vector<std::string> &aux_labels)
{
  // Per pole–residue (pₖ, rₖ) at port p, with the signed symmetric factorization
  // M_proj = Σⱼ λⱼ uⱼuⱼᵀ:
  // for each kept singular vector j, allocate one aux state sₖⱼ. On the state [v; sₖⱼ],
  // the augmented pencil is
  //
  //   A(ω) = [ K      a·uⱼ ] + iω · [ C   0  ] − ω² · [ M   0 ]
  //          [ a·uⱼᵀ  -pₖ  ]        [ 0   -i ]        [ 0   0 ]
  //
  // (K/C/M are the n×n base blocks; the aux row/column carries K_va = K_avᵀ = a·uⱼ,
  // K_aa = -pₖ, C_aa = -i, M_aa = 0). Eliminating sₖⱼ gives the Schur contribution
  // -a²·uⱼuⱼᵀ·v/(ω - pₖ). Choosing a² = -i·rₖ·λⱼ makes the sum over j equal
  // +i·rₖ·M_proj·v/(ω - pₖ), the desired rₖ-residue contribution to the unaugmented
  // wave-port pencil i·kₙ(ω)·Mp_r·v. Both couplings are the same a·uⱼ (no conjugation),
  // so the augmented pencil stays complex-symmetric (not Hermitian — downstream
  // eigensolvers handle this fine).
  std::size_t n_aux_total = 0;
  for (const auto &blk : aux_blocks)
  {
    n_aux_total += blk.poles.size() * blk.weights.size();
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
      for (std::size_t j = 0; j < blk.weights.size(); j++)
      {
        std::complex<double> coupling =
            std::sqrt(std::complex<double>(0.0, -1.0) * rk * blk.weights[j]);
        aug.Kr(aux_row, aux_row) = -pk;
        aug.Cr(aux_row, aux_row) = std::complex<double>(0.0, -1.0);
        for (long i = 0; i < n_v; i++)
        {
          aug.Kr(i, aux_row) = coupling * blk.u_dirs[j](i);
          aug.Kr(aux_row, i) = coupling * blk.u_dirs[j](i);
        }
        const std::string prefix =
            blk.label.empty() ? fmt::format("waveport_{:d}", blk.port_idx) : blk.label;
        aux_labels.push_back(fmt::format("{}_p{:d}d{:d}", prefix, k, j));
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
  // The scan above assumes the wave-port rows form one contiguous block directly after
  // the lumped-port rows (AddWavePortModesForSynthesis runs immediately after
  // AddLumpedPortModesForSynthesis, before any HDM sample). Verify no waveport_ label
  // appears later, which would mean the block was split and undercounted. Row existence is
  // independent of Active: an inactive included port is still an unloaded synthesis node.
  for (long j = n_port_modes + n_waveport_rows; j < static_cast<long>(v_node_label.size());
       j++)
  {
    MFEM_VERIFY(v_node_label[j].rfind("waveport_", 0) != 0,
                "Wave-port basis rows are not contiguous after the lumped-port rows!");
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
  // FitWavePortDispersion does the per-port sampling, fitting and regime selection;
  // the Apply*FitCorrections helpers accumulate the top-left contributions;
  // BuildAugmentedPencil extends the pencil with the aux states from collected fit results.
  Eigen::MatrixXcd Kr_total_corr = Eigen::MatrixXcd::Zero(Kr.rows(), Kr.cols());
  Eigen::MatrixXcd Cr_total_corr = Eigen::MatrixXcd::Zero(Kr.rows(), Kr.cols());
  Eigen::MatrixXcd Mr_total_corr = Eigen::MatrixXcd::Zero(Kr.rows(), Kr.cols());
  std::vector<WavePortAuxBlock> aux_blocks_total;
  struct PendingPortLoad
  {
    std::string label;
    Eigen::MatrixXcd Kr_corr;
    Eigen::MatrixXcd Cr_corr;
    Eigen::MatrixXcd Mr_corr;
    std::vector<WavePortAuxBlock> aux_blocks;
  };
  std::vector<PendingPortLoad> pending_port_loads;
  if (!Mwp_p_r.empty() && sweep_omega_max > sweep_omega_min)
  {
    // Sample kₙ(ω)/M(ω) below at a tight, EigenTol-independent EVP tolerance so the fit is
    // not resolved below the port-mode accuracy floor (see WAVEPORT_SYNTHESIS_EIG_TOL).
    space_op.GetWavePortOp().SetSynthesisEigTol(WAVEPORT_SYNTHESIS_EIG_TOL,
                                                WAVEPORT_SYNTHESIS_KSP_TOL);
    for (auto &[port_idx, Mp_r] : Mwp_p_r)
    {
      auto fit = FitWavePortDispersion(port_idx, Mp_r);
      out.wave_port_fits.push_back(fit);

      // Inactive included ports are unloaded synthesis terminals. They need the fit for
      // rom-port-reference.csv, but their Robin termination must not be added to either the
      // loaded total pencil or the removable per-port load.
      const auto &port_data = space_op.GetWavePortOp().GetPort(port_idx);
      if (!port_data.active)
      {
        continue;
      }

      PendingPortLoad port_load;
      port_load.label = fmt::format("waveport_{:d}_re", port_idx);
      port_load.Kr_corr = Eigen::MatrixXcd::Zero(Kr.rows(), Kr.cols());
      port_load.Cr_corr = Eigen::MatrixXcd::Zero(Kr.rows(), Kr.cols());
      port_load.Mr_corr = Eigen::MatrixXcd::Zero(Kr.rows(), Kr.cols());
      ApplyPolynomialFitCorrections(fit, Mp_r, Kr_total_corr, Cr_total_corr, Mr_total_corr);
      ApplyPolynomialFitCorrections(fit, Mp_r, port_load.Kr_corr, port_load.Cr_corr,
                                    port_load.Mr_corr);
      if (fit.aux)
      {
        aux_blocks_total.push_back(*fit.aux);
        port_load.aux_blocks.push_back(*fit.aux);
      }
      pending_port_loads.push_back(std::move(port_load));
    }
  }

  // Wave-port modal correction W = Σ_ports (W_full − W_scalar), folded into synthesis. Per
  // active port, sample the recomputed reduced n×H vectors ŝ_full(ω), ŝ_scalar(ω) across
  // the band, build an orthonormal modal subspace Q (truncated SVD), and write
  // Wᵣ(ω)=Q·M(ω)·Qᵀ with M(ω)=g_full·a_full a_fullᵀ + g_scalar·a_scalar a_scalarᵀ,
  // aₓ=Qᴴŝₓ(ω), an r×r complex-symmetric matrix. Since Q(Qᴴŝ)(Qᴴŝ)ᵀQᵀ=ŝŝᵀ for ŝ∈range(Q),
  // this reproduces the eigenmode-path Wᵣ(ω) as Q resolves the mode-shape rotation (a
  // center-frozen rank-2 span cannot). Each unique entry M_pq(ω) is fit with
  // FitScalarDispersion, added to the loaded pencil and the matching removable per-port
  // load.
  if (!Mwp_p_r.empty() && sweep_omega_max > sweep_omega_min)
  {
    const double omega_c = 0.5 * (sweep_omega_min + sweep_omega_max);
    const int nr = static_cast<int>(Kr.rows());
    const auto build_omegas = SampleChebyshevLobatto(sweep_omega_min, sweep_omega_max,
                                                     WAVEPORT_SYNTHESIS_SUBSPACE_SAMPLES);
    auto sample_reduced =
        [this, nr](int port_idx, std::complex<double> w, Eigen::VectorXcd &sf,
                   Eigen::VectorXcd &ss) -> WavePortOperator::ModalCorrectionSample
    {
      auto smp = space_op.SampleModalCorrectionVectors(port_idx, w);
      sf = Eigen::VectorXcd::Zero(nr);
      ss = Eigen::VectorXcd::Zero(nr);
      if (smp.active)
      {
        ProjectVecInternal(space_op.GetComm(), V, *smp.s_full, sf, 0);
        ProjectVecInternal(space_op.GetComm(), V, *smp.s_scalar, ss, 0);
      }
      return smp;
    };
    for (int port_idx : space_op.GetModalCorrectionSynthesisPorts(omega_c))
    {
      // Stack reduced, unit-normalized samples as columns; truncated-SVD left vectors span
      // Q.
      Eigen::MatrixXcd cols(nr, 2 * static_cast<int>(build_omegas.size()));
      int ncol = 0;
      for (double wr : build_omegas)
      {
        Eigen::VectorXcd sf, ss;
        auto smp = sample_reduced(port_idx, std::complex<double>(wr, 0.0), sf, ss);
        if (!smp.active)
        {
          continue;
        }
        if (double n = sf.norm(); n > 0.0)
        {
          cols.col(ncol++) = sf / n;
        }
        if (double n = ss.norm(); n > 0.0)
        {
          cols.col(ncol++) = ss / n;
        }
      }
      if (ncol == 0)
      {
        continue;
      }
      cols.conservativeResize(nr, ncol);
      Eigen::BDCSVD<Eigen::MatrixXcd, Eigen::ComputeThinU> svd;
      svd.compute(cols);
      const auto &sv = svd.singularValues();
      int r = 0;
      while (r < sv.size() && sv(r) > waveport_synthesis_rank_tol * sv(0))
      {
        r++;
      }
      if (r == 0)
      {
        continue;
      }
      if (r > WAVEPORT_SYNTHESIS_SUBSPACE_RANK_MAX)
      {
        Mpi::Warning(" Wave port {:d} modal subspace rank {:d} capped at {:d}; synthesis "
                     "accuracy may be reduced\n",
                     port_idx, r, WAVEPORT_SYNTHESIS_SUBSPACE_RANK_MAX);
        r = WAVEPORT_SYNTHESIS_SUBSPACE_RANK_MAX;
      }
      const Eigen::MatrixXcd Q = svd.matrixU().leftCols(r);
      Mpi::Print(" Wave port {:d} modal-correction synthesis subspace rank {:d}\n",
                 port_idx, r);
      // M(ω), shared by a common-pole matrix-valued rational fit. Fitting every entry
      // independently is not covariant under rotations of the SVD basis Q and creates a
      // large nonminimal realization with nearly repeated poles. Instead obtain one AAA
      // denominator from the vectorized matrix residual, then solve a linear least-squares
      // problem for all polynomial and residue matrices using those common poles.
      std::map<std::pair<double, double>, Eigen::MatrixXcd> M_cache;
      auto eval_M = [&](std::complex<double> w) -> const Eigen::MatrixXcd &
      {
        const auto key = std::make_pair(w.real(), w.imag());
        auto it = M_cache.find(key);
        if (it == M_cache.end())
        {
          Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(r, r);
          Eigen::VectorXcd sf, ss;
          auto smp = sample_reduced(port_idx, w, sf, ss);
          if (smp.active)
          {
            const Eigen::VectorXcd af = Q.adjoint() * sf, as = Q.adjoint() * ss;
            M = smp.g_full * (af * af.transpose()) + smp.g_scalar * (as * as.transpose());
            M = 0.5 * (M + M.transpose()).eval();
          }
          it = M_cache.emplace(key, std::move(M)).first;
        }
        return it->second;
      };
      auto pack_matrix = [r](const Eigen::MatrixXcd &M)
      {
        Eigen::RowVectorXcd v(r * r);
        for (int i = 0; i < r; i++)
        {
          for (int j = 0; j < r; j++)
          {
            v(i * r + j) = M(i, j);
          }
        }
        return v;
      };
      auto unpack_matrix = [r](const Eigen::RowVectorXcd &v)
      {
        Eigen::MatrixXcd M(r, r);
        for (int i = 0; i < r; i++)
        {
          for (int j = 0; j < r; j++)
          {
            M(i, j) = v(i * r + j);
          }
        }
        return (0.5 * (M + M.transpose())).eval();
      };

      constexpr int n_fit = 30;
      constexpr int n_dense = 60;
      const auto fit_real = SampleChebyshevLobatto(sweep_omega_min, sweep_omega_max, n_fit);
      const auto dense_real =
          SampleChebyshevGauss(sweep_omega_min, sweep_omega_max, n_dense);
      Eigen::VectorXcd z(n_fit);
      Eigen::MatrixXcd Y(n_fit, r * r), poly_design(n_fit, 3);
      for (int i = 0; i < n_fit; i++)
      {
        const std::complex<double> w(fit_real[static_cast<std::size_t>(i)], 0.0);
        z(i) = w;
        Y.row(i) = pack_matrix(eval_M(w));
        poly_design(i, 0) = 1.0;
        poly_design(i, 1) = w;
        poly_design(i, 2) = w * w;
      }
      const Eigen::MatrixXcd poly_coeff = poly_design.colPivHouseholderQr().solve(Y);
      const Eigen::MatrixXcd residual = Y - poly_design * poly_coeff;
      auto aaa = utils::RunAAA(z, residual, waveport_synthesis_tol,
                               std::max<std::size_t>(waveport_synthesis_order_max, 1));
      auto denominator = utils::AAAToPoleResidue(aaa);
      const int n_poles = static_cast<int>(denominator.poles.size());

      Eigen::MatrixXcd design(n_fit, 3 + n_poles);
      design.leftCols(3) = poly_design;
      for (int i = 0; i < n_fit; i++)
      {
        for (int k = 0; k < n_poles; k++)
        {
          design(i, 3 + k) = 1.0 / (z(i) - denominator.poles(k));
        }
      }
      const Eigen::MatrixXcd coeff = design.colPivHouseholderQr().solve(Y);
      const Eigen::MatrixXcd P0 = unpack_matrix(coeff.row(0));
      const Eigen::MatrixXcd P1 = unpack_matrix(coeff.row(1));
      const Eigen::MatrixXcd P2 = unpack_matrix(coeff.row(2));

      const Eigen::MatrixXcd P0_full = Q * P0 * Q.transpose();
      const Eigen::MatrixXcd P1_full = Q * P1 * Q.transpose();
      const Eigen::MatrixXcd P2_full = Q * P2 * Q.transpose();
      Kr_total_corr += P0_full;
      Cr_total_corr += std::complex<double>(0.0, -1.0) * P1_full;
      Mr_total_corr -= P2_full;

      const auto load_label = fmt::format("waveport_{:d}_re", port_idx);
      auto pl = std::find_if(pending_port_loads.begin(), pending_port_loads.end(),
                             [&](const auto &p) { return p.label == load_label; });
      if (pl != pending_port_loads.end())
      {
        pl->Kr_corr += P0_full;
        pl->Cr_corr += std::complex<double>(0.0, -1.0) * P1_full;
        pl->Mr_corr -= P2_full;
      }

      const std::complex<double> imag_unit(0.0, 1.0);
      for (int k = 0; k < n_poles; k++)
      {
        const Eigen::MatrixXcd R = unpack_matrix(coeff.row(3 + k));
        const Eigen::MatrixXcd R_full = Q * R * Q.transpose();
        const double residue_norm = R_full.norm();
        const std::array<std::pair<Eigen::MatrixXd, std::complex<double>>, 2> parts = {
            std::make_pair(R_full.real(), std::complex<double>(0.0, -1.0)),
            std::make_pair(R_full.imag(), std::complex<double>(1.0, 0.0))};
        const char *suffix[2] = {"re", "im"};
        for (int part = 0; part < 2; part++)
        {
          if (residue_norm == 0.0 ||
              parts[part].first.norm() <= WAVEPORT_SYNTHESIS_MODAL_PART_TOL * residue_norm)
          {
            continue;
          }
          const Eigen::MatrixXcd Mp =
              imag_unit * parts[part].first.cast<std::complex<double>>();
          auto blk = MakeAuxBlock(
              fmt::format("waveport_{:d}_modal_p{:d}_{}", port_idx, k, suffix[part]), Mp,
              {denominator.poles(k)}, {parts[part].second}, waveport_synthesis_rank_tol);
          if (blk)
          {
            aux_blocks_total.push_back(*blk);
            if (pl != pending_port_loads.end())
            {
              pl->aux_blocks.push_back(*blk);
            }
          }
        }
      }

      double max_err = 0.0, max_ref = 0.0;
      for (const double wr : dense_real)
      {
        const std::complex<double> w(wr, 0.0);
        Eigen::MatrixXcd M_fit = P0 + w * P1 + w * w * P2;
        for (int k = 0; k < n_poles; k++)
        {
          M_fit += unpack_matrix(coeff.row(3 + k)) / (w - denominator.poles(k));
        }
        max_err = std::max(max_err, (M_fit - eval_M(w)).norm());
        max_ref = std::max(max_ref, eval_M(w).norm());
      }
      const double rel_err = (max_ref > 0.0) ? max_err / max_ref : 0.0;
      Mpi::Print(" Wave port {:d} common-pole modal synthesis: rank {:d}, {:d} pole{}, "
                 "matrix residual {:.3e} (tol {:.3e})\n",
                 port_idx, r, n_poles, n_poles == 1 ? "" : "s", rel_err,
                 waveport_synthesis_tol);
      if (rel_err > waveport_synthesis_tol)
      {
        Mpi::Warning("Wave port {:d} common-pole modal synthesis residual {:.3e} exceeds "
                     "AdaptiveTol={:.3e}!\n",
                     port_idx, rel_err, waveport_synthesis_tol);
      }
    }
  }

  // Other frequency-dependent BCs, folded into the same aug-pencil form (each contributes
  // i·f(ω)·M_proj·v with M_proj projected onto the basis and the imaginary slot carrying
  // the i). Only meaningful with a nonzero sweep band.
  if (sweep_omega_max > sweep_omega_min)
  {
    // 2nd-order farfield ABC: i·(0.5/ω)·M_ff. f(ω) = 0.5/ω is EXACTLY a single pole at ω=0
    // with residue 0.5 — inject it analytically (no polynomial part, no fit), as a frozen
    // aux block. This is the synthesis analogue of the NLEPS frozen-ABC seed.
    if (M_ff_ && M_ff_r.rows() == Kr.rows())
    {
      auto blk =
          MakeAuxBlock("farfield", M_ff_r, {std::complex<double>(0.0, 0.0)},
                       {std::complex<double>(0.5, 0.0)}, waveport_synthesis_rank_tol);
      if (blk)
      {
        aux_blocks_total.push_back(*blk);
        Mpi::Print(" Second-order farfield ABC: folded into synthesis as 1 pole at ω=0 "
                   "(residue 0.5)\n");
      }
    }
    // Surface conductivity, one group at a time: f_g(ω) = ω/Z_g(ω) (the i is the implicit
    // slot factor; EvaluateScalar returns i·ω/Z, so f = EvaluateScalar/i). Fit complex
    // poly + AAA via FitScalarDispersion.
    const auto &surf_op = space_op.GetSurfaceConductivityOp();
    for (std::size_t g = 0; g < Asig_g_.size(); g++)
    {
      if (!Asig_g_[g] || Asig_g_r[g].rows() != Kr.rows())
      {
        continue;
      }
      const auto label = fmt::format("surfsigma_{:d}", g);
      auto f = [&surf_op, g](std::complex<double> omega) -> std::complex<double>
      { return surf_op.EvaluateScalar(g, omega) / std::complex<double>(0.0, 1.0); };
      auto fit = FitScalarDispersion(label, Asig_g_r[g], f, /*allow_augment=*/true);
      ApplyComplexPolynomialFitCorrections(fit.alpha0c, fit.alpha1c, fit.alpha2c,
                                           Asig_g_r[g], Kr_total_corr, Cr_total_corr,
                                           Mr_total_corr);
      if (fit.aux)
      {
        aux_blocks_total.push_back(*fit.aux);
      }
    }
    // Rational surface impedance, one boundary at a time: the full Robin term is
    // g(iω)·M_b = i·f_b(ω)·M_b with f_b(ω) = g(iω)/i and g(s) = s·D(s)/N(s). The exact
    // long-division split (polynomial part into the pencil, AAA directly on the strictly
    // proper remainder) reproduces the rational coefficient to rounding.
    for (std::size_t b = 0; b < Arz_b_.size(); b++)
    {
      if (!Arz_b_[b] || Arz_b_r[b].rows() != Kr.rows())
      {
        continue;
      }
      const auto label = fmt::format("rationalz_{:d}", b);
      auto fit = FitRationalImpedanceDispersion(label, Arz_b_r[b], static_cast<int>(b));
      ApplyComplexPolynomialFitCorrections(fit.alpha0c, fit.alpha1c, fit.alpha2c,
                                           Arz_b_r[b], Kr_total_corr, Cr_total_corr,
                                           Mr_total_corr);
      if (fit.aux)
      {
        aux_blocks_total.push_back(*fit.aux);
      }
    }
  }

  // Polynomial-only matrices (basis dim n × n). The legacy matrices are loaded by the
  // matched port/reference realization. Per-port load matrices are emitted separately below
  // so downstream tools can remove internal port loads and add back only external loads
  // during matrix-level network assembly.
  Eigen::MatrixXcd Kr_total = Kr + Kr_total_corr;
  Eigen::MatrixXcd Mr_total = Mr + Mr_total_corr;
  Eigen::MatrixXcd Cr_total = Cr_total_corr;
  if (C)
  {
    Cr_total += Cr;
  }

  auto aug =
      BuildAugmentedPencil(Kr_total, Cr_total, Mr_total, aux_blocks_total, out.aux_labels);

  // v_d port-row scaling: extend with 1's for aux rows (no port-impedance scaling on
  // aux states — they're internal circuit nodes).
  auto unit_henry_inv = 1.0 / units.GetScaleFactor<Units::ValueType::INDUCTANCE>();
  auto unit_farad = units.GetScaleFactor<Units::ValueType::CAPACITANCE>();
  auto unit_ohm_inv = 1.0 / units.GetScaleFactor<Units::ValueType::IMPEDANCE>();

  auto normalize_augmented =
      [&](const AugmentedPencil &aug_in, std::unique_ptr<mat_t> &L_inv,
          std::unique_ptr<mat_t> &R_inv, std::unique_ptr<mat_t> &C_out)
  {
    const long n_v = Kr.rows();
    const long n_aug = aug_in.Kr.rows();
    Eigen::VectorXcd v_conc_aug = Eigen::VectorXcd::Ones(n_aug);
    for (long j = 0; j < n_v; j++)
    {
      v_conc_aug(j) = v_conc(j);
    }
    auto v_d_aug = v_conc_aug.asDiagonal();

    L_inv =
        std::make_unique<mat_t>((unit_henry_inv * v_d_aug * aug_in.Kr * v_d_aug).eval());
    C_out = std::make_unique<mat_t>((unit_farad * v_d_aug * aug_in.Mr * v_d_aug).eval());
    // Emit R⁻¹ whenever there's any dissipative contribution: lumped resistance, surface
    // conductivity, wave-port α₁, or aux-state damping.
    if (aug_in.Cr.cwiseAbs().maxCoeff() > 0.0)
    {
      R_inv =
          std::make_unique<mat_t>((unit_ohm_inv * v_d_aug * aug_in.Cr * v_d_aug).eval());
    }
  };

  normalize_augmented(aug, out.L_inv, out.R_inv, out.C);

  std::vector<std::string> total_labels = v_node_label;
  for (const auto &lab : out.aux_labels)
  {
    total_labels.push_back(lab);
  }

  auto make_zero_augmented = [](long n) -> AugmentedPencil
  {
    AugmentedPencil z;
    z.Kr = Eigen::MatrixXcd::Zero(n, n);
    z.Cr = Eigen::MatrixXcd::Zero(n, n);
    z.Mr = Eigen::MatrixXcd::Zero(n, n);
    return z;
  };

  auto embed_augmented = [&](const AugmentedPencil &local_aug,
                             const std::vector<std::string> &local_aux_labels)
  {
    const long n_v = Kr.rows();
    const long n_total = aug.Kr.rows();
    auto full_aug = make_zero_augmented(n_total);
    auto global_index = [&](long i) -> long
    {
      if (i < n_v)
      {
        return i;
      }
      const long aux_i = i - n_v;
      MFEM_VERIFY(aux_i < static_cast<long>(local_aux_labels.size()),
                  "Malformed port-load auxiliary label list!");
      const long gi = LabelIndex(total_labels, local_aux_labels[aux_i]);
      MFEM_VERIFY(gi >= 0, "Missing port-load auxiliary row in total PROM labels!");
      return gi;
    };
    for (long i = 0; i < local_aug.Kr.rows(); i++)
    {
      const long gi = global_index(i);
      for (long j = 0; j < local_aug.Kr.cols(); j++)
      {
        const long gj = global_index(j);
        full_aug.Kr(gi, gj) += local_aug.Kr(i, j);
        full_aug.Cr(gi, gj) += local_aug.Cr(i, j);
        full_aug.Mr(gi, gj) += local_aug.Mr(i, j);
      }
    }
    return full_aug;
  };

  for (const auto &pending : pending_port_loads)
  {
    std::vector<std::string> local_aux_labels;
    auto local_aug = BuildAugmentedPencil(pending.Kr_corr, pending.Cr_corr, pending.Mr_corr,
                                          pending.aux_blocks, local_aux_labels);
    auto full_aug = embed_augmented(local_aug, local_aux_labels);

    NormalizedMatrices::PortLoad load;
    load.label = pending.label;
    normalize_augmented(full_aug, load.L_inv, load.R_inv, load.C);
    out.port_loads.push_back(std::move(load));
  }

  // Lumped-port R/L/C boundary conditions are part of the legacy loaded pencil. Export
  // their terminal admittance as per-port load matrices so downstream tools can form the
  // connectable device by subtracting selected port loads from the legacy total matrices.
  const long n_total = aug.Kr.rows();
  for (const auto &[port_idx, port_data] : space_op.GetLumpedPortOp())
  {
    if (!port_data.active || !port_data.include_in_synthesis)
    {
      continue;
    }
    const auto label = fmt::format("port_{:d}_re", port_idx);
    const long row = LabelIndex(total_labels, label);
    MFEM_VERIFY(row >= 0, "Missing synthesized lumped-port row for port-load export!");

    NormalizedMatrices::PortLoad load;
    load.label = label;
    load.L_inv = std::make_unique<mat_t>(mat_t::Zero(n_total, n_total));
    load.C = std::make_unique<mat_t>(mat_t::Zero(n_total, n_total));
    if (std::abs(port_data.L) > 0.0)
    {
      (*load.L_inv)(row, row) += unit_henry_inv / port_data.L;
    }
    if (std::abs(port_data.C) > 0.0)
    {
      (*load.C)(row, row) += unit_farad * port_data.C;
    }
    if (std::abs(port_data.R) > 0.0)
    {
      load.R_inv = std::make_unique<mat_t>(mat_t::Zero(n_total, n_total));
      (*load.R_inv)(row, row) += unit_ohm_inv / port_data.R;
    }
    out.port_loads.push_back(std::move(load));
  }

  return out;
}

void RomOperator::PrintPortReferenceData(const Units &units, const fs::path &post_dir,
                                         const NormalizedMatrices &matrices) const
{
  if (sweep_omega_samples.empty())
  {
    return;
  }

  enum class RefType
  {
    Lumped,
    Wave
  };
  struct RefPort
  {
    RefType type;
    int port_idx;
    std::string label;
    const WavePortDispersionFit *wave_fit = nullptr;
  };

  std::vector<RefPort> refs;
  refs.reserve(NumSynthesisPortModes() + NumSynthesisWavePortModes());

  for (const auto &[port_idx, port_data] : space_op.GetLumpedPortOp())
  {
    if (!port_data.include_in_synthesis)
    {
      continue;
    }
    const auto label = fmt::format("port_{:d}_re", port_idx);
    if (LabelIndex(v_node_label, label) >= 0)
    {
      refs.push_back({RefType::Lumped, port_idx, label, nullptr});
    }
  }
  for (const auto &[port_idx, port_data] : space_op.GetWavePortOp())
  {
    if (!port_data.include_in_synthesis)
    {
      continue;
    }
    const int wp_idx = port_idx;
    const auto fit_it =
        std::find_if(matrices.wave_port_fits.begin(), matrices.wave_port_fits.end(),
                     [wp_idx](const auto &fit) { return fit.port_idx == wp_idx; });
    const auto label = fmt::format("waveport_{:d}_re", port_idx);
    if (fit_it != matrices.wave_port_fits.end() && LabelIndex(v_node_label, label) >= 0)
    {
      refs.push_back({RefType::Wave, port_idx, label, &(*fit_it)});
    }
  }
  if (refs.empty())
  {
    return;
  }

  const double unit_GHz =
      units.Dimensionalize<Units::ValueType::FREQUENCY>(1.0) / (2.0 * M_PI);
  const double unit_ohm_inv = 1.0 / units.GetScaleFactor<Units::ValueType::IMPEDANCE>();

  // Physical frequency scale s_phys = iω·ω0 with ω0 = unit_henry_inv/unit_ohm_inv, so
  // Y(s_phys) = L⁻¹/s_phys + R⁻¹ + s_phys·C reproduces the old scalar reference and also
  // carries the W modal correction and aux states baked into the per-port load pencil.
  const double unit_henry_inv = 1.0 / units.GetScaleFactor<Units::ValueType::INDUCTANCE>();
  const double omega0 = unit_henry_inv / unit_ohm_inv;

  std::vector<std::string> total_labels = v_node_label;
  for (const auto &lab : matrices.aux_labels)
  {
    total_labels.push_back(lab);
  }

  auto wave_y_ref = [&](const RefPort &ref, double omega) -> std::complex<double>
  {
    if (!(omega > 0.0))
    {
      return {0.0, 0.0};
    }
    const auto load_it =
        std::find_if(matrices.port_loads.begin(), matrices.port_loads.end(),
                     [&ref](const auto &pl) { return pl.label == ref.label; });
    MFEM_VERIFY(load_it != matrices.port_loads.end(),
                "Missing wave-port load pencil for port reference output!");
    const long phys = LabelIndex(total_labels, ref.label);
    MFEM_VERIFY(phys >= 0, "Missing wave-port physical row for port reference output!");

    // Full-size port admittance from the load pencil (already in physical units).
    const std::complex<double> s_phys = 1i * omega * omega0;
    Eigen::MatrixXcd Y = (*load_it->L_inv) / s_phys + s_phys * (*load_it->C);
    if (load_it->R_inv)
    {
      Y += *load_it->R_inv;
    }

    // Order the physical port row first, then any internal (imag/aux) rows carrying load,
    // and Schur-complement the internal rows onto the physical row for the terminal
    // admittance.
    std::vector<long> rows{phys};
    for (long i = 0; i < Y.rows(); i++)
    {
      if (i != phys &&
          (Y.row(i).cwiseAbs().maxCoeff() > 0.0 || Y.col(i).cwiseAbs().maxCoeff() > 0.0))
      {
        rows.push_back(i);
      }
    }
    const long n = static_cast<long>(rows.size());
    Eigen::MatrixXcd Ys(n, n);
    for (long i = 0; i < n; i++)
    {
      for (long j = 0; j < n; j++)
      {
        Ys(i, j) = Y(rows[i], rows[j]);
      }
    }
    std::complex<double> y_eff = Ys(0, 0);
    if (n > 1)
    {
      // The internal block spans a huge dynamic range (L⁻¹ ~ 1e14, C ~ 1e-14), so the
      // default rank threshold spuriously flags it as deficient; keep all pivots and guard
      // only against a genuinely singular solve producing non-finite entries.
      const long ni = n - 1;
      Eigen::FullPivLU<Eigen::MatrixXcd> lu(Ys.bottomRightCorner(ni, ni));
      lu.setThreshold(1.0e-300);
      const std::complex<double> corr =
          (Ys.block(0, 1, 1, ni) * lu.solve(Ys.block(1, 0, ni, 1)))(0, 0);
      if (std::isfinite(corr.real()) && std::isfinite(corr.imag()))
      {
        y_eff -= corr;
      }
    }
    return y_eff;
  };

  Mpi::Print(" Printing PROM port reference admittance to disk.\n");
  auto out = TableWithCSVFile(post_dir / "rom-port-reference.csv");
  out.table.col_options.float_precision = 17;
  out.table.reserve(sweep_omega_samples.size(), 1 + 4 * refs.size());
  out.table.insert("idx", "f (GHz)", -1, 0, std::size_t{12}, "");
  for (const auto &ref : refs)
  {
    const auto key = ref.label;
    out.table.insert(fmt::format("re_yref_{}", key),
                     fmt::format("Re{{Y_ref[{}]}} (S)", ref.label));
    out.table.insert(fmt::format("im_yref_{}", key),
                     fmt::format("Im{{Y_ref[{}]}} (S)", ref.label));
    out.table.insert(fmt::format("re_zref_{}", key),
                     fmt::format("Re{{Z_ref[{}]}} (Ohm)", ref.label));
    out.table.insert(fmt::format("im_zref_{}", key),
                     fmt::format("Im{{Z_ref[{}]}} (Ohm)", ref.label));
  }

  for (const auto omega : sweep_omega_samples)
  {
    out.table["idx"] << omega * unit_GHz;
    for (const auto &ref : refs)
    {
      std::complex<double> y_ref = 0.0;
      if (ref.type == RefType::Lumped)
      {
        // Use the same real reference resistance as Palace's own S-parameter measurement
        // (GetExcitationRefResistance: R for a resistive port, the unit internal reference
        // for a purely reactive R == 0 port). This keeps the table consistent with its
        // documented purpose — converting synthesized admittances to Palace/Kurokawa S —
        // after reactive lumped excitation referenced direct S to the real R_ref. The
        // physical (complex) RLC port load is exported separately in the rom-portload-*
        // matrices. For resistive ports this coincides with the characteristic impedance.
        const auto &port_data = space_op.GetLumpedPortOp().GetPort(ref.port_idx);
        y_ref = unit_ohm_inv / port_data.GetExcitationRefResistance();
      }
      else
      {
        y_ref = wave_y_ref(ref, omega);
      }
      const std::complex<double> z_ref =
          (std::abs(y_ref) > 0.0) ? (1.0 / y_ref) : std::complex<double>{0.0, 0.0};
      out.table[fmt::format("re_yref_{}", ref.label)] << y_ref.real();
      out.table[fmt::format("im_yref_{}", ref.label)] << y_ref.imag();
      out.table[fmt::format("re_zref_{}", ref.label)] << z_ref.real();
      out.table[fmt::format("im_zref_{}", ref.label)] << z_ref.imag();
    }
  }
  out.WriteFullTableTrunc();
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

  // Eigenvalue estimates of the synthesized system and their HDM eigenpair errors. The
  // QZ runs redundantly on every rank (replicated dense matrices); the HDM error
  // evaluation is collective (prolongation, operator application, norms), so it must
  // also run on every rank before the root-only output below.
  const double fmin_GHz =
      units.Dimensionalize<Units::ValueType::FREQUENCY>(sweep_omega_min) / (2.0 * M_PI);
  const double fmax_GHz =
      units.Dimensionalize<Units::ValueType::FREQUENCY>(sweep_omega_max) / (2.0 * M_PI);
  auto eigs = ComputeEigenvalueEstimates(*matrices.L_inv, matrices.R_inv.get(), *matrices.C,
                                         fmin_GHz, fmax_GHz, GetReducedDimension());
  ComputeEigenvalueEstimateErrors(units, eigs);

  // Diagnostic scattering reconstruction from the synthesized loaded pencil plus the
  // frequency-dependent wave-port input/output coupling vectors. A rotating hybrid mode
  // cannot in general be represented by a fixed terminal selector alone.
  std::vector<int> coupled_ports;
  for (const auto &[port_idx, port] : space_op.GetWavePortOp())
  {
    if (port.active && port.include_in_synthesis)
    {
      coupled_ports.push_back(port_idx);
    }
  }
  std::vector<std::vector<std::complex<double>>> coupled_s;
  if (!coupled_ports.empty())
  {
    const long nr = static_cast<long>(V.size());
    const long na = matrices.L_inv->rows();
    Eigen::VectorXd d = Eigen::VectorXd::Ones(na);
    const long n_lumped = static_cast<long>(NumSynthesisPortModes());
    long n_wave = 0;
    for (long j = n_lumped; j < static_cast<long>(v_node_label.size()) &&
                            v_node_label[j].rfind("waveport_", 0) == 0;
         j++)
    {
      n_wave++;
    }
    for (long j = 0; j < n_lumped + n_wave; j++)
    {
      d(j) = orth_R(j, j);
    }
    const double unit_ohm_inv = 1.0 / units.GetScaleFactor<Units::ValueType::IMPEDANCE>();
    const double unit_henry_inv =
        1.0 / units.GetScaleFactor<Units::ValueType::INDUCTANCE>();
    const double omega0 = unit_henry_inv / unit_ohm_inv;
    coupled_s.resize(sweep_omega_samples.size());
    for (std::size_t fi = 0; fi < sweep_omega_samples.size(); fi++)
    {
      const double omega = sweep_omega_samples[fi];
      std::vector<Eigen::VectorXcd> sv;
      for (int port_idx : coupled_ports)
      {
        auto s = space_op.GetWavePortModeVector(port_idx, omega);
        Eigen::VectorXcd p(nr);
        ProjectVecInternal(space_op.GetComm(), V, *s, p, 0);
        sv.push_back(std::move(p));
      }
      const std::complex<double> s_phys = 1i * omega * omega0;
      Eigen::MatrixXcd Y = (*matrices.L_inv) / s_phys + s_phys * (*matrices.C);
      if (matrices.R_inv)
      {
        Y += *matrices.R_inv;
      }
      auto &vals = coupled_s[fi];
      vals.resize(coupled_ports.size() * coupled_ports.size());
      for (std::size_t drive = 0; drive < coupled_ports.size(); drive++)
      {
        Eigen::VectorXcd rhs = Eigen::VectorXcd::Zero(na);
        rhs.head(nr) = -2.0 * unit_ohm_inv *
                       d.head(nr).cast<std::complex<double>>().cwiseProduct(sv[drive]);
        const Eigen::VectorXcd y = Y.fullPivLu().solve(rhs);
        for (std::size_t obs = 0; obs < coupled_ports.size(); obs++)
        {
          vals[obs * coupled_ports.size() + drive] =
              -(sv[obs].adjoint() *
                d.head(nr).cast<std::complex<double>>().cwiseProduct(y.head(nr)))(0, 0) -
              ((obs == drive) ? 1.0 : 0.0);
        }
      }
    }
  }

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
  if (!coupled_s.empty())
  {
    auto out = TableWithCSVFile(post_dir / "rom-coupled-S.csv");
    out.table.col_options.float_precision = 17;
    out.table.insert("f", "f (GHz)");
    for (std::size_t obs = 0; obs < coupled_ports.size(); obs++)
    {
      for (std::size_t drive = 0; drive < coupled_ports.size(); drive++)
      {
        const auto key = fmt::format("S[{:d}][{:d}]", obs + 1, drive + 1);
        out.table.insert(fmt::format("re_{}", key), fmt::format("Re{{{}}}", key));
        out.table.insert(fmt::format("im_{}", key), fmt::format("Im{{{}}}", key));
      }
    }
    const double unit_GHz =
        units.Dimensionalize<Units::ValueType::FREQUENCY>(1.0) / (2.0 * M_PI);
    for (std::size_t fi = 0; fi < sweep_omega_samples.size(); fi++)
    {
      out.table["f"] << sweep_omega_samples[fi] * unit_GHz;
      for (std::size_t obs = 0; obs < coupled_ports.size(); obs++)
      {
        for (std::size_t drive = 0; drive < coupled_ports.size(); drive++)
        {
          const auto key = fmt::format("S[{:d}][{:d}]", obs + 1, drive + 1);
          const auto value = coupled_s[fi][obs * coupled_ports.size() + drive];
          out.table[fmt::format("re_{}", key)] << value.real();
          out.table[fmt::format("im_{}", key)] << value.imag();
        }
      }
    }
    out.WriteFullTableTrunc();
  }
  auto print_table = [post_dir](const Eigen::MatrixXd &mat, std::string_view filename,
                                const std::vector<std::string> &table_labels)
  {
    MFEM_VERIFY((table_labels.size() == mat.cols()) && (table_labels.size() == mat.rows()),
                "Inconsistent PROM size!");

    auto out = TableWithCSVFile(post_dir / filename);
    out.table.col_options.float_precision = 17;
    for (long i = 0; i < mat.cols(); i++)
    {
      out.table.insert(table_labels[i], table_labels[i]);
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
  auto print_if_nonzero = [&](const Eigen::MatrixXd &mat, std::string_view filename,
                              const std::vector<std::string> &table_labels)
  {
    if (mat.cwiseAbs().maxCoeff() > 0.0)
    {
      print_table(mat, filename, table_labels);
    }
  };
  print_if_nonzero(inductance_L_inv->real(), "rom-Linv-re.csv", labels);
  print_if_nonzero(inductance_L_inv->imag(), "rom-Linv-im.csv", labels);
  print_if_nonzero(capacitance_C->real(), "rom-C-re.csv", labels);
  print_if_nonzero(capacitance_C->imag(), "rom-C-im.csv", labels);
  if (resistance_R_inv)
  {
    print_if_nonzero(resistance_R_inv->real(), "rom-Rinv-re.csv", labels);
    print_if_nonzero(resistance_R_inv->imag(), "rom-Rinv-im.csv", labels);
  }

  // Port-load matrices decompose the legacy loaded rom-* matrices into selectable
  // per-port terminations. Each file has the same row/column labels and dimension as the
  // total matrices, so downstream tools can form device matrices by subtracting any subset
  // of these loads from rom-* and can add back only the external loads after connecting
  // internal ports.
  for (const auto &load : matrices.port_loads)
  {
    const auto prefix = fmt::format("rom-portload-{}", load.label);
    print_if_nonzero(load.L_inv->real(), fmt::format("{}-Linv-re.csv", prefix), labels);
    print_if_nonzero(load.L_inv->imag(), fmt::format("{}-Linv-im.csv", prefix), labels);
    print_if_nonzero(load.C->real(), fmt::format("{}-C-re.csv", prefix), labels);
    print_if_nonzero(load.C->imag(), fmt::format("{}-C-im.csv", prefix), labels);
    if (load.R_inv)
    {
      print_if_nonzero(load.R_inv->real(), fmt::format("{}-Rinv-re.csv", prefix), labels);
      print_if_nonzero(load.R_inv->imag(), fmt::format("{}-Rinv-im.csv", prefix), labels);
    }
  }

  // Print orth-R. Don't divide by diagonal to keep state normalization info.
  // Pad with identity for aux states (regime 2): aux rows are synthetic circuit
  // nodes with no basis-orthogonality content; identity preserves the diagonal
  // form expected by downstream consumers.
  Eigen::MatrixXd orth_R_padded = Eigen::MatrixXd::Identity(labels.size(), labels.size());
  orth_R_padded.topLeftCorner(orth_R.rows(), orth_R.cols()) = orth_R;
  print_table(orth_R_padded, "rom-orthogonalization-matrix-R.csv", labels);

  PrintPortReferenceData(units, post_dir, matrices);

  // Write the synthesized-system eigenvalue estimates (computed above, before the
  // root-only guard) together with their HDM eigenpair errors.
  if (!eigs.empty())
  {
    auto out = TableWithCSVFile(post_dir / "rom-eigenvalues.csv");
    out.table.col_options.float_precision = 12;
    out.table.reserve(eigs.size(), 5);
    out.table.insert("re_f", "Re{f} (GHz)");
    out.table.insert("im_f", "Im{f} (GHz)");
    out.table.insert("Q", "Q");
    out.table.insert("err_bkwd", "Error (Bkwd.)");
    out.table.insert("err_abs", "Error (Abs.)");
    for (const auto &e : eigs)
    {
      out.table["re_f"] << e.freq_re_GHz;
      out.table["im_f"] << e.freq_im_GHz;
      out.table["Q"] << e.Q;
      out.table["err_bkwd"] << e.error_bkwd;
      out.table["err_abs"] << e.error_abs;
    }
    out.WriteFullTableTrunc();

    // Write the matching eigenvectors in the synthesized node space: one column per node
    // (basis + aux, same labels and order as the rom-Linv/Rinv/C matrix columns), one row
    // per mode (row k is the eigenvector of row k in rom-eigenvalues.csv). Unit 2-norm,
    // phase fixed so the largest-magnitude entry is real-positive. A consumer can check
    // (L⁻¹ + iωR⁻¹ − ω²C)·y ≈ 0 directly from the printed artifacts.
    auto print_eigvecs = [&](auto accessor, std::string_view filename)
    {
      auto vec_out = TableWithCSVFile(post_dir / filename);
      vec_out.table.col_options.float_precision = 17;
      vec_out.table.reserve(eigs.size(), labels.size());
      for (std::size_t i = 0; i < labels.size(); i++)
      {
        vec_out.table.insert(labels[i], labels[i]);
        auto &col = vec_out.table[static_cast<long>(i)];
        for (const auto &e : eigs)
        {
          // The estimates are computed on the augmented pencil, so the eigenvector
          // dimension always equals the label count (basis + aux).
          MFEM_VERIFY(e.eigvec.size() == static_cast<Eigen::Index>(labels.size()),
                      "Synthesized eigenvector dimension does not match the node count!");
          col << accessor(e.eigvec(static_cast<Eigen::Index>(i)));
        }
      }
      vec_out.WriteFullTableTrunc();
    };
    print_eigvecs([](std::complex<double> v) { return v.real(); },
                  "rom-eigenvectors-re.csv");
    print_eigvecs([](std::complex<double> v) { return v.imag(); },
                  "rom-eigenvectors-im.csv");
    Mpi::Print("\n Synthesized-system eigenvalue estimates ({:d} modes in [{:.3f}, "
               "{:.3f}] GHz):\n",
               eigs.size(), fmin_GHz, fmax_GHz);
    for (const auto &e : eigs)
    {
      Mpi::Print("   f = {:+.6e} {:+.6e}i GHz,  Q = {:.3e},  bkwd. error = {:.3e}\n",
                 e.freq_re_GHz, e.freq_im_GHz, e.Q, e.error_bkwd);
    }
  }
}

std::vector<RomOperator::EigenvalueEstimate> RomOperator::ComputeEigenvalueEstimates(
    const Eigen::MatrixXcd &L_inv, const Eigen::MatrixXcd *R_inv, const Eigen::MatrixXcd &C,
    double fmin_GHz, double fmax_GHz, int n_basis)
{
  // Solve the quadratic eigenvalue problem (L⁻¹ + iωR⁻¹ − ω²C)v = 0 via companion
  // linearization. SI matrices span ~28 orders of magnitude (L⁻¹ ~ 1e14, C ~ 1e-14), so
  // we nondimensionalize first: substitute ω = ω₀·ω' with ω₀ = √(‖K‖/‖C‖), then divide
  // through by ‖K‖, making both leading blocks O(1). Results are converted back to physical
  // units.
  const int n = static_cast<int>(L_inv.rows());
  MFEM_VERIFY(n > 0 && C.rows() == n, "Empty synthesis matrices in eigenvalue estimate!");

  const double norm_K = L_inv.cwiseAbs().maxCoeff();
  const double norm_C = C.cwiseAbs().maxCoeff();
  const double w0 = (norm_K > 0.0 && norm_C > 0.0) ? std::sqrt(norm_K / norm_C) : 1.0;
  const double scale = (norm_K > 0.0) ? norm_K : 1.0;

  Eigen::MatrixXcd Ks = L_inv / scale;
  Eigen::MatrixXcd Cs = (w0 * w0 / scale) * C;
  Eigen::MatrixXcd Gs = Eigen::MatrixXcd::Zero(n, n);
  if (R_inv)
  {
    Gs = (w0 / scale) * (*R_inv);
  }

  // First-companion linearization in s = iω: substituting s into the physical pencil
  // (L⁻¹ + iωR⁻¹ − ω²C) v = 0 gives P(s) = Ks + s·Gs + s²·Cs (all coefficient signs
  // positive), so with x = [v; s·v]:
  //   [0 I; -Ks -Gs] x = s [I 0; 0 Cs] x.
  // The eigenvalue is s = iω, converted back below via ω = s / i.
  Eigen::MatrixXcd A = Eigen::MatrixXcd::Zero(2 * n, 2 * n);
  Eigen::MatrixXcd B = Eigen::MatrixXcd::Zero(2 * n, 2 * n);
  A.topRightCorner(n, n) = Eigen::MatrixXcd::Identity(n, n);
  A.bottomLeftCorner(n, n) = -Ks;
  A.bottomRightCorner(n, n) = -Gs;
  B.topLeftCorner(n, n) = Eigen::MatrixXcd::Identity(n, n);
  B.bottomRightCorner(n, n) = Cs;

  // Dense generalized eigenvalue solve via LAPACK zggev (complex QZ). B may be singular
  // (augmented-state rows with zero capacitance → infinite eigenvalues), which a B⁻¹A
  // standard-EVP approach cannot handle; the QZ factorization returns inf for those
  // and finite values for physical modes.
  int N2 = 2 * n;
  Eigen::MatrixXcd A_col = A;  // zggev overwrites; column-major by default
  Eigen::MatrixXcd B_col = B;
  Eigen::VectorXcd alpha(N2), beta(N2);
  Eigen::MatrixXcd vl_dummy(1, 1), vr(N2, N2);
  int ldvl = 1, ldvr = N2;
  int lwork = 4 * N2;
  Eigen::VectorXcd work(lwork);
  Eigen::VectorXd rwork(8 * N2);
  int info = 0;
  char jobN = 'N', jobV = 'V';
  zggev_(&jobN, &jobV, &N2, A_col.data(), &N2, B_col.data(), &N2, alpha.data(), beta.data(),
         vl_dummy.data(), &ldvl, vr.data(), &ldvr, work.data(), &lwork, rwork.data(),
         &info);
  MFEM_VERIFY(info == 0, "zggev failed with info = " << info);

  // Eigenvalues are alpha/beta; infinite when |beta| ≈ 0.
  Eigen::VectorXcd s(N2);
  for (int k = 0; k < N2; k++)
  {
    s(k) = (std::abs(beta(k)) > 1.0e-300) ? alpha(k) / beta(k)
                                          : std::complex<double>(1.0e300, 0.0);
  }

  // Convert eigenvalues back to physical frequencies (GHz): the companion eigenvalue is
  // s = iω' (nondimensional), so ω = w0 · s / i, then f = ω / (2π·1e9). The pencil is
  // complex symmetric (not Hermitian), so roots do not come in conjugate pairs: each
  // physical mode contributes one root whose Im{f} is positive for decay, matching the
  // eigenmode solver's eig.csv convention directly. Filter to the trained band and
  // Q > 0.5: the augmented realization's aux states (zero capacitance rows,
  // cond(C) = ∞) produce spurious near-critically-damped roots at Q ≲ 0.5 that carry no
  // physical content..
  const std::complex<double> inv_i(0.0, -1.0);  // 1/i = -i
  std::vector<EigenvalueEstimate> modes;
  for (int k = 0; k < s.size(); k++)
  {
    if (!std::isfinite(s(k).real()) || !std::isfinite(s(k).imag()))
    {
      continue;
    }
    const std::complex<double> omega_phys = w0 * (s(k) * inv_i);
    const double f_re = omega_phys.real() / (2.0 * M_PI * 1.0e9);
    const double f_im = omega_phys.imag() / (2.0 * M_PI * 1.0e9);
    if (f_re < fmin_GHz || f_re > fmax_GHz)
    {
      continue;
    }
    // Match the eigenmode postprocessor's convention Q = |ω| / (2 |Im ω|) (full complex
    // magnitude in the numerator, not |Re ω|) so synthesized and eigensolver Q agree.
    const double abs_omega = std::abs(omega_phys);
    const double abs_omega_im = std::abs(omega_phys.imag());
    const double Q = (abs_omega_im > 1.0e-20) ? abs_omega / (2.0 * abs_omega_im)
                                              : std::numeric_limits<double>::infinity();
    if (Q <= SYNTHESIS_EIG_Q_MIN)
    {
      continue;
    }
    EigenvalueEstimate est;
    est.freq_re_GHz = f_re;
    est.freq_im_GHz = f_im;
    est.Q = Q;
    // The companion eigenvector is x = [v; s·v]; keep the "v" block (the synthesized
    // pencil's node coordinates: basis + aux rows), normalized to unit 2-norm with the
    // largest-magnitude entry rotated to be real and positive. The LAPACK phase is
    // arbitrary (and MPI-partition dependent), so fixing it makes the printed
    // eigenvectors reproducible across runs.
    est.eigvec = vr.col(k).head(n);
    const double vnorm = est.eigvec.norm();
    if (vnorm > 0.0)
    {
      est.eigvec /= vnorm;
      // Drop spurious aux-pole roots: the augmented realization's aux states inject roots
      // at their pole frequencies whose eigenvector lives almost entirely in the trailing
      // aux rows, whereas a physical resonance has substantial energy in the basis rows.
      if (n_basis > 0 && n_basis < n)
      {
        const double basis_frac = est.eigvec.head(n_basis).norm();  // eigvec is unit-norm
        if (basis_frac < SYNTHESIS_EIG_BASIS_FRAC_MIN)
        {
          continue;
        }
      }
      Eigen::Index i_max;
      est.eigvec.cwiseAbs().maxCoeff(&i_max);
      const std::complex<double> pivot = est.eigvec(i_max);
      if (std::abs(pivot) > 0.0)
      {
        est.eigvec *= std::abs(pivot) / pivot;
      }
    }
    modes.push_back(std::move(est));
  }
  std::sort(modes.begin(), modes.end(),
            [](const auto &a, const auto &b) { return a.freq_re_GHz < b.freq_re_GHz; });

  // No frequency-based deduplication: the complex eigenfrequency is not a unique mode
  // identifier (a genuinely degenerate eigenspace yields one QZ eigenvalue per
  // multiplicity, each with its own eigenvector, and all must be reported). The pencil is
  // complex-symmetric, so roots do not come in conjugate pairs and there is no redundant
  // numerical representation to remove here; if one is ever observed it needs an
  // eigenvector/rank-aware policy, not frequency clustering.
  return modes;
}

void RomOperator::ComputeEigenvalueEstimateErrors(
    const Units &units, std::vector<EigenvalueEstimate> &estimates) const
{
  // Evaluate each synthesized eigenpair on the HDM, following the eigensolver residual
  // conventions (cf. QuasiNewtonSolver::GetResidualNorm / GetBackwardScaling):
  //   error_abs  = ‖(K + iωC − ω²M + A2(ω)) u‖₂ / ‖u‖₂,
  //   error_bkwd = error_abs / (‖K‖₂ + |ω|‖C‖₂ + |ω|²‖M‖₂),
  // with ω the complex eigenfrequency in nondimensional angular units and u = V·y the
  // prolongation of the eigenvector's basis block. The frequency-dependent A2(ω) —
  // wave ports, second-order farfield ABC, surface conductivity — is evaluated exactly
  // at complex ω via SpaceOperator::GetExtraSystemMatrix, so this also validates the
  // rational dispersion fit content that the aux rows encode. The aux entries of the
  // eigenvector have no HDM image and are dropped from the prolongation; for a converged
  // physical mode their rational contribution is reproduced by the true A2(ω), while for
  // a spurious aux-dominated root the truncated u cannot satisfy the HDM equation and the
  // backward error is O(1). The wave-port modal correction W is included via
  // GetExtraSystemOperator (not GetExtraSystemMatrix) so the residual accounts for it.
  // Collective on the space communicator.
  if (estimates.empty() || V.empty())
  {
    return;
  }
  MPI_Comm comm = space_op.GetComm();
  const std::size_t n_basis = V.size();

  // The synthesized matrices are the reduced pencil conjugated by the diagonal port
  // scaling v_d (CalculateNormalizedPROMMatrices: L⁻¹ = h⁻¹·v_d·Kr·v_d etc., with
  // v_conc(j) = orth_R(j,j) on the port rows and 1 elsewhere). An eigenvector y of the
  // synthesized pencil therefore corresponds to reduced-basis coordinates
  // y_reduced = v_conc ∘ y.
  Eigen::VectorXd v_conc = Eigen::VectorXd::Ones(static_cast<long>(n_basis));
  {
    const long n_port_modes = static_cast<long>(NumSynthesisPortModes());
    long n_waveport_rows = 0;
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
    for (long j = 0; j < n_port_modes + n_waveport_rows; j++)
    {
      v_conc[j] = orth_R(j, j);
    }
  }

  // Spectral norms of the ω-independent operators (SLEPc/ARPACK use these for backward
  // scaling too). Computed once.
  const double normK = linalg::SpectralNorm(comm, *K, K->IsReal());
  const double normC = C ? linalg::SpectralNorm(comm, *C, C->IsReal()) : 0.0;
  const double normM = linalg::SpectralNorm(comm, *M, M->IsReal());

  ComplexVector u(K->Width()), res(K->Width());
  u.UseDevice(true);
  res.UseDevice(true);
  const double freq_to_omega_nd =
      2.0 * M_PI * units.Nondimensionalize<Units::ValueType::FREQUENCY>(1.0);
  for (auto &est : estimates)
  {
    // Physical complex frequency f (GHz) → nondimensional complex angular frequency.
    // Palace's driven/eigenmode time convention is e^{iωt}: a decaying mode (positive
    // Im{f} in eig.csv) has eigenfrequency ω = 2π(Re{f} + i·Im{f})·tc.
    const std::complex<double> omega(freq_to_omega_nd * est.freq_re_GHz,
                                     freq_to_omega_nd * est.freq_im_GHz);

    // Prolongate the basis block: u = Σ_j (v_conc_j y_j) V[j] (V real, y complex; the
    // v_conc factor undoes the port-row scaling baked into the synthesized pencil).
    u = 0.0;
    for (std::size_t j = 0; j < n_basis && j < static_cast<std::size_t>(est.eigvec.size());
         j++)
    {
      const std::complex<double> yj = v_conc[static_cast<long>(j)] * est.eigvec(j);
      linalg::AXPY(yj.real(), V[j], u.Real());
      linalg::AXPY(yj.imag(), V[j], u.Imag());
    }
    const double unorm = linalg::Norml2(comm, u);
    if (unorm == 0.0)
    {
      // Pure aux-state root: no basis content at all; flag with an O(1) sentinel since
      // the mode has no HDM representation.
      est.error_abs = est.error_bkwd = 1.0;
      continue;
    }

    // res = (K + iωC − ω²M) u.
    K->Mult(u, res);
    if (C)
    {
      C->AddMult(u, res, 1i * omega);
    }
    M->AddMult(u, res, -omega * omega);
    // The frequency-dependent boundary terms at complex ω (null when absent), including
    // the wave-port modal correction W (GetExtraSystemOperator, not GetExtraSystemMatrix).
    auto A2_omega = space_op.GetExtraSystemOperator(omega, Operator::DIAG_ZERO);
    if (A2_omega)
    {
      A2_omega->AddMult(u, res, 1.0);
    }

    est.error_abs = linalg::Norml2(comm, res) / unorm;
    const double t = std::abs(omega);
    est.error_bkwd = est.error_abs / (normK + t * normC + t * t * normM);
  }
}

}  // namespace palace
