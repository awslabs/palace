// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "romoperator.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <limits>
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

namespace palace
{

using namespace std::complex_literals;

namespace
{

constexpr auto ORTHOG_TOL = 1.0e-12;

// ---------------------------------------------------------------------------------------
// DtN structured-√ wave-port synthesis helpers. The target is the boundary multiplier
// q(s) = i·kₙ(ω) (s = iω), modeled as q(s) = √c₂·u·P(t)/Q(t) with u = s + shift,
// t = ω_s²/u² (so 1−t = (u²−ω_s²)/u² and q² = c₂·u²(1−t) = c₂(u²−ω_s²) = c₀+c₁s+c₂s²),
// then converted to an s-domain pole-residue form q(s) = poly0 + poly1·s + Σ R_k/(s−p_k).
// See prom-waveport-validation/dtn_realization_lock.py for the validated reference; all
// expressions below mirror it (verified to 1e-15 against the structured fit and the pencil
// realization).

// Roots of a polynomial given by ASCENDING complex coefficients (c[0] + c[1] x + ...), via
// the companion-matrix eigenvalues. Leading near-zero coefficients are trimmed.
inline std::vector<std::complex<double>> PolyRootsAscending(std::vector<std::complex<double>> c)
{
  while (c.size() > 1 && std::abs(c.back()) <= 1.0e-300)
  {
    c.pop_back();
  }
  const long n = static_cast<long>(c.size()) - 1;  // degree
  if (n <= 0)
  {
    return {};
  }
  // Monic companion matrix (eigenvalues = roots).
  Eigen::MatrixXcd C = Eigen::MatrixXcd::Zero(n, n);
  for (long i = 1; i < n; i++)
  {
    C(i, i - 1) = 1.0;
  }
  for (long i = 0; i < n; i++)
  {
    C(i, n - 1) = -c[i] / c[n];
  }
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(C, /*computeEigenvectors=*/false);
  std::vector<std::complex<double>> roots(es.eigenvalues().data(),
                                          es.eigenvalues().data() + n);
  return roots;
}

// Evaluate a polynomial (ascending complex coefficients) at z (Horner from the top).
inline std::complex<double> PolyEvalAscending(const std::vector<std::complex<double>> &c,
                                              std::complex<double> z)
{
  std::complex<double> y = 0.0;
  for (long k = static_cast<long>(c.size()) - 1; k >= 0; k--)
  {
    y = y * z + c[k];
  }
  return y;
}

// Convert the structured fit q(s)=√c₂·u·P(t)/Q(t), t=ω_s²/u², u=s+shift into an s-domain
// pole-residue form q(s)=poly0+poly1·s + Σ R_k/(s−p_k). Returns the poly part (ascending,
// expected degree ≤ 1) and the pole/residue lists. Mirrors to_s_partial_fractions() in the
// validated Python (verified pf_recon_err ~1e-15).
inline void DtnPartialFractions(std::complex<double> sqrt_c2, double shift, double omega_s2,
                                const Eigen::VectorXd &a, const Eigen::VectorXd &b,
                                std::vector<std::complex<double>> &poly_asc,
                                std::vector<std::complex<double>> &poles,
                                std::vector<std::complex<double>> &residues)
{
  const long deg = b.size() - 1;
  // Build numerator/denominator as polynomials in u:
  //   num_u = √c₂ · Σ_k a_k·ω_s²^k · u^(2(deg−k)+1)   (the leading u factor),
  //   den_u =       Σ_k b_k·ω_s²^k · u^(2(deg−k)).
  std::vector<std::complex<double>> num_u(2 * deg + 2, 0.0), den_u(2 * deg + 1, 0.0);
  double os2_pow = 1.0;
  std::vector<double> os2_pows(std::max(a.size(), b.size()));
  for (long k = 0; k < static_cast<long>(os2_pows.size()); k++)
  {
    os2_pows[k] = os2_pow;
    os2_pow *= omega_s2;
  }
  for (long k = 0; k < a.size(); k++)
  {
    num_u[2 * (deg - k) + 1] += sqrt_c2 * a(k) * os2_pows[k];
  }
  for (long k = 0; k < b.size(); k++)
  {
    den_u[2 * (deg - k)] += b(k) * os2_pows[k];
  }
  // Substitute u = s + shift: expand each u^p via the binomial into ascending s-powers.
  auto shift_poly = [shift](const std::vector<std::complex<double>> &cu)
  {
    std::vector<std::complex<double>> out(cu.size(), 0.0);
    for (std::size_t pu = 0; pu < cu.size(); pu++)
    {
      if (std::abs(cu[pu]) == 0.0)
      {
        continue;
      }
      // (s + shift)^pu = Σ_j C(pu,j)·shift^(pu−j)·s^j. Accumulate shift^(pu−j) from j=pu
      // downward (shift^0=1 at j=pu, then multiply by `shift` as j decreases) so the leading
      // s^pu term is always present even when shift==0 (the lossless case shift=c₁/(2c₂)≈0,
      // where (s+0)^pu = s^pu and only j=pu survives). Binomials C(pu,j) built the same way.
      double binom = 1.0;  // C(pu, pu)
      double shift_pow = 1.0;  // shift^(pu-j) at j=pu
      for (long j = static_cast<long>(pu); j >= 0; j--)
      {
        out[j] += cu[pu] * binom * shift_pow;
        // Step j -> j-1: C(pu,j-1) = C(pu,j)·j/(pu-j+1); shift^(pu-j) -> shift^(pu-j+1).
        binom *= static_cast<double>(j) / static_cast<double>(pu - j + 1);
        shift_pow *= shift;
      }
    }
    return out;
  };
  std::vector<std::complex<double>> num_s = shift_poly(num_u);
  std::vector<std::complex<double>> den_s = shift_poly(den_u);
  poles = PolyRootsAscending(den_s);
  // Proper part by polynomial long division num_s / den_s (ascending). Compute quotient and
  // remainder. Quotient degree = deg(num) − deg(den).
  // Work in descending order for stable synthetic division.
  std::vector<std::complex<double>> num_d(num_s.rbegin(), num_s.rend());
  std::vector<std::complex<double>> den_d(den_s.rbegin(), den_s.rend());
  while (num_d.size() > 1 && std::abs(num_d.front()) <= 1.0e-300)
  {
    num_d.erase(num_d.begin());
  }
  while (den_d.size() > 1 && std::abs(den_d.front()) <= 1.0e-300)
  {
    den_d.erase(den_d.begin());
  }
  std::vector<std::complex<double>> quo_d, rem_d(num_d);
  if (num_d.size() >= den_d.size())
  {
    const std::size_t nq = num_d.size() - den_d.size() + 1;
    quo_d.assign(nq, 0.0);
    for (std::size_t i = 0; i < nq; i++)
    {
      std::complex<double> coef = rem_d[i] / den_d[0];
      quo_d[i] = coef;
      for (std::size_t j = 0; j < den_d.size(); j++)
      {
        rem_d[i + j] -= coef * den_d[j];
      }
    }
  }
  // poly part (ascending).
  poly_asc.assign(quo_d.rbegin(), quo_d.rend());
  // Residues R_k = rem(p_k) / den'(p_k), with rem and den' as descending polynomials.
  std::vector<std::complex<double>> rem_asc(rem_d.rbegin(), rem_d.rend());
  std::vector<std::complex<double>> den_asc(den_d.rbegin(), den_d.rend());
  std::vector<std::complex<double>> dprime_asc(
      den_asc.size() > 1 ? den_asc.size() - 1 : 0, 0.0);
  for (long k = 1; k < static_cast<long>(den_asc.size()); k++)
  {
    dprime_asc[k - 1] = den_asc[k] * static_cast<double>(k);
  }
  residues.clear();
  residues.reserve(poles.size());
  for (const auto &p : poles)
  {
    residues.push_back(PolyEvalAscending(rem_asc, p) / PolyEvalAscending(dprime_asc, p));
  }
}

// Online evaluation of the DtN pole-residue model: q(iω) = poly0 + poly1·(iω) + Σ R_k/(iω−p_k).
inline std::complex<double> EvalDtnPoleResidue(std::complex<double> poly0,
                                               std::complex<double> poly1, double omega,
                                               const std::vector<std::complex<double>> &poles,
                                               const std::vector<std::complex<double>> &residues)
{
  const std::complex<double> s(0.0, omega);
  std::complex<double> q = poly0 + poly1 * s;
  for (std::size_t k = 0; k < poles.size(); k++)
  {
    q += residues[k] / (s - poles[k]);
  }
  return q;
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

  // Cache the ω-independent boundary masses for the other frequency-dependent BCs so they
  // can be folded into circuit synthesis (projected per basis growth in UpdatePROM, fit
  // in CalculateNormalizedPROMMatrices). On the imaginary slot, matching Mwp_p.
  //   - 2nd-order farfield ABC: M_ff with full term i·(0.5/ω)·M_ff.
  //   - Surface conductivity: one boundary mass per active attribute group, term
  //     (i·ω/Z_g(ω))·A_σ_g.
  M_ff_ = space_op.GetFarfieldExtraBoundaryMatrix<ComplexOperator>(Operator::DIAG_ZERO,
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

  // Capture sweep band and synthesis tolerance for later use in
  // CalculateNormalizedPROMMatrices (polynomial fit residual + bound). After
  // config::Nondimensionalize (utils/configfile.cpp), sample_f stores
  // 2π·f_nondim — i.e. ω in nondimensional units, not the linear frequency. Use it
  // directly. Default the synthesis tol to adaptive_tol when the user did not set one
  // explicitly.
  const auto &sample_f = iodata.solver.driven.sample_f;
  if (!sample_f.empty())
  {
    sweep_omega_samples = sample_f;
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
  auto P = space_op.GetPreconditionerMatrix<ComplexOperator>(
      1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i, omega + 0.0i);
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
  // Other frequency-dependent BC boundary masses folded into circuit synthesis: the
  // 2nd-order farfield ABC (M_ff_) and each surface-conductivity group (Asig_g_). Projected
  // like the wave-port masses so the synthesis path treats them uniformly.
  if (M_ff_)
  {
    M_ff_r.conservativeResize(dim_V_new, dim_V_new);
    ProjectMatInternal(comm, V, *M_ff_, M_ff_r, r, dim_V_old);
  }
  for (std::size_t g = 0; g < Asig_g_.size(); g++)
  {
    if (Asig_g_[g])
    {
      Asig_g_r[g].conservativeResize(dim_V_new, dim_V_new);
      ProjectMatInternal(comm, V, *Asig_g_[g], Asig_g_r[g], r, dim_V_old);
    }
  }
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
  if (V.empty())
  {
    u = 0.0;
    return;
  }

  Ar.resize(V.size(), V.size());
  RHSr.resize(V.size());

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
  // (a future ω-dependent non-wave-port BC added to GetExtraSystemMatrix would have no
  // factored operator and would be silently dropped). So on the first factored online solve
  // we additionally verify the factored Aᵣ contribution against the full HDM projection; on
  // any mismatch we latch other_A2_factored_ok = false and use the slow fallback for the
  // rest of the sweep.
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
      other_A2_self_checked = true;
      Eigen::MatrixXcd Ar_factored = Ar;
      Eigen::MatrixXcd Ar_hdm = Eigen::MatrixXcd::Zero(V.size(), V.size());
      A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO,
                                                          /*include_wave_ports=*/false);
      if (A2)
      {
        ProjectMatInternal(space_op.GetComm(), V, *A2, Ar_hdm, r, 0);
      }
      const double err = (Ar_factored - Ar_hdm).cwiseAbs().maxCoeff();
      const double ref = std::max(Ar_hdm.cwiseAbs().maxCoeff(), 1.0e-300);
      if (err / ref > 1.0e-9)
      {
        other_A2_factored_ok = false;
        Ar = Ar_hdm;  // Use the trusted HDM projection for this solve.
        Mpi::Warning("Factored online A2 (farfield ABC / surface conductivity) disagrees "
                     "with the full operator (rel. err {:.3e})!\n"
                     "Reverting to the per-frequency assembled A2 for the remaining sweep. "
                     "This indicates an ω-dependent boundary condition not covered by the "
                     "factored path.\n",
                     err / ref);
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
      ProjectMatInternal(space_op.GetComm(), V, *A2, Ar, r, 0);
    }
  }
  Ar += Kr;
  if (C)
  {
    Ar += (1i * omega) * Cr;
  }
  Ar += (-omega * omega) * Mr;
  // Wave-port contribution: A_wp(ω) = i·Σ_p kₙ,p(ω)·M^(p)_{μ⁻¹}. The Mwp_p HDM operators
  // are purely imaginary (boundary mass on the imaginary slot), so multiplying their
  // projection by the real scalar kₙ recovers the i·kₙ·M contribution. kₙ,p(ω) comes from
  // GetWavePortKn, which re-solves the per-port cross-section EVP at this ω. This is the
  // dominant ω-dependent online work (timed under WAVE_PORT), but it is ALSO REQUIRED for
  // correctness: GetWavePortKn → WavePortData::Initialize(ω) refreshes the modal
  // post-processing state (E0t/E0n, port_sr/port_si, Z_PV) that MeasureWavePorts uses to
  // project the field onto the port mode for S-parameters/power. An online kₙ SURROGATE
  // (the WavePortDispersion fit) would skip this refresh, leaving stale modal fields and
  // producing wrong S-parameters for dispersive non-TEM ports near cutoff (verified: |ΔS|
  // 0.25 with surrogate vs 8e-4 with the per-ω EVP on the adapter). The dispersion fit is
  // therefore used ONLY for the offline L⁻¹/R⁻¹/C circuit synthesis (where no modal field
  // post-processing is needed), not for the online sweep.
  {
    BlockTimer bt(Timer::WAVE_PORT);
    for (const auto &[port_idx, Mp_r] : Mwp_p_r)
    {
      const double kn = space_op.GetWavePortOp().GetWavePortKn(port_idx, omega);
      Ar += std::complex<double>(kn, 0.0) * Mp_r;
    }
  }

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
  // QR solve, for maximal stability. The small system is cheap to compute but can be
  // numerically poorly conditioned to due the splitting of HDM solutions into Re and Im
  // into separate columns.
  RHSr = Ar.fullPivHouseholderQr().solve(RHSr);
  ProlongatePROMSolution(V.size(), V, RHSr, u);
}

RomOperator::WavePortRegime RomOperator::SelectWavePortRegime(int port_idx, double rel_err,
                                                              bool meets_tol) const
{
  // AUTO: polynomial if residual meets tolerance, else the structured-√ DtN rational fit
  // (handles lossless and lossy non-TEM ports uniformly).
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
    case WavePortSynthesisRegime::DTN_RATIONAL:
      return WavePortRegime::DtnRational;
    case WavePortSynthesisRegime::AUTO:
    default:
      return meets_tol ? WavePortRegime::Polynomial : WavePortRegime::DtnRational;
  }
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
  // Both grids are Chebyshev–Gauss–Lobatto nodes on [w_lo, w_hi]: they cluster at the band
  // edges where polynomial/rational interpolation error peaks, giving a tighter max-error
  // bound per sample than a uniform grid. The dense (validation) grid is offset to the
  // Chebyshev interior nodes so its points are DISTINCT from the fit grid — the LSQ
  // residual is artificially small AT the fit points, so the residual must be measured
  // between them.
  constexpr int n_fit = 12;
  constexpr int n_dense = 24;
  const double w_lo = sweep_omega_min;
  const double w_hi = sweep_omega_max;
  const double w_mid = 0.5 * (w_lo + w_hi);
  const double w_half = 0.5 * (w_hi - w_lo);
  // Chebyshev–Gauss–Lobatto: ω_i = mid − half·cos(π i/(n−1)), i = 0..n−1 (endpoints incl.).
  auto sample_cgl = [w_mid, w_half](int n)
  {
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
  };
  // Chebyshev–Gauss (interior) nodes: ω_j = mid − half·cos(π(2j+1)/(2n)), strictly inside
  // (w_lo, w_hi) and interlacing the CGL fit nodes, so validation points never coincide
  // with fit points.
  auto sample_cg = [w_mid, w_half](int n)
  {
    std::vector<double> ws(n);
    for (int j = 0; j < n; j++)
    {
      ws[j] = w_mid - w_half * std::cos(M_PI * (2 * j + 1) / (2.0 * n));
    }
    return ws;
  };
  auto fit_omegas = sample_cgl(n_fit);
  auto dense_omegas = sample_cg(n_dense);

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
    // fit.aux is still empty here, so EvaluateWavePortKnFit returns the bare polynomial.
    const double poly = EvaluateWavePortKnFit(fit, w);
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
          port_idx, fit.rel_err_polynomial, waveport_synthesis_tol, fit.alpha0, fit.alpha1,
          fit.alpha2);
    }
    return fit;
  }

  if (fit.regime == WavePortRegime::DtnRational)
  {
    FitWavePortDtnRational(fit, Mp_r, fit_omegas, y_fit, dense_omegas);
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

  std::vector<std::complex<double>> poles, residues;
  poles.reserve(pr.poles.size());
  residues.reserve(pr.residues.size());
  for (long k = 0; k < pr.poles.size(); k++)
  {
    poles.push_back(pr.poles(k));
    residues.push_back(pr.residues(k));
  }

  WavePortAuxBlock blk;
  blk.port_idx = port_idx;
  AddAuxBlockDirections(blk, Mp_r, waveport_synthesis_rank_tol);
  blk.poles = std::move(poles);
  blk.residues = std::move(residues);
  // Attach the aux block now so EvaluateWavePortKnFit picks up the AAA poles for the
  // augmented-model residual check below (and so the stored fit is complete for the online
  // surrogate).
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
  const std::size_t rank_used = fit.aux->sigmas.size();
  const std::size_t aux_per_port = n_poles * rank_used;
  Mpi::Print(" Wave port {:d}: augmented synthesis residual {:.3e} → {:.3e} "
             "(tol {:.3e}, {:d} pole{} × rank-{:d} mass = +{:d} aux state{}, "
             "α₀={:.3e}, α₁={:.3e}, α₂={:.3e}, d={:.3e})\n",
             port_idx, fit.rel_err_polynomial, fit.rel_err_augmented,
             waveport_synthesis_tol, n_poles, n_poles == 1 ? "" : "s", rank_used,
             aux_per_port, aux_per_port == 1 ? "" : "s", fit.alpha0, fit.alpha1, fit.alpha2,
             aaa_d);
  return fit;
}

void RomOperator::FitWavePortDtnRational(WavePortDispersionFit &fit,
                                         const Eigen::MatrixXcd &Mp_r,
                                         const std::vector<double> &fit_omegas,
                                         const Eigen::VectorXd &y_fit,
                                         const std::vector<double> &dense_omegas) const
{
  using cplx = std::complex<double>;
  const int port_idx = fit.port_idx;
  const int n_fit = static_cast<int>(fit_omegas.size());

  // 1) Fit q(s)² = c₀ + c₁s + c₂s² from q(iω)² = −kₙ². Real-ω samples give a real kₙ for a
  // lossless port (c₁ ≈ 0); for a lossy port kₙ has an imaginary part and c₁ ≠ 0. Here y_fit
  // is the real propagating kₙ (driven path), so the fit is effectively lossless (c₁ ≈ 0);
  // the structure still admits c₁ for a future complex-kₙ sampling without code changes.
  Eigen::MatrixXcd Aq(n_fit, 3);
  Eigen::VectorXcd bq(n_fit);
  for (int i = 0; i < n_fit; i++)
  {
    const cplx s(0.0, fit_omegas[i]);
    Aq(i, 0) = 1.0;
    Aq(i, 1) = s;
    Aq(i, 2) = s * s;
    bq(i) = -cplx(y_fit(i), 0.0) * y_fit(i);  // (i·kₙ)² = −kₙ²
  }
  Eigen::VectorXcd cq = Aq.colPivHouseholderQr().solve(bq);
  const cplx c0 = cq(0), c1 = cq(1), c2 = cq(2);

  // 2) Complete the square: q² = c₂(u² − ω_s²), u = s + shift.
  const cplx shift_c = c1 / (2.0 * c2);
  const cplx omega_s2_c = (c1 * c1 / (4.0 * c2) - c0) / c2;
  // Driven (real-kₙ) path: c₁, c₂ are real to round-off, so shift/ω_s² are real. Use the real
  // parts for the structured fit and realization (the imaginary residue is fit noise here).
  const double shift = shift_c.real();
  const double omega_s2 = omega_s2_c.real();
  const cplx sqrt_c2 = std::sqrt(c2);

  // 3) Fit √(1−t) ≈ P(t)/Q(t), t = ω_s²/u², real coefficients, column-scaled LSQ. The model
  // order reuses waveport_synthesis_order_max (P,Q degree).
  const int degree = std::max<int>(static_cast<int>(waveport_synthesis_order_max), 1);
  Eigen::MatrixXcd Asq(n_fit, 2 * degree + 1);
  Eigen::VectorXcd ysq(n_fit);
  for (int i = 0; i < n_fit; i++)
  {
    const cplx s(0.0, fit_omegas[i]);
    const cplx u = s + shift;
    const cplx t = omega_s2 / (u * u);
    const cplx y = std::sqrt(1.0 - t);
    ysq(i) = y;
    cplx tk = 1.0;
    for (int k = 0; k <= degree; k++)
    {
      Asq(i, k) = tk;  // P columns
      tk *= t;
    }
    tk = t;
    for (int k = 1; k <= degree; k++)
    {
      Asq(i, degree + k) = -y * tk;  // Q columns (b0 = 1 fixed)
      tk *= t;
    }
  }
  // Column scaling for conditioning, then solve the real-embedded LSQ.
  Eigen::VectorXd colscale(2 * degree + 1);
  for (int j = 0; j < 2 * degree + 1; j++)
  {
    double nrm = std::sqrt(Asq.col(j).real().squaredNorm() + Asq.col(j).imag().squaredNorm());
    colscale(j) = std::max(nrm, 1.0e-300);
  }
  Eigen::MatrixXcd Asq_s = Asq;
  for (int j = 0; j < 2 * degree + 1; j++)
  {
    Asq_s.col(j) /= colscale(j);
  }
  Eigen::MatrixXd Ar(2 * n_fit, 2 * degree + 1);
  Ar.topRows(n_fit) = Asq_s.real();
  Ar.bottomRows(n_fit) = Asq_s.imag();
  Eigen::VectorXd br(2 * n_fit);
  br.head(n_fit) = ysq.real();
  br.tail(n_fit) = ysq.imag();
  Eigen::VectorXd xs = Ar.colPivHouseholderQr().solve(br);
  for (int j = 0; j < 2 * degree + 1; j++)
  {
    xs(j) /= colscale(j);
  }
  Eigen::VectorXd a_coef(degree + 1), b_coef(degree + 1);
  for (int k = 0; k <= degree; k++)
  {
    a_coef(k) = xs(k);
  }
  b_coef(0) = 1.0;
  for (int k = 1; k <= degree; k++)
  {
    b_coef(k) = xs(degree + k);
  }

  // 4) Convert q(s) = √c₂·u·P(t)/Q(t) to the s-domain pole-residue form.
  std::vector<cplx> poly_asc, poles, residues;
  DtnPartialFractions(sqrt_c2, shift, omega_s2, a_coef, b_coef, poly_asc, poles, residues);
  fit.dtn_poly0 = poly_asc.size() > 0 ? poly_asc[0] : cplx(0.0, 0.0);
  fit.dtn_poly1 = poly_asc.size() > 1 ? poly_asc[1] : cplx(0.0, 0.0);

  WavePortAuxBlock blk;
  blk.port_idx = port_idx;
  AddAuxBlockDirections(blk, Mp_r, waveport_synthesis_rank_tol);
  blk.poles = std::move(poles);
  blk.residues = std::move(residues);
  fit.dtn_aux = std::move(blk);

  // Dense-grid residual on the realized pole-residue model vs true kₙ (q = i·kₙ).
  double max_rel = 0.0, max_abs_truth = 0.0;
  for (double w : dense_omegas)
  {
    const double truth = space_op.GetWavePortOp().GetWavePortKn(port_idx, w);
    const cplx model = EvalDtnPoleResidue(fit.dtn_poly0, fit.dtn_poly1, w,
                                          fit.dtn_aux->poles, fit.dtn_aux->residues);
    // Compare to q = i·kₙ.
    max_abs_truth = std::max(max_abs_truth, std::abs(truth));
    max_rel = std::max(max_rel, std::abs(model - cplx(0.0, truth)));
  }
  fit.rel_err_dtn = (max_abs_truth > 0.0) ? max_rel / max_abs_truth : 0.0;

  const std::size_t n_poles = fit.dtn_aux->poles.size();
  const std::size_t rank_used = fit.dtn_aux->sigmas.size();
  double max_re_pole = 0.0;
  for (const auto &p : fit.dtn_aux->poles)
  {
    max_re_pole = std::max(max_re_pole, p.real());
  }
  Mpi::Print(" Wave port {:d}: DtN structured-√ synthesis residual {:.3e} → {:.3e} "
             "(tol {:.3e}, {:d} pole{} × rank-{:d} = +{:d} aux state{}, "
             "kₙ²={:.3e}{:+.3e}·s{:+.3e}·s², max Re(pole)={:.2e})\n",
             port_idx, fit.rel_err_polynomial, fit.rel_err_dtn, waveport_synthesis_tol,
             n_poles, n_poles == 1 ? "" : "s", rank_used, n_poles * rank_used,
             n_poles * rank_used == 1 ? "" : "s", c0.real(), c1.real(), c2.real(),
             max_re_pole);
  if (fit.rel_err_dtn > waveport_synthesis_tol)
  {
    Mpi::Warning("Wave port {:d}: DtN structured-√ synthesis residual {:.3e} exceeds "
                 "WavePortSynthesisTol={:.3e}; the online surrogate falls back to the exact "
                 "per-ω cross-section EVP for this port.\n",
                 port_idx, fit.rel_err_dtn, waveport_synthesis_tol);
  }
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

std::complex<double>
RomOperator::EvaluateWavePortMultiplierFit(const WavePortDispersionFit &fit, double omega)
{
  if (fit.regime == WavePortRegime::DtnRational)
  {
    // DtN multiplier q(iω) = poly0 + poly1·(iω) + Σ R_k/(iω − p_k). The online assembly does
    // Ar += scalar·Mp_r with Mp_r = i·M_proj, but the DtN contribution to Aᵣ is q(iω)·M_proj
    // (since q = i·kₙ and the boundary term is i·kₙ·M_proj). So the scalar multiplying Mp_r
    // is q(iω)/i = −i·q(iω), matching the −i factors in ApplyDtnRationalFitCorrections.
    const std::vector<std::complex<double>> empty;
    const auto &poles = fit.dtn_aux ? fit.dtn_aux->poles : empty;
    const auto &residues = fit.dtn_aux ? fit.dtn_aux->residues : empty;
    const std::complex<double> q =
        EvalDtnPoleResidue(fit.dtn_poly0, fit.dtn_poly1, omega, poles, residues);
    return std::complex<double>(0.0, -1.0) * q;
  }
  return std::complex<double>(EvaluateWavePortKnFit(fit, omega), 0.0);
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

void RomOperator::ApplyDtnRationalFitCorrections(const WavePortDispersionFit &fit,
                                                 const Eigen::MatrixXcd &Mp_r,
                                                 Eigen::MatrixXcd &Kr_corr,
                                                 Eigen::MatrixXcd &Cr_corr)
{
  // The DtN contribution to Aᵣ(ω) is q(iω)·M_proj with q(s)=poly0+poly1·s+Σ R_k/(s−p_k).
  // The proper-polynomial part: poly0·M_proj folds into Kr, and poly1·s·M_proj = poly1·iω·
  // M_proj folds into Cr (the iω·C channel). Here M_proj = (−i)·Mp_r (Mp_r = i·M_proj), so
  // Kr += poly0·(−i)·Mp_r and Cr += poly1·(−i)·Mp_r. The pole part Σ R_k/(s−p_k)·M_proj is
  // realized by the DtN aux states in BuildAugmentedPencil.
  const Eigen::MatrixXcd M_proj = (std::complex<double>(0.0, -1.0) * Mp_r).eval();
  Kr_corr += fit.dtn_poly0 * M_proj;
  Cr_corr += fit.dtn_poly1 * M_proj;
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

bool RomOperator::AddAuxBlockDirections(WavePortAuxBlock &blk,
                                        const Eigen::MatrixXcd &Mp_r, double rank_tol)
{
  Eigen::MatrixXd M_proj = Mp_r.imag().eval();
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(M_proj, Eigen::ComputeThinU);
  if (svd.singularValues().size() == 0)
  {
    return false;
  }
  const double sigma_max = svd.singularValues()(0);
  for (long j = 0; j < svd.singularValues().size(); j++)
  {
    const double s = svd.singularValues()(j);
    if (sigma_max > 0.0 && s / sigma_max > rank_tol)
    {
      blk.sigmas.push_back(s);
      blk.u_dirs.push_back(svd.matrixU().col(j));
    }
  }
  return !blk.sigmas.empty();
}

std::optional<RomOperator::WavePortAuxBlock>
RomOperator::MakeAuxBlock(int label_idx, const Eigen::MatrixXcd &Mp_r,
                          const std::vector<std::complex<double>> &poles,
                          const std::vector<std::complex<double>> &residues,
                          double rank_tol)
{
  // Build an aux block from a projected boundary mass Mp_r (purely imaginary = i·M_proj)
  // and a pole-residue list: SVD the real symmetric M_proj to get coupling directions, then
  // attach the poles.
  if (poles.empty())
  {
    return std::nullopt;
  }
  WavePortAuxBlock blk;
  blk.port_idx = label_idx;
  AddAuxBlockDirections(blk, Mp_r, rank_tol);
  if (blk.sigmas.empty())
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

  // Complex LSQ polynomial fit of f(ω) ≈ c0 + c1 ω + c2 ω² at order 2.
  Eigen::MatrixXcd Vand(n_fit, 3);
  Eigen::VectorXcd y(n_fit);
  for (int i = 0; i < n_fit; i++)
  {
    const double w = fit_omegas[i];
    Vand(i, 0) = 1.0;
    Vand(i, 1) = w;
    Vand(i, 2) = w * w;
    y(i) = f(std::complex<double>(w, 0.0));
  }
  Eigen::Vector3cd c = Vand.colPivHouseholderQr().solve(y);
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
  Eigen::VectorXcd z_aaa(n_fit), F_aaa(n_fit);
  for (int i = 0; i < n_fit; i++)
  {
    const double w = fit_omegas[i];
    z_aaa(i) = w;
    F_aaa(i) = y(i) - (c(0) + c(1) * w + c(2) * w * w);
  }
  const double aaa_tol_rel = (max_abs_truth > 0.0)
                                 ? waveport_synthesis_tol * max_abs_truth /
                                       std::max(F_aaa.cwiseAbs().maxCoeff(), 1.0e-300)
                                 : waveport_synthesis_tol;
  auto aaa = utils::RunAAA(z_aaa, F_aaa, aaa_tol_rel,
                           std::max<std::size_t>(waveport_synthesis_order_max, 1));
  auto pr = utils::AAAToPoleResidue(aaa);
  fit.alpha0c += pr.d;  // fold AAA asymptote into the constant term

  max_rel = 0.0;
  for (int i = 0; i < n_dense; i++)
  {
    const double w = dense_omegas[i];
    const std::complex<double> truth = f(std::complex<double>(w, 0.0));
    std::complex<double> aug = fit.alpha0c + fit.alpha1c * w + fit.alpha2c * w * w;
    for (long k = 0; k < pr.poles.size(); k++)
    {
      aug += pr.residues(k) / (std::complex<double>(w, 0.0) - pr.poles(k));
    }
    max_rel = std::max(max_rel, std::abs(aug - truth));
  }
  fit.rel_err_augmented = (max_abs_truth > 0.0) ? max_rel / max_abs_truth : 0.0;

  std::vector<std::complex<double>> poles, residues;
  for (long k = 0; k < pr.poles.size(); k++)
  {
    poles.push_back(pr.poles(k));
    residues.push_back(pr.residues(k));
  }
  fit.aux = MakeAuxBlock(-1, Mp_r, poles, residues, waveport_synthesis_rank_tol);
  fit.regime = WavePortRegime::Augmented;
  Mpi::Print(
      " {}: augmented synthesis residual {:.3e} → {:.3e} (tol {:.3e}, {:d} pole{})\n",
      label, fit.rel_err_polynomial, fit.rel_err_augmented, waveport_synthesis_tol,
      pr.poles.size(), pr.poles.size() == 1 ? "" : "s");
  return fit;
}

RomOperator::AugmentedPencil RomOperator::BuildAugmentedPencil(
    const Eigen::MatrixXcd &Kr_total, const Eigen::MatrixXcd &Cr_total,
    const Eigen::MatrixXcd &Mr_total, const std::vector<WavePortAuxBlock> &aux_blocks,
    const std::vector<WavePortAuxBlock> &dtn_aux_blocks,
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
  for (const auto &blk : dtn_aux_blocks)
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
        const std::string prefix =
            blk.label.empty() ? fmt::format("waveport_{:d}", blk.port_idx) : blk.label;
        aux_labels.push_back(fmt::format("{}_p{:d}d{:d}", prefix, k, j));
        aux_row++;
      }
    }
  }
  for (const auto &blk : dtn_aux_blocks)
  {
    for (std::size_t k = 0; k < blk.poles.size(); k++)
    {
      const std::complex<double> pk = blk.poles[k];
      const std::complex<double> rk = blk.residues[k];
      for (std::size_t j = 0; j < blk.sigmas.size(); j++)
      {
        // DtN s-domain pole term R_k/(s−p_k)·M_proj. The aux state has Kr_aa=−p_k, Cr_aa=1
        // so A_aa = −p_k + iω = s − p_k, and eliminating it gives −coupling²·uⱼuⱼᵀ/(s−p_k).
        // Matching needs coupling² = −R_k·σ_j ⇒ coupling = √(−R_k·σ_j) (complex symmetric,
        // reciprocal). Lossless ports give jω-axis p_k (undamped, pure reactance); lossy give
        // LHP p_k (damped). Differs from the AAA aux (ω-domain, Cr_aa=−i): this is the
        // s-domain complex-pole realization.
        const std::complex<double> coupling = std::sqrt(-rk * blk.sigmas[j]);
        aug.Kr(aux_row, aux_row) = -pk;
        aug.Cr(aux_row, aux_row) = std::complex<double>(1.0, 0.0);
        for (long i = 0; i < n_v; i++)
        {
          aug.Kr(i, aux_row) = coupling * blk.u_dirs[j](i);
          aug.Kr(aux_row, i) = coupling * blk.u_dirs[j](i);
        }
        const std::string prefix =
            blk.label.empty() ? fmt::format("waveport_{:d}", blk.port_idx) : blk.label;
        aux_labels.push_back(fmt::format("{}_dtn{:d}d{:d}", prefix, k, j));
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
  // FitWavePortDispersion does the per-port sampling, fitting and regime selection;
  // the Apply*FitCorrections helpers accumulate the top-left contributions;
  // BuildAugmentedPencil extends the pencil with the aux states from collected fit results.
  Eigen::MatrixXcd Kr_total_corr = Eigen::MatrixXcd::Zero(Kr.rows(), Kr.cols());
  Eigen::MatrixXcd Cr_total_corr = Eigen::MatrixXcd::Zero(Kr.rows(), Kr.cols());
  Eigen::MatrixXcd Mr_total_corr = Eigen::MatrixXcd::Zero(Kr.rows(), Kr.cols());
  std::vector<WavePortAuxBlock> aux_blocks_total;
  std::vector<WavePortAuxBlock> dtn_aux_blocks_total;
  struct PendingPortLoad
  {
    std::string label;
    Eigen::MatrixXcd Kr_corr;
    Eigen::MatrixXcd Cr_corr;
    Eigen::MatrixXcd Mr_corr;
    std::vector<WavePortAuxBlock> aux_blocks;
    std::vector<WavePortAuxBlock> dtn_aux_blocks;
  };
  std::vector<PendingPortLoad> pending_port_loads;
  if (!Mwp_p_r.empty() && sweep_omega_max > sweep_omega_min)
  {
    for (auto &[port_idx, Mp_r] : Mwp_p_r)
    {
      auto fit = FitWavePortDispersion(port_idx, Mp_r);
      PendingPortLoad port_load;
      port_load.label = fmt::format("waveport_{:d}_re", port_idx);
      port_load.Kr_corr = Eigen::MatrixXcd::Zero(Kr.rows(), Kr.cols());
      port_load.Cr_corr = Eigen::MatrixXcd::Zero(Kr.rows(), Kr.cols());
      port_load.Mr_corr = Eigen::MatrixXcd::Zero(Kr.rows(), Kr.cols());
      out.wave_port_fits.push_back(fit);
      if (fit.regime == WavePortRegime::DtnRational)
      {
        ApplyDtnRationalFitCorrections(fit, Mp_r, Kr_total_corr, Cr_total_corr);
        ApplyDtnRationalFitCorrections(fit, Mp_r, port_load.Kr_corr,
                                       port_load.Cr_corr);
      }
      else
      {
        ApplyPolynomialFitCorrections(fit, Mp_r, Kr_total_corr, Cr_total_corr,
                                      Mr_total_corr);
        ApplyPolynomialFitCorrections(fit, Mp_r, port_load.Kr_corr,
                                      port_load.Cr_corr, port_load.Mr_corr);
      }
      if (fit.aux)
      {
        aux_blocks_total.push_back(*fit.aux);
        port_load.aux_blocks.push_back(*fit.aux);
      }
      if (fit.dtn_aux)
      {
        dtn_aux_blocks_total.push_back(*fit.dtn_aux);
        port_load.dtn_aux_blocks.push_back(*fit.dtn_aux);
      }
      pending_port_loads.push_back(std::move(port_load));
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
          MakeAuxBlock(/*label_idx=*/0, M_ff_r, {std::complex<double>(0.0, 0.0)},
                       {std::complex<double>(0.5, 0.0)}, waveport_synthesis_rank_tol);
      if (blk)
      {
        blk->label = "farfield";
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
        fit.aux->label = label;
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

  auto aug = BuildAugmentedPencil(Kr_total, Cr_total, Mr_total, aux_blocks_total,
                                  dtn_aux_blocks_total, out.aux_labels);

  // v_d port-row scaling: extend with 1's for aux rows (no port-impedance scaling on
  // aux states — they're internal circuit nodes).
  auto unit_henry_inv = 1.0 / units.GetScaleFactor<Units::ValueType::INDUCTANCE>();
  auto unit_farad = units.GetScaleFactor<Units::ValueType::CAPACITANCE>();
  auto unit_ohm_inv = 1.0 / units.GetScaleFactor<Units::ValueType::IMPEDANCE>();

  auto normalize_augmented = [&](const AugmentedPencil &aug_in,
                                 std::unique_ptr<mat_t> &L_inv,
                                 std::unique_ptr<mat_t> &R_inv,
                                 std::unique_ptr<mat_t> &C_out)
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
    C_out =
        std::make_unique<mat_t>((unit_farad * v_d_aug * aug_in.Mr * v_d_aug).eval());
    // Emit R⁻¹ whenever there's any dissipative contribution: lumped resistance, surface
    // conductivity, wave-port α₁, or aux-state damping.
    if (aug_in.Cr.cwiseAbs().maxCoeff() > 0.0)
    {
      R_inv =
          std::make_unique<mat_t>((unit_ohm_inv * v_d_aug * aug_in.Cr * v_d_aug).eval());
    }
  };

  normalize_augmented(aug, out.L_inv, out.R_inv, out.C);

  auto label_index = [](const std::vector<std::string> &labels,
                        const std::string &target) -> long
  {
    auto it = std::find(labels.begin(), labels.end(), target);
    if (it == labels.end())
    {
      return -1;
    }
    return static_cast<long>(std::distance(labels.begin(), it));
  };

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
      const long gi = label_index(total_labels, local_aux_labels[aux_i]);
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
    auto local_aug = BuildAugmentedPencil(pending.Kr_corr, pending.Cr_corr,
                                          pending.Mr_corr, pending.aux_blocks,
                                          pending.dtn_aux_blocks, local_aux_labels);
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
    const long row = label_index(total_labels, label);
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

  auto label_index = [](const std::vector<std::string> &labels,
                        const std::string &target) -> long
  {
    auto it = std::find(labels.begin(), labels.end(), target);
    if (it == labels.end())
    {
      return -1;
    }
    return static_cast<long>(std::distance(labels.begin(), it));
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
    if (label_index(v_node_label, label) >= 0)
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
    if (fit_it != matrices.wave_port_fits.end() && label_index(v_node_label, label) >= 0)
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

  auto wave_y_ref = [&](const RefPort &ref, double omega) -> std::complex<double>
  {
    MFEM_VERIFY(ref.wave_fit != nullptr, "Missing wave-port fit for port reference output!");
    const auto Mp_it = Mwp_p_r.find(ref.port_idx);
    MFEM_VERIFY(Mp_it != Mwp_p_r.end(), "Missing wave-port boundary mass projection!");

    std::vector<long> rows;
    const long re = label_index(v_node_label, fmt::format("waveport_{:d}_re", ref.port_idx));
    const long im = label_index(v_node_label, fmt::format("waveport_{:d}_im", ref.port_idx));
    MFEM_VERIFY(re >= 0, "Missing wave-port real row for port reference output!");
    rows.push_back(re);
    if (im >= 0)
    {
      rows.push_back(im);
    }

    const auto scalar = EvaluateWavePortMultiplierFit(*ref.wave_fit, omega);
    Eigen::MatrixXcd A = Eigen::MatrixXcd::Zero(rows.size(), rows.size());
    for (std::size_t i = 0; i < rows.size(); i++)
    {
      for (std::size_t j = 0; j < rows.size(); j++)
      {
        A(static_cast<long>(i), static_cast<long>(j)) =
            scalar * orth_R(rows[i], rows[i]) * Mp_it->second(rows[i], rows[j]) *
            orth_R(rows[j], rows[j]);
      }
    }

    std::complex<double> A_eff = A(0, 0);
    if (rows.size() > 1)
    {
      const long ni = static_cast<long>(rows.size() - 1);
      const Eigen::MatrixXcd Aii = A.bottomRightCorner(ni, ni);
      Eigen::FullPivLU<Eigen::MatrixXcd> lu(Aii);
      if (lu.rank() == ni)
      {
        A_eff -=
            (A.block(0, 1, 1, ni) * lu.solve(A.block(1, 0, ni, 1)))(0, 0);
      }
    }

    return (omega > 0.0) ? (unit_ohm_inv * A_eff / (1i * omega))
                         : std::complex<double>{0.0, 0.0};
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
        const auto &port_data = space_op.GetLumpedPortOp().GetPort(ref.port_idx);
        y_ref = unit_ohm_inv / port_data.GetCharacteristicImpedance(omega);
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
      [post_dir](const Eigen::MatrixXd &mat, std::string_view filename,
                 const std::vector<std::string> &table_labels)
  {
    MFEM_VERIFY((table_labels.size() == mat.cols()) &&
                    (table_labels.size() == mat.rows()),
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
}

}  // namespace palace
