// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "waveportreducedmodel.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include "utils/communication.hpp"

extern "C"
{
  void zggev_(char *, char *, int *, std::complex<double> *, int *, std::complex<double> *,
              int *, std::complex<double> *, std::complex<double> *, std::complex<double> *,
              int *, std::complex<double> *, int *, std::complex<double> *, int *, double *,
              int *);
}

namespace palace
{

namespace
{

// Number of vectors copied at a time into the dense panels of LocalInnerProducts.
constexpr std::size_t DENSE_PANEL_SIZE = 32;

// Copy the local entries of the complex vectors get(first), ..., get(first + count - 1),
// each of local size m, into the columns of the dense matrix X.
template <typename GetVector>
void CopyColumns(std::size_t first, std::size_t count, int m, const GetVector &get,
                 Eigen::MatrixXcd &X)
{
  X.resize(m, static_cast<Eigen::Index>(count));
  for (std::size_t j = 0; j < count; j++)
  {
    const ComplexVector &x = get(first + j);
    MFEM_VERIFY(x.Size() == m,
                "Invalid vector size for the dense wave-port PROM products!");
    const double *xr = x.Real().HostRead();
    const double *xi = x.Imag().HostRead();
    for (int i = 0; i < m; i++)
    {
      X(i, static_cast<Eigen::Index>(j)) = std::complex<double>(xr[i], xi[i]);
    }
  }
}

// Local inner products X(i, j) = a(i)ᴴ b(j) of the complex vectors a(0), ..., a(p - 1) and
// b(0), ..., b(q - 1), each of local size m, as dense matrix products of panels of the
// vectors rather than one inner product per entry. The panels bound the temporary dense
// copies to two panels of vectors. For the Gram matrix of a single list of vectors (gram,
// with b = a), only the panels on and above the diagonal are computed.
template <typename GetA, typename GetB>
Eigen::MatrixXcd LocalInnerProducts(std::size_t p, const GetA &a, std::size_t q,
                                    const GetB &b, int m, bool gram)
{
  Eigen::MatrixXcd X(static_cast<Eigen::Index>(p), static_cast<Eigen::Index>(q));
  Eigen::MatrixXcd A, B;
  for (std::size_t i0 = 0; i0 < p; i0 += DENSE_PANEL_SIZE)
  {
    const std::size_t ni = std::min(DENSE_PANEL_SIZE, p - i0);
    CopyColumns(i0, ni, m, a, A);
    for (std::size_t j0 = gram ? i0 : 0; j0 < q; j0 += DENSE_PANEL_SIZE)
    {
      const std::size_t nj = std::min(DENSE_PANEL_SIZE, q - j0);
      auto block = X.block(static_cast<Eigen::Index>(i0), static_cast<Eigen::Index>(j0),
                           static_cast<Eigen::Index>(ni), static_cast<Eigen::Index>(nj));
      if (gram && j0 == i0)
      {
        block.noalias() = A.adjoint() * A;
        continue;
      }
      CopyColumns(j0, nj, m, b, B);
      block.noalias() = A.adjoint() * B;
      if (gram)
      {
        X.block(static_cast<Eigen::Index>(j0), static_cast<Eigen::Index>(i0),
                static_cast<Eigen::Index>(nj), static_cast<Eigen::Index>(ni)) =
            block.adjoint();
      }
    }
  }
  return X;
}

}  // namespace

WavePortReducedModel::WavePortReducedModel(
    int num_modes, double eig_tol, int nd_size, int h1_size, MPI_Comm solver_comm,
    bool wave_port, const mode_assembly::ModeOperatorModel &mode_op_model,
    const ComplexOperator &opB)
  : num_modes(num_modes), eig_tol(eig_tol), nd_size(nd_size), h1_size(h1_size),
    solver_comm(solver_comm), wave_port(wave_port), mode_op_model(mode_op_model), opB(opB)
{
}

void WavePortReducedModel::ConfigureTraining(std::size_t max_basis_size)
{
  basis_cap = std::max(max_basis_size, static_cast<std::size_t>(num_modes));
}

void WavePortReducedModel::Enable(bool enable, double adaptive_tol,
                                  const ComplexOperator *truth_operator,
                                  std::complex<double> truth_omega, double truth_sigma)
{
  enabled = enable;
  stats.offline_basis_rank = basis.size();
  const std::size_t half_rank = basis.size() / 2 + basis.size() % 2;
  MFEM_VERIFY(static_cast<std::size_t>(num_modes) <=
                  std::numeric_limits<std::size_t>::max() / 4,
              "Wave-port PROM basis headroom overflow!");
  const std::size_t headroom = std::max(4 * static_cast<std::size_t>(num_modes), half_rank);
  MFEM_VERIFY(basis.size() <= std::numeric_limits<std::size_t>::max() - headroom,
              "Wave-port PROM online basis capacity overflow!");
  basis_cap = std::max(basis.size() + headroom, static_cast<std::size_t>(num_modes));
  stats.online_basis_cap = basis_cap;
  tolerance = std::max(10.0 * eig_tol, std::min(1.0e-4, 0.01 * adaptive_tol));
  if (enabled && wave_port)
  {
    BuildAffineModel(truth_operator, truth_omega, truth_sigma);
  }
}

void WavePortReducedModel::RecordExactSolve()
{
  stats.exact_solves++;
}

void WavePortReducedModel::RecordReducedSolve()
{
  has_solution = true;
  stats.reduced_solves++;
}

void WavePortReducedModel::BuildAffineModel(const ComplexOperator *truth_operator,
                                            std::complex<double> truth_omega,
                                            double truth_sigma)
{
  affine_model_ready = false;
  components.clear();
  B_actions.clear();
  B_projection.resize(0, 0);
  if (!wave_port)
  {
    return;
  }

  for (const auto &source : mode_op_model.GetComponents())
  {
    components.push_back({&source, {}, {}});
  }

  UpdateProjection(0);
  affine_model_ready = solver_comm == MPI_COMM_NULL || !components.empty();
  if (solver_comm != MPI_COMM_NULL && truth_operator && truth_omega.imag() == 0.0)
  {
    const double discrepancy =
        ValidateAffineModel(*truth_operator, truth_omega.real(), truth_sigma);
    if (!(discrepancy <= 1.0e-9))
    {
      affine_model_ready = false;
      Mpi::Warning(solver_comm,
                   "Disabling affine wave-port reduced operator: component/action "
                   "discrepancy {:.3e} exceeds 1e-9.\n",
                   discrepancy);
    }
  }
}

void WavePortReducedModel::UpdateProjection(std::size_t old_basis_size)
{
  if (solver_comm == MPI_COMM_NULL)
  {
    return;
  }
  const std::size_t n = basis.size();
  const std::size_t n_families = 1 + components.size();
  MFEM_VERIFY(old_basis_size <= n, "Invalid old affine basis size!");
  MFEM_VERIFY(n == 0 || n <= static_cast<std::size_t>(std::numeric_limits<int>::max()) /
                                 (n_families * n),
              "Affine projection matrix exceeds MPI count range!");

  // Actions of the operators on the new basis vectors.
  auto update = [&](const ComplexOperator &op, std::vector<ComplexVector> &actions)
  {
    actions.resize(n);
    for (std::size_t j = old_basis_size; j < n; j++)
    {
      actions[j].SetSize(op.Height());
      actions[j].UseDevice(true);
      op.Mult(basis[j], actions[j]);
    }
  };
  update(opB, B_actions);
  for (auto &component : components)
  {
    update(*component.source->op, component.actions);
  }

  // Projections basis[i]ᴴ actions[j] of all operators, as dense products.
  const int m = (n > 0) ? basis[0].Size() : 0;
  auto basis_vector = [this](std::size_t j) -> const ComplexVector & { return basis[j]; };
  auto action = [this](std::size_t col) -> const ComplexVector & { return Action(col); };
  Eigen::MatrixXcd projection =
      LocalInnerProducts(n, basis_vector, n_families * n, action, m, false);
  if (projection.size() > 0)
  {
    Mpi::GlobalSum(static_cast<int>(projection.size()), projection.data(), solver_comm);
  }
  const auto n_cols = static_cast<Eigen::Index>(n);
  B_projection = projection.leftCols(n_cols);
  for (std::size_t q = 0; q < components.size(); q++)
  {
    components[q].projection =
        projection.middleCols(static_cast<Eigen::Index>(q + 1) * n_cols, n_cols);
  }
  BuildActionGram();
}

const ComplexVector &WavePortReducedModel::Action(std::size_t col) const
{
  // Column col of the actions of all operators, B first: B_actions, then the actions of
  // each component, each family indexed by the basis vectors.
  const std::size_t n = basis.size(), family = col / n, j = col % n;
  return (family == 0) ? B_actions[j] : components[family - 1].actions[j];
}

void WavePortReducedModel::BuildActionGram()
{
  if (solver_comm == MPI_COMM_NULL)
  {
    return;
  }
  const std::size_t n = basis.size();
  const std::size_t n_families = 1 + components.size();
  MFEM_VERIFY(n > 0 && B_actions.size() == n,
              "Invalid affine action data for Gram construction!");
  const std::size_t dim = n_families * n;
  MFEM_VERIFY(dim <= static_cast<std::size_t>(
                         std::sqrt(static_cast<double>(std::numeric_limits<int>::max()))),
              "Affine action Gram matrix exceeds MPI count range!");

  // Gram matrix W_rowᴴ W_col of the actions, as dense products.
  auto action = [this](std::size_t col) -> const ComplexVector & { return Action(col); };
  action_gram = LocalInnerProducts(dim, action, dim, action, B_actions[0].Size(), true);
  Mpi::GlobalSum(static_cast<int>(dim * dim), action_gram.data(), solver_comm);
  action_gram = 0.5 * (action_gram + action_gram.adjoint()).eval();
  gram_direct_self_checked = false;
  gram_residual_trusted = true;
}

double WavePortReducedModel::ValidateAffineModel(const ComplexOperator &truth_operator,
                                                 double omega, double sigma) const
{
  if (solver_comm == MPI_COMM_NULL || components.empty())
  {
    return 0.0;
  }
  ComplexVector x(truth_operator.Width()), exact(truth_operator.Height()),
      affine(truth_operator.Height()), component_action(truth_operator.Height());
  x.UseDevice(true);
  exact.UseDevice(true);
  affine.UseDevice(true);
  component_action.UseDevice(true);
  double discrepancy = 0.0;
  // Deterministic random probes measure operator-level agreement without the misleading
  // relative amplification obtained from an eigenvector for which A*x is nearly null.
  for (int probe = 0; probe < 2; probe++)
  {
    linalg::SetRandom(solver_comm, x, 9137 + probe);
    linalg::Normalize(solver_comm, x);
    truth_operator.Mult(x, exact);
    affine = 0.0;
    for (const auto &component : components)
    {
      component.source->op->Mult(x, component_action);
      affine.Add(mode_op_model.EvaluateCoefficient(*component.source, omega, sigma),
                 component_action);
    }
    affine.Add(-1.0, exact);
    const double relative = linalg::Norml2(solver_comm, affine) /
                            std::max(linalg::Norml2(solver_comm, exact), 1.0e-300);
    discrepancy = std::max(discrepancy, relative);
  }
  return discrepancy;
}

bool WavePortReducedModel::AddBasisVector(const ComplexVector &x)
{
  if (solver_comm == MPI_COMM_NULL)
  {
    return false;
  }
  if (basis.size() >= basis_cap)
  {
    if (!basis_cap_warned)
    {
      Mpi::Warning(solver_comm,
                   "Wave-port PROM basis reached its capacity ({:d}); exact solves remain "
                   "authoritative but further enrichment is disabled.\n",
                   basis_cap);
      basis_cap_warned = true;
    }
    return false;
  }

  ComplexVector q(x);
  const double norm0 = linalg::Norml2(solver_comm, q);
  if (!(norm0 > 0.0))
  {
    return false;
  }
  // Two-pass modified Gram-Schmidt is sufficient for the deliberately small basis and
  // avoids making the port ROM depend on the real-valued 3D PROM orthogonalization code.
  for (int pass = 0; pass < 2; pass++)
  {
    for (const auto &v : basis)
    {
      const std::complex<double> alpha = linalg::Dot(solver_comm, q, v);
      q.Add(-alpha, v);
    }
  }
  const double norm = linalg::Norml2(solver_comm, q);
  if (!(norm > 1.0e-10 * norm0))
  {
    return false;
  }
  q *= 1.0 / norm;
  basis.push_back(std::move(q));
  return true;
}

void WavePortReducedModel::ObserveExactEigenvectors(int num_converged,
                                                    const EigenvalueSolver &eigen,
                                                    const std::vector<int> &mode_perm)
{
  if (num_converged <= 0 || basis_cap == 0)
  {
    return;
  }
  const std::size_t old_basis_size = basis.size();
  const int n = std::min(num_modes, num_converged);
  ComplexVector x(nd_size + h1_size);
  x.UseDevice(true);
  for (int i = 0; i < n; i++)
  {
    eigen.GetEigenvector(mode_perm[i], x);
    AddBasisVector(x);
  }
  if (!components.empty() && basis.size() > old_basis_size)
  {
    UpdateProjection(old_basis_size);
  }
}

bool WavePortReducedModel::TrySolve(double omega, double sigma)
{
  const int n = static_cast<int>(basis.size());
  if (!affine_model_ready || n < num_modes || solver_comm == MPI_COMM_NULL)
  {
    return false;
  }

  Eigen::MatrixXcd Ar = Eigen::MatrixXcd::Zero(n, n);
  std::vector<std::complex<double>> coefficients;
  coefficients.reserve(components.size());
  for (const auto &component : components)
  {
    const auto coeff = mode_op_model.EvaluateCoefficient(*component.source, omega, sigma);
    coefficients.push_back(coeff);
    Ar.noalias() += coeff * component.projection;
  }
  return SolveFromGram(sigma, Ar, coefficients);
}

std::optional<double> WavePortReducedModel::EvaluateGramResidual(
    const Eigen::VectorXcd &y, std::complex<double> lambda,
    const std::vector<std::complex<double>> &coefficients)
{
  const long n = static_cast<long>(basis.size());
  const long n_families = 1 + static_cast<long>(components.size());
  const long dim = n * n_families;
  MFEM_VERIFY(y.size() == n && coefficients.size() == components.size() &&
                  action_gram.rows() == dim && action_gram.cols() == dim,
              "Invalid affine Gram residual dimensions!");

  Eigen::VectorXcd cb = Eigen::VectorXcd::Zero(dim);
  Eigen::VectorXcd ca = Eigen::VectorXcd::Zero(dim);
  cb.head(n) = y;
  for (long q = 0; q < static_cast<long>(coefficients.size()); q++)
  {
    ca.segment((q + 1) * n, n) = coefficients[q] * y;
  }
  const Eigen::VectorXcd cr = cb - lambda * ca;

  auto norm_from_gram = [&](const Eigen::VectorXcd &c) -> std::optional<double>
  {
    const std::complex<double> value = c.dot(action_gram * c);
    const double scale =
        std::max(1.0, c.squaredNorm() * action_gram.cwiseAbs().maxCoeff() * dim);
    const double roundoff = 256.0 * std::numeric_limits<double>::epsilon() * scale;
    if (!std::isfinite(value.real()) || !std::isfinite(value.imag()) ||
        std::abs(value.imag()) > roundoff || value.real() < -roundoff)
    {
      return std::nullopt;
    }
    return std::sqrt(std::max(value.real(), 0.0));
  };

  const auto norm_b = norm_from_gram(cb);
  const auto norm_a = norm_from_gram(ca);
  const auto norm_r = norm_from_gram(cr);
  if (!norm_b || !norm_a || !norm_r)
  {
    return std::nullopt;
  }
  const double denom = *norm_b + std::abs(lambda) * *norm_a;
  const double eta = *norm_r / std::max(denom, std::numeric_limits<double>::min());
  if (!std::isfinite(eta))
  {
    return std::nullopt;
  }
  return eta;
}

double WavePortReducedModel::EvaluateDirectResidual(
    const Eigen::VectorXcd &y, std::complex<double> lambda,
    const std::vector<std::complex<double>> &coefficients) const
{
  const int n = static_cast<int>(basis.size());
  MFEM_VERIFY(y.size() == n && coefficients.size() == components.size(),
              "Invalid affine direct residual dimensions!");
  ComplexVector ax(nd_size + h1_size), bx(nd_size + h1_size), residual(nd_size + h1_size);
  ax.UseDevice(true);
  bx.UseDevice(true);
  residual.UseDevice(true);
  ax = 0.0;
  bx = 0.0;
  for (int j = 0; j < n; j++)
  {
    bx.Add(y(j), B_actions[j]);
  }
  for (std::size_t q = 0; q < components.size(); q++)
  {
    for (int j = 0; j < n; j++)
    {
      ax.Add(coefficients[q] * y(j), components[q].actions[j]);
    }
  }
  residual = bx;
  residual.Add(-lambda, ax);
  const double denom =
      linalg::Norml2(solver_comm, bx) + std::abs(lambda) * linalg::Norml2(solver_comm, ax);
  return linalg::Norml2(solver_comm, residual) /
         std::max(denom, std::numeric_limits<double>::min());
}

bool WavePortReducedModel::SolveFromGram(
    double sigma, const Eigen::MatrixXcd &Ar,
    const std::vector<std::complex<double>> &coefficients)
{
  int n = static_cast<int>(basis.size());
  MFEM_VERIFY(Ar.rows() == n && Ar.cols() == n && B_projection.rows() == n &&
                  B_projection.cols() == n,
              "Invalid reduced affine mode-operator dimensions!");

  Eigen::MatrixXcd Bq = B_projection, Aq = Ar;
  Eigen::VectorXcd alpha(n), beta(n);
  Eigen::MatrixXcd vl_dummy(1, 1), vr(n, n);
  int ldvl = 1, ldvr = n, lwork = std::max(4 * n, 1), info = 0;
  Eigen::VectorXcd work(lwork);
  Eigen::VectorXd rwork(std::max(8 * n, 1));
  char job_n = 'N', job_v = 'V';
  zggev_(&job_n, &job_v, &n, Bq.data(), &n, Aq.data(), &n, alpha.data(), beta.data(),
         vl_dummy.data(), &ldvl, vr.data(), &ldvr, work.data(), &lwork, rwork.data(),
         &info);
  if (info != 0)
  {
    return false;
  }

  struct Candidate
  {
    std::complex<double> lambda;
    std::complex<double> kn;
    Eigen::VectorXcd coefficients;
    double residual;
    double distance;
  };
  std::vector<Candidate> candidates;
  candidates.reserve(n);
  const double kn_target = std::sqrt(-sigma);
  for (int k = 0; k < n; k++)
  {
    const double scale = std::max(std::abs(alpha(k)), std::abs(beta(k)));
    if (!(scale > 0.0) || std::abs(beta(k)) <= 1.0e-13 * scale ||
        std::abs(alpha(k)) <= 1.0e-13 * scale)
    {
      continue;
    }
    const std::complex<double> lambda = alpha(k) / beta(k);
    const std::complex<double> kn = std::sqrt(-sigma - 1.0 / lambda);
    if (!std::isfinite(lambda.real()) || !std::isfinite(lambda.imag()) ||
        !std::isfinite(kn.real()) || !std::isfinite(kn.imag()) ||
        std::abs(kn) <= 1.0e-10 * std::max(kn_target, 1.0))
    {
      continue;
    }

    Eigen::VectorXcd y = vr.col(k);
    const double ynorm = y.norm();
    if (!(ynorm > 0.0))
    {
      continue;
    }
    y /= ynorm;

    const auto gram_eta = EvaluateGramResidual(y, lambda, coefficients);
    double eta = gram_eta.value_or(std::numeric_limits<double>::infinity());
    const bool near_acceptance_threshold =
        eta >= 0.25 * tolerance && eta <= 4.0 * tolerance;
    const bool unchecked_acceptance_candidate =
        !gram_direct_self_checked && eta <= tolerance;
    const bool verify_direct = !gram_eta || !gram_residual_trusted ||
                               unchecked_acceptance_candidate || near_acceptance_threshold;
    if (verify_direct)
    {
      const double direct_eta = EvaluateDirectResidual(y, lambda, coefficients);
      if (gram_eta && std::isfinite(direct_eta))
      {
        // Normalize by the acceptance scale: relative error against a residual already at
        // roundoff is not meaningful, while disagreement comparable to tolerance can
        // change the accept/fallback decision.
        const double discrepancy =
            std::abs(*gram_eta - direct_eta) / std::max(direct_eta, tolerance);
        if (discrepancy > 5.0e-2 && gram_residual_trusted)
        {
          gram_residual_trusted = false;
          Mpi::Warning(solver_comm,
                       "Wave-port residual Gram verification disagrees with direct "
                       "evaluation (scaled discrepancy {:.3e}); using direct residuals "
                       "for the remaining basis lifetime.\n",
                       discrepancy);
        }
      }
      eta = direct_eta;
      if (gram_eta && (*gram_eta <= tolerance || near_acceptance_threshold))
      {
        gram_direct_self_checked = true;
      }
    }
    if (!std::isfinite(eta))
    {
      continue;
    }
    if (eta > tolerance)
    {
      continue;
    }
    candidates.push_back(
        {lambda, kn, std::move(y), eta, mode_assembly::TargetDistance(kn, kn_target)});
  }

  // Rank with the same complex target distance as the exact solve (ModeEigenSolver::Solve)
  // so that both select the same mode.
  std::ranges::stable_sort(candidates, [](const Candidate &a, const Candidate &b)
                           { return a.distance < b.distance; });
  if (candidates.size() < static_cast<std::size_t>(num_modes))
  {
    return false;
  }
  if (candidates.size() > static_cast<std::size_t>(num_modes))
  {
    const double gap = candidates[num_modes].distance - candidates[num_modes - 1].distance;
    if (gap <= 1.0e-10 * std::max(kn_target, 1.0))
    {
      return false;
    }
  }

  eigenvalues.clear();
  eigenvectors.clear();
  errors.clear();
  for (int i = 0; i < num_modes; i++)
  {
    ComplexVector x(nd_size + h1_size);
    x.UseDevice(true);
    x = 0.0;
    for (int j = 0; j < n; j++)
    {
      x.Add(candidates[i].coefficients(j), basis[j]);
    }
    const double xnorm = linalg::Norml2(solver_comm, x);
    if (!(xnorm > 0.0))
    {
      eigenvalues.clear();
      eigenvectors.clear();
      errors.clear();
      return false;
    }
    x *= 1.0 / xnorm;
    stats.worst_residual = std::max(stats.worst_residual, candidates[i].residual);
    eigenvalues.push_back(candidates[i].lambda);
    eigenvectors.push_back(std::move(x));
    errors.push_back(candidates[i].residual);
  }
  return true;
}

}  // namespace palace
