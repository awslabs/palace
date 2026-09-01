// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "magnetostaticsolver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <vector>
#include <mfem.hpp>
#include "fem/errorindicator.hpp"
#include "fem/mesh.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "models/curlcurloperator.hpp"
#include "models/postoperator.hpp"
#include "models/surfacecurlsolver.hpp"
#include "models/surfacecurrentoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace
{

// Invert a matrix in place, returning false and leaving it untouched if it cannot be
// inverted. The entries are checked first because mfem's LU factorization aborts rather
// than reporting failure on a singular matrix.
bool InvertOrMarkSingular(mfem::DenseMatrix &mat)
{
  const int n = mat.Height();
  if (n == 0)
  {
    return true;
  }
  // The entries are scanned explicitly rather than via DenseMatrix::MaxMaxNorm, whose
  // running maximum never updates on a NaN comparison and so hides NaNs behind a finite
  // norm.
  double norm = 0.0;
  for (int i = 0; i < n * n; i++)
  {
    const double entry = std::abs(mat.GetData()[i]);
    if (!std::isfinite(entry))
    {
      return false;
    }
    norm = std::max(norm, entry);
  }
  if (norm == 0.0)
  {
    return false;
  }
  // Reject a matrix whose LU pivots underflow relative to its scale, which is what an
  // uninvertible block looks like here (for example a flux loop that links no flux).
  mfem::DenseMatrix lu(mat);
  mfem::Array<int> ipiv(n);
  mfem::LUFactors factors(lu.GetData(), ipiv.GetData());
  if (!factors.Factor(n, std::numeric_limits<double>::epsilon() * norm))
  {
    return false;
  }
  mat.Invert();
  return true;
}

// Fill the reluctance Minv = M^-1 and the mutual (current-difference) matrix Mm from a
// fully populated inductance matrix M. If M is not usable every entry is left as NaN.
void ComputeMinvAndMm(const mfem::DenseMatrix &M, mfem::DenseMatrix &Minv,
                      mfem::DenseMatrix &Mm, bool valid)
{
  if (!valid)
  {
    return;
  }
  const int n = M.Height();
  mfem::DenseMatrix M_inv(M);
  if (InvertOrMarkSingular(M_inv))
  {
    Minv = M_inv;
  }
  for (int i = 0; i < n; i++)
  {
    Mm(i, i) = M(i, i);
    for (int j = 0; j < n; j++)
    {
      if (i != j)
      {
        Mm(i, j) = -M(i, j);
        Mm(i, i) += M(i, j);
      }
    }
  }
}

}  // namespace

std::pair<ErrorIndicator, long long int>
MagnetostaticSolver::Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  // Construct the system matrix defining the linear operator. Dirichlet boundaries are
  // handled eliminating the rows and columns of the system matrix for the corresponding
  // dofs.
  BlockTimer bt0(Timer::CONSTRUCT);
  CurlCurlOperator curlcurl_op(iodata, mesh);
  auto K = curlcurl_op.GetStiffnessMatrix();
  const auto &Curl = curlcurl_op.GetCurlMatrix();
  SaveMetadata(curlcurl_op.GetNDSpaces());

  // Set up the linear solver. Each inactive surface current port is treated during the
  // sweep either as Open (natural BC, no current across it) or Short (PEC, screening
  // current allowed). The mode is resolved per port: a port's own "InactiveMode" overrides
  // the global "/Solver/Magnetostatic/InactivePorts" default when set.
  const InactivePortMode global_mode = iodata.solver.magnetostatic.inactive_port_mode;
  auto port_is_short = [&](int port_idx)
  {
    const auto &port_data = iodata.boundaries.current.at(port_idx);
    return port_data.inactive_port_mode.value_or(global_mode) == InactivePortMode::SHORT;
  };

  // A single linear solver is reused across all excitation steps so that its solver
  // statistics are accumulated for the metadata. Each step swaps in the stiffness matrix
  // whose essential DOFs match that step's shorted inactive ports (an empty set for a step
  // with no shorted ports reuses the base operator K). The base operator is bound lazily so
  // that an all-short configuration does not pay for an unused preconditioner setup.
  KspSolver ksp(iodata, curlcurl_op.GetNDSpaces(), &curlcurl_op.GetH1Spaces());
  const Operator *bound_op = nullptr;
  auto set_operator = [&](const Operator &op)
  {
    if (bound_op != &op)
    {
      ksp.SetOperators(op, op);
      bound_op = &op;
    }
  };

  // Surface current source indices define the boundaries over which to compute the
  // inductance matrix.
  PostOperator<ProblemType::MAGNETOSTATIC> post_op(iodata, curlcurl_op);
  int n_current_steps = static_cast<int>(curlcurl_op.GetSurfaceCurrentOp().Size());
  int n_flux_steps = static_cast<int>(curlcurl_op.GetSurfaceFluxOp().Size());
  int n_step = n_current_steps + n_flux_steps;

  MFEM_VERIFY(n_step > 0, "No surface current boundaries or flux loops specified for "
                          "magnetostatic simulation!");
  MFEM_VERIFY(
      n_flux_steps == 0 || iodata.model.refinement.max_it == 0 ||
          !iodata.model.refinement.nonconformal,
      "Flux loop excitation is only supported with conformal adaptation or no adaptation!");

  // Source term and solution vector storage.
  Vector RHS(Curl.Width()), B(Curl.Height());
  std::vector<Vector> A(n_step);
  std::vector<double> I_inc(n_step);
  std::vector<double> Phi_inc(n_step);

  // Magnetic flux linked through each current port's element apertures during each flux
  // excitation, used to recover the current-flux mutual inductance in a mixed simulation.
  // This coupling is invisible to the cross-energy formula: a current step drives zero loop
  // flux while a flux step drives zero port current, so the two states are
  // energy-orthogonal and the mutual has to be measured directly.
  mfem::DenseMatrix linked_flux(n_current_steps, n_flux_steps);
  linked_flux = 0.0;

  // Initialize structures for storing and reducing the results of error estimation.
  // In 2D, the curl is scalar: use the L2 curl space and H1 spaces for the estimator.
  CurlFluxErrorEstimator<Vector> estimator(
      curlcurl_op.GetMaterialOp(),
      curlcurl_op.GetCurlSpace(),  // RT (3D) or L2 curl (2D)
      (curlcurl_op.GetMesh().Dimension() < 3) ? curlcurl_op.GetH1Spaces()
                                              : curlcurl_op.GetNDSpaces(),
      iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
      iodata.solver.linear.estimator_mg);
  ErrorIndicator indicator;

  // Unified loop over all excitation sources (current and flux loops).
  if (n_current_steps > 0 && n_flux_steps > 0)
  {
    Mpi::Print(
        "\nComputing magnetostatic fields for {:d} current source{} and {:d} flux loop{}\n",
        n_current_steps, (n_current_steps > 1) ? "s" : "", n_flux_steps,
        (n_flux_steps > 1) ? "s" : "");
  }
  else
  {
    Mpi::Print("\nComputing magnetostatic fields for {:d} source {}\n", n_step,
               (n_step > 1) ? "boundaries" : "boundary");
  }
  auto t0 = Timer::Now();

  // Pre-allocate boundary values vector for flux loop optimization
  Vector boundary_values;

  // Shorted-attribute set of each current-source step and its sorted, deduplicated key
  // (empty = base operator). Steps sharing a key share a cached screened operator.
  std::vector<mfem::Array<int>> short_attrs_per_step(n_current_steps);
  std::vector<std::vector<int>> step_keys(n_step);
  {
    int step = 0;
    for (const auto &[idx, data] : curlcurl_op.GetSurfaceCurrentOp())
    {
      auto &short_attrs = short_attrs_per_step[step];
      for (const auto &[other_idx, other_data] : curlcurl_op.GetSurfaceCurrentOp())
      {
        if (other_idx != idx && port_is_short(other_idx))
        {
          for (const auto &elem : other_data.elements)
          {
            short_attrs.Append(elem.source->GetAttrList());
          }
        }
      }
      std::set<int> unique_attr(short_attrs.begin(), short_attrs.end());
      step_keys[step++].assign(unique_attr.begin(), unique_attr.end());
    }
  }

  // Order steps by each key's first appearance so steps sharing an operator are adjacent
  // and the bound operator (and its preconditioner) is reused across the group;
  // already-grouped orderings are untouched. Results go to the canonical slot A[step], so
  // postprocessing is unchanged.
  std::map<std::vector<int>, int> key_first_seen;
  for (int step = 0; step < n_step; step++)
  {
    key_first_seen.emplace(step_keys[step], step);
  }
  std::vector<int> solve_order(n_step);
  std::iota(solve_order.begin(), solve_order.end(), 0);
  std::stable_sort(
      solve_order.begin(), solve_order.end(), [&](int a, int b)
      { return key_first_seen.at(step_keys[a]) < key_first_seen.at(step_keys[b]); });

  // Pass 1: solve each excitation in operator-grouped order.
  int solve_it = 0;
  for (int step : solve_order)
  {
    A[step].SetSize(RHS.Size());
    A[step].UseDevice(true);
    A[step] = 0.0;
    solve_it++;

    if (step < n_current_steps)
    {
      // Current source excitation
      const auto &[idx, data] = *std::next(curlcurl_op.GetSurfaceCurrentOp().begin(), step);
      Mpi::Print("\nIt {:d}/{:d}: Current Index = {:d} (elapsed time = {:.2e} s)\n",
                 solve_it, n_step, idx, Timer::Duration(Timer::Now() - t0).count());

      curlcurl_op.GetCurrentExcitationVector(idx, RHS);
      I_inc[step] = data.GetExcitationCurrent();
      Phi_inc[step] = 0.0;  // Zero flux for current sources

      const auto &short_attrs = short_attrs_per_step[step];
      if (short_attrs.Size() > 0)
      {
        // Fetch the cached stiffness matrix with the shorted inactive ports added as
        // essential (PEC) boundaries, and zero the excitation on those DOFs so DIAG_ONE
        // elimination injects no spurious values on edges shared with the active port.
        const Operator &K_step = curlcurl_op.GetScreenedStiffnessMatrix(short_attrs);
        curlcurl_op.ZeroEssentialTrueDofs(short_attrs, RHS);
        set_operator(K_step);
        ksp.Mult(RHS, A[step]);
      }
      else
      {
        // No inactive ports shorted for this step: all inactive ports are open.
        set_operator(*K);
        ksp.Mult(RHS, A[step]);
      }
    }
    else
    {
      // Flux loop excitation
      int flux_idx = step - n_current_steps;
      const auto &[idx, data] =
          *std::next(curlcurl_op.GetSurfaceFluxOp().begin(), flux_idx);
      Mpi::Print("\nIt {:d}/{:d}: FluxLoop Index = {:d} (elapsed time = {:.2e} s)\n",
                 solve_it, n_step, idx, Timer::Duration(Timer::Now() - t0).count());

      curlcurl_op.GetFluxExcitationVector(idx, RHS, post_op, &boundary_values);
      I_inc[step] = 0.0;                         // Zero current for flux loops
      Phi_inc[step] = data.GetExcitationFlux();  // Store prescribed flux

      // Solve 3D magnetostatic problem (flux loops use the base operator).
      set_operator(*K);
      ksp.Mult(RHS, A[step]);
    }
  }

  // Pass 2: postprocess in canonical step order so output is independent of solve order.
  for (int step = 0; step < n_step; step++)
  {
    Curl.Mult(A[step], B);

    // Flux verification for flux loops.
    if (step >= n_current_steps)
    {
      int flux_idx = step - n_current_steps;
      const auto &[idx, data] =
          *std::next(curlcurl_op.GetSurfaceFluxOp().begin(), flux_idx);
      auto &B_gf = post_op.GetBGridFunction().Real();
      B_gf.SetFromTrueDofs(B);

      // Create flux direction vector from data
      mfem::Vector flux_direction(const_cast<double *>(data.direction.data()),
                                  static_cast<int>(data.direction.size()));
      // Verify that the flux is correctly computed
      VerifyFluxThroughHoles(B_gf, data.hole_attributes, data.flux_amounts,
                             curlcurl_op.GetMesh(), curlcurl_op.GetMaterialOp(),
                             flux_direction, curlcurl_op.GetComm());

      // Measure the flux this excitation links through each current port's aperture. With
      // the current ports open (I_c = 0) the response is Φ_c = M_cf*M_ff^-1*Φ_f, which is
      // what makes the current-flux mutual inductance recoverable.
      int col = 0;
      for (const auto &current : curlcurl_op.GetSurfaceCurrentOp())
      {
        const auto &j_data = current.second;
        for (const auto &elem : j_data.elements)
        {
          MFEM_VERIFY(elem.aperture,
                      "Missing SurfaceCurrent aperture after configuration validation!");
          linked_flux(col, flux_idx) +=
              elem.current_fraction *
              ComputeFluxThroughSurface(B_gf, elem.aperture->attributes,
                                        curlcurl_op.GetMesh(), curlcurl_op.GetMaterialOp(),
                                        elem.aperture->direction, curlcurl_op.GetComm());
        }
        col++;
      }
    }

    // Energy calculation and error estimation.
    int terminal_idx =
        (step < n_current_steps)
            ? std::next(curlcurl_op.GetSurfaceCurrentOp().begin(), step)->first
            : std::next(curlcurl_op.GetSurfaceFluxOp().begin(), step - n_current_steps)
                  ->first;
    auto total_domain_energy = post_op.MeasureAndPrintAll(step, A[step], B, terminal_idx);
    estimator.AddErrorIndicator(B, total_domain_energy, indicator);
  }

  // Postprocess the inductance matrix from the computed field solutions.
  BlockTimer bt1(Timer::POSTPRO);
  SaveMetadata(ksp);
  PostprocessTerminals(post_op, curlcurl_op.GetSurfaceCurrentOp(),
                       curlcurl_op.GetSurfaceFluxOp(), A, I_inc, Phi_inc, linked_flux);
  post_op.MeasureFinalize(indicator);
  return {indicator, curlcurl_op.GlobalTrueVSize()};
}

void MagnetostaticSolver::PostprocessTerminals(
    PostOperator<ProblemType::MAGNETOSTATIC> &post_op,
    const SurfaceCurrentOperator &surf_j_op, const SurfaceFluxOperator &surf_flux_op,
    const std::vector<Vector> &A, const std::vector<double> &I_inc,
    const std::vector<double> &Phi_inc, const mfem::DenseMatrix &linked_flux) const
{
  // Postprocess the Maxwell inductance matrix. See p. 97 of the COMSOL AC/DC Module manual
  // for the associated formulas based on the magnetic field energy based on a current
  // excitation for each port:
  //                          M_ij = (A_j^T*K*A_i)/(I_i*I_j)
  // Alternatively, we could compute the resulting loop fluxes to get M directly as:
  //                         Φ_i = ∫ B ⋅ n_j dS
  // and M_ij = Φ_i/I_j. The energy formulation avoids having to locally integrate B =
  // ∇ x A.
  // If flux excitation is employed, inductance matrix is computed by first computing
  // the reluctance:
  //                          R_ij = (A_j^T*K*A_i)/(Φ_i*Φ_j)
  // and then M = R^-1.
  //
  // With both current and flux excitations present the normalized cross-energies form a
  // hybrid (mixed impedance/admittance) matrix H rather than an inductance matrix: each
  // column is normalized by its own excitation, so a current column carries units of H and
  // a flux column of 1/H. Writing c for the current ports and f for the flux ports, a
  // current step holds the flux loops at zero flux (their PEC surfaces are part of the base
  // essential set) while a flux step leaves the current ports open (I_c = 0), so
  //                  H_cc = M_cc - M_cf*M_ff^-1*M_fc      (screened, a Schur complement)
  //                  H_ff = M_ff^-1                        (flux-loop reluctance)
  // The off-diagonal cross-energy H_cf, however, vanishes identically: pairing a current
  // state with a flux state gives E = Σ_c I_c*Φ_c + Σ_f I_f*Φ_f, in which the first sum is
  // zero because the flux state carries no port current and the second is zero because the
  // current state links no loop flux. The two states are energy-orthogonal, so the mutual
  // inductance cannot be obtained from cross-energies and has to be measured directly. A
  // flux excitation responds with Φ_c = M_cf*M_ff^-1*Φ_f, so integrating B ⋅ n over each
  // current port's aperture during each flux step gives G = M_cf*M_ff^-1 and hence
  //                  M_ff = H_ff^-1
  //                  M_cf = G*M_ff
  //                  M_cc = H_cc + G*M_ff*G^T
  // recovering the bare inductance matrix in henries; only H_ff has to be inverted, exactly
  // as in the pure flux case. If H_ff is singular the flux rows and columns, and the
  // un-screening correction to M_cc, are not defined and the affected entries are reported
  // as NaN. Mixed configurations require a complete oriented aperture for every current
  // element, so terminal-M.csv never mixes bare and screened inductance entries.
  //
  // Reciprocity note for Short inactive ports: the reciprocal cross-energy formula only
  // holds when every excitation is solved with the same operator (same essential-DOF set).
  // A port shorted while inactive is added to the essential set on every step except the
  // one that excites it, so exciting a Short port uses a different operator than the others
  // and its mutual entries are not reciprocal (and reflect the screened field of a
  // different boundary-value problem, not a coupling). Ports that are Open when inactive,
  // by contrast, are never essential, so every excitation of an Open port shares one common
  // operator (the Short ports acting as fixed screens) and its mutual entries are
  // well-defined and reciprocal. We therefore report self-inductances for all ports, but
  // mutual inductances only between ports that are Open when inactive; other off-diagonals
  // are set to NaN and Minv/Mm are computed over the Open-Open sub-block only.
  int n_current = static_cast<int>(surf_j_op.Size());
  int n_flux = static_cast<int>(surf_flux_op.Size());
  int n = A.size();

  // Mark which columns have a well-defined reciprocal mutual (Open-when-inactive ports).
  // Flux loops always share one operator, so they are all reciprocal.
  const InactivePortMode global_mode = iodata.solver.magnetostatic.inactive_port_mode;
  std::vector<bool> reciprocal(n, true);
  {
    int col = 0;
    for (const auto &[idx, data] : surf_j_op)
    {
      const auto &port_cfg = iodata.boundaries.current.at(idx);
      // A single current port is never inactive during its own (only) solve, so its
      // self-inductance is reciprocal regardless of the chosen mode.
      reciprocal[col++] = n_current == 1 || port_cfg.inactive_port_mode.value_or(
                                                global_mode) == InactivePortMode::OPEN;
    }
    // Remaining columns are flux loops (reciprocal), already initialized to true.
  }
  auto reciprocal_pair = [&](int i, int j) { return reciprocal[i] && reciprocal[j]; };
  const double nan = std::numeric_limits<double>::quiet_NaN();

  // Allocate final result matrices
  mfem::DenseMatrix M(n), Minv(n), Mm(n);
  M = nan;
  Minv = nan;
  Mm = nan;

  // Compute cross-energy matrix and diagonals. Off-diagonals are only meaningful between
  // reciprocal (Open-Open) port pairs; other pairs are left as NaN.
  mfem::DenseMatrix cross_energy(n);
  cross_energy = nan;
  for (int i = 0; i < n; i++)
  {
    auto &A_gf = post_op.GetAGridFunction().Real();
    auto &H_gf = post_op.GetDomainPostOp().H;
    A_gf.SetFromTrueDofs(A[i]);
    post_op.GetDomainPostOp().M_mag->Mult(A_gf, H_gf);
    cross_energy(i, i) = linalg::Dot<Vector>(post_op.GetComm(), A_gf, H_gf);

    // Off-diagonal cross-energies (only for reciprocal Open-Open pairs).
    for (int j = i + 1; j < n; j++)
    {
      if (!reciprocal_pair(i, j))
      {
        continue;
      }
      A_gf.SetFromTrueDofs(A[j]);
      cross_energy(i, j) = cross_energy(j, i) =
          linalg::Dot<Vector>(post_op.GetComm(), A_gf, H_gf);
    }
  }

  if (n_flux == n)
  {
    // Pure flux case: compute reluctance first, then get inductance by inversion. All flux
    // loops are reciprocal, so this always spans the full matrix.
    for (int i = 0; i < n; i++)
    {
      Minv(i, i) = cross_energy(i, i) / (Phi_inc[i] * Phi_inc[i]);
      for (int j = i + 1; j < n; j++)
      {
        Minv(i, j) = Minv(j, i) = cross_energy(i, j) / (Phi_inc[i] * Phi_inc[j]);
      }
    }
    // Get inductance from reluctance: M = R^{-1}
    M = Minv;
    M.Invert();
  }
  else if (n_current > 0 && n_flux > 0)
  {
    // Mixed current/flux case: normalize each column by its own excitation to build the
    // hybrid matrix H, then convert the blocks to a single inductance matrix in henries.
    // Current columns come first, followed by flux columns. Only the diagonal blocks of H
    // are meaningful: the current-flux cross-energies vanish identically (the two
    // excitation types are energy-orthogonal), so the coupling comes from the measured
    // linked flux.
    auto excitation = [&](int i) { return (i < n_current) ? I_inc[i] : Phi_inc[i]; };
    mfem::DenseMatrix H(n);
    H = nan;
    for (int i = 0; i < n; i++)
    {
      H(i, i) = cross_energy(i, i) / (excitation(i) * excitation(i));
      for (int j = i + 1; j < n; j++)
      {
        // Skip the current-flux block, whose cross-energy is identically zero.
        if ((i < n_current) != (j < n_current))
        {
          continue;
        }
        H(i, j) = H(j, i) = cross_energy(i, j) / (excitation(i) * excitation(j));
      }
    }

    // G = M_cf*M_ff^-1 from the flux linked through each current port's aperture, which is
    // the only place the current-flux coupling is observable.
    mfem::DenseMatrix G(n_current, n_flux);
    for (int i = 0; i < n_current; i++)
    {
      for (int b = 0; b < n_flux; b++)
      {
        G(i, b) = linked_flux(i, b) / Phi_inc[n_current + b];
      }
    }
    // Invert the flux-flux reluctance block H_ff to get the flux-loop inductance block.
    mfem::DenseMatrix Mff(n_flux);
    for (int a = 0; a < n_flux; a++)
    {
      for (int b = 0; b < n_flux; b++)
      {
        Mff(a, b) = H(n_current + a, n_current + b);
      }
    }
    const bool flux_invertible = InvertOrMarkSingular(Mff);

    // M_ff = H_ff^-1.
    for (int a = 0; a < n_flux; a++)
    {
      for (int b = 0; b < n_flux; b++)
      {
        M(n_current + a, n_current + b) = flux_invertible ? Mff(a, b) : nan;
      }
    }
    // M_cf = G*M_ff, and its transpose M_fc.
    for (int i = 0; i < n_current; i++)
    {
      for (int b = 0; b < n_flux; b++)
      {
        double sum = 0.0;
        for (int a = 0; a < n_flux; a++)
        {
          sum += G(i, a) * Mff(a, b);
        }
        M(i, n_current + b) = M(n_current + b, i) = flux_invertible ? sum : nan;
      }
    }
    // M_cc = H_cc + G*M_ff*G^T un-screens the current-port block, which the zero-flux loops
    // would otherwise report as a Schur complement.
    for (int i = 0; i < n_current; i++)
    {
      for (int j = i; j < n_current; j++)
      {
        // correction = (G*M_ff*G^T)_ij, reusing the already-computed M_cf = G*M_ff row.
        double correction = 0.0;
        for (int b = 0; b < n_flux; b++)
        {
          correction += M(i, n_current + b) * G(j, b);
        }
        M(i, j) = M(j, i) = flux_invertible ? H(i, j) + correction : nan;
      }
    }

    ComputeMinvAndMm(M, Minv, Mm, flux_invertible);
  }
  else
  {
    // Pure current case: compute inductance directly. Self-inductances (diagonal) are
    // always defined; mutuals only for reciprocal Open-Open pairs, otherwise left as NaN.
    for (int i = 0; i < n; i++)
    {
      M(i, i) = cross_energy(i, i) / (I_inc[i] * I_inc[i]);
      for (int j = i + 1; j < n; j++)
      {
        if (!reciprocal_pair(i, j))
        {
          continue;
        }
        M(i, j) = M(j, i) = cross_energy(i, j) / (I_inc[i] * I_inc[j]);
      }
    }
    // Reluctance and mutual matrices are only well-defined over the reciprocal sub-block.
    // Extract the Open-Open sub-block, invert/transform it, and scatter back;
    // non-reciprocal rows/columns remain NaN.
    std::vector<int> recip_idx;
    for (int i = 0; i < n; i++)
    {
      if (reciprocal[i])
      {
        recip_idx.push_back(i);
      }
    }
    const int nr = static_cast<int>(recip_idx.size());
    if (nr > 0)
    {
      mfem::DenseMatrix M_sub(nr), Minv_sub(nr);
      for (int a = 0; a < nr; a++)
      {
        for (int b = 0; b < nr; b++)
        {
          M_sub(a, b) = M(recip_idx[a], recip_idx[b]);
        }
      }
      // Reluctance R = M^{-1} over the sub-block.
      Minv_sub = M_sub;
      Minv_sub.Invert();
      // Mutual inductance matrix Mm (current-difference form) over the sub-block.
      mfem::DenseMatrix Mm_sub(nr);
      for (int a = 0; a < nr; a++)
      {
        Mm_sub(a, a) = M_sub(a, a);
        for (int b = 0; b < nr; b++)
        {
          if (a != b)
          {
            Mm_sub(a, b) = -M_sub(a, b);
            Mm_sub(a, a) += M_sub(a, b);
          }
        }
      }
      for (int a = 0; a < nr; a++)
      {
        for (int b = 0; b < nr; b++)
        {
          Minv(recip_idx[a], recip_idx[b]) = Minv_sub(a, b);
          Mm(recip_idx[a], recip_idx[b]) = Mm_sub(a, b);
        }
      }
    }
  }

  if (n_flux == n)
  {
    // Compute Mm matrix from final M (full matrix in the flux case).
    for (int i = 0; i < n; i++)
    {
      Mm(i, i) = M(i, i);
      for (int j = 0; j < n; j++)
      {
        if (i != j)
        {
          Mm(i, j) = -M(i, j);
          Mm(i, i) += M(i, j);
        }
      }
    }
  }

  // Only root writes to disk
  if (!root)
  {
    return;
  }

  // Write matrix data using existing pattern
  auto PrintMatrix = [&surf_j_op, &surf_flux_op, this, n_current,
                      n_flux](const std::string &file, const std::string &name,
                              const std::string &unit, const mfem::DenseMatrix &mat,
                              double scale)
  {
    TableWithCSVFile output(post_dir / file);
    output.table.insert(Column("i", "i", 0, 0, 2, ""));

    auto AddTerminal = [&](int idx, int col_pos)
    {
      output.table.insert(fmt::format("i2{}", idx),
                          fmt::format("{}[i][{}] {}", name, idx, unit));
      output.table["i"] << idx;
      auto &col = output.table[fmt::format("i2{}", idx)];
      for (int i = 0; i < mat.Height(); i++)
      {
        col << mat(i, col_pos) * scale;
      }
    };

    // Add current sources
    int j = 0;
    for (const auto &[idx, data] : surf_j_op)
      AddTerminal(idx, j++);

    // Add flux loops
    for (const auto &[idx, data] : surf_flux_op)
      AddTerminal(idx, j++);

    output.WriteFullTableTrunc();
  };

  const double H = iodata.units.GetScaleFactor<Units::ValueType::INDUCTANCE>();
  PrintMatrix("terminal-M.csv", "M", "(H)", M, H);
  PrintMatrix("terminal-Minv.csv", "M⁻¹", "(1/H)", Minv, 1.0 / H);
  PrintMatrix("terminal-Mm.csv", "M_m", "(H)", Mm, H);

  // Write out a file with source current excitations.
  if (n_current > 0)
  {
    TableWithCSVFile terminal_I(post_dir / "terminal-I.csv");
    terminal_I.table.insert(Column("i", "i", 0, 0, 2, ""));
    terminal_I.table.insert("Iinc", "I_inc[i] (A)");
    int i = 0;
    for (const auto &[idx, data] : surf_j_op)
    {
      terminal_I.table["i"] << double(idx);
      terminal_I.table["Iinc"] << iodata.units.Dimensionalize<Units::ValueType::CURRENT>(
          I_inc[i]);
      i++;
    }
    terminal_I.WriteFullTableTrunc();
  }

  // Write out a file with flux loop excitations.
  if (n_flux > 0)
  {
    TableWithCSVFile terminal_Phi(post_dir / "terminal-Phi.csv");
    terminal_Phi.table.insert(Column("i", "i", 0, 0, 2, ""));
    terminal_Phi.table.insert("Phiinc", "Phi_inc[i] (Wb)");
    const double magnetic_flux_scale =
        iodata.units.GetScaleFactor<Units::ValueType::INDUCTANCE>() *
        iodata.units.GetScaleFactor<Units::ValueType::CURRENT>();
    int i = n_current;
    for (const auto &[idx, data] : surf_flux_op)
    {
      terminal_Phi.table["i"] << double(idx);
      terminal_Phi.table["Phiinc"] << Phi_inc[i] * magnetic_flux_scale;
      i++;
    }
    terminal_Phi.WriteFullTableTrunc();
  }
}

}  // namespace palace
