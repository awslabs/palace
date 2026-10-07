// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "drivensolver.hpp"

#include <complex>
#include <cstddef>
#include <iostream>
#include <numbers>
#include <Eigen/Dense>
#include <fmt/core.h>
#include <mfem.hpp>
#include "fem/errorindicator.hpp"
#include "fem/mesh.hpp"
#include "fem/substructure.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/floquetcorrection.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/drivensubstructure.hpp"
#include "models/floquetportoperator.hpp"
#include "models/lumpedportoperator.hpp"
#include "models/portexcitations.hpp"
#include "models/postoperator.hpp"
#include "models/romoperator.hpp"
#include "models/spaceoperator.hpp"
#include "models/surfacecurrentoperator.hpp"
#include "models/waveportoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"
#include "utils/timer.hpp"

namespace palace
{

using namespace std::complex_literals;

std::pair<ErrorIndicator, long long int>
DrivenSolver::Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  // Set up the spatial discretization and frequency sweep.
  BlockTimer bt0(Timer::CONSTRUCT);
  SpaceOperator space_op(iodata, mesh);
  const auto &port_excitations = space_op.GetPortExcitations();
  SaveMetadata(port_excitations);

  const auto &omega_sample = iodata.solver.driven.sample_f;

  bool adaptive = (iodata.solver.driven.adaptive_tol > 0.0);
  if (adaptive && omega_sample.size() <= iodata.solver.driven.prom_indices.size() &&
      !iodata.solver.driven.adaptive_circuit_synthesis)
  {
    Mpi::Warning("Adaptive frequency sweep requires > {} total frequency samples!\n"
                 "Reverting to uniform sweep!\n",
                 iodata.solver.driven.prom_indices.size());
    adaptive = false;
  }
  SaveMetadata(space_op.GetNDSpaces());
  Mpi::Print("\nComputing {}frequency response for:\n{}", adaptive ? "adaptive fast " : "",
             port_excitations.FmtLog());

  std::size_t restart = iodata.solver.driven.restart;
  if (restart != 1)
  {
    std::size_t max_iter = omega_sample.size() * space_op.GetPortExcitations().Size();
    MFEM_VERIFY(
        restart - 1 < max_iter,
        fmt::format("\"Restart\" ({}) is greater than the number of total samples ({})!",
                    restart, max_iter));

    Mpi::Print("\nRestarting from solve {}", iodata.solver.driven.restart);
  }

  // Main frequency sweep loop.
  if (iodata.solver.substructuring)
  {
    return {SweepSubstructured(space_op), space_op.GlobalTrueVSize()};
  }
  return {adaptive ? SweepAdaptive(space_op) : SweepUniform(space_op),
          space_op.GlobalTrueVSize()};
}

namespace
{

// The lumped ports on the environment side, and an abort for ports touching Γ (their
// excitation and voltage would depend on both sides).
std::vector<int> EnvironmentPorts(const IoData &iodata, const DrivenSubstructure &ds,
                                  std::map<int, std::vector<char>> &port_dofs)
{
  std::vector<int> env_ports;
  for (const auto &[idx, data] : iodata.boundaries.lumpedport)
  {
    std::vector<int> attrs;
    for (const auto &elem : data.elements)
    {
      attrs.insert(attrs.end(), elem.attributes.begin(), elem.attributes.end());
    }
    auto mark = ds.BoundaryTrueDofs(attrs);
    double n[3] = {0.0, 0.0, 0.0};  // environment interior, interface, region interior
    for (std::size_t i = 0; i < mark.size(); i++)
    {
      n[0] += mark[i] && ds.EnvironmentInterior()[i];
      n[1] += mark[i] && ds.Interface()[i];
      n[2] += mark[i] && ds.RegionInterior()[i];
    }
    Mpi::GlobalSum(3, n, ds.GetComm());
    MFEM_VERIFY(n[1] == 0.0 && (n[0] == 0.0 || n[2] == 0.0),
                "Lumped port " << idx
                               << " touches the substructuring interface: a saved driven "
                                  "model needs the lumped ports away from it!");
    if (n[0] > 0.0)
    {
      env_ports.push_back(idx);
      port_dofs[idx] = std::move(mark);
    }
  }
  return env_ports;
}

// The voltage functional of a lumped port as a true-DOF vector l (V = l^T E), evaluated
// with the port's own voltage postprocessing on the unit vectors of the true DOFs it
// touches.
ComplexVector VoltageFunctional(SpaceOperator &space_op, int idx,
                                const std::vector<char> &mark)
{
  auto &nd = space_op.GetNDSpace();
  const auto &fes = nd.Get();
  MPI_Comm comm = space_op.GetComm();
  const int nt = fes.GetTrueVSize();
  const HYPRE_BigInt tstart = fes.GetMyTDofOffset();
  std::vector<HYPRE_BigInt> mine, all;
  for (int i = 0; i < nt; i++)
  {
    if (mark[i])
    {
      mine.push_back(tstart + i);
    }
  }
  {
    const int nranks = Mpi::Size(comm);
    int nloc = static_cast<int>(mine.size());
    std::vector<int> cnt(nranks), disp(nranks, 0);
    MPI_Allgather(&nloc, 1, MPI_INT, cnt.data(), 1, MPI_INT, comm);
    for (int r = 1; r < nranks; r++)
    {
      disp[r] = disp[r - 1] + cnt[r - 1];
    }
    all.resize(disp[nranks - 1] + cnt[nranks - 1]);
    MPI_Allgatherv(mine.data(), nloc, HYPRE_MPI_BIG_INT, all.data(), cnt.data(),
                   disp.data(), HYPRE_MPI_BIG_INT, comm);
  }
  const auto &port = space_op.GetLumpedPortOp().GetPort(idx);
  GridFunction E(nd, true);
  E.Imag() = 0.0;
  Vector e(nt);
  ComplexVector l(nt);
  l = 0.0;
  for (HYPRE_BigInt q : all)
  {
    const HYPRE_BigInt i = q - tstart;
    const bool own = (i >= 0 && i < nt);
    e = 0.0;
    if (own)
    {
      e(static_cast<int>(i)) = 1.0;
    }
    E.Real().SetFromTrueDofs(e);
    const std::complex<double> V = port.GetVoltage(E);
    if (own)
    {
      l.Real()(static_cast<int>(i)) = V.real();
    }
  }
  return l;
}

// The index of each value of x in y (relative tolerance tol), -1 if none.
int FindValue(const std::vector<double> &y, double x, double tol = 1.0e-10)
{
  for (std::size_t j = 0; j < y.size(); j++)
  {
    if (std::abs(y[j] - x) <= tol * std::max(std::abs(x), std::abs(y[j])))
    {
      return static_cast<int>(j);
    }
  }
  return -1;
}

// Fingerprints match: the first nreal entries relative to themselves, then complex values
// (pairs of entries) relative to their magnitude (a real or imaginary part can be zero).
bool SameFingerprint(const std::vector<double> &a, const std::vector<double> &b, int nreal,
                     double tol = 1.0e-10)
{
  if (a.size() != b.size() || (a.size() - nreal) % 2 != 0)
  {
    return false;
  }
  for (std::size_t q = 0; q < a.size(); q++)
  {
    const bool complex = (q >= static_cast<std::size_t>(nreal));
    const std::complex<double> za(a[q], complex ? a[q + 1] : 0.0),
        zb(b[q], complex ? b[q + 1] : 0.0);
    if (std::abs(za - zb) > tol * std::max(std::abs(za), std::abs(zb)))
    {
      return false;
    }
    q += complex;
  }
  return true;
}

}  // namespace

ErrorIndicator DrivenSolver::SweepSubstructured(SpaceOperator &space_op) const
{
  // Exact per-frequency substructuring: at each frequency, condense the environment and
  // factor the region against it once, then solve all excitations (so the frequency loop
  // is the outer one). With a saved model, an offline sweep saves S_E(ω) and the
  // environment's source and port data per frequency, and an online sweep (at saved
  // frequencies) factors the region only.
  const auto &port_excitations = space_op.GetPortExcitations();
  const auto &omega_sample = iodata.solver.driven.sample_f;
  const auto &sub = *iodata.solver.substructuring;
  const bool online = (sub.mode == SubstructuringMode::ONLINE);
  const std::string &model_path = sub.save_model;
  const bool save = !online && !model_path.empty();
  MPI_Comm comm = space_op.GetComm();
  const bool root = Mpi::Root(comm);
  DrivenSubstructure ds(space_op, sub.region_attributes, sub.environment_attributes,
                        online);
  const int nG = ds.InterfaceSize();
  Mpi::Print("\nSubstructuring: |Γ| = {:d} interface unknowns\n", nG);
  const auto &Curl = space_op.GetCurlMatrix();
  std::vector<int> ex_idx;
  for (const auto &[idx, spec] : port_excitations)
  {
    ex_idx.push_back(idx);
  }
  const int n = static_cast<int>(ex_idx.size());
  std::vector<ComplexVector> rhs(n, ComplexVector(Curl.Width())), E;
  std::vector<const ComplexVector *> rhs_ptr(n);
  for (int k = 0; k < n; k++)
  {
    rhs[k].UseDevice(true);
    rhs_ptr[k] = &rhs[k];
  }

  // The lumped ports of the environment (for a saved model).
  std::map<int, std::vector<char>> port_dofs;
  std::vector<int> env_ports;
  if (save || online)
  {
    env_ports = EnvironmentPorts(iodata, ds, port_dofs);
  }

  // Online, the field is known in the region (and on Γ) only: energies of region domains,
  // without totals or participation ratios, and the environment's port voltages from the
  // model (without their power).
  config::DomainData domains = iodata.domains;
  if (online)
  {
    domains.postpro.partial = true;
    std::vector<int> skipped;
    for (auto it = domains.postpro.energy.begin(); it != domains.postpro.energy.end();)
    {
      const bool region =
          std::ranges::all_of(it->second.attributes,
                              [&](int a)
                              {
                                return std::ranges::find(sub.region_attributes, a) !=
                                       sub.region_attributes.end();
                              });
      if (!region)
      {
        skipped.push_back(it->first);
        it = domains.postpro.energy.erase(it);
      }
      else
      {
        ++it;
      }
    }
    Mpi::Warning("Online driven substructuring writes no total energies or participation "
                 "ratios{}{}!\n",
                 skipped.empty() ? std::string()
                                 : fmt::format(", no energies of the domain postprocessing "
                                               "outside the region ({})",
                                               fmt::join(skipped, ", ")),
                 env_ports.empty() ? std::string()
                                   : fmt::format(", and no power of the environment's "
                                                 "lumped ports ({})",
                                                 fmt::join(env_ports, ", ")));
  }
  PostOperator<ProblemType::DRIVEN> post_op(iodata.problem, iodata.solver, domains,
                                            iodata.boundaries, iodata.units, space_op);

  // The model: header and environment ports.
  DrivenSubstructureModel model;
  std::vector<ComplexVector> port_l;
  SignatureMap gamma_map;           // online: the saved interface on the current one
  std::vector<int> record, ex_col;  // online: frequency records, model columns
  if (save)
  {
    model.nG = nG;
    model.sig_w = SignatureWidth(space_op.GetNDSpace().Get());
    model.signatures =
        TrueDofSignatures(space_op.GetNDSpace().Get(), ds.InterfaceIndex(), nG);
    model.env_fp = ds.EnvironmentFingerprint(omega_sample[0]);
    for (int k = 0; k < n; k++)
    {
      space_op.GetExcitationVector(ex_idx[k], omega_sample[0], rhs[k]);
      const auto fp = ds.SourceFingerprint(rhs[k]);
      if (fp[0] != 0.0)
      {
        model.excitations.push_back(ex_idx[k]);
        model.exc_fp.insert(model.exc_fp.end(), fp.begin(), fp.end());
      }
    }
    model.ports = env_ports;
    model.omega = omega_sample;
    for (int p : env_ports)
    {
      port_l.push_back(VoltageFunctional(space_op, p, port_dofs.at(p)));
    }
    if (root)
    {
      model.WriteHeader(model_path);
    }
    Mpi::Print(" Saving the substructuring model to {} ({:d} frequencies, {:d} "
               "excitations with environment sources, {:d} environment ports)\n",
               model_path, omega_sample.size(), model.excitations.size(), env_ports.size());
  }
  else if (online)
  {
    model.ReadHeader(model_path, comm);
    MFEM_VERIFY(model.nG == nG,
                "Saved substructuring model interface size ("
                    << model.nG << ") does not match this run (" << nG
                    << "); the interface (Gamma) must be identical between the offline and "
                       "online runs.");
    gamma_map = MatchSignatureBasis(
        TrueDofSignatures(space_op.GetNDSpace().Get(), ds.InterfaceIndex(), nG),
        model.signatures, model.sig_w);
    MFEM_VERIFY(SameFingerprint(ds.EnvironmentFingerprint(model.omega[0]), model.env_fp, 2),
                "The environment differs from the one the saved substructuring model was "
                "condensed from (environment mesh, materials, boundary conditions, "
                "order or problem type changed): rerun in \"Offline\" mode to "
                "condense it again!");
    MFEM_VERIFY(env_ports == model.ports,
                "The environment's lumped ports differ from those of the saved "
                "substructuring model!");
    int nrec = root ? model.NumRecords(model_path) : 0;
    Mpi::Broadcast(1, &nrec, 0, comm);
    for (double omega : omega_sample)
    {
      const int j = FindValue(model.omega, omega);
      MFEM_VERIFY(
          j >= 0 && j < nrec,
          "Frequency " << iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(omega) /
                              (2 * std::numbers::pi)
                       << " GHz is not in the saved substructuring model: online driven "
                          "substructuring runs at saved frequencies only!");
      record.push_back(j);
    }
    // The model columns of the excitations with environment sources (-1: none).
    for (int k = 0; k < n; k++)
    {
      space_op.GetExcitationVector(ex_idx[k], model.omega[0], rhs[k]);
      const auto fp = ds.SourceFingerprint(rhs[k]);
      int col = -1;
      for (std::size_t e = 0; e < model.excitations.size(); e++)
      {
        if (model.excitations[e] == ex_idx[k] &&
            SameFingerprint(
                fp,
                {model.exc_fp.begin() + DrivenSubstructureModel::kSourceFp * e,
                 model.exc_fp.begin() + DrivenSubstructureModel::kSourceFp * (e + 1)},
                1))
        {
          col = static_cast<int>(e);
        }
      }
      MFEM_VERIFY(fp[0] == 0.0 || col >= 0,
                  "Excitation " << ex_idx[k]
                                << " has environment sources that differ from the saved "
                                   "substructuring model!");
      ex_col.push_back(col);
    }
  }
  const int np = static_cast<int>(env_ports.size()),
            ne = static_cast<int>(model.excitations.size());

  ComplexVector B(Curl.Height());
  B.UseDevice(true);
  auto t0 = Timer::Now();
  for (std::size_t omega_i = 0; omega_i < omega_sample.size(); omega_i++)
  {
    const double omega = omega_sample[omega_i];
    Mpi::Print("\nIt {:d}/{:d}: ω/2π = {:.3e} GHz (total elapsed time = {:.2e} s)\n",
               omega_i + 1, omega_sample.size(),
               iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(omega) /
                   (2 * std::numbers::pi),
               Timer::Duration(Timer::Now() - t0).count());
    std::vector<std::complex<double>> V(static_cast<std::size_t>(np) * n);  // port voltages
    {
      BlockTimer bt(Timer::KSP);
      for (int k = 0; k < n; k++)
      {
        space_op.GetExcitationVector(ex_idx[k], omega, rhs[k]);
      }
      if (!online)
      {
        ds.Condense(omega);
        ds.Solve(rhs_ptr, E);
      }
      else
      {
        // The saved record on the current interface (dual quantities: x_cur = M^-T x).
        DrivenSubstructureModel::Record rec;
        std::vector<std::complex<double>> S, g;
        if (root)
        {
          rec = model.ReadRecord(model_path, record[omega_i]);
          S = gamma_map.DualMatrix(rec.S);
          g.assign(static_cast<std::size_t>(nG) * n, 0.0);
          for (int k = 0; k < n; k++)
          {
            if (ex_col[k] >= 0)
            {
              const auto gk =
                  gamma_map.Dual(rec.g.data() + static_cast<std::size_t>(ex_col[k]) * nG);
              std::copy(gk.begin(), gk.end(),
                        g.begin() + static_cast<std::ptrdiff_t>(k) * nG);
            }
          }
          std::vector<std::complex<double>> h;
          for (int j = 0; j < np; j++)
          {
            const auto hj = gamma_map.Dual(rec.h.data() + static_cast<std::size_t>(j) * nG);
            h.insert(h.end(), hj.begin(), hj.end());
          }
          rec.h = std::move(h);
        }
        ds.Condense(omega, std::move(S));
        ds.Solve(rhs_ptr, E, &g);
        // The environment's port voltages, V_jk = c_jk + h_j^T u_Γ,k.
        if (root)
        {
          const auto &uG = ds.InterfaceSolution();
          for (int k = 0; k < n; k++)
          {
            for (int j = 0; j < np; j++)
            {
              std::complex<double> v =
                  (ex_col[k] >= 0) ? rec.c[static_cast<std::size_t>(ex_col[k]) * np + j]
                                   : 0.0;
              for (int a = 0; a < nG; a++)
              {
                v += rec.h[static_cast<std::size_t>(j) * nG + a] *
                     uG[static_cast<std::size_t>(k) * nG + a];
              }
              V[static_cast<std::size_t>(k) * np + j] = v;
            }
          }
        }
        Mpi::Broadcast(static_cast<int>(V.size()), V.data(), 0, comm);
      }
      if (save)
      {
        // The record of this frequency: the condensed voltage functionals of the
        // environment's ports, h_j = Reduce(l_j), and c_jk = l_j^T u_k - h_j^T u_Γ,k.
        std::vector<const ComplexVector *> L(np);
        for (int j = 0; j < np; j++)
        {
          L[j] = &port_l[j];
        }
        DrivenSubstructureModel::Record rec;
        rec.h = ds.CondenseEnvironment(L);
        std::vector<std::complex<double>> lu(static_cast<std::size_t>(np) * n);
        for (int k = 0; k < n; k++)
        {
          for (int j = 0; j < np; j++)
          {
            lu[static_cast<std::size_t>(k) * np + j] = {
                mfem::InnerProduct(port_l[j].Real(), E[k].Real()),
                mfem::InnerProduct(port_l[j].Real(), E[k].Imag())};
          }
        }
        Mpi::GlobalSum(static_cast<int>(lu.size()), lu.data(), comm);
        if (root)
        {
          rec.S = ds.Schur();
          const auto &gl = ds.EnvironmentSourceCondensation();
          const auto &uG = ds.InterfaceSolution();
          for (int e = 0; e < ne; e++)
          {
            const int k = static_cast<int>(std::ranges::find(ex_idx, model.excitations[e]) -
                                           ex_idx.begin());
            rec.g.insert(rec.g.end(), gl.begin() + static_cast<std::ptrdiff_t>(k) * nG,
                         gl.begin() + static_cast<std::ptrdiff_t>(k + 1) * nG);
            for (int j = 0; j < np; j++)
            {
              std::complex<double> c = lu[static_cast<std::size_t>(k) * np + j];
              for (int a = 0; a < nG; a++)
              {
                c -= rec.h[static_cast<std::size_t>(j) * nG + a] *
                     uG[static_cast<std::size_t>(k) * nG + a];
              }
              rec.c.push_back(c);
            }
          }
          model.AppendRecord(model_path, rec);
        }
      }
    }
    for (int k = 0; k < n; k++)
    {
      BlockTimer bt(Timer::POSTPRO);
      // B = -1/(iω) ∇ x E on the true dofs.
      Curl.Mult(E[k].Real(), B.Real());
      Curl.Mult(E[k].Imag(), B.Imag());
      B *= -1.0 / (1i * omega);
      if (online)
      {
        std::map<int, std::complex<double>> Vk;
        for (int j = 0; j < np; j++)
        {
          Vk[env_ports[j]] = V[static_cast<std::size_t>(k) * np + j];
        }
        post_op.SetLumpedPortVoltages(std::move(Vk));
      }
      post_op.MeasureAndPrintAll(ex_idx[k], int(omega_i), E[k], B, omega);
    }
  }
  // Substructuring computes no error estimate.
  ErrorIndicator indicator;
  post_op.MeasureFinalize(indicator);
  return indicator;
}

ErrorIndicator DrivenSolver::SweepUniform(SpaceOperator &space_op) const
{
  const auto &port_excitations = space_op.GetPortExcitations();
  const auto &omega_sample = iodata.solver.driven.sample_f;

  // Initialize postprocessing for measurement and printers.
  // Initialize write directory with default path; will be changed for multi-excitations.
  PostOperator<ProblemType::DRIVEN> post_op(iodata, space_op);

  // Construct the system matrices defining the linear operator. PEC boundaries are handled
  // simply by setting diagonal entries of the system matrix for the corresponding dofs.
  // Because the Dirichlet BC is always homogeneous, no special elimination is required on
  // the RHS. Assemble the linear system for the initial frequency (so we can call
  // KspSolver::SetOperators). Compute everything at the first frequency step.
  auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
  auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  const auto &Curl = space_op.GetCurlMatrix();

  // Set up the linear solver.
  // The operators are constructed for each frequency step and used to initialize the ksp.
  ComplexKspSolver ksp(iodata, space_op.GetNDSpaces(), &space_op.GetH1Spaces());

  // Set up RHS vector for the incident field at port boundaries, and the vector for the
  // first frequency step.
  ComplexVector RHS(Curl.Width()), E(Curl.Width()), B(Curl.Height());
  RHS.UseDevice(true);
  E.UseDevice(true);
  B.UseDevice(true);
  E = 0.0;
  B = 0.0;

  // Initialize structures for storing and reducing the results of error estimation.
  const bool is_2d = (space_op.GetNDSpace().Dimension() < 3);
  std::unique_ptr<TimeDependentFluxErrorEstimator<ComplexVector>> estimator_3d;
  std::unique_ptr<BoundaryModeFluxErrorEstimator<ComplexVector>> estimator_2d;
  if (is_2d)
  {
    estimator_2d = std::make_unique<BoundaryModeFluxErrorEstimator<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpaces(), space_op.GetRTSpaces(),
        space_op.GetCurlSpace(), space_op.GetH1Spaces(), iodata.solver.linear.estimator_tol,
        iodata.solver.linear.estimator_max_it, 0, iodata.solver.linear.estimator_mg);
  }
  else
  {
    estimator_3d = std::make_unique<TimeDependentFluxErrorEstimator<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpaces(), space_op.GetRTSpaces(),
        iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
        iodata.solver.linear.estimator_mg);
  }
  auto AddEstimate =
      [&](const ComplexVector &E, const ComplexVector &B, double Et, ErrorIndicator &ind)
  {
    if (is_2d)
      estimator_2d->AddErrorIndicator(E, B, Et, ind);
    else
      estimator_3d->AddErrorIndicator(E, B, Et, ind);
  };
  ErrorIndicator indicator;

  // If using Floquet BCs, a correction term (kp x E) needs to be added to the B field.
  std::unique_ptr<FloquetCorrSolver<ComplexVector>> floquet_corr;
  if (space_op.GetMaterialOp().HasWaveVector())
  {
    floquet_corr = std::make_unique<FloquetCorrSolver<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpace(), space_op.GetRTSpace(),
        iodata.solver.linear.tol, iodata.solver.linear.max_it, 0);
  }

  // Main excitation and frequency loop.
  auto t0 = Timer::Now();
  std::size_t excitation_counter = 0;
  const std::size_t excitation_restart_counter =
      ((iodata.solver.driven.restart - 1) / omega_sample.size()) + 1;
  const std::size_t freq_restart_idx =
      (iodata.solver.driven.restart - 1) % omega_sample.size();
  for (const auto &[excitation_idx, excitation_spec] : port_excitations)
  {
    if (++excitation_counter < excitation_restart_counter)
    {
      continue;
    }
    if (port_excitations.Size() > 1)
    {
      Mpi::Print("\nSweeping excitation index {:d} ({:d}/{:d}):\n", excitation_idx,
                 excitation_counter, port_excitations.Size());
    }
    // Frequency loop. Delay the potentially large visualization operators until the
    // first solution is ready, so they do not increase the solve/preconditioner peak.
    const std::size_t omega_start =
        (excitation_counter == excitation_restart_counter) ? freq_restart_idx : 0;
    for (std::size_t omega_i = omega_start; omega_i < omega_sample.size(); omega_i++)
    {
      auto omega = omega_sample[omega_i];
      // Assemble frequency dependent matrices and initialize operators in linear
      // solver.
      auto A2 = space_op.GetExtraSystemOperator(omega, Operator::DIAG_ZERO);
      auto A = space_op.GetSystemMatrix(1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i,
                                        K.get(), C.get(), M.get(), A2.get());
      auto P = space_op.GetPreconditionerMatrix<ComplexOperator>(
          1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i, omega);
      ksp.SetOperators(*A, *P);

      Mpi::Print(
          "\nIt {:d}/{:d}: ω/2π = {:.3e} GHz (total elapsed time = {:.2e} s{})\n",
          omega_i + 1, omega_sample.size(),
          iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(omega) /
              (2 * std::numbers::pi),
          Timer::Duration(Timer::Now() - t0).count(),
          (port_excitations.Size() > 1)
              ? fmt::format(", solve {:d}/{:d}",
                            1 + omega_i + (excitation_counter - 1) * omega_sample.size(),
                            omega_sample.size() * port_excitations.Size())
              : "");

      // Solve linear system.
      space_op.GetExcitationVector(excitation_idx, omega, RHS);

      Mpi::Print("\n");
      ksp.Mult(RHS, E);

      // Start Post-processing.
      BlockTimer bt0(Timer::POSTPRO);
      Mpi::Print(" Sol. ||E|| = {:.6e} (||RHS|| = {:.6e})\n",
                 linalg::Norml2(space_op.GetComm(), E),
                 linalg::Norml2(space_op.GetComm(), RHS));

      // Compute B = -1/(iω) ∇ x E on the true dofs.
      Curl.Mult(E.Real(), B.Real());
      Curl.Mult(E.Imag(), B.Imag());
      B *= -1.0 / (1i * omega);
      if (space_op.GetMaterialOp().HasWaveVector())
      {
        // Calculate B field correction for Floquet BCs: B += k_F(ω)/ω × E.
        // With k₀ = k_F_ref/ω_ref stored, k_F(ω)/ω = k_F_ref/ω_ref = k₀, so scale = 1.
        floquet_corr->AddMult(
            E, B,
            space_op.GetMaterialOp().HasFloquetFrequencyScaling() ? 1.0 : 1.0 / omega);
      }

      if (omega_i == omega_start)
      {
        // Switch ParaView subfolders once per excitation, after the first solve.
        post_op.InitializeParaviewDataCollection(excitation_idx);
      }
      auto total_domain_energy =
          post_op.MeasureAndPrintAll(excitation_idx, int(omega_i), E, B, omega);

      // Calculate and record the error indicators.
      Mpi::Print(" Updating solution error estimates\n");
      AddEstimate(E, B, total_domain_energy, indicator);
    }

    // Final postprocessing & printing.
    BlockTimer bt0(Timer::POSTPRO);
    SaveMetadata(ksp);
  }
  post_op.MeasureFinalize(indicator);
  return indicator;
}

ErrorIndicator DrivenSolver::SweepAdaptive(SpaceOperator &space_op) const
{
  const auto &port_excitations = space_op.GetPortExcitations();
  const auto &omega_sample = iodata.solver.driven.sample_f;
  // Initialize postprocessing for measurement and printers.
  // Initialize write directory with default path; will be changed for multi-excitations.
  PostOperator<ProblemType::DRIVEN> post_op(iodata, space_op);

  // Configure PROM parameters if not specified.
  double offline_tol = iodata.solver.driven.adaptive_tol;
  std::size_t convergence_memory = iodata.solver.driven.adaptive_memory;
  std::size_t max_size_per_excitation = iodata.solver.driven.adaptive_max_size;
  std::size_t nprom_indices = iodata.solver.driven.prom_indices.size();
  MFEM_VERIFY(max_size_per_excitation <= 0 || max_size_per_excitation >= nprom_indices,
              "Adaptive frequency sweep must sample at least " << nprom_indices
                                                               << " frequency points!");

  // Allocate negative curl matrix for postprocessing the B-field and vectors for the
  // high-dimensional field solution.
  const auto &Curl = space_op.GetCurlMatrix();
  ComplexVector E(Curl.Width()), Eh(Curl.Width()), B(Curl.Height());
  E.UseDevice(true);
  Eh.UseDevice(true);
  B.UseDevice(true);
  E = 0.0;
  Eh = 0.0;
  B = 0.0;

  // Initialize structures for storing and reducing the results of error estimation.
  const bool is_2d = (space_op.GetNDSpace().Dimension() < 3);
  std::unique_ptr<TimeDependentFluxErrorEstimator<ComplexVector>> estimator_3d;
  std::unique_ptr<BoundaryModeFluxErrorEstimator<ComplexVector>> estimator_2d;
  if (is_2d)
  {
    estimator_2d = std::make_unique<BoundaryModeFluxErrorEstimator<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpaces(), space_op.GetRTSpaces(),
        space_op.GetCurlSpace(), space_op.GetH1Spaces(), iodata.solver.linear.estimator_tol,
        iodata.solver.linear.estimator_max_it, 0, iodata.solver.linear.estimator_mg);
  }
  else
  {
    estimator_3d = std::make_unique<TimeDependentFluxErrorEstimator<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpaces(), space_op.GetRTSpaces(),
        iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
        iodata.solver.linear.estimator_mg);
  }
  auto AddEstimate =
      [&](const ComplexVector &E, const ComplexVector &B, double Et, ErrorIndicator &ind)
  {
    if (is_2d)
      estimator_2d->AddErrorIndicator(E, B, Et, ind);
    else
      estimator_3d->AddErrorIndicator(E, B, Et, ind);
  };
  ErrorIndicator indicator;

  // If using Floquet BCs, a correction term (kp x E) needs to be added to the B field.
  std::unique_ptr<FloquetCorrSolver<ComplexVector>> floquet_corr;
  if (space_op.GetMaterialOp().HasWaveVector())
  {
    floquet_corr = std::make_unique<FloquetCorrSolver<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpace(), space_op.GetRTSpace(),
        iodata.solver.linear.tol, iodata.solver.linear.max_it, 0);
  }

  // Configure the PROM operator which performs the parameter space sampling and basis
  // construction during the offline phase as well as the PROM solution during the online
  // phase.
  auto t0 = Timer::Now();
  const double unit_GHz = iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(1.0) /
                          (2 * std::numbers::pi);
  Mpi::Print("\nBeginning PROM construction offline phase:\n"
             " {:d} points for frequency sweep over [{:.3e}, {:.3e}] GHz\n",
             omega_sample.size(), omega_sample.front() * unit_GHz,
             omega_sample.back() * unit_GHz);
  RomOperator prom_op(iodata, space_op, max_size_per_excitation);
  space_op.GetWavePortOp().SetSuppressOutput(true);
  space_op.GetWavePortOp().ConfigureReducedModelTraining(
      max_size_per_excitation, port_excitations.Size(),
      iodata.solver.driven.adaptive_circuit_synthesis);

  // Add ports to PROM if we do synthesis.
  if (iodata.solver.driven.adaptive_circuit_synthesis)
  {
    prom_op.AddLumpedPortModesForSynthesis();
    if (space_op.GetWavePortOp().Size() > 0)
    {
      // Use the band center as the reference frequency for seeding wave-port modes.
      // The choice rescales the basis vector but does not change correctness.
      const double omega_ref = 0.5 * (omega_sample.front() + omega_sample.back());
      prom_op.AddWavePortModesForSynthesis(omega_ref);
    }
  }

  // Initialize the basis with samples from the top and bottom of the frequency
  // range of interest. Each call for an HDM solution adds the frequency sample to P_S and
  // removes it from P \ P_S. Timing for the HDM construction and solve is handled inside
  // of the RomOperator.
  auto UpdatePROM = [&](int excitation_idx, double omega, std::size_t sample_idx)
  {
    // Add the HDM solution to the PROM reduced basis.
    prom_op.UpdatePROM(E, fmt::format("sample_e{:d}_s{:d}", excitation_idx, sample_idx));
    prom_op.UpdateMRI(excitation_idx, omega, E);

    // Compute B = -1/(iω) ∇ x E on the true dofs, and set the internal GridFunctions in
    // PostOperator for energy postprocessing and error estimation.
    BlockTimer bt0(Timer::POSTPRO);
    Curl.Mult(E.Real(), B.Real());
    Curl.Mult(E.Imag(), B.Imag());
    B *= -1.0 / (1i * omega);
    if (space_op.GetMaterialOp().HasWaveVector())
    {
      // Calculate B field correction for Floquet BCs: B += k_F(ω)/ω × E.
      // With k₀ = k_F_ref/ω_ref stored, k_F(ω)/ω = k₀, so scale = 1.
      floquet_corr->AddMult(
          E, B, space_op.GetMaterialOp().HasFloquetFrequencyScaling() ? 1.0 : 1.0 / omega);
    }

    // Measure domain energies for the error indicator only. Don't exchange face_nbr_data,
    // unless printing paraview fields.
    auto total_domain_energy = post_op.MeasureDomainFieldEnergyOnly(E, B);
    AddEstimate(E, B, total_domain_energy, indicator);
  };

  // Loop excitations to add to PROM.
  //
  // Restart should not really be used for adaptive sweeps, but must work. Construct PROM in
  // the same way same regardless of restart for consistency. Don't shift excitation start.
  int excitation_counter = 0;
  for (const auto &[excitation_idx, excitation_spec] : port_excitations)
  {
    if (port_excitations.Size() > 1)
    {
      Mpi::Print("\nAdding excitation index {:d} ({:d}/{:d}):\n", excitation_idx,
                 ++excitation_counter, port_excitations.Size());
    }
    prom_op.SetExcitationIndex(excitation_idx);  // Pre-compute RHS1

    // Initialize PROM with explicit HDM samples, record the estimate but do not act on it.
    std::vector<double> max_errors;
    std::size_t counter_rom_sample = 0;
    for (auto i : iodata.solver.driven.prom_indices)
    {
      auto omega = omega_sample[i];
      prom_op.SolveHDM(excitation_idx, omega, E);
      prom_op.SolvePROM(excitation_idx, omega, Eh);
      linalg::AXPY(-1.0, E, Eh);
      max_errors.push_back(linalg::Norml2(space_op.GetComm(), Eh) /
                           linalg::Norml2(space_op.GetComm(), E));
      UpdatePROM(excitation_idx, omega, counter_rom_sample);
      counter_rom_sample++;
    }
    // The estimates associated to the end points are assumed inaccurate.
    max_errors[0] = std::numeric_limits<double>::infinity();
    max_errors[1] = std::numeric_limits<double>::infinity();
    auto memory = std::distance(max_errors.rbegin(),
                                std::find_if(max_errors.rbegin(), max_errors.rend(),
                                             [=](auto x) { return x > offline_tol; }));
    memory = std::max(0L, memory);  // Ensure memory >= 0 as it should be.

    // Greedy procedure for basis construction (offline phase). Basis is initialized with
    // solutions at frequency sweep endpoints and explicit sample frequencies.
    std::size_t it = max_errors.size();
    for (std::size_t it0 = it; it < max_size_per_excitation && memory < convergence_memory;
         it++)
    {
      // Compute the location of the maximum error in parameter domain (bounded by the
      // previous samples).
      double omega_star = prom_op.FindMaxError(excitation_idx)[0];

      // Sample HDM and add solution to basis.
      prom_op.SolveHDM(excitation_idx, omega_star, E);
      prom_op.SolvePROM(excitation_idx, omega_star, Eh);
      linalg::AXPY(-1.0, E, Eh);

      max_errors.push_back(linalg::Norml2(space_op.GetComm(), Eh) /
                           linalg::Norml2(space_op.GetComm(), E));
      memory = max_errors.back() < offline_tol ? memory + 1 : 0;

      Mpi::Print("\nGreedy iteration {:d} (n = {:d}): ω* = {:.3e} GHz ({:.3e}), error = "
                 "{:.3e}, memory = {:d}/{:d}\n",
                 it - it0 + 1, prom_op.GetReducedDimension(), omega_star * unit_GHz,
                 omega_star, max_errors.back(), memory, convergence_memory);
      UpdatePROM(excitation_idx, omega_star, counter_rom_sample);
      counter_rom_sample++;
    }
    Mpi::Print("\nAdaptive sampling{} {:d} frequency samples:\n"
               " n = {:d}, error = {:.3e}, tol = {:.3e}, memory = {:d}/{:d}\n",
               (it == max_size_per_excitation) ? " reached maximum" : " converged with", it,
               prom_op.GetReducedDimension(), max_errors.back(), offline_tol, memory,
               convergence_memory);
    utils::PrettyPrint(prom_op.GetSamplePoints(excitation_idx), unit_GHz,
                       " Sampled frequencies (GHz):");
    utils::PrettyPrint(max_errors, 1.0, " Sample errors:");
  }

  Mpi::Print(" Total offline phase elapsed time: {:.2e} s\n",
             Timer::Duration(Timer::Now() - t0).count());  // Timing on root

  // Circuit synthesis samples the exact cross-section EVP at a tightened tolerance so the
  // fitted pencil does not depend on the port-mode accuracy floor; keep those solves exact
  // (they also enrich the per-port reduced basis).
  if (iodata.solver.driven.adaptive_circuit_synthesis)
  {
    prom_op.PrintPROMMatrices(iodata.units, iodata.problem.output);
  }

  // Exact wave-port modes computed by the HDM samples (and any synthesis samples) have
  // trained a separate per-port reduced eigenspace. Enable guarded Rayleigh-Ritz evaluation
  // for the online output sweep; failed residual checks transparently fall back to the
  // exact port eigensolver and enrich the basis.
  space_op.GetWavePortOp().EnableReducedModel(iodata.solver.driven.adaptive_tol);
  prom_op.PrepareOnlineExcitations();
  post_op.ConfigureReducedPostprocessing(prom_op);

  // Main fast frequency sweep loop (online phase).
  Mpi::Print("\nBeginning fast frequency sweep online phase\n");
  space_op.GetWavePortOp().SetSuppressOutput(false);  // Disable output suppression
  auto solve_online_point = [&](int excitation_idx, std::size_t omega_i)
  {
    const auto omega = omega_sample[omega_i];
    // Refresh every port, including inactive observation-only ports, before either PROM
    // assembly or postprocessing. Repeated calls at the same frequency are cached.
    space_op.GetWavePortOp().PrepareFrequency(omega);
    Mpi::Print("\nIt {:d}/{:d}, excitation {:d}: ω/2π = {:.3e} GHz "
               "(total elapsed time = {:.2e} s)\n",
               omega_i + 1, omega_sample.size(), excitation_idx,
               iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(omega) /
                   (2 * std::numbers::pi),
               Timer::Duration(Timer::Now() - t0).count());

    // Assemble and solve the PROM linear system.
    prom_op.SolvePROM(excitation_idx, omega, E);
    Mpi::Print("\n");
    if (omega_i == 0 && post_op.WillWriteFields())
    {
      // Switch ParaView subfolders once per excitation. Delay visualization setup until
      // the first online solution is ready.
      post_op.InitializeParaviewDataCollection(excitation_idx);
    }

    if (post_op.HasReducedPostprocessing() &&
        !post_op.WillWriteFields(static_cast<int>(omega_i)))
    {
      post_op.MeasureAndPrintReduced(excitation_idx, int(omega_i), E, omega,
                                     prom_op.GetReducedSolution());
      return;
    }

    // Start full post-processing.
    BlockTimer bt0(Timer::POSTPRO);
    Mpi::Print(" Sol. ||E|| = {:.6e}\n", linalg::Norml2(space_op.GetComm(), E));

    // Compute B = -1/(iω) ∇ x E on the true dofs.
    Curl.Mult(E.Real(), B.Real());
    Curl.Mult(E.Imag(), B.Imag());
    B *= -1.0 / (1i * omega);
    if (space_op.GetMaterialOp().HasWaveVector())
    {
      // Calculate B field correction for Floquet BCs: B += k_F(ω)/ω × E.
      // With k₀ = k_F_ref/ω_ref stored, k_F(ω)/ω = k_F_ref/ω_ref = k₀, so scale = 1.
      floquet_corr->AddMult(
          E, B, space_op.GetMaterialOp().HasFloquetFrequencyScaling() ? 1.0 : 1.0 / omega);
    }
    post_op.MeasureAndPrintAll(excitation_idx, int(omega_i), E, B, omega);
  };

  if (!post_op.WillWriteFields())
  {
    // Port modes depend on frequency, not excitation. Frequency-major traversal keeps the
    // one-frequency modal cache hot for every excitation.
    for (std::size_t omega_i = 0; omega_i < omega_sample.size(); omega_i++)
    {
      for (const auto &[excitation_idx, excitation_spec] : port_excitations)
      {
        solve_online_point(excitation_idx, omega_i);
      }
    }
  }
  else
  {
    // Preserve excitation-major ordering for excitation-specific field collections.
    for (const auto &[excitation_idx, excitation_spec] : port_excitations)
    {
      if (port_excitations.Size() > 1)
      {
        Mpi::Print("\nSweeping excitation index {:d}:\n", excitation_idx);
      }
      for (std::size_t omega_i = 0; omega_i < omega_sample.size(); omega_i++)
      {
        solve_online_point(excitation_idx, omega_i);
      }
    }
  }

  // Final postprocessing & printing: no change to indicator since these are in PROM.
  {
    BlockTimer bt0(Timer::POSTPRO);
    SaveMetadata(prom_op.GetLinearSolver());
  }
  space_op.GetWavePortOp().PrintReducedModelStats();
  post_op.MeasureFinalize(indicator);
  return indicator;
}

}  // namespace palace
