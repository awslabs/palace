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
    return {SweepSubstructured(space_op, adaptive), space_op.GlobalTrueVSize()};
  }
  return {adaptive ? SweepAdaptive(space_op) : SweepUniform(space_op),
          space_op.GlobalTrueVSize()};
}

namespace
{

// The lumped ports on the environment side (and its wave ports, as -index), and an abort
// for lumped ports touching Γ (their excitation and voltage would depend on both sides).
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
  for (int idx : ds.EnvironmentWavePorts())
  {
    env_ports.push_back(-idx);
  }
  return env_ports;
}

// The functional of a wave port's overlap with its mode at ω as a true-DOF vector l (the
// overlap is l^T E): -conj(s), with s the port's n×H mode vector.
ComplexVector WaveOverlapFunctional(SpaceOperator &space_op, int idx, double omega)
{
  auto s = space_op.GetWavePortModeVector(idx, omega);
  ComplexVector l(s->Size());
  l.Real() = s->Real();
  l.Real() *= -1.0;
  l.Imag() = s->Imag();
  return l;
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

// The full S_E (n x n, column-major) of a record's lower triangle (by columns), mapped to
// the current interface as M^-T S_E M^-1 with map (rows[g]: row g of M^-T), on one
// triangle.
std::vector<std::complex<double>> FullSchur(const std::complex<double> *lower, int n,
                                            const SignatureMap *map)
{
  auto at = [&](int a, int b)
  {
    if (a < b)
    {
      std::swap(a, b);
    }
    return lower[static_cast<std::size_t>(b) * n -
                 static_cast<std::size_t>(b) * (b - 1) / 2 + (a - b)];
  };
  std::vector<std::complex<double>> S(static_cast<std::size_t>(n) * n);
  for (int j = 0; j < n; j++)
  {
    for (int i = j; i < n; i++)
    {
      std::complex<double> v = 0.0;
      if (map)
      {
        for (const auto &[a, ca] : map->rows[i])
        {
          for (const auto &[b, cb] : map->rows[j])
          {
            v += ca * cb * at(a, b);
          }
        }
      }
      else
      {
        v = at(i, j);
      }
      S[static_cast<std::size_t>(j) * n + i] = S[static_cast<std::size_t>(i) * n + j] = v;
    }
  }
  return S;
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

}  // namespace

ErrorIndicator DrivenSolver::SweepSubstructured(SpaceOperator &space_op,
                                                bool adaptive) const
{
  // Exact per-frequency substructuring: at each frequency, condense the environment and
  // factor the region against it once, then solve all excitations (so the frequency loop
  // is the outer one). An adaptive sweep condenses the environment at frequencies chosen
  // greedily for a rational model of its condensed data (BarycentricInterpolant), then
  // solves every frequency on the region against the model. With a saved model, an offline
  // sweep saves the condensed data per frequency, and an online sweep factors the region
  // only: at the saved frequencies, or anywhere in the band of a rational model.
  const auto &port_excitations = space_op.GetPortExcitations();
  const auto &omega_sample = iodata.solver.driven.sample_f;
  const auto &sub = *iodata.solver.substructuring;
  const bool online = (sub.mode == SubstructuringMode::ONLINE);
  adaptive = adaptive && !online;
  const bool from_model = online || adaptive;
  const std::string &model_path = sub.save_model;
  const bool save = !online && !model_path.empty();
  MPI_Comm comm = space_op.GetComm();
  const bool root = Mpi::Root(comm);
  auto ds = std::make_unique<DrivenSubstructure>(space_op, sub.region_attributes,
                                                 sub.environment_attributes, online);
  const int nG = ds->InterfaceSize();
  Mpi::Print("\nSubstructuring: |Γ| = {:d} interface unknowns\n", nG);
  const auto &Curl = space_op.GetCurlMatrix();
  const double unit_GHz = iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(1.0) /
                          (2 * std::numbers::pi);
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

  // The lumped ports of the environment (for a model).
  std::map<int, std::vector<char>> port_dofs;
  std::vector<int> env_ports;
  if (save || from_model)
  {
    env_ports = EnvironmentPorts(iodata, *ds, port_dofs);
  }

  // From a model, the field is known in the region (and on Γ) only: energies of region
  // domains, without totals or participation ratios, and the environment's port voltages
  // from the model (without their power).
  config::DomainData domains = iodata.domains;
  if (from_model)
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
    Mpi::Warning("{} driven substructuring writes no total energies or participation "
                 "ratios{}{}!\n",
                 online ? "Online" : "Adaptive",
                 skipped.empty() ? std::string()
                                 : fmt::format(", no energies of the domain postprocessing "
                                               "outside the region ({})",
                                               fmt::join(skipped, ", ")),
                 env_ports.empty() ? std::string()
                                   : ", and no power (or wave port voltage) of the "
                                     "environment's ports");
  }
  PostOperator<ProblemType::DRIVEN> post_op(iodata.problem, iodata.solver, domains,
                                            iodata.boundaries, iodata.units, space_op);
  if (from_model && post_op.WillWriteFields())
  {
    Mpi::Warning("{} driven substructuring writes the fields in the region and on the "
                 "interface only (zero in the environment)!\n",
                 online ? "Online" : "Adaptive");
  }

  // The model: header and environment ports.
  DrivenSubstructureModel model;
  std::vector<ComplexVector> port_l;
  SignatureMap gamma_map;           // online: the saved interface on the current one
  std::vector<int> record, ex_col;  // the records of an exact model, model columns
  std::vector<std::vector<std::complex<double>>> records;  // a rational model (root)
  if (save || adaptive)
  {
    model.nG = nG;
    model.sig_w = SignatureWidth(space_op.GetNDSpace().Get());
    model.signatures =
        TrueDofSignatures(space_op.GetNDSpace().Get(), ds->InterfaceIndex(), nG);
    model.env_fp = ds->EnvironmentFingerprint(omega_sample[0]);
    for (int k = 0; k < n; k++)
    {
      space_op.GetExcitationVector(ex_idx[k], omega_sample[0], rhs[k]);
      const auto fp = ds->SourceFingerprint(rhs[k]);
      ex_col.push_back(-1);
      if (fp[0] != 0.0)
      {
        ex_col.back() = static_cast<int>(model.excitations.size());
        model.excitations.push_back(ex_idx[k]);
        model.exc_fp.insert(model.exc_fp.end(), fp.begin(), fp.end());
      }
    }
    model.ports = env_ports;
    if (adaptive)
    {
      const std::size_t n_init = iodata.solver.driven.prom_indices.size();
      MFEM_VERIFY(iodata.solver.driven.adaptive_max_size >= n_init,
                  "Adaptive frequency sweep must sample at least " << n_init
                                                                   << " frequency points!");
      model.capacity = static_cast<int>(iodata.solver.driven.adaptive_max_size);
    }
    else
    {
      model.omega = omega_sample;
      model.capacity = static_cast<int>(omega_sample.size());
    }
    for (int p : env_ports)
    {
      port_l.push_back((p > 0) ? VoltageFunctional(space_op, p, port_dofs.at(p))
                               : ComplexVector());
    }
    if (save && root)
    {
      model.WriteHeader(model_path);
    }
    if (save)
    {
      Mpi::Print(" Saving the substructuring model to {} ({}{:d} frequencies, {:d} "
                 "excitations with environment sources, {:d} environment ports)\n",
                 model_path, adaptive ? "up to " : "", model.capacity,
                 model.excitations.size(), env_ports.size());
    }
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
        TrueDofSignatures(space_op.GetNDSpace().Get(), ds->InterfaceIndex(), nG),
        model.signatures, model.sig_w);
    MFEM_VERIFY(DrivenSubstructureModel::SameEnvironment(
                    ds->EnvironmentFingerprint(model.omega[0]), model.env_fp),
                "The environment differs from the one the saved substructuring model was "
                "condensed from (environment mesh, materials, boundary conditions, "
                "order or problem type changed): rerun in \"Offline\" mode to "
                "condense it again!");
    MFEM_VERIFY(env_ports == model.ports,
                "The environment's lumped ports differ from those of the saved "
                "substructuring model!");
    int nrec = root ? model.NumRecords(model_path) : 0;
    Mpi::Broadcast(1, &nrec, 0, comm);
    MFEM_VERIFY(nrec >= static_cast<int>(model.omega.size()),
                "Truncated substructuring model \"" << model_path << "\"!");
    if (!model.weights.empty())
    {
      // A rational model: any frequency between its samples.
      const auto [lo, hi] = std::ranges::minmax(model.omega);
      for (double omega : omega_sample)
      {
        MFEM_VERIFY(omega >= lo * (1.0 - 1.0e-10) && omega <= hi * (1.0 + 1.0e-10),
                    fmt::format("Frequency {:.6g} GHz is outside the band of the saved "
                                "substructuring model ({:.6g} to {:.6g} GHz)!",
                                omega * unit_GHz, lo * unit_GHz, hi * unit_GHz));
      }
      for (std::size_t j = 0; root && j < model.omega.size(); j++)
      {
        records.push_back(model.ReadFlatRecord(model_path, static_cast<int>(j)));
      }
      Mpi::Print(" Rational substructuring model: {:d} samples from {:.3e} to {:.3e} GHz\n",
                 model.omega.size(), lo * unit_GHz, hi * unit_GHz);
    }
    else
    {
      for (double omega : omega_sample)
      {
        const int j = FindValue(model.omega, omega);
        MFEM_VERIFY(j >= 0, fmt::format("Frequency {:.6g} GHz is not in the saved "
                                        "substructuring model: online driven "
                                        "substructuring runs at saved frequencies only "
                                        "(or anywhere in the band of a model saved by an "
                                        "adaptive sweep)!",
                                        omega * unit_GHz));
        record.push_back(j);
      }
    }
    // The model columns of the excitations with environment sources (-1: none). The modes
    // of the environment's wave ports, and so their sources, are reproducible to their
    // eigensolver tolerance over the spectral gap only.
    for (int k = 0; k < n; k++)
    {
      space_op.GetExcitationVector(ex_idx[k], model.omega[0], rhs[k]);
      const auto fp = ds->SourceFingerprint(rhs[k]);
      double tol = 1.0e-10;
      for (int p : port_excitations.excitations.at(ex_idx[k]).wave_port)
      {
        if (std::ranges::find(env_ports, -p) != env_ports.end())
        {
          tol = std::max(tol, 100.0 * iodata.boundaries.waveport.at(p).eig_tol);
        }
      }
      int col = -1;
      for (std::size_t e = 0; e < model.excitations.size(); e++)
      {
        if (model.excitations[e] == ex_idx[k] &&
            DrivenSubstructureModel::SameSource(
                fp.data(), model.exc_fp.data() + DrivenSubstructureModel::kSourceFp * e,
                tol))
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

  // The functionals l_j of the environment's ports at ω.
  auto port_functionals = [&](double omega)
  {
    std::vector<const ComplexVector *> L(np);
    for (int j = 0; j < np; j++)
    {
      if (env_ports[j] < 0)
      {
        port_l[j] = WaveOverlapFunctional(space_op, -env_ports[j], omega);
      }
      L[j] = &port_l[j];
    }
    return L;
  };
  // Offline: condense the environment at ω and solve all excitations, and the record of the
  // model (on rank 0), with the condensed functionals of the environment's ports,
  // h_j = Reduce(l_j), and c_jk = l_j^T u_k - h_j^T u_Γ,k.
  auto condense = [&](double omega, bool with_record)
  {
    for (int k = 0; k < n; k++)
    {
      space_op.GetExcitationVector(ex_idx[k], omega, rhs[k]);
    }
    ds->Condense(omega);
    ds->Solve(rhs_ptr, E);
    DrivenSubstructureModel::Record rec;
    if (!with_record)
    {
      return rec;
    }
    rec.h = ds->CondenseEnvironment(port_functionals(omega));
    std::vector<std::complex<double>> lu(static_cast<std::size_t>(np) * n);
    for (int k = 0; k < n; k++)
    {
      for (int j = 0; j < np; j++)
      {
        const auto &l = port_l[j];
        lu[static_cast<std::size_t>(k) * np + j] = {
            mfem::InnerProduct(l.Real(), E[k].Real()) -
                mfem::InnerProduct(l.Imag(), E[k].Imag()),
            mfem::InnerProduct(l.Real(), E[k].Imag()) +
                mfem::InnerProduct(l.Imag(), E[k].Real())};
      }
    }
    Mpi::GlobalSum(static_cast<int>(lu.size()), lu.data(), comm);
    if (root)
    {
      rec.S = ds->TakeSchur();
      const auto &gl = ds->EnvironmentSourceCondensation();
      const auto &uG = ds->InterfaceSolution();
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
    }
    return rec;
  };

  // The record alone, from the environment (not factoring the region): with the condensed
  // sources g_k and interior solutions x_k = A_EE^-1 b_E,k (Γ held at 0) of the excitations
  // with environment sources, c_jk = l_j^T x_k, since l_j^T u_k = l_j^T x_k + h_j^T u_Γ,k.
  auto condense_environment = [&](double omega)
  {
    std::vector<const ComplexVector *> B(ne);
    for (int e = 0; e < ne; e++)
    {
      const int k = static_cast<int>(std::ranges::find(ex_idx, model.excitations[e]) -
                                     ex_idx.begin());
      space_op.GetExcitationVector(ex_idx[k], omega, rhs[k]);
      B[e] = &rhs[k];
    }
    ds->Condense(omega, false);
    DrivenSubstructureModel::Record rec;
    rec.h = ds->CondenseEnvironment(port_functionals(omega));
    std::vector<ComplexVector> x;
    auto g = ds->CondenseEnvironment(B, &x);
    std::vector<std::complex<double>> c(static_cast<std::size_t>(np) * ne);
    for (int e = 0; e < ne; e++)
    {
      for (int j = 0; j < np; j++)
      {
        const auto &l = port_l[j];
        c[static_cast<std::size_t>(e) * np + j] = {
            mfem::InnerProduct(l.Real(), x[e].Real()) -
                mfem::InnerProduct(l.Imag(), x[e].Imag()),
            mfem::InnerProduct(l.Real(), x[e].Imag()) +
                mfem::InnerProduct(l.Imag(), x[e].Real())};
      }
    }
    Mpi::GlobalSum(static_cast<int>(c.size()), c.data(), comm);
    if (root)
    {
      rec.S = ds->TakeSchur();
      rec.g = std::move(g);
      rec.c = std::move(c);
    }
    return rec;
  };

  // Adaptive: the samples of the rational model, chosen greedily as the adaptive driven
  // solver chooses its own, with the records scaled part by part (by their norms at the
  // first sample) so that the error measure weighs them alike.
  BarycentricInterpolant fit;
  std::vector<double> part_scale;
  auto scale_record = [&](std::vector<std::complex<double>> &x, bool inverse)
  {
    const auto sz = model.PartSizes();
    std::size_t off = 0;
    if (part_scale.empty())
    {
      for (std::size_t p = 0; p < sz.size(); p++)
      {
        double nrm = 0.0;
        for (std::size_t q = off; q < off + sz[p]; q++)
        {
          nrm += std::norm(x[q]);
        }
        part_scale.push_back((nrm > 0.0) ? 1.0 / std::sqrt(nrm) : 1.0);
        off += sz[p];
      }
      off = 0;
    }
    for (std::size_t p = 0; p < sz.size(); p++)
    {
      const double c = inverse ? 1.0 / part_scale[p] : part_scale[p];
      for (std::size_t q = off; q < off + sz[p]; q++)
      {
        x[q] *= c;
      }
      off += sz[p];
    }
  };
  if (adaptive)
  {
    BlockTimer bt(Timer::KSP);
    const double tol = iodata.solver.driven.adaptive_tol;
    const std::size_t max_samples = iodata.solver.driven.adaptive_max_size,
                      memory_size = iodata.solver.driven.adaptive_memory;
    Mpi::Print("\nAdaptive sampling of the condensed environment (tol = {:.3e}, at most "
               "{:d} frequencies):\n",
               tol, max_samples);
    // The relative error of the rational model at a new sample, before adding it.
    auto sample = [&](double omega)
    {
      auto rec = condense_environment(omega);
      double err = std::numeric_limits<double>::infinity();
      if (root)
      {
        if (save)
        {
          model.AppendRecord(model_path, rec);
        }
        auto x = model.Flatten(rec);
        scale_record(x, false);
        if (fit.Samples().size() >= 2)
        {
          const auto y = fit.Evaluate(omega);
          double d = 0.0, m = 0.0;
          for (std::size_t q = 0; q < x.size(); q++)
          {
            d += std::norm(y[q] - x[q]);
            m += std::norm(x[q]);
          }
          err = std::sqrt(d / m);
        }
        fit.AddSample(omega, x);
        model.weights = fit.Weights();
      }
      Mpi::Broadcast(1, &err, 0, comm);
      model.omega.push_back(omega);
      if (save && root)
      {
        model.WriteHeader(model_path, false);
      }
      return err;
    };
    std::vector<double> errors;
    for (auto i : iodata.solver.driven.prom_indices)
    {
      errors.push_back(sample(omega_sample[i]));
    }
    // The errors of the end points are not estimates.
    errors[0] = errors[1] = std::numeric_limits<double>::infinity();
    std::size_t memory = 0;
    for (auto it = errors.rbegin(); it != errors.rend() && *it < tol; ++it)
    {
      memory++;
    }
    while (model.omega.size() < max_samples && memory < memory_size)
    {
      double omega = root ? fit.FindMaxError() : 0.0;
      Mpi::Broadcast(1, &omega, 0, comm);
      errors.push_back(sample(omega));
      memory = (errors.back() < tol) ? memory + 1 : 0;
      Mpi::Print(" Greedy iteration {:d}: ω* = {:.6e} GHz, error = {:.3e}, memory = "
                 "{:d}/{:d}\n",
                 model.omega.size() - iodata.solver.driven.prom_indices.size(),
                 omega * unit_GHz, errors.back(), memory, memory_size);
    }
    Mpi::Print(
        "\nAdaptive sampling{} {:d} frequency samples: error = {:.3e}, tol = {:.3e}, "
        "memory = {:d}/{:d}\n",
        (memory < memory_size) ? " reached maximum" : " converged with", model.omega.size(),
        errors.back(), tol, memory, memory_size);
    // Passivity between the samples (Im S_E positive semidefinite, §10.4): a violation
    // beyond the tolerance adds a sample where it is largest.
    double passivity = 0.0;
    while (true)
    {
      double at = 0.0;
      passivity = std::numeric_limits<double>::infinity();
      if (root)
      {
        auto z = model.omega;
        std::ranges::sort(z);
        for (std::size_t j = 0; j + 1 < z.size(); j++)
        {
          const double mid = 0.5 * (z[j] + z[j + 1]);
          const double p = DrivenSubstructureModel::Passivity(fit.Evaluate(mid).data(), nG);
          if (p < passivity)
          {
            passivity = p;
            at = mid;
          }
        }
      }
      Mpi::Broadcast(1, &passivity, 0, comm);
      Mpi::Broadcast(1, &at, 0, comm);
      if (passivity >= -tol || model.omega.size() >= max_samples)
      {
        break;
      }
      errors.push_back(sample(at));
      Mpi::Print(" Passivity sample: ω = {:.6e} GHz (least eigenvalue of Im S_E / ‖S_E‖ = "
                 "{:.3e}), error = {:.3e}\n",
                 at * unit_GHz, passivity, errors.back());
    }
    Mpi::Print(
        " Passivity between the samples: least eigenvalue of Im S_E / ‖S_E‖ = {:.3e}\n",
        passivity);
    if (passivity < -tol)
    {
      Mpi::Warning("The rational environment model is not passive between its samples "
                   "(least eigenvalue of Im S_E / ‖S_E‖ = {:.3e}): raise "
                   "\"AdaptiveMaxSamples\"!\n",
                   passivity);
    }
    utils::PrettyPrint(model.omega, unit_GHz, " Sampled frequencies (GHz):");
    utils::PrettyPrint(errors, 1.0, " Sample errors:");

    // The frequencies of the sweep from the model, on the region only.
    ds.reset();
    ds = std::make_unique<DrivenSubstructure>(space_op, sub.region_attributes,
                                              sub.environment_attributes, true);
  }

  // The model's record at a frequency of the sweep, on the current interface (on rank 0;
  // online, dual quantities transform as x_cur = M^-T x): S_E, the condensed environment
  // sources g of the excitations (|Γ| x n), and the environment ports' h and c.
  auto model_record = [&](double omega, std::size_t omega_i)
  {
    std::vector<std::complex<double>> x;  // the record, flat (S_E by its lower triangle)
    if (adaptive)
    {
      x = fit.Evaluate(omega);
      scale_record(x, true);
    }
    else if (!model.weights.empty())
    {
      const auto a =
          BarycentricInterpolant::Coefficients(model.omega, model.weights, omega);
      x.assign(records[0].size(), 0.0);
      for (std::size_t j = 0; j < records.size(); j++)
      {
        for (std::size_t q = 0; q < x.size(); q++)
        {
          x[q] += a[j] * records[j][q];
        }
      }
    }
    else
    {
      x = model.ReadFlatRecord(model_path, record[omega_i]);
    }
    const auto sz = model.PartSizes();
    DrivenSubstructureModel::Record out;
    out.S = FullSchur(x.data(), nG, online ? &gamma_map : nullptr);
    const auto *g = x.data() + sz[0], *h = g + sz[1], *c = h + sz[2];
    auto dual = [&](const std::complex<double> *y)
    {
      std::vector<std::complex<double>> z(y, y + nG);
      return online ? gamma_map.DualRows(z, 1) : z;
    };
    out.g.assign(static_cast<std::size_t>(nG) * n, 0.0);
    for (int k = 0; k < n; k++)
    {
      if (ex_col[k] >= 0)
      {
        const auto gk = dual(g + static_cast<std::size_t>(ex_col[k]) * nG);
        std::copy(gk.begin(), gk.end(),
                  out.g.begin() + static_cast<std::ptrdiff_t>(k) * nG);
      }
    }
    for (int j = 0; j < np; j++)
    {
      const auto hj = dual(h + static_cast<std::size_t>(j) * nG);
      out.h.insert(out.h.end(), hj.begin(), hj.end());
    }
    out.c.assign(c, c + sz[3]);
    return out;
  };

  // From the model at ω: the region condensed against S_E, and all excitations solved.
  auto solve_region = [&](double omega, DrivenSubstructureModel::Record &rec)
  {
    for (int k = 0; k < n; k++)
    {
      space_op.GetExcitationVector(ex_idx[k], omega, rhs[k]);
    }
    ds->Condense(omega, std::move(rec.S));
    ds->Solve(rhs_ptr, E, &rec.g);
  };

  // The environment's port values V_jk = c_jk + h_j^T u_Γ,k (replicated), from the model's
  // record and the interface solution on rank 0.
  auto port_values = [&](const DrivenSubstructureModel::Record &rec,
                         const std::vector<std::complex<double>> &uG)
  {
    std::vector<std::complex<double>> V(static_cast<std::size_t>(np) * n, 0.0);
    for (int k = 0; root && k < n; k++)
    {
      for (int j = 0; j < np; j++)
      {
        std::complex<double> v =
            (ex_col[k] >= 0) ? rec.c[static_cast<std::size_t>(ex_col[k]) * np + j] : 0.0;
        for (int a = 0; a < nG; a++)
        {
          v += rec.h[static_cast<std::size_t>(j) * nG + a] *
               uG[static_cast<std::size_t>(k) * nG + a];
        }
        V[static_cast<std::size_t>(k) * np + j] = v;
      }
    }
    Mpi::Broadcast(static_cast<int>(V.size()), V.data(), 0, comm);
    return V;
  };

  // From a rational model with an adaptive sweep, the region is reduced: Galerkin
  // projection onto its fields at frequencies chosen greedily, as in the adaptive driven
  // solver (the next one where the minimal rational interpolant of the fields has its least
  // denominator).
  const bool region_adaptive =
      from_model && iodata.solver.driven.adaptive_tol > 0.0 &&
      (adaptive || !model.weights.empty()) &&
      omega_sample.size() > iodata.solver.driven.prom_indices.size();
  if (online && iodata.solver.driven.adaptive_tol > 0.0 && model.weights.empty())
  {
    Mpi::Warning("An adaptive online sweep needs a model saved by an adaptive sweep: every "
                 "frequency is solved!\n");
  }

  // The environment model as Σ_i c_i(ω) x_i over records x_i in the saved interface basis
  // (on rank 0): the basis of the adaptive fit (scaled part by part), or the records of a
  // saved rational model.
  const auto &env_x = adaptive ? fit.Basis() : records;
  auto env_coef = [&](double omega)
  {
    return adaptive
               ? fit.BasisCoefficients(omega)
               : BarycentricInterpolant::Coefficients(model.omega, model.weights, omega);
  };
  auto part_factor = [&](int p) { return adaptive ? 1.0 / part_scale[p] : 1.0; };
  const auto part_size = model.PartSizes();
  const std::size_t g_off = part_size[0], h_off = g_off + part_size[1],
                    c_off = h_off + part_size[2];

  // Interface values in the saved basis, u_saved = M^-1 u (online; M^-T has rows rows[g]).
  auto to_saved = [&](const auto *x)
  {
    std::vector<std::complex<double>> y(nG, 0.0);
    for (int g = 0; g < nG; g++)
    {
      if (!online)
      {
        y[g] = x[g];
        continue;
      }
      for (const auto &[a, c] : gamma_map.rows[g])
      {
        y[a] += c * x[g];
      }
    }
    return y;
  };

  // The environment's reduced terms on rank 0, with W = M^-1 V_Γ: W^T S_i W and W^T g_i per
  // record (and S_i W, to extend them as the basis grows).
  std::vector<std::complex<double>> W;
  std::vector<std::vector<std::complex<double>>> SW(env_x.size()), RS(env_x.size()),
      RG(env_x.size());
  auto update_env_terms = [&]()
  {
    const int r = ds->ReducedDimension(), r0 = static_cast<int>(W.size()) / std::max(nG, 1);
    const auto &VG = ds->ReducedInterfaceBasis();
    for (int j = r0; j < r; j++)
    {
      const auto w = to_saved(VG.data() + static_cast<std::size_t>(j) * nG);
      W.insert(W.end(), w.begin(), w.end());
    }
    for (std::size_t i = 0; i < env_x.size(); i++)
    {
      const auto *S = env_x[i].data();
      SW[i].resize(static_cast<std::size_t>(nG) * r, 0.0);
      for (int j = r0; j < r; j++)
      {
        // S_i (symmetric, by its lower triangle) times column j of W.
        const auto *wj = W.data() + static_cast<std::size_t>(j) * nG;
        auto *y = SW[i].data() + static_cast<std::size_t>(j) * nG;
        for (int c = 0, q = 0; c < nG; c++)
        {
          for (int a = c; a < nG; a++, q++)
          {
            y[a] += S[q] * wj[c];
            if (a != c)
            {
              y[c] += S[q] * wj[a];
            }
          }
        }
      }
      std::vector<std::complex<double>> R(static_cast<std::size_t>(r) * r, 0.0),
          G(static_cast<std::size_t>(r) * ne, 0.0);
      for (int j = 0; j < r0; j++)
      {
        std::copy(RS[i].begin() + static_cast<std::ptrdiff_t>(j) * r0,
                  RS[i].begin() + static_cast<std::ptrdiff_t>(j + 1) * r0,
                  R.begin() + static_cast<std::ptrdiff_t>(j) * r);
        for (int e = 0; e < ne; e++)
        {
          G[static_cast<std::size_t>(e) * r + j] =
              RG[i][static_cast<std::size_t>(e) * r0 + j];
        }
      }
      for (int j = r0; j < r; j++)
      {
        for (int a = 0; a < r; a++)
        {
          std::complex<double> d = 0.0;
          for (int q = 0; q < nG; q++)
          {
            d += W[static_cast<std::size_t>(a) * nG + q] *
                 SW[i][static_cast<std::size_t>(j) * nG + q];
          }
          R[static_cast<std::size_t>(j) * r + a] = R[static_cast<std::size_t>(a) * r + j] =
              d;
        }
        for (int e = 0; e < ne; e++)
        {
          std::complex<double> d = 0.0;
          for (int q = 0; q < nG; q++)
          {
            d += W[static_cast<std::size_t>(j) * nG + q] *
                 env_x[i][g_off + static_cast<std::size_t>(e) * nG + q];
          }
          G[static_cast<std::size_t>(e) * r + j] = d;
        }
      }
      RS[i] = std::move(R);
      RG[i] = std::move(G);
    }
  };

  // At ω, the environment's reduced terms (rank 0), and its port values from the interface
  // solution (replicated).
  auto env_terms = [&](double omega, std::vector<std::complex<double>> &A_env,
                       std::vector<std::complex<double>> &b_env)
  {
    const int r = ds->ReducedDimension();
    A_env.assign(static_cast<std::size_t>(r) * r, 0.0);
    b_env.assign(static_cast<std::size_t>(r) * n, 0.0);
    if (!root)
    {
      return;
    }
    const auto c = env_coef(omega);
    for (std::size_t i = 0; i < env_x.size(); i++)
    {
      for (std::size_t q = 0; q < A_env.size(); q++)
      {
        A_env[q] += c[i] * part_factor(0) * RS[i][q];
      }
      for (int k = 0; k < n; k++)
      {
        for (int a = 0; ex_col[k] >= 0 && a < r; a++)
        {
          b_env[static_cast<std::size_t>(k) * r + a] +=
              c[i] * part_factor(1) * RG[i][static_cast<std::size_t>(ex_col[k]) * r + a];
        }
      }
    }
  };
  auto env_port_values = [&](double omega, const std::vector<std::complex<double>> &uG)
  {
    std::vector<std::complex<double>> V(static_cast<std::size_t>(np) * n, 0.0);
    if (root)
    {
      const auto c = env_coef(omega);
      std::vector<std::complex<double>> h(static_cast<std::size_t>(nG) * np, 0.0),
          cc(static_cast<std::size_t>(np) * ne, 0.0);
      for (std::size_t i = 0; i < env_x.size(); i++)
      {
        for (std::size_t q = 0; q < h.size(); q++)
        {
          h[q] += c[i] * part_factor(2) * env_x[i][h_off + q];
        }
        for (std::size_t q = 0; q < cc.size(); q++)
        {
          cc[q] += c[i] * part_factor(3) * env_x[i][c_off + q];
        }
      }
      for (int k = 0; k < n; k++)
      {
        const auto us = to_saved(uG.data() + static_cast<std::size_t>(k) * nG);
        for (int j = 0; j < np; j++)
        {
          std::complex<double> v =
              (ex_col[k] >= 0) ? cc[static_cast<std::size_t>(ex_col[k]) * np + j] : 0.0;
          for (int a = 0; a < nG; a++)
          {
            v += h[static_cast<std::size_t>(j) * nG + a] * us[a];
          }
          V[static_cast<std::size_t>(k) * np + j] = v;
        }
      }
    }
    Mpi::Broadcast(static_cast<int>(V.size()), V.data(), 0, comm);
    return V;
  };

  BarycentricInterpolant field(comm);  // the fields' minimal rational interpolant
  const int nt = static_cast<int>(Curl.Width());
  if (region_adaptive)
  {
    BlockTimer bt(Timer::KSP);
    const double tol = iodata.solver.driven.adaptive_tol;
    const std::size_t max_samples = iodata.solver.driven.adaptive_max_size,
                      memory_size = iodata.solver.driven.adaptive_memory;
    Mpi::Print(
        "\nAdaptive sampling of the region (tol = {:.3e}, at most {:d} frequencies):\n",
        tol, max_samples);
    // The relative error of the reduced region at a new sample, before adding it.
    auto sample = [&](double omega)
    {
      auto rec = root ? model_record(omega, 0) : DrivenSubstructureModel::Record();
      solve_region(omega, rec);
      double err = std::numeric_limits<double>::infinity();
      if (ds->ReducedDimension() > 0)
      {
        std::vector<std::complex<double>> A_env, b_env;
        std::vector<ComplexVector> Eh;
        env_terms(omega, A_env, b_env);
        ds->SolveReduced(omega, rhs_ptr, A_env, b_env, Eh);
        double d[2] = {0.0, 0.0};
        for (int k = 0; k < n; k++)
        {
          ComplexVector diff(Eh[k]);
          diff -= E[k];
          d[0] += std::pow(linalg::Norml2(comm, diff), 2);
          d[1] += std::pow(linalg::Norml2(comm, E[k]), 2);
        }
        err = std::sqrt(d[0] / d[1]);
      }
      ds->AddReducedBasis(E);
      if (root)
      {
        update_env_terms();
      }
      std::vector<std::complex<double>> x;
      x.reserve(static_cast<std::size_t>(n) * nt);
      for (int k = 0; k < n; k++)
      {
        const double *er = E[k].Real().HostRead(), *ei = E[k].Imag().HostRead();
        for (int i = 0; i < nt; i++)
        {
          x.emplace_back(er[i], ei[i]);
        }
      }
      field.AddSample(omega, x);
      return err;
    };
    std::vector<double> errors;
    for (auto i : iodata.solver.driven.prom_indices)
    {
      errors.push_back(sample(omega_sample[i]));
    }
    errors[0] = errors[1] = std::numeric_limits<double>::infinity();
    std::size_t memory = 0;
    for (auto it = errors.rbegin(); it != errors.rend() && *it < tol; ++it)
    {
      memory++;
    }
    while (field.Samples().size() < max_samples && memory < memory_size)
    {
      double omega = field.FindMaxError();
      Mpi::Broadcast(1, &omega, 0, comm);
      errors.push_back(sample(omega));
      memory = (errors.back() < tol) ? memory + 1 : 0;
      Mpi::Print(
          " Greedy iteration {:d} (n = {:d}): ω* = {:.6e} GHz, error = {:.3e}, memory "
          "= {:d}/{:d}\n",
          field.Samples().size() - iodata.solver.driven.prom_indices.size(),
          ds->ReducedDimension(), omega * unit_GHz, errors.back(), memory, memory_size);
    }
    Mpi::Print(
        "\nAdaptive sampling{} {:d} frequency samples: n = {:d}, error = {:.3e}, tol = "
        "{:.3e}, memory = {:d}/{:d}\n",
        (memory < memory_size) ? " reached maximum" : " converged with",
        field.Samples().size(), ds->ReducedDimension(), errors.back(), tol, memory,
        memory_size);
    utils::PrettyPrint(field.Samples(), unit_GHz, " Sampled frequencies (GHz):");
    utils::PrettyPrint(errors, 1.0, " Sample errors:");
  }

  ComplexVector B(Curl.Height());
  B.UseDevice(true);
  auto t0 = Timer::Now();
  for (std::size_t omega_i = 0; omega_i < omega_sample.size(); omega_i++)
  {
    const double omega = omega_sample[omega_i];
    Mpi::Print("\nIt {:d}/{:d}: ω/2π = {:.3e} GHz (total elapsed time = {:.2e} s)\n",
               omega_i + 1, omega_sample.size(), omega * unit_GHz,
               Timer::Duration(Timer::Now() - t0).count());
    std::vector<std::complex<double>> V;  // the environment's port values
    {
      BlockTimer bt(Timer::KSP);
      if (!from_model)
      {
        const auto rec = condense(omega, save);
        if (save && root)
        {
          model.AppendRecord(model_path, rec);
        }
      }
      else if (region_adaptive)
      {
        // The fields and interface solution from the reduced region.
        for (int k = 0; k < n; k++)
        {
          space_op.GetExcitationVector(ex_idx[k], omega, rhs[k]);
        }
        std::vector<std::complex<double>> A_env, b_env;
        env_terms(omega, A_env, b_env);
        ds->SolveReduced(omega, rhs_ptr, A_env, b_env, E);
        V = env_port_values(omega, ds->InterfaceSolution());
      }
      else
      {
        auto rec = root ? model_record(omega, omega_i) : DrivenSubstructureModel::Record();
        solve_region(omega, rec);
        V = port_values(rec, ds->InterfaceSolution());
      }
    }
    for (int k = 0; k < n; k++)
    {
      BlockTimer bt(Timer::POSTPRO);
      // B = -1/(iω) ∇ x E on the true dofs.
      Curl.Mult(E[k].Real(), B.Real());
      Curl.Mult(E[k].Imag(), B.Imag());
      B *= -1.0 / (1i * omega);
      if (from_model)
      {
        std::map<int, std::complex<double>> Vk, Sk;
        for (int j = 0; j < np; j++)
        {
          const int p = env_ports[j];
          ((p > 0) ? Vk[p] : Sk[-p]) = V[static_cast<std::size_t>(k) * np + j];
        }
        post_op.SetLumpedPortVoltages(std::move(Vk));
        post_op.SetWavePortOverlaps(std::move(Sk));
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
