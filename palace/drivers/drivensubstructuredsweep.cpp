// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "drivensubstructuredsweep.hpp"

#include <algorithm>
#include <complex>
#include <cstddef>
#include <limits>
#include <map>
#include <memory>
#include <numbers>
#include <string>
#include <vector>
#include <fmt/core.h>
#include <mfem.hpp>
#include "fem/errorindicator.hpp"
#include "fem/substructure.hpp"
#include "linalg/vector.hpp"
#include "models/drivensubstructure.hpp"
#include "models/lumpedportoperator.hpp"
#include "models/portexcitations.hpp"
#include "models/postoperator.hpp"
#include "models/spaceoperator.hpp"
#include "models/waveportoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"
#include "utils/timer.hpp"

namespace palace
{

using namespace std::complex_literals;

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

// The parameters of the environment's ports that enter S_E (see
// DrivenSubstructureModel::port_data).
std::vector<double> EnvironmentPortData(const IoData &iodata,
                                        const std::vector<int> &env_ports)
{
  std::vector<double> data;
  for (int p : env_ports)
  {
    if (p > 0)
    {
      const auto &lp = iodata.boundaries.lumpedport.at(p);
      data.insert(data.end(),
                  {lp.R, lp.L, lp.C, lp.Rs, lp.Ls, lp.Cs, lp.active ? 1.0 : 0.0});
    }
    else
    {
      const auto &wp = iodata.boundaries.waveport.at(-p);
      data.insert(data.end(), {static_cast<double>(wp.mode_idx), wp.d_offset,
                               wp.active ? 1.0 : 0.0, 0.0, 0.0, 0.0, 0.0});
    }
  }
  return data;
}

// The functional of a wave port's S-parameter at ω as a true-DOF vector l (S = l^T E), as
// in WavePortData::GetSParameter: s / R, the unconjugated reaction of the field with the
// port's n×H mode vector s relative to that of the mode itself, R.
ComplexVector WaveOverlapFunctional(SpaceOperator &space_op, int idx, double omega)
{
  auto &wp_op = space_op.GetWavePortOp();
  auto s = wp_op.GetModeVector(idx, omega, space_op.GetNDSpace(),
                               space_op.GetNDDbcTDofLists().back());
  ComplexVector l(s->Size());
  l.Real() = s->Real();
  l.Imag() = s->Imag();
  l *= 1.0 / wp_op.GetPort(idx).modal_reaction;
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

// The index of x in y (relative tolerance tol), -1 if none.
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

// The frequency sweep of driven substructuring. Exact per-frequency substructuring: at each
// frequency, condense the environment and factor the region against it once, then solve all
// excitations (so the frequency loop is the outer one). An adaptive sweep condenses the
// environment at frequencies chosen greedily for a rational model of its condensed data
// (BarycentricInterpolant), then solves every frequency on the region against the model.
// With a saved model, an offline sweep saves the condensed data per frequency, and an
// online sweep factors the region only: at the saved frequencies, or anywhere in the band
// of a rational model.
class DrivenSubstructuredSweep
{
public:
  DrivenSubstructuredSweep(const IoData &iodata, SpaceOperator &space_op,
                           bool adaptive_sweep);
  ErrorIndicator Run();

private:
  using Record = DrivenSubstructureModel::Record;

  const IoData &iodata;
  SpaceOperator &space_op;
  const PortExcitations &port_excitations;
  const std::vector<double> &omega_sample;
  const config::SubstructuringData &sub;
  const bool online, adaptive, from_model;
  const std::string &model_path;
  const bool save;
  MPI_Comm comm;
  const bool root;
  const Operator &Curl;
  const double unit_GHz;
  const int nt;
  DrivenSubstructureModel model;
  std::unique_ptr<DrivenSubstructure> ds;
  int nG = 0, n = 0;
  std::vector<int> ex_idx;
  std::vector<ComplexVector> rhs, E;
  std::vector<const ComplexVector *> rhs_ptr;
  // The lumped ports of the environment (for a model).
  std::map<int, std::vector<char>> port_dofs;
  std::vector<int> env_ports;
  // A restart continues the sweep in its order, all excitations at a frequency before the
  // next one: "Restart" - 1 results are written.
  std::size_t first_freq = 0, first_ex = 0;
  // From a model, the field is known in the region (and on Γ) only: energies of region
  // domains, without totals or participation ratios, surface fluxes and probes the region
  // holds, and the environment's port voltages from the model (without their power).
  config::DomainData domains;
  config::BoundaryData boundaries;
  std::unique_ptr<PostOperator<ProblemType::DRIVEN>> post_op;
  // The model: header and environment ports.
  std::vector<ComplexVector> port_l;
  SignatureMap gamma_map;           // online: the saved interface on the current one
  std::vector<int> record, ex_col;  // the records of an exact model, model columns
  std::vector<std::vector<std::complex<double>>> records;  // a rational model (root)
  int np = 0, ne = 0;
  std::size_t g_off = 0, h_off = 0, c_off = 0;  // the parts of a flat record
  // Adaptive: the samples of the rational model, chosen greedily as the adaptive driven
  // solver chooses its own, with the records scaled part by part (by their norms at the
  // first sample) so that the error measure weighs them alike.
  BarycentricInterpolant fit;
  std::vector<double> part_scale;
  // The environment's reduced terms on rank 0, with W = M^-1 V_Γ: W^T S_i W and W^T g_i per
  // record (and S_i W, to extend them as the basis grows).
  std::vector<std::complex<double>> W;
  std::vector<std::vector<std::complex<double>>> SW, RS, RG;
  BarycentricInterpolant field;  // the fields' minimal rational interpolant

  // The model of an offline sweep that saves one, or of an adaptive sweep, and the model
  // of an online sweep.
  void SetUpModel();
  void OpenModel();

  // The tolerance of an excitation's source fingerprint: the modes of the environment's
  // wave ports, and so their sources, are reproducible to their eigensolver tolerance over
  // the spectral gap only.
  double SourceTol(int ex) const;
  // The sources of a model's excitations lie away from Γ (as the lumped ports do).
  void CheckSource(int k) const;
  // The functionals l_j of the environment's ports at ω.
  std::vector<const ComplexVector *> PortFunctionals(double omega);
  // Offline: condense the environment at ω and solve all excitations, and the record of the
  // model (on rank 0), with the condensed functionals of the environment's ports,
  // h_j = Reduce(l_j), and c_jk = l_j^T u_k - h_j^T u_Γ,k.
  Record CondenseAndSolve(double omega, bool with_record);
  // The record alone, from the environment (not factoring the region): with the condensed
  // sources g_k and interior solutions x_k = A_EE^-1 b_E,k (Γ held at 0) of the excitations
  // with environment sources, c_jk = l_j^T x_k, since l_j^T u_k = l_j^T x_k + h_j^T u_Γ,k.
  Record EnvironmentRecord(double omega);
  void ScaleRecord(std::vector<std::complex<double>> &x, bool inverse);
  // Adaptive: sample the environment for the rational model, and check its passivity.
  void SampleEnvironment();
  // The model's record at a frequency of the sweep, on the current interface (on rank 0;
  // online, dual quantities transform as x_cur = M^-T x): S_E, the condensed environment
  // sources g of the excitations (|Γ| x n), and the environment ports' h and c.
  Record ModelRecord(double omega, std::size_t omega_i);
  // From the model at ω: the region condensed against S_E, and all excitations solved.
  void SolveRegion(double omega, Record &rec);
  // The environment's port values V_jk = c_jk + h_j^T u_Γ,k (replicated), from the model's
  // record and the interface solution on rank 0.
  std::vector<std::complex<double>>
  PortValues(const Record &rec, const std::vector<std::complex<double>> &uG) const;
  // The environment model as Σ_i c_i(ω) x_i over records x_i in the saved interface basis
  // (on rank 0): the basis of the adaptive fit (scaled part by part), or the records of a
  // saved rational model.
  const std::vector<std::vector<std::complex<double>>> &EnvRecords() const
  {
    return adaptive ? fit.Basis() : records;
  }
  std::vector<std::complex<double>> EnvCoefficients(double omega) const;
  double PartFactor(int p) const { return adaptive ? 1.0 / part_scale[p] : 1.0; }
  // Interface values in the saved basis, u_saved = M^-1 u (online; M^-T has rows rows[g]).
  template <typename T>
  std::vector<std::complex<double>> ToSaved(const T *x) const;
  void UpdateEnvTerms();
  // At ω, the environment's reduced terms (rank 0), and its port values from the interface
  // solution (replicated).
  void EnvTerms(double omega, std::vector<std::complex<double>> &A_env,
                std::vector<std::complex<double>> &b_env) const;
  std::vector<std::complex<double>>
  EnvPortValues(double omega, const std::vector<std::complex<double>> &uG) const;
  // From a rational model with an adaptive sweep: sample the region for its reduced model.
  void SampleRegion();
};

DrivenSubstructuredSweep::DrivenSubstructuredSweep(const IoData &iodata,
                                                   SpaceOperator &space_op,
                                                   bool adaptive_sweep)
  : iodata(iodata), space_op(space_op), port_excitations(space_op.GetPortExcitations()),
    omega_sample(iodata.solver.driven.sample_f), sub(*iodata.solver.substructuring),
    online(sub.mode == SubstructuringMode::ONLINE), adaptive(adaptive_sweep && !online),
    from_model(online || adaptive), model_path(sub.save_model),
    save(!online && !model_path.empty()), comm(space_op.GetComm()), root(Mpi::Root(comm)),
    Curl(space_op.GetCurlMatrix()),
    unit_GHz(iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(1.0) /
             (2 * std::numbers::pi)),
    nt(static_cast<int>(Curl.Width())), domains(iodata.domains),
    boundaries(iodata.boundaries), field(comm)
{
  if (online)
  {
    model.ReadHeader(model_path, comm);
  }
  ds = std::make_unique<DrivenSubstructure>(space_op, sub.region_attributes,
                                            sub.environment_attributes, online,
                                            online ? &model : nullptr);
  nG = ds->InterfaceSize();
  Mpi::Print("\nSubstructuring: |Γ| = {:d} interface unknowns\n", nG);
  for (const auto &[idx, spec] : port_excitations)
  {
    ex_idx.push_back(idx);
  }
  n = static_cast<int>(ex_idx.size());
  rhs.assign(n, ComplexVector(Curl.Width()));
  rhs_ptr.resize(n);
  for (int k = 0; k < n; k++)
  {
    rhs[k].UseDevice(true);
    rhs_ptr[k] = &rhs[k];
  }
  if (save || from_model)
  {
    env_ports = EnvironmentPorts(iodata, *ds, port_dofs);
  }
  const std::size_t restart = iodata.solver.driven.restart - 1;
  first_freq = restart / n;
  first_ex = restart % n;
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
    // On the mesh of the region, which lacks the environment's part of a surface, the
    // boundaries in the region are those of the model too.
    auto contains = [](const std::vector<int> &v, int a)
    { return std::ranges::find(v, a) != v.end(); };
    const auto &mesh = space_op.GetMesh().Get();
    auto bdr_in_region = BoundaryAttributesInRegion(mesh, sub.region_attributes);
    if (ds->RegionOnly())
    {
      std::erase_if(bdr_in_region,
                    [&](int a) { return !contains(model.region_bdr_attributes, a); });
    }
    std::vector<int> skipped_flux, skipped_probe;
    for (auto it = boundaries.postpro.flux.begin(); it != boundaries.postpro.flux.end();)
    {
      const auto &attrs = it->second.attributes;
      if (attrs.empty() ||
          !std::ranges::all_of(attrs, [&](int a) { return contains(bdr_in_region, a); }))
      {
        skipped_flux.push_back(it->first);
        it = boundaries.postpro.flux.erase(it);
      }
      else
      {
        ++it;
      }
    }
    for (auto it = domains.postpro.probe.begin(); it != domains.postpro.probe.end();)
    {
      if (!PointInRegion(mesh, it->second.center, sub.region_attributes))
      {
        skipped_probe.push_back(it->first);
        it = domains.postpro.probe.erase(it);
      }
      else
      {
        ++it;
      }
    }
    auto none = [](const std::vector<int> &idx, const char *what)
    {
      return idx.empty() ? std::string()
                         : fmt::format(", no {} ({})", what, fmt::join(idx, ", "));
    };
    Mpi::Warning("{} driven substructuring writes no total energies or participation "
                 "ratios{}{}{}{}!\n",
                 online ? "Online" : "Adaptive",
                 none(skipped, "energies of the domain postprocessing outside the region"),
                 none(skipped_flux, "surface fluxes outside the region"),
                 none(skipped_probe, "probes outside the region"),
                 env_ports.empty() ? std::string()
                                   : ", and no power (or wave port voltage) of the "
                                     "environment's ports");
  }
  post_op = std::make_unique<PostOperator<ProblemType::DRIVEN>>(
      iodata.problem, iodata.solver, domains, boundaries, iodata.units, space_op);
  if (from_model && post_op->WillWriteFields())
  {
    Mpi::Warning("{} driven substructuring writes the fields in the region and on the "
                 "interface only (zero in the environment)!\n",
                 online ? "Online" : "Adaptive");
  }
  if (save || adaptive)
  {
    SetUpModel();
  }
  else if (online)
  {
    OpenModel();
  }
  np = static_cast<int>(env_ports.size());
  ne = static_cast<int>(model.excitations.size());
  const auto part_size = model.PartSizes();
  g_off = part_size[0];
  h_off = g_off + part_size[1];
  c_off = h_off + part_size[2];
}

void DrivenSubstructuredSweep::SetUpModel()
{
  model.nG = nG;
  model.sig_w = SignatureWidth(space_op.GetNDSpace().Get());
  model.signatures =
      TrueDofSignatures(space_op.GetNDSpace().Get(), ds->InterfaceIndex(), nG);
  model.env_fp = ds->EnvironmentFingerprint(omega_sample[0]);
  for (int k = 0; k < n; k++)
  {
    space_op.GetExcitationVector(ex_idx[k], omega_sample[0], rhs[k]);
    CheckSource(k);
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
  model.port_data = EnvironmentPortData(iodata, env_ports);
  model.region_bdr_attributes =
      BoundaryAttributesInRegion(space_op.GetMesh().Get(), sub.region_attributes);
  model.length_scale = iodata.model.Lc;
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
  if (save)
  {
    Mpi::Print(" Saving the substructuring model to {} ({}{:d} frequencies, {:d} "
               "excitations with environment sources, {:d} environment ports)\n",
               model_path, adaptive ? "up to " : "", model.capacity,
               model.excitations.size(), env_ports.size());
    // The mesh of the region for online runs, with the environment's ports, and the
    // element moments of the environment and of its part in that mesh (in the header).
    std::vector<int> port_attrs;
    for (int p : env_ports)
    {
      if (p > 0)
      {
        for (const auto &elem : iodata.boundaries.lumpedport.at(p).elements)
        {
          port_attrs.insert(port_attrs.end(), elem.attributes.begin(),
                            elem.attributes.end());
        }
      }
      else
      {
        const auto &attrs = iodata.boundaries.waveport.at(-p).attributes;
        port_attrs.insert(port_attrs.end(), attrs.begin(), attrs.end());
      }
    }
    const auto region_mesh = RegionMeshPath(model_path);
    const std::vector<char> written = ds->WriteRegionMesh(
        region_mesh, port_attrs, iodata.units.GetMeshLengthRelativeScale());
    Mpi::Print(" Writing the mesh of the region for online runs to {}\n", region_mesh);
    model.env_moments = ds->EnvironmentMoments();
    model.region_moments = ds->EnvironmentMoments(written);
  }
  if (save && first_freq > 0)
  {
    // A restarted sweep continues the saved model, from its records before the restart.
    DrivenSubstructureModel saved;
    saved.ReadHeader(model_path, comm);
    double tol = 1.0e-10;
    for (int ex : model.excitations)
    {
      tol = std::max(tol, SourceTol(ex));
    }
    int nrec = root ? saved.NumRecords(model_path) : 0;
    Mpi::Broadcast(1, &nrec, 0, comm);
    MFEM_VERIFY(model.SameHeader(saved, tol) && nrec >= static_cast<int>(first_freq),
                "The substructuring model \""
                    << model_path
                    << "\" does not continue this sweep: a restarted offline sweep needs "
                       "the model it saved before the restart (in the same "
                       "configuration, on as many processes)!");
    if (root)
    {
      model.TruncateRecords(model_path, static_cast<int>(first_freq));
    }
    Mpi::Print(" Continuing the substructuring model from its first {:d} records\n",
               first_freq);
  }
  else if (save && root)
  {
    model.WriteHeader(model_path);
  }
}

void DrivenSubstructuredSweep::OpenModel()
{
  MFEM_VERIFY(model.nG == nG,
              "Saved substructuring model interface size ("
                  << model.nG << ") does not match this run (" << nG
                  << "); the interface (Gamma) must be identical between the offline and "
                     "online runs.");
  gamma_map = MatchSignatureBasis(
      TrueDofSignatures(space_op.GetNDSpace().Get(), ds->InterfaceIndex(), nG),
      model.signatures, model.sig_w);
  if (ds->RegionOnly())
  {
    Mpi::Print(" The mesh holds the region (the environment is that of the model)\n");
  }
  else
  {
    MFEM_VERIFY(DrivenSubstructureModel::SameEnvironment(
                    ds->EnvironmentFingerprint(model.omega[0]), model.env_fp),
                kEnvironmentDiffers);
  }
  MFEM_VERIFY(env_ports == model.ports,
              "The environment's lumped ports differ from those of the saved "
              "substructuring model!");
  MFEM_VERIFY(EnvironmentPortData(iodata, env_ports) == model.port_data,
              "The parameters of the environment's ports (the R, L and C of a lumped port, "
              "the mode of a wave port) differ from those of the saved substructuring "
              "model!");
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
  // The model columns of the excitations with environment sources (-1: none).
  for (int k = 0; k < n; k++)
  {
    space_op.GetExcitationVector(ex_idx[k], model.omega[0], rhs[k]);
    CheckSource(k);
    const auto fp = ds->SourceFingerprint(rhs[k]);
    const double tol = SourceTol(ex_idx[k]);
    int col = -1;
    for (std::size_t e = 0; e < model.excitations.size(); e++)
    {
      if (model.excitations[e] == ex_idx[k] &&
          DrivenSubstructureModel::SameSource(
              fp.data(), model.exc_fp.data() + DrivenSubstructureModel::kSourceFp * e, tol))
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

// The tolerance of an excitation's source fingerprint: the modes of the environment's
// wave ports, and so their sources, are reproducible to their eigensolver tolerance over
// the spectral gap only.
double DrivenSubstructuredSweep::SourceTol(int ex) const
{
  double tol = 1.0e-10;
  for (int p : port_excitations.excitations.at(ex).wave_port)
  {
    if (std::ranges::find(env_ports, -p) != env_ports.end())
    {
      tol = std::max(tol, 100.0 * iodata.boundaries.waveport.at(p).eig_tol);
    }
  }
  return tol;
}

// The sources of a model's excitations lie away from Γ (as the lumped ports do).
void DrivenSubstructuredSweep::CheckSource(int k) const
{
  MFEM_VERIFY(!ds->OnInterface(rhs[k]),
              "Excitation " << ex_idx[k]
                            << " has sources on the substructuring interface: a saved "
                               "driven model, or an adaptive sweep, needs the sources "
                               "(surface currents, current dipoles) away from it!");
}

// The functionals l_j of the environment's ports at ω.
std::vector<const ComplexVector *> DrivenSubstructuredSweep::PortFunctionals(double omega)
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
}

// Offline: condense the environment at ω and solve all excitations, and the record of the
// model (on rank 0), with the condensed functionals of the environment's ports,
// h_j = Reduce(l_j), and c_jk = l_j^T u_k - h_j^T u_Γ,k.
DrivenSubstructureModel::Record DrivenSubstructuredSweep::CondenseAndSolve(double omega,
                                                                           bool with_record)
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
  rec.h = ds->CondenseEnvironment(PortFunctionals(omega));
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
}

// The record alone, from the environment (not factoring the region): with the condensed
// sources g_k and interior solutions x_k = A_EE^-1 b_E,k (Γ held at 0) of the excitations
// with environment sources, c_jk = l_j^T x_k, since l_j^T u_k = l_j^T x_k + h_j^T u_Γ,k.
DrivenSubstructureModel::Record DrivenSubstructuredSweep::EnvironmentRecord(double omega)
{
  std::vector<const ComplexVector *> B(ne);
  for (int e = 0; e < ne; e++)
  {
    const int k =
        static_cast<int>(std::ranges::find(ex_idx, model.excitations[e]) - ex_idx.begin());
    space_op.GetExcitationVector(ex_idx[k], omega, rhs[k]);
    B[e] = &rhs[k];
  }
  ds->Condense(omega, false);
  DrivenSubstructureModel::Record rec;
  rec.h = ds->CondenseEnvironment(PortFunctionals(omega));
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
}

void DrivenSubstructuredSweep::ScaleRecord(std::vector<std::complex<double>> &x,
                                           bool inverse)
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
}

void DrivenSubstructuredSweep::SampleEnvironment()
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
    auto rec = EnvironmentRecord(omega);
    double err = std::numeric_limits<double>::infinity();
    if (root)
    {
      if (save)
      {
        model.AppendRecord(model_path, rec);
      }
      auto x = model.Flatten(rec);
      ScaleRecord(x, false);
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
  Mpi::Print("\nAdaptive sampling{} {:d} frequency samples: error = {:.3e}, tol = {:.3e}, "
             "memory = {:d}/{:d}\n",
             (memory < memory_size) ? " reached maximum" : " converged with",
             model.omega.size(), errors.back(), tol, memory, memory_size);
  // Passivity between the samples (Im S_E positive semidefinite), tested by Cholesky
  // factorizations: a violation beyond the tolerance adds a sample where it is largest.
  double passivity = 0.0;  // the least eigenvalue of the violations, relative (0: none)
  while (true)
  {
    double at = 0.0;
    passivity = 0.0;
    if (root)
    {
      auto z = model.omega;
      std::ranges::sort(z);
      for (std::size_t j = 0; j + 1 < z.size(); j++)
      {
        const double mid = 0.5 * (z[j] + z[j + 1]);
        const auto x = fit.Evaluate(mid);
        if (!DrivenSubstructureModel::Passive(x.data(), nG, tol))
        {
          const double p = DrivenSubstructureModel::Passivity(x.data(), nG);
          if (p < passivity)
          {
            passivity = p;
            at = mid;
          }
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
  if (passivity < -tol)
  {
    Mpi::Warning("The rational environment model is not passive between its samples "
                 "(least eigenvalue of Im S_E / ‖S_E‖ = {:.3e}): raise "
                 "\"AdaptiveMaxSamples\"!\n",
                 passivity);
  }
  else
  {
    Mpi::Print(" Passive between the samples: least eigenvalue of Im S_E / ‖S_E‖ above "
               "-{:.3e}\n",
               tol);
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
DrivenSubstructureModel::Record DrivenSubstructuredSweep::ModelRecord(double omega,
                                                                      std::size_t omega_i)
{
  std::vector<std::complex<double>> x;  // the record, flat (S_E by its lower triangle)
  if (adaptive)
  {
    x = fit.Evaluate(omega);
    ScaleRecord(x, true);
  }
  else if (!model.weights.empty())
  {
    const auto a = BarycentricInterpolant::Coefficients(model.omega, model.weights, omega);
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
      std::copy(gk.begin(), gk.end(), out.g.begin() + static_cast<std::ptrdiff_t>(k) * nG);
    }
  }
  for (int j = 0; j < np; j++)
  {
    const auto hj = dual(h + static_cast<std::size_t>(j) * nG);
    out.h.insert(out.h.end(), hj.begin(), hj.end());
  }
  out.c.assign(c, c + sz[3]);
  return out;
}

// From the model at ω: the region condensed against S_E, and all excitations solved.
void DrivenSubstructuredSweep::SolveRegion(double omega,
                                           DrivenSubstructureModel::Record &rec)
{
  for (int k = 0; k < n; k++)
  {
    space_op.GetExcitationVector(ex_idx[k], omega, rhs[k]);
  }
  ds->Condense(omega, std::move(rec.S));
  ds->Solve(rhs_ptr, E, &rec.g);
}

// The environment's port values V_jk = c_jk + h_j^T u_Γ,k (replicated), from the model's
// record and the interface solution on rank 0.
std::vector<std::complex<double>>
DrivenSubstructuredSweep::PortValues(const DrivenSubstructureModel::Record &rec,
                                     const std::vector<std::complex<double>> &uG) const
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
}

std::vector<std::complex<double>>
DrivenSubstructuredSweep::EnvCoefficients(double omega) const
{
  return adaptive ? fit.BasisCoefficients(omega)
                  : BarycentricInterpolant::Coefficients(model.omega, model.weights, omega);
}

template <typename T>
std::vector<std::complex<double>> DrivenSubstructuredSweep::ToSaved(const T *x) const
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
}

void DrivenSubstructuredSweep::UpdateEnvTerms()
{
  const int r = ds->ReducedDimension(), r0 = static_cast<int>(W.size()) / std::max(nG, 1);
  const auto &VG = ds->ReducedInterfaceBasis();
  for (int j = r0; j < r; j++)
  {
    const auto w = ToSaved(VG.data() + static_cast<std::size_t>(j) * nG);
    W.insert(W.end(), w.begin(), w.end());
  }
  for (std::size_t i = 0; i < EnvRecords().size(); i++)
  {
    const auto *S = EnvRecords()[i].data();
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
        R[static_cast<std::size_t>(j) * r + a] = R[static_cast<std::size_t>(a) * r + j] = d;
      }
      for (int e = 0; e < ne; e++)
      {
        std::complex<double> d = 0.0;
        for (int q = 0; q < nG; q++)
        {
          d += W[static_cast<std::size_t>(j) * nG + q] *
               EnvRecords()[i][g_off + static_cast<std::size_t>(e) * nG + q];
        }
        G[static_cast<std::size_t>(e) * r + j] = d;
      }
    }
    RS[i] = std::move(R);
    RG[i] = std::move(G);
  }
}

// At ω, the environment's reduced terms (rank 0), and its port values from the interface
// solution (replicated).
void DrivenSubstructuredSweep::EnvTerms(double omega,
                                        std::vector<std::complex<double>> &A_env,
                                        std::vector<std::complex<double>> &b_env) const
{
  const int r = ds->ReducedDimension();
  A_env.assign(static_cast<std::size_t>(r) * r, 0.0);
  b_env.assign(static_cast<std::size_t>(r) * n, 0.0);
  if (!root)
  {
    return;
  }
  const auto c = EnvCoefficients(omega);
  for (std::size_t i = 0; i < EnvRecords().size(); i++)
  {
    for (std::size_t q = 0; q < A_env.size(); q++)
    {
      A_env[q] += c[i] * PartFactor(0) * RS[i][q];
    }
    for (int k = 0; k < n; k++)
    {
      for (int a = 0; ex_col[k] >= 0 && a < r; a++)
      {
        b_env[static_cast<std::size_t>(k) * r + a] +=
            c[i] * PartFactor(1) * RG[i][static_cast<std::size_t>(ex_col[k]) * r + a];
      }
    }
  }
}

std::vector<std::complex<double>>
DrivenSubstructuredSweep::EnvPortValues(double omega,
                                        const std::vector<std::complex<double>> &uG) const
{
  std::vector<std::complex<double>> V(static_cast<std::size_t>(np) * n, 0.0);
  if (root)
  {
    const auto c = EnvCoefficients(omega);
    std::vector<std::complex<double>> h(static_cast<std::size_t>(nG) * np, 0.0),
        cc(static_cast<std::size_t>(np) * ne, 0.0);
    for (std::size_t i = 0; i < EnvRecords().size(); i++)
    {
      for (std::size_t q = 0; q < h.size(); q++)
      {
        h[q] += c[i] * PartFactor(2) * EnvRecords()[i][h_off + q];
      }
      for (std::size_t q = 0; q < cc.size(); q++)
      {
        cc[q] += c[i] * PartFactor(3) * EnvRecords()[i][c_off + q];
      }
    }
    for (int k = 0; k < n; k++)
    {
      const auto us = ToSaved(uG.data() + static_cast<std::size_t>(k) * nG);
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
}

void DrivenSubstructuredSweep::SampleRegion()
{
  BlockTimer bt(Timer::KSP);
  const double tol = iodata.solver.driven.adaptive_tol;
  const std::size_t max_samples = iodata.solver.driven.adaptive_max_size,
                    memory_size = iodata.solver.driven.adaptive_memory;
  Mpi::Print(
      "\nAdaptive sampling of the region (tol = {:.3e}, at most {:d} frequencies):\n", tol,
      max_samples);
  // The relative error of the reduced region at a new sample, before adding it.
  auto sample = [&](double omega)
  {
    auto rec = root ? ModelRecord(omega, 0) : DrivenSubstructureModel::Record();
    SolveRegion(omega, rec);
    double err = std::numeric_limits<double>::infinity();
    if (ds->ReducedDimension() > 0)
    {
      std::vector<std::complex<double>> A_env, b_env;
      std::vector<ComplexVector> Eh;
      EnvTerms(omega, A_env, b_env);
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
      UpdateEnvTerms();
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
    Mpi::Print(" Greedy iteration {:d} (n = {:d}): ω* = {:.6e} GHz, error = {:.3e}, memory "
               "= {:d}/{:d}\n",
               field.Samples().size() - iodata.solver.driven.prom_indices.size(),
               ds->ReducedDimension(), omega * unit_GHz, errors.back(), memory,
               memory_size);
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

ErrorIndicator DrivenSubstructuredSweep::Run()
{
  if (adaptive)
  {
    SampleEnvironment();
  }
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
  SW.assign(EnvRecords().size(), {});
  RS.assign(EnvRecords().size(), {});
  RG.assign(EnvRecords().size(), {});
  if (region_adaptive)
  {
    SampleRegion();
  }

  ComplexVector B(Curl.Height());
  B.UseDevice(true);
  auto t0 = Timer::Now();
  for (std::size_t omega_i = first_freq; omega_i < omega_sample.size(); omega_i++)
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
        const auto rec = CondenseAndSolve(omega, save);
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
        EnvTerms(omega, A_env, b_env);
        ds->SolveReduced(omega, rhs_ptr, A_env, b_env, E);
        V = EnvPortValues(omega, ds->InterfaceSolution());
      }
      else
      {
        auto rec = root ? ModelRecord(omega, omega_i) : DrivenSubstructureModel::Record();
        SolveRegion(omega, rec);
        V = PortValues(rec, ds->InterfaceSolution());
      }
    }
    for (int k = (omega_i == first_freq) ? static_cast<int>(first_ex) : 0; k < n; k++)
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
        post_op->SetLumpedPortVoltages(std::move(Vk));
        post_op->SetWavePortOverlaps(std::move(Sk));
      }
      post_op->MeasureAndPrintAll(ex_idx[k], int(omega_i), E[k], B, omega);
    }
  }
  // Substructuring computes no error estimate.
  ErrorIndicator indicator;
  post_op->MeasureFinalize(indicator);
  return indicator;
}

}  // namespace

ErrorIndicator SweepSubstructured(const IoData &iodata, SpaceOperator &space_op,
                                  bool adaptive)
{
  return DrivenSubstructuredSweep(iodata, space_op, adaptive).Run();
}

}  // namespace palace
