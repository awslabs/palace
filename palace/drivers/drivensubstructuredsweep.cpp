// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "drivensubstructuredsweep.hpp"

#include <algorithm>
#include <array>
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
#include "fem/substructure.hpp"
#include "fem/substructuremodel.hpp"
#include "linalg/blockcyclic.hpp"
#include "linalg/vector.hpp"
#include "models/drivensubstructure.hpp"
#include "models/environmentoperator.hpp"
#include "models/lumpedportoperator.hpp"
#include "models/portexcitations.hpp"
#include "models/spaceoperator.hpp"
#include "models/waveportoperator.hpp"
#include "utils/communication.hpp"
#include "utils/constants.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace
{

// The functional of a wave port's S-parameter at ω as a true-DOF vector l (S = l^T E), as
// in WavePortData::GetSParameter: s / (sqrt(2) R), the unconjugated reaction of the field
// with the port's n×H mode vector s relative to that of the incident wave, sqrt(2) times
// the mode (reaction R of the mode itself).
ComplexVector WaveOverlapFunctional(SpaceOperator &space_op, int idx, double omega)
{
  auto &wp_op = space_op.GetWavePortOp();
  auto s = wp_op.GetModeVector(idx, omega, space_op.GetNDSpace(),
                               space_op.GetNDDbcTDofLists().back());
  ComplexVector l(s->Size());
  l.Real() = s->Real();
  l.Imag() = s->Imag();
  l *= 1.0 /
       (electromagnetics::UnitPowerAmplitude(true) * wp_op.GetPort(idx).modal_reaction);
  return l;
}

// The voltage functional of a lumped port as a true-DOF vector l (V = l^T E, summed over
// the processes), from its voltage linear form.
ComplexVector VoltageFunctional(SpaceOperator &space_op, int idx)
{
  const auto &port = space_op.GetLumpedPortOp().GetPort(idx);
  const Vector v = port.GetVoltageFunctional(space_op.GetNDSpace().Get());
  ComplexVector l(v.Size());
  l.Real() = v;
  l.Imag() = 0.0;
  return l;
}

// The condensation of the environment of an offline driven substructuring, as the
// environment of the region (see EnvironmentOperator): exactly at each frequency of the
// sweep, or at frequencies chosen greedily for a rational model of its condensed data
// (BarycentricInterpolant). The region is then solved against it as online (with a saved
// model, which this saves: at the saved frequencies, or anywhere in the band of a rational
// model).
class DrivenSubstructuredSweep
{
public:
  DrivenSubstructuredSweep(const IoData &iodata, SpaceOperator &space_op,
                           bool adaptive_sweep);

  // The exact model, whose records are condensed (and saved) when the region's sweep
  // reaches their frequencies, by its record function, which holds the sweep (self).
  static std::unique_ptr<EnvironmentOperator>
  ExactModel(const std::shared_ptr<DrivenSubstructuredSweep> &self);

  // The adaptive sampling of the environment's rational model, which the adaptive sweep
  // then reduces the region against.
  std::unique_ptr<EnvironmentOperator> SampleModel();

private:
  using Record = DrivenSubstructureModel::Record;

  const IoData &iodata;
  SpaceOperator &space_op;
  const PortExcitations &port_excitations;
  const std::vector<double> &omega_sample;
  const config::SubstructuringData &sub;
  const bool adaptive;
  const std::string model_path;
  const bool save;
  MPI_Comm comm;
  const bool root;
  const double unit_GHz;
  DrivenSubstructureModel model;
  std::unique_ptr<DrivenSubstructure> ds;
  int nG = 0, n = 0;
  std::vector<int> ex_idx;
  std::vector<ComplexVector> rhs;
  // The ports of the environment, and their attributes (for the model).
  std::vector<int> env_ports, env_port_attributes;
  // A restart continues the sweep in its order, all excitations at a frequency before the
  // next one: the records of the frequencies before that of solve "Restart" are kept.
  std::size_t first_freq = 0;
  // The model: the functionals of the environment's ports, and the model column of each
  // excitation (-1: no sources in the environment).
  std::vector<ComplexVector> port_l;
  std::vector<int> ex_col;
  int np = 0, ne = 0;
  // Adaptive: the samples of the rational model, chosen greedily as the adaptive driven
  // solver chooses its own, with the records scaled part by part (by their norms at the
  // first sample) so that the error measure weighs them alike. A sample is distributed as
  // its record: each process holds its share of the lower triangle of S_E (see
  // TriangleShare), and root also the parts after it.
  BarycentricInterpolant fit;
  std::vector<double> part_scale;

  // The model of the sweep.
  void SetUpModel();

  // The tolerance of an excitation's source fingerprint: the modes of the environment's
  // wave ports, and so their sources, are reproducible to their eigensolver tolerance over
  // the spectral gap only.
  double SourceTol(int ex) const;
  // The sources of a model's excitations lie away from Γ (as the lumped ports do).
  void CheckSource(int k) const;
  // The functionals l_j of the environment's ports at ω.
  std::vector<const ComplexVector *> PortFunctionals(double omega);
  // The record of the environment at ω: the share of this process of S_E (see
  // TriangleShare), and on rank 0 the condensed functionals of the environment's ports,
  // h_j = Reduce(l_j), and with the condensed sources g_k and interior solutions
  // x_k = A_EE^-1 b_E,k (Γ held at 0) of the excitations with environment sources,
  // c_jk = l_j^T x_k, since l_j^T u_k = l_j^T x_k + h_j^T u_Γ,k.
  Record CondensedEnvironment(double omega);
  void ScaleRecord(std::vector<std::complex<double>> &x, bool inverse);
  // Adaptive: sample the environment for the rational model, and check its passivity.
  void SampleEnvironment();
};

DrivenSubstructuredSweep::DrivenSubstructuredSweep(const IoData &iodata,
                                                   SpaceOperator &space_op,
                                                   bool adaptive_sweep)
  : iodata(iodata), space_op(space_op), port_excitations(space_op.GetPortExcitations()),
    omega_sample(iodata.solver.driven.sample_f), sub(*iodata.solver.substructuring),
    adaptive(adaptive_sweep), model_path(substructuring::SaveModelPath(iodata)),
    save(!model_path.empty()), comm(space_op.GetComm()), root(Mpi::Root(comm)),
    unit_GHz(iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(1.0) /
             (2 * std::numbers::pi)),
    fit(comm)
{
  ds = std::make_unique<DrivenSubstructure>(space_op, sub.region_attributes,
                                            sub.environment_attributes, sub.direct_solver,
                                            iodata.solver.linear.sym_factorization);
  nG = ds->GetInterface().size;
  MFEM_VERIFY(nG > 0, "The substructuring interface holds no unknowns, as of a region "
                      "enclosed by PEC: the region does not couple to the environment, so "
                      "solve it without \"Substructuring\"!");
  Mpi::Print("\nSubstructuring: |Γ| = {:d} interface unknowns\n", nG);
  for (const auto &[idx, spec] : port_excitations)
  {
    ex_idx.push_back(idx);
  }
  n = static_cast<int>(ex_idx.size());
  rhs.assign(n, ComplexVector(space_op.GetNDSpace().GetTrueVSize()));
  for (int k = 0; k < n; k++)
  {
    rhs[k].UseDevice(true);
  }
  env_ports = EnvironmentPorts(iodata, *ds, env_port_attributes);
  first_freq = (iodata.solver.driven.restart - 1) / n;
  SetUpModel();
  np = static_cast<int>(env_ports.size());
  ne = static_cast<int>(model.excitations.size());
}

void DrivenSubstructuredSweep::SetUpModel()
{
  model.nG = nG;
  model.sig_w = SignatureWidth(space_op.GetNDSpace().Get());
  model.signatures =
      TrueDofSignatures(space_op.GetNDSpace().Get(), ds->GetInterface().index, nG);
  model.fp_centers = ds->FingerprintCenters();
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
  model.port_attributes = env_port_attributes;
  model.length_scale = iodata.model.Lc;
  if (adaptive)
  {
    const std::size_t n_init = iodata.solver.driven.prom_indices.size();
    MFEM_VERIFY(iodata.solver.driven.adaptive_max_size >= n_init,
                "Adaptive frequency sweep must sample at least " << n_init
                                                                 << " frequency points!");
  }
  else
  {
    model.omega = omega_sample;
  }
  for (int p : env_ports)
  {
    port_l.push_back((p > 0) ? VoltageFunctional(space_op, p) : ComplexVector());
  }
  if (save)
  {
    Mpi::Print(" Saving the substructuring model to {} ({}{:d} frequencies, {:d} "
               "excitations with environment sources, {:d} environment ports)\n",
               model_path, adaptive ? "up to " : "",
               adaptive ? iodata.solver.driven.adaptive_max_size : omega_sample.size(),
               model.excitations.size(), env_ports.size());
    // The mesh of the region for online runs, with the environment's ports, and the record
    // of the environment and of its part in that mesh.
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
    model.environment = EnvironmentRecord::Make(space_op.GetNDSpace().GetParMesh(),
                                                sub.environment_attributes,
                                                *iodata.resolved_config, written);
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

// The tolerance of an excitation's source fingerprint: the modes of the environment's
// wave ports, and so their sources, are reproducible to their eigensolver tolerance over
// the spectral gap only.
double DrivenSubstructuredSweep::SourceTol(int ex) const
{
  return SourceTolerance(iodata, port_excitations, env_ports, ex);
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

// The record alone, from the environment (not factoring the region): with the condensed
// sources g_k and interior solutions x_k = A_EE^-1 b_E,k (Γ held at 0) of the excitations
// with environment sources, c_jk = l_j^T x_k, since l_j^T u_k = l_j^T x_k + h_j^T u_Γ,k.
DrivenSubstructureModel::Record DrivenSubstructuredSweep::CondensedEnvironment(double omega)
{
  std::vector<const ComplexVector *> B(ne);
  for (int e = 0; e < ne; e++)
  {
    const int k =
        static_cast<int>(std::ranges::find(ex_idx, model.excitations[e]) - ex_idx.begin());
    space_op.GetExcitationVector(ex_idx[k], omega, rhs[k]);
    B[e] = &rhs[k];
  }
  ds->Condense(omega);
  DrivenSubstructureModel::Record rec;
  rec.h = ds->CondenseEnvironment(PortFunctionals(omega));
  std::vector<ComplexVector> x;
  auto g = ds->CondenseEnvironment(B, &x);
  std::vector<std::complex<double>> c(static_cast<std::size_t>(np) * ne);
  for (int e = 0; e < ne; e++)
  {
    for (int j = 0; j < np; j++)
    {
      c[static_cast<std::size_t>(e) * np + j] = port_l[j].TransposeDot(x[e]);
    }
  }
  Mpi::GlobalSum(static_cast<int>(c.size()), c.data(), comm);
  rec.S = ds->TakeSchur();
  if (root)
  {
    rec.g = std::move(g);
    rec.c = std::move(c);
  }
  return rec;
}

void DrivenSubstructuredSweep::ScaleRecord(std::vector<std::complex<double>> &x,
                                           bool inverse)
{
  // The parts of a sample on this process: its share of S_E, then on root g, h and c.
  const auto sz = model.PartSizes();
  const std::size_t rest = root ? sz[1] + sz[2] + sz[3] : 0;
  const std::array<std::size_t, 5> off = {0, x.size() - rest,
                                          x.size() - rest + (root ? sz[1] : 0),
                                          x.size() - (root ? sz[3] : 0), x.size()};
  if (part_scale.empty())
  {
    std::array<double, 4> nrm = {};
    for (std::size_t p = 0; p < nrm.size(); p++)
    {
      for (std::size_t q = off[p]; q < off[p + 1]; q++)
      {
        nrm[p] += std::norm(x[q]);
      }
    }
    Mpi::GlobalSum(static_cast<int>(nrm.size()), nrm.data(), comm);
    for (double v : nrm)
    {
      part_scale.push_back((v > 0.0) ? 1.0 / std::sqrt(v) : 1.0);
    }
  }
  for (std::size_t p = 0; p < part_scale.size(); p++)
  {
    const double c = inverse ? 1.0 / part_scale[p] : part_scale[p];
    for (std::size_t q = off[p]; q < off[p + 1]; q++)
    {
      x[q] *= c;
    }
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
  model.tol = tol;
  // The relative error of the rational model at a new sample, before adding it.
  auto sample = [&](double omega)
  {
    auto rec = CondensedEnvironment(omega);
    if (save)
    {
      model.AppendRecord(model_path, rec, comm);
    }
    auto x = DrivenSubstructureModel::Flatten(std::move(rec));
    ScaleRecord(x, false);
    double err = std::numeric_limits<double>::infinity();
    if (fit.Samples().size() >= 2)
    {
      const auto y = fit.Evaluate(omega);
      std::array<double, 2> dm = {0.0, 0.0};
      for (std::size_t q = 0; q < x.size(); q++)
      {
        dm[0] += std::norm(y[q] - x[q]);
        dm[1] += std::norm(x[q]);
      }
      Mpi::GlobalSum(2, dm.data(), comm);
      err = std::sqrt(dm[0] / dm[1]);
    }
    fit.AddSample(omega, x);
    model.weights = fit.Weights();
    model.omega.push_back(omega);
    // A saved model is one stopped here, until the sampling ends.
    model.converged = false;
    model.error = err;
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
  const bool converged = (memory >= memory_size);
  const double error = errors.back();
  // Passivity between the samples, tested by Cholesky factorizations: the least eigenvalue
  // of Im S_E at a midpoint is to be no lower than at the samples beside it, which is zero
  // for a passive environment and can be negative for one that is not, as with a
  // second-order absorbing boundary. A violation beyond the tolerance adds a sample where
  // it is largest. The tests run on a grid of the processes, from the shares of S_E.
  const BlockCyclicGrid grid(comm, nG);
  double passivity = 0.0;               // the least violation, relative (0: none)
  std::map<double, double> at_samples;  // the least eigenvalue at a sample, at most 0
  auto sample_passivity = [&](double omega)
  {
    auto [it, added] = at_samples.try_emplace(omega, 0.0);
    if (added)
    {
      const auto x = fit.Evaluate(omega);
      it->second = std::min(0.0, DrivenSubstructureModel::Passivity(x.data(), nG, grid));
    }
    return it->second;
  };
  while (true)
  {
    double at = 0.0;
    passivity = 0.0;
    auto z = model.omega;
    std::ranges::sort(z);
    for (std::size_t j = 0; j + 1 < z.size(); j++)
    {
      const double mid = 0.5 * (z[j] + z[j + 1]);
      const auto x = fit.Evaluate(mid);
      if (!DrivenSubstructureModel::Passive(x.data(), nG, tol, grid))
      {
        const double p = DrivenSubstructureModel::Passivity(x.data(), nG, grid) -
                         std::min(sample_passivity(z[j]), sample_passivity(z[j + 1]));
        if (p < passivity)
        {
          passivity = p;
          at = mid;
        }
      }
    }
    if (passivity >= -tol || model.omega.size() >= max_samples)
    {
      break;
    }
    errors.push_back(sample(at));
    Mpi::Print(" Passivity sample: ω = {:.6e} GHz (least eigenvalue of Im S_E / ‖S_E‖ "
               "{:.3e} below its samples'), error = {:.3e}\n",
               at * unit_GHz, passivity, errors.back());
  }
  if (passivity < -tol)
  {
    Mpi::Warning("The rational environment model is not passive between its samples "
                 "(least eigenvalue of Im S_E / ‖S_E‖ {:.3e} below the samples'): raise "
                 "\"AdaptiveMaxSamples\"!\n",
                 passivity);
  }
  else
  {
    Mpi::Print(" Passive between the samples: least eigenvalue of Im S_E / ‖S_E‖ above "
               "the samples' less {:.3e}\n",
               tol);
  }
  utils::PrettyPrint(model.omega, unit_GHz, " Sampled frequencies (GHz):");
  utils::PrettyPrint(errors, 1.0, " Sample errors:");
  model.converged = converged;
  model.error = error;
  model.passivity = passivity;
  if (save && root)
  {
    model.WriteHeader(model_path, false);
  }
}
std::unique_ptr<EnvironmentOperator> DrivenSubstructuredSweep::SampleModel()
{
  SampleEnvironment();
  // The region is factored without the environment's factorization, against the samples
  // as recorded, one at a time, each process with its part.
  ds.reset();
  return std::make_unique<EnvironmentOperator>(iodata, space_op, model,
                                               [this](std::size_t j)
                                               {
                                                 auto x = fit.Sample(j);
                                                 ScaleRecord(x, true);
                                                 return x;
                                               });
}

std::unique_ptr<EnvironmentOperator>
DrivenSubstructuredSweep::ExactModel(const std::shared_ptr<DrivenSubstructuredSweep> &self)
{
  return std::make_unique<EnvironmentOperator>(
      self->iodata, self->space_op, self->model,
      [self](std::size_t j)
      {
        BlockTimer bt(Timer::KSP);
        auto rec = self->CondensedEnvironment(self->omega_sample[j]);
        if (self->save)
        {
          self->model.AppendRecord(self->model_path, rec, self->comm);
        }
        return DrivenSubstructureModel::Flatten(std::move(rec));
      });
}

}  // namespace

std::unique_ptr<EnvironmentOperator> ExactEnvironmentModel(const IoData &iodata,
                                                           SpaceOperator &space_op)
{
  return DrivenSubstructuredSweep::ExactModel(
      std::make_shared<DrivenSubstructuredSweep>(iodata, space_op, false));
}

std::unique_ptr<EnvironmentOperator> SampleEnvironmentModel(const IoData &iodata,
                                                            SpaceOperator &space_op)
{
  return DrivenSubstructuredSweep(iodata, space_op, true).SampleModel();
}

}  // namespace palace
