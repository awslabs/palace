// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "condensedenvironment.hpp"

#include <algorithm>
#include <cstdlib>
#include <numbers>
#include <string>
#include <utility>
#include <fmt/core.h>
#include <fmt/ranges.h>
#include "fem/substructure.hpp"
#include "models/drivensubstructure.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

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

int FindValue(const std::vector<double> &y, double x, double tol)
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

void RestrictPostprocessingToRegion(const mfem::ParMesh &mesh,
                                    const config::SubstructuringData &sub,
                                    const DrivenSubstructure &ds,
                                    const DrivenSubstructureModel &model, bool env_ports,
                                    const char *mode, config::DomainData &domains,
                                    config::BoundaryData &boundaries)
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
  auto bdr_in_region = BoundaryAttributesInRegion(mesh, sub.region_attributes);
  if (ds.RegionOnly())
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
               mode, none(skipped, "energies of the domain postprocessing outside the region"),
               none(skipped_flux, "surface fluxes outside the region"),
               none(skipped_probe, "probes outside the region"),
               !env_ports ? std::string()
                          : ", and no power (or wave port voltage) of the "
                            "environment's ports");
}

// S_E as an operator on the true DOFs, by the dense rows of the local interface DOFs
// (shared with the record they come from): y_Γ = S_E x_Γ, with x_Γ gathered on every rank
// in interface order, and y = x on the pinned true DOFs which are not essential.
class CondensedEnvironment::SchurOperator : public ComplexOperator
{
public:
  SchurOperator(const CondensedEnvironment &env,
                std::shared_ptr<const std::vector<std::complex<double>>> rows)
    : ComplexOperator(env.space_op.GetNDSpace().GetTrueVSize()), env(env),
      rows(std::move(rows))
  {
  }

  void Mult(const ComplexVector &x, ComplexVector &y) const override
  {
    y = 0.0;
    AddMult(x, y, 1.0);
  }

  void AddMult(const ComplexVector &x, ComplexVector &y,
               const std::complex<double> a = 1.0) const override
  {
    const int m = static_cast<int>(env.gamma_local.size()), n = env.nG;
    const double *xr = x.Real().HostRead(), *xi = x.Imag().HostRead();
    std::vector<std::complex<double>> xl(m), xg(n);
    for (int l = 0; l < m; l++)
    {
      xl[l] = {xr[env.gamma_local[l]], xi[env.gamma_local[l]]};
    }
    MPI_Allgatherv(xl.data(), m, MPI_C_DOUBLE_COMPLEX, xg.data(), env.gamma_cnt.data(),
                   env.gamma_disp.data(), MPI_C_DOUBLE_COMPLEX, env.comm);
    double *yr = y.Real().HostReadWrite(), *yi = y.Imag().HostReadWrite();
    for (int l = 0; l < m; l++)
    {
      const auto *row = rows->data() + static_cast<std::size_t>(l) * n;
      std::complex<double> v = 0.0;
      for (int b = 0; b < n; b++)
      {
        v += row[b] * xg[b];
      }
      v *= a;
      yr[env.gamma_local[l]] += v.real();
      yi[env.gamma_local[l]] += v.imag();
    }
    for (int i : env.identity_tdofs)
    {
      const auto v = a * std::complex<double>(xr[i], xi[i]);
      yr[i] += v.real();
      yi[i] += v.imag();
    }
  }

private:
  const CondensedEnvironment &env;
  std::shared_ptr<const std::vector<std::complex<double>>> rows;
};

bool CondensedEnvironment::Enabled(const config::SolverData &solver)
{
  const char *env = std::getenv("PALACE_SUBSTRUCTURING_DRIVENSOLVER");
  return solver.substructuring &&
         solver.substructuring->mode == SubstructuringMode::ONLINE &&
         solver.driven.adaptive_tol <= 0.0 && env && std::string(env) == "1";
}

CondensedEnvironment::CondensedEnvironment(const IoData &iodata, SpaceOperator &space_op)
  : space_op(space_op), sub(*iodata.solver.substructuring), comm(space_op.GetComm()),
    model(std::make_unique<DrivenSubstructureModel>()), domains(iodata.domains),
    boundaries(iodata.boundaries)
{
  MFEM_VERIFY(!space_op.GetNDSpace().IsBroken(),
              "Online driven substructuring through the driven solver (prototype) needs an "
              "unbroken space!");
  MFEM_VERIFY(space_op.GetWavePortOp().Size() == 0,
              "Online driven substructuring through the driven solver (prototype) supports "
              "lumped ports only!");
  model->ReadHeader(sub.save_model, comm);
  MFEM_VERIFY(model->weights.empty(),
              "Online driven substructuring through the driven solver (prototype) needs an "
              "exact model (saved by a uniform sweep)!");
  ds = std::make_unique<DrivenSubstructure>(space_op, sub.region_attributes,
                                            sub.environment_attributes, true, model.get(),
                                            sub.direct_solver);
  if (ds->RegionOnly())
  {
    const auto &refinement = iodata.model.refinement;
    MFEM_VERIFY(refinement.uniform_ref_levels == 0 &&
                    refinement.ser_uniform_ref_levels == 0 &&
                    refinement.GetBoxes().empty() && refinement.GetSpheres().empty(),
                kRegionMeshRefined);
  }
  nG = ds->InterfaceSize();
  Mpi::Print("\nSubstructuring: |Γ| = {:d} interface unknowns (online, by the driven "
             "solver)\n",
             nG);
  std::map<int, std::vector<char>> port_dofs;
  env_ports = EnvironmentPorts(iodata, *ds, port_dofs);

  // The checks of the online sweep (DrivenSubstructuredSweep::OpenModel).
  const auto &fes = space_op.GetNDSpace().Get();
  MFEM_VERIFY(model->nG == nG,
              "Saved substructuring model interface size ("
                  << model->nG << ") does not match this run (" << nG
                  << "); the interface (Gamma) must be identical between the offline and "
                     "online runs.");
  gamma_map = MatchSignatureBasis(TrueDofSignatures(fes, ds->InterfaceIndex(), nG),
                                  model->signatures, model->sig_w);
  if (ds->RegionOnly())
  {
    Mpi::Print(" The mesh holds the region (the environment is that of the model)\n");
  }
  else
  {
    MFEM_VERIFY(DrivenSubstructureModel::SameEnvironment(
                    ds->EnvironmentFingerprint(model->omega[0]), model->env_fp),
                kEnvironmentDiffers);
  }
  MFEM_VERIFY(env_ports == model->ports,
              "The environment's lumped ports differ from those of the saved "
              "substructuring model!");
  MFEM_VERIFY(EnvironmentPortData(iodata, env_ports) == model->port_data,
              "The parameters of the environment's ports (the R, L and C of a lumped port, "
              "the mode of a wave port) differ from those of the saved substructuring "
              "model!");
  int nrec = Mpi::Root(comm) ? model->NumRecords(sub.save_model) : 0;
  Mpi::Broadcast(1, &nrec, 0, comm);
  MFEM_VERIFY(nrec >= static_cast<int>(model->omega.size()),
              "Truncated substructuring model \"" << sub.save_model << "\"!");
  const double unit_GHz = iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(1.0) /
                          (2 * std::numbers::pi);
  for (double omega : iodata.solver.driven.sample_f)
  {
    MFEM_VERIFY(FindValue(model->omega, omega) >= 0,
                fmt::format("Frequency {:.6g} GHz is not in the saved substructuring "
                            "model: online driven substructuring runs at saved "
                            "frequencies only!",
                            omega * unit_GHz));
  }
  // The model columns of the excitations with environment sources (-1: none), by their
  // source fingerprints (of lumped ports: exact, so to the default tolerance).
  ComplexVector rhs(fes.GetTrueVSize());
  rhs.UseDevice(true);
  for (const auto &[idx, spec] : space_op.GetPortExcitations())
  {
    space_op.GetExcitationVector(idx, model->omega[0], rhs);
    MFEM_VERIFY(!ds->OnInterface(rhs),
                "Excitation " << idx
                              << " has sources on the substructuring interface: a saved "
                                 "driven model needs the sources (surface currents, "
                                 "current dipoles) away from it!");
    const auto fp = ds->SourceFingerprint(rhs);
    int col = -1;
    for (std::size_t e = 0; e < model->excitations.size(); e++)
    {
      if (model->excitations[e] == idx &&
          DrivenSubstructureModel::SameSource(
              fp.data(), model->exc_fp.data() + DrivenSubstructureModel::kSourceFp * e))
      {
        col = static_cast<int>(e);
      }
    }
    MFEM_VERIFY(fp[0] == 0.0 || col >= 0,
                "Excitation " << idx
                              << " has environment sources that differ from the saved "
                                 "substructuring model!");
    ex_idx.push_back(idx);
    ex_col.push_back(col);
  }
  RestrictPostprocessingToRegion(space_op.GetMesh().Get(), sub, *ds, *model,
                                 !env_ports.empty(), "Online", domains, boundaries);
  if (!iodata.solver.driven.save_indices.empty())
  {
    Mpi::Warning(comm, "Online driven substructuring writes the fields in the region and on "
                       "the interface only (zero in the environment)!\n");
  }

  // The local true DOFs of Γ (in interface order, a contiguous range of it on each rank)
  // and of the region interior. The true DOFs touched by region elements other than these
  // are Dirichlet DOFs, essential for the region's operators; the others are pinned with
  // the identity.
  const int nt = fes.GetTrueVSize();
  const auto index = ds->InterfaceIndex();
  region_int = ds->RegionInterior();
  gamma_off = nG;
  for (int i = 0; i < nt; i++)
  {
    if (index[i] >= 0)
    {
      gamma_off = std::min(gamma_off, index[i]);
      gamma_local.push_back(i);
    }
  }
  for (std::size_t l = 0; l < gamma_local.size(); l++)
  {
    MFEM_VERIFY(index[gamma_local[l]] == gamma_off + static_cast<int>(l),
                "Unexpected interface order!");
  }
  gamma_cnt.resize(Mpi::Size(comm));
  gamma_disp.assign(Mpi::Size(comm), 0);
  const int m = static_cast<int>(gamma_local.size());
  MPI_Allgather(&m, 1, MPI_INT, gamma_cnt.data(), 1, MPI_INT, comm);
  for (std::size_t r = 1; r < gamma_cnt.size(); r++)
  {
    gamma_disp[r] = gamma_disp[r - 1] + gamma_cnt[r - 1];
  }
  mfem::Array<int> rm, em;
  {
    mfem::Array<int> ra(const_cast<int *>(sub.region_attributes.data()),
                        static_cast<int>(sub.region_attributes.size())),
        ea(const_cast<int *>(sub.environment_attributes.data()),
           static_cast<int>(sub.environment_attributes.size()));
    MarkInterfaceTrueDofs(const_cast<mfem::ParFiniteElementSpace &>(fes), ra, ea, rm, em);
  }
  std::vector<char> dbc(nt, 0);
  for (int i : space_op.GetNDDbcTDofLists().back())
  {
    dbc[i] = 1;
  }
  for (int i = 0; i < nt; i++)
  {
    if (dbc[i] && rm[i])
    {
      ess_tdofs.Append(i);
    }
    else if (!region_int[i] && index[i] < 0)
    {
      MFEM_VERIFY(dbc[i] || !rm[i], "A true DOF of the region is not classified!");
      identity_tdofs.Append(i);
    }
  }
}

CondensedEnvironment::~CondensedEnvironment() = default;

const std::vector<int> &CondensedEnvironment::RegionAttributes() const
{
  return sub.region_attributes;
}

void CondensedEnvironment::Load(double omega)
{
  const int j = FindValue(model->omega, omega);
  MFEM_VERIFY(j >= 0, "Frequency not in the saved substructuring model!");
  if (j == record)
  {
    return;
  }
  // The record on the current interface (as DrivenSubstructuredSweep::ModelRecord), on
  // rank 0: S_E in full (symmetric, so its row a is its column a), g, h and c.
  const bool root = Mpi::Root(comm);
  const int n = static_cast<int>(ex_idx.size()), np = static_cast<int>(env_ports.size()),
            ne = static_cast<int>(model->excitations.size());
  std::vector<std::complex<double>> S;
  g.assign(static_cast<std::size_t>(nG) * n, 0.0);
  h.assign(static_cast<std::size_t>(nG) * np, 0.0);
  c.assign(static_cast<std::size_t>(np) * ne, 0.0);
  if (root)
  {
    const auto x = model->ReadFlatRecord(sub.save_model, j);
    const auto sz = model->PartSizes();
    S = FullSchur(x.data(), nG, &gamma_map);
    const auto *gx = x.data() + sz[0], *hx = gx + sz[1], *cx = hx + sz[2];
    auto dual = [&](const std::complex<double> *y, std::complex<double> *out)
    {
      const auto z = gamma_map.DualRows(std::vector<std::complex<double>>(y, y + nG), 1);
      std::copy(z.begin(), z.end(), out);
    };
    for (int k = 0; k < n; k++)
    {
      if (ex_col[k] >= 0)
      {
        dual(gx + static_cast<std::size_t>(ex_col[k]) * nG,
             g.data() + static_cast<std::size_t>(k) * nG);
      }
    }
    for (int p = 0; p < np; p++)
    {
      dual(hx + static_cast<std::size_t>(p) * nG, h.data() + static_cast<std::size_t>(p) * nG);
    }
    std::copy(cx, cx + sz[3], c.begin());
  }
  Mpi::Broadcast(static_cast<int>(g.size()), g.data(), 0, comm);
  Mpi::Broadcast(static_cast<int>(h.size()), h.data(), 0, comm);
  Mpi::Broadcast(static_cast<int>(c.size()), c.data(), 0, comm);

  // The rows of S_E to the ranks of their interface DOFs (S_E is symmetric: row a is
  // column a).
  const int nranks = Mpi::Size(comm), m = static_cast<int>(gamma_local.size());
  std::vector<int> cnt(nranks), disp(nranks);
  for (int r = 0; r < nranks; r++)
  {
    cnt[r] = gamma_cnt[r] * nG;
    disp[r] = gamma_disp[r] * nG;
  }
  auto rows =
      std::make_shared<std::vector<std::complex<double>>>(static_cast<std::size_t>(m) * nG);
  MPI_Scatterv(S.data(), cnt.data(), disp.data(), MPI_C_DOUBLE_COMPLEX, rows->data(),
               m * nG, MPI_C_DOUBLE_COMPLEX, 0, comm);
  S_rows = std::move(rows);
  S_root = std::move(S);
  omega_record = omega;
  record = j;
}

void CondensedEnvironment::Factor()
{
  MFEM_VERIFY(record >= 0, "No S_E to factor the region against!");
  if (factored != record)
  {
    // Condense reads S_E (on rank 0), which is released after.
    ds->Condense(omega_record, std::exchange(S_root, {}));
    factored = record;
  }
}

class CondensedEnvironment::Preconditioner : public Solver<ComplexOperator>
{
public:
  explicit Preconditioner(CondensedEnvironment &env) : env(env) {}

  void SetOperator(const ComplexOperator &op) override
  {
    height = op.Height();
    width = op.Width();
    env.Factor();
  }

  void Mult(const ComplexVector &x, ComplexVector &y) const override
  {
    env.ds->SolveRegion({&x}, {&y});
  }

private:
  CondensedEnvironment &env;
};

std::unique_ptr<Solver<ComplexOperator>> CondensedEnvironment::GetPreconditioner()
{
  return std::make_unique<Preconditioner>(*this);
}

std::unique_ptr<ComplexOperator> CondensedEnvironment::GetSchurOperator(double omega)
{
  Load(omega);
  return std::make_unique<SchurOperator>(*this, S_rows);
}

void CondensedEnvironment::ApplyExcitation(int excitation_idx, double omega,
                                           ComplexVector &RHS)
{
  Load(omega);
  const int k = static_cast<int>(std::ranges::find(ex_idx, excitation_idx) - ex_idx.begin());
  MFEM_VERIFY(k < static_cast<int>(ex_idx.size()), "Unknown excitation!");
  double *br = RHS.Real().HostReadWrite(), *bi = RHS.Imag().HostReadWrite();
  for (int i = 0; i < RHS.Size(); i++)
  {
    if (!region_int[i])
    {
      br[i] = bi[i] = 0.0;
    }
  }
  const auto *gk = g.data() + static_cast<std::size_t>(k) * nG + gamma_off;
  for (std::size_t l = 0; l < gamma_local.size(); l++)
  {
    br[gamma_local[l]] = gk[l].real();
    bi[gamma_local[l]] = gk[l].imag();
  }
}

std::map<int, std::complex<double>>
CondensedEnvironment::PortVoltages(int excitation_idx, double omega, const ComplexVector &E)
{
  Load(omega);
  const int k = static_cast<int>(std::ranges::find(ex_idx, excitation_idx) - ex_idx.begin());
  MFEM_VERIFY(k < static_cast<int>(ex_idx.size()), "Unknown excitation!");
  const int np = static_cast<int>(env_ports.size());
  std::vector<std::complex<double>> V(np, 0.0);
  const double *er = E.Real().HostRead(), *ei = E.Imag().HostRead();
  for (int p = 0; p < np; p++)
  {
    const auto *hp = h.data() + static_cast<std::size_t>(p) * nG + gamma_off;
    for (std::size_t l = 0; l < gamma_local.size(); l++)
    {
      V[p] += hp[l] * std::complex<double>(er[gamma_local[l]], ei[gamma_local[l]]);
    }
  }
  Mpi::GlobalSum(np, V.data(), comm);
  std::map<int, std::complex<double>> voltages;
  for (int p = 0; p < np; p++)
  {
    voltages[env_ports[p]] =
        V[p] + ((ex_col[k] >= 0) ? c[static_cast<std::size_t>(ex_col[k]) * np + p] : 0.0);
  }
  return voltages;
}

}  // namespace palace
