// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "environmentoperator.hpp"

#include <algorithm>
#include <climits>
#include <cstddef>
#include <numbers>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include "models/portexcitations.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

EnvironmentOperator::EnvironmentOperator(const IoData &iodata, SpaceOperator &space_op)
  : space_op(space_op), comm(space_op.GetComm()), root(Mpi::Root(comm)),
    model_path(iodata.solver.substructuring->model_file),
    region_attrs(iodata.solver.substructuring->region_attributes)
{
  const auto &sub = *iodata.solver.substructuring;
  const double unit_GHz = iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(1.0) /
                          (2 * std::numbers::pi);

  // The mesh and the configuration of the environment are those of the model (the sources
  // of the excitations, which may be fewer, are checked one by one below), and the mesh of
  // the region is written refined.
  model.ReadHeader(model_path, comm);
  region_mesh = model.environment.CheckMesh(space_op.GetNDSpace().GetParMesh(),
                                            DrivenSubstructureModel::kEnvironmentTol);
  if (region_mesh)
  {
    const auto &refinement = iodata.model.refinement;
    MFEM_VERIFY(refinement.uniform_ref_levels == 0 &&
                    refinement.ser_uniform_ref_levels == 0 &&
                    refinement.GetBoxes().empty() && refinement.GetSpheres().empty(),
                kRegionMeshRefined);
  }
  model.environment.CheckConfig(*iodata.resolved_config, sub.environment_attributes,
                                region_mesh);
  ds = std::make_unique<DrivenSubstructure>(
      space_op, sub.region_attributes, sub.environment_attributes, true, &model,
      region_mesh, sub.direct_solver, iodata.solver.linear.sym_factorization);
  if (region_mesh)
  {
    // Γ is on the boundary of the mesh of the region, so a boundary term on it would be the
    // region's alone: in a full run it is of the side of the lower-numbered element (so
    // possibly in S_E), and a conducting sheet couples the two sides. (A model holds no
    // lumped port or surface current on Γ.)
    std::vector<int> attrs;
    for (const auto &data : iodata.boundaries.impedance)
    {
      attrs.insert(attrs.end(), data.attributes.begin(), data.attributes.end());
    }
    for (const auto &data : iodata.boundaries.conductivity)
    {
      attrs.insert(attrs.end(), data.attributes.begin(), data.attributes.end());
    }
    const std::vector<int> on_gamma = ds->InterfaceBoundaryAttributes(attrs);
    MFEM_VERIFY(on_gamma.empty(),
                "On the mesh of the region, boundary conditions on the substructuring "
                "interface (attributes "
                    << fmt::format("{}", fmt::join(on_gamma, ", "))
                    << ") would be held by the region alone: run online on the full mesh!");
  }
  nG = ds->InterfaceSize();
  Mpi::Print("\nSubstructuring: |Γ| = {:d} interface unknowns\n", nG);
  std::map<int, std::vector<char>> port_dofs;
  env_ports = EnvironmentPorts(iodata, *ds, port_dofs);

  // The pinned DOFs, and the rows of Γ of each process: its interface DOFs are a block of
  // the interface order (by owner, then local index).
  const auto &dbc = space_op.GetNDDbcTDofLists().back();
  std::vector<char> is_dbc(ds->Interface().size(), 0);
  const int *dl = dbc.HostRead();
  for (int q = 0; q < dbc.Size(); q++)
  {
    is_dbc[dl[q]] = 1;
  }
  for (std::size_t i = 0; i < ds->Interface().size(); i++)
  {
    if (ds->Interface()[i])
    {
      gamma_local.push_back(static_cast<int>(i));
    }
    else if (!ds->RegionInterior()[i] && !is_dbc[i])
    {
      pinned.Append(static_cast<int>(i));
    }
  }
  const int nranks = Mpi::Size(comm), nloc = static_cast<int>(gamma_local.size());
  gamma_cnt.assign(nranks, 0);
  gamma_disp.assign(nranks, 0);
  MPI_Allgather(&nloc, 1, MPI_INT, gamma_cnt.data(), 1, MPI_INT, comm);
  for (int r = 1; r < nranks; r++)
  {
    gamma_disp[r] = gamma_disp[r - 1] + gamma_cnt[r - 1];
  }
  MFEM_VERIFY(static_cast<long long>(nG) * nG <= INT_MAX,
              "The substructuring interface is too large to distribute its rows!");

  np = static_cast<int>(env_ports.size());
  ne = static_cast<int>(model.excitations.size());

  // The interface, environment and ports of the model.
  MFEM_VERIFY(model.nG == nG,
              "Saved substructuring model interface size ("
                  << model.nG << ") does not match this run (" << nG
                  << "); the interface (Gamma) must be identical between the offline and "
                     "online runs.");
  gamma_map = MatchSignatureBasis(
      TrueDofSignatures(space_op.GetNDSpace().Get(), ds->InterfaceIndex(), nG),
      model.signatures, model.sig_w);
  if (region_mesh)
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
  int nrec = root ? model.NumRecords(model_path) : 0;
  Mpi::Broadcast(1, &nrec, 0, comm);
  MFEM_VERIFY(nrec >= static_cast<int>(model.omega.size()),
              "Truncated substructuring model \"" << model_path << "\"!");

  // The frequencies of the sweep: anywhere between the samples of a rational model, or at
  // the frequencies of an exact one.
  const auto &omega_sample = iodata.solver.driven.sample_f;
  if (!model.weights.empty())
  {
    const auto [lo, hi] = std::ranges::minmax(model.omega);
    for (double omega : omega_sample)
    {
      MFEM_VERIFY(omega >= lo * (1.0 - 1.0e-10) && omega <= hi * (1.0 + 1.0e-10),
                  fmt::format("Frequency {:.6g} GHz is outside the band of the saved "
                              "substructuring model ({:.6g} to {:.6g} GHz)!",
                              omega * unit_GHz, lo * unit_GHz, hi * unit_GHz));
    }
    for (std::size_t j = 0; j < model.omega.size(); j++)
    {
      records.push_back(
          Distribute(root ? model.ReadFlatRecord(model_path, static_cast<int>(j))
                          : std::vector<std::complex<double>>()));
    }
    Mpi::Print(" Rational substructuring model: {:d} samples from {:.3e} to {:.3e} GHz\n",
               model.omega.size(), lo * unit_GHz, hi * unit_GHz);
    if (!model.converged)
    {
      Mpi::Warning("The saved rational substructuring model did not converge (its sampling "
                   "stopped at {:d} samples, with an error of {:.3e} for \"AdaptiveTol\" "
                   "{:.3e}): rerun the offline sweep (with a larger "
                   "\"AdaptiveMaxSamples\" if its sampling reached it)!\n",
                   model.omega.size(), model.error, model.tol);
    }
    if (model.passivity < -model.tol)
    {
      Mpi::Warning("The saved rational substructuring model is not passive between its "
                   "samples (least eigenvalue of Im S_E / ‖S_E‖ = {:.3e}): rerun the "
                   "offline sweep with a larger \"AdaptiveMaxSamples\"!\n",
                   model.passivity);
    }
  }
  else
  {
    for (double omega : omega_sample)
    {
      MFEM_VERIFY(FindValue(model.omega, omega) >= 0,
                  fmt::format("Frequency {:.6g} GHz is not in the saved substructuring "
                              "model: online driven substructuring runs at saved "
                              "frequencies only (or anywhere in the band of a model saved "
                              "by an adaptive sweep)!",
                              omega * unit_GHz));
    }
  }

  // The model columns of the excitations with environment sources, which must lie away
  // from Γ (as the lumped ports do).
  const auto &port_excitations = space_op.GetPortExcitations();
  ComplexVector rhs(space_op.GetNDSpace().GetTrueVSize());
  rhs.UseDevice(true);
  for (const auto &[idx, spec] : port_excitations)
  {
    space_op.GetExcitationVector(idx, model.omega[0], rhs);
    MFEM_VERIFY(!ds->OnInterface(rhs),
                "Excitation " << idx
                              << " has sources on the substructuring interface: a saved "
                                 "driven model, or an adaptive sweep, needs the sources "
                                 "(surface currents, current dipoles) away from it!");
    const auto fp = ds->SourceFingerprint(rhs);
    const double tol = SourceTolerance(iodata, port_excitations, env_ports, idx);
    int col = -1;
    for (std::size_t e = 0; e < model.excitations.size(); e++)
    {
      if (model.excitations[e] == idx &&
          DrivenSubstructureModel::SameSource(
              fp.data(), model.exc_fp.data() + DrivenSubstructureModel::kSourceFp * e, tol))
      {
        col = static_cast<int>(e);
      }
    }
    MFEM_VERIFY(fp[0] == 0.0 || col >= 0,
                "Excitation " << idx
                              << " has environment sources that differ from the saved "
                                 "substructuring model!");
    if (col >= 0)
    {
      ex_col[idx] = col;
    }
  }
}

EnvironmentOperator::~EnvironmentOperator() = default;

EnvironmentOperator::Rows
EnvironmentOperator::Distribute(const std::vector<std::complex<double>> &flat) const
{
  // On root, the record on the current interface: dual quantities transform as
  // x_cur = M^-T x.
  std::vector<std::complex<double>> S, all_g, all_h, c;
  if (root)
  {
    const auto sz = model.PartSizes();
    S = FullSchur(flat.data(), nG, &gamma_map);
    const auto *gp = flat.data() + sz[0], *hp = gp + sz[1], *cp = hp + sz[2];
    auto dual = [&](const std::complex<double> *y)
    { return gamma_map.DualRows(std::vector<std::complex<double>>(y, y + nG), 1); };
    for (int e = 0; e < ne; e++)
    {
      const auto ge = dual(gp + static_cast<std::size_t>(e) * nG);
      all_g.insert(all_g.end(), ge.begin(), ge.end());
    }
    for (int j = 0; j < np; j++)
    {
      const auto hj = dual(hp + static_cast<std::size_t>(j) * nG);
      all_h.insert(all_h.end(), hj.begin(), hj.end());
    }
    c.assign(cp, cp + sz[3]);
  }

  // The rows of Γ of each process.
  const auto type = mpi::DataType<std::complex<double>>();
  const int m = gamma_cnt[Mpi::Rank(comm)];
  auto scatter = [&](const std::complex<double> *x, int width)
  {
    std::vector<int> cnt(gamma_cnt.size()), disp(gamma_cnt.size());
    for (std::size_t q = 0; q < gamma_cnt.size(); q++)
    {
      cnt[q] = gamma_cnt[q] * width;
      disp[q] = gamma_disp[q] * width;
    }
    std::vector<std::complex<double>> y(static_cast<std::size_t>(m) * width);
    MPI_Scatterv(x, cnt.data(), disp.data(), type, y.data(), m * width, type, 0, comm);
    return y;
  };
  Rows rows;
  rows.S = scatter(S.data(), nG);
  for (int e = 0; e < ne; e++)
  {
    rows.g.push_back(
        scatter(root ? all_g.data() + static_cast<std::size_t>(e) * nG : nullptr, 1));
  }
  for (int j = 0; j < np; j++)
  {
    rows.h.push_back(
        scatter(root ? all_h.data() + static_cast<std::size_t>(j) * nG : nullptr, 1));
  }
  rows.c.resize(static_cast<std::size_t>(np) * ne);
  if (root)
  {
    rows.c = c;
  }
  Mpi::Broadcast(static_cast<int>(rows.c.size()), rows.c.data(), 0, comm);
  return rows;
}

std::vector<std::complex<double>> EnvironmentOperator::Coefficients(double omega) const
{
  return BarycentricInterpolant::Coefficients(model.omega, model.weights, omega);
}

void EnvironmentOperator::Prepare(double omega)
{
  if (omega == omega_cur)
  {
    return;
  }
  if (records.empty())
  {
    cur = Distribute(root ? model.ReadFlatRecord(model_path, FindValue(model.omega, omega))
                          : std::vector<std::complex<double>>());
  }
  else
  {
    // The combination of the records with their coefficients at ω.
    const auto a = Coefficients(omega);
    auto combine = [&](auto part)
    {
      std::vector<std::complex<double>> x(part(records[0]).size(), 0.0);
      for (std::size_t k = 0; k < records.size(); k++)
      {
        const auto &y = part(records[k]);
        for (std::size_t q = 0; q < x.size(); q++)
        {
          x[q] += a[k] * y[q];
        }
      }
      return x;
    };
    cur.S = combine([](const Rows &r) -> const auto & { return r.S; });
    cur.c = combine([](const Rows &r) -> const auto & { return r.c; });
    cur.g.resize(ne);
    for (int e = 0; e < ne; e++)
    {
      cur.g[e] = combine([e](const Rows &r) -> const auto & { return r.g[e]; });
    }
    cur.h.resize(np);
    for (int j = 0; j < np; j++)
    {
      cur.h[j] = combine([j](const Rows &r) -> const auto & { return r.h[j]; });
    }
  }
  omega_cur = omega;
}

std::vector<std::complex<double>>
EnvironmentOperator::LocalInterfaceValues(const ComplexVector &x) const
{
  const double *xr = x.Real().HostRead(), *xi = x.Imag().HostRead();
  std::vector<std::complex<double>> v(gamma_local.size());
  for (std::size_t a = 0; a < gamma_local.size(); a++)
  {
    v[a] = {xr[gamma_local[a]], xi[gamma_local[a]]};
  }
  return v;
}

void EnvironmentOperator::AddRows(const std::vector<std::complex<double>> &rows,
                                  const ComplexVector &x, ComplexVector &y) const
{
  const auto type = mpi::DataType<std::complex<double>>();
  const auto x_loc = LocalInterfaceValues(x);
  std::vector<std::complex<double>> xG(nG);
  MPI_Allgatherv(x_loc.data(), static_cast<int>(x_loc.size()), type, xG.data(),
                 gamma_cnt.data(), gamma_disp.data(), type, comm);
  double *yr = y.Real().HostReadWrite(), *yi = y.Imag().HostReadWrite();
  for (std::size_t a = 0; a < gamma_local.size(); a++)
  {
    const auto *row = rows.data() + a * nG;
    std::complex<double> v = 0.0;
    for (int b = 0; b < nG; b++)
    {
      v += row[b] * xG[b];
    }
    yr[gamma_local[a]] += v.real();
    yi[gamma_local[a]] += v.imag();
  }
}

void EnvironmentOperator::AddMult(double omega, const ComplexVector &x, ComplexVector &y)
{
  Prepare(omega);
  AddRows(cur.S, x, y);
}

void EnvironmentOperator::RecordAddMult(std::size_t k, const ComplexVector &x,
                                        ComplexVector &y) const
{
  AddRows(records[k].S, x, y);
}

std::array<std::unique_ptr<mfem::HypreParMatrix>, 2>
EnvironmentOperator::Assemble(double omega)
{
  Prepare(omega);
  const auto &fes = space_op.GetNDSpace().Get();
  const int nt = fes.GetTrueVSize();
  const auto &cols = ds->InterfaceTrueDofs();
  // The rows of Γ, and a zero diagonal entry in every other row, so that the eliminated
  // rows of a sum with this matrix can hold their diagonal.
  std::vector<int> I(nt + 1, 1);
  I[0] = 0;
  for (int i : gamma_local)
  {
    I[i + 1] = nG;
  }
  for (int i = 0; i < nt; i++)
  {
    I[i + 1] += I[i];
  }
  std::vector<HYPRE_BigInt> J(static_cast<std::size_t>(I[nt]));
  std::vector<double> data[2] = {std::vector<double>(J.size(), 0.0),
                                 std::vector<double>(J.size(), 0.0)};
  const HYPRE_BigInt tstart = fes.GetMyTDofOffset();
  for (int i = 0; i < nt; i++)
  {
    J[I[i]] = tstart + i;
  }
  for (std::size_t a = 0; a < gamma_local.size(); a++)
  {
    const auto *row = cur.S.data() + a * nG;
    const std::size_t at = I[gamma_local[a]];
    for (int b = 0; b < nG; b++)
    {
      J[at + b] = cols[b];
      data[0][at + b] = row[b].real();
      data[1][at + b] = row[b].imag();
    }
  }
  HYPRE_BigInt *offsets =
      const_cast<mfem::ParFiniteElementSpace &>(fes).GetTrueDofOffsets();
  std::array<std::unique_ptr<mfem::HypreParMatrix>, 2> parts;
  for (int q = 0; q < 2; q++)
  {
    parts[q] = std::make_unique<mfem::HypreParMatrix>(
        comm, nt, fes.GlobalTrueVSize(), fes.GlobalTrueVSize(), I.data(), J.data(),
        data[q].data(), offsets, offsets);
  }
  return parts;
}

bool EnvironmentOperator::AddExcitation(int excitation_idx, double omega,
                                        ComplexVector &RHS)
{
  const auto it = ex_col.find(excitation_idx);
  if (it == ex_col.end())
  {
    return false;
  }
  Prepare(omega);
  const auto &g = cur.g[it->second];
  double *br = RHS.Real().HostReadWrite(), *bi = RHS.Imag().HostReadWrite();
  for (std::size_t a = 0; a < gamma_local.size(); a++)
  {
    br[gamma_local[a]] += g[a].real();
    bi[gamma_local[a]] += g[a].imag();
  }
  return true;
}

std::vector<std::complex<double>>
EnvironmentOperator::AllInterfaceValues(const Vector &x) const
{
  const double *px = x.HostRead();
  std::vector<double> x_loc(gamma_local.size()), xG(nG);
  for (std::size_t a = 0; a < gamma_local.size(); a++)
  {
    x_loc[a] = px[gamma_local[a]];
  }
  MPI_Allgatherv(x_loc.data(), static_cast<int>(x_loc.size()), MPI_DOUBLE, xG.data(),
                 gamma_cnt.data(), gamma_disp.data(), MPI_DOUBLE, comm);
  return {xG.begin(), xG.end()};
}

void EnvironmentOperator::RecordAddMultPart(std::size_t k, bool imag, const Vector &x,
                                            Vector &y) const
{
  const auto xG = AllInterfaceValues(x);
  double *py = y.HostReadWrite();
  for (std::size_t a = 0; a < gamma_local.size(); a++)
  {
    const auto *row = records[k].S.data() + a * nG;
    double v = 0.0;
    for (int b = 0; b < nG; b++)
    {
      v += (imag ? row[b].imag() : row[b].real()) * xG[b].real();
    }
    py[gamma_local[a]] += v;
  }
}

std::vector<std::complex<double>>
EnvironmentOperator::RecordExcitationPairings(int excitation_idx, const Vector &v) const
{
  const auto it = ex_col.find(excitation_idx);
  if (it == ex_col.end())
  {
    return {};
  }
  const double *pv = v.HostRead();
  std::vector<std::complex<double>> p(records.size(), 0.0);
  for (std::size_t k = 0; k < records.size(); k++)
  {
    const auto &g = records[k].g[it->second];
    for (std::size_t a = 0; a < gamma_local.size(); a++)
    {
      p[k] += pv[gamma_local[a]] * g[a];
    }
  }
  Mpi::GlobalSum(static_cast<int>(p.size()), p.data(), comm);
  return p;
}

void EnvironmentOperator::PortValues(int excitation_idx, double omega,
                                     const ComplexVector &E,
                                     std::map<int, std::complex<double>> &V,
                                     std::map<int, std::complex<double>> &S)
{
  V.clear();
  S.clear();
  if (np == 0)
  {
    return;
  }
  Prepare(omega);
  // V_j = c_jk + h_j^T u_Γ, from the rows of Γ of each process.
  const auto e_loc = LocalInterfaceValues(E);
  std::vector<std::complex<double>> values(np, 0.0);
  for (int j = 0; j < np; j++)
  {
    for (std::size_t a = 0; a < e_loc.size(); a++)
    {
      values[j] += cur.h[j][a] * e_loc[a];
    }
  }
  Mpi::GlobalSum(np, values.data(), comm);
  const auto it = ex_col.find(excitation_idx);
  for (int j = 0; j < np; j++)
  {
    if (it != ex_col.end())
    {
      values[j] += cur.c[static_cast<std::size_t>(it->second) * np + j];
    }
    const int p = env_ports[j];
    ((p > 0) ? V[p] : S[-p]) = values[j];
  }
}

}  // namespace palace
