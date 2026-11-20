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
    for (std::size_t j = 0; root && j < model.omega.size(); j++)
    {
      records.push_back(model.ReadFlatRecord(model_path, static_cast<int>(j)));
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
  np = static_cast<int>(env_ports.size());
  ne = static_cast<int>(model.excitations.size());

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
}

EnvironmentOperator::~EnvironmentOperator() = default;

void EnvironmentOperator::Prepare(double omega)
{
  if (omega == omega_cur)
  {
    return;
  }
  // On root, the record at ω on the current interface: dual quantities transform as
  // x_cur = M^-T x.
  std::vector<std::complex<double>> S;
  std::vector<std::vector<std::complex<double>>> g(ne);
  if (root)
  {
    std::vector<std::complex<double>> x;  // the record, flat (S_E by its lower triangle)
    if (!model.weights.empty())
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
      x = model.ReadFlatRecord(model_path, FindValue(model.omega, omega));
    }
    const auto sz = model.PartSizes();
    S = FullSchur(x.data(), nG, &gamma_map);
    const auto *gp = x.data() + sz[0], *hp = gp + sz[1], *cp = hp + sz[2];
    auto dual = [&](const std::complex<double> *y)
    { return gamma_map.DualRows(std::vector<std::complex<double>>(y, y + nG), 1); };
    for (int e = 0; e < ne; e++)
    {
      g[e] = dual(gp + static_cast<std::size_t>(e) * nG);
    }
    h.assign(np, {});
    for (int j = 0; j < np; j++)
    {
      h[j] = dual(hp + static_cast<std::size_t>(j) * nG);
    }
    c.assign(cp, cp + sz[3]);
  }

  // The rows of Γ of each process.
  const auto type = mpi::DataType<std::complex<double>>();
  const int r = Mpi::Rank(comm), m = gamma_cnt[r];
  std::vector<int> cnt_S(gamma_cnt.size()), disp_S(gamma_cnt.size());
  for (std::size_t q = 0; q < gamma_cnt.size(); q++)
  {
    cnt_S[q] = gamma_cnt[q] * nG;
    disp_S[q] = gamma_disp[q] * nG;
  }
  S_loc.resize(static_cast<std::size_t>(m) * nG);
  MPI_Scatterv(S.data(), cnt_S.data(), disp_S.data(), type, S_loc.data(), m * nG, type, 0,
               comm);
  g_loc.assign(ne, std::vector<std::complex<double>>(m));
  for (int e = 0; e < ne; e++)
  {
    MPI_Scatterv(g[e].data(), gamma_cnt.data(), gamma_disp.data(), type, g_loc[e].data(), m,
                 type, 0, comm);
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

void EnvironmentOperator::AddMult(double omega, const ComplexVector &x, ComplexVector &y)
{
  Prepare(omega);
  const auto type = mpi::DataType<std::complex<double>>();
  const auto x_loc = LocalInterfaceValues(x);
  std::vector<std::complex<double>> xG(nG);
  MPI_Allgatherv(x_loc.data(), static_cast<int>(x_loc.size()), type, xG.data(),
                 gamma_cnt.data(), gamma_disp.data(), type, comm);
  double *yr = y.Real().HostReadWrite(), *yi = y.Imag().HostReadWrite();
  for (std::size_t a = 0; a < gamma_local.size(); a++)
  {
    const auto *row = S_loc.data() + a * nG;
    std::complex<double> s = 0.0;
    for (int b = 0; b < nG; b++)
    {
      s += row[b] * xG[b];
    }
    yr[gamma_local[a]] += s.real();
    yi[gamma_local[a]] += s.imag();
  }
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
    const auto *row = S_loc.data() + a * nG;
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
  for (int p = 0; p < 2; p++)
  {
    parts[p] = std::make_unique<mfem::HypreParMatrix>(
        comm, nt, fes.GlobalTrueVSize(), fes.GlobalTrueVSize(), I.data(), J.data(),
        data[p].data(), offsets, offsets);
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
  const auto &g = g_loc[it->second];
  double *br = RHS.Real().HostReadWrite(), *bi = RHS.Imag().HostReadWrite();
  for (std::size_t a = 0; a < gamma_local.size(); a++)
  {
    br[gamma_local[a]] += g[a].real();
    bi[gamma_local[a]] += g[a].imag();
  }
  return true;
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
  const auto type = mpi::DataType<std::complex<double>>();
  const auto e_loc = LocalInterfaceValues(E);
  std::vector<std::complex<double>> eG(root ? nG : 0), values(np, 0.0);
  MPI_Gatherv(e_loc.data(), static_cast<int>(e_loc.size()), type, eG.data(),
              gamma_cnt.data(), gamma_disp.data(), type, 0, comm);
  if (root)
  {
    const auto it = ex_col.find(excitation_idx);
    for (int j = 0; j < np; j++)
    {
      std::complex<double> v =
          (it != ex_col.end()) ? c[static_cast<std::size_t>(it->second) * np + j] : 0.0;
      for (int a = 0; a < nG; a++)
      {
        v += h[j][a] * eG[a];
      }
      values[j] = v;
    }
  }
  Mpi::Broadcast(np, values.data(), 0, comm);
  for (int j = 0; j < np; j++)
  {
    const int p = env_ports[j];
    ((p > 0) ? V[p] : S[-p]) = values[j];
  }
}

}  // namespace palace
