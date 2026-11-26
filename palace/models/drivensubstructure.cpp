// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "drivensubstructure.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <span>
#include <type_traits>
#include <utility>
#include <Eigen/Dense>
#include "fem/interpolator.hpp"
#include "fem/substructure.hpp"
#include "fem/substructuremodel.hpp"
#include "linalg/blockcyclic.hpp"
#include "linalg/lowertriangle.hpp"
#include "linalg/orthog.hpp"
#include "linalg/rap.hpp"
#include "linalg/schurfactory.hpp"
#include "models/portexcitations.hpp"
#include "models/romoperator.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

namespace
{

// Free a vector's storage (v = {} empties it but keeps its capacity).
template <typename V>
void Release(V &v)
{
  V().swap(v);
}

// The assembled real (imag = false) or imaginary part of an operator, if any, owned.
std::unique_ptr<mfem::HypreParMatrix> StealPart(const ComplexOperator *A, bool imag)
{
  const auto *part =
      A ? dynamic_cast<const ParOperator *>(imag ? A->Imag() : A->Real()) : nullptr;
  return part ? part->StealParallelAssemble() : nullptr;
}

}  // namespace

DrivenSubstructure::DrivenSubstructure(SpaceOperator &space_op,
                                       const std::vector<int> &region_attributes,
                                       const std::vector<int> &environment_attributes,
                                       LinearSolver solver, SymbolicFactorization ordering)
  : DrivenSubstructure(space_op, region_attributes, environment_attributes, nullptr, false,
                       solver, ordering)
{
}

DrivenSubstructure::DrivenSubstructure(SpaceOperator &space_op,
                                       const std::vector<int> &region_attributes,
                                       const std::vector<int> &environment_attributes,
                                       const DrivenSubstructureModel &model,
                                       bool region_mesh)
  : DrivenSubstructure(space_op, region_attributes, environment_attributes, &model,
                       region_mesh, LinearSolver::DEFAULT, SymbolicFactorization::DEFAULT)
{
}

DrivenSubstructure::DrivenSubstructure(SpaceOperator &space_op,
                                       const std::vector<int> &region_attributes,
                                       const std::vector<int> &environment_attributes,
                                       const DrivenSubstructureModel *saved,
                                       bool region_mesh, LinearSolver solver,
                                       SymbolicFactorization ordering)
  : space_op(space_op), solver(solver), ordering(ordering), online(saved != nullptr),
    region_attrs(region_attributes)
{
  env.attrs = environment_attributes;
  const auto &fes = space_op.GetNDSpace().Get();
  const auto &mesh = *fes.GetParMesh();
  CheckSubstructuringAttributes(mesh, region_attrs, env.attrs);
  fp_centers =
      saved ? saved->fp_centers
            : palace::FingerprintCenters(space_op.GetNDSpace().GetParMesh(), env.attrs);

  // Interface DOFs: on region and environment elements; environment (region) interior: on
  // environment (region) elements only; all without the Dirichlet DOFs.
  mfem::Array<int> rm, em;
  {
    mfem::Array<int> ra(region_attrs.data(), static_cast<int>(region_attrs.size())),
        ea(env.attrs.data(), static_cast<int>(env.attrs.size()));
    MarkInterfaceTrueDofs(const_cast<mfem::ParFiniteElementSpace &>(fes), ra, ea, rm, em);
  }
  const int nt = fes.GetTrueVSize();
  std::vector<char> dbc(nt, 0);
  const auto &dbc_list = space_op.GetNDDbcTDofLists().back();
  const int *dl = dbc_list.HostRead();
  for (int q = 0; q < dbc_list.Size(); q++)
  {
    dbc[dl[q]] = 1;
  }
  MPI_Comm comm = fes.GetComm();

  // A mesh of the region (online): Γ also holds the region's boundary DOFs at the points
  // of the saved interface DOFs, where the environment was cut off.
  std::vector<char> cut(nt, 0);
  // With the saved model, the mesh holds its environment or the part of it in the mesh of
  // the region, which the sweep tells by the model's record.
  region_only = saved && region_mesh;
  if (region_only)
  {
    std::vector<char> region_side(nt);
    for (int i = 0; i < nt; i++)
    {
      region_side[i] = rm[i] && !em[i] && !dbc[i];
    }
    cut = RegionMeshInterface(const_cast<mfem::ParFiniteElementSpace &>(fes), region_side,
                              {}, saved->signatures, saved->sig_w);
  }
  iface = SubstructureInterface(ClassifyInterfaceTrueDofs(rm, em, dbc, cut),
                                fes.GetMyTDofOffset(), comm);
  for (int i = 0; i < nt; i++)
  {
    if (!iface.gamma[i] && !iface.env_int[i])
    {
      env.pinned.Append(i);
    }
  }

  // Wave ports: each inside the region or the environment, away from Γ, so that its modal
  // terms (if active) border the environment's system or are the region's, and its overlap
  // is a functional of that side's field.
  for (const auto &[idx, data] : space_op.GetWavePortOp())
  {
    const auto &list = data.GetAttrList();
    auto mark = BoundaryTrueDofs(std::vector<int>(list.begin(), list.end()));
    double n[3] = {0.0, 0.0, 0.0};  // environment interior, interface, region interior
    for (int i = 0; i < nt; i++)
    {
      n[0] += mark[i] && iface.env_int[i];
      n[1] += mark[i] && iface.gamma[i];
      n[2] += mark[i] && iface.region_int[i];
    }
    Mpi::GlobalSum(3, n, comm);
    MFEM_VERIFY(
        n[1] == 0.0 && (n[0] == 0.0 || n[2] == 0.0),
        "Wave port " << idx
                     << " touches the substructuring interface: wave ports must lie "
                        "inside the region or the environment!");
    if (n[0] > 0.0)
    {
      env.wave_ports.push_back(idx);
    }
    wave_port_dofs[idx] = std::move(mark);
  }
}

std::vector<int> DrivenSubstructure::EnvironmentWavePorts() const
{
  return env.wave_ports;
}

std::vector<std::pair<int, WavePortOperator::ModalCorrectionTerm>>
DrivenSubstructure::ModalTerms(const Side &side, double omega)
{
  // The terms of the side's active wave ports (each term carries its port).
  std::vector<std::pair<int, WavePortOperator::ModalCorrectionTerm>> terms;
  if (side.wave_ports.empty())
  {
    return terms;
  }
  for (auto &term : space_op.GetWavePortOp().GetModalCorrectionTerms(
           omega, space_op.GetNDSpace(), space_op.GetNDDbcTDofLists().back()))
  {
    const int port = term.port;
    if (std::ranges::find(side.wave_ports, port) != side.wave_ports.end())
    {
      terms.emplace_back(port, std::move(term));
    }
  }
  return terms;
}

DrivenSubstructure::~DrivenSubstructure() = default;

std::vector<std::unique_ptr<mfem::HypreParMatrix>>
DrivenSubstructure::ExtraParts(const Side &side, double omega)
{
  std::unique_ptr<ComplexOperator> A2;
  {
    AssemblyDomains::Restriction restriction(space_op.GetAssemblyDomains(), side.attrs);
    A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO);
  }
  std::vector<std::unique_ptr<mfem::HypreParMatrix>> X(2);
  for (int p = 0; p < 2; p++)
  {
    X[p] = StealPart(A2.get(), p == 1);
    if (X[p])
    {
      X[p]->EliminateBC(side.pinned, Operator::DIAG_ZERO);
    }
  }
  return X;
}

namespace
{

// A(ω) = K + iω C - ω² M + A2(ω): the coefficients of the parts Kr, Ki, Cr, Ci, Mr, Mi, and
// of A2r, A2i.
std::array<std::complex<double>, 8> Coefficients(double omega)
{
  const double w2 = omega * omega;
  return {1.0, {0.0, 1.0}, {0.0, omega}, -omega, -w2, {0.0, -w2}, 1.0, {0.0, 1.0}};
}

}  // namespace

void DrivenSubstructure::Setup(Side &side, double omega)
{
#if defined(MFEM_USE_MUMPS)
  static_assert(std::is_same_v<MUMPS_INT, int>, "MUMPS_INT must be a 32-bit int!");
#endif
  // The factored system is the side's interior and Γ, numbered by rank (see SideSystem):
  // the pinned DOFs would be identity rows, and MUMPS's per-rank workspace grows with the
  // system size.
  const auto &fes = space_op.GetNDSpace().Get();
  MPI_Comm comm = fes.GetComm();
  const int nt = fes.GetTrueVSize();
  std::vector<char> keep(nt, 1);
  for (int i : side.pinned)
  {
    keep[i] = 0;
  }

  // The frequency-independent parts, assembled once and one at a time: the pattern of the
  // stiffness matrix on the system, which contains those of the other parts, coupling DOFs
  // of elements of the side.
  side.system.reset();
  for (int op = 0; op < 3; op++)
  {
    auto parts = SideParts(side, op);
    for (int p = 0; p < 2; p++)
    {
      const auto &X = parts[p];
      if (op == 0 && p == 0)
      {
        MFEM_VERIFY(X, "Missing stiffness matrix of a substructure operator!");
        side.system = std::make_unique<SideSystem>(
            std::vector<const mfem::HypreParMatrix *>{X.get()}, keep);
        side.system->parts.resize(6);
      }
      else if (X)
      {
        auto &part = side.system->parts[2 * op + p];
        part.assign(side.system->parts[0].size(), 0.0);
        side.system->AddPart(*X, 1.0, part);
      }
    }
  }
  auto &sys = *side.system;

  // The modal terms of the side's wave ports border the system, with an extra unknown per
  // term on the last rank, against the DOFs of the term's port. Its diagonal is -σ and its
  // entries √(gσ) s, with σ the largest diagonal entry of the side's operator at this
  // frequency, so that the bordering does not change the scale of the factored matrix.
  const auto terms = ModalTerms(side, omega);
  {
    const auto coef = Coefficients(omega);
    side.border = 0.0;
    for (std::size_t q = 0; q < sys.irn.size(); q++)
    {
      if (sys.irn[q] == sys.jcn[q])
      {
        std::complex<double> d = 0.0;
        for (int p = 0; p < 6; p++)
        {
          d += sys.parts[p].empty() ? 0.0 : coef[p] * sys.parts[p][q];
        }
        side.border = std::max(side.border, std::abs(d));
      }
    }
    Mpi::GlobalMax(1, &side.border, comm);
    side.border = (side.border > 0.0) ? side.border : 1.0;
  }
  const int nterm = static_cast<int>(terms.size());
  const bool last = (Mpi::Rank(comm) == Mpi::Size(comm) - 1);
  MFEM_VERIFY(sys.n_glob + nterm <= std::numeric_limits<int>::max(),
              "The system of the environment needs fewer than 2^31 unknowns (32-bit "
              "indices)!");
  side.term_ports.clear();
  side.term_dofs.assign(nterm, {});
  for (int t = 0; t < nterm; t++)
  {
    side.term_ports.push_back(terms[t].first);
    const auto &mark = wave_port_dofs.at(terms[t].first);
    const int row = static_cast<int>(sys.n_glob + t + 1);
    for (int i = 0; i < nt; i++)
    {
      if (mark[i] && sys.Index(i) >= 0)
      {
        side.term_dofs[t].push_back(i);
        sys.irn.push_back(row);
        sys.jcn.push_back(static_cast<int>(sys.Index(i) + 1));
      }
    }
    if (last)
    {
      sys.irn.push_back(row);
      sys.jcn.push_back(row);
      sys.rows.push_back(-1);
    }
  }
  side.val.resize(sys.irn.size());
  auto extra = ExtraParts(side, omega);
  side.extra = extra[0] || extra[1];
  Fill(side, omega, extra, terms, side.val);
}

std::array<std::unique_ptr<mfem::HypreParMatrix>, 2>
DrivenSubstructure::SideParts(const Side &side, int op) const
{
  std::unique_ptr<ComplexOperator> A;
  {
    AssemblyDomains::Restriction restriction(space_op.GetAssemblyDomains(), side.attrs);
    A = (op == 0)   ? space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO)
        : (op == 1) ? space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO)
                    : space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  }
  std::array<std::unique_ptr<mfem::HypreParMatrix>, 2> X;
  for (int p = 0; p < 2; p++)
  {
    X[p] = StealPart(A.get(), p == 1);
    if (X[p])
    {
      X[p]->EliminateBC(side.pinned, Operator::DIAG_ZERO);
    }
  }
  return X;
}

void DrivenSubstructure::Fill(
    const Side &side, double omega,
    const std::vector<std::unique_ptr<mfem::HypreParMatrix>> &extra,
    const std::vector<std::pair<int, WavePortOperator::ModalCorrectionTerm>> &terms,
    std::vector<std::complex<double>> &val) const
{
  const auto coef = Coefficients(omega);
  const auto &parts = side.system->parts;
  std::fill(val.begin(), val.end(), 0.0);
  for (int p = 0; p < 6; p++)
  {
    for (std::size_t q = 0; q < parts[p].size(); q++)
    {
      val[q] += coef[p] * parts[p][q];
    }
  }
  for (int p = 0; p < 2; p++)
  {
    if (extra[p])
    {
      side.system->AddPart(*extra[p], coef[6 + p], val);
    }
  }
  MFEM_VERIFY(terms.size() == side.term_dofs.size(),
              "The modal terms of the wave ports changed between frequencies!");
  const bool last = (Mpi::Rank(space_op.GetComm()) == Mpi::Size(space_op.GetComm()) - 1);
  std::size_t q = parts[0].size();
  for (std::size_t t = 0; t < terms.size(); t++)
  {
    MFEM_VERIFY(terms[t].first == side.term_ports[t],
                "The modal terms of the wave ports changed between frequencies!");
    const std::complex<double> c = std::sqrt(terms[t].second.g * side.border);
    const double *sr = terms[t].second.s->Real().HostRead(),
                 *si = terms[t].second.s->Imag().HostRead();
    for (int i : side.term_dofs[t])
    {
      val[q++] = c * std::complex<double>(sr[i], si[i]);
    }
    if (last)
    {
      val[q++] = -side.border;
    }
  }
}

void DrivenSubstructure::Factor(Side &side, double omega)
{
  if (!side.schur)
  {
    MPI_Comm comm = space_op.GetComm();
    gamma_grid = std::make_unique<BlockCyclicGrid>(comm, iface.size, Mpi::Size(comm));
    auto &sys = *side.system;
    const HYPRE_BigInt n_sys =
        sys.n_glob + static_cast<HYPRE_BigInt>(side.term_dofs.size());
    const int n_loc = static_cast<int>(sys.rows.size());
    side.schur = MakeSchurFactorization<std::complex<double>>(
        SchurSolverFor(solver, SchurKind::ENVIRONMENT), ordering, comm, n_sys, n_loc,
        {std::move(sys.irn), std::move(sys.jcn), std::move(side.val)},
        sys.Indices(iface.tdofs), sys.rows, false, true, gamma_grid.get());
    return;
  }

  // The entries at a new frequency, in place.
  auto extra = ExtraParts(side, omega);
  MFEM_VERIFY(side.extra == (extra[0] || extra[1]),
              "The frequency-dependent part of a substructure operator changed!");
  Fill(side, omega, extra, ModalTerms(side, omega), side.schur->Values());
  side.schur->Refactor();
}

MPI_Comm DrivenSubstructure::GetComm() const
{
  return space_op.GetComm();
}

void DrivenSubstructure::Condense(double omega)
{
  MFEM_VERIFY(!online, "An online substructure has no environment to condense!");
  if (!env.schur)
  {
    // The pattern first, so that the assembled operators are released before the
    // factorization.
    Setup(env, omega);
  }
  // Each process takes its share of S_E from the grid of the Schur complement, which MUMPS
  // hands over (and allocates again at the next factorization), or which the others
  // scatter from rank 0.
  Release(S);
  Factor(env, omega);
  S = ShareFromGrid(*gamma_grid, env.schur->ReleaseSchurBlocks(*gamma_grid));
}

namespace
{

// b on the masked local true DOFs (and on Γ with gamma), 0 elsewhere.
ComplexVector Masked(const ComplexVector &b, const std::vector<char> &mask,
                     const std::vector<char> *gamma = nullptr)
{
  const int nt = static_cast<int>(mask.size());
  ComplexVector y(nt);
  const double *br = b.Real().HostRead(), *bi = b.Imag().HostRead();
  double *yr = y.Real().HostWrite(), *yi = y.Imag().HostWrite();
  for (int i = 0; i < nt; i++)
  {
    const bool keep = mask[i] || (gamma && (*gamma)[i]);
    yr[i] = keep ? br[i] : 0.0;
    yi[i] = keep ? bi[i] : 0.0;
  }
  return y;
}

}  // namespace

std::vector<std::complex<double>>
DrivenSubstructure::CondenseEnvironment(const std::vector<const ComplexVector *> &b,
                                        std::vector<ComplexVector> *xo)
{
  MFEM_VERIFY(!online && env.schur, "CondenseEnvironment needs the environment factor!");
  std::vector<std::complex<double>> red;
  if (b.empty())
  {
    return red;  // (MUMPS rejects an empty batch)
  }
  std::vector<ComplexVector> x(b.size());
  std::vector<const ComplexVector *> X(b.size());
  std::vector<ComplexVector *> Xo(b.size());
  for (std::size_t k = 0; k < b.size(); k++)
  {
    x[k] = Masked(*b[k], iface.env_int, &iface.gamma);
    X[k] = &x[k];
    Xo[k] = &x[k];
  }
  env.schur->Reduce(X, red);
  // The reduction is paired with an expansion, of a zero interface solution: the interior
  // solutions.
  std::vector<std::complex<double>> zero(Mpi::Root(space_op.GetComm()) ? red.size() : 0,
                                         0.0);
  env.schur->Expand(zero, Xo);
  if (xo)
  {
    *xo = std::move(x);
  }
  return red;
}

void DrivenSubstructure::ExpandEnvironment(const ComplexVector &b, ComplexVector &E)
{
  MFEM_VERIFY(!online && env.schur, "ExpandEnvironment needs the environment factor!");
  // The reduction of b paired with the expansion of the values of E on Γ, gathered on root
  // in interface order.
  ComplexVector x = Masked(b, iface.env_int, &iface.gamma);
  std::vector<std::complex<double>> red;
  env.schur->Reduce({&x}, red);
  const int nt = static_cast<int>(iface.gamma.size());
  const double *er = E.Real().HostRead(), *ei = E.Imag().HostRead();
  std::vector<std::complex<double>> u_gamma;
  u_gamma.reserve(iface.count);
  for (int i = 0; i < nt; i++)
  {
    if (iface.gamma[i])
    {
      u_gamma.emplace_back(er[i], ei[i]);
    }
  }
  env.schur->Expand(iface.GatherRows(space_op.GetComm(), u_gamma, 1), {&x});
  const double *xr = x.Real().HostRead(), *xi = x.Imag().HostRead();
  double *Er = E.Real().HostReadWrite(), *Ei = E.Imag().HostReadWrite();
  for (int i = 0; i < nt; i++)
  {
    if (iface.env_int[i])
    {
      Er[i] = xr[i];
      Ei[i] = xi[i];
    }
  }
}

std::vector<char>
DrivenSubstructure::BoundaryTrueDofs(const std::vector<int> &bdr_attributes) const
{
  return palace::BoundaryTrueDofs(
      const_cast<mfem::ParFiniteElementSpace &>(space_op.GetNDSpace().Get()),
      bdr_attributes);
}

std::vector<int> DrivenSubstructure::InterfaceBoundaryAttributes(
    const std::vector<int> &bdr_attributes) const
{
  // On the local DOFs, also those other ranks own: 2 on Γ, 1 on the Dirichlet boundary
  // (which Γ excludes, and which a face on Γ can touch), else 0.
  auto &fes = const_cast<mfem::ParFiniteElementSpace &>(space_op.GetNDSpace().Get());
  const auto &mesh = *fes.GetParMesh();
  mfem::Vector t(fes.GetTrueVSize());
  t = 0.0;
  const auto &dbc_list = space_op.GetNDDbcTDofLists().back();
  const int *dl = dbc_list.HostRead();
  for (int q = 0; q < dbc_list.Size(); q++)
  {
    t(dl[q]) = 1.0;
  }
  for (int i = 0; i < t.Size(); i++)
  {
    if (iface.gamma[i])
    {
      t(i) = 2.0;
    }
  }
  mfem::ParGridFunction g(&fes);
  g.SetFromTrueDofs(t);
  const int amax = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  std::vector<int> found(amax + 1, 0);
  mfem::Array<int> dofs;
  const double *h_g = g.HostRead();
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    const int a = mesh.GetBdrAttribute(be);
    if (found[a] || std::ranges::find(bdr_attributes, a) == bdr_attributes.end())
    {
      continue;
    }
    fes.GetBdrElementDofs(be, dofs);
    bool all = true, any = false;
    for (int d : dofs)
    {
      const double v = h_g[(d >= 0) ? d : -1 - d];
      all = all && v > 0.5;
      any = any || v > 1.5;
    }
    found[a] = all && any;
  }
  MPI_Allreduce(MPI_IN_PLACE, found.data(), amax + 1, MPI_INT, MPI_MAX, fes.GetComm());
  std::vector<int> attrs;
  for (int a = 1; a <= amax; a++)
  {
    if (found[a])
    {
      attrs.push_back(a);
    }
  }
  return attrs;
}

std::vector<char>
DrivenSubstructure::WriteRegionMesh(const std::string &path,
                                    const std::vector<int> &port_attributes,
                                    double length_scale) const
{
  // The region, and the environment elements touching the DOFs of the environment's ports
  // or the region's Dirichlet DOFs (on Γ, set by boundary conditions of the environment).
  auto &fes = const_cast<mfem::ParFiniteElementSpace &>(space_op.GetNDSpace().Get());
  if (fes.GetParMesh()->Nonconforming())
  {
    // As for the static substructuring: the submesh of a nonconforming mesh also selects by
    // the attributes of the parents of refined elements, which the written environment
    // elements do not shift.
    Mpi::Warning(fes.GetComm(),
                 "No mesh of the region is written for a nonconforming mesh: online runs "
                 "need the full mesh!\n");
    return std::vector<char>(fes.GetParMesh()->GetNE(), 0);
  }
  std::vector<char> keep = BoundaryTrueDofs(port_attributes);
  mfem::Array<int> rm, em, ra, ea;
  for (int a : region_attrs)
  {
    ra.Append(a);
  }
  for (int a : env.attrs)
  {
    ea.Append(a);
  }
  MarkInterfaceTrueDofs(fes, ra, ea, rm, em);
  const auto &dbc_list = space_op.GetNDDbcTDofLists().back();
  const int *dl = dbc_list.HostRead();
  for (int q = 0; q < dbc_list.Size(); q++)
  {
    keep[dl[q]] = keep[dl[q]] || rm[dl[q]];
  }
  std::vector<char> written =
      palace::WriteRegionMesh(path, fes, region_attrs, env.attrs, keep, length_scale);
  Mpi::Print(" Writing the mesh of the region for online runs to {}\n", path);
  return written;
}

const Vector &DrivenSubstructure::FingerprintField(int k) const
{
  if (!fp_fields[k])
  {
    fp_fields[k] = palace::FingerprintField(space_op.GetNDSpace().Get(), k, fp_centers[k]);
  }
  return *fp_fields[k];
}

std::vector<double> DrivenSubstructure::EnvironmentFingerprint(double omega) const
{
  const auto &fes = space_op.GetNDSpace().Get();
  MPI_Comm comm = space_op.GetComm();
  const int nt = fes.GetTrueVSize();
  double counts[2] = {0.0, 0.0};
  for (int i = 0; i < nt; i++)
  {
    counts[0] += iface.env_int[i];
    counts[1] += iface.gamma[i];
  }
  Mpi::GlobalSum(2, counts, comm);
  std::vector<double> fp = {counts[0], counts[1]};
  // A_E(ω) r = (K + iω C - ω² M + A2(ω)) r, partially assembled on the environment.
  std::unique_ptr<ComplexOperator> K, C, M, A2;
  {
    AssemblyDomains::Restriction restriction(space_op.GetAssemblyDomains(), env.attrs);
    K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO);
  }
  ComplexVector r(nt), y(nt);
  r.UseDevice(true);
  y.UseDevice(true);
  for (int k = 0; k < 3; k++)
  {
    r.Real() = FingerprintField(k);
    r.Imag() = 0.0;
    y = 0.0;
    K->AddMult(r, y, 1.0);
    if (C)
    {
      C->AddMult(r, y, std::complex<double>(0.0, omega));
    }
    if (M)
    {
      M->AddMult(r, y, -omega * omega);
    }
    if (A2)
    {
      A2->AddMult(r, y, 1.0);
    }
    double v[2] = {mfem::InnerProduct(r.Real(), y.Real()),
                   mfem::InnerProduct(r.Real(), y.Imag())};
    Mpi::GlobalSum(2, v, comm);
    fp.push_back(v[0]);
    fp.push_back(v[1]);
  }
  return fp;
}

std::vector<double> DrivenSubstructure::SourceFingerprint(const ComplexVector &b) const
{
  const auto &fes = space_op.GetNDSpace().Get();
  const double *br = b.Real().HostRead(), *bi = b.Imag().HostRead();
  constexpr int nfp = DrivenSubstructureModel::kSourceFp;
  double v[nfp] = {};
  for (int i = 0; i < fes.GetTrueVSize(); i++)
  {
    v[0] = std::max(v[0], (iface.env_int[i] && (br[i] != 0.0 || bi[i] != 0.0)) ? 1.0 : 0.0);
  }
  for (int k = 0; k < 3; k++)
  {
    const double *rr = FingerprintField(k).HostRead();
    for (int i = 0; i < fes.GetTrueVSize(); i++)
    {
      if (iface.env_int[i])
      {
        v[1 + 2 * k] += rr[i] * br[i];
        v[2 + 2 * k] += rr[i] * bi[i];
      }
    }
  }
  Mpi::GlobalMax(1, v, space_op.GetComm());
  Mpi::GlobalSum(nfp - 1, v + 1, space_op.GetComm());
  return {v, v + nfp};
}

bool DrivenSubstructure::OnInterface(const ComplexVector &b) const
{
  const double *br = b.Real().HostRead(), *bi = b.Imag().HostRead();
  int on = 0;
  for (std::size_t i = 0; i < iface.gamma.size() && !on; i++)
  {
    on = iface.gamma[i] && (br[i] != 0.0 || bi[i] != 0.0);
  }
  Mpi::GlobalMax(1, &on, space_op.GetComm());
  return on != 0;
}

namespace
{

// The records of a driven model in its file (see substructuring::ModelFile): the interface
// records; the header, a record per field or group of fields after DHED [kVersion], and the
// record of the environment (later ones supersede them, as the frequencies of an adaptive
// sweep are added; a reader skips the records it does not know); and the condensed
// environment at each frequency, in the order of the frequencies. The header records:
// DFPC [the fingerprint centers, 9 doubles], DENF [env_fp], DEXC [excitations], DSFP
// [exc_fp], DPRT [ports], DPRA [port_attributes], DFRQ [omega], DWGT [weights, empty for an
// exact model] and DSMP [error, tol, passivity, converged (1 or 0)].
using ModelFile = substructuring::ModelFile;
using substructuring::InterfaceRecords;
constexpr std::uint32_t kHeaderTag = ModelFile::Tag("DHED"),
                        kCentersTag = ModelFile::Tag("DFPC"),
                        kEnvFpTag = ModelFile::Tag("DENF"),
                        kExcitationsTag = ModelFile::Tag("DEXC"),
                        kSourceFpTag = ModelFile::Tag("DSFP"),
                        kPortsTag = ModelFile::Tag("DPRT"),
                        kPortAttributesTag = ModelFile::Tag("DPRA"),
                        kFrequenciesTag = ModelFile::Tag("DFRQ"),
                        kWeightsTag = ModelFile::Tag("DWGT"),
                        kSamplingTag = ModelFile::Tag("DSMP"),
                        kEnvRecordTag = ModelFile::Tag("ENVR"),
                        kRecordTag = ModelFile::Tag("DREC");

// Append the entries of x as the record tag.
template <typename T>
void AppendArray(const std::string &path, std::uint32_t tag, const std::vector<T> &x)
{
  ModelFile::Append(path, MPI_COMM_SELF, tag, static_cast<long long>(sizeof(T) * x.size()),
                    [&](std::ostream &f)
                    {
                      f.write(reinterpret_cast<const char *>(x.data()),
                              static_cast<std::streamsize>(sizeof(T) * x.size()));
                    });
}

// The records of the frequencies of a model file.
std::vector<ModelFile::Record> FrequencyRecords(const ModelFile &file)
{
  std::vector<ModelFile::Record> records;
  for (const auto &r : file.records)
  {
    if (r.tag == kRecordTag)
    {
      records.push_back(r);
    }
  }
  return records;
}

}  // namespace

bool DrivenSubstructureModel::SameEnvironment(const std::vector<double> &a,
                                              const std::vector<double> &b)
{
  return FingerprintsMatch(a, b, 2, kEnvironmentTol, true);
}

bool DrivenSubstructureModel::SameSource(const double *a, const double *b, double tol)
{
  return a[0] == b[0] &&
         PairingsMatch({a + 1, kSourceFp - 1}, {b + 1, kSourceFp - 1}, tol, true);
}

std::vector<std::complex<double>> BarycentricInterpolant::Sample(std::size_t j) const
{
  const std::size_t m = z.size();
  std::vector<std::complex<double>> x(Q[0].size(), 0.0);
  for (std::size_t i = 0; i <= j; i++)
  {
    const auto r = R[j * m + i];
    for (std::size_t q = 0; q < x.size(); q++)
    {
      x[q] += r * Q[i][q];
    }
  }
  return x;
}

void BarycentricInterpolant::AddSample(double omega,
                                       const std::vector<std::complex<double>> &x)
{
  // Orthonormalize x against the previous samples (classical Gram-Schmidt, twice), for the
  // factor R of the samples. The snapshots {x_j, iω_j x_j} are [Q 0; 0 Q] [R; R iΩ], so
  // their MRI weights are the right singular vector of [R; R iΩ] of least singular value,
  // also when it is rank deficient: weights that annihilate it reproduce the samples'
  // rational function (the adaptive driven solver's indicator, which only places samples,
  // takes the least singular value above its tolerance instead).
  const std::size_t m = z.size(), n = x.size();
  MFEM_VERIFY(m == 0 || Q[0].size() == n, "Samples of different sizes!");
  std::vector<std::complex<double>> v(x), r(m + 1, 0.0), d(m);
  for (int pass = 0; pass < 2; pass++)
  {
    for (std::size_t i = 0; i < m; i++)
    {
      d[i] = 0.0;
      for (std::size_t q = 0; q < n; q++)
      {
        d[i] += std::conj(Q[i][q]) * v[q];
      }
    }
    Mpi::GlobalSum(static_cast<int>(m), d.data(), comm);
    for (std::size_t i = 0; i < m; i++)
    {
      for (std::size_t q = 0; q < n; q++)
      {
        v[q] -= d[i] * Q[i][q];
      }
      r[i] += d[i];
    }
  }
  double nv = 0.0;
  for (const auto &c : v)
  {
    nv += std::norm(c);
  }
  Mpi::GlobalSum(1, &nv, comm);
  r[m] = std::sqrt(nv);
  for (auto &c : v)
  {
    c = (nv > 0.0) ? c / r[m].real() : 0.0;
  }
  Q.push_back(std::move(v));
  std::vector<std::complex<double>> R1((m + 1) * (m + 1), 0.0);
  for (std::size_t j = 0; j < m; j++)
  {
    std::copy(R.begin() + j * m, R.begin() + (j + 1) * m, R1.begin() + j * (m + 1));
  }
  std::copy(r.begin(), r.end(), R1.begin() + m * (m + 1));
  R = std::move(R1);
  z.push_back(omega);

  const int k = static_cast<int>(m + 1);
  Eigen::MatrixXcd A(2 * k, k);
  for (int j = 0; j < k; j++)
  {
    for (int i = 0; i < k; i++)
    {
      A(i, j) = R[j * k + i];
      A(k + i, j) = R[j * k + i] * std::complex<double>(0.0, z[j]);
    }
  }
  Eigen::JacobiSVD<Eigen::MatrixXcd, Eigen::ComputeFullV> svd(A);
  const auto q = svd.matrixV().col(k - 1);
  w.assign(q.data(), q.data() + k);
}

std::vector<std::complex<double>> BarycentricInterpolant::Coefficients(
    const std::vector<double> &z, const std::vector<std::complex<double>> &w, double omega)
{
  std::vector<std::complex<double>> a(z.size(), 0.0);
  std::complex<double> sum = 0.0;
  for (std::size_t j = 0; j < z.size(); j++)
  {
    if (std::abs(omega - z[j]) <= 1.0e-14 * std::abs(z[j]))
    {
      std::fill(a.begin(), a.end(), 0.0);
      a[j] = 1.0;
      return a;
    }
    a[j] = w[j] / (omega - z[j]);
    sum += a[j];
  }
  for (auto &c : a)
  {
    c /= sum;
  }
  return a;
}

std::vector<std::complex<double>>
BarycentricInterpolant::BasisCoefficients(double omega) const
{
  // y = R a.
  const std::size_t m = z.size();
  const auto a = Coefficients(z, w, omega);
  std::vector<std::complex<double>> y(m, 0.0);
  for (std::size_t j = 0; j < m; j++)
  {
    for (std::size_t i = 0; i <= j; i++)
    {
      y[i] += R[j * m + i] * a[j];
    }
  }
  return y;
}

std::vector<std::complex<double>> BarycentricInterpolant::Evaluate(double omega) const
{
  const auto y = BasisCoefficients(omega);
  std::vector<std::complex<double>> x(Q.empty() ? 0 : Q[0].size(), 0.0);
  for (std::size_t i = 0; i < Q.size(); i++)
  {
    for (std::size_t q = 0; q < x.size(); q++)
    {
      x[q] += Q[i][q] * y[i];
    }
  }
  return x;
}

namespace
{

// Im S_E on a grid, from the share of each process of the lower triangle of S_E, and
// ‖S_E‖_F. Collective.
std::pair<std::vector<double>, double> ImaginaryPart(const std::complex<double> *share,
                                                     int n, const BlockCyclicGrid &grid)
{
  MPI_Comm comm = grid.GetComm();
  const auto range = TriangleShare(n, comm);
  const long long e0 = range[0];
  double f = 0.0;
  auto [a, b] = TriangleEntry(n, e0);
  for (long long e = e0; e < range[1]; e++)
  {
    f += ((a == b) ? 1.0 : 2.0) * std::norm(share[e - e0]);
    if (++a == n)
    {
      a = ++b;
    }
  }
  Mpi::GlobalSum(1, &f, comm);
  std::vector<double> A(grid.LocalSize(), 0.0);
  RedistributeTriangle<double>(
      comm, n,
      [share, e0](long long e, long long m, double *x)
      {
        for (long long k = 0; k < m; k++)
        {
          x[k] = share[e - e0 + k].imag();
        }
      },
      {}, [&grid](int g, int h) { return grid.Owner(g, h); },
      [&grid](int g, int h) { return grid.LocalIndex(g, h); }, A);
  return {std::move(A), std::sqrt(f)};
}

}  // namespace

double DrivenSubstructureModel::Passivity(const std::complex<double> *share, int n,
                                          const BlockCyclicGrid &grid)
{
  auto [A, f] = ImaginaryPart(share, n, grid);
  const double lambda = LeastEigenvalue(grid, A);
  return (f > 0.0) ? lambda / f : 0.0;
}

bool DrivenSubstructureModel::Passive(const std::complex<double> *share, int n, double tol,
                                      const BlockCyclicGrid &grid)
{
  auto [A, f] = ImaginaryPart(share, n, grid);
  const int rank = Mpi::Rank(grid.GetComm());
  for (int g = 0; grid.InGrid() && g < n; g++)
  {
    if (grid.Owner(g, g) == rank)
    {
      A[grid.LocalIndex(g, g)] += tol * f;
    }
  }
  return PositiveDefinite(grid, A);
}

double BarycentricInterpolant::FindMaxError() const
{
  return MinimizeBarycentricDenominator(z, w.data());
}

std::array<std::size_t, 4> DrivenSubstructureModel::PartSizes() const
{
  const std::size_t n = nG, ne = excitations.size(), np = ports.size();
  return {n * (n + 1) / 2, n * ne, n * np, np * ne};
}

std::vector<std::complex<double>> DrivenSubstructureModel::Flatten(Record &&r)
{
  std::vector<std::complex<double>> x = std::move(r.S);
  x.insert(x.end(), r.g.begin(), r.g.end());
  x.insert(x.end(), r.h.begin(), r.h.end());
  x.insert(x.end(), r.c.begin(), r.c.end());
  r = Record();
  return x;
}

std::size_t DrivenSubstructureModel::RecordBytes() const
{
  const auto sz = PartSizes();
  return sizeof(std::complex<double>) * (sz[0] + sz[1] + sz[2] + sz[3]);
}

void DrivenSubstructureModel::WriteHeader(const std::string &path, bool create) const
{
  MFEM_VERIFY(weights.empty() || weights.size() == omega.size(),
              "Invalid substructuring model header!");
  if (create)
  {
    ModelFile::Create(path, MPI_COMM_SELF);
    InterfaceRecords{nG, sig_w, signatures, length_scale}.Append(path, MPI_COMM_SELF);
  }
  std::vector<double> centers;
  for (const auto &c : fp_centers)
  {
    centers.insert(centers.end(), c.begin(), c.end());
  }
  AppendArray(path, kHeaderTag, std::vector<int>{kVersion});
  AppendArray(path, kCentersTag, centers);
  AppendArray(path, kEnvFpTag, env_fp);
  AppendArray(path, kExcitationsTag, excitations);
  AppendArray(path, kSourceFpTag, exc_fp);
  AppendArray(path, kPortsTag, ports);
  AppendArray(path, kPortAttributesTag, port_attributes);
  AppendArray(path, kFrequenciesTag, omega);
  AppendArray(path, kWeightsTag, weights);
  AppendArray(path, kSamplingTag,
              std::vector<double>{error, tol, passivity, converged ? 1.0 : 0.0});
  const std::vector<char> record = environment.Serialize();
  ModelFile::Append(
      path, MPI_COMM_SELF, kEnvRecordTag, static_cast<long long>(record.size()),
      [&](std::ostream &f)
      { f.write(record.data(), static_cast<std::streamsize>(record.size())); });
}

void DrivenSubstructureModel::AppendRecord(const std::string &path, const Record &r,
                                           MPI_Comm comm) const
{
  // S_E by the columns of its lower triangle, from the shares of the processes, then g, h
  // and c.
  using T = std::complex<double>;
  const auto size = TriangleShareSizes(nG, comm);
  ModelFile::Append(
      path, comm, kRecordTag, static_cast<long long>(RecordBytes()),
      [&](std::ostream &f)
      {
        auto write = [&f](const T *x, std::size_t n)
        {
          f.write(reinterpret_cast<const char *>(x),
                  static_cast<std::streamsize>(sizeof(T) * n));
        };
        ToRoot<T>(comm, r.S, size, kTriangleChunk, write);
        write(r.g.data(), r.g.size());
        write(r.h.data(), r.h.size());
        write(r.c.data(), r.c.size());
      },
      [&]() { ToRoot<T>(comm, r.S, size, kTriangleChunk, {}); });
}

void DrivenSubstructureModel::ReadHeader(const std::string &path, MPI_Comm comm)
{
  const auto file = ModelFile::Read(path, comm);
  MFEM_VERIFY(file.version > 0, "Cannot read the driven substructuring model \""
                                    << path
                                    << "\" (run offline first, without \"ModelFile\")!");
  auto incomplete = [&](bool ok)
  {
    MFEM_VERIFY(ok, "The saved substructuring model \""
                        << path
                        << "\" is incomplete: rerun offline (without \"ModelFile\")!");
  };
  // The entries of the last record with the tag (a whole number of them).
  auto read = [&](std::uint32_t tag, auto &v)
  {
    using T = typename std::decay_t<decltype(v)>::value_type;
    const auto *r = file.Find(tag);
    incomplete(r && r->size % static_cast<long long>(sizeof(T)) == 0);
    const std::vector<char> p = ModelFile::ReadPayload(path, *r, comm);
    incomplete(p.size() == static_cast<std::size_t>(r->size));
    v.resize(p.size() / sizeof(T));
    std::memcpy(v.data(), p.data(), p.size());
  };
  std::vector<int> head;
  read(kHeaderTag, head);
  MFEM_VERIFY(head.size() == 1 && head[0] == kVersion,
              "The saved substructuring model \""
                  << path
                  << "\" was written by another version of Palace: rerun offline (without "
                     "\"ModelFile\")!");
  std::vector<double> centers, sampling;
  read(kCentersTag, centers);
  incomplete(centers.size() == 9);
  for (int k = 0; k < 3; k++)
  {
    std::copy_n(centers.begin() + 3 * k, 3, fp_centers[k].begin());
  }
  read(kEnvFpTag, env_fp);
  read(kExcitationsTag, excitations);
  read(kSourceFpTag, exc_fp);
  incomplete(exc_fp.size() == kSourceFp * excitations.size());
  read(kPortsTag, ports);
  read(kPortAttributesTag, port_attributes);
  read(kFrequenciesTag, omega);
  read(kWeightsTag, weights);
  incomplete(!omega.empty() && (weights.empty() || weights.size() == omega.size()));
  read(kSamplingTag, sampling);
  incomplete(sampling.size() == 4);
  error = sampling[0];
  tol = sampling[1];
  passivity = sampling[2];
  converged = sampling[3] != 0.0;
  InterfaceRecords iface;
  incomplete(iface.Read(path, file, comm));
  nG = iface.n;
  sig_w = iface.width;
  signatures = std::move(iface.signatures);
  length_scale = iface.length;
  const auto *e = file.Find(kEnvRecordTag);
  incomplete(e);
  const auto record =
      EnvironmentRecord::Deserialize(ModelFile::ReadPayload(path, *e, comm));
  incomplete(record.has_value());
  environment = *record;
}

int DrivenSubstructureModel::NumRecords(const std::string &path) const
{
  return static_cast<int>(FrequencyRecords(ModelFile::Read(path, MPI_COMM_SELF)).size());
}

void DrivenSubstructureModel::TruncateRecords(const std::string &path, int n) const
{
  // Before the tag and size of the next record, if any, else at the end of the last
  // complete record (dropping an incomplete one).
  const auto file = ModelFile::Read(path, MPI_COMM_SELF);
  const auto records = FrequencyRecords(file);
  const auto end = (static_cast<std::size_t>(n) < records.size())
                       ? records[n].pos - static_cast<long long>(sizeof(std::uint32_t) +
                                                                 sizeof(long long))
                       : file.records.back().pos + file.records.back().size;
  std::filesystem::resize_file(path, end);
}

bool DrivenSubstructureModel::SameHeader(const DrivenSubstructureModel &saved,
                                         double source_tol) const
{
  auto same = [](const std::vector<std::array<double, 6>> &a,
                 const std::vector<std::array<double, 6>> &b)
  {
    return std::ranges::equal(a, b, [](const auto &x, const auto &y)
                              { return MomentsMatch(x, y, kEnvironmentTol); });
  };
  const EnvironmentRecord &env = environment, &saved_env = saved.environment;
  if (saved.nG != nG || saved.sig_w != sig_w || saved.signatures != signatures ||
      saved.excitations != excitations || saved.exc_fp.size() != exc_fp.size() ||
      saved.ports != ports || saved.port_attributes != port_attributes ||
      saved.omega != omega || saved.weights != weights ||
      saved.length_scale != length_scale || !SameEnvironment(saved.env_fp, env_fp) ||
      saved_env.attributes != env.attributes ||
      saved_env.bdr_attributes != env.bdr_attributes ||
      saved_env.config_hash != env.config_hash ||
      saved_env.excitations_hash != env.excitations_hash ||
      !same(saved_env.moments, env.moments) ||
      !same(saved_env.region_moments, env.region_moments) ||
      !same(saved_env.bdr_moments, env.bdr_moments))
  {
    return false;
  }
  for (std::size_t e = 0; e < excitations.size(); e++)
  {
    if (!SameSource(exc_fp.data() + kSourceFp * e, saved.exc_fp.data() + kSourceFp * e,
                    source_tol))
    {
      return false;
    }
  }
  return true;
}

std::vector<long long> DrivenSubstructureModel::RecordPositions(const std::string &path,
                                                                MPI_Comm comm) const
{
  std::vector<long long> pos;
  for (const auto &r : FrequencyRecords(ModelFile::Read(path, comm)))
  {
    MFEM_VERIFY(r.size == static_cast<long long>(RecordBytes()),
                "Truncated substructuring model \"" << path << "\"!");
    pos.push_back(r.pos);
  }
  return pos;
}

std::vector<int> EnvironmentPorts(const IoData &iodata, const DrivenSubstructure &ds,
                                  std::vector<int> &port_attributes)
{
  std::vector<int> env_ports;
  port_attributes.clear();
  auto add = [&port_attributes](std::vector<int> attrs)
  {
    std::ranges::sort(attrs);
    port_attributes.push_back(static_cast<int>(attrs.size()));
    port_attributes.insert(port_attributes.end(), attrs.begin(), attrs.end());
  };
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
      n[0] += mark[i] && ds.GetInterface().env_int[i];
      n[1] += mark[i] && ds.GetInterface().gamma[i];
      n[2] += mark[i] && ds.GetInterface().region_int[i];
    }
    Mpi::GlobalSum(3, n, ds.GetComm());
    MFEM_VERIFY(n[1] == 0.0 && (n[0] == 0.0 || n[2] == 0.0),
                "Lumped port " << idx
                               << " touches the substructuring interface: a saved driven "
                                  "model needs the lumped ports away from it!");
    if (n[0] > 0.0)
    {
      env_ports.push_back(idx);
      add(attrs);
    }
  }
  for (int idx : ds.EnvironmentWavePorts())
  {
    env_ports.push_back(-idx);
    add(iodata.boundaries.waveport.at(idx).attributes);
  }
  return env_ports;
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

double SourceTolerance(const IoData &iodata, const PortExcitations &port_excitations,
                       const std::vector<int> &env_ports, int excitation_idx)
{
  double tol = 1.0e-10;
  for (int p : port_excitations.excitations.at(excitation_idx).wave_port)
  {
    if (std::ranges::find(env_ports, -p) != env_ports.end())
    {
      tol = std::max(tol, 100.0 * iodata.boundaries.waveport.at(p).eig_tol);
    }
  }
  return tol;
}

}  // namespace palace
