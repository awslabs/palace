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
#include "linalg/directschur.hpp"
#include "linalg/lowertriangle.hpp"
#include "linalg/mumpsschur.hpp"
#include "linalg/orthog.hpp"
#include "linalg/rap.hpp"
#include "models/portexcitations.hpp"
#include "models/romoperator.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

extern "C"
{
  void zgetrf_(const int *, const int *, std::complex<double> *, const int *, int *, int *);
  void zgetrs_(const char *, const int *, const int *, const std::complex<double> *,
               const int *, const int *, std::complex<double> *, const int *, int *);
  void dpotrf_(const char *, const int *, double *, const int *, int *);
  void dsyevr_(const char *, const char *, const char *, const int *, double *, const int *,
               const double *, const double *, const int *, const int *, const double *,
               int *, double *, double *, const int *, int *, double *, const int *, int *,
               const int *, int *);
}

namespace palace
{

namespace
{

#if defined(MFEM_USE_MUMPS)
// BLR tolerance of the factorizations, at round-off: MUMPS's BLR factorization is faster
// than its full-rank one here, at no loss of accuracy.
constexpr double kBlrTol = 1.0e-14;
#endif

using Factorization = SchurFactorization<std::complex<double>>;

// The sparse direct solver of the coupled system: that given, or the preferences.
LinearSolver CoupledSolver(LinearSolver type)
{
  if (type == LinearSolver::STRUMPACK || type == LinearSolver::MUMPS ||
      type == LinearSolver::SUPERLU)
  {
    return type;
  }
#if defined(MFEM_USE_STRUMPACK)
  return LinearSolver::STRUMPACK;
#elif defined(MFEM_USE_MUMPS)
  return LinearSolver::MUMPS;
#else
  return LinearSolver::SUPERLU;
#endif
}

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

LinearSolver SchurSolverType(LinearSolver type)
{
  if (type == LinearSolver::STRUMPACK || type == LinearSolver::MUMPS ||
      type == LinearSolver::SUPERLU)
  {
    return type;
  }
#if defined(MFEM_USE_MUMPS)
  return LinearSolver::MUMPS;
#elif defined(MFEM_USE_STRUMPACK)
  return LinearSolver::STRUMPACK;
#else
  return LinearSolver::SUPERLU;
#endif
}

std::unique_ptr<SchurFactorization<std::complex<double>>>
MakeSchurFactorization(LinearSolver type, SymbolicFactorization ordering, MPI_Comm comm,
                       HYPRE_BigInt n_glob,
                       SchurFactorization<std::complex<double>>::Coo &&coo,
                       const std::vector<HYPRE_BigInt> &schur_vars,
                       const std::vector<int> &rows, const BlockCyclicGrid *grid)
{
  const int n_loc = static_cast<int>(rows.size());
#if defined(MFEM_USE_MUMPS)
  if (type == LinearSolver::MUMPS)
  {
    return grid ? std::make_unique<MumpsSchurSolverT<std::complex<double>>>(
                      comm, n_glob, n_loc, std::move(coo), schur_vars, kBlrTol, true, false,
                      rows, ordering, grid->BlockSize(), grid->NumProcessRows(),
                      grid->NumProcessCols())
                : std::make_unique<MumpsSchurSolverT<std::complex<double>>>(
                      comm, n_glob, n_loc, std::move(coo), schur_vars, kBlrTol, true, false,
                      rows, ordering);
  }
#else
  MFEM_CONTRACT_VAR(grid);
#endif
  return std::make_unique<DirectSchurSolver>(comm, n_glob, n_loc, std::move(coo),
                                             schur_vars, rows, type, ordering);
}

DrivenSubstructure::DrivenSubstructure(SpaceOperator &space_op,
                                       const std::vector<int> &region_attributes,
                                       const std::vector<int> &environment_attributes,
                                       bool online, const DrivenSubstructureModel *saved,
                                       bool region_mesh, LinearSolver solver,
                                       SymbolicFactorization ordering)
  : space_op(space_op), solver(solver), ordering(ordering), online(online)
{
  env.attrs = environment_attributes;
  region.attrs = region_attributes;
  const auto &fes = space_op.GetNDSpace().Get();
  const auto &mesh = *fes.GetParMesh();
  CheckSubstructuringAttributes(mesh, region.attrs, env.attrs);
  fp_centers =
      saved ? saved->fp_centers
            : palace::FingerprintCenters(space_op.GetNDSpace().GetParMesh(), env.attrs);

  // Interface DOFs: on region and environment elements; environment (region) interior: on
  // environment (region) elements only; all without the Dirichlet DOFs.
  mfem::Array<int> rm, em;
  {
    mfem::Array<int> ra(region.attrs.data(), static_cast<int>(region.attrs.size())),
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
  auto classes = ClassifyInterfaceTrueDofs(rm, em, dbc, cut);
  is_gamma = std::move(classes.gamma);
  is_env_int = std::move(classes.env_int);
  is_region_int = std::move(classes.region_int);
  for (int i = 0; i < nt; i++)
  {
    if (!is_gamma[i] && !is_env_int[i])
    {
      env.pinned.Append(i);
    }
    if (!is_gamma[i] && !is_region_int[i])
    {
      region.pinned.Append(i);
    }
  }

  // Interface order: by owner rank, then local index.
  std::vector<HYPRE_BigInt> mine;
  const HYPRE_BigInt tstart = fes.GetMyTDofOffset();
  for (int i = 0; i < nt; i++)
  {
    if (is_gamma[i])
    {
      mine.push_back(tstart + i);
    }
  }
  const int nranks = Mpi::Size(comm), nloc = static_cast<int>(mine.size());
  gamma_cnt.assign(nranks, 0);
  gamma_disp.assign(nranks, 0);
  MPI_Allgather(&nloc, 1, MPI_INT, gamma_cnt.data(), 1, MPI_INT, comm);
  for (int r = 1; r < nranks; r++)
  {
    gamma_disp[r] = gamma_disp[r - 1] + gamma_cnt[r - 1];
  }
  gamma_tdofs.resize(gamma_disp[nranks - 1] + gamma_cnt[nranks - 1]);
  MPI_Allgatherv(mine.data(), nloc, HYPRE_MPI_BIG_INT, gamma_tdofs.data(), gamma_cnt.data(),
                 gamma_disp.data(), HYPRE_MPI_BIG_INT, comm);

  // Wave ports: each inside the region or the environment, away from Γ, so that its modal
  // terms (if active) border one side's system and its overlap is a functional of that
  // side's field.
  for (const auto &[idx, data] : space_op.GetWavePortOp())
  {
    const auto &list = data.GetAttrList();
    auto mark = BoundaryTrueDofs(std::vector<int>(list.begin(), list.end()));
    double n[3] = {0.0, 0.0, 0.0};  // environment interior, interface, region interior
    for (int i = 0; i < nt; i++)
    {
      n[0] += mark[i] && is_env_int[i];
      n[1] += mark[i] && is_gamma[i];
      n[2] += mark[i] && is_region_int[i];
    }
    Mpi::GlobalSum(3, n, comm);
    MFEM_VERIFY(
        n[1] == 0.0 && (n[0] == 0.0 || n[2] == 0.0),
        "Wave port " << idx
                     << " touches the substructuring interface: wave ports must lie "
                        "inside the region or the environment!");
    ((n[0] > 0.0) ? env : region).wave_ports.push_back(idx);
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
  // The terms of all active wave ports, each assigned to the port that holds its support.
  std::vector<std::pair<int, WavePortOperator::ModalCorrectionTerm>> terms;
  if (wave_port_dofs.empty())
  {
    return terms;
  }
  for (auto &term : space_op.GetWavePortOp().GetModalCorrectionTerms(
           omega, space_op.GetNDSpace(), space_op.GetNDDbcTDofLists().back()))
  {
    const double *sr = term.s->Real().HostRead(), *si = term.s->Imag().HostRead();
    std::vector<double> count;  // per port: support inside, outside
    for (const auto &[idx, mark] : wave_port_dofs)
    {
      double c[2] = {0.0, 0.0};
      for (std::size_t i = 0; i < mark.size(); i++)
      {
        c[mark[i] ? 0 : 1] += (sr[i] != 0.0 || si[i] != 0.0);
      }
      count.insert(count.end(), c, c + 2);
    }
    Mpi::GlobalSum(static_cast<int>(count.size()), count.data(), space_op.GetComm());
    int idx = -1;
    for (std::size_t p = 0; p < wave_port_dofs.size(); p++)
    {
      if (count[2 * p] > 0.0 && count[2 * p + 1] == 0.0)
      {
        idx = std::next(wave_port_dofs.begin(), static_cast<std::ptrdiff_t>(p))->first;
      }
    }
    MFEM_VERIFY(idx >= 0, "A wave-port modal term outside of its port!");
    if (std::ranges::find(side.wave_ports, idx) != side.wave_ports.end())
    {
      terms.emplace_back(idx, std::move(term));
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

std::vector<HYPRE_BigInt> OwnerValues(MPI_Comm comm, HYPRE_BigInt tstart,
                                      const std::vector<HYPRE_BigInt> &sys,
                                      const std::vector<int> &cols)
{
  const int nranks = Mpi::Size(comm);
  const HYPRE_BigInt tend = tstart + static_cast<HYPRE_BigInt>(sys.size());
  std::vector<HYPRE_BigInt> starts(nranks);
  MPI_Allgather(&tstart, 1, HYPRE_MPI_BIG_INT, starts.data(), 1, HYPRE_MPI_BIG_INT, comm);
  auto owner = [&](HYPRE_BigInt j)
  {
    return static_cast<int>(std::upper_bound(starts.begin(), starts.end(), j) -
                            starts.begin()) -
           1;
  };
  std::vector<std::vector<HYPRE_BigInt>> req(nranks);
  for (int c : cols)
  {
    const HYPRE_BigInt j = c - 1;
    if (j < tstart || j >= tend)
    {
      req[owner(j)].push_back(j);
    }
  }
  std::vector<int> scnt(nranks), sdisp(nranks, 0), rcnt(nranks), rdisp(nranks, 0);
  std::vector<HYPRE_BigInt> sbuf;
  for (int r = 0; r < nranks; r++)
  {
    std::sort(req[r].begin(), req[r].end());
    req[r].erase(std::unique(req[r].begin(), req[r].end()), req[r].end());
    scnt[r] = static_cast<int>(req[r].size());
    sdisp[r] = static_cast<int>(sbuf.size());
    sbuf.insert(sbuf.end(), req[r].begin(), req[r].end());
  }
  MPI_Alltoall(scnt.data(), 1, MPI_INT, rcnt.data(), 1, MPI_INT, comm);
  for (int r = 1; r < nranks; r++)
  {
    rdisp[r] = rdisp[r - 1] + rcnt[r - 1];
  }
  std::vector<HYPRE_BigInt> rbuf(rdisp[nranks - 1] + rcnt[nranks - 1]);
  MPI_Alltoallv(sbuf.data(), scnt.data(), sdisp.data(), HYPRE_MPI_BIG_INT, rbuf.data(),
                rcnt.data(), rdisp.data(), HYPRE_MPI_BIG_INT, comm);
  for (auto &j : rbuf)
  {
    j = sys[j - tstart];
  }
  std::vector<HYPRE_BigInt> ans(sbuf.size());
  MPI_Alltoallv(rbuf.data(), rcnt.data(), rdisp.data(), HYPRE_MPI_BIG_INT, ans.data(),
                scnt.data(), sdisp.data(), HYPRE_MPI_BIG_INT, comm);
  std::vector<HYPRE_BigInt> v(cols.size());
  for (std::size_t q = 0; q < cols.size(); q++)
  {
    const HYPRE_BigInt j = cols[q] - 1;
    if (j >= tstart && j < tend)
    {
      v[q] = sys[j - tstart];
    }
    else
    {
      const int r = owner(j);
      v[q] = ans[sdisp[r] +
                 (std::lower_bound(req[r].begin(), req[r].end(), j) - req[r].begin())];
    }
  }
  return v;
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
  // The factored system is the side's interior and Γ, numbered by rank: the pinned DOFs
  // would be identity rows, and MUMPS's per-rank workspace grows with the system size.
  const auto &fes = space_op.GetNDSpace().Get();
  MPI_Comm comm = fes.GetComm();
  const int nt = fes.GetTrueVSize();
  const HYPRE_BigInt tstart = fes.GetMyTDofOffset();
  std::vector<HYPRE_BigInt> sys(nt, 0);
  for (int i : side.pinned)
  {
    sys[i] = -1;
  }
  side.rows.clear();
  for (int i = 0; i < nt; i++)
  {
    if (sys[i] == 0)
    {
      side.rows.push_back(i);
    }
  }
  HYPRE_BigInt nloc = static_cast<HYPRE_BigInt>(side.rows.size()), off = 0;
  MPI_Exscan(&nloc, &off, 1, HYPRE_MPI_BIG_INT, MPI_SUM, comm);
  off = Mpi::Root(comm) ? 0 : off;
  MPI_Allreduce(&nloc, &side.n_sys, 1, HYPRE_MPI_BIG_INT, MPI_SUM, comm);
  for (std::size_t q = 0; q < side.rows.size(); q++)
  {
    sys[side.rows[q]] = off + static_cast<HYPRE_BigInt>(q);
  }
  std::vector<HYPRE_BigInt> mine;
  for (int i = 0; i < nt; i++)
  {
    if (is_gamma[i])
    {
      mine.push_back(sys[i]);
    }
  }
  side.gamma_sys.resize(gamma_tdofs.size());
  MPI_Allgatherv(mine.data(), static_cast<int>(mine.size()), HYPRE_MPI_BIG_INT,
                 side.gamma_sys.data(), gamma_cnt.data(), gamma_disp.data(),
                 HYPRE_MPI_BIG_INT, comm);

  // The frequency-independent parts, assembled once and one at a time: the pattern of the
  // stiffness matrix on the system, which contains those of the other parts, coupling DOFs
  // of elements of the side.
  side.parts.assign(6, {});
  for (int op = 0; op < 3; op++)
  {
    auto parts = SideParts(side, op);
    for (int p = 0; p < 2; p++)
    {
      const auto &X = parts[p];
      if (!X)
      {
        MFEM_VERIFY(op > 0 || p > 0,
                    "Missing stiffness matrix of a substructure operator!");
        continue;
      }
      auto &part = side.parts[2 * op + p];
      if (op == 0 && p == 0)
      {
        // The pattern without the (zero) entries of pinned rows and columns.
        std::vector<int> irn, jcn, row_ptr;
        std::vector<std::vector<double>> vals;
        LowerTrianglePattern({X.get()}, irn, jcn, vals, row_ptr);
        const auto col_sys = OwnerValues(comm, tstart, sys, jcn);
        side.row_ptr.assign(nt + 1, 0);
        side.jcn.clear();
        side.irn_sys.clear();
        side.jcn_sys.clear();
        part.clear();
        for (int i = 0; i < nt; i++)
        {
          for (int q = row_ptr[i]; sys[i] >= 0 && q < row_ptr[i + 1]; q++)
          {
            if (col_sys[q] >= 0)
            {
              side.jcn.push_back(jcn[q]);
              side.irn_sys.push_back(static_cast<int>(sys[i] + 1));
              side.jcn_sys.push_back(static_cast<int>(col_sys[q] + 1));
              part.push_back(vals[0][q]);
            }
          }
          side.row_ptr[i + 1] = static_cast<int>(side.jcn.size());
        }
      }
      else
      {
        part.assign(side.jcn.size(), 0.0);
        AddToPattern(*X, 1.0, side.row_ptr, side.jcn, part);
      }
    }
  }
  // The modal terms of the side's wave ports border the system, with an extra unknown per
  // term on the last rank, against the DOFs of the term's port. Its diagonal is -σ and its
  // entries √(gσ) s, with σ the largest diagonal entry of the side's operator at this
  // frequency, so that the bordering does not change the scale of the factored matrix.
  const auto terms = ModalTerms(side, omega);
  {
    const auto coef = Coefficients(omega);
    side.border = 0.0;
    for (int i = 0; i < nt; i++)
    {
      for (int q = side.row_ptr[i]; q < side.row_ptr[i + 1]; q++)
      {
        if (side.jcn[q] == tstart + i + 1)
        {
          std::complex<double> d = 0.0;
          for (int p = 0; p < 6; p++)
          {
            d += side.parts[p].empty() ? 0.0 : coef[p] * side.parts[p][q];
          }
          side.border = std::max(side.border, std::abs(d));
        }
      }
    }
    Mpi::GlobalMax(1, &side.border, comm);
    side.border = (side.border > 0.0) ? side.border : 1.0;
  }
  const int nterm = static_cast<int>(terms.size());
  const bool last = (Mpi::Rank(comm) == Mpi::Size(comm) - 1);
  side.n_dofs = side.n_sys;
  side.n_sys += nterm;
  side.term_ports.clear();
  side.term_dofs.assign(nterm, {});
  for (int t = 0; t < nterm; t++)
  {
    side.term_ports.push_back(terms[t].first);
    const auto &mark = wave_port_dofs.at(terms[t].first);
    for (int i = 0; i < nt; i++)
    {
      if (mark[i] && sys[i] >= 0)
      {
        side.term_dofs[t].push_back(i);
        side.irn_sys.push_back(static_cast<int>(side.n_dofs + t + 1));
        side.jcn_sys.push_back(static_cast<int>(sys[i] + 1));
      }
    }
    if (last)
    {
      side.irn_sys.push_back(static_cast<int>(side.n_dofs + t + 1));
      side.jcn_sys.push_back(static_cast<int>(side.n_dofs + t + 1));
      side.rows.push_back(-1);
    }
  }
  side.val.resize(side.irn_sys.size());
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
  std::fill(val.begin(), val.end(), 0.0);
  for (int p = 0; p < 6; p++)
  {
    for (std::size_t q = 0; q < side.parts[p].size(); q++)
    {
      val[q] += coef[p] * side.parts[p][q];
    }
  }
  for (int p = 0; p < 2; p++)
  {
    if (extra[p])
    {
      AddToPattern(*extra[p], coef[6 + p], side.row_ptr, side.jcn, val);
    }
  }
  MFEM_VERIFY(terms.size() == side.term_dofs.size(),
              "The modal terms of the wave ports changed between frequencies!");
  const bool last = (Mpi::Rank(space_op.GetComm()) == Mpi::Size(space_op.GetComm()) - 1);
  std::size_t q = side.jcn.size();
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
    side.schur = MakeSchurFactorization(
        SchurSolverType(solver), ordering, space_op.GetComm(), side.n_sys,
        {std::move(side.irn_sys), std::move(side.jcn_sys), std::move(side.val)},
        side.gamma_sys, side.rows);
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

std::vector<int> DrivenSubstructure::InterfaceIndex() const
{
  const int nt = static_cast<int>(is_gamma.size());
  std::vector<int> index(nt, -1);
  const int off = gamma_disp[Mpi::Rank(space_op.GetComm())];
  for (int i = 0, g = 0; i < nt; i++)
  {
    if (is_gamma[i])
    {
      index[i] = off + g++;
    }
  }
  return index;
}

void DrivenSubstructure::Condense(double omega, bool region_too)
{
  MFEM_VERIFY(!online, "An online substructure is condensed with a given S_E!");
  if (!env.schur)
  {
    // The patterns first, so that the assembled operators are released before the
    // factorizations.
    Setup(env, omega);
  }
  const bool set_up = region_too && !coupled;
  if (set_up)
  {
    Setup(region, omega);
  }
  // S_E on rank 0 moves out of MUMPS's Schur buffer without a copy (the buffer is allocated
  // again at the next factorization).
  Release(S);
  Factor(env, omega);
  S = env.schur->ReleaseSchur();
  if (region_too)
  {
    FactorCoupled(omega, S, false, set_up);
  }
  solvable = region_too;
}

void DrivenSubstructure::Condense(double omega, std::vector<std::complex<double>> &&S_env)
{
  MFEM_VERIFY(online, "An offline substructure condenses its own environment!");
  FactorCoupled(omega, S_env, true, false);
  solvable = true;
}

void DrivenSubstructure::FactorCoupled(double omega,
                                       std::vector<std::complex<double>> &S_env,
                                       bool release, bool set_up)
{
  MPI_Comm comm = space_op.GetComm();
  const bool root = Mpi::Root(comm);
  const int nG = InterfaceSize();
  MFEM_VERIFY(!root || S_env.size() == static_cast<std::size_t>(nG) * nG,
              "Wrong size of a given S_E!");
  // The region with S_E added on Γ, factored as one system: the region's entries, then on
  // rank 0 those of S_E's lower triangle, in the rows of the interface DOFs (whose system
  // indices increase with the interface order).
  auto add_schur = [&](auto &&add)
  {
    for (int q = 0; q < nG && root; q++)
    {
      for (int p = q; p < nG; p++)
      {
        add(region.gamma_sys[p], region.gamma_sys[q],
            S_env[p + static_cast<std::size_t>(nG) * q]);
      }
    }
  };
  if (!coupled)
  {
    if (!set_up)
    {
      Setup(region, omega);
    }
    Factorization::Coo coo{std::move(region.irn_sys), std::move(region.jcn_sys),
                           std::move(region.val)};
    coupled_region_entries = coo.val.size();
    add_schur(
        [&](HYPRE_BigInt a, HYPRE_BigInt b, std::complex<double> v)
        {
          coo.irn.push_back(static_cast<int>(a + 1));
          coo.jcn.push_back(static_cast<int>(b + 1));
          coo.val.push_back(v);
        });
    if (release)
    {
      Release(S_env);
    }
    coupled = MakeSchurFactorization(CoupledSolver(solver), ordering, comm, region.n_sys,
                                     std::move(coo), {}, region.rows);
  }
  else
  {
    auto &val = coupled->Values();
    auto extra = ExtraParts(region, omega);
    MFEM_VERIFY(region.extra == (extra[0] || extra[1]),
                "The frequency-dependent part of a substructure operator changed!");
    Fill(region, omega, extra, ModalTerms(region, omega), val);
    std::size_t k = coupled_region_entries;
    add_schur([&](HYPRE_BigInt, HYPRE_BigInt, std::complex<double> v) { val[k++] = v; });
    if (release)
    {
      Release(S_env);
    }
    coupled->Refactor();
  }
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

void DrivenSubstructure::Solve(const std::vector<const ComplexVector *> &rhs,
                               std::vector<ComplexVector> &u,
                               const std::vector<std::complex<double>> *g_env)
{
  MFEM_VERIFY(solvable, "Condense (with the region) must be called before Solve!");
  MFEM_VERIFY(online == (g_env != nullptr),
              "The environment's source condensation is given online only!");
  if (online)
  {
    SolveCoupled(rhs, u, *g_env);
    return;
  }

  // The environment's sources condensed onto Γ with the interface loads b_Γ,
  // b_Γ - A_ΓE A_EE^-1 b_E, the region and Γ solved with them as online, and the
  // environment interior from the interface solution, A_EE^-1 (b_E - A_EΓ u_Γ).
  const int n = static_cast<int>(rhs.size()), nt = space_op.GetNDSpace().GetTrueVSize();
  std::vector<ComplexVector> w(n);
  std::vector<const ComplexVector *> W(n);
  std::vector<ComplexVector *> Wo(n);
  for (int k = 0; k < n; k++)
  {
    w[k] = Masked(*rhs[k], is_env_int, &is_gamma);
    W[k] = &w[k];
    Wo[k] = &w[k];
  }
  env.schur->Reduce(W, g_last);
  SolveCoupled(rhs, u, g_last);
  env.schur->Expand(u_last, Wo);
  for (int k = 0; k < n; k++)
  {
    double *ur = u[k].Real().HostReadWrite(), *ui = u[k].Imag().HostReadWrite();
    const double *wr = w[k].Real().HostRead(), *wi = w[k].Imag().HostRead();
    for (int i = 0; i < nt; i++)
    {
      if (is_env_int[i])
      {
        ur[i] = wr[i];
        ui[i] = wi[i];
      }
    }
  }
}

void DrivenSubstructure::SolveCoupled(const std::vector<const ComplexVector *> &rhs,
                                      std::vector<ComplexVector> &u,
                                      const std::vector<std::complex<double>> &g_env)
{
  // The region's sources on its interior and the environment's condensation on Γ (from rank
  // 0, by the ranks of the interface DOFs), solved for the region and Γ at once.
  MPI_Comm comm = space_op.GetComm();
  const int n = static_cast<int>(rhs.size()), nt = static_cast<int>(is_gamma.size());
  const int nG = InterfaceSize(), rank = Mpi::Rank(comm), nranks = Mpi::Size(comm);
  const bool root = (rank == 0);
  MFEM_VERIFY(!root || g_env.size() == static_cast<std::size_t>(nG) * n,
              "Wrong size of the environment's source condensation!");
  std::vector<int> cnt(nranks), disp(nranks);
  for (int r = 0; r < nranks; r++)
  {
    cnt[r] = gamma_cnt[r] * n;
    disp[r] = gamma_disp[r] * n;
  }
  const int m = gamma_cnt[rank];
  std::vector<std::complex<double>> blocks(root ? static_cast<std::size_t>(nG) * n : 0),
      loc(static_cast<std::size_t>(m) * n);
  // Rank-major blocks of the interface rows, each column-major.
  auto block = [&](int r, int l, int k)
  {
    return static_cast<std::size_t>(disp[r]) + l +
           static_cast<std::size_t>(gamma_cnt[r]) * k;
  };
  for (int r = 0; r < nranks && root; r++)
  {
    for (int k = 0; k < n; k++)
    {
      for (int l = 0; l < gamma_cnt[r]; l++)
      {
        blocks[block(r, l, k)] =
            g_env[gamma_disp[r] + l + static_cast<std::size_t>(nG) * k];
      }
    }
  }
  MPI_Scatterv(blocks.data(), cnt.data(), disp.data(), MPI_C_DOUBLE_COMPLEX, loc.data(),
               m * n, MPI_C_DOUBLE_COMPLEX, 0, comm);
  std::vector<ComplexVector> y(n);
  std::vector<const ComplexVector *> Y(n);
  std::vector<ComplexVector *> Yo(n);
  for (int k = 0; k < n; k++)
  {
    y[k] = Masked(*rhs[k], is_region_int);
    double *yr = y[k].Real().HostReadWrite(), *yi = y[k].Imag().HostReadWrite();
    for (int i = 0, l = 0; i < nt; i++)
    {
      if (is_gamma[i])
      {
        const auto v = loc[l++ + static_cast<std::size_t>(m) * k];
        yr[i] = v.real();
        yi[i] = v.imag();
      }
    }
    Y[k] = &y[k];
    Yo[k] = &y[k];
  }
  coupled->SolveInternal(Y, Yo);

  // The solution on the region's interior and Γ, and its interface part on rank 0.
  u.resize(n);
  for (int k = 0; k < n; k++)
  {
    u[k].SetSize(nt);
    u[k].UseDevice(true);
    u[k] = y[k];
    const double *yr = y[k].Real().HostRead(), *yi = y[k].Imag().HostRead();
    for (int i = 0, l = 0; i < nt; i++)
    {
      if (is_gamma[i])
      {
        loc[l++ + static_cast<std::size_t>(m) * k] = {yr[i], yi[i]};
      }
    }
  }
  MPI_Gatherv(loc.data(), m * n, MPI_C_DOUBLE_COMPLEX, blocks.data(), cnt.data(),
              disp.data(), MPI_C_DOUBLE_COMPLEX, 0, comm);
  if (root)
  {
    g_last = g_env;
    u_last.resize(static_cast<std::size_t>(nG) * n);
    for (int r = 0; r < nranks; r++)
    {
      for (int k = 0; k < n; k++)
      {
        for (int l = 0; l < gamma_cnt[r]; l++)
        {
          u_last[gamma_disp[r] + l + static_cast<std::size_t>(nG) * k] =
              blocks[block(r, l, k)];
        }
      }
    }
  }
}

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
    x[k] = Masked(*b[k], is_env_int, &is_gamma);
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
    if (is_gamma[i])
    {
      t(i) = 2.0;
    }
  }
  mfem::ParGridFunction g(&fes);
  g.SetFromTrueDofs(t);
  const int amax = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  std::vector<int> found(amax + 1, 0);
  mfem::Array<int> dofs;
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
      const double v = g((d >= 0) ? d : -1 - d);
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
  for (int a : region.attrs)
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
      palace::WriteRegionMesh(path, fes, region.attrs, env.attrs, keep, length_scale);
  Mpi::Print(" Writing the mesh of the region for online runs to {}\n", path);
  return written;
}

const Vector &DrivenSubstructure::FingerprintField(int k) const
{
  if (fp_fields[k].Size() == 0)
  {
    fp_fields[k] = palace::FingerprintField(space_op.GetNDSpace().Get(), k, fp_centers[k]);
  }
  return fp_fields[k];
}

std::vector<double> DrivenSubstructure::EnvironmentFingerprint(double omega) const
{
  const auto &fes = space_op.GetNDSpace().Get();
  MPI_Comm comm = space_op.GetComm();
  const int nt = fes.GetTrueVSize();
  double counts[2] = {0.0, 0.0};
  for (int i = 0; i < nt; i++)
  {
    counts[0] += is_env_int[i];
    counts[1] += is_gamma[i];
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
    v[0] = std::max(v[0], (is_env_int[i] && (br[i] != 0.0 || bi[i] != 0.0)) ? 1.0 : 0.0);
  }
  for (int k = 0; k < 3; k++)
  {
    const double *rr = FingerprintField(k).HostRead();
    for (int i = 0; i < fes.GetTrueVSize(); i++)
    {
      if (is_env_int[i])
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
  for (std::size_t i = 0; i < is_gamma.size() && !on; i++)
  {
    on = is_gamma[i] && (br[i] != 0.0 || bi[i] != 0.0);
  }
  Mpi::GlobalMax(1, &on, space_op.GetComm());
  return on != 0;
}

namespace
{

// The records of a driven model in its file (see substructuring::ModelFile): the interface
// records, the header and the record of the environment (later ones supersede them, as the
// frequencies of an adaptive sweep are added) and the condensed environment at each
// frequency, in the order of the frequencies.
using ModelFile = substructuring::ModelFile;
using substructuring::InterfaceRecords;
constexpr std::uint32_t kHeaderTag = ModelFile::Tag("DHED"),
                        kEnvRecordTag = ModelFile::Tag("ENVR"),
                        kRecordTag = ModelFile::Tag("DREC");

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

double DrivenSubstructureModel::Passivity(const std::complex<double> *lower, int n)
{
  std::vector<double> A(static_cast<std::size_t>(n) * n, 0.0);
  double f = 0.0;
  std::size_t q = 0;
  for (int j = 0; j < n; j++)
  {
    for (int i = j; i < n; i++, q++)
    {
      A[static_cast<std::size_t>(j) * n + i] = lower[q].imag();
      f += ((i == j) ? 1.0 : 2.0) * std::norm(lower[q]);
    }
  }
  // The least eigenvalue only (LAPACK dsyevr). W needs n entries even for one eigenvalue
  // (dstebz stores the whole bracket first).
  const int one = 1;
  const double zero = 0.0;
  int found = 0, info = 0, lwork = -1, liwork = -1, iwq = 0;
  double wq = 0.0, z = 0.0;
  std::vector<double> w(n);
  std::vector<int> isuppz(2);
  dsyevr_("N", "I", "L", &n, A.data(), &n, &zero, &zero, &one, &one, &zero, &found,
          w.data(), &z, &one, isuppz.data(), &wq, &lwork, &iwq, &liwork, &info);
  lwork = static_cast<int>(wq);
  liwork = iwq;
  std::vector<double> work(lwork);
  std::vector<int> iwork(liwork);
  dsyevr_("N", "I", "L", &n, A.data(), &n, &zero, &zero, &one, &one, &zero, &found,
          w.data(), &z, &one, isuppz.data(), work.data(), &lwork, iwork.data(), &liwork,
          &info);
  MFEM_VERIFY(info == 0 && found == 1,
              "Eigenvalue of Im S_E failed: info = " << info << "!");
  return (f > 0.0) ? w[0] / std::sqrt(f) : 0.0;
}

bool DrivenSubstructureModel::Passive(const std::complex<double> *lower, int n, double tol)
{
  std::vector<double> A(static_cast<std::size_t>(n) * n, 0.0);
  double f = 0.0;
  std::size_t q = 0;
  for (int j = 0; j < n; j++)
  {
    for (int i = j; i < n; i++, q++)
    {
      A[static_cast<std::size_t>(j) * n + i] = lower[q].imag();
      f += ((i == j) ? 1.0 : 2.0) * std::norm(lower[q]);
    }
  }
  for (int j = 0; j < n; j++)
  {
    A[static_cast<std::size_t>(j) * n + j] += tol * std::sqrt(f);
  }
  int info = 0;
  dpotrf_("L", &n, A.data(), &n, &info);
  return info == 0;
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

std::vector<std::complex<double>> DrivenSubstructureModel::Flatten(const Record &r) const
{
  // S_E is symmetric: its lower triangle, by columns.
  std::vector<std::complex<double>> x;
  const auto sz = PartSizes();
  x.reserve(sz[0] + sz[1] + sz[2] + sz[3]);
  for (int j = 0; j < nG; j++)
  {
    for (int i = j; i < nG; i++)
    {
      x.push_back(r.S[static_cast<std::size_t>(j) * nG + i]);
    }
  }
  x.insert(x.end(), r.g.begin(), r.g.end());
  x.insert(x.end(), r.h.begin(), r.h.end());
  x.insert(x.end(), r.c.begin(), r.c.end());
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
  const int head[8] = {kVersion,
                       static_cast<int>(env_fp.size()),
                       static_cast<int>(excitations.size()),
                       static_cast<int>(ports.size()),
                       static_cast<int>(omega.size()),
                       weights.empty() ? 0 : 1,
                       static_cast<int>(region_bdr_attributes.size()),
                       converged ? 1 : 0};
  const double scalars[3] = {error, tol, passivity};
  const long long size =
      sizeof(head) + sizeof(scalars) + sizeof(fp_centers) +
      sizeof(double) * (env_fp.size() + exc_fp.size() + omega.size()) +
      sizeof(int) * (excitations.size() + ports.size() + region_bdr_attributes.size()) +
      sizeof(std::complex<double>) * weights.size();
  ModelFile::Append(path, MPI_COMM_SELF, kHeaderTag, size,
                    [&](std::ostream &f)
                    {
                      auto write = [&f](const auto *x, std::size_t n)
                      {
                        f.write(reinterpret_cast<const char *>(x),
                                static_cast<std::streamsize>(sizeof(*x) * n));
                      };
                      write(head, 8);
                      write(scalars, 3);
                      for (const auto &c : fp_centers)
                      {
                        write(c.data(), 3);
                      }
                      write(env_fp.data(), env_fp.size());
                      write(excitations.data(), excitations.size());
                      write(exc_fp.data(), exc_fp.size());
                      write(ports.data(), ports.size());
                      write(region_bdr_attributes.data(), region_bdr_attributes.size());
                      write(omega.data(), omega.size());
                      write(weights.data(), weights.size());
                    });
  const std::vector<char> record = environment.Serialize();
  ModelFile::Append(
      path, MPI_COMM_SELF, kEnvRecordTag, static_cast<long long>(record.size()),
      [&](std::ostream &f)
      { f.write(record.data(), static_cast<std::streamsize>(record.size())); });
}

void DrivenSubstructureModel::AppendRecord(const std::string &path, const Record &r) const
{
  const auto x = Flatten(r);
  ModelFile::Append(
      path, MPI_COMM_SELF, kRecordTag,
      static_cast<long long>(sizeof(std::complex<double>) * x.size()),
      [&](std::ostream &f)
      {
        f.write(reinterpret_cast<const char *>(x.data()),
                static_cast<std::streamsize>(sizeof(std::complex<double>) * x.size()));
      });
}

void DrivenSubstructureModel::ReadHeader(const std::string &path, MPI_Comm comm)
{
  const auto file = ModelFile::Read(path, comm);
  MFEM_VERIFY(file.version > 0,
              "Cannot read the driven substructuring model \""
                  << path << "\" (run in \"Offline\" mode with \"ModelFile\" first)!");
  auto incomplete = [&](bool ok)
  {
    MFEM_VERIFY(ok, "The saved substructuring model \""
                        << path << "\" is incomplete: rerun in \"Offline\" mode!");
  };
  const auto *h = file.Find(kHeaderTag);
  incomplete(h);
  const std::vector<char> p = ModelFile::ReadPayload(path, *h, comm);
  // The next n entries of the payload, if it holds them (n counted in floating point, so
  // that any size is checked before it is allocated).
  std::size_t at = 0;
  auto read = [&](auto &v, double n)
  {
    using T = typename std::decay_t<decltype(v)>::value_type;
    incomplete(n >= 0.0 && n * sizeof(T) <= static_cast<double>(p.size() - at));
    v.resize(static_cast<std::size_t>(n));
    std::memcpy(v.data(), p.data() + at, sizeof(T) * v.size());
    at += sizeof(T) * v.size();
  };
  std::vector<int> head;
  read(head, 1);
  MFEM_VERIFY(head[0] == kVersion, "The saved substructuring model \""
                                       << path
                                       << "\" was written by another version of Palace: "
                                          "rerun in \"Offline\" mode!");
  at = 0;  // the version, then the whole head
  read(head, 8);
  std::vector<double> scalars;
  read(scalars, 3);
  error = scalars[0];
  tol = scalars[1];
  passivity = scalars[2];
  converged = head[7] != 0;
  std::vector<double> centers;
  read(centers, 9);
  for (int k = 0; k < 3; k++)
  {
    std::copy_n(centers.begin() + 3 * k, 3, fp_centers[k].begin());
  }
  read(env_fp, head[1]);
  read(excitations, head[2]);
  read(exc_fp, static_cast<double>(kSourceFp) * head[2]);
  read(ports, head[3]);
  read(region_bdr_attributes, head[6]);
  read(omega, head[4]);
  read(weights, head[5] ? head[4] : 0);
  incomplete(head[4] > 0 && at == p.size());
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
      saved.ports != ports || saved.region_bdr_attributes != region_bdr_attributes ||
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

std::vector<std::complex<double>>
DrivenSubstructureModel::ReadFlatRecord(const std::string &path, int j) const
{
  const auto records = FrequencyRecords(ModelFile::Read(path, MPI_COMM_SELF));
  MFEM_VERIFY(static_cast<std::size_t>(j) < records.size() &&
                  records[j].size == static_cast<long long>(RecordBytes()),
              "Truncated substructuring model \"" << path << "\"!");
  const std::vector<char> p = ModelFile::ReadPayload(path, records[j], MPI_COMM_SELF);
  std::vector<std::complex<double>> x(RecordBytes() / sizeof(std::complex<double>));
  MFEM_VERIFY(p.size() == RecordBytes(),
              "Truncated substructuring model \"" << path << "\"!");
  std::memcpy(x.data(), p.data(), p.size());
  return x;
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
  auto read_row = [&](int a, std::complex<double> *x)
  {
    for (int b = 0; b < n; b++)
    {
      x[b] = at(a, b);
    }
  };
  if (map)
  {
    return map->DualMatrix<std::complex<double>>(0, n, read_row);
  }
  std::vector<std::complex<double>> S(static_cast<std::size_t>(n) * n);
  for (int a = 0; a < n; a++)
  {
    read_row(a, S.data() + static_cast<std::size_t>(a) * n);
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

void RegionPostprocessing(const std::vector<int> &region_attributes,
                          const mfem::ParMesh &mesh,
                          const std::vector<int> *region_bdr_attributes, bool env_ports,
                          const char *what, config::DomainData &domains,
                          config::BoundaryData &boundaries)
{
  domains.postpro.partial = true;
  std::vector<int> skipped;
  for (auto it = domains.postpro.energy.begin(); it != domains.postpro.energy.end();)
  {
    const bool region = std::ranges::all_of(
        it->second.attributes, [&](int a)
        { return std::ranges::find(region_attributes, a) != region_attributes.end(); });
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
  auto bdr_in_region = BoundaryAttributesInRegion(mesh, region_attributes);
  if (region_bdr_attributes)
  {
    std::erase_if(bdr_in_region,
                  [&](int a) { return !contains(*region_bdr_attributes, a); });
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
  // A probe is interpolated in the element it is found in, which for a probe on Γ can be
  // one of the environment.
  const auto probes_in_region =
      fem::ProbesInDomains(domains.postpro.probe, mesh, region_attributes);
  for (auto it = domains.postpro.probe.begin(); it != domains.postpro.probe.end();)
  {
    if (!contains(probes_in_region, it->first))
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
               what,
               none(skipped, "energies of the domain postprocessing outside the region"),
               none(skipped_flux, "surface fluxes outside the region"),
               none(skipped_probe, "probes outside the region"),
               !env_ports ? std::string()
                          : ", and no power (or wave port voltage) of the "
                            "environment's ports");
}

}  // namespace palace
