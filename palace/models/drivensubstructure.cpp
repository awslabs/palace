// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "drivensubstructure.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <utility>
#include <Eigen/Dense>
#include "fem/substructure.hpp"
#include "linalg/mumpsschur.hpp"
#include "linalg/rap.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"

extern "C"
{
  void zsytrf_(const char *, const int *, std::complex<double> *, const int *, int *,
               std::complex<double> *, const int *, int *);
  void zsytrs_(const char *, const int *, const int *, const std::complex<double> *,
               const int *, const int *, std::complex<double> *, const int *, int *);
  void dsyevr_(const char *, const char *, const char *, const int *, double *, const int *,
               const double *, const double *, const int *, const int *, const double *,
               int *, double *, double *, const int *, int *, double *, const int *, int *,
               const int *, int *);
}

namespace palace
{

#if !defined(MFEM_USE_MUMPS)
template <typename T>
class MumpsSchurSolverT
{
};
#endif

namespace
{

// BLR tolerance of the factorizations, at round-off: MUMPS's BLR factorization of these
// complex symmetric operators is several times faster than its full-rank one even
// without compression, at no loss of accuracy.
constexpr double kBlrTol = 1.0e-14;

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
                                       bool online)
  : space_op(space_op), online(online)
{
  env.attrs = environment_attributes;
  region.attrs = region_attributes;
  const auto &fes = space_op.GetNDSpace().Get();
  const auto &mesh = *fes.GetParMesh();
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    const int a = mesh.GetAttribute(e);
    MFEM_VERIFY(std::ranges::find(region.attrs, a) != region.attrs.end() ||
                    std::ranges::find(env.attrs, a) != env.attrs.end(),
                "Domain attribute " << a
                                    << " is neither in the region nor in the "
                                       "environment!");
  }

  // Interface DOFs: on region and environment elements; environment interior: on
  // environment elements only; region-free: on region elements (all without the Dirichlet
  // DOFs).
  mfem::Array<int> rm, em, im;
  {
    mfem::Array<int> ra(region.attrs.data(), static_cast<int>(region.attrs.size())),
        ea(env.attrs.data(), static_cast<int>(env.attrs.size()));
    MarkInterfaceTrueDofs(const_cast<mfem::ParFiniteElementSpace &>(fes), ra, ea, rm, em,
                          im);
  }
  const int nt = fes.GetTrueVSize();
  std::vector<char> dbc(nt, 0);
  for (int d : space_op.GetNDDbcTDofLists().back())
  {
    dbc[d] = 1;
  }
  is_env_int.assign(nt, 0);
  is_gamma.assign(nt, 0);
  is_region_int.assign(nt, 0);
  for (int i = 0; i < nt; i++)
  {
    is_gamma[i] = !dbc[i] && rm[i] && em[i];
    is_env_int[i] = !dbc[i] && em[i] && !rm[i];
    is_region_int[i] = !dbc[i] && rm[i] && !em[i];
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
  MPI_Comm comm = fes.GetComm();
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
  // terms border one side's system.
  for (const auto &[idx, data] : space_op.GetWavePortOp())
  {
    MFEM_VERIFY(data.active, "Driven substructuring does not support inactive wave ports "
                             "(wave port "
                                 << idx << ")!");
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
  for (auto &term : space_op.GetModalCorrectionTerms(omega))
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
    SpaceOperator::AssemblyRestriction restriction(space_op, side.attrs);
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

// For each global column (1-based) of cols, the value of sys at that true DOF on its owner:
// a request to each owner for its distinct columns.
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
    std::unique_ptr<ComplexOperator> A;
    {
      SpaceOperator::AssemblyRestriction restriction(space_op, side.attrs);
      A = (op == 0)   ? space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO)
          : (op == 1) ? space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO)
                      : space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    }
    for (int p = 0; p < 2; p++)
    {
      auto X = StealPart(A.get(), p == 1);
      if (!X)
      {
        MFEM_VERIFY(op > 0 || p > 0,
                    "Missing stiffness matrix of a substructure operator!");
        continue;
      }
      X->EliminateBC(side.pinned, Operator::DIAG_ZERO);
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
#endif
}

void DrivenSubstructure::Fill(
    const Side &side, double omega,
    const std::vector<std::unique_ptr<mfem::HypreParMatrix>> &extra,
    const std::vector<std::pair<int, WavePortOperator::ModalCorrectionTerm>> &terms,
    std::vector<std::complex<double>> &val) const
{
#if defined(MFEM_USE_MUMPS)
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
#endif
}

void DrivenSubstructure::Factor(Side &side, double omega)
{
#if defined(MFEM_USE_MUMPS)
  using Solver = MumpsSchurSolverT<std::complex<double>>;
  if (!side.schur)
  {
    Solver::Coo coo{std::move(side.irn_sys), std::move(side.jcn_sys), std::move(side.val)};
    side.irn_sys = {};
    side.jcn_sys = {};
    side.val = {};
    side.schur = std::make_unique<Solver>(
        space_op.GetComm(), side.n_sys, static_cast<int>(side.rows.size()), std::move(coo),
        side.gamma_sys, kBlrTol, 0, true, false, side.rows);
    return;
  }

  // The entries at a new frequency, in place.
  auto extra = ExtraParts(side, omega);
  MFEM_VERIFY(side.extra == (extra[0] || extra[1]),
              "The frequency-dependent part of a substructure operator changed!");
  Fill(side, omega, extra, ModalTerms(side, omega), side.schur->Values());
  side.schur->Refactor();
#endif
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

void DrivenSubstructure::Condense(double omega)
{
  MFEM_VERIFY(!online, "An online substructure is condensed with a given S_E!");
#if !defined(MFEM_USE_MUMPS)
  MFEM_ABORT("Driven substructuring requires MUMPS!");
#else
  if (!env.schur)
  {
    // Both patterns first, so that the assembled operators are released before the
    // factorizations.
    Setup(env, omega);
    Setup(region, omega);
  }
  Factor(env, omega);
  Factor(region, omega);
  if (Mpi::Root(space_op.GetComm()))
  {
    S = env.schur->Schur();
    S.resize(static_cast<std::size_t>(InterfaceSize()) * InterfaceSize());
  }
  FactorInterface();
#endif
}

void DrivenSubstructure::Condense(double omega, std::vector<std::complex<double>> &&S_env)
{
  MFEM_VERIFY(online, "An offline substructure condenses its own environment!");
#if !defined(MFEM_USE_MUMPS)
  MFEM_ABORT("Driven substructuring requires MUMPS!");
#else
  if (!region.schur)
  {
    Setup(region, omega);
  }
  Factor(region, omega);
  if (Mpi::Root(space_op.GetComm()))
  {
    MFEM_VERIFY(S_env.size() == static_cast<std::size_t>(InterfaceSize()) * InterfaceSize(),
                "Wrong size of a given S_E!");
    S = std::move(S_env);
  }
  FactorInterface();
#endif
}

void DrivenSubstructure::FactorInterface()
{
  // The interface system S_R + S_E, factored on rank 0 (complex symmetric).
#if defined(MFEM_USE_MUMPS)
  if (Mpi::Root(space_op.GetComm()))
  {
    const int nG = InterfaceSize();
    T = region.schur->Schur();
    T.resize(static_cast<std::size_t>(nG) * nG);
    for (std::size_t k = 0; k < T.size(); k++)
    {
      T[k] += S[k];
    }
    T_piv.resize(nG);
    int info = 0, lwork = -1;
    std::complex<double> wq;
    zsytrf_("L", &nG, T.data(), &nG, T_piv.data(), &wq, &lwork, &info);
    lwork = std::max(1, static_cast<int>(wq.real()));
    std::vector<std::complex<double>> work(lwork);
    zsytrf_("L", &nG, T.data(), &nG, T_piv.data(), work.data(), &lwork, &info);
    MFEM_VERIFY(info == 0, "Factorization of the interface system failed: info = " << info);
  }
#endif
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
#if !defined(MFEM_USE_MUMPS)
  MFEM_ABORT("Driven substructuring requires MUMPS!");
#else
  MFEM_VERIFY(region.schur && (online || env.schur),
              "Condense must be called before Solve!");
  MFEM_VERIFY(online == (g_env != nullptr),
              "The environment's source condensation is given online only!");
  MPI_Comm comm = space_op.GetComm();
  const int n = static_cast<int>(rhs.size()), nt = space_op.GetNDSpace().GetTrueVSize();
  const int nG = InterfaceSize();

  // Condensation of the sources onto Γ: the environment's with the interface loads b_Γ,
  // b_Γ - A_ΓE A_EE^-1 b_E (offline, or given online), and the region's, -A_ΓR A_RR^-1 b_R.
  std::vector<ComplexVector> w(online ? 0 : n), y(n);
  std::vector<const ComplexVector *> W(w.size()), Y(n);
  for (int k = 0; k < n; k++)
  {
    if (!online)
    {
      w[k] = Masked(*rhs[k], is_env_int, &is_gamma);
      W[k] = &w[k];
    }
    y[k] = Masked(*rhs[k], is_region_int);
    Y[k] = &y[k];
  }
  std::vector<std::complex<double>> rR;
  if (!online)
  {
    env.schur->Reduce(W, g_last);
  }
  else if (Mpi::Root(comm))
  {
    MFEM_VERIFY(g_env->size() == static_cast<std::size_t>(nG) * n,
                "Wrong size of the environment's source condensation!");
    g_last = *g_env;
  }
  region.schur->Reduce(Y, rR);
  if (Mpi::Root(comm))
  {
    u_last = g_last;
    for (std::size_t q = 0; q < u_last.size(); q++)
    {
      u_last[q] += rR[q];
    }
    if (nG > 0)
    {
      int info = 0;
      zsytrs_("L", &nG, &n, T.data(), &nG, T_piv.data(), u_last.data(), &nG, &info);
      MFEM_VERIFY(info == 0, "Solve of the interface system failed: info = " << info);
    }
  }

  // The interiors from the interface solution, A_II^-1 (b_I - A_IΓ u_Γ) (the interface
  // values come out of either expansion; the Dirichlet DOFs are 0).
  std::vector<ComplexVector *> Wo(w.size()), Yo(n);
  for (int k = 0; k < n; k++)
  {
    if (!online)
    {
      Wo[k] = &w[k];
    }
    Yo[k] = &y[k];
  }
  if (!online)
  {
    env.schur->Expand(u_last, Wo);
  }
  region.schur->Expand(u_last, Yo);
  u.resize(n);
  for (int k = 0; k < n; k++)
  {
    u[k].SetSize(nt);
    u[k].UseDevice(true);
    double *ur = u[k].Real().HostWrite(), *ui = u[k].Imag().HostWrite();
    const double *yr = y[k].Real().HostRead(), *yi = y[k].Imag().HostRead();
    const double *wr = online ? nullptr : w[k].Real().HostRead();
    const double *wi = online ? nullptr : w[k].Imag().HostRead();
    for (int i = 0; i < nt; i++)
    {
      const bool from_region = is_region_int[i] || (online && is_gamma[i]);
      const bool from_env = !online && (is_env_int[i] || is_gamma[i]);
      ur[i] = from_region ? yr[i] : (from_env ? wr[i] : 0.0);
      ui[i] = from_region ? yi[i] : (from_env ? wi[i] : 0.0);
    }
  }
#endif
}

void DrivenSubstructure::AddReducedBasis(const std::vector<ComplexVector> &u)
{
  MPI_Comm comm = space_op.GetComm();
  const int nt = static_cast<int>(is_gamma.size());
  if (red_ops.empty())
  {
    // The region's frequency-independent parts, as in Setup.
    red_ops.resize(6);
    for (int op = 0; op < 3; op++)
    {
      std::unique_ptr<ComplexOperator> A;
      {
        SpaceOperator::AssemblyRestriction restriction(space_op, region.attrs);
        A = (op == 0)   ? space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO)
            : (op == 1) ? space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO)
                        : space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
      }
      for (int p = 0; p < 2; p++)
      {
        red_ops[2 * op + p] = StealPart(A.get(), p == 1);
        if (red_ops[2 * op + p])
        {
          red_ops[2 * op + p]->EliminateBC(region.pinned, Operator::DIAG_ZERO);
        }
      }
    }
    red_parts.assign(6, {});
  }

  // Orthonormalize the real and imaginary parts (classical Gram-Schmidt, twice), dropping
  // the ones already in the span.
  const std::size_t r0 = red_basis.size();
  for (const auto &uk : u)
  {
    for (const Vector *part : {&uk.Real(), &uk.Imag()})
    {
      Vector v(*part);
      const double nrm0 = linalg::Norml2(comm, v);
      if (nrm0 == 0.0)
      {
        continue;
      }
      for (int pass = 0; pass < 2; pass++)
      {
        std::vector<double> d(red_basis.size());
        for (std::size_t i = 0; i < red_basis.size(); i++)
        {
          d[i] = mfem::InnerProduct(red_basis[i], v);
        }
        Mpi::GlobalSum(static_cast<int>(d.size()), d.data(), comm);
        for (std::size_t i = 0; i < red_basis.size(); i++)
        {
          v.Add(-d[i], red_basis[i]);
        }
      }
      const double nrm = linalg::Norml2(comm, v);
      if (nrm > 1.0e-10 * nrm0)
      {
        v *= 1.0 / nrm;
        red_basis.push_back(std::move(v));
      }
    }
  }

  // The projections of the parts on the new columns (the parts are symmetric).
  const std::size_t r = red_basis.size();
  for (int p = 0; p < 6; p++)
  {
    if (!red_ops[p])
    {
      continue;
    }
    std::vector<double> P(r * r, 0.0);
    for (std::size_t j = 0; j < r0; j++)
    {
      std::copy(red_parts[p].begin() + j * r0, red_parts[p].begin() + (j + 1) * r0,
                P.begin() + j * r);
    }
    Vector w(nt);
    for (std::size_t j = r0; j < r; j++)
    {
      red_ops[p]->Mult(red_basis[j], w);
      std::vector<double> d(j + 1);
      for (std::size_t i = 0; i <= j; i++)
      {
        d[i] = mfem::InnerProduct(red_basis[i], w);
      }
      Mpi::GlobalSum(static_cast<int>(d.size()), d.data(), comm);
      for (std::size_t i = 0; i <= j; i++)
      {
        P[j * r + i] = P[i * r + j] = d[i];
      }
    }
    red_parts[p] = std::move(P);
  }

  // The interface rows of the new columns on rank 0, in interface order.
  const int nG = InterfaceSize();
  for (std::size_t j = r0; j < r; j++)
  {
    std::vector<double> mine;
    for (int i = 0; i < nt; i++)
    {
      if (is_gamma[i])
      {
        mine.push_back(red_basis[j](i));
      }
    }
    const bool root = Mpi::Root(comm);
    red_gamma.resize(root ? (j + 1) * nG : 0);
    MPI_Gatherv(mine.data(), static_cast<int>(mine.size()), MPI_DOUBLE,
                root ? red_gamma.data() + j * nG : nullptr, gamma_cnt.data(),
                gamma_disp.data(), MPI_DOUBLE, 0, comm);
  }
}

void DrivenSubstructure::SolveReduced(double omega,
                                      const std::vector<const ComplexVector *> &rhs,
                                      const std::vector<std::complex<double>> &A_env,
                                      const std::vector<std::complex<double>> &b_env,
                                      std::vector<ComplexVector> &u)
{
  MPI_Comm comm = space_op.GetComm();
  const int nt = static_cast<int>(is_gamma.size()), r = ReducedDimension(),
            n = static_cast<int>(rhs.size()), nG = InterfaceSize();
  const auto coef = Coefficients(omega);

  // V^T A_R(ω) V: the parts, A2(ω) and the wave-port terms of the region, and V^T b for
  // the region interior of the right-hand sides.
  std::vector<std::complex<double>> A(static_cast<std::size_t>(r) * r, 0.0),
      b(static_cast<std::size_t>(r) * n, 0.0);
  for (int p = 0; p < 6; p++)
  {
    for (std::size_t q = 0; q < red_parts[p].size(); q++)
    {
      A[q] += coef[p] * red_parts[p][q];
    }
  }
  std::vector<std::complex<double>> loc(static_cast<std::size_t>(r) * r, 0.0);
  const auto extra = ExtraParts(region, omega);
  Vector w(nt);
  for (int p = 0; p < 2; p++)
  {
    for (int j = 0; extra[p] && j < r; j++)
    {
      extra[p]->Mult(red_basis[j], w);
      for (int i = 0; i < r; i++)
      {
        loc[static_cast<std::size_t>(j) * r + i] +=
            coef[6 + p] * mfem::InnerProduct(red_basis[i], w);
      }
    }
  }
  for (const auto &[idx, term] : ModalTerms(region, omega))
  {
    // g (V^T s)(V^T s)^T, with the local parts of V^T s summed below.
    std::vector<std::complex<double>> z(r);
    for (int i = 0; i < r; i++)
    {
      z[i] = {mfem::InnerProduct(red_basis[i], term.s->Real()),
              mfem::InnerProduct(red_basis[i], term.s->Imag())};
    }
    Mpi::GlobalSum(r, z.data(), comm);
    if (Mpi::Root(comm))
    {
      for (int j = 0; j < r; j++)
      {
        for (int i = 0; i < r; i++)
        {
          loc[static_cast<std::size_t>(j) * r + i] += term.g * z[i] * z[j];
        }
      }
    }
  }
  for (int k = 0; k < n; k++)
  {
    const auto y = Masked(*rhs[k], is_region_int);
    for (int i = 0; i < r; i++)
    {
      b[static_cast<std::size_t>(k) * r + i] = {mfem::InnerProduct(red_basis[i], y.Real()),
                                                mfem::InnerProduct(red_basis[i], y.Imag())};
    }
  }
  Mpi::GlobalSum(static_cast<int>(loc.size()), loc.data(), comm);
  Mpi::GlobalSum(static_cast<int>(b.size()), b.data(), comm);

  // The reduced system with the environment's terms, on rank 0 (complex symmetric).
  if (Mpi::Root(comm))
  {
    for (std::size_t q = 0; q < A.size(); q++)
    {
      A[q] += loc[q] + A_env[q];
    }
    for (std::size_t q = 0; q < b.size(); q++)
    {
      b[q] += b_env[q];
    }
    std::vector<int> piv(r);
    int info = 0, lwork = -1;
    std::complex<double> wq;
    zsytrf_("L", &r, A.data(), &r, piv.data(), &wq, &lwork, &info);
    lwork = std::max(1, static_cast<int>(wq.real()));
    std::vector<std::complex<double>> work(lwork);
    zsytrf_("L", &r, A.data(), &r, piv.data(), work.data(), &lwork, &info);
    MFEM_VERIFY(info == 0, "Factorization of the reduced system failed: info = " << info);
    zsytrs_("L", &r, &n, A.data(), &r, piv.data(), b.data(), &r, &info);
    MFEM_VERIFY(info == 0, "Solve of the reduced system failed: info = " << info);
  }
  Mpi::Broadcast(static_cast<int>(b.size()), b.data(), 0, comm);

  // u = V y, and the interface solution on rank 0.
  u.resize(n);
  for (int k = 0; k < n; k++)
  {
    u[k].SetSize(nt);
    u[k].UseDevice(true);
    u[k] = 0.0;
    for (int i = 0; i < r; i++)
    {
      const auto y = b[static_cast<std::size_t>(k) * r + i];
      u[k].Real().Add(y.real(), red_basis[i]);
      u[k].Imag().Add(y.imag(), red_basis[i]);
    }
  }
  if (Mpi::Root(comm))
  {
    u_last.assign(static_cast<std::size_t>(nG) * n, 0.0);
    for (int k = 0; k < n; k++)
    {
      for (int i = 0; i < r; i++)
      {
        const auto y = b[static_cast<std::size_t>(k) * r + i];
        for (int a = 0; a < nG; a++)
        {
          u_last[static_cast<std::size_t>(k) * nG + a] +=
              y * red_gamma[static_cast<std::size_t>(i) * nG + a];
        }
      }
    }
  }
}

std::vector<std::complex<double>>
DrivenSubstructure::CondenseEnvironment(const std::vector<const ComplexVector *> &b)
{
  MFEM_VERIFY(!online && env.schur, "CondenseEnvironment needs the environment factor!");
  std::vector<std::complex<double>> red;
#if defined(MFEM_USE_MUMPS)
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
  // The reduction is paired with an expansion (here of a zero interface solution).
  std::vector<std::complex<double>> zero(Mpi::Root(space_op.GetComm()) ? red.size() : 0,
                                         0.0);
  env.schur->Expand(zero, Xo);
#endif
  return red;
}

std::vector<char>
DrivenSubstructure::BoundaryTrueDofs(const std::vector<int> &bdr_attributes) const
{
  // Mark the local DOFs of the boundary elements, then the true DOFs they reach (|P|^T: a
  // true DOF shared between ranks is marked on its owner).
  auto &fes = const_cast<mfem::ParFiniteElementSpace &>(space_op.GetNDSpace().Get());
  const auto &mesh = *fes.GetParMesh();
  mfem::Vector lmark(fes.GetVSize()), tmark(fes.GetTrueVSize());
  lmark = 0.0;
  mfem::Array<int> dofs;
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    if (std::ranges::find(bdr_attributes, mesh.GetBdrAttribute(be)) != bdr_attributes.end())
    {
      fes.GetBdrElementDofs(be, dofs);
      for (int d : dofs)
      {
        lmark(d >= 0 ? d : -1 - d) = 1.0;
      }
    }
  }
  fes.Dof_TrueDof_Matrix()->AbsMultTranspose(1.0, lmark, 0.0, tmark);
  std::vector<char> mark(fes.GetTrueVSize(), 0);
  for (int i = 0; i < fes.GetTrueVSize(); i++)
  {
    mark[i] = tmark(i) > 0.0;
  }
  return mark;
}

namespace
{

// Fixed polynomial fields of the fingerprints (true DOFs).
Vector FingerprintField(const mfem::ParFiniteElementSpace &fespace, int k)
{
  auto &fes = const_cast<mfem::ParFiniteElementSpace &>(fespace);
  mfem::ParGridFunction gf(&fes);
  mfem::VectorFunctionCoefficient c(
      3,
      [k](const mfem::Vector &x, mfem::Vector &v)
      {
        const double X = x(0), Y = x(1), Z = x(2);
        const double f[3][3] = {{Z, X, Y}, {Y * Y, Z * Z, X * X}, {X * Y, Y * Z, Z * X}};
        v.SetSize(3);
        for (int d = 0; d < 3; d++)
        {
          v(d) = f[k][d];
        }
      });
  gf.ProjectCoefficient(c);
  Vector r(fes.GetTrueVSize());
  gf.GetTrueDofs(r);
  return r;
}

}  // namespace

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
    SpaceOperator::AssemblyRestriction restriction(space_op, env.attrs);
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
    r.Real() = FingerprintField(fes, k);
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
  double v[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  for (int i = 0; i < fes.GetTrueVSize(); i++)
  {
    v[0] = std::max(v[0], (is_env_int[i] && (br[i] != 0.0 || bi[i] != 0.0)) ? 1.0 : 0.0);
  }
  for (int k = 0; k < 3; k++)
  {
    const Vector r = FingerprintField(fes, k);
    const double *rr = r.HostRead();
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
  Mpi::GlobalSum(6, v + 1, space_op.GetComm());
  return {v, v + 7};
}

namespace
{

constexpr int kDrivenModelMagic = 0x31565244;  // "DRV1" (and later versions)

template <typename T>
void ReadVec(std::ifstream &f, std::vector<T> &v, std::size_t n)
{
  v.resize(n);
  f.read(reinterpret_cast<char *>(v.data()), static_cast<std::streamsize>(sizeof(T) * n));
}

}  // namespace

bool DrivenSubstructureModel::SameEnvironment(const std::vector<double> &a,
                                              const std::vector<double> &b)
{
  constexpr double tol = 1.0e-10;
  if (a.size() != b.size() || a.size() < 2 || a.size() % 2 != 0)
  {
    return false;
  }
  for (std::size_t q = 0; q < a.size(); q += (q < 2) ? 1 : 2)
  {
    const std::complex<double> za(a[q], (q < 2) ? 0.0 : a[q + 1]),
        zb(b[q], (q < 2) ? 0.0 : b[q + 1]);
    if (std::abs(za - zb) > tol * std::max(std::abs(za), std::abs(zb)))
    {
      return false;
    }
  }
  return true;
}

bool DrivenSubstructureModel::SameSource(const double *a, const double *b)
{
  constexpr double tol = 1.0e-10;
  if (a[0] != b[0])
  {
    return false;
  }
  double scale = 0.0, diff = 0.0;
  for (int q = 1; q < kSourceFp; q += 2)
  {
    const std::complex<double> za(a[q], a[q + 1]), zb(b[q], b[q + 1]);
    scale = std::max({scale, std::abs(za), std::abs(zb)});
    diff = std::max(diff, std::abs(za - zb));
  }
  return diff <= tol * scale;
}

void BarycentricInterpolant::AddSample(double omega,
                                       const std::vector<std::complex<double>> &x)
{
  // Orthonormalize x against the previous samples (classical Gram-Schmidt, twice), for the
  // factor R of the samples. The snapshots {x_j, iω_j x_j} are [Q 0; 0 Q] [R; R iΩ], so
  // their MRI weights are the right singular vector of [R; R iΩ] of least singular value.
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
  for (int j = 0, q = 0; j < n; j++)
  {
    for (int i = j; i < n; i++, q++)
    {
      A[static_cast<std::size_t>(j) * n + i] = lower[q].imag();
      f += ((i == j) ? 1.0 : 2.0) * std::norm(lower[q]);
    }
  }
  // The least eigenvalue only (LAPACK dsyevr).
  const int one = 1;
  const double zero = 0.0;
  int found = 0, info = 0, lwork = -1, liwork = -1, iwq = 0;
  double w = 0.0, wq = 0.0, z = 0.0;
  std::vector<int> isuppz(2);
  dsyevr_("N", "I", "L", &n, A.data(), &n, &zero, &zero, &one, &one, &zero, &found, &w, &z,
          &one, isuppz.data(), &wq, &lwork, &iwq, &liwork, &info);
  lwork = static_cast<int>(wq);
  liwork = iwq;
  std::vector<double> work(lwork);
  std::vector<int> iwork(liwork);
  dsyevr_("N", "I", "L", &n, A.data(), &n, &zero, &zero, &one, &one, &zero, &found, &w, &z,
          &one, isuppz.data(), work.data(), &lwork, iwork.data(), &liwork, &info);
  MFEM_VERIFY(info == 0 && found == 1, "Eigenvalue of Im S_E failed: info = " << info);
  return (f > 0.0) ? w / std::sqrt(f) : 0.0;
}

double BarycentricInterpolant::FindMaxError() const
{
  MFEM_VERIFY(z.size() >= 2, "The next sample needs two samples to bound the band!");
  const auto [lo, hi] = std::ranges::minmax(z);
  constexpr int n = 1000000;
  double best = lo, dmin = std::numeric_limits<double>::infinity();
  for (int i = 1; i < n; i++)
  {
    const double omega = lo + (hi - lo) * i / n;
    std::complex<double> d = 0.0;
    for (std::size_t j = 0; j < z.size(); j++)
    {
      d += w[j] / (omega - z[j]);
    }
    if (std::abs(d) < dmin)
    {
      dmin = std::abs(d);
      best = omega;
    }
  }
  return best;
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

DrivenSubstructureModel::Record
DrivenSubstructureModel::Unflatten(const std::vector<std::complex<double>> &x) const
{
  const auto sz = PartSizes();
  MFEM_VERIFY(x.size() == sz[0] + sz[1] + sz[2] + sz[3], "Wrong size of a record!");
  Record r;
  const std::size_t n = nG;
  r.S.resize(n * n);
  for (std::size_t c = 0, q = 0; c < n; c++)
  {
    for (std::size_t i = c; i < n; i++, q++)
    {
      r.S[c * n + i] = r.S[i * n + c] = x[q];
    }
  }
  auto it = x.begin() + static_cast<std::ptrdiff_t>(sz[0]);
  r.g.assign(it, it + static_cast<std::ptrdiff_t>(sz[1]));
  it += static_cast<std::ptrdiff_t>(sz[1]);
  r.h.assign(it, it + static_cast<std::ptrdiff_t>(sz[2]));
  it += static_cast<std::ptrdiff_t>(sz[2]);
  r.c.assign(it, it + static_cast<std::ptrdiff_t>(sz[3]));
  return r;
}

std::size_t DrivenSubstructureModel::HeaderBytes() const
{
  // Version 1: 8 ints and the frequencies; version 2: 10 ints, and room for capacity
  // frequencies and weights.
  const std::size_t common = sizeof(double) * signatures.size() +
                             sizeof(double) * env_fp.size() +
                             sizeof(int) * excitations.size() +
                             sizeof(double) * exc_fp.size() + sizeof(int) * ports.size();
  return (version == 1) ? sizeof(int) * 8 + common + sizeof(double) * omega.size()
                        : sizeof(int) * 10 + common +
                              (sizeof(double) + sizeof(std::complex<double>)) * capacity;
}

std::size_t DrivenSubstructureModel::RecordBytes() const
{
  const auto sz = PartSizes();
  return sizeof(std::complex<double>) * (sz[0] + sz[1] + sz[2] + sz[3]);
}

void DrivenSubstructureModel::WriteHeader(const std::string &path, bool create) const
{
  MFEM_VERIFY(version == 2 && capacity >= static_cast<int>(omega.size()) &&
                  (weights.empty() || weights.size() == omega.size()),
              "Invalid substructuring model header!");
  std::fstream f(path, create ? (std::ios::binary | std::ios::out | std::ios::trunc)
                              : (std::ios::binary | std::ios::in | std::ios::out));
  MFEM_VERIFY(f.good(), "Cannot write the substructuring model \"" << path << "\"!");
  const int head[10] = {kDrivenModelMagic,
                        2,
                        nG,
                        sig_w,
                        static_cast<int>(env_fp.size()),
                        static_cast<int>(excitations.size()),
                        static_cast<int>(ports.size()),
                        capacity,
                        static_cast<int>(omega.size()),
                        weights.empty() ? 0 : 1};
  std::vector<double> om(omega);
  std::vector<std::complex<double>> wt(weights);
  om.resize(capacity, 0.0);
  wt.resize(capacity, 0.0);
  f.write(reinterpret_cast<const char *>(head), sizeof(head));
  for (const auto *v : {&signatures, &env_fp})
  {
    f.write(reinterpret_cast<const char *>(v->data()),
            static_cast<std::streamsize>(sizeof(double) * v->size()));
  }
  f.write(reinterpret_cast<const char *>(excitations.data()),
          static_cast<std::streamsize>(sizeof(int) * excitations.size()));
  f.write(reinterpret_cast<const char *>(exc_fp.data()),
          static_cast<std::streamsize>(sizeof(double) * exc_fp.size()));
  f.write(reinterpret_cast<const char *>(ports.data()),
          static_cast<std::streamsize>(sizeof(int) * ports.size()));
  f.write(reinterpret_cast<const char *>(om.data()),
          static_cast<std::streamsize>(sizeof(double) * om.size()));
  f.write(reinterpret_cast<const char *>(wt.data()),
          static_cast<std::streamsize>(sizeof(std::complex<double>) * wt.size()));
  MFEM_VERIFY(f.good(), "Cannot write the substructuring model \"" << path << "\"!");
}

void DrivenSubstructureModel::AppendRecord(const std::string &path, const Record &r) const
{
  std::ofstream f(path, std::ios::binary | std::ios::app);
  MFEM_VERIFY(f.good(), "Cannot write the substructuring model \"" << path << "\"!");
  const auto x = Flatten(r);
  f.write(reinterpret_cast<const char *>(x.data()),
          static_cast<std::streamsize>(sizeof(std::complex<double>) * x.size()));
}

void DrivenSubstructureModel::ReadHeader(const std::string &path, MPI_Comm comm)
{
  int head[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  std::ifstream f;
  if (Mpi::Root(comm))
  {
    f.open(path, std::ios::binary);
    if (f.good())
    {
      f.read(reinterpret_cast<char *>(head), sizeof(int) * 2);
      f.read(reinterpret_cast<char *>(head + 2), sizeof(int) * ((head[1] == 1) ? 6 : 8));
    }
  }
  Mpi::Broadcast(10, head, 0, comm);
  MFEM_VERIFY(head[0] == kDrivenModelMagic && (head[1] == 1 || head[1] == 2),
              "Cannot read the driven substructuring model \""
                  << path << "\" (run in \"Offline\" mode with \"SaveModel\" first)!");
  version = head[1];
  nG = head[2];
  sig_w = head[3];
  auto read = [&](auto &v, std::size_t n)
  {
    if (Mpi::Root(comm))
    {
      ReadVec(f, v, n);
    }
    else
    {
      v.resize(n);
    }
    Mpi::Broadcast(static_cast<int>(n), v.data(), 0, comm);
  };
  read(signatures, static_cast<std::size_t>(nG) * sig_w);
  read(env_fp, head[4]);
  read(excitations, head[5]);
  read(exc_fp, kSourceFp * static_cast<std::size_t>(head[5]));
  read(ports, head[6]);
  capacity = head[7];
  read(omega, capacity);
  weights.clear();
  if (version == 2)
  {
    read(weights, capacity);
    omega.resize(head[8]);
    weights.resize(head[9] ? head[8] : 0);
  }
}

int DrivenSubstructureModel::NumRecords(const std::string &path) const
{
  std::ifstream f(path, std::ios::binary | std::ios::ate);
  if (!f.good())
  {
    return 0;
  }
  const auto size = static_cast<std::size_t>(f.tellg());
  return (size < HeaderBytes()) ? 0
                                : static_cast<int>((size - HeaderBytes()) / RecordBytes());
}

std::vector<std::complex<double>>
DrivenSubstructureModel::ReadFlatRecord(const std::string &path, int j) const
{
  std::ifstream f(path, std::ios::binary);
  f.seekg(static_cast<std::streamoff>(HeaderBytes() + RecordBytes() * j));
  std::vector<std::complex<double>> x;
  ReadVec(f, x, RecordBytes() / sizeof(std::complex<double>));
  MFEM_VERIFY(f.good(), "Truncated substructuring model \"" << path << "\"!");
  return x;
}

}  // namespace palace
