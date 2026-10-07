// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "drivensubstructure.hpp"

#include <algorithm>
#include <array>
#include <utility>
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
                                       const std::vector<int> &environment_attributes)
  : space_op(space_op)
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
  // The frequency-independent parts, assembled once and one at a time: the pattern of the
  // stiffness matrix (with the diagonal of the pinned DOFs: a side-restricted operator
  // stores no entries on the DOFs no element of its side touches), which contains those
  // of the other parts, coupling DOFs of elements of the side.
  side.parts.assign(6, {});
  HYPRE_BigInt row0 = 0;
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
        std::vector<std::vector<double>> vals;
        LowerTrianglePattern({X.get()}, side.pinned, side.irn, side.jcn, vals,
                             side.row_ptr);
        part = std::move(vals[0]);
        row0 = X->GetRowStarts()[0];
      }
      else
      {
        part.assign(side.irn.size(), 0.0);
        AddToPattern(*X, 1.0, side.row_ptr, side.jcn, part);
      }
    }
  }
  side.unit.clear();
  for (int i : side.pinned)
  {
    const auto first = side.jcn.begin() + side.row_ptr[i],
               last = side.jcn.begin() + side.row_ptr[i + 1];
    side.unit.push_back(static_cast<int>(
        std::lower_bound(first, last, static_cast<int>(row0 + i + 1)) - side.jcn.begin()));
  }
  side.val.resize(side.irn.size());
  auto extra = ExtraParts(side, omega);
  side.extra = extra[0] || extra[1];
  Fill(side, omega, extra, side.jcn, side.val);
#endif
}

void DrivenSubstructure::Fill(
    const Side &side, double omega,
    const std::vector<std::unique_ptr<mfem::HypreParMatrix>> &extra,
    const std::vector<int> &jcn, std::vector<std::complex<double>> &val) const
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
  for (int q : side.unit)
  {
    val[q] += 1.0;
  }
  for (int p = 0; p < 2; p++)
  {
    if (extra[p])
    {
      AddToPattern(*extra[p], coef[6 + p], side.row_ptr, jcn, val);
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
    const auto &fes = space_op.GetNDSpace().Get();
    Solver::Coo coo{std::move(side.irn), std::move(side.jcn), std::move(side.val)};
    side.irn = {};
    side.jcn = {};
    side.val = {};
    side.schur =
        std::make_unique<Solver>(fes.GetComm(), fes.GlobalTrueVSize(), fes.GetTrueVSize(),
                                 std::move(coo), gamma_tdofs, kBlrTol, false, true, false);
    return;
  }

  // The entries at a new frequency, in place.
  auto extra = ExtraParts(side, omega);
  MFEM_VERIFY(side.extra == (extra[0] || extra[1]),
              "The frequency-dependent part of a substructure operator changed!");
  Fill(side, omega, extra, side.schur->Columns(), side.schur->Values());
  side.schur->Refactor();
#endif
}

void DrivenSubstructure::Condense(double omega)
{
#if !defined(MFEM_USE_MUMPS)
  MFEM_ABORT("Driven substructuring requires MUMPS!");
#else
  MPI_Comm comm = space_op.GetComm();
  const int nG = InterfaceSize();
  if (!env.schur)
  {
    // Both patterns first, so that the assembled operators are released before the
    // factorizations.
    Setup(env, omega);
    Setup(region, omega);
  }
  Factor(env, omega);
  Factor(region, omega);

  // The interface system S_R + S_E, factored on rank 0 (complex symmetric).
  if (Mpi::Root(comm))
  {
    S = env.schur->Schur();
    T = region.schur->Schur();
    S.resize(static_cast<std::size_t>(nG) * nG);
    T.resize(static_cast<std::size_t>(nG) * nG);
  }
  if (Mpi::Root(comm))
  {
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

void DrivenSubstructure::Solve(const std::vector<const ComplexVector *> &rhs,
                               std::vector<ComplexVector> &u)
{
#if !defined(MFEM_USE_MUMPS)
  MFEM_ABORT("Driven substructuring requires MUMPS!");
#else
  MFEM_VERIFY(env.schur && region.schur, "Condense must be called before Solve!");
  MPI_Comm comm = space_op.GetComm();
  const int n = static_cast<int>(rhs.size()), nt = space_op.GetNDSpace().GetTrueVSize();
  const int nG = InterfaceSize();
  auto masked = [&](const ComplexVector &b, const std::vector<char> &mask, bool gamma)
  {
    ComplexVector y(nt);
    const double *br = b.Real().HostRead(), *bi = b.Imag().HostRead();
    double *yr = y.Real().HostWrite(), *yi = y.Imag().HostWrite();
    for (int i = 0; i < nt; i++)
    {
      const bool keep = mask[i] || (gamma && is_gamma[i]);
      yr[i] = keep ? br[i] : 0.0;
      yi[i] = keep ? bi[i] : 0.0;
    }
    return y;
  };

  // Condensation of the sources onto Γ: the environment's with the interface loads b_Γ,
  // b_Γ - A_ΓE A_EE^-1 b_E, and the region's, -A_ΓR A_RR^-1 b_R.
  std::vector<ComplexVector> w(n), y(n);
  std::vector<const ComplexVector *> W(n), Y(n);
  for (int k = 0; k < n; k++)
  {
    w[k] = masked(*rhs[k], is_env_int, true);
    y[k] = masked(*rhs[k], is_region_int, false);
    W[k] = &w[k];
    Y[k] = &y[k];
  }
  std::vector<std::complex<double>> rG, rR;
  env.schur->Reduce(W, rG);
  region.schur->Reduce(Y, rR);
  if (Mpi::Root(comm) && nG > 0)
  {
    for (std::size_t q = 0; q < rG.size(); q++)
    {
      rG[q] += rR[q];
    }
    int info = 0;
    zsytrs_("L", &nG, &n, T.data(), &nG, T_piv.data(), rG.data(), &nG, &info);
    MFEM_VERIFY(info == 0, "Solve of the interface system failed: info = " << info);
  }

  // Both interiors from the interface solution, A_II^-1 (b_I - A_IΓ u_Γ).
  std::vector<ComplexVector *> Wo(n), Yo(n);
  for (int k = 0; k < n; k++)
  {
    Wo[k] = &w[k];
    Yo[k] = &y[k];
  }
  env.schur->Expand(rG, Wo);
  region.schur->Expand(rG, Yo);
  u.resize(n);
  for (int k = 0; k < n; k++)
  {
    u[k].SetSize(nt);
    u[k].UseDevice(true);
    double *ur = u[k].Real().HostWrite(), *ui = u[k].Imag().HostWrite();
    const double *wr = w[k].Real().HostRead(), *wi = w[k].Imag().HostRead();
    const double *yr = y[k].Real().HostRead(), *yi = y[k].Imag().HostRead();
    for (int i = 0; i < nt; i++)
    {
      // The interface values come out of either expansion; the Dirichlet DOFs are 0.
      const bool region = is_region_int[i];
      ur[i] = (is_env_int[i] || is_gamma[i]) ? wr[i] : (region ? yr[i] : 0.0);
      ui[i] = (is_env_int[i] || is_gamma[i]) ? wi[i] : (region ? yi[i] : 0.0);
    }
  }
#endif
}

}  // namespace palace
