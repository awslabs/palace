// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "drivensubstructure.hpp"

#include <algorithm>
#include <array>
#include <fstream>
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

std::vector<std::complex<double>>
DrivenSubstructure::CondenseEnvironment(const std::vector<const ComplexVector *> &b)
{
  MFEM_VERIFY(!online && env.schur, "CondenseEnvironment needs the environment factor!");
  std::vector<std::complex<double>> red;
#if defined(MFEM_USE_MUMPS)
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

constexpr int kDrivenModelMagic = 0x31565244;  // "DRV1"

template <typename T>
void WriteVec(std::ofstream &f, const std::vector<T> &v)
{
  f.write(reinterpret_cast<const char *>(v.data()),
          static_cast<std::streamsize>(sizeof(T) * v.size()));
}

template <typename T>
void ReadVec(std::ifstream &f, std::vector<T> &v, std::size_t n)
{
  v.resize(n);
  f.read(reinterpret_cast<char *>(v.data()), static_cast<std::streamsize>(sizeof(T) * n));
}

}  // namespace

std::size_t DrivenSubstructureModel::HeaderBytes() const
{
  return sizeof(int) * 8 + sizeof(double) * signatures.size() +
         sizeof(double) * env_fp.size() + sizeof(int) * excitations.size() +
         sizeof(double) * exc_fp.size() + sizeof(int) * ports.size() +
         sizeof(double) * omega.size();
}

std::size_t DrivenSubstructureModel::RecordBytes() const
{
  const std::size_t n = nG, ne = excitations.size(), np = ports.size();
  return sizeof(std::complex<double>) * (n * (n + 1) / 2 + n * ne + n * np + np * ne);
}

void DrivenSubstructureModel::WriteHeader(const std::string &path) const
{
  std::ofstream f(path, std::ios::binary | std::ios::trunc);
  MFEM_VERIFY(f.good(), "Cannot write the substructuring model \"" << path << "\"!");
  const int head[8] = {kDrivenModelMagic,
                       1,
                       nG,
                       sig_w,
                       static_cast<int>(env_fp.size()),
                       static_cast<int>(excitations.size()),
                       static_cast<int>(ports.size()),
                       static_cast<int>(omega.size())};
  f.write(reinterpret_cast<const char *>(head), sizeof(head));
  WriteVec(f, signatures);
  WriteVec(f, env_fp);
  WriteVec(f, excitations);
  WriteVec(f, exc_fp);
  WriteVec(f, ports);
  WriteVec(f, omega);
}

void DrivenSubstructureModel::AppendRecord(const std::string &path, const Record &r) const
{
  // S_E is symmetric: its lower triangle, by columns.
  std::vector<std::complex<double>> lower;
  lower.reserve(static_cast<std::size_t>(nG) * (nG + 1) / 2);
  for (int j = 0; j < nG; j++)
  {
    for (int i = j; i < nG; i++)
    {
      lower.push_back(r.S[static_cast<std::size_t>(j) * nG + i]);
    }
  }
  std::ofstream f(path, std::ios::binary | std::ios::app);
  MFEM_VERIFY(f.good(), "Cannot write the substructuring model \"" << path << "\"!");
  WriteVec(f, lower);
  WriteVec(f, r.g);
  WriteVec(f, r.h);
  WriteVec(f, r.c);
}

void DrivenSubstructureModel::ReadHeader(const std::string &path, MPI_Comm comm)
{
  int head[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  std::ifstream f;
  if (Mpi::Root(comm))
  {
    f.open(path, std::ios::binary);
    if (f.good())
    {
      f.read(reinterpret_cast<char *>(head), sizeof(head));
    }
  }
  Mpi::Broadcast(8, head, 0, comm);
  MFEM_VERIFY(head[0] == kDrivenModelMagic && head[1] == 1,
              "Cannot read the driven substructuring model \""
                  << path << "\" (run in \"Offline\" mode with \"SaveModel\" first)!");
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
  read(omega, head[7]);
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

DrivenSubstructureModel::Record DrivenSubstructureModel::ReadRecord(const std::string &path,
                                                                    int j) const
{
  std::ifstream f(path, std::ios::binary);
  f.seekg(static_cast<std::streamoff>(HeaderBytes() + RecordBytes() * j));
  std::vector<std::complex<double>> lower;
  Record r;
  const std::size_t n = nG, ne = excitations.size(), np = ports.size();
  ReadVec(f, lower, n * (n + 1) / 2);
  ReadVec(f, r.g, n * ne);
  ReadVec(f, r.h, n * np);
  ReadVec(f, r.c, np * ne);
  MFEM_VERIFY(f.good(), "Truncated substructuring model \"" << path << "\"!");
  r.S.resize(n * n);
  for (std::size_t c = 0, q = 0; c < n; c++)
  {
    for (std::size_t i = c; i < n; i++, q++)
    {
      r.S[c * n + i] = r.S[i * n + c] = lower[q];
    }
  }
  return r;
}

}  // namespace palace
