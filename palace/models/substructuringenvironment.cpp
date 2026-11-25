// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "substructuringsolverimpl.hpp"

#include <cmath>
#include <span>
#include "linalg/lowertriangle.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace
{

// A row of w doubles, so that the counts of a row exchange are rows (n w entries can
// overflow an int).
class RowType
{
public:
  explicit RowType(int w)
  {
    MPI_Type_contiguous(w, MPI_DOUBLE, &type);
    MPI_Type_commit(&type);
  }
  ~RowType() { MPI_Type_free(&type); }
  RowType(const RowType &) = delete;
  RowType &operator=(const RowType &) = delete;
  MPI_Datatype type;
};

}  // namespace

void SubstructuringSolver::Impl::ComputeSourceModes(const std::vector<int> &ids,
                                                    const std::vector<Source> &src)
{
  EnsureEnv(true);
  MPI_Comm comm = parent_fes.GetComm();
  const int K = static_cast<int>(ids.size());
  src_ids = ids;
  src_fp.assign(2 * K, 0.0);
  std::vector<mfem::Vector> b(K), W(K, mfem::Vector(nt)), q(K, mfem::Vector(nt));
  std::vector<const mfem::Vector *> X(K);
  std::vector<mfem::Vector *> Y(K);
  for (int k = 0; k < K; k++)
  {
    std::array<double, 3> fp;
    b[k] = EnvSource(src[k], fp);
    src_fp[2 * k] = fp[1];
    src_fp[2 * k + 1] = fp[2];
    W[k] = 0.0;
    X[k] = &b[k];
    Y[k] = &W[k];
  }
  if (K > 0)
  {
    ApplyAeeInvMulti(X, Y);
  }
  mfem::Vector t(nt);
  auto rows = [&](std::vector<double> &out, int k, const mfem::Vector &v)
  {
    const double *vv = v.HostRead();
    for (int i = 0; i < nt; i++)
    {
      if (iface.gamma[i])
      {
        out[static_cast<std::size_t>(iface.index[i] - iface.offset) * K + k] = vv[i];
      }
    }
  };
  SG_rows.assign(static_cast<std::size_t>(iface.count) * K, 0.0);
  std::vector<mfem::Vector> KW(K, mfem::Vector(nt)), MD(K, mfem::Vector(nt)),
      D(K, mfem::Vector(nt));
  // By condensations: h_k = Z^T q_k for q_k = (A_env - D_env) W_k - b^E_k, and with S^Kin
  // h^kin_k = Z^T M_kin^E D_k for D_k = W_k - a_k. The London mass is in A_env, the sheet
  // mass not: the energies c take the sheet loads MD, the kinetic parts the loads MK.
  const bool kin = have_skin;
  std::vector<mfem::Vector> MK(kin && MlondonEnv() ? K : 0, mfem::Vector(nt));
  std::vector<const mfem::Vector *> V(kin ? 2 * K : K);
  for (int k = 0; k < K; k++)
  {
    Aenv().Mult(W[k], t);
    t.Neg();
    t += b[k];  // g_k = b^E_k|_Gamma - (A_env W_k)|_Gamma
    rows(SG_rows, k, t);
    Aenv().Mult(W[k], q[k]);
    Denv().Mult(W[k], t);
    q[k] -= t;
    q[k] -= b[k];
    V[k] = &q[k];
    D[k] = W[k];
    if (src[k].HasA())
    {
      D[k] -= *src[k].a;
    }
    MD[k] = 0.0;
    if (has_sheets)
    {
      MsheetEnv().Mult(D[k], MD[k]);
    }
    if (kin)
    {
      V[K + k] = &MD[k];
      if (!MK.empty())
      {
        KineticEnvMult(D[k], MK[k]);
        V[K + k] = &MK[k];
      }
    }
  }
  std::vector<double> hrows;
  CondenseRows(V, hrows);
  const int nv = static_cast<int>(V.size());
  SH_rows.assign(static_cast<std::size_t>(iface.count) * K, 0.0);
  SHkin_rows.assign(kin ? static_cast<std::size_t>(iface.count) * K : 0, 0.0);
  for (int r = 0; r < iface.count; r++)
  {
    for (int k = 0; k < K; k++)
    {
      const std::size_t q0 = static_cast<std::size_t>(r) * nv,
                        rk = static_cast<std::size_t>(r) * K + k;
      SH_rows[rk] = hrows[q0 + k];
      if (kin)
      {
        SHkin_rows[rk] = hrows[q0 + K + k];
      }
    }
  }
  SC.assign(static_cast<std::size_t>(K) * K, 0.0);
  SCkin.assign(kin ? static_cast<std::size_t>(K) * K : 0, 0.0);
  std::vector<mfem::Vector> JE(K, mfem::Vector(nt));
  for (int k = 0; k < K; k++)
  {
    Kenv().Mult(W[k], KW[k]);
    JE[k] = 0.0;
    if (src[k].HasJ())
    {
      const double *J = src[k].J->HostRead();
      double *je = JE[k].HostWrite();
      for (int i = 0; i < nt; i++)
      {
        je[i] = iface.env_int[i] ? J[i] : 0.0;
      }
    }
  }
  SB.assign(static_cast<std::size_t>(K) * K, 0.0);
  for (int i = 0; i < K; i++)
  {
    for (int j = 0; j < K; j++)
    {
      const std::size_t ij = static_cast<std::size_t>(i) * K + j;
      const double sheets = D[i] * MD[j];
      SC[ij] = (W[i] * KW[j]) + sheets - (JE[i] * W[j]) - (JE[j] * W[i]);
      SB[ij] = b[i] * W[j];
      if (kin)
      {
        SCkin[ij] = MK.empty() ? sheets : D[i] * MK[j];
      }
    }
  }
  if (K > 0)
  {
    MPI_Allreduce(MPI_IN_PLACE, SC.data(), K * K, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(MPI_IN_PLACE, SB.data(), K * K, MPI_DOUBLE, MPI_SUM, comm);
    if (kin)
    {
      MPI_Allreduce(MPI_IN_PLACE, SCkin.data(), K * K, MPI_DOUBLE, MPI_SUM, comm);
    }
  }
  src_ready = true;
}

bool SubstructuringSolver::Impl::SourceModesMatch(const std::vector<int> &ids,
                                                  const std::vector<Source> &src,
                                                  std::vector<int> &col,
                                                  std::vector<char> &need)
{
  col.assign(ids.size(), -1);
  need.assign(ids.size(), 0);
  // On the mesh of the region, whether each source has DOFs in it: one without them is of
  // the environment, whose modes the model must hold.
  std::vector<double> reach(region_only ? ids.size() : 0, 0.0);
  for (std::size_t j = 0; j < reach.size(); j++)
  {
    reach[j] = std::max(src[j].HasA() ? src[j].a->Normlinf() : 0.0,
                        src[j].HasJ() ? src[j].J->Normlinf() : 0.0);
  }
  if (!reach.empty())
  {
    MPI_Allreduce(MPI_IN_PLACE, reach.data(), static_cast<int>(reach.size()), MPI_DOUBLE,
                  MPI_MAX, parent_fes.GetComm());
  }
  bool ok = true;
  for (std::size_t j = 0; j < ids.size(); j++)
  {
    std::array<double, 3> fp;
    (void)EnvSource(src[j], fp);  // collective: always evaluated
    for (std::size_t k = 0; src_ready && k < src_ids.size(); k++)
    {
      if (src_ids[k] == ids[j])
      {
        col[j] = static_cast<int>(k);
      }
    }
    if ((region_only || model_sources.contains(ids[j])) && col[j] >= 0)
    {
      // The model holds the environment and its sources, also those without a generator
      // here (left out of the mesh of the region).
      need[j] = 1;
      continue;
    }
    MFEM_VERIFY(!region_only || reach[j] > 0.0,
                (ids[j] > 0 ? "Excitation " : "The aperture of surface current ")
                    << std::abs(ids[j])
                    << " has no DOFs in the mesh of the region and is not in the saved "
                       "substructuring model: run on the full mesh!");
    if (fp[0] == 0.0)
    {
      col[j] = -1;
      continue;
    }
    need[j] = 1;
    if (col[j] < 0)
    {
      ok = false;
      continue;
    }
    ok = ok && PairingsMatch(std::span<const double>(fp).subspan(1),
                             std::span<const double>(src_fp).subspan(2 * col[j], 2),
                             kFingerprintTol);
  }
  return ok;
}

void SubstructuringSolver::Impl::EnsureEnv(bool condense)
{
  MFEM_VERIFY(!region_only,
              "Substructuring on the mesh of the region cannot solve the environment (an "
              "excitation not in the saved model, or general Dirichlet data): run on the "
              "full mesh!");
  if (env_schur && (env_condenses || !condense))
  {
    return;
  }
  BlockTimer bt(Timer::KSP_SETUP);
  // The environment's side: its interior and, to condense onto it, Gamma as the Schur
  // variables.
  std::vector<char> keep(nt);
  for (int i = 0; i < nt; i++)
  {
    keep[i] = iface.env_int[i] || (condense && iface.gamma[i]);
  }
  SideSystem side({&Aenv()}, keep);
  const LinearSolver type =
      SchurSolverFor(iodata.solver.substructuring->direct_solver, SchurKind::ENVIRONMENT);
  env_factorization = SchurFactorizationName(type);
  const auto gamma = condense ? side.Indices(iface.tdofs) : std::vector<HYPRE_BigInt>{};
  MFEM_VERIFY(gamma_grid, "The environment is factored before the interface is numbered!");
  env_schur.reset();
  env_schur = MakeSchurFactorization<double>(
      type, iodata.solver.linear.sym_factorization, parent_fes.GetComm(), side.n_glob,
      static_cast<int>(side.rows.size()),
      {std::move(side.irn), std::move(side.jcn), std::move(side.parts[0])}, gamma,
      side.rows, true, false, condense ? gamma_grid.get() : nullptr);
  env_condenses = condense;
}

void SubstructuringSolver::Impl::RecoverEnvInterior(std::vector<mfem::Vector> &u,
                                                    const std::vector<mfem::Vector> *f)
{
  if (region_only)
  {
    Mpi::Warning(parent_fes.GetComm(),
                 "Substructuring on the mesh of the region writes the fields in the "
                 "region and on the interface only!\n");
    return;
  }
  EnsureEnv();
  const int n = static_cast<int>(u.size());
  std::vector<mfem::Vector> rhs(n, mfem::Vector(nt));
  std::vector<const mfem::Vector *> X(n);
  std::vector<mfem::Vector *> Y(n);
  for (int k = 0; k < n; k++)
  {
    Aenv().Mult(u[k], rhs[k]);
    rhs[k].Neg();
    if (f)
    {
      rhs[k] += (*f)[k];
    }
    X[k] = &rhs[k];
    Y[k] = &rhs[k];  // solve in place
  }
  if (n > 0)
  {
    ApplyAeeInvMulti(X, Y);
  }
  for (int k = 0; k < n; k++)
  {
    const double *r = rhs[k].HostRead();
    double *uu = u[k].HostReadWrite();
    for (int i = 0; i < nt; i++)
    {
      if (iface.env_int[i])
      {
        uu[i] = r[i];
      }
    }
  }
}

mfem::DenseMatrix SubstructuringSolver::Impl::EnergyGram(const std::vector<mfem::Vector> &u)
{
  const int n = static_cast<int>(u.size());
  mfem::DenseMatrix E(n);
  mfem::Vector Ku(nt), t(nt);
  for (int j = 0; j < n; j++)
  {
    K_region_e->Mult(u[j], Ku);
    Kenv().Mult(u[j], t);
    Ku += t;
    for (int i = 0; i < n; i++)
    {
      E(i, j) = u[i] * Ku;
    }
  }
  if (n > 0)
  {
    MPI_Allreduce(MPI_IN_PLACE, E.GetData(), n * n, MPI_DOUBLE, MPI_SUM,
                  parent_fes.GetComm());
  }
  return E;
}

std::array<double, 2>
SubstructuringSolver::Impl::LiftFingerprint(const mfem::Vector &x,
                                            const mfem::Vector &t) const
{
  const mfem::Vector &r = CachedFingerprintField(0);
  const double rAr = EnvironmentFingerprint()[2];  // r^T A_env r
  double loc[2] = {x * t, r * t}, glob[2];
  MPI_Allreduce(loc, glob, 2, MPI_DOUBLE, MPI_SUM, parent_fes.GetComm());
  return {glob[0], glob[0] + 2.0 * glob[1] + rAr};
}

void SubstructuringSolver::Impl::ComputeModes(const std::vector<int> &ids,
                                              const std::vector<mfem::Vector> &x)
{
  EnsureEnv();
  MPI_Comm comm = parent_fes.GetComm();
  const int K = static_cast<int>(ids.size());
  mode_ids = ids;
  std::vector<mfem::Vector> t(K, mfem::Vector(nt)), w(K, mfem::Vector(nt));
  std::vector<int> touches(K, 0);
  mode_fp.assign(2 * K, 0.0);
  {
    std::vector<const mfem::Vector *> X(K);
    for (int k = 0; k < K; k++)
    {
      X[k] = &x[k];
    }
    t = EnvApply(X);
  }
  for (int k = 0; k < K; k++)
  {
    const auto fp = LiftFingerprint(x[k], t[k]);
    mode_fp[2 * k] = fp[0];
    mode_fp[2 * k + 1] = fp[1];
    const double *tt = t[k].HostRead();
    for (int i = 0; i < nt && !touches[k]; i++)
    {
      touches[k] = (iface.env_int[i] && tt[i] != 0.0);
    }
  }
  if (K > 0)
  {
    MPI_Allreduce(MPI_IN_PLACE, touches.data(), K, MPI_INT, MPI_MAX, comm);
  }
  // w_k = A_EE^-1 t_k|_E, for the lifts that touch the environment interior (zero for the
  // others).
  {
    std::vector<const mfem::Vector *> X;
    std::vector<mfem::Vector *> Y;
    for (int k = 0; k < K; k++)
    {
      w[k] = 0.0;
      if (touches[k])
      {
        X.push_back(&t[k]);
        Y.push_back(&w[k]);
      }
    }
    if (!X.empty())
    {
      ApplyAeeInvMulti(X, Y);
    }
  }
  // g_k = (t_k - A_env w_k)|_Gamma and C_kl = x_k . t_l - t_k . w_l.
  mfem::Vector tmp(nt);
  G_rows.assign(static_cast<std::size_t>(iface.count) * K, 0.0);
  for (int k = 0; k < K; k++)
  {
    Aenv().Mult(w[k], tmp);
    const double *tk = t[k].HostRead(), *tt = tmp.HostRead();
    for (int i = 0; i < nt; i++)
    {
      if (iface.gamma[i])
      {
        G_rows[static_cast<std::size_t>(iface.index[i] - iface.offset) * K + k] =
            tk[i] - tt[i];
      }
    }
  }
  Cmode.assign(static_cast<std::size_t>(K) * K, 0.0);
  for (int k = 0; k < K; k++)
  {
    for (int l = 0; l < K; l++)
    {
      Cmode[static_cast<std::size_t>(k) * K + l] = (x[k] * t[l]) - (t[k] * w[l]);
    }
  }
  if (K > 0)
  {
    MPI_Allreduce(MPI_IN_PLACE, Cmode.data(), K * K, MPI_DOUBLE, MPI_SUM, comm);
  }
  modes_ready = true;
}

bool SubstructuringSolver::Impl::ModesMatch(const std::vector<int> &ids,
                                            const std::vector<mfem::Vector> &xd,
                                            std::vector<int> &col,
                                            std::vector<char> &need) const
{
  const std::size_t n = ids.size();
  col.assign(n, -1);
  need.assign(n, 0);
  // Whether each excitation has DOFs in the mesh (on the mesh of the region, one without
  // them is of the environment) and its lift reaches the environment (A_env x != 0), and,
  // on the full mesh, the fingerprint of its lift.
  constexpr std::size_t kBlock = 4;  // bounds the memory of the full vectors
  std::vector<std::array<double, 2>> fps(n);
  std::vector<double> reach(2 * n, 0.0);
  for (std::size_t j0 = 0; j0 < n; j0 += kBlock)
  {
    std::vector<mfem::Vector> x;
    std::vector<const mfem::Vector *> X;
    for (std::size_t j = j0; j < std::min(n, j0 + kBlock); j++)
    {
      x.push_back(DirichletField(xd[j]));
    }
    for (const auto &v : x)
    {
      X.push_back(&v);
    }
    const std::vector<mfem::Vector> t = EnvApply(X);
    for (std::size_t j = j0; j < j0 + x.size(); j++)
    {
      reach[j] = t[j - j0].Normlinf();
      reach[n + j] = xd[j].Normlinf();
      if (!region_only)
      {
        fps[j] = LiftFingerprint(x[j - j0], t[j - j0]);
      }
    }
  }
  if (n > 0)
  {
    MPI_Allreduce(MPI_IN_PLACE, reach.data(), static_cast<int>(2 * n), MPI_DOUBLE, MPI_MAX,
                  parent_fes.GetComm());
  }
  bool ok = true;
  for (std::size_t j = 0; j < n; j++)
  {
    for (std::size_t k = 0; modes_ready && k < mode_ids.size(); k++)
    {
      col[j] = (mode_ids[k] == ids[j]) ? static_cast<int>(k) : col[j];
    }
    if (region_only && col[j] >= 0)
    {
      continue;  // the environment, and so its modes, are those of the model
    }
    MFEM_VERIFY(!region_only || reach[n + j] > 0.0,
                "Excitation " << ids[j]
                              << " has no DOFs in the mesh of the region and is not in the "
                                 "saved substructuring model: run on the full mesh!");
    if (reach[j] == 0.0)
    {
      col[j] = -1;  // zero modes
      continue;
    }
    need[j] = 1;
    if (col[j] < 0 || region_only)
    {
      ok = false;
      continue;
    }
    ok = ok &&
         FingerprintsMatch(fps[j], std::span<const double>(mode_fp).subspan(2 * col[j], 2),
                           0, kFingerprintTol);
  }
  return ok;
}

std::vector<double> SubstructuringSolver::Impl::GatherRows(const std::vector<double> &rows,
                                                           int w) const
{
  MPI_Comm comm = parent_fes.GetComm();
  const bool root = (Mpi::Rank(comm) == 0);
  std::vector<double> full;
  if (w == 0)
  {
    return full;
  }
  const std::vector<int> &cnt = iface.counts, &disp = iface.offsets;
  if (root)
  {
    full.assign(static_cast<std::size_t>(iface.size) * w, 0.0);
  }
  const RowType row(w);
  MPI_Gatherv(rows.data(), iface.count, row.type, root ? full.data() : nullptr, cnt.data(),
              disp.data(), row.type, 0, comm);
  return full;
}

void SubstructuringSolver::Impl::ScatterRows(const std::vector<double> &full, int w,
                                             std::vector<double> &rows) const
{
  rows.assign(static_cast<std::size_t>(iface.count) * w, 0.0);
  if (w == 0)
  {
    return;
  }
  const std::vector<int> &cnt = iface.counts, &disp = iface.offsets;
  MPI_Comm comm = parent_fes.GetComm();
  const RowType row(w);
  MPI_Scatterv(Mpi::Rank(comm) == 0 ? full.data() : nullptr, cnt.data(), disp.data(),
               row.type, rows.data(), iface.count, row.type, 0, comm);
}

void SubstructuringSolver::Impl::ApplyAeeInvMulti(
    const std::vector<const mfem::Vector *> &X, const std::vector<mfem::Vector *> &Y) const
{
  MFEM_ASSERT(X.size() == Y.size(), "ApplyAeeInvMulti size mismatch!");
  BlockTimer bt(Timer::KSP);
  MFEM_VERIFY(env_schur, "Environment solver used before EnsureEnv()!");
  // The internal solve: Gamma held at 0, and 0 off the environment's side.
  env_schur->SolveInternal(X, Y);
}

void SubstructuringSolver::Impl::CondenseRows(const std::vector<const mfem::Vector *> &V,
                                              std::vector<double> &rows) const
{
  MFEM_ASSERT(env_condenses, "The environment does not condense onto Gamma!");
  // The condensations (rank 0, nG x |V|, column-major), to the owners of their rows.
  const int nG = iface.size, nv = static_cast<int>(V.size());
  std::vector<double> g, full;
  env_schur->Condense(V, g);
  if (Mpi::Root(parent_fes.GetComm()))
  {
    full.resize(static_cast<std::size_t>(nG) * nv);
    for (int k = 0; k < nv; k++)
    {
      for (int q = 0; q < nG; q++)
      {
        full[static_cast<std::size_t>(q) * nv + k] =
            g[q + static_cast<std::size_t>(nG) * k];
      }
    }
  }
  ScatterRows(full, nv, rows);
}

void SubstructuringSolver::Impl::MaterializeWithEnergy()
{
  const int nG = iface.size, nloc = iface.count, off = iface.offset;
  // Owned interface column -> local true DOF.
  std::vector<int> col_to_dof(nloc, -1);
  for (int i = 0; i < nt; i++)
  {
    if (iface.gamma[i])
    {
      col_to_dof[iface.index[i] - off] = i;
    }
  }
  // S_E from the environment's factorization, on gamma_grid. The energy operator S^K = S_E
  // - P^T D P, for the lift p = e_c - y with y = A_EE^-1 (A_env e_c)|_E, takes the
  // condensation P^T D p onto Gamma, and the kinetic one S^Kin = P^T M_kin^E P that of
  // M_kin^E p, in batches of kMaterializeBlock lifts whose condensations (by rows) are
  // added to the grid; the thin coupling is G = A_env E_Gamma.
  S_grid = env_schur->ReleaseSchurBlocks(*gamma_grid);
  std::unique_ptr<mfem::HypreParMatrix> E = AssembleSelection(iface.gamma);
  std::unique_ptr<mfem::HypreParMatrix> G(mfem::ParMult(&Aenv(), E.get()));
  E.reset();
  const int B = kMaterializeBlock;
  const bool kin = EnvKinetic();
  std::vector<Vector> t(B, Vector(nt)), v(kin ? 2 * B : B, Vector(nt));
  std::vector<const mfem::Vector *> X(B), V;
  std::vector<mfem::Vector *> Y(B);
  Vector ec(nloc);
  SK_grid = S_grid;
  have_sk = true;
  SKin_grid.assign(kin ? gamma_grid->LocalSize() : 0, 0.0);
  have_skin = kin;
  const mfem::HypreParMatrix &Dm = Denv();
  std::vector<double> gr, panel;
  for (int c0 = 0; c0 < nG; c0 += B)
  {
    const int nb = std::min(B, nG - c0);
    X.resize(nb);
    Y.resize(nb);
    V.resize(kin ? 2 * nb : nb);
    for (int k = 0; k < nb; k++)
    {
      const int c = c0 + k;
      ec = 0.0;
      if (c >= off && c < off + nloc)
      {
        ec(c - off) = 1.0;
      }
      // A_env e_c, whose env-interior part is the solve RHS.
      G->Mult(ec, t[k]);
      X[k] = &t[k];
      Y[k] = &t[k];  // solve in place
    }
    ApplyAeeInvMulti(X, Y);
    for (int k = 0; k < nb; k++)
    {
      const int c = c0 + k;
      // The lift p = e_c - y (t[k] holds y), its regularization load D p and its kinetic
      // load M_kin^E p.
      t[k].Neg();
      if (c >= off && c < off + nloc)
      {
        t[k].HostReadWrite()[col_to_dof[c - off]] += 1.0;
      }
      Dm.Mult(t[k], v[k]);
      V[k] = &v[k];
      if (kin)
      {
        KineticEnvMult(t[k], v[B + k]);
        V[nb + k] = &v[B + k];
      }
    }
    CondenseRows(V, gr);
    // The columns of the batch of S^K (less the condensations of D p) and of S^Kin.
    const std::size_t nv = V.size();
    panel.resize(static_cast<std::size_t>(nloc) * nb);
    for (int part = 0; part < (kin ? 2 : 1); part++)
    {
      for (int r = 0; r < nloc; r++)
      {
        for (int k = 0; k < nb; k++)
        {
          const double v = gr[r * nv + part * nb + k];
          panel[static_cast<std::size_t>(r) * nb + k] = (part == 0) ? -v : v;
        }
      }
      gamma_grid->AddFromRows(panel, off, nloc, c0, nb, (part == 0) ? SK_grid : SKin_grid);
    }
  }
}

}  // namespace palace
