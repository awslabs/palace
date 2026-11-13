// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "substructuringsolverimpl.hpp"

#include <cmath>
#include "utils/timer.hpp"

namespace palace
{

void SubstructuringSolver::Impl::ComputeSourceModes(const std::vector<int> &ids,
                                                    const std::vector<Source> &src)
{
  EnsureEnv();
  MPI_Comm comm = parent_fes.GetComm();
  const int K = static_cast<int>(ids.size());
  src_ids = ids;
  src_fp.assign(2 * K, 0.0);
  std::vector<mfem::Vector> b(K), W(K, mfem::Vector(nt)), q(K, mfem::Vector(nt)),
      z(K, mfem::Vector(nt));
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
      if (is_gamma[i])
      {
        out[static_cast<std::size_t>(gamma_global[i] - gamma_off) * K + k] = vv[i];
      }
    }
  };
  SG_rows.assign(static_cast<std::size_t>(gamma_nloc) * K, 0.0);
  for (int k = 0; k < K; k++)
  {
    Aenv().Mult(W[k], t);
    t.Neg();
    t += b[k];  // g_k = b^E_k|_Gamma - (A_env W_k)|_Gamma
    rows(SG_rows, k, t);
    // q_k = (A_env - D_env) W_k - b^E_k, then z_k = A_EE^-1 q_k|_E.
    Aenv().Mult(W[k], q[k]);
    Denv().Mult(W[k], t);
    q[k] -= t;
    q[k] -= b[k];
    z[k] = 0.0;
    X[k] = &q[k];
    Y[k] = &z[k];
  }
  if (K > 0)
  {
    ApplyAeeInvMulti(X, Y);
  }
  SH_rows.assign(static_cast<std::size_t>(gamma_nloc) * K, 0.0);
  for (int k = 0; k < K; k++)
  {
    Aenv().Mult(z[k], t);
    t.Neg();
    t += q[k];  // h_k = q_k|_Gamma - (A_env z_k)|_Gamma
    rows(SH_rows, k, t);
  }
  SC.assign(static_cast<std::size_t>(K) * K, 0.0);
  std::vector<mfem::Vector> KW(K, mfem::Vector(nt)), MD(K, mfem::Vector(nt)),
      D(K, mfem::Vector(nt));
  std::vector<mfem::Vector> JE(K, mfem::Vector(nt));
  for (int k = 0; k < K; k++)
  {
    Kenv().Mult(W[k], KW[k]);
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
    JE[k] = 0.0;
    if (src[k].HasJ())
    {
      const double *J = src[k].J->HostRead();
      double *je = JE[k].HostWrite();
      for (int i = 0; i < nt; i++)
      {
        je[i] = is_env_int[i] ? J[i] : 0.0;
      }
    }
  }
  SB.assign(static_cast<std::size_t>(K) * K, 0.0);
  for (int i = 0; i < K; i++)
  {
    for (int j = 0; j < K; j++)
    {
      const std::size_t ij = static_cast<std::size_t>(i) * K + j;
      SC[ij] = (W[i] * KW[j]) + (D[i] * MD[j]) - (JE[i] * W[j]) - (JE[j] * W[i]);
      SB[ij] = b[i] * W[j];
    }
  }
  if (K > 0)
  {
    MPI_Allreduce(MPI_IN_PLACE, SC.data(), K * K, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(MPI_IN_PLACE, SB.data(), K * K, MPI_DOUBLE, MPI_SUM, comm);
  }
  src_ready = true;
  src_saved = false;
}

bool SubstructuringSolver::Impl::SourceModesMatch(const std::vector<int> &ids,
                                                  const std::vector<Source> &src,
                                                  std::vector<int> &col,
                                                  std::vector<char> &need)
{
  col.assign(ids.size(), -1);
  need.assign(ids.size(), 0);
  bool ok = true;
  for (std::size_t j = 0; j < ids.size(); j++)
  {
    std::array<double, 3> fp;
    (void)EnvSource(src[j], fp);  // collective: always evaluated
    if (fp[0] == 0.0)
    {
      continue;
    }
    need[j] = 1;
    for (std::size_t k = 0; src_ready && k < src_ids.size(); k++)
    {
      if (src_ids[k] == ids[j])
      {
        col[j] = static_cast<int>(k);
      }
    }
    if (col[j] < 0)
    {
      ok = false;
      continue;
    }
    if (region_only)
    {
      continue;  // the environment, and so its sources, are those of the model
    }
    // Relative to the larger pairing (either may vanish, up to rounding).
    const double scale =
        std::max({std::abs(fp[1]), std::abs(fp[2]), std::abs(src_fp[2 * col[j]]),
                  std::abs(src_fp[2 * col[j] + 1])});
    for (int q = 0; q < 2; q++)
    {
      if (std::abs(fp[1 + q] - src_fp[2 * col[j] + q]) > kFingerprintTol * scale + 1e-300)
      {
        ok = false;
      }
    }
  }
  return ok;
}

void SubstructuringSolver::Impl::MaterializePartial()
{
  BlockTimer bt(Timer::KSP_SETUP);
  mfem::Array<int> other;
  for (int i = 0; i < nt; i++)
  {
    if (!is_env_int[i] && !is_gamma[i])
    {
      other.Append(i);
    }
  }
  env_schur = substructuring::MakePartialFactorization(
      iodata, *EliminateToIdentity(Aenv(), other), InterfaceTrueDofs(), env_partial);
  env_built = (env_schur != nullptr);
}

void SubstructuringSolver::Impl::EnsureEnv()
{
  MFEM_VERIFY(!region_only,
              "Substructuring on the mesh of the region cannot solve the environment (an "
              "excitation not in the saved model, or general Dirichlet data): run on the "
              "full mesh!");
  if (!env_built)
  {
    BlockTimer bt(Timer::KSP_SETUP);
    FactorEnvironmentInterior();
    env_schur = std::make_unique<BackSolveSchurSolver>(
        Aenv(), InterfaceTrueDofs(),
        [this](const std::vector<const mfem::Vector *> &X,
               const std::vector<mfem::Vector *> &Y)
        {
          env_lu_parent->Mult(X, Y);  // the factor of [A_EE, 0; 0, I]
          MaskToEnvInterior(Y);
        },
        kMaterializeBlock);
    env_built = true;
  }
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
      if (is_env_int[i])
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
      touches[k] = (is_env_int[i] && tt[i] != 0.0);
    }
  }
  if (K > 0)
  {
    MPI_Allreduce(MPI_IN_PLACE, touches.data(), K, MPI_INT, MPI_MAX, comm);
  }
  auto solve_touching = [&](std::vector<mfem::Vector> &in, std::vector<mfem::Vector> &out)
  {
    std::vector<const mfem::Vector *> X;
    std::vector<mfem::Vector *> Y;
    for (int k = 0; k < K; k++)
    {
      out[k] = 0.0;
      if (touches[k])
      {
        X.push_back(&in[k]);
        Y.push_back(&out[k]);
      }
    }
    if (!X.empty())
    {
      ApplyAeeInvMulti(X, Y);
    }
  };
  solve_touching(t, w);
  // g = (a - A_env b)|_Gamma, C_kl = x_k . a_l - t_k . b_l  for (a, b) = (t, w) [solve]
  // or (v, z) [energy correction].
  mfem::Vector tmp(nt);
  auto couplings = [&](const std::vector<mfem::Vector> &a,
                       const std::vector<mfem::Vector> &b, std::vector<double> &Grows,
                       std::vector<double> &C)
  {
    Grows.assign(static_cast<std::size_t>(gamma_nloc) * K, 0.0);
    for (int k = 0; k < K; k++)
    {
      Aenv().Mult(b[k], tmp);
      const double *ak = a[k].HostRead(), *tt = tmp.HostRead();
      for (int i = 0; i < nt; i++)
      {
        if (is_gamma[i])
        {
          Grows[static_cast<std::size_t>(gamma_global[i] - gamma_off) * K + k] =
              ak[i] - tt[i];
        }
      }
    }
    C.assign(static_cast<std::size_t>(K) * K, 0.0);
    for (int k = 0; k < K; k++)
    {
      for (int l = 0; l < K; l++)
      {
        C[static_cast<std::size_t>(k) * K + l] = (x[k] * a[l]) - (t[k] * b[l]);
      }
    }
    if (K > 0)
    {
      MPI_Allreduce(MPI_IN_PLACE, C.data(), K * K, MPI_DOUBLE, MPI_SUM, comm);
    }
  };
  couplings(t, w, G_rows, Cmode);
  GK_rows.clear();
  have_gk = EnergyDiffers();
  if (have_gk)
  {
    // P^T D P blocks: p_k = x_k - w_k (the A-harmonic lift), v_k = D p_k,
    // z_k = A_EE^-1 v_k|_E; subtract them from the A-blocks.
    std::vector<mfem::Vector> v(K, mfem::Vector(nt)), z(K, mfem::Vector(nt));
    for (int k = 0; k < K; k++)
    {
      tmp = x[k];
      tmp -= w[k];
      Denv().Mult(tmp, v[k]);
      touches[k] = 1;  // v_k is generally nonzero in E
    }
    solve_touching(v, z);
    std::vector<double> GD, CD;
    couplings(v, z, GD, CD);
    GK_rows = G_rows;
    for (std::size_t q = 0; q < GK_rows.size(); q++)
    {
      GK_rows[q] -= GD[q];
    }
    for (std::size_t q = 0; q < Cmode.size(); q++)
    {
      Cmode[q] -= CD[q];
    }
  }
  modes_ready = true;
  modes_saved = false;
}

bool SubstructuringSolver::Impl::ModesMatch(const std::vector<int> &ids,
                                            const std::vector<mfem::Vector> &xd,
                                            std::vector<int> &col) const
{
  col.assign(ids.size(), -1);
  bool ok = modes_ready;
  if (region_only)
  {
    // The environment is that of the model: the modes are matched by id.
    for (std::size_t j = 0; j < ids.size(); j++)
    {
      for (std::size_t k = 0; k < mode_ids.size(); k++)
      {
        col[j] = (mode_ids[k] == ids[j]) ? static_cast<int>(k) : col[j];
      }
      ok = ok && col[j] >= 0;
    }
    return ok;
  }
  constexpr std::size_t kBlock = 4;  // bounds the memory of the full vectors
  std::vector<std::array<double, 2>> fps(xd.size());
  for (std::size_t j0 = 0; j0 < xd.size(); j0 += kBlock)
  {
    std::vector<mfem::Vector> x;
    std::vector<const mfem::Vector *> X;
    for (std::size_t j = j0; j < std::min(xd.size(), j0 + kBlock); j++)
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
      fps[j] = LiftFingerprint(x[j - j0], t[j - j0]);
    }
  }
  for (std::size_t j = 0; j < ids.size(); j++)
  {
    for (std::size_t k = 0; k < mode_ids.size(); k++)
    {
      if (mode_ids[k] == ids[j])
      {
        col[j] = static_cast<int>(k);
      }
    }
    const auto &fp = fps[j];
    if (col[j] < 0)
    {
      ok = false;
      continue;
    }
    for (int q = 0; q < 2; q++)
    {
      const double ref = mode_fp[2 * col[j] + q];
      if (std::abs(fp[q] - ref) >
          kFingerprintTol * std::max(std::abs(ref), std::abs(fp[q])) + 1e-300)
      {
        ok = false;
      }
    }
  }
  return ok;
}

void SubstructuringSolver::Impl::RowLayout(int w, std::vector<int> &cnt,
                                           std::vector<int> &disp) const
{
  MPI_Comm comm = parent_fes.GetComm();
  const int nranks = Mpi::Size(comm);
  cnt.assign(nranks, 0);
  disp.assign(nranks, 0);
  const int mine = gamma_nloc * w;
  MPI_Allgather(&mine, 1, MPI_INT, cnt.data(), 1, MPI_INT, comm);
  for (int r = 1; r < nranks; r++)
  {
    disp[r] = disp[r - 1] + cnt[r - 1];
  }
}

std::vector<double> SubstructuringSolver::Impl::GatherRows(const std::vector<double> &rows,
                                                           int w) const
{
  std::vector<int> cnt, disp;
  RowLayout(w, cnt, disp);
  MPI_Comm comm = parent_fes.GetComm();
  const bool root = (Mpi::Rank(comm) == 0);
  std::vector<double> full;
  if (root)
  {
    full.assign(static_cast<std::size_t>(nG_global) * w, 0.0);
  }
  MPI_Gatherv(rows.data(), gamma_nloc * w, MPI_DOUBLE, root ? full.data() : nullptr,
              cnt.data(), disp.data(), MPI_DOUBLE, 0, comm);
  return full;
}

void SubstructuringSolver::Impl::ScatterRows(const std::vector<double> &full, int w,
                                             std::vector<double> &rows) const
{
  std::vector<int> cnt, disp;
  RowLayout(w, cnt, disp);
  MPI_Comm comm = parent_fes.GetComm();
  rows.assign(static_cast<std::size_t>(gamma_nloc) * w, 0.0);
  MPI_Scatterv(Mpi::Rank(comm) == 0 ? full.data() : nullptr, cnt.data(), disp.data(),
               MPI_DOUBLE, rows.data(), gamma_nloc * w, MPI_DOUBLE, 0, comm);
}

void SubstructuringSolver::Impl::FactorEnvironmentInterior()
{
  non_env_int.SetSize(0);
  for (int i = 0; i < nt; i++)
  {
    if (!is_env_int[i])
    {
      non_env_int.Append(i);
    }
  }
  env_A_ee_parent = EliminateToIdentity(Aenv(), non_env_int);
  env_lu_parent = std::make_unique<substructuring::BlockedDirectSolver>(
      substructuring::MakeDirectSolver(iodata, parent_fes.GetComm()));
  env_lu_parent->SetOperator(*env_A_ee_parent);
}

void SubstructuringSolver::Impl::ApplyAeeInvMulti(
    const std::vector<const mfem::Vector *> &X, const std::vector<mfem::Vector *> &Y) const
{
  MFEM_ASSERT(X.size() == Y.size(), "ApplyAeeInvMulti size mismatch!");
  BlockTimer bt(Timer::KSP);
  MFEM_VERIFY(env_built, "Environment solver used before EnsureEnv()!");
  // The internal solve (Gamma held at 0), masked to E.
  env_schur->SolveInternal(X, Y);
  MaskToEnvInterior(Y);
}

void SubstructuringSolver::Impl::MaskToEnvInterior(
    const std::vector<mfem::Vector *> &Y) const
{
  for (auto *y : Y)
  {
    double *yy = y->HostReadWrite();
    for (int i = 0; i < nt; i++)
    {
      if (!is_env_int[i])
      {
        yy[i] = 0.0;
      }
    }
  }
}

void SubstructuringSolver::Impl::ApplyAeeInv(const mfem::Vector &x, mfem::Vector &y) const
{
  ApplyAeeInvMulti({&x}, {&y});
}

void SubstructuringSolver::Impl::MaterializeWithEnergy()
{
  const int nG = nG_global, nloc = gamma_nloc, off = gamma_off;
  // Owned interface column -> local true DOF.
  std::vector<int> col_to_dof(nloc, -1);
  for (int i = 0; i < nt; i++)
  {
    if (is_gamma[i])
    {
      col_to_dof[gamma_global[i] - off] = i;
    }
  }
  // The first solve of a direct factor fixes its block size at kMaterializeBlock; the
  // thin couplings are G = A_env E_Gamma and R = E_Gamma^T A_env. Each rank stores its
  // own rows.
  std::unique_ptr<mfem::HypreParMatrix> E = AssembleSelection(is_gamma);
  std::unique_ptr<mfem::HypreParMatrix> G(mfem::ParMult(&Aenv(), E.get()));
  std::unique_ptr<mfem::HypreParMatrix> Et(E->Transpose());
  std::unique_ptr<mfem::HypreParMatrix> R(mfem::ParMult(Et.get(), &Aenv()));
  E.reset();
  Et.reset();
  const int B = kMaterializeBlock;
  std::vector<Vector> t(B, Vector(nt));
  std::vector<const mfem::Vector *> X(B);
  std::vector<mfem::Vector *> Y(B);
  Vector ec(nloc), z(nloc);
  std::vector<Vector> v(B, Vector(nt));
  SK_rows.assign(static_cast<std::size_t>(nloc) * nG, 0.0);
  have_sk = true;
  const mfem::HypreParMatrix &Dm = Denv();
  std::vector<double> agg(static_cast<std::size_t>(B) * nloc);  // A_GG e_c, own rows
  for (int c0 = 0; c0 < nG; c0 += B)
  {
    const int nb = std::min(B, nG - c0);
    X.resize(nb);
    Y.resize(nb);
    for (int k = 0; k < nb; k++)
    {
      const int c = c0 + k;
      ec = 0.0;
      if (c >= off && c < off + nloc)
      {
        ec(c - off) = 1.0;
      }
      // A_env e_c: its env-interior part is the solve RHS, its Gamma part is A_GG e_c.
      G->Mult(ec, t[k]);
      const double *tk = t[k].HostRead();
      for (int r = 0; r < nloc; r++)
      {
        agg[static_cast<std::size_t>(k) * nloc + r] = tk[col_to_dof[r]];
      }
      X[k] = &t[k];
      Y[k] = &t[k];  // solve in place
    }
    ApplyAeeInvMulti(X, Y);
    for (int k = 0; k < nb; k++)
    {
      R->Mult(t[k], z);  // (A_env y)|_Gamma = A_GE A_EE^-1 A_EG e_c, this rank's rows
      const double *zz = z.HostRead();
      const int c = c0 + k;
      for (int r = 0; r < nloc; r++)
      {
        S_rows[static_cast<std::size_t>(r) * nG + c] =
            agg[static_cast<std::size_t>(k) * nloc + r] - zz[r];
      }
    }
    // Energy correction P^T D P e_c: p = e_c - y (the A-harmonic lift of e_c; t[k]
    // holds y), v = D p, then S^K col = S col - (v|_Gamma - R A_EE^-1 v|_E).
    for (int k = 0; k < nb; k++)
    {
      const int c = c0 + k;
      t[k].Neg();
      if (c >= off && c < off + nloc)
      {
        t[k].HostReadWrite()[col_to_dof[c - off]] += 1.0;
      }
      Dm.Mult(t[k], v[k]);
      const double *vk = v[k].HostRead();
      for (int r = 0; r < nloc; r++)
      {
        agg[static_cast<std::size_t>(k) * nloc + r] = vk[col_to_dof[r]];
      }
      X[k] = &v[k];
      Y[k] = &v[k];
    }
    ApplyAeeInvMulti(X, Y);
    for (int k = 0; k < nb; k++)
    {
      R->Mult(v[k], z);
      const double *zz = z.HostRead();
      const int c = c0 + k;
      for (int r = 0; r < nloc; r++)
      {
        const std::size_t q = static_cast<std::size_t>(r) * nG + c;
        SK_rows[q] = S_rows[q] - (agg[static_cast<std::size_t>(k) * nloc + r] - zz[r]);
      }
    }
  }
}

}  // namespace palace
