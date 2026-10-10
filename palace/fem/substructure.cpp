// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "substructure.hpp"

#include <algorithm>
#include <array>
#include <cmath>

namespace palace
{

int MarkInterfaceTrueDofs(mfem::ParFiniteElementSpace &parent_fespace,
                          const mfem::Array<int> &region_attrs,
                          const mfem::Array<int> &environment_attrs,
                          mfem::Array<int> &region_marker, mfem::Array<int> &env_marker,
                          mfem::Array<int> &interface_marker)
{
  auto in = [](const mfem::Array<int> &s, int a)
  {
    for (int x : s)
    {
      if (x == a)
      {
        return true;
      }
    }
    return false;
  };

  // L-vector attribute markers over element DOFs.
  mfem::Vector reg_l(parent_fespace.GetVSize()), env_l(parent_fespace.GetVSize());
  reg_l = 0.0;
  env_l = 0.0;
  mfem::Array<int> dofs;
  auto *mesh = parent_fespace.GetParMesh();
  for (int e = 0; e < mesh->GetNE(); e++)
  {
    const int attr = mesh->GetAttribute(e);
    const bool is_reg = in(region_attrs, attr), is_env = in(environment_attrs, attr);
    if (!is_reg && !is_env)
    {
      continue;
    }
    parent_fespace.GetElementDofs(e, dofs);
    for (int i = 0; i < dofs.Size(); i++)
    {
      const int d = dofs[i] >= 0 ? dofs[i] : -1 - dofs[i];
      if (is_reg)
      {
        reg_l(d) = 1.0;
      }
      else
      {
        env_l(d) = 1.0;
      }
    }
  }

  // Reduce to true DOFs with cross-rank accumulation via |P|^T.
  const mfem::HypreParMatrix *P = parent_fespace.Dof_TrueDof_Matrix();
  const int nt = parent_fespace.GetTrueVSize();
  mfem::Vector reg_t(nt), env_t(nt);
  P->AbsMultTranspose(1.0, reg_l, 0.0, reg_t);
  P->AbsMultTranspose(1.0, env_l, 0.0, env_t);

  region_marker.SetSize(nt);
  env_marker.SetSize(nt);
  interface_marker.SetSize(nt);
  int local_iface = 0;
  for (int i = 0; i < nt; i++)
  {
    region_marker[i] = (reg_t(i) > 0.5) ? 1 : 0;
    env_marker[i] = (env_t(i) > 0.5) ? 1 : 0;
    interface_marker[i] = (region_marker[i] && env_marker[i]) ? 1 : 0;
    local_iface += interface_marker[i];
  }
  int global_iface = 0;
  MPI_Allreduce(&local_iface, &global_iface, 1, MPI_INT, MPI_SUM, parent_fespace.GetComm());
  return global_iface;
}

SignatureMap MatchSignatureBasis(const std::vector<double> &current,
                                 const std::vector<double> &saved, int width)
{
  MFEM_VERIFY(width == 12, "MatchSignatureBasis expects Nédélec signatures!");
  MFEM_VERIFY(current.size() == saved.size(),
              "Saved and current DOF sets of different sizes cannot be matched!");
  const int n = static_cast<int>(current.size() / width);
  double scale = 1.0;
  for (double v : saved)
  {
    scale = std::max(scale, std::abs(v));
  }
  const double tol = 1.0e-8 * scale;
  auto row = [width](const std::vector<double> &sig, int g) { return &sig[g * width]; };

  // Signed one-to-one matches.
  SignatureMap map;
  map.rows.assign(n, {});
  std::vector<char> used(n, 0), matched(n, 0);
  for (int g = 0; g < n; g++)
  {
    const double *c = row(current, g);
    int best = -1;
    double bd = tol * tol, bs = 1.0;
    for (int s = 0; s < n; s++)
    {
      const double *v = row(saved, s);
      double dp = 0.0, dm = 0.0;
      for (int q = 0; q < width; q++)
      {
        dp += (c[q] - v[q]) * (c[q] - v[q]);
        dm += (c[q] + v[q]) * (c[q] + v[q]);
      }
      if (!used[s] && std::min(dp, dm) <= bd)
      {
        best = s;
        bd = std::min(dp, dm);
        bs = (dp <= dm) ? 1.0 : -1.0;
      }
    }
    if (best >= 0)
    {
      used[best] = matched[g] = 1;
      map.rows[g] = {{best, bs}};
    }
  }

  // The others with all DOFs at their point (for second-order Nédélec face DOFs on
  // tetrahedra, the two tangents at the face center, one of which may have matched): the
  // block M_g with current = M_g saved, M_g = C S^T (S S^T)^-1 for the signature blocks.
  // A point-tangent DOF u(p).t has the signature (t_b, p_a t_b), so p = Q t / |t|^2 with
  // Q_ab = p_a t_b.
  bool ok = true;
  double worst = 0.0;  // largest block residual
  if (std::find(matched.begin(), matched.end(), 0) != matched.end())
  {
    auto points = [&](const std::vector<double> &sig)
    {
      std::vector<std::array<double, 3>> p(n, {0.0, 0.0, 0.0});
      for (int g = 0; g < n; g++)
      {
        const double *v = row(sig, g);
        const double t2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        for (int a = 0; a < 3 && t2 > 0.0; a++)
        {
          for (int b = 0; b < 3; b++)
          {
            p[g][a] += v[3 + 3 * a + b] * v[b] / t2;
          }
        }
      }
      return p;
    };
    const auto pc = points(current), ps = points(saved);
    auto near = [tol](const std::array<double, 3> &p, const std::array<double, 3> &q)
    { return std::hypot(p[0] - q[0], p[1] - q[1], p[2] - q[2]) <= tol; };
    std::vector<char> done(n, 0);
    for (int g0 = 0; g0 < n && ok; g0++)
    {
      if (matched[g0] || done[g0])
      {
        continue;
      }
      std::vector<int> C, Sv;
      for (int g = 0; g < n; g++)
      {
        if (near(pc[g], pc[g0]))
        {
          C.push_back(g);
        }
        if (near(ps[g], pc[g0]))
        {
          Sv.push_back(g);
        }
      }
      const int k = static_cast<int>(C.size());
      if (k != static_cast<int>(Sv.size()))
      {
        ok = false;
        break;
      }
      mfem::DenseMatrix Cm(k, width), Sm(k, width), SSt(k, k), CSt(k, k), Mg(k, k),
          R(k, width);
      for (int i = 0; i < k; i++)
      {
        for (int q = 0; q < width; q++)
        {
          Cm(i, q) = row(current, C[i])[q];
          Sm(i, q) = row(saved, Sv[i])[q];
        }
      }
      mfem::MultABt(Sm, Sm, SSt);
      mfem::MultABt(Cm, Sm, CSt);
      SSt.Invert();
      mfem::Mult(CSt, SSt, Mg);
      mfem::Mult(Mg, Sm, R);
      R -= Cm;
      worst = std::max(worst, R.MaxMaxNorm());
      ok = (R.MaxMaxNorm() <= tol);  // also false for a singular block (NaN)
      // Row i of M_g^-T: column i of M_g^-1.
      mfem::DenseMatrix Mi(Mg);
      Mi.Invert();
      for (int i = 0; i < k; i++)
      {
        map.rows[C[i]].clear();
        for (int j = 0; j < k; j++)
        {
          map.rows[C[i]].push_back({Sv[j], Mi(j, i)});
        }
        done[C[i]] = 1;
      }
    }
  }
  // A bijection: every current DOF mapped, every saved DOF used.
  std::fill(used.begin(), used.end(), 0);
  for (const auto &r : map.rows)
  {
    ok = ok && !r.empty();
    for (const auto &[s, c] : r)
    {
      used[s] = 1;
    }
  }
  ok = ok && std::find(used.begin(), used.end(), 0) == used.end();
  MFEM_VERIFY(ok, "Online interface DOFs do not match the saved model (max mismatch "
                      << worst << "); the interface Gamma must be identical.");
  return map;
}

}  // namespace palace
