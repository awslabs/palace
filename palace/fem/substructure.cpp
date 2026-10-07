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

namespace
{

bool IsNedelec(const mfem::ParFiniteElementSpace &fespace)
{
  return dynamic_cast<const mfem::ND_FECollection *>(fespace.FEColl()) != nullptr;
}

}  // namespace

int SignatureWidth(const mfem::ParFiniteElementSpace &fespace)
{
  return IsNedelec(fespace) ? 12 : 3;
}

std::vector<double> TrueDofSignatures(const mfem::ParFiniteElementSpace &fespace,
                                      const std::vector<int> &index, int n)
{
  const int w = SignatureWidth(fespace), nt = fespace.GetTrueVSize();
  std::vector<double> loc(static_cast<std::size_t>(n) * w, 0.0),
      glob(static_cast<std::size_t>(n) * w, 0.0);
  auto &fes = const_cast<mfem::ParFiniteElementSpace &>(fespace);
  mfem::ParGridFunction gf(&fes);
  mfem::Vector td(nt);
  auto stamp = [&](int slot)
  {
    gf.GetTrueDofs(td);
    for (int i = 0; i < nt; i++)
    {
      if (index[i] >= 0)
      {
        loc[static_cast<std::size_t>(index[i]) * w + slot] = td(i);
      }
    }
  };
  if (!IsNedelec(fespace))
  {
    for (int d = 0; d < fes.GetParMesh()->SpaceDimension(); d++)
    {
      mfem::FunctionCoefficient xc([d](const mfem::Vector &x) { return x(d); });
      gf.ProjectCoefficient(xc);
      stamp(d);
    }
  }
  else
  {
    for (int b = 0; b < 3; b++)
    {
      mfem::Vector e(3);
      e = 0.0;
      e(b) = 1.0;
      mfem::VectorConstantCoefficient ec(e);
      gf.ProjectCoefficient(ec);
      stamp(b);
    }
    for (int a = 0; a < 3; a++)
    {
      for (int b = 0; b < 3; b++)
      {
        mfem::VectorFunctionCoefficient xc(3,
                                           [a, b](const mfem::Vector &x, mfem::Vector &v)
                                           {
                                             v = 0.0;
                                             v(b) = x(a);
                                           });
        gf.ProjectCoefficient(xc);
        stamp(3 + a * 3 + b);
      }
    }
  }
  MPI_Allreduce(loc.data(), glob.data(), n * w, MPI_DOUBLE, MPI_SUM, fespace.GetComm());
  return glob;
}

void MatchSignatures(const std::vector<double> &current, const std::vector<double> &saved,
                     int width, bool signed_match, std::vector<int> &perm,
                     std::vector<double> &sgn)
{
  const int n = static_cast<int>(current.size() / std::max(width, 1));
  MFEM_VERIFY(current.size() == saved.size(),
              "Saved and current DOF sets of different sizes cannot be matched!");
  perm.assign(n, 0);
  sgn.assign(n, 1.0);
  double worst = 0.0, scale = 1.0;
  for (double v : saved)
  {
    scale = std::max(scale, std::abs(v));
  }
  for (int g = 0; g < n; g++)
  {
    int best = 0;
    double bs = 1.0, bd = 1e300;
    for (int s = 0; s < n; s++)
    {
      double dp = 0.0, dm = 0.0;
      for (int d = 0; d < width; d++)
      {
        const double a = current[static_cast<std::size_t>(g) * width + d];
        const double b = saved[static_cast<std::size_t>(s) * width + d];
        dp += (a - b) * (a - b);
        if (signed_match)
        {
          dm += (a + b) * (a + b);
        }
      }
      if (dp < bd)
      {
        bd = dp;
        best = s;
        bs = 1.0;
      }
      if (signed_match && dm < bd)
      {
        bd = dm;
        best = s;
        bs = -1.0;
      }
    }
    perm[g] = best;
    sgn[g] = bs;
    worst = std::max(worst, bd);
  }
  std::vector<char> used(n, 0);
  bool bijective = true;
  for (int g = 0; g < n; g++)
  {
    bijective = bijective && !used[perm[g]];
    used[perm[g]] = 1;
  }
  MFEM_VERIFY(bijective && std::sqrt(worst) < 1e-8 * scale,
              "Online interface DOFs do not match the saved model (max mismatch "
                  << std::sqrt(worst) << "); the interface Gamma must be identical.");
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

  // The point of a point-tangent DOF u(p).t, whose signature is (t_b, p_a t_b):
  // p = Q t / |t|^2 with Q_ab = p_a t_b.
  auto point = [&](const std::vector<double> &sig, int g)
  {
    const double *v = row(sig, g);
    std::array<double, 3> p = {0.0, 0.0, 0.0};
    const double t2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    for (int a = 0; a < 3 && t2 > 0.0; a++)
    {
      for (int b = 0; b < 3; b++)
      {
        p[a] += v[3 + 3 * a + b] * v[b] / t2;
      }
    }
    return p;
  };
  auto near = [tol](const std::array<double, 3> &p, const std::array<double, 3> &q)
  { return std::hypot(p[0] - q[0], p[1] - q[1], p[2] - q[2]) <= tol; };

  // Signed one-to-one matches.
  SignatureMap map;
  map.rows.assign(n, {});
  std::vector<int> match(n, -1), matched_by(n, -1);
  for (int g = 0; g < n; g++)
  {
    double bd = tol;
    for (int s = 0; s < n; s++)
    {
      for (double sg : {1.0, -1.0})
      {
        double d = 0.0;
        for (int q = 0; q < width; q++)
        {
          const double e = row(current, g)[q] - sg * row(saved, s)[q];
          d += e * e;
        }
        if (std::sqrt(d) <= bd && matched_by[s] < 0)
        {
          bd = std::sqrt(d);
          if (match[g] >= 0)
          {
            matched_by[match[g]] = -1;
          }
          match[g] = s;
          matched_by[s] = g;
          map.rows[g] = {{s, sg}};
        }
      }
    }
  }

  // The others with all DOFs at their point (for second-order Nédélec face DOFs on
  // tetrahedra, the two tangents at the face center, one of which may have matched): the
  // block M_g with current = M_g saved, M_g = C S^T (S S^T)^-1 for the signature blocks.
  std::vector<char> done(n, 0);
  double worst = 0.0;
  bool ok = true;
  for (int g0 = 0; g0 < n && ok; g0++)
  {
    if (match[g0] >= 0 || done[g0])
    {
      continue;
    }
    const auto p = point(current, g0);
    std::vector<int> C, Sv;
    for (int g = 0; g < n; g++)
    {
      if (near(point(current, g), p))
      {
        C.push_back(g);
      }
    }
    for (int s = 0; s < n; s++)
    {
      if (near(point(saved, s), p))
      {
        Sv.push_back(s);
      }
    }
    if (C.size() != Sv.size())
    {
      ok = false;
      break;
    }
    const int k = static_cast<int>(C.size());
    mfem::DenseMatrix Cm(k, width), Sm(k, width), SSt(k, k), CSt(k, k), Mg(k, k);
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
    mfem::DenseMatrix R(k, width);
    mfem::Mult(Mg, Sm, R);
    R -= Cm;
    worst = std::max(worst, R.MaxMaxNorm());
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
  // Every saved DOF is used once.
  std::vector<int> uses(n, 0);
  for (const auto &r : map.rows)
  {
    ok = ok && !r.empty();
    for (const auto &[s, c] : r)
    {
      uses[s]++;
    }
  }
  for (int g = 0; g < n && ok; g++)
  {
    ok = (uses[g] >= 1);
  }
  MFEM_VERIFY(ok && worst <= tol,
              "Online interface DOFs do not match the saved model (max mismatch "
                  << worst << "); the interface Gamma must be identical.");
  return map;
}

}  // namespace palace
