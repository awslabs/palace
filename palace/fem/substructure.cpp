// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "substructure.hpp"

namespace palace
{

bool BuildSubMeshDofMap(const mfem::ParFiniteElementSpace &sub_fespace,
                        const mfem::ParFiniteElementSpace &parent_fespace,
                        const mfem::Array<int> &parent_element_ids,
                        std::vector<int> &par_dof, std::vector<double> &sign)
{
  // Signed dofs are encoded as raw >= 0 -> (dof = raw, +1) and raw < 0 -> (-1-raw, -1).
  auto decode = [](int raw, int &dof, double &s)
  {
    if (raw >= 0)
    {
      dof = raw;
      s = 1.0;
    }
    else
    {
      dof = -1 - raw;
      s = -1.0;
    }
  };

  par_dof.assign(sub_fespace.GetVSize(), -1);
  sign.assign(sub_fespace.GetVSize(), 0.0);
  std::vector<char> is_set(sub_fespace.GetVSize(), 0);
  mfem::Array<int> sd, pd;
  bool ok = true;
  for (int e = 0; e < sub_fespace.GetNE(); e++)
  {
    sub_fespace.GetElementDofs(e, sd);
    parent_fespace.GetElementDofs(parent_element_ids[e], pd);
    if (sd.Size() != pd.Size())
    {
      return false;
    }
    for (int l = 0; l < sd.Size(); l++)
    {
      int ds, dp;
      double ss, sp;
      decode(sd[l], ds, ss);
      decode(pd[l], dp, sp);
      const double rel = ss * sp;
      if (!is_set[ds])
      {
        par_dof[ds] = dp;
        sign[ds] = rel;
        is_set[ds] = 1;
      }
      else if (par_dof[ds] != dp || sign[ds] != rel)
      {
        // A submesh dof mapping inconsistently to two parent dofs (or with two signs)
        // indicates a nonconforming interface or active DOF transformations.
        ok = false;
      }
    }
  }
  return ok;
}

Substructure::Substructure(mfem::ParFiniteElementSpace &parent_fespace,
                           const mfem::Array<int> &domain_attrs,
                           mfem::FiniteElementCollection &fec)
{
  auto &parent_mesh = *parent_fespace.GetParMesh();
  submesh = std::make_unique<mfem::ParSubMesh>(
      mfem::ParSubMesh::CreateFromDomain(parent_mesh, domain_attrs));
  fespace = std::make_unique<mfem::ParFiniteElementSpace>(submesh.get(), &fec);
  mfem::Array<int> emap = submesh->GetParentElementIDMap();
  map_ok = BuildSubMeshDofMap(*fespace, parent_fespace, emap, par_dof, sign);
}

void MarkParentDofOwnership(const std::vector<const Substructure *> &subs,
                            int n_parent_dofs, std::vector<int> &owner)
{
  owner.assign(n_parent_dofs, 0);
  for (std::size_t k = 0; k < subs.size(); k++)
  {
    const int bit = 1 << k;
    for (int pd : subs[k]->GetParentDof())
    {
      if (pd >= 0)
      {
        owner[pd] |= bit;
      }
    }
  }
}

int BuildInterfaceIndex(const std::vector<int> &owner, std::vector<int> &gamma_index)
{
  auto popcount = [](int m)
  {
    int c = 0;
    while (m)
    {
      c += m & 1;
      m >>= 1;
    }
    return c;
  };
  gamma_index.assign(owner.size(), -1);
  int nG = 0;
  for (std::size_t p = 0; p < owner.size(); p++)
  {
    if (popcount(owner[p]) > 1)
    {
      gamma_index[p] = nG++;
    }
  }
  return nG;
}

DtNBoundaryOperator::DtNBoundaryOperator(const Substructure &environment,
                                         const std::vector<int> &gamma_index,
                                         const mfem::SparseMatrix &A_env,
                                         const mfem::Vector &f_env)
{
  const auto &par = environment.GetParentDof();
  const auto &sgn = environment.GetSign();
  const int n = A_env.Height();

  // Classify environment DOFs: interface (compact gamma index) vs interior (compact index).
  std::vector<int> g_of(n, -1), i_of(n, -1);
  int nG = 0, nI = 0;
  for (int i = 0; i < n; i++)
  {
    if (gamma_index[par[i]] >= 0)
    {
      g_of[i] = nG++;
    }
    else
    {
      i_of[i] = nI++;
    }
  }
  // Map each local interface DOF to its global gamma index and orientation sign.
  std::vector<int> gam_glob(nG);
  std::vector<double> gam_sign(nG);
  for (int i = 0; i < n; i++)
  {
    if (g_of[i] >= 0)
    {
      gam_glob[g_of[i]] = gamma_index[par[i]];
      gam_sign[g_of[i]] = sgn[i];
    }
  }

  // Unsigned local blocks (interior signs cancel in the Schur complement).
  mfem::DenseMatrix Aii(nI), AiG(nI, nG), AGi(nG, nI), AGG(nG);
  Aii = 0.0;
  AiG = 0.0;
  AGi = 0.0;
  AGG = 0.0;
  mfem::Vector fi(nI), fG(nG);
  fi = 0.0;
  fG = 0.0;
  for (int i = 0; i < n; i++)
  {
    const int *cols = A_env.GetRowColumns(i);
    const double *vals = A_env.GetRowEntries(i);
    const bool ig = g_of[i] >= 0;
    for (int k = 0; k < A_env.RowSize(i); k++)
    {
      const int j = cols[k];
      const bool jg = g_of[j] >= 0;
      const double v = vals[k];
      if (ig && jg)
      {
        AGG(g_of[i], g_of[j]) += v;
      }
      else if (ig && !jg)
      {
        AGi(g_of[i], i_of[j]) += v;
      }
      else if (!ig && jg)
      {
        AiG(i_of[i], g_of[j]) += v;
      }
      else
      {
        Aii(i_of[i], i_of[j]) += v;
      }
    }
    if (ig)
    {
      fG(g_of[i]) += f_env(i);
    }
    else
    {
      fi(i_of[i]) += f_env(i);
    }
  }

  // Raw Schur: S = AGG - AGi Aii^-1 AiG ; g = fG - AGi Aii^-1 fi.
  mfem::DenseMatrix Aii_inv(Aii);
  if (nI)
  {
    Aii_inv.Invert();
  }
  mfem::DenseMatrix tmp(nG, nI);
  if (nI)
  {
    mfem::Mult(AGi, Aii_inv, tmp);
  }
  mfem::DenseMatrix rawS(AGG);
  if (nI)
  {
    mfem::DenseMatrix t2(nG, nG);
    mfem::Mult(tmp, AiG, t2);
    rawS -= t2;
  }
  mfem::Vector rawg(fG);
  if (nI)
  {
    mfem::Vector t3(nG);
    tmp.Mult(fi, t3);
    rawg -= t3;
  }

  // Reindex into global gamma order and apply interface orientation signs:
  //   S_E(a,b) = s_a s_b rawS ,  g_E(a) = s_a rawg .
  S_E.SetSize(nG);
  S_E = 0.0;
  g_E.SetSize(nG);
  g_E = 0.0;
  for (int a = 0; a < nG; a++)
  {
    g_E(gam_glob[a]) = gam_sign[a] * rawg(a);
    for (int b = 0; b < nG; b++)
    {
      S_E(gam_glob[a], gam_glob[b]) = gam_sign[a] * gam_sign[b] * rawS(a, b);
    }
  }
}

}  // namespace palace
