// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "substructure.hpp"

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

}  // namespace palace
