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

}  // namespace palace
