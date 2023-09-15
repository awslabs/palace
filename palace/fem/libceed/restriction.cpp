// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "restriction.hpp"

namespace palace::ceed
{

namespace internal
{

std::unordered_map<RestrKey, CeedElemRestriction, RestrHash> restr_map;

}  // namespace internal

namespace
{

void InitLexicoRestr(const mfem::FiniteElementSpace &fespace,
                     const std::vector<int> &indices, bool use_bdr, Ceed ceed,
                     CeedElemRestriction *restr)
{
  const std::size_t ne = indices.size();
  const mfem::FiniteElement &fe =
      use_bdr ? *fespace.GetBE(indices[0]) : *fespace.GetFE(indices[0]);
  const int P = fe.GetDof();
  const mfem::TensorBasisElement *tfe = dynamic_cast<const mfem::TensorBasisElement *>(&fe);
  const mfem::Array<int> &dof_map = tfe->GetDofMap();
  CeedInt compstride =
      (fespace.GetOrdering() == mfem::Ordering::byVDIM) ? 1 : fespace.GetNDofs();
  const int stride = (compstride == 1) ? fespace.GetVDim() : 1;
  mfem::Array<int> tp_el_dof(ne * P), dofs;
  mfem::Array<bool> tp_el_orients(ne * P);
  bool use_el_orients = false;
  mfem::DofTransformation dof_trans;

  for (std::size_t i = 0; i < ne; i++)
  {
    // No need to handle DofTransformation for tensor-product elements.
    const int elem_index = indices[i];
    if (use_bdr)
    {
      fespace.GetBdrElementDofs(elem_index, dofs, dof_trans);
    }
    else
    {
      fespace.GetElementDofs(elem_index, dofs, dof_trans);
    }
    MFEM_VERIFY(!dof_trans.GetDofTransformation(),
                "Unexpected DofTransformation for lexicographic element "
                "restriction.");
    for (int j = 0; j < P; j++)
    {
      const int sdid = dof_map[j];  // signed
      const int did = (sdid >= 0) ? sdid : -1 - sdid;
      const int sgid = dofs[did];  // signed
      const int gid = (sgid >= 0) ? sgid : -1 - sgid;
      tp_el_dof[j + P * i] = stride * gid;
      tp_el_orients[j + P * i] = (sgid >= 0 && sdid < 0) || (sgid < 0 && sdid >= 0);
      use_el_orients = use_el_orients || tp_el_orients[j + P * i];
    }
  }

  if (use_el_orients)
  {
    PalaceCeedCall(ceed, CeedElemRestrictionCreateOriented(
                             ceed, ne, P, fespace.GetVDim(), compstride,
                             fespace.GetVDim() * fespace.GetNDofs(), CEED_MEM_HOST,
                             CEED_COPY_VALUES, tp_el_dof.GetData(), tp_el_orients.GetData(),
                             restr));
  }
  else
  {
    PalaceCeedCall(ceed, CeedElemRestrictionCreate(
                             ceed, ne, P, fespace.GetVDim(), compstride,
                             fespace.GetVDim() * fespace.GetNDofs(), CEED_MEM_HOST,
                             CEED_COPY_VALUES, tp_el_dof.GetData(), restr));
  }
}

void InitNativeRestr(const mfem::FiniteElementSpace &fespace,
                     const std::vector<int> &indices, bool use_bdr, bool has_dof_trans,
                     bool is_interp_range, Ceed ceed, CeedElemRestriction *restr)
{
  const std::size_t ne = indices.size();
  const mfem::FiniteElement &fe =
      use_bdr ? *fespace.GetBE(indices[0]) : *fespace.GetFE(indices[0]);
  const int P = fe.GetDof();
  CeedInt compstride =
      (fespace.GetOrdering() == mfem::Ordering::byVDIM) ? 1 : fespace.GetNDofs();
  const int stride = (compstride == 1) ? fespace.GetVDim() : 1;
  mfem::Array<int> tp_el_dof(ne * P), dofs;
  mfem::Array<bool> tp_el_orients;
  mfem::Array<int8_t> tp_el_curl_orients;
  bool use_el_orients = false;
  mfem::DofTransformation dof_trans;
  mfem::Vector el_trans_j;
  if (!has_dof_trans)
  {
    tp_el_orients.SetSize(ne * P);
  }
  else
  {
    tp_el_curl_orients.SetSize(ne * P * 3, 0);
    el_trans_j.SetSize(P);
  }

  for (std::size_t i = 0; i < ne; i++)
  {
    const auto e = indices[i];
    if (use_bdr)
    {
      fespace.GetBdrElementDofs(e, dofs, dof_trans);
    }
    else
    {
      fespace.GetElementDofs(e, dofs, dof_trans);
    }
    if (!has_dof_trans)
    {
      for (int j = 0; j < P; j++)
      {
        const int sgid = dofs[j];  // signed
        const int gid = (sgid >= 0) ? sgid : -1 - sgid;
        tp_el_dof[j + P * i] = stride * gid;
        tp_el_orients[j + P * i] = (sgid < 0);
        use_el_orients = use_el_orients || tp_el_orients[j + P * i];
      }
    }
    else
    {
      for (int j = 0; j < P; j++)
      {
        const int sgid = dofs[j];  // signed
        const int gid = (sgid >= 0) ? sgid : -1 - sgid;
        tp_el_dof[j + P * i] = stride * gid;

        // Fill column j of element tridiagonal matrix tp_el_curl_orients.
        el_trans_j = 0.0;
        el_trans_j(j) = 1.0;
        if (is_interp_range)
        {
          dof_trans.InvTransformDual(el_trans_j);
        }
        else
        {
          dof_trans.InvTransformPrimal(el_trans_j);
        }
        double sign_j = (sgid < 0) ? -1.0 : 1.0;
        tp_el_curl_orients[3 * (j + 0 + P * i) + 1] =
            static_cast<int8_t>(sign_j * el_trans_j(j + 0));
        if (j > 0)
        {
          tp_el_curl_orients[3 * (j - 1 + P * i) + 2] =
              static_cast<int8_t>(sign_j * el_trans_j(j - 1));
        }
        if (j < P - 1)
        {
          tp_el_curl_orients[3 * (j + 1 + P * i) + 0] =
              static_cast<int8_t>(sign_j * el_trans_j(j + 1));
        }
#ifdef MFEM_DEBUG
        int nnz = 0;
        for (int k = 0; k < P; k++)
        {
          if (k < j - 1 && k > j + 1 && el_trans_j(k) != 0.0)
          {
            nnz++;
          }
        }
        MFEM_ASSERT(nnz == 0, "Element transformation matrix is not tridiagonal at column "
                                  << j << " (nnz = " << nnz << ")!");
#endif
      }
    }
  }

  if (tp_el_curl_orients.Size())
  {
    PalaceCeedCall(ceed, CeedElemRestrictionCreateCurlOriented(
                             ceed, ne, P, fespace.GetVDim(), compstride,
                             fespace.GetVDim() * fespace.GetNDofs(), CEED_MEM_HOST,
                             CEED_COPY_VALUES, tp_el_dof.GetData(),
                             tp_el_curl_orients.GetData(), restr));
  }
  else if (use_el_orients)
  {
    PalaceCeedCall(ceed, CeedElemRestrictionCreateOriented(
                             ceed, ne, P, fespace.GetVDim(), compstride,
                             fespace.GetVDim() * fespace.GetNDofs(), CEED_MEM_HOST,
                             CEED_COPY_VALUES, tp_el_dof.GetData(), tp_el_orients.GetData(),
                             restr));
  }
  else
  {
    PalaceCeedCall(ceed, CeedElemRestrictionCreate(
                             ceed, ne, P, fespace.GetVDim(), compstride,
                             fespace.GetVDim() * fespace.GetNDofs(), CEED_MEM_HOST,
                             CEED_COPY_VALUES, tp_el_dof.GetData(), restr));
  }
}

}  // namespace

void InitRestriction(const mfem::FiniteElementSpace &fespace,
                     const std::vector<int> &indices, bool use_bdr, bool is_interp,
                     bool is_range, Ceed ceed, CeedElemRestriction *restr)
{
  // Check for fespace -> restriction in hash table.
  // The restriction for an interpolator range space is slightly different as
  // the output is a primal vector instead of a dual vector, and lexicographic
  // ordering is never used (no use of tensor-product basis).
  const mfem::FiniteElement &fe =
      use_bdr ? *fespace.GetBE(indices[0]) : *fespace.GetFE(indices[0]);
  const int ncomp = fespace.GetVDim();
  mfem::Array<int> dofs;
  mfem::DofTransformation dof_trans;
  if (use_bdr)
  {
    fespace.GetBdrElementDofs(indices[0], dofs, dof_trans);
  }
  else
  {
    fespace.GetElementDofs(indices[0], dofs, dof_trans);
  }
  const bool has_dof_trans = dof_trans.GetDofTransformation() && !dof_trans.IsEmpty();
  const bool unique_range_restr = (is_interp && is_range && has_dof_trans);
  internal::RestrKey key(ceed, (void *)&fe, ncomp, unique_range_restr);

  // Initialize or retrieve key values.
  auto restr_itr = internal::restr_map.find(key);
  if (restr_itr == internal::restr_map.end())
  {
    const mfem::TensorBasisElement *tfe =
        dynamic_cast<const mfem::TensorBasisElement *>(&fe);
    const bool vector = fe.GetRangeType() == mfem::FiniteElement::VECTOR;
    const bool lexico = (tfe && tfe->GetDofMap().Size() > 0 && !vector && !is_interp);
    if (lexico)
    {
      // Lexicographic ordering using dof_map.
      InitLexicoRestr(fespace, indices, use_bdr, ceed, restr);
    }
    else
    {
      // Native ordering.
      InitNativeRestr(fespace, indices, use_bdr, has_dof_trans, is_interp && is_range, ceed,
                      restr);
    }
    PalacePragmaOmp(critical(InitRestriction))
    {
      internal::restr_map[key] = *restr;
    }
  }
  else
  {
    *restr = restr_itr->second;
  }
}

}  // namespace palace::ceed
