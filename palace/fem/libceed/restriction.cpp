// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "restriction.hpp"

#include <mfem.hpp>
#include "utils/omp.hpp"

namespace palace::ceed
{

namespace
{

const mfem::FiniteElement *GetTraceElement(const mfem::FiniteElementSpace &fespace,
                                           const std::vector<int> &indices)
{
  int elem_id, face_info;
  fespace.GetMesh()->GetBdrElementAdjacentElement(indices[0], elem_id, face_info);
  mfem::Geometry::Type face_geom = fespace.GetMesh()->GetBdrElementGeometry(indices[0]);
  return fespace.GetTraceElement(elem_id, face_geom);
};

mfem::Array<int> GetFaceDofsFromAdjacentElement(const mfem::FiniteElementSpace &fespace,
                                                mfem::DofTransformation &dof_trans,
                                                const int P, const int e)
{
  // Get coordinates of face dofs.
  int elem_id, face_info;
  fespace.GetMesh()->GetBdrElementAdjacentElement(e, elem_id, face_info);
  mfem::Geometry::Type face_geom = fespace.GetMesh()->GetBdrElementGeometry(e);
  face_info = fespace.GetMesh()->EncodeFaceInfo(
      fespace.GetMesh()->DecodeFaceInfoLocalIndex(face_info),
      mfem::Geometry::GetInverseOrientation(
          face_geom, fespace.GetMesh()->DecodeFaceInfoOrientation(face_info)));
  mfem::IntegrationPointTransformation Loc1;
  fespace.GetMesh()->GetLocalFaceTransformation(fespace.GetMesh()->GetBdrElementType(e),
                                                fespace.GetMesh()->GetElementType(elem_id),
                                                Loc1.Transf, face_info);
  const mfem::FiniteElement *face_el = fespace.GetTraceElement(elem_id, face_geom);
  MFEM_VERIFY(dynamic_cast<const mfem::NodalFiniteElement *>(face_el),
              "Mesh requires nodal Finite Element.");
  mfem::IntegrationRule face_ir(face_el->GetDof());
  Loc1.Transf.ElementNo = elem_id;
  Loc1.Transf.mesh = fespace.GetMesh();
  Loc1.Transf.ElementType = mfem::ElementTransformation::ELEMENT;
  Loc1.Transform(face_el->GetNodes(), face_ir);
  mfem::DenseMatrix face_pm;
  fespace.GetMesh()->GetNodes()->GetVectorValues(Loc1.Transf, face_ir, face_pm);

  // Get coordinates of element dofs.
  mfem::DenseMatrix elem_pm;
  const mfem::FiniteElement *fe_elem = fespace.GetFE(elem_id);
  mfem::IsoparametricTransformation T;
  fespace.GetMesh()->GetElementTransformation(elem_id, &T);
  T.Transform(fe_elem->GetNodes(), elem_pm);

  // Find the dofs.
  double tol = 1E-5;
  mfem::Array<int> elem_dofs, dofs(P);
  fespace.GetElementDofs(elem_id, elem_dofs, dof_trans);
  for (int l = 0; l < P; l++)
  {
    double norm2_f = 0.0;
    for (int m = 0; m < face_pm.Height(); m++)
    {
      norm2_f += face_pm(m, l) * face_pm(m, l);
    }

    bool found_match = false;
    MFEM_CONTRACT_VAR(found_match);  // silence unused warning
    for (int m = 0; m < elem_pm.Width(); m++)
    {
      double norm2_e = 0.0;
      for (int n = 0; n < elem_pm.Height(); n++)
      {
        norm2_e += elem_pm(n, m) * elem_pm(n, m);
      }
      double relative_tol = tol * std::max(std::max(norm2_f, norm2_e), 1.0E-6);
      double diff = 0.0;
      for (int o = 0; o < elem_pm.Height(); o++)
      {
        diff += std::fabs(elem_pm(o, m) - face_pm(o, l));
      }
      if (diff <= relative_tol)
      {
        dofs[l] = elem_dofs[m];
        found_match = true;
        break;
      }
    }

    MFEM_ASSERT(found_match,
                [&]()
                {
                  std::stringstream msg;
                  msg << "l " << l << '\n';
                  msg << "elem_dofs\n";
                  for (auto x : elem_dofs)
                    msg << x << ' ';

                  msg << "\ndofs\n";
                  for (auto x : dofs)
                    msg << x << ' ';
                  msg << '\n';
                  return msg.str();
                }());
  }

  return dofs;
};

void InitLexicoRestr(const mfem::FiniteElementSpace &fespace,
                     const std::vector<int> &indices, bool use_bdr, Ceed ceed,
                     CeedElemRestriction *restr)
{
  const std::size_t num_elem = indices.size();
  const mfem::FiniteElement *fe;
  bool face_flg = false;
  if (!use_bdr)
  {
    fe = fespace.GetFE(indices[0]);
  }
  else
  {
    fe = fespace.GetBE(indices[0]);
    if (!fe)
    {
      fe = GetTraceElement(fespace, indices);
      face_flg = true;
    }
  }
  const int P = fe->GetDof();
  const mfem::TensorBasisElement *tfe = dynamic_cast<const mfem::TensorBasisElement *>(fe);
  const mfem::Array<int> &dof_map = tfe->GetDofMap();
  const bool dof_map_is_identity = dof_map.Size() == 0;
  const CeedInt comp_stride =
      (fespace.GetVDim() == 1 || fespace.GetOrdering() == mfem::Ordering::byVDIM)
          ? 1
          : fespace.GetNDofs();
  const int stride =
      (fespace.GetOrdering() == mfem::Ordering::byVDIM) ? fespace.GetVDim() : 1;
  mfem::Array<int> tp_el_dof(num_elem * P);
  mfem::Array<bool> tp_el_orients(num_elem * P);
  int use_el_orients = 0;

  PalacePragmaOmp(parallel reduction(+ : use_el_orients))
  {
    mfem::Array<int> dofs;
    mfem::DofTransformation dof_trans;
    bool use_el_orients_loc = false;

    PalacePragmaOmp(for schedule(static))
    for (std::size_t i = 0; i < num_elem; i++)
    {
      // No need to handle DofTransformation for tensor-product elements.
      const int e = indices[i];
      if (use_bdr)
      {
        if (!face_flg)
        {
          fespace.GetBdrElementDofs(e, dofs, dof_trans);
        }
        else
        {
          dofs = GetFaceDofsFromAdjacentElement(fespace, dof_trans, P, e);
        }
      }
      else
      {
        fespace.GetElementDofs(e, dofs, dof_trans);
      }
      MFEM_VERIFY(!dof_trans.GetDofTransformation(),
                  "Unexpected DofTransformation for lexicographic element "
                  "restriction.");
      for (int j = 0; j < P; j++)
      {
        const int sdid = dof_map_is_identity ? j : dof_map[j];  // signed
        const int did = (sdid >= 0) ? sdid : -1 - sdid;
        const int sgid = dofs[did];  // signed
        const int gid = (sgid >= 0) ? sgid : -1 - sgid;
        tp_el_dof[j + P * i] = stride * gid;
        tp_el_orients[j + P * i] = (sgid >= 0 && sdid < 0) || (sgid < 0 && sdid >= 0);
        use_el_orients_loc = use_el_orients_loc || tp_el_orients[j + P * i];
      }
    }
    use_el_orients += use_el_orients_loc;
  }

  if (use_el_orients)
  {
    PalaceCeedCall(ceed, CeedElemRestrictionCreateOriented(
                             ceed, num_elem, P, fespace.GetVDim(), comp_stride,
                             fespace.GetVDim() * fespace.GetNDofs(), CEED_MEM_HOST,
                             CEED_COPY_VALUES, tp_el_dof.GetData(), tp_el_orients.GetData(),
                             restr));
  }
  else
  {
    PalaceCeedCall(ceed, CeedElemRestrictionCreate(
                             ceed, num_elem, P, fespace.GetVDim(), comp_stride,
                             fespace.GetVDim() * fespace.GetNDofs(), CEED_MEM_HOST,
                             CEED_COPY_VALUES, tp_el_dof.GetData(), restr));
  }
}

void InitNativeRestr(const mfem::FiniteElementSpace &fespace,
                     const std::vector<int> &indices, bool use_bdr, bool is_interp_range,
                     Ceed ceed, CeedElemRestriction *restr)
{
  const std::size_t num_elem = indices.size();
  const mfem::FiniteElement *fe;
  bool face_flg = false;
  if (!use_bdr)
  {
    fe = fespace.GetFE(indices[0]);
  }
  else
  {
    fe = fespace.GetBE(indices[0]);
    if (!fe)
    {
      fe = GetTraceElement(fespace, indices);
      face_flg = true;
    }
  }
  const int P = fe->GetDof();
  const CeedInt comp_stride =
      (fespace.GetVDim() == 1 || fespace.GetOrdering() == mfem::Ordering::byVDIM)
          ? 1
          : fespace.GetNDofs();
  const int stride =
      (fespace.GetOrdering() == mfem::Ordering::byVDIM) ? fespace.GetVDim() : 1;
  const bool has_dof_trans = [&]()
  {
    if (fespace.GetMesh()->Dimension() < 3)
    {
      return false;
    }
    const auto geom = fe->GetGeomType();
    const auto *dof_trans = fespace.FEColl()->DofTransformationForGeometry(geom);
    return (dof_trans && !dof_trans->IsIdentity());
  }();
  mfem::Array<int> tp_el_dof(num_elem * P);
  mfem::Array<bool> tp_el_orients;
  mfem::Array<int8_t> tp_el_curl_orients;
  if (!has_dof_trans)
  {
    tp_el_orients.SetSize(num_elem * P);
  }
  else
  {
    tp_el_curl_orients.SetSize(num_elem * P * 3, 0);
  }
  int use_el_orients = 0;

  PalacePragmaOmp(parallel reduction(+ : use_el_orients))
  {
    mfem::Array<int> dofs;
    mfem::DofTransformation dof_trans;
    mfem::Vector el_trans_j;
    if (has_dof_trans)
    {
      el_trans_j.SetSize(P);
      el_trans_j = 0.0;
    }
    bool use_el_orients_loc = false;

    PalacePragmaOmp(for schedule(static))
    for (std::size_t i = 0; i < num_elem; i++)
    {
      const auto e = indices[i];
      if (use_bdr)
      {
        if (!face_flg)
        {
          fespace.GetBdrElementDofs(e, dofs, dof_trans);
        }
        else
        {
          dofs = GetFaceDofsFromAdjacentElement(fespace, dof_trans, P, e);
        }
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
          use_el_orients_loc = use_el_orients_loc || tp_el_orients[j + P * i];
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
              static_cast<int8_t>(sign_j * el_trans_j(j));
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

#if defined(MFEM_DEBUG)
          // Debug check that transformation is actually tridiagonal.
          int nnz = 0;
          for (int k = 0; k < P; k++)
          {
            if ((k < j - 1 || k > j + 1) && el_trans_j(k) != 0.0)
            {
              nnz++;
            }
          }
          MFEM_ASSERT(nnz == 0,
                      "Element transformation matrix is not tridiagonal at column "
                          << j << " (nnz = " << nnz << ")!");
#endif

          // Zero out column vector for next iteration.
          el_trans_j(j) = 0.0;
          if (j > 0)
          {
            el_trans_j(j - 1) = 0.0;
          }
          if (j < P - 1)
          {
            el_trans_j(j + 1) = 0.0;
          }
        }
      }
    }
    use_el_orients += use_el_orients_loc;
  }

  if (has_dof_trans)
  {
    PalaceCeedCall(ceed, CeedElemRestrictionCreateCurlOriented(
                             ceed, num_elem, P, fespace.GetVDim(), comp_stride,
                             fespace.GetVDim() * fespace.GetNDofs(), CEED_MEM_HOST,
                             CEED_COPY_VALUES, tp_el_dof.GetData(),
                             tp_el_curl_orients.GetData(), restr));
  }
  else if (use_el_orients)
  {
    PalaceCeedCall(ceed, CeedElemRestrictionCreateOriented(
                             ceed, num_elem, P, fespace.GetVDim(), comp_stride,
                             fespace.GetVDim() * fespace.GetNDofs(), CEED_MEM_HOST,
                             CEED_COPY_VALUES, tp_el_dof.GetData(), tp_el_orients.GetData(),
                             restr));
  }
  else
  {
    PalaceCeedCall(ceed, CeedElemRestrictionCreate(
                             ceed, num_elem, P, fespace.GetVDim(), comp_stride,
                             fespace.GetVDim() * fespace.GetNDofs(), CEED_MEM_HOST,
                             CEED_COPY_VALUES, tp_el_dof.GetData(), restr));
  }
}

}  // namespace

void InitRestriction(const mfem::FiniteElementSpace &fespace,
                     const std::vector<int> &indices, bool use_bdr, bool is_interp,
                     bool is_interp_range, Ceed ceed, CeedElemRestriction *restr)
{
  MFEM_ASSERT(!indices.empty(), "Empty element index set for libCEED element restriction!");
  if constexpr (false)
  {
    std::cout << "New element restriction (" << ceed << ", " << &fespace << ", "
              << indices[0] << ", " << use_bdr << ", " << is_interp << ", "
              << is_interp_range << ")\n";
  }
  const mfem::FiniteElement *fe;
  if (!use_bdr)
  {
    fe = fespace.GetFE(indices[0]);
  }
  else
  {
    fe = fespace.GetBE(indices[0]);
    if (!fe)
    {
      fe = GetTraceElement(fespace, indices);
    }
  }
  const mfem::TensorBasisElement *tfe = dynamic_cast<const mfem::TensorBasisElement *>(fe);
  const bool vector = fe->GetRangeType() == mfem::FiniteElement::VECTOR;
  const bool lexico = (tfe && !vector && !is_interp);
  if (lexico)
  {
    // Lexicographic ordering using dof_map.
    InitLexicoRestr(fespace, indices, use_bdr, ceed, restr);
  }
  else
  {
    // Native ordering.
    InitNativeRestr(fespace, indices, use_bdr, is_interp_range, ceed, restr);
  }
}

}  // namespace palace::ceed
