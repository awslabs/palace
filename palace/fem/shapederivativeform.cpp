// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "shapederivativeform.hpp"

#include "fem/fespace.hpp"
#include "fem/libceed/ceed.hpp"
#include "fem/mesh.hpp"
#include "utils/omp.hpp"

namespace palace
{

ShapeDerivativeOperator::ShapeDerivativeOperator(
    const FiniteElementSpace &fespace, const mfem::ParFiniteElementSpace &nodal_fespace)
  : fespace(fespace), nodal_fespace(nodal_fespace),
    op(nodal_fespace.GetVSize(), fespace.GetVSize()),
    u_vec(ceed::internal::NumCeeds(), nullptr), w_vec(ceed::internal::NumCeeds(), nullptr),
    u_l(fespace.GetVSize()), w_l(fespace.GetVSize()), g_l(nodal_fespace.GetVSize())
{
  u_l.UseDevice(true);
  w_l.UseDevice(true);
  g_l.UseDevice(true);
  u_l = 0.0;
  w_l = 0.0;
}

ShapeDerivativeOperator::~ShapeDerivativeOperator()
{
  PalacePragmaOmp(parallel if (u_vec.size() > 1))
  {
    const int id = utils::GetThreadNum();
    Ceed ceed = ceed::internal::GetCeedObjects()[id];
    PalaceCeedCall(ceed, CeedVectorDestroy(&u_vec[id]));
    PalaceCeedCall(ceed, CeedVectorDestroy(&w_vec[id]));
  }
}

void ShapeDerivativeOperator::Mult(const Vector &u, const Vector &w, Vector &g) const
{
  fespace.GetProlongationMatrix()->Mult(u, u_l);
  fespace.GetProlongationMatrix()->Mult(w, w_l);
  g_l = 0.0;

  Ceed ceed;
  PalaceCeedCallBackend(CeedOperatorGetCeed(op[0], &ceed));
  const CeedMemType mem = ceed::GetUsableMemType(ceed);
  auto *g_data = g_l.ReadWrite(mem == CEED_MEM_DEVICE);
  PalacePragmaOmp(parallel if (op.Size() > 1))
  {
    const int id = utils::GetThreadNum();
    Ceed ceed;
    PalaceCeedCallBackend(CeedOperatorGetCeed(op[id], &ceed));
    CeedInt num_sub;
    PalaceCeedCall(ceed, CeedOperatorCompositeGetNumSub(op[id], &num_sub));
    if (num_sub > 0)
    {
      // Re-pointing the passive inputs updates their state, so the sub-operators read the
      // new fields.
      ceed::InitCeedVector(u_l, ceed, &u_vec[id], false);
      ceed::InitCeedVector(w_l, ceed, &w_vec[id], false);
      CeedVector g_vec;
      PalaceCeedCall(ceed, CeedVectorCreate(ceed, g_l.Size(), &g_vec));
      PalaceCeedCall(ceed, CeedVectorSetArray(g_vec, mem, CEED_USE_POINTER, g_data));
      PalaceCeedCall(ceed, CeedOperatorApplyAdd(op[id], CEED_VECTOR_NONE, g_vec,
                                                CEED_REQUEST_IMMEDIATE));
      PalaceCeedCall(ceed, CeedVectorTakeArray(g_vec, mem, nullptr));
      PalaceCeedCall(ceed, CeedVectorDestroy(&g_vec));
    }
  }

  // Sum the contributions to the nodes shared between processes.
  g.SetSize(nodal_fespace.GetTrueVSize());
  nodal_fespace.GetProlongationMatrix()->MultTranspose(g_l, g);
}

void GetComplexShapeGradient(const ShapeDerivativeOperator &dR,
                             const ShapeDerivativeOperator &dI, const ComplexVector &x,
                             Vector &gr, Vector &gi)
{
  // With x = xr + i xi, xᵀ (R + i I) x = xrᵀ R xr - xiᵀ R xi - xrᵀ I xi - xiᵀ I xr
  //                                    + i (xrᵀ R xi + xiᵀ R xr + xrᵀ I xr - xiᵀ I xi).
  const auto &xr = x.Real(), &xi = x.Imag();
  const int size = dR.GetNodalSpace().GetTrueVSize();
  for (auto *g : {&gr, &gi})
  {
    g->SetSize(size);
    g->UseDevice(true);
    *g = 0.0;
  }
  Vector t;
  auto Add = [&t](const ShapeDerivativeOperator &op, const Vector &u, const Vector &w,
                  double a, Vector &g)
  {
    op.Mult(u, w, t);
    g.Add(a, t);
  };
  Add(dR, xr, xr, 1.0, gr);
  Add(dR, xi, xi, -1.0, gr);
  Add(dI, xi, xr, -1.0, gr);
  Add(dI, xr, xi, -1.0, gr);
  Add(dR, xi, xr, 1.0, gi);
  Add(dR, xr, xi, 1.0, gi);
  Add(dI, xr, xr, 1.0, gi);
  Add(dI, xi, xi, -1.0, gi);
}

std::unique_ptr<ShapeDerivativeOperator> ShapeDerivativeForm::Assemble() const
{
  const auto &mesh = fespace.GetMesh();
  const auto *nodes = mesh.Get().GetNodes();
  const auto *nodal_fespace =
      nodes ? dynamic_cast<const mfem::ParFiniteElementSpace *>(nodes->FESpace()) : nullptr;
  MFEM_VERIFY(nodal_fespace &&
                  dynamic_cast<const mfem::H1_FECollection *>(nodal_fespace->FEColl()),
              "Shape derivatives require a mesh with continuous nodes (periodic meshes are "
              "not supported)!");
  MFEM_VERIFY(fespace.GetVDim() == 1 &&
                  (dynamic_cast<const mfem::H1_FECollection *>(&fespace.GetFEColl()) ||
                   dynamic_cast<const mfem::ND_FECollection *>(&fespace.GetFEColl())),
              "Shape derivatives are only available for scalar H1 and H(curl) spaces!");

  std::unique_ptr<ShapeDerivativeOperator> sd_op(
      new ShapeDerivativeOperator(fespace, *nodal_fespace));
  PalacePragmaOmp(parallel if (ceed::internal::NumCeeds() > 1))
  {
    const int id = utils::GetThreadNum();
    Ceed ceed = ceed::internal::GetCeedObjects()[id];
    ceed::InitCeedVector(sd_op->u_l, ceed, &sd_op->u_vec[id]);
    ceed::InitCeedVector(sd_op->w_l, ceed, &sd_op->w_vec[id]);
    for (const auto &[geom, data] : mesh.GetCeedGeomFactorData(ceed))
    {
      const int dim = mfem::Geometry::Dimension[geom];
      const auto &integs = (dim != mesh.Dimension()) ? boundary_integs : domain_integs;
      if (dim < mesh.Dimension() - 1 || integs.empty())
      {
        continue;
      }
      CeedElemRestriction restr = fespace.GetCeedElemRestriction(ceed, geom, data);
      CeedBasis basis = fespace.GetCeedBasis(ceed, geom);
      CeedElemRestriction mesh_restr = FiniteElementSpace::BuildCeedElemRestriction(
          *nodal_fespace, ceed, geom, data.indices);
      CeedBasis mesh_basis = FiniteElementSpace::BuildCeedBasis(*nodal_fespace, ceed, geom);
      for (const auto &integ : integs)
      {
        CeedOperator sub_op;
        integ->Assemble(ceed, restr, basis, sd_op->u_vec[id], sd_op->w_vec[id], mesh_restr,
                        mesh_basis, data.geom_data, data.geom_data_restr, &sub_op);
        sd_op->op.AddSubOperator(sub_op);  // Sub-operator owned by ceed::Operator
      }

      // The restriction and basis are owned by the sub-operators.
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&mesh_restr));
      PalaceCeedCall(ceed, CeedBasisDestroy(&mesh_basis));
    }
  }
  sd_op->op.Finalize();
  return sd_op;
}

}  // namespace palace
