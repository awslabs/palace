// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "integrator.hpp"

#include "fem/libceed/integrator.hpp"
#include "utils/communication.hpp"
#include "utils/omp.hpp"

namespace palace
{

namespace fem
{

int DefaultIntegrationOrder::Get(const mfem::IsoparametricTransformation &T)
{
  return 2 * p_trial + (q_order_jac ? T.OrderW() : 0) +
         (T.GetFE()->Space() == mfem::FunctionSpace::Pk ? q_order_extra_pk
                                                        : q_order_extra_qk);
}

int DefaultIntegrationOrder::Get(const mfem::ElementTransformation &T)
{
#if defined(MFEM_DEBUG)
  const auto *T_iso = dynamic_cast<const mfem::IsoparametricTransformation *>(&T);
  MFEM_VERIFY(
      T_iso,
      "Unexpected non-isoparametric element transformation to calculate quadrature order!");
  return Get(*T_iso);
#else
  return Get(static_cast<const mfem::IsoparametricTransformation &>(T));
#endif
}

int DefaultIntegrationOrder::Get(const mfem::Mesh &mesh, mfem::Geometry::Type geom)
{
  MFEM_VERIFY(mesh.GetNodes(), "The mesh has no nodal FE space!");
  mfem::IsoparametricTransformation T;
  T.SetFE(mesh.GetNodalFESpace()->FEColl()->FiniteElementForGeometry(geom));
  return Get(T);
}

double IntegrateFunction(
    const mfem::ParMesh &mesh, const mfem::Array<int> &marker, bool bdr,
    mfem::Coefficient &Q,
    std::function<int(const mfem::ElementTransformation &)> GetQuadratureOrder)
{
  double sum = IntegrateFunctionLocal(mesh, marker, bdr, Q, GetQuadratureOrder);
  Mpi::GlobalSum(1, &sum, mesh.GetComm());
  return sum;
}

double IntegrateFunctionLocal(
    const mfem::ParMesh &mesh, const mfem::Array<int> &marker, bool bdr,
    mfem::Coefficient &Q,
    std::function<int(const mfem::ElementTransformation &)> GetQuadratureOrder)
{
  auto ElementIntegral = [&Q, &GetQuadratureOrder](mfem::ElementTransformation &T)
  {
    double sum = 0.0;
    const mfem::IntegrationRule &ir =
        mfem::IntRules.Get(T.GetGeometryType(), GetQuadratureOrder(T));
    for (int j = 0; j < ir.GetNPoints(); j++)
    {
      const mfem::IntegrationPoint &ip = ir.IntPoint(j);
      T.SetIntPoint(&ip);
      sum += Q.Eval(T, ip) * ip.weight * T.Weight();
    }
    return sum;
  };
  double sum = 0.0;
  PalacePragmaOmp(parallel reduction(+ : sum))
  {
    mfem::IsoparametricTransformation T;
    if (bdr)
    {
      PalacePragmaOmp(for schedule(static))
      for (int i = 0; i < mesh.GetNBE(); i++)
      {
        if (!marker[mesh.GetBdrAttribute(i) - 1])
        {
          continue;
        }
        mesh.GetBdrElementTransformation(i, &T);
        sum += ElementIntegral(T);
      }
    }
    else
    {
      PalacePragmaOmp(for schedule(static))
      for (int i = 0; i < mesh.GetNE(); i++)
      {
        if (!marker[mesh.GetAttribute(i) - 1])
        {
          continue;
        }
        mesh.GetElementTransformation(i, &T);
        sum += ElementIntegral(T);
      }
    }
  }
  return sum;
}

}  // namespace fem

void DiscreteInterpolator::Assemble(Ceed ceed, CeedElemRestriction trial_restr,
                                    CeedElemRestriction test_restr, CeedBasis interp_basis,
                                    CeedOperator *op, CeedOperator *op_t)
{
  // Interpolators do not use an integration rule to map between the test and trial spaces.
  ceed::AssembleCeedInterpolator(ceed, trial_restr, test_restr, interp_basis, op, op_t);
}

void VectorFEBoundaryLFIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &fe,
                                                          mfem::ElementTransformation &T,
                                                          mfem::Vector &elvect)
{
  const int dof = fe.GetDof();
  const int dim = fe.GetDim();
  const int q_order = fem::DefaultIntegrationOrder::Get(T);
  const mfem::IntegrationRule &ir = mfem::IntRules.Get(fe.GetGeomType(), q_order);
  f_hat.SetSize(dim);
  vshape.SetSize(dof, dim);
  elvect.SetSize(dof);
  elvect = 0.0;

  for (int i = 0; i < ir.GetNPoints(); i++)
  {
    const mfem::IntegrationPoint &ip = ir.IntPoint(i);
    T.SetIntPoint(&ip);
    fe.CalcVShape(ip, vshape);

    Q.Eval(f_loc, T, ip);
    T.InverseJacobian().Mult(f_loc, f_hat);
    f_hat *= ip.weight * T.Weight();
    vshape.AddMult(f_hat, elvect);
  }
}

void BoundaryLFIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &fe,
                                                  mfem::ElementTransformation &T,
                                                  mfem::Vector &elvect)
{
  const int dof = fe.GetDof();
  const int q_order = fem::DefaultIntegrationOrder::Get(T);
  const mfem::IntegrationRule &ir = mfem::IntRules.Get(fe.GetGeomType(), q_order);
  shape.SetSize(dof);
  elvect.SetSize(dof);
  elvect = 0.0;

  for (int i = 0; i < ir.GetNPoints(); i++)
  {
    const mfem::IntegrationPoint &ip = ir.IntPoint(i);
    T.SetIntPoint(&ip);
    fe.CalcShape(ip, shape);

    double val = ip.weight * T.Weight() * Q.Eval(T, ip);
    elvect.Add(val, shape);
  }
}

}  // namespace palace
