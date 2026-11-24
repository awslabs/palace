// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "integrator.hpp"

#include <algorithm>
#include <initializer_list>
#include "fem/libceed/integrator.hpp"
#include "models/materialoperator.hpp"

namespace palace
{

namespace fem
{

int DefaultIntegrationOrder::Get(const mfem::IsoparametricTransformation &T)
{
  MFEM_VERIFY(
      p_trial >= 0,
      "DefaultIntegrationOrder not initialized, call IoData::CheckConfiguration()!");
  return 2 * p_trial + (q_order_jac ? T.OrderW() : 0) +
         (T.GetFE()->Space() == mfem::FunctionSpace::Pk ? q_order_extra_pk
                                                        : q_order_extra_qk);
}

int DefaultIntegrationOrder::Get(const mfem::ElementTransformation &T)
{
  const auto *T_iso = dynamic_cast<const mfem::IsoparametricTransformation *>(&T);
  MFEM_VERIFY(
      T_iso,
      "Unexpected non-isoparametric element transformation to calculate quadrature order!");
  return Get(*T_iso);
}

int DefaultIntegrationOrder::Get(const mfem::Mesh &mesh, mfem::Geometry::Type geom)
{
  MFEM_VERIFY(mesh.GetNodes(), "The mesh has no nodal FE space!");
  const auto *fec = mesh.GetNodalFESpace()->FEColl();
  const mfem::FiniteElement *fe = fec->FiniteElementForGeometry(geom);
  if (!fe)
  {
    // Periodic meshes use a discontinuous nodal space so that identified vertices can
    // retain distinct physical coordinates. Such collections only provide lower-
    // dimensional mesh geometry through their trace finite elements.
    fe = fec->TraceFiniteElementForGeometry(geom);
  }
  MFEM_VERIFY(fe, "Unable to get mesh finite element for geometry " << geom << "!");
  mfem::IsoparametricTransformation T;
  T.SetFE(fe);
  return Get(T);
}

}  // namespace fem

namespace
{

// Union of the (libCEED, 1-based) attributes on which the given coefficients are not
// identically zero. A null coefficient is the identity on every attribute, in which case
// there is no restricted attribute set.
std::optional<std::vector<int>>
NonzeroCoefficientAttributes(std::initializer_list<const MaterialPropertyCoefficient *> Q)
{
  std::vector<int> attr_list;
  for (const auto *Q_i : Q)
  {
    if (!Q_i)
    {
      return std::nullopt;
    }
    Q_i->AddNonzeroAttributes(attr_list);
  }
  std::sort(attr_list.begin(), attr_list.end());
  attr_list.erase(std::unique(attr_list.begin(), attr_list.end()), attr_list.end());
  return attr_list;
}

}  // namespace

// The forms of these integrators are symmetric for symmetric coefficients, and their
// transposes are the same integrators with the transposed coefficients.
#define PALACE_DEFINE_SYMMETRIC_FORM_INTEGRATOR(Integrator)                \
  bool Integrator::IsSymmetric() const                                     \
  {                                                                        \
    return !Q || Q->IsSymmetric();                                         \
  }                                                                        \
  std::unique_ptr<BilinearFormIntegrator> Integrator::GetTranspose() const \
  {                                                                        \
    auto integ = std::make_unique<Integrator>(*this);                      \
    integ->transpose = !transpose;                                         \
    return integ;                                                          \
  }
#define PALACE_DEFINE_SYMMETRIC_FORM_PAIR_INTEGRATOR(Integrator)           \
  bool Integrator::IsSymmetric() const                                     \
  {                                                                        \
    return (!Q || Q->IsSymmetric()) && (!Q_mass || Q_mass->IsSymmetric()); \
  }                                                                        \
  std::unique_ptr<BilinearFormIntegrator> Integrator::GetTranspose() const \
  {                                                                        \
    auto integ = std::make_unique<Integrator>(*this);                      \
    integ->transpose = !transpose;                                         \
    integ->transpose_mass = !transpose_mass;                               \
    return integ;                                                          \
  }

PALACE_DEFINE_SYMMETRIC_FORM_INTEGRATOR(MassIntegrator)
PALACE_DEFINE_SYMMETRIC_FORM_INTEGRATOR(VectorFEMassIntegrator)
PALACE_DEFINE_SYMMETRIC_FORM_INTEGRATOR(RotatedVectorFEMassIntegrator)
PALACE_DEFINE_SYMMETRIC_FORM_INTEGRATOR(DiffusionIntegrator)
PALACE_DEFINE_SYMMETRIC_FORM_INTEGRATOR(CurlCurlIntegrator)
PALACE_DEFINE_SYMMETRIC_FORM_INTEGRATOR(DivDivIntegrator)
PALACE_DEFINE_SYMMETRIC_FORM_PAIR_INTEGRATOR(DiffusionMassIntegrator)
PALACE_DEFINE_SYMMETRIC_FORM_PAIR_INTEGRATOR(CurlCurlMassIntegrator)
PALACE_DEFINE_SYMMETRIC_FORM_PAIR_INTEGRATOR(DivDivMassIntegrator)

#undef PALACE_DEFINE_SYMMETRIC_FORM_INTEGRATOR
#undef PALACE_DEFINE_SYMMETRIC_FORM_PAIR_INTEGRATOR

std::optional<std::vector<int>>
BilinearFormIntegrator::GetNonzeroCoefficientAttributes() const
{
  return NonzeroCoefficientAttributes({Q});
}

std::optional<std::vector<int>>
DiffusionMassIntegrator::GetNonzeroCoefficientAttributes() const
{
  return NonzeroCoefficientAttributes({Q, Q_mass});
}

std::optional<std::vector<int>>
CurlCurlMassIntegrator::GetNonzeroCoefficientAttributes() const
{
  return NonzeroCoefficientAttributes({Q, Q_mass});
}

std::optional<std::vector<int>>
DivDivMassIntegrator::GetNonzeroCoefficientAttributes() const
{
  return NonzeroCoefficientAttributes({Q, Q_mass});
}

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
