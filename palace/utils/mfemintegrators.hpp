// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MFEM_INTEGRATORS_HPP
#define PALACE_MFEM_INTEGRATORS_HPP

#include <mfem.hpp>

#include "quadrules.hpp"

namespace palace
{

// Similar to MFEM's VectorFEBoundaryTangentLFIntegrator for ND spaces, but instead of
// computing (n x f, v), this just computes (f, v). Also eliminates the a and b quadrature
// parameters and uses utils::GetDefaultRule instead.
class VectorFEBoundaryLFIntegrator : public mfem::LinearFormIntegrator
{
private:
  mfem::VectorCoefficient &Q;
  mfem::DenseMatrix vshape;
  mfem::Vector f_loc, f_hat;

public:
  VectorFEBoundaryLFIntegrator(mfem::VectorCoefficient &QG) : Q(QG), f_loc(QG.GetVDim()) {}

  void AssembleRHSElementVect(const mfem::FiniteElement &fe,
                              mfem::ElementTransformation &Tr,
                              mfem::Vector &elvect) override
  {
    const int dof = fe.GetDof();
    const int dim = fe.GetDim();
    const mfem::IntegrationRule *ir =
        (IntRule != nullptr) ? IntRule : utils::GetDefaultRule(fe, Tr);
    vshape.SetSize(dof, dim);
    elvect.SetSize(dof);
    elvect = 0.0;
    f_hat.SetSize(dim);

    for (int i = 0; i < ir->GetNPoints(); i++)
    {
      const mfem::IntegrationPoint &ip = ir->IntPoint(i);
      Tr.SetIntPoint(&ip);
      fe.CalcVShape(ip, vshape);

      Q.Eval(f_loc, Tr, ip);

      Tr.InverseJacobian().Mult(f_loc, f_hat);
      f_hat *= ip.weight * Tr.Weight();
      vshape.AddMult(f_hat, elvect);
    }
  }
};

// Similar to MFEM's BoundaryLFIntegrator for H1 spaces, but eliminates the a and b
// quadrature parameters and uses utils::GetDefaultRule instead.
class BoundaryLFIntegrator : public mfem::LinearFormIntegrator
{
private:
  mfem::Coefficient &Q;
  mfem::Vector shape;

public:
  BoundaryLFIntegrator(mfem::Coefficient &QG) : Q(QG) {}

  void AssembleRHSElementVect(const mfem::FiniteElement &fe,
                              mfem::ElementTransformation &Tr,
                              mfem::Vector &elvect) override
  {
    const int dof = fe.GetDof();
    const mfem::IntegrationRule *ir =
        (IntRule != nullptr) ? IntRule : utils::GetDefaultRule(fe, Tr);
    shape.SetSize(dof);
    elvect.SetSize(dof);
    elvect = 0.0;

    for (int i = 0; i < ir->GetNPoints(); i++)
    {
      const mfem::IntegrationPoint &ip = ir->IntPoint(i);
      Tr.SetIntPoint(&ip);
      fe.CalcShape(ip, shape);

      double val = ip.weight * Tr.Weight() * Q.Eval(Tr, ip);
      add(elvect, val, shape, elvect);
    }
  }
};

using VectorFEDomainLFIntegrator = VectorFEBoundaryLFIntegrator;
using DomainLFIntegrator = BoundaryLFIntegrator;

}  // namespace palace

#endif  // PALACE_MFEM_INTEGRATORS_HPP
