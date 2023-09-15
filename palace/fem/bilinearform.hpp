// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_BILINEARFORM_HPP
#define PALACE_FEM_BILINEARFORM_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/integrator.hpp"
#include "linalg/operator.hpp"

namespace palace
{

// XX TODO WIP: NEED TO HANDLE SYMMETRIC FOR LIBCEED OPERATORS (Symmetric bilinear form
// creates
//              a symmetric ceed::Operator with MultTranspose = Mult?)

//
// XX TODO WIP
//
class BilinearForm
{
protected:
  // Domain and range finite element spaces.
  const mfem::FiniteElementSpace &trial_fespace, &test_fespace;

  // Integration order for quadrature rules, defaults to p_trial + p_test + w, where p_test
  // and p_trial are the test and trial space basis function orders and w is the geometry
  // order.
  mutable int q_order;

  // List of domain and boundary integrators making up the bilinear form.
  std::vector<std::unique_ptr<BilinearFormIntegrator>> domain_integs, boundary_integs;

public:
  BilinearForm(const mfem::ParFiniteElementSpace &trial_fespace,
               const mfem::ParFiniteElementSpace &test_fespace, int q_order = -1)
    : trial_fespace(trial_fespace), test_fespace(test_fespace), q_order(q_order)
  {
  }

  void AddDomainIntegrator(std::unique_ptr<BilinearFormIntegrator> &&bfi)
  {
    domain_integs.push_back(std::move(bfi));
  }

  void AddBoundaryIntegrator(std::unique_ptr<BilinearFormIntegrator> &&bfi)
  {
    boundary_integs.push_back(std::move(bfi));
  }

  virtual std::unique_ptr<Operator> Assemble() const;

  virtual std::unique_ptr<mfem::SparseMatrix> FullAssemble(bool skip_zeros) const;
};

class DiscreteLinearOperator : public BilinearForm
{
public:
  using BilinearForm::BilinearForm;

  std::unique_ptr<Operator> Assemble() const override;

  std::unique_ptr<mfem::SparseMatrix> FullAssemble(bool skip_zeros) const override;
};

}  // namespace palace

#endif  // PALACE_FEM_BILINEARFORM_HPP
