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

//
// This class implements bilinear and mixed bilinear forms based on integrators assembled
// using the libCEED library. Assembly in the form of a partially assembled operator or
// fully assembled sparse matrix is available.
//
class BilinearForm
{
protected:
  // Domain and range finite element spaces.
  const mfem::FiniteElementSpace &trial_fespace, &test_fespace;

  // List of domain and boundary integrators making up the bilinear form.
  std::vector<std::unique_ptr<BilinearFormIntegrator>> domain_integs, boundary_integs;

  // Integration order for quadrature rules, defaults to p_trial + p_test + w, where p_test
  // and p_trial are the test and trial space basis function orders and w is the geometry
  // order.
  mutable int q_order;

public:
  BilinearForm(const mfem::ParFiniteElementSpace &trial_fespace,
               const mfem::ParFiniteElementSpace &test_fespace, int q_order = -1)
    : trial_fespace(trial_fespace), test_fespace(test_fespace), q_order(q_order)
  {
  }
  BilinearForm(const mfem::ParFiniteElementSpace &fespace, int q_order = -1)
    : BilinearForm(fespace, fespace, q_order)
  {
  }

  const auto &GetTrialSpace() const { return trial_fespace; }
  const auto &GetTestSpace() const { return test_fespace; }

  // MFEM's RT_FECollection actually returns order + 1 for GetOrder() for historical
  // reasons.
  auto GetMaxElementOrder() const
  {
    const auto &trial_fec = *trial_fespace.FEColl();
    const auto &test_fec = *test_fespace.FEColl();
    return std::max(
        dynamic_cast<const mfem::RT_FECollection *>(&trial_fec) ? trial_fec.GetOrder() - 1
                                                                : trial_fec.GetOrder(),
        dynamic_cast<const mfem::RT_FECollection *>(&test_fec) ? test_fec.GetOrder() - 1
                                                               : test_fec.GetOrder());
  }

  void AddDomainIntegrator(std::unique_ptr<BilinearFormIntegrator> &&bfi)
  {
    domain_integs.push_back(std::move(bfi));
  }

  void AddBoundaryIntegrator(std::unique_ptr<BilinearFormIntegrator> &&bfi)
  {
    boundary_integs.push_back(std::move(bfi));
  }

  std::unique_ptr<Operator> Assemble(int pa_order_threshold, bool skip_zeros) const
  {
    return (GetMaxElementOrder() >= pa_order_threshold) ? Assemble()
                                                        : FullAssemble(skip_zeros);
  }

  std::unique_ptr<Operator> Assemble() const;

  std::unique_ptr<mfem::SparseMatrix> FullAssemble(bool skip_zeros) const;
};

// Discrete linear operators map primal vectors to primal vectors for interpolation between
// spaces.
class DiscreteLinearOperator
{
private:
  BilinearForm a;

public:
  DiscreteLinearOperator(const mfem::ParFiniteElementSpace &trial_fespace,
                         const mfem::ParFiniteElementSpace &test_fespace)
    : a(trial_fespace, test_fespace)
  {
  }

  void AddDomainInterpolator(std::unique_ptr<DiscreteInterpolator> &&di)
  {
    a.AddDomainIntegrator(std::move(di));
  }

  std::unique_ptr<Operator> Assemble(int pa_order_threshold, bool skip_zeros) const
  {
    return (a.GetMaxElementOrder() >= pa_order_threshold) ? Assemble()
                                                          : FullAssemble(skip_zeros);
  }

  std::unique_ptr<Operator> Assemble() const;

  std::unique_ptr<mfem::SparseMatrix> FullAssemble(bool skip_zeros) const;
};

}  // namespace palace

#endif  // PALACE_FEM_BILINEARFORM_HPP
