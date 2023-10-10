// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_BILINEARFORM_HPP
#define PALACE_FEM_BILINEARFORM_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/integrator.hpp"
#include "fem/libceed/operator.hpp"

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
  const mfem::ParFiniteElementSpace &trial_fespace, &test_fespace;

  // List of domain and boundary integrators making up the bilinear form.
  std::vector<std::unique_ptr<BilinearFormIntegrator>> domain_integs, boundary_integs;

  // Integration order for quadrature rules is calculated as p_trial + p_test + w + q_extra,
  // where p_test and p_trial are the test and trial space basis function orders and w is
  // the geometry order.
  int q_extra_pk, q_extra_qk;

public:
  BilinearForm(const mfem::ParFiniteElementSpace &trial_fespace,
               const mfem::ParFiniteElementSpace &test_fespace, int q_extra_pk,
               int q_extra_qk)
    : trial_fespace(trial_fespace), test_fespace(test_fespace), q_extra_pk(q_extra_pk),
      q_extra_qk(q_extra_qk)
  {
  }
  BilinearForm(const mfem::ParFiniteElementSpace &trial_fespace,
               const mfem::ParFiniteElementSpace &test_fespace, int q_extra = 0)
    : BilinearForm(trial_fespace, test_fespace, q_extra, q_extra)
  {
  }
  BilinearForm(const mfem::ParFiniteElementSpace &fespace, int q_extra_pk, int q_extra_qk)
    : BilinearForm(fespace, fespace, q_extra_pk, q_extra_qk)
  {
  }
  BilinearForm(const mfem::ParFiniteElementSpace &fespace, int q_extra = 0)
    : BilinearForm(fespace, fespace, q_extra, q_extra)
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

  template <typename T, typename... U>
  void AddDomainIntegrator(U &&...args)
  {
    domain_integs.push_back(std::make_unique<T>(std::forward<U>(args)...));
  }

  template <typename T, typename... U>
  void AddBoundaryIntegrator(U &&...args)
  {
    boundary_integs.push_back(std::make_unique<T>(std::forward<U>(args)...));
  }

  std::unique_ptr<Operator> Assemble(int pa_order_threshold, bool skip_zeros) const
  {
    if (GetMaxElementOrder() >= pa_order_threshold)
    {
      return Assemble();
    }
    else
    {
      return FullAssemble(skip_zeros);
    }
  }

  std::unique_ptr<mfem::SparseMatrix> FullAssemble(bool skip_zeros) const
  {
    return FullAssemble(*Assemble(), skip_zeros);
  }

  std::unique_ptr<ceed::Operator> Assemble() const;

  static std::unique_ptr<mfem::SparseMatrix> FullAssemble(const ceed::Operator &op,
                                                          bool skip_zeros);
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

  const auto &GetTrialSpace() const { return a.GetTrialSpace(); }
  const auto &GetTestSpace() const { return a.GetTestSpace(); }

  template <typename T, typename... U>
  void AddDomainInterpolator(U &&...args)
  {
    a.AddDomainIntegrator<T>(std::forward<U>(args)...);
  }

  std::unique_ptr<Operator> Assemble(int pa_order_threshold, bool skip_zeros) const
  {
    if (a.GetMaxElementOrder() >= pa_order_threshold)
    {
      return Assemble();
    }
    else
    {
      return FullAssemble(skip_zeros);
    }
  }

  std::unique_ptr<mfem::SparseMatrix> FullAssemble(bool skip_zeros) const
  {
    return FullAssemble(*Assemble(), skip_zeros);
  }

  std::unique_ptr<ceed::Operator> Assemble() const;

  static std::unique_ptr<mfem::SparseMatrix> FullAssemble(const ceed::Operator &op,
                                                          bool skip_zeros);
};

}  // namespace palace

#endif  // PALACE_FEM_BILINEARFORM_HPP
