// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_BILINEARFORM_HPP
#define PALACE_FEM_BILINEARFORM_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/fespace.hpp"
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
  const FiniteElementSpace &trial_fespace, &test_fespace;

  // List of domain and boundary integrators making up the bilinear form.
  std::vector<std::unique_ptr<BilinearFormIntegrator>> domain_integs, boundary_integs;

public:
  // Order above which to use partial assembly vs. full.
  inline static int pa_order_threshold = 1;

public:
  BilinearForm(const FiniteElementSpace &trial_fespace,
               const FiniteElementSpace &test_fespace)
    : trial_fespace(trial_fespace), test_fespace(test_fespace)
  {
  }
  BilinearForm(const FiniteElementSpace &fespace) : BilinearForm(fespace, fespace) {}

  const auto &GetTrialSpace() const { return trial_fespace; }
  const auto &GetTestSpace() const { return test_fespace; }

  // Returns order such that the miniumum for all element types is 1. MFEM's RT_FECollection
  // actually already returns order + 1 for GetOrder() for historical reasons.
  auto GetMaxElementOrder() const
  {
    const auto &trial_fec = *trial_fespace.FEColl();
    const auto &test_fec = *test_fespace.FEColl();
    return std::max(
        dynamic_cast<const mfem::L2_FECollection *>(&trial_fec) ? trial_fec.GetOrder() + 1
                                                                : trial_fec.GetOrder(),
        dynamic_cast<const mfem::L2_FECollection *>(&test_fec) ? test_fec.GetOrder() + 1
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

  std::unique_ptr<Operator> Assemble(bool skip_zeros) const
  {
    if (GetMaxElementOrder() >= pa_order_threshold)
    {
      return PartialAssemble();
    }
    else
    {
      return FullAssemble(skip_zeros);
    }
  }

  std::unique_ptr<ceed::Operator> PartialAssemble() const;

  std::unique_ptr<mfem::SparseMatrix> FullAssemble(bool skip_zeros) const
  {
    return FullAssemble(*PartialAssemble(), skip_zeros);
  }

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
  DiscreteLinearOperator(const FiniteElementSpace &trial_fespace,
                         const FiniteElementSpace &test_fespace)
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

  std::unique_ptr<Operator> Assemble(bool skip_zeros) const
  {
    if (a.GetMaxElementOrder() >= a.pa_order_threshold)
    {
      return PartialAssemble();
    }
    else
    {
      return FullAssemble(skip_zeros);
    }
  }

  std::unique_ptr<ceed::Operator> PartialAssemble() const;

  std::unique_ptr<mfem::SparseMatrix> FullAssemble(bool skip_zeros) const
  {
    return FullAssemble(*a.PartialAssemble(), skip_zeros);
  }

  static std::unique_ptr<mfem::SparseMatrix> FullAssemble(const ceed::Operator &op,
                                                          bool skip_zeros);
};

}  // namespace palace

#endif  // PALACE_FEM_BILINEARFORM_HPP
