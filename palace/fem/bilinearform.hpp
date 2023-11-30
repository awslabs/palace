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

  std::unique_ptr<ceed::Operator> PartialAssemble() const;

  std::unique_ptr<mfem::SparseMatrix> FullAssemble(bool skip_zeros) const
  {
    return FullAssemble(*PartialAssemble(), skip_zeros, false);
  }

  std::unique_ptr<mfem::SparseMatrix> FullAssemble(const ceed::Operator &op,
                                                   bool skip_zeros) const
  {
    return FullAssemble(op, skip_zeros, false);
  }

  static std::unique_ptr<mfem::SparseMatrix> FullAssemble(const ceed::Operator &op,
                                                          bool skip_zeros, bool set);

  std::unique_ptr<Operator> Assemble(bool skip_zeros) const
  {
    return Assemble(*this, skip_zeros);
  }

  template <typename T>
  static std::unique_ptr<Operator> Assemble(const T &a, bool skip_zeros)
  {
    // Returns order such that the miniumum for all element types is 1. MFEM's
    // RT_FECollection actually already returns order + 1 for GetOrder() for historical
    // reasons.
    const auto &trial_fec = a.GetTrialSpace().GetFEColl();
    const auto &test_fec = a.GetTestSpace().GetFEColl();
    int max_order = std::max(
        dynamic_cast<const mfem::L2_FECollection *>(&trial_fec) ? trial_fec.GetOrder() + 1
                                                                : trial_fec.GetOrder(),
        dynamic_cast<const mfem::L2_FECollection *>(&test_fec) ? test_fec.GetOrder() + 1
                                                               : test_fec.GetOrder());
    if (max_order >= pa_order_threshold)
    {
      return a.PartialAssemble();
    }
    else
    {
      return a.FullAssemble(skip_zeros);
    }
  }
};

// Discrete linear operators map primal vectors to primal vectors for interpolation between
// spaces.
class DiscreteLinearOperator
{
private:
  // Domain and range finite element spaces.
  const FiniteElementSpace &trial_fespace, &test_fespace;

  // List of domain interpolators making up the discrete linear operator.
  std::vector<std::unique_ptr<DiscreteInterpolator>> domain_interps;

public:
  DiscreteLinearOperator(const FiniteElementSpace &trial_fespace,
                         const FiniteElementSpace &test_fespace)
    : trial_fespace(trial_fespace), test_fespace(test_fespace)
  {
  }

  const auto &GetTrialSpace() const { return trial_fespace; }
  const auto &GetTestSpace() const { return test_fespace; }

  template <typename T, typename... U>
  void AddDomainInterpolator(U &&...args)
  {
    domain_interps.push_back(std::make_unique<T>(std::forward<U>(args)...));
  }

  std::unique_ptr<ceed::Operator> PartialAssemble() const;

  std::unique_ptr<mfem::SparseMatrix> FullAssemble(bool skip_zeros) const
  {
    return BilinearForm::FullAssemble(*PartialAssemble(), skip_zeros, true);
  }

  std::unique_ptr<mfem::SparseMatrix> FullAssemble(const ceed::Operator &op,
                                                   bool skip_zeros) const
  {
    return BilinearForm::FullAssemble(op, skip_zeros, true);
  }

  std::unique_ptr<Operator> Assemble(bool skip_zeros) const
  {
    return BilinearForm::Assemble(*this, skip_zeros);
  }
};

}  // namespace palace

#endif  // PALACE_FEM_BILINEARFORM_HPP
