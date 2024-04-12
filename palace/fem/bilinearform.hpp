// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_BILINEARFORM_HPP
#define PALACE_FEM_BILINEARFORM_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/integrator.hpp"
#include "fem/libceed/operator.hpp"
#include "linalg/hypre.hpp"

namespace palace
{

class FiniteElementSpace;
class FiniteElementSpaceHierarchy;

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

  std::unique_ptr<ceed::Operator>
  PartialAssemble(const FiniteElementSpace &trial_fespace,
                  const FiniteElementSpace &test_fespace) const;

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

  void AssembleQuadratureData();

  std::unique_ptr<ceed::Operator> PartialAssemble() const
  {
    return PartialAssemble(GetTrialSpace(), GetTestSpace());
  }

  std::unique_ptr<hypre::HypreCSRMatrix> FullAssemble(bool skip_zeros) const
  {
    return FullAssemble(*PartialAssemble(), skip_zeros, false);
  }

  static std::unique_ptr<hypre::HypreCSRMatrix> FullAssemble(const ceed::Operator &op,
                                                             bool skip_zeros)
  {
    return FullAssemble(op, skip_zeros, false);
  }

  static std::unique_ptr<hypre::HypreCSRMatrix> FullAssemble(const ceed::Operator &op,
                                                             bool skip_zeros, bool set);

  std::unique_ptr<Operator> Assemble(bool skip_zeros) const;

  std::vector<std::unique_ptr<Operator>>
  Assemble(const FiniteElementSpaceHierarchy &fespaces, bool skip_zeros,
           std::size_t l0 = 0) const;
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

  std::unique_ptr<hypre::HypreCSRMatrix> FullAssemble(bool skip_zeros) const
  {
    return BilinearForm::FullAssemble(*PartialAssemble(), skip_zeros, true);
  }

  static std::unique_ptr<hypre::HypreCSRMatrix> FullAssemble(const ceed::Operator &op,
                                                             bool skip_zeros)
  {
    return BilinearForm::FullAssemble(op, skip_zeros, true);
  }
};

}  // namespace palace

#endif  // PALACE_FEM_BILINEARFORM_HPP
