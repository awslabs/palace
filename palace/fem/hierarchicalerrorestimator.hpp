// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_HIERARCHICALERRORESTIMATOR_HPP
#define PALACE_FEM_HIERARCHICALERRORESTIMATOR_HPP

#include <set>
#include <vector>
#include <mfem.hpp>

namespace palace::fem::hierarchical
{

// Sparse coefficients of one complete conforming patch direction in a combined standard
// plus singular true-DOF numbering. The owning estimator is responsible for orientation and
// constraint transformations before constructing this column.
struct SparseColumn
{
  std::vector<int> dofs;
  std::vector<double> values;
};

// One uneliminated element or boundary-facet contribution. Several contributions may share
// the same support element (for example domain curl/mass plus impedance facets).
struct LocalOperatorContribution
{
  int support_element = -1;
  std::vector<int> dofs;
  mfem::DenseMatrix matrix;
  mfem::Vector rhs;
};

// Assemble b-Ax without constructing A. Essential equations are set to zero after all local
// contributions are scattered.
mfem::Vector AssembleResidual(int combined_size,
                              const std::vector<LocalOperatorContribution> &contributions,
                              const mfem::Vector &injected,
                              const std::vector<bool> &essential);

// Assemble B^T A B and B^T r over a complete support union. Returns the number of local
// matrix entries touched, a deterministic work proxy.
long long
AssembleRestrictedOperator(const std::vector<LocalOperatorContribution> &contributions,
                           const std::set<int> &support_elements,
                           const std::vector<SparseColumn> &basis,
                           const mfem::Vector &residual, mfem::DenseMatrix &restricted,
                           mfem::Vector &restricted_rhs);

// Evaluate x^T A x by local contributions. This is also used for support-local scalar line
// searches and element indicators.
double Energy(const std::vector<LocalOperatorContribution> &contributions,
              const mfem::Vector &vector);

}  // namespace palace::fem::hierarchical

#endif  // PALACE_FEM_HIERARCHICALERRORESTIMATOR_HPP
