// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_ERROR_INDICATORS_HPP
#define PALACE_ERROR_INDICATORS_HPP

#include <mfem.hpp>

namespace palace
{

// Storage for error estimation results from the solve. Required in the AMR loop. An error
// indicator is non-negative, whilst an error estimate is signed.
struct ErrorIndicators
{
  ErrorIndicators(int ndof) : ndof(ndof) {}
  ErrorIndicators() = delete;
  ErrorIndicators(const ErrorIndicators &) = default;
  ErrorIndicators(ErrorIndicators &&) = default;
  ErrorIndicators &operator=(const ErrorIndicators &) = default;
  ErrorIndicators &operator=(ErrorIndicators &&) = default;
  ~ErrorIndicators() = default;

  // Elemental localized error indicators. Used for marking elements for
  // refinement and coarsening.
  mfem::Vector local_error_indicators;
  // Global error indicator. Used for driving AMR and diagnostics. This
  // combines the local error indicators across all ranks. This number is the
  // same on all ranks.
  double global_error_indicator = 0;
  // Number of global dof in the mesh.
  int ndof;
};

// Operator for performing reduction of a vector of local indicators into a
// global running total for use in adaptation.
class ErrorReductionOperator
{
  // number of samples. mutability required to ensure operation.
  mutable int n = 0;

public:
  void operator()(ErrorIndicators &e, mfem::Vector &&i) const;
};

}  // namespace palace

#endif  // PALACE_ERROR_INDICATORS_HPP
