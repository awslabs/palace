// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_ERROR_INDICATORS_HPP
#define PALACE_UTILS_ERROR_INDICATORS_HPP

#include "linalg/vector.hpp"

namespace palace
{

// Storage for error estimation results from the solve. Required in the AMR loop. An error
// indicator is non-negative, while an error estimate is signed.
class ErrorIndicators
{
public:
  ErrorIndicators(int ndof) : ndof(ndof) {}
  ErrorIndicators() = delete;
  ErrorIndicators(const ErrorIndicators &) = default;
  ErrorIndicators(ErrorIndicators &&) = default;
  ErrorIndicators &operator=(const ErrorIndicators &) = default;
  ErrorIndicators &operator=(ErrorIndicators &&) = default;
  ~ErrorIndicators() = default;

  // Elemental localized error indicators. Used for marking elements for
  // refinement and coarsening.
  Vector local_error_indicators;
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
  // Reduce a vector indicators, i, computed with norm p, into the combined
  // error indicator e.
  void operator()(ErrorIndicators &e, const Vector &i, double p = 2) const;
};

}  // namespace palace

#endif  // PALACE_UTILS_ERROR_INDICATORS_HPP
