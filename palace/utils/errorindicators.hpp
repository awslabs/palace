// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_ERROR_INDICATORS_HPP
#define PALACE_UTILS_ERROR_INDICATORS_HPP

#include "linalg/vector.hpp"
#include "utils/communication.hpp"

namespace palace
{
class ErrorReductionOperator;

// Unnormalized error indicators and the normalization factor.
struct IndicatorsAndNormalization
{
  Vector indicators;
  double normalization;
};

// Storage for error estimation results from the solve. Required in the AMR loop. This is
// richer than the IndicatorsAndNormalization because it stores derived quantities, and a
// communicator for use in adaptation.
class ErrorIndicators
{
public:
  // Construct an Error indicator from an initial set of indicators
  explicit ErrorIndicators(const IndicatorsAndNormalization &indicators,
                           int global_true_v_size, MPI_Comm comm);
  // Construct an empty ErrorIndicators.
  explicit ErrorIndicators(int global_true_v_size, MPI_Comm comm)
    : global_true_v_size(global_true_v_size), comm(comm)
  {
  }
  ErrorIndicators() = delete;
  ErrorIndicators(const ErrorIndicators &) = default;
  ErrorIndicators(ErrorIndicators &&) = default;
  ErrorIndicators &operator=(const ErrorIndicators &) = default;
  ErrorIndicators &operator=(ErrorIndicators &&) = default;
  ~ErrorIndicators() = default;

  // Return the average normalized local error indicators.
  const auto &GetLocalErrorIndicators() const { return local_error_indicators; }
  // Return the global error indicator.
  auto GetGlobalErrorIndicator() const { return global_error_indicator; }
  // Return the largest normalized local error indicator.
  auto GetMaxErrorIndicator() const { return max; }
  // Return the smallest normalized local error indicator.
  auto GetMinErrorIndicator() const { return min; }
  // Return the mean normalized local error indicator.
  auto GetMeanErrorIndicator() const { return mean; }
  // Return the normalization constant for the absolute error.
  auto GetNormalization() const { return normalization; }
  // The communicator used in any reductions over processors.
  const MPI_Comm &GetComm() { return comm; }
  // Return the global number of true dofs associated with this set of error indicators.
  auto GlobalTrueVSize() const { return global_true_v_size; }

protected:
  friend class ErrorReductionOperator;
  // Elemental localized error indicators. Used for marking elements for
  // refinement and coarsening.
  Vector local_error_indicators;
  // Global error indicator. Used for driving AMR and diagnostics.
  double global_error_indicator = 0;
  // Global number of true dof in the finite element solution for which this indicator is
  // calculated.
  int global_true_v_size;
  // Communicator used in calculation of the global error indicator.
  MPI_Comm comm;
  // Statistics, updated simultaneously with the global error indicator.
  double min, max, mean;
  // Mean normalization constant.
  double normalization;
};

// Operator for performing reduction of a vector of local indicators into a
// global running total for use in adaptation.
class ErrorReductionOperator
{
  // Number of samples. Mutability required to guarantee operation.
  mutable int n = 0;

public:
  // Reduce a vector indicators, i, computed with norm p, into the combined
  // error indicator e.
  void operator()(ErrorIndicators &e, const IndicatorsAndNormalization &i,
                  double p = 2) const;
  // Resets the internal counter for number of samples.
  void Reset() { n = 0; }
};

}  // namespace palace

#endif  // PALACE_UTILS_ERROR_INDICATORS_HPP
