// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_ERROR_INDICATORS_HPP
#define PALACE_UTILS_ERROR_INDICATORS_HPP

#include "linalg/vector.hpp"
#include "utils/communication.hpp"

namespace palace
{
// Storage for error estimation results from the solve. Required in the AMR loop. This is
// richer than the IndicatorsAndNormalization because it stores derived quantities, and a
// communicator for use in adaptation.
class ErrorIndicators
{
public:
  ErrorIndicators(Vector local, double normalization)
    : local(std::move(local)), normalization(normalization), n(1)
  {
  }

  ErrorIndicators() = default;
  ErrorIndicators(const ErrorIndicators &) = default;
  ErrorIndicators(ErrorIndicators &&) = default;
  ErrorIndicators &operator=(const ErrorIndicators &) = default;
  ErrorIndicators &operator=(ErrorIndicators &&) = default;
  ~ErrorIndicators() = default;

  // Return the local error indicators.
  const auto &GetLocalErrorIndicators() const { return local; }
  // Return the global error indicator.
  inline auto GetGlobalErrorIndicator(MPI_Comm comm) const
  {
    constexpr int p = 2;
    double global_error_indicator =
        std::transform_reduce(local.begin(), local.end(), 0.0, std::plus(),
                              [](auto val) { return std::pow(val, p); });
    Mpi::GlobalSum(1, &global_error_indicator, comm);
    return std::pow(global_error_indicator, 1.0 / p);
  }
  // Return the largest local error indicator.
  inline auto GetMaxErrorIndicator(MPI_Comm comm) const
  {
    double max = local.Max();
    Mpi::GlobalMax(1, &max, comm);
    return max;
  }
  // Return the smallest local error indicator.
  inline auto GetMinErrorIndicator(MPI_Comm comm) const
  {
    double min = local.Min();
    Mpi::GlobalMin(1, &min, comm);
    return min;
  }
  // Return the mean local error indicator.
  inline auto GetMeanErrorIndicator(MPI_Comm comm) const
  {
    int size = local.Size();
    auto global_error_indicator = GetGlobalErrorIndicator(comm);
    return global_error_indicator / size;
  }
  // Return the normalization constant for the absolute error.
  inline auto GetNormalization() const { return normalization; }
  // Add a set of indicators to the running totals.
  void AddIndicators(const ErrorIndicators &indicators);
  // Reset a running total of error indicators ready for computing a new running average.
  inline void Reset()
  {
    n = 0;
    local = Vector();
    normalization = 0;
  }

protected:
  // Elemental localized error indicators. Used for marking elements for
  // refinement and coarsening.
  Vector local;
  // Normalization constant.
  double normalization;
  // Number of samples. Mutability required to guarantee operation.
  int n = 0;
};

}  // namespace palace

#endif  // PALACE_UTILS_ERROR_INDICATORS_HPP
