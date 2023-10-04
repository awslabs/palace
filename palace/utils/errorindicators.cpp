// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "errorindicators.hpp"

#include "utils/communication.hpp"

namespace palace
{

ErrorIndicators::ErrorIndicators(const IndicatorsAndNormalization &indicators,
                                 int global_true_v_size, MPI_Comm comm)
  : local_error_indicators(indicators.indicators), global_true_v_size(global_true_v_size),
    comm(comm), mean_normalization(indicators.normalization)
{
  // Compute the global indicator across all processors.
  constexpr int p = 2;
  global_error_indicator =
      std::transform_reduce(local_error_indicators.begin(), local_error_indicators.end(),
                            0.0, std::plus(), [](auto val) { return std::pow(val, p); });
  Mpi::GlobalSum(1, &global_error_indicator, comm);
  global_error_indicator = std::pow(global_error_indicator, 1.0 / p);
  min = local_error_indicators.Min();
  max = local_error_indicators.Max();
  int size = local_error_indicators.Size();
  Mpi::GlobalMin(1, &min, comm);
  Mpi::GlobalMax(1, &max, comm);
  Mpi::GlobalSum(1, &size, comm);
  mean = global_error_indicator / size;
}

void ErrorIndicators::Reset()
{
  n = 0;
  local_error_indicators = Vector();
  global_error_indicator = std::numeric_limits<double>::max();
  mean_normalization = 0.0;
}

void ErrorIndicators::AddEstimates(const Vector &indicators, double normalization)
{
  // Compute the global indicator across all processors.
  constexpr int p = 2;
  double candidate_global_error_indicator =
      std::transform_reduce(indicators.begin(), indicators.end(), 0.0, std::plus(),
                            [](auto val) { return std::pow(val, p); });
  Mpi::GlobalSum(1, &candidate_global_error_indicator, comm);
  candidate_global_error_indicator = std::pow(candidate_global_error_indicator, 1.0 / p);

  // Update the global error indicator and local error indicators to the new mean normalized
  // values. The average local indicator is used rather than the indicator for the maximum
  // error to drive the adaptation, to account for a local error that might be marginally
  // important to many solves, rather than only large in one solve.
  if (n > 0)
  {
    global_error_indicator =
        (n * global_error_indicator + candidate_global_error_indicator) / (n + 1);
    mean_normalization = (n * mean_normalization + normalization) / (n + 1);
    auto running_average = [this](const auto &xbar, const auto &x)
    { return (xbar * n + x) / (n + 1); };
    MFEM_VERIFY(local_error_indicators.Size() == indicators.Size(),
                "Local error indicator vectors mismatch.");
    // Combine these error indicators into the current average.
    std::transform(local_error_indicators.begin(), local_error_indicators.end(),
                   indicators.begin(), local_error_indicators.begin(), running_average);
  }
  else
  {
    global_error_indicator = candidate_global_error_indicator;
    local_error_indicators = indicators;
    mean_normalization = normalization;
  }

  // Assumes that the global error indicator is already reduced across processors.
  min = local_error_indicators.Min();
  max = local_error_indicators.Max();
  int size = local_error_indicators.Size();
  Mpi::GlobalMin(1, &min, comm);
  Mpi::GlobalMax(1, &max, comm);
  Mpi::GlobalSum(1, &size, comm);
  mean = global_error_indicator / size;

  // Another sample has been added, increment for the running average lambda.
  n++;
}

}  // namespace palace
