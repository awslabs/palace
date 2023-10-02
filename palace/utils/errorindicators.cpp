// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "errorindicators.hpp"

#include "utils/communication.hpp"

namespace palace
{

ErrorIndicators::ErrorIndicators(const IndicatorsAndNormalization &indicators,
                                 int global_true_v_size, MPI_Comm comm)
  : local_error_indicators(indicators.indicators), global_true_v_size(global_true_v_size),
    comm(comm), normalization(indicators.normalization)
{
  // Compute the global indicator across all processors.
  constexpr int p = 2;
  global_error_indicator =
      std::transform_reduce(local_error_indicators.begin(), local_error_indicators.end(),
                            0.0, std::plus(), [&p](auto val) { return std::pow(val, p); });
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

void ErrorReductionOperator::operator()(ErrorIndicators &ebar,
                                        const IndicatorsAndNormalization &ind,
                                        double p) const
{
  if (n == 0)
  {
    // No direct reduction necessary.
    ebar = ErrorIndicators(ind, ebar.global_true_v_size, ebar.comm);
    n++;
    return;
  }

  const auto comm = ebar.GetComm();

  // Compute the global indicator across all processors.
  double candidate_global_error_indicator =
      std::transform_reduce(ind.indicators.begin(), ind.indicators.end(), 0.0, std::plus(),
                            [&p](auto val) { return std::pow(val, p); });
  Mpi::GlobalSum(1, &candidate_global_error_indicator, comm);
  candidate_global_error_indicator = std::pow(candidate_global_error_indicator, 1.0 / p);

  // Update the global error indicator and local error indicators to the new mean normalized
  // values. The average local indicator is used rather than the indicator for the maximum
  // error to drive the adaptation, to account for a local error that might be marginally
  // important to many solves, rather than only large in one solve.

  // TODO: Could alternatively consider the maximum.
  ebar.global_error_indicator =
      (n * ebar.global_error_indicator + candidate_global_error_indicator) / (n + 1);
  auto running_average = [this](const auto &xbar, const auto &x)
  { return (xbar * n + x) / (n + 1); };

  MFEM_VERIFY(ebar.local_error_indicators.Size() == ind.indicators.Size(),
              "Local error indicator vectors mismatch.");
  // Combine these error indicators into the current average.
  std::transform(ebar.local_error_indicators.begin(), ebar.local_error_indicators.end(),
                 ind.indicators.begin(), ebar.local_error_indicators.begin(),
                 running_average);

  // Assumes that the global error indicator is already reduced across processors.
  ebar.min = ebar.local_error_indicators.Min();
  ebar.max = ebar.local_error_indicators.Max();
  int size = ebar.local_error_indicators.Size();
  Mpi::GlobalMin(1, &ebar.min, comm);
  Mpi::GlobalMax(1, &ebar.max, comm);
  Mpi::GlobalSum(1, &size, comm);
  ebar.mean = ebar.global_error_indicator / size;

  // Another sample has been added, increment for the running average lambda.
  n++;
}

}  // namespace palace
