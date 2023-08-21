// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "errorindicators.hpp"

#include "utils/communication.hpp"

namespace palace
{

void ErrorReductionOperator::operator()(ErrorIndicators &ebar, const Vector &ind,
                                        double p) const
{
  // Compute the global indicator across all processors.
  auto comm = Mpi::World();
  double candidate_global_error_indicator =
      std::transform_reduce(ind.begin(), ind.end(), 0.0, std::plus(),
                            [&p](auto val) { return std::pow(val, p); });

  Mpi::GlobalSum(1, &candidate_global_error_indicator, comm);

  candidate_global_error_indicator = std::pow(candidate_global_error_indicator, 1.0 / p);

  ebar.global_error_indicator =
      std::max(ebar.global_error_indicator, candidate_global_error_indicator);

  // update the average local indicator. Using running average update rather
  // than sum and final division to maintain validity at all times.
  auto running_average = [this](const auto &xbar, const auto &x)
  { return (xbar * n + x) / (n + 1); };

  if (n > 0)
  {
    MFEM_VERIFY(ebar.local_error_indicators.Size() == ind.Size(),
                "Local error indicator vectors mismatch.");
    // Combine these error indicators into the current average.
    std::transform(ebar.local_error_indicators.begin(), ebar.local_error_indicators.end(),
                   ind.begin(), ebar.local_error_indicators.begin(), running_average);
  }
  else
  {
    // This is the first sample, just steal the data.
    ebar.local_error_indicators = ind;
  }

  // Another sample has been added, increment for the running average lambda.
  n++;
}

}  // namespace palace
