// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "errorindicators.hpp"

#include "utils/communication.hpp"
#include "utils/quadrules.hpp"

namespace palace
{

void ErrorReductionOperator::operator()(ErrorIndicators &ebar, mfem::Vector &&ind) const
{

  // Compute the global indicator across all processors
  auto comm = Mpi::World();
  double candidate_global_error_indicator = std::accumulate(ind.begin(), ind.end(), 0.0);
  Mpi::GlobalSum(1, &candidate_global_error_indicator, comm);

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
    ebar.local_error_indicators = std::move(ind);
  }

  // Another sample has been added, increment for the running average lambda.
  ++n;
}

}  // namespace palace
