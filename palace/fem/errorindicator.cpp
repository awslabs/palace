// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "errorindicator.hpp"
#include <fmt/os.h>
#include "utils/communication.hpp"

namespace palace
{
void ErrorIndicator::AddIndicator(const ErrorIndicator &error_indicators)
{
  // The average local indicator is used rather than the indicator for the maximum
  // error to drive the adaptation, to account for a local error that might be marginally
  // important to many solves, rather than only large in one solve.
  if (n > 0)
  {
    normalization =
        (n * normalization + error_indicators.normalization * error_indicators.n) /
        (n + error_indicators.n);

    // The local indicators must be squared before combining, so that the global error
    // calculation is valid:
    //    E = sqrt(1/N sum_N sum_K eta_{kn}^2)
    // from which it follows that:
    //    E^2 = 1/N sum_N sum_K eta_{kn}^2
    //        = 1/N sum_N E_n
    // Namely the average of the global error indicators included in the reduction.
    // Squaring both sides means the summation can be rearranged, and then the local error
    // indicators become:
    //    e_K = sqrt(1/N sum_N eta_{Kn}^2)
    auto running_average = [this, error_indicators](const auto &xbar, const auto &x)
    {
      return std::sqrt((xbar * xbar * n + x * x * error_indicators.n) /
                       (n + error_indicators.n));
    };
    MFEM_VERIFY(local.Size() == error_indicators.local.Size(),
                "Local error indicator vectors mismatch.");
    // Combine these error indicators into the current average.
    std::transform(local.begin(), local.end(), error_indicators.local.begin(),
                   local.begin(), running_average);

    // More samples have been added, update for the running average lambda.
    n += error_indicators.n;
  }
  else
  {
    // This indicator was empty, just steal.
    (*this) = error_indicators;
  }
}
}  // namespace palace
