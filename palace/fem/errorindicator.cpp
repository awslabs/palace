// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "errorindicator.hpp"

#include <mfem/general/forall.hpp>

namespace palace
{

void ErrorIndicator::AddIndicator(const Vector &indicator, double weight)
{
  if (n == 0.0)
  {
    local = indicator;
    n = weight;
    return;
  }

  // The average local indicator is used rather than the indicator for the maximum
  // error to drive the adaptation, to account for a local error that might be marginally
  // important to many solves, rather than only large in one solve.
  MFEM_ASSERT(local.Size() == indicator.Size(),
              "Unexpected size mismatch for ErrorIndicator::AddIndicator!");

  // The local indicators must be squared before combining, so that the global error
  // calculation is valid:
  //                            E = √(1/N ∑ₙ ∑ₖ ηₖₙ²)
  // from which it follows that:
  //                            E² = 1/N ∑ₙ ∑ₖ ηₖₙ²
  //                               = 1/N ∑ₙ Eₙ²
  // Namely the average of the global error indicators included in the reduction.
  // Squaring both sides means the summation can be rearranged, and then the local error
  // indicators become:
  //                            eₖ = √(1/N ∑ₙ ηₖₙ²)
  // An indicator of weight w counts as w samples (N is then the sum of the weights).
  const bool use_dev = local.UseDevice() || indicator.UseDevice();
  const int N = local.Size();
  const double Dn = n, Dw = weight;
  const auto *DI = indicator.Read();
  auto *DL = local.ReadWrite();
  mfem::forall_switch(
      use_dev, N, [=] MFEM_HOST_DEVICE(int i)
      { DL[i] = std::sqrt((DL[i] * DL[i] * Dn + Dw * DI[i] * DI[i]) / (Dn + Dw)); });

  // More samples have been added, update for the running average.
  n += weight;
}

}  // namespace palace
