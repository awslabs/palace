// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "errorindicator.hpp"

#include <mfem/general/forall.hpp>

namespace palace
{

void ErrorIndicator::AddIndicator(const ErrorIndicator &indicator)
{
  if (n == 0)
  {
    local = indicator.local;
    normalization = indicator.normalization;
    n = indicator.n;
    return;
  }

  // The average local indicator is used rather than the indicator for the maximum
  // error to drive the adaptation, to account for a local error that might be marginally
  // important to many solves, rather than only large in one solve.
  MFEM_ASSERT(local.Size() == indicator.local.Size(),
              "Unexpected size mismatch for ErrorIndicator::AddIndicator!");

  // The local indicators must be squared before combining, so that the global error
  // calculation is valid:
  //                            E = √(1/N ∑ₙ ∑ₖ ηₖₙ²)
  // from which it follows that:
  //                            E² = 1/N ∑ₙ ∑ₖ ηₖₙ²
  //                               = 1/N ∑ₙ Eₙ
  // Namely the average of the global error indicators included in the reduction.
  // Squaring both sides means the summation can be rearranged, and then the local error
  // indicators become:
  //                            eₖ = √(1/N ∑ₙ ηₖₙ²)
  const int N = local.Size();
  const auto *DIL = indicator.local.Read();
  auto *DL = local.ReadWrite();
  mfem::forall(N,
               [=] MFEM_HOST_DEVICE(int i)
               {
                 DL[i] = std::sqrt((DL[i] * DL[i] * n + DIL[i] * DIL[i] * indicator.n) /
                                   (n + indicator.n));
               });

  // More samples have been added, update for the running average.
  normalization =
      (normalization * n + indicator.normalization * indicator.n) / (n + indicator.n);
  n += indicator.n;
}

}  // namespace palace
