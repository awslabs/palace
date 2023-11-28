// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_COEFFICIENT_HPP
#define PALACE_LIBCEED_COEFFICIENT_HPP

#include <mfem.hpp>
#include "fem/libceed/ceed.hpp"

#include "fem/qfunctions/types_qf.h"

namespace palace::ceed
{

// XX TODO WIP

inline MatCoeffContext1 PopulateCoefficientContext1()
{
  MatCoeffContext1 ctx = {0};
  ctx.mat_coeff[0][0] = 1.0;
  return ctx;
}

inline MatCoeffContext2 PopulateCoefficientContext2()
{
  MatCoeffContext2 ctx = {0};
  ctx.mat_coeff[0][0] = 1.0;
  ctx.mat_coeff[0][1] = 0.0;
  ctx.mat_coeff[0][2] = 1.0;
  return ctx;
}

inline MatCoeffContext3 PopulateCoefficientContext3()
{
  MatCoeffContext3 ctx = {0};
  ctx.mat_coeff[0][0] = 1.0;
  ctx.mat_coeff[0][1] = 0.0;
  ctx.mat_coeff[0][2] = 0.0;
  ctx.mat_coeff[0][3] = 1.0;
  ctx.mat_coeff[0][4] = 0.0;
  ctx.mat_coeff[0][5] = 1.0;
  return ctx;
}

}  // namespace palace::ceed

#endif  // PALACE_LIBCEED_COEFFICIENT_HPP
