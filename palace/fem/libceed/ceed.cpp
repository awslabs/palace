// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "ceed.hpp"

#include "fem/libceed/utils.hpp"

namespace palace::ceed
{

CeedGeomFactorData_private::~CeedGeomFactorData_private()
{
  PalaceCeedCall(ceed, CeedVectorDestroy(&wdetJ_vec));
  PalaceCeedCall(ceed, CeedVectorDestroy(&adjJt_vec));
  PalaceCeedCall(ceed, CeedVectorDestroy(&J_vec));
  PalaceCeedCall(ceed, CeedVectorDestroy(&attr_vec));
  PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&wdetJ_restr));
  PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&adjJt_restr));
  PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&J_restr));
  PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&attr_restr));
}

}  // namespace palace::ceed
