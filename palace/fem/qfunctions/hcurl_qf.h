// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_QF_H
#define PALACE_LIBCEED_HCURL_QF_H

// libCEED QFunctions for H(curl) operators (Piola transformation u = adj(J)^T / det(J) ̂u).
// in[0] is geometry quadrature data, shape [ncomp=2+space_dim*dim, Q]
// in[1] is active vector, shape [qcomp=dim, ncomp=1, Q]
// out[0] is active vector, shape [qcomp=dim, ncomp=1, Q]

// Build functions assemble the quadrature point data.

#include "21/hcurl_21_qf.h"
#include "21/hcurl_build_21_qf.h"
#include "22/hcurl_22_qf.h"
#include "22/hcurl_build_22_qf.h"
#include "32/hcurl_32_qf.h"
#include "32/hcurl_build_32_qf.h"
#include "33/hcurl_33_qf.h"
#include "33/hcurl_build_33_qf.h"

#endif  // PALACE_LIBCEED_HCURL_QF_H
