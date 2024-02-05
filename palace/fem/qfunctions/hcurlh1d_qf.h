// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_H1D_QF_H
#define PALACE_LIBCEED_HCURL_H1D_QF_H

// libCEED QFunctions for mixed H(curl)-(H1)ᵈ operators (Piola transformation u =
// adj(J)^T / det(J) ̂u and u = ̂u)
// in[0] is geometry quadrature data, shape [ncomp=2+space_dim*dim, Q]
// in[1] is active vector, shape [qcomp=dim, ncomp=1, Q]
// out[0] is active vector, shape [ncomp=space_dim, Q]

// Build functions assemble the quadrature point data.

#include "21/hcurlh1d_21_qf.h"
#include "21/hcurlh1d_build_21_qf.h"
#include "22/hcurlh1d_22_qf.h"
#include "22/hcurlh1d_build_22_qf.h"
#include "32/hcurlh1d_32_qf.h"
#include "32/hcurlh1d_build_32_qf.h"
#include "33/hcurlh1d_33_qf.h"
#include "33/hcurlh1d_build_33_qf.h"

#endif  // PALACE_LIBCEED_HCURL_H1D_QF_H
