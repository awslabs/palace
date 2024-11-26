// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_HDIV_QF_H
#define PALACE_LIBCEED_HCURL_HDIV_QF_H

// libCEED QFunctions for mixed H(curl)-H(div) operators (Piola transformations u =
// adj(J)^T / det(J) ̂u and u = J / det(J) ̂u).
// Note: J / det(J) = adj(adj(J)^T / det(J))^T
// in[0] is geometry quadrature data, shape [ncomp=2+space_dim*dim, Q]
// in[1] is active vector, shape [qcomp=dim, ncomp=1, Q]
// out[0] is active vector, shape [qcomp=dim, ncomp=1, Q]

// Build functions assemble the quadrature point data.

#include "21/hcurlhdiv_21_qf.h"
#include "21/hcurlhdiv_build_21_qf.h"
#include "22/hcurlhdiv_22_qf.h"
#include "22/hcurlhdiv_build_22_qf.h"
#include "32/hcurlhdiv_32_qf.h"
#include "32/hcurlhdiv_build_32_qf.h"
#include "33/hcurlhdiv_33_qf.h"
#include "33/hcurlhdiv_build_33_qf.h"

#endif  // PALACE_LIBCEED_HCURL_HDIV_QF_H
