// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_HDIV_ERROR_QF_H
#define PALACE_LIBCEED_HCURL_HDIV_ERROR_QF_H

// libCEED QFunctions for computing errors between two functions, one in H(curl) and one in
// H(div) (Piola transformations u = adj(J)^T / det(J) ̂u and u = J / det(J) ̂u).
// Note: J / det(J) = adj(adj(J)^T / det(J))^T
// in[0] is geometry quadrature data, shape [ncomp=2+space_dim*dim, Q]
// in[1] is active vector 1, shape [qcomp=dim, ncomp=1, Q]
// in[2] is active vector 2, shape [qcomp=dim, ncomp=1, Q]
// out[0] is active vector, shape [ncomp=1, Q]

// Only for the square Jacobian case where dim = space_dim.

#include "22/hcurlhdiv_error_22_qf.h"
#include "33/hcurlhdiv_error_33_qf.h"

#endif  // PALACE_LIBCEED_HCURL_HDIV_ERROR_QF_H
