// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_L2_QF_H
#define PALACE_LIBCEED_HCURL_L2_QF_H

// libCEED QFunctions for mixed H(curl)-L2 operators coupling an H(curl) function on a 2D
// manifold in 3D with its curl (Piola transformations u = adj(J)^T / det(J) ̂u and
// u = 1 / det(J) ̂u). The scalar curl is used as a vector along the element normal.
// in[0] is geometry quadrature data, shape [ncomp=2+space_dim*dim, Q]
// in[1] is active vector, shape [qcomp=dim or 1, ncomp=1, Q]
// out[0] is active vector, shape [qcomp=1 or dim, ncomp=1, Q]

#include "32/hcurll2_32_qf.h"

#endif  // PALACE_LIBCEED_HCURL_L2_QF_H
