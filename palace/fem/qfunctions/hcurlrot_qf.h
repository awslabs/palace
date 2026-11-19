// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_ROT_QF_H
#define PALACE_LIBCEED_HCURL_ROT_QF_H

// libCEED QFunctions for H(curl) operators on 2D manifolds in 3D acting on the tangential
// vector rotated by 90° about the element normal (n x u = J / det(J) R ̂u for the Piola
// transformation u = adj(J)^T / det(J) ̂u and the reference rotation R).
// in[0] is geometry quadrature data, shape [ncomp=2+space_dim*dim, Q]
// in[1] is active vector, shape [qcomp=dim, ncomp=1, Q]
// out[0] is active vector, shape [qcomp=dim, ncomp=1, Q]

// Build functions assemble the quadrature point data.

#include "32/hcurlrot_32_qf.h"
#include "32/hcurlrot_build_32_qf.h"

#endif  // PALACE_LIBCEED_HCURL_ROT_QF_H
