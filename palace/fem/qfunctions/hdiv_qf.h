// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HDIV_QF_H
#define PALACE_LIBCEED_HDIV_QF_H

// libCEED QFunctions for H(div) operators (Piola transformation u = J / det(J) ̂u).
// Note: J / det(J) = adj(adj(J)^T / det(J))^T
// in[0] is geometry quadrature data, shape [ncomp=2+space_dim*dim, Q]
// in[1] is active vector, shape [qcomp=dim, ncomp=1, Q]
// out[0] is active vector, shape [qcomp=dim, ncomp=1, Q]

// Build functions assemble the quadrature point data.

#include "21/hdiv_21_qf.h"
#include "21/hdiv_build_21_qf.h"
#include "22/hdiv_22_qf.h"
#include "22/hdiv_build_22_qf.h"
#include "32/hdiv_32_qf.h"
#include "32/hdiv_build_32_qf.h"
#include "33/hdiv_33_qf.h"
#include "33/hdiv_build_33_qf.h"

#endif  // PALACE_LIBCEED_HDIV_QF_H
