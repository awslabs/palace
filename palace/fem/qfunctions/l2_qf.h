// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_L2_QF_H
#define PALACE_LIBCEED_L2_QF_H

// libCEED QFunctions for L2 operators (Piola transformation u = 1 / det(J) ̂u).
// in[0] is geometry quadrature data, shape [ncomp=2+space_dim*dim, Q]
// in[1] is quadrature weights, shape [Q]
// in[2] is active vector, shape [ncomp=vdim, Q]
// out[0] is active vector, shape [ncomp=vdim, Q]

// Build functions assemble the quadrature point data.

#include "1/l2_1_qf.h"
#include "1/l2_build_1_qf.h"
#include "2/l2_2_qf.h"
#include "2/l2_build_2_qf.h"
#include "3/l2_3_qf.h"
#include "3/l2_build_3_qf.h"

#endif  // PALACE_LIBCEED_L2_QF_H
