// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_L2H1_QF_H
#define PALACE_LIBCEED_L2H1_QF_H

// libCEED QFunctions for mixed L2 (INTEGRAL map, u = 1 / det(J) ̂u) trial and H1 (VALUE
// map) test operators.
// in[0] is geometry quadrature data, shape [ncomp=2+space_dim*dim, Q]
// in[1] is quadrature weights, shape [Q]
// in[2] is active vector, shape [ncomp=1, Q]
// out[0] is active vector, shape [ncomp=1, Q]

// Build functions assemble the quadrature point data.

#include "1/l2h1_1_qf.h"
#include "1/l2h1_build_1_qf.h"

#endif  // PALACE_LIBCEED_L2H1_QF_H
