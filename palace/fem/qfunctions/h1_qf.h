// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_H1_QF_H
#define PALACE_LIBCEED_H1_QF_H

// libCEED QFunctions for H1 operators (Piola transformation u = ̂u).
// in[0] is geometry quadrature data, shape [ncomp=2+space_dim*dim, Q]
// in[1] is active vector, shape [ncomp=vdim, Q]
// out[0] is active vector, shape [ncomp=vdim, Q]

// Build functions assemble the quadrature point data.

#include "1/h1_1_qf.h"
#include "1/h1_build_1_qf.h"
#include "2/h1_2_qf.h"
#include "2/h1_build_2_qf.h"
#include "3/h1_3_qf.h"
#include "3/h1_build_3_qf.h"

#endif  // PALACE_LIBCEED_H1_QF_H
